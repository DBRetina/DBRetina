import json
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import networkx as nx
from matplotlib.patches import Patch

edges_f, labels_f, out_f, cutoff, title = sys.argv[1], sys.argv[2], sys.argv[3], float(sys.argv[4]), sys.argv[5]
distmode = sys.argv[6] if len(sys.argv) > 6 else "linear"   # "linear" (0-100 metric) or "inv" (odds ratio)
cat = {k[4:]: v for k, v in json.load(open(labels_f)).items()}
CATCOL = {"autism": "#c0392b", "neurodev": "#2980b9", "psychiatric": "#8e44ad", "syndromic": "#e67e22",
          "control": "#7f8c8d", "autoimmune": "#c0392b", "comorbidity": "#16a085",
          "alzheimer": "#c0392b", "neurodegenerative": "#2980b9", "vascular": "#8e44ad", "metabolic": "#e67e22",
          "type_2_diabetes": "#8e44ad", "type_1_diabetes": "#2980b9", "complication": "#f1c40f"}
CATNAME = {"autism": "autism (seed)", "neurodev": "neurodevelopmental", "psychiatric": "psychiatric",
           "syndromic": "monogenic syndrome", "control": "non-neuro control",
           "autoimmune": "autoimmune", "comorbidity": "comorbidity",
           "alzheimer": "Alzheimer's (seed)", "neurodegenerative": "neurodegenerative", "vascular": "vascular",
           "metabolic": "metabolic", "type_2_diabetes": "type 2 diabetes (seed)",
           "type_1_diabetes": "type 1 diabetes", "complication": "diabetes complication"}
LABEL = {"asd": "autism", "intellectual_disability": "intellectual disability", "adhd": "ADHD",
         "dev_delay": "developmental delay", "schizophrenia": "schizophrenia", "bipolar": "bipolar",
         "depression": "depression", "anxiety": "anxiety", "ocd": "OCD", "tourette": "Tourette",
         "rett": "Rett", "fragile_x": "Fragile X", "tuberous_sclerosis": "tuberous sclerosis",
         "down_syndrome": "Down syndrome", "breast_cancer": "breast cancer", "ra": "rheumatoid arthritis",
         "cad": "coronary artery disease", "epilepsy": "epilepsy", "lupus": "lupus", "ms": "MS", "t1d": "type 1 diabetes",
         "uc": "ulcerative colitis", "crohn": "Crohn's", "psoriasis": "psoriasis", "celiac": "celiac",
         "atherosclerosis": "atherosclerosis", "coronary_artery_disease": "coronary artery disease",
         "myocardial_infarction": "myocardial infarction", "major_depression": "major depression", "parkinson": "Parkinson's",
         "alzheimer": "Alzheimer's", "t2d": "type 2 diabetes", "als": "ALS", "huntington": "Huntington's",
         "dementia": "dementia", "ftd": "frontotemporal dementia", "lewy_body": "Lewy body", "tauopathies": "tauopathies",
         "prion": "prion disease", "ischemic_stroke": "ischemic stroke", "cerebral_infarction": "cerebral infarction",
         "vascular_dementia": "vascular dementia", "obesity": "obesity", "hypertension": "hypertension",
         "metabolic_syndrome": "metabolic syndrome", "nafld": "fatty liver (NAFLD)", "dyslipidemia": "dyslipidemia",
         "mi": "myocardial infarction", "diabetic_nephropathy": "diabetic nephropathy",
         "diabetic_retinopathy": "diabetic retinopathy"}

alled = {}
for i, line in enumerate(open(edges_f)):
    if i == 0:
        continue
    a, b, w = line.rstrip("\n").split("\t")
    alled[(a[4:], b[4:])] = float(w)

G = nx.Graph()
G.add_nodes_from(cat)
for (a, b), w in alled.items():
    if w >= cutoff:
        G.add_edge(a, b, weight=w)
# no floating nodes: give any disconnected node its single strongest edge
strongest = {}
for (a, b), w in alled.items():
    for x, y in ((a, b), (b, a)):
        if w > strongest.get(x, (0, None))[0]:
            strongest[x] = (w, y)
for n in G.nodes():
    if G.degree(n) == 0 and n in strongest:
        w, m = strongest[n]
        G.add_edge(n, m, weight=w)

for u, v in G.edges():
    w = G[u][v]["weight"]
    G[u][v]["dist"] = (1.0 / (w + 1e-6)) if distmode == "inv" else (100.0 - w)   # similar -> close
pos = nx.kamada_kawai_layout(G, weight="dist")

# de-overlap: push apart any nodes closer than min_d, so labels never collide
import numpy as np
rng = np.random.default_rng(0)
nodes = list(G.nodes())
P = np.array([pos[n] for n in nodes], dtype=float)
min_d = 0.34
for _ in range(300):
    moved = False
    for i in range(len(nodes)):
        for j in range(i + 1, len(nodes)):
            d = P[i] - P[j]
            dist = float(np.hypot(*d))
            if dist < min_d:
                if dist < 1e-9:
                    d = rng.standard_normal(2) * 1e-3
                    dist = float(np.hypot(*d))
                push = (min_d - dist) / 2.0
                u = d / dist
                P[i] += u * push
                P[j] -= u * push
                moved = True
    if not moved:
        break
pos = {n: P[k] for k, n in enumerate(nodes)}

fig, ax = plt.subplots(figsize=(12, 9.5))
ws = [G[u][v]["weight"] for u, v in G.edges()]
wmin, wmax = min(ws), max(ws)
for u, v in G.edges():
    f = (G[u][v]["weight"] - wmin) / (wmax - wmin + 1e-9)
    ax.plot([pos[u][0], pos[v][0]], [pos[u][1], pos[v][1]], "-", color="#9aa0a6",
            lw=0.8 + 4.5 * f, alpha=0.25 + 0.55 * f, zorder=1, solid_capstyle="round")
for d in G.nodes():
    ax.scatter(*pos[d], s=460, color=CATCOL.get(cat[d], "#777"), edgecolors="white", linewidths=2.0, zorder=2)
    ax.annotate(LABEL.get(d, d), (pos[d][0], pos[d][1] - 0.055), ha="center", va="top", fontsize=10,
                fontweight="bold", color="#2c3e50", zorder=3,
                path_effects=[pe.withStroke(linewidth=3, foreground="white")])
cats = list(dict.fromkeys(cat.values()))
ax.legend(handles=[Patch(facecolor=CATCOL[c], label=CATNAME.get(c, c)) for c in cats],
          fontsize=10.5, loc="upper left", framealpha=.96)
ax.set_title(title, fontsize=14, fontweight="bold")
ax.margins(0.16)
ax.axis("off")
plt.tight_layout()
plt.savefig(out_f, dpi=140, bbox_inches="tight")
print(f"wrote {out_f}: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges, "
      f"{nx.number_connected_components(G)} components at cutoff {cutoff}")
