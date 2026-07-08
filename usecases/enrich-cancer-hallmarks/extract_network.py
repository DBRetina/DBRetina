"""Extract the cancer-enrichment subgraph for the interactive tutorial figure.

Nodes: the 15 cancer pathways + the recovered oncogenic-signaling pathways.
Edges: their Ochiai similarities, thresholded so the graph shows structure, not a
hairball. Also records, per node, how many of its FULL-graph neighbours are cancer
pathways (the guilt-by-association evidence the holdout demo shows). Writes
network.json, which the EnrichNetwork component embeds.
"""
import glob
import json
import re
from collections import defaultdict

import networkx as nx
import pyarrow.parquet as pq

PW = "kegg_DBRetina_pairwise"
NBR_CUT = 20.0     # the neighbourhood cutoff enrich uses
GRAPH_CUT = 28.0   # only draw edges this strong, so the layout reads clearly

nt = pq.read_table(f"{PW}/names.parquet").to_pydict()
id2name = dict(zip(nt["group_id"], nt["group_name"]))
name2id = {v: k for k, v in id2name.items()}
cancer = {g for g in id2name if re.search(r"cancer|carcinoma|leukemia|melanoma|glioma", id2name[g])}

# recovered nodes: the clearly-oncogenic top hits (immune-receptor ones left out to keep it clean)
ONCO = ["focal_adhesion", "erbb_signaling_pathway", "mtor_signaling_pathway", "apoptosis",
        "vegf_signaling_pathway", "insulin_signaling_pathway", "neurotrophin_signaling_pathway",
        "regulation_of_actin_cytoskeleton"]
recovered = [name2id["kegg_" + k] for k in ONCO if "kegg_" + k in name2id]
sel = set(cancer) | set(recovered)

full_nbr = defaultdict(set)
cand = defaultdict(list)
for f in glob.glob(f"{PW}/data/*.parquet"):
    t = pq.read_table(f, columns=["group_1_id", "group_2_id", "ochiai"]).to_pydict()
    for a, b, o in zip(t["group_1_id"], t["group_2_id"], t["ochiai"]):
        if o >= NBR_CUT:
            full_nbr[a].add(b); full_nbr[b].add(a)
        if a in sel and b in sel:
            cand[a].append((o, b)); cand[b].append((o, a))
# draw only each node's K strongest links -> a readable backbone, not a hairball.
# (The "surrounded by cancer" evidence uses the FULL neighbourhood, below, not this.)
K = 2
edge_map = {}
for g in sel:
    for o, h in sorted(cand[g], reverse=True)[:K]:
        edge_map[(min(g, h), max(g, h))] = o
sub_edges = [(a, b, o) for (a, b), o in edge_map.items()]

G = nx.Graph()
G.add_nodes_from(sel)
for a, b, o in sub_edges:
    G.add_edge(a, b, weight=o / 8.0)
pos = nx.spring_layout(G, seed=3, k=1.15, iterations=500)
xs = [p[0] for p in pos.values()]
ys = [p[1] for p in pos.values()]
sx = lambda x: round(12 + (x - min(xs)) / (max(xs) - min(xs)) * 216, 1)
sy = lambda y: round(10 + (y - min(ys)) / (max(ys) - min(ys)) * 100, 1)

ABBR = {"non small cell lung cancer": "NSCLC", "small cell lung cancer": "SCLC",
        "acute myeloid leukemia": "AML", "chronic myeloid leukemia": "CML",
        "renal cell carcinoma": "renal cell carcinoma", "basal cell carcinoma": "basal cell carcinoma",
        "pathways in cancer": "pan-cancer", "endometrial cancer": "endometrial",
        "colorectal cancer": "colorectal", "pancreatic cancer": "pancreatic",
        "prostate cancer": "prostate", "thyroid cancer": "thyroid", "bladder cancer": "bladder"}


def short(nm):
    s = nm.replace("kegg_", "").replace("_", " ")
    if s in ABBR:
        return ABBR[s]
    s = (s.replace("signaling pathway", "signaling").replace("regulation of ", "")
           .replace(" cancer", "").replace(" mediated", ""))
    return s


# held-out (leave-one-out) rank of each pathway, from validate.py's scores.tsv,
# for the holdout caption ("enrich still ranks it #N")
rank = {}
try:
    srows = []
    with open("scores.tsv") as fh:
        next(fh)
        for line in fh:
            nm, _isc, sc = line.rstrip("\n").split("\t")
            srows.append((float(sc), nm))
    srows.sort(reverse=True)
    rank = {nm: i + 1 for i, (_sc, nm) in enumerate(srows)}
    n_ranked = len(srows)
except FileNotFoundError:
    n_ranked = 0

nodes = []
for g in sel:
    is_c = g in cancer
    cancer_nbrs = len(full_nbr[g] & cancer) - (1 if is_c else 0)
    nodes.append({"id": str(g), "label": short(id2name[g]),
                  "color": "#c0392b" if is_c else "#16a085",  # cancer red, recovered teal
                  "hold": 1 if is_c else 0,                   # the button cycles cancer pathways
                  "x": sx(pos[g][0]), "y": sy(pos[g][1]),
                  "cancerNbrs": cancer_nbrs, "deg": len(full_nbr[g]),
                  "looRank": rank.get(id2name[g], 0)})
edges = [{"s": str(a), "t": str(b), "w": round(o, 1)} for a, b, o in sub_edges]
json.dump({"nodes": nodes, "edges": edges, "totalRanked": n_ranked}, open("network.json", "w"), indent=1)
print(f"nodes {len(nodes)} ({len(cancer)} cancer, {len(recovered)} recovered) | "
      f"edges {len(edges)} (top-2 backbone, {2*len(edges)/len(nodes):.1f} per node)")
