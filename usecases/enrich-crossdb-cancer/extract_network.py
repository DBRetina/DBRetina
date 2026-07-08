"""Extract the cross-database subgraph for the interactive tutorial figure.

Nodes: the 15 KEGG cancer pathways (the seed) + the top cancer pathways that
Reactome and WikiPathways curated independently and that the KEGG seed recovers.
Coloured by the database that curated each. The holdout button cycles through the
independently-curated ones, showing each is surrounded by cancer pathways (from
every database) and so is ranked back near the top by the KEGG-only seed.
Writes network.json, which the EnrichNetwork component embeds.
"""
import glob
import json
import re
from collections import defaultdict

import networkx as nx
import pyarrow.parquet as pq

PW = "c2cp_DBRetina_pairwise"
NBR_CUT = 20.0
COLOR = {"kegg": "#c0392b", "reactome": "#2c7fb8", "wp": "#2ca25f"}

nt = pq.read_table(f"{PW}/names.parquet").to_pydict()
id2name = dict(zip(nt["group_id"], nt["group_name"]))
name2id = {v: k for k, v in id2name.items()}


def db(nm):
    m = re.match(r"(kegg|reactome|biocarta|pid|wp)_", nm)
    return m.group(1) if m else "other"


CANCER = re.compile(r"cancer|carcinoma|leukemia|melanoma|glioma")
cancer = {g for g in id2name if CANCER.search(id2name[g])}
kegg_cancer = [g for g in cancer if db(id2name[g]) == "kegg"]

# recovered independent cancer nodes: the top ones per database (by enrichment for
# the KEGG seed), from validate.py's crossdb_recovered.tsv
by_db = {"reactome": [], "wp": []}
with open("crossdb_recovered.tsv") as fh:
    next(fh)
    for line in fh:
        nm, d, hits, pv = line.rstrip("\n").split("\t")
        if d in by_db and nm in name2id:
            by_db[d].append(name2id[nm])
recovered = by_db["reactome"][:4] + by_db["wp"][:6]   # keep it as readable as the KEGG graph
sel = set(kegg_cancer) | set(recovered)

# LOO rank of each pathway for the KEGG seed, from validate.py's scores.tsv
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

full_nbr = defaultdict(set)
cand = defaultdict(list)
for f in glob.glob(f"{PW}/data/*.parquet"):
    t = pq.read_table(f, columns=["group_1_id", "group_2_id", "ochiai"]).to_pydict()
    for a, b, o in zip(t["group_1_id"], t["group_2_id"], t["ochiai"]):
        if o >= NBR_CUT:
            full_nbr[a].add(b); full_nbr[b].add(a)
        if a in sel and b in sel:
            cand[a].append((o, b)); cand[b].append((o, a))
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
pos = nx.spring_layout(G, seed=5, k=1.5, iterations=600)
xs = [p[0] for p in pos.values()]
ys = [p[1] for p in pos.values()]
sx = lambda x: round(12 + (x - min(xs)) / (max(xs) - min(xs)) * 216, 1)
sy = lambda y: round(10 + (y - min(ys)) / (max(ys) - min(ys)) * 100, 1)


NAME_MAP = {
    "reactome_aberrant_regulation_of_mitotic_g1_s_transition_in_cancer_due_to_rb1_defects": "RB1 loss (cancer)",
    "reactome_loss_of_function_of_smad2_3_in_cancer": "SMAD2/3 loss (cancer)",
    "reactome_constitutive_signaling_by_ligand_responsive_egfr_cancer_variants": "EGFR variants (cancer)",
    "reactome_signaling_by_tgf_beta_receptor_complex_in_cancer": "TGF-beta (cancer)",
    "reactome_pi3k_akt_signaling_in_cancer": "PI3K-AKT (cancer)",
    "reactome_signaling_by_egfr_in_cancer": "EGFR (cancer)",
    "wp_chromosomal_and_microsatellite_instability_in_colorectal_cancer": "colorectal (MSI)",
    "wp_epithelial_to_mesenchymal_transition_in_colorectal_cancer": "EMT colorectal",
    "wp_pancreatic_adenocarcinoma_pathway": "pancreatic",
    "wp_mirna_regulation_of_prostate_cancer_signaling_pathways": "prostate (miRNA)",
}


def short(nm):
    if nm in NAME_MAP:
        return NAME_MAP[nm]
    s = nm
    for pre in ("kegg_", "reactome_", "biocarta_", "pid_", "wp_"):
        if s.startswith(pre):
            s = s[len(pre):]
    s = s.replace("_", " ")
    s = re.sub(r"signaling by (.+?) in cancer", r"\1 (cancer)", s)
    s = (s.replace("pi3k akt signaling in cancer", "PI3K-AKT (cancer)")
           .replace("nonsmall cell lung cancer", "NSCLC")
           .replace("head and neck squamous cell carcinoma", "head & neck")
           .replace("mirna regulation of prostate cancer signaling pathways", "prostate (miRNA)")
           .replace(" cancer pathway", "").replace(" cancer", "").replace(" pathway", ""))
    s = s[:24].strip()
    for gg in ("egfr", "erbb2", "alk", "wnt", "pi3k", "akt", "tgf", "rb1", "notch1", "smad"):
        s = re.sub(rf"\b{gg}\b", gg.upper(), s)
    return s


nodes = []
for g in sel:
    d = db(id2name[g])
    cancer_nbrs = len(full_nbr[g] & cancer) - (1 if g in cancer else 0)
    nodes.append({"id": str(g), "label": short(id2name[g]), "color": COLOR[d],
                  "hold": 0 if d == "kegg" else 1,   # cycle the independently-curated ones
                  "x": sx(pos[g][0]), "y": sy(pos[g][1]),
                  "cancerNbrs": cancer_nbrs, "deg": len(full_nbr[g]),
                  "looRank": rank.get(id2name[g], 0)})
edges = [{"s": str(a), "t": str(b), "w": round(o, 1)} for a, b, o in sub_edges]
json.dump({"nodes": nodes, "edges": edges, "totalRanked": n_ranked}, open("network.json", "w"), indent=1)
print(f"nodes {len(nodes)} ({len(kegg_cancer)} KEGG, {len(recovered)} recovered from Reactome/WP) | "
      f"edges {len(edges)} (top-2 backbone)")
