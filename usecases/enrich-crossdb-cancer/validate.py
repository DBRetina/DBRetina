"""Cross-database replication test for `DBRetina enrich`.

Seed enrich with KEGG's cancer pathways only, then ask how well it recovers the
cancer pathways that Reactome and WikiPathways curated independently (the seed
never saw them). Reported against a null of random seeds of the same size.

Run from this directory after run.sh has built c2cp_DBRetina_pairwise/.
"""
import glob
import json
import math
import os
import random
import re
import statistics as st
from collections import defaultdict

import pyarrow.parquet as pq
from scipy.stats import hypergeom

random.seed(0)
PW = "c2cp_DBRetina_pairwise"
CUTOFF = 20.0
N_PERM = 2000

nt = pq.read_table(f"{PW}/names.parquet").to_pydict()
id2name = dict(zip(nt["group_id"], nt["group_name"]))
all_ids = list(id2name)
N = len(all_ids)
lab = json.load(open("data/labels.json"))
kegg_cancer = {g for g in all_ids if lab[id2name[g]]["db"] == "kegg" and lab[id2name[g]]["cancer"]}
ext_cancer = {g for g in all_ids if lab[id2name[g]]["db"] in ("reactome", "wp") and lab[id2name[g]]["cancer"]}

# neighbor graph at ochiai >= CUTOFF
neighbors = defaultdict(set)
for f in glob.glob(f"{PW}/data/*.parquet"):
    t = pq.read_table(f, columns=["group_1_id", "group_2_id", "ochiai"]).to_pydict()
    for a, b, o in zip(t["group_1_id"], t["group_2_id"], t["ochiai"]):
        if o >= CUTOFF:
            neighbors[a].add(b)
            neighbors[b].add(a)
deg = {g: len(neighbors[g]) for g in all_ids}
non_kegg = [g for g in all_ids if lab[id2name[g]]["db"] != "kegg" and deg[g] > 0]


def score(g, module):
    d = deg[g]
    if d == 0:
        return 0.0
    h = len(neighbors[g] & module)
    if h == 0:
        return 0.0
    return -math.log10(max(float(hypergeom.sf(h - 1, N, len(module), d)), 1e-300))


def auc_for(seed, positives):
    """Mann-Whitney AUC (average ranks, ties at 0.5): positives vs the rest, over non-KEGG."""
    scored = sorted((score(g, seed), g in positives) for g in non_kegg)
    pos_n = sum(1 for _, p in scored if p)
    neg_n = len(scored) - pos_n
    if pos_n == 0 or neg_n == 0:
        return 0.5
    ranks = [0.0] * len(scored)
    i = 0
    while i < len(scored):
        j = i
        while j < len(scored) and scored[j][0] == scored[i][0]:
            j += 1
        for k in range(i, j):
            ranks[k] = (i + 1 + j) / 2.0
        i = j
    spr = sum(ranks[k] for k in range(len(scored)) if scored[k][1])
    return (spr - pos_n * (pos_n + 1) / 2.0) / (pos_n * neg_n)


real = auc_for(kegg_cancer, ext_cancer)

# per-pathway scores (all non-KEGG pathways), for the rank figure
with open("scores.tsv", "w") as fh:
    fh.write("pathway\tis_cancer\tscore\n")
    for g in non_kegg:
        fh.write(f"{id2name[g]}\t{int(g in ext_cancer)}\t{score(g, kegg_cancer):.4f}\n")

null = sorted(auc_for(set(random.sample(all_ids, len(kegg_cancer))), ext_cancer) for _ in range(N_PERM))
mu, sd = st.mean(null), st.pstdev(null)
p = (sum(a >= real for a in null) + 1) / (len(null) + 1)
z = (real - mu) / sd if sd else float("inf")

print(f"compendium: {N} pathways from 5 databases | edges at ochiai>={CUTOFF:.0f}%")
print(f"seed: {len(kegg_cancer)} KEGG cancer pathways | held-out: {len(ext_cancer)} Reactome+WP cancer pathways")
print("\nCROSS-DATABASE RECOVERY")
print(f"  real AUC          = {real:.3f}")
print(f"  null AUC ({N_PERM}x)  = {mu:.3f} +/- {sd:.3f}")
print(f"  {z:.1f} SD above null   permutation p = {p:.2e}")

# top non-cancer-named cross-DB hits, for the literature check + the figure
rows = []
for g in non_kegg:
    if lab[id2name[g]]["cancer"]:
        continue
    h = len(neighbors[g] & kegg_cancer)
    if h == 0:
        continue
    rows.append((float(hypergeom.sf(h - 1, N, len(kegg_cancer), deg[g])), id2name[g], lab[id2name[g]]["db"]))
rows.sort()
with open("crossdb_top.tsv", "w") as fh:
    fh.write("pathway\tdatabase\tenrichment_p\n")
    for pv, nm, d in rows[:25]:
        fh.write(f"{nm}\t{d}\t{pv:.3g}\n")

# the independently-curated cancer pathways, ranked by enrichment for the KEGG seed
# (this is the recovery itself: known cancer drivers, from other databases, pulled up)
rec = []
for g in ext_cancer:
    if deg[g] == 0:
        continue
    h = len(neighbors[g] & kegg_cancer)
    pv = float(hypergeom.sf(h - 1, N, len(kegg_cancer), deg[g])) if h else 1.0
    rec.append((pv, id2name[g], lab[id2name[g]]["db"], h))
rec.sort()
with open("crossdb_recovered.tsv", "w") as fh:
    fh.write("pathway\tdatabase\thits\tenrichment_p\n")
    for pv, nm, d, h in rec:
        fh.write(f"{nm}\t{d}\t{h}\t{pv:.3g}\n")
with open("validation.txt", "w") as fh:
    for k, v in [("substrate", f"{N} C2:CP pathways"), ("cutoff", f"ochiai >= {CUTOFF:.0f}%"),
                 ("seed", f"{len(kegg_cancer)} KEGG cancer"), ("heldout", f"{len(ext_cancer)} Reactome+WP cancer"),
                 ("real_auc", f"{real:.4f}"), ("null_auc_mean", f"{mu:.4f}"), ("null_auc_sd", f"{sd:.4f}"),
                 ("z", f"{z:.2f}"), ("perm_p", f"{p:.2e}"), ("n_perm", str(N_PERM))]:
        fh.write(f"{k}\t{v}\n")
with open("null_aucs.txt", "w") as fh:
    fh.write("\n".join(f"{a:.5f}" for a in null) + f"\n#real\t{real:.5f}\n")
print("\ntop cross-database hits -> crossdb_top.tsv")
