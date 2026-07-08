"""Held-out validation of `DBRetina enrich` on the KEGG pathway compendium.

Question: told only the names of specific cancers, does guilt-by-association
recover the cancer-pathway class from gene overlap alone? We leave each cancer
pathway out, score it against the other fourteen, and measure how well the
enrichment score separates cancer from non-cancer pathways (AUC), against a
size-matched random-module null.

Run from this directory after run.sh has built kegg_DBRetina_pairwise/.
"""
import os
import re
import glob
import math
import random
import statistics as st
from collections import defaultdict

import pyarrow.parquet as pq
from scipy.stats import hypergeom

random.seed(0)
PW = "kegg_DBRetina_pairwise"
CUTOFF = 20.0          # ochiai % for a neighbor edge (sparse, informative neighborhoods)
N_PERM = 2000          # random modules for the null

# names <-> ids
nt = pq.read_table(os.path.join(PW, "names.parquet")).to_pydict()
id2name = dict(zip(nt["group_id"], nt["group_name"]))
all_ids = list(id2name)
N = len(all_ids)

# neighbor graph at ochiai >= CUTOFF
neighbors = defaultdict(set)
for f in glob.glob(os.path.join(PW, "data", "*.parquet")):
    t = pq.read_table(f, columns=["group_1_id", "group_2_id", "ochiai"]).to_pydict()
    for a, b, o in zip(t["group_1_id"], t["group_2_id"], t["ochiai"]):
        if o >= CUTOFF:
            neighbors[a].add(b)
            neighbors[b].add(a)
deg = {g: len(neighbors[g]) for g in all_ids}
active = [g for g in all_ids if deg[g] > 0]

# objective cancer module: KEGG cancer-type pathways, defined by name only
cancer = {g for g in all_ids if re.search(r"cancer|carcinoma|leukemia|melanoma|glioma", id2name[g])}
M = len(cancer)


def score(g, module):
    """enrich's score: -log10 hypergeometric sf of module hits among g's neighbors."""
    d = deg[g]
    if d == 0:
        return 0.0
    hits = len(neighbors[g] & module)
    if hits == 0:
        return 0.0
    p = float(hypergeom.sf(hits - 1, N, len(module), d))
    return -math.log10(max(p, 1e-300))


def recovery_auc(members):
    """AUC of (is-member) vs enrichment score, leave-one-out for the members."""
    pos, neg = [], []
    for g in active:
        mod = (members - {g}) if g in members else members
        (pos if g in members else neg).append(score(g, mod))
    if not pos or not neg:
        return 0.5
    wins = ties = 0
    for a in pos:
        for b in neg:
            if a > b:
                wins += 1
            elif a == b:
                ties += 1
    return (wins + 0.5 * ties) / (len(pos) * len(neg))


real = recovery_auc(cancer)

# per-pathway leave-one-out scores, for the "where do cancer pathways rank" figure
with open("scores.tsv", "w") as fh:
    fh.write("pathway\tis_cancer\tscore\n")
    for g in active:
        mod = (cancer - {g}) if g in cancer else cancer
        fh.write(f"{id2name[g]}\t{int(g in cancer)}\t{score(g, mod):.4f}\n")

null = sorted(recovery_auc(set(random.sample(all_ids, M))) for _ in range(N_PERM))
mu, sd = st.mean(null), st.pstdev(null)
p = (sum(a >= real for a in null) + 1) / (len(null) + 1)
z = (real - mu) / sd if sd else float("inf")

print(f"KEGG pathways: {N} | edges at ochiai>={CUTOFF:.0f}% | active groups: {len(active)}")
print(f"cancer module (by name): {M} pathways")
print("\nHELD-OUT RECOVERY (leave-one-out)")
print(f"  real AUC          = {real:.3f}")
print(f"  null AUC ({N_PERM}x)  = {mu:.3f} +/- {sd:.3f}   (95th pct {null[int(0.95 * len(null))]:.3f})")
print(f"  {z:.1f} SD above null   permutation p = {p:.2e}"
      + ("   (0 of %d random modules reached it)" % N_PERM if p <= 1.0 / N_PERM else ""))

with open("validation.txt", "w") as fh:
    fh.write(f"substrate\t{N} KEGG pathways\n")
    fh.write(f"cutoff\tochiai >= {CUTOFF:.0f}%\n")
    fh.write(f"cancer_module\t{M} pathways\n")
    fh.write(f"real_auc\t{real:.4f}\n")
    fh.write(f"null_auc_mean\t{mu:.4f}\n")
    fh.write(f"null_auc_sd\t{sd:.4f}\n")
    fh.write(f"z\t{z:.2f}\n")
    fh.write(f"perm_p\t{p:.2e}\n")
    fh.write(f"n_perm\t{N_PERM}\n")

# stash the null distribution for the figure
with open("null_aucs.txt", "w") as fh:
    fh.write("\n".join(f"{a:.5f}" for a in null) + "\n")
    fh.write(f"#real\t{real:.5f}\n")
