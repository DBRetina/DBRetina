#!/usr/bin/env python3
"""Two benchmarks that validate DBRetina's similarity metric.

Reads a combined GMT with pathways (names prefixed kegg_/reactome_/wp_) and diseases (dis_).
Uses ../integrative-network/data/multidb.gmt by default (set GMT to override).

Test 1: cross-database pathway matching. For pathways carrying the same name across two
databases, is the correct counterpart DBRetina's top match? (precision@k)
Test 2: gene prediction. Does network-weighted voting beat a gene-frequency baseline
(leave-one-disease-out AUC)?
Writes figures/pathway_matching_proof.png.
"""
import os
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy.stats import wilcoxon

GMT = os.environ.get("GMT", "data/diseases_and_pathways.gmt")
STOP = set("signaling signalling pathway pathways of the in and to a by response system "
           "disease regulation via mediated".split())

sets, db, dis = {}, {}, {}
for line in open(GMT):
    p = line.rstrip("\n").split("\t")
    n, genes = p[0], set(g.upper() for g in p[2:] if g.strip())
    for pre in ("kegg_", "reactome_", "wp_"):
        if n.startswith(pre):
            sets[n] = genes
            db[n] = pre.rstrip("_")
    if n.startswith("dis_"):
        dis[n[4:]] = genes


def ochiai(a, b):
    o = len(a & b)
    return o / np.sqrt(len(a) * len(b)) if o else 0.0


# ---- Test 1: cross-database pathway matching ----
paths = list(sets)
tok = lambda n: frozenset(t for t in re.sub(r"^(kegg|reactome|wp)_", "", n).split("_")
                          if t not in STOP and len(t) > 2)
T = {n: tok(n) for n in paths}
gt = [(a, b) for i, a in enumerate(paths) for b in paths[i + 1:]
      if db[a] != db[b] and T[a] and T[a] == T[b]]
ranks, same = [], []
for a, b in gt:
    for x, y in ((a, b), (b, a)):
        cands = sorted(((ochiai(sets[x], sets[c]), c) for c in paths if db[c] == db[y] and c != x),
                       reverse=True)
        ranks.append([c for _, c in cands].index(y) + 1)
        same.append(ochiai(sets[x], sets[y]) * 100)
ranks = np.array(ranks)
print(f"cross-DB pathway matching: {len(gt)} same-process pairs")
print(f"  precision@1 = {np.mean(ranks == 1) * 100:.0f}%   precision@5 = {np.mean(ranks <= 5) * 100:.0f}%")

# ---- Test 2: gene prediction vs frequency baseline ----
diseases = list(dis)
sim = {(a, b): ochiai(dis[a], dis[b]) for a in diseases for b in diseases if a != b}


def auc(y, s):
    y, s = np.array(y), np.array(s)
    pos, neg = s[y == 1], s[y == 0]
    return (np.sum(pos[:, None] > neg[None, :]) + 0.5 * np.sum(pos[:, None] == neg[None, :])) / (len(pos) * len(neg))


gains = []
for D in diseases:
    others = [d for d in diseases if d != D]
    universe = set().union(*[dis[d] for d in others])
    net = {g: 0.0 for g in universe}
    freq = {g: 0 for g in universe}
    for d in others:
        for g in dis[d]:
            net[g] += sim[(D, d)]
            freq[g] += 1
    genes = list(universe)
    y = [1 if g in dis[D] else 0 for g in genes]
    gains.append(auc(y, [net[g] for g in genes]) - auc(y, [freq[g] for g in genes]))
gains = np.array(gains)
_, p = wilcoxon(gains)
print(f"gene prediction: network beats frequency in {np.sum(gains > 0)}/{len(gains)} diseases, "
      f"mean gain {gains.mean():+.3f} AUC, Wilcoxon p = {p:.1e}")

# ---- figure for test 1 ----
rng = np.random.default_rng(0)
rand = [ochiai(sets[paths[i]], sets[paths[j]]) * 100
        for i, j in ((rng.integers(len(paths)), rng.integers(len(paths))) for _ in range(4000))
        if db[paths[i]] != db[paths[j]]]
sns.set_theme(style="whitegrid")
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))
ks = range(1, 11)
ax1.plot(list(ks), [np.mean(ranks <= k) * 100 for k in ks], "o-", color="#2e86c1", lw=2.5, ms=7)
for k in (1, 5):
    ax1.annotate(f"{np.mean(ranks <= k) * 100:.0f}%", (k, np.mean(ranks <= k) * 100),
                 textcoords="offset points", xytext=(0, 10), ha="center", fontweight="bold")
ax1.set_xlabel("top-k cross-database matches considered")
ax1.set_ylabel("% of same-process pathways found (precision@k)")
ax1.set_ylim(0, 100)
ax1.set_title(f"DBRetina matches the same pathway\nacross databases ({len(gt)} name-matched pairs)",
              fontsize=12, fontweight="bold")
ax2.hist(rand, bins=30, color="#95a5a6", alpha=.8, label="random cross-DB pairs", density=True)
ax2.hist(same, bins=30, color="#2e86c1", alpha=.7, label="same-process pairs", density=True)
ax2.set_xlabel("ochiai similarity (%)")
ax2.set_ylabel("density")
ax2.legend()
ax2.set_title(f"Same-process pairs: median {np.median(same):.0f}%\nrandom pairs: median {np.median(rand):.0f}%",
              fontsize=12, fontweight="bold")
plt.tight_layout()
os.makedirs("figures", exist_ok=True)
plt.savefig("figures/pathway_matching_proof.png", dpi=140, bbox_inches="tight")
print("wrote figures/pathway_matching_proof.png")
