#!/usr/bin/env python3
"""Show that a size-aware metric recovers Alzheimer's real neurodegeneration neighbours.

Reads data/panel.gmt and data/labels.json, ranks Alzheimer's neighbours by ochiai (size-biased) and
by odds ratio (size-aware), and writes figures/alzheimer_metric.png.
"""
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib.patches import Patch

sets = {}
for line in open("data/panel.gmt"):
    p = line.rstrip("\n").split("\t")
    sets[p[0][4:]] = set(g.upper() for g in p[2:] if g.strip())
cat = {k[4:]: v for k, v in json.load(open("data/labels.json")).items()}
U = len(set().union(*sets.values()))


def ochiai(a, b):
    o = len(sets[a] & sets[b])
    return o / np.sqrt(len(sets[a]) * len(sets[b])) * 100 if o else 0.0


def odds(a, b):
    i = len(sets[a] & sets[b])
    ab = len(sets[a]) - i
    ba = len(sets[b]) - i
    d = U - i - ab - ba
    return (i * d) / (ab * ba) if ab * ba > 0 else float("inf")


GRP = {"neurodegenerative": "neurodegenerative / dementia", "vascular": "neurodegenerative / dementia",
       "metabolic": "metabolic", "control": "non-neuro control"}
COL = {"neurodegenerative / dementia": "#2980b9", "metabolic": "#e67e22", "non-neuro control": "#7f8c8d"}
LAB = {"parkinson": "Parkinson's", "als": "ALS", "huntington": "Huntington's", "dementia": "dementia",
       "ftd": "frontotemporal dementia", "lewy_body": "Lewy body", "tauopathies": "tauopathies",
       "prion": "prion disease", "ischemic_stroke": "ischemic stroke", "cerebral_infarction": "cerebral infarction",
       "vascular_dementia": "vascular dementia", "t2d": "type 2 diabetes", "obesity": "obesity",
       "metabolic_syndrome": "metabolic syndrome", "breast_cancer": "breast cancer", "ra": "rheumatoid arthritis"}
others = [d for d in sets if d != "alzheimer"]

sns.set_theme(style="whitegrid")
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 7))


def panel(ax, scorefn, title, xlabel):
    rows = sorted(((scorefn("alzheimer", d), d) for d in others), reverse=True)
    y = list(range(len(rows)))[::-1]
    ax.barh(y, [v for v, _ in rows], color=[COL[GRP[cat[d]]] for _, d in rows], alpha=.9)
    ax.set_yticks(y)
    ax.set_yticklabels([LAB.get(d, d) for _, d in rows], fontsize=9.5)
    ax.set_xlabel(xlabel, fontsize=10.5)
    ax.set_title(title, fontsize=12.5, fontweight="bold")


panel(ax1, ochiai, "Raw overlap (ochiai):\nfooled by gene-set size", "similarity to Alzheimer's (ochiai %)")
panel(ax2, odds, "Size-aware (odds ratio):\nthe real neurodegeneration signal", "enrichment over chance (odds ratio)")
ax2.axvline(1, ls=":", color="#7f8c8d")
fig.legend(handles=[Patch(facecolor=COL[k], label=k) for k in COL], fontsize=10.5, loc="lower center",
           ncol=3, bbox_to_anchor=(0.5, -0.02))
plt.tight_layout(rect=[0, 0.04, 1, 1])
plt.savefig("figures/alzheimer_metric.png", dpi=140, bbox_inches="tight")
print("wrote figures/alzheimer_metric.png")
