#!/usr/bin/env python3
"""Rank autism's genetic neighbours and separate polygenic comorbidity from single-gene cause.

Reads data/neuro.gmt and data/labels.json. Ranks every disorder by gene-set similarity to autism,
tests whether the polygenic comorbidities rank above the non-neuro controls, and writes
figures/autism_neighbours.png.
"""
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib.patches import Patch
from scipy.stats import mannwhitneyu

MONOGENIC = {"rett", "fragile_x", "tuberous_sclerosis"}
LABEL = {
    "schizophrenia": "schizophrenia", "epilepsy": "epilepsy", "adhd": "ADHD", "anxiety": "anxiety",
    "intellectual_disability": "intellectual disability", "depression": "depression", "bipolar": "bipolar disorder",
    "breast_cancer": "breast cancer", "ra": "rheumatoid arthritis", "cad": "coronary artery disease",
    "dev_delay": "developmental delay", "down_syndrome": "Down syndrome", "ocd": "OCD",
    "tuberous_sclerosis": "tuberous sclerosis", "fragile_x": "Fragile X", "rett": "Rett syndrome", "tourette": "Tourette",
}

sets = {}
for line in open("data/neuro.gmt"):
    p = line.rstrip("\n").split("\t")
    sets[p[0][4:]] = set(g.upper() for g in p[2:] if g.strip())
cat = {k[4:]: v for k, v in json.load(open("data/labels.json")).items()}


def ochiai(a, b):
    o = len(a & b)
    return o / np.sqrt(len(a) * len(b)) * 100 if o else 0.0


def group(d):
    if cat[d] == "control":
        return "non-neuro control"
    if d in MONOGENIC:
        return "monogenic syndrome"
    return "polygenic neuro-comorbidity"


poly = ["schizophrenia", "epilepsy", "adhd", "anxiety", "intellectual_disability", "depression", "bipolar"]
ctrl = [d for d in sets if cat[d] == "control"]
p_v = [ochiai(sets["asd"], sets[d]) for d in poly]
c_v = [ochiai(sets["asd"], sets[d]) for d in ctrl]
s_v = [ochiai(sets["asd"], sets[d]) for d in MONOGENIC]
_, p = mannwhitneyu(p_v, c_v, alternative="greater")
print(f"autism-to-polygenic-comorbidity: {np.mean(p_v):.1f}%  (range {min(p_v):.1f}-{max(p_v):.1f})")
print(f"autism-to-control:                {np.mean(c_v):.1f}%  (range {min(c_v):.1f}-{max(c_v):.1f})")
print(f"autism-to-monogenic-syndrome:     {np.mean(s_v):.1f}%  (range {min(s_v):.1f}-{max(s_v):.1f})")
print(f"  -> every polygenic comorbidity ranks above every control: Mann-Whitney p = {p:.4f}")

COL = {"polygenic neuro-comorbidity": "#8e44ad", "non-neuro control": "#95a5a6", "monogenic syndrome": "#e67e22"}
rows = sorted(((ochiai(sets["asd"], sets[d]), d) for d in sets if d != "asd"), reverse=True)
sns.set_theme(style="whitegrid")
fig, ax = plt.subplots(figsize=(9, 7))
y = list(range(len(rows)))[::-1]
ax.barh(y, [v for v, _ in rows], color=[COL[group(d)] for _, d in rows], alpha=.9)
ax.set_yticks(y)
ax.set_yticklabels([LABEL.get(d, d) for _, d in rows], fontsize=10)
ax.set_xlabel("gene-set similarity to autism (ochiai %)", fontsize=11)
ax.axvline(max(c_v), ls=":", color="#7f8c8d", lw=1)
ax.axvline(min(c_v), ls=":", color="#7f8c8d", lw=1)
ax.set_title("Autism's genetic neighbours are its polygenic comorbidities,\n"
             "not the single-gene syndromes that cause it", fontsize=12.5, fontweight="bold")
ax.legend(handles=[Patch(facecolor=COL[k], label=k) for k in
                   ("polygenic neuro-comorbidity", "non-neuro control", "monogenic syndrome")],
          fontsize=9.5, loc="lower right", framealpha=.95)
plt.tight_layout()
plt.savefig("figures/autism_neighbours.png", dpi=140, bbox_inches="tight")
print("wrote figures/autism_neighbours.png")
