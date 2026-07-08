#!/usr/bin/env python3
"""Show that type 1 diabetes leans autoimmune and type 2 leans metabolic.

Reads data/panel.gmt and data/labels.json, compares each diabetes type's mean similarity to the
autoimmune and metabolic groups, and writes figures/diabetes_split.png.
"""
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

sets = {}
for line in open("data/panel.gmt"):
    p = line.rstrip("\n").split("\t")
    sets[p[0][4:]] = set(g.upper() for g in p[2:] if g.strip())
cat = {k[4:]: v for k, v in json.load(open("data/labels.json")).items()}


def ochiai(a, b):
    o = len(a & b)
    return o / np.sqrt(len(a) * len(b)) * 100 if o else 0.0


auto = [d for d in sets if cat[d] == "autoimmune"]
metab = [d for d in sets if cat[d] == "metabolic"]
mean = lambda seed, grp: np.mean([ochiai(sets[seed], sets[d]) for d in grp])
t1 = [mean("t1d", auto), mean("t1d", metab)]
t2 = [mean("t2d", auto), mean("t2d", metab)]
print(f"type 1 diabetes -> autoimmune {t1[0]:.0f}%, metabolic {t1[1]:.0f}%")
print(f"type 2 diabetes -> autoimmune {t2[0]:.0f}%, metabolic {t2[1]:.0f}%")

sns.set_theme(style="whitegrid")
fig, ax = plt.subplots(figsize=(8, 6))
x = np.arange(2)
w = 0.36
b1 = ax.bar(x - w / 2, [t1[0], t2[0]], w, label="to autoimmune diseases", color="#c0392b", alpha=.9)
b2 = ax.bar(x + w / 2, [t1[1], t2[1]], w, label="to metabolic diseases", color="#e67e22", alpha=.9)
for b in list(b1) + list(b2):
    ax.annotate(f"{b.get_height():.0f}%", (b.get_x() + b.get_width() / 2, b.get_height()),
                ha="center", va="bottom", fontsize=11, fontweight="bold")
ax.set_xticks(x)
ax.set_xticklabels(["type 1 diabetes\n(autoimmune)", "type 2 diabetes\n(metabolic)"], fontsize=11.5)
ax.set_ylabel("mean gene-set similarity (ochiai %)", fontsize=11)
ax.set_ylim(0, 50)
ax.legend(fontsize=11)
ax.set_title("Same name, different genetics:\ntype 1 diabetes leans autoimmune, type 2 leans metabolic",
             fontsize=13, fontweight="bold")
plt.tight_layout()
plt.savefig("figures/diabetes_split.png", dpi=140, bbox_inches="tight")
print("wrote figures/diabetes_split.png")
