#!/usr/bin/env python3
"""Draw the rheumatoid-arthritis mechanism/treatment figure from the DBRetina pairwise output.

Expects multidb_DBRetina_pairwise.tsv (from run.sh). Writes figures/multidb_ra.png.
Pulls rheumatoid arthritis's neighbours and splits them into pathways (by source prefix) and drugs.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Patch

RA = "dis_rheumatoid_arthritis"
SRC = {"kegg_": ("KEGG", "#e67e22"), "reactome_": ("Reactome", "#28b463"), "wp_": ("WikiPathways", "#2e86c1")}

paths, drugs = [], []
for line in open("multidb_DBRetina_pairwise.tsv"):
    if line.startswith("#"):
        continue
    p = line.rstrip("\n").split("\t")
    if p[0] == "group_1_ID":
        continue
    other = p[3] if p[2] == RA else (p[2] if p[3] == RA else None)
    if not other:
        continue
    containment, shared = float(p[5]), int(p[4])
    if other.startswith("drug_"):
        drugs.append((other[5:].title(), containment))
    else:
        for pref, (label, colour) in SRC.items():
            if other.startswith(pref):
                name = other[len(pref):].replace("_", " ")
                name = (name[:34] + "...") if len(name) > 37 else name
                paths.append((name + f" ({label[:2]})", shared, colour))
                break

# top pathways by shared genes (gene-rich, RA-relevant), top drugs by containment
paths = sorted(paths, key=lambda x: -x[1])[:7][::-1]
drugs = sorted(drugs, key=lambda x: -x[1])[:7][::-1]

sns.set_theme(style="whitegrid")
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.5))
ax1.barh(range(len(paths)), [p[1] for p in paths], color=[p[2] for p in paths], alpha=.9)
ax1.set_yticks(range(len(paths))); ax1.set_yticklabels([p[0] for p in paths], fontsize=9)
ax1.set_xlabel("RA-associated genes in the pathway", fontsize=10)
ax1.set_title("Mechanism: pathways that connect to RA", fontsize=12, fontweight="bold")
ax1.legend(handles=[Patch(color=c, label=l) for l, c in
                    [("Reactome", "#28b463"), ("WikiPathways", "#2e86c1"), ("KEGG", "#e67e22")]],
           loc="lower right", frameon=True, fontsize=9)
ax2.barh(range(len(drugs)), [d[1] for d in drugs], color="#c0392b", alpha=.9)
ax2.set_yticks(range(len(drugs))); ax2.set_yticklabels([d[0] for d in drugs], fontsize=9)
ax2.set_xlabel("% of the drug's targets that are RA genes (containment)", fontsize=10)
ax2.set_title("Treatment: drugs that connect to RA", fontsize=12, fontweight="bold")
ax2.set_xlim(0, 100)
fig.suptitle("What connects to rheumatoid arthritis across 5 databases", fontsize=14, fontweight="bold", y=1.03)
plt.tight_layout()
plt.savefig("figures/multidb_ra.png", dpi=140, bbox_inches="tight")
print("wrote figures/multidb_ra.png")
