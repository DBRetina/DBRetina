#!/usr/bin/env python3
"""Draw the drug clustermap (coloured by mechanism) from the DBRetina pairwise output.

Expects drugs_DBRetina_pairwise.tsv (from run.sh). Writes figures/drug_clustermap.png.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.patches import Patch

CLASS = {
    "imatinib": "BCR-ABL", "dasatinib": "BCR-ABL", "nilotinib": "BCR-ABL",
    "bosutinib": "BCR-ABL", "ponatinib": "BCR-ABL",
    "gefitinib": "EGFR", "erlotinib": "EGFR", "afatinib": "EGFR", "lapatinib": "EGFR",
    "sorafenib": "VEGFR/multi", "sunitinib": "VEGFR/multi", "pazopanib": "VEGFR/multi",
    "regorafenib": "VEGFR/multi", "vandetanib": "VEGFR/multi", "cabozantinib": "VEGFR/multi",
    "vemurafenib": "BRAF", "dabrafenib": "BRAF",
    "crizotinib": "other", "trametinib": "other", "ruxolitinib": "other", "ibrutinib": "other",
    "vorinostat": "other", "tamoxifen": "other", "methotrexate": "other",
}
PALETTE = {"BCR-ABL": "#e8743b", "EGFR": "#2e86c1", "VEGFR/multi": "#28b463",
           "BRAF": "#8e44ad", "other": "#bdc3c7"}

drugs = sorted(CLASS)
idx = {d: i for i, d in enumerate(drugs)}
M = np.eye(len(drugs)) * 100
for line in open("drugs_DBRetina_pairwise.tsv"):
    if line.startswith("#"):
        continue
    p = line.rstrip("\n").split("\t")
    if p[0] == "group_1_ID" or p[2] not in idx or p[3] not in idx:
        continue
    M[idx[p[2]], idx[p[3]]] = M[idx[p[3]], idx[p[2]]] = float(p[6])

df = pd.DataFrame(M, index=[d.title() for d in drugs], columns=[d.title() for d in drugs])
rc = pd.Series({d.title(): PALETTE[CLASS[d]] for d in drugs}, name="")
sns.set_theme(style="white")
g = sns.clustermap(df, cmap="mako_r", annot=True, fmt=".0f", annot_kws={"size": 6},
                   row_colors=rc, col_colors=rc, linewidths=.4, linecolor="white",
                   figsize=(11, 10), cbar_kws={"label": "Ochiai similarity (%)"}, vmin=0, vmax=60)
g.ax_heatmap.set_xlabel(""); g.ax_heatmap.set_ylabel("")
g.ax_heatmap.legend(handles=[Patch(color=c, label=l) for l, c in PALETTE.items()],
                    loc="upper left", bbox_to_anchor=(1.22, 1.28), frameon=False, title="Mechanism")
g.figure.suptitle("Cancer drugs grouped by the gene targets they share", y=1.02, fontsize=14, fontweight="bold")
g.savefig("figures/drug_clustermap.png", dpi=140, bbox_inches="tight")
print("wrote figures/drug_clustermap.png")
