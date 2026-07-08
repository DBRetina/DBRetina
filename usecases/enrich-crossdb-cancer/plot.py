"""Figures for the cross-database enrich use case.

  figures/crossdb_recovered.png -- the independently-curated cancer pathways that
                                   the KEGG-only seed pulls to the top, by database.
  figures/crossdb_recovery.png  -- cross-database recovery AUC vs the random-seed null.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

DB_COLOR = {"wp": "#2ca25f", "reactome": "#2c7fb8", "biocarta": "#e6801a",
            "pid": "#7a4fa3", "other": "#95a5a6"}
DB_LABEL = {"wp": "WikiPathways", "reactome": "Reactome", "biocarta": "BioCarta",
            "pid": "PID", "other": "other"}
RED = "#c0392b"


import math
import re as _re


def short(name):
    for pre in ("wp_", "reactome_", "biocarta_", "pid_"):
        if name.startswith(pre):
            name = name[len(pre):]
    s = name.replace("_", " ")
    s = _re.sub(r"signaling by (.+?) in cancer", r"\1 signaling (cancer)", s)
    s = _re.sub(r"(pi3k akt) signaling in cancer", r"\1 signaling (cancer)", s)
    s = (s.replace("nonsmall cell lung cancer", "non-small cell lung cancer")
           .replace(" pathway", "").replace(" signaling pathways", " signaling")
           .replace("head and neck squamous cell carcinoma", "head & neck squamous carcinoma")
           .replace("epithelial to mesenchymal transition in ", "EMT in ")
           .replace("mirna regulation of prostate cancer signaling pathways", "prostate cancer (miRNA)"))
    s = s[:44].strip().title()
    for g in ("Egfr", "Erbb2", "Erbb3", "Alk", "Pi3K", "Pi3k", "Akt", "Wnt", "Mapk", "Emt", "Rb1"):
        s = s.replace(g, g.upper())
    return s


# --- Figure 1: recovered independent cancer pathways, top from EACH database ---
by_db = {"wp": [], "reactome": []}
with open("crossdb_recovered.tsv") as fh:
    next(fh)
    for line in fh:
        nm, db, hits, pv = line.rstrip("\n").split("\t")
        if db in by_db:
            by_db[db].append((nm, db, float(pv)))
recs = by_db["wp"][:8] + by_db["reactome"][:5]
recs.sort(key=lambda r: r[2])   # most significant first
recs = recs[::-1]               # so barh puts the most significant at the top
fig, ax = plt.subplots(figsize=(9.6, 6.4))
ys = range(len(recs))
ax.barh(list(ys), [-math.log10(p) for _, _d, p in recs],
        color=[DB_COLOR[d] for _n, d, _p in recs], edgecolor="white", height=0.72)
ax.set_yticks(list(ys))
ax.set_yticklabels([short(n) for n, _d, _p in recs], fontsize=9.5)
ax.set_xlabel("enrichment for the KEGG cancer seed  ( -log10 p )", fontsize=11)
ax.set_title("A KEGG-only cancer seed recovers cancers curated by other databases",
             fontsize=12.5, weight="bold")
seen = {}
for _n, d, _p in recs:
    seen[d] = DB_COLOR[d]
handles = [plt.Line2D([0], [0], marker="s", ls="", ms=10, mfc=c, mec=c, label=DB_LABEL[d])
           for d, c in seen.items()]
ax.legend(handles=handles, loc="lower right", fontsize=10, frameon=False, title="curated by")
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
fig.tight_layout()
fig.savefig("figures/crossdb_recovered.png", dpi=150, bbox_inches="tight")
plt.close(fig)

# --- Figure 2: where the independent cancer pathways rank among all non-KEGG ---
rows = []
with open("scores.tsv") as fh:
    next(fh)
    for line in fh:
        nm, isc, sc = line.rstrip("\n").split("\t")
        rows.append((float(sc), int(isc)))
rows.sort(key=lambda r: -r[0])            # rank 1 = strongest enrichment for the KEGG seed
n = len(rows)
cancer_ranks = [i + 1 for i, (_s, isc) in enumerate(rows) if isc]
val = {}
for line in open("validation.txt"):
    k, v = line.rstrip("\n").split("\t")
    val[k] = v
top10 = sum(1 for r in cancer_ranks if r <= n * 0.10)
median_rank = sorted(cancer_ranks)[len(cancer_ranks) // 2]

fig, ax = plt.subplots(figsize=(9.6, 3.2))
ax.axvspan(1, n * 0.10, color="#f2d7d3", zorder=0)                       # the top 10% band
ax.hlines(0, 1, n, color="#e6e6e6", lw=10, zorder=1)
for r in cancer_ranks:
    ax.vlines(r, -0.5, 0.5, color=RED, lw=1.1, alpha=0.85, zorder=2)     # each independent cancer pathway
ax.set_xlim(-n * 0.02, n * 1.02)
ax.set_ylim(-1.9, 1.7)
ax.set_yticks([])
ax.annotate(f"{top10} of the {len(cancer_ranks)} independent cancer pathways\nrank in the top 10% (shaded); "
            f"median rank {median_rank} of {n}",
            (n * 0.05, 0.7), (n * 0.30, 1.25), fontsize=10.5, weight="bold", color=RED,
            arrowprops=dict(arrowstyle="->", color=RED))
ax.text(n, -1.5, f"AUC {float(val['real_auc']):.2f}   (random seeds {float(val['null_auc_mean']):.2f}, p < 5e-4)",
        ha="right", fontsize=9.5, color="#444", style="italic")
ax.set_xlabel(f"rank by enrichment for the KEGG cancer seed   (1 = strongest, of {n} pathways from the other databases)",
              fontsize=10)
ax.set_title("The cancer pathways from other databases rank near the top", fontsize=12.5, weight="bold")
for s in ("top", "right", "left"):
    ax.spines[s].set_visible(False)
fig.tight_layout()
fig.savefig("figures/crossdb_recovery.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print("wrote figures/crossdb_recovered.png and figures/crossdb_recovery.png")
