"""Figures for the enrich / Hallmarks-of-Cancer use case.

  figures/enrich_hallmarks.png -- the pathways enrich flags, with the evidence
                                  (how many of each one's neighbours are cancer)
                                  and the Hallmark each embodies.
  figures/enrich_recovery.png  -- where the cancer pathways rank when each is
                                  held out and scored against the others.
"""
import math

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

GREY = "#95a5a6"
RED = "#c0392b"

# each recovered pathway -> the Hallmark of Cancer it embodies (Hanahan & Weinberg, 2011)
HALLMARK = {
    "erbb_signaling_pathway": ("Proliferative signaling", "#2c7fb8"),
    "mtor_signaling_pathway": ("Proliferative signaling", "#2c7fb8"),
    "insulin_signaling_pathway": ("Proliferative signaling", "#2c7fb8"),
    "neurotrophin_signaling_pathway": ("Proliferative signaling", "#2c7fb8"),
    "apoptosis": ("Resisting cell death", "#7a4fa3"),
    "vegf_signaling_pathway": ("Inducing angiogenesis", "#2ca25f"),
    "focal_adhesion": ("Invasion & metastasis", "#e6801a"),
    "regulation_of_actin_cytoskeleton": ("Invasion & metastasis", "#e6801a"),
    "b_cell_receptor_signaling_pathway": ("Shared immune-receptor kinases", GREY),
    "t_cell_receptor_signaling_pathway": ("Shared immune-receptor kinases", GREY),
    "toll_like_receptor_signaling_pathway": ("Shared immune-receptor kinases", GREY),
    "fc_epsilon_ri_signaling_pathway": ("Shared immune-receptor kinases", GREY),
    "fc_gamma_r_mediated_phagocytosis": ("Shared immune-receptor kinases", GREY),
}


def short(name):
    s = name.replace("kegg_", "").replace("_", " ")
    s = (s.replace("signaling pathway", "signaling").replace(" receptor", " rec.")
           .replace("regulation of ", "").replace(" mediated", ""))
    return s.title()


# --- Figure 1: the flagged pathways, with the guilt-by-association evidence ---
recs = []
with open("enrich_top.tsv") as fh:
    next(fh)
    for line in fh:
        nm, hits, deg, pv, _in = line.rstrip("\n").split("\t")
        key = nm.replace("kegg_", "")
        if key in HALLMARK:
            recs.append((nm, int(hits), int(deg), float(pv), HALLMARK[key]))
recs = recs[:11][::-1]

fig, ax = plt.subplots(figsize=(10.5, 6.4))
ys = list(range(len(recs)))
vals = [-math.log10(p) for _n, _h, _d, p, _hm in recs]
ax.barh(ys, vals, color=[hm[1] for *_r, hm in recs], edgecolor="white", height=0.7)
for i, (nm, hits, deg, pv, _hm) in enumerate(recs):
    ax.text(-math.log10(pv) + 0.15, i, f"{hits} of its {deg} neighbours are cancer pathways",
            va="center", fontsize=8.3, color="#444")
ax.set_yticks(ys)
ax.set_yticklabels([short(n) for n, *_r in recs], fontsize=10)
ax.set_xlim(0, max(vals) * 1.6)
ax.set_xlabel("strength of the cancer signal in its neighbourhood   (higher = less likely by chance)", fontsize=10.5)
ax.set_title("Given only the cancer pathways, enrich returns cancer's machinery",
             fontsize=13, weight="bold")
seen = {}
for *_r, hm in recs:
    seen[hm[0]] = hm[1]
handles = [plt.Line2D([0], [0], marker="s", ls="", ms=10, mfc=c, mec=c, label=l) for l, c in seen.items()]
ax.legend(handles=handles, loc="lower right", fontsize=9, frameon=False, title="Hallmark of Cancer")
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
fig.tight_layout()
fig.savefig("figures/enrich_hallmarks.png", dpi=150, bbox_inches="tight")
plt.close(fig)

# --- Figure 2: where the cancer pathways rank (held-out) ---
rows = []
with open("scores.tsv") as fh:
    next(fh)
    for line in fh:
        nm, isc, sc = line.rstrip("\n").split("\t")
        rows.append((float(sc), int(isc)))
rows.sort(key=lambda r: -r[0])          # rank 1 = strongest enrichment
n = len(rows)
cancer_ranks = [i + 1 for i, (_s, isc) in enumerate(rows) if isc]
val = {}
for line in open("validation.txt"):
    k, v = line.rstrip("\n").split("\t")
    val[k] = v

fig, ax = plt.subplots(figsize=(9.5, 2.9))
ax.hlines(0, 1, n, color="#e6e6e6", lw=10, zorder=1)                    # every pathway's rank
for r in cancer_ranks:
    ax.vlines(r, -0.5, 0.5, color=RED, lw=1.8, zorder=2)               # the 15 cancer pathways
ax.set_xlim(-3, n + 3)
ax.set_ylim(-1.8, 1.6)
ax.set_yticks([])
worst = max(cancer_ranks)
ax.annotate(f"all {len(cancer_ranks)} cancer pathways land in the top {worst} of {n}",
            (worst, 0.7), (worst + 12, 1.15), fontsize=11, weight="bold", color=RED,
            arrowprops=dict(arrowstyle="->", color=RED))
ax.text(n, -1.35, f"AUC {float(val['real_auc']):.2f}   (random modules {float(val['null_auc_mean']):.2f}, "
        f"p < 5e-4)", ha="right", fontsize=9.5, color="#444", style="italic")
ax.set_xlabel("rank by cancer-neighbourhood enrichment   (1 = strongest, each pathway held out in turn)",
              fontsize=10.5)
ax.set_title("Where the cancer pathways rank among all 154: at the very top", fontsize=12.5, weight="bold")
for s in ("top", "right", "left"):
    ax.spines[s].set_visible(False)
fig.tight_layout()
fig.savefig("figures/enrich_recovery.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print("wrote figures/enrich_hallmarks.png and figures/enrich_recovery.png")
