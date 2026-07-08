# Recovering the Hallmarks of Cancer with `enrich`

Given only the names of specific cancers, DBRetina's `enrich` command recovers the oncogenic
signaling network from shared genes alone: guilt-by-association over the KEGG pathway similarity graph.

## Data

- `data/kegg.gmt` — 186 KEGG pathways (MSigDB `c2.cp.kegg`, gene symbols).
- `data/cancer_module.txt` — the 15 cancer-type pathways, defined from names only (written by `build_module.py`).
- `data/labels.json` — cancer / other class labels.

## Run

```bash
bash run.sh
```

It indexes the GMT, builds the pairwise, defines the cancer module, runs `enrich`, then validates and plots.

## What it produces

- `enrich_top.tsv` — the ranked candidates (top: focal adhesion, ErbB, mTOR, apoptosis, VEGF, ...).
- `validation.txt` — held-out recovery AUC and the permutation null.
- `figures/enrich_hallmarks.png` — top hits by significance, coloured by Hallmark of Cancer.
- `figures/enrich_recovery.png` — recovery AUC vs the random-module null.

## Result

Held-out (leave-one-out) recovery of the cancer class: **AUC 0.986**, versus a size-matched random-module
null of 0.48 +/- 0.10 (**4.9 SD above null, permutation p < 5e-4**; 0 of 2000 random modules reached it).
The top non-cancer hits are the Hallmarks of Cancer (Hanahan and Weinberg, 2011): proliferative signaling
(ErbB/EGFR, mTOR), resisting cell death (apoptosis), inducing angiogenesis (VEGF), and activating invasion
(focal adhesion).
