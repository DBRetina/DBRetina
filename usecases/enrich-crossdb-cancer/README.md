# Cancer's shared genes replicate across databases (`enrich`)

Does cancer's shared-gene structure survive independent curation? Seeded with KEGG's 15 cancer
pathways, DBRetina `enrich` recovers the cancer pathways that Reactome and WikiPathways curated
separately, out of the 3090-pathway MSigDB C2:CP compendium.

## Data

- `data/c2cp.gmt` — MSigDB C2:CP: 3090 pathways from KEGG, Reactome, WikiPathways, BioCarta, PID.
- `data/kegg_cancer_seed.txt` — the 15 KEGG cancer pathways (the seed), by name (written by `build_labels.py`).
- `data/independent_cancer.txt` — the 64 Reactome + WikiPathways cancer pathways (the held-out truth).

## Run

```bash
bash run.sh
```

## What it produces

- `enrich_top.tsv` — pathways ranked by enrichment for the KEGG cancer seed.
- `crossdb_recovered.tsv` — the independent cancer pathways, ranked by that enrichment.
- `validation.txt` — cross-database recovery AUC and the random-seed null.
- `figures/crossdb_recovered.png` — recovered cancer pathways, coloured by source database.
- `figures/crossdb_recovery.png` — recovery AUC vs the random-seed null.

## Result

Cross-database recovery of the independently-curated cancer pathways: **AUC 0.728**, versus a
random-seed null of 0.53 +/- 0.05 (**4.1 SD above null, permutation p < 5e-4**; 0 of 2000 random
seeds reached it). The recovered pathways are the known cancers and drivers, each curated by a
different team: bladder, breast, melanoma, lung, pancreatic, and colorectal cancer from
WikiPathways, and EGFR, PI3K-AKT, TGF-beta, and RB1-loss oncogenic signaling from Reactome. The
shared-gene structure of cancer crosses curation boundaries, so it is biology, not a database
artifact.
