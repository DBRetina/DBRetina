---
title: pairwise
description: Compare every pair of groups and compute their similarity.
---

Compare all groups in an index against each other. `pairwise` is the heart of DBRetina: it produces the table of
similarities that `query`, `cluster`, `export`, and every other analysis reads.

## When to use it

Right after `index`. Run it once per index (or once per metric/threshold you care about) and reuse the output.

## Synopsis

```bash
DBRetina pairwise -i <index> -m <metric> -c <cutoff> [-t <threads>] [--pvalue] [--legacy-output]
```

## Inputs

| Option | Meaning | Required |
| --- | --- | --- |
| `-i, --index-prefix` | index prefix (from `index`) | yes |
| `-m, --metric` | metric the cutoff applies to: `containment`, `ochiai`, or `jaccard` | no (default keeps all) |
| `-c, --cutoff` | drop pairs whose metric is below this (0–100) | no (default `0`) |
| `-t, --threads` | number of cores to use | no (default `1`) |
| `--pvalue` | also compute a hypergeometric p-value per pair | no |
| `--fdr` | apply Benjamini–Hochberg FDR correction (needs `--pvalue`) | no |
| `--fdr-alpha` | FDR significance threshold (default `0.05`) | no |
| `--legacy-output` | also write the flat `.tsv` and `.dbrp` files | no |
| `--no-progress` / `--debug` | quiet the progress bar / print diagnostics | no |

## Outputs

By default `pairwise` writes a **Parquet directory**, `<index>_DBRetina_pairwise/`, which is what the analysis
commands read. It also writes `<index>_DBRetina_featuresNo.tsv` (feature counts per group) and two distribution
plots.

Each pair has these columns:

```text
group_1_ID  group_2_ID  group_1_name  group_2_name  shared_features
containment  ochiai  jaccard  csi  dice  odds_ratio
```

Add `--pvalue` for a `pvalue` column, `--fdr` for a `<index>_DBRetina_pairwise_fdr.tsv`, and `--legacy-output`
for a flat `.tsv` plus the binary `.dbrp`.

## Example

```bash
DBRetina pairwise -i example -m ochiai -c 0
```

Use all cores and keep only strong pairs:

```bash
DBRetina pairwise -i pathways -m ochiai -c 50 -t 16
```

## How it works

At index time each group's genes are hashed, and identical genes across groups collapse to a shared colour, so
`pairwise` never re-intersects gene lists at run time: it walks that colour index and, for every pair of groups
that share at least one gene, accumulates the shared-gene count directly. Pairs that share nothing never appear.
For each shared pair it computes the similarity metrics from three numbers, the shared count and the two group
sizes (|A∩B|, |A|, |B|), plus the feature-universe size for the size-aware ones. The exact formulas are in the
[metrics reference](/DBRetina/reference/metrics/); note that `containment` is symmetric (it divides by the smaller
set), and `csi` here is the squared cosine, not the graph-theoretic index the name suggests.

`-c` keeps only the pairs whose chosen metric (`-m`, one of containment, ochiai, or jaccard) clears the cutoff,
though every metric is still written for the pairs that pass. `--pvalue` adds a one-sided hypergeometric
over-representation p-value per pair (the right tail of Fisher's exact test on the same 2x2 table as the odds
ratio), and `--fdr` adds Benjamini-Hochberg q-values on top.

## Notes

- `pairwise` runs in parallel, so pass `-t` to use more cores. On a 16-core node it compares a 31.8-million-pair
  disease set in about 34 seconds.
- All similarity metrics are computed for every pair regardless of `-m`; the metric you pass only decides which
  one the `-c` cutoff filters on.
- The p-value is *lower is better*. A p-value threshold keeps the pairs at or below it, opposite to the
  similarity metrics.

## Next

Filter the result with [query](/DBRetina/commands/query/), group it with
[cluster](/DBRetina/commands/cluster/), or draw it with [export](/DBRetina/commands/export/).
