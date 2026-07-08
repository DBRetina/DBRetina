---
title: bipartite
description: Build a bipartite similarity graph between two sets of groups.
---

Compare one set of groups against another, rather than everything against everything. `bipartite` builds the
connections between two named sides, useful for asking "which of these diseases resemble which of those?"

## When to use it

When you have two lists of groups and want the cross-similarities between them, not the full all-pairs matrix.

## Synopsis

```bash
DBRetina bipartite -p <pairwise> (--group1 <f> --group2 <f> | --gmt1 <f> --gmt2 <f>) -m <metric> [-c <cutoff>] -o <prefix>
```

## Inputs

| Option | Meaning | Required |
| --- | --- | --- |
| `-p, --pairwise` | the pairwise output (Parquet dir, `.dbrp`, or `.tsv`) | yes |
| `--group1` / `--group2` | single-column files naming the two sides | one pair |
| `--gmt1` / `--gmt2` | GMT files defining the two sides instead | one pair |
| `-m, --metric` | metric for the connections | yes |
| `-c, --cutoff` | keep connections above this (0–100) | no |
| `--no-plot` | skip the plots | no |
| `--no-1-1` | drop one-to-one mappings | no |
| `-o, --output` | output prefix | yes |

## Outputs

| File | What it is |
| --- | --- |
| `<prefix>_bipartite_pairwise.tsv` | the cross-connections with their metrics |
| `<prefix>_similarity_metrics_histogram.png` / `.json` | distribution of the connection strengths |

## Example

```bash
DBRetina bipartite -p example_DBRetina_pairwise --group1 grp1.txt --group2 grp2.txt -m ochiai -c 0 -o example_bp
```

```text
group_1         group_2      containment  ochiai  jaccard
breast_cancer   melanoma     33.3         30.9    18.2
ovarian_cancer  lung_cancer  16.7         16.7    9.1
```

## Notes

Rendering the plot to a static image needs the optional `kaleido` package. If it isn't installed, `bipartite`
still writes the data and warns instead of failing; pass `--no-plot` to skip plotting entirely.

## Next

Filter the connections with [query](/DBRetina/commands/query/), or look at the shared genes behind a connection
with [shared-genes](/DBRetina/commands/lookups/).
