---
title: query
description: Pull specific groups, clusters, or strong pairs out of a pairwise result.
---

Filter a pairwise result down to what you care about, a handful of named groups, the members of a cluster, or
just the pairs above a threshold.

## When to use it

After `pairwise`, when you want a focused slice of the full similarity table rather than the whole thing.

## Synopsis

```bash
DBRetina query -p <pairwise> [-g <groups.txt> | --clusters-file <f> --cluster-ids <ids>] [-m <metric> -c <cutoff>] [--extend] [-o <prefix>] [--tsv | --tsv-only | --inplace]
```

## Inputs

| Option | Meaning | Required |
| --- | --- | --- |
| `-p, --pairwise` | the pairwise output, a Parquet directory, `.dbrp`, or `.tsv` | yes |
| `-g, --groups-file` | single-column file of group names to keep | no |
| `--clusters-file` | a clusters file from `cluster` | no |
| `--cluster-ids` | comma-separated cluster IDs to pull from that file | no |
| `-m, --metric` | metric to filter on | no |
| `-c, --cutoff` | keep pairs at/above this (or at/below, for `pvalue`) | no |
| `--extend` | also include every group linked to your selection | no |
| `-o, --output` | output prefix; a Parquet input writes `<prefix>_DBRetina_pairwise/` | yes, unless `--inplace` |
| `--tsv` | also write a legacy `.tsv` beside the directory | no |
| `--tsv-only` | write only the legacy `.tsv`, no directory | no |
| `--inplace` | filter the input directory in place, no new output | no |

## Outputs

When the input is a Parquet directory, `query` writes a new Parquet pairwise directory,
`<prefix>_DBRetina_pairwise/`. That directory is itself a complete pairwise, so you can run
`cluster`, `export`, or `query` on it again, with the metrics kept at full precision.

`--tsv` writes that directory and a `<prefix>.tsv` beside it; `--tsv-only` writes just the
`.tsv` (the earlier behaviour); `--inplace` rewrites the input directory in place and takes no
`-o`. A legacy `.tsv` or `.dbrp` input writes a `<prefix>.tsv`. P-values are not carried into
the filtered directory, so use `--tsv-only` when you need to keep them.

## Example

Keep everything involving breast cancer, and everything linked to it, as a new pairwise you can
carry on analysing:

```bash
DBRetina query -p example_DBRetina_pairwise -g groups.txt --extend -o breast_neighbourhood
# breast_neighbourhood_DBRetina_pairwise/  is a full pairwise
DBRetina cluster -p breast_neighbourhood_DBRetina_pairwise -m ochiai -c 50 --community -o bn
```

Add `--tsv` (or `--tsv-only`) to also get the matching pairs as a table:

```text
group_1_name    group_2_name    shared_features  containment  ochiai  jaccard
breast_cancer   ovarian_cancer  4                66.7         61.7    44.4
breast_cancer   melanoma        2                33.3         30.9    18.2
breast_cancer   lung_cancer     1                16.7         15.4    8.3
```

## Notes

- `-p` accepts all three pairwise forms, the Parquet directory, the `.dbrp`, or a flat `.tsv`.
- The default output is a Parquet directory, so a filtered slice is a first-class pairwise you can
  keep chaining. `--tsv-only` gives the earlier flat-`.tsv` behaviour.
- `--inplace` overwrites the input directory with the filtered pairs and keeps no copy of the
  original, so reach for it only when you no longer need the full result.
- `--extend` starts from your named groups and pulls in everything connected to them, which is useful
  for building a neighbourhood around a few seeds.
- With `-m pvalue`, the cutoff keeps the *most* significant pairs (p-value at or below the cutoff).

## Next

Cluster the slice with [cluster](/DBRetina/commands/cluster/) or draw it with
[export](/DBRetina/commands/export/).
