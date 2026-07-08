---
title: Quick lookups
description: Four one-line commands for reading a pairwise result and its index.
---

Four small commands answer the questions you ask most often while exploring a result. Each takes the pairwise
output with `-d` (a Parquet directory or a `.tsv`), and prints to the screen unless you give it `-o`.

## search: find groups by name

Find groups whose name contains a pattern.

```bash
DBRetina search -d <pairwise> <pattern> [-o <file>]
```

```bash
DBRetina search -d example_DBRetina_pairwise "cancer"
```

```text
[INFO] Found 3 groups matching "cancer"
group_id  group_name
2         breast_cancer
6         lung_cancer
3         ovarian_cancer
```

## neighbors: nearest groups to one group

Show the top co-occurring groups for a target, sorted by a metric.

```bash
DBRetina neighbors -d <pairwise> <group> -m <metric> -c <cutoff> [-n <count>] [-o <file>]
```

```bash
DBRetina neighbors -d example_DBRetina_pairwise "breast_cancer" -m ochiai -c 0
```

```text
neighbor         ochiai  jaccard  shared_features
ovarian_cancer   61.7    44.4     4
melanoma         30.9    18.2     2
lung_cancer      15.4    8.3      1
```

`-n` caps the number of neighbours (default 20).

## shared-genes: features common to two groups

List the features two groups have in common. This one also needs the index (`-i`) for the gene data.

```bash
DBRetina shared-genes -d <pairwise> -i <index> <group_a> <group_b> [-o <file>]
```

```bash
DBRetina shared-genes -d example_DBRetina_pairwise -i example "breast_cancer" "ovarian_cancer"
```

```text
[INFO] 4 shared features between "breast_cancer" and "ovarian_cancer"
feature
brca1
brca2
palb2
tp53
```

## gene-search: groups containing a feature

Find every group that contains a given feature. The search is case-insensitive and needs the index (`-i`).

```bash
DBRetina gene-search -d <pairwise> -i <index> <feature> [-o <file>]
```

```bash
DBRetina gene-search -d example_DBRetina_pairwise -i example "TP53"
```

```text
[INFO] TP53 found in 4 groups
group_name
breast_cancer
lung_cancer
melanoma
ovarian_cancer
```

## Next

For a feature-by-feature breakdown across several groups at once, use
[geneinfo](/DBRetina/commands/geneinfo/).
