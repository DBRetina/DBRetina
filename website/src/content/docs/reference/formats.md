---
title: File formats
description: The input formats DBRetina reads and the files it writes.
---

## Input

DBRetina reads your data in one of two plain-text formats.

### GMT

Tab-delimited, no header. Each row is one group: name, a description column (often unused), then the features.

```text
breast_cancer  na  BRCA1  BRCA2  TP53  PTEN
lung_cancer    na  EGFR   KRAS   TP53  ALK
```

### Association (ASC)

A two-column TSV **with** a header, one feature per row.

```text
group          feature
breast_cancer  BRCA1
breast_cancer  BRCA2
lung_cancer    EGFR
```

:::note
Names are lowercased on import, double quotes are stripped, and `|` is reserved. GMT input must be uncompressed , 
gunzip a `.gmt.gz` first.
:::

## The index

`index` writes three files per prefix:

| File | Contents |
| --- | --- |
| `<prefix>.dbri` | the binary index, the file every command reads |
| `<prefix>_raw.json` | your groups and their features, in plain text |
| `<prefix>_hashes.json` | the same, with features hashed (used internally) |

## The pairwise output

`pairwise` writes a **Parquet directory**, `<index>_DBRetina_pairwise/`:

| Entry | Contents |
| --- | --- |
| `manifest.json` | row and partition counts, and whether a p-value column is present |
| `names.parquet` | the group-name lookup |
| `group_index.parquet` | per-group index data |
| `data/part_*.parquet` | the pairwise rows, partitioned |
| `statistics.json` | summary statistics of the run |

This directory is what the analysis commands expect for `-p` (or `-d`). Add `--legacy-output` to also get a flat
`.tsv` and the binary `.dbrp`, which some commands still accept as input.

The pairwise columns are:

```text
group_1_ID  group_2_ID  group_1_name  group_2_name  shared_features
containment  ochiai  jaccard  csi  dice  odds_ratio  [pvalue]
```
