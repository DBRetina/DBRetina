---
title: geneinfo
description: List which groups each feature belongs to.
---

Turn the question around: instead of "what features does this group have?", ask "which groups contain this
feature?". `geneinfo` inverts an index and reports, for every feature in the groups you name, the groups it
appears in.

## When to use it

When you want to see how features are shared across a set of groups, for example, which genes drive the overlap
inside a cluster.

## Synopsis

```bash
DBRetina geneinfo -i <index> (-g <groups.txt> | --clusters-file <f> --cluster-ids <ids>) -o <prefix>
```

## Inputs

| Option | Meaning | Required |
| --- | --- | --- |
| `-i, --index-prefix` | the index to read | yes |
| `-g, --groups-file` | single-column file of group names to inspect | one of these |
| `--clusters-file` + `--cluster-ids` | a clusters file and the IDs to inspect | one of these |
| `-o, --output` | output prefix | yes |

## Outputs

`<prefix>_feature_to_groups.tsv`, one row per feature:

```text
feature  supergroups_count  supergroups
brca2    2                  breast_cancer|ovarian_cancer
tp53     2                  breast_cancer|ovarian_cancer
brca1    2                  breast_cancer|ovarian_cancer
palb2    2                  breast_cancer|ovarian_cancer
chek2    1                  breast_cancer
rad51c   1                  ovarian_cancer
```

## Example

```bash
DBRetina geneinfo -i example -g groups.txt -o example_genes
```

Here `groups.txt` lists `breast_cancer` and `ovarian_cancer`; the output shows that BRCA1, BRCA2, TP53, and PALB2
are shared between them, which is why the two come out so similar.

## Next

Find the groups behind a single feature with [gene-search](/DBRetina/commands/lookups/), or the shared features
of a pair with [shared-genes](/DBRetina/commands/lookups/).
