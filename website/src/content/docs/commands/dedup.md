---
title: dedup
description: Collapse near-duplicate groups into one representative each.
---

Find groups that are near-duplicates of each other and keep one representative from each set. `dedup` uses ochiai
similarity to decide what counts as a duplicate.

## When to use it

When a collection has redundant entries, the same pathway under two names, the same disease from two databases , 
and you want a non-redundant subset.

## Synopsis

```bash
DBRetina dedup -i <index> -p <pairwise> -c <cutoff> -o <prefix>
```

## Inputs

| Option | Meaning | Required |
| --- | --- | --- |
| `-i, --index-prefix` | the index the pairwise came from | yes |
| `-p, --pairwise` | the pairwise output (Parquet dir, `.dbrp`, or `.tsv`) | yes |
| `-c, --cutoff` | ochiai similarity above which two groups are duplicates (0–100) | yes |
| `-o, --output` | output prefix | yes |

## Outputs

`<prefix>_deduplicated_groups.txt`, the list of representative groups that survive deduplication.

## Example

Measure how redundant the WikiPathways collection is (733 pathways), treating two pathways with ochiai ≥ 80 as
duplicates:

```bash
DBRetina dedup -i wikipath -p wikipath_DBRetina_pairwise -c 80 -o wikipath_dedup
```

```text
[INFO] Number of gene sets after deduplication: 718
[INFO] Number of removed gene sets due to deduplication: 15
```

Only 15 of the 733 pathways are near-duplicates, so WikiPathways carries little internal redundancy. The
[database redundancy example](/DBRetina/examples/redundancy/) works through this in full.

## How it works

`dedup` works on the ochiai metric (there is no `-m`). It builds a graph joining every pair of groups with ochiai
at or above the cutoff, takes the connected components (so A~B and B~C group together even when A and C are not
directly similar), and keeps one representative per component while dropping the rest. The representative is the
most connected group, the one with the most such edges, ties broken by the largest gene set. Groups that match
nothing at the cutoff are kept untouched. The output, `{prefix}_deduplicated_groups.txt`, is the list of groups
that survive, one per line.

## Notes

A higher cutoff is stricter, only very-similar groups collapse. Lower it to merge more aggressively.

## Next

See where the duplicates come from with [cluster](/DBRetina/commands/cluster/), or reduce a database to a
covering subset with [setcov](/DBRetina/commands/setcov/).
