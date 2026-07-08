---
title: setcov
description: Reduce a database to a small covering subset with a set-cover algorithm.
---

Pick a small set of groups that together cover the whole collection. `setcov` applies a set-cover algorithm so
you can shrink a redundant database down to a representative core.

## When to use it

When you want the *smallest* useful subset of a database, a covering set, rather than just removing duplicates.

## Synopsis

```bash
DBRetina setcov -i <index> [--modularity <c>] [--dedup <c>] [--stop-cov <pct>] -o <prefix>
```

## Inputs

| Option | Meaning | Required |
| --- | --- | --- |
| `-i, --index-prefix` | the index to reduce | yes |
| `--modularity` | containment cutoff for the internal modularity step (default `80`) | no |
| `--dedup` | ochiai cutoff for the internal dedup step (default `100`) | no |
| `--stop-cov` | stop once this percentage of items is covered | no |
| `-o, --output` | output prefix | yes |

`setcov` reads only the index; it runs the modularity, deduplication, and community steps internally.

## Outputs

| File | What it is |
| --- | --- |
| `<prefix>_new.gmt` | the covering subset as a GMT |
| `<prefix>_original.gmt` | the full collection as a GMT |
| `<prefix>_associations.tsv` | the reduced associations |
| `<prefix>_groups_metadata.tsv` | metadata for the groups |
| `<prefix>_remaining_groups_metadata.tsv` | groups kept in the cover |
| `<prefix>_removed_groups_metadata.tsv` | groups dropped |
| `<prefix>_item_to_GPI_CSI.tsv` | per-item GPI/CSI scores |

## Example

```bash
DBRetina setcov -i pathways_multidb --modularity 50 --stop-cov 90 -o pathways_cover
```

```text
[SUCCESS] Group modularities computed with containment cutoff 50.0
[SUCCESS] Deduplication completed
[SUCCESS] Set coverage process completed
[SUCCESS] The new association file exported to pathways_cover_associations.tsv
```

:::caution
The internal community step runs at the `--modularity` cutoff. On sparse data the default of `80` can leave no
edges; if `setcov` reports the data is too sparse, lower `--modularity` (for example to `30`–`50`).
:::

## How it works

`setcov` reduces a database to a smaller set of groups that still covers almost all the features, reading the index
(`-i`) and the pairwise built from the same prefix, in two stages.

First it removes near-exact duplicates: groups linked at very high ochiai (`--dedup`, default 100) are collapsed to
a single representative each (the most connected, tie-broken by the largest set).

Then it runs a static-priority greedy set cover over the survivors. Rather than re-ranking by marginal gain at each
step, it orders the candidates once, the cleanest and most informative first: lowest containment
[`modularity`](/DBRetina/commands/modularity/) (least redundant), then highest cluster specificity, then largest
set. It walks that fixed order, keeps any group that still adds features not yet covered, and stops once
`--stop-cov` percent of all features are covered. The reduced set is written as `{prefix}_new.gmt`, alongside
per-group metadata recording why each group was kept or dropped.

## Next

Build a fresh index from the covering GMT with [index](/DBRetina/commands/), or compare the kept groups with
[pairwise](/DBRetina/commands/pairwise/).
