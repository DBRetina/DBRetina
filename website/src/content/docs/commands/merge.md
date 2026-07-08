---
title: merge
description: Merge two DBRetina indexes into one.
---

Combine two indexes into a single one. `merge` fuses them in feature space: a feature that appears in both
indexes is recognised as the same feature, so the merged index behaves exactly as if you had indexed all the
groups together from the start.

## When to use it

When your groups live in separate indexes and you want to compare them against each other. Merge them, then run
`pairwise` on the result.

:::tip
If you have the original GMT files, you can also build a combined index directly with
`index -g a.gmt -g b.gmt -o merged`. Use `merge` when you already have the `.dbri` files and not the sources.
:::

## Synopsis

```bash
DBRetina merge -a <indexA.dbri> -b <indexB.dbri> -o <merged.dbri>
```

## Inputs

| Option | Meaning | Required |
| --- | --- | --- |
| `-a, --index-a` | first (base) `.dbri` index | yes |
| `-b, --index-b` | second `.dbri` index to merge in | yes |
| `-o, --output` | output `.dbri` path | yes |

Index A is the base; every group from index B is added to it.

## Outputs

A single `.dbri` index (plus its `_raw.json` / `_hashes.json`) holding every group from both inputs.

## Example

```bash
DBRetina merge -a diseases_omim.dbri -b diseases_disgenet.dbri -o diseases_merged.dbri
```

```text
[SUCCESS] Total groups: 175 (27 from A + 148 from B)
```

## Notes

Group names must be **unique across both indexes**. If the same name appears in A and B, the merge stops with an
error rather than silently combining two different groups. Prefix your group names by source (as MSigDB does with
`KEGG_…`, `REACTOME_…`) to keep them distinct.

## Next

Run [pairwise](/DBRetina/commands/pairwise/) on the merged index to compare groups across the two sources.
