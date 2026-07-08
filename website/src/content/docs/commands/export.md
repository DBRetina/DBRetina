---
title: export
description: Export the similarities as a distance matrix, dendrogram, heatmap, and Newick tree.
---

Turn the pairwise result into the standard shapes for downstream tools and figures: a distance matrix, a
hierarchical-clustering dendrogram, a heatmap, and a Newick tree.

## When to use it

When you want a picture of the whole collection, or a matrix/tree to hand to another program.

## Synopsis

```bash
DBRetina export -p <pairwise> -m <metric> [--newick] [-l names|ids] [--linkage <method>] -o <prefix>
```

## Inputs

| Option | Meaning | Required |
| --- | --- | --- |
| `-p, --pairwise` | the pairwise output (Parquet dir, `.dbrp`, or `.tsv`) | yes |
| `-m, --metric` | metric to build the matrix from | yes |
| `--newick` | also write a Newick tree | no |
| `-l, --labels` | label rows/columns by `names` (default) or `ids` | no |
| `--linkage` | linkage for the dendrogram: `single`, `complete`, `average`, `weighted`, `centroid`, `median`, `ward` (default `ward`) | no |
| `-o, --output` | output prefix | yes |

## Outputs

| File | What it is |
| --- | --- |
| `<prefix>_distmat.tsv` | the matrix (similarity values, not distances, despite the name) |
| `<prefix>_distmat.pkl` | the same matrix as a pickled DataFrame |
| `<prefix>_dendrogram.png` | hierarchical-clustering dendrogram |
| `<prefix>_heatmap.png` | clustered heatmap |
| `<prefix>.newick` | Newick tree (with `--newick`) |

## Example

```bash
DBRetina export -p example_DBRetina_pairwise -m ochiai --newick -o example_tree
```

```text
[INFO] Writing newick tree to example_tree.newick
[SUCCESS] Done.
```

![Dendrogram exported from the example set](/DBRetina/img/quickstart_dendrogram.png)

## How it works

`export` turns one metric's pairwise values into a square, symmetric similarity matrix (self-similarity on the
diagonal set to 100) and writes it as `{prefix}_distmat.tsv` and `.pkl`, alongside a clustered
`{prefix}_heatmap.png`.

With `--newick` it also builds a tree, and the distance it clusters on is deliberately not the raw
`100 - similarity`. Each group is described by its whole row of dissimilarities to every other group, and the
distance between two groups is the **Euclidean distance between those two rows**, so groups that relate to the rest
of the panel in the same way sit close together. That distance feeds hierarchical clustering
(`scipy.cluster.hierarchy.linkage`, the method set by `--linkage`, default **ward**), and the tree is written as
`{prefix}.newick` and drawn as a circular `{prefix}_dendrogram.png`.

## Notes

All seven linkage methods work; `ward` is a good default for compact, evenly-sized clusters. The matrix file
holds the raw similarity values; read the name as "the matrix used for distance-based clustering", not as a
distance matrix.

## Next

Open the Newick tree in any tree viewer, or cluster the same network with
[cluster](/DBRetina/commands/cluster/).
