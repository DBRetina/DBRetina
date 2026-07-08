---
title: cluster
description: Group similar sets together by clustering the similarity graph.
---

Turn the pairwise similarities into groups. `cluster` treats the pairwise result as a graph (groups are nodes,
strong similarities are edges) and partitions it into clusters of related sets.

## When to use it

After `pairwise`, when you want to find families of similar groups rather than inspect individual pairs.

## Synopsis

```bash
DBRetina cluster -p <pairwise> -m <metric> -c <cutoff> [--community] [--resolution <r>] -o <prefix>
```

## Inputs

| Option | Meaning | Required |
| --- | --- | --- |
| `-p, --pairwise` | the pairwise output (Parquet dir, `.dbrp`, or `.tsv`) | yes |
| `-m, --metric` | metric used for the edges | yes |
| `-c, --cutoff` | minimum metric value for an edge (0–100) | yes |
| `--community` | use Leiden community detection; without it, connected components | no |
| `--resolution` | Leiden resolution, higher gives more, smaller clusters | no |
| `--node-weight-transform` | `log2`, `linear`, or `sqrt` weighting of group sizes | no |
| `-o, --output` | output prefix | yes |

## Outputs

| File | What it is |
| --- | --- |
| `<prefix>_clusters.tsv` | each cluster, its size, and its members (pipe-separated) |
| `<prefix>_clusters_histogram.png` | distribution of cluster sizes |
| `<prefix>_clusters_bubbles.png` | bubble plot of the clusters |

## Example

```bash
DBRetina cluster -p example_DBRetina_pairwise -m ochiai -c 20 --community -o example_clusters
```

```text
[INFO] Total number of clustered supergroups: 3
[INFO] number of clusters: 1
```

```text title="example_clusters_clusters.tsv"
cluster_id  cluster_size  cluster_members
1           3             breast_cancer|ovarian_cancer|melanoma
```

## How it works

`cluster` builds a graph in which every pair of groups whose similarity clears the cutoff (`-m`/`-c`) is an edge,
weighted by that similarity.

By default it returns the **connected components** of that graph: any groups linked, directly or through a chain of
above-cutoff edges, land in the same cluster (single-linkage at the fixed threshold). Groups with no passing edge
are dropped.

With `--community` it instead runs the **Leiden** algorithm (the Constant Potts Model partition), which carves out
tightly interlinked communities rather than merely connected ones. `--resolution` sets the resolution (higher
gives more, smaller communities), edges are weighted by similarity, and each group's node size is its gene-set
size transformed by `--node-weight-transform` (log2 by default) so large sets do not dominate. Output is
`{prefix}_clusters.tsv` with a cluster id, size, and the `|`-joined members, plus size-distribution figures.

## Notes

- **Connected components** (the default) put two groups in the same cluster if a chain of edges links them.
  **`--community`** (Leiden) finds tighter, more modular groups and is usually what you want for a dense network.
- Raise `--resolution` for more and smaller clusters; lower it for fewer and larger.
- The cutoff controls how the graph is built; a higher `-c` keeps only strong edges, giving smaller, tighter
  clusters.

## Next

Pull the members of a cluster with [query](/DBRetina/commands/query/) `--clusters-file`, or look up their genes
with [geneinfo](/DBRetina/commands/geneinfo/).
