---
title: graph
description: Export the similarity network as node and edge files for graph tools.
---

Write the similarity network as plain node and edge tables you can load into Cytoscape, Gephi, or any graph
library.

## When to use it

When you want to visualise or analyse the network in an external graph tool rather than inside DBRetina.

## Synopsis

```bash
DBRetina graph -i <index> -p <pairwise> -m <metric> -c <cutoff> [--intra-targets <f>] [--inter-targets <f>] -o <prefix>
```

## Inputs

| Option | Meaning | Required |
| --- | --- | --- |
| `-i, --index-prefix` | the index the pairwise came from | yes |
| `-p, --pairwise` | the pairwise output (Parquet dir, `.dbrp`, or `.tsv`) | yes |
| `-m, --metric` | metric used for the edges | yes |
| `-c, --cutoff` | minimum metric value for an edge (0–100) | no (default `0`) |
| `--intra-targets` | comma-separated TSVs; keep edges *within* these group sets | no |
| `--inter-targets` | comma-separated TSVs; keep edges *between* these group sets | no |
| `--include-isolates` | keep nodes that have no edges | no |
| `--visualize` | open an interactive view (needs the optional `dash` package) | no |
| `-o, --output` | output prefix | yes |

## Outputs

| File | What it is |
| --- | --- |
| `<prefix>_edges.tsv` | the edges and their metric values |
| `<prefix>_nodes.tsv` | the nodes |

## Example

```bash
DBRetina graph -i kegg -p kegg_DBRetina_pairwise -m ochiai -c 30 -o kegg_graph
```

```text
[INFO] Total number of edges: 368
[INFO] total number of nodes: 122
[SUCCESS] Done!
```

## Notes

Targets are optional. With neither `--intra-targets` nor `--inter-targets`, `graph` exports the whole network
above the cutoff. Use the target options to restrict the export to edges within or between particular sets of
groups.

## Next

Load the node and edge files into your graph tool, or cluster the same network with
[cluster](/DBRetina/commands/cluster/).
