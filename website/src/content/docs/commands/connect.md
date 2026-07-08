---
title: connect
description: Find the shortest path between two groups through the similarity network.
---

Find how two groups are linked through the network, even when they do not overlap directly. `connect` returns the
shortest path between them and the metric at each step, so you can see the intermediaries that connect them (for
example a disease reaching a drug through a shared pathway).

## When to use it

After `pairwise`, when two groups have little or no direct similarity but you want to know whether, and how, the
network connects them.

## Synopsis

```bash
DBRetina connect -d <pairwise> GROUP_A GROUP_B -m <metric> -c <cutoff> [-o <file>]
```

## Inputs

| Option | Meaning | Required |
| --- | --- | --- |
| `-d, --data` | the pairwise output (Parquet dir or `.tsv`) | yes |
| `GROUP_A GROUP_B` | the two group names to connect | yes |
| `-m, --metric` | metric used for the edges | yes |
| `-c, --cutoff` | minimum metric value for an edge (0–100) | yes |
| `-o, --output` | output file (default: stdout) | no |

## Outputs

One row per step: the step number, the group at that step, and the edge metric leading into it.

## Example

```bash
DBRetina connect -d diseases_DBRetina_pairwise dis_ra dis_atherosclerosis -m containment -c 40
```

```text
step  group              containment
0     dis_ra             -
1     dis_atherosclerosis 63.6
```

## How it works

`connect` treats the pairwise table as a graph: every pair of groups whose similarity clears the cutoff
(`-m`/`-c`) becomes an undirected edge, and every group is a vertex (unconnected groups included). It then finds
the shortest path between the two named groups by **hop count**, a breadth-first search for the fewest
intermediate groups that chain them together, each link a real similarity above the cutoff. The search is
unweighted: it minimises the number of steps, not the total similarity along the way, so the metric value shown
for each hop is reported, not optimised. If no chain of above-cutoff links joins them, the groups are reported as
not connected at that threshold. Unlike the neighbourhood commands, `connect` takes `-m`/`-c` only, not
`--min-shared`.

## Notes

The path is the shortest by number of hops (an unweighted shortest path), so it minimises intermediaries rather
than maximising edge strength. `-c` is required: without a cutoff the graph would include every edge, which on a
large network is both slow and uninformative. If the two groups are not connected at the cutoff, `connect` says
so and stops cleanly.

## Next

Rank the groups near one of them with [module](/DBRetina/commands/module/), or find related groups by shared
neighbourhood with [enrich](/DBRetina/commands/enrich/).
