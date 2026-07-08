---
title: modularity
description: Score how modular each group is within the similarity network.
---

Give each group a modularity score describing how it sits in the similarity network, whether its overlap is
concentrated in a few strong partners or spread thinly across many.

## When to use it

When you want a per-group number summarising its connectivity, rather than a clustering of the whole network.

## Synopsis

```bash
DBRetina modularity -i <index> -p <pairwise> -c <cutoff> -o <prefix>
```

## Inputs

| Option | Meaning | Required |
| --- | --- | --- |
| `-i, --index-prefix` | the index the pairwise came from | yes |
| `-p, --pairwise` | the pairwise output (Parquet dir, `.dbrp`, or `.tsv`) | yes |
| `-c, --cutoff` | containment cutoff defining an edge (0–100) | yes |
| `-o, --output` | output prefix | yes |

## Outputs

`<prefix>_modularity.tsv`, one row per group:

```text
gene_set                                            fragmentation  heterogeneity  modularity
reactome_opioid_signalling                          -13            55             42
reactome_g_beta_gamma_signalling_through_pi3kgamma  -53            6              47
reactome_nonhomologous_end_joining_nhej             -35            11             24
pid_dna_pk_pathway                                  -4             3              1
```

Modularity is the magnitude of `fragmentation + heterogeneity`.

## How it works

`modularity` is not graph (Newman) modularity; it is a per-group **containment-directionality** score. Working
only on the containment metric (there is no `-m`), it takes every pair of groups whose containment clears `-c` and
asks which group is the larger. For each such pair the smaller, contained group gets a point of *fragmentation* and
the larger, containing group a point of *heterogeneity*. A group's score is then the absolute difference,
`|heterogeneity - fragmentation|`: a group that is almost always the small subset of bigger groups, or almost
always the big set absorbing smaller ones, scores high, while a group that plays both roles equally scores near
zero. It flags the lopsided groups, the redundant subsets and the broad umbrellas, which are the ones worth
pruning with `dedup` or `setcov`. Output is `{prefix}_modularity.tsv` with each group's `fragmentation`,
`heterogeneity`, and `modularity`.

## Notes

`modularity` handles extreme cutoffs gracefully: if the cutoff is so high that no pairs qualify, every group is
reported with a modularity of 0 rather than erroring.

## Next

Cluster the same network with [cluster](/DBRetina/commands/cluster/), or reduce it to a covering subset with
[setcov](/DBRetina/commands/setcov/).
