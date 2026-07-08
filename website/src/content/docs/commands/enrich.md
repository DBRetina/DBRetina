---
title: enrich
description: Rank groups by enrichment of a reference module among their neighbours.
---

Rank every group by how much its neighbourhood overlaps a reference module. Given a module (a list of groups),
`enrich` scores each group by the hypergeometric enrichment of module members among its neighbours, so groups
that sit inside the module's neighbourhood rank first. This is a guilt-by-association ranking.

## When to use it

After building a module (with [module](/DBRetina/commands/module/)), to find the groups across the whole network
most associated with it, for example the diseases or drugs closest to a disease's mechanism.

## Synopsis

```bash
DBRetina enrich -d <pairwise> --module <groups.txt> -m <metric> -c <cutoff> [--min-shared K] [--exclude-module] [-n N] [-o <file>]
```

## Inputs

| Option | Meaning | Required |
| --- | --- | --- |
| `-d, --data` | the pairwise output (Parquet dir or `.tsv`) | yes |
| `--module` | single-column file of reference group names | yes |
| `-m, --metric` | metric used for the neighbour edges | yes |
| `-c, --cutoff` | minimum metric value for an edge (0–100) | yes |
| `--min-shared` | drop neighbour edges sharing fewer than K features | no |
| `--exclude-module` | drop the module's own members, keep only new candidates | no |
| `-n, --top` | number of top groups to show (0 = all) | no |
| `-o, --output` | output file (default: stdout) | no |

## Outputs

`group`, `hits` (module neighbours), `neighbors` (total), `enrichment_p`, `in_module`.

## Example

```bash
DBRetina enrich -d diseases_DBRetina_pairwise --module ra_module.txt -m ochiai -c 20 --exclude-module
```

## How it works

`enrich` runs one hypergeometric test per group, scoring how enriched that group's neighbourhood is for the
module.

**The neighbourhood.** Two groups are neighbours when their similarity on the chosen metric clears the cutoff
(`-m`/`-c`) and they share at least `--min-shared` features. That turns the pairwise table into a graph. For a
group *g*, its neighbourhood is the set of groups joined to it by such an edge, and its degree *d* is how many
there are.

**The test.** Let *N* be the number of groups in the pairwise (the universe), *M* the module size, and *k* the
number of *g*'s neighbours that are module members. Under the null that *g*'s *d* neighbours are a random draw,
without replacement, from the *N* groups (of which *M* are in the module), *k* follows a hypergeometric
distribution, and `enrich` reports the upper-tail probability of seeing at least *k* module members among them:

```text
enrichment_p  =  P(X >= k),   X ~ Hypergeometric(population = N, successes = M, draws = d)
```

A small *p* means *g* has more module members among its neighbours than a random group of the same degree would,
which is the guilt-by-association signal. `hits` in the output is *k*, `neighbors` is *d*.

**The ranking.** Groups are sorted by `enrichment_p` ascending, ties broken by more hits. `--exclude-module` drops
the module's own members so only new candidates remain; `-n` limits the output to the top *N*.

**What the null does and does not assume.** Because the draw is over all *N* groups, the test corrects for a
group's degree: a highly connected group needs proportionally more module neighbours to score. It does *not*
correct for module members being unusually well or poorly connected, so a very broad module can pull in large,
hub-like groups. The neighbourhood is a modelling choice: a low cutoff makes the graph dense and the test
permissive, while a higher cutoff or `--min-shared` keeps neighbourhoods specific. And because a separate test is
run for every group, `enrichment_p` is **uncorrected** for multiple comparisons: use it to rank, and validate a
shortlist with a permutation test (as in the [cancer](/DBRetina/examples/enrich-cancer-hallmarks/) and
[autoimmune](/DBRetina/examples/autoimmune-diseases/) case studies) rather than reading it as a calibrated
significance.

## Next

Trace how a top-ranked group connects back to the seed with [connect](/DBRetina/commands/connect/).
