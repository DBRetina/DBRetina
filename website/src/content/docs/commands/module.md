---
title: module
description: Extract a group's characteristic module, the groups most similar to it.
---

Pull out a group's characteristic module: the groups most similar to it, as a plain list you can feed to other
commands. It is `neighbors` reduced to just the names, with a `--min-shared` filter to drop the noisy small-set
matches.

## When to use it

When you want a group's neighbourhood as a reusable set, for example to hand to `geneinfo`, `query -g`, or
`enrich --module`.

## Synopsis

```bash
DBRetina module -d <pairwise> GROUP -m <metric> -c <cutoff> [--min-shared K] [-n N] [-o <file>]
```

## Inputs

| Option | Meaning | Required |
| --- | --- | --- |
| `-d, --data` | the pairwise output (Parquet dir or `.tsv`) | yes |
| `GROUP` | the seed group | yes |
| `-m, --metric` | metric to rank neighbours by | yes |
| `-c, --cutoff` | minimum metric value (0–100) | yes |
| `--min-shared` | drop neighbours sharing fewer than K features | no (default 0) |
| `-n, --top` | keep only the top N neighbours (0 = all) | no |
| `-o, --output` | output file (default: stdout) | no |

## Outputs

One group name per line: the seed's neighbours passing the cutoff and `--min-shared`.

## Example

```bash
DBRetina module -d diseases_DBRetina_pairwise dis_ra -m ochiai -c 15 --min-shared 20 -n 6
```

```text
dis_lupus
dis_t1d
dis_atherosclerosis
dis_ms
dis_breast_cancer
dis_uc
```

Rheumatoid arthritis's nearest diseases are the other autoimmune ones.

## How it works

`module` is a direct-neighbour query, not a graph traversal. It returns every group joined to the seed by a
single edge that clears the cutoff (`-m`/`-c`) and shares at least `--min-shared` features, sorted by that
similarity (strongest first, or smallest p first for `-m pvalue`) and optionally cut to the top `-n`. It is
exactly one hop: every group in the output is similar to the seed itself, not to the seed's neighbours. The result
is a plain list of group names, one per line, ready to feed into `enrich --module`, `query -g`, or `geneinfo -g`.

## Notes

`--min-shared` matters most with containment, where a small gene set trivially reaches 100% and would otherwise
dominate the top of the list. Requiring a minimum number of shared features keeps the module substantial.

## Next

Rank the wider network against this module with [enrich](/DBRetina/commands/enrich/), or list the shared genes of
its members with [geneinfo](/DBRetina/commands/geneinfo/).
