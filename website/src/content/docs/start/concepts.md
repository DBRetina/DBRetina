---
title: Core concepts
description: Groups, features, similarity metrics, and the files DBRetina reads and writes.
---

A few terms and file types show up everywhere in DBRetina. This page defines them once.

## Groups and features

A **group** (also called a supergroup) is an entry in your database: a disease, a drug, a pathway. Each group is
described by a set of **features**: the genes for a disease, the targets for a drug, the members of a pathway.

DBRetina never looks at feature order or count within a group beyond membership: a group is a *set* of features.
Two groups are similar when they share many features.

:::caution
Group and feature names are lowercased when you build an index, and the pipe character `|` is reserved as an
internal separator, so it cannot appear in a name. Double quotes are stripped.
:::

## Similarity metrics

For a pair of groups with `s1` and `s2` features and `o` shared, DBRetina reports several metrics. All are scaled
to a **0–100 percentage**.

| Metric | Definition | Reads as |
| --- | --- | --- |
| `shared_features` | `o` | how many features the two groups have in common |
| `containment` | `o / min(s1, s2)` | how much of the *smaller* group sits inside the larger |
| `ochiai` | `o / sqrt(s1 · s2)` | geometric-mean overlap, size-balanced |
| `jaccard` | `o / (s1 + s2 − o)` | overlap over the combined set |
| `dice` | `2o / (s1 + s2)` | harmonic-style overlap |
| `odds_ratio` | enrichment vs. chance | how surprising the overlap is |

You choose which metric a command filters or sorts on with `-m`, and the threshold with `-c`. Containment is a
good default for "is this group contained in that one"; ochiai and jaccard are good for symmetric similarity.

Optionally, `pairwise --pvalue` adds a hypergeometric **p-value** testing whether an overlap is larger than
chance would predict. Unlike the similarity metrics, a p-value is *lower is better*, and a p-value cutoff keeps
the pairs at or below it.

## The index (`.dbri`)

`DBRetina index` turns your input into three files:

- **`<name>.dbri`**: the binary index every other command reads.
- **`<name>_raw.json`**: your groups and their features in plain text, for you to inspect.
- **`<name>_hashes.json`**: the same, with features hashed, used internally.

Build the index once, then reuse it across `pairwise`, `geneinfo`, `merge`, and the lookups.

## The pairwise output

`DBRetina pairwise` writes its results to a directory, **`<index>_DBRetina_pairwise/`**, stored as Parquet. This
is the default and it is what the analysis commands (`query`, `cluster`, `export`, `dedup`, `modularity`,
`graph`, `bipartite`, and the lookups) expect for their `-p`/`-d` argument.

If you want a single flat table instead, for example to open in a spreadsheet, add `--legacy-output` and
`pairwise` will also write the classic `.tsv` and `.dbrp` files.

## Chaining commands

Most work is a short pipeline: build an **index**, run **pairwise**, then feed the pairwise output into whichever
analysis you need. The [Command chaining](/DBRetina/reference/chaining/) reference shows how the pieces fit
together.
