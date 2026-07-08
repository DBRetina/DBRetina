---
title: index
description: Build a DBRetina index (.dbri) from your gene sets.
---

Build an index from your groups and their features. This is the first step in every workflow. Every command that
follows reads the index it produces.

## When to use it

Run `index` once per dataset. After that, reuse the `.dbri` it writes for `pairwise`, `geneinfo`, `merge`, and
the lookup commands.

## Synopsis

```bash
DBRetina index (-g <file.gmt> | -a <file.tsv>) -o <prefix>
```

## Inputs

| Option | Meaning | Required |
| --- | --- | --- |
| `-g, --gmt` | GMT file(s). Repeat `-g` for several files. | one of `-g`/`-a` |
| `-a, --asc` | Association TSV file(s): two columns, `group` and `feature`, with a header. Repeat `-a` for several. | one of `-g`/`-a` |
| `-o, --output` | Output prefix for the index files. | yes |
| `--no-progress` | Turn off the progress bar. | no |
| `--debug` | Print per-phase diagnostics. | no |

You give **either** GMT files **or** association files, not both in one call.

### GMT format

A tab-delimited, headerless file. Each row is one group: the first column is the group name, the second is a
description (often unused), and the rest are its features.

```text
breast_cancer  na  BRCA1  BRCA2  TP53  PTEN
```

### Association (ASC) format

A two-column TSV **with** a header. Each row is one feature and the group it belongs to.

```text
group          feature
breast_cancer  BRCA1
breast_cancer  BRCA2
```

## Outputs

| File | What it is |
| --- | --- |
| `<prefix>.dbri` | the binary index every command reads |
| `<prefix>_raw.json` | your groups and features in plain text, for inspection |
| `<prefix>_hashes.json` | the same with hashed features, used internally |

## Example

```bash
DBRetina index -g example.gmt -o example
```

```text
[SUCCESS] File(s) have been indexed.
[SUCCESS] Index written to example.dbri
✓ index done · total 0.0s · peak memory 0.26 GiB
```

Index several GMT files into one database by repeating `-g`. This is how you combine sources:

```bash
DBRetina index -g kegg.gmt -g reactome.gmt -g wikipathways.gmt -o pathways
```

:::caution
Group and feature names are lowercased, double quotes are stripped, and the pipe character `|` is not allowed
(it is a reserved separator). GMT input must be plain text; gunzip a `.gmt.gz` file first.
:::

## How it works

`index` reads your GMT (or association file), lowercases every group and gene name, and hashes each gene to a
64-bit FNV-1a value, so groups are compared by hashed gene identity rather than by string. Genes are deduplicated
within a group, so a group's size is its count of distinct genes, and a gene that appears in several groups is
stored once as a shared "colour" pointing at all of them. That colour index is what lets `pairwise` count shared
genes without re-reading gene lists. The index also records the feature universe, the number of distinct genes
across the whole database, which the size-aware metrics and the p-value later use as their population. Everything
is written to a single binary `.dbri`. Two rules matter: the `|` character is reserved as an internal delimiter and
rejected in names, and groups that share a name are merged into one (with a warning).

## Next

Compare every pair with [pairwise](/DBRetina/commands/pairwise/), or grow an index with
[merge](/DBRetina/commands/merge/).
