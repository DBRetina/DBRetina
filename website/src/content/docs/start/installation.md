---
title: Installation
description: Install DBRetina with pip, or build it from source in a conda environment.
---

DBRetina runs on `x86_64` Linux and requires **Python 3.11 or newer** (3.12 recommended). Prebuilt wheels are
published for CPython 3.12 and 3.13, so most people can install it in one line.

## With pip

```bash
pip install DBRetina
```

That pulls a prebuilt wheel, no compiler needed. Check it worked:

```bash
DBRetina --help
```

If you prefer an isolated environment:

```bash
conda create -n dbretina python=3.12
conda activate dbretina
pip install DBRetina
```

## Building from source

The commands above install the published wheel. Build from source only if you are developing DBRetina or need a
platform without a wheel. The C++ core links against Arrow and Parquet, so the build needs a compiler toolchain
and those libraries. The reproducible way to get them is the bundled conda environment:

```bash
git clone https://github.com/DBRetina/DBRetina
cd DBRetina
conda env create -f environment.yml
conda activate dbretina
pip install '.[all]'
```

The `[all]` extra installs the optional feature groups (statistics helpers and the export integrations) alongside
the core.

:::note
Building from source without conda's Arrow/Parquet will fail at `find_package(Arrow REQUIRED)`. If you manage
those libraries yourself, point CMake at them with `export CMAKE_PREFIX_PATH="$CONDA_PREFIX"` before
`pip install`.
:::

## What you need as input

DBRetina reads your data as either of two plain-text formats, a **GMT** file or a two-column **association
(ASC)** file. Both are covered on the [index](/DBRetina/commands/) page. If you just want to try the tool, the
[Quickstart](/DBRetina/start/quickstart/) ships a tiny example file so you don't have to prepare anything.
