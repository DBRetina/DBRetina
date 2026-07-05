# DBRetina installation

## System requirements

The current version of DBRetina targets `x86_64` Linux and requires **Python 3.11 or newer** (3.12 recommended; prebuilt wheels are published for CPython 3.12 and 3.13). Docker support is planned.

## Installation

To install DBRetina, type in the terminal:

```sh
pip install DBRetina
```

Alternatively, if you're working in a conda environment, you can type:

```sh
conda create -n dbretina python=3.12
conda activate dbretina
pip install DBRetina
```

### Building from source

The commands above install the published wheel. To build from source you need the C++ toolchain and Arrow/Parquet; create the environment from `environment.yml`:

```sh
conda env create -f environment.yml
conda activate dbretina
pip install '.[all]'
```
