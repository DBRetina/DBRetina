# Interactome

```
Usage: DBRetina interactome [OPTIONS]

  Construct a features-interactome.

  Detailed description:

  For a groups pairwise file, construct an interactome between the features of
  each group and all other features in the pairwise file.

Options:
  -i, --index-prefix TEXT   index file prefix  [required]
  -p, --pairwise PATH       pairwise TSV file  [required]
  --graphml                 export interactome as graphml file
  --gexf                    export interactome as gexf file
  -o, --output-prefix TEXT  output file prefix  [required]
  --help                    Show this message and exit.
```

## Command arguments

<span class="cmd"> -i, --index-prefix TEXT   Index file prefix  [required] </span>

This is the user-defined prefix that was used in the indexing step as an output prefix.

<span class="cmd">  -p, --pairwise PATH  the pairwise TSV file  [required] </span>

The original or a filtered pairwise TSV file.

<span class="cmd"> --graphml </span>

Export interactome as graphml file.

<span class="cmd"> --gexf </span>

Export interactome as gexf file.

<span class="cmd"> -o, --output TEXT    output prefix  [required] </span>

The output files prefix.


<hr class="fancy-hr">


## Output files format

<span class="cmd"> {output_prefix}_interactome.tsv </span>

A pairwise TSV file with the interactome connections.

<span class="cmd"> {output_prefix}_interactome.garphml </span>

An interactome graph in graphml format.

<span class="cmd"> {output_prefix}_interactome.gexf </span>

An interactome graph in gexf format.