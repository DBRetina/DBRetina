# Export

The Export command converts the pairwise TSV file to a dissimilarity matrix and (optionally) a newick-format file.

```
Usage: DBRetina export [OPTIONS]

  Export to distance matrix and tree formats.

  Export a pairwise TSV file into a distance matrix, newick-format file and
  circular dendrogram.

Options:
  -p, --pairwise PATH  pairwise TSV file  [required]
  -m, --metric TEXT    select from ['containment', 'ochiai', 'jaccard',
                       'pvalue']  [required]
  --newick             Convert the distance matrix to newick tree format
  -l, --labels TEXT    select from ['ids', 'names']  [default: names]
  --linkage TEXT       select from ['single', 'complete', 'average',
                       'weighted', 'centroid', 'median', 'ward']  [default:
                       ward]
  -o, --output TEXT    output prefix  [required]
  --help               Show this message and exit.
```

## Command arguments


<span class="cmd"> -p, --pairwise PATH  pairwise TSV file  [required] </span>

The original or a filtered pairwise TSV file.

<span class="cmd"> -m, --metric TEXT    select from ['containment', 'ochiai', 'jaccard', 'pvalue']  [required] </span>

The similarity metric to be utilized in the dissimilarity matrix.

<span class="cmd"> --newick              Convert the dissimilarity matrix to newick tree format </span>

This will convert the dissimilarity matrix to a newick tree format.

<span class="cmd"> -l, --labels TEXT    select from ['ids', 'names']  [default: names] </span>

The labels to be used in all outputs (plots, newick, dissimilarity matrix). Default is names. The option option is the IDs which are the IDs found in the [`{index_prefix}.namesMap`](dbretina_index.md){:target="_blank"} file.

<span class="cmd"> --linkage TEXT       select from ['single', 'complete', 'average', 'weighted', 'centroid', 'median', 'ward']  [default: ward] </span>

<span class="cmd"> -o, --output TEXT    output prefix  [required] </span>

The output files prefix.


---


## Output files format

<span class="cmd"> {output_prefix}_distmat.tsv </span>

A labels tab-separated similarity matrix file.

<span class="cmd"> {output_prefix}_distmat.pkl </span>

A labels pickle similarity matrix file in Python pickle format.

<span class="cmd"> {output_prefix}.newick </span>

A newick tree file that can be used in tree visualization tools.

<span class="cmd"> {output_prefix}_dendrogram.png </span>

A circular dendrogram plot of the tree.

<span class="cmd"> {output_prefix}_heatmap.png </span>

A heatmap plot of the similarity matrix.