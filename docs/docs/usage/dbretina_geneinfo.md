# Querying a DBRetina index

Extract gene information from a DBRetina index with a set of groups (provided as a single-column file or cluster IDs in a DBRetina cluster file).

```
Usage: DBRetina genenet [OPTIONS]

  Construct a Gene Network.

  Detailed description:

  For a groups pairwise file, construct an gene network between the features of
  each group and all other features in the pairwise file.

Options:
  -i, --index-prefix TEXT   index file prefix  [required]
  -p, --pairwise PATH       pairwise TSV file  [required]
  --graphml                 export genenet as graphml file
  --gexf                    export genenet as gexf file
  -o, --output-prefix TEXT  output file prefix  [required]
  --help                    Show this message and exit.
```

## Command arguments


<span class="cmd"> -i, --index-prefix TEXT  Index file prefix  [required] </span>

This is the user-defined prefix that was used in the indexing step.

<span class="cmd"> -g, --groups-file PATH    single-column supergroups file </span>

This will filter out all pairwise similarities that are between supergroups that are not in the provided groups file. The groups file is a single-column file that contains the names of the supergroups to be included in the filtering.

<span class="cmd"> --cluster-ids TEXT        comma-separated list of cluster IDs </span>

The cluster IDs selected from the clusters file. This argument is only used if the clusters file is not provided.

<span class="cmd"> -o, --output TEXT        output file prefix  [required] </span>

The output prefix.

---

## Output files format


<!-- TODO: Implement later -->
<!-- <span class="cmd"> {output_prefix}_features_count_per_group.tsv </span>

A TSV file containing two columns, the first column is the supergroup name and the second column is the number of features that are contained in that supergroup, and the third column is a PIPE-separated list of supergroups from the user geneinfo query that are associated with the feature in the first column. -->

<span class="cmd"> {output_prefix}_feature_to_groups.tsv </span>

A three-columns TSV file containing the feature name, the number of supergroups associated with that feature, and a PIPE-separated list of supergroups associated with that feature.

<!-- TODO: Fix later -->
<!-- <span class="cmd"> {output_prefix}_features_count_per_group_histogram.png </span>

A histogram plot showing the distribution of the number of features per supergroup. -->