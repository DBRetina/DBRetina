# Indexing (Index the input data files)

The Indexing process in DBRetina primarily focuses on creating an index structure for the input entities and their associated features.
This structure is utilized for calculating pairwise distances between input entities using the "pairwise" command, as well as querying one or more features to determine their associated entities with the "query" command.

```
Usage: DBRetina index [OPTIONS]

  Index the input data files.

Options:
  -a, --asc TEXT     associations file col1: gene_set, col2: single gene. 1st
                     line is header.
  -g, --gmt TEXT     GMT file(s)
  -o, --output TEXT  output file prefix  [required]
  --help             Show this message and exit.
```

## Command arguments

!!! warning

    - The text of all input files is automatically converted to lowercase.
    - All double quotes are removed.
    
!!! danger

    Pipe character `|` can't be used in the input data.

<span class="cmd"> -a, --asc TEXT     associations file(s) col1: gene_set, col2: single gene. 1st line is header.</span>

The "Association File" is a two-column TSV (tab-separated values) file with an included header. The first column denotes groups, while the second column indicates the features associated with each respective group. Each row signifies a single feature and its corresponding group.

**Example of an Association File with features as features:**

```tsv
Disease       feature
Breast Cancer BRCA1
Breast Cancer BRCA2
Lung Cancer   EGFR
Lung Cancer   KRAS
```
<br>

<span class="cmd">-g, --gmt TEXT     GMT file(s).</span>

The "GMT File" is a tab-delimited headerless file that contains gene sets. Each row denotes a single gene set, while the first column indicates the name of the gene set. The second column contains a description of the gene set, while the remaining columns contain the genes that belong to the gene set. <a target=target="_blank" href="https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29">Click here to read more about GMT format.</a>

<br>

<span class="cmd">-o, --output TEXT  output file prefix  [required]</span>

The output file prefix is used to name the output files. The output files are explained in detail in the next section.


<br>

!!! example

    === "Multiple GMT files"

        ```bash
        DBRetina index -g gmt_file1.gmt -g gmt_file2.gmt -o idx_example
        ```

    === "Multiple association files"

        ```bash
        DBRetina index -a asc_file1.tsv -a asc_file2.tsv -o idx_example
        ```

<br>

!!! warning
      <b>You can't use combination of gmt files and association files, the command accepts only one type of input.</b>


---

## Output files format

The output featurerated by the Index command consists of two JSON files (private and public) and a set of index files. These files are explained in detail below:

### Primary output files

<span class="cmd"> {prefix}_raw.json </span>

This JSON file contains supergroups and their related features in plain text. This JSON file is prepared for the user to understand the final input data that is used for indexing.

<span class="cmd"> {prefix}_hashes.json </span>

This is another JSON file that contains supergroups and their related features, but the features are hashed for indexing. This JSON file is used internally by DBRetina for indexing.

??? info end "Advanced Output (For developers)"
    `{output_prefix}_groupID_to_featureCount.bin`: binary file that contains the number of features for each supergroup.

    `{output_prefix}_groupID_to_featureCount.tsv`: tab-delimited file that contains the number of features for each supergroup.

    `{output_prefix}_color_count.bin`: binary file that contains the colors count.

    `{output_prefix}.phmap`: The parallel-hash-map binary index file contains features and their associated colors.

    `{output_prefix}.namesMap`: pipe-delimited file of supergroup ID and its original name with total number of supergroups in header.

    `{output_prefix}_color_to_sources.bin`: binary file that contains the colors and their associated supergroup IDs.

    `{output_prefix}.extra`: text file that contains metadata for inter-command communication.
