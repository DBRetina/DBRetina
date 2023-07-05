# 1. Indexing (Index the input data files)

The Indexing process in DBRetina primarily focuses on featurerating an index structure for the input entities and their associated features. This structure is utilized for calculating pairwise distances between input entities using the "pairwise" command, as well as querying one or more features to determine their associated entities with the "query" command.

```
Usage: DBRetina index [OPTIONS]

  Index the input data files.

Options:
  -a, --asc PATH     associations file  [required]
  -n, --names PATH   names file
  -o, --output TEXT  output file prefix  [required]
  --help             Show this message and exit.
```

## 1.1 Command arguments

!!! warning

    - The text of all input files is automatically converted to lowercase.
    - All double quotes are removed.
    
!!! danger

    Pipe character `|` can't be used in the input data.

<span style="color:orange;">**-a, --asc PATH     associations file  [required]**</span>

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


<span style="color:orange;">** -n, --names PATH      names file [optional]**</span>

The "Names File" is another two-column TSV file ^^==(without header)==^^. The first column represents the group names from the Association File, and the second column specifies the corresponding supergroup or alias name. The Names File serves to allow the grouping of related groups under a broader category or supergroup. This file is optional; if not provided, the groups from the Association File will be treated as supergroups. Also, it can also provide partial grouping of the association file.

**Example of a Names File:**

```tsv
Breast Cancer Cancer
Lung Cancer   Cancer
```

<br>

<span style="color:orange;">** -n, --names PATH      names file [optional]**</span>

Out output files of the index will carry the provided index. 
XX add somthing about feature to groups then refer to the output file format section

---


## 1.2 Output files format

The output featurerated by the Index command consists of two JSON files (private and public) and a set of index files. These files are explained in detail below:


<span style="color:orange;">** 1.2.1 {prefix}_raw.json **</span>

This JSON file contains supergroups and their related features in plain text.

<span style="color:orange;">** 1.2.2 {prefix}_hashes.json **</span>

This is another JSON file that contains supergroups and their related features, but the features are hashed for indexing.

<span style="color:orange;">** 1.2.3 {prefix}_* **</span>

The generated index files contains binary and text files that hold the index information.