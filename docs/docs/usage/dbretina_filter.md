# 3. Filter

The Filter command in DBRetina is designed to filter out the pairwise TSV file. The command requires the **full path** of the pairwise TSV file.

```

Usage: DBRetina filter [OPTIONS]

  Filter a pairwise file.

  Detailed description:

      Filter a pairwise file by distance cutoff and/or a set of groups
      (provided as a single-column file or cluster IDs in a DBRetina cluster
      file).

  Examples:

      1- distance cutoff only              | dbretina filter -p pairwise.tsv
      -d ochiai -c 0.5 -o filtered.tsv

      2- distance cutoff and groups file   | dbretina filter -p pairwise.tsv
      -d min_cont -c 0.5 -g groups.tsv -o filtered.tsv

      3- distance cutoff and a cluster IDs | dbretina filter -p pairwise.tsv
      -d max_cont -c 0.5 --clusters-file clusters.tsv --clusters-id 8 -o
      filtered.tsv

      4- groups file only                  | dbretina filter -p pairwise.tsv
      -g groups.tsv -o filtered.tsv

      5- cluster file with cluster IDs     | dbretina filter -p pairwise.tsv
      --clusters-file clusters.tsv --clusters-id 8 -o filtered.tsv

Options:
  -p, --pairwise PATH       the pairwise TSV file  [required]
  -g, --groups-file PATH    single-column supergroups file
  --clusters-file PATH      DBRetina clusters file
  --cluster-ids TEXT        comma-separated list of cluster IDs
  -d, --dist-type TEXT      select from ['min_cont', 'avg_cont', 'max_cont',
                            'ochiai', 'jaccard']  [default: NA]
  -c, --cutoff FLOAT RANGE  filter out distances < cutoff  [default: 0.0;
                            0<=x<=100]
  -o, --output TEXT         output file prefix  [required]
  --help                    Show this message and exit.
```

## 3.1 Command arguments

### 3.1.1 Filtering by distance's cutoff

<span style="color:orange;">** -c, --cutoff FLOAT RANGE  filter out distances < cutoff  [default: 0.0; 0<=x<=100] **</span>

This will filter out all pairwise distances that are below the cutoff value.

<span style="color:orange;">** -d, --dist-type TEXT      select from ['min_cont', 'avg_cont', 'max_cont', 'ochiai', 'jaccard']  [default: NA] **</span>

The distance metric to apply the cutoff on. 

### 3.1.2 Filtering by supergroups

<span style="color:orange;">** -g, --groups-file PATH    single-column supergroups file **</span>

This will filter out all pairwise distances that are between supergroups that are not in the provided groups file. The groups file is a single-column file that contains the names of the supergroups to be included in the filtering.

### 3.1.3 Filtering by clusters

This will filter out all pairwise distances that are between clusters that are not in the provided clusters file.

<span style="color:orange;">** --clusters-file PATH      DBRetina clusters file **</span>

The clusters file is a DBRetina clusters file that contains the cluster IDs to be included in the filtering.

<span style="color:orange;">** --cluster-ids TEXT        comma-separated list of cluster IDs **</span>

The cluster IDs selected from the clusters file. This argument is only used if the clusters file is not provided.

---

## 3.2 Output files format

<span style="color:orange;">** {output_prefix}.tsv **</span>

Filtered version of the pairwise TSV file.




