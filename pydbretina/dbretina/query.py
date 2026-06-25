#!/usr/bin/python
# -*- coding: utf-8 -*-
from __future__ import division
import sys
import _dbretina_internal as dbretina_internal
import click
from dbretina.click_context import cli
from dbretina.validators import validate_metric
import subprocess
import os
import dbretina.dbretina_doc_url as dbretina_doc


# Metrics where a SMALLER value is better. p-value is the only one: a smaller p
# is more significant, so a p-value cutoff keeps pairs with value <= cutoff,
# whereas every similarity metric (containment/ochiai/jaccard/csi/dice/
# odds_ratio) is "higher is better" and keeps value >= cutoff. Filtering p-value
# with >= inverted the result -- it returned the LEAST significant pairs and
# silently dropped the significant ones (issue 068). This mirrors
# pairwise_store.LOWER_IS_BETTER and metric_passes_cutoff() in
# src/DBRetinaPairwise.cpp; the three definitions must stay in sync. It is kept
# local here so the awk-only query path needs no pyarrow/duckdb import.
LOWER_IS_BETTER = frozenset({"pvalue"})


def cutoff_operator(metric):
    """Comparison operator for a metric's cutoff: ``<=`` for p-value (lower is
    better), ``>=`` for the similarity metrics (higher is better)."""
    return "<=" if metric in LOWER_IS_BETTER else ">="


def _read_group_name_set(groups_file):
    """Read a single-column groups file into a lowercased name set.

    Mirrors the awk ``BEGIN`` block exactly:
    ``gsub(/"/, "", $1); id_map[tolower($1)]=1`` -- take the first tab field,
    strip ALL double quotes, lowercase. Blank lines are ignored.
    """
    names = set()
    with open(groups_file) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            first = line.split("\t", 1)[0].replace('"', "")
            if first:
                names.add(first.lower())
    return names


def _format_store_row(row, gid1, gid2, names_map, has_pvalue):
    """Format one store DataFrame row to the canonical query output line.

    Matches the cutoff-only store branch byte-for-byte: ids, names from the
    names map, integer shared_features, then the six similarity metrics at
    ``%.1f`` and (when present) the raw-string pvalue. Names are joined via the
    store names map (the parquet rows carry only ids), exactly like the awk path
    which prints the source TSV's own name columns ($3/$4)."""
    fields = [str(gid1), str(gid2),
              names_map.get(gid1, ""), names_map.get(gid2, ""),
              str(int(row['shared_features'])),
              f"{row['containment']:.1f}", f"{row['ochiai']:.1f}",
              f"{row['jaccard']:.1f}", f"{row['csi']:.1f}",
              f"{row['dice']:.1f}", f"{row['odds_ratio']:.1f}"]
    if has_pvalue:
        fields.append(str(row['pvalue']))
    return '\t'.join(fields)


def is_awk_available():
    try:
        subprocess.run(["awk"], stdin=subprocess.DEVNULL,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        return True
    except FileNotFoundError:
        return False


def execute_bash_command(command):
    try:
        process = subprocess.Popen(
            command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, text=True
        )
        output, error = process.communicate()

        if process.returncode == 0:
            return True
        print("Command execution failed.", file=sys.stderr)
        print(f"Error:\n{error}", file=sys.stderr)
        return False

    except Exception as e:
        print(f"An error occurred: {e}", file=sys.stderr)
        return False


def get_command():
    _sys_argv = sys.argv
    for i in range(len(_sys_argv)):
        if os.path.isfile(_sys_argv[i]):
            _sys_argv[i] = os.path.abspath(_sys_argv[i])
    return "DBRetina " + " ".join(_sys_argv[1:])


def increment_version(output):
    version = 1
    base_output = output.rsplit('.', 1)[0]
    if not os.path.isfile(output):
        return output
    while True:
        output_version = f"{base_output}_{version}.tsv"
        if not os.path.isfile(output_version):
            return output_version
        version += 1


def validate_numbers(ctx, param, value):
    if not len(value):
        return []
    try:
        return [int(num) for num in value.split(',')]
    except ValueError as e:
        raise click.BadParameter(
            'Numbers must be a comma-separated list of integers') from e


def path_to_absolute_path(ctx, param, value):
    return value if value == "NA" else os.path.abspath(value)

def check_cutoff_value(ctx, param, value):
    # if value not == -1 and not between 0 to 100
    if value != -1 and (value < 0 or value > 100):
        raise click.BadParameter(
            'cutoff must be between 0 and 100 or -1 for no cutoff.')
    else:
        return value
        


def _classify_pairwise_input(pairwise_file):
    """Classify a query -p input as 'tsv', 'store', or 'dbri'.

    -p accepts the pairwise TSV, but pairwise emits several sibling forms next
    to it (a parquet directory, a .dbrp binary, and the .dbri index). Without
    this guard the text header-copy below open()s the path unconditionally and
    crashes with a raw traceback on the non-TSV forms (issue 020):
      - parquet directory -> IsADirectoryError
      - .dbrp / .dbri     -> UnicodeDecodeError

    Returns:
      'tsv'   -> a readable text TSV: use the existing awk/TSV code paths.
      'store' -> a parquet directory, or a .dbrp with a sibling parquet store:
                 query via dbretina.compat.open_pairwise (cutoff-only).
      'dbri'  -> a DBRetina index, not a pairwise file: caller emits a clean
                 error (an index can't be queried as a pairwise file).
    """
    # A .dbri index is not a pairwise file. Detect by extension or magic bytes.
    if pairwise_file.endswith(".dbri"):
        return "dbri"
    if os.path.isfile(pairwise_file):
        try:
            with open(pairwise_file, "rb") as fh:
                magic = fh.read(4)
        except OSError:
            magic = b""
        if magic == b"DBRI":
            return "dbri"

    # A directory input is a parquet pairwise dir -> route through the store.
    if os.path.isdir(pairwise_file):
        return "store"

    # A regular file: is it decodable text (the TSV) or a binary (.dbrp)?
    try:
        with open(pairwise_file, "r", encoding="utf-8") as fh:
            fh.read(4096)
        return "tsv"
    except (UnicodeDecodeError, IsADirectoryError):
        return "store"


@cli.command(name="query", epilog = dbretina_doc.doc_url("query") , help_priority=3)
@click.option('-p', '--pairwise', 'pairwise_file', callback=path_to_absolute_path, required=True, type=click.Path(exists=True), help="the pairwise parquet dir (or legacy .dbrp/.tsv)")
@click.option('-g', '--groups-file', "groups_file", callback=path_to_absolute_path, required=False, default="NA", type=click.Path(exists=False), help="single-column supergroups file")
@click.option('--clusters-file', "clusters_file", callback=path_to_absolute_path, required=False, default="NA", type=click.Path(exists=False), help="DBRetina clusters file")
@click.option('--cluster-ids', "cluster_ids", callback=validate_numbers, required=False, default="", help="comma-separated list of cluster IDs")
@click.option('-m', '--metric', "metric", required=False, default="NA", type=click.STRING, callback=validate_metric, help="select from ['containment', 'ochiai', 'jaccard', 'pvalue']")
@click.option('-c', '--cutoff', callback=check_cutoff_value, required=False, default=-1, type=click.FLOAT, help="filter out similarities < cutoff (for -m pvalue, keeps pairs with pvalue <= cutoff)")
@click.option('--extend', "extend", is_flag=True, default=False, show_default=True, help="include all supergroups that are linked to the given supergroups.")
@click.option('-o', '--output', "output_file", required=True, type=click.STRING, help="output file prefix")
@click.pass_context
def main(ctx, pairwise_file, groups_file, metric, cutoff, output_file, clusters_file, cluster_ids, extend):
    # sourcery skip: low-code-quality
    """Query a pairwise file.

Detailed description:

    Query a pairwise file by similarity cutoff and/or a set of groups (provided as a single-column file or cluster IDs in a DBRetina cluster file).
    """
    
    # NOTE: cutoff range is validated by the check_cutoff_value callback, which
    # correctly treats the default -1 as "no cutoff". A second inline guard used
    # to live here as `if cutoff > 100 or cutoff < 0` -- it fired on the default
    # -1 and aborted groups-only / clusters-only / no-args queries with a
    # misleading "cutoff must be between 0 and 100" (issue 021). Removed so those
    # modes fall through to their proper branches / accurate errors below.

    # Extend must be used only when clusters file or groups file is provided
    if extend and groups_file == "NA" and clusters_file == "NA":
        ctx.obj.ERROR("DBRetina's query command requires a groups_file or clusters_file if --extend is provided.")


    # check if not any option is provided for filteration
    if metric == "NA" and cutoff == -1 and groups_file == "NA" and clusters_file == "NA":
        ctx.obj.ERROR(
            "DBRetina's query command requires at least one option to query the pairwise file.")

    # if clusters_file then must be cluster_id
    if clusters_file != "NA" and not len(cluster_ids):
        ctx.obj.ERROR(
            "DBRetina's query command requires cluster_id(s) if clusters_file is provided.")
    elif clusters_file == "NA" and len(cluster_ids):
        ctx.obj.ERROR(
            "DBRetina's query command requires a clusters_file if cluster_id(s) is provided.")

    # can't query by groups_file and clusters_file at the same time
    if groups_file != "NA" and clusters_file != "NA":
        ctx.obj.ERROR(
            "DBRetina's query command can't query by groups_file and clusters_file at the same time.")

    # if metric is provided then cutoff must be provided
    if metric != "NA" and cutoff == -1:
        ctx.obj.ERROR(
            "DBRetina's query command requires a cutoff if metric is provided.")
    elif metric == "NA" and cutoff != -1:
        ctx.obj.ERROR(
            "DBRetina's query command requires a metric if cutoff is provided.")

    if not is_awk_available():
        ctx.obj.ERROR(
            "DBRetina's query command requires awk to be installed and available in the PATH.")

    # Detect the -p input form. -p is the pairwise TSV, but pairwise emits sibling
    # forms (parquet dir, .dbrp, .dbri) that used to crash with a raw traceback
    # (issue 020). Route them coherently instead.
    input_kind = _classify_pairwise_input(pairwise_file)
    if input_kind == "dbri":
        ctx.obj.ERROR(
            f"'{pairwise_file}' is a DBRetina index (.dbri), not a pairwise file. "
            "Pass the pairwise TSV (or its parquet directory / .dbrp) to -p.")
    # A groups/clusters filter over a parquet store input (PLAN-094 Step-1b): the
    # awk path below needs the text TSV, but the producer also emits a parquet
    # directory, so port the SAME group/cluster/extend selection to the store
    # (see _query_store_groups). A bare pre-existing .tsv stays on the awk path.
    store_for_filter = None
    if input_kind == "store" and (groups_file != "NA" or clusters_file != "NA"):
        from dbretina.compat import open_pairwise
        try:
            store_for_filter = open_pairwise(pairwise_file)
        except Exception:
            store_for_filter = None
        if store_for_filter is None:
            ctx.obj.ERROR(
                f"'{pairwise_file}' is a pairwise directory without a usable "
                f"parquet store (no manifest.json); pass the pairwise TSV, its "
                f"parquet directory, or the .dbrp to -p.")

    metric_to_col = {
        "containment": 5,
        "ochiai": 6,
        "jaccard": 7,
        "csi": 8,
        "dice": 9,
        "odds_ratio": 10,
        "pvalue": 11,
    }
    
    # check if pvalue (format-aware: covers the .tsv / parquet-dir / .dbrp forms
    # up-front so a -m pvalue request on a pvalue-less dataset fails with the same
    # clean error on every route, instead of crashing in (or silently skipping)
    # the store/.dbrp branches below (issues 043/053).
    from dbretina.compat import pairwise_has_pvalue
    input_has_pvalue = pairwise_has_pvalue(pairwise_file)
    if metric == "pvalue" and not input_has_pvalue:
        ctx.obj.ERROR("pvalue not found in pairwise file!")

    if metric in metric_to_col:
        # +1 because awk is 1-indexed
        awk_column = metric_to_col[metric] + 1
        # p-value is "lower is better", so its cutoff keeps value <= cutoff;
        # similarity metrics keep value >= cutoff (issue 068).
        awk_op = cutoff_operator(metric)
    elif metric != "NA":
        ctx.obj.ERROR(f"DBRetina's query command doesn't support the metric {metric}.")

    # check if output_file already exist
    output_file += ".tsv"
    if os.path.exists(output_file):
        ctx.obj.WARNING(f"Output file {output_file} already exists, overwriting ...")

    if groups_file != "NA" and not os.path.exists(groups_file):
        ctx.obj.ERROR(f"Groups file {groups_file} doesn't exist.")

    # The output's first non-comment line must be the canonical column header,
    # identical across input forms. The .tsv form copies the source's full text
    # header (comments + that header row). A parquet directory / .dbrp has no text
    # header to copy, so we seed the #command line and then write the SAME column
    # header explicitly -- otherwise the store/.dbrp output was headerless while
    # the data rows matched, an inconsistency across forms (issue 047).
    _COLUMN_HEADER = (
        "group_1_ID\tgroup_2_ID\tgroup_1_name\tgroup_2_name\tshared_features\t"
        "containment\tochiai\tjaccard\tcsi\tdice\todds_ratio"
    )
    if input_kind == "tsv":
        with (open(pairwise_file) as f, open(output_file, 'w') as w):
            for line in f:
                if line.startswith("#"):
                    w.write(line)
                else:
                    w.write(f"#command: {get_command()}\n")
                    w.write(line)
                    break
    else:
        with open(output_file, 'w') as w:
            w.write(f"#command: {get_command()}\n")
            header = _COLUMN_HEADER + ("\tpvalue" if input_has_pvalue else "")
            w.write(header + "\n")

    ctx.obj.INFO(
        f"Querying the pairwise matrix on the {metric} column with a cutoff of {cutoff} and groups file {groups_file}."
    )

    _tmp_file = ".DBRetina.tmp.group"

    all_ids = set()
    if clusters_file != "NA":
        groups_file = _tmp_file
        with open(clusters_file) as f, open(_tmp_file, 'w') as W:

            # skip comments
            while True:
                pos = f.tell()
                line = f.readline()
                if not line.startswith('#'):
                    f.seek(pos)
                    break

            next(f)
            for line in f:
                line = line.strip().split('\t')
                cluster_id = int(line[0])
                if cluster_id in cluster_ids:
                    all_ids.add(cluster_id)
                    W.write(line[2].replace('|', '\n') + '\n')

    unfound_ids = set(cluster_ids).difference(all_ids)
    if len(unfound_ids):
        ctx.obj.WARNING(
            f"Couldn't find the following cluster IDs: {unfound_ids}")


    extended_ids_list = ".DBRetina_extended_ids_list"

    with open(extended_ids_list, 'w') as f:
        f.write("")

    extended_supergroups_file = f"{output_file.replace('.tsv','')}_extended_supergroups.txt"

    # ── Parquet-store group/cluster/extend filtering (PLAN-094 Step-1b) ──
    # Reproduces the awk selection below over the parquet store so the legacy
    # .tsv is no longer required for a groups/clusters query. A bare pre-existing
    # .tsv (input_kind == "tsv") still takes the awk path further down.
    if store_for_filter is not None:
        store = store_for_filter
        names_map = store.get_names_map()                 # id -> name (lowercased)
        name_to_ids = {}                                  # name -> {ids}
        for _gid, _name in names_map.items():
            name_to_ids.setdefault(_name, set()).add(int(_gid))
        has_pvalue = store.has_pvalue

        # The working id-set, from the (lowercased) group-name set. Names not in
        # the index simply never match -- drop them, exactly as awk would.
        working_names = _read_group_name_set(groups_file)
        working_ids = set()
        for nm in working_names:
            working_ids.update(name_to_ids.get(nm, ()))

        # Pull the passing pairs once (cutoff applied in SQL when present).
        cols = ["group_1_id", "group_2_id", "shared_features",
                "containment", "ochiai", "jaccard", "csi", "dice", "odds_ratio"]
        if has_pvalue:
            cols.append("pvalue")
        if cutoff != -1:
            df = store.to_pandas(metric=metric, cutoff=cutoff, columns=cols)
        else:
            df = store.to_pandas(columns=cols)
        id1 = df["group_1_id"].astype("int64")
        id2 = df["group_2_id"].astype("int64")

        if extend:
            # EITHER endpoint in the set (cutoff already applied above) -> collect
            # BOTH endpoints of those pairs as the extended id-set.
            mask = id1.isin(working_ids) | id2.isin(working_ids)
            extended_ids = set(id1[mask]).union(set(id2[mask]))
            # Write the extended supergroup NAMES (sorted unique), same file the
            # awk path produces, then the final filter uses this set.
            ext_names = sorted({names_map.get(i, "") for i in extended_ids} - {""})
            with open(extended_supergroups_file, "w") as ef:
                for nm in ext_names:
                    ef.write(nm + "\n")
            working_ids = extended_ids

        # Final filter: BOTH endpoints in the (possibly extended) id-set.
        both = id1.isin(working_ids) & id2.isin(working_ids)
        sel = df[both]
        with open(output_file, "a") as out:
            for _, row in sel.iterrows():
                out.write(_format_store_row(
                    row, int(row["group_1_id"]), int(row["group_2_id"]),
                    names_map, has_pvalue) + "\n")
        store.close()

    elif extend:
        # Pairs with EITHER endpoint in the group set, then collect BOTH endpoints
        # as the extended set. The metric cutoff is applied only when one was given;
        # without -m/-c (metric==NA, cutoff==-1) extend over all such pairs (issue
        # 095 -- awk_column/awk_op are unset then; matches the parquet-store extend
        # path, which likewise skips the cutoff here).
        either_in = f"""grep '^[^#;]' {pairwise_file} | tail -n+2 | LC_ALL=C awk -F'\t' 'BEGIN {{ while ( getline < "{groups_file}" ) {{ gsub(/"/, "", $1); id_map[tolower($1)]=1 }} }} {{ if ( (tolower($3) in id_map) || (tolower($4) in id_map) ) {{ print $0 }} }}'"""
        if cutoff != -1:
            collect = f"""awk -F'\t' '{{if (${awk_column} {awk_op} {cutoff}) {{ print $3 >> "{extended_ids_list}"; print $4 >> "{extended_ids_list}"}}}}'"""
        else:
            collect = f"""awk -F'\t' '{{ print $3 >> "{extended_ids_list}"; print $4 >> "{extended_ids_list}" }}'"""
        awk_script = either_in + " | " + collect

        result = execute_bash_command(awk_script)
        bash_script = f"""sort -u {extended_ids_list} -o {extended_supergroups_file}"""
        result = execute_bash_command(bash_script)
        groups_file = extended_supergroups_file

    # filter by both cutoff and groups
    if store_for_filter is not None:
        pass  # handled above in the parquet-store branch
    elif cutoff != -1 and groups_file != "NA":
        awk_script = f"""grep '^[^#;]' {pairwise_file} | tail -n+2 | LC_ALL=C awk -F'\t' 'BEGIN {{ while ( getline < "{groups_file}" ) {{ gsub(/"/, "", $1); id_map[tolower($1)]=1 }} }} {{ if ( (tolower($3) in id_map) && (tolower($4) in id_map) ) {{ print $0 }} }}' | awk -F'\t' '{{if (${awk_column} {awk_op} {cutoff}) print $0}}' >> {output_file}"""
        result = execute_bash_command(awk_script)

    elif cutoff != -1:
        ctx.obj.INFO(
            f"Querying the pairwise matrix on the {metric} column with a cutoff of {cutoff}.")
        # Try Parquet/PairwiseStore first
        try:
            from dbretina.compat import open_pairwise
            store = open_pairwise(pairwise_file)
        except Exception:
            store = None

        if store is not None:
            ctx.obj.INFO("using Parquet pairwise data via PairwiseStore")
            if metric == "pvalue" and not store.has_pvalue:
                store.close()
                ctx.obj.ERROR("pvalue not found in pairwise file!")
            names_map = store.get_names_map()
            cols = ["group_1_id", "group_2_id", "shared_features",
                    "containment", "ochiai", "jaccard", "csi", "dice", "odds_ratio"]
            if store.has_pvalue:
                cols.append("pvalue")
            df = store.to_pandas(metric=metric, cutoff=cutoff, columns=cols)
            with open(output_file, 'a') as out:
                for _, row in df.iterrows():
                    gid1 = int(row['group_1_id'])
                    gid2 = int(row['group_2_id'])
                    fields = [str(gid1), str(gid2),
                             names_map.get(gid1, ""), names_map.get(gid2, ""),
                             str(int(row['shared_features'])),
                             f"{row['containment']:.1f}", f"{row['ochiai']:.1f}",
                             f"{row['jaccard']:.1f}", f"{row['csi']:.1f}",
                             f"{row['dice']:.1f}", f"{row['odds_ratio']:.1f}"]
                    if store.has_pvalue:
                        fields.append(str(row['pvalue']))
                    out.write('\t'.join(fields) + '\n')
            store.close()
        else:
            # Fallback: existing .dbrp / TSV code. Resolve the canonical .dbrp
            # via resolve_dbrp_path (an existing *file* or None) instead of the
            # raw str.replace, which no-op'd for a directory/.dbrp form and then
            # leaked the directory path into the binary reader -> "Invalid .dbrp
            # file (bad magic bytes)" (issue 052).
            from dbretina.compat import resolve_dbrp_path
            dbrp_path = resolve_dbrp_path(pairwise_file)
            metric_name_to_id = {
                "containment": 0, "ochiai": 1, "jaccard": 2, "csi": 3,
                "dice": 4, "odds_ratio": 5, "pvalue": 6,
            }
            if dbrp_path is None and os.path.isdir(pairwise_file):
                # A directory with no usable parquet store (open_pairwise
                # returned None: no manifest.json) and no sibling .dbrp; the awk
                # TSV fallback below would also fail on a directory.
                ctx.obj.ERROR(
                    f"'{pairwise_file}' is a pairwise directory without a usable "
                    f"parquet store (no manifest.json) and no sibling .dbrp; pass "
                    f"the pairwise TSV, its parquet directory, or the .dbrp to -p.")
            if dbrp_path is not None and metric in metric_name_to_id:
                mid = metric_name_to_id[metric]
                records = dbretina_internal.dbrp_filter_pairs(dbrp_path, mid, cutoff)
                with open(output_file, 'a') as out:
                    for rec in records:
                        fields = [str(rec['group_1_id']), str(rec['group_2_id']),
                                 rec['group_1_name'], rec['group_2_name'],
                                 str(rec['shared_features']),
                                 f"{rec['containment']:.1f}", f"{rec['ochiai']:.1f}",
                                 f"{rec['jaccard']:.1f}", f"{rec['csi']:.1f}",
                                 f"{rec['dice']:.1f}", f"{rec['odds_ratio']:.1f}"]
                        if 'pvalue' in rec:
                            fields.append(str(rec['pvalue']))
                        out.write('\t'.join(fields) + '\n')
            else:
                command = f"grep '^[^#;]' {pairwise_file} | tail -n+2 | LC_ALL=C awk -F'\t' '{{if (${awk_column} {awk_op} {cutoff}) print $0}}' >> {output_file}"
                result = execute_bash_command(command)

    elif groups_file != "NA":
        ctx.obj.INFO(f"Querying by groups file {groups_file}\nPlease wait...")
        awk_script = f"""grep '^[^#;]' {pairwise_file} | tail -n+2 |LC_ALL=C awk -F'\t' 'BEGIN {{ while ( getline < "{groups_file}" ) {{ gsub(/"/, "", $1); id_map[tolower($1)]=1 }} }} {{ if ( (tolower($3) in id_map) && (tolower($4) in id_map) ) {{ print $0 }} }}' >> {output_file}"""
        result = execute_bash_command(awk_script)


    # if _tmp_file exists, remove it
    if os.path.exists(_tmp_file):
        os.remove(_tmp_file)
        
    if os.path.exists(extended_ids_list):
        os.remove(extended_ids_list)

    ctx.obj.SUCCESS("Done.")
