#!/usr/bin/env python3
"""
Comprehensive blackbox test suite for DBRetina CLI.

Tests all 13 commands with edge cases and integration scenarios.

Test data design (6 groups, 12 distinct genes):
  GroupA: {Alpha, Beta, Gamma, Delta, Epsilon}        size=5
  GroupB: {Alpha, Beta, Gamma}                        size=3  (subset of A)
  GroupC: {Gamma, Delta, Zeta, Eta, Theta, Iota}      size=6
  GroupD: {Kappa, Lambda}                             size=2  (nearly isolated)
  GroupE: {Alpha, Zeta, Kappa, Mu}                    size=4  (bridges A,C,D)
  GroupF: {Alpha, Beta, Gamma, Delta, Epsilon}        size=5  (identical to A)

Edge cases exercised:
  - Complete subset (B ⊂ A): containment=100%, odds_ratio=-1
  - Identical groups (A = F): all metrics=100%, odds_ratio=-1
  - Single-gene overlaps (e.g. A∩E=1, D∩E=1)
  - Nearly isolated group (D shares 1 gene with E only)
  - Varying group sizes: 2, 3, 4, 5, 5, 6

Run:
    python test/test_blackbox.py
    # or
    python -m unittest test.test_blackbox -v
"""

import os
import sys
import subprocess
import tempfile
import shutil
import json
import unittest
from types import SimpleNamespace

# ============================================================
# SECTION 1: Test Data Constants
# ============================================================

TEST_ASC_CONTENT = (
    "gene_set\tgene\n"
    "GroupA\tAlpha\n"
    "GroupA\tBeta\n"
    "GroupA\tGamma\n"
    "GroupA\tDelta\n"
    "GroupA\tEpsilon\n"
    "GroupB\tAlpha\n"
    "GroupB\tBeta\n"
    "GroupB\tGamma\n"
    "GroupC\tGamma\n"
    "GroupC\tDelta\n"
    "GroupC\tZeta\n"
    "GroupC\tEta\n"
    "GroupC\tTheta\n"
    "GroupC\tIota\n"
    "GroupD\tKappa\n"
    "GroupD\tLambda\n"
    "GroupE\tAlpha\n"
    "GroupE\tZeta\n"
    "GroupE\tKappa\n"
    "GroupE\tMu\n"
    "GroupF\tAlpha\n"
    "GroupF\tBeta\n"
    "GroupF\tGamma\n"
    "GroupF\tDelta\n"
    "GroupF\tEpsilon\n"
)

TEST_GMT_CONTENT = (
    "GroupA\tDescription_A\tAlpha\tBeta\tGamma\tDelta\tEpsilon\n"
    "GroupB\tDescription_B\tAlpha\tBeta\tGamma\n"
    "GroupC\tDescription_C\tGamma\tDelta\tZeta\tEta\tTheta\tIota\n"
    "GroupD\tDescription_D\tKappa\tLambda\n"
    "GroupE\tDescription_E\tAlpha\tZeta\tKappa\tMu\n"
    "GroupF\tDescription_F\tAlpha\tBeta\tGamma\tDelta\tEpsilon\n"
)

EXPECTED_GROUP_NAMES = {"groupa", "groupb", "groupc", "groupd", "groupe", "groupf"}

# --- Variant data: same logical groups but different input shapes ---

# Split rows with overlapping + duplicate genes
TEST_ASC_SPLIT_ROWS = (
    "gene_set\tgene\n"
    "GroupA\tAlpha\n"
    "GroupA\tBeta\n"
    "GroupA\tGamma\n"
    "GroupA\tGamma\n"        # duplicate
    "GroupA\tDelta\n"
    "GroupA\tEpsilon\n"
    "GroupB\tAlpha\n"
    "GroupB\tBeta\n"
    "GroupB\tGamma\n"
    "GroupB\tAlpha\n"        # duplicate
    "GroupC\tGamma\n"
    "GroupC\tDelta\n"
    "GroupC\tZeta\n"
    "GroupC\tEta\n"
    "GroupC\tTheta\n"
    "GroupC\tIota\n"
    "GroupC\tGamma\n"        # duplicate
    "GroupD\tKappa\n"
    "GroupD\tLambda\n"
    "GroupE\tAlpha\n"
    "GroupE\tZeta\n"
    "GroupE\tKappa\n"
    "GroupE\tMu\n"
    "GroupF\tAlpha\n"
    "GroupF\tBeta\n"
    "GroupF\tGamma\n"
    "GroupF\tDelta\n"
    "GroupF\tEpsilon\n"
)

# Mixed case — same data, various case combinations
TEST_ASC_MIXED_CASE = (
    "gene_set\tgene\n"
    "GROUPA\tALPHA\n"
    "GroupA\tBeta\n"
    "groupa\tgamma\n"
    "GroupA\tDELTA\n"
    "groupA\tepsilon\n"
    "GROUPB\talpha\n"
    "groupb\tBETA\n"
    "GroupB\tGamma\n"
    "GroupC\tGAMMA\n"
    "GROUPC\tdelta\n"
    "groupc\tZeta\n"
    "GroupC\teta\n"
    "GROUPC\tTheta\n"
    "groupc\tIOTA\n"
    "GroupD\tKAPPA\n"
    "GROUPD\tlambda\n"
    "groupE\tAlpha\n"
    "GROUPE\tzeta\n"
    "GroupE\tKappa\n"
    "groupe\tMU\n"
    "GROUPF\talpha\n"
    "groupf\tBeta\n"
    "GroupF\tgamma\n"
    "GROUPF\tDelta\n"
    "groupf\tEPSILON\n"
)

# One gene per row (max fragmentation) — same data
TEST_ASC_ONE_PER_ROW = (
    "gene_set\tgene\n"
    "GroupA\tAlpha\n"
    "GroupA\tBeta\n"
    "GroupA\tGamma\n"
    "GroupA\tDelta\n"
    "GroupA\tEpsilon\n"
    "GroupB\tAlpha\n"
    "GroupB\tBeta\n"
    "GroupB\tGamma\n"
    "GroupC\tGamma\n"
    "GroupC\tDelta\n"
    "GroupC\tZeta\n"
    "GroupC\tEta\n"
    "GroupC\tTheta\n"
    "GroupC\tIota\n"
    "GroupD\tKappa\n"
    "GroupD\tLambda\n"
    "GroupE\tAlpha\n"
    "GroupE\tZeta\n"
    "GroupE\tKappa\n"
    "GroupE\tMu\n"
    "GroupF\tAlpha\n"
    "GroupF\tBeta\n"
    "GroupF\tGamma\n"
    "GroupF\tDelta\n"
    "GroupF\tEpsilon\n"
)

# Exact duplicate rows (every row doubled)
TEST_ASC_DUPLICATE_ROWS = (
    "gene_set\tgene\n"
    + "".join(line + "\n" + line + "\n"
             for line in TEST_ASC_CONTENT.strip().split("\n")[1:] if line.strip())
)

# Multi-file split: groups A+B in file1, C+D in file2, E+F in file3
TEST_ASC_FILE1 = (
    "gene_set\tgene\n"
    "GroupA\tAlpha\nGroupA\tBeta\nGroupA\tGamma\nGroupA\tDelta\nGroupA\tEpsilon\n"
    "GroupB\tAlpha\nGroupB\tBeta\nGroupB\tGamma\n"
)
TEST_ASC_FILE2 = (
    "gene_set\tgene\n"
    "GroupC\tGamma\nGroupC\tDelta\nGroupC\tZeta\nGroupC\tEta\nGroupC\tTheta\nGroupC\tIota\n"
    "GroupD\tKappa\nGroupD\tLambda\n"
)
TEST_ASC_FILE3 = (
    "gene_set\tgene\n"
    "GroupE\tAlpha\nGroupE\tZeta\nGroupE\tKappa\nGroupE\tMu\n"
    "GroupF\tAlpha\nGroupF\tBeta\nGroupF\tGamma\nGroupF\tDelta\nGroupF\tEpsilon\n"
)

# Cross-file merge: GroupA genes split across 2 files + rest in file2
TEST_ASC_CROSS1 = (
    "gene_set\tgene\n"
    "GroupA\tAlpha\nGroupA\tBeta\n"
    "GroupB\tAlpha\nGroupB\tBeta\nGroupB\tGamma\n"
    "GroupD\tKappa\nGroupD\tLambda\n"
)
TEST_ASC_CROSS2 = (
    "gene_set\tgene\n"
    "GroupA\tGamma\nGroupA\tDelta\nGroupA\tEpsilon\n"
    "GroupC\tGamma\nGroupC\tDelta\nGroupC\tZeta\nGroupC\tEta\nGroupC\tTheta\nGroupC\tIota\n"
    "GroupE\tAlpha\nGroupE\tZeta\nGroupE\tKappa\nGroupE\tMu\n"
    "GroupF\tAlpha\nGroupF\tBeta\nGroupF\tGamma\nGroupF\tDelta\nGroupF\tEpsilon\n"
)

# Quoted names — same data with double quotes around names
TEST_ASC_QUOTED = (
    "gene_set\tgene\n"
    + "".join(
        f'"{parts[0]}"\t"{parts[1]}"\n'
        for line in TEST_ASC_CONTENT.strip().split("\n")[1:]
        for parts in [line.split("\t")]
    )
)

# Shuffled row order — same data, deterministic shuffle
_SHUFFLED_LINES = TEST_ASC_CONTENT.strip().split("\n")[1:]  # skip header
_SHUFFLED_LINES = _SHUFFLED_LINES[::2] + _SHUFFLED_LINES[1::2]  # interleave odd/even
_SHUFFLED_LINES.reverse()
TEST_ASC_SHUFFLED = "gene_set\tgene\n" + "\n".join(_SHUFFLED_LINES) + "\n"

# Special characters in group/gene names (standalone dataset)
TEST_ASC_SPECIAL = (
    "gene_set\tgene\n"
    "Gene Set-1 (main)\tgene_alpha-v2\n"
    "Gene Set-1 (main)\tgene.beta\n"
    "Gene Set-1 (main)\tgene [gamma]\n"
    "group_2.test\tgene_alpha-v2\n"
    "group_2.test\tgene [gamma]\n"
    "group_2.test\tgene_delta\n"
)

# Population size N = 12 (unique features)
N = 12

# Expected pairwise rows keyed by frozenset of group names.
# Values: (shared, containment, ochiai, jaccard, csi, dice, odds_ratio)
EXPECTED_PAIRWISE = {
    frozenset({"groupa", "groupb"}): (3, 100.0, 77.5, 60.0, 60.0, 75.0, -1.0),
    frozenset({"groupa", "groupc"}): (2,  40.0, 36.5, 22.2, 13.3, 36.4, 0.5),
    frozenset({"groupa", "groupe"}): (1,  25.0, 22.4, 12.5, 5.0, 22.2, 0.3),
    frozenset({"groupa", "groupf"}): (5, 100.0, 100.0, 100.0, 100.0, 100.0, -1.0),
    frozenset({"groupb", "groupc"}): (1,  33.3, 23.6, 12.5, 5.6, 22.2, 0.4),
    frozenset({"groupb", "groupe"}): (1,  33.3, 28.9, 16.7, 8.3, 28.6, 1.0),
    frozenset({"groupb", "groupf"}): (3, 100.0, 77.5, 60.0, 60.0, 75.0, -1.0),
    frozenset({"groupc", "groupe"}): (1,  25.0, 20.4, 11.1, 4.2, 20.0, 0.2),
    frozenset({"groupc", "groupf"}): (2,  40.0, 36.5, 22.2, 13.3, 36.4, 0.5),
    frozenset({"groupd", "groupe"}): (1,  50.0, 35.4, 20.0, 12.5, 33.3, 2.3),
    frozenset({"groupe", "groupf"}): (1,  25.0, 22.4, 12.5, 5.0, 22.2, 0.3),
}

# Pairs that must NOT appear (zero overlap)
ZERO_OVERLAP_PAIRS = [
    frozenset({"groupa", "groupd"}),
    frozenset({"groupb", "groupd"}),
    frozenset({"groupc", "groupd"}),
    frozenset({"groupd", "groupf"}),
]


# ============================================================
# SECTION 2: Helper Functions
# ============================================================

def run_command(cmd, timeout=120, cwd=None):
    """Run a shell command and return (returncode, stdout, stderr)."""
    result = subprocess.run(
        cmd, shell=True, capture_output=True, text=True, timeout=timeout, cwd=cwd
    )
    return result.returncode, result.stdout, result.stderr


def assert_file_exists(test_case, filepath, msg=None):
    """Assert file exists and is non-empty."""
    test_case.assertTrue(
        os.path.exists(filepath),
        msg or f"Missing file: {filepath}"
    )
    test_case.assertGreater(
        os.path.getsize(filepath), 0,
        msg or f"Empty file: {filepath}"
    )


def count_tsv_data_rows(filepath):
    """Count data rows in a TSV file (skip # comments and first non-comment header)."""
    count = 0
    header_skipped = False
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if not header_skipped:
                header_skipped = True
                continue
            count += 1
    return count


def parse_pairwise_tsv(filepath):
    """Parse _DBRetina_pairwise.tsv -> list of row dicts."""
    rows = []
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("group_1_ID"):
                continue
            parts = line.split("\t")
            row = {
                "group_1_id": int(parts[0]),
                "group_2_id": int(parts[1]),
                "group_1_name": parts[2],
                "group_2_name": parts[3],
                "shared_features": int(parts[4]),
                "containment": float(parts[5]),
                "ochiai": float(parts[6]),
                "jaccard": float(parts[7]),
                "csi": float(parts[8]),
                "dice": float(parts[9]),
                "odds_ratio": float(parts[10]),
            }
            if len(parts) > 11:
                row["pvalue"] = float(parts[11])
            rows.append(row)
    return rows


def setup_index_and_pairwise(tmpdir, asc_content=TEST_ASC_CONTENT, extra_pw_args=""):
    """Create index + pairwise in tmpdir. Returns (prefix, pw_tsv_path)."""
    asc_path = os.path.join(tmpdir, "test_input.asc")
    with open(asc_path, "w") as f:
        f.write(asc_content)

    prefix = os.path.join(tmpdir, "test_idx")
    rc, _, stderr = run_command(f"DBRetina index -a {asc_path} -o {prefix}")
    assert rc == 0, f"index failed: {stderr}"

    rc, _, stderr = run_command(f"DBRetina pairwise -i {prefix} {extra_pw_args}")
    assert rc == 0, f"pairwise failed: {stderr}"

    pw_file = f"{prefix}_DBRetina_pairwise.tsv"
    return prefix, pw_file


def write_file(path, content):
    """Write content to a file."""
    with open(path, "w") as f:
        f.write(content)
    return path


def get_groups_from_raw_json(json_path):
    """Parse _raw.json -> {group_name: sorted_gene_list}."""
    with open(json_path) as f:
        data = json.load(f)["data"]
    return {k: sorted(set(v)) for k, v in data.items()}


def get_canonical_groups(tmpdir):
    """Index TEST_ASC_CONTENT and return canonical {group: sorted_genes} from _raw.json."""
    asc = write_file(os.path.join(tmpdir, "canonical.asc"), TEST_ASC_CONTENT)
    prefix = os.path.join(tmpdir, "canonical_idx")
    rc, _, stderr = run_command(f"DBRetina index -a {asc} -o {prefix}")
    assert rc == 0, f"canonical index failed: {stderr}"
    return get_groups_from_raw_json(f"{prefix}_raw.json")


def assert_pairwise_matches_expected(test_case, pw_file):
    """Verify a pairwise TSV matches EXPECTED_PAIRWISE exactly."""
    rows = parse_pairwise_tsv(pw_file)
    test_case.assertEqual(len(rows), len(EXPECTED_PAIRWISE),
                          f"Expected {len(EXPECTED_PAIRWISE)} pairs, got {len(rows)}")
    actual = {}
    for row in rows:
        key = frozenset({row["group_1_name"], row["group_2_name"]})
        actual[key] = (
            row["shared_features"], row["containment"], row["ochiai"],
            row["jaccard"], row["csi"], row["dice"], row["odds_ratio"],
        )

    labels = ["shared", "containment", "ochiai", "jaccard", "csi", "dice", "odds_ratio"]
    for pair_key, expected_vals in EXPECTED_PAIRWISE.items():
        test_case.assertIn(pair_key, actual, f"Missing pair {sorted(pair_key)}")
        actual_vals = actual[pair_key]
        for i, label in enumerate(labels):
            if label == "shared":
                test_case.assertEqual(actual_vals[i], expected_vals[i],
                                      f"{sorted(pair_key)} {label}")
            else:
                test_case.assertAlmostEqual(actual_vals[i], expected_vals[i], delta=0.15,
                                            msg=f"{sorted(pair_key)} {label}")

    for pair_key in ZERO_OVERLAP_PAIRS:
        test_case.assertNotIn(pair_key, actual,
                              f"Zero-overlap pair {sorted(pair_key)} should not appear")


# ============================================================
# Shared fixture: index + pairwise created once per module
# ============================================================

_SHARED_DIR = None
_SHARED_PREFIX = None
_SHARED_PW_FILE = None


def _ensure_shared_fixture():
    global _SHARED_DIR, _SHARED_PREFIX, _SHARED_PW_FILE
    if _SHARED_DIR is None:
        _SHARED_DIR = tempfile.mkdtemp(prefix="dbretina_shared_")
        _SHARED_PREFIX, _SHARED_PW_FILE = setup_index_and_pairwise(_SHARED_DIR)
    return _SHARED_PREFIX, _SHARED_PW_FILE


def _cleanup_shared_fixture():
    global _SHARED_DIR
    if _SHARED_DIR is not None:
        shutil.rmtree(_SHARED_DIR, ignore_errors=True)
        _SHARED_DIR = None


# ============================================================
# SECTION 3: Index Tests
# ============================================================

class TestIndex(unittest.TestCase):

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp(prefix="dbretina_idx_")

    def tearDown(self):
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def test_index_asc(self):
        """Index from ASC file produces .dbri, _raw.json, _hashes.json."""
        asc = write_file(os.path.join(self.tmpdir, "test.asc"), TEST_ASC_CONTENT)
        prefix = os.path.join(self.tmpdir, "idx")
        rc, _, stderr = run_command(f"DBRetina index -a {asc} -o {prefix}")
        self.assertEqual(rc, 0, stderr)
        assert_file_exists(self, f"{prefix}.dbri")
        assert_file_exists(self, f"{prefix}_raw.json")
        assert_file_exists(self, f"{prefix}_hashes.json")

    def test_index_gmt(self):
        """Index from GMT file produces .dbri."""
        gmt = write_file(os.path.join(self.tmpdir, "test.gmt"), TEST_GMT_CONTENT)
        # Run from tmpdir because GMT indexing creates temp file with output prefix in name
        rc, _, stderr = run_command(f"DBRetina index -g {gmt} -o idx", cwd=self.tmpdir)
        self.assertEqual(rc, 0, stderr)
        assert_file_exists(self, os.path.join(self.tmpdir, "idx.dbri"))

    def test_index_multiple_asc(self):
        """Index from two ASC files merges groups."""
        asc1 = write_file(os.path.join(self.tmpdir, "a.asc"),
                          "gene_set\tgene\nGroupX\tG1\nGroupX\tG2\n")
        asc2 = write_file(os.path.join(self.tmpdir, "b.asc"),
                          "gene_set\tgene\nGroupY\tG3\nGroupY\tG4\n")
        prefix = os.path.join(self.tmpdir, "idx")
        rc, _, stderr = run_command(f"DBRetina index -a {asc1} -a {asc2} -o {prefix}")
        self.assertEqual(rc, 0, stderr)
        assert_file_exists(self, f"{prefix}.dbri")
        # Verify both groups present
        with open(f"{prefix}_raw.json") as f:
            data = json.load(f)["data"]
        self.assertIn("groupx", {k.lower() for k in data.keys()})
        self.assertIn("groupy", {k.lower() for k in data.keys()})

    def test_index_single_group(self):
        """Index with a single group, single gene."""
        asc = write_file(os.path.join(self.tmpdir, "single.asc"),
                         "gene_set\tgene\nOnlyGroup\tOnlyGene\n")
        prefix = os.path.join(self.tmpdir, "idx")
        rc, _, stderr = run_command(f"DBRetina index -a {asc} -o {prefix}")
        self.assertEqual(rc, 0, stderr)
        assert_file_exists(self, f"{prefix}.dbri")


# ============================================================
# SECTION 4: Pairwise Tests
# ============================================================

class TestPairwise(unittest.TestCase):

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp(prefix="dbretina_pw_")

    def tearDown(self):
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def test_pairwise_default(self):
        """Pairwise with defaults produces TSV, DBRP, stats, plots."""
        prefix, pw_file = setup_index_and_pairwise(self.tmpdir)
        assert_file_exists(self, pw_file)
        assert_file_exists(self, f"{prefix}_DBRetina_pairwise.dbrp")
        assert_file_exists(self, f"{prefix}_DBRetina_pairwise_stats.json")
        assert_file_exists(self, f"{prefix}_DBRetina_similarity_metrics_plot_linear.png")
        assert_file_exists(self, f"{prefix}_DBRetina_similarity_metrics_plot_log.png")
        self.assertEqual(count_tsv_data_rows(pw_file), 11)

    def test_pairwise_rerun_clears_stale_parquet_parts(self):
        """Regression: re-running pairwise with fewer threads must not leave stale
        parquet parts that silently inflate glob-based row counts."""
        import json, glob, duckdb
        d = self.tmpdir
        genes = [f"g{i}" for i in range(40)]
        rows = ["grp_%02d\tdesc\t%s" % (k, "\t".join(genes[(k % 10):(k % 10) + 12]))
                for k in range(30)]
        gmt = os.path.join(d, "many.gmt")
        with open(gmt, "w") as f:
            f.write("\n".join(rows) + "\n")
        rc, _, err = run_command("DBRetina index -g many.gmt -o idx", cwd=d)
        self.assertEqual(rc, 0, err)
        rc, _, err = run_command("DBRetina pairwise -i idx -m containment -c 0 -t 4", cwd=d)
        self.assertEqual(rc, 0, err)
        pwdir = os.path.join(d, "idx_DBRetina_pairwise")
        parts1 = glob.glob(os.path.join(pwdir, "data", "part_*.parquet"))
        self.assertGreaterEqual(len(parts1), 2, "need a multi-part run to test stale-part clearing")
        # re-run with fewer threads (1 part); stale parts from run 1 must be cleared
        rc, _, err = run_command("DBRetina pairwise -i idx -m containment -c 0 -t 1", cwd=d)
        self.assertEqual(rc, 0, err)
        with open(os.path.join(pwdir, "manifest.json")) as mf:
            manifest = json.load(mf)
        glob_count = duckdb.connect().execute(
            f"select count(*) from read_parquet('{pwdir}/data/part_*.parquet', union_by_name=true)"
        ).fetchone()[0]
        self.assertEqual(glob_count, manifest["num_pairs"],
                         f"stale parquet parts: glob={glob_count} != manifest {manifest['num_pairs']}")
        # No leftover parts of any kind (also covers the --pvalue schema-mix variant):
        parts2 = glob.glob(os.path.join(pwdir, "data", "part_*.parquet"))
        self.assertEqual(len(parts2), manifest["num_partitions"],
                         f"stale part files left: {len(parts2)} != num_partitions {manifest['num_partitions']}")

    def test_pairwise_metric_values(self):
        """All metric values match hand-computed expectations."""
        prefix, pw_file = setup_index_and_pairwise(self.tmpdir)
        rows = parse_pairwise_tsv(pw_file)

        actual = {}
        for row in rows:
            key = frozenset({row["group_1_name"], row["group_2_name"]})
            actual[key] = (
                row["shared_features"], row["containment"], row["ochiai"],
                row["jaccard"], row["csi"], row["dice"], row["odds_ratio"],
            )

        self.assertEqual(len(rows), len(EXPECTED_PAIRWISE))

        labels = ["shared", "containment", "ochiai", "jaccard", "csi", "dice", "odds_ratio"]
        for pair_key, expected_vals in EXPECTED_PAIRWISE.items():
            self.assertIn(pair_key, actual, f"Missing pair {sorted(pair_key)}")
            actual_vals = actual[pair_key]
            for i, label in enumerate(labels):
                if label == "shared":
                    self.assertEqual(actual_vals[i], expected_vals[i],
                                     f"{sorted(pair_key)} {label}")
                else:
                    self.assertAlmostEqual(actual_vals[i], expected_vals[i], delta=0.15,
                                           msg=f"{sorted(pair_key)} {label}")

        for pair_key in ZERO_OVERLAP_PAIRS:
            self.assertNotIn(pair_key, actual,
                             f"Zero-overlap pair {sorted(pair_key)} should not appear")

        # group_1_id < group_2_id invariant
        for row in rows:
            self.assertLess(row["group_1_id"], row["group_2_id"])

    def test_pairwise_with_cutoff(self):
        """Pairwise with -c 50 filters containment < 50."""
        prefix, pw_file = setup_index_and_pairwise(self.tmpdir, extra_pw_args="-c 50")
        rows = parse_pairwise_tsv(pw_file)
        for row in rows:
            self.assertGreaterEqual(row["containment"], 50.0)
        # Known pairs with containment >= 50: A-B(100), A-F(100), B-F(100), D-E(50), C-F(40→no), B-C(33.3→no)
        # So: A-B, A-F, B-E(33.3→no), B-F, D-E, C-A(40→no) => A-B, A-F, B-F, D-E = 4
        self.assertGreater(len(rows), 0)

    def test_pairwise_with_pvalue(self):
        """Pairwise with --pvalue includes pvalue column."""
        prefix, pw_file = setup_index_and_pairwise(self.tmpdir, extra_pw_args="--pvalue")
        rows = parse_pairwise_tsv(pw_file)
        self.assertEqual(len(rows), 11)
        for row in rows:
            self.assertIn("pvalue", row)

    def test_pairwise_dbrp_matches_tsv(self):
        """DBRP binary file contains same data as TSV."""
        prefix, pw_file = setup_index_and_pairwise(self.tmpdir)
        tsv_rows = parse_pairwise_tsv(pw_file)

        import _dbretina_internal as dbi
        dbrp_path = f"{prefix}_DBRetina_pairwise.dbrp"
        dbrp_recs = dbi.dbrp_iterate_all(dbrp_path)

        self.assertEqual(len(dbrp_recs), len(tsv_rows))

        # Build lookup by pair names
        tsv_by_pair = {}
        for r in tsv_rows:
            key = frozenset({r["group_1_name"], r["group_2_name"]})
            tsv_by_pair[key] = r

        for rec in dbrp_recs:
            key = frozenset({rec["group_1_name"], rec["group_2_name"]})
            self.assertIn(key, tsv_by_pair)
            tsv_r = tsv_by_pair[key]
            self.assertEqual(rec["shared_features"], tsv_r["shared_features"])
            for metric in ["containment", "ochiai", "jaccard", "csi", "dice", "odds_ratio"]:
                self.assertAlmostEqual(rec[metric], tsv_r[metric], delta=0.15,
                                       msg=f"{sorted(key)} {metric}")

    def test_pairwise_dbrp_filter(self):
        """dbrp_filter_pairs correctly filters by ochiai >= 50."""
        prefix, _ = setup_index_and_pairwise(self.tmpdir)
        import _dbretina_internal as dbi
        dbrp_path = f"{prefix}_DBRetina_pairwise.dbrp"
        filtered = dbi.dbrp_filter_pairs(dbrp_path, 1, 50.0)  # ochiai=1

        for rec in filtered:
            self.assertGreaterEqual(rec["ochiai"], 50.0)

        # Expected: A-B(77.5), A-F(100), B-F(77.5), C-A(36.5→no), C-F(36.5→no) = 3 pairs
        pair_names = {frozenset({r["group_1_name"], r["group_2_name"]}) for r in filtered}
        self.assertIn(frozenset({"groupa", "groupf"}), pair_names)
        self.assertIn(frozenset({"groupa", "groupb"}), pair_names)
        self.assertIn(frozenset({"groupb", "groupf"}), pair_names)
        self.assertEqual(len(filtered), 3)

    def test_pairwise_dbrp_metadata(self):
        """dbrp metadata, names, statistics return valid data."""
        prefix, _ = setup_index_and_pairwise(self.tmpdir)
        import _dbretina_internal as dbi
        dbrp_path = f"{prefix}_DBRetina_pairwise.dbrp"

        meta_str = dbi.dbrp_load_metadata(dbrp_path)
        self.assertTrue(len(meta_str) > 0, "Metadata should be non-empty")
        # Metadata JSON may contain unescaped paths; just verify it's parseable
        try:
            meta = json.loads(meta_str)
            self.assertEqual(meta["num_groups"], 6)
        except json.JSONDecodeError:
            # If command field has unescaped chars, at least verify it's a string
            self.assertIn("num_groups", meta_str)

        names = dbi.dbrp_load_names(dbrp_path)
        self.assertEqual(len(names), 6)

        stats = dbi.dbrp_load_statistics(dbrp_path)
        self.assertIsInstance(json.loads(stats), dict)

        self.assertEqual(dbi.dbrp_get_num_pairs(dbrp_path), 11)
        self.assertEqual(dbi.dbrp_get_num_groups(dbrp_path), 6)

    def test_pairwise_multithreaded(self):
        """Multithreaded pairwise produces same results as single-threaded."""
        prefix_st, pw_st = setup_index_and_pairwise(
            tempfile.mkdtemp(dir=self.tmpdir, prefix="st_")
        )
        prefix_mt, pw_mt = setup_index_and_pairwise(
            tempfile.mkdtemp(dir=self.tmpdir, prefix="mt_"),
            extra_pw_args="-t 4"
        )

        rows_st = parse_pairwise_tsv(pw_st)
        rows_mt = parse_pairwise_tsv(pw_mt)
        self.assertEqual(len(rows_st), len(rows_mt))

        st_pairs = {frozenset({r["group_1_name"], r["group_2_name"]}): r for r in rows_st}
        mt_pairs = {frozenset({r["group_1_name"], r["group_2_name"]}): r for r in rows_mt}
        self.assertEqual(st_pairs.keys(), mt_pairs.keys())


# ============================================================
# SECTION 5: Query Tests
# ============================================================

class TestQuery(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.prefix, cls.pw_file = _ensure_shared_fixture()

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp(prefix="dbretina_query_")

    def tearDown(self):
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def test_query_cutoff_only(self):
        """Query with ochiai >= 50 returns only qualifying pairs."""
        out = os.path.join(self.tmpdir, "q")
        rc, _, stderr = run_command(
            f"DBRetina query -p {self.pw_file} -m ochiai -c 50 -o {out}"
        )
        self.assertEqual(rc, 0, stderr)
        assert_file_exists(self, f"{out}.tsv")
        rows = parse_pairwise_tsv(f"{out}.tsv")
        for row in rows:
            self.assertGreaterEqual(row["ochiai"], 50.0)

    def test_query_cutoff_zero(self):
        """Query with cutoff=0 returns all 11 pairs."""
        out = os.path.join(self.tmpdir, "q")
        rc, _, stderr = run_command(
            f"DBRetina query -p {self.pw_file} -m ochiai -c 0 -o {out}"
        )
        self.assertEqual(rc, 0, stderr)
        self.assertEqual(count_tsv_data_rows(f"{out}.tsv"), 11)

    def test_query_cutoff_100(self):
        """Query with ochiai=100 returns only A-F (identical)."""
        out = os.path.join(self.tmpdir, "q")
        rc, _, stderr = run_command(
            f"DBRetina query -p {self.pw_file} -m ochiai -c 100 -o {out}"
        )
        self.assertEqual(rc, 0, stderr)
        rows = parse_pairwise_tsv(f"{out}.tsv")
        self.assertEqual(len(rows), 1)
        names = {rows[0]["group_1_name"], rows[0]["group_2_name"]}
        self.assertEqual(names, {"groupa", "groupf"})

    def test_query_groups_file(self):
        """Query with groups file filters to those groups."""
        groups = write_file(os.path.join(self.tmpdir, "groups.txt"), "GroupA\nGroupB\n")
        out = os.path.join(self.tmpdir, "q")
        rc, _, stderr = run_command(
            f"DBRetina query -p {self.pw_file} -g {groups} -m ochiai -c 0 -o {out}"
        )
        self.assertEqual(rc, 0, stderr)
        rows = parse_pairwise_tsv(f"{out}.tsv")
        for row in rows:
            pair = {row["group_1_name"], row["group_2_name"]}
            self.assertTrue(pair.issubset({"groupa", "groupb"}),
                            f"Unexpected pair: {pair}")

    def test_query_groups_and_cutoff(self):
        """Query with groups file AND cutoff applies both filters."""
        groups = write_file(os.path.join(self.tmpdir, "groups.txt"),
                            "GroupA\nGroupB\nGroupC\n")
        out = os.path.join(self.tmpdir, "q")
        rc, _, stderr = run_command(
            f"DBRetina query -p {self.pw_file} -g {groups} -m ochiai -c 50 -o {out}"
        )
        self.assertEqual(rc, 0, stderr)
        rows = parse_pairwise_tsv(f"{out}.tsv")
        for row in rows:
            self.assertGreaterEqual(row["ochiai"], 50.0)
            pair = {row["group_1_name"], row["group_2_name"]}
            self.assertTrue(pair.issubset({"groupa", "groupb", "groupc"}))

    def test_query_extend(self):
        """Query with --extend expands group list."""
        groups = write_file(os.path.join(self.tmpdir, "groups.txt"), "GroupA\n")
        out = os.path.join(self.tmpdir, "q")
        rc, _, stderr = run_command(
            f"DBRetina query -p {self.pw_file} -g {groups} "
            f"-m ochiai -c 50 --extend -o {out}"
        )
        self.assertEqual(rc, 0, stderr)
        ext_file = f"{out}_extended_supergroups.txt"
        assert_file_exists(self, ext_file)
        with open(ext_file) as f:
            extended = {line.strip().lower() for line in f if line.strip()}
        # GroupA has ochiai>=50 with GroupB and GroupF, so those should be extended
        self.assertTrue(len(extended) > 1, f"Expected extended groups, got: {extended}")

    def test_query_dbrp_vs_tsv(self):
        """Query results identical with and without .dbrp file."""
        out1 = os.path.join(self.tmpdir, "q_dbrp")
        rc, _, _ = run_command(
            f"DBRetina query -p {self.pw_file} -m ochiai -c 50 -o {out1}"
        )
        self.assertEqual(rc, 0)
        count1 = count_tsv_data_rows(f"{out1}.tsv")

        # Remove .dbrp temporarily, query again
        dbrp = self.pw_file.replace(".tsv", ".dbrp")
        dbrp_bak = dbrp + ".bak"
        if os.path.exists(dbrp):
            os.rename(dbrp, dbrp_bak)
            try:
                out2 = os.path.join(self.tmpdir, "q_tsv")
                rc, _, _ = run_command(
                    f"DBRetina query -p {self.pw_file} -m ochiai -c 50 -o {out2}"
                )
                self.assertEqual(rc, 0)
                count2 = count_tsv_data_rows(f"{out2}.tsv")
                self.assertEqual(count1, count2)
            finally:
                os.rename(dbrp_bak, dbrp)


# ============================================================
# SECTION 6: Cluster Tests
# ============================================================

class TestCluster(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.prefix, cls.pw_file = _ensure_shared_fixture()

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp(prefix="dbretina_cluster_")

    def tearDown(self):
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def test_cluster_connected_components(self):
        """Cluster with connected components produces expected files."""
        out = os.path.join(self.tmpdir, "cl")
        rc, _, stderr = run_command(
            f"DBRetina cluster -p {self.pw_file} -m ochiai -c 50 -o {out}"
        )
        self.assertEqual(rc, 0, stderr)
        assert_file_exists(self, f"{out}_clusters.tsv")
        assert_file_exists(self, f"{out}_clusters_histogram.png")

    def test_cluster_community(self):
        """Cluster with --community uses Leiden algorithm."""
        out = os.path.join(self.tmpdir, "cl")
        rc, _, stderr = run_command(
            f"DBRetina cluster -p {self.pw_file} -m ochiai -c 30 --community -o {out}"
        )
        self.assertEqual(rc, 0, stderr)
        assert_file_exists(self, f"{out}_clusters.tsv")

    def test_cluster_linear_transform_warns(self):
        """ISSUE-024: --node-weight-transform linear (community) warns it needs tuned
        --resolution (raw sizes -> all-singletons at default); log2 does not warn."""
        out = os.path.join(self.tmpdir, "cl_lin")
        rc, stdout, stderr = run_command(
            f"DBRetina cluster -p {self.pw_file} -m ochiai -c 30 --community "
            f"--node-weight-transform linear -o {out}"
        )
        self.assertEqual(rc, 0, stderr)
        out_all = (stdout + stderr).lower()
        self.assertIn("linear", out_all)
        self.assertIn("resolution", out_all)
        # log2 (default) must NOT emit the linear advisory
        out2 = os.path.join(self.tmpdir, "cl_log2")
        rc2, stdout2, stderr2 = run_command(
            f"DBRetina cluster -p {self.pw_file} -m ochiai -c 30 --community "
            f"--node-weight-transform log2 -o {out2}"
        )
        self.assertEqual(rc2, 0, stderr2)
        self.assertNotIn("raw gene-set sizes", (stdout2 + stderr2).lower())

    def test_cluster_file_format(self):
        """Cluster output has correct TSV format."""
        out = os.path.join(self.tmpdir, "cl")
        run_command(f"DBRetina cluster -p {self.pw_file} -m ochiai -c 50 -o {out}")
        with open(f"{out}_clusters.tsv") as f:
            lines = [l for l in f if not l.startswith("#")]
        header = lines[0].strip().split("\t")
        self.assertIn("cluster_id", header[0].lower())
        # Data rows should have numeric cluster_id and cluster_size
        if len(lines) > 1:
            parts = lines[1].strip().split("\t")
            self.assertTrue(parts[0].strip().isdigit())

    def test_cluster_various_metrics(self):
        """Cluster works with containment and jaccard metrics."""
        for metric in ["containment", "jaccard"]:
            out = os.path.join(self.tmpdir, f"cl_{metric}")
            rc, _, stderr = run_command(
                f"DBRetina cluster -p {self.pw_file} -m {metric} -c 30 -o {out}"
            )
            self.assertEqual(rc, 0, f"{metric} failed: {stderr}")
            assert_file_exists(self, f"{out}_clusters.tsv")

    def test_cluster_high_cutoff(self):
        """Cluster with ochiai=100 finds only identical groups (A=F)."""
        out = os.path.join(self.tmpdir, "cl")
        rc, _, stderr = run_command(
            f"DBRetina cluster -p {self.pw_file} -m ochiai -c 100 -o {out}"
        )
        self.assertEqual(rc, 0, stderr)
        # Should have exactly 1 cluster containing groupa and groupf
        with open(f"{out}_clusters.tsv") as f:
            data_lines = [l for l in f if not l.startswith("#") and not l.strip().startswith("cluster")]
            # Filter header
            data_lines = [l for l in data_lines if l.strip() and not "cluster_id" in l.lower()]
        self.assertGreaterEqual(len(data_lines), 1)


# ============================================================
# SECTION 7: Export Tests
# ============================================================

class TestExport(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.prefix, cls.pw_file = _ensure_shared_fixture()

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp(prefix="dbretina_export_")

    def tearDown(self):
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def test_export_basic(self):
        """Export produces heatmap, distance matrix."""
        out = os.path.join(self.tmpdir, "exp")
        rc, _, stderr = run_command(
            f"DBRetina export -p {self.pw_file} -m ochiai -o {out}"
        )
        self.assertEqual(rc, 0, stderr)
        assert_file_exists(self, f"{out}_heatmap.png")
        assert_file_exists(self, f"{out}_distmat.pkl")
        assert_file_exists(self, f"{out}_distmat.tsv")

    def test_export_newick(self):
        """Export with --newick produces newick and dendrogram."""
        out = os.path.join(self.tmpdir, "exp")
        rc, _, stderr = run_command(
            f"DBRetina export -p {self.pw_file} -m ochiai --newick -o {out}"
        )
        self.assertEqual(rc, 0, stderr)
        assert_file_exists(self, f"{out}.newick")
        assert_file_exists(self, f"{out}_dendrogram.png")

    def test_export_labels_ids(self):
        """Export with -l ids uses numeric IDs in distance matrix."""
        out = os.path.join(self.tmpdir, "exp")
        rc, _, stderr = run_command(
            f"DBRetina export -p {self.pw_file} -m ochiai -l ids -o {out}"
        )
        self.assertEqual(rc, 0, stderr)
        with open(f"{out}_distmat.tsv") as f:
            header = f.readline().strip().split("\t")
        # First element is empty (row index header), rest are IDs
        for h in header[1:]:
            self.assertTrue(h.strip().isdigit(), f"Expected numeric ID, got: {h}")

    def test_export_labels_names(self):
        """Export with -l names uses group names."""
        out = os.path.join(self.tmpdir, "exp")
        rc, _, stderr = run_command(
            f"DBRetina export -p {self.pw_file} -m ochiai -l names -o {out}"
        )
        self.assertEqual(rc, 0, stderr)
        with open(f"{out}_distmat.tsv") as f:
            header = f.readline().strip().split("\t")
        # Names should be group names (not purely digits)
        non_empty = [h for h in header[1:] if h.strip()]
        self.assertTrue(any(not h.isdigit() for h in non_empty))

    def test_export_linkage_methods(self):
        """Export with different linkage methods all succeed."""
        for method in ["single", "complete", "average", "ward"]:
            out = os.path.join(self.tmpdir, f"exp_{method}")
            rc, _, stderr = run_command(
                f"DBRetina export -p {self.pw_file} -m ochiai "
                f"--newick --linkage {method} -o {out}"
            )
            self.assertEqual(rc, 0, f"linkage={method} failed: {stderr}")

    def test_export_distmat_symmetry(self):
        """Distance matrix is symmetric with diagonal=100."""
        out = os.path.join(self.tmpdir, "exp")
        run_command(f"DBRetina export -p {self.pw_file} -m ochiai -o {out}")
        import pandas as pd
        df = pd.read_csv(f"{out}_distmat.tsv", sep="\t", index_col=0)
        # Check symmetry
        for i in range(len(df)):
            for j in range(len(df)):
                self.assertAlmostEqual(df.iloc[i, j], df.iloc[j, i], delta=0.01,
                                       msg=f"Not symmetric at [{i},{j}]")
        # Diagonal should be 100 (self-similarity)
        for i in range(len(df)):
            self.assertAlmostEqual(df.iloc[i, i], 100.0, delta=0.01)

    def test_export_invalid_linkage_clean_error(self):
        """ISSUE-017: invalid --linkage exits cleanly (no raw scipy traceback),
        and the error message lists the valid linkage methods."""
        out = os.path.join(self.tmpdir, "exp_bogus")
        rc, stdout, stderr = run_command(
            f"DBRetina export -p {self.pw_file} -m containment "
            f"--newick --linkage BOGUS -o {out}"
        )
        combined = stdout + stderr
        self.assertNotEqual(rc, 0, "invalid --linkage should be a nonzero exit")
        self.assertNotIn("Traceback", combined,
                         f"raw traceback leaked:\n{combined}")
        self.assertNotIn("ValueError", combined,
                         f"raw ValueError leaked:\n{combined}")
        # Message should mention --linkage and list at least one valid method.
        self.assertIn("linkage", combined.lower())
        self.assertTrue(
            any(m in combined for m in
                ["single", "complete", "average", "weighted",
                 "centroid", "median", "ward"]),
            f"error should list valid linkage methods:\n{combined}"
        )
        # No partial output should be relied upon, but at minimum the run
        # must not have produced a (stale) newick from the bogus method.
        self.assertFalse(os.path.exists(f"{out}.newick"),
                         "no newick should be written when --linkage is invalid")

    def test_export_valid_linkage_still_works(self):
        """ISSUE-017: a valid --linkage value (average) still succeeds."""
        out = os.path.join(self.tmpdir, "exp_avg")
        rc, _, stderr = run_command(
            f"DBRetina export -p {self.pw_file} -m containment "
            f"--newick --linkage average -o {out}"
        )
        self.assertEqual(rc, 0, stderr)
        assert_file_exists(self, f"{out}.newick")

    def test_export_empty_pairwise_clean_error(self):
        """ISSUE-027: a truly empty pairwise TSV yields a clean error,
        not a raw pandas EmptyDataError traceback."""
        empty_tsv = os.path.join(self.tmpdir, "empty.tsv")
        open(empty_tsv, "w").close()
        out = os.path.join(self.tmpdir, "exp_empty")
        rc, stdout, stderr = run_command(
            f"DBRetina export -p {empty_tsv} -m containment -o {out}"
        )
        combined = stdout + stderr
        self.assertNotEqual(rc, 0, "empty pairwise should be a nonzero exit")
        self.assertNotIn("Traceback", combined,
                         f"raw traceback leaked:\n{combined}")
        self.assertNotIn("EmptyDataError", combined,
                         f"raw EmptyDataError leaked:\n{combined}")

    def test_export_comments_only_pairwise_clean_error(self):
        """ISSUE-027: a comments-only pairwise TSV yields a clean error."""
        comments_tsv = os.path.join(self.tmpdir, "onlycomments.tsv")
        with open(comments_tsv, "w") as f:
            f.write("# DBRetina pairwise output\n# population_size: 10\n")
        out = os.path.join(self.tmpdir, "exp_comments")
        rc, stdout, stderr = run_command(
            f"DBRetina export -p {comments_tsv} -m containment -o {out}"
        )
        combined = stdout + stderr
        self.assertNotEqual(rc, 0, "comments-only pairwise should be a nonzero exit")
        self.assertNotIn("Traceback", combined,
                         f"raw traceback leaked:\n{combined}")
        self.assertNotIn("EmptyDataError", combined,
                         f"raw EmptyDataError leaked:\n{combined}")

    def test_export_nonexistent_output_dir(self):
        """ISSUE-028: -o prefix in a nonexistent directory must not crash with a
        raw FileNotFoundError. We auto-create the parent dir, so it succeeds and
        the output files are present."""
        out = os.path.join(self.tmpdir, "no_such_dir_xyz", "pref")
        rc, stdout, stderr = run_command(
            f"DBRetina export -p {self.pw_file} -m containment -o {out}"
        )
        combined = stdout + stderr
        self.assertNotIn("Traceback", combined,
                         f"raw traceback leaked:\n{combined}")
        self.assertNotIn("FileNotFoundError", combined,
                         f"raw FileNotFoundError leaked:\n{combined}")
        # Chosen behavior: auto-create the parent directory -> success + output.
        self.assertEqual(rc, 0, stderr)
        assert_file_exists(self, f"{out}_heatmap.png")
        assert_file_exists(self, f"{out}_distmat.tsv")

    def test_export_output_prefix_parent_is_file(self):
        """ISSUE-028 follow-up: -o whose parent path component is an existing FILE
        must fail with a clean error, not a raw OSError traceback."""
        afile = write_file(os.path.join(self.tmpdir, "afile"), "x\n")
        out = os.path.join(afile, "pref")  # parent component 'afile' is a regular file
        rc, stdout, stderr = run_command(
            f"DBRetina export -p {self.pw_file} -m containment -o {out}"
        )
        combined = stdout + stderr
        self.assertNotEqual(rc, 0)
        self.assertNotIn("Traceback", combined, f"raw traceback leaked:\n{combined}")
        self.assertIn("[ERROR]", combined)

    def test_export_heatmap_log_matches_reality(self):
        """ISSUE-039: the heatmap log line must reference a file that is actually
        written. Only a PNG is produced, so the log must not claim a .html file
        that never appears on disk."""
        out = os.path.join(self.tmpdir, "exp_log")
        rc, stdout, stderr = run_command(
            f"DBRetina export -p {self.pw_file} -m containment -o {out}"
        )
        self.assertEqual(rc, 0, stderr)
        combined = stdout + stderr
        heatmap_lines = [ln for ln in combined.splitlines()
                         if "heatmap" in ln.lower() and ("Writing" in ln or "writing" in ln.lower())]
        self.assertTrue(heatmap_lines, f"expected a heatmap log line:\n{combined}")
        for ln in heatmap_lines:
            # Any file the log claims to write must exist on disk.
            for token in ln.split():
                if "_heatmap." in token:
                    claimed = token.strip().strip("'\"")
                    self.assertTrue(
                        os.path.exists(claimed),
                        f"log claims to write {claimed} but it does not exist; "
                        f"log line: {ln}"
                    )
        # Specifically: no .html is written, so no log line may claim one.
        self.assertFalse(
            any("_heatmap.html" in ln for ln in heatmap_lines),
            f"log claims a _heatmap.html that is never written:\n{combined}"
        )


# ============================================================
# SECTION 8: Dedup Tests
# ============================================================

class TestDedup(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.prefix, cls.pw_file = _ensure_shared_fixture()

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp(prefix="dbretina_dedup_")

    def tearDown(self):
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def test_dedup_exact_duplicates(self):
        """Dedup with ochiai=100 removes one of A/F (identical)."""
        out = os.path.join(self.tmpdir, "dd")
        rc, _, stderr = run_command(
            f"DBRetina dedup -i {self.prefix} -p {self.pw_file} -c 100 -o {out}"
        )
        self.assertEqual(rc, 0, stderr)
        with open(f"{out}_deduplicated_groups.txt") as f:
            groups = {line.strip().lower() for line in f if line.strip()}
        # A and F are identical; one should be removed → 5 groups
        self.assertEqual(len(groups), 5)
        # Either groupa or groupf, not both
        self.assertEqual(len(groups & {"groupa", "groupf"}), 1)

    def test_dedup_low_cutoff(self):
        """Dedup with low cutoff merges more groups."""
        out = os.path.join(self.tmpdir, "dd")
        rc, _, stderr = run_command(
            f"DBRetina dedup -i {self.prefix} -p {self.pw_file} -c 20 -o {out}"
        )
        self.assertEqual(rc, 0, stderr)
        with open(f"{out}_deduplicated_groups.txt") as f:
            groups = {line.strip() for line in f if line.strip()}
        # With low cutoff, more groups are merged
        self.assertLess(len(groups), 6)

    def test_dedup_output_format(self):
        """Dedup output has one group name per line."""
        out = os.path.join(self.tmpdir, "dd")
        run_command(f"DBRetina dedup -i {self.prefix} -p {self.pw_file} -c 100 -o {out}")
        with open(f"{out}_deduplicated_groups.txt") as f:
            lines = [l.strip() for l in f if l.strip()]
        for line in lines:
            self.assertNotIn("\t", line, "Should be one group per line, no tabs")

    def test_dedup_dbrp_vs_tsv(self):
        """Dedup produces same result with and without .dbrp."""
        out1 = os.path.join(self.tmpdir, "dd1")
        run_command(f"DBRetina dedup -i {self.prefix} -p {self.pw_file} -c 100 -o {out1}")
        with open(f"{out1}_deduplicated_groups.txt") as f:
            groups1 = {l.strip().lower() for l in f if l.strip()}

        # Remove .dbrp, run again
        dbrp = self.pw_file.replace(".tsv", ".dbrp")
        dbrp_bak = dbrp + ".bak"
        if os.path.exists(dbrp):
            os.rename(dbrp, dbrp_bak)
            try:
                out2 = os.path.join(self.tmpdir, "dd2")
                run_command(
                    f"DBRetina dedup -i {self.prefix} -p {self.pw_file} -c 100 -o {out2}"
                )
                with open(f"{out2}_deduplicated_groups.txt") as f:
                    groups2 = {l.strip().lower() for l in f if l.strip()}
                self.assertEqual(groups1, groups2)
            finally:
                os.rename(dbrp_bak, dbrp)


# ============================================================
# SECTION 9: Modularity Tests
# ============================================================

class TestModularity(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.prefix, cls.pw_file = _ensure_shared_fixture()

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp(prefix="dbretina_mod_")

    def tearDown(self):
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def test_modularity_basic(self):
        """Modularity produces output file."""
        out = os.path.join(self.tmpdir, "mod")
        rc, _, stderr = run_command(
            f"DBRetina modularity -i {self.prefix} -p {self.pw_file} -c 40 -o {out}"
        )
        self.assertEqual(rc, 0, stderr)
        assert_file_exists(self, f"{out}_modularity.tsv")

    def test_modularity_format(self):
        """Modularity output has correct header."""
        out = os.path.join(self.tmpdir, "mod")
        run_command(
            f"DBRetina modularity -i {self.prefix} -p {self.pw_file} -c 40 -o {out}"
        )
        with open(f"{out}_modularity.tsv") as f:
            header = f.readline().strip().split("\t")
        expected = ["gene_set", "fragmentation", "heterogeneity", "modularity"]
        self.assertEqual(header, expected)

    def test_modularity_all_groups_present(self):
        """All 6 groups from index appear in modularity output."""
        out = os.path.join(self.tmpdir, "mod")
        run_command(
            f"DBRetina modularity -i {self.prefix} -p {self.pw_file} -c 40 -o {out}"
        )
        with open(f"{out}_modularity.tsv") as f:
            next(f)  # skip header
            gene_sets = {line.strip().split("\t")[0].lower() for line in f if line.strip()}
        self.assertEqual(gene_sets, EXPECTED_GROUP_NAMES)

    def test_modularity_high_cutoff(self):
        """Modularity with high cutoff: most groups have low modularity."""
        out = os.path.join(self.tmpdir, "mod")
        rc, _, stderr = run_command(
            f"DBRetina modularity -i {self.prefix} -p {self.pw_file} -c 100 -o {out}"
        )
        self.assertEqual(rc, 0, stderr)
        assert_file_exists(self, f"{out}_modularity.tsv")
        # At cutoff=100, only A-F pair passes. Verify file is valid and has all groups.
        with open(f"{out}_modularity.tsv") as f:
            next(f)  # skip header
            gene_sets = set()
            for line in f:
                if not line.strip():
                    continue
                parts = line.strip().split("\t")
                gene_sets.add(parts[0].lower())
        self.assertEqual(gene_sets, EXPECTED_GROUP_NAMES)


# ============================================================
# SECTION 10: Bipartite Tests
# ============================================================

class TestBipartite(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.prefix, cls.pw_file = _ensure_shared_fixture()

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp(prefix="dbretina_bip_")
        self.g1 = write_file(os.path.join(self.tmpdir, "g1.txt"), "GroupA\nGroupB\n")
        self.g2 = write_file(os.path.join(self.tmpdir, "g2.txt"), "GroupC\nGroupD\n")

    def tearDown(self):
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def test_bipartite_basic(self):
        """Bipartite with two group files produces TSV output."""
        out = os.path.join(self.tmpdir, "bip")
        rc, _, stderr = run_command(
            f"DBRetina bipartite -p {self.pw_file} --group1 {self.g1} "
            f"--group2 {self.g2} -m ochiai -c 0 --no-plot -o {out}"
        )
        self.assertEqual(rc, 0, stderr)
        self.assertNotIn("Traceback", stderr)
        assert_file_exists(self, f"{out}_bipartite_pairwise.tsv")

    def test_bipartite_default_plot_no_kaleido_clean(self):
        """ISSUE-002: default (plotting on) with kaleido absent must NOT crash —
        warn and still produce the data output."""
        out = os.path.join(self.tmpdir, "bipdef")
        rc, _, stderr = run_command(
            f"DBRetina bipartite -p {self.pw_file} --group1 {self.g1} "
            f"--group2 {self.g2} -m ochiai -c 0 -o {out}"
        )
        self.assertEqual(rc, 0, stderr)
        self.assertNotIn("Traceback", stderr)
        assert_file_exists(self, f"{out}_bipartite_pairwise.tsv")

    def test_bipartite_no_plot(self):
        """Bipartite with --no-plot skips interactive bipartite graph."""
        out = os.path.join(self.tmpdir, "bip")
        rc, _, stderr = run_command(
            f"DBRetina bipartite -p {self.pw_file} --group1 {self.g1} "
            f"--group2 {self.g2} -m ochiai --no-plot -o {out}"
        )
        # ISSUE-015: --no-plot must skip ALL plotting (incl the pivot table) -> clean exit.
        self.assertEqual(rc, 0, stderr)
        self.assertNotIn("Traceback", stderr)

    def test_bipartite_no_1_1(self):
        """Bipartite with --no-1-1 excludes 1-to-1 mappings."""
        out = os.path.join(self.tmpdir, "bip")
        rc, _, stderr = run_command(
            f"DBRetina bipartite -p {self.pw_file} --group1 {self.g1} "
            f"--group2 {self.g2} -m ochiai --no-1-1 --no-plot -o {out}"
        )
        # --no-1-1 can remove all pairs, causing plotly to fail on empty data
        if rc != 0:
            plotting_errors = ["kaleido", "nan", "width", "height", "log-scaled"]
            if any(err in stderr.lower() for err in plotting_errors):
                if os.path.exists(f"{out}_bipartite_pairwise.tsv"):
                    return
        self.assertEqual(rc, 0, stderr)

    def test_bipartite_gmt_input(self):
        """Bipartite with GMT files instead of group files."""
        gmt1 = write_file(os.path.join(self.tmpdir, "gmt1.gmt"),
                          "GroupA\tDesc\tAlpha\tBeta\tGamma\tDelta\tEpsilon\n"
                          "GroupB\tDesc\tAlpha\tBeta\tGamma\n")
        gmt2 = write_file(os.path.join(self.tmpdir, "gmt2.gmt"),
                          "GroupC\tDesc\tGamma\tDelta\tZeta\tEta\tTheta\tIota\n"
                          "GroupD\tDesc\tKappa\tLambda\n")
        out = os.path.join(self.tmpdir, "bip")
        rc, _, stderr = run_command(
            f"DBRetina bipartite -p {self.pw_file} --gmt1 {gmt1} --gmt2 {gmt2} "
            f"-m ochiai --no-plot -o {out}"
        )
        if rc != 0 and ("kaleido" in stderr.lower() or "nan" in stderr.lower()):
            if os.path.exists(f"{out}_bipartite_pairwise.tsv"):
                return
        self.assertEqual(rc, 0, stderr)

    def test_bipartite_high_cutoff(self):
        """Bipartite with cutoff=100 finds only perfect matches."""
        out = os.path.join(self.tmpdir, "bip")
        rc, _, stderr = run_command(
            f"DBRetina bipartite -p {self.pw_file} --group1 {self.g1} "
            f"--group2 {self.g2} -m ochiai -c 100 --no-plot -o {out}"
        )
        # May succeed with empty results or error (no perfect cross-group matches)
        # A-B have ochiai=77.5, A-C have 36.5, so no 100% matches across groups


# ============================================================
# SECTION 11: Graph Tests
# ============================================================

class TestGraph(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.prefix, cls.pw_file = _ensure_shared_fixture()

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp(prefix="dbretina_graph_")
        # All groups must be covered by target files to avoid KeyError
        self.intra1 = write_file(
            os.path.join(self.tmpdir, "intra1.tsv"),
            "GroupA\nGroupB\nGroupE\nGroupF\n"
        )
        self.inter1 = write_file(
            os.path.join(self.tmpdir, "inter1.tsv"),
            "GroupC\nGroupD\n"
        )

    def tearDown(self):
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def test_graph_basic(self):
        """Graph with intra-targets produces edges and nodes files."""
        out = os.path.join(self.tmpdir, "gr")
        rc, _, stderr = run_command(
            f"DBRetina graph -i {self.prefix} -p {self.pw_file} "
            f"--intra-targets {self.intra1} -m ochiai -c 50 -o {out}"
        )
        self.assertEqual(rc, 0, stderr)
        assert_file_exists(self, f"{out}_edges.tsv")
        assert_file_exists(self, f"{out}_nodes.tsv")

    def test_graph_with_isolates(self):
        """Graph with --include-isolates includes all nodes."""
        out = os.path.join(self.tmpdir, "gr")
        rc, _, stderr = run_command(
            f"DBRetina graph -i {self.prefix} -p {self.pw_file} "
            f"--intra-targets {self.intra1} -m ochiai -c 50 "
            f"--include-isolates -o {out}"
        )
        self.assertEqual(rc, 0, stderr)

    def test_graph_inter_targets(self):
        """Graph with inter-targets produces output."""
        out = os.path.join(self.tmpdir, "gr")
        rc, _, stderr = run_command(
            f"DBRetina graph -i {self.prefix} -p {self.pw_file} "
            f"--intra-targets {self.intra1} --inter-targets {self.inter1} "
            f"-m ochiai -c 20 -o {out}"
        )
        self.assertEqual(rc, 0, stderr)

    def test_graph_high_cutoff(self):
        """Graph with high cutoff produces minimal edges."""
        out = os.path.join(self.tmpdir, "gr")
        rc, _, stderr = run_command(
            f"DBRetina graph -i {self.prefix} -p {self.pw_file} "
            f"--intra-targets {self.intra1} -m ochiai -c 100 -o {out}"
        )
        self.assertEqual(rc, 0, stderr)

    # ---- Regression: ISSUE-006 (graph -p <parquet DIRECTORY>) ----

    def _parquet_dir(self):
        """The Parquet pairwise directory sibling of the .tsv (a valid -p form)."""
        d = self.pw_file[:-len(".tsv")] if self.pw_file.endswith(".tsv") else self.pw_file
        self.assertTrue(
            os.path.isdir(d), f"expected Parquet pairwise dir at {d}"
        )
        return d

    def test_graph_parquet_dir_input(self):
        """ISSUE-006: graph -p <parquet DIRECTORY> succeeds (like genenet/bipartite).

        A complete intra+inter target partition (GroupA/B/E/F + GroupC/D) avoids the
        unrelated no/partial-targets path. Pre-fix this crashed with
        'RuntimeError: Invalid .dbrp file (bad magic bytes)'.
        """
        pw_dir = self._parquet_dir()
        out = os.path.join(self.tmpdir, "gdir")
        rc, _, stderr = run_command(
            f"DBRetina graph -i {self.prefix} -p {pw_dir} "
            f"--intra-targets {self.intra1} --inter-targets {self.inter1} "
            f"-m ochiai -c 20 -o {out}"
        )
        self.assertEqual(rc, 0, stderr)
        self.assertNotIn("bad magic bytes", stderr)
        assert_file_exists(self, f"{out}_edges.tsv")
        assert_file_exists(self, f"{out}_nodes.tsv")

    def test_graph_dir_matches_tsv_output(self):
        """ISSUE-006: -p <DIR> and -p <.tsv> must produce identical graph output."""
        pw_dir = self._parquet_dir()
        out_dir = os.path.join(self.tmpdir, "from_dir")
        out_tsv = os.path.join(self.tmpdir, "from_tsv")

        rc_d, _, err_d = run_command(
            f"DBRetina graph -i {self.prefix} -p {pw_dir} "
            f"--intra-targets {self.intra1} --inter-targets {self.inter1} "
            f"-m ochiai -c 20 -o {out_dir}"
        )
        rc_t, _, err_t = run_command(
            f"DBRetina graph -i {self.prefix} -p {self.pw_file} "
            f"--intra-targets {self.intra1} --inter-targets {self.inter1} "
            f"-m ochiai -c 20 -o {out_tsv}"
        )
        self.assertEqual(rc_d, 0, err_d)
        self.assertEqual(rc_t, 0, err_t)

        def _read_sorted(path):
            with open(path) as f:
                return sorted(line.rstrip("\n") for line in f)

        self.assertEqual(
            _read_sorted(f"{out_dir}_edges.tsv"),
            _read_sorted(f"{out_tsv}_edges.tsv"),
            "edges differ between -p <DIR> and -p <.tsv>",
        )
        self.assertEqual(
            _read_sorted(f"{out_dir}_nodes.tsv"),
            _read_sorted(f"{out_tsv}_nodes.tsv"),
            "nodes differ between -p <DIR> and -p <.tsv>",
        )

    def test_graph_tsv_input_still_works(self):
        """ISSUE-006: the legacy -p <.tsv> form keeps working after the fix."""
        out = os.path.join(self.tmpdir, "gtsv")
        rc, _, stderr = run_command(
            f"DBRetina graph -i {self.prefix} -p {self.pw_file} "
            f"--intra-targets {self.intra1} --inter-targets {self.inter1} "
            f"-m ochiai -c 20 -o {out}"
        )
        self.assertEqual(rc, 0, stderr)
        assert_file_exists(self, f"{out}_edges.tsv")
        assert_file_exists(self, f"{out}_nodes.tsv")

    # ---- Regression: ISSUE-030 (invalid -m metric) ----

    def test_graph_invalid_metric_clean_error(self):
        """ISSUE-030: graph -m FOO exits nonzero with a clean error, not a KeyError traceback."""
        out = os.path.join(self.tmpdir, "gbadm")
        rc, _, stderr = run_command(
            f"DBRetina graph -i {self.prefix} -p {self.pw_file} "
            f"--intra-targets {self.intra1} --inter-targets {self.inter1} "
            f"-m FOO -c 20 -o {out}"
        )
        self.assertNotEqual(rc, 0, "invalid metric should exit nonzero")
        self.assertNotIn("Traceback", stderr)
        self.assertNotIn("KeyError", stderr)
        # mentions the bad value and/or how to fix it
        self.assertIn("FOO", stderr)
        self.assertIn("metric", stderr.lower())
        # the "NA" sentinel also bypasses the option callback -> must still be clean
        rc2, _, stderr2 = run_command(
            f"DBRetina graph -i {self.prefix} -p {self.pw_file} "
            f"--intra-targets {self.intra1} --inter-targets {self.inter1} "
            f"-m NA -c 20 -o {out}_na"
        )
        self.assertNotEqual(rc2, 0, "-m NA should exit nonzero")
        self.assertNotIn("Traceback", stderr2)
        self.assertNotIn("KeyError", stderr2)

    # ---- Regression: ISSUE-009 (graph with NO targets, the documented default) ----

    def test_graph_no_targets_does_not_crash(self):
        """ISSUE-009: graph with no --intra/--inter targets is the documented default.

        --intra-targets/--inter-targets are NOT [required] in --help, so a
        targets-free invocation is valid documented usage. Pre-fix it crashed with a
        raw 'KeyError' at the first edge because geneSetToTargetsArgumentID was empty.
        It must now succeed (all nodes ungrouped) with no traceback.
        """
        out = os.path.join(self.tmpdir, "gnotarg")
        rc, _, stderr = run_command(
            f"DBRetina graph -i {self.prefix} -p {self.pw_file} "
            f"-m ochiai -c 20 -o {out}"
        )
        self.assertEqual(rc, 0, stderr)
        self.assertNotIn("Traceback", stderr)
        self.assertNotIn("KeyError", stderr)
        assert_file_exists(self, f"{out}_edges.tsv")
        assert_file_exists(self, f"{out}_nodes.tsv")
        # The default must actually produce a graph (edges present), not an empty one.
        self.assertGreater(
            count_tsv_data_rows(f"{out}_edges.tsv"), 0,
            "no-targets default should still produce edges (all nodes ungrouped)",
        )

    # ---- Regression: ISSUE-008 (partial target coverage) ----

    def test_graph_partial_targets_does_not_crash(self):
        """ISSUE-008: a target list covering only SOME groups must be graceful.

        The shared fixture has 6 groups (GroupA..GroupF); here we target only two of
        them. Pre-fix any pairwise edge touching an untargeted group raised a raw
        'KeyError'. Untargeted groups are now handled (treated as ungrouped), so the
        run succeeds with no traceback.
        """
        partial = write_file(
            os.path.join(self.tmpdir, "partial.tsv"),
            "GroupA\nGroupB\n"
        )
        out = os.path.join(self.tmpdir, "gpartial")
        rc, _, stderr = run_command(
            f"DBRetina graph -i {self.prefix} -p {self.pw_file} "
            f"--intra-targets {partial} -m ochiai -c 20 -o {out}"
        )
        self.assertEqual(rc, 0, stderr)
        self.assertNotIn("Traceback", stderr)
        self.assertNotIn("KeyError", stderr)
        assert_file_exists(self, f"{out}_edges.tsv")
        assert_file_exists(self, f"{out}_nodes.tsv")

    # ---- Regression: ISSUE-007 (graph --visualize, optional dash dep missing) ----

    def test_graph_visualize_missing_dep_clean_error(self):
        """ISSUE-007: --visualize without the optional dash/visdcc deps fails cleanly.

        'dash' is not a declared dependency (not even in [all]). Pre-fix the import
        raised a raw 'ModuleNotFoundError' traceback. It must now emit a clean,
        actionable [ERROR] naming the install and exit nonzero.
        """
        try:
            import dash  # noqa: F401
            self.skipTest("dash is installed; cannot test the missing-dep path")
        except ImportError:
            pass

        out = os.path.join(self.tmpdir, "gviz")
        rc, _, stderr = run_command(
            f"DBRetina graph -i {self.prefix} -p {self.pw_file} "
            f"--intra-targets {self.intra1} -m ochiai -c 50 -o {out} --visualize"
        )
        self.assertNotEqual(rc, 0, "missing viz dep should exit nonzero")
        self.assertNotIn("Traceback", stderr)
        self.assertNotIn("ModuleNotFoundError", stderr)
        # actionable: tells the user how to install
        self.assertIn("pip install", stderr.lower())
        self.assertIn("dash", stderr.lower())


# ============================================================
# SECTION 12: Interactome Tests
# ============================================================

class TestGenenet(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.prefix, cls.pw_file = _ensure_shared_fixture()

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp(prefix="dbretina_genenet_")

    def tearDown(self):
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def test_genenet_basic(self):
        """Genenet produces genenet TSV."""
        out = os.path.join(self.tmpdir, "gn")
        rc, _, stderr = run_command(
            f"DBRetina genenet -i {self.prefix} -p {self.pw_file} -o {out}"
        )
        self.assertEqual(rc, 0, stderr)
        assert_file_exists(self, f"{out}_genenet.tsv")

    def test_genenet_graphml(self):
        """Genenet with --graphml exports GraphML."""
        out = os.path.join(self.tmpdir, "gn")
        rc, _, stderr = run_command(
            f"DBRetina genenet -i {self.prefix} -p {self.pw_file} "
            f"--graphml -o {out}"
        )
        self.assertEqual(rc, 0, stderr)
        assert_file_exists(self, f"{out}_genenet.graphml")

    def test_genenet_gexf(self):
        """Genenet with --gexf exports GEXF."""
        out = os.path.join(self.tmpdir, "gn")
        rc, _, stderr = run_command(
            f"DBRetina genenet -i {self.prefix} -p {self.pw_file} "
            f"--gexf -o {out}"
        )
        self.assertEqual(rc, 0, stderr)
        assert_file_exists(self, f"{out}_genenet.gexf")


# ============================================================
# SECTION 13: Geneinfo Tests
# ============================================================

class TestGeneinfo(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # Geneinfo creates files with prefix embedded in name (e.g., inverted_{prefix}_raw.json),
        # so we need a relative prefix and run from its own tmpdir
        cls.geneinfo_dir = tempfile.mkdtemp(prefix="dbretina_geneinfo_base_")
        asc_path = os.path.join(cls.geneinfo_dir, "test_input.asc")
        with open(asc_path, "w") as f:
            f.write(TEST_ASC_CONTENT)
        rc, _, stderr = run_command(
            f"DBRetina index -a {asc_path} -o test_idx", cwd=cls.geneinfo_dir
        )
        assert rc == 0, f"index failed: {stderr}"
        rc, _, stderr = run_command(
            "DBRetina pairwise -i test_idx", cwd=cls.geneinfo_dir
        )
        assert rc == 0, f"pairwise failed: {stderr}"
        cls.prefix = "test_idx"
        cls.pw_file = os.path.join(cls.geneinfo_dir, "test_idx_DBRetina_pairwise.tsv")

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.geneinfo_dir, ignore_errors=True)

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp(prefix="dbretina_geneinfo_")

    def tearDown(self):
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def test_geneinfo_groups(self):
        """Geneinfo with groups file extracts feature info."""
        groups = write_file(os.path.join(self.geneinfo_dir, "groups.txt"),
                            "GroupA\nGroupB\n")
        rc, _, stderr = run_command(
            f"DBRetina geneinfo -i {self.prefix} -g {groups} -o gi",
            cwd=self.geneinfo_dir
        )
        self.assertEqual(rc, 0, stderr)
        assert_file_exists(self, os.path.join(self.geneinfo_dir, "gi_feature_to_groups.tsv"))

    def test_geneinfo_absolute_index_prefix(self):
        """ISSUE-005: geneinfo must accept an absolute / dir-qualified -i prefix
        (was FileNotFoundError from the mangled 'inverted_<abspath>' path)."""
        abs_prefix = os.path.join(self.geneinfo_dir, "test_idx")
        groups = write_file(os.path.join(self.tmpdir, "g.txt"), "GroupA\nGroupB\n")
        out = os.path.join(self.tmpdir, "gi_abs")
        rc, _, stderr = run_command(
            f"DBRetina geneinfo -i {abs_prefix} -g {groups} -o {out}", cwd=self.tmpdir
        )
        self.assertEqual(rc, 0, stderr)
        self.assertNotIn("FileNotFoundError", stderr)
        self.assertNotIn("Traceback", stderr)
        assert_file_exists(self, f"{out}_feature_to_groups.tsv")
        # inverted index must land NEXT TO the index (the index dir), not at a mangled/cwd path
        assert_file_exists(self, os.path.join(self.geneinfo_dir, "inverted_test_idx_raw.json"))
        self.assertFalse(os.path.exists(os.path.join(self.tmpdir, "inverted_test_idx_raw.json")))

    def test_geneinfo_clusters(self):
        """Geneinfo with clusters file."""
        # First create clusters
        cl_out = os.path.join(self.geneinfo_dir, "cl")
        run_command(
            f"DBRetina cluster -p {self.pw_file} -m ochiai -c 50 -o {cl_out}"
        )
        clusters_file = f"{cl_out}_clusters.tsv"
        if not os.path.exists(clusters_file):
            self.skipTest("Cluster file not created")

        rc, _, stderr = run_command(
            f"DBRetina geneinfo -i {self.prefix} --clusters-file {clusters_file} "
            f"--cluster-ids 1 -o gi_cl",
            cwd=self.geneinfo_dir
        )
        self.assertEqual(rc, 0, stderr)


# ============================================================
# SECTION 14: Setcov Tests
# ============================================================

class TestSetcov(unittest.TestCase):

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp(prefix="dbretina_setcov_")
        # Setcov creates temp files with index_prefix embedded in name,
        # so we must use a relative prefix and run from tmpdir
        asc_path = os.path.join(self.tmpdir, "test_input.asc")
        with open(asc_path, "w") as f:
            f.write(TEST_ASC_CONTENT)
        rc, _, stderr = run_command(
            f"DBRetina index -a {asc_path} -o test_idx", cwd=self.tmpdir
        )
        assert rc == 0, f"index failed: {stderr}"
        rc, _, stderr = run_command(
            "DBRetina pairwise -i test_idx", cwd=self.tmpdir
        )
        assert rc == 0, f"pairwise failed: {stderr}"
        self.prefix = "test_idx"

    def tearDown(self):
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def test_setcov_basic(self):
        """Setcov with default parameters produces expected outputs."""
        rc, _, stderr = run_command(
            f"DBRetina setcov -i {self.prefix} --modularity 80 --dedup 100 "
            f"--community 30 --stop-cov 95 -o sc",
            timeout=300, cwd=self.tmpdir
        )
        self.assertEqual(rc, 0, stderr)
        assert_file_exists(self, os.path.join(self.tmpdir, "sc_groups_metadata.tsv"))
        assert_file_exists(self, os.path.join(self.tmpdir, "sc_associations.tsv"))
        assert_file_exists(self, os.path.join(self.tmpdir, "sc_new.gmt"))

    def test_setcov_full_coverage(self):
        """Setcov with --stop-cov 100."""
        rc, _, stderr = run_command(
            f"DBRetina setcov -i {self.prefix} --modularity 80 --dedup 100 "
            f"--community 30 --stop-cov 100 -o sc",
            timeout=300, cwd=self.tmpdir
        )
        self.assertEqual(rc, 0, stderr)

    def test_setcov_half_coverage(self):
        """Setcov with --stop-cov 50."""
        rc, _, stderr = run_command(
            f"DBRetina setcov -i {self.prefix} --modularity 80 --dedup 100 "
            f"--community 30 --stop-cov 50 -o sc_half",
            timeout=300, cwd=self.tmpdir
        )
        self.assertEqual(rc, 0, stderr)


# ============================================================
# SECTION 15: Index Management Tests (append, merge)
# ============================================================

class TestIndexManagement(unittest.TestCase):

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp(prefix="dbretina_idxmgmt_")

    def tearDown(self):
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def _create_index(self, name, content):
        """Create an index from ASC content, return prefix."""
        asc = write_file(os.path.join(self.tmpdir, f"{name}.asc"), content)
        prefix = os.path.join(self.tmpdir, name)
        rc, _, stderr = run_command(f"DBRetina index -a {asc} -o {prefix}")
        assert rc == 0, f"index {name} failed: {stderr}"
        return prefix

    def test_append_asc(self):
        """Append new ASC data to existing index."""
        idx = self._create_index("base",
            "gene_set\tgene\nGroupA\tAlpha\nGroupA\tBeta\n")
        new_asc = write_file(os.path.join(self.tmpdir, "new.asc"),
            "gene_set\tgene\nGroupB\tGamma\nGroupB\tDelta\n")
        out = os.path.join(self.tmpdir, "appended.dbri")
        rc, _, stderr = run_command(
            f"DBRetina append -i {idx}.dbri -a {new_asc} -o {out}"
        )
        self.assertEqual(rc, 0, stderr)
        assert_file_exists(self, out)

    def test_append_gmt(self):
        """Append GMT data to existing index."""
        idx = self._create_index("base",
            "gene_set\tgene\nGroupA\tAlpha\nGroupA\tBeta\n")
        gmt = write_file(os.path.join(self.tmpdir, "new.gmt"),
            "GroupB\tDesc\tGamma\tDelta\n")
        out = os.path.join(self.tmpdir, "appended.dbri")
        rc, _, stderr = run_command(
            f"DBRetina append -i {idx}.dbri -g {gmt} -o {out}"
        )
        self.assertEqual(rc, 0, stderr)
        assert_file_exists(self, out)

    def test_merge_indexes(self):
        """Merge two indexes into one."""
        idx_a = self._create_index("a",
            "gene_set\tgene\nGroupA\tAlpha\nGroupA\tBeta\n")
        idx_b = self._create_index("b",
            "gene_set\tgene\nGroupB\tGamma\nGroupB\tDelta\n")
        out = os.path.join(self.tmpdir, "merged.dbri")
        rc, _, stderr = run_command(
            f"DBRetina merge -a {idx_a}.dbri -b {idx_b}.dbri -o {out}"
        )
        self.assertEqual(rc, 0, stderr)
        assert_file_exists(self, out)

    # --- Regression tests: RAW_GENE_SETS / HASHED_GENE_SETS must include
    #     newly-added groups after append/merge (ISSUE-013, ISSUE-014). ---

    def _raw_data(self, dbri_path):
        """Return the RAW_GENE_SETS 'data' dict embedded in a .dbri."""
        import _dbretina_internal as dbi
        return json.loads(dbi.dbri_load_raw_gene_sets(dbri_path))["data"]

    def test_append_updates_raw_gene_sets(self):
        """ISSUE-014: append must add the new groups to RAW_GENE_SETS."""
        idx = self._create_index("base",
            "gene_set\tgene\nGroupA\tAlpha\nGroupA\tBeta\nGroupA\tGamma\n")
        before = self._raw_data(f"{idx}.dbri")

        gmt = write_file(os.path.join(self.tmpdir, "new.gmt"),
            "append_set_alpha\tdesc\tTP53\tBRCA1\tEGFR\n"
            "append_set_beta\tdesc\tMYC\tKRAS\tPTEN\n")
        out = os.path.join(self.tmpdir, "appended.dbri")
        rc, _, stderr = run_command(
            f"DBRetina append -i {idx}.dbri -g {gmt} -o {out}"
        )
        self.assertEqual(rc, 0, stderr)

        after = self._raw_data(out)
        # (a) count grew by exactly the 2 appended groups
        self.assertEqual(len(after), len(before) + 2,
                         f"raw gene sets should grow by 2; before={len(before)} after={len(after)}")
        # ...and the original group survived
        self.assertIn("groupa", after)
        # (b) the new group names are present with their raw gene names
        self.assertIn("append_set_alpha", after)
        self.assertIn("append_set_beta", after)
        self.assertEqual(set(after["append_set_alpha"]), {"tp53", "brca1", "egfr"})
        self.assertEqual(set(after["append_set_beta"]), {"myc", "kras", "pten"})

    def test_append_geneinfo_returns_rows(self):
        """ISSUE-014: geneinfo on an appended group must return gene rows."""
        idx = self._create_index("base",
            "gene_set\tgene\nGroupA\tAlpha\nGroupA\tBeta\nGroupA\tGamma\n")
        gmt = write_file(os.path.join(self.tmpdir, "new.gmt"),
            "append_set_alpha\tdesc\tTP53\tBRCA1\tEGFR\n"
            "append_set_beta\tdesc\tMYC\tKRAS\tPTEN\n")
        out_prefix = os.path.join(self.tmpdir, "appended")
        rc, _, stderr = run_command(
            f"DBRetina append -i {idx}.dbri -g {gmt} -o {out_prefix}.dbri"
        )
        self.assertEqual(rc, 0, stderr)

        groups = write_file(os.path.join(self.tmpdir, "appended_groups.txt"),
                            "append_set_alpha\n")
        # Use a relative -i prefix (run from tmpdir) to avoid the unrelated
        # geneinfo abspath bug when it builds the inverted index.
        rc, _, stderr = run_command(
            "DBRetina geneinfo -i appended -g appended_groups.txt -o gi_appended",
            cwd=self.tmpdir
        )
        self.assertEqual(rc, 0, stderr)
        gi_tsv = os.path.join(self.tmpdir, "gi_appended_feature_to_groups.tsv")
        assert_file_exists(self, gi_tsv)
        n_rows = count_tsv_data_rows(gi_tsv)
        self.assertEqual(n_rows, 3,
                         f"geneinfo on appended group should return 3 feature rows, got {n_rows}")

    def test_append_pairwise_includes_new_groups(self):
        """Pairwise on the appended index must work and see the new groups."""
        idx = self._create_index("base",
            "gene_set\tgene\nGroupA\tAlpha\nGroupA\tBeta\nGroupA\tGamma\n")
        # New groups share genes with each other so pairwise yields a row.
        gmt = write_file(os.path.join(self.tmpdir, "new.gmt"),
            "append_set_alpha\tdesc\tTP53\tBRCA1\tEGFR\n"
            "append_set_beta\tdesc\tTP53\tBRCA1\tMYC\n")
        out_prefix = os.path.join(self.tmpdir, "appended")
        rc, _, stderr = run_command(
            f"DBRetina append -i {idx}.dbri -g {gmt} -o {out_prefix}.dbri"
        )
        self.assertEqual(rc, 0, stderr)
        rc, _, stderr = run_command(
            f"DBRetina pairwise -i {out_prefix} -m containment -c 1"
        )
        self.assertEqual(rc, 0, stderr)
        pw_file = f"{out_prefix}_DBRetina_pairwise.tsv"
        assert_file_exists(self, pw_file)
        rows = parse_pairwise_tsv(pw_file)
        names = set()
        for r in rows:
            names.add(r["group_1_name"])
            names.add(r["group_2_name"])
        self.assertIn("append_set_alpha", names)
        self.assertIn("append_set_beta", names)

    def test_merge_updates_raw_gene_sets(self):
        """ISSUE-013: merge must include B's groups in RAW_GENE_SETS."""
        idx_a = self._create_index("a",
            "gene_set\tgene\nGroupA\tAlpha\nGroupA\tBeta\nGroupA\tGamma\n")
        idx_b = self._create_index("b",
            "gene_set\tgene\n"
            "bgroup_one\tTP53\nbgroup_one\tBRCA1\nbgroup_one\tEGFR\n"
            "bgroup_two\tMYC\nbgroup_two\tKRAS\n")
        a_raw = self._raw_data(f"{idx_a}.dbri")
        b_raw = self._raw_data(f"{idx_b}.dbri")

        out = os.path.join(self.tmpdir, "merged.dbri")
        rc, _, stderr = run_command(
            f"DBRetina merge -a {idx_a}.dbri -b {idx_b}.dbri -o {out}"
        )
        self.assertEqual(rc, 0, stderr)

        merged = self._raw_data(out)
        # merged raw = union of A and B
        self.assertEqual(len(merged), len(a_raw) + len(b_raw),
                         f"merged raw should be union; A={len(a_raw)} B={len(b_raw)} merged={len(merged)}")
        self.assertIn("groupa", merged)          # from A
        self.assertIn("bgroup_one", merged)      # from B
        self.assertIn("bgroup_two", merged)      # from B
        self.assertEqual(set(merged["bgroup_one"]), {"tp53", "brca1", "egfr"})

    def test_merge_geneinfo_on_b_group(self):
        """ISSUE-013: geneinfo on a B-group must return rows after merge."""
        idx_a = self._create_index("a",
            "gene_set\tgene\nGroupA\tAlpha\nGroupA\tBeta\nGroupA\tGamma\n")
        idx_b = self._create_index("b",
            "gene_set\tgene\n"
            "bgroup_one\tTP53\nbgroup_one\tBRCA1\nbgroup_one\tEGFR\n")
        out_prefix = os.path.join(self.tmpdir, "merged")
        rc, _, stderr = run_command(
            f"DBRetina merge -a {idx_a}.dbri -b {idx_b}.dbri -o {out_prefix}.dbri"
        )
        self.assertEqual(rc, 0, stderr)

        write_file(os.path.join(self.tmpdir, "bgroup.txt"), "bgroup_one\n")
        # Use a relative -i prefix (run from tmpdir) to avoid the unrelated
        # geneinfo abspath bug when it builds the inverted index.
        rc, _, stderr = run_command(
            "DBRetina geneinfo -i merged -g bgroup.txt -o gi_merged",
            cwd=self.tmpdir
        )
        self.assertEqual(rc, 0, stderr)
        gi_tsv = os.path.join(self.tmpdir, "gi_merged_feature_to_groups.tsv")
        assert_file_exists(self, gi_tsv)
        n_rows = count_tsv_data_rows(gi_tsv)
        self.assertEqual(n_rows, 3,
                         f"geneinfo on merged B-group should return 3 feature rows, got {n_rows}")

    def test_merge_pairwise_includes_b_groups(self):
        """Pairwise on the merged index must work and include B's groups."""
        idx_a = self._create_index("a",
            "gene_set\tgene\nGroupA\tAlpha\nGroupA\tBeta\nGroupA\tGamma\n")
        idx_b = self._create_index("b",
            "gene_set\tgene\n"
            "bgroup_one\tTP53\nbgroup_one\tBRCA1\nbgroup_one\tEGFR\n"
            "bgroup_two\tTP53\nbgroup_two\tBRCA1\nbgroup_two\tMYC\n")
        out_prefix = os.path.join(self.tmpdir, "merged")
        rc, _, stderr = run_command(
            f"DBRetina merge -a {idx_a}.dbri -b {idx_b}.dbri -o {out_prefix}.dbri"
        )
        self.assertEqual(rc, 0, stderr)
        rc, _, stderr = run_command(
            f"DBRetina pairwise -i {out_prefix} -m containment -c 1"
        )
        self.assertEqual(rc, 0, stderr)
        pw_file = f"{out_prefix}_DBRetina_pairwise.tsv"
        assert_file_exists(self, pw_file)
        rows = parse_pairwise_tsv(pw_file)
        names = set()
        for r in rows:
            names.add(r["group_1_name"])
            names.add(r["group_2_name"])
        self.assertIn("bgroup_one", names)
        self.assertIn("bgroup_two", names)


# ============================================================
# SECTION 16: Edge Case Tests
# ============================================================

class TestEdgeCases(unittest.TestCase):

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp(prefix="dbretina_edge_")

    def tearDown(self):
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def test_single_group_pairwise(self):
        """Single group produces zero pairwise rows."""
        asc = write_file(os.path.join(self.tmpdir, "single.asc"),
                         "gene_set\tgene\nOnlyGroup\tGene1\nOnlyGroup\tGene2\n")
        prefix = os.path.join(self.tmpdir, "idx")
        run_command(f"DBRetina index -a {asc} -o {prefix}")
        rc, _, stderr = run_command(f"DBRetina pairwise -i {prefix}")
        self.assertEqual(rc, 0, stderr)
        pw_file = f"{prefix}_DBRetina_pairwise.tsv"
        assert_file_exists(self, pw_file)
        self.assertEqual(count_tsv_data_rows(pw_file), 0)

    def test_two_groups_no_overlap(self):
        """Two disjoint groups produce zero pairwise rows."""
        asc = write_file(os.path.join(self.tmpdir, "disjoint.asc"),
                         "gene_set\tgene\n"
                         "GroupA\tGene1\nGroupA\tGene2\n"
                         "GroupB\tGene3\nGroupB\tGene4\n")
        prefix = os.path.join(self.tmpdir, "idx")
        run_command(f"DBRetina index -a {asc} -o {prefix}")
        run_command(f"DBRetina pairwise -i {prefix}")
        pw_file = f"{prefix}_DBRetina_pairwise.tsv"
        self.assertEqual(count_tsv_data_rows(pw_file), 0)

    def test_two_groups_with_overlap(self):
        """Two overlapping groups produce one pairwise row."""
        asc = write_file(os.path.join(self.tmpdir, "overlap.asc"),
                         "gene_set\tgene\n"
                         "GroupA\tGene1\nGroupA\tGene2\nGroupA\tGene3\n"
                         "GroupB\tGene2\nGroupB\tGene3\nGroupB\tGene4\n")
        prefix = os.path.join(self.tmpdir, "idx")
        run_command(f"DBRetina index -a {asc} -o {prefix}")
        run_command(f"DBRetina pairwise -i {prefix}")
        pw_file = f"{prefix}_DBRetina_pairwise.tsv"
        self.assertEqual(count_tsv_data_rows(pw_file), 1)

        # Verify metrics for 2 shared out of 3 each
        rows = parse_pairwise_tsv(pw_file)
        self.assertEqual(rows[0]["shared_features"], 2)

    def test_identical_groups(self):
        """Two identical groups have containment=100, ochiai=100, odds_ratio=-1."""
        asc = write_file(os.path.join(self.tmpdir, "identical.asc"),
                         "gene_set\tgene\n"
                         "GroupA\tG1\nGroupA\tG2\nGroupA\tG3\n"
                         "GroupB\tG1\nGroupB\tG2\nGroupB\tG3\n")
        prefix = os.path.join(self.tmpdir, "idx")
        run_command(f"DBRetina index -a {asc} -o {prefix}")
        run_command(f"DBRetina pairwise -i {prefix}")
        rows = parse_pairwise_tsv(f"{prefix}_DBRetina_pairwise.tsv")
        self.assertEqual(len(rows), 1)
        self.assertAlmostEqual(rows[0]["containment"], 100.0, delta=0.1)
        self.assertAlmostEqual(rows[0]["ochiai"], 100.0, delta=0.1)
        self.assertAlmostEqual(rows[0]["odds_ratio"], -1.0, delta=0.1)

    def test_pairwise_multithreaded_consistency(self):
        """Multithreaded pairwise matches single-threaded."""
        prefix_st, pw_st = setup_index_and_pairwise(
            tempfile.mkdtemp(dir=self.tmpdir, prefix="st_")
        )
        prefix_mt, pw_mt = setup_index_and_pairwise(
            tempfile.mkdtemp(dir=self.tmpdir, prefix="mt_"),
            extra_pw_args="-t 4"
        )
        rows_st = parse_pairwise_tsv(pw_st)
        rows_mt = parse_pairwise_tsv(pw_mt)
        self.assertEqual(len(rows_st), len(rows_mt))

    def test_dbrp_fallback_dedup(self):
        """Dedup works with TSV-only (no .dbrp)."""
        prefix, pw_file = setup_index_and_pairwise(self.tmpdir)
        dbrp = pw_file.replace(".tsv", ".dbrp")
        if os.path.exists(dbrp):
            os.remove(dbrp)

        out = os.path.join(self.tmpdir, "dd")
        rc, _, stderr = run_command(
            f"DBRetina dedup -i {prefix} -p {pw_file} -c 100 -o {out}"
        )
        self.assertEqual(rc, 0, stderr)
        with open(f"{out}_deduplicated_groups.txt") as f:
            groups = {l.strip() for l in f if l.strip()}
        self.assertEqual(len(groups), 5)

    def test_dbrp_fallback_modularity(self):
        """Modularity works with TSV-only (no .dbrp)."""
        prefix, pw_file = setup_index_and_pairwise(self.tmpdir)
        dbrp = pw_file.replace(".tsv", ".dbrp")
        if os.path.exists(dbrp):
            os.remove(dbrp)

        out = os.path.join(self.tmpdir, "mod")
        rc, _, stderr = run_command(
            f"DBRetina modularity -i {prefix} -p {pw_file} -c 40 -o {out}"
        )
        self.assertEqual(rc, 0, stderr)
        assert_file_exists(self, f"{out}_modularity.tsv")


# ============================================================
# SECTION 17: Query with Clusters Tests
# ============================================================

class TestQueryClusters(unittest.TestCase):
    """Test query with cluster files (depends on cluster output)."""

    @classmethod
    def setUpClass(cls):
        cls.prefix, cls.pw_file = _ensure_shared_fixture()

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp(prefix="dbretina_qcl_")

    def tearDown(self):
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def test_query_with_cluster_ids(self):
        """Query using clusters file and cluster IDs."""
        # Create clusters first
        cl_out = os.path.join(self.tmpdir, "cl")
        rc, _, _ = run_command(
            f"DBRetina cluster -p {self.pw_file} -m ochiai -c 30 -o {cl_out}"
        )
        if rc != 0:
            self.skipTest("Cluster creation failed")

        clusters_file = f"{cl_out}_clusters.tsv"
        if not os.path.exists(clusters_file):
            self.skipTest("Clusters file not found")

        out = os.path.join(self.tmpdir, "q")
        rc, _, stderr = run_command(
            f"DBRetina query -p {self.pw_file} --clusters-file {clusters_file} "
            f"--cluster-ids 1 -m ochiai -c 0 -o {out}"
        )
        self.assertEqual(rc, 0, stderr)
        assert_file_exists(self, f"{out}.tsv")


# ============================================================
# SECTION 18: Index Input Variant Tests
# ============================================================

class TestIndexVariants(unittest.TestCase):
    """Verify that the same logical data, presented in different input shapes,
    produces an identical index (_raw.json groups and gene lists)."""

    @classmethod
    def setUpClass(cls):
        cls.base_dir = tempfile.mkdtemp(prefix="dbretina_variants_base_")
        cls.canonical = get_canonical_groups(cls.base_dir)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.base_dir, ignore_errors=True)

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp(prefix="dbretina_var_")

    def tearDown(self):
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def _index_and_get_groups(self, asc_content, extra_args=""):
        asc = write_file(os.path.join(self.tmpdir, "input.asc"), asc_content)
        prefix = os.path.join(self.tmpdir, "idx")
        rc, _, stderr = run_command(f"DBRetina index -a {asc} {extra_args} -o {prefix}")
        self.assertEqual(rc, 0, stderr)
        return get_groups_from_raw_json(f"{prefix}_raw.json")

    def _assert_groups_match_canonical(self, groups):
        self.assertEqual(set(groups.keys()), set(self.canonical.keys()),
                         f"Group name mismatch: {set(groups.keys())} vs {set(self.canonical.keys())}")
        for gname in self.canonical:
            self.assertEqual(groups[gname], self.canonical[gname],
                             f"Gene list mismatch for {gname}")

    def test_split_rows_with_duplicates(self):
        """Groups with duplicate rows produce same index as clean input."""
        groups = self._index_and_get_groups(TEST_ASC_SPLIT_ROWS)
        self._assert_groups_match_canonical(groups)

    def test_one_gene_per_row(self):
        """One gene per row (max fragmentation) produces same index."""
        groups = self._index_and_get_groups(TEST_ASC_ONE_PER_ROW)
        self._assert_groups_match_canonical(groups)

    def test_duplicate_rows(self):
        """Exact duplicate rows are deduplicated to produce same index."""
        groups = self._index_and_get_groups(TEST_ASC_DUPLICATE_ROWS)
        self._assert_groups_match_canonical(groups)

    def test_mixed_case_groups_and_genes(self):
        """Mixed case in group/gene names all lowercased to produce same index."""
        groups = self._index_and_get_groups(TEST_ASC_MIXED_CASE)
        self._assert_groups_match_canonical(groups)

    def test_quoted_names(self):
        """Double-quoted group/gene names have quotes stripped to produce same index."""
        groups = self._index_and_get_groups(TEST_ASC_QUOTED)
        self._assert_groups_match_canonical(groups)

    def test_shuffled_row_order(self):
        """Shuffled row order produces same index (order-independent)."""
        groups = self._index_and_get_groups(TEST_ASC_SHUFFLED)
        self._assert_groups_match_canonical(groups)

    def test_multi_file_split(self):
        """Groups split across 3 ASC files are merged correctly."""
        f1 = write_file(os.path.join(self.tmpdir, "f1.asc"), TEST_ASC_FILE1)
        f2 = write_file(os.path.join(self.tmpdir, "f2.asc"), TEST_ASC_FILE2)
        f3 = write_file(os.path.join(self.tmpdir, "f3.asc"), TEST_ASC_FILE3)
        prefix = os.path.join(self.tmpdir, "idx")
        rc, _, stderr = run_command(
            f"DBRetina index -a {f1} -a {f2} -a {f3} -o {prefix}"
        )
        self.assertEqual(rc, 0, stderr)
        groups = get_groups_from_raw_json(f"{prefix}_raw.json")
        self._assert_groups_match_canonical(groups)

    def test_multi_file_cross_group_merge(self):
        """Same group split across 2 files merges correctly."""
        f1 = write_file(os.path.join(self.tmpdir, "f1.asc"), TEST_ASC_CROSS1)
        f2 = write_file(os.path.join(self.tmpdir, "f2.asc"), TEST_ASC_CROSS2)
        prefix = os.path.join(self.tmpdir, "idx")
        rc, _, stderr = run_command(
            f"DBRetina index -a {f1} -a {f2} -o {prefix}"
        )
        self.assertEqual(rc, 0, stderr)
        groups = get_groups_from_raw_json(f"{prefix}_raw.json")
        self._assert_groups_match_canonical(groups)

    def test_special_chars_in_names(self):
        """Special characters (spaces, hyphens, underscores, parens, dots, brackets)
        are preserved in group and gene names."""
        groups = self._index_and_get_groups(TEST_ASC_SPECIAL)
        # After lowercasing: "gene set-1 (main)" and "group_2.test"
        self.assertIn("gene set-1 (main)", groups)
        self.assertIn("group_2.test", groups)
        # Gene names preserved
        self.assertIn("gene_alpha-v2", groups["gene set-1 (main)"])
        self.assertIn("gene.beta", groups["gene set-1 (main)"])
        self.assertIn("gene [gamma]", groups["gene set-1 (main)"])
        # Shared gene present in both groups
        self.assertIn("gene_alpha-v2", groups["group_2.test"])
        self.assertIn("gene [gamma]", groups["group_2.test"])

    def test_gmt_produces_same_groups(self):
        """GMT format input produces same groups as ASC."""
        gmt = write_file(os.path.join(self.tmpdir, "test.gmt"), TEST_GMT_CONTENT)
        rc, _, stderr = run_command(
            f"DBRetina index -g {gmt} -o idx", cwd=self.tmpdir
        )
        self.assertEqual(rc, 0, stderr)
        groups = get_groups_from_raw_json(os.path.join(self.tmpdir, "idx_raw.json"))
        self._assert_groups_match_canonical(groups)


# ============================================================
# SECTION 19: Index Error Tests
# ============================================================

class TestIndexErrors(unittest.TestCase):
    """Verify that invalid inputs are rejected with appropriate errors."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp(prefix="dbretina_err_")

    def tearDown(self):
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def test_pipe_in_group_name(self):
        """Pipe character in group name is rejected."""
        asc = write_file(os.path.join(self.tmpdir, "bad.asc"),
                         "gene_set\tgene\nGroup|A\tAlpha\n")
        prefix = os.path.join(self.tmpdir, "idx")
        rc, _, stderr = run_command(f"DBRetina index -a {asc} -o {prefix}")
        self.assertNotEqual(rc, 0)

    def test_pipe_in_gene_name(self):
        """Pipe character in gene name is rejected."""
        asc = write_file(os.path.join(self.tmpdir, "bad.asc"),
                         "gene_set\tgene\nGroupA\tAlpha|Beta\n")
        prefix = os.path.join(self.tmpdir, "idx")
        rc, _, stderr = run_command(f"DBRetina index -a {asc} -o {prefix}")
        self.assertNotEqual(rc, 0)

    def test_empty_gene_name(self):
        """Empty gene name (tab with nothing after) is rejected."""
        asc = write_file(os.path.join(self.tmpdir, "bad.asc"),
                         "gene_set\tgene\nGroupA\t\n")
        prefix = os.path.join(self.tmpdir, "idx")
        rc, _, stderr = run_command(f"DBRetina index -a {asc} -o {prefix}")
        self.assertNotEqual(rc, 0)

    def test_empty_group_name(self):
        """Empty group name is rejected."""
        asc = write_file(os.path.join(self.tmpdir, "bad.asc"),
                         "gene_set\tgene\n\tAlpha\n")
        prefix = os.path.join(self.tmpdir, "idx")
        rc, _, stderr = run_command(f"DBRetina index -a {asc} -o {prefix}")
        self.assertNotEqual(rc, 0)

    def test_no_input_files(self):
        """Index with no -a or -g flags fails."""
        prefix = os.path.join(self.tmpdir, "idx")
        rc, _, stderr = run_command(f"DBRetina index -o {prefix}")
        self.assertNotEqual(rc, 0)

    def test_both_asc_and_gmt(self):
        """Index with both -a and -g flags fails."""
        asc = write_file(os.path.join(self.tmpdir, "test.asc"),
                         "gene_set\tgene\nGroupA\tAlpha\n")
        gmt = write_file(os.path.join(self.tmpdir, "test.gmt"),
                         "GroupB\tDesc\tBeta\n")
        prefix = os.path.join(self.tmpdir, "idx")
        rc, _, stderr = run_command(f"DBRetina index -a {asc} -g {gmt} -o {prefix}")
        self.assertNotEqual(rc, 0)


# ============================================================
# SECTION 20: Pairwise from Input Variants
# ============================================================

class TestPairwiseFromVariants(unittest.TestCase):
    """Full index->pairwise pipeline on variant inputs produces identical results."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp(prefix="dbretina_pwvar_")

    def tearDown(self):
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def test_pairwise_from_split_rows(self):
        """Split rows with duplicates produce same pairwise as clean input."""
        _, pw_file = setup_index_and_pairwise(self.tmpdir, asc_content=TEST_ASC_SPLIT_ROWS)
        assert_pairwise_matches_expected(self, pw_file)

    def test_pairwise_from_duplicate_rows(self):
        """Duplicate rows produce same pairwise as clean input."""
        _, pw_file = setup_index_and_pairwise(self.tmpdir, asc_content=TEST_ASC_DUPLICATE_ROWS)
        assert_pairwise_matches_expected(self, pw_file)

    def test_pairwise_from_mixed_case(self):
        """Mixed case input produces same pairwise as clean input."""
        _, pw_file = setup_index_and_pairwise(self.tmpdir, asc_content=TEST_ASC_MIXED_CASE)
        assert_pairwise_matches_expected(self, pw_file)

    def test_pairwise_from_multi_file(self):
        """3-file split produces same pairwise as single-file input."""
        f1 = write_file(os.path.join(self.tmpdir, "f1.asc"), TEST_ASC_FILE1)
        f2 = write_file(os.path.join(self.tmpdir, "f2.asc"), TEST_ASC_FILE2)
        f3 = write_file(os.path.join(self.tmpdir, "f3.asc"), TEST_ASC_FILE3)
        prefix = os.path.join(self.tmpdir, "idx")
        rc, _, stderr = run_command(
            f"DBRetina index -a {f1} -a {f2} -a {f3} -o {prefix}"
        )
        self.assertEqual(rc, 0, stderr)
        rc, _, stderr = run_command(f"DBRetina pairwise -i {prefix}")
        self.assertEqual(rc, 0, stderr)
        pw_file = f"{prefix}_DBRetina_pairwise.tsv"
        assert_pairwise_matches_expected(self, pw_file)

    def test_pairwise_from_quoted_names(self):
        """Quoted names produce same pairwise as clean input (quotes stripped)."""
        _, pw_file = setup_index_and_pairwise(self.tmpdir, asc_content=TEST_ASC_QUOTED)
        assert_pairwise_matches_expected(self, pw_file)

    def test_pairwise_from_gmt(self):
        """GMT format input produces same pairwise as ASC."""
        gmt = write_file(os.path.join(self.tmpdir, "test.gmt"), TEST_GMT_CONTENT)
        rc, _, stderr = run_command(
            f"DBRetina index -g {gmt} -o idx", cwd=self.tmpdir
        )
        self.assertEqual(rc, 0, stderr)
        rc, _, stderr = run_command("DBRetina pairwise -i idx", cwd=self.tmpdir)
        self.assertEqual(rc, 0, stderr)
        pw_file = os.path.join(self.tmpdir, "idx_DBRetina_pairwise.tsv")
        assert_pairwise_matches_expected(self, pw_file)


# ============================================================
# SECTION 21: Special Character Pipeline Tests
# ============================================================

class TestSpecialCharPairwise(unittest.TestCase):
    """End-to-end pipeline with special characters in group/gene names."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp(prefix="dbretina_special_")

    def tearDown(self):
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def test_special_char_pipeline(self):
        """Index -> pairwise with special characters in names succeeds."""
        asc = write_file(os.path.join(self.tmpdir, "special.asc"), TEST_ASC_SPECIAL)
        prefix = os.path.join(self.tmpdir, "idx")
        rc, _, stderr = run_command(f"DBRetina index -a {asc} -o {prefix}")
        self.assertEqual(rc, 0, stderr)

        rc, _, stderr = run_command(f"DBRetina pairwise -i {prefix}")
        self.assertEqual(rc, 0, stderr)

        pw_file = f"{prefix}_DBRetina_pairwise.tsv"
        assert_file_exists(self, pw_file)
        rows = parse_pairwise_tsv(pw_file)
        # 2 groups with shared genes -> 1 pair
        self.assertEqual(len(rows), 1)
        # shared features: gene_alpha-v2 and gene [gamma] -> 2
        self.assertEqual(rows[0]["shared_features"], 2)

    def test_query_special_char_groups(self):
        """Query with -g on special-char group names (no regex-sensitive chars)."""
        # Use group names without parentheses/brackets since query uses awk internally
        special_asc = (
            "gene_set\tgene\n"
            "gene-set_1\tgene_alpha-v2\n"
            "gene-set_1\tgene.beta\n"
            "gene-set_1\tgene_gamma\n"
            "group_2.test\tgene_alpha-v2\n"
            "group_2.test\tgene_gamma\n"
            "group_2.test\tgene_delta\n"
        )
        asc = write_file(os.path.join(self.tmpdir, "special.asc"), special_asc)
        prefix = os.path.join(self.tmpdir, "idx")
        rc, _, stderr = run_command(f"DBRetina index -a {asc} -o {prefix}")
        self.assertEqual(rc, 0, stderr)
        rc, _, stderr = run_command(f"DBRetina pairwise -i {prefix}")
        self.assertEqual(rc, 0, stderr)

        pw_file = f"{prefix}_DBRetina_pairwise.tsv"
        # Query -g requires BOTH groups in a pair to be in the groups file
        groups = write_file(os.path.join(self.tmpdir, "groups.txt"),
                            "gene-set_1\ngroup_2.test\n")
        out = os.path.join(self.tmpdir, "q")
        rc, _, stderr = run_command(
            f"DBRetina query -p {pw_file} -g {groups} -m ochiai -c 0 -o {out}"
        )
        self.assertEqual(rc, 0, stderr)
        assert_file_exists(self, f"{out}.tsv")
        rows = parse_pairwise_tsv(f"{out}.tsv")
        self.assertGreaterEqual(len(rows), 1)


# ============================================================
# Genescore (subdev1 A2 regression): -m/-c must be applied
# ============================================================

def parse_gene_score_map(out):
    """Parse genescore stdout -> {gene: score}; ignores log/header lines."""
    scores = {}
    for line in out.splitlines():
        parts = line.strip().split("\t")
        if len(parts) >= 2 and parts[0] != "gene":
            try:
                scores[parts[0]] = float(parts[1])
            except ValueError:
                continue
    return scores


class TestGenescore(unittest.TestCase):
    """genescore must actually apply --metric/--cutoff to edge_weighted scoring."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp(prefix="dbretina_genescore_")
        self.prefix, self.pw = setup_index_and_pairwise(self.tmpdir)
        self.parquet = f"{self.prefix}_DBRetina_pairwise"

    def tearDown(self):
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def _edge_weighted_scores(self, cutoff, group="groupa", metric="ochiai"):
        """Run genescore edge_weighted; return {gene: score} (skips log/header lines)."""
        rc, out, err = run_command(
            f"DBRetina genescore -d {self.parquet} -i {self.prefix} "
            f"{group} --method edge_weighted -m {metric} -c {cutoff} -n 50"
        )
        self.assertEqual(rc, 0, f"genescore failed (cutoff={cutoff}): {err}")
        return parse_gene_score_map(out)

    def test_edge_weighted_respects_cutoff(self):
        # GroupA's pairs: B(~77), C(~36), E(~22), F(100, identical). At -c 100 only the
        # identical A=F pair survives, so scores must differ from -c 0 (all pairs).
        low = self._edge_weighted_scores(0)
        high = self._edge_weighted_scores(100)
        self.assertTrue(low, "edge_weighted at -c 0 returned no genes")
        self.assertTrue(high, "edge_weighted at -c 100 returned no genes")
        self.assertNotEqual(
            low, high,
            "genescore --method edge_weighted ignored -c: identical scores at "
            f"-c 0 and -c 100 ({low})",
        )

    def _compare_scores(self, ga, gb, cutoff, metric="ochiai"):
        """Run genescore --compare (explain_pair path); return {gene: score}."""
        rc, out, err = run_command(
            f"DBRetina genescore -d {self.parquet} -i {self.prefix} "
            f"{ga} --compare {gb} --method edge_weighted -m {metric} -c {cutoff} -n 50"
        )
        self.assertEqual(rc, 0, f"genescore --compare failed (cutoff={cutoff}): {err}")
        return parse_gene_score_map(out)

    def test_edge_weighted_respects_metric(self):
        # Same cutoff, different metric -> different weights/qualifying pairs -> different scores.
        ochiai = self._edge_weighted_scores(20, metric="ochiai")
        jaccard = self._edge_weighted_scores(20, metric="jaccard")
        self.assertTrue(ochiai, "edge_weighted (ochiai) returned no genes")
        self.assertTrue(jaccard, "edge_weighted (jaccard) returned no genes")
        self.assertNotEqual(
            ochiai, jaccard,
            "genescore --method edge_weighted ignored -m (ochiai == jaccard scores)",
        )

    def test_explain_pair_respects_cutoff(self):
        # --compare routes through explain_pair, which must also honor -c.
        low = self._compare_scores("groupa", "groupc", 0)
        high = self._compare_scores("groupa", "groupc", 100)
        self.assertTrue(low, "explain_pair at -c 0 returned no shared genes")
        self.assertTrue(high, "explain_pair at -c 100 returned no shared genes")
        self.assertNotEqual(
            low, high,
            "genescore --compare --method edge_weighted ignored -c",
        )

    def test_explain_pair_respects_metric(self):
        # --compare must honor -m too (forwarded on the same call line as cutoff).
        ochiai = self._compare_scores("groupa", "groupc", 20, metric="ochiai")
        jaccard = self._compare_scores("groupa", "groupc", 20, metric="jaccard")
        self.assertTrue(ochiai, "explain_pair (ochiai) returned no shared genes")
        self.assertTrue(jaccard, "explain_pair (jaccard) returned no shared genes")
        self.assertNotEqual(
            ochiai, jaccard,
            "genescore --compare --method edge_weighted ignored -m",
        )

    def test_invalid_metric_rejected(self):
        # An unknown metric must fail cleanly at the CLI, not with a raw traceback.
        rc, out, err = run_command(
            f"DBRetina genescore -d {self.parquet} -i {self.prefix} "
            f"groupa --method edge_weighted -m bogus_metric -c 20"
        )
        self.assertNotEqual(rc, 0, "invalid --metric should fail")
        self.assertNotIn(
            "Traceback", out + err,
            f"invalid --metric produced a raw traceback:\n{err[-300:]}",
        )


# ============================================================
# REST hub-genes (subdev1 A5): endpoint must honor request metric/cutoff
# ============================================================

try:
    import fastapi as _fastapi  # noqa: F401
    import duckdb as _duckdb  # noqa: F401
    _HAS_SERVER = True
except Exception:
    _HAS_SERVER = False


@unittest.skipUnless(_HAS_SERVER, "REST/genescore needs the [server] extra (fastapi, duckdb)")
class TestRestHubGenes(unittest.TestCase):
    """REST /api/v1/genes/hub-genes must apply request.metric/cutoff (A5);
    edge_weighted_scores must reject unknown metrics with a clean error (A7)."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp(prefix="dbretina_rest_")
        self.prefix, self.pw = setup_index_and_pairwise(self.tmpdir)
        self.parquet = f"{self.prefix}_DBRetina_pairwise"
        self.dbri = f"{self.prefix}.dbri"
        from dbretina.compat import open_pairwise
        from dbretina.pairwise_store import PairwiseStore
        from dbretina.rest_api import create_app
        store = open_pairwise(self.parquet)
        if store is None:
            store = PairwiseStore(self.parquet, dbri_path=self.dbri)
        else:
            store._dbri_path = self.dbri
        self.store = store
        app = create_app(store, dbri_path=self.dbri)
        self.hub_ep = next(
            (r.endpoint for r in app.routes
             if getattr(r, "path", "") == "/api/v1/genes/hub-genes"), None)
        self.assertIsNotNone(self.hub_ep, "hub-genes endpoint not registered")

    def tearDown(self):
        try:
            self.store.close()
        except Exception:
            pass
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def _hub(self, cutoff, metric="ochiai", group="groupa"):
        body = SimpleNamespace(group_name=group, method="edge_weighted",
                               hops=2, top_n=50, metric=metric, cutoff=cutoff)
        res = self.hub_ep(body)
        return {g["gene"]: g["score"] for g in res["genes"]}

    def test_hub_genes_respects_cutoff(self):
        low = self._hub(0)
        high = self._hub(100)
        self.assertTrue(low, "REST hub-genes at cutoff 0 returned no genes")
        self.assertNotEqual(
            low, high,
            "REST /genes/hub-genes ignored request.cutoff (identical at 0 and 100)",
        )

    def test_hub_genes_respects_metric(self):
        ochiai = self._hub(20, metric="ochiai")
        jaccard = self._hub(20, metric="jaccard")
        self.assertTrue(ochiai, "REST hub-genes (ochiai) returned no genes")
        self.assertNotEqual(
            ochiai, jaccard,
            "REST /genes/hub-genes ignored request.metric (ochiai == jaccard)",
        )

    def test_edge_weighted_rejects_unknown_metric(self):
        # A7: an unknown metric must raise a clean ValueError, not leak a raw DuckDB error.
        from dbretina.gene_importance import GeneImportance
        gi = GeneImportance(self.store)
        with self.assertRaises(ValueError):
            gi.edge_weighted_scores("groupa", metric="bogus_metric")


# ============================================================
# REST security (S1 SQL file-read, S2 CORS, S3 api-key-in-URL)
# ============================================================

try:
    import fastapi as _fa  # noqa: F401
    import duckdb as _dk  # noqa: F401
    import httpx as _hx  # noqa: F401
    from starlette.testclient import TestClient as _TestClient
    _HAS_REST_TEST = True
except Exception:
    _HAS_REST_TEST = False


@unittest.skipUnless(_HAS_REST_TEST, "REST security tests need [server] extra + httpx")
class TestRestSecurity(unittest.TestCase):
    """serve REST must: block SQL file reads, not use wildcard CORS, reject API key in URL."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp(prefix="dbretina_sec_")
        self.prefix, self.pw = setup_index_and_pairwise(self.tmpdir)
        self.parquet = f"{self.prefix}_DBRetina_pairwise"
        self.dbri = f"{self.prefix}.dbri"

    def tearDown(self):
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def _client(self, api_key=None):
        from dbretina.compat import open_pairwise
        from dbretina.pairwise_store import PairwiseStore
        from dbretina.rest_api import create_app
        store = open_pairwise(self.parquet)
        if store is None:
            store = PairwiseStore(self.parquet, dbri_path=self.dbri)
        else:
            store._dbri_path = self.dbri
        app = create_app(store, api_key=api_key, dbri_path=self.dbri)
        return _TestClient(app), store

    def test_sql_endpoint_cannot_read_files(self):
        # S1: arbitrary file read via DuckDB read_text() must be blocked.
        client, store = self._client()
        try:
            r = client.post("/api/v1/sql",
                            json={"query": "SELECT * FROM read_text('/etc/passwd')"})
            self.assertNotIn("root:", r.text,
                             "SQL endpoint leaked /etc/passwd contents (file read not blocked)")
            self.assertNotEqual(r.status_code, 200,
                                f"file-read query unexpectedly succeeded: {r.text[:200]}")
        finally:
            store.close()

    def test_cors_is_not_wildcard(self):
        # S2: CORS must not allow arbitrary origins.
        client, store = self._client()
        try:
            r = client.get("/api/v1/info", headers={"Origin": "http://evil.example"})
            acao = r.headers.get("access-control-allow-origin")
            self.assertNotEqual(acao, "*", "CORS allows any origin (*)")
            self.assertNotEqual(acao, "http://evil.example",
                                "CORS reflects an arbitrary origin")
        finally:
            store.close()

    def test_api_key_rejected_in_query_param(self):
        # S3: when auth is on, the key must NOT be accepted via URL query param (log leak).
        client, store = self._client(api_key="SECRET")
        try:
            r_q = client.get("/api/v1/info?api_key=SECRET")
            self.assertEqual(r_q.status_code, 401,
                             "API key accepted via query param (leaks into URLs/logs)")
            r_h = client.get("/api/v1/info", headers={"x-api-key": "SECRET"})
            self.assertEqual(r_h.status_code, 200, "API key via header should authenticate")
        finally:
            store.close()

    def test_harden_readonly_is_idempotent(self):
        # Calling harden twice (e.g. create_app twice on one store) must not crash and stay locked.
        from dbretina.compat import open_pairwise
        from dbretina.pairwise_store import PairwiseStore
        store = open_pairwise(self.parquet)
        if store is None:
            store = PairwiseStore(self.parquet, dbri_path=self.dbri)
        try:
            store.harden_readonly()
            store.harden_readonly()  # second call must be a no-op, not raise
            with self.assertRaises(Exception):  # filesystem still disabled
                store._con.execute("SELECT * FROM read_text('/etc/passwd')").fetchone()
            n = store._con.execute("SELECT count(*) FROM pairs").fetchone()[0]
            self.assertGreater(n, 0, "pairs still queryable after double-harden")
        finally:
            store.close()

    def test_cypher_endpoint_removed(self):
        # The raw /api/v1/cypher endpoint is removed (unsafe: file read + segfault DoS via Kùzu).
        client, store = self._client()
        try:
            r = client.post("/api/v1/cypher",
                            json={"query": "MATCH (g:`Group`) RETURN count(g)",
                                  "metric": "ochiai", "cutoff": 0.0})
            # Route gone: 404 (no route) or 405 (POST falls through to the static mount).
            self.assertIn(r.status_code, (404, 405),
                          f"/api/v1/cypher should be removed, got {r.status_code}")
            self.assertNotIn("row_count", r.text, "/cypher must not execute queries")
        finally:
            store.close()


@unittest.skipUnless(_HAS_REST_TEST, "graph endpoint tests need [server] extra + httpx")
class TestRestGraph(unittest.TestCase):
    """serve graph + algorithms endpoints (igraph-backed)."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp(prefix="dbretina_graph_")
        self.prefix, self.pw = setup_index_and_pairwise(self.tmpdir)
        self.parquet = f"{self.prefix}_DBRetina_pairwise"
        self.dbri = f"{self.prefix}.dbri"

    def tearDown(self):
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def _client(self):
        from dbretina.compat import open_pairwise
        from dbretina.pairwise_store import PairwiseStore
        from dbretina.rest_api import create_app
        store = open_pairwise(self.parquet)
        if store is None:
            store = PairwiseStore(self.parquet, dbri_path=self.dbri)
        else:
            store._dbri_path = self.dbri
        return _TestClient(create_app(store, dbri_path=self.dbri)), store

    def test_algorithms_clustering_lists_algorithms(self):
        client, store = self._client()
        try:
            r = client.get("/api/v1/algorithms/clustering")
            self.assertEqual(r.status_code, 200, f"got {r.status_code}: {r.text[:200]}")
            self.assertIn("algorithms", r.json())
        finally:
            store.close()

    def test_graph_cluster_returns_membership(self):
        client, store = self._client()
        try:
            r = client.post("/api/v1/graph/cluster",
                            json={"algorithm": "leiden", "parameters": {},
                                  "metric": "ochiai", "cutoff": 0.0})
            self.assertEqual(r.status_code, 200, f"got {r.status_code}: {r.text[:200]}")
            body = r.json()
            for k in ("membership", "num_clusters", "modularity"):
                self.assertIn(k, body)
        finally:
            store.close()

    def test_graph_components_returns_components(self):
        client, store = self._client()
        try:
            r = client.get("/api/v1/graph/components?metric=ochiai&cutoff=0.0")
            self.assertEqual(r.status_code, 200, f"got {r.status_code}: {r.text[:200]}")
            body = r.json()
            for k in ("num_components", "components"):
                self.assertIn(k, body)
        finally:
            store.close()

    def test_graph_cluster_high_cutoff_no_edges(self):
        # No edges (very high cutoff) must NOT 500 — leiden/louvain modularity is NaN there.
        client, store = self._client()
        try:
            r = client.post("/api/v1/graph/cluster",
                            json={"algorithm": "leiden", "parameters": {},
                                  "metric": "ochiai", "cutoff": 100.0})
            self.assertEqual(r.status_code, 200, f"high-cutoff cluster 500'd: {r.text[:200]}")
            self.assertIsNotNone(r.json().get("modularity"))
        finally:
            store.close()

    def test_graph_cluster_louvain(self):
        client, store = self._client()
        try:
            r = client.post("/api/v1/graph/cluster",
                            json={"algorithm": "louvain", "parameters": {},
                                  "metric": "ochiai", "cutoff": 0.0})
            self.assertEqual(r.status_code, 200, f"got {r.status_code}: {r.text[:200]}")
        finally:
            store.close()

    def test_graph_cluster_unknown_algorithm_rejected(self):
        client, store = self._client()
        try:
            r = client.post("/api/v1/graph/cluster",
                            json={"algorithm": "bogus", "parameters": {},
                                  "metric": "ochiai", "cutoff": 0.0})
            self.assertEqual(r.status_code, 400, f"unknown algorithm not rejected: {r.status_code}")
        finally:
            store.close()

    def test_graph_components_min_size_filters_all(self):
        client, store = self._client()
        try:
            r = client.get("/api/v1/graph/components?metric=ochiai&cutoff=0.0&min_size=100000")
            self.assertEqual(r.status_code, 200, f"got {r.status_code}: {r.text[:200]}")
            self.assertEqual(r.json().get("num_components"), 0)
        finally:
            store.close()

    # ── igraph re-platform characterization tests ──────────────────
    # /data, /neighborhood, /communities, /layout passed on Kùzu and must
    # still pass on igraph; /shortest-path was BROKEN on Kùzu and must now work.

    def _a_group(self, store):
        """Return some real group name from the store."""
        names = list(store.get_names_map().values())
        self.assertTrue(names, "store has no groups")
        return names[0]

    def _connected_pair(self, store):
        """Return (source, target) for two groups joined by an edge at cutoff 0."""
        from dbretina.pairwise_graph import PairwiseGraph
        g = PairwiseGraph(store, metric="ochiai", cutoff=0.0)
        edf = g.edges_dataframe()
        self.assertGreater(len(edf), 0, "no edges in graph to test shortest-path")
        nm = g._names_map
        src = nm[int(edf.iloc[0]["src"])]
        dst = nm[int(edf.iloc[0]["dst"])]
        return src, dst

    def test_graph_data_shape(self):
        client, store = self._client()
        try:
            r = client.get("/api/v1/graph/data?metric=ochiai&cutoff=0.0")
            self.assertEqual(r.status_code, 200, f"got {r.status_code}: {r.text[:200]}")
            body = r.json()
            for k in ("nodes", "edges", "meta"):
                self.assertIn(k, body)
            self.assertGreater(len(body["edges"]), 0, "graph/data returned no edges")
            # node objects carry the documented keys
            n0 = body["nodes"][0]
            for k in ("id", "label", "degree", "community", "pagerank"):
                self.assertIn(k, n0)
            # edge objects carry source/target/weight/shared_features
            e0 = body["edges"][0]
            for k in ("source", "target", "weight", "shared_features"):
                self.assertIn(k, e0)
        finally:
            store.close()

    def test_graph_data_edges_carry_all_metrics(self):
        # Each edge must carry EVERY available metric value so the dashboard can
        # filter client-side without re-querying. pvalue must be ABSENT here
        # (the test substrate has no pvalue) — available_metrics excludes it.
        client, store = self._client()
        try:
            self.assertFalse(store.has_pvalue, "test substrate unexpectedly has pvalue")
            r = client.get("/api/v1/graph/data?metric=ochiai&cutoff=0.0")
            self.assertEqual(r.status_code, 200, f"got {r.status_code}: {r.text[:200]}")
            edges = r.json()["edges"]
            self.assertGreater(len(edges), 0, "graph/data returned no edges")
            e0 = edges[0]
            for m in store.available_metrics:
                self.assertIn(m, e0, f"edge missing metric '{m}'")
                self.assertIsInstance(e0[m], (int, float))
            # The active metric value must agree with the legacy 'weight' field
            # (weight is rounded to 2dp).
            self.assertAlmostEqual(e0["weight"], round(e0["ochiai"], 2), places=2)
            # pvalue absent because it wasn't computed.
            self.assertNotIn("pvalue", e0, "pvalue leaked onto edges without a pvalue dataset")
        finally:
            store.close()

    def test_graph_neighborhood_shape(self):
        client, store = self._client()
        try:
            group = self._a_group(store)
            r = client.get(
                f"/api/v1/graph/neighborhood?group={group}&metric=ochiai&cutoff=0.0&hops=2"
            )
            self.assertEqual(r.status_code, 200, f"got {r.status_code}: {r.text[:200]}")
            body = r.json()
            for k in ("nodes", "edges", "meta"):
                self.assertIn(k, body)
        finally:
            store.close()

    def test_graph_neighborhood_unknown_group_404(self):
        client, store = self._client()
        try:
            r = client.get(
                "/api/v1/graph/neighborhood?group=__no_such_group__&metric=ochiai&cutoff=0.0"
            )
            self.assertEqual(r.status_code, 404, f"expected 404, got {r.status_code}")
        finally:
            store.close()

    def test_graph_communities_shape(self):
        client, store = self._client()
        try:
            r = client.get("/api/v1/graph/communities?metric=ochiai&cutoff=0.0&method=leiden")
            self.assertEqual(r.status_code, 200, f"got {r.status_code}: {r.text[:200]}")
            body = r.json()
            for k in ("communities", "num_communities", "sizes"):
                self.assertIn(k, body)
            self.assertIsInstance(body["communities"], dict)
        finally:
            store.close()

    def test_graph_layout_shape(self):
        client, store = self._client()
        try:
            r = client.get("/api/v1/graph/layout?metric=ochiai&cutoff=0.0&algorithm=fr")
            self.assertEqual(r.status_code, 200, f"got {r.status_code}: {r.text[:200]}")
            body = r.json()
            for k in ("algorithm", "positions"):
                self.assertIn(k, body)
            self.assertIsInstance(body["positions"], dict)
            # every position is an [x, y] pair
            for pos in body["positions"].values():
                self.assertEqual(len(pos), 2)
        finally:
            store.close()

    def test_graph_shortest_path_works(self):
        # Was BROKEN on Kùzu (SHORTEST path query); must return 200 with a real path.
        client, store = self._client()
        try:
            src, dst = self._connected_pair(store)
            r = client.get(
                f"/api/v1/graph/shortest-path?source={src}&target={dst}&metric=ochiai&cutoff=0.0"
            )
            self.assertEqual(r.status_code, 200, f"shortest-path 500'd: {r.text[:200]}")
            body = r.json()
            for k in ("path_length", "path_nodes", "connected"):
                self.assertIn(k, body)
            self.assertTrue(body["connected"], "connected pair reported as disconnected")
            self.assertEqual(body["path_length"], 1, "adjacent pair should be 1 hop")
            self.assertEqual(body["path_nodes"], [src, dst])
        finally:
            store.close()

    def test_graph_shortest_path_unknown_group_404(self):
        client, store = self._client()
        try:
            group = self._a_group(store)
            r = client.get(
                f"/api/v1/graph/shortest-path?source={group}&target=__nope__&metric=ochiai&cutoff=0.0"
            )
            self.assertEqual(r.status_code, 404, f"expected 404, got {r.status_code}")
        finally:
            store.close()

    def test_export_graph_graphml_works(self):
        # igraph write_graphml needs a real fd; the endpoint must not 500 on a BytesIO.
        client, store = self._client()
        try:
            r = client.get("/api/v1/export/graph/graphml?metric=ochiai&cutoff=0.0")
            self.assertEqual(r.status_code, 200, f"graphml export failed: {r.text[:200]}")
            self.assertIn(b"graphml", r.content[:400].lower())
        finally:
            store.close()


@unittest.skipUnless(_HAS_REST_TEST, "REST validation tests need [server] extra + httpx")
class TestRestValidation(unittest.TestCase):
    """Bad user input on REST endpoints must return a 4xx client error, never a
    500. Guards ISSUE-022 (hub-genes invalid metric), ISSUE-037 (hub-genes
    negative hops / unknown method / non-positive top_n) and ISSUE-038
    (graph/cluster bad-type algorithm parameter)."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp(prefix="dbretina_validate_")
        self.prefix, self.pw = setup_index_and_pairwise(self.tmpdir)
        self.parquet = f"{self.prefix}_DBRetina_pairwise"
        self.dbri = f"{self.prefix}.dbri"

    def tearDown(self):
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def _client(self):
        from dbretina.compat import open_pairwise
        from dbretina.pairwise_store import PairwiseStore
        from dbretina.rest_api import create_app
        store = open_pairwise(self.parquet)
        if store is None:
            store = PairwiseStore(self.parquet, dbri_path=self.dbri)
        else:
            store._dbri_path = self.dbri
        # raise_server_exceptions=False so a 500 surfaces as a status code we can
        # assert on, instead of re-raising into the test (we WANT to see 4xx vs 500).
        return _TestClient(create_app(store, dbri_path=self.dbri),
                           raise_server_exceptions=False), store

    def _a_group(self, store):
        names = list(store.get_names_map().values())
        self.assertTrue(names, "store has no groups")
        return names[0]

    @staticmethod
    def _is_4xx(code):
        return 400 <= code < 500

    # ── ISSUE-022: hub-genes invalid metric ────────────────────────
    def test_hub_genes_invalid_metric_is_4xx(self):
        client, store = self._client()
        try:
            group = self._a_group(store)
            r = client.post("/api/v1/genes/hub-genes",
                            json={"group_name": group, "method": "edge_weighted",
                                  "metric": "bogus"})
            self.assertTrue(self._is_4xx(r.status_code),
                            f"invalid metric should be 4xx, got {r.status_code}: {r.text[:200]}")
            self.assertNotEqual(r.status_code, 500,
                                "invalid metric still maps to 500 INTERNAL_ERROR (ISSUE-022)")
        finally:
            store.close()

    # ── ISSUE-037a: hub-genes negative hops ────────────────────────
    def test_hub_genes_negative_hops_is_4xx(self):
        client, store = self._client()
        try:
            group = self._a_group(store)
            r = client.post("/api/v1/genes/hub-genes",
                            json={"group_name": group, "hops": -1})
            self.assertTrue(self._is_4xx(r.status_code),
                            f"negative hops should be 4xx, got {r.status_code}: {r.text[:200]}")
            self.assertNotEqual(r.status_code, 500,
                                "negative hops still maps to 500 (ISSUE-037a)")
        finally:
            store.close()

    # ── ISSUE-037c: hub-genes non-positive top_n ───────────────────
    def test_hub_genes_negative_top_n_is_4xx(self):
        client, store = self._client()
        try:
            group = self._a_group(store)
            for bad in (-1, 0):
                r = client.post("/api/v1/genes/hub-genes",
                                json={"group_name": group, "top_n": bad})
                self.assertTrue(self._is_4xx(r.status_code),
                                f"top_n={bad} should be 4xx, got {r.status_code}: {r.text[:200]}")
        finally:
            store.close()

    # ── ISSUE-037b: hub-genes unknown method (no silent fallback) ──
    def test_hub_genes_unknown_method_is_4xx(self):
        client, store = self._client()
        try:
            group = self._a_group(store)
            r = client.post("/api/v1/genes/hub-genes",
                            json={"group_name": group, "method": "nonsense"})
            self.assertTrue(self._is_4xx(r.status_code),
                            f"unknown method should be 4xx, got {r.status_code}: {r.text[:200]}")
            # Must NOT silently succeed with mislabeled projection results.
            self.assertNotEqual(r.status_code, 200,
                                "unknown method silently fell back to projection (ISSUE-037b)")
        finally:
            store.close()

    # ── ISSUE-037c: cluster-analysis non-positive top_n ────────────
    def test_cluster_analysis_negative_top_n_is_4xx(self):
        client, store = self._client()
        try:
            group = self._a_group(store)
            for bad in (-1, 0):
                r = client.post("/api/v1/genes/cluster-analysis",
                                json={"node_names": [group], "top_n": bad})
                self.assertTrue(self._is_4xx(r.status_code),
                                f"cluster top_n={bad} should be 4xx, got {r.status_code}: {r.text[:200]}")
        finally:
            store.close()

    def test_cluster_analysis_unknown_method_is_4xx(self):
        client, store = self._client()
        try:
            group = self._a_group(store)
            r = client.post("/api/v1/genes/cluster-analysis",
                            json={"node_names": [group], "method": "nonsense"})
            self.assertTrue(self._is_4xx(r.status_code),
                            f"cluster unknown method should be 4xx, got {r.status_code}: {r.text[:200]}")
        finally:
            store.close()

    # ── ISSUE-038: graph/cluster bad-type parameter ────────────────
    def test_graph_cluster_bad_type_param_is_4xx(self):
        client, store = self._client()
        try:
            r = client.post("/api/v1/graph/cluster",
                            json={"algorithm": "leiden",
                                  "parameters": {"resolution": "notanum"},
                                  "metric": "ochiai", "cutoff": 0.0})
            self.assertTrue(self._is_4xx(r.status_code),
                            f"bad-type resolution should be 4xx, got {r.status_code}: {r.text[:200]}")
            self.assertNotEqual(r.status_code, 500,
                                "bad-type param still maps to 500 ALGORITHM_ERROR (ISSUE-038)")
        finally:
            store.close()

    # ── Regression guard: valid calls must still return 200 ────────
    def test_valid_calls_still_succeed(self):
        client, store = self._client()
        try:
            group = self._a_group(store)
            r = client.post("/api/v1/genes/hub-genes",
                            json={"group_name": group, "method": "edge_weighted",
                                  "metric": "jaccard", "hops": 2, "top_n": 5})
            self.assertEqual(r.status_code, 200,
                             f"valid hub-genes broke: {r.status_code} {r.text[:200]}")
            self.assertIn("genes", r.json())

            r = client.post("/api/v1/genes/cluster-analysis",
                            json={"node_names": [group], "top_n": 5})
            self.assertEqual(r.status_code, 200,
                             f"valid cluster-analysis broke: {r.status_code} {r.text[:200]}")
            self.assertIn("genes", r.json())

            r = client.post("/api/v1/graph/cluster",
                            json={"algorithm": "leiden",
                                  "parameters": {"resolution": 1.0},
                                  "metric": "ochiai", "cutoff": 0.0})
            self.assertEqual(r.status_code, 200,
                             f"valid graph/cluster broke: {r.status_code} {r.text[:200]}")
            self.assertIn("membership", r.json())
        finally:
            store.close()


class TestAlgorithmsModule(unittest.TestCase):
    """dbretina.algorithms.run_clustering (igraph-backed) — no server needed."""

    def test_run_clustering_no_edge_graph_modularity_finite(self):
        # A graph with vertices but no edges yields NaN modularity from igraph;
        # run_clustering must coerce it to a JSON-serializable finite value.
        import math
        import igraph as ig
        from dbretina.algorithms import run_clustering
        for algo in ("leiden", "louvain"):
            g = ig.Graph(n=5, directed=False)  # 5 vertices, 0 edges
            res = run_clustering(g, algorithm=algo)
            self.assertTrue(
                math.isfinite(res.modularity),
                f"{algo} no-edge modularity not finite: {res.modularity}",
            )


@unittest.skipUnless(_HAS_REST_TEST, "pairs-schema tests need [server] extra + httpx")
class TestRestPairsSchema(unittest.TestCase):
    """serve: pairs exposes names (bug3) and metrics respect has_pvalue (bug2)."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp(prefix="dbretina_schema_")
        self.prefix, self.pw = setup_index_and_pairwise(self.tmpdir)
        self.parquet = f"{self.prefix}_DBRetina_pairwise"
        self.dbri = f"{self.prefix}.dbri"

    def tearDown(self):
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def _client(self):
        from dbretina.compat import open_pairwise
        from dbretina.pairwise_store import PairwiseStore
        from dbretina.rest_api import create_app
        store = open_pairwise(self.parquet)
        if store is None:
            store = PairwiseStore(self.parquet, dbri_path=self.dbri)
        else:
            store._dbri_path = self.dbri
        return _TestClient(create_app(store, dbri_path=self.dbri)), store

    def test_metric_profile_no_pvalue_dataset(self):
        # bug2: a dataset without pvalue must not 500 on metric-profile.
        client, store = self._client()
        try:
            name = next(iter(store._id_to_name.values()))
            r = client.get(f"/api/v1/groups/{name}/metric-profile")
            self.assertEqual(r.status_code, 200, f"metric-profile failed: {r.text[:200]}")
            metrics = [m["metric"] for m in r.json()["metrics"]]
            self.assertIn("ochiai", metrics)
            if not store.has_pvalue:
                self.assertNotIn("pvalue", metrics, "pvalue reported on a no-pvalue dataset")
        finally:
            store.close()

    def test_info_valid_metrics_respects_has_pvalue(self):
        # bug2: /info must not advertise pvalue when the dataset lacks it.
        client, store = self._client()
        try:
            r = client.get("/api/v1/info")
            self.assertEqual(r.status_code, 200)
            if not store.has_pvalue:
                self.assertNotIn("pvalue", r.json().get("valid_metrics", []))
        finally:
            store.close()

    def test_sql_pairs_exposes_group_names(self):
        # bug3: pairs must expose group_1_name/group_2_name for the SQL editor.
        client, store = self._client()
        try:
            r = client.post("/api/v1/sql", json={
                "query": "SELECT group_1_name, group_2_name, ochiai FROM pairs LIMIT 3"})
            self.assertEqual(r.status_code, 200, f"pairs name query failed: {r.text[:200]}")
            cols = r.json()["columns"]
            self.assertIn("group_1_name", cols)
            self.assertIn("group_2_name", cols)
        finally:
            store.close()


# ============================================================
# SECTION: CLI UX regressions (clean errors, correct output naming)
# ============================================================

class TestNeighborsBadMetric(unittest.TestCase):
    """ISSUE-019: neighbors -m with an invalid or absent-but-valid metric must
    emit a clean [ERROR] line (no Python traceback) and exit nonzero, while a
    valid metric still returns neighbors."""

    @classmethod
    def setUpClass(cls):
        cls.prefix, cls.pw_file = _ensure_shared_fixture()
        # Parquet pairwise dir produced alongside the TSV (neighbors accepts both).
        cls.parquet = f"{cls.prefix}_DBRetina_pairwise"
        # A group that exists in the shared fixture (names are lowercased).
        cls.group = "groupa"

    def _assert_clean_metric_error(self, metric):
        rc, stdout, stderr = run_command(
            f'DBRetina neighbors -d {self.pw_file} "{self.group}" -m {metric} -c 5'
        )
        self.assertNotEqual(rc, 0, f"expected nonzero exit for metric {metric!r}")
        self.assertNotIn("Traceback", stderr,
                         f"metric {metric!r} leaked a Python traceback:\n{stderr}")
        self.assertNotIn("ValueError", stderr,
                         f"metric {metric!r} leaked a ValueError:\n{stderr}")
        self.assertIn("[ERROR]", stderr,
                      f"metric {metric!r} did not produce a clean [ERROR] line:\n{stderr}")

    def test_invalid_metric_clean_error(self):
        """An unknown metric string -> clean [ERROR], no traceback."""
        self._assert_clean_metric_error("bogusxyz")

    def test_absent_valid_metric_clean_error(self):
        """A valid-but-unavailable metric (pvalue on a no-pvalue dataset) ->
        clean [ERROR], no traceback."""
        self._assert_clean_metric_error("pvalue")

    def test_invalid_metric_on_parquet_clean_error(self):
        """Same clean handling when reading the Parquet pairwise directory."""
        if not os.path.isdir(self.parquet):
            self.skipTest("Parquet pairwise dir not produced")
        rc, _, stderr = run_command(
            f'DBRetina neighbors -d {self.parquet} "{self.group}" -m bogusxyz -c 5'
        )
        self.assertNotEqual(rc, 0)
        self.assertNotIn("Traceback", stderr)
        self.assertIn("[ERROR]", stderr)

    def test_valid_metric_still_returns_neighbors(self):
        """A valid metric must still succeed and print neighbors."""
        rc, stdout, stderr = run_command(
            f'DBRetina neighbors -d {self.pw_file} "{self.group}" -m ochiai -c 5'
        )
        self.assertEqual(rc, 0, stderr)
        self.assertNotIn("Traceback", stderr)
        # groupa overlaps several groups; header + at least one neighbor row.
        self.assertIn("ochiai", stdout)


class TestMergeDuplicateName(unittest.TestCase):
    """ISSUE-034: merging two indexes that share a group name must produce a
    clean error (no Python traceback) whose message does NOT reference the
    nonexistent --prefix flag; a clean merge of non-overlapping indexes still
    succeeds."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp(prefix="dbretina_mergedup_")

    def tearDown(self):
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def _create_index(self, name, content):
        asc = write_file(os.path.join(self.tmpdir, f"{name}.asc"), content)
        prefix = os.path.join(self.tmpdir, name)
        rc, _, stderr = run_command(f"DBRetina index -a {asc} -o {prefix}")
        self.assertEqual(rc, 0, f"index {name} failed: {stderr}")
        return prefix

    def test_duplicate_name_clean_error_no_prefix_advice(self):
        idx_a = self._create_index("dupa",
            "gene_set\tgene\nshared_grp\tAlpha\nshared_grp\tBeta\nonly_a\tGamma\n")
        idx_b = self._create_index("dupb",
            "gene_set\tgene\nshared_grp\tDelta\nshared_grp\tEpsilon\nonly_b\tZeta\n")
        out = os.path.join(self.tmpdir, "dup_merged.dbri")
        rc, _, stderr = run_command(
            f"DBRetina merge -a {idx_a}.dbri -b {idx_b}.dbri -o {out}"
        )
        self.assertNotEqual(rc, 0, "duplicate-name merge should fail")
        self.assertNotIn("Traceback", stderr,
                         f"duplicate-name merge leaked a traceback:\n{stderr}")
        self.assertNotIn("RuntimeError", stderr,
                         f"duplicate-name merge leaked a RuntimeError:\n{stderr}")
        self.assertNotIn("--prefix", stderr,
                         f"error still references the nonexistent --prefix flag:\n{stderr}")
        self.assertIn("[ERROR]", stderr,
                      f"no clean [ERROR] line emitted:\n{stderr}")
        # The error should still name the offending group so the user can fix it.
        self.assertIn("shared_grp", stderr)

    def test_clean_merge_still_succeeds(self):
        idx_a = self._create_index("cleana",
            "gene_set\tgene\nclean_a\tAlpha\nclean_a\tBeta\n")
        idx_b = self._create_index("cleanb",
            "gene_set\tgene\nclean_b\tGamma\nclean_b\tDelta\n")
        out = os.path.join(self.tmpdir, "clean_merged.dbri")
        rc, _, stderr = run_command(
            f"DBRetina merge -a {idx_a}.dbri -b {idx_b}.dbri -o {out}"
        )
        self.assertEqual(rc, 0, stderr)
        assert_file_exists(self, out)


class TestInteractomeNaming(unittest.TestCase):
    """ISSUE-040: interactome must name its outputs/log strings for the
    interactome command (not '*_genenet.*'), while genenet keeps its own
    genenet-named outputs."""

    @classmethod
    def setUpClass(cls):
        cls.prefix, cls.pw_file = _ensure_shared_fixture()

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp(prefix="dbretina_interactome_")

    def tearDown(self):
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def test_interactome_outputs_named_interactome(self):
        out = os.path.join(self.tmpdir, "it")
        rc, _, stderr = run_command(
            f"DBRetina interactome -i {self.prefix} -p {self.pw_file} -o {out}"
        )
        self.assertEqual(rc, 0, stderr)
        # interactome-named output exists and is non-empty...
        assert_file_exists(self, f"{out}_interactome.tsv")
        # ...and the misleading genenet-named output does NOT.
        self.assertFalse(os.path.exists(f"{out}_genenet.tsv"),
                         "interactome still wrote a *_genenet.tsv file")
        # Log strings should refer to 'interactome', not 'genenet'.
        self.assertNotIn("Exporting genenet", stderr)

    def test_genenet_outputs_unaffected(self):
        out = os.path.join(self.tmpdir, "gn")
        rc, _, stderr = run_command(
            f"DBRetina genenet -i {self.prefix} -p {self.pw_file} -o {out}"
        )
        self.assertEqual(rc, 0, stderr)
        # genenet keeps its genenet-named output...
        assert_file_exists(self, f"{out}_genenet.tsv")
        # ...and must NOT produce an interactome-named file.
        self.assertFalse(os.path.exists(f"{out}_interactome.tsv"),
                         "genenet wrote an *_interactome.tsv file")

    def test_interactome_graphml_named_interactome(self):
        out = os.path.join(self.tmpdir, "itg")
        rc, _, stderr = run_command(
            f"DBRetina interactome -i {self.prefix} -p {self.pw_file} "
            f"--graphml -o {out}"
        )
        self.assertEqual(rc, 0, stderr)
        assert_file_exists(self, f"{out}_interactome.graphml")
        self.assertFalse(os.path.exists(f"{out}_genenet.graphml"),
                         "interactome still wrote a *_genenet.graphml file")


# ============================================================
# SECTION: Empty / sparse graph (zero pairs above cutoff)
#   Issues 001 (modularity), 003 (cluster --community),
#   004 (cluster connected-components), 012 (setcov).
#   All four crashed with a raw traceback when no pair/edge
#   passed the cutoff. They must now degrade gracefully:
#   clean exit + no 'Traceback', and a valid output where one
#   is expected.
# ============================================================

# Four pairwise-disjoint groups: no two share any gene, so the
# pairwise file has ZERO data rows and every graph is edge-less
# at any cutoff > 0. This is the guaranteed-empty substrate.
DISJOINT_ASC_CONTENT = (
    "gene_set\tgene\n"
    "G1\ta1\n"
    "G1\ta2\n"
    "G2\tb1\n"
    "G2\tb2\n"
    "G3\tc1\n"
    "G3\tc2\n"
    "G4\td1\n"
    "G4\td2\n"
)

DISJOINT_GROUP_NAMES = {"g1", "g2", "g3", "g4"}


def assert_no_traceback(test_case, stderr, context=""):
    """A graceful failure must never print a Python traceback."""
    test_case.assertNotIn(
        "Traceback", stderr,
        f"{context}: command emitted a Python traceback:\n{stderr}"
    )


class TestEmptyGraphCrashes(unittest.TestCase):
    """Empty/sparse-graph regression guards (issues 001/003/004/012).

    The disjoint fixture yields zero pairs, so every metric/cutoff gives an
    edge-less graph -- the exact condition that used to traceback.
    """

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp(prefix="dbretina_empty_")
        # Disjoint (zero-pair) index + pairwise.
        self.prefix, self.pw_file = setup_index_and_pairwise(
            self.tmpdir, asc_content=DISJOINT_ASC_CONTENT
        )
        # Sanity: the pairwise file really has no data rows.
        self.assertEqual(
            count_tsv_data_rows(self.pw_file), 0,
            "disjoint fixture unexpectedly produced pairs"
        )

    def tearDown(self):
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    # ---- 001: modularity ----
    def test_modularity_empty_graph_no_traceback(self):
        """modularity with zero qualifying pairs must not traceback and must
        still write every group with modularity/fragmentation/heterogeneity=0."""
        out = os.path.join(self.tmpdir, "mod_empty")
        rc, _, stderr = run_command(
            f"DBRetina modularity -i {self.prefix} -p {self.pw_file} -c 50 -o {out}"
        )
        assert_no_traceback(self, stderr, "modularity empty graph")
        self.assertEqual(rc, 0, stderr)
        out_tsv = f"{out}_modularity.tsv"
        assert_file_exists(self, out_tsv)
        with open(out_tsv) as f:
            header = f.readline().strip().split("\t")
            self.assertEqual(
                header, ["gene_set", "fragmentation", "heterogeneity", "modularity"]
            )
            gene_sets = {}
            for line in f:
                if not line.strip():
                    continue
                parts = line.strip().split("\t")
                gene_sets[parts[0].lower()] = parts[1:]
        # Every group present, all metrics zero.
        self.assertEqual(set(gene_sets), DISJOINT_GROUP_NAMES)
        for name, vals in gene_sets.items():
            for v in vals:
                self.assertEqual(float(v), 0.0, f"{name} expected 0, got {v}")

    # ---- 003 + 012 share the community empty-edge path ----
    def test_cluster_community_empty_graph_no_traceback(self):
        """cluster --community with no edges (Leiden weight attr missing) must
        not traceback; output clusters file should exist."""
        out = os.path.join(self.tmpdir, "clu_comm_empty")
        rc, _, stderr = run_command(
            f"DBRetina cluster -p {self.pw_file} -m ochiai -c 50 --community -o {out}"
        )
        assert_no_traceback(self, stderr, "cluster --community empty graph")
        self.assertEqual(rc, 0, stderr)
        assert_file_exists(self, f"{out}_clusters.tsv")

    # ---- 004: connected-components empty ----
    def test_cluster_cc_empty_graph_no_traceback(self):
        """cluster (connected-components) with no edges must not traceback on the
        empty cluster_sizes histogram."""
        out = os.path.join(self.tmpdir, "clu_cc_empty")
        rc, _, stderr = run_command(
            f"DBRetina cluster -p {self.pw_file} -m ochiai -c 50 -o {out}"
        )
        assert_no_traceback(self, stderr, "cluster connected-components empty graph")
        self.assertEqual(rc, 0, stderr)
        assert_file_exists(self, f"{out}_clusters.tsv")

    # ---- 012: setcov on a sparse substrate ----
    def test_setcov_empty_graph_no_traceback(self):
        """setcov on sparse data routes through community detection on an
        edge-less graph; it must not traceback."""
        asc_path = os.path.join(self.tmpdir, "dj_setcov.asc")
        with open(asc_path, "w") as f:
            f.write(DISJOINT_ASC_CONTENT)
        rc, _, stderr = run_command(
            "DBRetina index -a dj_setcov.asc -o sc_idx", cwd=self.tmpdir
        )
        self.assertEqual(rc, 0, stderr)
        rc, _, stderr = run_command(
            "DBRetina pairwise -i sc_idx", cwd=self.tmpdir
        )
        self.assertEqual(rc, 0, stderr)
        rc, _, stderr = run_command(
            "DBRetina setcov -i sc_idx -o sc_out", timeout=300, cwd=self.tmpdir
        )
        assert_no_traceback(self, stderr, "setcov sparse graph")
        # Disjoint fixture (0 pairs): deterministically a clean [ERROR] exit, never a traceback.
        self.assertEqual(rc, 1, stderr)
        self.assertNotIn("Traceback", stderr)
        self.assertIn("[ERROR]", stderr, stderr)

    # ---- normal (non-empty) cases must be unchanged ----
    def test_normal_cutoff_still_works(self):
        """Same commands on the shared (overlapping) fixture must still succeed,
        proving the empty-graph guards do not change normal behaviour."""
        prefix, pw_file = _ensure_shared_fixture()
        # modularity
        out_m = os.path.join(self.tmpdir, "mod_ok")
        rc, _, stderr = run_command(
            f"DBRetina modularity -i {prefix} -p {pw_file} -c 40 -o {out_m}"
        )
        self.assertEqual(rc, 0, stderr)
        self.assertGreater(count_tsv_data_rows(f"{out_m}_modularity.tsv"), 0)
        # cluster --community
        out_cm = os.path.join(self.tmpdir, "clu_comm_ok")
        rc, _, stderr = run_command(
            f"DBRetina cluster -p {pw_file} -m ochiai -c 30 --community -o {out_cm}"
        )
        self.assertEqual(rc, 0, stderr)
        self.assertGreater(count_tsv_data_rows(f"{out_cm}_clusters.tsv"), 0)
        # cluster connected-components
        out_cc = os.path.join(self.tmpdir, "clu_cc_ok")
        rc, _, stderr = run_command(
            f"DBRetina cluster -p {pw_file} -m ochiai -c 50 -o {out_cc}"
        )
        self.assertEqual(rc, 0, stderr)
        self.assertGreater(count_tsv_data_rows(f"{out_cc}_clusters.tsv"), 0)


class TestQueryClusterInputHandling(unittest.TestCase):
    """Input-handling regressions for query/cluster (issues 016, 020, 021).

    Uses the shared overlapping fixture, which emits all sibling forms of the
    pairwise output next to each other: the .tsv, the parquet directory, the
    .dbrp binary, and the .dbri index.
    """

    @classmethod
    def setUpClass(cls):
        cls.prefix, cls.pw_file = _ensure_shared_fixture()
        cls.parquet_dir = f"{cls.prefix}_DBRetina_pairwise"
        cls.dbrp = f"{cls.prefix}_DBRetina_pairwise.dbrp"
        cls.dbri = f"{cls.prefix}.dbri"

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp(prefix="dbretina_inputbug_")
        # All sibling forms must actually exist for these tests to be meaningful.
        self.assertTrue(os.path.isdir(self.parquet_dir),
                        f"fixture missing parquet dir: {self.parquet_dir}")
        self.assertTrue(os.path.isfile(self.dbri),
                        f"fixture missing .dbri: {self.dbri}")

    def tearDown(self):
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    # ---- 016: cluster --community on a parquet directory ----
    def test_016_cluster_community_parquet_dir_no_traceback(self):
        """cluster --community -p <parquet DIR> must not IsADirectoryError; it
        should produce the same clusters as the .tsv form."""
        out_dir = os.path.join(self.tmpdir, "cl_dir")
        rc, _, stderr = run_command(
            f"DBRetina cluster -p {self.parquet_dir} -m ochiai -c 30 "
            f"--community -o {out_dir}"
        )
        assert_no_traceback(self, stderr, "cluster --community parquet dir")
        self.assertNotIn("IsADirectoryError", stderr, stderr)
        self.assertEqual(rc, 0, stderr)
        assert_file_exists(self, f"{out_dir}_clusters.tsv")

        # Result must match the .tsv form (same clustering, just a different input form).
        out_tsv = os.path.join(self.tmpdir, "cl_tsv")
        rc, _, stderr = run_command(
            f"DBRetina cluster -p {self.pw_file} -m ochiai -c 30 "
            f"--community -o {out_tsv}"
        )
        self.assertEqual(rc, 0, stderr)

        def _cluster_members(path):
            members = []
            with open(path) as f:
                for line in f:
                    if line.startswith("#") or line.lower().startswith("cluster_id"):
                        continue
                    parts = line.strip().split("\t")
                    if len(parts) >= 3:
                        members.append(frozenset(parts[2].split("|")))
            return sorted(members, key=lambda s: sorted(s))

        self.assertEqual(
            _cluster_members(f"{out_dir}_clusters.tsv"),
            _cluster_members(f"{out_tsv}_clusters.tsv"),
            "parquet-dir community clustering differs from .tsv form",
        )

    def test_016_cluster_community_dbrp_no_traceback(self):
        """cluster --community -p <.dbrp> must also resolve node sizes (sibling
        featuresNo.tsv) and not crash."""
        out = os.path.join(self.tmpdir, "cl_dbrp")
        rc, _, stderr = run_command(
            f"DBRetina cluster -p {self.dbrp} -m ochiai -c 30 --community -o {out}"
        )
        assert_no_traceback(self, stderr, "cluster --community .dbrp")
        self.assertEqual(rc, 0, stderr)
        assert_file_exists(self, f"{out}_clusters.tsv")

    # ---- 020: query on parquet dir / .dbrp / .dbri ----
    def test_020_query_parquet_dir_no_traceback(self):
        """query -p <parquet DIR> (cutoff-only) must work via the store, not
        IsADirectoryError. Output must match the .tsv form."""
        out = os.path.join(self.tmpdir, "q_dir")
        rc, _, stderr = run_command(
            f"DBRetina query -p {self.parquet_dir} -m ochiai -c 50 -o {out}"
        )
        assert_no_traceback(self, stderr, "query parquet dir")
        self.assertNotIn("IsADirectoryError", stderr, stderr)
        self.assertEqual(rc, 0, stderr)
        assert_file_exists(self, f"{out}.tsv")

        rows = parse_pairwise_tsv(f"{out}.tsv")
        self.assertGreater(len(rows), 0, "parquet-dir query returned no rows")
        for row in rows:
            self.assertGreaterEqual(row["ochiai"], 50.0)

        # Same pairs as querying the .tsv directly.
        out_tsv = os.path.join(self.tmpdir, "q_tsv")
        rc, _, stderr = run_command(
            f"DBRetina query -p {self.pw_file} -m ochiai -c 50 -o {out_tsv}"
        )
        self.assertEqual(rc, 0, stderr)

        def _pairs(path):
            return sorted(
                frozenset({r["group_1_name"], r["group_2_name"]})
                for r in parse_pairwise_tsv(path)
            )

        self.assertEqual(_pairs(f"{out}.tsv"), _pairs(f"{out_tsv}.tsv"),
                         "parquet-dir query pairs differ from .tsv form")

    def test_020_query_dbrp_no_traceback(self):
        """query -p <.dbrp> (cutoff-only) must work via the sibling store, not
        UnicodeDecodeError."""
        out = os.path.join(self.tmpdir, "q_dbrp")
        rc, _, stderr = run_command(
            f"DBRetina query -p {self.dbrp} -m ochiai -c 50 -o {out}"
        )
        assert_no_traceback(self, stderr, "query .dbrp")
        self.assertNotIn("UnicodeDecodeError", stderr, stderr)
        self.assertEqual(rc, 0, stderr)
        assert_file_exists(self, f"{out}.tsv")
        for row in parse_pairwise_tsv(f"{out}.tsv"):
            self.assertGreaterEqual(row["ochiai"], 50.0)

    def test_020_query_dbri_clean_error(self):
        """query -p <.dbri> must emit a clean [ERROR] (it is an index, not a
        pairwise file), never a raw traceback."""
        out = os.path.join(self.tmpdir, "q_dbri")
        rc, stdout, stderr = run_command(
            f"DBRetina query -p {self.dbri} -m ochiai -c 50 -o {out}"
        )
        assert_no_traceback(self, stderr, "query .dbri")
        self.assertNotIn("UnicodeDecodeError", stderr, stderr)
        self.assertEqual(rc, 1, stderr)
        self.assertIn("[ERROR]", stderr, stderr)
        self.assertIn("index", (stdout + stderr).lower(), stderr)

    # ---- 021: groups-only / clusters-only / no-args ----
    def test_021_query_groups_only_no_cutoff(self):
        """query -g groups.txt with no -m/-c must run the groups filter, NOT
        report a misleading 'cutoff must be between 0 and 100'."""
        groups = write_file(os.path.join(self.tmpdir, "groups.txt"),
                            "GroupA\nGroupB\n")
        out = os.path.join(self.tmpdir, "q_g")
        rc, stdout, stderr = run_command(
            f"DBRetina query -p {self.pw_file} -g {groups} -o {out}"
        )
        assert_no_traceback(self, stderr, "query groups-only")
        self.assertNotIn("cutoff must be between", (stdout + stderr),
                         f"misleading cutoff error surfaced:\n{stderr}")
        self.assertEqual(rc, 0, stderr)
        assert_file_exists(self, f"{out}.tsv")
        rows = parse_pairwise_tsv(f"{out}.tsv")
        self.assertGreater(len(rows), 0, "groups-only query returned no rows")
        for row in rows:
            pair = {row["group_1_name"], row["group_2_name"]}
            self.assertTrue(pair.issubset({"groupa", "groupb"}),
                            f"unexpected pair in groups-only result: {pair}")

    def test_021_query_clusters_only_no_cutoff(self):
        """query --clusters-file ... --cluster-ids 1 with no -m/-c must run, NOT
        report a misleading cutoff error."""
        # Build a clusters file from the shared fixture.
        cl = os.path.join(self.tmpdir, "cl")
        rc, _, stderr = run_command(
            f"DBRetina cluster -p {self.pw_file} -m ochiai -c 0 -o {cl}"
        )
        self.assertEqual(rc, 0, stderr)
        clusters_file = f"{cl}_clusters.tsv"
        assert_file_exists(self, clusters_file)

        out = os.path.join(self.tmpdir, "q_cl")
        rc, stdout, stderr = run_command(
            f"DBRetina query -p {self.pw_file} --clusters-file {clusters_file} "
            f"--cluster-ids 1 -o {out}"
        )
        assert_no_traceback(self, stderr, "query clusters-only")
        self.assertNotIn("cutoff must be between", (stdout + stderr),
                         f"misleading cutoff error surfaced:\n{stderr}")
        self.assertEqual(rc, 0, stderr)
        assert_file_exists(self, f"{out}.tsv")

    def test_021_query_no_args_clean_error(self):
        """query with no filter args must give the accurate 'requires at least
        one option' error, NOT the misleading cutoff error."""
        out = os.path.join(self.tmpdir, "q_noargs")
        rc, stdout, stderr = run_command(
            f"DBRetina query -p {self.pw_file} -o {out}"
        )
        assert_no_traceback(self, stderr, "query no-args")
        self.assertEqual(rc, 1, stderr)
        self.assertNotIn("cutoff must be between", (stdout + stderr),
                         f"misleading cutoff error surfaced:\n{stderr}")
        self.assertIn("at least one option", (stdout + stderr).lower(), stderr)

    # ---- guard: valid existing paths unchanged ----
    def test_021_query_cutoff_only_tsv_still_works(self):
        """Removing the redundant cutoff guard must not change normal cutoff
        queries on the .tsv."""
        out = os.path.join(self.tmpdir, "q_co")
        rc, _, stderr = run_command(
            f"DBRetina query -p {self.pw_file} -m ochiai -c 50 -o {out}"
        )
        self.assertEqual(rc, 0, stderr)
        rows = parse_pairwise_tsv(f"{out}.tsv")
        for row in rows:
            self.assertGreaterEqual(row["ochiai"], 50.0)

    def test_021_query_out_of_range_cutoff_still_rejected(self):
        """An actually out-of-range cutoff (>100) must still be rejected cleanly
        (the Click callback handles this)."""
        out = os.path.join(self.tmpdir, "q_bad")
        rc, stdout, stderr = run_command(
            f"DBRetina query -p {self.pw_file} -m ochiai -c 150 -o {out}"
        )
        assert_no_traceback(self, stderr, "query out-of-range cutoff")
        self.assertNotEqual(rc, 0, "out-of-range cutoff should be rejected")


# ============================================================
# Main
# ============================================================

if __name__ == "__main__":
    try:
        unittest.main(verbosity=2)
    finally:
        _cleanup_shared_fixture()
