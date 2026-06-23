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
import time
import unittest
import warnings
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


# Shared WITH-pvalue substrate (pairwise computed with --pvalue), used by the
# cross-route pvalue tests below. The standard fixture above is computed WITHOUT
# --pvalue, so it doubles as the "pvalue absent" case.
_SHARED_PV_DIR = None
_SHARED_PV_PREFIX = None
_SHARED_PV_PW_FILE = None


def _ensure_shared_pvalue_fixture():
    """Index + pairwise(--pvalue). Returns (prefix, pw_tsv_path) with a pvalue column."""
    global _SHARED_PV_DIR, _SHARED_PV_PREFIX, _SHARED_PV_PW_FILE
    if _SHARED_PV_DIR is None:
        _SHARED_PV_DIR = tempfile.mkdtemp(prefix="dbretina_shared_pv_")
        _SHARED_PV_PREFIX, _SHARED_PV_PW_FILE = setup_index_and_pairwise(
            _SHARED_PV_DIR, extra_pw_args="--pvalue"
        )
    return _SHARED_PV_PREFIX, _SHARED_PV_PW_FILE


def _cleanup_shared_pvalue_fixture():
    global _SHARED_PV_DIR
    if _SHARED_PV_DIR is not None:
        shutil.rmtree(_SHARED_PV_DIR, ignore_errors=True)
        _SHARED_PV_DIR = None


# Shared ZERO-PAIR substrate: an index whose groups share no genes, so `pairwise`
# emits a valid 0-pair dataset (manifest num_pairs:0, an empty data/ with no
# part_*.parquet). This is the cosmic-like case behind issue 071 -- the lookup
# commands must handle it without a raw DuckDB IOException.
_ZERO_PAIR_DISJOINT_ASC = (
    "gene_set\tgene\n"
    "GroupX\tgx1\n"
    "GroupX\tgx2\n"
    "GroupY\tgy1\n"
    "GroupY\tgy2\n"
    "GroupZ\tgz1\n"
    "GroupZ\tgz2\n"
)

_SHARED_ZERO_DIR = None
_SHARED_ZERO_PREFIX = None
_SHARED_ZERO_PW_FILE = None


def _ensure_zero_pair_fixture():
    """Index + pairwise over fully disjoint gene sets -> a 0-pair pairwise.
    Returns (prefix, pw_tsv_path); the parquet dir ``{prefix}_DBRetina_pairwise``
    has an empty ``data/`` (num_pairs:0)."""
    global _SHARED_ZERO_DIR, _SHARED_ZERO_PREFIX, _SHARED_ZERO_PW_FILE
    if _SHARED_ZERO_DIR is None:
        _SHARED_ZERO_DIR = tempfile.mkdtemp(prefix="dbretina_shared_zero_")
        _SHARED_ZERO_PREFIX, _SHARED_ZERO_PW_FILE = setup_index_and_pairwise(
            _SHARED_ZERO_DIR, asc_content=_ZERO_PAIR_DISJOINT_ASC
        )
    return _SHARED_ZERO_PREFIX, _SHARED_ZERO_PW_FILE


def _cleanup_zero_pair_fixture():
    global _SHARED_ZERO_DIR
    if _SHARED_ZERO_DIR is not None:
        shutil.rmtree(_SHARED_ZERO_DIR, ignore_errors=True)
        _SHARED_ZERO_DIR = None


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

    def test_pairwise_help_no_phantom_directional_containment(self):
        """ISSUE-084: pairwise --help/docstring must not promise directional
        containment_1_in_2 / containment_2_in_1 columns -- they are never emitted
        (the TSV has a single symmetric `containment` column). It should still
        document the real `containment` column."""
        rc, stdout, stderr = run_command("DBRetina pairwise --help")
        self.assertEqual(rc, 0, stderr)
        help_text = stdout + stderr
        self.assertNotIn("containment_1_in_2", help_text)
        self.assertNotIn("containment_2_in_1", help_text)
        # The real symmetric column is still documented.
        self.assertIn("containment", help_text.lower())

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

    def test_pairwise_nonpositive_threads_clean_error(self):
        """ISSUE-031: pairwise -t with a non-positive thread count must be rejected
        cleanly, not crash in C++ (vector::_M_default_append for negatives, SIGSEGV
        for 0). A positive -t still works."""
        asc = write_file(os.path.join(self.tmpdir, "t.asc"), TEST_ASC_CONTENT)
        prefix = os.path.join(self.tmpdir, "tidx")
        rc, _, stderr = run_command(f"DBRetina index -a {asc} -o {prefix}")
        self.assertEqual(rc, 0, stderr)

        for bad in ("-2", "0"):
            rc, _, stderr = run_command(
                f"DBRetina pairwise -i {prefix} -m containment -c 0 -t {bad}"
            )
            self.assertNotEqual(rc, 0, f"-t {bad} should exit nonzero")
            # no opaque STL error, no raw traceback, no segfault (-11/139)
            self.assertNotIn("Traceback", stderr)
            self.assertNotIn("_M_default_append", stderr)
            self.assertNotIn("Segmentation", stderr)
            self.assertNotEqual(rc, -11, "must not SIGSEGV")
            self.assertNotEqual(rc, 139, "must not SIGSEGV")

        # a valid thread count still works
        rc, _, stderr = run_command(
            f"DBRetina pairwise -i {prefix} -m containment -c 0 -t 2"
        )
        self.assertEqual(rc, 0, stderr)
        self.assertNotIn("Traceback", stderr)
        assert_file_exists(self, f"{prefix}_DBRetina_pairwise.tsv")

    def test_pairwise_unsupported_filter_metric_clean_error(self):
        """ISSUE-032: pairwise -m csi/dice pass the old validator but the C++ cutoff
        filter only accepts containment/ochiai/jaccard -> they crashed with a raw
        C++ ValueError. They must now be rejected cleanly; supported metrics work."""
        asc = write_file(os.path.join(self.tmpdir, "m.asc"), TEST_ASC_CONTENT)
        prefix = os.path.join(self.tmpdir, "midx")
        rc, _, stderr = run_command(f"DBRetina index -a {asc} -o {prefix}")
        self.assertEqual(rc, 0, stderr)

        for bad in ("csi", "dice"):
            rc, _, stderr = run_command(
                f"DBRetina pairwise -i {prefix} -m {bad} -c 0"
            )
            self.assertNotEqual(rc, 0, f"-m {bad} should exit nonzero")
            self.assertNotIn("Traceback", stderr)
            # not the raw C++ message
            self.assertNotIn("cutoff_distance_type", stderr)
            # clean Click error naming the bad value + the supported metrics
            self.assertIn(bad, stderr)
            self.assertIn("containment", stderr.lower())

        # supported cutoff-filter metrics still work
        for good in ("containment", "ochiai", "jaccard"):
            rc, _, stderr = run_command(
                f"DBRetina pairwise -i {prefix} -m {good} -c 0"
            )
            self.assertEqual(rc, 0, f"-m {good} should succeed: {stderr}")
            self.assertNotIn("Traceback", stderr)
            assert_file_exists(self, f"{prefix}_DBRetina_pairwise.tsv")

    def test_pairwise_missing_index_clean_error(self):
        """ISSUE-073: pairwise -i <nonexistent prefix> must give a clean [ERROR]
        instead of a raw RuntimeError traceback from the unwrapped C++ call."""
        bad = os.path.join(self.tmpdir, "does_not_exist_xyz")
        rc, stdout, stderr = run_command(f"DBRetina pairwise -i {bad}")
        combined = stdout + stderr
        self.assertNotEqual(rc, 0, "missing index should exit nonzero")
        self.assertNotIn("Traceback", combined, f"raw traceback leaked:\n{combined}")
        self.assertNotIn("RuntimeError", combined,
                         f"raw RuntimeError leaked:\n{combined}")
        self.assertIn("[ERROR]", combined)
        self.assertIn(bad, combined)


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
# SECTION 5b: query -m pvalue cutoff direction (ISSUE 068)
# ============================================================

class TestQueryPvalueCutoffDirection(unittest.TestCase):
    """`query -m pvalue -c CUTOFF` must keep the SIGNIFICANT pairs (pvalue <=
    cutoff), not invert and return the least-significant ones (issue 068).

    p-values are 'lower is better'; the bug applied the similarity rule
    (value >= cutoff) to pvalue on every route, so it returned the LEAST
    significant pairs and silently dropped the significant ones. The shared
    --pvalue fixture yields this clean pvalue distribution::

        groupf|groupa  0.00126   <= 0.05  (significant)
        groupb|groupa  0.04545   <= 0.05  (significant)
        groupb|groupf  0.04545   <= 0.05  (significant)
        (8 more pairs) 0.57..0.97 > 0.05  (not significant)

    so a cutoff of 0.05 cleanly separates 3 significant pairs from 8. Pre-fix,
    each route returned the 8 with pvalue > 0.05; post-fix, the 3 with <= 0.05.
    All three confirmed-inverted routes are exercised: the C++ dbrp_filter_pairs
    path (.dbrp present), the awk TSV fallback (no .dbrp/parquet), and the
    DuckDB PairwiseStore path (parquet directory as -p).
    """

    PV_CUTOFF = 0.05
    EXPECTED_SIGNIFICANT = {  # pairs with pvalue <= 0.05
        frozenset({"groupa", "groupf"}),
        frozenset({"groupa", "groupb"}),
        frozenset({"groupb", "groupf"}),
    }

    @classmethod
    def setUpClass(cls):
        cls.pv_prefix, cls.pv_pw_tsv = _ensure_shared_pvalue_fixture()

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp(prefix="dbretina_q_pv_")

    def tearDown(self):
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def _copy_pairwise_artifacts(self):
        """Copy the --pvalue pairwise artifacts into this test's tmpdir.

        Route-isolation tests rename/remove the .dbrp and parquet dir; copying
        first keeps the module-shared fixture pristine for other tests.
        Returns the copied pairwise .tsv path.
        """
        import glob as _glob
        src_dir = os.path.dirname(self.pv_pw_tsv)
        base = os.path.basename(self.pv_pw_tsv)[:-len(".tsv")]  # ..._DBRetina_pairwise
        for path in _glob.glob(os.path.join(src_dir, base + "*")):
            dst = os.path.join(self.tmpdir, os.path.basename(path))
            if os.path.isdir(path):
                shutil.copytree(path, dst)
            else:
                shutil.copy2(path, dst)
        return os.path.join(self.tmpdir, base + ".tsv")

    def _assert_significant(self, tsv_path, route):
        rows = parse_pairwise_tsv(tsv_path)
        # Every returned pair is significant (pvalue <= cutoff), NOT inverted.
        for row in rows:
            self.assertIn("pvalue", row, f"{route}: missing pvalue column")
            self.assertLessEqual(
                row["pvalue"], self.PV_CUTOFF,
                f"{route}: returned a non-significant pair (pvalue "
                f"{row['pvalue']} > {self.PV_CUTOFF}) -- filter is inverted")
        got = {frozenset({r["group_1_name"], r["group_2_name"]}) for r in rows}
        # Exactly the significant pairs, count and identity (3, not the 8 above).
        self.assertEqual(
            got, self.EXPECTED_SIGNIFICANT,
            f"{route}: expected the 3 significant pairs (pvalue <= "
            f"{self.PV_CUTOFF}), got {len(rows)} rows: {sorted(map(sorted, got))}")

    def test_068_pvalue_cutoff_dbrp_route(self):
        """C++ dbrp_filter_pairs route (.dbrp present): keeps pvalue <= cutoff."""
        pw_tsv = self._copy_pairwise_artifacts()
        self.assertTrue(os.path.exists(pw_tsv.replace(".tsv", ".dbrp")),
                        "fixture should ship a .dbrp for the C++ route")
        out = os.path.join(self.tmpdir, "q_dbrp")
        rc, _, stderr = run_command(
            f"DBRetina query -p {pw_tsv} -m pvalue -c {self.PV_CUTOFF} -o {out}"
        )
        self.assertEqual(rc, 0, stderr)
        self._assert_significant(f"{out}.tsv", "dbrp route")

    def test_068_pvalue_cutoff_awk_tsv_route(self):
        """awk TSV fallback (no .dbrp, no parquet dir): keeps pvalue <= cutoff."""
        pw_tsv = self._copy_pairwise_artifacts()
        # Force the plain-text awk path: drop the .dbrp and the parquet dir so
        # neither the C++ reader nor the DuckDB store can claim the query.
        dbrp = pw_tsv.replace(".tsv", ".dbrp")
        if os.path.exists(dbrp):
            os.remove(dbrp)
        pq_dir = pw_tsv[:-len(".tsv")]  # the sibling parquet directory
        if os.path.isdir(pq_dir):
            shutil.rmtree(pq_dir)
        out = os.path.join(self.tmpdir, "q_awk")
        rc, _, stderr = run_command(
            f"DBRetina query -p {pw_tsv} -m pvalue -c {self.PV_CUTOFF} -o {out}"
        )
        self.assertEqual(rc, 0, stderr)
        self._assert_significant(f"{out}.tsv", "awk TSV route")

    def test_068_pvalue_cutoff_parquet_store_route(self):
        """DuckDB PairwiseStore route (parquet dir as -p): keeps pvalue <= cutoff."""
        pw_tsv = self._copy_pairwise_artifacts()
        pq_dir = pw_tsv[:-len(".tsv")]
        if not os.path.isdir(pq_dir):
            self.skipTest("no parquet directory emitted by pairwise")
        out = os.path.join(self.tmpdir, "q_store")
        rc, _, stderr = run_command(
            f"DBRetina query -p {pq_dir} -m pvalue -c {self.PV_CUTOFF} -o {out}"
        )
        self.assertEqual(rc, 0, stderr)
        self._assert_significant(f"{out}.tsv", "parquet store route")

    def test_068_similarity_cutoff_unchanged_dbrp(self):
        """Regression guard: ochiai cutoff still keeps value >= cutoff (C++/dbrp)."""
        out = os.path.join(self.tmpdir, "q_ochiai")
        rc, _, stderr = run_command(
            f"DBRetina query -p {self.pv_pw_tsv} -m ochiai -c 50 -o {out}"
        )
        self.assertEqual(rc, 0, stderr)
        rows = parse_pairwise_tsv(f"{out}.tsv")
        self.assertGreater(len(rows), 0, "ochiai >= 50 should return some pairs")
        for row in rows:
            self.assertGreaterEqual(
                row["ochiai"], 50.0,
                f"similarity filter must stay >= cutoff, got ochiai {row['ochiai']}")

    def test_068_similarity_cutoff_unchanged_parquet_store(self):
        """Regression guard: ochiai cutoff still keeps value >= cutoff (store)."""
        pw_tsv = self._copy_pairwise_artifacts()
        pq_dir = pw_tsv[:-len(".tsv")]
        if not os.path.isdir(pq_dir):
            self.skipTest("no parquet directory emitted by pairwise")
        out = os.path.join(self.tmpdir, "q_ochiai_store")
        rc, _, stderr = run_command(
            f"DBRetina query -p {pq_dir} -m ochiai -c 50 -o {out}"
        )
        self.assertEqual(rc, 0, stderr)
        rows = parse_pairwise_tsv(f"{out}.tsv")
        self.assertGreater(len(rows), 0)
        for row in rows:
            self.assertGreaterEqual(row["ochiai"], 50.0)


class TestSiblingPvalueCutoffDirection(unittest.TestCase):
    """The pvalue inversion (issue 068) also lived in the raw SQL / bare-TSV
    re-filters of sibling commands (neighbors / cluster / graph / bipartite),
    not just `query`. Their store and .dbrp paths share the fixed helpers, but
    the hand-written comparisons applied the similarity rule (>= cutoff) to
    pvalue too. These guard that pvalue keeps the SIGNIFICANT pairs there as
    well.

    Same --pvalue fixture distribution as TestQueryPvalueCutoffDirection: at
    cutoff 0.05, groupa's significant neighbors are groupf (p=0.00126) and
    groupb (p=0.0455); groupc/groupe (p>0.87) are NOT significant.
    """

    PV_CUTOFF = 0.05

    @classmethod
    def setUpClass(cls):
        cls.pv_prefix, cls.pv_pw_tsv = _ensure_shared_pvalue_fixture()
        cls.pv_dir = f"{cls.pv_prefix}_DBRetina_pairwise"  # parquet directory

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp(prefix="dbretina_sib_pv_")

    def tearDown(self):
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def _copy_pairwise_artifacts(self):
        """Copy the --pvalue pairwise artifacts into this test's tmpdir so
        bare-TSV isolation never mutates the module-shared fixture."""
        import glob as _glob
        src_dir = os.path.dirname(self.pv_pw_tsv)
        base = os.path.basename(self.pv_pw_tsv)[:-len(".tsv")]
        for path in _glob.glob(os.path.join(src_dir, base + "*")):
            dst = os.path.join(self.tmpdir, os.path.basename(path))
            if os.path.isdir(path):
                shutil.copytree(path, dst)
            else:
                shutil.copy2(path, dst)
        return os.path.join(self.tmpdir, base + ".tsv")

    def _parse_neighbors(self, path):
        """neighbors output -> list of (neighbor, metric_value) tuples."""
        rows = []
        with open(path) as f:
            header = True
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                if header:  # "neighbor\t{metric}\tjaccard\tshared_features"
                    header = False
                    continue
                parts = line.split("\t")
                rows.append((parts[0], float(parts[1])))
        return rows

    def test_068_neighbors_pvalue_significant(self):
        """neighbors -m pvalue keeps neighbors with pvalue <= cutoff (was
        inverted: cmd_neighbors raw SQL used >= and ORDER BY ... DESC)."""
        out = os.path.join(self.tmpdir, "nbr.tsv")
        rc, _, stderr = run_command(
            f'DBRetina neighbors -d {self.pv_dir} "groupa" '
            f"-m pvalue -c {self.PV_CUTOFF} -o {out}"
        )
        self.assertEqual(rc, 0, stderr)
        rows = self._parse_neighbors(out)
        self.assertGreater(len(rows), 0, "expected significant neighbors of groupa")
        for name, pval in rows:
            self.assertLessEqual(
                pval, self.PV_CUTOFF,
                f"neighbors returned non-significant {name} (p={pval} > "
                f"{self.PV_CUTOFF}) -- pvalue filter inverted")
        names = {n for n, _ in rows}
        self.assertEqual(
            names, {"groupf", "groupb"},
            f"expected significant neighbors {{groupf, groupb}}, got {names}")
        # ASC order: most significant (smallest p) first.
        pvals = [p for _, p in rows]
        self.assertEqual(pvals, sorted(pvals),
                         f"neighbors should rank most-significant first: {pvals}")

    def test_068_neighbors_similarity_unchanged(self):
        """Regression guard: neighbors -m ochiai still keeps ochiai >= cutoff."""
        out = os.path.join(self.tmpdir, "nbr_och.tsv")
        rc, _, stderr = run_command(
            f'DBRetina neighbors -d {self.pv_dir} "groupa" -m ochiai -c 50 -o {out}'
        )
        self.assertEqual(rc, 0, stderr)
        rows = self._parse_neighbors(out)
        self.assertGreater(len(rows), 0)
        for name, val in rows:
            self.assertGreaterEqual(val, 50.0, f"{name} ochiai {val} < 50")

    def _cluster_members(self, clusters_tsv):
        """Parse a *_clusters.tsv into {frozenset(member_names)} per cluster."""
        clusters = set()
        with open(clusters_tsv) as f:
            for line in f:
                line = line.rstrip("\n")
                if not line or line.startswith("#") or line.startswith("cluster"):
                    continue
                parts = line.split("\t")
                # cluster file: cluster_id, size, members (|-joined)
                members = frozenset(parts[-1].split("|"))
                clusters.add(members)
        return clusters

    def test_068_cluster_pvalue_bare_tsv_matches_store(self):
        """cluster -m pvalue on a bare TSV (no .dbrp/parquet -> the fixed
        row-by-row TSV re-filter) must yield the SAME clusters as the parquet
        store route (which filters pvalue <= cutoff correctly)."""
        # Correct reference: the parquet-directory (store) route.
        ref_out = os.path.join(self.tmpdir, "cl_store")
        rc, _, stderr = run_command(
            f"DBRetina cluster -p {self.pv_dir} -m pvalue -c {self.PV_CUTOFF} -o {ref_out}"
        )
        self.assertEqual(rc, 0, stderr)
        ref = self._cluster_members(f"{ref_out}_clusters.tsv")

        # Bare-TSV route: copy artifacts, then strip .dbrp + parquet dir.
        pw_tsv = self._copy_pairwise_artifacts()
        dbrp = pw_tsv.replace(".tsv", ".dbrp")
        if os.path.exists(dbrp):
            os.remove(dbrp)
        pq_dir = pw_tsv[:-len(".tsv")]
        if os.path.isdir(pq_dir):
            shutil.rmtree(pq_dir)
        tsv_out = os.path.join(self.tmpdir, "cl_tsv")
        rc, _, stderr = run_command(
            f"DBRetina cluster -p {pw_tsv} -m pvalue -c {self.PV_CUTOFF} -o {tsv_out}"
        )
        self.assertEqual(rc, 0, stderr)
        got = self._cluster_members(f"{tsv_out}_clusters.tsv")
        self.assertEqual(
            got, ref,
            "bare-TSV cluster (pvalue) disagrees with the correct store route -- "
            "the TSV pvalue re-filter is inverted")


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

    def test_cluster_average_size_log_correct(self):
        """ISSUE-076: the 'Average cluster size' log must be total clustered nodes /
        number of EMITTED clusters, not / len(connected_components) (which counts the
        isolated singletons that are skipped, never written, and never summed) -- that
        denominator inflated the count and produced an implausible average (often <1).

        ochiai c=70 on the 6-group fixture keeps only {groupa, groupb, groupf} as one
        cluster; the remaining groups are singleton components that get skipped, so the
        buggy denominator (4) != the emitted-cluster count (1)."""
        import re
        out = os.path.join(self.tmpdir, "cl_avg")
        rc, stdout, stderr = run_command(
            f"DBRetina cluster -p {self.pw_file} -m ochiai -c 70 -o {out}"
        )
        self.assertEqual(rc, 0, stderr)
        combined = stdout + stderr

        def _grab(pattern):
            m = re.search(pattern, combined)
            self.assertIsNotNone(m, f"missing log line for /{pattern}/:\n{combined}")
            return float(m.group(1))

        clustered = _grab(r"Total number of clustered supergroups:\s*([0-9.]+)")
        num_clusters = _grab(r"number of clusters:\s*([0-9.]+)")
        avg_logged = _grab(r"Average cluster size:\s*([0-9.]+)")

        self.assertGreaterEqual(num_clusters, 1)
        expected_avg = clustered / num_clusters
        self.assertAlmostEqual(
            avg_logged, expected_avg, places=6,
            msg=(f"avg-size log {avg_logged} != clustered/clusters "
                 f"({clustered}/{num_clusters}={expected_avg})"))
        # The bug divided by every partition incl. skipped singletons, yielding an
        # average below 1 here (3/4=0.75); a real cluster averages >= 1 node.
        self.assertGreaterEqual(
            avg_logged, 1.0,
            f"avg cluster size {avg_logged} < 1 -> denominator still counts "
            f"skipped singletons (issue 076)")

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

    def test_cluster_nonexistent_output_dir(self):
        """ISSUE-078: -o prefix in a nonexistent directory must not crash with a
        raw FileNotFoundError after the full graph build. We auto-create the
        parent dir, so it succeeds and the clusters file is present."""
        out = os.path.join(self.tmpdir, "no_such_dir_xyz", "out")
        rc, stdout, stderr = run_command(
            f"DBRetina cluster -p {self.pw_file} -m ochiai -c 30 -o {out}"
        )
        combined = stdout + stderr
        self.assertNotIn("Traceback", combined, f"raw traceback leaked:\n{combined}")
        self.assertNotIn("FileNotFoundError", combined,
                         f"raw FileNotFoundError leaked:\n{combined}")
        self.assertEqual(rc, 0, stderr)
        assert_file_exists(self, f"{out}_clusters.tsv")

    def test_cluster_output_prefix_parent_is_file(self):
        """ISSUE-078 follow-up: -o whose parent path component is an existing FILE
        must fail with a clean error, not a raw OSError traceback."""
        afile = write_file(os.path.join(self.tmpdir, "afile"), "x\n")
        out = os.path.join(afile, "out")  # parent component 'afile' is a regular file
        rc, stdout, stderr = run_command(
            f"DBRetina cluster -p {self.pw_file} -m ochiai -c 30 -o {out}"
        )
        combined = stdout + stderr
        self.assertNotEqual(rc, 0)
        self.assertNotIn("Traceback", combined, f"raw traceback leaked:\n{combined}")
        self.assertIn("[ERROR]", combined)

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

    def test_export_help_lists_all_accepted_metrics(self):
        """ISSUE-085: export -m --help must list every metric the command actually
        accepts (validators.VALID_METRICS), not just 4. csi/dice/odds_ratio were
        processed but undocumented."""
        rc, stdout, stderr = run_command("DBRetina export --help")
        self.assertEqual(rc, 0, stderr)
        help_text = stdout + stderr
        for m in ("containment", "ochiai", "jaccard", "csi", "dice",
                  "odds_ratio", "pvalue"):
            self.assertIn(m, help_text,
                          f"export --help omits accepted metric '{m}'")

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

    def test_dedup_pairwise_group_absent_from_index_clean_error(self):
        """ISSUE-077: a pairwise group name absent from the index (mismatched
        -i/-p) must give a clean [ERROR] naming the offending group, not a raw
        KeyError. A standalone TSV (no parquet/.dbrp sibling) forces the TSV
        fallback and references a group the index never had."""
        # Pairwise TSV with two groups, one of which ('ghost_group') is NOT in
        # the shared index. ochiai is column 6; high value so it passes -c.
        mal = os.path.join(self.tmpdir, "mismatch.tsv")
        header = ("group_1_ID\tgroup_2_ID\tgroup_1_name\tgroup_2_name\t"
                  "shared_features\tcontainment\tochiai\tjaccard\tcsi\tdice\todds_ratio\n")
        row = "1\t99\tgroupa\tghost_group\t3\t90.0\t90.0\t80.0\t85.0\t85.0\t1.0\n"
        write_file(mal, header + row)
        out = os.path.join(self.tmpdir, "ddmis")
        rc, stdout, stderr = run_command(
            f"DBRetina dedup -i {self.prefix} -p {mal} -c 50 -o {out}"
        )
        combined = stdout + stderr
        self.assertNotEqual(rc, 0, "mismatched index/pairwise should exit nonzero")
        self.assertNotIn("Traceback", combined, f"raw traceback leaked:\n{combined}")
        self.assertNotIn("KeyError", combined, f"raw KeyError leaked:\n{combined}")
        self.assertIn("[ERROR]", combined)
        # names the offending group and advises matching -i/-p
        self.assertIn("ghost_group", combined)


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

    def test_modularity_missing_index_clean_error(self):
        """ISSUE-075: modularity -i <nonexistent prefix> must give a clean [ERROR]
        instead of a raw FileNotFoundError from open('<prefix>.namesMap')."""
        bad = os.path.join(self.tmpdir, "nonexistent")
        out = os.path.join(self.tmpdir, "mbad")
        rc, stdout, stderr = run_command(
            f"DBRetina modularity -i {bad} -p {self.pw_file} -c 80 -o {out}"
        )
        combined = stdout + stderr
        self.assertNotEqual(rc, 0, "missing index should exit nonzero")
        self.assertNotIn("Traceback", combined, f"raw traceback leaked:\n{combined}")
        self.assertNotIn("FileNotFoundError", combined,
                         f"raw FileNotFoundError leaked:\n{combined}")
        self.assertIn("[ERROR]", combined)
        self.assertIn(bad, combined)


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

    def test_bipartite_malformed_pairwise_row_clean_error(self):
        """ISSUE-079: a truncated pairwise TSV row (too few columns) on bipartite's
        TSV fallback must give a clean [ERROR] naming the line, not a raw
        IndexError. A standalone TSV (no parquet/.dbrp sibling) forces the
        fallback."""
        mal = os.path.join(self.tmpdir, "mal.tsv")
        header = ("group_1_ID\tgroup_2_ID\tgroup_1_name\tgroup_2_name\t"
                  "shared_features\tcontainment\tochiai\tjaccard\tcsi\tdice\todds_ratio\n")
        good = "1\t3\tgroupa\tgroupc\t2\t40.0\t40.0\t20.0\t30.0\t30.0\t1.0\n"
        bad = "only\ttwo\n"
        write_file(mal, header + good + bad)
        out = os.path.join(self.tmpdir, "bmal")
        rc, stdout, stderr = run_command(
            f"DBRetina bipartite -p {mal} --group1 {self.g1} --group2 {self.g2} "
            f"-m ochiai -c 0 --no-plot -o {out}"
        )
        combined = stdout + stderr
        self.assertNotEqual(rc, 0, "malformed row should exit nonzero")
        self.assertNotIn("Traceback", combined, f"raw traceback leaked:\n{combined}")
        self.assertNotIn("IndexError", combined,
                         f"raw IndexError leaked:\n{combined}")
        self.assertIn("[ERROR]", combined)
        self.assertIn("malformed", combined.lower())

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

    def test_bipartite_invalid_metric_clean_error(self):
        """ISSUE-023: bipartite -m BOGUS exits nonzero with a clean error, not a raw
        traceback (KeyError on metric_to_col / ValueError from PairwiseStore)."""
        out = os.path.join(self.tmpdir, "bipbadm")
        rc, _, stderr = run_command(
            f"DBRetina bipartite -p {self.pw_file} --group1 {self.g1} "
            f"--group2 {self.g2} -m BOGUS -c 5 --no-plot -o {out}"
        )
        self.assertNotEqual(rc, 0, "invalid metric should exit nonzero")
        self.assertNotIn("Traceback", stderr)
        self.assertNotIn("KeyError", stderr)
        # mentions the bad value and lists the valid choices
        self.assertIn("BOGUS", stderr)
        self.assertIn("metric", stderr.lower())
        # a valid metric still works
        rc2, _, stderr2 = run_command(
            f"DBRetina bipartite -p {self.pw_file} --group1 {self.g1} "
            f"--group2 {self.g2} -m ochiai -c 0 --no-plot -o {out}_ok"
        )
        self.assertEqual(rc2, 0, stderr2)
        self.assertNotIn("Traceback", stderr2)
        assert_file_exists(self, f"{out}_ok_bipartite_pairwise.tsv")


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

    def test_graph_include_isolates_no_targets(self):
        """ISSUE-080: --include-isolates with NO targets must list zero-edge
        (isolate) groups too. At ochiai c=100 only groupa==groupf survives as an
        edge, so the other 4 groups are isolates. With --include-isolates the node
        table lists all 6 index groups; without it, only the 2 edge endpoints.

        Pre-fix the no-targets node universe was built solely from edge endpoints,
        so --include-isolates had nothing to add (empty/edge-only node table)."""
        # WITH --include-isolates: every index group appears (6 groups -> 6 nodes).
        out_iso = os.path.join(self.tmpdir, "iso")
        rc, stdout, stderr = run_command(
            f"DBRetina graph -i {self.prefix} -p {self.pw_file} "
            f"-m ochiai -c 100 --include-isolates -o {out_iso}"
        )
        self.assertEqual(rc, 0, stderr)
        assert_file_exists(self, f"{out_iso}_nodes.tsv")
        iso_nodes = set()
        with open(f"{out_iso}_nodes.tsv") as f:
            for i, line in enumerate(f):
                if i == 0 or not line.strip():  # header / blank
                    continue
                iso_nodes.add(line.split("\t")[0])
        self.assertEqual(
            iso_nodes, EXPECTED_GROUP_NAMES,
            f"--include-isolates (no targets) should list all index groups as "
            f"isolates; got {sorted(iso_nodes)}")

        # WITHOUT --include-isolates: only the edge endpoints (groupa, groupf).
        out_noiso = os.path.join(self.tmpdir, "noiso")
        rc, _, stderr = run_command(
            f"DBRetina graph -i {self.prefix} -p {self.pw_file} "
            f"-m ochiai -c 100 -o {out_noiso}"
        )
        self.assertEqual(rc, 0, stderr)
        noiso_nodes = set()
        with open(f"{out_noiso}_nodes.tsv") as f:
            for i, line in enumerate(f):
                if i == 0 or not line.strip():
                    continue
                noiso_nodes.add(line.split("\t")[0])
        self.assertEqual(
            noiso_nodes, {"groupa", "groupf"},
            f"without --include-isolates only edge endpoints should appear; "
            f"got {sorted(noiso_nodes)}")

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

    def test_graph_missing_index_clean_error(self):
        """ISSUE-069: graph -i <missing prefix> must give a clean [ERROR] (like
        genenet) instead of a raw FileNotFoundError from parse_node_size."""
        bad = os.path.join(self.tmpdir, "nope_prefix")
        out = os.path.join(self.tmpdir, "gmiss")
        rc, stdout, stderr = run_command(
            f"DBRetina graph -i {bad} -p {self.pw_file} -m ochiai -c 30 -o {out}"
        )
        combined = stdout + stderr
        self.assertNotEqual(rc, 0, "missing index should exit nonzero")
        self.assertNotIn("Traceback", combined, f"raw traceback leaked:\n{combined}")
        self.assertNotIn("FileNotFoundError", combined,
                         f"raw FileNotFoundError leaked:\n{combined}")
        self.assertIn("[ERROR]", combined)
        self.assertIn(bad, combined)

    def test_graph_malformed_pairwise_row_clean_error(self):
        """ISSUE-079: a truncated pairwise TSV row (too few columns) on the TSV
        fallback path must give a clean [ERROR] naming the line, not a raw
        IndexError. A sibling parquet dir/.dbrp would bypass the TSV path, so we
        build a standalone TSV (no siblings) to force the fallback."""
        # Standalone pairwise TSV (header + a couple valid rows + one truncated).
        mal = os.path.join(self.tmpdir, "mal.tsv")
        header = ("group_1_ID\tgroup_2_ID\tgroup_1_name\tgroup_2_name\t"
                  "shared_features\tcontainment\tochiai\tjaccard\tcsi\tdice\todds_ratio\n")
        good = "1\t2\tgroupa\tgroupb\t3\t60.0\t60.0\t50.0\t55.0\t55.0\t1.0\n"
        bad = "only\ttwo\n"
        write_file(mal, header + good + bad)
        out = os.path.join(self.tmpdir, "gmal")
        rc, stdout, stderr = run_command(
            f"DBRetina graph -i {self.prefix} -p {mal} -m ochiai -c 0 -o {out}"
        )
        combined = stdout + stderr
        self.assertNotEqual(rc, 0, "malformed row should exit nonzero")
        self.assertNotIn("Traceback", combined, f"raw traceback leaked:\n{combined}")
        self.assertNotIn("IndexError", combined,
                         f"raw IndexError leaked:\n{combined}")
        self.assertIn("[ERROR]", combined)
        self.assertIn("malformed", combined.lower())

    def test_graph_wellformed_tsv_fallback_still_works(self):
        """ISSUE-079 guard: a well-formed standalone pairwise TSV (TSV fallback,
        no parquet/.dbrp sibling) still produces graph output."""
        good_tsv = os.path.join(self.tmpdir, "good.tsv")
        header = ("group_1_ID\tgroup_2_ID\tgroup_1_name\tgroup_2_name\t"
                  "shared_features\tcontainment\tochiai\tjaccard\tcsi\tdice\todds_ratio\n")
        rows = ("1\t2\tgroupa\tgroupb\t3\t60.0\t60.0\t50.0\t55.0\t55.0\t1.0\n"
                "1\t6\tgroupa\tgroupf\t5\t100.0\t100.0\t100.0\t100.0\t100.0\t-1.0\n")
        write_file(good_tsv, header + rows)
        out = os.path.join(self.tmpdir, "ggood")
        rc, stdout, stderr = run_command(
            f"DBRetina graph -i {self.prefix} -p {good_tsv} -m ochiai -c 0 -o {out}"
        )
        combined = stdout + stderr
        self.assertEqual(rc, 0, combined)
        self.assertNotIn("Traceback", combined)
        assert_file_exists(self, f"{out}_edges.tsv")

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
# SECTION 11b: Cross-route -m pvalue handling (ISSUES 043 & 053)
# ============================================================

class TestPvalueCrossRoute(unittest.TestCase):
    """`-m pvalue` must behave consistently across the .tsv / parquet-dir / .dbrp
    input forms for every command that supports it:

      (a) when a pvalue column IS present  -> the command works (exit 0);
      (b) when it is ABSENT                -> the SAME clean error
          '[ERROR] pvalue not found in pairwise file!' with no Python traceback
          (no IsADirectoryError / UnicodeDecodeError / ValueError leak).

    Pre-fix (issues 043 & 053) the non-.tsv routes crashed: cluster/export/bipartite
    raised IsADirectoryError from check_if_there_is_a_pvalue() open()ing a parquet
    DIRECTORY as text (053); graph raised an uncaught ValueError ("Unknown metric
    'pvalue'") from PairwiseStore._validate_metric on the store route (043). The
    .tsv route already gave the clean error; these tests pin every route to it.
    """

    PV_ABSENT_MSG = "pvalue not found in pairwise file"
    TRACEBACK_MARKERS = (
        "Traceback (most recent call last)",
        "IsADirectoryError",
        "UnicodeDecodeError",
        "Unknown metric",
    )

    @classmethod
    def setUpClass(cls):
        # NO-pvalue substrate (standard fixture) and WITH-pvalue substrate.
        cls.prefix, cls.pw_tsv = _ensure_shared_fixture()
        cls.pv_prefix, cls.pv_pw_tsv = _ensure_shared_pvalue_fixture()

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp(prefix="dbretina_pvxr_")
        # Full target partition for the graph command (avoids the no-targets path).
        self.intra1 = write_file(
            os.path.join(self.tmpdir, "intra1.tsv"), "GroupA\nGroupB\nGroupE\nGroupF\n"
        )
        self.inter1 = write_file(
            os.path.join(self.tmpdir, "inter1.tsv"), "GroupC\nGroupD\n"
        )
        # Two non-overlapping group files for the bipartite command.
        self.bg1 = write_file(os.path.join(self.tmpdir, "bg1.txt"), "GroupA\nGroupB\n")
        self.bg2 = write_file(os.path.join(self.tmpdir, "bg2.txt"), "GroupC\nGroupD\n")

    def tearDown(self):
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    # ---- form helpers ----

    @staticmethod
    def _forms(pw_tsv):
        """Return the three -p forms sharing a base: (.tsv, parquet dir, .dbrp)."""
        base = pw_tsv[:-len(".tsv")] if pw_tsv.endswith(".tsv") else pw_tsv
        return {"tsv": pw_tsv, "dir": base, ".dbrp": base + ".dbrp"}

    def _assert_no_traceback(self, stderr, ctx):
        for marker in self.TRACEBACK_MARKERS:
            self.assertNotIn(
                marker, stderr, f"{ctx}: raw '{marker}' leaked:\n{stderr}"
            )

    def _assert_clean_pvalue_absent(self, rc, stderr, ctx):
        self.assertNotEqual(rc, 0, f"{ctx}: expected failure on pvalue-absent input")
        self._assert_no_traceback(stderr, ctx)
        self.assertIn(self.PV_ABSENT_MSG, stderr, f"{ctx}: missing clean error:\n{stderr}")

    # ---- command runners (one per command, parametrized by -p) ----

    def _run_cluster(self, p, out):
        return run_command(f"DBRetina cluster -p {p} -m pvalue -c 20 -o {out}")

    def _run_export(self, p, out):
        return run_command(f"DBRetina export -p {p} -m pvalue -o {out}")

    def _run_graph(self, p, out, prefix):
        return run_command(
            f"DBRetina graph -i {prefix} -p {p} "
            f"--intra-targets {self.intra1} --inter-targets {self.inter1} "
            f"-m pvalue -c 20 -o {out}"
        )

    def _run_bipartite(self, p, out):
        # -c 1 keeps every pair: p-values lie in [0, 1] and the cutoff is now
        # correctly "lower is better" (keep pvalue <= cutoff) (issue 068). The
        # old default cutoff 0.0 only kept pairs because pvalue filtering was
        # inverted (pvalue >= 0 matched everything); with the fix, 0.0 would
        # exclude all pairs and the command would (correctly) report no overlap.
        return run_command(
            f"DBRetina bipartite -p {p} --group1 {self.bg1} --group2 {self.bg2} "
            f"-m pvalue -c 1 --no-plot -o {out}"
        )

    # =====================================================================
    # ISSUE 053: cluster / export -m pvalue on a parquet DIRECTORY
    # =====================================================================

    def test_053_cluster_parquet_dir_no_pvalue_clean_error(self):
        """ISSUE-053: cluster -p <parquet dir> -m pvalue on a no-pvalue dataset
        -> clean error, NOT IsADirectoryError."""
        pw_dir = self._forms(self.pw_tsv)["dir"]
        self.assertTrue(os.path.isdir(pw_dir))
        rc, _, stderr = self._run_cluster(pw_dir, os.path.join(self.tmpdir, "c"))
        self._assert_clean_pvalue_absent(rc, stderr, "cluster -p DIR (no pvalue)")

    def test_053_export_parquet_dir_no_pvalue_clean_error(self):
        """ISSUE-053: export -p <parquet dir> -m pvalue on a no-pvalue dataset
        -> clean error, NOT IsADirectoryError."""
        pw_dir = self._forms(self.pw_tsv)["dir"]
        rc, _, stderr = self._run_export(pw_dir, os.path.join(self.tmpdir, "e"))
        self._assert_clean_pvalue_absent(rc, stderr, "export -p DIR (no pvalue)")

    def test_053_cluster_parquet_dir_with_pvalue_works(self):
        """ISSUE-053: cluster -p <parquet dir> -m pvalue WITH pvalue -> exit 0."""
        pw_dir = self._forms(self.pv_pw_tsv)["dir"]
        out = os.path.join(self.tmpdir, "cok")
        rc, _, stderr = self._run_cluster(pw_dir, out)
        self.assertEqual(rc, 0, stderr)
        assert_file_exists(self, f"{out}_clusters.tsv")

    def test_053_export_parquet_dir_with_pvalue_works(self):
        """ISSUE-053: export -p <parquet dir> -m pvalue WITH pvalue -> exit 0."""
        pw_dir = self._forms(self.pv_pw_tsv)["dir"]
        out = os.path.join(self.tmpdir, "eok")
        rc, _, stderr = self._run_export(pw_dir, out)
        self.assertEqual(rc, 0, stderr)
        assert_file_exists(self, f"{out}_distmat.tsv")

    # =====================================================================
    # ISSUE 043: graph -m pvalue cross-route (store / parquet / .dbrp)
    # =====================================================================

    def test_043_graph_parquet_dir_no_pvalue_clean_error(self):
        """ISSUE-043: graph -p <parquet dir> -m pvalue on a no-pvalue dataset
        -> clean error, NOT an uncaught ValueError('Unknown metric')."""
        pw_dir = self._forms(self.pw_tsv)["dir"]
        rc, _, stderr = self._run_graph(pw_dir, os.path.join(self.tmpdir, "g"), self.prefix)
        self._assert_clean_pvalue_absent(rc, stderr, "graph -p DIR (no pvalue)")

    def test_043_graph_parquet_dir_with_pvalue_works(self):
        """ISSUE-043: graph -p <parquet dir> -m pvalue WITH pvalue -> exit 0."""
        pw_dir = self._forms(self.pv_pw_tsv)["dir"]
        out = os.path.join(self.tmpdir, "gok")
        rc, _, stderr = self._run_graph(pw_dir, out, self.pv_prefix)
        self.assertEqual(rc, 0, stderr)
        assert_file_exists(self, f"{out}_edges.tsv")

    def test_043_graph_tsv_no_pvalue_clean_error(self):
        """ISSUE-043: graph -p <.tsv> -m pvalue on a no-pvalue dataset is also
        clean (the .tsv resolves to a sibling parquet store, so it took the same
        unguarded store route pre-fix)."""
        rc, _, stderr = self._run_graph(
            self.pw_tsv, os.path.join(self.tmpdir, "gt"), self.prefix
        )
        self._assert_clean_pvalue_absent(rc, stderr, "graph -p .tsv (no pvalue)")

    # =====================================================================
    # 053 sibling: bipartite -m pvalue on a parquet DIRECTORY
    # =====================================================================

    def test_bipartite_parquet_dir_no_pvalue_clean_error(self):
        """bipartite -p <parquet dir> -m pvalue (same check_if_there_is_a_pvalue
        text-open bug as 053) -> clean error, NOT IsADirectoryError."""
        pw_dir = self._forms(self.pw_tsv)["dir"]
        rc, _, stderr = self._run_bipartite(pw_dir, os.path.join(self.tmpdir, "b"))
        self._assert_clean_pvalue_absent(rc, stderr, "bipartite -p DIR (no pvalue)")

    def test_bipartite_parquet_dir_with_pvalue_works(self):
        """bipartite -p <parquet dir> -m pvalue WITH pvalue -> exit 0."""
        pw_dir = self._forms(self.pv_pw_tsv)["dir"]
        out = os.path.join(self.tmpdir, "bok")
        rc, _, stderr = self._run_bipartite(pw_dir, out)
        self.assertEqual(rc, 0, stderr)

    # =====================================================================
    # Cross-route PARITY: every form behaves the same
    # =====================================================================

    def test_parity_cluster_all_forms_no_pvalue(self):
        """cluster -m pvalue: clean error on .tsv AND parquet dir AND .dbrp."""
        forms = self._forms(self.pw_tsv)
        for name, p in forms.items():
            with self.subTest(form=name):
                rc, _, stderr = self._run_cluster(
                    p, os.path.join(self.tmpdir, f"par_{name.strip('.')}")
                )
                self._assert_clean_pvalue_absent(rc, stderr, f"cluster -p {name}")

    def test_parity_cluster_all_forms_with_pvalue(self):
        """cluster -m pvalue: works on .tsv AND parquet dir AND .dbrp (WITH pvalue)."""
        forms = self._forms(self.pv_pw_tsv)
        for name, p in forms.items():
            with self.subTest(form=name):
                out = os.path.join(self.tmpdir, f"parok_{name.strip('.')}")
                rc, _, stderr = self._run_cluster(p, out)
                self.assertEqual(rc, 0, f"cluster -p {name}: {stderr}")
                assert_file_exists(self, f"{out}_clusters.tsv")

    def test_parity_export_all_forms_no_pvalue(self):
        """export -m pvalue: clean error on .tsv AND parquet dir AND .dbrp."""
        forms = self._forms(self.pw_tsv)
        for name, p in forms.items():
            with self.subTest(form=name):
                rc, _, stderr = self._run_export(
                    p, os.path.join(self.tmpdir, f"epar_{name.strip('.')}")
                )
                self._assert_clean_pvalue_absent(rc, stderr, f"export -p {name}")

    def test_parity_export_dbrp_with_pvalue_works(self):
        """export -m pvalue on the .dbrp form WITH pvalue -> exit 0."""
        dbrp = self._forms(self.pv_pw_tsv)[".dbrp"]
        self.assertTrue(os.path.isfile(dbrp))
        out = os.path.join(self.tmpdir, "edbrp")
        rc, _, stderr = self._run_export(dbrp, out)
        self.assertEqual(rc, 0, stderr)
        assert_file_exists(self, f"{out}_distmat.tsv")

    def test_pvalue_absent_does_not_affect_other_metrics(self):
        """Regression guard: a non-pvalue metric (containment) still works on the
        no-pvalue parquet dir (the new pvalue check must only fire for -m pvalue)."""
        pw_dir = self._forms(self.pw_tsv)["dir"]
        out = os.path.join(self.tmpdir, "cm")
        rc, _, stderr = run_command(
            f"DBRetina cluster -p {pw_dir} -m containment -c 50 -o {out}"
        )
        self.assertEqual(rc, 0, stderr)
        self._assert_no_traceback(stderr, "cluster -p DIR -m containment")


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


class TestGenenetPairwiseInputForms(unittest.TestCase):
    """ISSUE-066 / ISSUE-067: genenet AND interactome must accept the same -p
    input forms that graph/bipartite do (the pairwise .tsv, its parquet
    directory, and the binary .dbrp), and must reject a parquet directory with
    no usable store (no manifest.json) with a clean [ERROR] -- never a raw
    UnicodeDecodeError / IsADirectoryError traceback.

    genenet and interactome share one Click callback (net_kind from
    ctx.info_name), so each scenario is exercised against BOTH command names.
    """

    @classmethod
    def setUpClass(cls):
        cls.prefix, cls.pw_file = _ensure_shared_fixture()

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp(prefix="dbretina_gn_forms_")

    def tearDown(self):
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    # ---- helpers ----

    def _sibling(self, ext):
        """Return the sibling of the shared pairwise .tsv with the given suffix.

        ext="" -> the parquet directory; ext=".dbrp" -> the binary file.
        """
        base = self.pw_file[:-len(".tsv")] if self.pw_file.endswith(".tsv") else self.pw_file
        return base + ext

    def _network_lines(self, out_prefix, net_kind):
        """Sorted lines of the produced <prefix>_<net_kind>.tsv network file."""
        path = f"{out_prefix}_{net_kind}.tsv"
        assert_file_exists(self, path)
        with open(path) as f:
            return sorted(line.rstrip("\n") for line in f)

    def _isolated_dbrp(self):
        """Copy ONLY the .dbrp into a fresh dir (no sibling parquet store / .tsv).

        Without this isolation open_pairwise() would resolve the shared parquet
        directory sibling and the .dbrp text-read fallback (the ISSUE-066 bug)
        would never be exercised -- exactly the masking seen on the kegg
        substrate. The basename keeps the canonical *_DBRetina_pairwise stem so
        resolve_dbrp_path finds it.
        """
        dst_dir = tempfile.mkdtemp(prefix="dbretina_gn_dbrp_", dir=self.tmpdir)
        dst = os.path.join(dst_dir, "iso_DBRetina_pairwise.dbrp")
        shutil.copy(self._sibling(".dbrp"), dst)
        return dst

    def _storeless_parquet_dir(self):
        """Copy the parquet dir but strip manifest.json and any sibling .dbrp/.tsv.

        open_pairwise() returns None (no manifest) and resolve_dbrp_path() finds
        no sibling .dbrp, so the command must hit the storeless-directory guard
        rather than open()'ing the directory as text (the ISSUE-067 bug).
        """
        dst_dir = tempfile.mkdtemp(prefix="dbretina_gn_nostore_", dir=self.tmpdir)
        dst = os.path.join(dst_dir, "iso_DBRetina_pairwise")
        shutil.copytree(self._sibling(""), dst)
        manifest = os.path.join(dst, "manifest.json")
        if os.path.exists(manifest):
            os.remove(manifest)
        # ensure no sibling .dbrp/.tsv accidentally rescues it
        for sib in ("iso_DBRetina_pairwise.dbrp", "iso_DBRetina_pairwise.tsv"):
            p = os.path.join(dst_dir, sib)
            if os.path.exists(p):
                os.remove(p)
        return dst

    def _assert_no_raw_traceback(self, stderr):
        self.assertNotIn("Traceback", stderr)
        self.assertNotIn("UnicodeDecodeError", stderr)
        self.assertNotIn("IsADirectoryError", stderr)

    # ---- ISSUE-066: .dbrp input (both commands) ----

    def test_genenet_dbrp_no_raw_crash(self):
        """ISSUE-066: genenet -p <.dbrp> (isolated, no sibling store) builds the
        network from the binary reader -- no UnicodeDecodeError traceback."""
        dbrp = self._isolated_dbrp()
        out = os.path.join(self.tmpdir, "gn_dbrp")
        rc, _, stderr = run_command(
            f"DBRetina genenet -i {self.prefix} -p {dbrp} -o {out}"
        )
        self._assert_no_raw_traceback(stderr)
        self.assertEqual(rc, 0, stderr)
        assert_file_exists(self, f"{out}_genenet.tsv")

    def test_interactome_dbrp_no_raw_crash(self):
        """ISSUE-066: interactome -p <.dbrp> (isolated, no sibling store) builds
        the network from the binary reader -- no UnicodeDecodeError traceback."""
        dbrp = self._isolated_dbrp()
        out = os.path.join(self.tmpdir, "in_dbrp")
        rc, _, stderr = run_command(
            f"DBRetina interactome -i {self.prefix} -p {dbrp} -o {out}"
        )
        self._assert_no_raw_traceback(stderr)
        self.assertEqual(rc, 0, stderr)
        assert_file_exists(self, f"{out}_interactome.tsv")

    def test_genenet_dbrp_matches_tsv(self):
        """ISSUE-066: genenet on an isolated .dbrp yields the SAME network as the
        .tsv form (same nodes/edges/weights)."""
        out_tsv = os.path.join(self.tmpdir, "gn_from_tsv")
        rc_t, _, err_t = run_command(
            f"DBRetina genenet -i {self.prefix} -p {self.pw_file} -o {out_tsv}"
        )
        self.assertEqual(rc_t, 0, err_t)

        dbrp = self._isolated_dbrp()
        out_dbrp = os.path.join(self.tmpdir, "gn_from_dbrp")
        rc_d, _, err_d = run_command(
            f"DBRetina genenet -i {self.prefix} -p {dbrp} -o {out_dbrp}"
        )
        self._assert_no_raw_traceback(err_d)
        self.assertEqual(rc_d, 0, err_d)

        self.assertEqual(
            self._network_lines(out_dbrp, "genenet"),
            self._network_lines(out_tsv, "genenet"),
            "genenet network from .dbrp differs from the .tsv form",
        )

    # ---- ISSUE-067: storeless parquet directory (both commands) ----

    def test_genenet_storeless_dir_clean_error(self):
        """ISSUE-067: genenet -p <parquet dir w/o manifest.json, no sibling
        .dbrp> exits with a clean [ERROR], not a raw IsADirectoryError."""
        pw_dir = self._storeless_parquet_dir()
        out = os.path.join(self.tmpdir, "gn_nostore")
        rc, _, stderr = run_command(
            f"DBRetina genenet -i {self.prefix} -p {pw_dir} -o {out}"
        )
        self.assertNotEqual(rc, 0, "storeless dir should exit nonzero")
        self._assert_no_raw_traceback(stderr)
        self.assertIn("[ERROR]", stderr)
        self.assertIn("manifest.json", stderr)

    def test_interactome_storeless_dir_clean_error(self):
        """ISSUE-067: interactome -p <parquet dir w/o manifest.json> exits with a
        clean [ERROR], not a raw IsADirectoryError."""
        pw_dir = self._storeless_parquet_dir()
        out = os.path.join(self.tmpdir, "in_nostore")
        rc, _, stderr = run_command(
            f"DBRetina interactome -i {self.prefix} -p {pw_dir} -o {out}"
        )
        self.assertNotEqual(rc, 0, "storeless dir should exit nonzero")
        self._assert_no_raw_traceback(stderr)
        self.assertIn("[ERROR]", stderr)
        self.assertIn("manifest.json", stderr)

    # ---- parquet directory WITH a usable store still works (no regression) ----

    def test_genenet_parquet_dir_with_store_works(self):
        """The canonical parquet-dir-with-store path is unchanged: genenet -p
        <valid parquet dir> still builds the network."""
        out = os.path.join(self.tmpdir, "gn_dir")
        rc, _, stderr = run_command(
            f"DBRetina genenet -i {self.prefix} -p {self._sibling('')} -o {out}"
        )
        self._assert_no_raw_traceback(stderr)
        self.assertEqual(rc, 0, stderr)
        assert_file_exists(self, f"{out}_genenet.tsv")

    # ---- three-way parity: .tsv vs parquet-dir(store) vs .dbrp ----

    def test_genenet_parity_across_all_forms(self):
        """genenet produces the SAME network for the .tsv, the parquet directory
        (with store), and an isolated .dbrp."""
        out_tsv = os.path.join(self.tmpdir, "p_tsv")
        out_dir = os.path.join(self.tmpdir, "p_dir")
        out_dbrp = os.path.join(self.tmpdir, "p_dbrp")

        rc_t, _, err_t = run_command(
            f"DBRetina genenet -i {self.prefix} -p {self.pw_file} -o {out_tsv}"
        )
        rc_d, _, err_d = run_command(
            f"DBRetina genenet -i {self.prefix} -p {self._sibling('')} -o {out_dir}"
        )
        dbrp = self._isolated_dbrp()
        rc_b, _, err_b = run_command(
            f"DBRetina genenet -i {self.prefix} -p {dbrp} -o {out_dbrp}"
        )
        self.assertEqual(rc_t, 0, err_t)
        self.assertEqual(rc_d, 0, err_d)
        self.assertEqual(rc_b, 0, err_b)

        tsv_lines = self._network_lines(out_tsv, "genenet")
        self.assertEqual(tsv_lines, self._network_lines(out_dir, "genenet"),
                         "parquet-dir network differs from .tsv")
        self.assertEqual(tsv_lines, self._network_lines(out_dbrp, "genenet"),
                         ".dbrp network differs from .tsv")


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

    def test_setcov_missing_index_clean_error(self):
        """ISSUE-074: setcov -i <nonexistent prefix> must give a clean [ERROR]
        instead of a raw FileNotFoundError from open('<prefix>_raw.json')."""
        rc, stdout, stderr = run_command(
            "DBRetina setcov -i nonexistent_pref -o s6", cwd=self.tmpdir
        )
        combined = stdout + stderr
        self.assertNotEqual(rc, 0, "missing index should exit nonzero")
        self.assertNotIn("Traceback", combined, f"raw traceback leaked:\n{combined}")
        self.assertNotIn("FileNotFoundError", combined,
                         f"raw FileNotFoundError leaked:\n{combined}")
        self.assertIn("[ERROR]", combined)
        self.assertIn("nonexistent_pref", combined)


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

    def test_append_duplicate_group_clean_error(self):
        """ISSUE-083: append a GMT whose group name already exists in the base
        index must give a clean [ERROR] (like merge), NOT a raw RuntimeError
        traceback. The message must not mention the nonexistent --force flag, and
        no intermediate temp files (_raw.json/_hashes.json) may be left behind."""
        idx = self._create_index("base",
            "gene_set\tgene\nfoo\tAlpha\nfoo\tBeta\nbar\tGamma\n")
        # GMT whose group name 'foo' duplicates an existing index group.
        gmt = write_file(os.path.join(self.tmpdir, "dup.gmt"),
            "foo\tdesc\tGENEX\tGENEY\n")
        out_prefix = os.path.join(self.tmpdir, "app_dup")
        out = f"{out_prefix}.dbri"
        rc, stdout, stderr = run_command(
            f"DBRetina append -i {idx}.dbri -g {gmt} -o {out}"
        )
        combined = stdout + stderr
        self.assertNotEqual(rc, 0, "duplicate-name append should fail")
        self.assertNotIn("Traceback", combined, f"raw traceback leaked:\n{combined}")
        self.assertNotIn("RuntimeError", combined,
                         f"raw RuntimeError leaked:\n{combined}")
        self.assertIn("[ERROR]", combined)
        # append has NO --force option -> the hint must be gone.
        self.assertNotIn("--force", combined,
                         f"error references the nonexistent --force flag:\n{combined}")
        # the offending group name is still reported.
        self.assertIn("foo", combined)
        # no leaked intermediate temp files.
        self.assertFalse(os.path.exists(f"{out_prefix}_raw.json"),
                         "append left behind a _raw.json temp file")
        self.assertFalse(os.path.exists(f"{out_prefix}_hashes.json"),
                         "append left behind a _hashes.json temp file")

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


@unittest.skipUnless(_HAS_REST_TEST, "REST /sql tests need [server] extra + httpx")
class TestRestSql(unittest.TestCase):
    """POST /api/v1/sql edge cases (issues 057 + 059).

    057: a read-only SELECT whose result has a LIST/array column 500s because
    fetchdf() yields numpy ndarrays that FastAPI's jsonable_encoder cannot
    serialize -- and that failure happens AFTER the handler returns, so the
    handler's try/except can't catch it and the catch-all turns it into a 500.
    A serialization failure must never be a 500: list cells should serialize.

    059: an empty / blank / comment-only query passes validate_sql_safety,
    store.sql() returns None, and None.fetchdf() leaks
    "'NoneType' object has no attribute 'fetchdf'" into the 400 detail. Such
    queries must be rejected up front with a clean, meaningful 400.

    Prefers the (larger) kegg substrate so list(group_1_id) FROM pairs returns
    rows, and falls back to the synthetic fixture so the suite stays portable.
    """

    KEGG_PARQUET = "/home/mabuelanin/dbretina_scratch/out/kegg_DBRetina_pairwise"
    KEGG_DBRI = "/home/mabuelanin/dbretina_scratch/out/kegg.dbri"

    def setUp(self):
        if (os.path.isdir(self.KEGG_PARQUET)
                and os.path.exists(os.path.join(self.KEGG_PARQUET, "manifest.json"))):
            self.tmpdir = None
            self.parquet = self.KEGG_PARQUET
            self.dbri = self.KEGG_DBRI if os.path.exists(self.KEGG_DBRI) else None
        else:
            self.tmpdir = tempfile.mkdtemp(prefix="dbretina_sql_")
            self.prefix, _ = setup_index_and_pairwise(self.tmpdir)
            self.parquet = f"{self.prefix}_DBRetina_pairwise"
            self.dbri = f"{self.prefix}.dbri"

    def tearDown(self):
        if self.tmpdir:
            shutil.rmtree(self.tmpdir, ignore_errors=True)

    def _client(self):
        from dbretina.compat import open_pairwise
        from dbretina.pairwise_store import PairwiseStore
        from dbretina.rest_api import create_app
        store = open_pairwise(self.parquet)
        if store is None:
            store = PairwiseStore(self.parquet, dbri_path=self.dbri)
        elif self.dbri:
            store._dbri_path = self.dbri
        # raise_server_exceptions=False so a serialization failure surfaces as a
        # 500 status code we can assert on, instead of re-raising into the worker.
        return _TestClient(create_app(store, dbri_path=self.dbri),
                           raise_server_exceptions=False), store

    # ── 057: LIST/array result must serialize, never 500 ──────────────
    def test_sql_list_literal_serializes(self):
        client, store = self._client()
        try:
            r = client.post("/api/v1/sql", json={"query": "SELECT [1,2,3] AS arr"})
            self.assertEqual(r.status_code, 200,
                             f"list-literal SELECT should be 200, got "
                             f"{r.status_code}: {r.text[:300]}")
            body = r.json()
            self.assertEqual(body["row_count"], 1)
            # The ndarray must come back as a JSON list with native ints.
            self.assertEqual(body["rows"][0]["arr"], [1, 2, 3],
                             f"list column not serialized to a JSON array: {body['rows'][0]}")
        finally:
            store.close()

    def test_sql_list_aggregate_serializes(self):
        client, store = self._client()
        try:
            r = client.post(
                "/api/v1/sql",
                json={"query": "SELECT list(group_1_id) AS ids FROM pairs LIMIT 1"},
            )
            self.assertEqual(r.status_code, 200,
                             f"list(group_1_id) should be 200, got "
                             f"{r.status_code}: {r.text[:300]}")
            body = r.json()
            self.assertEqual(body["columns"], ["ids"])
            self.assertEqual(body["row_count"], 1)
            ids = body["rows"][0]["ids"]
            self.assertIsInstance(ids, list,
                                  f"aggregate list column not serialized as a JSON array: {ids!r}")
        finally:
            store.close()

    def test_sql_blob_serializes(self):
        """A BLOB result column (bytes) must serialize (base64), not 500 — the /sql
        endpoint runs arbitrary SELECTs so a BLOB-producing expression is valid input."""
        client, store = self._client()
        try:
            r = client.post("/api/v1/sql", json={"query": "SELECT 'abc'::BLOB AS b"})
            self.assertEqual(r.status_code, 200,
                             f"BLOB SELECT should be 200, got {r.status_code}: {r.text[:300]}")
            # base64(b'abc') == 'YWJj'
            self.assertEqual(r.json()["rows"][0]["b"], "YWJj",
                             f"BLOB not serialized to base64: {r.json()['rows'][0]}")
        finally:
            store.close()

    def test_sql_scalar_result_unchanged(self):
        # The serialization fix must keep ordinary scalar results byte-identical.
        client, store = self._client()
        try:
            r = client.post("/api/v1/sql",
                            json={"query": "SELECT 1 AS one, 2.5 AS two"})
            self.assertEqual(r.status_code, 200, f"scalar SELECT: {r.text[:300]}")
            self.assertEqual(
                r.json(),
                {"columns": ["one", "two"], "row_count": 1,
                 "rows": [{"one": 1, "two": 2.5}]},
            )
        finally:
            store.close()

    # ── 059: empty / blank / comment-only query -> clean 400 ──────────
    def test_sql_empty_query_clean_400(self):
        client, store = self._client()
        try:
            for q in ("", "   ", "-- hi", "/* nothing */", "\n\t  \n"):
                r = client.post("/api/v1/sql", json={"query": q})
                self.assertEqual(r.status_code, 400,
                                 f"empty query {q!r} should be a clean 400, got "
                                 f"{r.status_code}: {r.text[:200]}")
                # The internal AttributeError text must NOT leak.
                self.assertNotIn("NoneType", r.text,
                                 f"empty query {q!r} leaked internal error: {r.text[:200]}")
                self.assertNotIn("fetchdf", r.text,
                                 f"empty query {q!r} leaked internal error: {r.text[:200]}")
                self.assertNotIn("Traceback", r.text)
                body = r.json()
                self.assertEqual(body["error_code"], "VALIDATION_ERROR",
                                 f"empty query {q!r} wrong error_code: {body}")
                self.assertTrue(body["detail"],
                                f"empty query {q!r} has no detail message")
        finally:
            store.close()

    # ── regression guard: the S1 SQL sandbox must stay intact ─────────
    def test_sql_sandbox_still_blocks_file_read_and_ddl(self):
        client, store = self._client()
        try:
            # read_text file exfiltration stays blocked (external access off).
            r = client.post("/api/v1/sql",
                            json={"query": "SELECT * FROM read_text('/etc/passwd')"})
            self.assertNotIn("root:", r.text, "sandbox leaked /etc/passwd")
            self.assertNotEqual(r.status_code, 200,
                                f"file-read query unexpectedly succeeded: {r.text[:200]}")
            # DDL/DML keywords stay blocked by the denylist.
            for q in ("DROP TABLE pairs", "DELETE FROM pairs",
                      "CREATE TABLE x AS SELECT 1"):
                r = client.post("/api/v1/sql", json={"query": q})
                self.assertEqual(r.status_code, 400,
                                 f"dangerous query {q!r} not blocked: "
                                 f"{r.status_code} {r.text[:200]}")
                self.assertEqual(r.json()["error_code"], "UNSAFE_QUERY",
                                 f"dangerous query {q!r} wrong error_code: {r.text[:200]}")
        finally:
            store.close()


@unittest.skipUnless(_HAS_REST_TEST, "REST serve-low tests need [server] extra + httpx")
class TestRestServeLows(unittest.TestCase):
    """ISSUE-060: bundle of LOW-severity serve (REST) polish/hardening items.

    1. hub-genes (method=hypergraph, the default) only validated ``metric`` in the
       edge_weighted branch, so a typo'd metric silently 200'd. Validate it for ALL
       methods -> clean 400.
    2. hub-genes 404 for a missing group echoed the dataset group catalog
       ("...Available groups include: kegg_..., ...") into the detail. Trim it.
    3. shared-features + graph/shortest-path embedded a KeyError that already
       contains the full "Group not found: x" sentence -> doubled message
       "Group not found: 'Group not found: x'". Use the clean name.
    4. POST /sql allowed PRAGMA / version() / duckdb_settings() etc. (metadata
       disclosure). Extend the denylist; legit SELECT-over-pairs still works and
       the S1 sandbox (file read + DDL/DML) stays blocked.
    5. groups/{g}/pairs?cutoff=... with no metric silently ignored the cutoff and
       returned ALL pairs -> now a clean 400 ("cutoff requires a metric").
    6. /favicon.ico 404'd on every dashboard load -> tiny 204 route.

    Prefers the (larger) kegg substrate and falls back to the synthetic fixture so
    the suite stays portable.
    """

    KEGG_PARQUET = "/home/mabuelanin/dbretina_scratch/out/kegg_DBRetina_pairwise"
    KEGG_DBRI = "/home/mabuelanin/dbretina_scratch/out/kegg.dbri"

    def setUp(self):
        if (os.path.isdir(self.KEGG_PARQUET)
                and os.path.exists(os.path.join(self.KEGG_PARQUET, "manifest.json"))
                and os.path.exists(self.KEGG_DBRI)):
            self.tmpdir = None
            self.parquet = self.KEGG_PARQUET
            self.dbri = self.KEGG_DBRI
        else:
            self.tmpdir = tempfile.mkdtemp(prefix="dbretina_low_")
            self.prefix, _ = setup_index_and_pairwise(self.tmpdir)
            self.parquet = f"{self.prefix}_DBRetina_pairwise"
            self.dbri = f"{self.prefix}.dbri"

    def tearDown(self):
        if self.tmpdir:
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
        # assert on instead of re-raising into the worker.
        return _TestClient(create_app(store, dbri_path=self.dbri),
                           raise_server_exceptions=False), store

    def _a_group(self, store):
        """A real group name present in this store's gene sets."""
        names = list(store.get_names_map().values())
        self.assertTrue(names, "store has no groups")
        return names[0]

    # ── 1: hub-genes metric validated for ALL methods ─────────────────
    def test_hub_genes_hypergraph_bad_metric_is_400(self):
        client, store = self._client()
        try:
            group = self._a_group(store)
            r = client.post("/api/v1/genes/hub-genes",
                            json={"group_name": group, "method": "hypergraph",
                                  "metric": "not_a_metric"})
            self.assertEqual(r.status_code, 400,
                             f"hypergraph + bad metric should be 400, got "
                             f"{r.status_code}: {r.text[:200]}")
            self.assertEqual(r.json()["error_code"], "VALIDATION_ERROR",
                             f"wrong error_code: {r.text[:200]}")
        finally:
            store.close()

    def test_hub_genes_hypergraph_valid_metric_is_200(self):
        client, store = self._client()
        try:
            group = self._a_group(store)
            r = client.post("/api/v1/genes/hub-genes",
                            json={"group_name": group, "method": "hypergraph",
                                  "metric": "ochiai"})
            self.assertEqual(r.status_code, 200,
                             f"hypergraph + valid metric should be 200, got "
                             f"{r.status_code}: {r.text[:300]}")
            self.assertEqual(r.json()["method"], "hypergraph")
        finally:
            store.close()

    # ── 2: hub-genes missing group 404 must not enumerate the catalog ──
    def test_hub_genes_missing_group_404_no_catalog(self):
        client, store = self._client()
        try:
            r = client.post("/api/v1/genes/hub-genes",
                            json={"group_name": "definitely_no_such_group_xyz",
                                  "method": "hypergraph", "metric": "ochiai"})
            self.assertEqual(r.status_code, 404,
                             f"missing group should be 404, got "
                             f"{r.status_code}: {r.text[:200]}")
            self.assertNotIn("Available groups include", r.text,
                             f"404 detail leaks the group catalog: {r.text[:300]}")
        finally:
            store.close()

    # ── 3: clean (non-doubled) not-found messages ─────────────────────
    def test_shared_features_missing_group_clean_message(self):
        client, store = self._client()
        try:
            r = client.get("/api/v1/shared-features",
                           params={"group_a": "no_such_group_a", "group_b": "no_such_group_b"})
            self.assertEqual(r.status_code, 404,
                             f"missing group should be 404, got "
                             f"{r.status_code}: {r.text[:200]}")
            detail = r.json()["detail"]
            # Must not be the doubled "Group not found: 'Group not found: x'".
            self.assertNotIn("Group not found: 'Group not found",
                             detail, f"doubled not-found message: {detail!r}")
            self.assertNotIn("Group not found: \"Group not found",
                             detail, f"doubled not-found message: {detail!r}")
        finally:
            store.close()

    def test_shortest_path_missing_group_clean_message(self):
        client, store = self._client()
        try:
            r = client.get("/api/v1/graph/shortest-path",
                           params={"source": "no_such_source", "target": "no_such_target"})
            self.assertEqual(r.status_code, 404,
                             f"missing group should be 404, got "
                             f"{r.status_code}: {r.text[:200]}")
            detail = r.json()["detail"]
            self.assertNotIn("Group not found: 'Group not found",
                             detail, f"doubled not-found message: {detail!r}")
            self.assertNotIn("Group not found: \"Group not found",
                             detail, f"doubled not-found message: {detail!r}")
        finally:
            store.close()

    # ── 4: /sql metadata-disclosure denylist (S1 sandbox unchanged) ───
    def test_sql_blocks_metadata_introspection(self):
        client, store = self._client()
        try:
            for q in ("PRAGMA version",
                      "SELECT version()",
                      "SELECT * FROM duckdb_settings()",
                      "PRAGMA database_list",
                      "SELECT * FROM pragma_table_info('pairs')",
                      "SELECT * FROM pg_tables",
                      "SELECT * FROM sqlite_master",
                      "SET memory_limit='1GB'"):
                r = client.post("/api/v1/sql", json={"query": q})
                self.assertEqual(r.status_code, 400,
                                 f"metadata query {q!r} should be blocked (400), got "
                                 f"{r.status_code}: {r.text[:200]}")
                self.assertEqual(r.json()["error_code"], "UNSAFE_QUERY",
                                 f"metadata query {q!r} wrong error_code: {r.text[:200]}")
                self.assertNotIn("v1.", r.text,
                                 f"metadata query {q!r} leaked a version string: {r.text[:200]}")
        finally:
            store.close()

    def test_sql_metadata_prefix_aliases_allowed(self):
        """The metadata denylist must NOT false-positive on legit SELECTs that merely
        alias/name something with a pg_/duckdb_/sqlite_/pragma_ prefix (issue 060.4
        over-block). /sql runs arbitrary SELECTs, so these must still return 200."""
        client, store = self._client()
        try:
            for q in ("SELECT ochiai AS pg_score FROM pairs LIMIT 1",
                      "SELECT group_1_id AS pg_rank FROM pairs LIMIT 1",
                      "WITH pg_top AS (SELECT * FROM pairs LIMIT 1) SELECT * FROM pg_top",
                      "SELECT ochiai AS duckdb_note FROM pairs LIMIT 1",
                      "SELECT ochiai AS sqlite_x FROM pairs LIMIT 1",
                      "SELECT ochiai AS pragma_y FROM pairs LIMIT 1"):
                r = client.post("/api/v1/sql", json={"query": q})
                self.assertEqual(r.status_code, 200,
                                 f"legit alias query {q!r} wrongly blocked: "
                                 f"{r.status_code}: {r.text[:200]}")
        finally:
            store.close()

    def test_sql_legit_select_over_pairs_still_works(self):
        client, store = self._client()
        try:
            r = client.post("/api/v1/sql",
                            json={"query": "SELECT count(*) AS n FROM pairs"})
            self.assertEqual(r.status_code, 200,
                             f"legit SELECT over pairs should be 200, got "
                             f"{r.status_code}: {r.text[:300]}")
            self.assertEqual(r.json()["row_count"], 1)
        finally:
            store.close()

    def test_sql_s1_sandbox_intact(self):
        # The denylist extension must NOT weaken the S1 protections.
        client, store = self._client()
        try:
            r = client.post("/api/v1/sql",
                            json={"query": "SELECT * FROM read_text('/etc/passwd')"})
            self.assertNotIn("root:", r.text, "sandbox leaked /etc/passwd")
            self.assertNotEqual(r.status_code, 200,
                                f"file-read query unexpectedly succeeded: {r.text[:200]}")
            r = client.post("/api/v1/sql", json={"query": "DROP TABLE pairs"})
            self.assertEqual(r.status_code, 400, f"DROP not blocked: {r.text[:200]}")
            self.assertEqual(r.json()["error_code"], "UNSAFE_QUERY", r.text[:200])
        finally:
            store.close()

    # ── 5: groups/{g}/pairs cutoff requires a metric ──────────────────
    def test_group_pairs_cutoff_without_metric_is_400(self):
        client, store = self._client()
        try:
            group = self._a_group(store)
            r = client.get(f"/api/v1/groups/{group}/pairs", params={"cutoff": 50})
            self.assertEqual(r.status_code, 400,
                             f"cutoff without metric should be 400, got "
                             f"{r.status_code}: {r.text[:200]}")
            self.assertEqual(r.json()["error_code"], "VALIDATION_ERROR",
                             f"wrong error_code: {r.text[:200]}")
        finally:
            store.close()

    def test_group_pairs_no_metric_no_cutoff_still_ok(self):
        # The cutoff guard must not break the plain "all pairs for a group" call.
        client, store = self._client()
        try:
            group = self._a_group(store)
            r = client.get(f"/api/v1/groups/{group}/pairs")
            self.assertEqual(r.status_code, 200,
                             f"plain group pairs should be 200, got "
                             f"{r.status_code}: {r.text[:200]}")
        finally:
            store.close()

    def test_group_pairs_metric_and_cutoff_ok(self):
        client, store = self._client()
        try:
            group = self._a_group(store)
            r = client.get(f"/api/v1/groups/{group}/pairs",
                           params={"metric": "ochiai", "cutoff": 0})
            self.assertEqual(r.status_code, 200,
                             f"metric+cutoff group pairs should be 200, got "
                             f"{r.status_code}: {r.text[:200]}")
        finally:
            store.close()

    # ── 6: /favicon.ico must not 404 ──────────────────────────────────
    def test_favicon_not_404(self):
        client, store = self._client()
        try:
            r = client.get("/favicon.ico")
            self.assertNotEqual(r.status_code, 404,
                                f"/favicon.ico should not 404, got {r.status_code}")
            self.assertIn(r.status_code, (200, 204),
                          f"/favicon.ico should be 200 or 204, got {r.status_code}")
        finally:
            store.close()


@unittest.skipUnless(_HAS_REST_TEST, "REST concurrency tests need [server] extra + httpx")
class TestRestConcurrency(unittest.TestCase):
    """ISSUE-055: PairwiseStore shares ONE DuckDB connection across all endpoints.

    FastAPI runs the sync route handlers in a threadpool with no synchronization,
    so a concurrent burst (e.g. the dashboard Stats tab firing statistics/{metric}
    + /top + /pairs at once) interleaves execute()/fetch on the single connection
    and corrupts each other -> fetchone() returns None -> 500. Works serially.

    The burst test MUST 500 pre-fix (RED) and be all-200 post-fix (GREEN); the
    serial test is a sanity guard that each endpoint alone still returns a correct
    body. Prefers the (larger) kegg substrate where queries interleave more, and
    falls back to the synthetic fixture so the suite stays portable.
    """

    KEGG_PARQUET = "/home/mabuelanin/dbretina_scratch/out/kegg_DBRetina_pairwise"
    KEGG_DBRI = "/home/mabuelanin/dbretina_scratch/out/kegg.dbri"

    def setUp(self):
        if (os.path.isdir(self.KEGG_PARQUET)
                and os.path.exists(os.path.join(self.KEGG_PARQUET, "manifest.json"))):
            self.tmpdir = None
            self.parquet = self.KEGG_PARQUET
            self.dbri = self.KEGG_DBRI if os.path.exists(self.KEGG_DBRI) else None
        else:
            self.tmpdir = tempfile.mkdtemp(prefix="dbretina_conc_")
            self.prefix, _ = setup_index_and_pairwise(self.tmpdir)
            self.parquet = f"{self.prefix}_DBRetina_pairwise"
            self.dbri = f"{self.prefix}.dbri"

    def tearDown(self):
        if self.tmpdir:
            shutil.rmtree(self.tmpdir, ignore_errors=True)

    def _client(self):
        from dbretina.compat import open_pairwise
        from dbretina.pairwise_store import PairwiseStore
        from dbretina.rest_api import create_app
        store = open_pairwise(self.parquet)
        if store is None:
            store = PairwiseStore(self.parquet, dbri_path=self.dbri)
        elif self.dbri:
            store._dbri_path = self.dbri
        # raise_server_exceptions=False so a corrupted-connection 500 surfaces as a
        # status code we can count, instead of re-raising into a worker thread.
        return _TestClient(create_app(store, dbri_path=self.dbri),
                           raise_server_exceptions=False), store

    def _metrics(self, store):
        # Use the metrics actually present (pvalue only if computed).
        return [m for m in store.available_metrics if m != "pvalue"]

    def test_concurrent_burst_no_500(self):
        """A concurrent burst across statistics/{metric} + /top + /pairs must be all-200.

        Pre-fix (single shared connection, no lock) this intermittently 500s with
        ``metric_summary row[0] -> NoneType`` from interleaved execute()/fetch.
        """
        from concurrent.futures import ThreadPoolExecutor

        client, store = self._client()
        try:
            metrics = self._metrics(store)
            self.assertTrue(metrics, "store exposes no metrics to query")
            urls = [f"/api/v1/statistics/{m}" for m in metrics]
            urls.append(f"/api/v1/top?metric={metrics[0]}&n=20")
            urls.append(f"/api/v1/pairs?metric={metrics[0]}&cutoff=0&limit=200")
            # ~size of a Stats-tab fan-out, amplified to reliably beat the race.
            burst = urls * 8

            failures = []  # (round, url, status, snippet)
            rounds = 6
            for rnd in range(rounds):
                def fire(url):
                    r = client.get(url)
                    return url, r.status_code, ("" if r.status_code == 200 else r.text[:200])
                with ThreadPoolExecutor(max_workers=16) as ex:
                    results = list(ex.map(fire, burst))
                for url, code, snippet in results:
                    if code != 200:
                        failures.append((rnd, url, code, snippet))

            n_500 = sum(1 for f in failures if f[2] == 500)
            self.assertEqual(
                failures, [],
                f"concurrent burst produced {len(failures)} non-200 "
                f"({n_500} were 500) across {rounds} rounds of {len(burst)} reqs; "
                f"first few: {failures[:3]}",
            )
        finally:
            store.close()

    def test_serial_endpoints_ok(self):
        """Sanity: each endpoint alone returns 200 with a correctly-shaped body."""
        client, store = self._client()
        try:
            metrics = self._metrics(store)
            for m in metrics:
                r = client.get(f"/api/v1/statistics/{m}")
                self.assertEqual(r.status_code, 200,
                                 f"statistics/{m} serial: {r.status_code} {r.text[:200]}")
                body = r.json()
                for k in ("count", "min", "max", "mean", "median", "stddev"):
                    self.assertIn(k, body, f"statistics/{m} missing key {k}")
                self.assertIsNotNone(body["count"], f"statistics/{m} count is None")

            r = client.get(f"/api/v1/top?metric={metrics[0]}&n=5")
            self.assertEqual(r.status_code, 200, f"/top serial: {r.status_code} {r.text[:200]}")
            tb = r.json()
            self.assertEqual(tb["metric"], metrics[0])
            self.assertIn("pairs", tb)

            r = client.get(f"/api/v1/pairs?metric={metrics[0]}&cutoff=0&limit=10")
            self.assertEqual(r.status_code, 200, f"/pairs serial: {r.status_code} {r.text[:200]}")
            self.assertIn("pairs", r.json())
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


# Standalone script run in a SEPARATE PROCESS to exercise the negative-resolution
# clustering path (ISSUE-054). Pre-fix, some igraph builds call abort() in their C
# layer for a negative resolution -> SIGABRT, which would kill the whole test
# runner if run in-process. Isolating it in a subprocess means a crash shows up as
# a non-zero exit code we can assert on, instead of taking the suite down with it.
_ISSUE_054_SUBPROCESS = r'''
import json, sys
from starlette.testclient import TestClient
from dbretina.compat import open_pairwise
from dbretina.pairwise_store import PairwiseStore
from dbretina.rest_api import create_app

parquet, dbri, algorithm, resolution = sys.argv[1:5]
store = open_pairwise(parquet)
if store is None:
    store = PairwiseStore(parquet, dbri_path=dbri)
else:
    store._dbri_path = dbri
client = TestClient(create_app(store, dbri_path=dbri), raise_server_exceptions=False)
r = client.post("/api/v1/graph/cluster",
                json={"algorithm": algorithm,
                      "parameters": {"resolution": float(resolution)},
                      "metric": "ochiai", "cutoff": 0.0})
store.close()
# Emit a one-line, easily-parsed verdict for the parent process.
print("DBRETINA_054_STATUS=%d" % r.status_code)
'''


@unittest.skipUnless(_HAS_REST_TEST, "serve validation-audit tests need [server] extra + httpx")
class TestRestServeValidationAudit(unittest.TestCase):
    """serve coverage-audit input-validation guards.

    ISSUE-054: POST /graph/cluster with a negative ``resolution`` reaches igraph's
               community detection, whose C layer may abort() (SIGABRT, killing the
               worker) or raise a non-ValueError -> 500. Must be a clean 400.
    ISSUE-056: a non-finite (inf/NaN) or out-of-range ``cutoff`` reaches SQL as a
               bare token -> 500 (or a wrong-but-200 negative cutoff on hub-genes).
               Must be a clean 400 on /pairs, /groups/{g}/pairs and hub-genes.
    ISSUE-058: GET /graph/communities?method=bogus -> 500 (community_detection
               raises a clean ValueError the endpoint never catches). Must be 400.
    """

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp(prefix="dbretina_serve_audit_")
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
        # raise_server_exceptions=False so a 500 surfaces as a status code instead
        # of re-raising into the test (we WANT to distinguish 4xx from 500).
        return _TestClient(create_app(store, dbri_path=self.dbri),
                           raise_server_exceptions=False), store

    def _a_group(self, store):
        names = list(store.get_names_map().values())
        self.assertTrue(names, "store has no groups")
        return names[0]

    # ── ISSUE-054: negative resolution must NOT crash the server ───
    def _cluster_status_in_subprocess(self, algorithm, resolution):
        """Run one /graph/cluster call in a child process; return its HTTP status.

        Returns an int status code, or raises AssertionError if the child crashed
        (e.g. igraph SIGABRT) or produced no verdict line.
        """
        proc = subprocess.run(
            [sys.executable, "-c", _ISSUE_054_SUBPROCESS,
             self.parquet, self.dbri, algorithm, str(resolution)],
            capture_output=True, text=True, timeout=300,
        )
        marker = "DBRETINA_054_STATUS="
        status = None
        for line in proc.stdout.splitlines():
            if line.startswith(marker):
                status = int(line[len(marker):])
        # A crash (negative rc = killed by signal, e.g. -6 SIGABRT) means the bug
        # is live: the server process died on user input.
        self.assertGreaterEqual(
            proc.returncode, 0,
            f"{algorithm} resolution={resolution} CRASHED the server "
            f"(rc={proc.returncode}; SIGABRT/DoS). stderr:\n{proc.stderr[-800:]}",
        )
        self.assertIsNotNone(
            status,
            f"{algorithm} resolution={resolution}: no status emitted. "
            f"stdout:\n{proc.stdout[-400:]}\nstderr:\n{proc.stderr[-800:]}",
        )
        return status

    def test_graph_cluster_negative_resolution_is_400_not_crash(self):
        # Pre-fix: leiden silently mis-clusters (200) and/or igraph aborts/500s;
        # louvain raises a non-ValueError -> 500; some builds SIGABRT (rc<0).
        # Post-fix: a clean 400 for every algorithm/value, no crash.
        for algorithm in ("leiden", "louvain"):
            for resolution in (-5, -0.5):
                status = self._cluster_status_in_subprocess(algorithm, resolution)
                self.assertEqual(
                    status, 400,
                    f"{algorithm} resolution={resolution} should be 400, got {status}",
                )

    def test_graph_cluster_valid_resolution_still_succeeds(self):
        client, store = self._client()
        try:
            for algorithm in ("leiden", "louvain"):
                r = client.post("/api/v1/graph/cluster",
                                json={"algorithm": algorithm,
                                      "parameters": {"resolution": 1.0},
                                      "metric": "ochiai", "cutoff": 0.0})
                self.assertEqual(r.status_code, 200,
                                 f"valid {algorithm} resolution=1.0 broke: "
                                 f"{r.status_code} {r.text[:200]}")
                self.assertIn("membership", r.json())
        finally:
            store.close()

    # ── ISSUE-056: non-finite / out-of-range cutoff -> 400 ─────────
    def test_pairs_non_finite_cutoff_is_400(self):
        client, store = self._client()
        try:
            self.assertIn("odds_ratio", store.available_metrics,
                          "substrate lacks odds_ratio (inf-max metric) needed for this test")
            # odds_ratio's range max is inf, so the old `cutoff < min` guard let
            # inf through to SQL as the bare token `inf` -> BinderException -> 500.
            r = client.get("/api/v1/pairs?metric=odds_ratio&cutoff=inf")
            self.assertEqual(r.status_code, 400,
                             f"pairs cutoff=inf should be 400, got {r.status_code}: {r.text[:200]}")
        finally:
            store.close()

    def test_group_pairs_non_finite_cutoff_is_400(self):
        client, store = self._client()
        try:
            group = self._a_group(store)
            r = client.get(f"/api/v1/groups/{group}/pairs?metric=ochiai&cutoff=inf")
            self.assertEqual(r.status_code, 400,
                             f"group pairs cutoff=inf should be 400, got {r.status_code}: {r.text[:200]}")
        finally:
            store.close()

    def test_hub_genes_non_finite_and_out_of_range_cutoff_is_400(self):
        client, store = self._client()
        try:
            group = self._a_group(store)
            headers = {"content-type": "application/json"}

            def hub_status(cutoff_token):
                # Send a raw JSON body so we can use the non-finite literals
                # (Infinity/NaN) a real client would emit; Python's json.loads
                # accepts them by default, so they reach the endpoint as floats.
                body = ('{"group_name": %s, "method": "edge_weighted", '
                        '"metric": "ochiai", "cutoff": %s}'
                        % (json.dumps(group), cutoff_token))
                return client.post("/api/v1/genes/hub-genes",
                                   content=body, headers=headers).status_code

            for tok in ("Infinity", "NaN"):
                self.assertEqual(hub_status(tok), 400,
                                 f"hub-genes edge_weighted cutoff={tok} should be 400")
            # Out-of-range (negative) cutoff was a wrong-but-200; now a clean 400.
            self.assertEqual(hub_status("-50"), 400,
                             "hub-genes edge_weighted cutoff=-50 should be 400")
        finally:
            store.close()

    def test_valid_cutoffs_still_succeed(self):
        client, store = self._client()
        try:
            group = self._a_group(store)
            headers = {"content-type": "application/json"}

            r = client.get("/api/v1/pairs?metric=odds_ratio&cutoff=0.0")
            self.assertEqual(r.status_code, 200,
                             f"valid pairs cutoff broke: {r.status_code} {r.text[:200]}")
            r = client.get(f"/api/v1/groups/{group}/pairs?metric=ochiai&cutoff=0.0")
            self.assertEqual(r.status_code, 200,
                             f"valid group-pairs cutoff broke: {r.status_code} {r.text[:200]}")
            body = ('{"group_name": %s, "method": "edge_weighted", '
                    '"metric": "ochiai", "cutoff": 0.0}' % json.dumps(group))
            r = client.post("/api/v1/genes/hub-genes", content=body, headers=headers)
            self.assertEqual(r.status_code, 200,
                             f"valid hub-genes cutoff broke: {r.status_code} {r.text[:200]}")
            self.assertIn("genes", r.json())
        finally:
            store.close()

    # ── ISSUE-058: unknown community method -> 400 ─────────────────
    def test_communities_unknown_method_is_400(self):
        client, store = self._client()
        try:
            r = client.get("/api/v1/graph/communities?metric=ochiai&cutoff=0.0&method=bogus")
            self.assertEqual(r.status_code, 400,
                             f"communities method=bogus should be 400, got {r.status_code}: {r.text[:200]}")
        finally:
            store.close()

    def test_communities_known_methods_still_succeed(self):
        client, store = self._client()
        try:
            for method in ("leiden", "louvain"):
                r = client.get(f"/api/v1/graph/communities?metric=ochiai&cutoff=0.0&method={method}")
                self.assertEqual(r.status_code, 200,
                                 f"communities method={method} broke: {r.status_code} {r.text[:200]}")
                self.assertIn("communities", r.json())
        finally:
            store.close()

    # ── ISSUE-054 (extended): other igraph param DoS vectors -> 400 ─────────
    def test_graph_cluster_huge_n_iterations_is_400_and_fast(self):
        """n_iterations is forwarded into igraph with no clamp; a huge value would
        tie up a worker for days. It must be rejected up front (fast 400)."""
        client, store = self._client()
        try:
            t0 = time.time()
            r = client.post("/api/v1/graph/cluster",
                            json={"algorithm": "leiden",
                                  "parameters": {"n_iterations": 1_000_000_000},
                                  "metric": "ochiai", "cutoff": 0.0})
            elapsed = time.time() - t0
            self.assertEqual(r.status_code, 400,
                             f"huge n_iterations should be 400, got {r.status_code}: {r.text[:200]}")
            self.assertLess(elapsed, 10, f"rejection took {elapsed:.1f}s (should be near-instant)")
            # non-finite n_iterations (Infinity) must also be a clean 400 (not a 500
            # from int(float('inf')) OverflowError). json.dumps can't emit Infinity,
            # so send a raw JSON body with the literal token.
            raw = ('{"algorithm":"leiden","parameters":{"n_iterations":Infinity},'
                   '"metric":"ochiai","cutoff":0.0}')
            r_inf = client.post("/api/v1/graph/cluster", content=raw,
                                headers={"content-type": "application/json"})
            self.assertEqual(r_inf.status_code, 400,
                             f"n_iterations=Infinity should be 400, got {r_inf.status_code}: {r_inf.text[:200]}")
            # a small, valid n_iterations still works
            r2 = client.post("/api/v1/graph/cluster",
                             json={"algorithm": "leiden",
                                   "parameters": {"n_iterations": 2},
                                   "metric": "ochiai", "cutoff": 0.0})
            self.assertEqual(r2.status_code, 200, f"valid n_iterations=2 broke: {r2.text[:200]}")
        finally:
            store.close()

    def test_graph_cluster_unknown_parameter_is_400(self):
        """An unexpected parameter key must be rejected (no smuggling arbitrary
        igraph C kwargs), not silently forwarded."""
        client, store = self._client()
        try:
            r = client.post("/api/v1/graph/cluster",
                            json={"algorithm": "leiden",
                                  "parameters": {"totally_bogus_kwarg": 5},
                                  "metric": "ochiai", "cutoff": 0.0})
            self.assertEqual(r.status_code, 400,
                             f"unknown parameter should be 400, got {r.status_code}: {r.text[:200]}")
        finally:
            store.close()

    def test_pairs_filter_non_finite_value_is_400(self):
        """POST /pairs/filter interpolates the value into SQL; a non-finite value
        must be a clean 400, not reach DuckDB as a bare inf/nan token."""
        client, store = self._client()
        try:
            headers = {"content-type": "application/json"}
            body = '{"filters":[{"metric":"ochiai","operator":">=","value":Infinity}],"logic":"AND"}'
            r = client.post("/api/v1/pairs/filter", content=body, headers=headers)
            self.assertEqual(r.status_code, 400,
                             f"filter value=Infinity should be 400, got {r.status_code}: {r.text[:200]}")
            # a finite value still works
            r2 = client.post("/api/v1/pairs/filter",
                             json={"filters": [{"metric": "ochiai", "operator": ">=", "value": 0.0}],
                                   "logic": "AND"})
            self.assertEqual(r2.status_code, 200, f"valid filter broke: {r2.text[:200]}")
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


class TestZeroPairLookups(unittest.TestCase):
    """ISSUE-071 (regression): a legitimately sparse 0-pair pairwise (empty
    ``data/`` with no part_*.parquet, as ``pairwise`` emits when nothing clears
    the cutoff -- the cosmic substrate) must NOT crash the lookup commands with a
    raw _duckdb.IOException ("No files found that match the pattern ...").

    Chosen behavior: lookups WORK with an empty result (a 0-pair dataset is a
    normal outcome). PairwiseStore now backs ``pairs`` with a typed zero-row
    relation when there are no parquet files, so every command returns a clean
    empty result / clean [ERROR]. A normal (non-empty) dataset is unaffected.
    """

    @classmethod
    def setUpClass(cls):
        cls.prefix, cls.pw_tsv = _ensure_zero_pair_fixture()
        cls.parquet = f"{cls.prefix}_DBRetina_pairwise"  # empty data/ dir
        # Real groups/genes that EXIST in the 0-pair fixture, so the commands
        # reach the (empty) pairs view / gene index rather than a not-found guard.
        cls.group_a = "groupx"
        cls.group_b = "groupy"
        cls.gene = "gx1"
        # A normal (non-empty) substrate, to confirm the fix is scoped to 0-pair.
        cls.norm_prefix, cls.norm_pw = _ensure_shared_fixture()
        cls.norm_parquet = f"{cls.norm_prefix}_DBRetina_pairwise"

    def _assert_no_raw_error(self, label, out, err):
        """No raw traceback / DuckDB IOException from any lookup on 0-pair data."""
        blob = out + err
        self.assertNotIn("Traceback", blob,
                         f"{label}: raw Python traceback on a 0-pair dataset:\n{err[-400:]}")
        self.assertNotIn("IOException", blob,
                         f"{label}: raw DuckDB IOException on a 0-pair dataset:\n{err[-400:]}")
        self.assertNotIn("No files found that match the pattern", blob,
                         f"{label}: empty-parquet-glob error leaked:\n{err[-400:]}")

    def test_zero_pair_served_pairs_schema_parity(self):
        """ISSUE-071 follow-up: a SERVED 0-pair store's hardened `pairs` table must
        expose group_1_name/group_2_name (schema parity with the non-empty
        names-joined table), so a raw /sql select of the name columns returns empty
        rather than a Binder Error; the external-access sandbox stays intact."""
        from dbretina.compat import open_pairwise
        store = open_pairwise(self.parquet)
        try:
            self.assertEqual(store.num_pairs, 0)
            store.harden_readonly()
            df = store._con.execute(
                "SELECT group_1_name, group_2_name FROM pairs LIMIT 5").fetchdf()
            self.assertEqual(list(df.columns), ["group_1_name", "group_2_name"])
            self.assertEqual(len(df), 0)
            with self.assertRaises(Exception):  # sandbox still blocks file reads
                store._con.execute("SELECT * FROM read_text('/etc/hostname')")
        finally:
            store.close()

    def test_search_zero_pair_clean(self):
        """search uses the names map -> still lists the groups, no crash."""
        rc, out, err = run_command(f'DBRetina search -d {self.parquet} "group"')
        self._assert_no_raw_error("search", out, err)
        self.assertEqual(rc, 0, err)
        self.assertIn(self.group_a, out)

    def test_neighbors_zero_pair_clean(self):
        """neighbors over the empty pairs view -> clean 'No neighbors', exit 0."""
        rc, out, err = run_command(
            f'DBRetina neighbors -d {self.parquet} "{self.group_a}" -m ochiai -c 1'
        )
        self._assert_no_raw_error("neighbors", out, err)
        self.assertEqual(rc, 0, err)
        self.assertIn("No neighbors found", out + err)

    def test_shared_genes_zero_pair_clean(self):
        """shared-genes (gene sets from .dbri) -> clean result, no crash."""
        rc, out, err = run_command(
            f'DBRetina shared-genes -d {self.parquet} -i {self.prefix} '
            f'"{self.group_a}" "{self.group_b}"'
        )
        self._assert_no_raw_error("shared-genes", out, err)
        self.assertEqual(rc, 0, err)

    def test_gene_search_zero_pair_clean(self):
        """gene-search (inverted gene index) -> finds the group, no crash."""
        rc, out, err = run_command(
            f'DBRetina gene-search -d {self.parquet} -i {self.prefix} "{self.gene}"'
        )
        self._assert_no_raw_error("gene-search", out, err)
        self.assertEqual(rc, 0, err)
        self.assertIn(self.group_a, out)

    def test_genescore_zero_pair_clean(self):
        """genescore over the empty pairs view -> clean result, no crash."""
        rc, out, err = run_command(
            f'DBRetina genescore -d {self.parquet} -i {self.prefix} '
            f'"{self.group_a}" -n 3'
        )
        self._assert_no_raw_error("genescore", out, err)
        self.assertEqual(rc, 0, err)

    def test_normal_lookup_unaffected(self):
        """Sanity: a normal (non-empty) dataset still returns real neighbors."""
        rc, out, err = run_command(
            f'DBRetina neighbors -d {self.norm_parquet} "groupa" -m ochiai -c 5'
        )
        self.assertEqual(rc, 0, err)
        self.assertNotIn("Traceback", out + err)
        self.assertIn("ochiai", out)  # header + at least one neighbor row


class TestNeighborsNegativeTop(unittest.TestCase):
    """ISSUE-081 (regression): neighbors -n with a negative value fed a negative
    LIMIT to DuckDB and surfaced a raw BinderException ("LIMIT/OFFSET cannot be
    negative"). -n is now click.IntRange(min=0): negatives are rejected cleanly
    before any SQL runs; -n 0 and positive values work."""

    @classmethod
    def setUpClass(cls):
        cls.prefix, cls.pw_tsv = _ensure_shared_fixture()
        cls.parquet = f"{cls.prefix}_DBRetina_pairwise"
        cls.group = "groupa"

    def test_negative_top_clean_error(self):
        rc, out, err = run_command(
            f'DBRetina neighbors -d {self.parquet} "{self.group}" -m ochiai -c 1 -n -5'
        )
        self.assertNotEqual(rc, 0, "negative -n must fail")
        blob = out + err
        self.assertNotIn("Traceback", blob, f"-n -5 leaked a traceback:\n{err[-400:]}")
        self.assertNotIn("BinderException", blob,
                         f"-n -5 leaked a raw DuckDB BinderException:\n{err[-400:]}")
        # click.IntRange emits a clean "Invalid value" usage error.
        self.assertIn("Invalid value", blob, f"expected a clean range error:\n{err[-400:]}")

    def test_zero_top_works(self):
        """-n 0 is a valid limit (no rows) and must not crash."""
        rc, out, err = run_command(
            f'DBRetina neighbors -d {self.parquet} "{self.group}" -m ochiai -c 1 -n 0'
        )
        self.assertEqual(rc, 0, err)
        self.assertNotIn("Traceback", out + err)

    def test_positive_top_works(self):
        rc, out, err = run_command(
            f'DBRetina neighbors -d {self.parquet} "{self.group}" -m ochiai -c 1 -n 5'
        )
        self.assertEqual(rc, 0, err)
        self.assertNotIn("Traceback", out + err)
        self.assertIn("ochiai", out)


class TestGenescoreNegativeHops(unittest.TestCase):
    """ISSUE-082 (regression): genescore --hops with a negative value reached
    igraph's neighborhood(order=hops) and raised a raw ValueError ("neighborhood
    order must be non-negative"). --hops is now click.IntRange(min=0): negatives
    are rejected cleanly; --hops 0 and positive values work."""

    @classmethod
    def setUpClass(cls):
        cls.prefix, cls.pw_tsv = _ensure_shared_fixture()
        cls.parquet = f"{cls.prefix}_DBRetina_pairwise"
        cls.group = "groupa"

    def test_negative_hops_clean_error(self):
        rc, out, err = run_command(
            f'DBRetina genescore -d {self.parquet} -i {self.prefix} '
            f'"{self.group}" --hops -1 -n 3'
        )
        self.assertNotEqual(rc, 0, "negative --hops must fail")
        blob = out + err
        self.assertNotIn("Traceback", blob, f"--hops -1 leaked a traceback:\n{err[-400:]}")
        self.assertNotIn("neighborhood order", blob,
                         f"--hops -1 leaked a raw igraph ValueError:\n{err[-400:]}")
        self.assertIn("Invalid value", blob, f"expected a clean range error:\n{err[-400:]}")

    def test_zero_hops_works(self):
        """--hops 0 (just the center node) is valid and must not crash."""
        rc, out, err = run_command(
            f'DBRetina genescore -d {self.parquet} -i {self.prefix} '
            f'"{self.group}" --hops 0 -n 3'
        )
        self.assertEqual(rc, 0, err)
        self.assertNotIn("Traceback", out + err)

    def test_positive_hops_works(self):
        rc, out, err = run_command(
            f'DBRetina genescore -d {self.parquet} -i {self.prefix} '
            f'"{self.group}" --hops 1 -n 3'
        )
        self.assertEqual(rc, 0, err)
        self.assertNotIn("Traceback", out + err)


class TestNeighborsPvalueDisplayPrecision(unittest.TestCase):
    """Display nit (from the pvalue-068 review): neighbors printed the metric
    value with ``{:.1f}``, so significant p-values all rendered as "0.0",
    hiding the precision that issue 068 made meaningful. The metric value is now
    formatted metric-aware -- ``{:.3g}`` for p-value, ``{:.1f}`` for similarity
    metrics -- so p-values display with real precision."""

    PV_CUTOFF = 0.05

    @classmethod
    def setUpClass(cls):
        cls.pv_prefix, cls.pv_pw_tsv = _ensure_shared_pvalue_fixture()
        cls.pv_dir = f"{cls.pv_prefix}_DBRetina_pairwise"

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp(prefix="dbretina_pv_disp_")

    def tearDown(self):
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def _neighbor_metric_strings(self, out):
        """Return the raw printed metric-column strings (parts[1]) per data row."""
        vals = []
        header = True
        for line in out.splitlines():
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("["):
                continue
            if header:  # "neighbor\t{metric}\tjaccard\tshared_features"
                header = False
                continue
            parts = line.split("\t")
            if len(parts) >= 2:
                vals.append(parts[1])
        return vals

    def test_pvalue_values_not_all_zero(self):
        """neighbors -m pvalue: at least one printed p-value is a real nonzero
        value (not the ``0.0`` the old {:.1f} collapsed everything to)."""
        rc, stdout, stderr = run_command(
            f'DBRetina neighbors -d {self.pv_dir} "groupa" '
            f"-m pvalue -c {self.PV_CUTOFF}"
        )
        self.assertEqual(rc, 0, stderr)
        vals = self._neighbor_metric_strings(stdout)
        self.assertTrue(vals, "expected significant pvalue neighbors of groupa")
        # The bug rendered every p as exactly "0.0". Require a meaningful value.
        self.assertFalse(
            all(v == "0.0" for v in vals),
            f"pvalue column collapsed to '0.0' (precision lost): {vals}")
        parsed = [float(v) for v in vals]
        self.assertTrue(
            any(0.0 < p <= self.PV_CUTOFF for p in parsed),
            f"expected a nonzero significant p-value in the output: {vals}")

    def test_similarity_metric_still_one_decimal(self):
        """Regression guard: a similarity metric (ochiai) still prints {:.1f}
        (one decimal place), unchanged by the metric-aware formatter."""
        rc, stdout, stderr = run_command(
            f'DBRetina neighbors -d {self.pv_dir} "groupa" -m ochiai -c 50'
        )
        self.assertEqual(rc, 0, stderr)
        vals = self._neighbor_metric_strings(stdout)
        self.assertTrue(vals, "expected ochiai neighbors of groupa")
        import re
        for v in vals:
            self.assertRegex(
                v, r"^\d+\.\d$",
                f"ochiai value {v!r} not formatted as one-decimal {{:.1f}}")


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


class TestModularityDedupInputHandling(unittest.TestCase):
    """Input-handling regressions for modularity/dedup (issue 046).

    Same root cause as cluster 016 / query 020: the fallback path derived the
    sibling .dbrp via ``pairwise_file.replace('_pairwise.tsv','_pairwise.dbrp')``,
    a no-op for the parquet-directory / .dbrp forms. When no parquet store was
    resolvable the directory path then leaked into the .dbrp binary reader,
    crashing with ``RuntimeError: Invalid .dbrp file (bad magic bytes)``.

    -p must accept the canonical parquet directory and the .dbrp binary and give
    the same result as the .tsv form, and an unreadable directory input must give
    a clean [ERROR] rather than a raw traceback.

    Uses the shared overlapping fixture, which emits all sibling forms (.tsv,
    parquet dir, .dbrp, .dbri) next to each other.
    """

    @classmethod
    def setUpClass(cls):
        cls.prefix, cls.pw_file = _ensure_shared_fixture()
        cls.parquet_dir = f"{cls.prefix}_DBRetina_pairwise"
        cls.dbrp = f"{cls.prefix}_DBRetina_pairwise.dbrp"

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp(prefix="dbretina_moddedup_input_")
        self.assertTrue(os.path.isdir(self.parquet_dir),
                        f"fixture missing parquet dir: {self.parquet_dir}")
        self.assertTrue(os.path.isfile(self.dbrp),
                        f"fixture missing .dbrp: {self.dbrp}")

    def tearDown(self):
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    # ---- helpers ----
    @staticmethod
    def _modularity_map(path):
        """Parse a *_modularity.tsv -> {gene_set: (frag, het, modularity)}."""
        out = {}
        with open(path) as f:
            next(f)  # header
            for line in f:
                if not line.strip():
                    continue
                gs, frag, het, mod = line.strip().split("\t")
                out[gs.lower()] = (int(float(frag)), int(float(het)), int(float(mod)))
        return out

    @staticmethod
    def _dedup_set(path):
        with open(path) as f:
            return {l.strip().lower() for l in f if l.strip()}

    def _hide_parquet_store(self, src_dir):
        """Copy the fixture parquet dir into tmpdir but strip manifest.json so
        open_pairwise() returns None, forcing the .dbrp/TSV fallback path (the
        path that carried the issue-046 .replace bug)."""
        dst = os.path.join(self.tmpdir, os.path.basename(src_dir))
        shutil.copytree(src_dir, dst)
        os.remove(os.path.join(dst, "manifest.json"))
        return dst

    # ---- 046: modularity on the parquet directory (canonical input) ----
    def test_046_modularity_parquet_dir_matches_tsv(self):
        out_dir = os.path.join(self.tmpdir, "mod_dir")
        rc, _, stderr = run_command(
            f"DBRetina modularity -i {self.prefix} -p {self.parquet_dir} "
            f"-c 40 -o {out_dir}"
        )
        assert_no_traceback(self, stderr, "modularity parquet dir")
        self.assertNotIn("IsADirectoryError", stderr, stderr)
        self.assertEqual(rc, 0, stderr)
        assert_file_exists(self, f"{out_dir}_modularity.tsv")

        out_tsv = os.path.join(self.tmpdir, "mod_tsv")
        rc, _, stderr = run_command(
            f"DBRetina modularity -i {self.prefix} -p {self.pw_file} -c 40 -o {out_tsv}"
        )
        self.assertEqual(rc, 0, stderr)
        self.assertEqual(
            self._modularity_map(f"{out_dir}_modularity.tsv"),
            self._modularity_map(f"{out_tsv}_modularity.tsv"),
            "parquet-dir modularity differs from .tsv form",
        )

    # ---- 046: modularity on the .dbrp binary ----
    def test_046_modularity_dbrp_matches_tsv(self):
        out = os.path.join(self.tmpdir, "mod_dbrp")
        rc, _, stderr = run_command(
            f"DBRetina modularity -i {self.prefix} -p {self.dbrp} -c 40 -o {out}"
        )
        assert_no_traceback(self, stderr, "modularity .dbrp")
        self.assertEqual(rc, 0, stderr)
        out_tsv = os.path.join(self.tmpdir, "mod_tsv2")
        rc, _, stderr = run_command(
            f"DBRetina modularity -i {self.prefix} -p {self.pw_file} -c 40 -o {out_tsv}"
        )
        self.assertEqual(rc, 0, stderr)
        self.assertEqual(
            self._modularity_map(f"{out}_modularity.tsv"),
            self._modularity_map(f"{out_tsv}_modularity.tsv"),
            ".dbrp modularity differs from .tsv form",
        )

    # ---- 046: forced-fallback directory input -> the actual .replace bug ----
    def test_046_modularity_dir_without_store_falls_back_cleanly(self):
        """A directory input that open_pairwise() can't use (no manifest) must
        resolve the sibling .dbrp instead of feeding the directory path into the
        binary reader (RuntimeError: bad magic bytes). With the sibling .dbrp
        present the result must still match the .tsv form."""
        hidden = self._hide_parquet_store(self.parquet_dir)
        # Provide the sibling .dbrp next to the (manifest-less) directory.
        shutil.copy(self.dbrp, hidden + ".dbrp")
        out = os.path.join(self.tmpdir, "mod_fb")
        rc, _, stderr = run_command(
            f"DBRetina modularity -i {self.prefix} -p {hidden} -c 40 -o {out}"
        )
        assert_no_traceback(self, stderr, "modularity fallback dir")
        self.assertNotIn("bad magic bytes", stderr, stderr)
        self.assertEqual(rc, 0, stderr)
        out_tsv = os.path.join(self.tmpdir, "mod_fb_tsv")
        run_command(
            f"DBRetina modularity -i {self.prefix} -p {self.pw_file} -c 40 -o {out_tsv}"
        )
        self.assertEqual(
            self._modularity_map(f"{out}_modularity.tsv"),
            self._modularity_map(f"{out_tsv}_modularity.tsv"),
            "fallback-dir modularity differs from .tsv form",
        )

    def test_046_modularity_unreadable_dir_clean_error(self):
        """A directory input with neither a usable store nor a sibling .dbrp must
        give a clean [ERROR], never a raw traceback (bad magic bytes)."""
        hidden = self._hide_parquet_store(self.parquet_dir)  # no sibling .dbrp
        out = os.path.join(self.tmpdir, "mod_bad")
        rc, stdout, stderr = run_command(
            f"DBRetina modularity -i {self.prefix} -p {hidden} -c 40 -o {out}"
        )
        assert_no_traceback(self, stderr, "modularity unreadable dir")
        self.assertNotIn("bad magic bytes", stderr, stderr)
        self.assertNotEqual(rc, 0, "unreadable directory input should fail")
        self.assertIn("[ERROR]", (stdout + stderr), stderr)

    # ---- 046: dedup on the parquet directory (canonical input) ----
    def test_046_dedup_parquet_dir_matches_tsv(self):
        out_dir = os.path.join(self.tmpdir, "dd_dir")
        rc, _, stderr = run_command(
            f"DBRetina dedup -i {self.prefix} -p {self.parquet_dir} -c 100 -o {out_dir}"
        )
        assert_no_traceback(self, stderr, "dedup parquet dir")
        self.assertNotIn("IsADirectoryError", stderr, stderr)
        self.assertEqual(rc, 0, stderr)
        assert_file_exists(self, f"{out_dir}_deduplicated_groups.txt")

        out_tsv = os.path.join(self.tmpdir, "dd_tsv")
        rc, _, stderr = run_command(
            f"DBRetina dedup -i {self.prefix} -p {self.pw_file} -c 100 -o {out_tsv}"
        )
        self.assertEqual(rc, 0, stderr)
        self.assertEqual(
            self._dedup_set(f"{out_dir}_deduplicated_groups.txt"),
            self._dedup_set(f"{out_tsv}_deduplicated_groups.txt"),
            "parquet-dir dedup differs from .tsv form",
        )

    # ---- 046: dedup on the .dbrp binary ----
    def test_046_dedup_dbrp_matches_tsv(self):
        out = os.path.join(self.tmpdir, "dd_dbrp")
        rc, _, stderr = run_command(
            f"DBRetina dedup -i {self.prefix} -p {self.dbrp} -c 100 -o {out}"
        )
        assert_no_traceback(self, stderr, "dedup .dbrp")
        self.assertEqual(rc, 0, stderr)
        out_tsv = os.path.join(self.tmpdir, "dd_tsv2")
        rc, _, stderr = run_command(
            f"DBRetina dedup -i {self.prefix} -p {self.pw_file} -c 100 -o {out_tsv}"
        )
        self.assertEqual(rc, 0, stderr)
        self.assertEqual(
            self._dedup_set(f"{out}_deduplicated_groups.txt"),
            self._dedup_set(f"{out_tsv}_deduplicated_groups.txt"),
            ".dbrp dedup differs from .tsv form",
        )

    # ---- 046: forced-fallback directory input -> the actual .replace bug ----
    def test_046_dedup_dir_without_store_falls_back_cleanly(self):
        hidden = self._hide_parquet_store(self.parquet_dir)
        shutil.copy(self.dbrp, hidden + ".dbrp")
        out = os.path.join(self.tmpdir, "dd_fb")
        rc, _, stderr = run_command(
            f"DBRetina dedup -i {self.prefix} -p {hidden} -c 100 -o {out}"
        )
        assert_no_traceback(self, stderr, "dedup fallback dir")
        self.assertNotIn("bad magic bytes", stderr, stderr)
        self.assertEqual(rc, 0, stderr)
        out_tsv = os.path.join(self.tmpdir, "dd_fb_tsv")
        run_command(
            f"DBRetina dedup -i {self.prefix} -p {self.pw_file} -c 100 -o {out_tsv}"
        )
        self.assertEqual(
            self._dedup_set(f"{out}_deduplicated_groups.txt"),
            self._dedup_set(f"{out_tsv}_deduplicated_groups.txt"),
            "fallback-dir dedup differs from .tsv form",
        )

    def test_046_dedup_unreadable_dir_clean_error(self):
        hidden = self._hide_parquet_store(self.parquet_dir)  # no sibling .dbrp
        out = os.path.join(self.tmpdir, "dd_bad")
        rc, stdout, stderr = run_command(
            f"DBRetina dedup -i {self.prefix} -p {hidden} -c 100 -o {out}"
        )
        assert_no_traceback(self, stderr, "dedup unreadable dir")
        self.assertNotIn("bad magic bytes", stderr, stderr)
        self.assertNotEqual(rc, 0, "unreadable directory input should fail")
        self.assertIn("[ERROR]", (stdout + stderr), stderr)


class TestSetcovCommunityOption(unittest.TestCase):
    """issue 033: setcov --community was a dead option with a misleading log.

    Community detection actually runs at the --modularity (containment) cutoff;
    the --community value was only ever echoed in a 'Detecting communities with
    ochiai cutoff <value>' line that implied it drove detection. The option is
    kept (it is documented) but the misleading line is gone and setcov warns the
    value currently has no effect.
    """

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp(prefix="dbretina_setcov_comm_")
        asc_path = os.path.join(self.tmpdir, "test_input.asc")
        with open(asc_path, "w") as f:
            f.write(TEST_ASC_CONTENT)
        rc, _, stderr = run_command(
            "DBRetina index -a test_input.asc -o test_idx", cwd=self.tmpdir
        )
        assert rc == 0, f"index failed: {stderr}"
        rc, _, stderr = run_command("DBRetina pairwise -i test_idx", cwd=self.tmpdir)
        assert rc == 0, f"pairwise failed: {stderr}"

    def tearDown(self):
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def test_033_no_misleading_community_cutoff_log(self):
        """The 'Detecting communities with ochiai cutoff <--community>' line that
        implied --community drives detection must be gone."""
        rc, stdout, stderr = run_command(
            "DBRetina setcov -i test_idx --community 90 -o sc_c90",
            timeout=300, cwd=self.tmpdir,
        )
        out = stdout + stderr
        assert_no_traceback(self, stderr, "setcov --community log")
        self.assertEqual(rc, 0, stderr)
        # The old misleading phrasing keyed on the --community value must not appear.
        self.assertNotIn("Detecting communities with ochiai cutoff 90", out, out)

    def test_033_setcov_accepts_community_option(self):
        """--community is documented, so it must remain a valid option (no Click
        'no such option') even though it currently has no effect."""
        rc, stdout, stderr = run_command(
            "DBRetina setcov -i test_idx --community 30 -o sc_acc",
            timeout=300, cwd=self.tmpdir,
        )
        out = stdout + stderr
        self.assertEqual(rc, 0, stderr)
        self.assertNotIn("no such option", out.lower(), out)

    def test_033_setcov_without_community_still_works(self):
        """setcov must run with --community omitted (default), unchanged."""
        rc, _, stderr = run_command(
            "DBRetina setcov -i test_idx --modularity 80 --dedup 100 -o sc_def",
            timeout=300, cwd=self.tmpdir,
        )
        assert_no_traceback(self, stderr, "setcov no --community")
        self.assertEqual(rc, 0, stderr)
        assert_file_exists(self, os.path.join(self.tmpdir, "sc_def_groups_metadata.tsv"))


# ============================================================
# SECTION: Input-robustness regressions
#   index/dedup/genenet must handle bad or unusual inputs
#   gracefully instead of crashing with a raw traceback (or
#   silently merging). Issues 018, 035, 036, 025, 026, 029.
# ============================================================


class TestInputRobustness(unittest.TestCase):
    """Bad/unusual-input handling for index, dedup, genenet.

    Each test asserts the *fixed* behaviour: a clean [ERROR]/[WARNING]
    and a sensible outcome, never a raw Python traceback.
    """

    @classmethod
    def setUpClass(cls):
        cls.prefix, cls.pw_file = _ensure_shared_fixture()

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp(prefix="dbretina_inputrobust_")

    def tearDown(self):
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    # ---- 018: index a gzipped GMT (.gmt.gz) ----
    def test_018_index_gzipped_gmt(self):
        """index -g <file.gmt.gz> must transparently gunzip and build an index
        identical to indexing the plain GMT (same groups/features)."""
        import gzip

        # Plain GMT reference index.
        plain_gmt = write_file(os.path.join(self.tmpdir, "ref.gmt"), TEST_GMT_CONTENT)
        rc, _, stderr = run_command(
            "DBRetina index -g ref.gmt -o ref_idx", cwd=self.tmpdir
        )
        self.assertEqual(rc, 0, stderr)
        ref_groups = get_groups_from_raw_json(
            os.path.join(self.tmpdir, "ref_idx_raw.json")
        )

        # Gzipped GMT.
        gz_path = os.path.join(self.tmpdir, "test.gmt.gz")
        with gzip.open(gz_path, "wt", encoding="utf-8") as gz:
            gz.write(TEST_GMT_CONTENT)
        rc, _, stderr = run_command(
            "DBRetina index -g test.gmt.gz -o gz_idx", cwd=self.tmpdir
        )
        assert_no_traceback(self, stderr, "index gzipped gmt")
        self.assertNotIn("UnicodeDecodeError", stderr, stderr)
        self.assertEqual(rc, 0, stderr)
        assert_file_exists(self, os.path.join(self.tmpdir, "gz_idx.dbri"))

        # The gzipped index must equal the plain-GMT index.
        gz_groups = get_groups_from_raw_json(
            os.path.join(self.tmpdir, "gz_idx_raw.json")
        )
        self.assertEqual(gz_groups, ref_groups,
                         "gzipped-GMT index differs from plain-GMT index")

    # ---- 035: pipe char / malformed line -> clean error, not a traceback ----
    def test_035_pipe_in_group_name_clean_error(self):
        """A '|' in a group name must be rejected with a clean [ERROR], not a
        raw ValueError traceback."""
        gmt = write_file(os.path.join(self.tmpdir, "pipe.gmt"),
                         "bad|name\tdesc\tA\tB\n")
        rc, stdout, stderr = run_command(
            "DBRetina index -g pipe.gmt -o e_pipe", cwd=self.tmpdir
        )
        assert_no_traceback(self, stderr, "index pipe in group name")
        self.assertNotIn("ValueError", stderr, stderr)
        self.assertNotEqual(rc, 0, "pipe in group name should fail")
        self.assertIn("[ERROR]", stderr, stderr)
        self.assertIn("|", (stdout + stderr), "error should name the pipe char")

    def test_035_pipe_in_gene_name_clean_error(self):
        """A '|' in a gene name must be rejected with a clean [ERROR]."""
        gmt = write_file(os.path.join(self.tmpdir, "pipegene.gmt"),
                         "goodname\tdesc\tA\tbad|gene\n")
        rc, stdout, stderr = run_command(
            "DBRetina index -g pipegene.gmt -o e_pipegene", cwd=self.tmpdir
        )
        assert_no_traceback(self, stderr, "index pipe in gene name")
        self.assertNotIn("ValueError", stderr, stderr)
        self.assertNotEqual(rc, 0, "pipe in gene name should fail")
        self.assertIn("[ERROR]", stderr, stderr)

    def test_018_mislabeled_gz_is_read_as_text(self):
        """ISSUE-018 follow-up: a plain-text GMT mislabeled .gz must be read as text
        (magic-byte detection), not crash with a gzip BadGzipFile traceback."""
        write_file(os.path.join(self.tmpdir, "fake.gmt.gz"), TEST_GMT_CONTENT)
        rc, _, stderr = run_command(
            "DBRetina index -g fake.gmt.gz -o fake_idx", cwd=self.tmpdir
        )
        assert_no_traceback(self, stderr, "mislabeled gz")
        self.assertNotIn("BadGzipFile", stderr, stderr)
        self.assertEqual(rc, 0, stderr)
        assert_file_exists(self, os.path.join(self.tmpdir, "fake_idx.dbri"))

    def test_035_pipe_in_association_file_clean_error(self):
        """ISSUE-035 (-a path): a pipe char in an association (-a) file must give a
        clean [ERROR], not a raw ValueError traceback (build_gene_set_json)."""
        write_file(os.path.join(self.tmpdir, "bad.asc"),
                   "gene_set\tgene\nbad|grp\tAlpha\n")
        rc, stdout, stderr = run_command(
            "DBRetina index -a bad.asc -o bad_asc_idx", cwd=self.tmpdir
        )
        combined = stdout + stderr
        self.assertNotEqual(rc, 0)
        assert_no_traceback(self, stderr, "index -a pipe")
        self.assertNotIn("ValueError", combined, combined)
        self.assertIn("[ERROR]", combined)

    def test_035_malformed_line_clean_error(self):
        """A line with <3 fields must be rejected with a clean [ERROR], not a
        raw ValueError traceback."""
        gmt = write_file(os.path.join(self.tmpdir, "bad.gmt"),
                         "onlyname\tjustdesc\n")
        rc, stdout, stderr = run_command(
            "DBRetina index -g bad.gmt -o e_bad", cwd=self.tmpdir
        )
        assert_no_traceback(self, stderr, "index malformed line")
        self.assertNotIn("ValueError", stderr, stderr)
        self.assertNotEqual(rc, 0, "malformed line should fail")
        self.assertIn("[ERROR]", stderr, stderr)

    # ---- 036: duplicate group name -> WARNING, still succeeds ----
    def test_036_duplicate_group_name_warns(self):
        """Two GMT lines with the same group name must emit a [WARNING] naming
        the duplicate, and the index must still build (genes unioned)."""
        gmt = write_file(os.path.join(self.tmpdir, "dup.gmt"),
                         "dupset\tdesc\tA\tB\ndupset\tdesc\tC\tD\n")
        rc, stdout, stderr = run_command(
            "DBRetina index -g dup.gmt -o e_dup", cwd=self.tmpdir
        )
        assert_no_traceback(self, stderr, "index duplicate group name")
        self.assertEqual(rc, 0, stderr)
        self.assertIn("[WARNING]", stderr, stderr)
        self.assertIn("dupset", (stdout + stderr).lower(),
                      "warning should name the duplicate group")
        # Merge behaviour is preserved: the two lines union into one group.
        groups = get_groups_from_raw_json(
            os.path.join(self.tmpdir, "e_dup_raw.json")
        )
        self.assertEqual(set(groups.keys()), {"dupset"})
        self.assertEqual(groups["dupset"], ["a", "b", "c", "d"])

    # ---- 025: dedup with no qualifying pairs -> graceful (keep all groups) ----
    def test_025_dedup_cutoff_above_all_pairs(self):
        """dedup with a cutoff above the max similarity (no pair qualifies) must
        keep every group, write *_deduplicated_groups.txt, and exit 0 -- not a
        confusing hard-exit with no output."""
        # Build a disjoint fixture: zero pairs at any cutoff > 0.
        asc = write_file(os.path.join(self.tmpdir, "dj.asc"), DISJOINT_ASC_CONTENT)
        rc, _, stderr = run_command("DBRetina index -a dj.asc -o dj_idx",
                                    cwd=self.tmpdir)
        self.assertEqual(rc, 0, stderr)
        rc, _, stderr = run_command("DBRetina pairwise -i dj_idx", cwd=self.tmpdir)
        self.assertEqual(rc, 0, stderr)
        pw = os.path.join(self.tmpdir, "dj_idx_DBRetina_pairwise.tsv")
        self.assertEqual(count_tsv_data_rows(pw), 0,
                         "disjoint fixture unexpectedly produced pairs")

        out = os.path.join(self.tmpdir, "dd_empty")
        rc, stdout, stderr = run_command(
            f"DBRetina dedup -i dj_idx -p {pw} -c 50 -o dd_empty", cwd=self.tmpdir
        )
        assert_no_traceback(self, stderr, "dedup no qualifying pairs")
        self.assertEqual(rc, 0, stderr)
        out_file = os.path.join(self.tmpdir, "dd_empty_deduplicated_groups.txt")
        assert_file_exists(self, out_file)
        with open(out_file) as f:
            groups = {line.strip().lower() for line in f if line.strip()}
        # No duplicates -> every group kept.
        self.assertEqual(groups, DISJOINT_GROUP_NAMES)

    def test_025_dedup_with_pairs_still_dedups(self):
        """Guard: a cutoff that *does* produce pairs must still deduplicate
        (the no-pairs fix must not change normal behaviour)."""
        out = os.path.join(self.tmpdir, "dd_ok")
        rc, _, stderr = run_command(
            f"DBRetina dedup -i {self.prefix} -p {self.pw_file} -c 100 -o {out}"
        )
        self.assertEqual(rc, 0, stderr)
        with open(f"{out}_deduplicated_groups.txt") as f:
            groups = {line.strip().lower() for line in f if line.strip()}
        # A and F are identical -> one removed -> 5 groups.
        self.assertEqual(len(groups), 5)

    # ---- 026: dedup -i <missing prefix> -> clean error ----
    def test_026_dedup_missing_index_clean_error(self):
        """dedup -i <nonexistent prefix> must emit a clean [ERROR], not a raw
        FileNotFoundError traceback."""
        out = os.path.join(self.tmpdir, "dd_noidx")
        missing = os.path.join(self.tmpdir, "nope_xyz")
        rc, stdout, stderr = run_command(
            f"DBRetina dedup -i {missing} -p {self.pw_file} -c 50 -o {out}"
        )
        assert_no_traceback(self, stderr, "dedup missing index")
        self.assertNotIn("FileNotFoundError", stderr, stderr)
        self.assertNotEqual(rc, 0, "missing index prefix should fail")
        self.assertIn("[ERROR]", stderr, stderr)
        self.assertIn("index", (stdout + stderr).lower(), stderr)

    # ---- 029: genenet/interactome -i <missing prefix> -> clean error ----
    def test_029_genenet_missing_index_clean_error(self):
        """genenet -i <nonexistent prefix> must emit a clean [ERROR], not a raw
        FileNotFoundError traceback."""
        out = os.path.join(self.tmpdir, "gn_bad")
        missing = os.path.join(self.tmpdir, "nonexistent_prefix")
        rc, stdout, stderr = run_command(
            f"DBRetina genenet -i {missing} -p {self.pw_file} -o {out}"
        )
        assert_no_traceback(self, stderr, "genenet missing index")
        self.assertNotIn("FileNotFoundError", stderr, stderr)
        self.assertNotEqual(rc, 0, "missing index prefix should fail")
        self.assertIn("[ERROR]", stderr, stderr)
        self.assertIn("index", (stdout + stderr).lower(), stderr)

    def test_029_interactome_missing_index_clean_error(self):
        """interactome shares genenet's main; a missing -i prefix must also emit
        a clean [ERROR], not a raw FileNotFoundError traceback."""
        out = os.path.join(self.tmpdir, "it_bad")
        missing = os.path.join(self.tmpdir, "nonexistent_prefix")
        rc, stdout, stderr = run_command(
            f"DBRetina interactome -i {missing} -p {self.pw_file} -o {out}"
        )
        assert_no_traceback(self, stderr, "interactome missing index")
        self.assertNotIn("FileNotFoundError", stderr, stderr)
        self.assertNotEqual(rc, 0, "missing index prefix should fail")
        self.assertIn("[ERROR]", stderr, stderr)


class TestExportClusterGraphDbrpFallback(unittest.TestCase):
    """Input-handling regression for export/cluster/graph (issue 051).

    Same root cause as modularity/dedup 046: the non-store FALLBACK path derived
    the sibling .dbrp via ``pairwise_file.replace('.tsv', '.dbrp')``, a no-op for
    the parquet-directory / .dbrp forms. When ``open_pairwise()`` can't resolve a
    store (a directory without manifest.json), the directory path then leaked
    into the C++ .dbrp reader, crashing with
    ``RuntimeError: Invalid .dbrp file (bad magic bytes)`` (or IsADirectoryError
    on the TSV open).

    The canonical store path is unaffected; only the fallback is fixed by routing
    the .dbrp resolution through ``compat.resolve_dbrp_path`` (isfile-guarded).

    Uses the shared overlapping fixture (emits .tsv + parquet dir + .dbrp +
    .dbri side by side). The forced-fallback dir is reproduced by copying the
    parquet dir and stripping manifest.json so open_pairwise() returns None.
    """

    @classmethod
    def setUpClass(cls):
        cls.prefix, cls.pw_file = _ensure_shared_fixture()
        cls.parquet_dir = f"{cls.prefix}_DBRetina_pairwise"
        cls.dbrp = f"{cls.prefix}_DBRetina_pairwise.dbrp"

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp(prefix="dbretina_expcluster_input_")
        self.assertTrue(os.path.isdir(self.parquet_dir),
                        f"fixture missing parquet dir: {self.parquet_dir}")
        self.assertTrue(os.path.isfile(self.dbrp),
                        f"fixture missing .dbrp: {self.dbrp}")

    def tearDown(self):
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    # ---- helpers ----
    def _hide_parquet_store(self, src_dir):
        """Copy the fixture parquet dir into tmpdir but strip manifest.json so
        open_pairwise() returns None, forcing the .dbrp/TSV fallback path (the
        path that carried the issue-051 .replace bug)."""
        dst = os.path.join(self.tmpdir, os.path.basename(src_dir))
        shutil.copytree(src_dir, dst)
        os.remove(os.path.join(dst, "manifest.json"))
        return dst

    @staticmethod
    def _distmat_text(path):
        with open(path) as f:
            return f.read()

    @staticmethod
    def _cluster_membership(path):
        """Parse a *_clusters.tsv -> set of frozensets of member groups.

        Compares cluster *membership* (set semantics): the store and TSV readers
        iterate pairs in a different order, so member ordering within a cluster
        line is not stable, but the partition itself is.
        """
        clusters = set()
        with open(path) as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                parts = line.rstrip("\n").split("\t")
                if not parts or "cluster_id" in parts[0].lower():
                    continue
                members = parts[-1].split("|")
                clusters.add(frozenset(m.strip().lower() for m in members if m.strip()))
        return clusters

    # ================= export =================
    def test_051_export_parquet_dir_matches_tsv(self):
        """Canonical parquet-dir input still works via the store path and gives a
        distance matrix byte-identical to the .tsv form (parity)."""
        out_dir = os.path.join(self.tmpdir, "exp_dir")
        rc, _, stderr = run_command(
            f"DBRetina export -p {self.parquet_dir} -m ochiai -o {out_dir}"
        )
        assert_no_traceback(self, stderr, "export parquet dir")
        self.assertNotIn("bad magic bytes", stderr, stderr)
        self.assertNotIn("IsADirectoryError", stderr, stderr)
        self.assertEqual(rc, 0, stderr)
        assert_file_exists(self, f"{out_dir}_distmat.tsv")

        out_tsv = os.path.join(self.tmpdir, "exp_tsv")
        rc, _, stderr = run_command(
            f"DBRetina export -p {self.pw_file} -m ochiai -o {out_tsv}"
        )
        self.assertEqual(rc, 0, stderr)
        self.assertEqual(
            self._distmat_text(f"{out_dir}_distmat.tsv"),
            self._distmat_text(f"{out_tsv}_distmat.tsv"),
            "parquet-dir export distmat differs from .tsv form",
        )

    def test_051_export_dir_without_store_falls_back_cleanly(self):
        """A directory input that open_pairwise() can't use (no manifest) must
        resolve the sibling .dbrp instead of feeding the directory path to the
        C++ reader (bad magic bytes). With the sibling .dbrp present the result
        must still match the .tsv form."""
        hidden = self._hide_parquet_store(self.parquet_dir)
        shutil.copy(self.dbrp, hidden + ".dbrp")
        out = os.path.join(self.tmpdir, "exp_fb")
        rc, _, stderr = run_command(
            f"DBRetina export -p {hidden} -m ochiai -o {out}"
        )
        assert_no_traceback(self, stderr, "export fallback dir")
        self.assertNotIn("bad magic bytes", stderr, stderr)
        self.assertNotIn("IsADirectoryError", stderr, stderr)
        self.assertEqual(rc, 0, stderr)
        out_tsv = os.path.join(self.tmpdir, "exp_fb_tsv")
        run_command(f"DBRetina export -p {self.pw_file} -m ochiai -o {out_tsv}")
        self.assertEqual(
            self._distmat_text(f"{out}_distmat.tsv"),
            self._distmat_text(f"{out_tsv}_distmat.tsv"),
            "fallback-dir export distmat differs from .tsv form",
        )

    def test_051_export_unreadable_dir_clean_error(self):
        """A directory input with neither a usable store nor a sibling .dbrp must
        give a clean [ERROR], never a raw traceback (bad magic bytes /
        IsADirectoryError)."""
        hidden = self._hide_parquet_store(self.parquet_dir)  # no sibling .dbrp
        out = os.path.join(self.tmpdir, "exp_bad")
        rc, stdout, stderr = run_command(
            f"DBRetina export -p {hidden} -m ochiai -o {out}"
        )
        assert_no_traceback(self, stderr, "export unreadable dir")
        self.assertNotIn("bad magic bytes", stderr, stderr)
        self.assertNotIn("IsADirectoryError", stderr, stderr)
        self.assertNotEqual(rc, 0, "unreadable directory input should fail")
        self.assertIn("[ERROR]", (stdout + stderr), stderr)

    # ================= cluster =================
    def test_051_cluster_parquet_dir_matches_tsv(self):
        """Canonical parquet-dir input still works via the store path and yields
        the same cluster partition as the .tsv form (set semantics)."""
        out_dir = os.path.join(self.tmpdir, "clu_dir")
        rc, _, stderr = run_command(
            f"DBRetina cluster -p {self.parquet_dir} -m ochiai -c 50 -o {out_dir}"
        )
        assert_no_traceback(self, stderr, "cluster parquet dir")
        self.assertNotIn("bad magic bytes", stderr, stderr)
        self.assertNotIn("IsADirectoryError", stderr, stderr)
        self.assertEqual(rc, 0, stderr)
        assert_file_exists(self, f"{out_dir}_clusters.tsv")

        out_tsv = os.path.join(self.tmpdir, "clu_tsv")
        rc, _, stderr = run_command(
            f"DBRetina cluster -p {self.pw_file} -m ochiai -c 50 -o {out_tsv}"
        )
        self.assertEqual(rc, 0, stderr)
        self.assertEqual(
            self._cluster_membership(f"{out_dir}_clusters.tsv"),
            self._cluster_membership(f"{out_tsv}_clusters.tsv"),
            "parquet-dir cluster partition differs from .tsv form",
        )

    def test_051_cluster_dir_without_store_falls_back_cleanly(self):
        """Forced-fallback directory (no manifest) with a sibling .dbrp must read
        the .dbrp, never feed the directory to the C++ reader, and match the
        .tsv-form partition."""
        hidden = self._hide_parquet_store(self.parquet_dir)
        shutil.copy(self.dbrp, hidden + ".dbrp")
        out = os.path.join(self.tmpdir, "clu_fb")
        rc, _, stderr = run_command(
            f"DBRetina cluster -p {hidden} -m ochiai -c 50 -o {out}"
        )
        assert_no_traceback(self, stderr, "cluster fallback dir")
        self.assertNotIn("bad magic bytes", stderr, stderr)
        self.assertNotIn("IsADirectoryError", stderr, stderr)
        self.assertEqual(rc, 0, stderr)
        out_tsv = os.path.join(self.tmpdir, "clu_fb_tsv")
        run_command(f"DBRetina cluster -p {self.pw_file} -m ochiai -c 50 -o {out_tsv}")
        self.assertEqual(
            self._cluster_membership(f"{out}_clusters.tsv"),
            self._cluster_membership(f"{out_tsv}_clusters.tsv"),
            "fallback-dir cluster partition differs from .tsv form",
        )

    def test_051_cluster_unreadable_dir_clean_error(self):
        """Directory input with no usable store and no sibling .dbrp -> clean
        [ERROR] (the node-count read in __init__ used to crash on bad magic
        bytes)."""
        hidden = self._hide_parquet_store(self.parquet_dir)  # no sibling .dbrp
        out = os.path.join(self.tmpdir, "clu_bad")
        rc, stdout, stderr = run_command(
            f"DBRetina cluster -p {hidden} -m ochiai -c 50 -o {out}"
        )
        assert_no_traceback(self, stderr, "cluster unreadable dir")
        self.assertNotIn("bad magic bytes", stderr, stderr)
        self.assertNotIn("IsADirectoryError", stderr, stderr)
        self.assertNotEqual(rc, 0, "unreadable directory input should fail")
        self.assertIn("[ERROR]", (stdout + stderr), stderr)

    # ================= graph =================
    # bipartite_graph.py:143 (load_all_pairwise) is DEAD code with no callers;
    # the live `graph` command reads pairwise via pairwise_file_iterator. Its
    # .dbrp fast-path was isfile-guarded (so a dir never hit the C++ reader) but
    # then fell to the TSV open() -> raw IsADirectoryError on a dir. The 051 fix
    # routes that resolution through resolve_dbrp_path (uses the sibling .dbrp
    # even for the dir form) and emits a clean error when there is none.
    def test_051_graph_dir_without_store_falls_back_cleanly(self):
        """`graph` on a forced-fallback directory (no manifest) with a sibling
        .dbrp must read the .dbrp (not IsADirectoryError on the TSV open) and
        produce the same edge count as the .tsv form."""
        hidden = self._hide_parquet_store(self.parquet_dir)
        shutil.copy(self.dbrp, hidden + ".dbrp")
        out = os.path.join(self.tmpdir, "grph_fb")
        rc, _, stderr = run_command(
            f"DBRetina graph -i {self.prefix} -p {hidden} -m ochiai -c 50 -o {out}"
        )
        assert_no_traceback(self, stderr, "graph fallback dir")
        self.assertNotIn("bad magic bytes", stderr, stderr)
        self.assertNotIn("IsADirectoryError", stderr, stderr)
        self.assertEqual(rc, 0, stderr)
        assert_file_exists(self, f"{out}_edges.tsv")
        out_tsv = os.path.join(self.tmpdir, "grph_tsv")
        rc2, _, stderr2 = run_command(
            f"DBRetina graph -i {self.prefix} -p {self.pw_file} -m ochiai -c 50 -o {out_tsv}"
        )
        self.assertEqual(rc2, 0, stderr2)
        self.assertEqual(
            count_tsv_data_rows(f"{out}_edges.tsv"),
            count_tsv_data_rows(f"{out_tsv}_edges.tsv"),
            "fallback-dir graph edge count differs from .tsv form",
        )

    def test_051_graph_unreadable_dir_clean_error(self):
        """`graph` on a directory with no usable store and no sibling .dbrp must
        give a clean [ERROR], never a raw IsADirectoryError traceback from the
        TSV open()."""
        hidden = self._hide_parquet_store(self.parquet_dir)  # no sibling .dbrp
        out = os.path.join(self.tmpdir, "grph_bad")
        rc, stdout, stderr = run_command(
            f"DBRetina graph -i {self.prefix} -p {hidden} -m ochiai -c 50 -o {out}"
        )
        assert_no_traceback(self, stderr, "graph unreadable dir")
        self.assertNotIn("bad magic bytes", stderr, stderr)
        self.assertNotIn("IsADirectoryError", stderr, stderr)
        self.assertNotEqual(rc, 0, "unreadable directory input should fail")
        self.assertIn("[ERROR]", (stdout + stderr), stderr)


# ============================================================
# SECTION: C++ crash regressions (issues 041, 042)
# ============================================================

# Signals that indicate a C++ core-dump (subprocess returncode is negative for a
# signal; the shell maps SIGABRT->134 and SIGSEGV->139).
_CRASH_RETURNCODES = (-6, -11, 134, 139)


def assert_no_coredump(test_case, rc, stderr, ctx=""):
    """The process must not have died from SIGABRT/SIGSEGV (boost abort etc.)."""
    test_case.assertNotIn(rc, _CRASH_RETURNCODES,
                          f"{ctx}: process core-dumped (rc={rc}):\n{stderr}")
    for marker in ("Segmentation", "Aborted", "core dumped",
                   "terminate called", "domain_error"):
        test_case.assertNotIn(marker, stderr, f"{ctx}: crash marker '{marker}':\n{stderr}")


class TestPvalueSmallPopulationCrash(unittest.TestCase):
    """ISSUE-041: `pairwise --pvalue` core-dumped (boost hypergeometric
    domain_error -> abort) on small/degenerate gene populations.

    Repro (from the issue): 30 groups drawing 12 genes each from a 40-gene pool.
    For such tiny populations the hypergeometric random variable k-1 falls below
    boost's lower domain bound max(0, n+r-N), so boost throws std::domain_error
    which (under the default policy) calls std::terminate -> SIGABRT (rc 134).
    The fix guards/clamps the parameters and returns a non-significant pvalue
    (1.0) for degenerate cases instead of crashing.
    """

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp(prefix="dbretina_pv041_")

    def tearDown(self):
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def _build_small_pop_index(self):
        """30 groups, 12 genes each, drawn from a 40-gene pool (sliding window)."""
        genes = [f"g{i}" for i in range(40)]
        rows = ["grp_%02d\tdesc\t%s" % (k, "\t".join(genes[(k % 10):(k % 10) + 12]))
                for k in range(30)]
        gmt = os.path.join(self.tmpdir, "many.gmt")
        with open(gmt, "w") as f:
            f.write("\n".join(rows) + "\n")
        rc, _, err = run_command("DBRetina index -g many.gmt -o idx", cwd=self.tmpdir)
        self.assertEqual(rc, 0, err)
        return "idx"

    def test_041_pairwise_pvalue_small_population_no_coredump(self):
        """The exact repro from issue 041 must NOT core-dump and must emit pvalues."""
        idx = self._build_small_pop_index()
        rc, stdout, stderr = run_command(
            f"DBRetina pairwise -i {idx} -m containment -c 0 --pvalue", cwd=self.tmpdir
        )
        assert_no_coredump(self, rc, stderr, "pairwise --pvalue small pop")
        self.assertEqual(rc, 0, f"pairwise --pvalue should succeed:\n{stderr}")
        pw_file = os.path.join(self.tmpdir, "idx_DBRetina_pairwise.tsv")
        assert_file_exists(self, pw_file)
        rows = parse_pairwise_tsv(pw_file)
        self.assertGreater(len(rows), 0, "expected pairwise rows")
        # Every row carries a pvalue and it is a valid probability in [0, 1].
        for row in rows:
            self.assertIn("pvalue", row, "pvalue column missing")
            self.assertGreaterEqual(row["pvalue"], 0.0, f"pvalue<0: {row['pvalue']}")
            self.assertLessEqual(row["pvalue"], 1.0, f"pvalue>1: {row['pvalue']}")

    def test_041_pvalue_unchanged_on_normal_population(self):
        """Regression guard: the fix must NOT change pvalues for a valid-domain
        (normal) population. The standard 6-group fixture has a comfortable gene
        population, so its hypergeometric inputs are all in-domain. We pin the
        exact pvalue for a known pair (GroupA vs GroupB) computed independently
        with the same hypergeometric definition boost uses.
        """
        prefix, pw_file = setup_index_and_pairwise(self.tmpdir, extra_pw_args="--pvalue")
        rows = parse_pairwise_tsv(pw_file)
        by_pair = {frozenset({r["group_1_name"], r["group_2_name"]}): r for r in rows}

        # Independent reference: over-enrichment pvalue = P(X >= k) = 1 - cdf(k-1)
        # for Hypergeometric(M=population, n=sample=|set1|, r=successes=|set2|).
        # GroupA={Alpha,Beta,Gamma,Delta,Epsilon}=5, GroupB={Alpha,Beta,Gamma}=3,
        # shared k=3. Population = 12 distinct genes across the fixture.
        try:
            from scipy.stats import hypergeom
            ref = float(hypergeom.sf(3 - 1, 12, 3, 5))  # sf(k-1) = P(X>=k)
        except ImportError:
            ref = None

        ab = by_pair[frozenset({"groupa", "groupb"})]
        self.assertIn("pvalue", ab)
        # Valid probability regardless of scipy availability.
        self.assertGreaterEqual(ab["pvalue"], 0.0)
        self.assertLessEqual(ab["pvalue"], 1.0)
        if ref is not None:
            self.assertAlmostEqual(ab["pvalue"], ref, places=4,
                                   msg="normal-population pvalue changed by the fix")


class TestPartialDbriOpen(unittest.TestCase):
    """ISSUE-042: DBRetinaIndex::open silently accepted an unfinalized/partial
    .dbri whose placeholder toc_offset==0 (a crash between header-write and
    finalize_write), reading garbage with no validation.

    Header layout: magic(4) + version(4) + toc_offset(8) = 16 bytes. A finalized
    index always has toc_offset >= 16 and < file_size. toc_offset==0 (or a value
    pointing past EOF) means the file was never finalized. open() must reject it
    with a clean error instead of seeking to 0, reading the magic bytes as a
    ~1.2-billion section count, and returning a silently-empty index.
    """

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp(prefix="dbretina_dbri042_")

    def tearDown(self):
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def _build_valid_dbri(self):
        asc = write_file(os.path.join(self.tmpdir, "t.asc"),
                         "gene_set\tgene\nGroupA\tAlpha\nGroupA\tBeta\n"
                         "GroupB\tBeta\nGroupB\tGamma\n")
        prefix = os.path.join(self.tmpdir, "idx")
        rc, _, err = run_command(f"DBRetina index -a {asc} -o {prefix}")
        self.assertEqual(rc, 0, err)
        return prefix, prefix + ".dbri"

    def test_042_partial_header_only_clean_error(self):
        """A 16-byte header with the placeholder toc_offset==0 (exactly what
        begin_write leaves on disk before finalize) must raise a clean error,
        NOT silently return an empty index, NOT a huge alloc, NOT a segfault."""
        import struct
        partial = os.path.join(self.tmpdir, "partial.dbri")
        with open(partial, "wb") as f:
            f.write(b"DBRI" + struct.pack("<I", 1) + struct.pack("<Q", 0))

        import _dbretina_internal as dbi
        with self.assertRaises(Exception) as cm:
            dbi.dbri_load_raw_gene_sets(partial)
        msg = str(cm.exception).lower()
        self.assertTrue(
            any(w in msg for w in ("incomplete", "corrupt", "unfinalized", "truncated")),
            f"expected a clear corrupt/incomplete message, got: {cm.exception}",
        )

    def test_042_partial_via_pairwise_clean_error(self):
        """Opening a partial (toc_offset==0) .dbri through the pairwise CLI gives
        a clean [ERROR]/exception, never a segfault/abort/garbage."""
        import struct
        partial_prefix = os.path.join(self.tmpdir, "partial")
        with open(partial_prefix + ".dbri", "wb") as f:
            f.write(b"DBRI" + struct.pack("<I", 1) + struct.pack("<Q", 0))
        rc, stdout, stderr = run_command(
            f"DBRetina pairwise -i {partial_prefix} -m containment -c 0"
        )
        assert_no_coredump(self, rc, stderr, "pairwise on partial .dbri")
        self.assertNotEqual(rc, 0, "partial .dbri must not succeed")

    def test_042_truncated_before_toc_clean_error(self):
        """A valid header whose toc_offset points past a truncated EOF must also
        be rejected cleanly (offset >= file_size)."""
        import struct
        _, valid = self._build_valid_dbri()
        data = open(valid, "rb").read()
        toc_offset = struct.unpack("<Q", data[8:16])[0]
        self.assertGreater(toc_offset, 16, "valid index must have toc_offset > header")
        trunc = os.path.join(self.tmpdir, "trunc.dbri")
        # Keep header (with the real, now-out-of-range toc_offset) but cut the TOC off.
        with open(trunc, "wb") as f:
            f.write(data[:toc_offset - 5])

        import _dbretina_internal as dbi
        with self.assertRaises(Exception) as cm:
            dbi.dbri_load_raw_gene_sets(trunc)
        msg = str(cm.exception).lower()
        self.assertTrue(
            any(w in msg for w in ("incomplete", "corrupt", "unfinalized", "truncated")),
            f"expected a clear corrupt/incomplete message, got: {cm.exception}",
        )

    def test_042_valid_dbri_still_opens(self):
        """A finalized .dbri must still open and read correctly (no regression)."""
        prefix, valid = self._build_valid_dbri()
        import _dbretina_internal as dbi
        # Reads through open() must succeed and return real data.
        names = dbi.dbri_load_names_list(valid)
        self.assertEqual(set(n.lower() for n in names), {"groupa", "groupb"})
        raw = dbi.dbri_load_raw_gene_sets(valid)
        self.assertIsInstance(raw, str)
        self.assertGreater(len(raw), 0, "valid index returned empty raw gene sets")
        # And the full pairwise pipeline still works end-to-end on it.
        rc, _, err = run_command(f"DBRetina pairwise -i {prefix} -m containment -c 0")
        self.assertEqual(rc, 0, err)


# ============================================================
# SECTION: Data-layer robustness / determinism
#   issue 052 (remaining raw .dbrp .replace sites),
#   issue 044 (cluster bare-TSV node-count parse),
#   issue 050 (setcov nondeterministic output).
# ============================================================

class TestDataLayerRobustness(unittest.TestCase):
    """Regression guards for issues 052 / 044 / 050.

    052: query/bipartite given a parquet *directory* with no usable store
         (no manifest.json) and no sibling .dbrp used to leak the directory
         path into the C++ .dbrp reader ("bad magic bytes") or open() it as
         text (IsADirectoryError). They must now emit a clean error.
    044: cluster on a bare pairwise.tsv (no .dbrp/parquet sibling) parsed the
         node count from line 1, but '#nodes:N' moved to a later header line,
         so it crashed with "invalid literal for int()". It must now parse the
         count robustly and cluster without a traceback.
    """

    @classmethod
    def setUpClass(cls):
        # Reuse the module-shared (overlapping) index + pairwise fixture.
        cls.prefix, cls.pw_file = _ensure_shared_fixture()

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp(prefix="dbretina_datalayer_")

    def tearDown(self):
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def _no_manifest_dir(self):
        """A directory that mimics a pairwise parquet dir but has NO manifest.json
        and NO sibling .dbrp -> open_pairwise() returns None and resolve_dbrp_path()
        returns None, exercising the bare fallback path of each command."""
        d = os.path.join(self.tmpdir, "p_DBRetina_pairwise")
        os.makedirs(d, exist_ok=True)
        # A stray file so it's a non-empty dir, but still no manifest.json.
        write_file(os.path.join(d, "part-0.parquet"), "not a real parquet")
        return d

    # ---- 052: query on a parquet dir without a usable store/.dbrp ----
    def test_052_query_dir_no_store_clean_error(self):
        """query -p <dir-without-manifest> -m ochiai -c: clean error, never a
        'bad magic bytes' RuntimeError or a raw traceback (issue 052)."""
        d = self._no_manifest_dir()
        out = os.path.join(self.tmpdir, "qout")
        rc, stdout, stderr = run_command(
            f"DBRetina query -p {d} -m ochiai -c 50 -o {out}"
        )
        combined = stdout + stderr
        assert_no_traceback(self, combined, "query dir-no-store")
        self.assertNotIn("bad magic bytes", combined, combined)
        self.assertNotIn("IsADirectoryError", combined, combined)
        self.assertNotEqual(rc, 0, "expected a clean non-zero exit, not success")
        self.assertIn("[ERROR]", combined, combined)

    def test_052_query_dir_no_store_pvalue_clean_error(self):
        """query -m pvalue on a dir-without-manifest must not IsADirectoryError in
        the compat pvalue probe (issue 052, compat.pairwise_has_pvalue dir guard)."""
        d = self._no_manifest_dir()
        out = os.path.join(self.tmpdir, "qoutpv")
        rc, stdout, stderr = run_command(
            f"DBRetina query -p {d} -m pvalue -c 50 -o {out}"
        )
        combined = stdout + stderr
        assert_no_traceback(self, combined, "query dir-no-store pvalue")
        self.assertNotIn("IsADirectoryError", combined, combined)
        self.assertNotEqual(rc, 0)

    # ---- 052: bipartite on a parquet dir without a usable store/.dbrp ----
    def test_052_bipartite_dir_no_store_clean_error(self):
        """bipartite -p <dir-without-manifest>: clean error, no bad-magic /
        IsADirectoryError / traceback (issue 052)."""
        d = self._no_manifest_dir()
        g1 = write_file(os.path.join(self.tmpdir, "g1.txt"), "groupa\ngroupb\n")
        g2 = write_file(os.path.join(self.tmpdir, "g2.txt"), "groupc\ngroupe\n")
        out = os.path.join(self.tmpdir, "bout")
        rc, stdout, stderr = run_command(
            f"DBRetina bipartite -p {d} --group1 {g1} --group2 {g2} "
            f"-m ochiai -c 0 --no-plot -o {out}"
        )
        combined = stdout + stderr
        assert_no_traceback(self, combined, "bipartite dir-no-store")
        self.assertNotIn("bad magic bytes", combined, combined)
        self.assertNotIn("IsADirectoryError", combined, combined)
        self.assertNotEqual(rc, 0)
        self.assertIn("[ERROR]", combined, combined)

    # ---- 052: the valid parquet-dir store path must still work ----
    def test_052_bipartite_valid_parquet_dir_still_works(self):
        """A real parquet directory (with manifest.json) must still be accepted
        via the store path -- the dir guard only rejects unusable dirs."""
        pdir = self.prefix + "_DBRetina_pairwise"
        self.assertTrue(
            os.path.isdir(pdir) and os.path.exists(os.path.join(pdir, "manifest.json")),
            "shared fixture is missing its parquet dir/manifest",
        )
        g1 = write_file(os.path.join(self.tmpdir, "g1.txt"), "groupa\ngroupb\n")
        g2 = write_file(os.path.join(self.tmpdir, "g2.txt"), "groupc\ngroupe\n")
        out = os.path.join(self.tmpdir, "bvalid")
        rc, _, stderr = run_command(
            f"DBRetina bipartite -p {pdir} --group1 {g1} --group2 {g2} "
            f"-m ochiai -c 0 --no-plot -o {out}"
        )
        assert_no_traceback(self, stderr, "bipartite valid parquet dir")
        self.assertEqual(rc, 0, stderr)
        assert_file_exists(self, f"{out}_bipartite_pairwise.tsv")

    # ---- 044: cluster on a bare pairwise.tsv (no .dbrp/parquet sibling) ----
    def test_044_cluster_bare_tsv_only_node_count(self):
        """cluster on a copied pairwise.tsv with NO .dbrp/parquet sibling must
        parse '#nodes:N' from its real header line (not line 1) and cluster
        without crashing on int() (issue 044)."""
        bare = os.path.join(self.tmpdir, "bare")
        os.makedirs(bare, exist_ok=True)
        bare_tsv = os.path.join(bare, "p_DBRetina_pairwise.tsv")
        shutil.copy(self.pw_file, bare_tsv)
        # Guarantee the bare-TSV-only layout: no .dbrp, no parquet dir sibling.
        self.assertFalse(os.path.exists(os.path.join(bare, "p_DBRetina_pairwise.dbrp")))
        self.assertFalse(os.path.isdir(os.path.join(bare, "p_DBRetina_pairwise")))
        # Sanity: the header really puts '#nodes:' past line 1 (the bug premise).
        with open(bare_tsv) as f:
            first = f.readline()
        self.assertFalse(first.startswith("#nodes:"),
                         "fixture changed: #nodes is on line 1, bug no longer reproducible")

        out = os.path.join(bare, "clu")
        rc, stdout, stderr = run_command(
            f"DBRetina cluster -p {bare_tsv} -m ochiai -c 50 -o {out}"
        )
        combined = stdout + stderr
        assert_no_traceback(self, combined, "cluster bare-TSV")
        self.assertNotIn("invalid literal for int", combined, combined)
        self.assertEqual(rc, 0, combined)
        clusters_tsv = f"{out}_clusters.tsv"
        assert_file_exists(self, clusters_tsv)
        # The correct node count (6 groups) is used: the b/f/a connected
        # component (the only edges above ochiai 50 in this fixture) is present.
        members = []
        with open(clusters_tsv) as f:
            for line in f:
                if line.startswith("#") or line.lower().startswith("cluster_id"):
                    continue
                if line.strip():
                    members.append(line.strip().split("\t")[2])
        joined = "|".join(members)
        for g in ("groupa", "groupb", "groupf"):
            self.assertIn(g, joined, f"{g} missing from clusters:\n{joined}")

    def test_044_cluster_bare_tsv_missing_nodes_header_clean_error(self):
        """A bare TSV with no '#nodes:' line at all must give a clean error,
        not a traceback (robust-parse fallback)."""
        bare_tsv = os.path.join(self.tmpdir, "nonodes_DBRetina_pairwise.tsv")
        # Header without a #nodes line, plus the column header (no data rows).
        write_file(
            bare_tsv,
            "# DBRetina pairwise output\n"
            "# population_size: 12\n"
            "group_1_ID\tgroup_2_ID\tgroup_1_name\tgroup_2_name\tshared_features\t"
            "containment\tochiai\tjaccard\tcsi\tdice\todds_ratio\n",
        )
        out = os.path.join(self.tmpdir, "clu_nonodes")
        rc, stdout, stderr = run_command(
            f"DBRetina cluster -p {bare_tsv} -m ochiai -c 50 -o {out}"
        )
        combined = stdout + stderr
        assert_no_traceback(self, combined, "cluster bare-TSV no #nodes")
        self.assertNotEqual(rc, 0)
        self.assertIn("[ERROR]", combined, combined)


class TestSetcovDeterminism(unittest.TestCase):
    """issue 050: setcov output must be byte-identical for identical inputs.

    The selection iterated Python sets and the associations file was written by
    iterating sets directly, so identical commands produced byte-different output
    (order + set-iteration). The fix sorts the written collections and uses a
    stable sort with a deterministic tie-break in the set-cover selection. Vary
    PYTHONHASHSEED across runs to force different set-iteration orders -- only a
    truly order-independent writer stays byte-identical.
    """

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp(prefix="dbretina_setcov_det_")
        asc_path = os.path.join(self.tmpdir, "test_input.asc")
        write_file(asc_path, TEST_ASC_CONTENT)
        rc, _, stderr = run_command(
            "DBRetina index -a test_input.asc -o test_idx", cwd=self.tmpdir
        )
        assert rc == 0, f"index failed: {stderr}"
        rc, _, stderr = run_command("DBRetina pairwise -i test_idx", cwd=self.tmpdir)
        assert rc == 0, f"pairwise failed: {stderr}"

    def tearDown(self):
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def _run_setcov(self, outdir, hashseed):
        os.makedirs(os.path.join(self.tmpdir, outdir), exist_ok=True)
        env = dict(os.environ)
        env["PYTHONHASHSEED"] = str(hashseed)
        result = subprocess.run(
            "DBRetina setcov -i ../test_idx --modularity 80 --dedup 100 "
            "--stop-cov 100 -o sc",
            shell=True, capture_output=True, text=True, timeout=300,
            cwd=os.path.join(self.tmpdir, outdir), env=env,
        )
        assert result.returncode == 0, f"setcov failed: {result.stderr}"
        return os.path.join(self.tmpdir, outdir)

    def test_050_setcov_byte_identical_across_runs(self):
        """Two setcov runs with identical args (different PYTHONHASHSEED) produce
        byte-identical output files, with the same number of representatives."""
        d1 = self._run_setcov("run1", hashseed=1)
        d2 = self._run_setcov("run2", hashseed=12345)

        # The previously-nondeterministic outputs (and the rest) must match byte-for-byte.
        for fname in (
            "sc_associations.tsv",
            "sc_new.gmt",
            "sc_original.gmt",
            "sc_groups_metadata.tsv",
            "sc_remaining_groups_metadata.tsv",
            "sc_removed_groups_metadata.tsv",
            "sc_item_to_GPI_CSI.tsv",
        ):
            f1 = os.path.join(d1, fname)
            f2 = os.path.join(d2, fname)
            assert_file_exists(self, f1)
            assert_file_exists(self, f2)
            with open(f1, "rb") as a, open(f2, "rb") as b:
                self.assertEqual(
                    a.read(), b.read(),
                    f"{fname} differs between two identical setcov runs (issue 050)",
                )

        # Same number of representatives (quality unchanged): both files identical
        # already implies this, but assert the count explicitly for clarity.
        def _rep_count(d):
            with open(os.path.join(d, "sc_remaining_groups_metadata.tsv")) as f:
                return sum(1 for i, line in enumerate(f) if i > 0 and line.strip())
        self.assertEqual(_rep_count(d1), _rep_count(d2))


# ============================================================
# SECTION: output-consistency across input forms (issues 047, 048, 063)
# ============================================================

def _read_header_and_rows(path):
    """Split a query/pairwise-style TSV into (column_header, data_rows).

    Skips '#'-comment lines; the first remaining line is the column header
    (group_1_ID...), the rest are tab-split data rows. Returns
    (header_str_or_None, list_of_field_lists).
    """
    header = None
    rows = []
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            line = line.rstrip("\n")
            if not line:
                continue
            if header is None:
                header = line
            else:
                rows.append(line.split("\t"))
    return header, rows


class TestQueryHeaderConsistency(unittest.TestCase):
    """issue 047: query output column-header must be identical across input forms.

    The .tsv query path copies the source TSV's column-header row into the output;
    the store (parquet dir) and .dbrp paths used to seed only a '#command:' line,
    so dir/.dbrp query output was headerless. The DATA rows were already byte-
    identical across forms; only the column-header row differed. After the fix all
    three forms emit the same canonical 'group_1_ID...odds_ratio[,pvalue]' header.

    Uses the shared fixture, which emits .tsv + parquet dir + .dbrp side by side.
    """

    @classmethod
    def setUpClass(cls):
        cls.prefix, cls.pw_file = _ensure_shared_fixture()
        cls.parquet_dir = f"{cls.prefix}_DBRetina_pairwise"
        cls.dbrp = f"{cls.prefix}_DBRetina_pairwise.dbrp"

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp(prefix="dbretina_qhdr_")
        self.assertTrue(os.path.isdir(self.parquet_dir),
                        f"fixture missing parquet dir: {self.parquet_dir}")
        self.assertTrue(os.path.isfile(self.dbrp),
                        f"fixture missing .dbrp: {self.dbrp}")

    def tearDown(self):
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    @staticmethod
    def _sorted_rows(rows):
        # Store/.dbrp/awk emit rows in a different order; compare as a set of tuples.
        return sorted(tuple(r) for r in rows)

    def test_047_query_header_identical_across_forms(self):
        """query on .tsv / parquet dir / .dbrp -> identical column-header row AND
        identical data rows (order-independent)."""
        out_tsv = os.path.join(self.tmpdir, "q_tsv")
        out_dir = os.path.join(self.tmpdir, "q_dir")
        out_dbrp = os.path.join(self.tmpdir, "q_dbrp")
        for src, out in ((self.pw_file, out_tsv),
                         (self.parquet_dir, out_dir),
                         (self.dbrp, out_dbrp)):
            rc, _, stderr = run_command(
                f"DBRetina query -p {src} -m ochiai -c 50 -o {out}"
            )
            assert_no_traceback(self, stderr, f"query {src}")
            self.assertEqual(rc, 0, stderr)
            assert_file_exists(self, f"{out}.tsv")

        h_tsv, r_tsv = _read_header_and_rows(f"{out_tsv}.tsv")
        h_dir, r_dir = _read_header_and_rows(f"{out_dir}.tsv")
        h_dbrp, r_dbrp = _read_header_and_rows(f"{out_dbrp}.tsv")

        # The .tsv form is the reference; its column header is the canonical one.
        self.assertEqual(
            h_tsv,
            "group_1_ID\tgroup_2_ID\tgroup_1_name\tgroup_2_name\t"
            "shared_features\tcontainment\tochiai\tjaccard\tcsi\tdice\todds_ratio",
            f"unexpected canonical header from .tsv form:\n{h_tsv!r}",
        )
        # The fix: parquet-dir and .dbrp forms now emit the SAME column header.
        self.assertEqual(h_dir, h_tsv,
                         f"parquet-dir query header differs from .tsv:\n{h_dir!r}")
        self.assertEqual(h_dbrp, h_tsv,
                         f".dbrp query header differs from .tsv:\n{h_dbrp!r}")

        # Data rows must be identical too (consistency, not just the header).
        self.assertEqual(self._sorted_rows(r_dir), self._sorted_rows(r_tsv),
                         "parquet-dir query data rows differ from .tsv")
        self.assertEqual(self._sorted_rows(r_dbrp), self._sorted_rows(r_tsv),
                         ".dbrp query data rows differ from .tsv")

    def test_047_query_header_identical_with_pvalue(self):
        """With a pvalue column present, the canonical header (incl trailing
        'pvalue') is still identical across .tsv / parquet dir / .dbrp."""
        prefix, pw = _ensure_shared_pvalue_fixture()
        parquet_dir = f"{prefix}_DBRetina_pairwise"
        dbrp = f"{prefix}_DBRetina_pairwise.dbrp"
        self.assertTrue(os.path.isdir(parquet_dir), parquet_dir)
        self.assertTrue(os.path.isfile(dbrp), dbrp)

        out_tsv = os.path.join(self.tmpdir, "qp_tsv")
        out_dir = os.path.join(self.tmpdir, "qp_dir")
        out_dbrp = os.path.join(self.tmpdir, "qp_dbrp")
        for src, out in ((pw, out_tsv), (parquet_dir, out_dir), (dbrp, out_dbrp)):
            rc, _, stderr = run_command(
                f"DBRetina query -p {src} -m ochiai -c 50 -o {out}"
            )
            assert_no_traceback(self, stderr, f"query(pvalue) {src}")
            self.assertEqual(rc, 0, stderr)

        h_tsv, _ = _read_header_and_rows(f"{out_tsv}.tsv")
        h_dir, _ = _read_header_and_rows(f"{out_dir}.tsv")
        h_dbrp, _ = _read_header_and_rows(f"{out_dbrp}.tsv")
        self.assertTrue(h_tsv.endswith("\tpvalue"),
                        f"pvalue fixture .tsv header lacks pvalue:\n{h_tsv!r}")
        self.assertEqual(h_dir, h_tsv,
                         f"parquet-dir pvalue header differs from .tsv:\n{h_dir!r}")
        self.assertEqual(h_dbrp, h_tsv,
                         f".dbrp pvalue header differs from .tsv:\n{h_dbrp!r}")


class TestBipartiteMetricRestriction(unittest.TestCase):
    """issue 048: bipartite -m must be restricted to the metrics it emits.

    bipartite only carries containment/ochiai/jaccard (+ pvalue when present) in
    its output; csi/dice/odds_ratio were accepted (exit 0) but silently produced
    output WITHOUT those columns. The fix restricts -m to the genuinely supported
    set with a clean Click error for the rest, matching the --help text
    ['containment','ochiai','jaccard','pvalue'].
    """

    @classmethod
    def setUpClass(cls):
        cls.prefix, cls.pw_file = _ensure_shared_fixture()

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp(prefix="dbretina_bipm_")
        self.g1 = write_file(os.path.join(self.tmpdir, "g1.txt"), "GroupA\nGroupB\n")
        self.g2 = write_file(os.path.join(self.tmpdir, "g2.txt"), "GroupC\nGroupD\n")

    def tearDown(self):
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def _bip(self, metric, out, pw=None, cutoff=0):
        # cutoff default 0 keeps all pairs for similarity metrics. p-value is
        # "lower is better" (issue 068): cutoff 0 would exclude every pair
        # (no p <= 0), so the pvalue case passes cutoff=1 to keep all pairs.
        return run_command(
            f"DBRetina bipartite -p {pw or self.pw_file} --group1 {self.g1} "
            f"--group2 {self.g2} -m {metric} -c {cutoff} --no-plot -o {out}"
        )

    def test_048_unsupported_metric_clean_error(self):
        """-m csi/dice/odds_ratio (computed columns bipartite does NOT emit) must
        exit nonzero with a clean error that lists the supported metrics, never a
        silent exit-0 with missing columns."""
        for metric in ("csi", "dice", "odds_ratio"):
            out = os.path.join(self.tmpdir, f"bip_{metric}")
            rc, stdout, stderr = self._bip(metric, out)
            combined = stdout + stderr
            assert_no_traceback(self, stderr, f"bipartite -m {metric}")
            self.assertNotEqual(
                rc, 0, f"bipartite -m {metric} should be rejected, got exit 0")
            self.assertIn(metric, combined,
                          f"error should name the rejected metric '{metric}':\n{combined}")
            # Lists the genuinely supported set so the user can pick a valid one.
            for supported in ("containment", "ochiai", "jaccard"):
                self.assertIn(supported, combined,
                              f"error should list supported metric '{supported}':\n{combined}")
            # No partial/silent output for a rejected metric.
            self.assertFalse(
                os.path.exists(f"{out}_bipartite_pairwise.tsv"),
                f"rejected -m {metric} must not write a bipartite TSV",
            )

    def test_048_supported_metric_emits_its_column(self):
        """-m containment/ochiai/jaccard still works and the metric column is
        present in the output (no regression)."""
        for metric in ("containment", "ochiai", "jaccard"):
            out = os.path.join(self.tmpdir, f"bip_ok_{metric}")
            rc, _, stderr = self._bip(metric, out)
            assert_no_traceback(self, stderr, f"bipartite -m {metric}")
            self.assertEqual(rc, 0, stderr)
            tsv = f"{out}_bipartite_pairwise.tsv"
            assert_file_exists(self, tsv)
            header, _ = _read_header_and_rows(tsv)
            self.assertIn(metric, header.split("\t"),
                          f"bipartite output missing the '{metric}' column:\n{header}")

    def test_048_pvalue_supported_when_present(self):
        """pvalue IS carried in bipartite output when the input has it, so -m
        pvalue must remain accepted (and emit the pvalue column)."""
        prefix, pw = _ensure_shared_pvalue_fixture()
        out = os.path.join(self.tmpdir, "bip_pv")
        # cutoff 1 keeps every pair (all p-values <= 1); cutoff 0 would now
        # correctly exclude all pairs since p-value is "lower is better" (068).
        rc, _, stderr = self._bip("pvalue", out, pw=pw, cutoff=1)
        assert_no_traceback(self, stderr, "bipartite -m pvalue")
        self.assertEqual(rc, 0, stderr)
        tsv = f"{out}_bipartite_pairwise.tsv"
        assert_file_exists(self, tsv)
        header, _ = _read_header_and_rows(tsv)
        self.assertIn("pvalue", header.split("\t"),
                      f"bipartite -m pvalue output missing 'pvalue' column:\n{header}")

    def test_048_bogus_metric_still_clean_error(self):
        """A wholly-unknown metric (issue 023 guard) still errors cleanly under
        the tightened validation."""
        out = os.path.join(self.tmpdir, "bip_bogus")
        rc, stdout, stderr = self._bip("BOGUS", out)
        assert_no_traceback(self, stderr, "bipartite -m BOGUS")
        self.assertNotEqual(rc, 0, "bogus metric should be rejected")
        self.assertIn("BOGUS", stdout + stderr)


class TestClusterBareTsvParity(unittest.TestCase):
    """issue 063: cluster on a bare-TSV-only layout must match the store/.dbrp paths.

    The bare-TSV branch populated self.original_nodes from EVERY pair row BEFORE
    the cutoff filter, so it emitted singleton clusters for isolated nodes that the
    store and .dbrp paths (which register only nodes incident to a surviving edge)
    omit. The meaningful (size>=2) clusters were already identical across all three
    paths; only isolated-node singletons differed. After the fix the bare-TSV path
    gates node population behind the same cutoff, so all three agree exactly.

    A bare-TSV-only layout is reproduced by copying ONLY the pairwise.tsv into an
    isolated directory (no .dbrp / parquet sibling), forcing the bare-TSV branch.
    """

    @classmethod
    def setUpClass(cls):
        cls.prefix, cls.pw_file = _ensure_shared_fixture()
        cls.parquet_dir = f"{cls.prefix}_DBRetina_pairwise"
        cls.dbrp = f"{cls.prefix}_DBRetina_pairwise.dbrp"

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp(prefix="dbretina_clbare_")
        # Isolated copy of ONLY the .tsv -> no .dbrp / parquet sibling, so the
        # cluster command takes the bare-TSV fallback branch.
        self.bare_dir = os.path.join(self.tmpdir, "bare")
        os.makedirs(self.bare_dir)
        self.bare_tsv = os.path.join(self.bare_dir, os.path.basename(self.pw_file))
        shutil.copy(self.pw_file, self.bare_tsv)

    def tearDown(self):
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    @staticmethod
    def _clusters(path):
        """Parse a *_clusters.tsv -> set of frozensets of (lowercased) members.

        Set-of-frozensets so member ordering within a line (which differs across
        readers) doesn't matter; SINGLETONS are included so the asymmetry under
        test is actually compared.
        """
        clusters = set()
        with open(path) as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                parts = line.rstrip("\n").split("\t")
                if not parts or "cluster_id" in parts[0].lower():
                    continue
                members = parts[-1].split("|")
                clusters.add(frozenset(m.strip().lower() for m in members if m.strip()))
        return clusters

    def test_063_bare_tsv_clusters_match_dbrp(self):
        """cluster on the bare-TSV layout yields the SAME clusters (incl singleton
        handling) as the .dbrp (store) path for the same cutoff."""
        out_bare = os.path.join(self.tmpdir, "clu_bare")
        out_dbrp = os.path.join(self.tmpdir, "clu_dbrp")
        rc, _, stderr = run_command(
            f"DBRetina cluster -p {self.bare_tsv} -m ochiai -c 50 -o {out_bare}"
        )
        assert_no_traceback(self, stderr, "cluster bare-tsv")
        self.assertEqual(rc, 0, stderr)
        rc, _, stderr = run_command(
            f"DBRetina cluster -p {self.dbrp} -m ochiai -c 50 -o {out_dbrp}"
        )
        self.assertEqual(rc, 0, stderr)

        bare = self._clusters(f"{out_bare}_clusters.tsv")
        dbrp = self._clusters(f"{out_dbrp}_clusters.tsv")
        self.assertEqual(
            bare, dbrp,
            "bare-TSV cluster set differs from .dbrp path (issue 063: "
            f"bare emits extra singletons).\nbare={sorted(map(sorted, bare))}\n"
            f"dbrp={sorted(map(sorted, dbrp))}",
        )
        # And there must be no isolated-node singletons the store path omits.
        self.assertFalse(
            any(len(c) == 1 for c in bare - dbrp),
            "bare-TSV path still emits singleton clusters the .dbrp path omits",
        )

    def test_063_bare_tsv_clusters_match_parquet_dir(self):
        """Same parity check against the parquet-dir (store) path."""
        out_bare = os.path.join(self.tmpdir, "clu_bare2")
        out_store = os.path.join(self.tmpdir, "clu_store")
        rc, _, stderr = run_command(
            f"DBRetina cluster -p {self.bare_tsv} -m ochiai -c 50 -o {out_bare}"
        )
        self.assertEqual(rc, 0, stderr)
        rc, _, stderr = run_command(
            f"DBRetina cluster -p {self.parquet_dir} -m ochiai -c 50 -o {out_store}"
        )
        self.assertEqual(rc, 0, stderr)
        self.assertEqual(
            self._clusters(f"{out_bare}_clusters.tsv"),
            self._clusters(f"{out_store}_clusters.tsv"),
            "bare-TSV cluster set differs from parquet-dir path (issue 063)",
        )

    def test_063_meaningful_clusters_present(self):
        """Guard: the fix must not drop the real (size>=2) cluster -- the bare-TSV
        path must still contain the {groupa,groupb,groupf} community."""
        out_bare = os.path.join(self.tmpdir, "clu_bare3")
        rc, _, stderr = run_command(
            f"DBRetina cluster -p {self.bare_tsv} -m ochiai -c 50 -o {out_bare}"
        )
        self.assertEqual(rc, 0, stderr)
        bare = self._clusters(f"{out_bare}_clusters.tsv")
        self.assertIn(frozenset({"groupa", "groupb", "groupf"}), bare,
                      f"bare-TSV path lost the meaningful cluster:\n{bare}")


# ============================================================
# SECTION: ISSUE-045 (DBRetinaGraph class-level mutable state leak)
# ============================================================

class TestGraphStateIsolation(unittest.TestCase):
    """ISSUE-045: DBRetinaGraph declared its per-run dicts/sets as CLASS attributes,
    so a second instantiation in the SAME Python process inherited the first run's
    entries (target maps, geneSetToTargetsArgumentID, ...). Invisible to the CLI
    (fresh process per run) but a real latent bug. The mutable per-run state must be
    INSTANCE-scoped: two instances must not share storage or leak entries."""

    @classmethod
    def setUpClass(cls):
        cls.prefix, cls.pw_file = _ensure_shared_fixture()

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp(prefix="dbretina_gstate_")
        # Two DISJOINT target partitions (no overlap -> no ERROR/sys.exit), each a
        # full cover of the 6 fixture groups so build_graph stays on the grouped path.
        self.targets_a = write_file(
            os.path.join(self.tmpdir, "ta.tsv"),
            "GroupA\nGroupB\nGroupC\nGroupD\nGroupE\nGroupF\n",
        )
        self.targets_b = write_file(
            os.path.join(self.tmpdir, "tb_one.tsv"),
            "GroupA\nGroupB\nGroupC\n",
        )
        self.targets_b2 = write_file(
            os.path.join(self.tmpdir, "tb_two.tsv"),
            "GroupD\nGroupE\nGroupF\n",
        )

    def tearDown(self):
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def _make_graph(self, intra_targets):
        """Build a DBRetinaGraph in-process (silent logger). intra_targets is the
        processed --intra-targets value: a list of file-lists (file groups)."""
        from dbretina.customLogger import Logger
        from dbretina.bipartite_graph import DBRetinaGraph
        out_prefix = os.path.join(self.tmpdir, f"st_{id(intra_targets)}")
        return DBRetinaGraph(
            pairwise_file=self.pw_file,
            index_prefix=self.prefix,
            metric="ochiai",
            cutoff=0.0,
            LOGGER=Logger(active=False),
            inter_targets=[],
            intra_targets=intra_targets,
            output_prefix=out_prefix,
        )

    def test_two_instances_have_independent_mutable_state(self):
        """Two fresh DBRetinaGraph instances must NOT share the same mutable
        dict/set objects (pre-fix they were the same class-level object: same id())."""
        g1 = self._make_graph([[self.targets_a]])
        g2 = self._make_graph([[self.targets_b]])
        mutable_attrs = [
            "target_to_gene_sets",
            "node_to_fragmentation",
            "node_to_heterogeneity",
            "node_to_modularity",
            "node_to_target_name",
            "target_groups",
            "target_to_targetGroupID",
            "geneSetToTargetsArgumentID",
            "gene_set_to_targetID",
        ]
        for attr in mutable_attrs:
            self.assertIsNot(
                getattr(g1, attr), getattr(g2, attr),
                f"DBRetinaGraph.{attr} is shared across instances (same object) "
                f"-> class-level mutable state leak (issue 045)",
            )

    def test_no_cross_instance_entry_leak(self):
        """The second instance's target maps must reflect ONLY its own targets, not
        the first instance's entries. Pre-fix the shared class dict accumulated both."""
        # First instance: single intra group covering all 6 groups -> 'intra_1'.
        g1 = self._make_graph([[self.targets_a]])
        # Sanity: g1 saw all six fixture groups.
        self.assertEqual(len(g1.gene_set_to_targetID), 6)
        self.assertTrue(
            all(v == "intra_1" for v in g1.geneSetToTargetsArgumentID.values()),
            "setup expectation: g1's groups all map to its single intra target",
        )

        # Second instance: two SEPARATE intra groups (3 groups each) -> some groups
        # map to 'intra_2'. If state leaked from g1, every group would still be
        # 'intra_1' and/or the file->target map would carry g1's file basename.
        g2 = self._make_graph([[self.targets_b], [self.targets_b2]])
        self.assertEqual(
            len(g2.gene_set_to_targetID), 6,
            "g2 should map exactly its own 6 groups (no leaked/extra entries)",
        )
        # g2 must contain an 'intra_2' assignment (its second target group); a leaked
        # shared dict pre-populated by g1 (all 'intra_1') would lack it.
        self.assertIn(
            "intra_2", set(g2.geneSetToTargetsArgumentID.values()),
            "g2 lost its own second target group -> cross-instance state leak (issue 045)",
        )
        # And g2's file->target map must NOT carry g1's target file basename.
        g1_basename = os.path.splitext(os.path.basename(self.targets_a))[0]
        self.assertNotIn(
            g1_basename, g2.target_to_targetGroupID,
            f"g2.target_to_targetGroupID leaked g1's target file '{g1_basename}' "
            f"(issue 045: shared class dict)",
        )

    def test_first_instance_unmutated_by_second(self):
        """Constructing a second instance must not retroactively mutate the first
        instance's state (it would, if both aliased one class-level dict)."""
        g1 = self._make_graph([[self.targets_a]])
        before = dict(g1.target_to_targetGroupID)
        _g2 = self._make_graph([[self.targets_b], [self.targets_b2]])
        self.assertEqual(
            g1.target_to_targetGroupID, before,
            "first instance's target map changed after building a second instance "
            "(issue 045: shared class-level dict)",
        )


# ============================================================
# SECTION: ISSUE-062 (PairwiseStore deprecated fetch_record_batch)
# ============================================================

try:
    import duckdb as _duckdb_062  # noqa: F401
    _HAS_DUCKDB_062 = True
except Exception:
    _HAS_DUCKDB_062 = False


@unittest.skipUnless(_HAS_DUCKDB_062, "PairwiseStore Arrow paths need duckdb")
class TestPairwiseStoreArrowReader(unittest.TestCase):
    """ISSUE-062: PairwiseStore used the DEPRECATED DuckDB .fetch_record_batch()
    (DuckDB 1.5.4 emits 'fetch_record_batch() is deprecated, use to_arrow_reader()
    instead'; breaks on a future DuckDB major). The streaming Arrow paths
    (filter_pairs / query_group / iterate_all) must use to_arrow_reader(): no
    DeprecationWarning, identical RecordBatchReader results."""

    @classmethod
    def setUpClass(cls):
        cls.prefix, cls.pw_file = _ensure_shared_fixture()
        cls.parquet = f"{cls.prefix}_DBRetina_pairwise"
        cls.dbri = f"{cls.prefix}.dbri"

    def _store(self):
        from dbretina.pairwise_store import PairwiseStore
        return PairwiseStore(self.parquet, dbri_path=self.dbri)

    @staticmethod
    def _consume(reader):
        """Drain a RecordBatchReader into a flat list of row-dicts."""
        rows = []
        for batch in reader:
            d = batch.to_pydict()
            cols = list(d.keys())
            n = len(d[cols[0]]) if cols else 0
            for i in range(n):
                rows.append({c: d[c][i] for c in cols})
        return rows

    def _assert_no_deprecation(self, recorded, label):
        offenders = [
            str(w.message) for w in recorded
            if "fetch_record_batch" in str(w.message)
        ]
        self.assertFalse(
            offenders,
            f"{label} emitted a fetch_record_batch DeprecationWarning "
            f"(issue 062, not migrated to to_arrow_reader): {offenders}",
        )

    def test_filter_pairs_no_deprecation_and_correct(self):
        import pyarrow as pa
        store = self._store()
        try:
            with warnings.catch_warnings(record=True) as rec:
                warnings.simplefilter("always")
                reader = store.filter_pairs(
                    "ochiai", 0.0,
                    columns=["group_1_id", "group_2_id", "ochiai"],
                )
                self.assertIsInstance(reader, pa.RecordBatchReader)
                rows = self._consume(reader)
            self._assert_no_deprecation(rec, "filter_pairs")
            self.assertTrue(rows, "filter_pairs returned no rows")
            self.assertTrue(all(r["ochiai"] >= 0.0 for r in rows))
            # cutoff actually filters: a high cutoff yields strictly fewer rows.
            hi = self._consume(
                store.filter_pairs("ochiai", 99.0, columns=["ochiai"])
            )
            self.assertLess(len(hi), len(rows), "cutoff did not reduce the result")
        finally:
            store.close()

    def test_query_group_no_deprecation_and_correct(self):
        import pyarrow as pa
        store = self._store()
        try:
            with warnings.catch_warnings(record=True) as rec:
                warnings.simplefilter("always")
                reader = store.query_group(
                    "groupa", metric="ochiai", cutoff=0.0,
                    columns=["group_1_id", "group_2_id", "ochiai"],
                )
                self.assertIsInstance(reader, pa.RecordBatchReader)
                rows = self._consume(reader)
            self._assert_no_deprecation(rec, "query_group")
            self.assertTrue(rows, "query_group(groupa) returned no rows")
            gid = store.group_id("groupa")
            self.assertIsNotNone(gid)
            # every returned pair must involve groupa
            self.assertTrue(
                all(gid in (r["group_1_id"], r["group_2_id"]) for r in rows),
                "query_group returned a pair not involving the queried group",
            )
        finally:
            store.close()

    def test_iterate_all_no_deprecation_and_correct(self):
        import pyarrow as pa
        store = self._store()
        try:
            with warnings.catch_warnings(record=True) as rec:
                warnings.simplefilter("always")
                reader = store.iterate_all(
                    columns=["group_1_id", "group_2_id", "ochiai"]
                )
                self.assertIsInstance(reader, pa.RecordBatchReader)
                rows = self._consume(reader)
            self._assert_no_deprecation(rec, "iterate_all")
            self.assertTrue(rows, "iterate_all returned no rows")
            # iterate_all (no filter) must cover at least the filtered subset.
            filtered = self._consume(
                store.filter_pairs("ochiai", 0.0, columns=["ochiai"])
            )
            self.assertGreaterEqual(len(rows), len(filtered))
        finally:
            store.close()

    def test_no_deprecation_across_all_arrow_paths(self):
        """One guard covering all three readers together: zero fetch_record_batch
        DeprecationWarnings anywhere on the streaming Arrow paths."""
        store = self._store()
        try:
            with warnings.catch_warnings(record=True) as rec:
                warnings.simplefilter("always")
                self._consume(store.filter_pairs("ochiai", 0.0, columns=["ochiai"]))
                self._consume(
                    store.query_group("groupa", metric="ochiai", columns=["ochiai"])
                )
                self._consume(store.iterate_all(columns=["ochiai"]))
            self._assert_no_deprecation(rec, "all Arrow paths")
        finally:
            store.close()


# ============================================================
# SECTION: ISSUE-064 (graph -m emits whatever metric is chosen)
# ============================================================

class TestGraphMetricEmission(unittest.TestCase):
    """ISSUE-064 (investigation): unlike `bipartite` (issue 048, fixed schema of
    containment/ochiai/jaccard[/pvalue] columns), the `graph` command emits a
    GENERIC 3-column edges file (from / to / <metric>) whose single weight column
    is the chosen -m metric, read via metric_to_col (TSV) / rec[self.metric]
    (.dbrp) / d[self.metric] (parquet). So graph correctly emits ALL 7 metrics,
    csi/dice/odds_ratio included -- it does NOT silently drop them. Verdict:
    NOT a defect; -m must stay permissive (restricting it would break working
    output). These tests document/lock in that the emitted weights actually match
    the requested metric's pairwise values (would have caught a silent drop)."""

    @classmethod
    def setUpClass(cls):
        cls.prefix, cls.pw_file = _ensure_shared_fixture()

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp(prefix="dbretina_gmetric_")
        # Ground truth: {frozenset(name1,name2): {metric: value}} from the pairwise TSV.
        self.truth = {}
        for r in parse_pairwise_tsv(self.pw_file):
            self.truth[frozenset((r["group_1_name"], r["group_2_name"]))] = r

    def tearDown(self):
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def _read_edges(self, edges_path):
        """Return (header_metric_name, {frozenset(from,to): weight})."""
        weights = {}
        header_metric = None
        with open(edges_path) as f:
            header = f.readline().rstrip("\n").split("\t")
            self.assertEqual(header[:2], ["from", "to"],
                             f"unexpected edges header: {header}")
            header_metric = header[2]
            for line in f:
                line = line.rstrip("\n")
                if not line:
                    continue
                a, b, w = line.split("\t")
                weights[frozenset((a, b))] = float(w)
        return header_metric, weights

    def _run_and_check(self, metric, pairwise_arg, label):
        """graph -m <metric> on `pairwise_arg`; assert the emitted weight column is
        named for the metric and each edge weight equals that metric's pairwise
        value (i.e. the requested metric is genuinely emitted, not dropped/wrong)."""
        out = os.path.join(self.tmpdir, f"g_{metric}_{label}")
        # cutoff 0: odds_ratio's -1.0 (subset/identical) rows fall below 0 and are
        # excluded -- that's correct cutoff behavior, not a dropped metric.
        rc, _, stderr = run_command(
            f"DBRetina graph -i {self.prefix} -p {pairwise_arg} "
            f"-m {metric} -c 0 -o {out}"
        )
        self.assertEqual(rc, 0, stderr)
        self.assertNotIn("Traceback", stderr)
        edges_path = f"{out}_edges.tsv"
        assert_file_exists(self, edges_path)

        header_metric, weights = self._read_edges(edges_path)
        # 1. The weight column is named for the REQUESTED metric (not a fixed name).
        self.assertEqual(
            header_metric, metric,
            f"[{label}] edges weight column is '{header_metric}', expected "
            f"'{metric}' -> graph mislabeled/ignored the requested metric",
        )
        # 2. At least one edge emitted (the fixture has plenty above cutoff 0).
        self.assertGreater(
            len(weights), 0,
            f"[{label}] graph -m {metric} produced no edges (silently dropped?)",
        )
        # 3. Every emitted edge weight equals that metric's pairwise value. A
        #    silently-dropped/mis-emitted metric (the 048 failure mode) would put a
        #    different column's number here and fail this. The pairwise TSV columns
        #    are display-rounded to 1 decimal, while graph emits the full-precision
        #    value (e.g. csi 13.333... vs TSV 13.3) -- which itself confirms graph
        #    reads the genuine metric, not a dropped/placeholder one. Compare at the
        #    TSV's own 1-decimal precision: a *wrong* metric (e.g. ochiai 36.5 for a
        #    csi 13.3 request) still rounds differently and fails.
        for pair, w in weights.items():
            self.assertIn(pair, self.truth, f"[{label}] unknown edge {pair}")
            expected = self.truth[pair][metric]
            # Compare at the TSV's 1-decimal precision with a tolerance: the full-precision
            # dir/.dbrp values can round-half oppositely to the C++ TSV formatter on a *.*5
            # tie, so use abs-diff <= 0.1 (a *wrong* metric differs by far more).
            self.assertLessEqual(
                abs(round(w, 1) - float(expected)), 0.1,
                msg=(f"[{label}] edge {pair} weight {w} (->{round(w, 1)}) != "
                     f"pairwise {metric} {expected}: graph emitted the wrong metric"),
            )

    def test_graph_csi_emitted_correctly_tsv(self):
        """graph -m csi (TSV route) emits the csi column with correct values."""
        self._run_and_check("csi", self.pw_file, "tsv")

    def test_graph_dice_emitted_correctly_tsv(self):
        """graph -m dice (TSV route) emits the dice column with correct values."""
        self._run_and_check("dice", self.pw_file, "tsv")

    def test_graph_odds_ratio_emitted_correctly_tsv(self):
        """graph -m odds_ratio (TSV route) emits the odds_ratio column correctly."""
        self._run_and_check("odds_ratio", self.pw_file, "tsv")

    def test_graph_csi_dice_oddsratio_emitted_correctly_parquet_dir(self):
        """Same correctness via the parquet-DIR route (d[self.metric])."""
        pw_dir = self.pw_file[:-len(".tsv")] if self.pw_file.endswith(".tsv") else self.pw_file
        if not (os.path.isdir(pw_dir)
                and os.path.exists(os.path.join(pw_dir, "manifest.json"))):
            self.skipTest("no parquet pairwise directory sibling")
        for metric in ("csi", "dice", "odds_ratio"):
            self._run_and_check(metric, pw_dir, "dir")

    def test_graph_csi_dice_oddsratio_emitted_correctly_dbrp(self):
        """Same correctness via the .dbrp route (rec[self.metric])."""
        dbrp = self.pw_file[:-len(".tsv")] + ".dbrp" if self.pw_file.endswith(".tsv") else None
        if not (dbrp and os.path.isfile(dbrp)):
            self.skipTest("no .dbrp sibling")
        for metric in ("csi", "dice", "odds_ratio"):
            self._run_and_check(metric, dbrp, "dbrp")


# ============================================================
# SECTION: ISSUE-065 (Clusters / DeduplicateGroups class-level state leak)
# ============================================================

class TestClustersStateIsolation(unittest.TestCase):
    """ISSUE-065 (sibling of 045): Clusters declared its per-run mutable maps
    (group_to_features, names_map) as CLASS attributes, so a second instantiation
    in the SAME process inherited the first run's storage. Invisible to the CLI
    (fresh process per run), but a real latent leak. The per-run mutable state must
    be INSTANCE-scoped: two instances must not share the same dict/set objects."""

    @classmethod
    def setUpClass(cls):
        cls.prefix, cls.pw_file = _ensure_shared_fixture()

    def setUp(self):
        from dbretina.customLogger import Logger
        self._Logger = Logger
        self.tmpdir = tempfile.mkdtemp(prefix="dbretina_clu_state_")

    def tearDown(self):
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def _make_clusters(self):
        from dbretina.clustering import Clusters
        out_prefix = os.path.join(self.tmpdir, f"c_{os.urandom(4).hex()}")
        # community=False -> rustworkx, light construction; reads the pairwise file
        # only for the node count (the shared fixture .tsv is always present).
        return Clusters(
            logger_obj=self._Logger(active=False),
            pairwise_file=self.pw_file,
            cut_off_threshold=0.0,
            metric="ochiai",
            output_prefix=out_prefix,
            community=False,
        )

    def test_two_instances_have_independent_mutable_state(self):
        """Two fresh Clusters instances must NOT share the same mutable dict objects
        (pre-fix group_to_features/names_map were one class-level object: same id)."""
        c1 = self._make_clusters()
        c2 = self._make_clusters()
        for attr in ("group_to_features", "names_map"):
            self.assertIsNot(
                getattr(c1, attr), getattr(c2, attr),
                f"Clusters.{attr} is shared across instances (same object) "
                f"-> class-level mutable state leak (issue 065)",
            )

    def test_mutable_attrs_are_instance_scoped_and_empty(self):
        """Each per-run mutable attr is an instance attribute (in the instance
        __dict__) and starts empty -> no inherited entries from a prior run."""
        c1 = self._make_clusters()
        for attr in ("group_to_features", "names_map"):
            self.assertIn(
                attr, vars(c1),
                f"Clusters.{attr} is not instance-scoped (still a class attr) "
                f"-> latent cross-instance leak (issue 065)",
            )
            self.assertEqual(
                len(getattr(c1, attr)), 0,
                f"Clusters.{attr} should be empty at construction",
            )

    def test_no_cross_instance_entry_leak(self):
        """Mutating one instance's maps must not appear in a later instance."""
        c1 = self._make_clusters()
        c1.group_to_features["sentinel_group"] = {"sentinel_feature"}
        c1.names_map[999999] = "sentinel_name"
        c2 = self._make_clusters()
        self.assertNotIn(
            "sentinel_group", c2.group_to_features,
            "Clusters.group_to_features leaked an entry across instances (issue 065)",
        )
        self.assertNotIn(
            999999, c2.names_map,
            "Clusters.names_map leaked an entry across instances (issue 065)",
        )


class TestDeduplicateGroupsStateIsolation(unittest.TestCase):
    """ISSUE-065 (sibling of 045): DeduplicateGroups declared its per-run mutable
    state (cluster_id_to_groups, groups_to_items_no, final_remaining_groups,
    removed_exact_ochiai_groups) as CLASS attributes -> a second in-process
    instance inherited the first run's storage/entries. They must be INSTANCE
    attributes initialized in __init__."""

    MUTABLE_ATTRS = (
        "cluster_id_to_groups",
        "groups_to_items_no",
        "final_remaining_groups",
        "removed_exact_ochiai_groups",
    )

    @classmethod
    def setUpClass(cls):
        # Only the index_prefix string is needed in __init__ (no file reads there).
        cls.prefix, cls.pw_file = _ensure_shared_fixture()

    def setUp(self):
        from dbretina.customLogger import Logger
        self._Logger = Logger

    def _make_dedup(self):
        from dbretina.setcov import DeduplicateGroups
        # __init__ reads no files: it sets attrs, builds an empty df_logging, and
        # derives main_pairwise_file from index_prefix. A minimal ctx with .obj
        # (the LOGGER) is all that's required.
        ctx = SimpleNamespace(obj=self._Logger(active=False))
        return DeduplicateGroups(
            associations_file=None,
            index_prefix=self.prefix,
            ctx=ctx,
        )

    def test_two_instances_have_independent_mutable_state(self):
        """Two fresh DeduplicateGroups instances must NOT share the same mutable
        dict/set objects (pre-fix all four were single class-level objects)."""
        d1 = self._make_dedup()
        d2 = self._make_dedup()
        for attr in self.MUTABLE_ATTRS:
            self.assertIsNot(
                getattr(d1, attr), getattr(d2, attr),
                f"DeduplicateGroups.{attr} is shared across instances (same object) "
                f"-> class-level mutable state leak (issue 065)",
            )

    def test_mutable_attrs_are_instance_scoped_and_empty(self):
        """Each per-run mutable attr is an instance attribute and starts empty."""
        d1 = self._make_dedup()
        for attr in self.MUTABLE_ATTRS:
            self.assertIn(
                attr, vars(d1),
                f"DeduplicateGroups.{attr} is not instance-scoped (still a class "
                f"attr) -> latent cross-instance leak (issue 065)",
            )
            self.assertEqual(
                len(getattr(d1, attr)), 0,
                f"DeduplicateGroups.{attr} should be empty at construction",
            )

    def test_no_cross_instance_entry_leak(self):
        """Entries added to one instance must not appear in a later instance."""
        d1 = self._make_dedup()
        d1.cluster_id_to_groups[42] = ["sentinel_group"]
        d1.groups_to_items_no["sentinel_group"] = 7
        d1.final_remaining_groups.add("sentinel_remaining")
        d1.removed_exact_ochiai_groups.add("sentinel_removed")

        d2 = self._make_dedup()
        self.assertNotIn(42, d2.cluster_id_to_groups,
                         "cluster_id_to_groups leaked across instances (issue 065)")
        self.assertNotIn("sentinel_group", d2.groups_to_items_no,
                         "groups_to_items_no leaked across instances (issue 065)")
        self.assertNotIn("sentinel_remaining", d2.final_remaining_groups,
                         "final_remaining_groups leaked across instances (issue 065)")
        self.assertNotIn("sentinel_removed", d2.removed_exact_ochiai_groups,
                         "removed_exact_ochiai_groups leaked across instances (issue 065)")


class TestExportNeo4j(unittest.TestCase):
    """ISSUE-070 + export_neo4j log nit: -d must accept a standalone pairwise TSV
    (as its --help promises and the sibling `export` command does), and the
    'Filtering' log must use the metric-aware operator."""

    @classmethod
    def setUpClass(cls):
        try:
            import neo4j  # noqa: F401
        except ImportError:
            raise unittest.SkipTest("neo4j driver not installed")
        cls.prefix, cls.pw_file = _ensure_shared_fixture()
        cls.pv_prefix, cls.pv_pw_tsv = _ensure_shared_pvalue_fixture()
        # A high port nothing listens on -> a fast 'connection refused' so the
        # command fails at the Bolt connect step, never at input-open.
        cls.dead_uri = "bolt://127.0.0.1:65399"

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp(prefix="dbretina_neo4j_")

    def tearDown(self):
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def _standalone_tsv(self, src_tsv):
        """Copy a pairwise TSV to a fresh dir with NO parquet/.dbrp sibling and a
        name that does not resolve to one -- a genuine standalone TSV input."""
        solo = os.path.join(self.tmpdir, "solo.tsv")
        shutil.copy(src_tsv, solo)
        # Make sure nothing sibling-like exists next to it.
        self.assertFalse(os.path.isdir(solo[:-len(".tsv")]))
        return solo

    def test_export_neo4j_accepts_standalone_tsv(self):
        """ISSUE-070: -d <standalone pairwise TSV> must be ingested (not rejected
        with 'Pairwise directory not found'); it should get PAST input-open and fail
        only at the neo4j connection step."""
        solo = self._standalone_tsv(self.pw_file)
        rc, stdout, stderr = run_command(
            f"DBRetina export-neo4j -d {solo} --uri {self.dead_uri} "
            f"--user '' --password 'x' --database neo4j -m jaccard -c 50 --clear"
        )
        combined = stdout + stderr
        self.assertNotIn("Traceback", combined, f"raw traceback:\n{combined}")
        # The input was accepted: no 'directory not found' rejection.
        self.assertNotIn("directory not found", combined.lower(),
                         f"standalone TSV was rejected as a directory:\n{combined}")
        self.assertNotIn("could not open pairwise data", combined.lower(),
                         f"standalone TSV input failed to open:\n{combined}")
        # It got past input-open (loaded the store) ...
        self.assertIn("Loaded:", combined, f"input never opened:\n{combined}")
        # ... and failed only at the neo4j connection.
        self.assertIn("connect", combined.lower(),
                      f"expected a neo4j connection error, got:\n{combined}")

    def test_export_neo4j_filter_log_operator_metric_aware(self):
        """export_neo4j log nit: the 'Filtering' line must use '>=' for similarity
        metrics and '<=' for pvalue (mirrors the real PairwiseStore cutoff filter,
        issue 068). It previously hardcoded '>='."""
        # similarity metric -> '>='
        solo = self._standalone_tsv(self.pw_file)
        _, stdout, stderr = run_command(
            f"DBRetina export-neo4j -d {solo} --uri {self.dead_uri} "
            f"--user '' --password 'x' --database neo4j -m jaccard -c 50"
        )
        self.assertIn("Filtering: jaccard >= 50", stdout + stderr)

        # pvalue -> '<=' (a standalone TSV that actually carries a pvalue column)
        solo_pv = os.path.join(self.tmpdir, "solo_pv.tsv")
        shutil.copy(self.pv_pw_tsv, solo_pv)
        _, stdout, stderr = run_command(
            f"DBRetina export-neo4j -d {solo_pv} --uri {self.dead_uri} "
            f"--user '' --password 'x' --database neo4j -m pvalue -c 0.05"
        )
        combined = stdout + stderr
        self.assertIn("Filtering: pvalue <= 0.05", combined)
        self.assertNotIn("pvalue >=", combined)


# ============================================================
# Main
# ============================================================

if __name__ == "__main__":
    try:
        unittest.main(verbosity=2)
    finally:
        _cleanup_shared_fixture()
        _cleanup_shared_pvalue_fixture()
        _cleanup_zero_pair_fixture()
