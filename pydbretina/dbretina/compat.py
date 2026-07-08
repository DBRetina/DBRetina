"""Format auto-detection for DBRetina pairwise data.

Given a TSV path (legacy CLI interface), auto-detect and prefer Parquet
directory if available, falling back to existing TSV/.dbrp behavior.
"""

import os
import pathlib
from typing import Optional

from .pairwise_store import PairwiseStore


def open_pairwise(pairwise_path: str) -> Optional[PairwiseStore]:
    """Try to open a PairwiseStore from a pairwise TSV or directory path.

    Auto-detection logic:
      1. If path is a directory with manifest.json → use directly
      2. If path is a TSV file → look for sibling ``_pairwise/`` directory
         (e.g. ``foo_DBRetina_pairwise.tsv`` → ``foo_DBRetina_pairwise/``,
          or ``foo_pairwise.tsv`` → strip ``_DBRetina_pairwise`` → ``foo_pairwise/``)
      3. Returns None if no Parquet directory found

    Args:
        pairwise_path: Path to pairwise TSV file or Parquet directory.

    Returns:
        PairwiseStore if Parquet data found, None otherwise.
    """
    p = pathlib.Path(pairwise_path)

    # Case 1: direct Parquet directory
    if p.is_dir() and (p / "manifest.json").exists():
        return PairwiseStore(str(p))

    # Case 2: TSV file → derive Parquet directory
    if p.is_file() or p.suffix == ".tsv":
        # foo_DBRetina_pairwise.tsv → foo_pairwise/
        stem = p.stem  # e.g. "foo_DBRetina_pairwise"
        parent = p.parent

        # Try: replace _DBRetina_pairwise with _pairwise
        if "_DBRetina_pairwise" in stem:
            parquet_name = stem.replace("_DBRetina_pairwise", "_pairwise")
            parquet_dir = parent / parquet_name
            if parquet_dir.is_dir() and (parquet_dir / "manifest.json").exists():
                return PairwiseStore(str(parquet_dir))

        # Try: just strip .tsv extension as directory name
        parquet_dir = parent / stem
        if parquet_dir.is_dir() and (parquet_dir / "manifest.json").exists():
            return PairwiseStore(str(parquet_dir))

    return None


def resolve_dbrp_path(pairwise_path: str) -> Optional[str]:
    """Resolve the canonical sibling ``*_DBRetina_pairwise.dbrp`` for a -p input.

    -p may be the pairwise TSV, the parquet directory, or the .dbrp binary, all
    sharing the base ``<prefix>_DBRetina_pairwise``. The fallback (non-store)
    readers in modularity/dedup derived the .dbrp via
    ``pairwise_file.replace('_pairwise.tsv', '_pairwise.dbrp')``, a no-op for the
    directory/.dbrp forms; the directory path then leaked into the .dbrp binary
    reader and crashed with "Invalid .dbrp file (bad magic bytes)" (issue 046).

    Strip a trailing ``.tsv``/``.dbrp`` and append ``.dbrp``. Returns the path
    only if it is an existing *file* (never a directory), else None so the caller
    can fall through to TSV reading or emit a clean error.
    """
    base = pairwise_path
    for ext in (".tsv", ".dbrp"):
        if base.endswith(ext):
            base = base[: -len(ext)]
            break
    dbrp = base + ".dbrp"
    if os.path.isfile(dbrp):
        return dbrp
    return None


# .dbrp header layout (see include/DBRetinaPairwise.hpp / src/DBRetinaPairwise.cpp::begin_write):
#   4  magic "DBRP" | 4 version | 8 toc_offset | 8 num_pairs | 4 num_groups
#   then a 1-byte metric_flags bitfield at offset 28, where bit 6 (0x40) == PVALUE.
_DBRP_MAGIC = b"DBRP"
_DBRP_METRIC_FLAGS_OFFSET = 28
_DBRP_PVALUE_BIT = 0x40  # DBRPMetric::PVALUE == 6


def _dbrp_has_pvalue(dbrp_path: str) -> bool:
    """Read the pvalue bit from a .dbrp header (O(1), no full materialization).

    Reads the fixed-offset ``metric_flags`` byte rather than scanning records, so
    it is cheap on large pairwise files. Falls back to probing the decoded records
    (``'pvalue'`` is only present per-record when the metric was computed) if the
    header is unexpectedly short or lacks the magic bytes.
    """
    try:
        with open(dbrp_path, "rb") as fh:
            head = fh.read(_DBRP_METRIC_FLAGS_OFFSET + 1)
        if len(head) > _DBRP_METRIC_FLAGS_OFFSET and head[:4] == _DBRP_MAGIC:
            return bool(head[_DBRP_METRIC_FLAGS_OFFSET] & _DBRP_PVALUE_BIT)
    except OSError:
        pass
    # Fallback: ask the C++ reader. Records carry a 'pvalue' key only when stored.
    import _dbretina_internal as dbretina_internal

    records = dbretina_internal.dbrp_filter_pairs(dbrp_path, 0, 0.0)
    return len(records) > 0 and "pvalue" in records[0]


def pairwise_has_pvalue(pairwise_path: str) -> bool:
    """Format-aware test for whether a pairwise input carries a pvalue column.

    ``-p`` may be the pairwise TSV, the parquet directory, or the .dbrp binary.
    Each command that supports ``-m pvalue`` must reject a pvalue request on a
    dataset that lacks it with the SAME clean error, on EVERY input form. The old
    per-module ``check_if_there_is_a_pvalue`` open()ed the path as text, which
    crashed with ``IsADirectoryError`` on a parquet directory and
    ``UnicodeDecodeError`` on a .dbrp before the format was resolved (issues
    043/053). This helper resolves the format first:

      - parquet directory / .tsv with a sibling parquet store -> ``PairwiseStore.has_pvalue``
      - .dbrp (or a -p form whose canonical .dbrp sibling exists) -> the .dbrp metric flags
      - plain .tsv -> scan the first non-comment line for a ``pvalue`` column

    Returns True iff a pvalue column is present. Callers emit the canonical
    ``"pvalue not found in pairwise file!"`` error when this is False.
    """
    # 1. Parquet store (the directory form, or a .tsv with a sibling parquet dir).
    try:
        store = open_pairwise(pairwise_path)
    except Exception:
        store = None
    if store is not None:
        try:
            return store.has_pvalue
        finally:
            store.close()

    # 2. .dbrp binary (a .dbrp passed directly, or a -p form whose .dbrp sibling exists).
    dbrp_path = resolve_dbrp_path(pairwise_path)
    if dbrp_path is not None:
        return _dbrp_has_pvalue(dbrp_path)

    # 3. Plain text TSV: scan the first non-comment (header) line.
    #
    # A directory with neither a manifest.json (step 1) nor a sibling .dbrp
    # (step 2) reaches here; opening it as text would raise IsADirectoryError
    # (issue 052). It cannot carry a pvalue column, so report False cleanly and
    # let the caller fall through to its own format handling / error.
    if os.path.isdir(pairwise_path):
        return False
    with open(pairwise_path) as fh:
        for line in fh:
            if not line.startswith("#"):
                return "pvalue" in line
    return False
