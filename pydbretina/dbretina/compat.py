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
