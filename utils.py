"""Data integrity and validation utilities for HIV-1 Molecular Explorer."""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import List

import pandas as pd

REQUIRED_COLUMNS = ["genbank_id", "year", "title", "country", "sequence"]

# Standard IUPAC nucleotide ambiguity codes for DNA/RNA.
IUPAC_DNA_PATTERN = re.compile(
    r"^[ATCGUNRYKMSWBDHVatcgunrykmswbdhv\-\.]+$",
    re.IGNORECASE,
)


@dataclass
class IntegrityCheckResult:
    """Structured result from a data integrity check."""

    passed: bool
    messages: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        return {
            "passed": self.passed,
            "messages": self.messages,
            "warnings": self.warnings,
        }


def is_valid_iupac_sequence(sequence: str) -> bool:
    """Return True if *sequence* contains only valid IUPAC DNA characters."""
    if not sequence or not isinstance(sequence, str):
        return False
    cleaned = sequence.strip().replace(" ", "").replace("\n", "")
    if not cleaned:
        return False
    return bool(IUPAC_DNA_PATTERN.match(cleaned))


def run_data_integrity_check(df: pd.DataFrame) -> IntegrityCheckResult:
    """
    Verify that a fetched sequence DataFrame meets BioVibe requirements.

    Checks:
      1. DataFrame is non-empty.
      2. Required columns exist.
      3. Sequences contain valid standard IUPAC DNA characters.
    """
    result = IntegrityCheckResult(passed=True)

    if df is None:
        result.passed = False
        result.messages.append("DataFrame is None — no data was loaded.")
        return result

    if not isinstance(df, pd.DataFrame):
        result.passed = False
        result.messages.append("Input is not a pandas DataFrame.")
        return result

    if df.empty:
        result.passed = False
        result.messages.append("DataFrame is empty — fetch sequences before running checks.")
        return result

    result.messages.append(f"Row count: {len(df)}")

    missing_cols = [col for col in REQUIRED_COLUMNS if col not in df.columns]
    if missing_cols:
        result.passed = False
        result.messages.append(f"Missing required columns: {', '.join(missing_cols)}")
        return result

    result.messages.append("All required columns present.")

    null_counts = df[REQUIRED_COLUMNS].isnull().sum()
    for col, count in null_counts.items():
        if count > 0:
            result.warnings.append(f"Column '{col}' has {count} null value(s).")

    empty_ids = df["genbank_id"].astype(str).str.strip().eq("").sum()
    if empty_ids > 0:
        result.warnings.append(f"{empty_ids} row(s) have empty GenBank IDs.")

    invalid_rows: List[int] = []
    for idx, seq in enumerate(df["sequence"].astype(str)):
        if not is_valid_iupac_sequence(seq):
            invalid_rows.append(idx)

    if invalid_rows:
        result.passed = False
        preview = invalid_rows[:10]
        suffix = "..." if len(invalid_rows) > 10 else ""
        result.messages.append(
            f"Invalid IUPAC sequences at row index(es): {preview}{suffix} "
            f"({len(invalid_rows)} total)."
        )
    else:
        result.messages.append("All sequences contain valid IUPAC DNA characters.")

    dup_ids = df["genbank_id"].duplicated().sum()
    if dup_ids > 0:
        result.warnings.append(f"{dup_ids} duplicate GenBank ID(s) detected.")

    short_seqs = (df["sequence"].astype(str).str.len() < 50).sum()
    if short_seqs > 0:
        result.warnings.append(
            f"{short_seqs} sequence(s) are shorter than 50 nt — may be partial fragments."
        )

    if result.passed and not result.warnings:
        result.messages.append("Data integrity check passed with no warnings.")
    elif result.passed:
        result.messages.append("Data integrity check passed with warnings (see below).")

    return result
