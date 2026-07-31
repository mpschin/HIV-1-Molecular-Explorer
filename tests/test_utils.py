"""Unit tests for HIV-1 Molecular Explorer data integrity utilities."""

import pandas as pd
import pytest

from utils import is_valid_iupac_sequence, run_data_integrity_check


@pytest.fixture
def valid_df():
    return pd.DataFrame(
        {
            "genbank_id": ["AB123", "CD456"],
            "year": ["2020", "2021"],
            "title": ["HIV-1 isolate A", "HIV-1 isolate B"],
            "country": ["Uganda", "Kenya"],
            "sequence": ["ATCGATCGNN", "GGCCAAUUTT"],
        }
    )


def test_valid_iupac():
    assert is_valid_iupac_sequence("ATCGN")
    assert is_valid_iupac_sequence("atcgn")
    assert not is_valid_iupac_sequence("ATCGX")
    assert not is_valid_iupac_sequence("")


def test_integrity_check_passes(valid_df):
    result = run_data_integrity_check(valid_df)
    assert result.passed
    assert any("passed" in m.lower() for m in result.messages)


def test_integrity_check_empty():
    result = run_data_integrity_check(pd.DataFrame())
    assert not result.passed


def test_integrity_check_missing_columns(valid_df):
    df = valid_df.drop(columns=["country"])
    result = run_data_integrity_check(df)
    assert not result.passed
    assert any("country" in m for m in result.messages)


def test_integrity_check_invalid_sequence(valid_df):
    df = valid_df.copy()
    df.loc[0, "sequence"] = "ATCG123"
    result = run_data_integrity_check(df)
    assert not result.passed
