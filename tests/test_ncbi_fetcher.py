"""Tests for NCBI fetcher country filtering."""

import pytest

from ncbi_fetcher import (
    _country_matches_filter,
    _country_search_term,
    build_search_query,
)


def test_build_query_uses_title_not_country_field():
    q = build_search_query("gag", 2016, 2026, "United States")
    assert '"USA"[Title]' in q
    assert "[Country]" not in q


def test_build_query_any_country():
    q = build_search_query("gag", 2016, 2026, "Any")
    assert "[Title]" not in q
    assert "[Country]" not in q


def test_country_search_term_usa_alias():
    assert _country_search_term("United States") == "USA"
    assert _country_search_term("Uganda") == "Uganda"


@pytest.mark.parametrize(
    "geo,requested,expected",
    [
        ("USA", "United States", True),
        ("Uganda", "United States", False),
        ("Uganda", "Uganda", True),
        ("United Kingdom", "United Kingdom", True),
        ("England", "United Kingdom", True),
        ("Botswana", "United States", False),
        ("Unknown", "United States", False),
    ],
)
def test_country_matches_filter(geo, requested, expected):
    assert _country_matches_filter(geo, requested) is expected
