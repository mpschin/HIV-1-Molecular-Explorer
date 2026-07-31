"""NCBI Entrez API query, parsing, and CSV export for HIV-1 sequences."""

from __future__ import annotations

import io
import re
import time
import xml.etree.ElementTree as ET
from typing import Optional, Tuple

import pandas as pd
from Bio import Entrez, SeqIO

# NCBI rate-limit: max 3 requests per second without an API key.
_ENTREZ_DELAY = 0.4

COUNTRY_PATTERNS = [
    re.compile(r"\b(?:country|isolation\s+country|geo_loc_name)\s*[:=]\s*([A-Za-z\s\-]+)", re.I),
    re.compile(r"\bfrom\s+([A-Z][a-z]+(?:\s+[A-Z][a-z]+)*)\b"),
]

YEAR_PATTERNS = [
    re.compile(r"\b(19|20)\d{2}\b"),
    re.compile(r"\b(\d{4})/(\d{2})/(\d{2})\b"),
]


# NCBI [Country] indexes submitter/author country, NOT sample isolation site.
# Use [Title] (e.g. "from USA") and post-fetch geo_loc_name filtering instead.
COUNTRY_SEARCH_TERMS: dict[str, str] = {
    "United States": "USA",
    "United Kingdom": "United Kingdom",
}

# Accepted geo_loc_name values per user-facing country label (lowercase).
COUNTRY_GEO_ALIASES: dict[str, set[str]] = {
    "United States": {"usa", "us", "united states", "u.s.a.", "u.s."},
    "United Kingdom": {"united kingdom", "uk", "u.k.", "england", "scotland", "wales", "northern ireland"},
    "South Africa": {"south africa"},
    "India": {"india"},
    "Brazil": {"brazil"},
    "China": {"china"},
    "Thailand": {"thailand"},
    "Uganda": {"uganda"},
    "Kenya": {"kenya"},
    "France": {"france"},
    "Germany": {"germany"},
    "Australia": {"australia"},
}


class NCBIFetchError(Exception):
    """Raised when NCBI Entrez operations fail."""


def configure_entrez(email: str, api_key: Optional[str] = None) -> None:
    """Set Entrez credentials required by NCBI."""
    if not email or not email.strip():
        raise NCBIFetchError("A valid email address is required for NCBI Entrez.")
    Entrez.email = email.strip()
    Entrez.api_key = api_key


def _normalize_geo(value: str) -> str:
    """Normalize a geo_loc_name / country string for comparison."""
    return value.strip().lower().split(":")[0].strip()


def _country_search_term(country: str) -> str:
    """Map a display country name to the best NCBI Title search term."""
    country = country.strip()
    return COUNTRY_SEARCH_TERMS.get(country, country)


def _country_matches_filter(record_country: str, requested_country: str) -> bool:
    """Return True if parsed isolation country matches the user's filter."""
    if not requested_country or requested_country.lower() in ("any", "all", ""):
        return True

    geo = _normalize_geo(record_country)
    if geo in ("unknown", ""):
        return False

    aliases = COUNTRY_GEO_ALIASES.get(requested_country)
    if aliases:
        return geo in aliases

    # Custom / unlisted country: case-insensitive substring match.
    requested = _normalize_geo(requested_country)
    return requested in geo or geo in requested


def build_search_query(
    gene: str,
    year_start: int,
    year_end: int,
    country: str = "Any",
) -> str:
    """
    Build an Entrez nucleotide search query for HIV-1 sequences.

    Template: HIV-1[Organism] AND <gene>[Gene] AND <year_range>[PDAT]
    """
    gene = gene.strip()
    if not gene:
        raise NCBIFetchError("Target gene cannot be empty.")

    query_parts = [
        "HIV-1[Organism]",
        f"{gene}[Gene]",
        f"{year_start}:{year_end}[PDAT]",
    ]

    country = country.strip()
    if country and country.lower() not in ("any", "all", ""):
        search_term = _country_search_term(country)
        query_parts.append(f'"{search_term}"[Title]')

    return " AND ".join(query_parts)


def _parse_year_from_text(text: str) -> str:
    """Extract a 4-digit year from free text, or return 'Unknown'."""
    if not text:
        return "Unknown"
    for pattern in YEAR_PATTERNS:
        match = pattern.search(text)
        if match:
            year_str = match.group(0)
            digits = re.search(r"\b(19|20)\d{2}\b", year_str)
            if digits:
                return digits.group(0)
    return "Unknown"


def _parse_country_from_text(text: str) -> str:
    """Extract country from free text, or return 'Unknown'."""
    if not text:
        return "Unknown"
    for pattern in COUNTRY_PATTERNS:
        match = pattern.search(text)
        if match:
            return match.group(1).strip().title()
    return "Unknown"


def _extract_geo_country(record) -> str:
    """Pull sample isolation country from GenBank source feature qualifiers."""
    for feature in record.features:
        if feature.type == "source":
            quals = dict(feature.qualifiers)
            for key in ("geo_loc_name", "country"):
                if key in quals:
                    raw = quals[key][0]
                    country = raw.split(":")[0].strip()
                    if country:
                        return country.title()
    return "Unknown"


def _extract_year(record) -> str:
    """Extract isolation year from record metadata."""
    if hasattr(record, "annotations"):
        date = record.annotations.get("date") or record.annotations.get("submission_date")
        if date:
            year = _parse_year_from_text(str(date))
            if year != "Unknown":
                return year

    for feature in record.features:
        if feature.type == "source":
            for key in ("collection_date", "isolation_date", "date"):
                if key in feature.qualifiers:
                    year = _parse_year_from_text(feature.qualifiers[key][0])
                    if year != "Unknown":
                        return year

    title = record.description if hasattr(record, "description") else ""
    return _parse_year_from_text(title)


def _safe_efetch(
    id_list: list,
    db: str = "nucleotide",
    rettype: str = "gb",
    retmode: str = "text",
    max_retries: int = 3,
) -> str:
    """Fetch records from NCBI with retry logic for rate limits."""
    last_error: Optional[Exception] = None
    for attempt in range(max_retries):
        try:
            time.sleep(_ENTREZ_DELAY)
            handle = Entrez.efetch(
                db=db,
                id=id_list,
                rettype=rettype,
                retmode=retmode,
            )
            data = handle.read()
            handle.close()
            if isinstance(data, bytes):
                return data.decode("utf-8", errors="replace")
            return data
        except Exception as exc:
            last_error = exc
            err_msg = str(exc).lower()
            if "429" in err_msg or "rate" in err_msg or "too many" in err_msg:
                wait = _ENTREZ_DELAY * (2 ** (attempt + 1))
                time.sleep(wait)
                continue
            raise NCBIFetchError(f"NCBI efetch failed: {exc}") from exc

    raise NCBIFetchError(
        f"NCBI efetch failed after {max_retries} retries: {last_error}"
    )


def _safe_esearch(
    query: str,
    db: str = "nucleotide",
    retmax: int = 50,
    max_retries: int = 3,
) -> Tuple[list, int]:
    """Search NCBI with retry logic; returns (id_list, total_count)."""
    last_error: Optional[Exception] = None
    for attempt in range(max_retries):
        try:
            time.sleep(_ENTREZ_DELAY)
            handle = Entrez.esearch(
                db=db,
                term=query,
                retmax=retmax,
                sort="relevance",
            )
            record = Entrez.read(handle)
            handle.close()
            id_list = record.get("IdList", [])
            total = int(record.get("Count", 0))
            return id_list, total
        except Exception as exc:
            last_error = exc
            err_msg = str(exc).lower()
            if "429" in err_msg or "rate" in err_msg or "too many" in err_msg:
                wait = _ENTREZ_DELAY * (2 ** (attempt + 1))
                time.sleep(wait)
                continue
            raise NCBIFetchError(f"NCBI esearch failed: {exc}") from exc

    raise NCBIFetchError(
        f"NCBI esearch failed after {max_retries} retries: {last_error}"
    )


def fetch_hiv1_sequences(
    email: str,
    gene: str,
    year_start: int,
    year_end: int,
    country: str = "Any",
    min_length: int = 1,
    max_length: int = 10_000,
    max_records: int = 50,
    progress_callback=None,
) -> pd.DataFrame:
    """
    Query NCBI Nucleotide database and return a cleaned DataFrame.

    Columns: genbank_id, year, title, country, sequence
    """
    configure_entrez(email)
    query = build_search_query(gene, year_start, year_end, country)

    if progress_callback:
        progress_callback(0.1, f"Searching NCBI: {query}")

    id_list, total_count = _safe_esearch(query, retmax=max_records)

    if not id_list:
        return pd.DataFrame(columns=["genbank_id", "year", "title", "country", "sequence"])

    if progress_callback:
        progress_callback(
            0.25,
            f"Found {total_count} hit(s); fetching {len(id_list)} record(s)...",
        )

    raw_gb = _safe_efetch(id_list)

    if progress_callback:
        progress_callback(0.5, "Parsing GenBank records...")

    records = list(SeqIO.parse(io.StringIO(raw_gb), "genbank"))
    rows = []

    for i, record in enumerate(records):
        sequence = str(record.seq).upper().replace("U", "T")
        seq_len = len(sequence)

        if seq_len < min_length or seq_len > max_length:
            continue

        genbank_id = record.id or "Unknown"
        title = record.description or "Unknown"
        year = _extract_year(record)
        country_val = _extract_geo_country(record)

        if country_val == "Unknown":
            country_val = _parse_country_from_text(title)

        if year == "Unknown":
            year = _parse_year_from_text(title)

        if not _country_matches_filter(country_val, country):
            continue

        rows.append(
            {
                "genbank_id": genbank_id,
                "year": year,
                "title": title,
                "country": country_val,
                "sequence": sequence,
                "length": seq_len,
            }
        )

        if progress_callback and records:
            progress_callback(0.5 + 0.4 * (i + 1) / len(records), "Parsing records...")

    df = pd.DataFrame(rows)

    if not df.empty:
        df = df.drop_duplicates(subset=["genbank_id"], keep="first")
        df = df[["genbank_id", "year", "title", "country", "sequence"]].reset_index(drop=True)

    if progress_callback:
        progress_callback(1.0, f"Done — {len(df)} sequence(s) after filtering.")

    return df


def dataframe_to_csv(df: pd.DataFrame) -> bytes:
    """Convert DataFrame to CSV bytes for download."""
    return df.to_csv(index=False).encode("utf-8")
