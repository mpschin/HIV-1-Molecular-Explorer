"""
HIV-1 Molecular Explorer
An interactive bioinformatics app for molecular epidemiology and phylogenetics.
"""

import streamlit as st
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from datetime import datetime
import time
import io
import re

# ─── Page Config ─────────────────────────────────────────────────────────────
st.set_page_config(
    page_title="HIV-1 Molecular Explorer",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

# ─── Biopython imports ────────────────────────────────────────────────────────
try:
    from Bio import Entrez, SeqIO
    from Bio import pairwise2
    from Bio.pairwise2 import format_alignment
    BIOPYTHON_OK = True
except ImportError:
    BIOPYTHON_OK = False

# ─── Session State Init ───────────────────────────────────────────────────────
if "df" not in st.session_state:
    st.session_state.df = None
if "fetch_done" not in st.session_state:
    st.session_state.fetch_done = False
if "alignment_result" not in st.session_state:
    st.session_state.alignment_result = None

# ─── Custom CSS ───────────────────────────────────────────────────────────────
st.markdown("""
<style>
    /* Header */
    .main-header {
        font-size: 2.2rem;
        font-weight: 800;
        color: #E8EBF0;
        letter-spacing: -0.5px;
        margin-bottom: 0.1rem;
    }
    .main-header span {
        color: #00C2CB;
    }
    .sub-header {
        font-size: 0.95rem;
        color: #6E7891;
        margin-bottom: 1.5rem;
    }
    /* Section headers */
    .section-title {
        font-size: 0.75rem;
        font-weight: 700;
        letter-spacing: 1.5px;
        color: #6E7891;
        text-transform: uppercase;
        border-bottom: 1px solid #252D40;
        padding-bottom: 0.4rem;
        margin-bottom: 1rem;
    }
    /* Badge styles */
    .badge-pass {
        background: #0D3B26;
        color: #00E676;
        padding: 2px 10px;
        border-radius: 4px;
        font-size: 0.75rem;
        font-weight: 700;
        letter-spacing: 0.5px;
    }
    .badge-fail {
        background: #3B1A1A;
        color: #FF5252;
        padding: 2px 10px;
        border-radius: 4px;
        font-size: 0.75rem;
        font-weight: 700;
        letter-spacing: 0.5px;
    }
    .check-row {
        display: flex;
        align-items: center;
        gap: 12px;
        padding: 6px 0;
        border-bottom: 1px solid #1E2435;
        font-size: 0.88rem;
        color: #B0BAD0;
    }
    /* Alignment block */
    .alignment-block {
        font-family: 'Courier New', monospace;
        font-size: 0.78rem;
        line-height: 1.6;
        background: #111827;
        border: 1px solid #252D40;
        border-radius: 6px;
        padding: 1rem 1.2rem;
        overflow-x: auto;
        white-space: pre;
        color: #B0BAD0;
    }
    .align-id {
        color: #00C2CB;
        font-weight: 700;
    }
    .consensus-label {
        color: #FFD740;
        font-weight: 700;
    }
    /* Metric card override */
    div[data-testid="stMetric"] {
        background: #1A1F2E;
        border: 1px solid #252D40;
        border-radius: 8px;
        padding: 0.8rem 1rem;
    }
    /* Sidebar */
    section[data-testid="stSidebar"] {
        background: #12161F;
    }
</style>
""", unsafe_allow_html=True)

# ─── Header ───────────────────────────────────────────────────────────────────
st.markdown('<div class="main-header"><span>NCBI</span> HIV-1 Explorer & Aligner</div>', unsafe_allow_html=True)
st.markdown('<div class="sub-header">Molecular epidemiology · Sequence retrieval · Multiple alignment · Phylogenetic visualization</div>', unsafe_allow_html=True)

# ─── Sidebar ──────────────────────────────────────────────────────────────────
with st.sidebar:
    st.markdown("## Controls")
    st.divider()

    # ── Connection ─────────────────────────────────────────────────────────
    ncbi_email = st.text_input(
        "NCBI Entrez Email",
        placeholder="your.email@example.com",
        help="Required by NCBI policy to identify your requests.",
    )

    st.divider()

    # ── Date & Volume ──────────────────────────────────────────────────────
    current_year = datetime.now().year
    st.markdown("**Date Range**")
    year_range = st.slider(
        "Year range",
        min_value=1990,
        max_value=current_year,
        value=(current_year - 10, current_year),
        label_visibility="collapsed",
    )

    max_seq = st.slider(
        "Max Sequences",
        min_value=10,
        max_value=300,
        value=150,
        step=10,
        help="Maximum number of sequences to retrieve from NCBI.",
    )

    st.divider()

    # ── Sequence Length ────────────────────────────────────────────────────
    st.markdown("**Sequence Length (nt)**")
    seq_len_range = st.slider(
        "Sequence length",
        min_value=1,
        max_value=10000,
        value=(200, 10000),
        step=50,
        label_visibility="collapsed",
        help="Filter sequences by nucleotide length. Applied both in the NCBI query and after fetching.",
    )
    st.caption(f"{seq_len_range[0]:,} – {seq_len_range[1]:,} nt")

    st.divider()

    # ── Gene Targets ───────────────────────────────────────────────────────
    st.markdown("**Gene Targets**")
    GENE_OPTIONS = ["Any gene", "gag", "pol", "env", "vif", "vpr", "tat", "rev", "vpu", "nef"]
    selected_genes = st.multiselect(
        "Genes",
        options=GENE_OPTIONS[1:],
        default=[],
        label_visibility="collapsed",
        placeholder="Leave empty to search all genes",
        help="Filter by HIV-1 gene region. Multiple selections use OR logic.",
    )

    st.divider()

    # ── Country Filter ─────────────────────────────────────────────────────
    st.markdown("**Country**")
    country_filter = st.text_input(
        "Country",
        value="",
        placeholder="e.g. USA, Kenya, Brazil",
        label_visibility="collapsed",
        help="Filter by country of isolation. Separate multiple countries with commas.",
    )
    country_list = [c.strip() for c in country_filter.split(",") if c.strip()]

    st.divider()
    fetch_btn = st.button("🔬 FETCH NCBI DATA", type="primary", use_container_width=True)

    if st.session_state.fetch_done and st.session_state.df is not None:
        st.success(f"✓ {len(st.session_state.df)} sequences loaded")

    # Active-filter summary
    active = []
    if selected_genes:
        active.append(f"Gene: {', '.join(selected_genes)}")
    if country_list:
        active.append(f"Country: {', '.join(country_list)}")
    active.append(f"Length: {seq_len_range[0]:,}–{seq_len_range[1]:,} nt")
    if active:
        st.caption("**Active filters:** " + " · ".join(active))


# ─── Helper Functions ─────────────────────────────────────────────────────────

def extract_year(record):
    """Extract isolation year from a GenBank record."""
    # Try feature qualifiers
    for feature in record.features:
        if feature.type == "source":
            for key in ("collection_date", "isolation_date", "note"):
                val = feature.qualifiers.get(key, [""])[0]
                m = re.search(r"\b(19|20)\d{2}\b", val)
                if m:
                    return int(m.group())
    # Fall back to record description
    m = re.search(r"\b(19|20)\d{2}\b", record.description)
    if m:
        return int(m.group())
    return "Unknown"


def extract_country(record):
    """
    Extract country/location from a GenBank record source feature.
    Priority order:
      1. /geo_loc_name qualifier (INSDC standard since ~2023)
      2. /country qualifier (legacy)
      3. Record description / title (free-text fallback)
    """
    for feature in record.features:
        if feature.type == "source":
            for key in ("geo_loc_name", "country"):
                vals = feature.qualifiers.get(key, [])
                if vals:
                    # Format can be "Country: Region" — take only the country part
                    return vals[0].split(":")[0].strip()

    # Fallback: scan description for known country names
    desc = record.description.lower()
    # Common HIV-1 surveillance countries (not exhaustive, catches the most frequent)
    COUNTRY_HINTS = [
        "Nigeria", "Kenya", "Uganda", "Tanzania", "Zimbabwe", "South Africa",
        "Cameroon", "Ethiopia", "Mozambique", "Rwanda", "Malawi", "Zambia",
        "USA", "United States", "Brazil", "Argentina", "Mexico", "Colombia",
        "UK", "United Kingdom", "France", "Germany", "Netherlands", "Spain",
        "China", "India", "Thailand", "Cambodia", "Vietnam", "Japan",
        "Russia", "Ukraine",
    ]
    for hint in COUNTRY_HINTS:
        if hint.lower() in desc:
            return hint

    return "Unknown"


def country_matches(record_country: str, record_desc: str, filter_terms: list) -> bool:
    """
    Return True if any of the filter_terms matches the record's country.
    Matching is case-insensitive substring; also checks the description
    so records with "Unknown" qualifier but a country in the title still pass.
    """
    haystack = f"{record_country} {record_desc}".lower()
    return any(term.lower() in haystack for term in filter_terms)


def extract_gene(record):
    """Extract gene name(s) from CDS/gene features of a GenBank record."""
    genes_found = []
    for feature in record.features:
        if feature.type in ("CDS", "gene"):
            gene_val = feature.qualifiers.get("gene", feature.qualifiers.get("product", []))
            for g in gene_val:
                g_clean = g.lower().split()[0]  # e.g. "gag protein" → "gag"
                if g_clean not in genes_found:
                    genes_found.append(g_clean)
    return ", ".join(genes_found) if genes_found else "Unknown"


def build_ncbi_query(start_year, end_year, seq_len_min, seq_len_max, genes, countries):
    """
    Build an NCBI Entrez query string from search parameters.

    Country: [Country] field is unreliable for nucleotide records, so we add
    each country as an [All Fields] free-text term instead. This biases NCBI
    to return country-relevant records. A second post-fetch pass (using the
    /geo_loc_name and /country source qualifiers) confirms the match.
    """
    parts = [
        '("HIV-1"[Organism] OR "Human immunodeficiency virus 1"[Organism])',
        f'{start_year}:{end_year}[PDAT]',
        f'{seq_len_min}:{seq_len_max}[SLEN]',
    ]

    # Gene filter — OR across selected genes
    if genes:
        gene_terms = " OR ".join(f'"{g}"[Gene]' for g in genes)
        parts.append(f"({gene_terms})")

    # Country filter — free-text OR so NCBI pre-biases results toward these locations
    if countries:
        country_terms = " OR ".join(f'"{c}"[All Fields]' for c in countries)
        parts.append(f"({country_terms})")

    return " AND ".join(parts)


def fetch_hiv_sequences(email, start_year, end_year, max_seqs,
                        seq_len_range=(1, 10000), genes=None, countries=None):
    """Fetch HIV-1 sequences from NCBI Nucleotide database."""
    Entrez.email = email
    genes = genes or []
    countries = countries or []

    query = build_ncbi_query(
        start_year, end_year,
        seq_len_range[0], seq_len_range[1],
        genes, countries,
    )

    # When country filtering is active, fetch a larger pool so that post-fetch
    # filtering still yields enough results. Many GenBank records lack explicit
    # /country or /geo_loc_name qualifiers, so we need headroom to filter from.
    fetch_size = max_seqs * 4 if countries else max_seqs

    records = []

    with st.status("Querying NCBI Entrez...", expanded=True) as status:
        st.write(f"🔍 Query: `{query[:140]}{'…' if len(query) > 140 else ''}`")
        if countries:
            st.write(f"🌍 Country post-filter active: **{', '.join(countries)}** "
                     f"(fetching up to {fetch_size} records to filter from)")
        try:
            handle = Entrez.esearch(db="nucleotide", term=query, retmax=fetch_size)
            search_results = Entrez.read(handle)
            handle.close()
        except Exception as e:
            st.error(f"Search failed: {e}")
            status.update(label="Search failed", state="error")
            return None

        id_list = search_results.get("IdList", [])
        if not id_list:
            st.warning("No sequences found. Try broadening your filters.")
            status.update(label="No results", state="error")
            return None

        st.write(f"📋 Found {len(id_list)} records. Fetching full GenBank entries…")

        # Fetch in one batch (NCBI allows up to ~500 IDs per request)
        try:
            handle = Entrez.efetch(
                db="nucleotide",
                id=",".join(id_list),
                rettype="gb",
                retmode="text",
            )
            raw = handle.read()
            handle.close()
        except Exception as e:
            st.error(f"Fetch failed: {e}")
            status.update(label="Fetch failed", state="error")
            return None

        st.write("🧬 Parsing GenBank records…")
        gb_records = list(SeqIO.parse(io.StringIO(raw), "genbank"))

        progress = st.progress(0)
        total = len(gb_records)

        for i, rec in enumerate(gb_records):
            # Stop once we have enough records
            if len(records) >= max_seqs:
                progress.progress(1.0)
                break

            seq_str = str(rec.seq).upper()
            seq_len = len(seq_str)

            # Skip empty / too-short sequences
            if not seq_str or seq_len < 20:
                progress.progress((i + 1) / total)
                continue

            # Enforce sequence length range
            if not (seq_len_range[0] <= seq_len <= seq_len_range[1]):
                progress.progress((i + 1) / total)
                continue

            year = extract_year(rec)
            if year != "Unknown":
                try:
                    yr_int = int(year)
                    if not (start_year <= yr_int <= end_year):
                        progress.progress((i + 1) / total)
                        continue
                except ValueError:
                    pass

            country = extract_country(rec)

            # Post-fetch country filter: use lenient matcher that also checks
            # the record description so records with missing qualifiers still match.
            if countries and not country_matches(country, rec.description, countries):
                progress.progress((i + 1) / total)
                continue

            gene = extract_gene(rec)

            records.append({
                "GenBank ID": rec.id,
                "Title": rec.description[:60] + ("…" if len(rec.description) > 60 else ""),
                "Year": year,
                "Country": country,
                "Gene": gene,
                "Length (nt)": seq_len,
                "Sequence (DNA/RNA)": seq_str[:40] + "…",
                "_full_seq": seq_str,
            })
            progress.progress((i + 1) / total)

        status.update(label=f"✅ Retrieved {len(records)} valid sequences.", state="complete")

    if not records:
        st.warning(
            "No sequences matched the country filter after fetching. "
            "Try checking the spelling or using a broader search — "
            "some records store the location differently in NCBI (e.g. 'USA' vs 'United States')."
        )
        return None

    return pd.DataFrame(records)


def simple_consensus(sequences):
    """Compute majority-vote consensus from a list of equal-length strings."""
    if not sequences:
        return ""
    length = max(len(s) for s in sequences)
    padded = [s.ljust(length, "-") for s in sequences]
    consensus = []
    for i in range(length):
        col = [s[i] for s in padded]
        base = max(set(col), key=col.count)
        consensus.append(base if col.count(base) > len(col) // 2 else "·")
    return "".join(consensus)


def run_msa(seqs_dict, max_seqs=10):
    """
    Run pairwise-progressive alignment on the first `max_seqs` sequences.
    Returns list of (id, aligned_seq) tuples.
    """
    ids = list(seqs_dict.keys())[:max_seqs]
    seqs = [seqs_dict[i] for i in ids]

    if len(seqs) < 2:
        return [(ids[0], seqs[0])], len(seqs[0])

    # Use shortest sequence as anchor; align all others against it pairwise
    anchor_idx = min(range(len(seqs)), key=lambda i: len(seqs[i]))
    anchor = seqs[anchor_idx]
    aligned = {}

    for idx, (sid, seq) in enumerate(zip(ids, seqs)):
        if idx == anchor_idx:
            aligned[sid] = seq
            continue
        # pairwise2 global alignment with match/mismatch/gap params
        alns = pairwise2.align.globalms(anchor, seq, 2, -1, -2, -0.5)
        if alns:
            _, aligned_seq, *_ = alns[0]
            aligned[sid] = aligned_seq
        else:
            aligned[sid] = seq

    result = [(sid, aligned[sid]) for sid in ids if sid in aligned]
    aln_len = max(len(s) for _, s in result) if result else 0
    return result, aln_len


def render_alignment_html(aligned_seqs, aln_len, max_cols=80):
    """Render aligned sequences as styled HTML block with consensus."""
    lines = []
    consensus = simple_consensus([s for _, s in aligned_seqs])

    for sid, seq in aligned_seqs:
        padded = seq.ljust(aln_len, "-")[:aln_len]
        label = f'<span class="align-id">{sid[:14]:<14}</span>'
        seq_colored = ""
        for base in padded[:max_cols]:
            if base == "-":
                seq_colored += f'<span style="color:#4A5568">-</span>'
            elif base == "A":
                seq_colored += f'<span style="color:#F6AD55">{base}</span>'
            elif base == "T":
                seq_colored += f'<span style="color:#68D391">{base}</span>'
            elif base == "G":
                seq_colored += f'<span style="color:#63B3ED">{base}</span>'
            elif base == "C":
                seq_colored += f'<span style="color:#FC8181">{base}</span>'
            else:
                seq_colored += f'<span style="color:#9AA5C0">{base}</span>'
        if aln_len > max_cols:
            seq_colored += f'<span style="color:#4A5568">…</span>'
        lines.append(f"{label}  {seq_colored}")

    # Consensus line
    cons_display = consensus[:max_cols]
    if aln_len > max_cols:
        cons_display += "…"
    lines.append(f'<span class="consensus-label">{"Consensus":<14}</span>  <span style="color:#9AA5C0">{cons_display}</span>')

    body = "\n".join(lines)
    return f'<div class="alignment-block">{body}</div>'


# ─── Fetch Logic ──────────────────────────────────────────────────────────────
if fetch_btn:
    if not ncbi_email or "@" not in ncbi_email:
        st.error("Please enter a valid email address in the sidebar (required by NCBI).")
    elif not BIOPYTHON_OK:
        st.error("Biopython is not installed. Please restart and try again.")
    else:
        df = fetch_hiv_sequences(
            ncbi_email,
            year_range[0], year_range[1],
            max_seq,
            seq_len_range=seq_len_range,
            genes=selected_genes,
            countries=country_list,
        )
        if df is not None and len(df) > 0:
            st.session_state.df = df
            st.session_state.fetch_done = True
            st.session_state.alignment_result = None
            st.rerun()

# ─── Main Content ─────────────────────────────────────────────────────────────
df = st.session_state.df

if df is None:
    # Empty state
    st.info("👈 Enter your NCBI email, configure parameters in the sidebar, then click **FETCH NCBI DATA** to begin.")

    col1, col2, col3 = st.columns(3)
    with col1:
        st.metric("Data Source", "NCBI Nucleotide")
    with col2:
        st.metric("Target Pathogen", "HIV-1")
    with col3:
        st.metric("Alignment Engine", "Biopython pairwise2")

else:
    # ── Tab Layout ─────────────────────────────────────────────────────────
    tab1, tab2, tab3 = st.tabs(["📊 Data Overview", "🧬 Alignment", "🔎 Integrity Check"])

    # ── Tab 1: Data Overview ───────────────────────────────────────────────
    with tab1:
        st.markdown('<div class="section-title">Data Overview</div>', unsafe_allow_html=True)

        col_info, col_dl = st.columns([3, 1])
        with col_info:
            st.caption(f"**{len(df)}** sequences fetched.")
        with col_dl:
            display_cols = ["GenBank ID", "Title", "Year", "Country", "Gene", "Length (nt)", "Sequence (DNA/RNA)"]
            csv_bytes = df[display_cols].to_csv(index=False).encode()
            st.download_button(
                "⬇ Download CSV",
                data=csv_bytes,
                file_name="hiv1_sequences.csv",
                mime="text/csv",
                use_container_width=True,
            )

        st.dataframe(
            df[display_cols],
            use_container_width=True,
            height=280,
            hide_index=True,
        )

        st.divider()

        # Charts
        chart_col1, chart_col2 = st.columns(2)

        with chart_col1:
            st.markdown('<div class="section-title">Top Countries (Past 10 Yrs)</div>', unsafe_allow_html=True)
            country_df = df[df["Country"] != "Unknown"]["Country"].value_counts().head(10).reset_index()
            country_df.columns = ["Country", "Count"]
            fig_bar = px.bar(
                country_df,
                x="Count",
                y="Country",
                orientation="h",
                color="Count",
                color_continuous_scale=[[0, "#1A3A4A"], [0.5, "#007A84"], [1.0, "#00C2CB"]],
                template="plotly_dark",
                labels={"Count": "Sequences", "Country": ""},
            )
            fig_bar.update_layout(
                plot_bgcolor="#1A1F2E",
                paper_bgcolor="#1A1F2E",
                margin=dict(l=0, r=10, t=10, b=10),
                coloraxis_showscale=False,
                yaxis=dict(autorange="reversed"),
                height=280,
            )
            fig_bar.update_traces(marker_line_width=0)
            st.plotly_chart(fig_bar, use_container_width=True)

        with chart_col2:
            st.markdown('<div class="section-title">Sequence Timeline</div>', unsafe_allow_html=True)
            yr_df = df[df["Year"] != "Unknown"].copy()
            yr_df["Year"] = pd.to_numeric(yr_df["Year"], errors="coerce").dropna()
            if not yr_df.empty:
                timeline = yr_df.groupby("Year").size().reset_index(name="Count")
                fig_line = px.area(
                    timeline,
                    x="Year",
                    y="Count",
                    template="plotly_dark",
                    labels={"Count": "Sequences", "Year": "Year"},
                    color_discrete_sequence=["#00C2CB"],
                )
                fig_line.update_traces(
                    fill="tozeroy",
                    line_color="#00C2CB",
                    fillcolor="rgba(0,194,203,0.18)",
                )
                fig_line.update_layout(
                    plot_bgcolor="#1A1F2E",
                    paper_bgcolor="#1A1F2E",
                    margin=dict(l=0, r=10, t=10, b=10),
                    height=280,
                    xaxis=dict(tickformat="d"),
                )
                st.plotly_chart(fig_line, use_container_width=True)
            else:
                st.info("No year data available for timeline.")

    # ── Tab 2: Alignment ───────────────────────────────────────────────────
    with tab2:
        st.markdown('<div class="section-title">Integrity Check & Alignment</div>', unsafe_allow_html=True)

        valid_seqs = df[df["_full_seq"].str.len() > 20].copy()

        col_a, col_b = st.columns([3, 1])
        with col_a:
            num_seqs = st.slider(
                "Sequences to align (fewer = faster)",
                min_value=2,
                max_value=min(20, len(valid_seqs)),
                value=min(6, len(valid_seqs)),
                help="Multiple sequence alignment using Biopython pairwise2.",
            )
        with col_b:
            st.write("")
            run_btn = st.button("▶ RUN SEQUENCE ALIGNMENT", type="primary", use_container_width=True)

        if run_btn:
            subset = valid_seqs.head(num_seqs)
            seqs_dict = dict(zip(subset["GenBank ID"], subset["_full_seq"]))
            with st.spinner(f"Aligning {num_seqs} sequences…"):
                try:
                    result, aln_len = run_msa(seqs_dict, max_seqs=num_seqs)
                    st.session_state.alignment_result = (result, aln_len)
                except Exception as e:
                    st.error(f"Alignment error: {e}")

        if st.session_state.alignment_result:
            result, aln_len = st.session_state.alignment_result
            st.caption(f"Alignment Length: **{aln_len:,} bp** · {len(result)} sequences")
            st.markdown(
                render_alignment_html(result, aln_len),
                unsafe_allow_html=True,
            )
            # Download alignment
            aln_text = "\n".join([f">{sid}\n{seq}" for sid, seq in result])
            st.download_button(
                "⬇ Download FASTA Alignment",
                data=aln_text.encode(),
                file_name="hiv1_alignment.fasta",
                mime="text/plain",
            )
        else:
            st.info("Configure the number of sequences above, then click **RUN SEQUENCE ALIGNMENT**.")

    # ── Tab 3: Data Integrity ──────────────────────────────────────────────
    with tab3:
        st.markdown('<div class="section-title">Data Integrity Dashboard</div>', unsafe_allow_html=True)

        def badge(passed):
            cls = "badge-pass" if passed else "badge-fail"
            label = "PASS" if passed else "FAIL"
            return f'<span class="{cls}">{label}</span>'

        checks = []

        # Check 1: Null sequences
        null_seqs = df["_full_seq"].isna().sum() + (df["_full_seq"] == "").sum()
        null_pass = (null_seqs == 0)
        checks.append((null_pass, "Null Sequence Check", f"{null_seqs} null/empty sequences found"))

        # Check 2: Year logic
        current_year = datetime.now().year
        yr_numeric = pd.to_numeric(df["Year"], errors="coerce")
        bad_years = yr_numeric[(yr_numeric < 1980) | (yr_numeric > current_year)].count()
        year_pass = (bad_years == 0)
        checks.append((year_pass, "Year Logic Check", f"{bad_years} records outside 1980–{current_year}"))

        # Check 3: GenBank ID format
        id_valid = df["GenBank ID"].str.match(r"^[A-Z]{1,2}\d{5,8}(\.\d+)?$", na=False)
        bad_ids = (~id_valid).sum()
        id_pass = (bad_ids == 0)
        checks.append((id_pass, "GenBank ID Format", f"{bad_ids} IDs with unexpected format"))

        # Check 4: Country coverage
        unknown_countries = (df["Country"] == "Unknown").sum()
        country_pct = unknown_countries / len(df) * 100
        country_pass = (country_pct < 30)
        checks.append((country_pass, "Country Coverage", f"{unknown_countries} ({country_pct:.0f}%) have unknown country"))

        # Check 5: Sequence length diversity
        seq_lens = df["_full_seq"].str.len()
        min_len, max_len = seq_lens.min(), seq_lens.max()
        length_pass = (min_len >= 20)
        checks.append((length_pass, "Minimum Sequence Length", f"Shortest: {min_len} bp, Longest: {max_len} bp"))

        for passed, name, detail in checks:
            st.markdown(
                f'<div class="check-row">{badge(passed)}<span><b>{name}</b> — {detail}</span></div>',
                unsafe_allow_html=True,
            )

        st.divider()

        # Summary metrics
        passed_count = sum(1 for p, _, _ in checks if p)
        m1, m2, m3, m4 = st.columns(4)
        m1.metric("Total Records", len(df))
        m2.metric("Valid Sequences", int((df["_full_seq"].str.len() > 20).sum()))
        m3.metric("Integrity Checks", f"{passed_count}/{len(checks)}")
        m4.metric("Countries Represented", int((df["Country"] != "Unknown").sum() and df[df["Country"] != "Unknown"]["Country"].nunique()))

        st.divider()

        # Distribution plots
        st.markdown('<div class="section-title">Sequence Length Distribution</div>', unsafe_allow_html=True)
        fig_hist = px.histogram(
            x=df["_full_seq"].str.len(),
            nbins=40,
            template="plotly_dark",
            color_discrete_sequence=["#00C2CB"],
            labels={"x": "Sequence Length (bp)", "y": "Count"},
        )
        fig_hist.update_layout(
            plot_bgcolor="#1A1F2E",
            paper_bgcolor="#1A1F2E",
            margin=dict(l=0, r=10, t=10, b=10),
            height=220,
        )
        st.plotly_chart(fig_hist, use_container_width=True)
