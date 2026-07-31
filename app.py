"""
HIV-1 Molecular Explorer
Main Streamlit UI and application layout.
"""

from __future__ import annotations

import io
from datetime import datetime

import matplotlib.pyplot as plt
import pandas as pd
import streamlit as st

from aligner import (
    AlignmentError,
    align_sequences,
    build_phylogenetic_tree,
    plot_tree_matplotlib,
    plot_tree_plotly,
    run_full_phylogenetics,
)
from ncbi_fetcher import NCBIFetchError, dataframe_to_csv, fetch_hiv1_sequences
from utils import run_data_integrity_check

# ---------------------------------------------------------------------------
# Page config
# ---------------------------------------------------------------------------
st.set_page_config(
    page_title="HIV-1 Molecular Explorer",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

CURRENT_YEAR = datetime.now().year
DEFAULT_YEAR_START = CURRENT_YEAR - 10

COUNTRY_OPTIONS = [
    "Any",
    "United States",
    "South Africa",
    "India",
    "Brazil",
    "China",
    "Thailand",
    "Uganda",
    "Kenya",
    "United Kingdom",
    "France",
    "Germany",
    "Australia",
    "Custom",
]

GENE_OPTIONS = ["gag", "env", "pol", "Custom"]


def _init_session_state() -> None:
    defaults = {
        "df": None,
        "msa": None,
        "tree": None,
        "fasta": None,
        "newick": None,
        "integrity_result": None,
        "fetch_query": "",
    }
    for key, val in defaults.items():
        if key not in st.session_state:
            st.session_state[key] = val


def _progress_bar(placeholder, fraction: float, message: str) -> None:
    placeholder.progress(min(max(fraction, 0.0), 1.0), text=message)


# ---------------------------------------------------------------------------
# Sidebar
# ---------------------------------------------------------------------------
def render_sidebar() -> dict:
    st.sidebar.header("🔬 Search Configuration")

    email = st.sidebar.text_input(
        "Email Address *",
        value="researcher@example.com",
        help="NCBI Entrez requires a contact email for API access.",
    )

    gene_choice = st.sidebar.selectbox(
        "Target Gene",
        GENE_OPTIONS,
        help=(
            "HIV-1 gene region to search. "
            "**gag** — structural core proteins; "
            "**env** — envelope glycoproteins; "
            "**pol** — reverse transcriptase / integrase."
        ),
    )
    if gene_choice == "Custom":
        gene = st.sidebar.text_input("Custom gene name", value="vif")
    else:
        gene = gene_choice

    country_mode = st.sidebar.selectbox(
        "Country / Region",
        COUNTRY_OPTIONS,
        help=(
            "Filters by **sample isolation country** (geo_loc_name in GenBank), "
            "not the submitting lab's country."
        ),
    )
    if country_mode == "Custom":
        country = st.sidebar.text_input("Enter country name", value="")
    else:
        country = country_mode

    st.sidebar.subheader("Year Range")
    year_range = st.sidebar.slider(
        "Publication / isolation year",
        min_value=1990,
        max_value=CURRENT_YEAR,
        value=(DEFAULT_YEAR_START, CURRENT_YEAR),
        help="Filter records by NCBI publication date (PDAT).",
    )

    st.sidebar.subheader("Sequence Length Filter")
    length_range = st.sidebar.slider(
        "Length (nucleotides)",
        min_value=1,
        max_value=10_000,
        value=(500, 5_000),
        step=50,
        help="Keep only sequences within this length range after fetching.",
    )

    max_records = st.sidebar.number_input(
        "Max Records Limit",
        min_value=5,
        max_value=500,
        value=50,
        step=5,
        help="Maximum number of GenBank records to retrieve (capped at 500 for speed).",
    )

    tree_method = st.sidebar.radio(
        "Tree Method",
        ["nj", "upgma"],
        format_func=lambda x: "Neighbor-Joining (NJ)" if x == "nj" else "UPGMA",
        help=(
            "**Neighbor-Joining** — standard distance-based method for molecular epidemiology. "
            "**UPGMA** — assumes molecular clock; faster but less accurate for divergent sequences."
        ),
    )

    with st.sidebar.expander("ℹ️ Glossary"):
        st.markdown(
            """
            - **PDAT** — Publication date filter in NCBI Entrez.
            - **IUPAC codes** — Standard nucleotide ambiguity letters (e.g., N = any base).
            - **Jukes-Cantor** — Simple nucleotide substitution model correcting for multiple hits.
            - **Newick (.nwk)** — Standard text format for phylogenetic trees.
            - **Profile HMM** — Hidden Markov Model profile for sensitive remote homology detection.
            - **mBed clustering** — Fast approximate phylogenetic clustering for large datasets.
            """
        )

    return {
        "email": email,
        "gene": gene,
        "country": country,
        "year_start": year_range[0],
        "year_end": year_range[1],
        "min_length": length_range[0],
        "max_length": length_range[1],
        "max_records": int(max_records),
        "tree_method": tree_method,
    }


# ---------------------------------------------------------------------------
# Main panels
# ---------------------------------------------------------------------------
def render_header() -> None:
    st.title("🧬 BioVibe")
    st.subheader("HIV-1 NCBI Phylogenetics Explorer")
    st.markdown(
        "Search, collect, align, and visualize **HIV-1** genomic sequences "
        "directly from the NCBI Nucleotide database. Built for molecular epidemiologists."
    )


def render_fetch_panel(config: dict) -> None:
    st.header("📡 Data Acquisition")

    col1, col2 = st.columns([1, 3])
    with col1:
        fetch_btn = st.button("🔍 Fetch Sequences from NCBI", type="primary", use_container_width=True)
    with col2:
        if st.session_state.fetch_query:
            st.caption(f"Last query: `{st.session_state.fetch_query}`")

    if fetch_btn:
        if not config["email"] or "@" not in config["email"]:
            st.error("Please provide a valid email address for NCBI Entrez.")
            return

        progress_slot = st.empty()
        status_slot = st.empty()

        def on_progress(fraction, message):
            _progress_bar(progress_slot, fraction, message)
            status_slot.info(message)

        try:
            with st.spinner("Contacting NCBI Entrez..."):
                from ncbi_fetcher import build_search_query

                st.session_state.fetch_query = build_search_query(
                    config["gene"],
                    config["year_start"],
                    config["year_end"],
                    config["country"],
                )
                df = fetch_hiv1_sequences(
                    email=config["email"],
                    gene=config["gene"],
                    year_start=config["year_start"],
                    year_end=config["year_end"],
                    country=config["country"],
                    min_length=config["min_length"],
                    max_length=config["max_length"],
                    max_records=config["max_records"],
                    progress_callback=on_progress,
                )
            st.session_state.df = df
            st.session_state.msa = None
            st.session_state.tree = None
            st.session_state.fasta = None
            st.session_state.newick = None
            st.session_state.integrity_result = None

            progress_slot.empty()
            if df.empty:
                status_slot.warning("No sequences matched your filters. Try broadening the search.")
            else:
                status_slot.success(f"Fetched **{len(df)}** sequence(s) successfully.")

        except NCBIFetchError as exc:
            progress_slot.empty()
            status_slot.error(f"NCBI fetch error: {exc}")
        except Exception as exc:
            progress_slot.empty()
            status_slot.error(f"Unexpected error during fetch: {exc}")

    df = st.session_state.df
    if df is not None and not df.empty:
        st.dataframe(
            df[["genbank_id", "year", "country", "title"]].assign(
                length=df["sequence"].str.len()
            ),
            use_container_width=True,
            hide_index=True,
        )
        st.download_button(
            label="⬇️ Download CSV",
            data=dataframe_to_csv(df),
            file_name="biovibe_hiv1_sequences.csv",
            mime="text/csv",
        )
    elif df is not None and df.empty:
        st.info("No records to display.")


def render_integrity_panel() -> None:
    st.header("✅ Data Integrity Check")

    if st.button("Run Data Integrity Check", use_container_width=False):
        result = run_data_integrity_check(st.session_state.df)
        st.session_state.integrity_result = result

    result = st.session_state.integrity_result
    if result is None:
        st.caption("Fetch data first, then run the integrity check.")
        return

    if result.passed and not result.warnings:
        st.success("✔ All integrity checks passed.")
    elif result.passed:
        st.success("✔ Integrity checks passed with warnings.")
    else:
        st.error("✘ Integrity check failed.")

    with st.expander("Check details", expanded=not result.passed):
        for msg in result.messages:
            st.write(f"• {msg}")
        for warn in result.warnings:
            st.warning(warn)


def render_alignment_panel(config: dict) -> None:
    st.header("🧬 Alignment & Phylogenetics")

    df = st.session_state.df
    if df is None or df.empty:
        st.info("Fetch sequences before running alignment.")
        return

    run_btn = st.button("⚙️ Run Alignment & Build Tree", type="primary")

    if run_btn:
        progress_slot = st.empty()
        status_slot = st.empty()

        def on_progress(fraction, message):
            _progress_bar(progress_slot, fraction, message)
            status_slot.info(message)

        try:
            with st.spinner("Aligning sequences and constructing phylogenetic tree..."):
                msa, tree, fasta, newick = run_full_phylogenetics(
                    df,
                    tree_method=config["tree_method"],
                    progress_callback=on_progress,
                )
            st.session_state.msa = msa
            st.session_state.tree = tree
            st.session_state.fasta = fasta
            st.session_state.newick = newick
            progress_slot.empty()
            status_slot.success("Alignment and tree construction complete.")
        except AlignmentError as exc:
            progress_slot.empty()
            status_slot.error(f"Alignment error: {exc}")
        except Exception as exc:
            progress_slot.empty()
            status_slot.error(f"Unexpected error: {exc}")

    if st.session_state.tree is not None:
        tab_plotly, tab_mpl, tab_downloads = st.tabs(
            ["Interactive Tree (Plotly)", "Static Tree (Matplotlib)", "Downloads"]
        )

        with tab_plotly:
            fig_plotly = plot_tree_plotly(st.session_state.tree)
            if fig_plotly is not None:
                st.plotly_chart(fig_plotly, use_container_width=True)
            else:
                st.warning("Plotly unavailable — use the Matplotlib tab.")

        with tab_mpl:
            fig = plot_tree_matplotlib(st.session_state.tree)
            st.pyplot(fig)
            plt.close(fig)

        with tab_downloads:
            dl_col1, dl_col2 = st.columns(2)
            with dl_col1:
                st.download_button(
                    "⬇️ Download Aligned FASTA",
                    data=st.session_state.fasta,
                    file_name="biovibe_aligned.fasta",
                    mime="text/plain",
                )
            with dl_col2:
                st.download_button(
                    "⬇️ Download Newick Tree (.nwk)",
                    data=st.session_state.newick,
                    file_name="biovibe_tree.nwk",
                    mime="text/plain",
                )

            with st.expander("Preview Newick"):
                st.code(st.session_state.newick, language=None)


def render_summary_panel(config: dict) -> None:
    df = st.session_state.df
    if df is None or df.empty:
        return

    st.header("📊 Summary Statistics")
    c1, c2, c3, c4 = st.columns(4)
    c1.metric("Sequences", len(df))
    c2.metric("Countries", df["country"].nunique())
    c3.metric("Median length (nt)", int(df["sequence"].str.len().median()))
    c4.metric("Year range", f"{df['year'].min()}–{df['year'].max()}")

    with st.expander("Country breakdown"):
        country_counts = df["country"].value_counts().reset_index()
        country_counts.columns = ["Country", "Count"]
        st.bar_chart(country_counts, x="Country", y="Count")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------
def main() -> None:
    _init_session_state()
    config = render_sidebar()
    render_header()

    tab_data, tab_integrity, tab_phylo = st.tabs(
        ["📡 Fetch Data", "✅ Integrity", "🌳 Phylogenetics"]
    )

    with tab_data:
        render_fetch_panel(config)
        render_summary_panel(config)

    with tab_integrity:
        render_integrity_panel()

    with tab_phylo:
        render_alignment_panel(config)


if __name__ == "__main__":
    main()
