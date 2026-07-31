"""Sequence alignment engine and Newick tree generator for HIV-1 Molecular Explorer."""

from __future__ import annotations

import io
import math
from collections import Counter
from typing import Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import Phylo, SeqIO
from Bio.Align import MultipleSeqAlignment, PairwiseAligner
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceMatrix, DistanceTreeConstructor
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

IUPAC_COMPLEMENT = {
    "A": "T", "T": "A", "C": "G", "G": "C",
    "N": "N", "R": "Y", "Y": "R", "S": "S",
    "W": "W", "K": "M", "M": "K", "B": "V",
    "V": "B", "D": "H", "H": "D",
}


class AlignmentError(Exception):
    """Raised when alignment or tree construction fails."""


def _make_label(row: pd.Series) -> str:
    """Format leaf label as GenBankID|Country|Year."""
    gid = str(row.get("genbank_id", "Unknown"))
    country = str(row.get("country", "Unknown"))
    year = str(row.get("year", "Unknown"))
    return f"{gid}|{country}|{year}"


def _kmer_set(sequence: str, k: int = 5) -> Counter:
    """Return k-mer counts for a sequence."""
    seq = sequence.upper().replace("U", "T")
    if len(seq) < k:
        return Counter({seq: 1})
    return Counter(seq[i : i + k] for i in range(len(seq) - k + 1))


def _kmer_distance(seq_a: str, seq_b: str, k: int = 5) -> float:
    """Jaccard distance between k-mer profiles."""
    ka = _kmer_set(seq_a, k)
    kb = _kmer_set(seq_b, k)
    if not ka and not kb:
        return 0.0
    intersection = sum((ka & kb).values())
    union = sum((ka | kb).values())
    if union == 0:
        return 1.0
    return 1.0 - (intersection / union)


def _pairwise_identity(seq_a: str, seq_b: str) -> float:
    """Fraction of identical positions (ungapped, same length assumed after pad)."""
    matches = sum(1 for a, b in zip(seq_a, seq_b) if a == b and a != "-")
    valid = sum(1 for a, b in zip(seq_a, seq_b) if a != "-" and b != "-")
    return matches / valid if valid else 0.0


def _progressive_pairwise_align(sequences: List[str]) -> List[str]:
    """
    Lightweight progressive alignment via guide tree from k-mer distances.

    For small datasets uses Bio.Align PairwiseAligner; merges progressively.
    """
    n = len(sequences)
    if n == 0:
        return []
    if n == 1:
        return [sequences[0]]

    if n > 100:
        max_len = max(len(s) for s in sequences)
        return [s.ljust(max_len, "-") for s in sequences]

    dist_matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            d = _kmer_distance(sequences[i], sequences[j])
            dist_matrix[i, j] = d
            dist_matrix[j, i] = d

    aligned: Dict[int, str] = {i: sequences[i] for i in range(n)}
    active = set(range(n))

    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -2
    aligner.extend_gap_score = -0.5

    while len(active) > 1:
        best_pair: Optional[Tuple[int, int]] = None
        best_dist = float("inf")
        active_list = sorted(active)
        for ii, i in enumerate(active_list):
            for j in active_list[ii + 1 :]:
                if dist_matrix[i, j] < best_dist:
                    best_dist = dist_matrix[i, j]
                    best_pair = (i, j)

        if best_pair is None:
            break

        i, j = best_pair
        try:
            result = aligner.align(aligned[i], aligned[j])[0]
            merged = str(result[0]) if hasattr(result, "__getitem__") else str(result)
        except Exception:
            max_len = max(len(aligned[i]), len(aligned[j]))
            merged = aligned[i].ljust(max_len, "-") + aligned[j].ljust(max_len, "-")

        aligned[i] = merged
        active.discard(j)

    remaining = next(iter(active))
    return [aligned[remaining]] if n == 1 else _realign_all(sequences, aligner, max_records=50)


def _realign_all(
    sequences: List[str],
    aligner: PairwiseAligner,
    max_records: int = 50,
) -> List[str]:
    """Align all sequences to the longest sequence as reference."""
    if len(sequences) > max_records:
        max_len = max(len(s) for s in sequences)
        return [s.ljust(max_len, "-") for s in sequences]

    ref_idx = max(range(len(sequences)), key=lambda i: len(sequences[i]))
    ref = sequences[ref_idx]
    aligned = [None] * len(sequences)
    aligned[ref_idx] = ref

    for i, seq in enumerate(sequences):
        if i == ref_idx:
            continue
        try:
            result = aligner.align(ref, seq)[0]
            if len(result[0]) >= len(result[1]):
                aligned[i] = str(result[1])
                if len(result[0]) > len(ref):
                    ref = str(result[0])
                    aligned[ref_idx] = ref
            else:
                aligned[i] = str(result[0])
        except Exception:
            aligned[i] = seq.ljust(len(ref), "-")

    max_len = max(len(s) for s in aligned if s)
    return [s.ljust(max_len, "-") if s else "-" * max_len for s in aligned]


def align_sequences(df: pd.DataFrame, progress_callback=None) -> Tuple[MultipleSeqAlignment, pd.DataFrame]:
    """
    Align sequences from a DataFrame and return MultipleSeqAlignment plus annotated df.

    Uses progressive k-mer guide-tree alignment with Bio.Align fallback.
    """
    if df is None or df.empty:
        raise AlignmentError("No sequences to align.")

    if progress_callback:
        progress_callback(0.1, "Preparing sequences...")

    labels = [_make_label(row) for _, row in df.iterrows()]
    sequences = [str(s).upper().replace("U", "T") for s in df["sequence"]]

    n = len(sequences)
    if n > 100:
        if progress_callback:
            progress_callback(0.3, f"Large dataset ({n} seqs) — using pad-only fallback.")
        max_len = max(len(s) for s in sequences)
        aligned_seqs = [s.ljust(max_len, "-") for s in sequences]
    else:
        if progress_callback:
            progress_callback(0.3, "Running progressive alignment...")
        aligner = PairwiseAligner()
        aligner.mode = "global"
        aligner.match_score = 2
        aligner.mismatch_score = -1
        aligner.open_gap_score = -2
        aligner.extend_gap_score = -0.5
        aligned_seqs = _realign_all(sequences, aligner)

    if progress_callback:
        progress_callback(0.7, "Building alignment object...")

    records = [
        SeqRecord(Seq(seq), id=label, description="")
        for label, seq in zip(labels, aligned_seqs)
    ]
    msa = MultipleSeqAlignment(records)

    annotated = df.copy()
    annotated["label"] = labels

    if progress_callback:
        progress_callback(1.0, f"Aligned {n} sequence(s).")

    return msa, annotated


def _build_distance_matrix(msa: MultipleSeqAlignment) -> DistanceMatrix:
    """Build a Jukes-Cantor distance matrix from aligned sequences."""
    calculator = DistanceCalculator("identity")
    dm = calculator.get_distance(msa)

    names = list(dm.names)
    n = len(names)
    matrix = []
    for i in range(n):
        row = []
        for j in range(i + 1):
            if i == j:
                row.append(0.0)
            else:
                p = min(float(dm[i, j]), 0.749)
                jc = -0.75 * math.log(max(1e-10, 1.0 - (4.0 / 3.0) * p))
                row.append(max(0.0, jc))
        matrix.append(row)

    return DistanceMatrix(names, matrix)


def build_phylogenetic_tree(
    msa: MultipleSeqAlignment,
    method: str = "nj",
    progress_callback=None,
):
    """
    Construct a Neighbor-Joining (default) or UPGMA tree from alignment.

    Leaf labels should already be formatted as GenBankID|Country|Year.
    """
    if msa is None or len(msa) == 0:
        raise AlignmentError("Empty alignment — cannot build tree.")

    if progress_callback:
        progress_callback(0.2, "Computing distance matrix (Jukes-Cantor)...")

    dm = _build_distance_matrix(msa)

    if progress_callback:
        progress_callback(0.6, f"Building tree ({method.upper()})...")

    constructor = DistanceTreeConstructor()
    if method.lower() == "upgma":
        tree = constructor.upgma(dm)
    else:
        tree = constructor.nj(dm)

    if progress_callback:
        progress_callback(1.0, "Tree construction complete.")

    return tree


def tree_to_newick(tree) -> str:
    """Serialize a Bio.Phylo tree to Newick format."""
    handle = io.StringIO()
    Phylo.write(tree, handle, "newick")
    return handle.getvalue().strip()


def alignment_to_fasta(msa: MultipleSeqAlignment) -> str:
    """Export MultipleSeqAlignment as FASTA text."""
    handle = io.StringIO()
    SeqIO.write(msa, handle, "fasta")
    return handle.getvalue()


def plot_tree_matplotlib(tree, figsize: Tuple[int, int] = (12, 8)):
    """Render phylogenetic tree with matplotlib; returns figure."""
    fig = plt.figure(figsize=figsize, dpi=100)
    ax = fig.add_subplot(1, 1, 1)
    Phylo.draw(tree, axes=ax, do_show=False)
    ax.set_title("HIV-1 Phylogenetic Tree (Neighbor-Joining)", fontsize=12, fontweight="bold")
    ax.set_xlabel("Branch length (Jukes-Cantor distance)")
    fig.tight_layout()
    return fig


def plot_tree_plotly(tree):
    """
    Build an interactive Plotly phylogram from a Bio.Phylo tree.

    Returns a plotly Figure or None if plotly is unavailable.
    """
    try:
        import plotly.graph_objects as go
    except ImportError:
        return None

    depths = tree.depths(unit_branch_lengths=True)
    max_depth = max(depths.values()) if depths else 1.0
    y_positions: Dict[str, float] = {}
    leaves = tree.get_terminals()

    for i, leaf in enumerate(leaves):
        y_positions[leaf] = float(i)

    def _assign_internal(node, y_lo, y_hi):
        if node.is_terminal():
            return y_positions[node]
        children = node.clades
        if not children:
            return (y_lo + y_hi) / 2
        child_ys = []
        span = (y_hi - y_lo) / len(children)
        for idx, child in enumerate(children):
            cy_lo = y_lo + idx * span
            cy_hi = cy_lo + span
            child_ys.append(_assign_internal(child, cy_lo, cy_hi))
        y_positions[node] = sum(child_ys) / len(child_ys)
        return y_positions[node]

    if tree.root:
        _assign_internal(tree.root, 0, len(leaves))

    x_lines, y_lines = [], []

    def _add_edges(node, x_base=0.0):
        if node.is_terminal():
            return
        for child in node.clades:
            x_child = x_base + (depths.get(child, 0) - depths.get(node, 0))
            y_node = y_positions.get(node, 0)
            y_child = y_positions.get(child, 0)
            x_lines.extend([x_base, x_child, None])
            y_lines.extend([y_node, y_child, None])
            _add_edges(child, x_child)

    if tree.root:
        _add_edges(tree.root)

    fig = go.Figure()

    fig.add_trace(
        go.Scatter(
            x=x_lines,
            y=y_lines,
            mode="lines",
            line=dict(color="steelblue", width=1.5),
            hoverinfo="skip",
            showlegend=False,
        )
    )

    leaf_x = [depths.get(leaf, 0) for leaf in leaves]
    leaf_y = [y_positions.get(leaf, 0) for leaf in leaves]
    leaf_labels = [leaf.name or "Unknown" for leaf in leaves]

    fig.add_trace(
        go.Scatter(
            x=leaf_x,
            y=leaf_y,
            mode="markers+text",
            text=leaf_labels,
            textposition="middle right",
            marker=dict(size=6, color="crimson"),
            hovertext=leaf_labels,
            hoverinfo="text",
            showlegend=False,
        )
    )

    fig.update_layout(
        title="HIV-1 Phylogenetic Tree (Interactive NJ)",
        xaxis_title="Branch length (Jukes-Cantor)",
        yaxis=dict(showticklabels=False, showgrid=False),
        height=max(400, len(leaves) * 30),
        margin=dict(l=20, r=200, t=50, b=40),
        template="plotly_white",
    )

    return fig


def run_full_phylogenetics(
    df: pd.DataFrame,
    tree_method: str = "nj",
    progress_callback=None,
) -> Tuple[MultipleSeqAlignment, object, str, str]:
    """
    End-to-end: align sequences, build tree, return (msa, tree, fasta, newick).
    """
    def _cb(fraction, msg):
        if progress_callback:
            progress_callback(fraction * 0.6, msg)

    msa, _ = align_sequences(df, progress_callback=_cb)

    def _tree_cb(fraction, msg):
        if progress_callback:
            progress_callback(0.6 + fraction * 0.4, msg)

    tree = build_phylogenetic_tree(msa, method=tree_method, progress_callback=_tree_cb)
    fasta = alignment_to_fasta(msa)
    newick = tree_to_newick(tree)

    return msa, tree, fasta, newick
