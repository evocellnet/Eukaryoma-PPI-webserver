"""Optional external protein-association scores (coabundance, cofractionation,
phyloprofiling), combined with the AF3 pool data into one "universe" table
covering every possible pair among the website's proteins.

Each source ships as a dense protein x protein score matrix (CSV or
parquet). Row/column labels come in different formats:

- coabundance: FASTA-header style, e.g. "XP_004340666.2,4-aminobutyrate...".
- cofractionation: bare NCBI accession with no "XP_" prefix, e.g.
  "004340666.2", and some labels are semicolon-joined groups of
  indistinguishable proteins, e.g. "004340769.1;004346401.1" -- the group's
  row/column applies to every member id.
- phyloprofiling: FASTA-header-style accession plus a trailing "_<taxon id>"
  the other sources don't have, e.g. "XP_004349908.1_595528". Also, unlike
  the other two, this matrix is *not* symmetric (matrix[A, B] can differ
  from matrix[B, A]); we always take the value at (protein_lo, protein_hi)
  in alphabetical order, matching the analysis this app's methodology was
  validated against.

All normalize to the website's id format (see normalize_external_id).
Phyloprofiling had no file for a while; every function here still treats a
missing file as "this source has no data" rather than an error.
"""

import csv

import numpy as np
import pandas as pd

from eukaryoma_ppi.config import (
    COABUNDANCE_FILE,
    COFRACTIONATION_FILE,
    PHYLOPROFILING_FILE,
    UNIVERSE_INDEX_FILE,
)

# display name -> (file path, score column name, symmetrize-zero-gaps flag)
# phyloprofiling's matrix has each pair's score written in only one of the
# two (A,B)/(B,A) directions, leaving the other at the default 0 -- treated
# as a gap to fill from the mirror, not a real "no signal" measurement (see
# fill_zero_gaps_from_mirror). The other two sources are already symmetric
# correlation matrices where a literal 0.0 is a real value, so leave them as-is.
SOURCES = {
    "coabundance": (COABUNDANCE_FILE, "coabundance_score", False),
    "cofractionation": (COFRACTIONATION_FILE, "cofractionation_score", False),
    "phyloprofiling": (PHYLOPROFILING_FILE, "phyloprofiling_score", True),
}
SCORE_COLUMNS = [col for _path, col, _sym in SOURCES.values()]

# All four "is this pair backed by this source" columns, AF3 pooling included,
# in the order the recap page and unified score should present them.
ALL_SOURCE_LABELS = {
    "corrected_chain_pair_iptm": "AF3 pool structure",
    "coabundance_score": "Coabundance",
    "cofractionation_score": "Cofractionation",
    "phyloprofiling_score": "Phyloprofiling",
}


def normalize_external_id(raw_label):
    """"XP_004340666.2,annotation..." -> "xp0043406662"
    "004340666.2" -> "xp0043406662"
    "XP_004349908.1_595528" -> "xp0043499081" (trailing "_<taxon id>" stripped)
    """
    s = raw_label.split(",")[0].strip()
    parts = s.split("_")
    if len(parts) >= 3 and parts[0].lower() == "xp":
        s = "_".join(parts[:2])
    s = s.replace(".", "").replace("_", "").lower()
    if not s.startswith("xp"):
        s = "xp" + s
    return s


def _label_to_member_ids(label):
    """A matrix label may be a ";"-joined group of indistinguishable proteins."""
    raw_id_part = label.split(",")[0]
    return [normalize_external_id(member) for member in raw_id_part.split(";")]


def load_correlation_matrix(path):
    """Read a square protein x protein score matrix (CSV or parquet) as-is
    (original labels)."""
    if path.suffix == ".parquet":
        return pd.read_parquet(path)
    with open(path, newline="") as fh:
        header = next(csv.reader(fh))
    labels = header[1:]
    df = pd.read_csv(path, index_col=0, low_memory=False)
    df.index = labels  # guard against pandas mangling duplicate/odd label text
    return df


def fill_zero_gaps_from_mirror(sub):
    """For a matrix where a pair's score may have been written in only one of
    the (A,B)/(B,A) directions -- the other left at a default 0, not a real
    "no signal" measurement -- fill each 0 from its mirror position whenever
    the mirror is non-zero. Where both sides are non-zero they're expected to
    already agree (verified for phyloprofiling: 0 conflicts); this doesn't
    check that, it just prefers whichever side is non-zero.
    """
    return np.where(sub == 0, sub.T, sub)


def build_source_pairs(path, website_ids, score_col, symmetrize_zero_gaps=False):
    """One row per pair among website_ids with that source's score, or None if
    the source file doesn't exist. Pairs where either protein isn't covered by
    this source are simply absent (not NaN rows) -- the caller left-joins.

    covered is sorted alphabetically before extraction so that, for sources
    whose matrix isn't symmetric (phyloprofiling), we deterministically take
    matrix[protein_lo, protein_hi] rather than whichever of the two array
    positions happened to come first.
    """
    if not path.exists():
        return None

    matrix = load_correlation_matrix(path)

    # Map each website id to the matrix row/column that covers it (a group
    # label may cover several website ids at once).
    label_for_id = {}
    for label in matrix.index:
        for member_id in _label_to_member_ids(label):
            label_for_id.setdefault(member_id, label)

    covered = sorted(wid for wid in website_ids if wid in label_for_id)
    matrix_labels = [label_for_id[wid] for wid in covered]
    sub = matrix.reindex(index=matrix_labels, columns=matrix_labels).to_numpy()
    if symmetrize_zero_gaps:
        sub = fill_zero_gaps_from_mirror(sub)

    n = len(covered)
    iu = np.triu_indices(n, k=1)
    covered_arr = np.array(covered)

    return pd.DataFrame({"protein_a": covered_arr[iu[0]], "protein_b": covered_arr[iu[1]], score_col: sub[iu]})


def build_all_pairs_universe(website_ids):
    """Every possible pair among website_ids, as the (protein_a < protein_b) frame."""
    ids = np.array(sorted(website_ids))
    n = len(ids)
    iu = np.triu_indices(n, k=1)
    return pd.DataFrame({"protein_a": ids[iu[0]], "protein_b": ids[iu[1]]})


def compute_unified_score(df, score_columns):
    """Combined rank across whichever score_columns are present for a pair.

    Placeholder aggregation, easy to swap out later: each column is turned
    into a percentile rank (0-1, higher = better, computed only over the
    pairs that have a value in that column), then a pair's unified_score is
    the mean of its available per-column percentile ranks.
    """
    available = [c for c in score_columns if c in df.columns]
    percentiles = pd.DataFrame({c: df[c].rank(ascending=True, pct=True) for c in available})
    return percentiles.mean(axis=1, skipna=True)


def build_universe(website_ids, af3_pairs_df):
    """The full universe table: every pair among website_ids, with a score
    column per available source (AF3 pool iptm included) plus unified_score.
    """
    from eukaryoma_ppi.index import collapse_best_per_pair

    universe = build_all_pairs_universe(website_ids)

    af3_best = collapse_best_per_pair(af3_pairs_df)[
        ["protein_a", "protein_b", "pool", "chain_a", "chain_b", "corrected_chain_pair_iptm"]
    ]
    universe = universe.merge(af3_best, on=["protein_a", "protein_b"], how="left")

    for source_name, (path, score_col, symmetrize_zero_gaps) in SOURCES.items():
        source_pairs = build_source_pairs(path, website_ids, score_col, symmetrize_zero_gaps)
        if source_pairs is None:
            universe[score_col] = float("nan")
            continue
        universe = universe.merge(source_pairs, on=["protein_a", "protein_b"], how="left")

    score_cols = ["corrected_chain_pair_iptm"] + SCORE_COLUMNS
    universe["unified_score"] = compute_unified_score(universe, score_cols)
    return universe


def add_presence_columns(df):
    """Add n_sources_present (int) and sources_present (str label) columns,
    describing which of ALL_SOURCE_LABELS' sources back each pair. Fully
    vectorized (no per-row Python loop) so it stays fast at ~2.3M rows: each
    row's presence pattern is packed into a bitmask, then only the handful of
    *distinct* bitmasks actually observed (at most 2**4) are turned into
    "A + B + C" labels and mapped back.
    """
    labels = {col: name for col, name in ALL_SOURCE_LABELS.items() if col in df.columns}
    names = list(labels.values())
    present = np.column_stack([df[col].notna().to_numpy() for col in labels])

    df = df.copy()
    df["n_sources_present"] = present.sum(axis=1)

    weights = 1 << np.arange(len(names))
    bitmask = present.astype(np.int64) @ weights
    label_for_bitmask = {
        b: (" + ".join(name for i, name in enumerate(names) if b & (1 << i)) or "none") for b in np.unique(bitmask)
    }
    df["sources_present"] = pd.Series(bitmask).map(label_for_bitmask).to_numpy()
    return df


def source_presence_summary(universe_df):
    """Cross-source coverage of the pair universe, for the recap page.

    Returns (pattern_counts, n_sources_counts):
    - pattern_counts: one row per distinct combination of sources present
      (e.g. "AF3 pool structure + Cofractionation"), with the pair count,
      sorted by count descending.
    - n_sources_counts: pair count by how many of the (up to 4) sources are
      present, indexed 0..len(labels).
    """
    df = add_presence_columns(universe_df)
    n_sources_counts = df["n_sources_present"].value_counts().sort_index()
    pattern_counts = (
        df["sources_present"]
        .value_counts()
        .rename_axis("sources_present")
        .reset_index(name="n_pairs")
        .sort_values("n_pairs", ascending=False)
        .reset_index(drop=True)
    )
    return pattern_counts, n_sources_counts


def save_universe(universe_df):
    UNIVERSE_INDEX_FILE.parent.mkdir(parents=True, exist_ok=True)
    universe_df.to_parquet(UNIVERSE_INDEX_FILE, index=False)


def load_universe():
    return pd.read_parquet(UNIVERSE_INDEX_FILE)
