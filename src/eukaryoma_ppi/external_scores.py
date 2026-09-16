"""Optional external protein-association scores (coabundance, cofractionation,
phyloprofiling), combined with the AF3 pool data into one "universe" table
covering every possible pair among the website's proteins.

Each source ships as a dense protein x protein correlation matrix CSV. Row/
column labels come in different formats:

- coabundance: FASTA-header style, e.g. "XP_004340666.2,4-aminobutyrate...".
- cofractionation: bare NCBI accession with no "XP_" prefix, e.g.
  "004340666.2", and some labels are semicolon-joined groups of
  indistinguishable proteins, e.g. "004340769.1;004346401.1" -- the group's
  row/column applies to every member id.

Both normalize to the website's id format (see normalize_external_id).
Phyloprofiling has no file yet; every function here treats a missing file as
"this source has no data" rather than an error.
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

# display name -> (file path, score column name)
SOURCES = {
    "coabundance": (COABUNDANCE_FILE, "coabundance_score"),
    "cofractionation": (COFRACTIONATION_FILE, "cofractionation_score"),
    "phyloprofiling": (PHYLOPROFILING_FILE, "phyloprofiling_score"),
}
SCORE_COLUMNS = [col for _path, col in SOURCES.values()]

# All four "is this pair backed by this source" columns, AF3 pooling included,
# in the order the recap page and unified score should present them.
ALL_SOURCE_LABELS = {
    "corrected_chain_pair_iptm": "AF3 pool structure",
    "coabundance_score": "Coabundance",
    "cofractionation_score": "Cofractionation",
    "phyloprofiling_score": "Phyloprofiling",
}


def normalize_external_id(raw_label):
    """"XP_004340666.2,annotation..." or "004340666.2" -> "xp0043406662"."""
    s = raw_label.split(",")[0].strip()
    s = s.replace(".", "").replace("_", "").lower()
    if not s.startswith("xp"):
        s = "xp" + s
    return s


def _label_to_member_ids(label):
    """A matrix label may be a ";"-joined group of indistinguishable proteins."""
    raw_id_part = label.split(",")[0]
    return [normalize_external_id(member) for member in raw_id_part.split(";")]


def load_correlation_matrix(path):
    """Read a square protein x protein correlation matrix CSV as-is (original labels)."""
    with open(path, newline="") as fh:
        header = next(csv.reader(fh))
    labels = header[1:]
    df = pd.read_csv(path, index_col=0, low_memory=False)
    df.index = labels  # guard against pandas mangling duplicate/odd label text
    return df


def build_source_pairs(path, website_ids, score_col):
    """One row per pair among website_ids with that source's score, or None if
    the source file doesn't exist. Pairs where either protein isn't covered by
    this source are simply absent (not NaN rows) -- the caller left-joins.
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

    covered = [wid for wid in website_ids if wid in label_for_id]
    matrix_labels = [label_for_id[wid] for wid in covered]
    sub = matrix.reindex(index=matrix_labels, columns=matrix_labels).to_numpy()

    n = len(covered)
    iu = np.triu_indices(n, k=1)
    covered_arr = np.array(covered)
    left, right = covered_arr[iu[0]], covered_arr[iu[1]]
    is_left_smaller = left <= right
    protein_lo = np.where(is_left_smaller, left, right)
    protein_hi = np.where(is_left_smaller, right, left)

    return pd.DataFrame({"protein_a": protein_lo, "protein_b": protein_hi, score_col: sub[iu]})


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

    for source_name, (path, score_col) in SOURCES.items():
        source_pairs = build_source_pairs(path, website_ids, score_col)
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
