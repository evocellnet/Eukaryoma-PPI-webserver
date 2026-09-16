import streamlit as st

from eukaryoma_ppi import annotations, external_scores, ui
from eukaryoma_ppi.config import ANNOTATIONS_INDEX_FILE, UNIVERSE_INDEX_FILE

st.set_page_config(page_title="Unified Ranking - Eukaryoma PPI", page_icon="🧬", layout="wide")
st.title("Unified Ranking")
st.write(
    "Every pair with a score from at least one source, ranked by a unified score. For now the unified "
    "score is the mean percentile rank across whichever sources have data for that pair (see "
    "`eukaryoma_ppi.external_scores.compute_unified_score`) -- a placeholder that can be swapped for a "
    "different combination method later. Select a row to view its AF3 structure, if one has been predicted."
)

if not UNIVERSE_INDEX_FILE.exists():
    st.error("No universe index found. Run `python scripts/build_data.py` first (see README.md).")
    st.stop()


@st.cache_data
def get_universe_for_ranking():
    universe_df = external_scores.load_universe()
    if ANNOTATIONS_INDEX_FILE.exists():
        universe_df = annotations.attach_annotations(universe_df)
    else:
        universe_df["annotation_a"] = ""
        universe_df["annotation_b"] = ""
    return external_scores.add_presence_columns(universe_df)


@st.cache_data
def get_pattern_counts(_universe_df):
    pattern_counts, _n_sources_counts = external_scores.source_presence_summary(_universe_df)
    return pattern_counts


universe_df = get_universe_for_ranking()
st.caption(f"{len(universe_df):,} possible pairs among all website proteins.")

max_sources = int(universe_df["n_sources_present"].max())
pattern_counts = get_pattern_counts(universe_df)
counts_by_pattern = pattern_counts.set_index("sources_present")["n_pairs"]

st.subheader("Filter by data source coverage")
filter_mode = st.radio(
    "Filter by data source coverage",
    ["All pairs", "Minimum number of sources", "Exact combination of sources"],
    horizontal=True,
    label_visibility="collapsed",
)

if filter_mode == "Minimum number of sources":
    min_sources = st.selectbox("At least this many sources present", list(range(1, max_sources + 1)))
    universe_df = universe_df[universe_df["n_sources_present"] >= min_sources]
elif filter_mode == "Exact combination of sources":
    combo_options = pattern_counts["sources_present"].tolist()
    selected_combos = st.multiselect(
        "Sources present",
        combo_options,
        default=combo_options,
        format_func=lambda combo: f"{combo} ({counts_by_pattern[combo]:,})",
    )
    universe_df = universe_df[universe_df["sources_present"].isin(selected_combos)]

if universe_df.empty:
    st.info("No pairs match this filter.")
    st.stop()

if filter_mode != "All pairs":
    st.caption(f"{len(universe_df):,} pairs match this filter.")

score_cols = list(external_scores.ALL_SOURCE_LABELS.items())

ui.render_ranked_pair_table(
    universe_df,
    sort_col="unified_score",
    column_order=["unified_score", "protein_a", "annotation_a", "protein_b", "annotation_b"]
    + [col for col, _ in score_cols],
    column_config={
        "unified_score": st.column_config.NumberColumn("Unified score", format="%.3f"),
        "protein_a": st.column_config.TextColumn("Protein A"),
        "annotation_a": st.column_config.TextColumn("Protein A annotation"),
        "protein_b": st.column_config.TextColumn("Protein B"),
        "annotation_b": st.column_config.TextColumn("Protein B annotation"),
        "corrected_chain_pair_iptm": st.column_config.NumberColumn("AF3 pool ipTM", format="%.3f"),
        "coabundance_score": st.column_config.NumberColumn("Coabundance score", format="%.3f"),
        "cofractionation_score": st.column_config.NumberColumn("Cofractionation score", format="%.3f"),
        "phyloprofiling_score": st.column_config.NumberColumn("Phyloprofiling score", format="%.3f"),
    },
    score_fields=[("unified_score", "Unified score")] + score_cols,
    key_prefix="unified_",
)
