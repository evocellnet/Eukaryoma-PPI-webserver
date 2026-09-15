import streamlit as st

from eukaryoma_ppi import annotations, external_scores, ui
from eukaryoma_ppi.config import ANNOTATIONS_INDEX_FILE, UNIVERSE_INDEX_FILE

st.set_page_config(page_title="Coabundance - Eukaryoma PPI", page_icon="🧬", layout="wide")
st.title("Coabundance")
st.write(
    "Protein pairs ranked by coabundance correlation score. Other sources' scores are shown when "
    "available; select a row to view its AF3 structure, if one has been predicted."
)

if not UNIVERSE_INDEX_FILE.exists():
    st.error("No universe index found. Run `python scripts/build_data.py` first (see README.md).")
    st.stop()

SCORE_COL = "coabundance_score"


@st.cache_data
def get_coabundance_pairs():
    universe_df = external_scores.load_universe()
    df = universe_df[universe_df[SCORE_COL].notna()].copy()
    if ANNOTATIONS_INDEX_FILE.exists():
        df = annotations.attach_annotations(df)
    else:
        df["annotation_a"] = ""
        df["annotation_b"] = ""
    return df


scored_pairs = get_coabundance_pairs()
if scored_pairs.empty:
    st.info("No coabundance data available.")
    st.stop()

st.caption(f"{len(scored_pairs):,} pairs have a coabundance score.")

other_scores = [
    (col, label)
    for col, label in external_scores.ALL_SOURCE_LABELS.items()
    if col != "corrected_chain_pair_iptm"
]

ui.render_ranked_pair_table(
    scored_pairs,
    sort_col=SCORE_COL,
    column_order=[SCORE_COL, "protein_a", "annotation_a", "protein_b", "annotation_b", "corrected_chain_pair_iptm"]
    + [col for col, _ in other_scores if col != SCORE_COL],
    column_config={
        SCORE_COL: st.column_config.NumberColumn("Coabundance score", format="%.3f"),
        "protein_a": st.column_config.TextColumn("Protein A"),
        "annotation_a": st.column_config.TextColumn("Protein A annotation"),
        "protein_b": st.column_config.TextColumn("Protein B"),
        "annotation_b": st.column_config.TextColumn("Protein B annotation"),
        "corrected_chain_pair_iptm": st.column_config.NumberColumn("AF3 pool ipTM", format="%.3f"),
        "cofractionation_score": st.column_config.NumberColumn("Cofractionation score", format="%.3f"),
        "phyloprofiling_score": st.column_config.NumberColumn("Phyloprofiling score", format="%.3f"),
    },
    score_fields=[(SCORE_COL, "Coabundance score")]
    + [("corrected_chain_pair_iptm", "AF3 pool ipTM")]
    + [(col, label) for col, label in other_scores if col != SCORE_COL],
    key_prefix="coab_",
)
