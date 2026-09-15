import streamlit as st

from eukaryoma_ppi import annotations, index, ui
from eukaryoma_ppi.config import ANNOTATIONS_INDEX_FILE, PAIRS_INDEX_FILE, POOLS_INDEX_FILE

st.set_page_config(page_title="Browse Pairs - Eukaryoma PPI", page_icon="🧬", layout="wide")
st.title("Browse Protein Pairs")
st.write("All AF3-pooled pairs, ranked by predicted interaction score. Select a row to view its structure.")

if not PAIRS_INDEX_FILE.exists() or not POOLS_INDEX_FILE.exists():
    st.error("No data index found. Run `python scripts/build_data.py` first (see README.md).")
    st.stop()


@st.cache_data
def get_ranked_pairs():
    pairs_df = index.load_pairs_index()
    if ANNOTATIONS_INDEX_FILE.exists():
        pairs_df = annotations.attach_annotations(pairs_df)
    else:
        pairs_df["annotation_a"] = ""
        pairs_df["annotation_b"] = ""
    return pairs_df


ranked_pairs = get_ranked_pairs()

ui.render_ranked_pair_table(
    ranked_pairs,
    sort_col="corrected_chain_pair_iptm",
    column_order=["corrected_chain_pair_iptm", "protein_a", "annotation_a", "protein_b", "annotation_b"],
    column_config={
        "corrected_chain_pair_iptm": st.column_config.NumberColumn("Corrected ipTM", format="%.3f"),
        "protein_a": st.column_config.TextColumn("Protein A"),
        "annotation_a": st.column_config.TextColumn("Protein A annotation"),
        "protein_b": st.column_config.TextColumn("Protein B"),
        "annotation_b": st.column_config.TextColumn("Protein B annotation"),
    },
    score_fields=[("corrected_chain_pair_iptm", "Predicted interaction score (corrected ipTM)")],
    key_prefix="browse_",
)
