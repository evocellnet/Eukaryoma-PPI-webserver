import streamlit as st

from eukaryoma_ppi import annotations, index, ui
from eukaryoma_ppi.config import ANNOTATIONS_INDEX_FILE, PAIRS_INDEX_FILE, POOLS_INDEX_FILE

st.set_page_config(page_title="Browse Pairs - Eukaryoma PPI", page_icon="🧬", layout="wide")
st.title("Browse Protein Pairs")
st.write("All predicted pairs, ranked by predicted interaction score. Select a row to view its structure.")

if not PAIRS_INDEX_FILE.exists() or not POOLS_INDEX_FILE.exists():
    st.error("No data index found. Run `python scripts/build_data.py` first (see README.md).")
    st.stop()


@st.cache_data
def get_ranked_pairs():
    pairs_df = index.load_pairs_index()

    if ANNOTATIONS_INDEX_FILE.exists():
        protein_annotations = annotations.load_annotations()[["protein_id", "annotation"]]
        pairs_df = pairs_df.merge(
            protein_annotations.rename(columns={"protein_id": "protein_a", "annotation": "annotation_a"}),
            on="protein_a",
            how="left",
        ).merge(
            protein_annotations.rename(columns={"protein_id": "protein_b", "annotation": "annotation_b"}),
            on="protein_b",
            how="left",
        )
    else:
        pairs_df["annotation_a"] = ""
        pairs_df["annotation_b"] = ""

    return pairs_df.sort_values("corrected_chain_pair_iptm", ascending=False).reset_index(drop=True)


ranked_pairs = get_ranked_pairs()

top_n = st.number_input(
    "Show top N pairs",
    min_value=10,
    max_value=len(ranked_pairs),
    value=min(500, len(ranked_pairs)),
    step=10,
)
displayed = ranked_pairs.head(top_n)

event = st.dataframe(
    displayed,
    column_order=["protein_a", "annotation_a", "protein_b", "annotation_b", "corrected_chain_pair_iptm"],
    column_config={
        "protein_a": st.column_config.TextColumn("Protein A"),
        "annotation_a": st.column_config.TextColumn("Protein A annotation"),
        "protein_b": st.column_config.TextColumn("Protein B"),
        "annotation_b": st.column_config.TextColumn("Protein B annotation"),
        "corrected_chain_pair_iptm": st.column_config.NumberColumn("Corrected ipTM", format="%.3f"),
    },
    hide_index=True,
    width="stretch",
    on_select="rerun",
    selection_mode="single-row",
)

selected_rows = event.selection.rows
if not selected_rows:
    st.info("Select a row above to view that pair's structure.")
    st.stop()

row = displayed.iloc[selected_rows[0]]
st.divider()
st.subheader(f"{row['protein_a']} &mdash; {row['protein_b']}")
ui.render_pair_detail(row, row["protein_a"], row["protein_b"], key_prefix="browse_")
