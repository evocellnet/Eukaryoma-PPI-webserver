import streamlit as st

from eukaryoma_ppi import index, ui
from eukaryoma_ppi.config import PAIRS_INDEX_FILE, POOLS_INDEX_FILE

st.set_page_config(page_title="Browse Pairs - Eukaryoma PPI", page_icon="🧬", layout="wide")
st.title("Browse Protein Pairs")
st.write("All predicted pairs, ranked by predicted interaction score. Select a row to view its structure.")

if not PAIRS_INDEX_FILE.exists() or not POOLS_INDEX_FILE.exists():
    st.error("No data index found. Run `python scripts/build_data.py` first (see README.md).")
    st.stop()


@st.cache_data
def get_ranked_pairs():
    pairs_df = index.load_pairs_index()
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
    column_order=["protein_a", "protein_b", "corrected_chain_pair_iptm"],
    column_config={
        "protein_a": st.column_config.TextColumn("Protein A"),
        "protein_b": st.column_config.TextColumn("Protein B"),
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
