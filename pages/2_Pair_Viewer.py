import streamlit as st

from eukaryoma_ppi import index, ui
from eukaryoma_ppi.config import PAIRS_INDEX_FILE, POOLS_INDEX_FILE

st.set_page_config(page_title="Pair Viewer - Eukaryoma PPI", page_icon="🧬", layout="wide")
st.title("Protein Pair Viewer")

if not PAIRS_INDEX_FILE.exists() or not POOLS_INDEX_FILE.exists():
    st.error("No data index found. Run `python scripts/build_data.py` first (see README.md).")
    st.stop()


@st.cache_data
def get_pairs_index():
    return index.load_pairs_index()


@st.cache_data
def get_protein_options(_pairs_df):
    return sorted(set(_pairs_df["protein_a"]) | set(_pairs_df["protein_b"]))


pairs_df = get_pairs_index()
protein_options = get_protein_options(pairs_df)

col_a, col_b = st.columns(2)
protein_a = col_a.selectbox("Protein A", protein_options, index=None, placeholder="Select or type an ID")
protein_b = col_b.selectbox("Protein B", protein_options, index=None, placeholder="Select or type an ID")

if not protein_a or not protein_b:
    st.info("Pick two proteins to look up their predicted interaction.")
    st.stop()

if protein_a == protein_b:
    st.warning("Pick two different proteins.")
    st.stop()

matches = index.lookup_pair(pairs_df, protein_a, protein_b)

if matches.empty:
    st.warning(f"No pool jointly predicted **{protein_a}** and **{protein_b}** together.")
    st.stop()

if len(matches) > 1:
    st.caption(f"This pair co-occurs in {len(matches)} pools; pick one below.")
    pool_choice = st.selectbox("Pool", matches["pool"].tolist())
    row = matches[matches["pool"] == pool_choice].iloc[0]
else:
    row = matches.iloc[0]

ui.render_pair_detail(row, protein_a, protein_b)
