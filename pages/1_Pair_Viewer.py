import streamlit.components.v1 as components
import streamlit as st

from eukaryoma_ppi import index, structures, viewer
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


@st.cache_data
def get_pdb_text(pool, chain_a, chain_b):
    return structures.extract_chains_as_pdb(pool, [chain_a, chain_b])


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

color_by = st.radio("Color by", ["chain", "plddt"], horizontal=True)

# row["chain_a"]/["chain_b"] follow the pool's own chain order, which may not
# match the order the user picked proteins in -- resolve by protein identity.
chain_for = {row["protein_a"]: row["chain_a"], row["protein_b"]: row["chain_b"]}
chain_a, chain_b = chain_for[protein_a], chain_for[protein_b]

pdb_text = get_pdb_text(row["pool"], chain_a, chain_b)
html = viewer.render_pair(pdb_text, chain_a, chain_b, color_by=color_by)

st.caption(f"Pool **{row['pool']}** &mdash; chain {chain_a} = {protein_a}, chain {chain_b} = {protein_b}")
components.html(html, height=580)
