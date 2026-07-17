import streamlit as st

from eukaryoma_ppi import index
from eukaryoma_ppi.config import PAIRS_INDEX_FILE, POOLS_INDEX_FILE

st.set_page_config(page_title="Eukaryoma PPI webserver", page_icon="🧬", layout="wide")

st.title("Eukaryoma PPI webserver")
st.write(
    "Browse pairwise protein-protein interactions extracted from pooled AlphaFold3 "
    "structure predictions. Each pool jointly predicts several proteins at once; this "
    "tool pulls out the two chains for any specific pair and renders their predicted "
    "3D interaction."
)

if not PAIRS_INDEX_FILE.exists() or not POOLS_INDEX_FILE.exists():
    st.error(
        "No data index found. Run `python scripts/build_data.py` first to decompress the "
        "pools and build the pair index (see README.md)."
    )
    st.stop()

pools_df = index.load_pools_index()
pairs_df = index.load_pairs_index()

n_unique_proteins = len(set(pairs_df["protein_a"]) | set(pairs_df["protein_b"]))
n_unique_pairs = len(pairs_df.drop_duplicates(subset=["protein_a", "protein_b"]))

col1, col2, col3 = st.columns(3)
col1.metric("Pools", f"{len(pools_df):,}")
col2.metric("Unique proteins", f"{n_unique_proteins:,}")
col3.metric("Protein pairs", f"{n_unique_pairs:,}")

st.info(
    "Use **Browse Pairs** in the sidebar to rank pairs by predicted interaction score, "
    "or **Pair Viewer** to look up a specific pair directly."
)
