import streamlit as st

from eukaryoma_ppi import external_scores, index, ui
from eukaryoma_ppi.config import PAIRS_INDEX_FILE, POOLS_INDEX_FILE, UNIVERSE_INDEX_FILE

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
def get_universe_for_pair_viewer():
    return external_scores.load_universe()


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

score_fields = [("corrected_chain_pair_iptm", "AF3 pool ipTM (this pool)")]

if UNIVERSE_INDEX_FILE.exists():
    protein_lo, protein_hi = min(protein_a, protein_b), max(protein_a, protein_b)
    universe_df = get_universe_for_pair_viewer()
    universe_match = universe_df[
        (universe_df["protein_a"] == protein_lo) & (universe_df["protein_b"] == protein_hi)
    ]
    if not universe_match.empty:
        universe_row = universe_match.iloc[0]
        for col, label in external_scores.ALL_SOURCE_LABELS.items():
            if col != "corrected_chain_pair_iptm":
                row[col] = universe_row[col]
                score_fields.append((col, label))

ui.render_scores(row, score_fields)
ui.render_true_positive_badges(row)
ui.render_structure_if_available(row, protein_a, protein_b)
