"""Shared Streamlit rendering used by more than one page."""

import pandas as pd
import streamlit as st
import streamlit.components.v1 as components

from eukaryoma_ppi import structures, viewer


@st.cache_data
def get_pdb_text(pool, chain_a, chain_b):
    return structures.extract_chains_as_pdb(pool, [chain_a, chain_b])


def render_pair_detail(row, protein_a, protein_b, key_prefix=""):
    """Render the score, color-by control, and 3D viewer for one pair row."""
    score = row.get("corrected_chain_pair_iptm")
    st.metric(
        "Predicted interaction score (corrected ipTM)",
        f"{score:.3f}" if pd.notna(score) else "n/a",
    )

    color_by = st.radio("Color by", ["chain", "plddt"], horizontal=True, key=f"{key_prefix}color_by")

    # row["chain_a"]/["chain_b"] follow the pool's own chain order, which may
    # not match the order the two proteins were passed in -- resolve by identity.
    chain_for = {row["protein_a"]: row["chain_a"], row["protein_b"]: row["chain_b"]}
    chain_a, chain_b = chain_for[protein_a], chain_for[protein_b]

    pdb_text = get_pdb_text(row["pool"], chain_a, chain_b)
    html = viewer.render_pair(pdb_text, chain_a, chain_b, color_by=color_by)

    st.caption(f"Pool **{row['pool']}** &mdash; chain {chain_a} = {protein_a}, chain {chain_b} = {protein_b}")
    components.html(html, height=580)
