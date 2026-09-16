"""Shared Streamlit rendering used by more than one page."""

import pandas as pd
import streamlit as st
import streamlit.components.v1 as components

from eukaryoma_ppi import structures, viewer

# Row highlight for pairs flagged true-positive by CORUM/Marcotte complex
# co-membership (see eukaryoma_ppi.complex_annotations). Both present wins.
TP_COLUMNS = ["corum_tp", "marcotte_tp"]
TP_ROW_STYLES = {
    (True, True): "background-color: #4a2f66; color: #E4C9FF",  # both
    (True, False): "background-color: #1e5c3a; color: #7CFFB2",  # corum only
    (False, True): "background-color: #1e4a6b; color: #8FD3FF",  # marcotte only
}
TP_LEGEND = (
    "Highlighted rows are pairs flagged as a known true-positive interaction by complex "
    "co-membership: \U0001F7E9 CORUM, \U0001F7E6 Marcotte, \U0001F7EA both."
)


@st.cache_data
def get_pdb_text(pool, chain_a, chain_b):
    return structures.extract_chains_as_pdb(pool, [chain_a, chain_b])


def style_true_positive_rows(df):
    """Highlight rows flagged true-positive by CORUM/Marcotte, if those
    columns are present; otherwise return df unchanged.
    """
    if not any(col in df.columns for col in TP_COLUMNS):
        return df

    def highlight(row):
        key = (bool(row.get("corum_tp", False)), bool(row.get("marcotte_tp", False)))
        style = TP_ROW_STYLES.get(key, "")
        return [style] * len(row)

    return df.style.apply(highlight, axis=1)


def render_true_positive_badges(row):
    """Small caption noting which annotation source(s) flag this pair, if any."""
    corum, marcotte = bool(row.get("corum_tp", False)), bool(row.get("marcotte_tp", False))
    if corum and marcotte:
        st.success("Known true-positive interaction: flagged by **both CORUM and Marcotte**.")
    elif corum:
        st.success("Known true-positive interaction: flagged by **CORUM** complex co-membership.")
    elif marcotte:
        st.success("Known true-positive interaction: flagged by **Marcotte** complex co-membership.")


def render_scores(row, score_fields):
    """One st.metric per (column_name, display_label) in score_fields."""
    cols = st.columns(len(score_fields))
    for col, (score_col, label) in zip(cols, score_fields):
        value = row.get(score_col)
        col.metric(label, f"{value:.3f}" if pd.notna(value) else "n/a")


def render_structure_if_available(row, protein_a, protein_b, key_prefix=""):
    """3D viewer if this pair has an AF3 pool structure, else a plain notice."""
    pool = row.get("pool")
    if pd.isna(pool):
        st.info("No AF3 structure has been predicted for this pair.")
        return

    color_by = st.radio("Color by", ["chain", "plddt"], horizontal=True, key=f"{key_prefix}color_by")

    # row["chain_a"]/["chain_b"] follow the pool's own chain order, which may
    # not match the order the two proteins were passed in -- resolve by identity.
    chain_for = {row["protein_a"]: row["chain_a"], row["protein_b"]: row["chain_b"]}
    chain_a, chain_b = chain_for[protein_a], chain_for[protein_b]

    pdb_text = get_pdb_text(pool, chain_a, chain_b)
    html = viewer.render_pair(pdb_text, chain_a, chain_b, color_by=color_by)

    st.caption(f"Pool **{pool}** &mdash; chain {chain_a} = {protein_a}, chain {chain_b} = {protein_b}")
    components.html(html, height=580)


def render_pair_detail(row, protein_a, protein_b, key_prefix=""):
    """Score + color-by control + 3D viewer for one AF3 pool pair row."""
    render_scores(row, [("corrected_chain_pair_iptm", "Predicted interaction score (corrected ipTM)")])
    render_true_positive_badges(row)
    render_structure_if_available(row, protein_a, protein_b, key_prefix=key_prefix)


def render_ranked_pair_table(
    df, sort_col, column_order, column_config, score_fields, default_top_n=500, key_prefix="", sort_ascending=False
):
    """Sortable/selectable pair table; selecting a row shows its scores and,
    when available, its 3D structure.

    df must have protein_a/protein_b columns, and pool/chain_a/chain_b for
    rows that have a predicted structure (NaN otherwise).
    score_fields: [(column_name, display_label), ...] shown as metrics.
    """
    if df.empty:
        st.info("No pairs match this filter.")
        return

    ranked = df.sort_values(sort_col, ascending=sort_ascending).reset_index(drop=True)
    top_n = st.number_input(
        "Show top N pairs",
        min_value=min(10, len(ranked)),
        max_value=len(ranked),
        value=min(default_top_n, len(ranked)),
        step=10,
        key=f"{key_prefix}topn",
    )
    displayed = ranked.head(top_n)

    if any(col in displayed.columns for col in TP_COLUMNS):
        st.caption(TP_LEGEND)

    event = st.dataframe(
        style_true_positive_rows(displayed),
        column_order=column_order,
        column_config=column_config,
        hide_index=True,
        width="stretch",
        on_select="rerun",
        selection_mode="single-row",
        key=f"{key_prefix}table",
    )

    selected_rows = event.selection.rows
    if not selected_rows:
        st.info("Select a row above to view that pair's scores and structure.")
        return

    row = displayed.iloc[selected_rows[0]]
    st.divider()
    st.subheader(f"{row['protein_a']} — {row['protein_b']}")
    render_scores(row, score_fields)
    render_true_positive_badges(row)
    render_structure_if_available(row, row["protein_a"], row["protein_b"], key_prefix=key_prefix)
