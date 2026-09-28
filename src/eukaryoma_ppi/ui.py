"""Shared Streamlit rendering used by more than one page."""

import pandas as pd
import streamlit as st
import streamlit.components.v1 as components

from eukaryoma_ppi import external_scores, structures, viewer

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

# Unified score combination mode, shared across every page via
# st.session_state -- set once by unified_score_mode_selector (Unified
# Ranking page), read everywhere else by get_unified_score_mode(). This is
# what keeps a mode picked on one page in sync on every other page without
# each page needing its own copy of the control.
#
# Two keys, not one: Streamlit clears a *widget's* session_state entry
# whenever that widget isn't instantiated on a script run -- which happens
# on every page that isn't Unified Ranking, since only that page calls the
# selectbox below. So the widget's own key (UNIFIED_SCORE_MODE_KEY) cannot
# be trusted to survive navigation; UNIFIED_SCORE_MODE_PERSIST_KEY is a
# plain (non-widget) session_state entry that does, and is what every page
# actually reads.
UNIFIED_SCORE_MODE_KEY = "unified_score_mode"
UNIFIED_SCORE_MODE_PERSIST_KEY = "unified_score_mode_persisted"
UNIFIED_SCORE_MODES = {
    "weightless": "Weightless (equal-weighted mean rank)",
    "weighted": "Weighted (default, by score distribution)",
}
# User-edited per-source weights for "weighted" mode, same persistence
# pattern as the mode above: a plain session_state entry (survives
# navigation) plus one widget key per source slider (doesn't).
UNIFIED_SCORE_WEIGHTS_PERSIST_KEY = "unified_score_weights_persisted"


def get_unified_score_mode():
    """Active unified-score mode. Defaults to "weightless" (the original,
    equal-weighted behaviour) until the user picks "weighted" on the Unified
    Ranking page.
    """
    return st.session_state.get(UNIFIED_SCORE_MODE_PERSIST_KEY, "weightless")


def _persist_unified_score_mode():
    st.session_state[UNIFIED_SCORE_MODE_PERSIST_KEY] = st.session_state[UNIFIED_SCORE_MODE_KEY]


def unified_score_mode_selector():
    """Selectbox for the unified-score combination mode. Meant to be shown
    once, on the Unified Ranking page -- every other page just reads the
    choice back via get_unified_score_mode().

    Seeds the widget's own session_state entry from the persisted mode
    before instantiating it (rather than passing `index=`), so the widget
    doesn't just spring back to the default the first time this page is
    revisited (see the module-level note above).
    """
    if UNIFIED_SCORE_MODE_KEY not in st.session_state:
        st.session_state[UNIFIED_SCORE_MODE_KEY] = get_unified_score_mode()
    st.selectbox(
        "Unified score type",
        list(UNIFIED_SCORE_MODES.keys()),
        format_func=lambda mode: UNIFIED_SCORE_MODES[mode],
        key=UNIFIED_SCORE_MODE_KEY,
        on_change=_persist_unified_score_mode,
    )
    return st.session_state[UNIFIED_SCORE_MODE_KEY]


def get_active_unified_score_weights(default_weights):
    """Weights to use in "weighted" mode: the user's last edited sliders, if
    they've touched them this session, otherwise `default_weights` as-is
    (computed fresh from the data by the caller, since that's cheap -- just
    a per-column std -- compared to the rank/combine step these feed into).
    """
    return st.session_state.get(UNIFIED_SCORE_WEIGHTS_PERSIST_KEY, default_weights)


def unified_score_weight_sliders(default_weights, labels):
    """One slider per source, defaulting to `default_weights` (or the user's
    own values from earlier this session). Meant to be shown once, on the
    Unified Ranking page, right under unified_score_mode_selector -- other
    pages read the result back via get_active_unified_score_weights().

    Weights don't need to sum to 1 -- compute_unified_score renormalizes
    per row over whichever sources that pair actually has -- but the
    sliders start at values that do, and a "reset" button gets back there.

    Returns the resulting {column: weight} dict.
    """
    active = get_active_unified_score_weights(default_weights)
    columns = st.columns(len(default_weights))
    weights = {}
    for col_widget, (col, default) in zip(columns, default_weights.items()):
        widget_key = f"unified_score_weight__{col}"
        if widget_key not in st.session_state:
            st.session_state[widget_key] = active.get(col, default)
        weights[col] = col_widget.slider(labels.get(col, col), 0.0, 1.0, step=0.01, key=widget_key)
    st.session_state[UNIFIED_SCORE_WEIGHTS_PERSIST_KEY] = weights

    if st.button("Reset to computed defaults"):
        for col in default_weights:
            st.session_state.pop(f"unified_score_weight__{col}", None)
        st.session_state.pop(UNIFIED_SCORE_WEIGHTS_PERSIST_KEY, None)
        st.rerun()

    return weights


def recompute_unified_score(df, score_columns, mode, weights=None):
    """Apply `mode` to df's unified_score column.

    Pure function of (df, mode, weights) so callers can wrap it in
    st.cache_data keyed on those alone (df should come from an already
    st.cache_data'd loader and be passed with a leading-underscore parameter
    name so Streamlit doesn't re-hash the whole ~2.3M-row table on every
    call) -- this way the rank/weight recompute (sub-second, but not free)
    only happens once per distinct (mode, weights), not on every widget
    interaction.

    weights: explicit {column: weight} for "weighted" mode (e.g. the user's
    edited sliders); if None, computed fresh from df's own distribution.
    Ignored in "weightless" mode.

    Returns (df, weights_used): weights_used is None in "weightless" mode.
    """
    if mode != "weighted":
        return df, None
    if weights is None:
        weights = external_scores.compute_default_weights(df, score_columns)
    df = df.copy()
    df["unified_score"] = external_scores.compute_unified_score(df, score_columns, weights=weights)
    return df, weights


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
    st.download_button(
        "Download structure (.pdb)",
        data=pdb_text,
        file_name=f"{pool}_{protein_a}_{protein_b}.pdb",
        mime="chemical/x-pdb",
        key=f"{key_prefix}download_structure",
    )


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
