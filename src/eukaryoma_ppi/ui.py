"""Shared Streamlit rendering used by more than one page."""

import altair as alt
import pandas as pd
import streamlit as st
import streamlit.components.v1 as components

from eukaryoma_ppi import complex_annotations, external_scores, external_structures, structures, viewer, yeast_orthologs
from eukaryoma_ppi.config import YEAST_ORTHOLOGS_INDEX_FILE

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

# Cross-page navigation: a page that already knows which pair the user wants
# (e.g. a clicked scatter point) can jump straight to Pair Viewer with that
# pair pre-selected, via a plain session_state handoff -- Pair Viewer's own
# selectbox keys get cleared on navigation away (the same widget-state gotcha
# unified_score_mode works around above), so this key is deliberately not a
# widget's own key, and Pair Viewer reads+clears it exactly once at the top
# of its script, before creating its selectboxes.
PENDING_PAIR_KEY = "pending_pair_selection"


def request_pair_view(protein_a, protein_b):
    """Ask Pair Viewer to preselect (protein_a, protein_b) and navigate
    there. Call this as the direct result of a widget interaction (e.g.
    inside an `if st.button(...):` block), not unconditionally every rerun.
    """
    st.session_state[PENDING_PAIR_KEY] = (protein_a, protein_b)
    st.switch_page("pages/2_Pair_Viewer.py")


def consume_pending_pair():
    """Pop and return a pending (protein_a, protein_b) pair set by
    request_pair_view, or None. Pair Viewer calls this once, before creating
    its own protein selectboxes, to seed them.
    """
    return st.session_state.pop(PENDING_PAIR_KEY, None)

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


@st.cache_data
def get_single_chain_pdb_text(pool, chain_id):
    return structures.extract_chains_as_pdb(pool, [chain_id])


@st.cache_data
def get_superposed_pdb_texts(pool_chain_pairs):
    return structures.superpose_chains_as_pdb(pool_chain_pairs)


MAX_STRUCTURE_COMPARISON = 8


def render_structure_comparison(protein_id, protein_pools, key_prefix=""):
    """One protein, several pools -- compare its predicted structure across
    them, either superposed in one view or side by side. protein_pools: a
    DataFrame with "pool"/"chain" columns, one row per pool this protein has
    an AF3 structure in (already deduplicated by caller).
    """
    if len(protein_pools) < 2:
        st.info("This protein has an AF3 structure in fewer than 2 pools, so there's nothing to compare.")
        return

    pool_options = protein_pools["pool"].tolist()
    default = pool_options[: min(4, len(pool_options))]
    selected_pools = st.multiselect(
        f"Pools to compare (up to {MAX_STRUCTURE_COMPARISON})",
        pool_options,
        default=default,
        key=f"{key_prefix}compare_pools",
    )
    if len(selected_pools) > MAX_STRUCTURE_COMPARISON:
        st.warning(f"Comparing this many structures at once gets slow to render -- showing the first {MAX_STRUCTURE_COMPARISON} selected.")
        selected_pools = selected_pools[:MAX_STRUCTURE_COMPARISON]
    if len(selected_pools) < 2:
        st.info("Select at least 2 pools to compare.")
        return

    chain_by_pool = dict(zip(protein_pools["pool"], protein_pools["chain"]))
    pool_chain_pairs = tuple((pool, chain_by_pool[pool]) for pool in selected_pools)

    mode = st.radio(
        "Comparison mode",
        ["Superimpose (aligned overlay)", "Side by side"],
        horizontal=True,
        key=f"{key_prefix}compare_mode",
    )

    if mode == "Superimpose (aligned overlay)":
        pdb_texts = get_superposed_pdb_texts(pool_chain_pairs)
        ok = [(pool, pdb) for (pool, _chain), pdb in zip(pool_chain_pairs, pdb_texts) if pdb is not None]
        skipped = [pool for (pool, _chain), pdb in zip(pool_chain_pairs, pdb_texts) if pdb is None]
        if skipped:
            st.warning(
                f"Skipped (different residue count than the first selected pool, so they can't be aligned "
                f"1:1): {', '.join(skipped)}."
            )
        if len(ok) < 2:
            st.info("Not enough compatible structures left to overlay.")
            return

        legend = " &nbsp; ".join(
            f'<span style="color:{viewer.COMPARISON_COLORS[i % len(viewer.COMPARISON_COLORS)]}">●</span> {pool}'
            for i, (pool, _pdb) in enumerate(ok)
        )
        st.markdown(legend, unsafe_allow_html=True)
        html = viewer.render_overlay([pdb for _pool, pdb in ok])
        components.html(html, height=580)
    else:
        cols = st.columns(min(len(pool_chain_pairs), 4))
        for i, (pool, chain) in enumerate(pool_chain_pairs):
            pdb_text = get_single_chain_pdb_text(pool, chain)
            color = viewer.COMPARISON_COLORS[i % len(viewer.COMPARISON_COLORS)]
            html = viewer.render_single(pdb_text, color)
            with cols[i % len(cols)]:
                st.caption(f"Pool **{pool}**")
                components.html(html, height=340)


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
    """3D viewer if this pair has an AF3 pool structure, else a plain notice.

    Also offers, right alongside it: an optional reference structure
    (human/yeast ortholog, aligned onto the matching chain and shown in the
    same view) and optional interface-contact highlighting.
    """
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

    reference_selection = render_reference_structure_selector(protein_a, protein_b, key_prefix=key_prefix)
    reference_pdb_text = None
    if reference_selection is not None:
        kind, accession, ref_protein_id, raw_reference_pdb = reference_selection
        target_chain_id = chain_for[ref_protein_id]
        aligned = get_aligned_reference(pdb_text, target_chain_id, raw_reference_pdb)
        if aligned is None:
            st.warning("Could not align this reference onto the pair (no usable sequence overlap).")
        else:
            reference_pdb_text, n_core, rmsd = aligned
            st.caption(
                f"Reference ({kind} {accession}, shown in green) aligned onto chain {target_chain_id} "
                f"({ref_protein_id}) -- {n_core}-residue conserved core, RMSD {rmsd:.2f} Å."
            )

    show_contacts = st.checkbox("Highlight contact residues (interface)", key=f"{key_prefix}show_contacts")
    contact_resi_a = contact_resi_b = None
    contacts_df = pd.DataFrame()
    if show_contacts:
        cutoff = st.slider(
            "Contact cutoff (Cα-Cα distance, Å)", 4.0, 15.0, 8.0, step=0.5, key=f"{key_prefix}contact_cutoff"
        )
        contacts_df = get_contacts(pdb_text, chain_a, chain_b, cutoff)
        if contacts_df.empty:
            st.caption("No contacts found at this cutoff.")
        else:
            contact_resi_a = sorted(contacts_df["resi_a"].unique().tolist())
            contact_resi_b = sorted(contacts_df["resi_b"].unique().tolist())
            st.caption(
                f"{len(contacts_df)} contacting residue pairs ({len(contact_resi_a)} residues in {protein_a}, "
                f"{len(contact_resi_b)} in {protein_b}) highlighted in gold."
            )

    html = viewer.render_pair(
        pdb_text,
        chain_a,
        chain_b,
        color_by=color_by,
        reference_pdb_text=reference_pdb_text,
        contact_resi_a=contact_resi_a,
        contact_resi_b=contact_resi_b,
    )

    st.caption(f"Pool **{pool}** &mdash; chain {chain_a} = {protein_a}, chain {chain_b} = {protein_b}")
    components.html(html, height=580)
    st.download_button(
        "Download structure (.pdb)",
        data=pdb_text,
        file_name=f"{pool}_{protein_a}_{protein_b}.pdb",
        mime="chemical/x-pdb",
        key=f"{key_prefix}download_structure",
    )

    if not contacts_df.empty:
        st.altair_chart(render_contact_position_chart(contacts_df, protein_a, protein_b), width="stretch")


@st.cache_data
def get_human_uniprot_orthologs():
    """{protein_id: {uniprot_accession, ...}} via CORUM/Marcotte's own
    "uniprot" column -- see complex_annotations.human_uniprot_orthologs.
    """
    return complex_annotations.human_uniprot_orthologs()


@st.cache_data
def get_yeast_orthologs():
    """{protein_id: {yeast_accession, ...}} -- see
    yeast_orthologs.yeast_orthologs_map. Empty if scripts/build_data.py
    hasn't resolved/decompressed any yeast reference structures.
    """
    if not YEAST_ORTHOLOGS_INDEX_FILE.exists():
        return {}
    return yeast_orthologs.yeast_orthologs_map()


@st.cache_data(show_spinner="Fetching reference structure from AlphaFold DB...")
def get_alphafold_reference_pdb(accession):
    return external_structures.fetch_alphafold_pdb(accession)


@st.cache_data
def get_yeast_reference_pdb(accession):
    return structures.get_yeast_structure_pdb(accession)


@st.cache_data
def get_aligned_reference(pair_pdb_text, target_chain_id, reference_pdb_text):
    return structures.align_reference_onto_chain(pair_pdb_text, target_chain_id, reference_pdb_text)


@st.cache_data
def get_contacts(pdb_text, chain_a_id, chain_b_id, cutoff):
    return structures.find_contact_residues(pdb_text, chain_a_id, chain_b_id, cutoff)


def render_reference_structure_selector(protein_a, protein_b, key_prefix=""):
    """Selectbox for an optional reference structure: either protein's
    ortholog's own solo prediction, human (fetched live from AlphaFold DB,
    by the UniProt accession CORUM/Marcotte cross-reference it to) or yeast
    (a locally pre-decompressed AF3 monomer prediction, by the UniProt
    accession a Capsaspora<->yeast HOG cross-reference resolves to -- see
    eukaryoma_ppi.yeast_orthologs). Returns None if neither protein has a
    usable accession of either kind, or the user picked "None"; otherwise
    (kind, accession, protein_id, reference_pdb_text) for the caller to
    align onto protein_id's chain in the pair and merge into one view.
    """
    human_orthologs = get_human_uniprot_orthologs()
    yeast_orthologs_by_protein = get_yeast_orthologs()
    options = {}
    for label, protein_id in [("Protein A", protein_a), ("Protein B", protein_b)]:
        for accession in sorted(human_orthologs.get(protein_id, ())):
            options[f"{label} ({protein_id}) — human {accession}"] = ("human", accession, protein_id)
        for accession in sorted(yeast_orthologs_by_protein.get(protein_id, ())):
            options[f"{label} ({protein_id}) — yeast {accession}"] = ("yeast", accession, protein_id)
    if not options:
        return None

    st.caption(
        "Optional reference: an ortholog's own solo prediction, either human (via CORUM/Marcotte's complex "
        "membership, fetched live from **AlphaFold DB**) or yeast (via a Capsaspora-yeast orthologous group, "
        "from a locally predicted AF3 monomer) -- shown aligned onto the matching chain above, in green."
    )
    choice = st.selectbox("Reference structure", ["None"] + list(options.keys()), key=f"{key_prefix}reference_structure")
    if choice == "None":
        return None

    kind, accession, ref_protein_id = options[choice]
    if kind == "human":
        reference_pdb = get_alphafold_reference_pdb(accession)
        if reference_pdb is None:
            st.warning(f"Could not fetch the AlphaFold DB structure for {accession} (not available, or a network error).")
            return None
    else:
        reference_pdb = get_yeast_reference_pdb(accession)
        if reference_pdb is None:
            st.warning(f"No local structure file for yeast ortholog {accession}.")
            return None

    return kind, accession, ref_protein_id, reference_pdb


def render_contact_position_chart(contacts_df, protein_a_label, protein_b_label):
    """A simple sequence/position view of the interface: one row per
    protein, a tick per residue position that has a contact.
    """
    a_positions = (
        contacts_df[["resi_a", "resn_a"]]
        .drop_duplicates()
        .rename(columns={"resi_a": "position", "resn_a": "resn"})
    )
    a_positions["protein"] = protein_a_label
    b_positions = (
        contacts_df[["resi_b", "resn_b"]]
        .drop_duplicates()
        .rename(columns={"resi_b": "position", "resn_b": "resn"})
    )
    b_positions["protein"] = protein_b_label
    combined = pd.concat([a_positions, b_positions], ignore_index=True)

    return (
        alt.Chart(combined)
        .mark_tick(thickness=3, size=20)
        .encode(
            x=alt.X("position:Q", title="Residue position"),
            y=alt.Y("protein:N", title=None, axis=alt.Axis(labelLimit=200), scale=alt.Scale(paddingOuter=0.6)),
            color=alt.Color("protein:N", legend=None, scale=alt.Scale(range=viewer.CHAIN_COLORS)),
            tooltip=["protein", "position", "resn"],
        )
        .properties(height=140, title="Contact residue positions")
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
