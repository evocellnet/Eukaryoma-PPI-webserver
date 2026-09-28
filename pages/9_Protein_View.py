import numpy as np
import pandas as pd
import altair as alt
import streamlit as st

from eukaryoma_ppi import annotations, eggnog, external_scores, go_enrichment, ui
from eukaryoma_ppi.config import ANNOTATIONS_INDEX_FILE, EGGNOG_INDEX_FILE, UNIVERSE_INDEX_FILE

st.set_page_config(page_title="Protein View - Eukaryoma PPI", page_icon="🧬", layout="wide")
st.title("Protein View")
st.write(
    "Everything the site knows about a single protein: every interaction it's part of (sortable by any "
    "score), links to external databases, and its eggNOG-mapper functional annotation."
)

if not UNIVERSE_INDEX_FILE.exists():
    st.error("No universe index found. Run `python scripts/build_data.py` first (see README.md).")
    st.stop()

SCORE_COLUMNS = list(external_scores.ALL_SOURCE_LABELS.items()) + [("unified_score", "Unified score")]


@st.cache_data
def get_universe_for_protein_view():
    df = external_scores.load_universe()
    if ANNOTATIONS_INDEX_FILE.exists():
        df = annotations.attach_annotations(df)
    else:
        df["annotation_a"] = ""
        df["annotation_b"] = ""
    return df


@st.cache_data
def get_scored_universe(_universe_df, mode, weights):
    return ui.recompute_unified_score(_universe_df, external_scores.UNIFIED_SCORE_COLUMNS, mode, weights=weights)[0]


@st.cache_data
def get_annotation_map():
    if not ANNOTATIONS_INDEX_FILE.exists():
        return {}
    df = annotations.load_annotations()
    return dict(zip(df["protein_id"], df["annotation"]))


@st.cache_data
def get_eggnog_map():
    if not EGGNOG_INDEX_FILE.exists():
        return {}
    return eggnog.load_eggnog_annotations().set_index("protein_id").to_dict("index")


_mode = ui.get_unified_score_mode()
_weights = None
if _mode == "weighted":
    _default_weights = external_scores.compute_default_weights(
        get_universe_for_protein_view(), external_scores.UNIFIED_SCORE_COLUMNS
    )
    _weights = ui.get_active_unified_score_weights(_default_weights)
universe_df = get_scored_universe(get_universe_for_protein_view(), _mode, _weights)
st.caption(f"Unified score: **{ui.UNIFIED_SCORE_MODES[_mode]}** (change this on the Unified Ranking page).")
annotation_map = get_annotation_map()
eggnog_map = get_eggnog_map()
protein_options = sorted(set(universe_df["protein_a"]) | set(universe_df["protein_b"]))


def format_protein_option(pid):
    annotation = annotation_map.get(pid)
    return f"{pid} — {annotation}" if annotation else pid


protein_id = st.selectbox(
    "Select a protein",
    protein_options,
    index=None,
    placeholder="Search by id or annotation",
    format_func=format_protein_option,
)

if not protein_id:
    st.info("Pick a protein to see everything the site knows about it.")
    st.stop()

st.header(protein_id)
if annotation_map.get(protein_id):
    st.caption(f"OMA/FASTA description: {annotation_map[protein_id]}")

st.divider()
st.subheader("External links & eggNOG-mapper annotation")
st.caption(
    "The functional annotation below comes from **eggNOG-mapper** (ortholog matching against reference "
    "databases) -- a different source from the OMA/FASTA description above, shown separately so the two "
    "aren't confused."
)

eggnog_row = eggnog_map.get(protein_id)
if eggnog_row is None:
    st.info("No eggNOG-mapper annotation available for this protein.")
else:
    accession, kind = eggnog_row["ortholog_accession"], eggnog_row["ortholog_kind"]
    if accession:
        link_cols = st.columns(3)
        link_cols[0].link_button(
            "NCBI Protein" if kind == "refseq" else "NCBI Protein (search)", eggnog.ncbi_protein_url(accession)
        )
        uniprot_href, uniprot_direct = eggnog.uniprot_url(accession, kind)
        link_cols[1].link_button("UniProt" if uniprot_direct else "UniProt (search)", uniprot_href)
        af_href, af_direct = eggnog.alphafold_url(accession, kind)
        link_cols[2].link_button("AlphaFold DB" if af_direct else "AlphaFold DB (search)", af_href)
        if kind == "refseq":
            st.caption(
                "eggNOG's ortholog match is this protein's own NCBI RefSeq entry, so UniProt/AlphaFold DB "
                "links are text searches rather than guaranteed direct hits."
            )
    else:
        st.caption("eggNOG has no ortholog identifier for this protein to link out with.")

    with st.expander("eggNOG-mapper details", expanded=True):
        if eggnog_row.get("Description"):
            st.write(f"**Description:** {eggnog_row['Description']}")
        if eggnog_row.get("Preferred_name"):
            st.write(f"**Preferred gene name:** {eggnog_row['Preferred_name']}")
        if eggnog_row.get("COG_category") not in (None, "-"):
            st.write(f"**COG category:** {eggnog_row['COG_category']}")
        if eggnog_row.get("eggNOG_OGs") not in (None, "-"):
            st.write(f"**eggNOG orthologous groups:** {eggnog_row['eggNOG_OGs']}")

        gos = eggnog.as_list(eggnog_row.get("GOs"))
        if gos:
            shown = gos[:40]
            links = " · ".join(f"[{g}]({eggnog.quickgo_url(g)})" for g in shown)
            st.markdown(f"**GO terms ({len(gos)}):** {links}" + (f" *(+{len(gos) - 40} more)*" if len(gos) > 40 else ""))

        kegg_ko = eggnog.as_list(eggnog_row.get("KEGG_ko"))
        if kegg_ko:
            links = " · ".join(f"[{k}]({eggnog.kegg_url(k)})" for k in kegg_ko)
            st.markdown(f"**KEGG orthologs:** {links}")

        kegg_pathway = eggnog.as_list(eggnog_row.get("KEGG_Pathway"))
        if kegg_pathway:
            links = " · ".join(f"[{k}]({eggnog.kegg_url(k)})" for k in kegg_pathway)
            st.markdown(f"**KEGG pathways:** {links}")

        pfams = eggnog.as_list(eggnog_row.get("PFAMs"))
        if pfams:
            st.write("**PFAM domains:** " + ", ".join(pfams))

st.divider()
st.header("Interactions involving this protein")

is_a = universe_df["protein_a"] == protein_id
protein_pairs = universe_df[is_a | (universe_df["protein_b"] == protein_id)].copy()
protein_pairs["partner"] = np.where(protein_pairs["protein_a"] == protein_id, protein_pairs["protein_b"], protein_pairs["protein_a"])
protein_pairs["partner_annotation"] = protein_pairs["partner"].map(annotation_map).fillna("")
st.caption(f"{len(protein_pairs):,} pairs involve this protein across any data source.")

sort_options = {label: col for col, label in SCORE_COLUMNS}
sort_label = st.selectbox("Sort interactions by", list(sort_options.keys()), index=len(sort_options) - 1)
sort_col = sort_options[sort_label]

ranked = protein_pairs.sort_values(sort_col, ascending=False, na_position="last").reset_index(drop=True)
top_n = st.number_input(
    "Show top N interactions",
    min_value=min(10, len(ranked)),
    max_value=len(ranked),
    value=min(500, len(ranked)),
    step=10,
)
displayed = ranked.head(top_n)

if any(col in displayed.columns for col in ui.TP_COLUMNS):
    st.caption(ui.TP_LEGEND)
st.caption("Select rows to add those partners to the interactome below (in addition to the auto-included ones).")

column_order = [sort_col, "partner", "partner_annotation"] + [c for c, _ in SCORE_COLUMNS if c != sort_col]
column_config = {
    sort_col: st.column_config.NumberColumn(sort_label, format="%.3f"),
    "partner": st.column_config.TextColumn("Partner"),
    "partner_annotation": st.column_config.TextColumn("Partner annotation"),
}
for col, label in SCORE_COLUMNS:
    if col != sort_col:
        column_config[col] = st.column_config.NumberColumn(label, format="%.3f")

table_event = st.dataframe(
    ui.style_true_positive_rows(displayed),
    column_order=column_order,
    column_config=column_config,
    hide_index=True,
    width="stretch",
    on_select="rerun",
    selection_mode="multi-row",
    key="protein_interactions_table",
)
selected_rows = table_event.selection.rows
manually_added = set(displayed.iloc[selected_rows]["partner"]) if selected_rows else set()

if len(selected_rows) == 1:
    row = displayed.iloc[selected_rows[0]]
    st.divider()
    st.subheader(f"{protein_id} — {row['partner']}")
    ui.render_scores(row, SCORE_COLUMNS)
    ui.render_true_positive_badges(row)
    ui.render_structure_if_available(row, protein_id, row["partner"], key_prefix="protview_")

st.divider()
st.header("Interactome")
st.write(
    "Starts small (top partners by unified score above the threshold); select rows in the table above to "
    "pull in additional proteins regardless of score."
)

unified_threshold = st.slider(
    "Auto-include partners with unified score above", 0.0, 1.0, 0.9, step=0.01, key="network_threshold"
)
auto_included = set(protein_pairs.loc[protein_pairs["unified_score"] >= unified_threshold, "partner"])
node_ids = [protein_id] + sorted(auto_included | manually_added)

if len(node_ids) < 2:
    st.info("No partners meet the threshold yet -- lower it, or select rows in the table above.")
else:
    n = len(node_ids)
    angles = np.linspace(0, 2 * np.pi, n - 1, endpoint=False)
    positions = {protein_id: (0.0, 0.0)}
    for angle, pid in zip(angles, node_ids[1:]):
        positions[pid] = (np.cos(angle), np.sin(angle))

    nodes_df = pd.DataFrame(
        [
            {
                "protein_id": pid,
                "x": positions[pid][0],
                "y": positions[pid][1],
                "annotation": annotation_map.get(pid, ""),
                "role": "center" if pid == protein_id else "partner",
            }
            for pid in node_ids
        ]
    )

    # universe_df's protein_a/protein_b are already canonically (lo, hi)
    # ordered, so isin() on both columns finds every edge among node_ids
    # without needing to generate/query individual (a, b) pairs by hand.
    node_set = set(node_ids)
    candidate_edges = universe_df[universe_df["protein_a"].isin(node_set) & universe_df["protein_b"].isin(node_set)]
    has_tp_columns = all(col in candidate_edges.columns for col in ui.TP_COLUMNS)
    edge_cols = ["protein_a", "protein_b", "unified_score"] + (ui.TP_COLUMNS if has_tp_columns else [])
    edges_df = candidate_edges[edge_cols].dropna(subset=["unified_score"]).copy()
    edges_df = edges_df.rename(columns={"protein_a": "a", "protein_b": "b"})
    edges_df["x"] = edges_df["a"].map(lambda p: positions[p][0])
    edges_df["y"] = edges_df["a"].map(lambda p: positions[p][1])
    edges_df["x2"] = edges_df["b"].map(lambda p: positions[p][0])
    edges_df["y2"] = edges_df["b"].map(lambda p: positions[p][1])

    if has_tp_columns:
        def _tp_label(row):
            corum, marcotte = bool(row["corum_tp"]), bool(row["marcotte_tp"])
            if corum and marcotte:
                return "Both"
            if corum:
                return "CORUM"
            if marcotte:
                return "Marcotte"
            return "None"

        edges_df["true_positive"] = edges_df.apply(_tp_label, axis=1)

    layers = []
    if not edges_df.empty:
        edge_encoding = dict(
            x=alt.X("x:Q", axis=None, scale=alt.Scale(domain=[-1.3, 1.3])),
            y=alt.Y("y:Q", axis=None, scale=alt.Scale(domain=[-1.3, 1.3])),
            x2=alt.X2("x2:Q"),
            y2=alt.Y2("y2:Q"),
            strokeWidth=alt.StrokeWidth("unified_score:Q", legend=None, scale=alt.Scale(range=[1, 6])),
        )
        if has_tp_columns:
            edge_encoding["color"] = alt.Color(
                "true_positive:N",
                title="Known true positive",
                scale=alt.Scale(
                    domain=["None", "CORUM", "Marcotte", "Both"],
                    range=["#999999", "#2ca02c", "#3d6fd6", "#9467bd"],
                ),
                legend=None,
            )
            edge_encoding["tooltip"] = ["a", "b", alt.Tooltip("unified_score:Q", format=".3f"), "true_positive"]
        else:
            edge_encoding["color"] = alt.value("#666666")
            edge_encoding["tooltip"] = ["a", "b", alt.Tooltip("unified_score:Q", format=".3f")]
        layers.append(alt.Chart(edges_df).mark_rule().encode(**edge_encoding))
    layers.append(
        alt.Chart(nodes_df)
        .mark_circle(size=500, stroke="#1a1a1a", strokeWidth=1)
        .encode(
            x=alt.X("x:Q", axis=None, scale=alt.Scale(domain=[-1.3, 1.3])),
            y=alt.Y("y:Q", axis=None, scale=alt.Scale(domain=[-1.3, 1.3])),
            color=alt.Color("role:N", legend=None, scale=alt.Scale(domain=["center", "partner"], range=["#F58518", "#4C78A8"])),
            tooltip=["protein_id", "annotation"],
        )
    )
    layers.append(
        alt.Chart(nodes_df)
        .mark_text(dy=-16, fontSize=11, fontWeight="bold", color="#1a1a1a")
        .encode(x="x:Q", y="y:Q", text="protein_id:N")
    )

    # Fixed light background regardless of the app's light/dark theme --
    # these marks use literal colors that don't adapt to Streamlit's theme,
    # so a dark app background left the (default black) text unreadable.
    network_chart = (
        alt.layer(*layers)
        .properties(width=600, height=600, background="#FFFFFF")
        .configure_view(strokeWidth=0)
    )
    st.altair_chart(network_chart, theme=None)
    caption = f"{len(node_ids)} proteins, {len(edges_df)} known interactions among them shown as edges."
    if has_tp_columns:
        caption += " Edge color marks known true-positive interactions: green CORUM, blue Marcotte, purple both."
    st.caption(caption)

    st.subheader("GO term enrichment of this interactome")
    st.write(
        "Which GO terms are over-represented among these proteins, compared to every eggNOG-annotated "
        "website protein as background (hypergeometric test, Benjamini-Hochberg FDR-corrected). GO term "
        "*names* aren't shown -- this site doesn't bundle a GO ontology file -- click a term to look it up "
        "on QuickGO."
    )

    @st.cache_data
    def get_go_enrichment(foreground_ids, _background_ids):
        return go_enrichment.enrichment(set(foreground_ids), eggnog.load_eggnog_annotations(), _background_ids)

    website_ids_set = set(protein_options)
    enrichment_df = get_go_enrichment(tuple(sorted(node_ids)), website_ids_set)
    fdr_threshold = st.slider("FDR threshold", 0.0, 1.0, 0.05, step=0.01, key="go_fdr_threshold")
    significant = enrichment_df[enrichment_df["fdr"] <= fdr_threshold] if not enrichment_df.empty else enrichment_df

    if significant.empty:
        st.info(f"No GO terms pass FDR ≤ {fdr_threshold:.2f} for this set of {len(node_ids)} proteins.")
    else:
        st.caption(
            f"{len(significant):,} GO term(s) enriched at FDR ≤ {fdr_threshold:.2f} "
            f"({len(node_ids)} proteins vs. {len(website_ids_set):,} website background)."
        )
        display_df = significant.head(50).copy()
        display_df["go_url"] = display_df["go_id"].map(eggnog.quickgo_url)
        st.dataframe(
            display_df[["go_url", "n_foreground", "n_background", "pvalue", "fdr"]],
            hide_index=True,
            width="stretch",
            column_config={
                "go_url": st.column_config.LinkColumn("GO term", display_text=r"QuickGO/term/(.*)"),
                "n_foreground": st.column_config.NumberColumn("In this set"),
                "n_background": st.column_config.NumberColumn("In website background"),
                "pvalue": st.column_config.NumberColumn("p-value", format="%.2e"),
                "fdr": st.column_config.NumberColumn("FDR", format="%.2e"),
            },
        )

    with st.expander("Just list the GO terms present (no statistics)"):
        freq_df = go_enrichment.term_frequency(node_ids, eggnog.load_eggnog_annotations())
        if freq_df.empty:
            st.info("None of these proteins have an eggNOG GO annotation.")
        else:
            freq_df = freq_df.copy()
            freq_df["go_url"] = freq_df["go_id"].map(eggnog.quickgo_url)
            st.dataframe(
                freq_df[["go_url", "n_proteins", "fraction"]],
                hide_index=True,
                width="stretch",
                column_config={
                    "go_url": st.column_config.LinkColumn("GO term", display_text=r"QuickGO/term/(.*)"),
                    "n_proteins": st.column_config.NumberColumn("Proteins with this term"),
                    "fraction": st.column_config.NumberColumn("Fraction of this set", format="percent"),
                },
            )

st.divider()
st.header("Compare structures across pools")
st.write(
    "This protein may have been predicted in several different pools, each time paired with a different "
    "partner. Compare those predictions of *this protein alone* to see how consistent its predicted fold is "
    "across pooling contexts."
)

protein_chain = np.where(protein_pairs["protein_a"] == protein_id, protein_pairs["chain_a"], protein_pairs["chain_b"])
protein_pools = (
    protein_pairs.assign(chain=protein_chain)
    .dropna(subset=["pool", "chain"])
    .drop_duplicates(subset=["pool"])[["pool", "chain"]]
    .sort_values("pool")
    .reset_index(drop=True)
)
ui.render_structure_comparison(protein_id, protein_pools, key_prefix="protview_compare_")
