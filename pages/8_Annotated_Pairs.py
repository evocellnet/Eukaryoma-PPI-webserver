import altair as alt
import numpy as np
import pandas as pd
import streamlit as st

from eukaryoma_ppi import annotations, external_scores, tp_analysis, ui
from eukaryoma_ppi.config import ANNOTATIONS_INDEX_FILE, TRUE_POSITIVE_INDEX_FILE, UNIVERSE_INDEX_FILE

st.set_page_config(page_title="Annotated Pairs - Eukaryoma PPI", page_icon="🧬", layout="wide")
st.title("Annotated Pairs")
st.write(
    "How well do the association scores agree with CORUM/Marcotte complex co-membership -- treated here as "
    "known true-positive interactions? Use this page to check score quality and to spot pairs where the "
    "score and the annotation disagree. For what CORUM/Marcotte actually are and how many complexes/pairs "
    "they contribute, see the **Annotation Sources** page."
)

if not UNIVERSE_INDEX_FILE.exists() or not TRUE_POSITIVE_INDEX_FILE.exists():
    st.error("No universe/true-positive index found. Run `python scripts/build_data.py` first (see README.md).")
    st.stop()


@st.cache_data
def get_universe_for_annotations():
    universe_df = external_scores.load_universe()
    if ANNOTATIONS_INDEX_FILE.exists():
        universe_df = annotations.attach_annotations(universe_df)
    else:
        universe_df["annotation_a"] = ""
        universe_df["annotation_b"] = ""
    return universe_df


@st.cache_data
def get_scored_universe(_universe_df, mode, weights):
    return ui.recompute_unified_score(_universe_df, external_scores.UNIFIED_SCORE_COLUMNS, mode, weights=weights)[0]


_mode = ui.get_unified_score_mode()
_weights = None
if _mode == "weighted":
    _default_weights = external_scores.compute_default_weights(
        get_universe_for_annotations(), external_scores.UNIFIED_SCORE_COLUMNS
    )
    _weights = ui.get_active_unified_score_weights(_default_weights)
universe_df = get_scored_universe(get_universe_for_annotations(), _mode, _weights)
st.caption(f"Unified score: **{ui.UNIFIED_SCORE_MODES[_mode]}** (change this on the Unified Ranking page).")

SCORE_COLUMNS = list(external_scores.ALL_SOURCE_LABELS.items()) + [("unified_score", "Unified score")]
SCORE_COLUMNS = [(col, label) for col, label in SCORE_COLUMNS if universe_df[col].notna().any()]

ground_truth_mode = st.selectbox("Ground truth", list(tp_analysis.GROUND_TRUTH_MODES.keys()), index=2)
ground_truth = tp_analysis.ground_truth_series(universe_df, ground_truth_mode)
st.caption(f"{int(ground_truth.sum()):,} of {len(universe_df):,} pairs are true-positive under this definition.")

st.divider()
st.header("True-positive rate by score quantile")
st.write(
    "Pairs are split into equal-sized bins by each score (bin 1 = lowest scores, bin N = highest), and this "
    "plots what fraction of pairs in each bin are true-positive. A score that tracks real interactions well "
    "should show a rising line."
)
n_bins = st.slider("Number of quantile bins", min_value=4, max_value=20, value=10)

quantile_df = tp_analysis.quantile_tp_rates_all_scores(universe_df, SCORE_COLUMNS, ground_truth, n_bins)
if quantile_df.empty:
    st.info("No scores available to compute quantile true-positive rates.")
else:
    chart = (
        alt.Chart(quantile_df)
        .mark_line(point=True)
        .encode(
            x=alt.X("quantile:O", title=f"Score quantile (1=lowest, {n_bins}=highest)"),
            y=alt.Y("tp_rate:Q", title="True-positive rate", axis=alt.Axis(format="%")),
            color=alt.Color("score:N", title="Score"),
            tooltip=["score", "quantile", alt.Tooltip("tp_rate:Q", format=".1%"), "n_pairs"],
        )
        .properties(height=350)
    )
    st.altair_chart(chart, width="stretch")

st.divider()
st.header("Explore score vs. annotation")
st.write(
    "Each dot is an *annotated* pair (CORUM/Marcotte/both), jittered by category, plotted against the grey "
    "violin -- the not-annotated population's full score density, not a sample, so its shape is always "
    "accurate even where no dots are drawn. Not-annotated pairs aren't plotted as individual dots: there can "
    "be millions of them, and a small random sample of dots would misrepresent where they actually sit -- "
    "drag a box over the violin (or the dot area) to pull up the actual not-annotated pairs in that score "
    "range below. Click a single dot to select just that pair."
)

explore_score_col, explore_score_label = st.selectbox(
    "Score to explore", SCORE_COLUMNS, format_func=lambda pair: pair[1], key="explore_score"
)


@st.cache_data
def get_baseline_density(_universe_df, score_col, n_bins=40):
    return tp_analysis.baseline_score_density(_universe_df, score_col, n_bins)


density_df = get_baseline_density(universe_df, explore_score_col)
# Renamed to match the scatter's x-field: the brush selection is shared
# between the violin and the scatter (see combined_chart below) by field
# name, so a box dragged on the violin (whose own x-field would otherwise
# be "score_bin") reports under the same key the code below reads.
if not density_df.empty:
    density_df = density_df.rename(columns={"score_bin": explore_score_col})

MAX_PLOTTED = 4000  # annotated pairs only now -- CORUM/Marcotte/both rarely exceed this for any one score


@st.cache_data
def get_plot_sample(_universe_df, score_col, max_plotted, seed=0):
    valid = _universe_df[score_col].notna()
    plot_df = _universe_df.loc[valid, ["protein_a", "protein_b", score_col]].copy()
    plot_df["tp_category"] = tp_analysis.tp_category(_universe_df.loc[valid])
    plot_df = plot_df[plot_df["tp_category"] != "Not annotated"]

    if len(plot_df) > max_plotted:
        plot_df = plot_df.sample(max_plotted, random_state=seed)

    rng = np.random.default_rng(seed)
    plot_df["jitter"] = rng.uniform(-0.4, 0.4, size=len(plot_df))
    return plot_df.reset_index(drop=True)


plot_df = get_plot_sample(universe_df, explore_score_col, MAX_PLOTTED)
n_annotated = int((tp_analysis.tp_category(universe_df) != "Not annotated").sum())
st.caption(f"Showing {len(plot_df):,} of {n_annotated:,} annotated pairs with a {explore_score_label} score.")

brush = alt.selection_interval(name="brush", encodings=["x", "y"])
point_select = alt.selection_point(name="point_select", fields=["protein_a", "protein_b"], on="click", empty=False)
category_order = ["CORUM", "Marcotte", "Both"]
scatter = (
    alt.Chart(plot_df)
    .mark_circle(size=40, opacity=0.5)
    .encode(
        x=alt.X(f"{explore_score_col}:Q", title=explore_score_label),
        y=alt.Y("tp_category:N", title="Annotation", sort=category_order),
        yOffset="jitter:Q",
        color=alt.Color("tp_category:N", title="Annotation", sort=category_order),
        tooltip=["protein_a", "protein_b", alt.Tooltip(f"{explore_score_col}:Q", format=".3f"), "tp_category"],
    )
    .add_params(brush, point_select)
    .properties(height=320)
)

if density_df.empty:
    combined_chart = scatter
else:
    max_density = float(density_df["density"].max()) * 1.15
    density_scale = alt.Scale(domain=[-max_density, max_density])
    violin_top = (
        alt.Chart(density_df)
        .mark_area(color="#888888", opacity=0.7, interpolate="monotone")
        .encode(x=alt.X(f"{explore_score_col}:Q", title=None), y=alt.Y("density:Q", axis=None, scale=density_scale))
    )
    violin_bottom = (
        alt.Chart(density_df)
        .mark_area(color="#888888", opacity=0.7, interpolate="monotone")
        .encode(x=alt.X(f"{explore_score_col}:Q", title=None), y=alt.Y("neg_density:Q", axis=None, scale=density_scale))
    )
    violin = (
        (violin_top + violin_bottom)
        .add_params(brush)
        .properties(height=70, title="Not-annotated baseline (score density, full population)")
    )
    combined_chart = alt.vconcat(violin, scatter).resolve_scale(x="shared")

event = st.altair_chart(
    combined_chart, on_select="rerun", selection_mode=["brush", "point_select"], key="annotation_scatter"
)

# A point selection reports as a list of dicts (one per selected mark,
# with the "fields" values for that point) -- unlike an interval/brush
# selection, which reports as {field: [range-or-categories]}.
point_sel = event.selection.get("point_select", []) if event.selection else []
if point_sel:
    sel_a, sel_b = point_sel[0]["protein_a"], point_sel[0]["protein_b"]
    st.success(f"Selected pair: **{sel_a}** — **{sel_b}**")
    st.code(f"{sel_a}\t{sel_b}", language=None)
    if st.button("Open in Pair Viewer", key="open_pair_viewer_from_scatter"):
        ui.request_pair_view(sel_a, sel_b)

brush_sel = event.selection.get("brush", {}) if event.selection else {}
x_range = brush_sel.get(explore_score_col)
categories = brush_sel.get("tp_category")

if not x_range:
    st.info("Drag a box on the plot above to list the pairs in that region.")
else:
    mask = universe_df[explore_score_col].between(min(x_range), max(x_range)).to_numpy()
    if categories:
        mask = mask & np.isin(tp_analysis.tp_category(universe_df), categories)
    selected_df = universe_df[mask]
    st.write(f"**{len(selected_df):,} pairs** in the selected region.")
    ui.render_ranked_pair_table(
        selected_df,
        sort_col=explore_score_col,
        column_order=[explore_score_col, "protein_a", "annotation_a", "protein_b", "annotation_b"],
        column_config={
            explore_score_col: st.column_config.NumberColumn(explore_score_label, format="%.3f"),
            "protein_a": st.column_config.TextColumn("Protein A"),
            "annotation_a": st.column_config.TextColumn("Protein A annotation"),
            "protein_b": st.column_config.TextColumn("Protein B"),
            "annotation_b": st.column_config.TextColumn("Protein B annotation"),
        },
        score_fields=[(explore_score_col, explore_score_label)],
        default_top_n=min(500, len(selected_df)) or 10,
        key_prefix="brush_",
    )

st.divider()
st.header("Find candidates directly")
st.write(
    "A more direct version of the plot above: set a threshold and list every pair on the surprising side of "
    "it, without needing to brush-select."
)

st.subheader("High score, not annotated")
st.caption("Possible false negatives in the annotation -- a strong score but no known complex co-membership.")
score_values = universe_df[explore_score_col].dropna()
default_high = float(score_values.quantile(0.95))
high_thresh = st.slider(
    "Minimum score", float(score_values.min()), float(score_values.max()), default_high, key="high_thresh"
)
high_candidates = universe_df[(universe_df[explore_score_col] >= high_thresh) & ~ground_truth]
st.caption(f"{len(high_candidates):,} pairs match.")
ui.render_ranked_pair_table(
    high_candidates,
    sort_col=explore_score_col,
    column_order=[explore_score_col, "protein_a", "annotation_a", "protein_b", "annotation_b"],
    column_config={
        explore_score_col: st.column_config.NumberColumn(explore_score_label, format="%.3f"),
        "protein_a": st.column_config.TextColumn("Protein A"),
        "annotation_a": st.column_config.TextColumn("Protein A annotation"),
        "protein_b": st.column_config.TextColumn("Protein B"),
        "annotation_b": st.column_config.TextColumn("Protein B annotation"),
    },
    score_fields=[(explore_score_col, explore_score_label)],
    default_top_n=min(500, len(high_candidates)) or 10,
    key_prefix="highcand_",
)

st.divider()
st.subheader("Low score, annotated true positive")
st.caption("Known interactions the score ranks poorly -- worth a second look at either.")
default_low = float(score_values.quantile(0.05))
low_thresh = st.slider(
    "Maximum score", float(score_values.min()), float(score_values.max()), default_low, key="low_thresh"
)
low_candidates = universe_df[(universe_df[explore_score_col] <= low_thresh) & ground_truth]
st.caption(f"{len(low_candidates):,} pairs match.")
ui.render_ranked_pair_table(
    low_candidates,
    sort_col=explore_score_col,
    column_order=[explore_score_col, "protein_a", "annotation_a", "protein_b", "annotation_b"],
    column_config={
        explore_score_col: st.column_config.NumberColumn(explore_score_label, format="%.3f"),
        "protein_a": st.column_config.TextColumn("Protein A"),
        "annotation_a": st.column_config.TextColumn("Protein A annotation"),
        "protein_b": st.column_config.TextColumn("Protein B"),
        "annotation_b": st.column_config.TextColumn("Protein B annotation"),
    },
    score_fields=[(explore_score_col, explore_score_label)],
    default_top_n=min(500, len(low_candidates)) or 10,
    sort_ascending=True,
    key_prefix="lowcand_",
)
