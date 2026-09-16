import altair as alt
import numpy as np
import pandas as pd
import streamlit as st

from eukaryoma_ppi import annotations, external_scores, tp_analysis, ui
from eukaryoma_ppi.config import ANNOTATIONS_INDEX_FILE, TRUE_POSITIVE_INDEX_FILE, UNIVERSE_INDEX_FILE

st.set_page_config(page_title="Annotations - Eukaryoma PPI", page_icon="🧬", layout="wide")
st.title("Annotations")
st.write(
    "How well do the association scores agree with CORUM/Marcotte complex co-membership -- treated here as "
    "known true-positive interactions? Use this page to check score quality and to spot pairs where the "
    "score and the annotation disagree."
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


universe_df = get_universe_for_annotations()

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
    "Each point is a pair, jittered by annotation category. Drag a box over a region (e.g. high scores in "
    "the \"Not annotated\" row) to list those pairs below -- candidates for false negatives in the "
    "annotation, or, in the true-positive rows at low score, cases the score under-ranks."
)

explore_score_col, explore_score_label = st.selectbox(
    "Score to explore", SCORE_COLUMNS, format_func=lambda pair: pair[1], key="explore_score"
)

MAX_PLOTTED = 4000


@st.cache_data
def get_plot_sample(_universe_df, score_col, max_plotted, seed=0):
    valid = _universe_df[score_col].notna()
    plot_df = _universe_df.loc[valid, ["protein_a", "protein_b", score_col]].copy()
    plot_df["tp_category"] = tp_analysis.tp_category(_universe_df.loc[valid])

    if len(plot_df) > max_plotted:
        annotated = plot_df[plot_df["tp_category"] != "Not annotated"]
        not_annotated = plot_df[plot_df["tp_category"] == "Not annotated"]
        budget = max(max_plotted - len(annotated), 0)
        if len(not_annotated) > budget:
            not_annotated = not_annotated.sample(budget, random_state=seed)
        plot_df = pd.concat([annotated, not_annotated])

    rng = np.random.default_rng(seed)
    plot_df["jitter"] = rng.uniform(-0.4, 0.4, size=len(plot_df))
    return plot_df.reset_index(drop=True)


plot_df = get_plot_sample(universe_df, explore_score_col, MAX_PLOTTED)
st.caption(f"Showing {len(plot_df):,} of {int(universe_df[explore_score_col].notna().sum()):,} pairs with a {explore_score_label} score.")

brush = alt.selection_interval(name="brush", encodings=["x", "y"])
category_order = ["Not annotated", "CORUM", "Marcotte", "Both"]
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
    .add_params(brush)
    .properties(height=320)
)

event = st.altair_chart(scatter, on_select="rerun", selection_mode="brush", key="annotation_scatter")

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

col_high, col_low = st.columns(2)

with col_high:
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

with col_low:
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
