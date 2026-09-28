import altair as alt
import streamlit as st

from eukaryoma_ppi import complex_annotations, index
from eukaryoma_ppi.config import PAIRS_INDEX_FILE, TRUE_POSITIVE_INDEX_FILE

st.set_page_config(page_title="Annotation Sources - Eukaryoma PPI", page_icon="🧬", layout="wide")
st.title("Annotation Sources")
st.write(
    "Where the site's \"true-positive\" interaction annotations come from: CORUM and Marcotte protein "
    "complexes, mapped from human onto Capsaspora orthologs (OMA/FastOMA). Two Capsaspora proteins are "
    "treated as a known true interaction if they co-occur as members of the same complex/category. For how "
    "well the association scores actually agree with these annotations, see the **Annotated Pairs** page."
)

if not PAIRS_INDEX_FILE.exists() or not TRUE_POSITIVE_INDEX_FILE.exists():
    st.error("No pair/true-positive index found. Run `python scripts/build_data.py` first (see README.md).")
    st.stop()

SOURCE_LABELS = {"corum_tp": "CORUM", "marcotte_tp": "Marcotte"}
SOURCE_DESCRIPTIONS = {
    "corum_tp": (
        "Curated, literature-backed human protein complexes. Grouped by `corumID` -- each distinct complex "
        "is its own group, shown below by its actual complex name (e.g. \"Ribosome, cytoplasmic\")."
    ),
    "marcotte_tp": (
        "Coarse functional modules from the Marcotte lab's complex map, grouped by `category` -- a handful "
        "of broad modules (e.g. one \"ribosome\" group covering the whole ribosome), not individual named "
        "complexes. A finer-grained `category-specific` column exists in the source file but isn't used here."
    ),
}
LABEL_COLUMNS = {"corum_tp": "category", "marcotte_tp": None}


@st.cache_data
def get_website_ids():
    pairs_df = index.load_pairs_index()
    return set(pairs_df["protein_a"]) | set(pairs_df["protein_b"])


@st.cache_data
def get_group_stats(source_key):
    path, group_col = complex_annotations.SOURCES[source_key]
    return complex_annotations.group_stats(path, group_col, get_website_ids(), label_col=LABEL_COLUMNS[source_key])


@st.cache_data
def get_total_groups(source_key):
    path, group_col = complex_annotations.SOURCES[source_key]
    return complex_annotations.total_groups_in_file(path, group_col)


@st.cache_data
def get_true_positive_flags():
    return complex_annotations.load_true_positive_flags()


website_ids = get_website_ids()
flags = get_true_positive_flags()

st.divider()
st.header("Overview")
n_either = int(flags[complex_annotations.FLAG_COLUMNS].any(axis=1).sum())
n_both = int(flags[complex_annotations.FLAG_COLUMNS].all(axis=1).sum())
metric_cols = st.columns(len(complex_annotations.FLAG_COLUMNS) + 2)
for col_widget, flag_col in zip(metric_cols, complex_annotations.FLAG_COLUMNS):
    col_widget.metric(f"{SOURCE_LABELS[flag_col]} true-positive pairs", f"{int(flags[flag_col].sum()):,}")
metric_cols[-2].metric("Either source", f"{n_either:,}")
metric_cols[-1].metric("Both sources", f"{n_both:,}")
st.caption(f"Out of {len(website_ids):,} website proteins ({len(website_ids) * (len(website_ids) - 1) // 2:,} possible pairs).")

for flag_col in complex_annotations.FLAG_COLUMNS:
    label = SOURCE_LABELS[flag_col]
    st.divider()
    st.header(label)
    st.write(SOURCE_DESCRIPTIONS[flag_col])

    stats_df = get_group_stats(flag_col)
    total_groups = get_total_groups(flag_col)
    usable = stats_df[stats_df["n_website_members"] >= 2]

    metric_cols = st.columns(4)
    metric_cols[0].metric(f"{label} complexes (total)", f"{total_groups:,}")
    metric_cols[1].metric("With ≥1 Capsaspora ortholog", f"{len(stats_df):,}")
    metric_cols[2].metric("Usable for pairs (≥2 website proteins)", f"{len(usable):,}")
    metric_cols[3].metric("True-positive pairs contributed", f"{int(flags[flag_col].sum()):,}")

    if usable.empty:
        st.info("No complex in this source has 2+ website proteins, so it contributes no true-positive pairs.")
        continue

    st.caption(
        f"Average complex size, among the {len(usable):,} usable complexes: "
        f"{usable['n_website_members'].mean():.1f} website proteins. "
        f"Median pairs contributed per complex: {usable['n_pairs'].median():.0f}."
    )

    left, right = st.columns(2)
    with left:
        if len(usable) >= 15:
            chart = (
                alt.Chart(usable)
                .mark_bar()
                .encode(
                    x=alt.X("n_website_members:Q", bin=alt.Bin(maxbins=30), title="Website proteins per complex"),
                    y=alt.Y("count():Q", title="Number of complexes"),
                )
                .properties(height=280, title=f"{label}: complex size distribution")
            )
            st.altair_chart(chart, width="stretch")
        else:
            st.write(f"Only {len(usable)} usable complexes -- listing all of them instead of a histogram:")
            st.dataframe(
                usable.sort_values("n_pairs", ascending=False)[["label", "n_website_members", "n_pairs"]],
                hide_index=True,
                width="stretch",
                column_config={
                    "label": st.column_config.TextColumn(f"{label} group"),
                    "n_website_members": st.column_config.NumberColumn("Website proteins"),
                    "n_pairs": st.column_config.NumberColumn("Pairs contributed"),
                },
            )

    with right:
        top10 = usable.sort_values("n_pairs", ascending=False).head(10)
        chart = (
            alt.Chart(top10)
            .mark_bar()
            .encode(
                x=alt.X("n_pairs:Q", title="Pairs contributed"),
                y=alt.Y("label:N", sort="-x", title=f"{label} group"),
                tooltip=["label", "n_website_members", "n_pairs"],
            )
            .properties(height=280, title=f"{label}: largest complexes by pairs contributed")
        )
        st.altair_chart(chart, width="stretch")
