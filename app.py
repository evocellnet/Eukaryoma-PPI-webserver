import streamlit as st

from eukaryoma_ppi import external_scores, index
from eukaryoma_ppi.config import PAIRS_INDEX_FILE, POOLS_INDEX_FILE, UNIVERSE_INDEX_FILE

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

st.divider()
st.header("Pages")
st.write("Every page in the sidebar, and what it's for:")

PAGES = [
    ("pages/1_Browse_Pairs.py", "🔍", "Browse Pairs", "Rank AF3 pool pairs by predicted interaction score (ipTM)."),
    (
        "pages/2_Pair_Viewer.py",
        "🔎",
        "Pair Viewer",
        "Look up one specific pair by protein id, or by matching text anywhere in its annotation.",
    ),
    (
        "pages/3_Coabundance.py",
        "📊",
        "Coabundance",
        "Rank pairs by co-abundance correlation across samples.",
    ),
    (
        "pages/4_Cofractionation.py",
        "📊",
        "Cofractionation",
        "Rank pairs by co-fractionation correlation.",
    ),
    (
        "pages/5_Phyloprofiling.py",
        "📊",
        "Phyloprofiling",
        "Rank pairs by phylogenetic co-occurrence (HogProf).",
    ),
    (
        "pages/6_Unified_Ranking.py",
        "🏆",
        "Unified Ranking",
        "Combine every source into one score -- equal-weighted by default, or weighted (with "
        "adjustable per-source sliders) by each score's own distribution.",
    ),
    (
        "pages/7_Annotation_Sources.py",
        "📚",
        "Annotation Sources",
        "What the CORUM/Marcotte true-positive annotations actually cover: complex/category counts, "
        "sizes, and how many pairs each one contributes.",
    ),
    (
        "pages/8_Annotated_Pairs.py",
        "🎯",
        "Annotated Pairs",
        "Check how well each score agrees with the CORUM/Marcotte annotations, and find pairs where "
        "they disagree.",
    ),
    (
        "pages/9_Protein_View.py",
        "🧬",
        "Protein View",
        "Pick one protein: every interaction it's in, external database links, its eggNOG annotation, "
        "an interactome graph with GO-term enrichment, and structure comparison across pools.",
    ),
]

page_cols = st.columns(3)
for i, (path, icon, title, description) in enumerate(PAGES):
    with page_cols[i % 3]:
        st.page_link(path, label=f"{icon} {title}")
        st.caption(description)

st.divider()
st.header("Data source coverage")

@st.cache_data
def get_universe_for_recap():
    return external_scores.load_universe()


@st.cache_data
def get_source_presence_summary(_universe_df):
    return external_scores.source_presence_summary(_universe_df)


if not UNIVERSE_INDEX_FILE.exists():
    st.warning("No universe index found. Run `python scripts/build_data.py` to include external data sources.")
else:
    universe_df = get_universe_for_recap()
    st.write(
        f"Beyond the AF3 pool structures above, {len(universe_df):,} possible pairs among the same "
        "proteins have been cross-referenced against optional external association-score sources "
        "(coabundance, cofractionation, phylogenetic profiling). Coverage differs by source -- see below."
    )

    pattern_counts, n_sources_counts = get_source_presence_summary(universe_df)
    # Max sources actually observed on any pair == how many sources currently
    # have data at all (a source with no file/no coverage never contributes).
    max_sources_present = int(n_sources_counts.index.max())
    all_sources_count = int(n_sources_counts.get(max_sources_present, 0))

    view = st.selectbox("View coverage by", ["Number of sources present", "Exact combination of sources"])
    if view == "Number of sources present":
        st.bar_chart(n_sources_counts.rename("n_pairs"), x_label="Sources present", y_label="Pairs")
        st.caption(f"{all_sources_count:,} pairs are backed by all {max_sources_present} currently available source(s).")
    else:
        st.dataframe(
            pattern_counts,
            column_config={
                "sources_present": st.column_config.TextColumn("Sources present"),
                "n_pairs": st.column_config.NumberColumn("Pairs", format="%d"),
            },
            hide_index=True,
            width="stretch",
        )
