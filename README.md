# Eukaryoma-PPI-webserver

Streamlit webserver for browsing pairwise protein-protein interactions extracted
from pooled AlphaFold3 predictions. Each pool jointly predicts several proteins
at once (a "pool"); this tool extracts and renders the two chains for any
specific protein pair contained in a pool. Loosely modelled on
[mutfunc](https://github.com/jurgjn/mutfunc)'s local-lookup + 3D-viewer UX.

Ways to pick a pair: **Browse Pairs** ranks AF3 pool pairs by predicted
interaction score; **Pair Viewer** looks up a specific pair by protein id *or*
by matching text in its annotation (e.g. searching "kinase" finds every
protein whose description mentions it); **Coabundance**/**Cofractionation**/
**Phyloprofiling** rank pairs by each external association-score source;
**Unified Ranking** combines all sources into one table; and **Annotations**
checks all of the above against known true-positive interactions. Selecting a
row anywhere shows that pair's AF3 structure when one has been predicted, with
a download button for the extracted PDB. Pairs known to be true positives
(see below) are highlighted in every table.

**Protein View** flips this around: pick one protein of interest and see
every interaction it's part of, external database links, its eggNOG-mapper
functional annotation, and a small interactome graph centered on it (see
below).

## Data layout

The webserver expects a `data/` directory (kept out of git, since it's large)
with:

```
data/
├── report_file.tsv        # <pool_name>\t<protein_1>_<protein_2>_..._<protein_n>
├── recap_set0_pairs.tsv    # per-pair AF3 interaction scores (see below)
├── OMAfiltered_..._annotated.fa  # protein annotations (see below)
├── other_data_sources/      # optional external association-score matrices (see below)
│   ├── coabundance/latest_coabundance_matrix.csv
│   ├── cofractionation/latest_cofrac_matrix.csv
│   ├── phyloprofiling/Capsaspora_hogprof_interactions_with_zeros.parquet
│   └── annotations/          # true-positive complex annotations (see below)
│       ├── corum/corum_annotations.tsv
│       └── marcotte/marcotte_annotations.txt
│   └── eggnogg_annotations/MICH_Capsaspora_owczarzaki_Schultz_A.tsv  # eggNOG-mapper output (see below)
├── pools/                   # <pool_name>.fcz, one pooled AF3 prediction per pool
├── structures/               # generated: <pool_name>.cif, decompressed by build_data.py
└── index/                     # generated: pools.parquet, pairs.parquet,
                                #            protein_annotations.parquet, eggnog_annotations.parquet,
                                #            universe_scores.parquet, true_positive_pairs.parquet
```

By default the app looks for `data/` as a sibling of this repo checkout
(`../../data`, matching the current on-disk layout). Override with the
`EUKARYOMA_DATA_DIR` environment variable if your layout differs (e.g. `/data`
inside Docker).

The protein order in `report_file.tsv` matches the chain order (A, B, C, ...)
in the pool's decompressed mmCIF structure.

`recap_set0_pairs.tsv` has one row per (pair, AF3 sample) -- each pool is
predicted with several seeds/samples, so a given pair appears several times
with slightly different `corrected_chain_pair_iptm` values.
`scripts/build_data.py` averages these into a single score per (pool, pair)
and stores it in `data/index/pairs.parquet`, which both pair pages display.
Override its path with `EUKARYOMA_RECAP_FILE`.

Protein annotations come from the FASTA headers in
`OMAfiltered_MSfiltered_MICH_Capsaspora_owczarzaki_Schultz_A_annotated.fa`,
e.g. `>XP_004340666.2,4-aminobutyrate:2-oxoglutarate transaminase
activity,Aminotran_3`. The first field is an NCBI-style id that
`eukaryoma_ppi.annotations.fasta_id_to_protein_id` converts to the website's
id format by dropping `.`/`_` and lowercasing (`XP_004340666.2` ->
`xp0043406662`); the remaining comma-separated fields are that protein's
annotations. `scripts/build_data.py` parses this into
`data/index/protein_annotations.parquet` and reports any website protein id
missing a FASTA match (currently: none, 2145/2145 matched). Override its path
with `EUKARYOMA_FASTA_FILE`.

## External data sources and the universe table

Beyond the AF3 pool structures, four optional sources give an association
score for arbitrary protein pairs (not just ones AF3 happened to pool
together): **coabundance**, **cofractionation**, and **phylogenetic
profiling** (HogProf) matrices, plus the CORUM/Marcotte true-positive
annotations (see below). Each of the three score matrices ships as a dense
protein x protein matrix, but with different formats and id conventions:

- coabundance (CSV): row/column labels are FASTA-header style, e.g.
  `"XP_004340666.2,4-aminobutyrate..."`.
- cofractionation (CSV): labels are a bare NCBI accession with no `XP_`
  prefix, e.g. `"004340666.2"`, and some labels are `;`-joined groups of
  proteins the experiment couldn't distinguish, e.g.
  `"004340769.1;004346401.1"` -- the group's row/column score applies to
  every member id.
- phyloprofiling / HogProf (parquet): labels are FASTA-header-style
  accessions with an extra trailing `_<taxon id>`, e.g.
  `"XP_004349908.1_595528"`. A `0.0` score is mostly a real measurement (no
  phylogenetic co-evolution signal detected) -- 97% of covered pairs score
  exactly `0.0` -- but the source matrix also writes each pair's score in
  only *one* of the `(A, B)`/`(B, A)` directions, leaving the other at that
  same default `0.0`, i.e. indistinguishable from a real zero without
  cross-checking the mirror. `build_source_pairs` (`symmetrize_zero_gaps`)
  fills each `0.0` from its mirror position whenever the mirror is non-zero
  (13,229 pairs affected among the website's proteins as of the current
  file; verified no pair had disagreeing non-zero values on both sides, so
  this never overwrites real data) -- the other two sources are already
  symmetric correlation matrices where a literal `0.0` is a genuine
  measurement, so they're left alone.

`eukaryoma_ppi.external_scores.normalize_external_id` converts all three id
formats to the website's id format, the same way as the FASTA annotations.
Coverage isn't complete: coabundance covers 2139/2145 website proteins,
cofractionation all 2145/2145, phyloprofiling 2003/2145.

`scripts/build_data.py` builds **`data/index/universe_scores.parquet`**: one
row for every possible pair among the website's proteins (~2.3M for 2145
proteins, since these are dense matrices, not just the 138,295 pairs AF3
happened to pool), with a score column per source (`NaN` where a source
doesn't cover that pair) plus a `unified_score` combining whichever sources
have data for that pair (`eukaryoma_ppi.external_scores.compute_unified_score`
-- currently the mean percentile rank across available sources; a
placeholder documented as easy to swap for a different combination method
later). The home page's "Data source coverage" section summarizes how many
pairs each source (and combination of sources) covers.

Override source file paths with `EUKARYOMA_COABUNDANCE_FILE`,
`EUKARYOMA_COFRACTIONATION_FILE`, `EUKARYOMA_PHYLOPROFILING_FILE`. Any of the
three can be swapped for an updated matrix (CSV or parquet, same row/column
format) at the configured path -- just re-run `scripts/build_data.py` to
pick it up everywhere it's used.

## True-positive annotations (CORUM / Marcotte)

**CORUM** and **Marcotte** list human protein complexes together with each
member's Capsaspora ortholog(s) (`;`-joined when a human gene has several
paralogs, `N/A` when it has none), e.g.:

```
corumID  category                          ...  Capsaspora
4        Multisubunit ACTR coactivator...       XP_004345483.1_595528;XP_004364839.1_595528
```

Two Capsaspora proteins are treated as a **known true-positive interaction**
if they co-occur as members of the same group: CORUM groups by `corumID`
(one complex), Marcotte by `category` (a coarser functional module). Ortholog
ids carry a trailing `_<taxon id>` the website's ids don't have, so
`eukaryoma_ppi.complex_annotations` strips that before normalizing them the
same way as the FASTA annotations. Parsing follows the approach in
`/Users/spascare/data/work/beltrao/combine_assoc_scores/Eukaryoma_PPI_analysis/src/cofrac_coab9partial_tmp_test.py`,
except a pair counts as true-positive if it shares membership in *any* group
(that script keeps only one category per protein and can miss some
co-membership pairs for proteins in more than one complex).

`scripts/build_data.py` builds `data/index/true_positive_pairs.parquet`
(`corum_tp`/`marcotte_tp` booleans, restricted to website proteins: 7,185 /
2,996 pairs respectively, 8,935 flagged by either, 1,246 by both) and merges
it into both `pairs.parquet` and `universe_scores.parquet`. Every pair table
in the app highlights true-positive rows (green/blue/purple for
CORUM/Marcotte/both); the **Annotations** page plots true-positive rate by
score quantile for each score, and a jittered scatter of true-positive pairs
per score plotted against a grey violin of the not-annotated population's
score distribution (`eukaryoma_ppi.tp_analysis.baseline_score_density` -- a
histogram computed over the full not-annotated population, not a sample,
since that stays fast even at millions of rows) as a baseline for whether
true-positive scores are actually higher, not just relative to each other.
Both the scatter and the violin support brush-select or threshold-filter for
pairs where a score and the annotation disagree (high score without
annotation, or low score despite it) -- candidates for annotation false
negatives or under-ranked real interactions. Override file paths with
`EUKARYOMA_CORUM_FILE` / `EUKARYOMA_MARCOTTE_FILE`.

## Protein View and eggNOG-mapper annotation

**Protein View** flips the browsing model from interaction-centric to
protein-centric: pick one protein and see everything the site knows about it
in one place -- every pair it appears in (sortable/filterable by any score,
same true-positive highlighting as elsewhere), direct links out to UniProt,
NCBI Protein and the AlphaFold DB, its eggNOG-mapper functional annotation,
and a small interactome graph centered on it.

The external-database links and the functional annotation (GO terms, KEGG
orthologs/pathways, PFAM domains, COG category, eggNOG orthologous groups)
come from **eggNOG-mapper** output (`eukaryoma_ppi.eggnog`), parsed into
`data/index/eggnog_annotations.parquet` by `scripts/build_data.py`. This is a
different source from the OMA/FASTA description used elsewhere in the app
(`protein_annotations.parquet`), so the UI always labels it "eggNOG" to keep
the two apart. Override the input file with `EUKARYOMA_EGGNOG_FILE`.

eggNOG's `seed_ortholog` column gives an accession to link out with, but for
~99% of Capsaspora proteins that "ortholog" is just the protein's own
existing NCBI RefSeq entry (a self-hit, since Capsaspora is already in
eggNOG's reference databases) rather than a genuine cross-species UniProt
match. `eukaryoma_ppi.eggnog` tells the two apart by accession shape
(RefSeq: `XP_`/`NP_` with an underscore; UniProt: no underscore) and links
accordingly: a RefSeq accession gets a direct NCBI Protein link but only
*search* links for UniProt/AlphaFold DB (no direct id to look up there),
while a genuine UniProt accession gets direct links to all three.

The interactome graph starts small: a slider sets a unified-score cutoff and
only partners scoring at or above it are auto-included as nodes (edge
thickness encodes unified score). Selecting rows in the interactions table
above adds those specific partners to the graph regardless of score, for
digging into a lower-confidence interaction of interest without lowering the
threshold for everyone else. Layout is a plain circular placement (the
protein of interest at the center, partners spaced evenly around it) computed
by hand with numpy -- no networkx/graphviz dependency, keeping the app
pip-installable per `COLLABORATOR_SETUP.md`.

## Downloading structures and plots

Every 3D structure viewer (Pair Viewer, Browse Pairs, the per-source ranking
pages, Unified Ranking, Protein View) has a **Download structure (.pdb)**
button right below it, serving the same extracted two-chain PDB the viewer
renders -- useful for figures or re-analysis outside the browser.

Every chart in the app (score plots, the Annotations scatter/violin, the
Protein View interactome) is a Vega-Embed chart under the hood, which ships
its own **"..." menu** in the top-right corner of the chart with "Save as
SVG"/"Save as PNG" built in -- no separate download button is needed for
plots.

## The `.fcz` format and `bin/foldcomp`

The pool structures are stored with [foldcomp](https://github.com/steineggerlab/foldcomp)
compression, but as a **custom internal build**: these files use the magic
number `FCZC`, not the public release's `FCMP`, so neither the `foldcomp` PyPI
package nor the official GitHub release binaries can decompress them. A
working binary is bundled at [bin/foldcomp](bin/foldcomp) (macOS universal
binary) and used only by `scripts/build_data.py`, via subprocess -- the
Streamlit app itself never calls foldcomp. If you need to run data conversion
on Linux/HPC, point `FOLDCOMP_BIN` at an equivalent Linux build of the same
fork (e.g. whatever produced [src/foldcomp_cmd.sh](../../src/foldcomp_cmd.sh)'s
output on the cluster).

## Setup

```bash
cd repo/Eukaryoma-PPI-webserver
source venv/bin/activate          # venv already exists in this repo
pip install -r requirements.txt
pip install -e .
```

## Convert the data (one-time, local, outside Docker)

```bash
python scripts/build_data.py
```

This decompresses every pool in `data/pools/*.fcz` to `data/structures/*.cif`
via `bin/foldcomp`, builds `data/index/pools.parquet` and
`data/index/pairs.parquet` from `report_file.tsv` (restricted to pools that
actually have a structure on disk), parses protein annotations and CORUM/
Marcotte true-positive flags, and builds `data/index/universe_scores.parquet`
from whichever external data sources are present (see above), printing a
coverage summary. Re-run with `--force` to redo decompression, or
`--skip-decompress` to only rebuild the indexes.

## Run locally

```bash
streamlit run app.py
```

## Run with Docker

The Docker image only serves already-converted data (`data/structures/` +
`data/index/`) -- it does not bundle or run `foldcomp` at all, since
`bin/foldcomp` is a macOS binary and won't run in a Linux container. Run
`scripts/build_data.py` on the host first, then:

```bash
docker build -t eukaryoma-ppi-webserver .
docker run -p 8501:8501 -v /absolute/path/to/data:/data eukaryoma-ppi-webserver
```

Then open http://localhost:8501.
