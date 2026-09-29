# Eukaryoma-PPI-webserver

Streamlit webserver for browsing pairwise protein-protein interactions extracted
from pooled AlphaFold3 predictions. Each pool jointly predicts several proteins
at once (a "pool"); this tool extracts and renders the two chains for any
specific protein pair contained in a pool. Loosely modelled on
[mutfunc](https://github.com/jurgjn/mutfunc)'s local-lookup + 3D-viewer UX.

The home page (`app.py`) lists every page with a one-line description and a
direct link to it -- a good starting point for getting oriented.

Ways to pick a pair: **Browse Pairs** ranks AF3 pool pairs by predicted
interaction score; **Pair Viewer** looks up a specific pair by protein id *or*
by matching text in its annotation (e.g. searching "kinase" finds every
protein whose description mentions it); **Coabundance**/**Cofractionation**/
**Phyloprofiling** rank pairs by each external association-score source;
**Unified Ranking** combines all sources into one table; **Annotation
Sources** describes where the CORUM/Marcotte true-positive annotations come
from and how many complexes/pairs each contributes; and **Annotated Pairs**
checks the scores against those annotations. Selecting a row anywhere shows
that pair's AF3 structure when one has been predicted, with a download
button for the extracted PDB and an optional live-fetched reference
structure from AlphaFold DB. Pairs known to be true positives (see below)
are highlighted in every table.

**Protein View** flips this around: pick one protein of interest and see
every interaction it's part of, external database links, its eggNOG-mapper
functional annotation, an interactome graph centered on it with GO-term
enrichment of whatever's currently in the graph, and a structure comparison
across every pool the protein was predicted in (see below).

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
│   ├── eggnogg_annotations/MICH_Capsaspora_owczarzaki_Schultz_A.tsv  # eggNOG-mapper output (see below)
│   ├── yeast_hog2gene/*.pkl        # Capsaspora<->yeast HOG map (see below)
│   └── yeast_structures/<accession>.fcz  # foldcomp-compressed yeast AF3 monomers (see below)
├── pools/                   # <pool_name>.fcz, one pooled AF3 prediction per pool
├── structures/               # generated: <pool_name>.cif, decompressed by build_data.py
│   └── yeast/<accession>.pdb  # generated: only the referenced yeast structures, decompressed
└── index/                     # generated: pools.parquet, pairs.parquet,
                                #            protein_annotations.parquet, eggnog_annotations.parquet,
                                #            universe_scores.parquet, true_positive_pairs.parquet,
                                #            yeast_orthologs.parquet
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
have data for that pair (`eukaryoma_ppi.external_scores.compute_unified_score`).
The home page's "Data source coverage" section summarizes how many pairs
each source (and combination of sources) covers.

### Unified score: weightless vs. weighted

Each source's score is turned into a percentile rank (0-1, computed only
over the pairs it covers), and a pair's `unified_score` is the mean of its
available per-column ranks -- this is the persisted, on-disk default
("weightless": every available source counts equally).

The **Unified Ranking** page can switch to **weighted** instead: each
source is weighted by its own raw (pre-rank) standard deviation, normalized
to sum to 1 (`eukaryoma_ppi.external_scores.compute_default_weights`) -- a
source whose values barely vary can't discriminate between pairs (e.g. an
AF3 ipTM that's ~0.95 for nearly everyone), so it counts for less. Rows
missing a given source have that source's weight dropped and the rest
renormalized, the same way the weightless mean already handles missing
columns.

Switching to "weighted" also reveals one slider per source, seeded from the
computed defaults -- drag any of them to override that source's weight (a
"Reset to computed defaults" button clears all the overrides). Weights
don't need to sum to 1; `compute_unified_score` renormalizes per pair over
whichever sources that pair actually has, so only their *relative*
proportions matter.

This choice is a **session-wide setting**, not a per-page one: picking
"weighted" (and any slider overrides) on the Unified Ranking page changes
what every other page shows too (Annotation Sources, Annotated Pairs,
Protein View), via `st.session_state` keys every page reads
(`eukaryoma_ppi.ui.get_unified_score_mode` /
`get_active_unified_score_weights`). Recomputing the weighted rank over all
~2.3M pairs takes under a second, so each page just redoes it in memory
(cached per mode/weights, so repeat page visits are instant) rather than
needing a second column on disk.

Getting this to actually survive page navigation needs a small workaround:
Streamlit clears a *widget's* session_state entry whenever that widget
isn't instantiated on a script run -- true for the mode selectbox and every
weight slider on every page except Unified Ranking. Each is therefore kept
in a second, plain session_state key that isn't tied to any widget, and the
widget is reseeded from it right before being (re)created.

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
CORUM/Marcotte/both), including edges in Protein View's interactome graph.
Override file paths with `EUKARYOMA_CORUM_FILE` / `EUKARYOMA_MARCOTTE_FILE`.

### Annotation Sources: what's actually in these files

CORUM's 2,655 complexes and Marcotte's 7 categories are very different in
kind, and the **Annotation Sources** page exists to make that visible before
trusting the true-positive counts above: `eukaryoma_ppi.complex_annotations.
group_stats` reports, per complex/category, how many Capsaspora orthologs it
lists, how many of those are website proteins, and how many pairs it alone
contributes (`n_website_members choose 2` -- before deduplicating against
pairs shared with other groups, so these per-group counts don't sum exactly
to the final `corum_tp`/`marcotte_tp` totals above).

Only 760 of CORUM's 1,880 orthology-mapped complexes have 2+ website
proteins (median 1 pair contributed -- most complexes are small); Marcotte
groups by the coarse `category` column (not the finer `category-specific`),
leaving only 6 categories with any Capsaspora ortholog and 5 usable for
pairs -- one of them, a single ~66-protein "ribosome" group, alone accounts
for 2,145 of Marcotte's 2,996 true-positive pairs. The page shows a size-
distribution histogram plus the largest complexes by pairs contributed for
CORUM, and the full (tiny) group table for Marcotte.

### Annotated Pairs: scores vs. annotations

The **Annotated Pairs** page checks the association scores against the
CORUM/Marcotte annotations described above: true-positive rate by score
quantile for each score, and a jittered scatter of true-positive pairs per
score plotted against a grey violin of the not-annotated population's score
distribution (`eukaryoma_ppi.tp_analysis.baseline_score_density` -- a
histogram computed over the full not-annotated population, not a sample,
since that stays fast even at millions of rows) as a baseline for whether
true-positive scores are actually higher, not just relative to each other.
Both the scatter and the violin support brush-select or threshold-filter for
pairs where a score and the annotation disagree (high score without
annotation, or low score despite it) -- candidates for annotation false
negatives or under-ranked real interactions.

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
only partners scoring at or above it are auto-included as nodes (edge color
marks CORUM/Marcotte true positives, thickness encodes unified score).
Selecting rows in the interactions table above adds those specific partners
to the graph regardless of score, for digging into a lower-confidence
interaction of interest without lowering the threshold for everyone else.
Layout is a plain circular placement (the protein of interest at the center,
partners spaced evenly around it) computed by hand with numpy -- no
networkx/graphviz dependency, keeping the app pip-installable per
`COLLABORATOR_SETUP.md`.

### GO term enrichment of an interactome

Right under the interactome graph, `eukaryoma_ppi.go_enrichment` tests
whether any GO term is over-represented among the graph's current proteins,
against every eggNOG-annotated website protein as background: a plain
hypergeometric test per term (`scipy`/`goatools`-free -- just Python's
`math.comb`) with Benjamini-Hochberg FDR correction across all terms tested.
An adjustable FDR threshold filters the results table; an expander below it
always shows the simpler fallback -- just the GO terms present in the set
and how many proteins carry each, no statistics -- for when nothing reaches
significance (small interactomes often won't) or a plain recap is more
useful than a test. Terms are linked out to QuickGO for their name/
definition, since this app doesn't bundle a GO ontology (OBO) file to look
that up locally.

### Comparing structures across pools

Below that, **"Compare structures across pools"** answers a different
question: not how a protein interacts, but how *consistently AF3 predicts
its own fold* across the different pools it was pooled into (each pool has
a different partner, so a different context). Pick 2-8 of the pools this
protein has a structure in, and either superpose them (Biopython's
`Superimposer`, CA-atom RMSD alignment onto the first selected pool, then
all rendered together in one py3Dmol view, one color per pool) or view them
side by side, unaligned. A pool is skipped from superposition (with a
warning, not a crash) if its chain has a different CA count than the
reference pool's -- `eukaryoma_ppi.structures.superpose_chains_as_pdb`.

### Reference structures: human (AlphaFold DB) and yeast (local AF3)

Wherever a pair's AF3 structure is shown (every pair page, Unified Ranking,
Protein View), an optional **reference structure** dropdown appears
alongside it when either protein has a human and/or yeast ortholog --
comparing this app's AF3 *pair* prediction against an ortholog's own solo
structure, predicted completely independently.

**Human**, via a live fetch from AlphaFold DB: the accession comes from
CORUM/Marcotte's own `uniprot` column, not eggNOG's `seed_ortholog` --
eggNOG's best-hit search spans many reference proteomes and is only very
rarely an actual human UniProt accession for this species, whereas
CORUM/Marcotte already curate the human gene each Capsaspora ortholog
corresponds to (`eukaryoma_ppi.complex_annotations.human_uniprot_orthologs`;
~750 of the website's 2,145 proteins have at least one, some several when
they're the ortholog of multiple human paralogs). `eukaryoma_ppi.
external_structures` calls AFDB's prediction API for the current `pdbUrl`
rather than guessing the file's model-version suffix directly
(`AF-<accession>-F1-model_v<N>.pdb` -- `N` changes release to release, so a
hardcoded version 404s eventually, as the original v4 guess did). Any fetch
failure (protein has no such ortholog, network error, accession not
actually in AFDB) is silent/graceful -- a nice-to-have next to the pair's
own structure, not something the rest of the page depends on.

**Yeast**, from AF3 monomer predictions run on the lab's cluster and
foldcomp-compressed locally (`data/other_data_sources/yeast_structures/
<accession>.fcz`, one per lowercase UniProt accession). The
Capsaspora-to-yeast ortholog call comes from a HOG (orthologous group) map
covering Capsaspora plus five other species
(`data/other_data_sources/yeast_hog2gene/*.pkl`, keyed by NCBI taxon id --
192875 for Capsaspora, 4932 for *S. cerevisiae*): two genes sharing a HOG
under those two taxon keys are treated as orthologs
(`eukaryoma_ppi.yeast_orthologs.capsaspora_to_yeast_oma_ids`). The yeast
side of that map uses OMA's own internal per-species ids ("YEAST02308"),
not UniProt accessions, so each one needs resolving via OMA's REST API --
`resolve_all_yeast_uniprot_accessions` does this once at build time
(network-bound, a few minutes even with concurrency, retrying transient
502s from OMA's public server), caching the result to
`data/index/yeast_orthologs.parquet` (943 of 2,145 website proteins covered
as of the current files -- notably better than the human mapping, since
this is a direct Capsaspora<->yeast call rather than routing through human
complexes).

`scripts/build_data.py` then decompresses only the *referenced* yeast
accessions (not the full ~6,000-protein yeast set -- most of it would never
be looked up) into `data/structures/yeast/<accession>.pdb`, going through
mmCIF and Biopython's own `PDBIO` rather than decompressing straight to PDB
like the pools do: the foldcomp CLI's own PDB writer turned out to
misformat negative coordinates for these particular monomer files (a fixed-
column parser reads through corrupted data past the first negative x
coordinate -- verified directly by inspecting a raw decompressed file), while
its mmCIF output is whitespace-delimited and unaffected. This keeps the
Streamlit app itself foldcomp-free at runtime for yeast structures too, the
same invariant the pools already rely on (the container "never calls
foldcomp itself" -- see the Docker section below) -- `structures.
get_yeast_structure_pdb` is a plain file read. Override paths with
`EUKARYOMA_YEAST_HOG2GENE_FILE` / `EUKARYOMA_YEAST_STRUCTURES_DIR`; skip
this step (e.g. if the yeast files aren't available yet) with
`--skip-yeast-orthologs`, or force re-resolution with
`--force-yeast-orthologs`.

## Downloading structures and plots

Every 3D structure viewer (Pair Viewer, Browse Pairs, the per-source ranking
pages, Unified Ranking, Protein View) has a **Download structure (.pdb)**
button right below it, serving the same extracted two-chain PDB the viewer
renders -- useful for figures or re-analysis outside the browser.

Every chart in the app (score plots, the Annotation Sources/Annotated Pairs
histograms and scatter/violin, the Protein View interactome) is a Vega-Embed
chart under the hood, which ships
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
