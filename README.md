# Eukaryoma-PPI-webserver

Streamlit webserver for browsing pairwise protein-protein interactions extracted
from pooled AlphaFold3 predictions. Each pool jointly predicts several proteins
at once (a "pool"); this tool extracts and renders the two chains for any
specific protein pair contained in a pool. Loosely modelled on
[mutfunc](https://github.com/jurgjn/mutfunc)'s local-lookup + 3D-viewer UX.

Ways to pick a pair: **Browse Pairs** ranks AF3 pool pairs by predicted
interaction score; **Pair Viewer** looks up a specific pair by protein id;
**Coabundance**/**Cofractionation**/**Phyloprofiling** rank pairs by each
external association-score source; and **Unified Ranking** combines all
sources into one table. Selecting a row anywhere shows that pair's AF3
structure when one has been predicted.

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
│   └── phyloprofiling/latest_phyloprofiling_matrix.csv   # not provided yet
├── pools/                   # <pool_name>.fcz, one pooled AF3 prediction per pool
├── structures/               # generated: <pool_name>.cif, decompressed by build_data.py
└── index/                     # generated: pools.parquet, pairs.parquet,
                                #            protein_annotations.parquet, universe_scores.parquet
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

Beyond the AF3 pool structures, three optional sources give an association
score for arbitrary protein pairs (not just ones AF3 happened to pool
together): **coabundance**, **cofractionation**, and **phylogenetic
profiling** (HogProf; no file yet). Each ships as a dense protein x protein
correlation matrix CSV, but with different id conventions:

- coabundance row/column labels are FASTA-header style, e.g.
  `"XP_004340666.2,4-aminobutyrate..."`.
- cofractionation labels are a bare NCBI accession with no `XP_` prefix, e.g.
  `"004340666.2"`, and some labels are `;`-joined groups of proteins the
  experiment couldn't distinguish, e.g. `"004340769.1;004346401.1"` -- the
  group's row/column score applies to every member id.

`eukaryoma_ppi.external_scores.normalize_external_id` converts both to the
website's id format, the same way as the FASTA annotations. Coverage isn't
complete: coabundance covers 2139/2145 website proteins, cofractionation
covers all 2145/2145.

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
`EUKARYOMA_COFRACTIONATION_FILE`, `EUKARYOMA_PHYLOPROFILING_FILE`. Drop a
phyloprofiling matrix at the configured path (same CSV format) and re-run
`scripts/build_data.py` to populate that source everywhere it's used.

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
actually have a structure on disk), parses protein annotations, and builds
`data/index/universe_scores.parquet` from whichever external data sources are
present (see above), printing a coverage summary. Re-run with `--force` to
redo decompression, or `--skip-decompress` to only rebuild the indexes.

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
