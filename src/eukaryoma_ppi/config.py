"""Central paths for the Eukaryoma PPI webserver.

All paths are overridable via environment variables so the same code works
for local development (data/ as a sibling of the repo) and inside Docker
(data/ mounted as a volume, e.g. at /data).
"""

import os
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]

# Local dev default: <project>/data, sibling of the repo checkout.
# Override with EUKARYOMA_DATA_DIR (e.g. "/data" inside the Docker container).
DATA_DIR = Path(os.environ.get("EUKARYOMA_DATA_DIR", REPO_ROOT.parent.parent / "data"))

REPORT_FILE = Path(os.environ.get("EUKARYOMA_REPORT_FILE", DATA_DIR / "report_file.tsv"))
POOLS_DIR = Path(os.environ.get("EUKARYOMA_POOLS_DIR", DATA_DIR / "pools"))

# Per-pair AlphaFold3 interaction scores (one row per pair per sample; see
# eukaryoma_ppi.index.read_recap_scores for how these are aggregated).
RECAP_FILE = Path(os.environ.get("EUKARYOMA_RECAP_FILE", DATA_DIR / "recap_set0_pairs.tsv"))

# Protein annotations (Pfam domains / descriptions) parsed from FASTA headers.
# See eukaryoma_ppi.annotations for the id transformation this file's header
# ids need to match the website's protein ids.
FASTA_FILE = Path(
    os.environ.get(
        "EUKARYOMA_FASTA_FILE",
        DATA_DIR / "OMAfiltered_MSfiltered_MICH_Capsaspora_owczarzaki_Schultz_A_annotated.fa",
    )
)

# Optional external association-score sources. Each is a dense protein x
# protein correlation matrix CSV; ids need normalizing to the website's id
# format (see eukaryoma_ppi.external_scores). Phyloprofiling has no file yet
# -- every code path treats a missing file as "source not available".
OTHER_DATA_SOURCES_DIR = Path(os.environ.get("EUKARYOMA_OTHER_DATA_SOURCES_DIR", DATA_DIR / "other_data_sources"))
COABUNDANCE_FILE = Path(
    os.environ.get(
        "EUKARYOMA_COABUNDANCE_FILE", OTHER_DATA_SOURCES_DIR / "coabundance" / "latest_coabundance_matrix.csv"
    )
)
COFRACTIONATION_FILE = Path(
    os.environ.get(
        "EUKARYOMA_COFRACTIONATION_FILE", OTHER_DATA_SOURCES_DIR / "cofractionation" / "latest_cofrac_matrix.csv"
    )
)
PHYLOPROFILING_FILE = Path(
    os.environ.get(
        "EUKARYOMA_PHYLOPROFILING_FILE",
        OTHER_DATA_SOURCES_DIR / "phyloprofiling" / "latest_phyloprofiling_matrix.csv",
    )
)

# "True positive" interaction annotations: two proteins co-occurring in the
# same complex/category are treated as a known-true interaction. See
# eukaryoma_ppi.complex_annotations for how these files are parsed.
CORUM_FILE = Path(
    os.environ.get(
        "EUKARYOMA_CORUM_FILE", OTHER_DATA_SOURCES_DIR / "annotations" / "corum" / "corum_annotations.tsv"
    )
)
MARCOTTE_FILE = Path(
    os.environ.get(
        "EUKARYOMA_MARCOTTE_FILE", OTHER_DATA_SOURCES_DIR / "annotations" / "marcotte" / "marcotte_annotations.txt"
    )
)

# Derived data, produced once by scripts/build_data.py and only ever *read*
# by the Streamlit app (the app never needs foldcomp at runtime).
STRUCTURES_DIR = Path(os.environ.get("EUKARYOMA_STRUCTURES_DIR", DATA_DIR / "structures"))
INDEX_DIR = Path(os.environ.get("EUKARYOMA_INDEX_DIR", DATA_DIR / "index"))
PAIRS_INDEX_FILE = INDEX_DIR / "pairs.parquet"
POOLS_INDEX_FILE = INDEX_DIR / "pools.parquet"
ANNOTATIONS_INDEX_FILE = INDEX_DIR / "protein_annotations.parquet"
# One row per possible pair among the website's proteins (~2.3M for 2145
# proteins), with a score column per source (NaN where that source doesn't
# cover the pair) plus a combined unified_score. See eukaryoma_ppi.external_scores.
UNIVERSE_INDEX_FILE = INDEX_DIR / "universe_scores.parquet"
# Sparse: only pairs flagged true-positive by CORUM and/or Marcotte complex
# co-membership. See eukaryoma_ppi.complex_annotations.
TRUE_POSITIVE_INDEX_FILE = INDEX_DIR / "true_positive_pairs.parquet"

# Only needed for scripts/build_data.py (local/HPC preprocessing), never at
# Streamlit runtime. This is a custom internal build of foldcomp that reads
# the lab's "FCZC" complex format -- the public PyPI/GitHub foldcomp release
# cannot decompress these files.
FOLDCOMP_BIN = Path(os.environ.get("FOLDCOMP_BIN", REPO_ROOT / "bin" / "foldcomp"))
