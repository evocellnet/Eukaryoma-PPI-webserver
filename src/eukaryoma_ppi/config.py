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

# Derived data, produced once by scripts/build_data.py and only ever *read*
# by the Streamlit app (the app never needs foldcomp at runtime).
STRUCTURES_DIR = Path(os.environ.get("EUKARYOMA_STRUCTURES_DIR", DATA_DIR / "structures"))
INDEX_DIR = Path(os.environ.get("EUKARYOMA_INDEX_DIR", DATA_DIR / "index"))
PAIRS_INDEX_FILE = INDEX_DIR / "pairs.parquet"
POOLS_INDEX_FILE = INDEX_DIR / "pools.parquet"

# Only needed for scripts/build_data.py (local/HPC preprocessing), never at
# Streamlit runtime. This is a custom internal build of foldcomp that reads
# the lab's "FCZC" complex format -- the public PyPI/GitHub foldcomp release
# cannot decompress these files.
FOLDCOMP_BIN = Path(os.environ.get("FOLDCOMP_BIN", REPO_ROOT / "bin" / "foldcomp"))
