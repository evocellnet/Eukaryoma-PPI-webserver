"""Protein annotations parsed from FASTA headers.

Each header looks like:

    >XP_004340666.2,4-aminobutyrate:2-oxoglutarate transaminase activity,Aminotran_3

The first comma-separated field is the protein id in NCBI format
(XP_004340666.2); the website uses a different id format
(xp0043406662) derived from it by dropping the "." and "_" and
lowercasing everything. The remaining comma-separated fields are
free-text/Pfam annotations for that protein.
"""

import pandas as pd

from eukaryoma_ppi.config import ANNOTATIONS_INDEX_FILE, FASTA_FILE


def fasta_id_to_protein_id(raw_id):
    """XP_004340666.2 -> xp0043406662."""
    return raw_id.replace(".", "").replace("_", "").lower()


def parse_fasta_annotations():
    """Parse FASTA_FILE headers into one row per protein."""
    records = []
    with open(FASTA_FILE) as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            fields = line[1:].rstrip("\n").split(",")
            raw_id, annotation_fields = fields[0], [f for f in fields[1:] if f]
            records.append(
                {
                    "protein_id": fasta_id_to_protein_id(raw_id),
                    "fasta_id": raw_id,
                    "annotation": "; ".join(annotation_fields),
                }
            )
    return pd.DataFrame.from_records(records)


def save_annotations(annotations_df):
    ANNOTATIONS_INDEX_FILE.parent.mkdir(parents=True, exist_ok=True)
    annotations_df.to_parquet(ANNOTATIONS_INDEX_FILE, index=False)


def load_annotations():
    return pd.read_parquet(ANNOTATIONS_INDEX_FILE)


def missing_protein_ids(protein_ids, annotations_df):
    """Website protein ids with no matching FASTA header, if any."""
    return set(protein_ids) - set(annotations_df["protein_id"])
