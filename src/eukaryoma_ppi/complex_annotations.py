""""True positive" interaction annotations from CORUM and Marcotte complex
membership: two proteins that co-occur in the same complex/category are
treated as a known-true interaction.

Both files list one row per (complex, human ortholog), with a "Capsaspora"
column giving that ortholog's Capsaspora protein id(s) -- ";"-joined when a
human gene has several Capsaspora paralogs, "N/A" when it has none. E.g.:

    corumID  category                        ...  Capsaspora
    4        Multisubunit ACTR coactivator...  XP_004345483.1_595528;XP_004364839.1_595528

Ids carry a trailing "_<taxon id>" the FASTA/website ids don't have, so they
need stripping before eukaryoma_ppi.external_scores.normalize_external_id
applies. Parsing approach follows
/Users/spascare/data/work/beltrao/combine_assoc_scores/Eukaryoma_PPI_analysis/src/cofrac_coab9partial_tmp_test.py,
except true-positive pairs are computed as the union of pairwise
combinations within *every* group a protein belongs to (that script instead
keeps only one category per protein, last-row-wins, which drops some
co-membership pairs when a protein belongs to more than one complex).

CORUM groups by "corumID" (one complex); Marcotte groups by "category" (a
coarse functional module, matching the reference script -- "category-specific"
would give finer-grained sub-groups if that's ever wanted instead).
"""

from itertools import combinations

import pandas as pd

from eukaryoma_ppi.config import CORUM_FILE, MARCOTTE_FILE, TRUE_POSITIVE_INDEX_FILE
from eukaryoma_ppi.external_scores import normalize_external_id

SOURCES = {
    "corum_tp": (CORUM_FILE, "corumID"),
    "marcotte_tp": (MARCOTTE_FILE, "category"),
}
FLAG_COLUMNS = list(SOURCES.keys())


def _member_website_ids(cell):
    """"XP_004345483.1_595528;XP_004364839.1_595528" -> ["xp0043454831", "xp0043648391"]."""
    if not cell or cell == "N/A":
        return []
    ids = []
    for raw in cell.split(";"):
        raw = raw.strip()
        if not raw or raw == "N/A":
            continue
        parts = raw.split("_")
        ncbi_id = "_".join(parts[:2]) if len(parts) >= 2 else raw
        ids.append(normalize_external_id(ncbi_id))
    return ids


def load_complex_groups(path, group_col, species_col="Capsaspora"):
    """{group_id: {protein_id, ...}} parsed from one CORUM/Marcotte-style file."""
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    groups = {}
    for _, row in df.iterrows():
        for protein_id in _member_website_ids(row[species_col]):
            groups.setdefault(row[group_col], set()).add(protein_id)
    return groups


def true_positive_pairs(path, group_col, website_ids):
    """{(protein_lo, protein_hi), ...}: every pair of website proteins sharing
    membership in some group. Annotation files cover the whole Capsaspora
    proteome, so members are restricted to website_ids -- otherwise this
    would also generate pairs for proteins outside the website's AF3-pooled
    protein set, which can never appear in any pair table here anyway.
    """
    groups = load_complex_groups(path, group_col)
    pairs = set()
    for members in groups.values():
        members = members & website_ids
        if len(members) < 2:
            continue
        pairs.update(combinations(sorted(members), 2))
    return pairs


def build_true_positive_flags(website_ids):
    """One row per pair flagged true-positive by CORUM and/or Marcotte, with
    a bool column per source (only True values are stored -- this table is
    sparse; callers left-join and fill missing with False).
    """
    frames = []
    for flag_col, (path, group_col) in SOURCES.items():
        if not path.exists():
            continue
        pairs = true_positive_pairs(path, group_col, website_ids)
        if not pairs:
            continue
        protein_a, protein_b = zip(*pairs)
        frames.append(pd.DataFrame({"protein_a": protein_a, "protein_b": protein_b, flag_col: True}))

    if not frames:
        return pd.DataFrame(columns=["protein_a", "protein_b", *FLAG_COLUMNS])

    merged = frames[0]
    for frame in frames[1:]:
        merged = merged.merge(frame, on=["protein_a", "protein_b"], how="outer")
    for flag_col in FLAG_COLUMNS:
        if flag_col not in merged.columns:
            merged[flag_col] = False
        merged[flag_col] = merged[flag_col].fillna(False)
    return merged


def save_true_positive_flags(df):
    TRUE_POSITIVE_INDEX_FILE.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(TRUE_POSITIVE_INDEX_FILE, index=False)


def load_true_positive_flags():
    return pd.read_parquet(TRUE_POSITIVE_INDEX_FILE)


def attach_true_positive_flags(df):
    """Left-join corum_tp/marcotte_tp onto any frame with protein_a/protein_b,
    filling non-matches with False (no complex evidence, not "unknown").
    """
    flags = load_true_positive_flags()
    merged = df.merge(flags, on=["protein_a", "protein_b"], how="left")
    for flag_col in FLAG_COLUMNS:
        merged[flag_col] = merged[flag_col].fillna(False)
    return merged
