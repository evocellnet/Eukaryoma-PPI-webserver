"""Build and query the pool/pair index derived from report_file.tsv.

report_file.tsv has one row per pooled AlphaFold3 prediction:
    <pool_name>\t<protein_1>_<protein_2>_..._<protein_n>

The protein order matches the chain order (A, B, C, ...) in the pool's
mmCIF structure. Only pools that actually have a decompressed structure
(data/structures/<pool>.cif) are indexed.
"""

import re
import warnings
from itertools import combinations

import pandas as pd

from eukaryoma_ppi.config import PAIRS_INDEX_FILE, POOLS_INDEX_FILE, RECAP_FILE, REPORT_FILE

# Matches the two quoted protein ids in a recap "ids" cell, e.g.
# "('xp0043429311', 'xp0043454971')".
RECAP_IDS_RE = re.compile(r"'([^']+)',\s*'([^']+)'")

def chain_letter(position):
    """0-indexed chain position -> mmCIF chain id (A, B, ..., Z, AA, AB, ...)."""
    letters = ""
    n = position
    while True:
        n, rem = divmod(n, 26)
        letters = chr(ord("A") + rem) + letters
        if n == 0:
            return letters
        n -= 1


def read_report(available_pools):
    """Read report_file.tsv, keeping only rows with a structure on disk.

    available_pools: set of pool names (no extension) that have a
    corresponding .cif in data/structures/.
    """
    rows = []
    with open(REPORT_FILE) as fh:
        for line in fh:
            pool, proteins = line.rstrip("\n").split("\t")
            if pool in available_pools:
                rows.append((pool, proteins.split("_")))
    return rows


def build_index(available_pools):
    """Build the pools and pairs index DataFrames from report_file.tsv."""
    rows = read_report(available_pools)

    pool_records = []
    pair_records = []
    for pool, proteins in rows:
        pool_records.append({"pool": pool, "n_proteins": len(proteins), "proteins": "_".join(proteins)})
        for (i, protein_a), (j, protein_b) in combinations(enumerate(proteins), 2):
            pair_records.append(
                {
                    "protein_a": protein_a,
                    "protein_b": protein_b,
                    "pool": pool,
                    "chain_a": chain_letter(i),
                    "chain_b": chain_letter(j),
                }
            )

    pools_df = pd.DataFrame.from_records(pool_records)
    pairs_df = pd.DataFrame.from_records(pair_records)
    pairs_df = attach_scores(pairs_df)
    return pools_df, pairs_df


def read_recap_scores():
    """Load per-pair AlphaFold3 interaction scores from RECAP_FILE.

    recap_set0_pairs.tsv has one row per (pair, sample) -- AF3 predicts each
    pool several times (multiple seeds/samples) and this file records every
    sample's score, so a given (pool, pair) shows up several times with
    slightly different corrected_chain_pair_iptm values. We average those
    into a single score per (pool, pair).
    """
    df = pd.read_csv(RECAP_FILE, sep="\t", usecols=["ids", "name", "corrected_chain_pair_iptm"])
    parsed = df["ids"].str.extract(RECAP_IDS_RE)

    lo = parsed[0].where(parsed[0] <= parsed[1], parsed[1])
    hi = parsed[0].where(parsed[0] > parsed[1], parsed[1])
    df["protein_lo"], df["protein_hi"] = lo, hi
    df = df.rename(columns={"name": "pool"})

    return (
        df.groupby(["pool", "protein_lo", "protein_hi"])["corrected_chain_pair_iptm"]
        .mean()
        .reset_index()
    )


def attach_scores(pairs_df):
    """Left-join the averaged AF3 interaction score onto the pairs index."""
    if not RECAP_FILE.exists():
        warnings.warn(f"RECAP_FILE not found at {RECAP_FILE}; pairs will have no iptm score.")
        pairs_df["corrected_chain_pair_iptm"] = float("nan")
        return pairs_df

    scores = read_recap_scores()
    protein_lo = pairs_df[["protein_a", "protein_b"]].min(axis=1)
    protein_hi = pairs_df[["protein_a", "protein_b"]].max(axis=1)
    merged = pairs_df.assign(protein_lo=protein_lo, protein_hi=protein_hi).merge(
        scores, on=["pool", "protein_lo", "protein_hi"], how="left"
    )
    return merged.drop(columns=["protein_lo", "protein_hi"])


def collapse_best_per_pair(pairs_df):
    """One row per (protein_a, protein_b) pair, keeping the pool occurrence
    with the highest corrected_chain_pair_iptm (a small number of pairs
    co-occur in more than one pool). Also normalizes protein_a/protein_b to
    (min, max) order, matching the convention used across all pair tables.
    Used to fold the per-pool AF3 pairs index into the universe table, which
    has exactly one row per pair.
    """
    df = pairs_df.copy()
    df["protein_a"], df["protein_b"] = (
        df[["protein_a", "protein_b"]].min(axis=1),
        df[["protein_a", "protein_b"]].max(axis=1),
    )
    df = df.sort_values("corrected_chain_pair_iptm", ascending=False)
    return df.drop_duplicates(subset=["protein_a", "protein_b"], keep="first").reset_index(drop=True)


def save_index(pools_df, pairs_df):
    POOLS_INDEX_FILE.parent.mkdir(parents=True, exist_ok=True)
    pools_df.to_parquet(POOLS_INDEX_FILE, index=False)
    pairs_df.to_parquet(PAIRS_INDEX_FILE, index=False)


def load_pairs_index():
    return pd.read_parquet(PAIRS_INDEX_FILE)


def load_pools_index():
    return pd.read_parquet(POOLS_INDEX_FILE)


def lookup_pair(pairs_df, protein_a, protein_b):
    """Return all pool occurrences of an unordered protein pair."""
    mask = ((pairs_df["protein_a"] == protein_a) & (pairs_df["protein_b"] == protein_b)) | (
        (pairs_df["protein_a"] == protein_b) & (pairs_df["protein_b"] == protein_a)
    )
    return pairs_df.loc[mask]
