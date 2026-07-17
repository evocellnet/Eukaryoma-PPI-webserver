"""Build and query the pool/pair index derived from report_file.tsv.

report_file.tsv has one row per pooled AlphaFold3 prediction:
    <pool_name>\t<protein_1>_<protein_2>_..._<protein_n>

The protein order matches the chain order (A, B, C, ...) in the pool's
mmCIF structure. Only pools that actually have a decompressed structure
(data/structures/<pool>.cif) are indexed.
"""

from itertools import combinations

import pandas as pd

from eukaryoma_ppi.config import PAIRS_INDEX_FILE, POOLS_INDEX_FILE, REPORT_FILE

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
    return pools_df, pairs_df


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
