#!/usr/bin/env python
"""One-time data conversion for the Eukaryoma PPI webserver.

Run this locally (never inside Docker) on a machine with the working
bin/foldcomp binary available:

    source venv/bin/activate
    python scripts/build_data.py

It:
  1. Decompresses every pool in data/pools/*.fcz to data/structures/<pool>.cif
     (skips pools already decompressed, unless --force is passed).
  2. Cross-references data/report_file.tsv against the pools that actually
     have a structure, and writes data/index/pools.parquet and
     data/index/pairs.parquet -- the pair-lookup table the Streamlit app
     queries at runtime.
  3. Parses protein annotations from the FASTA headers into
     data/index/protein_annotations.parquet, and reports whether every
     website protein id has a matching FASTA entry.
  3b. Parses eggNOG-mapper functional annotation (GO/KEGG/PFAM/ortholog
     ids) into data/index/eggnog_annotations.parquet.
  4. Parses CORUM/Marcotte complex co-membership into
     data/index/true_positive_pairs.parquet (corum_tp/marcotte_tp flags,
     merged into both pairs.parquet and universe_scores.parquet).
  5. Builds data/index/universe_scores.parquet -- every possible pair among
     the website's proteins, with a score column per optional external
     source (coabundance/cofractionation/phyloprofiling) when its file is
     present, plus the AF3 pool iptm and a combined unified_score.
  6. Resolves Capsaspora-to-yeast orthologs (network-bound, via OMA's REST
     API -- cached, so this only actually runs once) into
     data/index/yeast_orthologs.parquet, then decompresses just the
     referenced yeast AF3 monomer predictions into
     data/structures/yeast/<accession>.pdb for Protein View's reference-
     structure picker.

The Streamlit app only ever reads data/structures/ and data/index/; it does
not depend on foldcomp or the raw .fcz files.
"""

import argparse
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from eukaryoma_ppi import annotations, complex_annotations, eggnog, external_scores, index, structures, yeast_orthologs
from eukaryoma_ppi.config import (
    EGGNOG_FILE,
    FASTA_FILE,
    FOLDCOMP_BIN,
    POOLS_DIR,
    REPORT_FILE,
    STRUCTURES_DIR,
    YEAST_HOG2GENE_FILE,
    YEAST_ORTHOLOGS_INDEX_FILE,
    YEAST_STRUCTURES_DIR,
)
from eukaryoma_ppi.foldcomp_cli import FoldcompError, decompress_pool


def decompress_all_pools(force=False):
    STRUCTURES_DIR.mkdir(parents=True, exist_ok=True)
    fcz_files = sorted(POOLS_DIR.glob("*.fcz"))
    print(f"Found {len(fcz_files)} pool .fcz files in {POOLS_DIR}")

    available = set()
    failed = []
    start = time.time()
    for i, fcz_path in enumerate(fcz_files, 1):
        pool = fcz_path.stem
        out_path = STRUCTURES_DIR / f"{pool}.cif"
        if out_path.exists() and not force:
            available.add(pool)
            continue
        try:
            decompress_pool(fcz_path, out_path)
            available.add(pool)
        except FoldcompError as e:
            failed.append(pool)
            print(f"  [FAIL] {pool}: {e}", file=sys.stderr)
        if i % 250 == 0 or i == len(fcz_files):
            elapsed = time.time() - start
            print(f"  ...{i}/{len(fcz_files)} decompressed ({elapsed:.1f}s elapsed)")

    if failed:
        print(f"{len(failed)} pool(s) failed to decompress: {failed[:10]}{'...' if len(failed) > 10 else ''}")
    print(f"{len(available)} pool structures available in {STRUCTURES_DIR}")
    return available


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--force", action="store_true", help="Re-decompress pools even if a .cif already exists")
    parser.add_argument(
        "--skip-decompress",
        action="store_true",
        help="Skip decompression, only (re)build the index from data/structures/ already on disk",
    )
    parser.add_argument(
        "--skip-yeast-orthologs",
        action="store_true",
        help="Skip resolving/decompressing yeast reference structures (step 6 -- network-bound)",
    )
    parser.add_argument(
        "--force-yeast-orthologs",
        action="store_true",
        help="Re-resolve yeast orthologs via OMA's API even if yeast_orthologs.parquet already exists",
    )
    args = parser.parse_args()

    if not FOLDCOMP_BIN.exists() and not args.skip_decompress:
        sys.exit(
            f"foldcomp binary not found at {FOLDCOMP_BIN}. "
            "Set FOLDCOMP_BIN to point at a build that can read this lab's FCZ format."
        )
    if not REPORT_FILE.exists():
        sys.exit(f"report_file.tsv not found at {REPORT_FILE}. Set EUKARYOMA_DATA_DIR/EUKARYOMA_REPORT_FILE.")

    if args.skip_decompress:
        available = {p.stem for p in STRUCTURES_DIR.glob("*.cif")}
        print(f"Skipping decompression, found {len(available)} existing structures in {STRUCTURES_DIR}")
    else:
        available = decompress_all_pools(force=args.force)

    print("Building pool/pair index...")
    pools_df, pairs_df = index.build_index(available)
    print(f"Indexed {len(pools_df)} pools, {len(pairs_df)} protein pair occurrences.")

    if FASTA_FILE.exists():
        print("Parsing protein annotations from FASTA headers...")
        annotations_df = annotations.parse_fasta_annotations()
        annotations.save_annotations(annotations_df)
        website_ids = set(pairs_df["protein_a"]) | set(pairs_df["protein_b"])
        missing = annotations.missing_protein_ids(website_ids, annotations_df)
        if missing:
            print(
                f"  [WARN] {len(missing)}/{len(website_ids)} website protein ids have no FASTA annotation: "
                f"{sorted(missing)[:10]}{'...' if len(missing) > 10 else ''}"
            )
        else:
            print(f"  All {len(website_ids)} website protein ids matched a FASTA header ({len(annotations_df)} total).")
    else:
        print(f"  [WARN] FASTA_FILE not found at {FASTA_FILE}; skipping protein annotations.")

    if EGGNOG_FILE.exists():
        print("Parsing eggNOG-mapper functional annotation...")
        eggnog_df = eggnog.parse_eggnog_annotations()
        eggnog.save_eggnog_annotations(eggnog_df)
        website_ids = set(pairs_df["protein_a"]) | set(pairs_df["protein_b"])
        missing = eggnog.missing_protein_ids(website_ids, eggnog_df)
        if missing:
            print(
                f"  [WARN] {len(missing)}/{len(website_ids)} website protein ids have no eggNOG annotation: "
                f"{sorted(missing)[:10]}{'...' if len(missing) > 10 else ''}"
            )
        else:
            print(f"  All {len(website_ids)} website protein ids matched an eggNOG entry ({len(eggnog_df)} total).")
    else:
        print(f"  [WARN] EGGNOG_FILE not found at {EGGNOG_FILE}; skipping eggNOG annotations.")

    print("Parsing CORUM/Marcotte true-positive complex co-membership...")
    website_ids = set(pairs_df["protein_a"]) | set(pairs_df["protein_b"])
    for flag_col, (path, _group_col) in complex_annotations.SOURCES.items():
        print(f"  {flag_col}: {'found ' + str(path) if path.exists() else 'NOT FOUND, skipping'}")
    tp_flags_df = complex_annotations.build_true_positive_flags(website_ids)
    complex_annotations.save_true_positive_flags(tp_flags_df)
    for flag_col in complex_annotations.FLAG_COLUMNS:
        print(f"  {flag_col}: {int(tp_flags_df[flag_col].sum()):,} pairs")
    n_either = len(tp_flags_df[tp_flags_df[complex_annotations.FLAG_COLUMNS].any(axis=1)])
    n_both = len(tp_flags_df[tp_flags_df[complex_annotations.FLAG_COLUMNS].all(axis=1)])
    print(f"  either: {n_either:,} pairs, both: {n_both:,} pairs")

    pairs_df = complex_annotations.attach_true_positive_flags(pairs_df)
    index.save_index(pools_df, pairs_df)
    print("Saved data/index/pairs.parquet with true-positive flags attached.")

    print("Building the full pair universe (AF3 + external scores)...")
    website_ids = set(pairs_df["protein_a"]) | set(pairs_df["protein_b"])
    for source_name, (path, _col, _sym) in external_scores.SOURCES.items():
        print(f"  {source_name}: {'found ' + str(path) if path.exists() else 'NOT FOUND, skipping'}")
    universe_df = external_scores.build_universe(website_ids, pairs_df)
    universe_df = complex_annotations.attach_true_positive_flags(universe_df)
    external_scores.save_universe(universe_df)
    print(f"  Universe has {len(universe_df):,} possible pairs among {len(website_ids)} proteins.")

    pattern_counts, n_sources_counts = external_scores.source_presence_summary(universe_df)
    print("  Pairs by number of sources present:")
    for n, count in n_sources_counts.items():
        print(f"    {n} source(s): {count:,}")
    print("  Pairs by exact combination of sources present:")
    for _, row in pattern_counts.iterrows():
        print(f"    {row['sources_present']}: {row['n_pairs']:,}")

    if args.skip_yeast_orthologs:
        print("Skipping yeast ortholog resolution (--skip-yeast-orthologs).")
    elif not YEAST_HOG2GENE_FILE.exists() or not YEAST_STRUCTURES_DIR.exists():
        print(
            f"  [WARN] Yeast HOG map ({YEAST_HOG2GENE_FILE}) or structures dir ({YEAST_STRUCTURES_DIR}) not "
            "found; skipping yeast reference structures."
        )
    elif YEAST_ORTHOLOGS_INDEX_FILE.exists() and not args.force_yeast_orthologs:
        print(
            f"Yeast orthologs already resolved at {YEAST_ORTHOLOGS_INDEX_FILE} "
            "(pass --force-yeast-orthologs to re-resolve)."
        )
        yeast_df = yeast_orthologs.load_yeast_orthologs()
        decompressed, missing = structures.decompress_yeast_structures(set(yeast_df["yeast_accession"]))
        print(f"  {decompressed:,} yeast reference structure(s) available, {missing:,} missing a local .fcz.")
    else:
        print("Resolving Capsaspora-to-yeast orthologs via OMA's REST API (network-bound, a few minutes)...")

        def _progress(i, n):
            if i % 200 == 0 or i == n:
                print(f"    ...resolved {i}/{n} yeast OMA ids")

        yeast_df = yeast_orthologs.build_yeast_orthologs(website_ids, progress_callback=_progress)
        yeast_orthologs.save_yeast_orthologs(yeast_df)
        print(
            f"  {len(yeast_df):,} (protein, yeast ortholog) pairs resolved, covering "
            f"{yeast_df['protein_id'].nunique():,}/{len(website_ids):,} website proteins."
        )
        decompressed, missing = structures.decompress_yeast_structures(set(yeast_df["yeast_accession"]))
        print(f"  {decompressed:,} yeast reference structure(s) decompressed, {missing:,} missing a local .fcz.")


if __name__ == "__main__":
    main()
