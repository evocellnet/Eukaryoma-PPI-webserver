"""Capsaspora-to-yeast ortholog mapping, for Protein View's reference
structures -- but sourced from local AF3 monomer predictions instead of a
network fetch (see eukaryoma_ppi.structures for reading the decompressed
structures once scripts/build_data.py has produced them).

Source: a HOG (hierarchical orthologous group) -> gene map covering
Capsaspora plus five other species, pickled from an ancestral-genomes
OMA/FastOMA analysis (see config.YEAST_HOG2GENE_FILE). It's keyed by NCBI
taxon id; two genes sharing a HOG under Capsaspora's (192875) and yeast's
(4932) taxon keys are treated as orthologs. 1,375 HOGs contain both,
covering 1,896 of the website's 2,145 proteins (as of the current file) --
much better coverage than the human cross-reference in
eukaryoma_ppi.complex_annotations (~750 proteins), since this is a direct
Capsaspora<->yeast orthology call rather than routing through human
complexes.

The yeast side of the map uses OMA's own internal per-species ids (e.g.
"YEAST02308"), not the UniProt accessions the local yeast structure files
are named by, so each one needs resolving via OMA's REST API -- one-time,
at build time (a few minutes even with concurrency; never at Streamlit
runtime), cached to config.YEAST_ORTHOLOGS_INDEX_FILE.
"""

import pickle
import time
from concurrent.futures import ThreadPoolExecutor

import pandas as pd
import requests

from eukaryoma_ppi.config import YEAST_HOG2GENE_FILE, YEAST_ORTHOLOGS_INDEX_FILE
from eukaryoma_ppi.external_scores import normalize_external_id

CAPSASPORA_TAXON_ID = "192875"
YEAST_TAXON_ID = "4932"
OMA_PROTEIN_URL = "https://omabrowser.org/api/protein/{oma_id}/"
OMA_XREF_URL = "https://omabrowser.org/api/protein/{entry_nr}/xref/"


def load_hog_map():
    with open(YEAST_HOG2GENE_FILE, "rb") as fh:
        return pickle.load(fh)


def capsaspora_to_yeast_oma_ids(hog_map):
    """{protein_id: {yeast_oma_id, ...}}, one entry per Capsaspora gene that
    shares a HOG with at least one yeast gene. protein_id is normalized to
    the website's id format.
    """
    mapping = {}
    for species_dict in hog_map.values():
        cap_genes = species_dict.get(CAPSASPORA_TAXON_ID)
        yeast_genes = species_dict.get(YEAST_TAXON_ID)
        if not cap_genes or not yeast_genes:
            continue
        for raw in cap_genes:
            protein_id = normalize_external_id(raw)
            mapping.setdefault(protein_id, set()).update(yeast_genes)
    return mapping


def _resolve_one(oma_id, retries=3, backoff_seconds=1.5):
    """oma_id -> UniProt accession, or None. Retries a few times with
    backoff -- OMA's public server occasionally 502s under concurrent load
    (observed directly while building this), which is transient, not a
    real miss.
    """
    for attempt in range(retries):
        try:
            response = requests.get(OMA_PROTEIN_URL.format(oma_id=oma_id), timeout=20)
            if response.status_code == 200:
                entry_nr = response.json().get("entry_nr")
                xref_response = requests.get(OMA_XREF_URL.format(entry_nr=entry_nr), timeout=20)
                if xref_response.status_code == 200:
                    for item in xref_response.json():
                        xref = item.get("xref", "")
                        if item.get("source") == "UniProtKB/SwissProt" and xref.isalnum():
                            return oma_id, xref
                    return oma_id, None
        except requests.RequestException:
            pass
        time.sleep(backoff_seconds * (attempt + 1))
    return oma_id, None


def resolve_all_yeast_uniprot_accessions(oma_ids, max_workers=6, progress_callback=None):
    """{yeast_oma_id: uniprot_accession}, resolved via OMA's REST API.
    Network-bound and slow at ~1,800 ids (a few minutes even with
    concurrency) -- meant to run once at build time, cached via
    save_yeast_orthologs, never at Streamlit runtime.
    """
    oma_ids = list(oma_ids)
    results = {}
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        for i, (oma_id, accession) in enumerate(executor.map(_resolve_one, oma_ids), 1):
            if accession:
                results[oma_id] = accession
            if progress_callback:
                progress_callback(i, len(oma_ids))
    return results


def build_yeast_orthologs(website_ids, progress_callback=None):
    """[protein_id, yeast_accession] -- one row per (website protein,
    resolved yeast UniProt accession) pair.
    """
    hog_map = load_hog_map()
    cap_to_oma = capsaspora_to_yeast_oma_ids(hog_map)
    cap_to_oma = {protein_id: oma_ids for protein_id, oma_ids in cap_to_oma.items() if protein_id in website_ids}

    all_oma_ids = {oma_id for oma_ids in cap_to_oma.values() for oma_id in oma_ids}
    oma_to_uniprot = resolve_all_yeast_uniprot_accessions(all_oma_ids, progress_callback=progress_callback)

    rows = [
        {"protein_id": protein_id, "yeast_accession": oma_to_uniprot[oma_id]}
        for protein_id, oma_ids in cap_to_oma.items()
        for oma_id in oma_ids
        if oma_id in oma_to_uniprot
    ]
    return pd.DataFrame(rows, columns=["protein_id", "yeast_accession"]).drop_duplicates()


def save_yeast_orthologs(df):
    YEAST_ORTHOLOGS_INDEX_FILE.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(YEAST_ORTHOLOGS_INDEX_FILE, index=False)


def load_yeast_orthologs():
    return pd.read_parquet(YEAST_ORTHOLOGS_INDEX_FILE)


def yeast_orthologs_map():
    """{protein_id: {yeast_accession, ...}}, from the cached index."""
    df = load_yeast_orthologs()
    mapping = {}
    for protein_id, accession in zip(df["protein_id"], df["yeast_accession"]):
        mapping.setdefault(protein_id, set()).add(accession)
    return mapping
