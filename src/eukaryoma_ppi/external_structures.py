"""Fetch reference monomer structures from AlphaFold DB by UniProt
accession, for comparing a locally-predicted AF3 pool pair against each
protein's own solo AlphaFold DB prediction.

Network calls -- best-effort: any failure (not in AFDB, network error,
timeout) returns None rather than raising, since this is a nice-to-have
alongside a pair's own (locally predicted) structure, not data the rest of
the page depends on.
"""

import requests

AFDB_TIMEOUT_SECONDS = 10
AFDB_API_URL = "https://alphafold.ebi.ac.uk/api/prediction/{accession}"


def fetch_alphafold_pdb(accession):
    """Best-effort fetch of a UniProt accession's AlphaFold DB structure as
    PDB text, or None if it's unavailable for any reason.

    Goes through AFDB's prediction API rather than guessing the file URL's
    model-version suffix directly (AF-<accession>-F1-model_v<N>.pdb, where N
    has changed release to release -- v4 files 404 now that AFDB is on v6)
    -- the API returns the current pdbUrl for whatever version is live.
    """
    try:
        api_response = requests.get(AFDB_API_URL.format(accession=accession), timeout=AFDB_TIMEOUT_SECONDS)
        api_response.raise_for_status()
        predictions = api_response.json()
        if not predictions:
            return None
        pdb_url = predictions[0]["pdbUrl"]

        pdb_response = requests.get(pdb_url, timeout=AFDB_TIMEOUT_SECONDS)
        pdb_response.raise_for_status()
        return pdb_response.text
    except (requests.RequestException, KeyError, ValueError, IndexError):
        return None
