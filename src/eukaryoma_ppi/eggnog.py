"""eggNOG-mapper functional annotation: GO terms, KEGG, PFAM domains, and an
ortholog identifier usable for external database links (NCBI/UniProt/
AlphaFold DB). Distinct from the OMA/FASTA annotation in
eukaryoma_ppi.annotations -- always label this "eggNOG" in the UI so the two
aren't confused.

The file is standard emapper.py TSV output: a few "##" comment/header lines,
one "#query\t..." column-header line, then one row per protein. E.g.:

    #query  seed_ortholog  evalue  score  eggNOG_OGs  ...  GOs  ...  PFAMs
    XP_004339898.1  192875.XP_004339898.1  6.89e-195  541.0  ...

"seed_ortholog" is "<ortholog taxon id>.<accession>" -- for ~99% of proteins
here the ortholog is the protein's own existing NCBI RefSeq entry (self-hit,
since Capsaspora is already in eggNOG's reference databases), for the rest
it's a genuine UniProt accession from a different species. We tell the two
apart by whether the accession contains "_": RefSeq accessions always do
(XP_.../NP_...), UniProt accessions never do.
"""

import pandas as pd

from eukaryoma_ppi.config import EGGNOG_FILE, EGGNOG_INDEX_FILE
from eukaryoma_ppi.external_scores import normalize_external_id

MULTI_VALUE_COLUMNS = ["GOs", "KEGG_ko", "KEGG_Pathway", "KEGG_Module", "KEGG_Reaction", "KEGG_rclass", "PFAMs"]
NA_TOKEN = "-"


def _split_multi(value):
    if value == NA_TOKEN or not value:
        return []
    return value.split(",")


def _ortholog_link(seed_ortholog):
    """"192875.XP_004339898.1" -> ("XP_004339898.1", "refseq")
    "1392244.V5ENR7" -> ("V5ENR7", "uniprot")
    """
    if seed_ortholog == NA_TOKEN or "." not in seed_ortholog:
        return None, None
    accession = seed_ortholog.split(".", 1)[1]
    kind = "refseq" if "_" in accession else "uniprot"
    return accession, kind


def parse_eggnog_annotations():
    """Parse EGGNOG_FILE into one row per protein (website id format)."""
    header = None
    rows = []
    with open(EGGNOG_FILE) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("##"):
                continue
            if line.startswith("#"):
                header = line[1:].split("\t")
                continue
            if line:
                rows.append(line.split("\t"))

    df = pd.DataFrame(rows, columns=header)
    df["protein_id"] = df["query"].map(normalize_external_id)

    for col in MULTI_VALUE_COLUMNS:
        df[col] = df[col].map(_split_multi)

    ortholog = df["seed_ortholog"].map(_ortholog_link)
    df["ortholog_accession"] = ortholog.map(lambda t: t[0])
    df["ortholog_kind"] = ortholog.map(lambda t: t[1])

    df["Description"] = df["Description"].replace(NA_TOKEN, "")
    df["Preferred_name"] = df["Preferred_name"].replace(NA_TOKEN, "")

    return df


def save_eggnog_annotations(df):
    EGGNOG_INDEX_FILE.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(EGGNOG_INDEX_FILE, index=False)


def load_eggnog_annotations():
    return pd.read_parquet(EGGNOG_INDEX_FILE)


def as_list(value):
    """Multi-value columns (GOs, KEGG_ko, ...) come back from parquet as
    numpy arrays, not Python lists -- `array or []` raises "ambiguous truth
    value", so callers should use this instead of a plain `or []` fallback.
    """
    if value is None:
        return []
    return list(value)


def missing_protein_ids(protein_ids, eggnog_df):
    return set(protein_ids) - set(eggnog_df["protein_id"])


def ncbi_protein_url(accession):
    return f"https://www.ncbi.nlm.nih.gov/protein/{accession}"


def uniprot_url(accession, kind):
    if kind == "uniprot":
        return f"https://www.uniprot.org/uniprotkb/{accession}/entry", True
    return f"https://www.uniprot.org/uniprotkb?query={accession}", False


def alphafold_url(accession, kind):
    if kind == "uniprot":
        return f"https://alphafold.ebi.ac.uk/entry/{accession}", True
    return f"https://alphafold.ebi.ac.uk/search/text/{accession}", False


def quickgo_url(go_id):
    return f"https://www.ebi.ac.uk/QuickGO/term/{go_id}"


def kegg_url(kegg_id):
    return f"https://www.kegg.jp/entry/{kegg_id.split(':')[-1]}"
