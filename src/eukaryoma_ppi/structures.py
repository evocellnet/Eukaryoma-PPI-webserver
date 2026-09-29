"""Extract a protein pair's substructure from a cached, decompressed pool CIF.

The Streamlit app never touches foldcomp or the raw .fcz files: by the time
it runs, scripts/build_data.py has already decompressed every pool once into
data/structures/<pool>.cif (and, for the yeast reference structures used by
Protein View, every accession an ortholog mapping actually references into
data/structures/yeast/<accession>.pdb). This module just reads that already-
decompressed data back -- pulling the two relevant chains out of a pool's
mmCIF with Biopython, or a yeast monomer's PDB text as-is.
"""

import io
import tempfile
from pathlib import Path

from Bio.PDB import MMCIFParser, PDBIO, Select, Superimposer

from eukaryoma_ppi.config import STRUCTURES_DIR, YEAST_STRUCTURES_DECOMPRESSED_DIR, YEAST_STRUCTURES_DIR
from eukaryoma_ppi.foldcomp_cli import FoldcompError, decompress_to_file

_parser = MMCIFParser(QUIET=True)


class _ChainSelect(Select):
    def __init__(self, chain_ids):
        self.chain_ids = set(chain_ids)

    def accept_chain(self, chain):
        return chain.id in self.chain_ids


def pool_cif_path(pool):
    return STRUCTURES_DIR / f"{pool}.cif"


def extract_chains_as_pdb(pool, chain_ids):
    """Return a PDB-format string containing only the requested chains."""
    cif_path = pool_cif_path(pool)
    structure = _parser.get_structure(pool, str(cif_path))

    buf = io.StringIO()
    io_writer = PDBIO()
    io_writer.set_structure(structure)
    io_writer.save(buf, select=_ChainSelect(chain_ids))
    return buf.getvalue()


def _load_chain(pool, chain_id):
    cif_path = pool_cif_path(pool)
    structure = _parser.get_structure(pool, str(cif_path))
    return structure[0][chain_id]


def _ca_atoms(chain):
    return [residue["CA"] for residue in chain if "CA" in residue]


def superpose_chains_as_pdb(pool_chain_pairs):
    """One protein, predicted in several different pools (each time with a
    different partner) -- superpose just that protein's own chain from each
    pool onto the first one's coordinate frame, for comparing how
    consistently AF3 predicts its fold across pooling contexts.

    pool_chain_pairs: [(pool, chain_id), ...], at least one entry. Returns a
    list of PDB-format strings the same length and order, already
    superposed (the first is left as-is, it's the reference); an entry is
    None instead if that pool's chain has a different CA count than the
    reference's (e.g. a different construct boundary) -- Superimposer needs
    a 1:1 atom correspondence, so it can't be aligned.
    """
    chains = [_load_chain(pool, chain_id) for pool, chain_id in pool_chain_pairs]
    reference_atoms = _ca_atoms(chains[0])

    pdb_texts = []
    for chain in chains:
        atoms = _ca_atoms(chain)
        if len(atoms) != len(reference_atoms):
            pdb_texts.append(None)
            continue
        if chain is not chains[0]:
            superimposer = Superimposer()
            superimposer.set_atoms(reference_atoms, atoms)
            superimposer.apply(list(chain.get_atoms()))

        buf = io.StringIO()
        io_writer = PDBIO()
        io_writer.set_structure(chain)
        io_writer.save(buf)
        pdb_texts.append(buf.getvalue())
    return pdb_texts


def yeast_structure_pdb_path(accession):
    return YEAST_STRUCTURES_DECOMPRESSED_DIR / f"{accession.lower()}.pdb"


def get_yeast_structure_pdb(accession):
    """Pre-decompressed yeast AF3 monomer prediction, as PDB text -- None if
    scripts/build_data.py hasn't decompressed this accession (either no
    local .fcz for it, or it isn't referenced by any Capsaspora ortholog).
    Pure file read, no foldcomp involved at Streamlit runtime.
    """
    path = yeast_structure_pdb_path(accession)
    if not path.exists():
        return None
    return path.read_text()


def decompress_yeast_structures(accessions, force=False):
    """Build-time only: decompress each accession's compressed AF3 monomer
    prediction (YEAST_STRUCTURES_DIR/<accession>.fcz) to a plain .pdb file
    in YEAST_STRUCTURES_DECOMPRESSED_DIR, for get_yeast_structure_pdb to
    read later. Only decompresses the given accessions, not the full
    ~6,000-protein yeast set -- see config.YEAST_STRUCTURES_DECOMPRESSED_DIR
    for why (most of it would never be looked up).

    Goes through mmCIF and Biopython's own PDBIO rather than decompressing
    straight to PDB (like the pools do): the foldcomp CLI's own PDB writer
    misformats negative coordinates for these particular monomer files
    (verified directly -- e.g. resSeq/coordinate columns drift out of their
    fixed widths, corrupting every atom after the first negative x), while
    its mmCIF output is whitespace-delimited and unaffected, and Biopython's
    writer (already used for every pool structure) is known-correct.

    Returns (decompressed, missing): counts of accessions with/without a
    local .fcz file to decompress from.
    """
    YEAST_STRUCTURES_DECOMPRESSED_DIR.mkdir(parents=True, exist_ok=True)
    decompressed, missing = 0, 0
    for accession in accessions:
        out_path = yeast_structure_pdb_path(accession)
        if out_path.exists() and not force:
            decompressed += 1
            continue
        fcz_path = YEAST_STRUCTURES_DIR / f"{accession.lower()}.fcz"
        if not fcz_path.exists():
            missing += 1
            continue
        try:
            with tempfile.TemporaryDirectory() as tmp_dir:
                cif_path = Path(tmp_dir) / f"{accession}.cif"
                decompress_to_file(fcz_path, cif_path)
                structure = _parser.get_structure(accession, str(cif_path))
                io_writer = PDBIO()
                io_writer.set_structure(structure)
                io_writer.save(str(out_path))
            decompressed += 1
        except FoldcompError:
            missing += 1
    return decompressed, missing
