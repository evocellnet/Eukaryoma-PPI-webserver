"""Extract a protein pair's substructure from a cached, decompressed pool CIF.

The Streamlit app never touches foldcomp or the raw .fcz files: by the time
it runs, scripts/build_data.py has already decompressed every pool once into
data/structures/<pool>.cif. This module just pulls the two relevant chains
back out of that plain mmCIF text with Biopython.
"""

import io

from Bio.PDB import MMCIFParser, PDBIO, Select, Superimposer

from eukaryoma_ppi.config import STRUCTURES_DIR

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
