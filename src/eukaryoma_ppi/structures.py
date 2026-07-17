"""Extract a protein pair's substructure from a cached, decompressed pool CIF.

The Streamlit app never touches foldcomp or the raw .fcz files: by the time
it runs, scripts/build_data.py has already decompressed every pool once into
data/structures/<pool>.cif. This module just pulls the two relevant chains
back out of that plain mmCIF text with Biopython.
"""

import io

from Bio.PDB import MMCIFParser, PDBIO, Select

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
