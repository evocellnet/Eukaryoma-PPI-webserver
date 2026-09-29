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

import numpy as np
import pandas as pd
from Bio import Align
from Bio.Align import substitution_matrices
from Bio.PDB import MMCIFParser, PDBIO, PDBParser, Select, Superimposer
from Bio.SeqUtils import seq1

from eukaryoma_ppi.config import STRUCTURES_DIR, YEAST_STRUCTURES_DECOMPRESSED_DIR, YEAST_STRUCTURES_DIR
from eukaryoma_ppi.foldcomp_cli import FoldcompError, decompress_to_file

_parser = MMCIFParser(QUIET=True)
_pdb_parser = PDBParser(QUIET=True)


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


def _residues_with_ca(chain):
    return [r for r in chain if "CA" in r and r.id[0] == " "]


def _chain_one_letter_sequence(residues):
    return "".join(seq1(r.resname, undef_code="X") for r in residues)


def align_reference_onto_chain(pair_pdb_text, target_chain_id, reference_pdb_text, prune_fraction=0.3):
    """Sequence-guided structural superposition of an entire reference
    monomer structure (a human/yeast ortholog, from ui.render_reference_
    structure_picker) onto one chain of a pair's predicted structure.

    Unlike superpose_chains_as_pdb (repeat predictions of the *same*
    Capsaspora protein, so CA atom counts already match 1:1), a cross-
    species ortholog has a different sequence and residue count, so
    Superimposer can't be handed same-length CA lists directly -- a global
    pairwise alignment (BLOSUM62) finds the residue correspondence first,
    and only the aligned (non-gap) columns are superposed.

    prune_fraction: after an initial fit, the worst-fitting this fraction of
    aligned residue pairs are dropped and the fit redone once -- divergent
    orthologs often align well over a conserved core but poorly at
    peripheral loops/termini, and a single outlier-pruning pass gives a
    visually tighter core superposition (at the cost of a look that no
    longer reflects the full aligned length in its RMSD).

    Returns (aligned_reference_pdb_text, n_core_residues, core_rmsd), or
    None if there's no usable overlap (fewer than 3 aligned residues).
    """
    pair_structure = _pdb_parser.get_structure("pair", io.StringIO(pair_pdb_text))
    target_chain = pair_structure[0][target_chain_id]
    reference_structure = _pdb_parser.get_structure("reference", io.StringIO(reference_pdb_text))
    reference_chain = next(reference_structure[0].get_chains())

    target_residues = _residues_with_ca(target_chain)
    reference_residues = _residues_with_ca(reference_chain)
    target_seq = _chain_one_letter_sequence(target_residues)
    reference_seq = _chain_one_letter_sequence(reference_residues)
    if not target_seq or not reference_seq:
        return None

    aligner = Align.PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    aligner.mode = "global"
    alignment = aligner.align(target_seq, reference_seq)[0]

    target_atoms, reference_atoms = [], []
    for t_i, r_i in zip(*alignment.indices):
        if t_i == -1 or r_i == -1:
            continue
        target_atoms.append(target_residues[t_i]["CA"])
        reference_atoms.append(reference_residues[r_i]["CA"])
    if len(target_atoms) < 3:
        return None

    superimposer = Superimposer()
    superimposer.set_atoms(target_atoms, reference_atoms)

    if 0 < prune_fraction < 1 and len(target_atoms) >= 10:
        rotation, translation = superimposer.rotran
        reference_coords = np.array([a.coord for a in reference_atoms]) @ rotation + translation
        target_coords = np.array([a.coord for a in target_atoms])
        residuals = np.linalg.norm(reference_coords - target_coords, axis=1)
        keep_n = max(10, int(len(residuals) * (1 - prune_fraction)))
        keep_idx = np.argsort(residuals)[:keep_n]
        target_atoms = [target_atoms[i] for i in keep_idx]
        reference_atoms = [reference_atoms[i] for i in keep_idx]
        superimposer.set_atoms(target_atoms, reference_atoms)

    superimposer.apply(list(reference_structure.get_atoms()))

    buf = io.StringIO()
    io_writer = PDBIO()
    io_writer.set_structure(reference_structure)
    io_writer.save(buf)
    return buf.getvalue(), len(target_atoms), superimposer.rms


def find_contact_residues(pdb_text, chain_a_id, chain_b_id, cutoff=8.0):
    """Residue pairs from different chains of a pair's predicted structure
    whose CA atoms are within cutoff Angstroms -- a standard, simple proxy
    for "these residues are near the interface" (not a precise heavy-atom
    contact definition, but a good enough guide for browsing where an
    interaction happens, and cheap: only ~1 atom/residue to compare).

    Returns a DataFrame [resi_a, resn_a, resi_b, resn_b, distance], one row
    per contacting residue pair, sorted by distance ascending; empty if
    either chain has no standard residues or nothing is within cutoff.
    """
    structure = _pdb_parser.get_structure("pair", io.StringIO(pdb_text))
    chain_a = structure[0][chain_a_id]
    chain_b = structure[0][chain_b_id]

    columns = ["resi_a", "resn_a", "resi_b", "resn_b", "distance"]

    def ca_table(chain):
        residues = _residues_with_ca(chain)
        resi = np.array([r.id[1] for r in residues])
        resn = np.array([r.resname for r in residues])
        coords = np.array([r["CA"].coord for r in residues])
        return resi, resn, coords

    resi_a, resn_a, coords_a = ca_table(chain_a)
    resi_b, resn_b, coords_b = ca_table(chain_b)
    if len(coords_a) == 0 or len(coords_b) == 0:
        return pd.DataFrame(columns=columns)

    dists = np.linalg.norm(coords_a[:, None, :] - coords_b[None, :, :], axis=2)
    ai, bi = np.where(dists <= cutoff)
    if len(ai) == 0:
        return pd.DataFrame(columns=columns)

    records = [
        {
            "resi_a": int(resi_a[a]),
            "resn_a": resn_a[a],
            "resi_b": int(resi_b[b]),
            "resn_b": resn_b[b],
            "distance": float(dists[a, b]),
        }
        for a, b in zip(ai, bi)
    ]
    return pd.DataFrame(records).sort_values("distance").reset_index(drop=True)
