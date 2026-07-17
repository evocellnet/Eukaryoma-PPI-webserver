"""Thin subprocess wrapper around the `foldcomp` CLI binary.

We deliberately shell out to the binary instead of using the `foldcomp`
PyPI package: the pool .fcz files were written by a custom internal build
(magic number "FCZC") that neither the pip package nor the public
steineggerlab/foldcomp CLI releases (magic number "FCMP") can read. Only
the binary bundled at bin/foldcomp is known to decompress them correctly.
"""

import subprocess

from eukaryoma_ppi.config import FOLDCOMP_BIN


class FoldcompError(RuntimeError):
    pass


def decompress_pool(fcz_path, out_cif_path):
    """Decompress one pooled-complex .fcz file to mmCIF via the foldcomp CLI.

    The CLI infers the output format from the destination file extension,
    so `out_cif_path` must end in .cif.
    """
    if not str(out_cif_path).endswith(".cif"):
        raise ValueError("out_cif_path must end with .cif")

    result = subprocess.run(
        [str(FOLDCOMP_BIN), "decompress", str(fcz_path), str(out_cif_path)],
        capture_output=True,
        text=True,
    )
    if result.returncode != 0 or not out_cif_path.exists():
        raise FoldcompError(
            f"foldcomp decompress failed for {fcz_path}: "
            f"{result.stdout.strip()} {result.stderr.strip()}"
        )
