"""py3Dmol rendering helpers for a two-chain protein pair."""

import py3Dmol

CHAIN_COLORS = ["#4C78A8", "#F58518"]  # protein A, protein B

# Distinct colors for comparing several structures of the same protein at
# once (render_overlay/render_single) -- Vega's "category10" palette, reused
# here for visual consistency with the interactome graph's colors elsewhere.
COMPARISON_COLORS = [
    "#4C78A8", "#F58518", "#54A24B", "#E45756",
    "#B279A2", "#FF9DA6", "#9D755D", "#BAB0AC",
]

PLDDT_COLORSCHEME = {
    "prop": "b",
    "gradient": "roygb",
    "min": 50,
    "max": 90,
}


def render_pair(pdb_text, chain_a, chain_b, color_by="chain", width=760, height=560):
    """Build a py3Dmol view of a two-chain complex, returned as embeddable HTML."""
    view = py3Dmol.view(width=width, height=height)
    view.addModel(pdb_text, "pdb")

    if color_by == "plddt":
        view.setStyle({"chain": chain_a}, {"cartoon": {"colorscheme": PLDDT_COLORSCHEME}})
        view.setStyle({"chain": chain_b}, {"cartoon": {"colorscheme": PLDDT_COLORSCHEME}})
    else:
        view.setStyle({"chain": chain_a}, {"cartoon": {"color": CHAIN_COLORS[0]}})
        view.setStyle({"chain": chain_b}, {"cartoon": {"color": CHAIN_COLORS[1]}})

    view.zoomTo()
    view.spin(False)
    return view._make_html()


def render_overlay(pdb_texts, width=760, height=560):
    """Several already-superposed single-chain PDB texts in one view, each a
    distinct color -- for visually comparing how consistently AF3 predicts
    the same protein's fold across different pools (see
    eukaryoma_ppi.structures.superpose_chains_as_pdb).
    """
    view = py3Dmol.view(width=width, height=height)
    for i, pdb_text in enumerate(pdb_texts):
        view.addModel(pdb_text, "pdb")
        view.setStyle({"model": i}, {"cartoon": {"color": COMPARISON_COLORS[i % len(COMPARISON_COLORS)]}})
    view.zoomTo()
    view.spin(False)
    return view._make_html()


def render_single(pdb_text, color, width=320, height=320):
    """One structure, one flat color -- for side-by-side comparison panels
    (no alignment, unlike render_overlay)."""
    view = py3Dmol.view(width=width, height=height)
    view.addModel(pdb_text, "pdb")
    view.setStyle({"cartoon": {"color": color}})
    view.zoomTo()
    view.spin(False)
    return view._make_html()
