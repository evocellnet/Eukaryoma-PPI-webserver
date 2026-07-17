"""py3Dmol rendering helpers for a two-chain protein pair."""

import py3Dmol

CHAIN_COLORS = ["#4C78A8", "#F58518"]  # protein A, protein B

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
