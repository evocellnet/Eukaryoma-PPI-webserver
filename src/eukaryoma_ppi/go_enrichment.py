"""GO term enrichment for a foreground set of proteins (e.g. a Protein View
interactome) against the full eggNOG-annotated website background.

Implemented as a plain hypergeometric test + Benjamini-Hochberg FDR, using
only Python's stdlib `math.comb` -- no goatools/scipy dependency. This
doesn't know GO term *names* (that needs a GO OBO file this project doesn't
ship) so results are reported by GO id, linked out to QuickGO for the name.
"""

import math

import pandas as pd

from eukaryoma_ppi.eggnog import as_list


def _hypergeom_sf(k, K, n, N):
    """P(X >= k) for X ~ Hypergeometric(N, K, n): the chance of seeing at
    least k of the n foreground proteins carry this GO term, if n proteins
    were instead drawn at random (without replacement) from the N-protein
    background, K of which carry the term. Smaller = more enriched.
    """
    if k <= 0:
        return 1.0
    denom = math.comb(N, n)
    total = sum(math.comb(K, i) * math.comb(N - K, n - i) for i in range(k, min(n, K) + 1))
    return min(total / denom, 1.0)


def _benjamini_hochberg(pvalues):
    """FDR q-value per p-value (same order as input), by the standard BH
    step-up procedure: rank ascending, q = p * n / rank, then take a
    running minimum from the largest rank down so q is monotone.
    """
    n = len(pvalues)
    order = sorted(range(n), key=lambda i: pvalues[i])
    q = [0.0] * n
    running_min = 1.0
    for rank in range(n, 0, -1):
        i = order[rank - 1]
        running_min = min(running_min, pvalues[i] * n / rank)
        q[i] = running_min
    return q


def term_frequency(protein_ids, eggnog_df):
    """Simple recap (no statistics): how many of protein_ids carry each GO
    term, sorted by count descending. Returns [go_id, n_proteins, fraction].
    """
    protein_ids = set(protein_ids)
    annotated = eggnog_df[eggnog_df["protein_id"].isin(protein_ids)]
    if annotated.empty:
        return pd.DataFrame(columns=["go_id", "n_proteins", "fraction"])

    counts = {}
    for gos in annotated["GOs"]:
        for go_id in as_list(gos):
            counts[go_id] = counts.get(go_id, 0) + 1
    if not counts:
        return pd.DataFrame(columns=["go_id", "n_proteins", "fraction"])

    df = pd.DataFrame(sorted(counts.items(), key=lambda kv: -kv[1]), columns=["go_id", "n_proteins"])
    df["fraction"] = df["n_proteins"] / len(annotated)
    return df


def enrichment(foreground_ids, eggnog_df, background_ids, min_term_size=1):
    """Hypergeometric GO enrichment of foreground_ids against background_ids
    restricted to those with at least one eggNOG GO annotation (an
    unannotated protein can't meaningfully be "not enriched" for anything).

    background_ids must be the population foreground_ids was actually drawn
    from (e.g. every website protein) -- eggnog_df on its own covers the
    *whole* Capsaspora proteome (6,822 proteins), far more than the ~2,145
    ever selectable on the website, so using its rows directly as the
    background would understate every term's enrichment.

    min_term_size: skip GO terms covering fewer than this many background
    proteins (very rare terms make for unstable, uninterpretable p-values).

    Returns a DataFrame [go_id, n_foreground, n_background, foreground_size,
    background_size, pvalue, fdr], sorted by pvalue ascending -- only GO
    terms present in at least one foreground protein are tested at all.
    """
    eggnog_df = eggnog_df.assign(GOs=eggnog_df["GOs"].map(as_list))
    eggnog_df = eggnog_df[(eggnog_df["GOs"].map(len) > 0) & eggnog_df["protein_id"].isin(background_ids)]

    background_ids = set(eggnog_df["protein_id"])
    foreground_ids = set(foreground_ids) & background_ids
    N, n = len(background_ids), len(foreground_ids)
    if n == 0 or N == 0:
        return pd.DataFrame(columns=["go_id", "n_foreground", "n_background", "pvalue", "fdr"])

    term_members = {}
    for protein_id, gos in zip(eggnog_df["protein_id"], eggnog_df["GOs"]):
        for go_id in gos:
            term_members.setdefault(go_id, set()).add(protein_id)

    rows = []
    for go_id, members in term_members.items():
        K = len(members)
        if K < min_term_size:
            continue
        k = len(members & foreground_ids)
        if k == 0:
            continue
        rows.append({"go_id": go_id, "n_foreground": k, "n_background": K, "pvalue": _hypergeom_sf(k, K, n, N)})

    if not rows:
        return pd.DataFrame(columns=["go_id", "n_foreground", "n_background", "pvalue", "fdr"])

    result = pd.DataFrame(rows)
    result["fdr"] = _benjamini_hochberg(result["pvalue"].tolist())
    result["foreground_size"] = n
    result["background_size"] = N
    return result.sort_values("pvalue").reset_index(drop=True)
