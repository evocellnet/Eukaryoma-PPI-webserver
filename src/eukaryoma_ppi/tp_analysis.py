"""Analysis helpers for the Annotations page: how well each score predicts
the CORUM/Marcotte "true positive" complex co-membership annotation.
"""

import numpy as np
import pandas as pd

GROUND_TRUTH_MODES = {
    "CORUM only": lambda df: df["corum_tp"],
    "Marcotte only": lambda df: df["marcotte_tp"],
    "Either (CORUM or Marcotte)": lambda df: df["corum_tp"] | df["marcotte_tp"],
    "Both (CORUM and Marcotte)": lambda df: df["corum_tp"] & df["marcotte_tp"],
}


def ground_truth_series(df, mode):
    return GROUND_TRUTH_MODES[mode](df)


def quantile_bins(scores, n_bins):
    """1..n_bins, low score -> 1. Rank-based (not pd.qcut) so heavily-tied
    correlation scores near 1.0 don't collapse bin edges.
    """
    ranks = scores.rank(method="first", pct=True)
    return np.minimum((ranks * n_bins).astype(int), n_bins - 1) + 1


def quantile_tp_rates(df, score_col, ground_truth, n_bins=10):
    """TP rate per score quantile bin for one score column.

    Returns a DataFrame [quantile, tp_rate, n_pairs] restricted to rows where
    score_col is non-null; empty if there's no data for this score.
    """
    valid = df[score_col].notna()
    scores, gt = df.loc[valid, score_col], ground_truth.loc[valid]
    if len(scores) == 0:
        return pd.DataFrame(columns=["quantile", "tp_rate", "n_pairs"])

    bins = quantile_bins(scores, n_bins)
    return (
        pd.DataFrame({"quantile": bins, "tp": gt.to_numpy()})
        .groupby("quantile")["tp"]
        .agg(tp_rate="mean", n_pairs="count")
        .reset_index()
    )


def quantile_tp_rates_all_scores(df, score_columns, ground_truth, n_bins=10):
    """Long-format [score, quantile, tp_rate, n_pairs] across several score columns."""
    frames = []
    for score_col, label in score_columns:
        rates = quantile_tp_rates(df, score_col, ground_truth, n_bins)
        if rates.empty:
            continue
        rates["score"] = label
        frames.append(rates)
    if not frames:
        return pd.DataFrame(columns=["score", "quantile", "tp_rate", "n_pairs"])
    return pd.concat(frames, ignore_index=True)


def tp_category(df):
    """Categorical label per pair for coloring plots: which source(s) flag it."""
    corum, marcotte = df["corum_tp"], df["marcotte_tp"]
    return np.select(
        [corum & marcotte, corum, marcotte],
        ["Both", "CORUM", "Marcotte"],
        default="Not annotated",
    )


def baseline_score_density(df, score_col, n_bins=40):
    """Histogram-based score density for pairs not flagged true-positive by
    any source -- the baseline the true-positive rows should stand out
    against. A histogram is computed over the *full* not-annotated
    population (fast even at millions of rows -- no need to subsample);
    only the resulting n_bins-row curve is ever plotted.

    Returns a DataFrame [score_bin, density, neg_density] (the latter just
    -density, for mirroring into a violin shape); empty if there's no data.
    """
    not_annotated = ~(df["corum_tp"] | df["marcotte_tp"])
    scores = df.loc[not_annotated & df[score_col].notna(), score_col]
    if scores.empty:
        return pd.DataFrame(columns=["score_bin", "density", "neg_density"])

    counts, edges = np.histogram(scores, bins=n_bins, density=True)
    centers = (edges[:-1] + edges[1:]) / 2
    return pd.DataFrame({"score_bin": centers, "density": counts, "neg_density": -counts})
