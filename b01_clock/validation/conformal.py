"""Empirical (conformal) recalibration helpers shared by validation scripts.

The profile-LR certificate is miscalibrated on dense CRISPR data (systematic
late-shift, width too narrow). These helpers build binwise-conformal widths
from a benchmark's own internal-node residuals, so intervals can be widened
to honest marginal coverage. All quantiles exclude the evaluated point.

Bin = MLE decile cell (depth proxy). Width = standard conformal quantile of
the bin's residual margins (max(true-lo, hi-true)), with the test point
excluded from its own bin. See scripts/recalibration_conditional.py.
"""

from __future__ import annotations

from typing import Sequence, Tuple

import numpy as np

N_BINS = 9
ALPHA = 0.05


def bin_of(mle: float, n_bins: int = N_BINS) -> int:
    edges = np.linspace(0.0, 1.0, n_bins + 1)
    for i in range(n_bins):
        if edges[i] <= mle < edges[i + 1]:
            return i
    return n_bins - 1


def quantile_width(residuals: Sequence[float], alpha: float = ALPHA) -> float:
    """Standard conformal width from calibration residuals (test point excluded beforehand)."""
    s = np.sort(np.asarray(residuals, dtype=float))
    if s.size == 0:
        return 0.0
    k = min(int(np.ceil((s.size + 1) * (1 - alpha))), s.size)
    return float(s[k - 1])


def conditional_loo(rows, alpha: float = ALPHA) -> Tuple[float, float]:
    """True LOOCV coverage: each row's width from its bin excluding itself."""
    cov = 0
    per_bin = {b: [x for x in rows if bin_of(x["mle"]) == b] for b in range(N_BINS)}
    widths = []
    for x in rows:
        excl = [y for y in per_bin[bin_of(x["mle"])] if y is not x]
        w = quantile_width([y["resid"] for y in excl], alpha)
        widths.append(w)
        if x["mle"] - w <= x["true"] <= x["mle"] + w:
            cov += 1
    return cov / len(rows), float(np.mean(widths))


def split_widths(train, alpha: float = ALPHA):
    """Width per bin from a training set (split-conformal); test points disjoint."""
    return [quantile_width([x["resid"] for x in train if bin_of(x["mle"]) == b], alpha)
            for b in range(N_BINS)]


def calibrated_interval(bins: Sequence[float], mle: float):
    """[mle - w, mle + w] for a NEW point, purely function of its MLE (transfer)."""
    w = bins[bin_of(mle)]
    return mle - w, mle + w, w


def conditional_split(train, test, alpha: float = ALPHA):
    """Split-conformal: width from train bins only, evaluated on disjoint test rows."""
    bins = split_widths(train, alpha)
    cov = 0
    for x in test:
        lo, hi, _ = calibrated_interval(bins, x["mle"])
        if lo <= x["true"] <= hi:
            cov += 1
    mean_w = float(np.mean([calibrated_interval(bins, x["mle"])[2] for x in test])) if test else 0.0
    return cov / max(1, len(test)), mean_w, bins