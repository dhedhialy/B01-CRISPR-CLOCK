"""Leakage-free calibration / holdout splits."""

from __future__ import annotations

from typing import Dict, List, Sequence, Tuple

import numpy as np

from b01_clock.model.observations import CharacterMatrix, EDITED, UNEDITED, MISSING
from b01_clock.model.tree import LineageTree
from b01_clock.model.exposure import exposure_from_edit_fraction


def split_sites(
    n_sites: int,
    calib_frac: float = 0.4,
    seed: int = 0,
) -> Tuple[np.ndarray, np.ndarray]:
    """Pre-register site split before headline timing. Returns (calib_idx, analysis_idx)."""
    rng = np.random.default_rng(seed)
    perm = rng.permutation(n_sites)
    n_cal = max(1, int(round(calib_frac * n_sites)))
    calib = np.sort(perm[:n_cal])
    analysis = np.sort(perm[n_cal:])
    if len(analysis) == 0:
        analysis = calib.copy()
    return calib, analysis


def held_out_rate_bounds(
    tree: LineageTree,
    matrix: CharacterMatrix,
    calib_sites: Sequence[int],
    total_time: float,
    quantile: Tuple[float, float] = (0.1, 0.9),
    pad: float = 1.1,
    floor: float = 1e-3,
) -> Tuple[float, float]:
    """Estimate transferable global rate bounds from calibration sites only.

    Uses leaf edit fractions and dated total time: r̂ ≈ (-log(1-p̂)) / T.
    Returns a shared (lo, hi) box for the global-rate solver (leakage-free transfer).
    """
    leaves = tree.leaves()
    rates = []
    for s in calib_sites:
        edited = 0
        observed = 0
        for ℓ in leaves:
            st = int(matrix.leaf_states[ℓ][s])
            if st == MISSING:
                continue
            observed += 1
            if st == EDITED:
                edited += 1
        if observed == 0:
            continue
        p = edited / observed
        p = min(p, 1.0 - 1e-6)
        lam = exposure_from_edit_fraction(p)
        rates.append(lam / max(total_time, 1e-12))
    if not rates:
        return (floor, 10.0)
    lo = float(np.quantile(rates, quantile[0]) / pad)
    hi = float(np.quantile(rates, quantile[1]) * pad)
    lo = max(floor, lo)
    hi = max(lo * 1.01, hi)
    return (lo, hi)


def lock_split_record(
    calib: Sequence[int],
    analysis: Sequence[int],
    seed: int,
    note: str = "",
) -> dict:
    return {
        "calib_sites": list(map(int, calib)),
        "analysis_sites": list(map(int, analysis)),
        "seed": seed,
        "locked_before_headline": True,
        "note": note,
    }
