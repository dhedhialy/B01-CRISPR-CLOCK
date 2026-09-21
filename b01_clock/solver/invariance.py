"""Rate–time invariance guardrails against false absolute-time identification."""

from __future__ import annotations

from typing import Dict, Sequence, Tuple

import numpy as np

from b01_clock.model.observations import CharacterMatrix
from b01_clock.model.tree import LineageTree
from b01_clock.solver.likelihood import prepare_ancestral, tree_loglik


def check_rate_time_invariance(
    tree: LineageTree,
    matrix: CharacterMatrix,
    rates: np.ndarray,
    times: np.ndarray,
    scale: float = 2.0,
    atol: float = 1e-8,
) -> dict:
    """Likelihood must be unchanged under r → c r, t → t / c (same exposures)."""
    ancestral = prepare_ancestral(tree, matrix)
    ll0 = tree_loglik(tree, ancestral, rates, times, matrix=matrix)
    ll1 = tree_loglik(tree, ancestral, rates * scale, times / scale, matrix=matrix)
    ok = np.isfinite(ll0) and np.isfinite(ll1) and abs(ll0 - ll1) <= atol * (1 + abs(ll0))
    return {
        "ok": bool(ok),
        "ll_base": float(ll0),
        "ll_scaled": float(ll1),
        "abs_diff": float(abs(ll0 - ll1)),
        "scale": scale,
        "message": (
            "Rate–time invariance holds: absolute chronology is not identified from "
            "exposures alone."
            if ok
            else "Invariance violated — check implementation."
        ),
    }


def rate_time_transform_feasible(
    times: np.ndarray,
    rates: np.ndarray,
    rate_bounds,
    total_time: float,
    scale: float,
) -> dict:
    """Primary computational guardrail: after r→c r, t→t/c, do rate bounds + T still hold?

    If a proposed point-identification relies on a transform that exits the locked
    rate bounds or dated-anchor total duration, absolute time is not identified.
    """
    new_rates = rates * scale
    new_times = times / scale
    # Normalize bounds to per-site list
    if (
        isinstance(rate_bounds, tuple)
        and len(rate_bounds) == 2
        and not isinstance(rate_bounds[0], (tuple, list))
    ):
        lo, hi = float(rate_bounds[0]), float(rate_bounds[1])
        bounds_list = [(lo, hi)] * len(new_rates)
    else:
        bounds_list = list(rate_bounds)
    rates_ok = all(lo <= r <= hi for r, (lo, hi) in zip(new_rates, bounds_list))
    depth = float(np.sum(times))  # scale proxy
    new_depth = depth / scale
    feasible = bool(rates_ok and abs(scale - 1.0) <= 1e-12)
    feasible_without_anchor = bool(rates_ok and np.all(new_times >= -1e-12))
    return {
        "scale": scale,
        "rates_in_bounds": bool(rates_ok),
        "preserves_dated_anchor": bool(abs(new_depth - depth) <= 1e-9),
        "feasible_with_anchor_and_bounds": feasible,
        "feasible_without_anchor": feasible_without_anchor,
        "guardrail": (
            "False identification blocked: nontrivial rate–time rescaling exits the "
            "locked rate bounds and/or dated anchor."
            if not feasible
            else "Transform remains inside locked constraints."
        ),
    }
