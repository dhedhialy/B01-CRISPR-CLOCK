"""Likelihood and rate-bound helpers shared by solvers."""

from __future__ import annotations

from typing import Dict, Optional, Sequence, Tuple, Union

import numpy as np

from b01_clock.model.observations import CharacterMatrix, conservative_ancestral_states
from b01_clock.model.tree import LineageTree
from b01_clock.solver.pruning import (
    matrix_loglik_pruning,
    profile_loglik_global_rate,
    profile_loglik_per_site,
)

RateBounds = Union[Tuple[float, float], Sequence[Tuple[float, float]]]


def tree_loglik(
    tree: LineageTree,
    ancestral: Dict[str, np.ndarray],
    rates: np.ndarray,
    times: Optional[np.ndarray] = None,
    matrix: Optional[CharacterMatrix] = None,
) -> float:
    if matrix is not None:
        return matrix_loglik_pruning(tree, matrix, rates, times)
    from b01_clock.model.exposure import loglik_transition

    if times is not None:
        tree = tree.copy()
        tree.set_branch_times(times)
    total = 0.0
    for s in range(rates.shape[0]):
        r = float(rates[s])
        for p, c in tree.edges():
            ll = loglik_transition(
                int(ancestral[p][s]),
                int(ancestral[c][s]),
                r,
                tree.nodes[c].time_from_parent,
            )
            if not np.isfinite(ll):
                return -np.inf
            total += ll
    return float(total)


def _as_global_bound(rate_bounds: RateBounds) -> Tuple[float, float]:
    if isinstance(rate_bounds, tuple) and len(rate_bounds) == 2 and not isinstance(rate_bounds[0], tuple):
        return float(rate_bounds[0]), float(rate_bounds[1])  # type: ignore[arg-type]
    bounds = list(rate_bounds)  # type: ignore[arg-type]
    lo = min(b[0] for b in bounds)
    hi = max(b[1] for b in bounds)
    return float(lo), float(hi)


def profile_loglik_over_rates(
    tree: LineageTree,
    ancestral: Dict[str, np.ndarray],
    times: np.ndarray,
    rate_bounds: RateBounds,
    matrix: Optional[CharacterMatrix] = None,
    rate_mode: str = "global",
    site_weights: Optional[np.ndarray] = None,
) -> Tuple[float, np.ndarray]:
    """Maximize L over rates for fixed times.

    rate_mode:
      - \"global\": one rate in the bound (default; prevents false ID via per-site freedom)
      - \"per_site\": each site has its own box (requires tight/calibrated bounds)
    """
    if matrix is None:
        raise ValueError("matrix is required for profile likelihood")
    if rate_mode == "global":
        ll, r = profile_loglik_global_rate(
            tree, matrix, times, _as_global_bound(rate_bounds), site_weights=site_weights
        )
        w = np.ones(matrix.n_sites) if site_weights is None else site_weights
        return ll, r * np.asarray(w, dtype=float)
    return profile_loglik_per_site(tree, matrix, times, list(rate_bounds))  # type: ignore[arg-type]


def prepare_ancestral(tree: LineageTree, matrix: CharacterMatrix) -> Dict[str, np.ndarray]:
    if matrix.ancestral is not None:
        return matrix.ancestral
    return conservative_ancestral_states(tree, matrix)


def ultrametric_times_from_free(
    tree: LineageTree,
    free_internal: np.ndarray,
    total_time: float,
) -> np.ndarray:
    internals = [n for n in tree.internal_nodes() if n != tree.root]
    if len(free_internal) != len(internals):
        raise ValueError("free_internal length must match non-root internal nodes")
    ages = {tree.root: 0.0}
    for name, age in zip(internals, free_internal):
        ages[name] = float(np.clip(age, 0.0, total_time))
    for leaf in tree.leaves():
        ages[leaf] = total_time
    # Enforce age child >= age parent in parent-before-child order. A single pass
    # over an arbitrary edge order can RAISE a parent AFTER its child's branch time
    # was committed, so cumulative re-sums overflow total_time (observed mle > T).
    frontier = [tree.root]
    for u in frontier:
        for v in tree.nodes[u].children:
            if ages[v] < ages[u]:
                ages[v] = ages[u]
            frontier.append(v)
    times = [ages[c] - ages[p] for p, c in tree.edges()]
    return np.asarray(times, dtype=float)


def default_rate_bounds(n_sites: int, lo: float = 0.1, hi: float = 10.0) -> list[Tuple[float, float]]:
    return [(lo, hi)] * n_sites
