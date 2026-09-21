"""Felsenstein pruning likelihood for irreversible 2-state editing."""

from __future__ import annotations

from typing import Dict, Optional, Sequence, Tuple

import numpy as np
from scipy.optimize import minimize_scalar

from b01_clock.model.exposure import EPS
from b01_clock.model.observations import EDITED, MISSING, UNEDITED, CharacterMatrix
from b01_clock.model.tree import LineageTree


def _transition_log_scalar(parent: int, child: int, rt: float) -> float:
    if parent == EDITED:
        return 0.0 if child == EDITED else -np.inf
    if child == UNEDITED:
        return -rt
    if rt < 1e-8:
        return float(np.log(max(rt, EPS)))
    if rt > 30:
        return 0.0
    return float(np.log1p(-np.exp(-rt)))


def site_loglik_pruning(
    tree: LineageTree,
    leaf_states: Dict[str, int],
    rate: float,
) -> float:
    partial: Dict[str, np.ndarray] = {}

    def postorder(u: str) -> np.ndarray:
        if u in partial:
            return partial[u]
        if tree.nodes[u].is_leaf or not tree.nodes[u].children:
            st = leaf_states[u]
            if st == MISSING:
                out = np.array([0.0, 0.0])
            elif st == UNEDITED:
                out = np.array([0.0, -np.inf])
            else:
                out = np.array([-np.inf, 0.0])
            partial[u] = out
            return out
        acc0, acc1 = 0.0, 0.0
        for v in tree.nodes[u].children:
            child_ll = postorder(v)
            rt = max(rate, 0.0) * max(tree.nodes[v].time_from_parent, 0.0)
            terms0 = [
                _transition_log_scalar(UNEDITED, cs, rt) + child_ll[cs]
                for cs in (UNEDITED, EDITED)
            ]
            terms1 = [
                _transition_log_scalar(EDITED, cs, rt) + child_ll[cs]
                for cs in (UNEDITED, EDITED)
            ]
            acc0 += float(np.logaddexp(terms0[0], terms0[1]))
            acc1 += float(terms1[1]) if np.isfinite(terms1[1]) else -np.inf
        out = np.array([acc0, acc1])
        partial[u] = out
        return out

    return float(postorder(tree.root)[UNEDITED])


def _log_edit_prob(rt: np.ndarray) -> np.ndarray:
    """Stable log(1 - exp(-rt))."""
    out = np.empty_like(rt, dtype=float)
    small = rt < 1e-8
    large = rt > 30
    mid = ~small & ~large
    out[small] = np.log(np.maximum(rt[small], EPS))
    out[large] = 0.0
    out[mid] = np.log1p(-np.exp(-rt[mid]))
    return out


def matrix_loglik_vectorized(
    tree: LineageTree,
    leaf_array: np.ndarray,
    leaf_order: Sequence[str],
    rate: float,
) -> float:
    """Pruning over all sites at once for a single shared rate.

    leaf_array: shape (n_leaves, n_sites) with values in {-1,0,1}
    """
    n_sites = leaf_array.shape[1]
    # partial[node] -> (2, n_sites) log-likelihood given state
    partial: Dict[str, np.ndarray] = {}
    leaf_index = {ℓ: i for i, ℓ in enumerate(leaf_order)}

    def postorder(u: str) -> np.ndarray:
        if u in partial:
            return partial[u]
        if tree.nodes[u].is_leaf or not tree.nodes[u].children:
            row = leaf_array[leaf_index[u]]
            out = np.zeros((2, n_sites))
            # missing → 0,0 ; unedited → [0,-inf]; edited → [-inf,0]
            miss = row == MISSING
            un = row == UNEDITED
            ed = row == EDITED
            out[1, un] = -np.inf
            out[0, ed] = -np.inf
            # missing already 0
            _ = miss
            partial[u] = out
            return out

        acc = np.zeros((2, n_sites))
        for v in tree.nodes[u].children:
            child = postorder(v)
            rt = max(rate, 0.0) * max(tree.nodes[v].time_from_parent, 0.0)
            # from parent=0
            t00 = -rt + child[0]
            t01 = _log_edit_prob(np.full(n_sites, rt)) + child[1]
            # from parent=1: only child=1 with prob 1
            t11 = child[1].copy()
            t10 = np.full(n_sites, -np.inf)
            acc[0] += np.logaddexp(t00, t01)
            acc[1] += t11
            _ = t10
        partial[u] = acc
        return acc

    root = postorder(tree.root)
    return float(np.sum(root[0]))


def _leaf_array(matrix: CharacterMatrix, leaf_order: Sequence[str]) -> np.ndarray:
    return np.stack([matrix.leaf_states[ℓ] for ℓ in leaf_order], axis=0)


def matrix_loglik_pruning(
    tree: LineageTree,
    matrix: CharacterMatrix,
    rates: np.ndarray,
    branch_times: Optional[np.ndarray] = None,
) -> float:
    tree = tree.copy()
    if branch_times is not None:
        tree.set_branch_times(branch_times)
    leaves = tree.leaves()
    arr = _leaf_array(matrix, leaves)
    # If all rates equal, one vectorized call
    if np.allclose(rates, rates[0]):
        return matrix_loglik_vectorized(tree, arr, leaves, float(rates[0]))
    total = 0.0
    for s in range(matrix.n_sites):
        leaf_states = {ℓ: int(arr[i, s]) for i, ℓ in enumerate(leaves)}
        ll = site_loglik_pruning(tree, leaf_states, float(rates[s]))
        if not np.isfinite(ll):
            return -np.inf
        total += ll
    return float(total)


def profile_loglik_global_rate(
    tree: LineageTree,
    matrix: CharacterMatrix,
    branch_times: np.ndarray,
    rate_bound: Tuple[float, float],
    site_weights: Optional[np.ndarray] = None,
) -> Tuple[float, float]:
    tree = tree.copy()
    tree.set_branch_times(branch_times)
    lo, hi = rate_bound
    leaves = tree.leaves()
    arr = _leaf_array(matrix, leaves)
    w = None if site_weights is None else np.asarray(site_weights, dtype=float)

    if w is None or np.allclose(w, 1.0):
        def neg(r: float) -> float:
            ll = matrix_loglik_vectorized(tree, arr, leaves, r)
            return -ll if np.isfinite(ll) else 1e300
    else:
        def neg(r: float) -> float:
            rates = r * w
            ll = matrix_loglik_pruning(tree, matrix, rates, None)
            return -ll if np.isfinite(ll) else 1e300

    res = minimize_scalar(neg, bounds=(lo, hi), method="bounded", options={"xatol": 1e-10})
    return -float(res.fun), float(res.x)


def profile_loglik_per_site(
    tree: LineageTree,
    matrix: CharacterMatrix,
    branch_times: np.ndarray,
    rate_bounds: Sequence[Tuple[float, float]],
) -> Tuple[float, np.ndarray]:
    tree = tree.copy()
    tree.set_branch_times(branch_times)
    rates = np.zeros(matrix.n_sites)
    total = 0.0
    for s, (lo, hi) in enumerate(rate_bounds):
        leaf_states = {ℓ: int(matrix.leaf_states[ℓ][s]) for ℓ in tree.leaves()}
        if all(leaf_states[ℓ] == MISSING for ℓ in tree.leaves()):
            rates[s] = 0.5 * (lo + hi)
            continue

        def neg(r: float) -> float:
            ll = site_loglik_pruning(tree, leaf_states, r)
            return -ll if np.isfinite(ll) else 1e300

        res = minimize_scalar(neg, bounds=(lo, hi), method="bounded", options={"xatol": 1e-6})
        rates[s] = float(res.x)
        total += -float(res.fun)
    return total, rates
