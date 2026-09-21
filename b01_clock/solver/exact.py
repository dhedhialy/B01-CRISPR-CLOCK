"""Exact / certified small-tree feasible-set solver."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
from scipy.stats import chi2

from b01_clock.model.observations import CharacterMatrix
from b01_clock.model.tree import LineageTree
from b01_clock.solver.likelihood import (
    prepare_ancestral,
    profile_loglik_over_rates,
    ultrametric_times_from_free,
)


@dataclass
class ExactSolverResult:
    target: str
    lower: float
    upper: float
    mle: float
    status: str
    l_max: float
    threshold: float
    diagnostics: dict = field(default_factory=dict)
    profile_grid: Optional[np.ndarray] = None
    profile_ll: Optional[np.ndarray] = None


def _target_value(tree: LineageTree, times: np.ndarray, target: str) -> float:
    tree = tree.copy()
    tree.set_branch_times(times)
    if target.startswith("age:"):
        return tree.node_age(target.split(":", 1)[1])
    if target.startswith("branch:"):
        # branch:parent->child
        spec = target.split(":", 1)[1]
        p, c = spec.split("->")
        return tree.nodes[c].time_from_parent
    if target.startswith("order_gap:"):
        # order_gap:A,B → age(B) - age(A)
        a, b = target.split(":", 1)[1].split(",")
        return tree.node_age(b) - tree.node_age(a)
    raise ValueError(f"Unknown target {target}")


def exact_profile_bounds(
    tree: LineageTree,
    matrix: CharacterMatrix,
    rate_bounds: Sequence[Tuple[float, float]] | Tuple[float, float],
    total_time: float,
    targets: Sequence[str],
    confidence: float = 0.95,
    grid_size: int = 41,
    free_internal_names: Optional[Sequence[str]] = None,
    rate_mode: str = "global",
) -> List[ExactSolverResult]:
    """Compute sharp profile-likelihood bounds on small ultrametric trees.

    Free parameters = ages of non-root internal nodes in [0, T]; leaves fixed at T.
    Rates are profiled within bounds (primary guardrail vs false absolute-time ID).
    Default rate_mode=\"global\" profiles one rate in the bound.
    """
    ancestral = prepare_ancestral(tree, matrix)
    internals = list(free_internal_names) if free_internal_names is not None else [
        n for n in tree.internal_nodes() if n != tree.root
    ]
    d = len(internals)
    if d == 0:
        raise ValueError("Tree has no free internal nodes")
    if d > 3:
        raise ValueError("Exact grid solver supports ≤3 free internals; use scalable solver")

    alpha = 1.0 - confidence
    cutoff = 0.5 * float(chi2.ppf(1.0 - alpha, df=1))

    axis = np.linspace(0.0, total_time, grid_size)
    if d == 1:
        free_grid = axis.reshape(-1, 1)
    elif d == 2:
        a0, a1 = np.meshgrid(axis, axis, indexing="ij")
        free_grid = np.column_stack([a0.ravel(), a1.ravel()])
    else:
        pts = [(i, j, k) for i in axis for j in axis for k in axis]
        free_grid = np.asarray(pts)

    best_ll = -np.inf
    best_free = None
    lls = np.full(free_grid.shape[0], -np.inf)

    for i, free in enumerate(free_grid):
        times = ultrametric_times_from_free(tree, free, total_time)
        if np.any(times < -1e-12):
            continue
        ll, _ = profile_loglik_over_rates(
            tree, ancestral, times, rate_bounds, matrix=matrix, rate_mode=rate_mode
        )
        lls[i] = ll
        if ll > best_ll:
            best_ll = ll
            best_free = free.copy()

    if best_free is None or not np.isfinite(best_ll):
        return [
            ExactSolverResult(
                target=t,
                lower=0.0,
                upper=total_time,
                mle=float("nan"),
                status="infeasible",
                l_max=-np.inf,
                threshold=-np.inf,
                diagnostics={"reason": "no finite likelihood on grid"},
            )
            for t in targets
        ]

    threshold = best_ll - cutoff
    feasible_mask = lls >= threshold - 1e-9
    feasible_free = free_grid[feasible_mask]

    results = []
    for target in targets:
        values = []
        for free in feasible_free:
            times = ultrametric_times_from_free(tree, free, total_time)
            values.append(_target_value(tree, times, target))
        values = np.asarray(values, dtype=float)
        mle_times = ultrametric_times_from_free(tree, best_free, total_time)
        mle = _target_value(tree, mle_times, target)
        lo, hi = float(values.min()), float(values.max())
        width = hi - lo
        if width <= 1e-8:
            status = "point_identified"
        elif width >= total_time - 1e-6:
            status = "unidentified"
        else:
            status = "partially_identified"
        results.append(
            ExactSolverResult(
                target=target,
                lower=lo,
                upper=hi,
                mle=float(mle),
                status=status,
                l_max=float(best_ll),
                threshold=float(threshold),
                diagnostics={
                    "confidence": confidence,
                    "cutoff": cutoff,
                    "n_feasible_grid": int(feasible_mask.sum()),
                    "n_grid": int(len(lls)),
                    "free_internals": list(internals),
                    "total_time": total_time,
                    "rate_mode": rate_mode,
                    "certificate": (
                        f"Profile LR set {{θ: 2(Lmax-L(θ)) ≤ χ²_{{1,{confidence:.2f}}}}} "
                        f"with rates profiled in given bounds (mode={rate_mode}); "
                        f"bound attained on grid."
                    ),
                },
            )
        )
    return results


def verify_against_grid(
    tree: LineageTree,
    matrix: CharacterMatrix,
    rate_bounds: Sequence[Tuple[float, float]],
    total_time: float,
    target: str,
    result: ExactSolverResult,
    fine_grid: int = 81,
    tol: float = 0.05,
) -> dict:
    """Recompute bounds on a finer grid; check agreement within tol * T."""
    fine = exact_profile_bounds(
        tree,
        matrix,
        rate_bounds,
        total_time,
        targets=[target],
        confidence=float(result.diagnostics.get("confidence", 0.95)),
        grid_size=fine_grid,
    )[0]
    T = total_time
    ok = abs(fine.lower - result.lower) <= tol * T and abs(fine.upper - result.upper) <= tol * T
    return {
        "ok": ok,
        "coarse": (result.lower, result.upper),
        "fine": (fine.lower, fine.upper),
        "abs_diff": (abs(fine.lower - result.lower), abs(fine.upper - result.upper)),
    }
