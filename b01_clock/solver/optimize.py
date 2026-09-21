"""Scalable timing-certificate optimization (profile / bound formulation)."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
from scipy.optimize import differential_evolution, minimize, minimize_scalar
from scipy.stats import chi2

from b01_clock.model.observations import CharacterMatrix
from b01_clock.model.tree import LineageTree
from b01_clock.solver.likelihood import (
    prepare_ancestral,
    profile_loglik_over_rates,
    ultrametric_times_from_free,
    tree_loglik,
)


@dataclass
class TimingCertificate:
    target: str
    lower: float
    upper: float
    mle: float
    status: str
    l_max: float
    threshold: float
    solver_status: str
    diagnostics: dict = field(default_factory=dict)


def _pack_objective(
    tree: LineageTree,
    ancestral: Dict[str, np.ndarray],
    rate_bounds: Sequence[Tuple[float, float]],
    total_time: float,
    free_dim: int,
    matrix: CharacterMatrix,
):
    def neg_ll(free: np.ndarray) -> float:
        times = ultrametric_times_from_free(tree, free, total_time)
        ll, _ = profile_loglik_over_rates(
            tree, ancestral, times, rate_bounds, matrix=matrix
        )
        if not np.isfinite(ll):
            return 1e30
        return -ll

    return neg_ll


def _coordinate_descent(
    objective,
    x0: np.ndarray,
    bounds: Sequence[Tuple[float, float]],
    max_sweeps: int = 4,
    tol: float = 1e-2,
    seed: int = 0,
) -> Tuple[np.ndarray, float]:
    """Coordinate-wise bounded descent for high-dimensional internal-age problems.

    L-BFGS finite-difference gradients cost 2·d evals per iteration (unworkable at
    d≈400); one scalar bound per coordinate is ~30 evals regardless of d. Monotone:
    a coordinate is updated only when it strictly lowers the objective.
    """
    x = x0.copy()
    f = objective(x)
    for _ in range(max_sweeps):
        improved = False
        for k in range(len(x)):
            lo, hi = bounds[k]
            xc = x.copy()

            def f1(v: float) -> float:
                xc[k] = v
                return objective(xc)

            r = minimize_scalar(f1, bounds=(lo, hi), method="bounded",
                                options={"xatol": max(1e-4, (hi - lo) * 1e-3)})
            if r.fun < f - tol:
                x = x.copy()
                x[k] = float(r.x)
                f = r.fun
                improved = True
        if not improved:
            break
    return x, f


def scalable_timing_certificate(
    tree: LineageTree,
    matrix: CharacterMatrix,
    rate_bounds: Sequence[Tuple[float, float]],
    total_time: float,
    targets: Sequence[str],
    confidence: float = 0.95,
    n_starts: int = 4,
    seed: int = 0,
    max_sweeps: int = 3,
    n_restarts: int = 3,
    grid_n: int = 16,
    grid_sweeps: int = 2,
    progress_fn=None,
) -> List[TimingCertificate]:
    """Optimize profile likelihood; then bound each target on the LR confidence set.

    Uses multi-start local optimization over internal node ages. Returns certificates
    with solver diagnostics (not point estimates alone).
    `max_sweeps`/`n_restarts`/`grid_n`/`grid_sweeps` tune large-tree effort (defaults
    match the frozen protocol); `progress_fn(target)` is called after each target bound.
    """
    rng = np.random.default_rng(seed)
    ancestral = prepare_ancestral(tree, matrix)
    internals = [n for n in tree.internal_nodes() if n != tree.root]
    d = len(internals)
    bounds = [(0.0, total_time)] * d
    neg = _pack_objective(tree, ancestral, rate_bounds, total_time, d, matrix)

    best = None
    best_free = None
    # ponytail: coordinate descent for large trees (d>24); L-BFGS finite-difference
    # gradients cost O(d) evals/iter and stall at the ~400-internal real Dryad tree.
    if d > 24:
        x0 = np.clip(
            np.array([tree.nodes[n].age for n in internals]) + 0.0, 0.0, total_time
        )
        x_best, f_best = _coordinate_descent(neg, x0, bounds, max_sweeps=max_sweeps, seed=seed)
        for i in range(n_restarts):
            xj = np.clip(x_best + rng.normal(0.0, 0.03 * total_time, size=d), 0, total_time)
            xj, fj = _coordinate_descent(neg, xj, bounds, max_sweeps=max(1, max_sweeps - 1), seed=seed + 1)
            if fj < f_best:
                x_best, f_best = xj, fj
        best, best_free = f_best, x_best
    else:
        for i in range(n_starts):
            x0 = rng.uniform(0, total_time, size=d)
            res = minimize(neg, x0, method="L-BFGS-B", bounds=bounds)
            if best is None or res.fun < best:
                best = float(res.fun)
                best_free = res.x.copy()

    # global polish for small d
    if d <= 2:
        de = differential_evolution(neg, bounds, seed=seed, polish=True, popsize=15)
        if best is None or de.fun < best:
            best = float(de.fun)
            best_free = de.x.copy()

    assert best_free is not None and best is not None
    l_max = -best
    cutoff = 0.5 * float(chi2.ppf(confidence, df=1))  # note: ppf(1-alpha)=ppf(confidence)
    # chi2.ppf(confidence, 1) with confidence=0.95 is correct for 95% CI
    threshold = l_max - cutoff

    from b01_clock.solver.exact import _target_value

    mle_times = ultrametric_times_from_free(tree, best_free, total_time)
    certs = []
    # ponytail: DE over the full free-parameter space is only tractable up to a few
    # dozen internals. For large trees (real Dryad trees have ~400 leaves) bound each
    # target with multi-start local profile optimization around the MLE, which finds
    # the LR-set extent along plausible directions; upgrade to DE if miscoverage shows.
    use_de = d <= 24
    for target in targets:
        mle_val = _target_value(tree, mle_times, target)

        def objective_min(free: np.ndarray) -> float:
            times = ultrametric_times_from_free(tree, free, total_time)
            ll, _ = profile_loglik_over_rates(
                tree, ancestral, times, rate_bounds, matrix=matrix
            )
            if ll < threshold - 1e-6:
                return 1e6 + (threshold - ll)
            return _target_value(tree, times, target)

        def objective_max(free: np.ndarray) -> float:
            return -objective_min(free)

        if use_de:
            lo_res = differential_evolution(
                objective_min, bounds, seed=seed + 1, polish=True, popsize=12, atol=1e-6
            )
            hi_res = differential_evolution(
                objective_max, bounds, seed=seed + 2, polish=True, popsize=12, atol=1e-6
            )
            lo, hi = float(lo_res.fun), float(-hi_res.fun)
            bound_mode = "differential_evolution"
        else:
            if target.startswith("age:"):
                # Profile-grid over the focus age (robust): fix the target node's age
                # at each grid value, profile the rest by coordinate descent, collect
                # feasible focus values. The penalized full-space search below flattens
                # against the T-boundary on the 64-internal induced Dryad tree and
                # emitted spurious degenerate point certificates.
                k = internals.index(target.split(":", 1)[1])
                feas = []
                xwarm = best_free.copy()
                for f in np.linspace(0.0, total_time, grid_n):
                    x0g = np.clip(xwarm, 0.0, total_time)
                    x0g[k] = f
                    xs, fs = _coordinate_descent(neg, x0g, bounds, max_sweeps=grid_sweeps, seed=seed + 3)
                    xwarm = xs
                    if -fs >= threshold - 1e-6:
                        times_g = ultrametric_times_from_free(tree, xs, total_time)
                        feas.append(_target_value(tree, times_g, target))
                if feas:
                    lo, hi = min(feas), max(feas)
                    bound_mode = "profile_grid_target"
                else:
                    lo, hi = 0.0, total_time
                    bound_mode = "profile_grid_empty_unidentified"
            else:
                starts = [best_free.copy()]
                for _ in range(3):
                    j = np.clip(best_free + rng.normal(0.0, 0.03 * total_time, size=d), 0, total_time)
                    starts.append(j)
                lo_res, lo_best = None, None
                for x0 in starts:
                    xs, fs = _coordinate_descent(objective_min, x0, bounds, max_sweeps=grid_sweeps, seed=seed + 3)
                    if lo_best is None or fs < lo_best:
                        lo_best, lo_res = fs, xs
                hi_res, hi_best = None, None
                for x0 in starts:
                    xs, fs = _coordinate_descent(objective_max, x0, bounds, max_sweeps=grid_sweeps, seed=seed + 4)
                    if hi_best is None or fs < hi_best:
                        hi_best, hi_res = fs, xs
                lo, hi = lo_best, -hi_best
                bound_mode = "multi_start_lbfgs"

        lo = float(np.clip(lo, 0.0, total_time))
        hi = float(np.clip(hi, 0.0, total_time))
        if hi < lo:
            lo, hi = hi, lo
        width = hi - lo
        if width <= 1e-6 * total_time:
            status = "point_identified"
        elif width >= total_time * 0.99:
            status = "unidentified"
        else:
            status = "partially_identified"
        certs.append(
            TimingCertificate(
                target=target,
                lower=lo,
                upper=hi,
                mle=float(mle_val),
                status=status,
                l_max=float(l_max),
                threshold=float(threshold),
                solver_status="ok",
                diagnostics={
                    "confidence": confidence,
                    "cutoff": cutoff,
                    "n_free": d,
                    "internals": internals,
                    "bound_mode": bound_mode,
                    "certificate": (
                        "Bound from profile likelihood ratio set with rates constrained "
                        "to declared bounds; absolute time not inferred beyond those bounds."
                    ),
                },
            )
        )
        if progress_fn is not None:
            progress_fn(target, certs[-1])
    return certs
