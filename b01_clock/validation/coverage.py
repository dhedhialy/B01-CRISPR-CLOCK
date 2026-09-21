"""Coverage, interval width, temporal order, false-identification metrics."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Dict, List, Optional, Sequence

import numpy as np

from b01_clock.model.exposure import exposure_from_edit_fraction
from b01_clock.model.observations import EDITED, MISSING, UNEDITED, CharacterMatrix
from b01_clock.model.tree import LineageTree, small_scout_tree
from b01_clock.simulate.generator import draw_site_rates, simulate_recorder
from b01_clock.solver.exact import exact_profile_bounds


def per_site_rate_bounds(
    tree: LineageTree,
    matrix: CharacterMatrix,
    total_time: float,
    pad: float = 1.25,
    floor: float = 1e-3,
) -> list[tuple[float, float]]:
    """Empirical per-site rate boxes from leaf edit fractions (no truth leakage).

    r̂_s = -log(1 - p̂_s) / T for each site; box = [r̂_s/pad, r̂_s*pad].
    Site-wise profiling under these boxes correctly absorbs site heterogeneity,
    restoring profile-LR calibration that the global-rate model loses.
    """
    leaves = tree.leaves()
    bounds = []
    for s in range(matrix.n_sites):
        edited = observed = 0
        for ℓ in leaves:
            st = int(matrix.leaf_states[ℓ][s])
            if st == MISSING:
                continue
            observed += 1
            if st == EDITED:
                edited += 1
        if observed == 0:
            bounds.append((floor, max(floor * 2, 10.0)))
            continue
        p = min(edited / observed, 1.0 - 1e-6)
        r_hat = exposure_from_edit_fraction(p) / max(total_time, 1e-12)
        bounds.append((max(floor, r_hat / pad), max(floor * 1.01, r_hat * pad)))
    return bounds


@dataclass
class CoverageReport:
    confidence: float
    n_reps: int
    coverage: float
    mean_width_frac: float
    false_identification_rate: float
    abstention_rate: float
    temporal_order_accuracy: float
    details: list = field(default_factory=list)


def _true_target(tree: LineageTree, target: str) -> float:
    if target.startswith("age:"):
        return tree.node_age(target.split(":", 1)[1])
    if target.startswith("branch:"):
        p, c = target.split(":", 1)[1].split("->")
        return tree.nodes[c].time_from_parent
    if target.startswith("order_gap:"):
        a, b = target.split(":", 1)[1].split(",")
        return tree.node_age(b) - tree.node_age(a)
    raise ValueError(target)


def temporal_order_correct(tree: LineageTree, a: str, b: str, lo_gap: float, hi_gap: float) -> Optional[bool]:
    """If interval for age(b)-age(a) excludes 0, check whether sign matches truth."""
    if lo_gap > 0:
        pred = 1
    elif hi_gap < 0:
        pred = -1
    else:
        return None  # abstain on order
    truth = np.sign(tree.node_age(b) - tree.node_age(a))
    if truth == 0:
        return None
    return bool(pred == truth)


def evaluate_coverage(
    n_reps: int = 40,
    n_sites: int = 48,
    total_time: float = 1.0,
    heterogeneity: float = 0.4,
    confidence: float = 0.95,
    grid_size: int = 31,
    seed: int = 0,
    target: str = "age:A",
    rate_bound_pad: float = 1.25,
    rate_mode: str = "global",
) -> CoverageReport:
    rng = np.random.default_rng(seed)
    covered = 0
    widths = []
    false_id = 0
    abstain = 0
    order_correct = 0
    order_total = 0
    details = []

    for i in range(n_reps):
        tree = small_scout_tree(total_time)
        rates = draw_site_rates(n_sites, mean=0.9, heterogeneity=heterogeneity, seed=int(rng.integers(1e9)))
        truth = simulate_recorder(tree, rates, total_time=total_time, seed=int(rng.integers(1e9)))
        if rate_mode == "per_site":
            bounds = per_site_rate_bounds(truth.tree, truth.matrix, total_time, pad=rate_bound_pad)
        else:
            lo_r, hi_r = float(rates.min() / rate_bound_pad), float(rates.max() * rate_bound_pad)
            bounds = (lo_r, hi_r)
        res = exact_profile_bounds(
            truth.tree,
            truth.matrix,
            bounds,
            total_time,
            targets=[target],
            confidence=confidence,
            grid_size=grid_size,
            rate_mode=rate_mode,
        )
        r0 = res[0]
        true_val = _true_target(truth.tree, target)
        hit = r0.lower - 1e-9 <= true_val <= r0.upper + 1e-9
        covered += int(hit)
        width_frac = (r0.upper - r0.lower) / total_time
        widths.append(width_frac)
        if r0.status == "point_identified" and width_frac < 0.02:
            if abs(true_val - r0.mle) > 0.05 * total_time:
                false_id += 1
        if r0.status == "unidentified":
            abstain += 1
        if r0.lower > 0.35 * total_time:
            order_total += 1
            order_correct += int(true_val > 0.35 * total_time)
        details.append(
            {
                "hit": hit,
                "lower": r0.lower,
                "upper": r0.upper,
                "true": true_val,
                "status": r0.status,
                "width_frac": width_frac,
                "rate_mode": rate_mode,
            }
        )

    return CoverageReport(
        confidence=confidence,
        n_reps=n_reps,
        coverage=covered / max(n_reps, 1),
        mean_width_frac=float(np.mean(widths) if widths else np.nan),
        false_identification_rate=false_id / max(n_reps, 1),
        abstention_rate=abstain / max(n_reps, 1),
        temporal_order_accuracy=order_correct / max(order_total, 1),
        details=details,
    )
