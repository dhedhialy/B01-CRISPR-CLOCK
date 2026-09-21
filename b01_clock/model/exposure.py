"""Irreversible exposure model.

Site s on branch e:
    p_se = 1 - exp(-r_s * t_e)

Identified quantity is exposure λ = r * t. Absolute time is not identified
without rate constraints / anchors; this module never collapses that distinction.
"""

from __future__ import annotations

import numpy as np

EPS = 1e-12


def edit_prob(rate: float | np.ndarray, time: float | np.ndarray) -> np.ndarray:
    """p = 1 - exp(-r * t)."""
    r = np.asarray(rate, dtype=float)
    t = np.asarray(time, dtype=float)
    return 1.0 - np.exp(-r * t)


def survival_prob(rate: float | np.ndarray, time: float | np.ndarray) -> np.ndarray:
    """Probability of remaining unedited: exp(-r * t)."""
    r = np.asarray(rate, dtype=float)
    t = np.asarray(time, dtype=float)
    return np.exp(-r * t)


def path_exposure(rate: float, branch_times: np.ndarray) -> float:
    """Total exposure along a path with homogeneous site rate."""
    return float(rate * np.sum(branch_times))


def loglik_transition(
    parent_state: int,
    child_state: int,
    rate: float,
    time: float,
) -> float:
    """Log-likelihood of an irreversible 0/1 transition on one branch.

    States: 0 = unedited, 1 = edited, -1 = missing (ignored → 0 contrib).
    """
    if parent_state < 0 or child_state < 0:
        return 0.0
    if parent_state == 1 and child_state == 0:
        return -np.inf
    if parent_state == 1 and child_state == 1:
        return 0.0
    rt = max(rate, 0.0) * max(time, 0.0)
    if parent_state == 0 and child_state == 0:
        return -rt
    if parent_state == 0 and child_state == 1:
        # log(1 - exp(-rt)); stable for small and large rt
        if rt < 1e-8:
            return np.log(max(rt, EPS))
        if rt > 30:
            return 0.0
        return float(np.log1p(-np.exp(-rt)))
    raise ValueError(f"Invalid states parent={parent_state} child={child_state}")


def exposure_from_edit_fraction(p_hat: float) -> float:
    """MLE exposure from edit probability: -log(1 - p)."""
    p = float(np.clip(p_hat, 0.0, 1.0 - EPS))
    return float(-np.log1p(-p))


def time_bounds_from_exposure(
    exposure_lo: float,
    exposure_hi: float,
    rate_lo: float,
    rate_hi: float,
) -> tuple[float, float]:
    """Map exposure interval through rate bounds to a time interval.

    t = λ / r. With r ∈ [r_lo, r_hi] and λ ∈ [λ_lo, λ_hi]:
        t_lo = λ_lo / r_hi,  t_hi = λ_hi / r_lo
    (when bounds are positive). Preserves rate×time vs absolute-time distinction.
    """
    if rate_lo <= 0 or rate_hi <= 0:
        raise ValueError("Rate bounds must be strictly positive for time mapping")
    if rate_lo > rate_hi:
        raise ValueError("rate_lo > rate_hi")
    t_lo = exposure_lo / rate_hi
    t_hi = exposure_hi / rate_lo if rate_lo > 0 else np.inf
    return float(t_lo), float(t_hi)
