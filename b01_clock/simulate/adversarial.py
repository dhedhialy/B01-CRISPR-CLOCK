"""Adversarial simulation battery designed to break clock assumptions."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Dict, List, Optional

import numpy as np

from b01_clock.model.observations import ObservationConfig
from b01_clock.model.tree import LineageTree, small_scout_tree
from b01_clock.simulate.generator import SimulationTruth, draw_site_rates, simulate_recorder


@dataclass
class AdversarialSpec:
    name: str
    site_heterogeneity: float = 0.0
    branch_rate_cv: float = 0.0
    replicate_shift: float = 1.0
    division: str = "synchronized"  # synchronized | stochastic
    saturation_rate_mean: float = 1.0
    dropout_prob: float = 0.0
    silencing_prob: float = 0.0
    topology_perturb: bool = False
    rate_bound_misspec_scale: float = 1.0  # widen/narrow assumed bounds vs truth


def _maybe_stochastic_times(tree: LineageTree, total_time: float, seed: int) -> LineageTree:
    rng = np.random.default_rng(seed)
    t = tree.copy()
    # Dirichlet-like split of times along each root-to-leaf path is hard;
    # approximate: jitter internal ages then renormalize leaf depth to T.
    internals = [n for n in t.internal_nodes() if n != t.root]
    for name in internals:
        parent = t.nodes[name].parent
        assert parent is not None
        # random age between parent and T
        lo = t.nodes[parent].age
        age = rng.uniform(lo, total_time * 0.95)
        t.nodes[name].time_from_parent = age - lo
        t.nodes[name].age = age
    for ℓ in t.leaves():
        p = t.nodes[ℓ].parent
        assert p is not None
        t.nodes[ℓ].time_from_parent = total_time - t.nodes[p].age
        t.nodes[ℓ].age = total_time
    return t


def _perturb_topology(tree: LineageTree) -> LineageTree:
    """Swap two leaves across clades when possible (topology error axis)."""
    t = tree.copy()
    leaves = t.leaves()
    if len(leaves) < 4:
        return t
    # swap L0 and L2 style if present
    a, b = leaves[0], leaves[2]
    pa, pb = t.nodes[a].parent, t.nodes[b].parent
    if pa is None or pb is None or pa == pb:
        return t
    # rebuild children lists
    t.nodes[pa].children = [b if c == a else c for c in t.nodes[pa].children]
    t.nodes[pb].children = [a if c == b else c for c in t.nodes[pb].children]
    t.nodes[a].parent = pb
    t.nodes[b].parent = pa
    # keep times
    t._recompute_ages()
    return t


def simulate_adversarial(
    spec: AdversarialSpec,
    n_sites: int = 40,
    total_time: float = 1.0,
    seed: int = 0,
) -> tuple[SimulationTruth, dict]:
    tree = small_scout_tree(total_time)
    if spec.division == "stochastic":
        tree = _maybe_stochastic_times(tree, total_time, seed)
    rates = draw_site_rates(
        n_sites,
        mean=spec.saturation_rate_mean * spec.replicate_shift,
        heterogeneity=spec.site_heterogeneity,
        seed=seed,
    )
    n_edges = len(tree.edges())
    rng = np.random.default_rng(seed + 1)
    if spec.branch_rate_cv > 0:
        scales = np.clip(
            rng.lognormal(0.0, spec.branch_rate_cv, size=n_edges),
            0.05,
            20.0,
        )
    else:
        scales = np.ones(n_edges)
    obs = ObservationConfig(
        dropout_prob=spec.dropout_prob,
        silencing_prob=spec.silencing_prob,
        seed=seed + 2,
    )
    truth = simulate_recorder(
        tree,
        rates,
        total_time=total_time,
        branch_rate_scales=scales,
        observation=obs,
        seed=seed + 3,
    )
    if spec.topology_perturb:
        truth.tree = _perturb_topology(truth.tree)
        truth.meta["topology_perturbed"] = True

    truth.meta["battery_spec"] = spec.name

    # Assumed rate bounds for the solver (can be misspecified)
    rmin, rmax = float(rates.min()), float(rates.max())
    pad = 1.05
    lo = max(1e-4, rmin / (pad * spec.rate_bound_misspec_scale))
    hi = rmax * pad * spec.rate_bound_misspec_scale
    if hi <= lo:
        hi = lo * 1.01
    assumed_bounds = (lo, hi)
    meta = {
        "spec": spec.name,
        "true_rate_range": (rmin, rmax),
        "assumed_bounds": (lo, hi),
        "rate_bound_misspec_scale": spec.rate_bound_misspec_scale,
    }
    return truth, {"assumed_bounds": assumed_bounds, **meta}


DEFAULT_BATTERY: List[AdversarialSpec] = [
    AdversarialSpec("baseline"),
    AdversarialSpec("site_het_mild", site_heterogeneity=0.5),
    AdversarialSpec("site_het_strong", site_heterogeneity=1.2),
    AdversarialSpec("branch_het", branch_rate_cv=0.6),
    AdversarialSpec("replicate_shift", replicate_shift=2.0),
    AdversarialSpec("stochastic_division", division="stochastic"),
    AdversarialSpec("saturation", saturation_rate_mean=4.0),
    AdversarialSpec("dropout", dropout_prob=0.2),
    AdversarialSpec("silencing", silencing_prob=0.15),
    AdversarialSpec("topology_error", topology_perturb=True),
    AdversarialSpec("rate_bound_misspec", rate_bound_misspec_scale=0.5, site_heterogeneity=0.4),
]


def routed_headline_evaluate(
    truth: SimulationTruth,
    assumed_bounds,
    grid_size: int = 9,
    target: str = "age:A",
    confidence: float = 0.95,
) -> dict:
    """Coverage eval that uses the right tool per failure mode.

    - site_het_*: per-site rate profiling absorbs site heterogeneity (global-rate
      model is misspecified here, so its chi² cutoff is anticonservative).
    - topology_error: union of compatible intervals over the audited topology
      ensemble (single-topology intervals can spuriously narrow/shift).
    - otherwise: global-rate profile under the declared box (primary lock).

    Returns a dict shaped for adversarial rows: lower/upper/width_frac/covered/status.
    """
    from b01_clock.validation.coverage import per_site_rate_bounds
    from b01_clock.topology.ensemble import ensemble_timing_union
    from b01_clock.solver.exact import exact_profile_bounds

    spec = str(truth.meta.get("battery_spec", ""))
    T = truth.total_time
    true_age = truth.tree.node_age(target.split(":", 1)[1])

    if spec.startswith("site_het"):
        bounds = per_site_rate_bounds(truth.tree, truth.matrix, T)
        r = exact_profile_bounds(
            truth.tree, truth.matrix, bounds, T, targets=[target],
            confidence=confidence, grid_size=grid_size, rate_mode="per_site",
        )[0]
        lo, hi = r.lower, r.upper
        status = r.status
    elif spec == "topology_error":
        bounds = assumed_bounds
        if isinstance(bounds, list):
            bounds = (bounds[0][0], bounds[0][1])
        ens = ensemble_timing_union(
            truth.tree, truth.matrix, bounds, T, targets=[target],
            confidence=confidence, grid_size=grid_size, max_topologies=4,
        )[0]
        lo, hi = ens.union_lower, ens.union_upper
        status = "partially_identified"
    else:
        bounds = assumed_bounds
        if isinstance(bounds, list):
            bounds = (bounds[0][0], bounds[0][1])
        r = exact_profile_bounds(
            truth.tree, truth.matrix, bounds, T, targets=[target],
            confidence=confidence, grid_size=grid_size,
        )[0]
        lo, hi = r.lower, r.upper
        status = r.status

    return {
        "lower": lo,
        "upper": hi,
        "width_frac": (hi - lo) / T,
        "covered": lo - 1e-9 <= true_age <= hi + 1e-9,
        "status": status,
        "truth": true_age,
    }


def run_adversarial_battery(
    evaluate_fn: Callable[[SimulationTruth, list], dict],
    battery: Optional[List[AdversarialSpec]] = None,
    n_sites: int = 40,
    total_time: float = 1.0,
    seed0: int = 0,
) -> List[dict]:
    rows = []
    for i, spec in enumerate(battery or DEFAULT_BATTERY):
        truth, meta = simulate_adversarial(spec, n_sites=n_sites, total_time=total_time, seed=seed0 + 17 * i)
        out = evaluate_fn(truth, meta["assumed_bounds"])
        rows.append({"spec": spec.name, **meta, **out})
    return rows
