"""Per-site rate profiling restores coverage under site heterogeneity + routed battery."""

import numpy as np

from b01_clock.model.tree import small_scout_tree
from b01_clock.simulate.adversarial import (
    AdversarialSpec,
    routed_headline_evaluate,
    simulate_adversarial,
)
from b01_clock.simulate.generator import draw_site_rates, simulate_recorder
from b01_clock.solver.exact import exact_profile_bounds
from b01_clock.validation.coverage import evaluate_coverage, per_site_rate_bounds


def test_per_site_bounds_positive_and_ordered():
    tree = small_scout_tree(1.0)
    rates = draw_site_rates(10, mean=0.9, heterogeneity=0.8, seed=3)
    truth = simulate_recorder(tree, rates, total_time=1.0, seed=4)
    bounds = per_site_rate_bounds(truth.tree, truth.matrix, 1.0, pad=1.25)
    assert len(bounds) == 10
    for lo, hi in bounds:
        assert 0 < lo <= hi


def test_per_site_solver_covers_strong_heterogeneity():
    truth, _ = simulate_adversarial(
        AdversarialSpec("site_het_strong", site_heterogeneity=1.2), n_sites=12, seed=0
    )
    bounds = per_site_rate_bounds(truth.tree, truth.matrix, truth.total_time, pad=1.25)
    r = exact_profile_bounds(
        truth.tree, truth.matrix, bounds, truth.total_time,
        targets=["age:A"], grid_size=7, rate_mode="per_site",
    )[0]
    true_age = truth.tree.node_age("A")
    assert r.lower <= true_age <= r.upper


def test_routed_evaluate_topology_error_covers():
    truth, meta = simulate_adversarial(
        AdversarialSpec("topology_error", topology_perturb=True), n_sites=16, seed=17
    )
    out = routed_headline_evaluate(truth, meta["assumed_bounds"], grid_size=7)
    assert out["covered"]


def test_coverage_reports_nominal_under_per_site_mode():
    c = evaluate_coverage(
        n_reps=4, n_sites=10, heterogeneity=0.5, grid_size=7, seed=5, rate_mode="per_site"
    )
    assert c.coverage >= 0.5  # small n, but should not crater to 0
    assert not np.isnan(c.mean_width_frac)