"""Exact solver and analytic cases."""

import numpy as np

from b01_clock.model.tree import small_scout_tree
from b01_clock.simulate.generator import draw_site_rates, simulate_recorder
from b01_clock.solver.exact import exact_profile_bounds, verify_against_grid


def test_exact_solver_runs_and_covers_truth():
    T = 1.0
    tree = small_scout_tree(T)
    rates = draw_site_rates(40, mean=0.8, heterogeneity=0.15, seed=3)
    truth = simulate_recorder(tree, rates, total_time=T, seed=4)
    bounds = (float(rates.min() / 1.3), float(rates.max() * 1.3))
    res = exact_profile_bounds(
        truth.tree,
        truth.matrix,
        bounds,
        T,
        targets=["age:A"],
        confidence=0.95,
        grid_size=31,
    )[0]
    true_age = truth.tree.node_age("A")
    assert res.lower <= true_age <= res.upper
    assert 0 <= res.lower <= res.upper <= T
    assert (res.upper - res.lower) < T  # material narrowing vs [0,T]


def test_verify_against_finer_grid():
    T = 1.0
    tree = small_scout_tree(T)
    rates = draw_site_rates(30, mean=0.8, heterogeneity=0.1, seed=5)
    truth = simulate_recorder(tree, rates, total_time=T, seed=6)
    bounds = (float(rates.min() / 1.2), float(rates.max() * 1.2))
    res = exact_profile_bounds(
        truth.tree, truth.matrix, bounds, T, targets=["age:A"], grid_size=21
    )[0]
    v = verify_against_grid(truth.tree, truth.matrix, bounds, T, "age:A", res, fine_grid=41, tol=0.15)
    assert v["ok"]
