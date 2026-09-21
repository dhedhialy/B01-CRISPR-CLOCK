"""Rate–time invariance tests."""

import numpy as np

from b01_clock.model.tree import small_scout_tree
from b01_clock.simulate.generator import draw_site_rates, simulate_recorder
from b01_clock.solver.invariance import check_rate_time_invariance, rate_time_transform_feasible


def test_likelihood_rate_time_invariance():
    T = 1.0
    tree = small_scout_tree(T)
    rates = draw_site_rates(25, mean=1.0, heterogeneity=0.0, seed=1)
    truth = simulate_recorder(tree, rates, total_time=T, seed=2)
    out = check_rate_time_invariance(
        truth.tree, truth.matrix, rates, truth.branch_times, scale=3.0
    )
    assert out["ok"]


def test_guardrail_blocks_false_absolute_time():
    T = 1.0
    tree = small_scout_tree(T)
    times = tree.branch_times_vector()
    rates = np.ones(10)
    bounds = (0.5, 2.0)
    g = rate_time_transform_feasible(times, rates, bounds, T, scale=2.0)
    assert not g["feasible_with_anchor_and_bounds"]
    assert g["rates_in_bounds"]  # 2.0 still in [0.5,2.0]
