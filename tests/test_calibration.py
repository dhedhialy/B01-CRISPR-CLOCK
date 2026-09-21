"""Calibration split and bound transfer."""

from b01_clock.model.tree import small_scout_tree
from b01_clock.simulate.generator import draw_site_rates, simulate_recorder
from b01_clock.validation.calibration import held_out_rate_bounds, split_sites


def test_split_disjoint():
    calib, analysis = split_sites(50, calib_frac=0.4, seed=0)
    assert set(calib).isdisjoint(set(analysis))
    assert len(calib) + len(analysis) == 50


def test_held_out_bounds_positive():
    T = 1.0
    tree = small_scout_tree(T)
    rates = draw_site_rates(40, mean=1.0, heterogeneity=0.3, seed=8)
    truth = simulate_recorder(tree, rates, total_time=T, seed=9)
    calib, _ = split_sites(40, seed=10)
    bounds = held_out_rate_bounds(truth.tree, truth.matrix, calib, T)
    lo, hi = bounds
    assert 0 < lo < hi
