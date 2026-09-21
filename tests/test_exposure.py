"""Unit tests: exposure model."""

import numpy as np

from b01_clock.model.exposure import (
    edit_prob,
    survival_prob,
    loglik_transition,
    time_bounds_from_exposure,
)


def test_edit_survival_complement():
    r, t = 1.5, 0.7
    assert abs(edit_prob(r, t) + survival_prob(r, t) - 1.0) < 1e-12


def test_loglik_impossible_reverse():
    assert loglik_transition(1, 0, 1.0, 0.5) == -np.inf


def test_loglik_unedited():
    assert abs(loglik_transition(0, 0, 2.0, 0.5) - (-1.0)) < 1e-12


def test_time_bounds_mapping():
    lo, hi = time_bounds_from_exposure(1.0, 2.0, 0.5, 2.0)
    assert abs(lo - 0.5) < 1e-12  # 1/2
    assert abs(hi - 4.0) < 1e-12  # 2/0.5
