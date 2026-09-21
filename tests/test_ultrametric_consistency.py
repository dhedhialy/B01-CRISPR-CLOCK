"""Regression test: branch times re-accumulate to ages in [0, T] on any tree.

guards the ctree->copy() edges()-ordering bug (branch-time vector mis-paired
between an original tree and its BFS-rebuilt copy, producing reported node
ages > T) and the single-pass bump that raised a parent after its child's
branch time was committed.
"""

import numpy as np

from b01_clock.model.tree import small_scout_tree, LineageTree
from b01_clock.solver.likelihood import ultrametric_times_from_free


def _ages(tree, times):
    copy = tree.copy()
    copy.set_branch_times(times)
    return {n: copy.nodes[n].age for n in copy.nodes}


def test_ages_stay_within_total_time_on_deep_chain():
    T = 1.0
    tree = LineageTree("root")
    parent = "root"
    n = 0
    for depth in range(1, 12):
        node = f"n{depth}"
        tree.add_child(parent, node, 0.0)
        parent = node
    for depth, name in enumerate(["L0", "L1"]):
        tree.add_child(parent, name, 1.0)
        n += 1
    internals = [node for node in tree.internal_nodes() if node != tree.root]
    free = np.ones(len(internals)) * 0.5
    times = ultrametric_times_from_free(tree, free, T)
    ages = _ages(tree, times)
    assert max(ages.values()) <= T + 1e-9
    assert tree.copy().set_branch_times(times) or tree.is_ultrametric()
    assert max(ages.values()) == T  # all leaves at harvest time


def test_edge_order_identical_between_original_and_copy():
    T = 1.0
    tree = small_scout_tree(T)
    tree.add_child("A", "L5", T - 0.5)  # asymmetry beyond the 4-leaf default
    copy = tree.copy()
    assert tree.edges() == copy.edges()
    assert list(tree.edges()) == sorted(tree.edges(), key=lambda e: e[1])