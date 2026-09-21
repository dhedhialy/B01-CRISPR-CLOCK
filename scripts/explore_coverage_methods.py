#!/usr/bin/env python3
"""Exploration harness: can calibration be restored on the ConvexML benchmark?

Same canonical tree / site split / seeds as validate_convexml_ground_truth.py;
only the *method* varies. Cheap fidelity (fewer profile-grid focuses, single
inner sweep) to rank variants; winners get re-run at full fidelity.

Writes experiments/coverage_methods_exploration.json
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from b01_clock.data.adapters.convexml import load_convexml_example
from b01_clock.model.observations import CharacterMatrix, EDITED, UNEDITED, MISSING
from b01_clock.solver.exact import _target_value
from b01_clock.solver.likelihood import ultrametric_times_from_free
from b01_clock.solver.optimize import _coordinate_descent, prepare_ancestral
from b01_clock.solver.likelihood import profile_loglik_over_rates
from b01_clock.validation.calibration import held_out_rate_bounds, split_sites

from scipy.stats import chi2
from scipy.optimize import minimize_scalar

T = 1.0
OUT = ROOT / "experiments"
OUT.mkdir(exist_ok=True)
RNG = np.random.default_rng(13)
N_FOCUS = 8


def canonical_setup(rep: int = 0, keep_seed: int = 13, site_seed: int = 7):
    bundle = load_convexml_example(ROOT / "data_store" / "convexml", total_time=T, rep=rep)
    tree, matrix, meta = bundle.tree, bundle.matrix, bundle.meta
    calib, analysis = split_sites(matrix.n_sites, calib_frac=0.4, seed=site_seed)
    bounds = held_out_rate_bounds(tree, matrix, calib, T)
    ana = matrix.site_slice(analysis)
    rng = np.random.default_rng(keep_seed)
    keep = sorted(rng.choice(sorted(tree.leaves()), size=64, replace=False).tolist())
    ctree = tree.induced_subtree(keep)
    ana_mat = CharacterMatrix({l: ana.leaf_states[l] for l in keep}, ana.n_sites)
    return ctree, ana_mat, bounds


def target_nodes(ctree):
    cand = [u for u in ctree.internal_nodes() if u != ctree.root]
    ranked = sorted((ctree.node_age(u), u) for u in cand)
    n = len(ranked)
    qmap = {0.15: "shallow", 0.60: "mid", 0.90: "deep"}
    picks = {}
    for q, label in qmap.items():
        u = ranked[min(n - 1, int(q * n))][1]
        picks[label] = u
    u0, a0, b0 = ctree.basal_split()
    picks["basal"] = sorted([a0, b0])[0]
    return [picks["basal"], picks["shallow"], picks["deep"]]


def as_missing_is_unedited(matrix: CharacterMatrix) -> CharacterMatrix:
    states = {l: np.where(v == MISSING, UNEDITED, v).astype(int) for l, v in matrix.leaf_states.items()}
    m = CharacterMatrix(states, matrix.n_sites)
    m.ancestral = None
    return m


def per_site_bounds_from_calib(ctree, matrix, calib_sites):
    from b01_clock.model.exposure import exposure_from_edit_fraction
    lo, hi = [], []
    for s in calib_sites:
        edited = observed = 0
        for l in ctree.leaves():
            st = int(matrix.leaf_states[l][s])
            if st == MISSING:
                continue
            observed += 1
            if st == EDITED:
                edited += 1
        if observed == 0:
            lo.append(1e-4); hi.append(20.0); continue
        p = min(edited / observed, 1 - 1e-6)
        rhat = exposure_from_edit_fraction(p)
        lo.append(max(1e-4, 0.5 * rhat))
        hi.append(min(100.0, 2.0 * rhat))
    return list(zip(lo, hi))


def cert_fast(ctree, matrix, rate_bounds, targets, *, rate_mode="global",
              cutoff_mult=1.0):
    """Cheap-profile-grid certificate; returns (coverage, rows, notes)."""
    ancestral = prepare_ancestral(ctree, matrix)
    internals = [n for n in ctree.internal_nodes() if n != ctree.root]
    d = len(internals)
    bounds = [(0.0, T)] * d

    def neg(free):
        times = ultrametric_times_from_free(ctree, free, T)
        ll, _ = profile_loglik_over_rates(ctree, ancestral, times, rate_bounds,
                                          matrix=matrix, rate_mode=rate_mode)
        if not np.isfinite(ll):
            return 1e30
        return -ll

    x0 = np.clip(np.array([ctree.nodes[n].age for n in internals]) + 0.0, 0.0, T)
    xb, fb = _coordinate_descent(neg, x0, bounds, max_sweeps=2, seed=7)
    l_max = -fb
    cutoff = cutoff_mult * 0.5 * float(chi2.ppf(0.95, df=1))
    threshold = l_max - cutoff
    mle_times = ultrametric_times_from_free(ctree, xb, T)

    rows = []
    for target in targets:
        mle_val = _target_value(ctree, mle_times, target)
        if target.startswith("age:"):
            k = internals.index(target.split(":", 1)[1])
            feas = []
            xwarm = xb.copy()
            n_grid = N_FOCUS if rate_mode == "global" else 5
            for f in np.linspace(0.0, T, n_grid):
                x0g = np.clip(xwarm, 0.0, T)
                x0g[k] = f
                xs, fs = _coordinate_descent(neg, x0g, bounds, max_sweeps=1, seed=7)
                xwarm = xs
                if -fs >= threshold - 1e-6:
                    feas.append(_target_value(ctree, ultrametric_times_from_free(ctree, xs, T), target))
            if feas:
                lo, hi = min(feas), max(feas)
            else:
                lo, hi = 0.0, T
        else:
            lo, hi = 0.0, T
        lo = float(np.clip(lo, 0.0, T)); hi = float(np.clip(hi, 0.0, T))
        if hi < lo: lo, hi = hi, lo
        rows.append({"target": target, "true_age": ctree.node_age(target.split(':',1)[1]),
                     "mle": float(mle_val), "lower": lo, "upper": hi,
                     "covered": bool(lo <= ctree.node_age(target.split(':',1)[1]) <= hi)})
    covered = sum(r["covered"] for r in rows)
    return covered / len(rows), rows


def _is_below(tree, node, anc):
    while node is not None:
        if node == anc:
            return True
        node = tree.nodes[node].parent
    return False


def matrix_loglik_two_rate(tree, leaf_array, leaf_order, rate_of_edge):
    """Vectorized pruning with an explicit per-edge rate lookup."""
    n_sites = leaf_array.shape[1]
    partial = {}
    leaf_index = {l: i for i, l in enumerate(leaf_order)}

    def postorder(u):
        if u in partial:
            return partial[u]
        if tree.nodes[u].is_leaf or not tree.nodes[u].children:
            row = leaf_array[leaf_index[u]]
            out = np.zeros((2, n_sites))
            un = row == UNEDITED
            ed = row == EDITED
            out[1, un] = -np.inf
            out[0, ed] = -np.inf
            partial[u] = out
            return out
        acc = np.zeros((2, n_sites))
        for v in tree.nodes[u].children:
            child = postorder(v)
            rt = max(rate_of_edge[(u, v)], 0.0) * max(tree.nodes[v].time_from_parent, 0.0)
            if rt < 1e-8:
                edit_log = np.log(np.maximum(rt, 1e-12))
            elif rt > 30:
                edit_log = np.zeros(n_sites)
            else:
                edit_log = np.log1p(-np.exp(-rt))
            t00 = -rt + child[0]
            t01 = edit_log + child[1]
            acc[0] += np.logaddexp(t00, t01)
            acc[1] += child[1]
        partial[u] = acc
        return acc

    return float(np.sum(postorder(tree.root)[0]))


def loglik_two_rate(ctree, matrix, times, r_a, r_b, basal_children):
    """r_a for edges in clade basal_children[0], r_b for basal_children[1]; shared elsewhere."""
    tree = ctree.copy()
    tree.set_branch_times(times)
    arr = np.stack([matrix.leaf_states[l] for l in tree.leaves()], axis=0)
    lookup = {}
    for (p, c) in tree.edges():
        for side, r in zip(basal_children, (r_a, r_b)):
            if _is_below(tree, c, side):
                lookup[(p, c)] = r
                break
        else:
            lookup[(p, c)] = 0.5 * (r_a + r_b)
    return matrix_loglik_two_rate(tree, arr, tree.leaves(), lookup)


def profile_two_rates(ctree, matrix, times, rate_bounds):
    """maximize over r_a,r_b by coordinate ascent (fast, prunes per eval)."""
    lo, hi = rate_bounds
    basal_children = ctree.nodes[ctree.root].children
    if len(basal_children) < 2:
        u = ctree.nodes[ctree.root].children[0]
        while len(ctree.nodes[u].children) < 2:
            u = ctree.nodes[u].children[0]
        basal_children = ctree.nodes[u].children

    def f1(r_a):
        def tail(r_b):
            return loglik_two_rate(ctree, matrix, times, r_a, r_b, basal_children)
        res = minimize_scalar(lambda x: -tail(x), bounds=(lo, hi), method="bounded", options={"xatol": 1e-4})
        return -res.fun, float(res.x)
    ra = 0.5 * (lo + hi)
    for _ in range(3):
        ll, rb = f1(ra)
        def f2(ra2):
            return loglik_two_rate(ctree, matrix, times, ra2, rb, basal_children)
        res = minimize_scalar(lambda x: -f2(x), bounds=(lo, hi), method="bounded", options={"xatol": 1e-4})
        ra = float(res.x); ll = -res.fun
    br = (ra, rb)
    return ll, br, basal_children


def cert_two_rate(ctree, matrix, rate_bounds, targets):
    internals = [n for n in ctree.internal_nodes() if n != ctree.root]
    d = len(internals)
    bounds = [(0.0, T)] * d
    basal_children = None

    def neg(free):
        times = ultrametric_times_from_free(ctree, free, T)
        best, _, _ = profile_two_rates(ctree, matrix, times, rate_bounds)
        return -best if np.isfinite(best) else 1e30

    x0 = np.clip(np.array([ctree.nodes[n].age for n in internals]) + 0.0, 0.0, T)
    xb, fb = _coordinate_descent(neg, x0, bounds, max_sweeps=2, seed=7)
    l_max = -fb
    threshold = l_max - 0.5 * float(chi2.ppf(0.95, df=1))
    mle_times = ultrametric_times_from_free(ctree, xb, T)
    rows = []
    for target in targets:
        mle_val = _target_value(ctree, mle_times, target)
        k = internals.index(target.split(":", 1)[1])
        feas = []
        xwarm = xb.copy()
        for f in np.linspace(0.0, T, 5):
            x0g = np.clip(xwarm, 0.0, T)
            x0g[k] = f
            xs, fs = _coordinate_descent(neg, x0g, bounds, max_sweeps=1, seed=7)
            xwarm = xs
            if -fs >= threshold - 1e-6:
                feas.append(_target_value(ctree, ultrametric_times_from_free(ctree, xs, T), target))
        if feas:
            lo, hi = min(feas), max(feas)
        else:
            lo, hi = 0.0, T
        lo, hi = float(np.clip(lo, 0, T)), float(np.clip(hi, 0, T))
        rows.append({"target": target, "true_age": ctree.node_age(target.split(':',1)[1]),
                     "mle": float(mle_val), "lower": lo, "upper": hi,
                     "covered": bool(lo <= ctree.node_age(target.split(':',1)[1]) <= hi)})
    return sum(r["covered"] for r in rows) / len(rows), rows


def main():
    ctree, matrix, bounds = canonical_setup()
    matrix_uni = as_missing_is_unedited(matrix)
    picks = target_nodes(ctree)
    targets = [f"age:{u}" for u in picks]

    print(f"bounds canonical={tuple(round(x,3) for x in bounds)}  targets={targets}")
    results = {}

    def run(name, mat=None, rb=None, mode="global", mult=1.0):
        mat = matrix if mat is None else mat
        rb = bounds if rb is None else rb
        cov, rows = cert_fast(ctree, mat, rb, targets, rate_mode=mode, cutoff_mult=mult)
        results[name] = {"coverage": cov, "n": len(rows),
                         "mean_width": float(np.mean([r["upper"]-r["lower"] for r in rows])),
                         "mean_shift": float(np.mean([r["mle"]-r["true_age"] for r in rows])),
                         "rows": rows}
        print(f"{name:36s} cov={cov:.2f}  mean_width={results[name]['mean_width']:.2f}  mean_shift={results[name]['mean_shift']:+.2f}")

    run("A_baseline_global_canonical", rb=bounds)
    run("B_global_WIDE_box", rb=(0.001, 20.0))
    run("C_global_lowfloor_box", rb=(0.001, bounds[1] * 2))
    calib, _ = split_sites(matrix.n_sites, calib_frac=0.4, seed=7)
    rb_site = per_site_bounds_from_calib(ctree, matrix, calib)
    run("D_per_site_boxes", mode="per_site", rb=rb_site)
    run("E_missing_as_unedited_global", mat=matrix_uni)
    run("F_missing_as_unedited_per_site", mat=matrix_uni, mode="per_site", rb=rb_site)
    for m in (2.0, 4.0):
        run(f"G_baseline_cutoff_mult{m:.0f}", mult=m)
    ctree2, matrix2, bounds2 = canonical_setup(keep_seed=42)
    picks2 = target_nodes(ctree2)
    targets2 = [f"age:{u}" for u in picks2]

    def run_on(name, ct, mat, rb, tgts, mode="global", mult=1.0):
        cov, rows = cert_fast(ct, mat, rb, tgts, rate_mode=mode, cutoff_mult=mult)
        results[name] = {"coverage": cov, "n": len(rows),
                         "mean_width": float(np.mean([r["upper"]-r["lower"] for r in rows])),
                         "mean_shift": float(np.mean([r["mle"]-r["true_age"] for r in rows])),
                         "rows": rows}
        print(f"{name:36s} cov={cov:.2f}  mean_width={results[name]['mean_width']:.2f}  mean_shift={results[name]['mean_shift']:+.2f}")

    run_on("H_alt_draw_seed42_global", ctree2, matrix2, bounds2, targets2)

    cov, rows = cert_two_rate(ctree, matrix, bounds, [f"age:{u}" for u in picks])
    results["J_two_rate_clades"] = {"coverage": cov, "n": len(rows),
                                    "mean_width": float(np.mean([r["upper"]-r["lower"] for r in rows])),
                                    "mean_shift": float(np.mean([r["mle"]-r["true_age"] for r in rows])),
                                    "rows": rows}
    print(f"{'J_two_rate_clades':36s} cov={cov:.2f}  mean_width={results['J_two_rate_clades']['mean_width']:.2f}  "
          f"mean_shift={results['J_two_rate_clades']['mean_shift']:+.2f}")

    (OUT / "coverage_methods_exploration.json").write_text(json.dumps(results, indent=2, default=str))
    print("wrote experiments/coverage_methods_exploration.json")


if __name__ == "__main__":
    main()