#!/usr/bin/env python3
"""Empirical recalibration of certificate coverage on the ConvexML benchmark.

Question: the profile-LR certificate over-apes true splits (miscalibrated deep,
systematic late shift). Is the miscalibration LEARNABLE so that a recalibrated
interval restores nominal coverage? This is benchmark-as-it-should-be-used:
the deposit's newick lengths are known truth, so we can (i) measure the
residual structure and (ii) build a conformal width q from a training set of
internal nodes, then evaluate leave-one-out and cross-draw coverage.

Uses CHEAP certificates (5-grid focuses, single inner sweep) over ALL internal
nodes — a coverage/fidelity study, not the frozen headline numbers.

Writes experiments/recalibration_study.json
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scipy.stats import chi2

from b01_clock.data.adapters.convexml import load_convexml_example
from b01_clock.model.observations import CharacterMatrix
from b01_clock.solver.likelihood import (
    profile_loglik_over_rates, ultrametric_times_from_free,
)
from b01_clock.solver.exact import _target_value
from b01_clock.solver.optimize import _coordinate_descent, prepare_ancestral
from b01_clock.validation.calibration import held_out_rate_bounds, split_sites

T = 1.0
OUT = ROOT / "experiments"
OUT.mkdir(exist_ok=True)


def canonical_setup(keep_seed=13, site_seed=7):
    bundle = load_convexml_example(ROOT / "data_store" / "convexml", total_time=T, rep=0)
    tree, matrix, meta = bundle.tree, bundle.matrix, bundle.meta
    calib, analysis = split_sites(matrix.n_sites, calib_frac=0.4, seed=site_seed)
    rb = held_out_rate_bounds(tree, matrix, calib, T)
    ana = matrix.site_slice(analysis)
    rng = np.random.default_rng(keep_seed)
    keep = sorted(rng.choice(sorted(tree.leaves()), size=64, replace=False).tolist())
    ctree = tree.induced_subtree(keep)
    ana_mat = CharacterMatrix({l: ana.leaf_states[l] for l in keep}, ana.n_sites)
    return ctree, ana_mat, rb


def cheap_cert_all(ctree, matrix, rate_bounds):
    ancestral = prepare_ancestral(ctree, matrix)
    internals = [n for n in ctree.internal_nodes() if n != ctree.root]
    d = len(internals)
    bounds = [(0.0, T)] * d

    def neg(free):
        times = ultrametric_times_from_free(ctree, free, T)
        ll, _ = profile_loglik_over_rates(ctree, ancestral, times, rate_bounds,
                                          matrix=matrix, rate_mode="global")
        return -ll if np.isfinite(ll) else 1e30

    x0 = np.clip(np.array([ctree.nodes[n].age for n in internals]) + 0.0, 0.0, T)
    xb, fb = _coordinate_descent(neg, x0, bounds, max_sweeps=2, seed=7)
    l_max = -fb
    threshold = l_max - 0.5 * float(chi2.ppf(0.95, df=1))
    mle_times = ultrametric_times_from_free(ctree, xb, T)

    rows = []
    for u in internals:
        true = ctree.node_age(u)
        mle_val = _target_value(ctree, mle_times, f"age:{u}")
        k = internals.index(u)
        feas = []
        xwarm = xb.copy()
        for f in np.linspace(0.0, T, 5):
            x0g = np.clip(xwarm, 0.0, T)
            x0g[k] = f
            xs, fs = _coordinate_descent(neg, x0g, bounds, max_sweeps=1, seed=7)
            xwarm = xs
            if -fs >= threshold - 1e-6:
                feas.append(_target_value(ctree, ultrametric_times_from_free(ctree, xs, T), f"age:{u}"))
        if feas:
            lo, hi = float(min(feas)), float(max(feas))
        else:
            lo, hi = 0.0, T
        lo, hi = max(0.0, lo), min(T, hi)
        rows.append({
            "node": u, "true": float(true), "mle": float(mle_val), "lo": lo, "hi": hi,
            "profile_covered": bool(lo <= true <= hi),
            "resid": float(max(true - lo, hi - true)),  # width margin to cover
        })
        if (len(rows) % 5) == 0:
            print(f"  [{len(rows)} nodes done] mle={mle_val:.3f} true={true:.3f}", flush=True)
    return rows


def conformal_loo(rows, alpha=0.05):
    n = len(rows)
    r = np.array([x["resid"] for x in rows])
    t = np.array([x["true"] for x in rows])
    m = np.array([x["mle"] for x in rows])
    q_idx = int(np.ceil((n + 1) * (1 - alpha)))
    q_idx = min(q_idx, n)
    covered = []
    for i in range(n):
        r_train = np.sort(np.delete(r, i))
        q = r_train[q_idx - 1]
        covered.append(bool(m[i] - q <= t[i] <= m[i] + q))
    mean_w = float(np.mean(2 * np.sort(r)[q_idx - 1]))
    return float(np.mean(covered)), mean_w


def conformal_apply(train_rows, test_rows, alpha=0.05):
    r = np.sort([x["resid"] for x in train_rows])
    n = len(r)
    q_idx = min(int(np.ceil((n + 1) * (1 - alpha))), n)
    q = r[q_idx - 1]
    cov = 0
    for x in test_rows:
        if x["mle"] - q <= x["true"] <= x["mle"] + q:
            cov += 1
    return cov / max(1, len(test_rows)), 2 * q


def main():
    res = {}
    for name, seed in [("canonical_draw13", 13), ("alt_draw42", 42)]:
        ctree, matrix, rb = canonical_setup(keep_seed=seed)
        print(f"{name}: computing cheap certs over {len([n for n in ctree.internal_nodes() if n != ctree.root])} internals ...")
        rows = cheap_cert_all(ctree, matrix, rb)
        prof_cov = float(np.mean([x["profile_covered"] for x in rows]))
        loo_cov, loo_w = conformal_loo(rows)
        res[name] = {
            "n_nodes": len(rows),
            "profile_coverage": prof_cov,
            "conformal_loo_coverage": loo_cov,
            "conformal_width": loo_w,
            "mean_abs_bias": float(np.mean([abs(x["mle"] - x["true"]) for x in rows])),
            "mean_bias": float(np.mean([x["mle"] - x["true"] for x in rows])),
            "rows": rows,
        }
        print(f"  profile cov={prof_cov:.3f} | conformal LOOCV cov={loo_cov:.3f} width={loo_w:.3f}")

    # cross-draw: fit on one draw, apply to the other (external-benchmark transfer)
    a, b = res["canonical_draw13"], res["alt_draw42"]
    for fname, tname, train, test in [("13->42", "42", a, b), ("42->13", "13", b, a)]:
        cov, w = conformal_apply(train["rows"], test["rows"])
        res[f"transfer_{fname}"] = {"coverage": cov, "width": w}
        print(f"transfer {fname} applied from {train['n_nodes']} nodes: cov={cov:.3f} width={w:.3f}")

    (OUT / "recalibration_study.json").write_text(json.dumps(res, indent=2, default=str))
    print("wrote experiments/recalibration_study.json")


if __name__ == "__main__":
    main()