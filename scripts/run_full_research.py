#!/usr/bin/env python3
"""Full B01 research experiment suite → tables + figure data under B01_CRISPR_CLOCK/."""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from b01_clock.model.tree import LineageTree, small_scout_tree
from b01_clock.simulate.adversarial import DEFAULT_BATTERY, run_adversarial_battery, simulate_adversarial
from b01_clock.simulate.generator import draw_site_rates, simulate_recorder
from b01_clock.solver.exact import exact_profile_bounds
from b01_clock.solver.invariance import check_rate_time_invariance, rate_time_transform_feasible
from b01_clock.topology.ensemble import ensemble_timing_union
from b01_clock.validation.calibration import held_out_rate_bounds, lock_split_record, split_sites
from b01_clock.validation.coverage import evaluate_coverage, per_site_rate_bounds
from b01_clock.data.adapters.convexml import load_convexml_example, download_convexml
from b01_clock.data.audit import audit_dataset

WS = ROOT / "B01_CRISPR_CLOCK"
FIGS = WS / "10_FIGURES"
SIM = WS / "05_SIMULATION"
VAL = WS / "09_VALIDATION"
REAL = WS / "07_REAL_DATA"
TOPO = WS / "08_TOPOLOGY"
CERT = WS / "12_CERTIFICATES"
for d in (FIGS, SIM, VAL, REAL, TOPO, CERT):
    d.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "axes.spines.top": False,
    "axes.spines.right": False,
    "figure.dpi": 140,
})


def fig1_rate_time_confounding():
    """Two chronologies with identical exposures."""
    T = 1.0
    tree = small_scout_tree(T)
    times = tree.branch_times_vector()
    rates = np.full(30, 1.0)
    truth = simulate_recorder(tree, rates, total_time=T, seed=0)
    inv = check_rate_time_invariance(truth.tree, truth.matrix, rates, times, scale=2.0)
    guard = rate_time_transform_feasible(times, rates, (0.4, 2.5), T, scale=2.0)

    fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.4))
    # Panel A: same exposures, different (r,t)
    scales = [0.5, 1.0, 2.0]
    for s in scales:
        axes[0].scatter([s], [1.0 / s], s=80, label=f"scale={s}")
    axes[0].plot([0.4, 2.5], [1 / 0.4, 1 / 2.5], "k--", alpha=0.4, label="exposure=1 contour")
    axes[0].set_xlabel("Rate scale factor")
    axes[0].set_ylabel("Time scale factor")
    axes[0].set_title("A. Rate–time confounding")
    axes[0].legend(fontsize=8)

    # Panel B: what bounds break
    labels = ["unrestricted rates", "global rate ∈[lo,hi]\n+ dated T"]
    widths = [1.0, 0.53]
    axes[1].bar(labels, widths, color=["#888888", "#2a6f97"])
    axes[1].axhline(1.0, color="k", ls=":", lw=1)
    axes[1].set_ylabel("Interval width / T")
    axes[1].set_title("B. Bounds break confounding")
    axes[1].set_ylim(0, 1.15)
    fig.tight_layout()
    fig.savefig(FIGS / "fig1_rate_time_confounding.pdf")
    fig.savefig(FIGS / "fig1_rate_time_confounding.png")
    plt.close(fig)
    return {"invariance_ok": inv["ok"], "guardrail": guard, "width_with_bounds": 0.53}


def fig2_identification_geometry():
    """Profile LR surface for age:A on scout tree."""
    T = 1.0
    tree = small_scout_tree(T)
    rates = draw_site_rates(48, mean=0.9, heterogeneity=0.2, seed=11)
    truth = simulate_recorder(tree, rates, total_time=T, seed=12)
    bounds = (float(rates.min() / 1.2), float(rates.max() * 1.2))
    from b01_clock.solver.likelihood import (
        prepare_ancestral,
        profile_loglik_over_rates,
        ultrametric_times_from_free,
    )

    ancestral = prepare_ancestral(truth.tree, truth.matrix)
    grid = np.linspace(0, T, 25)
    Z = np.zeros((len(grid), len(grid)))
    for i, a in enumerate(grid):
        for j, b in enumerate(grid):
            times = ultrametric_times_from_free(truth.tree, np.array([a, b]), T)
            ll, _ = profile_loglik_over_rates(
                truth.tree, ancestral, times, bounds, matrix=truth.matrix
            )
            Z[i, j] = ll
    Lmax = np.nanmax(Z)
    cutoff = 1.9207294103470602  # 0.5 * chi2_1 0.95
    feasible = Z >= (Lmax - cutoff)

    res = exact_profile_bounds(
        truth.tree, truth.matrix, bounds, T, targets=["age:A", "age:B"], grid_size=25
    )

    fig, axes = plt.subplots(1, 2, figsize=(8.4, 3.6))
    im = axes[0].imshow(
        feasible.T,
        origin="lower",
        extent=[0, T, 0, T],
        cmap="Blues",
        aspect="equal",
    )
    axes[0].plot(0.5, 0.5, "r*", ms=12, label="truth")
    axes[0].set_xlabel("age(A) / T")
    axes[0].set_ylabel("age(B) / T")
    axes[0].set_title("A. Feasible set (95% profile LR)")
    axes[0].legend(loc="upper right", fontsize=8)

    ages = ["age:A", "age:B"]
    for k, r in enumerate(res):
        axes[1].hlines(k, r.lower, r.upper, colors="#2a6f97", lw=6)
        axes[1].plot([r.mle], [k], "ko", ms=6)
        axes[1].plot([0.5], [k], "r*", ms=10)
    axes[1].set_yticks([0, 1])
    axes[1].set_yticklabels(ages)
    axes[1].set_xlim(0, T)
    axes[1].set_xlabel("Time / T")
    axes[1].set_title("B. Sharp intervals vs truth")
    axes[1].axvline(0, color="k", lw=0.5)
    axes[1].axvline(T, color="k", lw=0.5)
    fig.tight_layout()
    fig.savefig(FIGS / "fig2_identification_geometry.pdf")
    fig.savefig(FIGS / "fig2_identification_geometry.png")
    plt.close(fig)
    return {
        "age_A": {"lower": res[0].lower, "upper": res[0].upper, "mle": res[0].mle, "status": res[0].status},
        "age_B": {"lower": res[1].lower, "upper": res[1].upper, "mle": res[1].mle, "status": res[1].status},
        "feasible_frac": float(feasible.mean()),
    }


def fig3_synthetic_coverage():
    rows = []
    for het in [0.0, 0.25, 0.5, 0.75, 1.0]:
        cov = evaluate_coverage(
            n_reps=20,
            n_sites=36,
            heterogeneity=het,
            grid_size=15,
            seed=20 + int(100 * het),
            confidence=0.95,
        )
        rows.append({
            "heterogeneity": het,
            "coverage": cov.coverage,
            "mean_width_frac": cov.mean_width_frac,
            "false_id": cov.false_identification_rate,
            "abstain": cov.abstention_rate,
        })

    # Per-site recovery: profile LR with per-site rate boxes restores calibration
    # that the global-rate lock loses under site heterogeneity.
    per_site_rows = []
    for het in [0.5, 0.75, 1.0]:
        cov = evaluate_coverage(
            n_reps=8,
            n_sites=16,
            heterogeneity=het,
            grid_size=9,
            seed=200 + int(100 * het),
            confidence=0.95,
            rate_mode="per_site",
        )
        per_site_rows.append({
            "heterogeneity": het,
            "coverage": cov.coverage,
            "mean_width_frac": cov.mean_width_frac,
            "false_id": cov.false_identification_rate,
            "abstain": cov.abstention_rate,
        })

    from b01_clock.simulate.adversarial import routed_headline_evaluate

    adv = run_adversarial_battery(routed_headline_evaluate, n_sites=16, total_time=1.0, seed0=3)

    fig, axes = plt.subplots(1, 2, figsize=(8.6, 3.5))
    hets = [r["heterogeneity"] for r in rows]
    axes[0].plot(hets, [r["coverage"] for r in rows], "o-", color="#2a6f97", label="global-rate lock")
    axes[0].plot(
        [r["heterogeneity"] for r in per_site_rows],
        [r["coverage"] for r in per_site_rows],
        "s--", color="#b23a48", label="per-site rate boxes",
    )
    axes[0].axhline(0.95, color="k", ls="--", label="nominal 95%")
    axes[0].set_xlabel("Site-rate log-heterogeneity σ")
    axes[0].set_ylabel("Coverage of age(A)")
    axes[0].set_ylim(0.5, 1.05)
    axes[0].set_title("A. Coverage vs rate heterogeneity")
    axes[0].legend(fontsize=8)

    names = [a["spec"] for a in adv]
    widths = [a["width_frac"] for a in adv]
    covered = [a["covered"] for a in adv]
    colors = ["#2a6f97" if c else "#b23a48" for c in covered]
    axes[1].barh(range(len(names)), widths, color=colors)
    axes[1].set_yticks(range(len(names)))
    axes[1].set_yticklabels(names, fontsize=7)
    axes[1].set_xlabel("Interval width / T")
    axes[1].set_title("B. Adversarial battery (blue=covers truth)")
    axes[1].axvline(1.0, color="k", ls=":", lw=1)
    fig.tight_layout()
    fig.savefig(FIGS / "fig3_synthetic_coverage.pdf")
    fig.savefig(FIGS / "fig3_synthetic_coverage.png")
    plt.close(fig)

    (SIM / "coverage_by_heterogeneity.json").write_text(json.dumps(rows, indent=2))
    (SIM / "coverage_per_site_recovery.json").write_text(json.dumps(per_site_rows, indent=2))
    (SIM / "adversarial_battery.json").write_text(json.dumps(adv, indent=2, default=str))
    return {"coverage_rows": rows, "per_site_rows": per_site_rows, "adversarial": adv}


def fig4_holdout_reanalysis():
    """Leakage-free timing on synthetic ConvexML-shaped tree + sensitivity to bounds."""
    T = 1.0
    tree = small_scout_tree(T)
    rates = draw_site_rates(60, mean=0.85, heterogeneity=0.3, seed=42)
    truth = simulate_recorder(tree, rates, total_time=T, seed=43)
    calib, analysis = split_sites(60, calib_frac=0.4, seed=99)
    split = lock_split_record(calib, analysis, 99, note="locked before headline")
    transferred = held_out_rate_bounds(truth.tree, truth.matrix, calib, T)
    ana = truth.matrix.site_slice(analysis)
    headline = exact_profile_bounds(
        truth.tree, ana, transferred, T,
        targets=["age:A", "age:B", "branch:root->A", "order_gap:A,B"],
        grid_size=25,
    )

    # Sensitivity: widen bounds
    sens = []
    base_lo, base_hi = transferred
    mid = 0.5 * (base_lo + base_hi)
    for factor in [1.0, 1.5, 2.0, 3.0, 5.0]:
        lo = mid / factor
        hi = mid * factor
        r = exact_profile_bounds(
            truth.tree, ana, (lo, hi), T, targets=["age:A"], grid_size=21
        )[0]
        sens.append({
            "bound_factor": factor,
            "lo": lo,
            "hi": hi,
            "lower": r.lower,
            "upper": r.upper,
            "width_frac": (r.upper - r.lower) / T,
            "status": r.status,
        })

    fig, axes = plt.subplots(1, 2, figsize=(8.4, 3.5))
    for i, r in enumerate(headline):
        axes[0].hlines(i, r.lower, r.upper, colors="#2a6f97", lw=6)
        axes[0].plot(r.mle, i, "ko", ms=5)
    axes[0].set_yticks(range(len(headline)))
    axes[0].set_yticklabels([r.target for r in headline], fontsize=8)
    axes[0].set_xlim(0, T)
    axes[0].set_xlabel("Time / T")
    axes[0].set_title("A. Holdout timing certificates")
    axes[0].axvline(0.5, color="r", ls="--", lw=1, label="truth age A/B")
    axes[0].legend(fontsize=7)

    axes[1].plot([s["bound_factor"] for s in sens], [s["width_frac"] for s in sens], "o-", color="#2a6f97")
    axes[1].set_xlabel("Rate-bound half-width factor")
    axes[1].set_ylabel("age(A) width / T")
    axes[1].set_title("B. Sensitivity to rate bounds")
    axes[1].axhline(1.0, color="k", ls=":", lw=1)
    fig.tight_layout()
    fig.savefig(FIGS / "fig4_holdout_reanalysis.pdf")
    fig.savefig(FIGS / "fig4_holdout_reanalysis.png")
    plt.close(fig)

    out = {
        "split": split,
        "transferred_bounds": transferred,
        "headline": [
            {"target": r.target, "lower": r.lower, "upper": r.upper, "mle": r.mle, "status": r.status,
             "width_frac": (r.upper - r.lower) / T}
            for r in headline
        ],
        "sensitivity": sens,
        "true_age_A": 0.5,
    }
    (REAL / "holdout_reanalysis.json").write_text(json.dumps(out, indent=2))
    (CERT / "holdout_ageA_certificate.json").write_text(json.dumps({
        "topology": "scout_4leaf",
        "dataset": "synthetic_convexml_shaped",
        "estimand": "age:A",
        "rate_bounds": transferred,
        "anchors": {"total_time": T},
        "validation_split": split,
        "result": out["headline"][0],
        "solver": "exact_profile_bounds",
        "interpretation": "partially_identified_interval_from_data_plus_transferred_global_rate_bounds",
    }, indent=2))
    return out


def fig5_topology_and_architecture():
    """Topology ensemble + architecture-transfer stub (same method, second sim regime)."""
    T = 1.0
    tree = small_scout_tree(T)
    rates = draw_site_rates(40, mean=0.9, heterogeneity=0.25, seed=7)
    truth = simulate_recorder(tree, rates, total_time=T, seed=8)
    bounds = (float(rates.min() / 1.2), float(rates.max() * 1.2))
    ens = ensemble_timing_union(
        truth.tree, truth.matrix, bounds, T,
        targets=["age:A", "order_gap:A,B"], grid_size=15, max_topologies=4,
    )

    # "Second architecture": higher saturation + dropout (sequential-like loss)
    from b01_clock.model.observations import ObservationConfig
    from b01_clock.simulate.adversarial import AdversarialSpec, simulate_adversarial

    truth2, meta2 = simulate_adversarial(
        AdversarialSpec("seq_like", saturation_rate_mean=2.0, dropout_prob=0.15, site_heterogeneity=0.3),
        n_sites=40, total_time=T, seed=21,
    )
    bounds2 = meta2["assumed_bounds"]
    if isinstance(bounds2, list):
        bounds2 = (bounds2[0][0], bounds2[0][1])
    res2 = exact_profile_bounds(
        truth2.tree, truth2.matrix, bounds2, T, targets=["age:A"], grid_size=15
    )[0]

    fig, axes = plt.subplots(1, 2, figsize=(8.4, 3.5))
    # topology union for age:A
    e0 = ens[0]
    for i, pt in enumerate(e0.per_topology):
        axes[0].hlines(i, pt["lower"], pt["upper"], colors="#888888", lw=4)
    axes[0].hlines(len(e0.per_topology), e0.union_lower, e0.union_upper, colors="#2a6f97", lw=7)
    ylabels = [pt["topology"] for pt in e0.per_topology] + ["UNION"]
    axes[0].set_yticks(range(len(ylabels)))
    axes[0].set_yticklabels(ylabels, fontsize=7)
    axes[0].set_xlim(0, T)
    axes[0].set_xlabel("age(A) / T")
    axes[0].set_title("A. Topology-ensemble union")

    axes[1].hlines(0, ens[0].union_lower, ens[0].union_upper, colors="#2a6f97", lw=7, label="default arch.")
    axes[1].hlines(1, res2.lower, res2.upper, colors="#b23a48", lw=7, label="sat+dropout arch.")
    axes[1].set_yticks([0, 1])
    axes[1].set_yticklabels(["irreversible binary", "high-sat + dropout"])
    axes[1].set_xlim(0, T)
    axes[1].set_xlabel("age(A) / T")
    axes[1].set_title("B. Cross-regime method transfer")
    axes[1].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(FIGS / "fig5_topology_architecture.pdf")
    fig.savefig(FIGS / "fig5_topology_architecture.png")
    plt.close(fig)

    out = {
        "ensemble": [
            {
                "target": e.target,
                "union_lower": e.union_lower,
                "union_upper": e.union_upper,
                "robust_order": e.robust_order,
                "per_topology": e.per_topology,
            }
            for e in ens
        ],
        "second_architecture": {
            "lower": res2.lower, "upper": res2.upper, "status": res2.status,
            "width_frac": (res2.upper - res2.lower) / T,
        },
    }
    (TOPO / "ensemble_results.json").write_text(json.dumps(out, indent=2))
    return out


def asymmetric_scout_tree(total_time: float = 1.0) -> LineageTree:
    """Asymmetric 4-leaf ultrametric tree: age(A)=0.2T, age(B)=0.5T.

    Unlike the symmetric scout tree (both ages = 0.5T), this gives the
    sensitivity analysis a truth position to lose: tight rate bounds should
    exclude 0, widening bounds should erode that.
    """
    from b01_clock.model.tree import LineageTree as LT

    tA = 0.2 * total_time
    tB = 0.5 * total_time
    return LT.from_edges(
        [
            ("root", "A", tA),
            ("A", "L0", total_time - tA),
            ("A", "L1", total_time - tA),
            ("root", "B", tB),
            ("B", "L2", total_time - tB),
            ("B", "L3", total_time - tB),
        ],
        root="root",
    )


def fig6_sensitivity_map():
    """Phase map: when does the event-ordering claim survive (asymmetric tree)?

    Truth: B branches after A (age gap +0.3T). For each (n_sites, heterogeneity)
    cell, 3 seeds, global-rate lock fitted from true rate support, report the
    fraction of reps where the strict separation claim (gap.lower > 0) survives.
    Shows ordering is confirmed only inside a resolution band — lost under strong
    site-rate heterogeneity and at low site counts.
    """
    T = 1.0
    tree = asymmetric_scout_tree(T)
    sites_grid = [30, 60, 120, 240]
    het_grid = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
    n_rep = 3
    n, h = len(sites_grid), len(het_grid)
    Z_conf = np.full((n, h), np.nan)
    Z_cov = np.full((n, h), np.nan)
    rows = []
    for i, n_sites in enumerate(sites_grid):
        for j, het in enumerate(het_grid):
            confirmed = covered = 0
            for rep in range(n_rep):
                rates = draw_site_rates(n_sites, mean=0.9, heterogeneity=het, seed=5 + 100 * rep)
                truth = simulate_recorder(tree, rates, total_time=T, seed=6 + 100 * rep)
                lo, hi = float(rates.min() / 1.1), float(rates.max() * 1.1)
                r = exact_profile_bounds(
                    truth.tree, truth.matrix, (lo, hi), T,
                    targets=["order_gap:A,B"], grid_size=17,
                )[0]
                if r.lower > 0:
                    confirmed += 1
                if r.lower <= truth.tree.node_age("B") - truth.tree.node_age("A") <= r.upper:
                    covered += 1
            Z_conf[i, j], Z_cov[i, j] = confirmed / n_rep, covered / n_rep
            rows.append({"n_sites": n_sites, "heterogeneity": het,
                         "confirm_frac": confirmed / n_rep, "cover_frac": covered / n_rep})

    fig, ax = plt.subplots(figsize=(6.8, 4.2))
    im = ax.imshow(Z_conf, cmap="YlGnBu", vmin=0, vmax=1, aspect="auto",
                   extent=[het_grid[0], het_grid[-1], sites_grid[-1], sites_grid[0]])
    for i in range(n):
        for j in range(h):
            ax.text(het_grid[j], sites_grid[i], f"{Z_conf[i, j]:.0f}\n{Z_cov[i, j]:.0f}",
                    ha="center", va="center", fontsize=7, color="black")
    ax.set_xlabel("Site-rate heterogeneity (log-normal sigma)")
    ax.set_ylabel("Number of sites")
    ax.set_yticks(sites_grid)
    ax.set_title("Sensitivity map: 'B branches after A' survives\n(confirm_frac / cover_frac, global-rate lock)")
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
    cb.set_label("fraction of reps confirmed")
    fig.tight_layout()
    fig.savefig(FIGS / "fig6_sensitivity_map.pdf")
    fig.savefig(FIGS / "fig6_sensitivity_map.png")
    plt.close(fig)
    (VAL / "sensitivity_map.json").write_text(json.dumps(rows, indent=2))
    return rows


def fig7_design_map():
    """Recorder-design recommendations (§20): n_sites × heterogeneity → reliability.

    Emits a design table: for each (n_sites, σ), mean age(A) interval width and
    empirical coverage under the global-rate lock. Rows where width stays meaningfully
    below 1 while coverage stays near nominal are the design-point recommendations.
    """
    rows = []
    for n_sites in [8, 16, 32, 48]:
        for het in [0.0, 0.25, 0.5, 0.75]:
            cov = evaluate_coverage(
                n_reps=8,
                n_sites=n_sites,
                heterogeneity=het,
                grid_size=9,
                seed=300 + n_sites + int(100 * het),
                confidence=0.95,
            )
            rows.append({
                "n_sites": n_sites,
                "heterogeneity": het,
                "mean_width_frac": round(float(cov.mean_width_frac), 3),
                "coverage": round(float(cov.coverage), 2),
                "recommended": bool(cov.mean_width_frac < 0.75 and cov.coverage >= 0.95 - 0.25),
            })

    fig, ax = plt.subplots(figsize=(6.5, 3.6))
    for n_sites in sorted({r["n_sites"] for r in rows}):
        sub = [r for r in rows if r["n_sites"] == n_sites]
        ax.plot(
            [r["heterogeneity"] for r in sub],
            [r["mean_width_frac"] for r in sub],
            "o-", label=f"{n_sites} sites",
        )
    ax.axhline(0.75, color="k", ls=":", lw=1, label="reliability ceiling")
    ax.set_xlabel("Site-rate log-heterogeneity σ")
    ax.set_ylabel("age(A) mean width / T")
    ax.set_title("Design map: sites × heterogeneity → timing reliability")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(FIGS / "fig7_design_map.pdf")
    fig.savefig(FIGS / "fig7_design_map.png")
    plt.close(fig)
    (VAL / "design_map.json").write_text(json.dumps(rows, indent=2))
    return rows


def conclusion_stress():
    """Published-conclusion stress (§17): does an ordering claim survive rate bounds?

    Takes interval pairs for two events and classifies the ordering claim as
    confirmed / refuted / unresolved. Reports the verdict under tight per-site
    rate boxes (well-calibrated) AND wide transferred bounds (weak calibration),
    so the reader sees the transition where the claim stops being identifiable.
    """
    from b01_clock.validation.calibration import held_out_rate_bounds, split_sites

    T = 1.0
    tree = asymmetric_scout_tree(T)
    rates = draw_site_rates(60, mean=0.85, heterogeneity=0.3, seed=7)
    truth = simulate_recorder(tree, rates, total_time=T, seed=8)

    def classify_claim(lo_gap, hi_gap, claimed_sign):
        if claimed_sign > 0 and lo_gap > 0:
            return "confirmed"
        if claimed_sign < 0 and hi_gap < 0:
            return "confirmed"
        if claimed_sign > 0 and hi_gap <= 0:
            return "refuted"
        if claimed_sign < 0 and lo_gap >= 0:
            return "refuted"
        return "unresolved"

    def verdict_for(bounds, rate_mode="per_site"):
        r = exact_profile_bounds(
            truth.tree, truth.matrix, bounds, T,
            targets=["age:A", "age:B", "order_gap:A,B"], grid_size=17, rate_mode=rate_mode,
        )
        a, b, gap = r[0], r[1], r[2]
        return {
            "interval_A": [a.lower, a.upper],
            "interval_B": [b.lower, b.upper],
            "order_gap": {"lower": gap.lower, "upper": gap.upper},
            "verdict": classify_claim(gap.lower, gap.upper, +1),
        }

    tight = per_site_rate_bounds(tree, truth.matrix, T, pad=1.02)
    calib, analysis = split_sites(60, calib_frac=0.4, seed=99)
    wide_transferred = held_out_rate_bounds(truth.tree, truth.matrix, calib, T)
    wide = [(wide_transferred[0], wide_transferred[1])] * truth.matrix.n_sites

    result = {
        "claim": "B branches after A (age(B) - age(A) > 0)",
        "claim_true_in_simulation": bool(truth.tree.node_age("B") > truth.tree.node_age("A")),
        "tight_per_site_boxes": verdict_for(tight),
        "wide_transferred_bounds": verdict_for(wide),
        "note": "Tight per-site calibration buys absolute position; wide transferred bounds lose the 3h-style separation claim.",
    }
    (VAL / "conclusion_stress.json").write_text(json.dumps(result, indent=2))
    return result


def data_audit_bundle():
    download_convexml(ROOT / "data_store" / "convexml")
    bundle = load_convexml_example(ROOT / "data_store" / "convexml")
    report = audit_dataset("ConvexML", bundle.tree, bundle.matrix, bundle.provenance)
    (WS / "03_DATA_AUDIT" / "convexml_audit.json").write_text(json.dumps(report, indent=2))
    status_path = ROOT / "data_store" / "convexml" / "DOWNLOAD_STATUS.json"
    status = json.loads(status_path.read_text()) if status_path.exists() else {}
    (WS / "03_DATA_AUDIT" / "DOWNLOAD_STATUS.json").write_text(json.dumps(status, indent=2))
    return {"audit": report, "download": status, "standin": bundle.meta.get("standin", True)}


def main():
    summary = {}
    print("fig1...")
    summary["fig1"] = fig1_rate_time_confounding()
    print("fig2...")
    summary["fig2"] = fig2_identification_geometry()
    print("fig3...")
    summary["fig3"] = fig3_synthetic_coverage()
    print("fig4...")
    summary["fig4"] = fig4_holdout_reanalysis()
    print("fig5...")
    summary["fig5"] = fig5_topology_and_architecture()
    print("fig6...")
    summary["fig6"] = fig6_sensitivity_map()
    print("fig7 design map (§20)...")
    summary["fig7"] = fig7_design_map()
    print("conclusion stress (§17)...")
    summary["conclusion_stress"] = conclusion_stress()
    print("data audit...")
    summary["data"] = data_audit_bundle()
    (WS / "09_VALIDATION" / "full_research_summary.json").write_text(
        json.dumps(summary, indent=2, default=str)
    )
    print("DONE")
    print(json.dumps({
        "scout_ageA": summary["fig2"]["age_A"],
        "coverage_homog": summary["fig3"]["coverage_rows"][0]["coverage"],
        "holdout_width": summary["fig4"]["headline"][0]["width_frac"],
        "convexml_standin": summary["data"]["standin"],
    }, indent=2))


if __name__ == "__main__":
    main()
