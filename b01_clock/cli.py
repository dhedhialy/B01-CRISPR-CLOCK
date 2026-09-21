"""Command-line entry points for simulation, fitting, validation, and external-benchmark analysis."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import yaml


def _load_config(path: str | None) -> dict:
    if path is None:
        return {}
    with open(path) as f:
        return yaml.safe_load(f) or {}


def cmd_scout(args: argparse.Namespace) -> int:
    from b01_clock.validation.scout import run_decisive_scout
    from b01_clock.data.adapters.convexml import download_convexml, load_convexml_example

    cfg = _load_config(args.config)
    out = Path(args.out or cfg.get("out_dir", "results/scout"))
    data_root = Path(cfg.get("convexml_dir", "data_store/convexml"))
    if args.download:
        download_convexml(data_root)
    bundle = load_convexml_example(data_root, total_time=cfg.get("total_time", 1.0))
    report = run_decisive_scout(
        out_dir=out,
        n_sites=cfg.get("n_sites", 60),
        total_time=cfg.get("total_time", 1.0),
        confidence=cfg.get("confidence", 0.95),
        coverage_reps=cfg.get("coverage_reps", 24),
        seed=cfg.get("seed", 7),
        convexml_matrix=bundle.matrix if not bundle.meta.get("standin") else None,
        convexml_tree=bundle.tree if not bundle.meta.get("standin") else None,
    )
    # Always also run holdout on bundle for wiring
    print(json.dumps({"verdict": report["verdict"], "pass": report["pass"], "out": str(out)}, indent=2))
    return 0 if report["pass"] else 2


def cmd_simulate(args: argparse.Namespace) -> int:
    from b01_clock.simulate.adversarial import run_adversarial_battery, DEFAULT_BATTERY, routed_headline_evaluate

    cfg = _load_config(args.config)
    total_time = cfg.get("total_time", 1.0)

    rows = run_adversarial_battery(
        routed_headline_evaluate,
        n_sites=cfg.get("n_sites", 40),
        total_time=total_time,
        seed0=cfg.get("seed", 0),
    )
    out = Path(args.out or "results/adversarial.json")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(rows, indent=2, default=str))
    print(f"Wrote {out} ({len(rows)} specs)")
    return 0


def cmd_fit(args: argparse.Namespace) -> int:
    from b01_clock.model.tree import small_scout_tree
    from b01_clock.simulate.generator import draw_site_rates, simulate_recorder
    from b01_clock.solver.exact import exact_profile_bounds
    from b01_clock.topology.ensemble import ensemble_timing_union

    cfg = _load_config(args.config)
    T = cfg.get("total_time", 1.0)
    tree = small_scout_tree(T)
    rates = draw_site_rates(cfg.get("n_sites", 50), mean=0.9, heterogeneity=0.3, seed=cfg.get("seed", 1))
    truth = simulate_recorder(tree, rates, total_time=T, seed=cfg.get("seed", 1) + 1)
    bounds = (float(rates.min() / 1.2), float(rates.max() * 1.2))
    exact = exact_profile_bounds(
        truth.tree, truth.matrix, bounds, T, targets=["age:A", "order_gap:A,B"], grid_size=41
    )
    ens = ensemble_timing_union(
        truth.tree, truth.matrix, bounds, T, targets=["age:A", "order_gap:A,B"], grid_size=25
    )
    out = {
        "exact": [r.__dict__ for r in exact],
        "ensemble": [e.__dict__ for e in ens],
    }
    path = Path(args.out or "results/fit.json")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=2, default=str))
    print(f"Wrote {path}")
    return 0


def cmd_audit(args: argparse.Namespace) -> int:
    from b01_clock.data.adapters.convexml import load_convexml_example
    from b01_clock.data.audit import audit_dataset

    bundle = load_convexml_example(args.data_root)
    report = audit_dataset("ConvexML", bundle.tree, bundle.matrix, bundle.provenance)
    path = Path(args.out or "results/audit_convexml.json")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2))
    return 0


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(prog="b01-clock", description="B01 CRISPR recorder timing certificates")
    sub = p.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("scout", help="Run first decisive scout")
    s.add_argument("--config", default="configs/scout.yaml")
    s.add_argument("--out", default=None)
    s.add_argument("--download", action="store_true")
    s.set_defaults(func=cmd_scout)

    s = sub.add_parser("simulate", help="Adversarial simulation battery")
    s.add_argument("--config", default="configs/default.yaml")
    s.add_argument("--out", default=None)
    s.set_defaults(func=cmd_simulate)

    s = sub.add_parser("fit", help="Fit timing certificates on synthetic scout tree")
    s.add_argument("--config", default="configs/default.yaml")
    s.add_argument("--out", default=None)
    s.set_defaults(func=cmd_fit)

    s = sub.add_parser("audit", help="Audit ConvexML/local dataset")
    s.add_argument("--data-root", default="data_store/convexml")
    s.add_argument("--out", default=None)
    s.set_defaults(func=cmd_audit)

    args = p.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
