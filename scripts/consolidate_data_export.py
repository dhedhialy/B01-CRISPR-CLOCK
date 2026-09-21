"""One-file consolidation of every B01 result + raw character data for the audit zip run.

Builds, in the user's Downloads dir:
  B01_consolidated_results.xlsx         - every derived result (benchmark certs,
                                          ground-truth points, recalibration
                                          symmetric + conditional, method
                                          exploration, real-data GESTALT run,
                                          plus a file manifest sheet)
  B01_raw_character_matrices_tidy.csv.gz - ALL raw site-level observations
                                          (256 ConvexML matrices + the GESTALT
                                          MIX matrix) as one tidy long table
                                          (tree_id, cell, site, state)

Reproduce:  python -m scripts.consolidate_data_export
"""

from __future__ import annotations

import csv
import gzip
import json
from pathlib import Path

from openpyxl import Workbook
from openpyxl.utils import get_column_letter

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "B01_CRISPR_CLOCK" / "07_REAL_DATA"
EXPERIMENTS = ROOT / "experiments"
DATA = ROOT / "data_store"
OUT_DIR = Path.home() / "Downloads"

XLSX_OUT = OUT_DIR / "B01_consolidated_results.xlsx"
LONG_OUT = OUT_DIR / "B01_raw_character_matrices_tidy.csv.gz"
WIDE_OUT = OUT_DIR / "B01_original_data.csv"

RESULT_FILES = [
    RESULTS / "convexml_protocol_certificate.json",
    RESULTS / "convexml_ground_truth_coverage.json",
    RESULTS / "convexml_protocol_reanalysis.json",
    RESULTS / "holdout_reanalysis.json",
    RESULTS / "gestalt_real_consistency.json",
    EXPERIMENTS / "recalibration_study.json",
    EXPERIMENTS / "recalibration_conditional.json",
    EXPERIMENTS / "coverage_methods_exploration.json",
]


def _records(rows, keys) -> "list[list]":
    return [[r.get(k, "") for k in keys] for r in rows]


def _scalar(v):
    return v if isinstance(v, (str, int, float, bool, type(None))) else str(v)


def add_sheet(wb, name, headers, rows):
    ws = wb.create_sheet(name)
    ws.append([_scalar(h) for h in headers])
    for r in rows:
        ws.append([_scalar(c) for c in r])
    for i, h in enumerate(headers, start=1):
        ws.column_dimensions[get_column_letter(i)].width = max(8, min(40, len(str(h)) + 2))
    return ws


def manifest_rows() -> "list[list]":
    rows = []
    srcs = [
        *RESULT_FILES,
        DATA / "convexml" / "convexml_dryad_trees.tgz",
        DATA / "convexml" / "extracted",
        DATA / "gestalt" / "GSE81713_fish_ADR1_PHYLIP_MIX_gte5_input.txt.gz",
        DATA / "gestalt" / "GSE81713_fish_ADR1_PHYLIP_MIX_gte5_output.newick.txt.gz",
        DATA / "gestalt" / "GSE81713_fish_ADR1_PHYLIP_MIX_gte5_json.txt.gz"
        if (DATA / "gestalt" / "GSE81713_fish_ADR1_PHYLIP_MIX_gte5_json.txt.gz").exists()
        else DATA / "gestalt" / "GSE81713_fish_ADR1_PHYLIP_MIX_gte5.json.gz",
        XLSX_OUT,
        LONG_OUT,
        WIDE_OUT,
    ]
    for p in srcs:
        if not p.exists():
            rows.append([p.name, str(p), "", "MISSING", ""])
            continue
        n = "dir" if p.is_dir() else f"{p.stat().st_size:,} B"
        rows.append([p.name, str(p), n, "ok", ""])
    return rows


def add_manifest(wb):
    if "manifest" in wb.sheetnames:
        del wb["manifest"]
    rows = manifest_rows()
    add_sheet(wb, "manifest",
              ["file", "path", "size", "status",
               "note/description"],
              [r + [""] for r in rows])


def add_benchmark(wb):
    cert = json.loads((RESULTS / "convexml_protocol_certificate.json").read_text())
    add_sheet(wb, "benchmark_est_single", ["field", "value"],
              [[k, v] for k, v in {
                  "estimand": cert.get("estimand"),
                  "dataset": cert.get("dataset"),
                  "dryad_doi": cert.get("dryad_doi"),
                  "lower": cert["result"]["lower"],
                  "upper": cert["result"]["upper"],
                  "mle": cert["result"]["mle"],
                  "status": cert["result"]["status"],
                  "width_frac": cert["result"]["width_frac"],
                  "rate_bounds": cert.get("rate_bounds"),
                  "total_time": cert.get("anchors", {}).get("total_time"),
                  "solver": cert.get("solver"),
              }.items()])

    gt = json.loads((RESULTS / "convexml_ground_truth_coverage.json").read_text())
    gk = ["target", "true_age", "lower", "upper", "mle", "status", "covered", "width_frac"]
    add_sheet(wb, "benchmark_ground_truth",
              gk + ["dataset"], [r + [gt.get("dataset")] for r in _records(gt["rows"], gk)])

    rean = json.loads((RESULTS / "convexml_protocol_reanalysis.json").read_text())
    ik = ["target", "lower", "upper", "mle", "status", "width_frac"]
    add_sheet(wb, "benchmark_reanalysis",
              ik + ["benchmark"], [r + ["convexml_reanalysis"] for r in _records(rean["headline"], ik)])

    hol = json.loads((RESULTS / "holdout_reanalysis.json").read_text())
    add_sheet(wb, "benchmark_holdout",
              ["target", "lower", "upper", "mle", "status", "width_frac", "regime"],
              [r + [None] for r in _records(hol["headline"], ik)])


def add_recalibration(wb):
    study = json.loads((EXPERIMENTS / "recalibration_study.json").read_text())

    rows = []
    sums = [["core_summary", "canonical_draw13", study["canonical_draw13"].get(x)]
                                for x in ["n_nodes", "profile_coverage", "conformal_loo_coverage", "conformal_width", "mean_bias", "mean_abs_bias"]]
    sums += [["core_summary", f"alt_draw42.{x}", study["alt_draw42"].get(x)]
                                for x in ["n_nodes", "profile_coverage", "conformal_loo_coverage", "conformal_width", "mean_bias", "mean_abs_bias"]]
    sums += [["core_summary", f"transfer_13_to_42.{x}", study["transfer_13->42"].get(x)]
                                for x in ["coverage", "width"]]
    sums += [["core_summary", f"transfer_42_to_13.{x}", study["transfer_42->13"].get(x)]
                                for x in ["coverage", "width"]]
    add_sheet(wb, "recalibration_symmetric_summary", ["group", "metric", "value"], sums)

    def node_rows(draw):
        k = ["node", "true", "mle", "lo", "hi", "profile_covered", "resid"]
        return _records(study[draw]["rows"], k)
    add_sheet(wb, "recalibration_draw13_rows",
              ["node", "true_age", "mle_age", "symmetric_lo", "symmetric_hi", "profile_covered", "resid"],
              node_rows("canonical_draw13"))
    add_sheet(wb, "recalibration_draw42_rows",
              ["node", "true_age", "mle_age", "symmetric_lo", "symmetric_hi", "profile_covered", "resid"],
              node_rows("alt_draw42"))

    cond = json.loads((EXPERIMENTS / "recalibration_conditional.json").read_text())
    rows = [["conditional_loo", "canonical_draw13", cond["canonical_draw13"]["conditional_loo_coverage"]],
            ["conditional_loo", "alt_draw42", cond["alt_draw42"]["conditional_loo_coverage"]],
            ["profile", "canonical_draw13", cond["canonical_draw13"]["profile_coverage"]],
            ["transfer", "13->42", cond["transfer_13->42"]["coverage"]],
            ["transfer", "42->13", cond["transfer_42->13"]["coverage"]],
            ["panel", "coverage", cond["canonical_panel_9_nodes_disjoint"]["coverage"]],
            ["panel", "mean_width", cond["canonical_panel_9_nodes_disjoint"]["mean_width"]],
            ["panel", "train_n_nodes", cond["canonical_panel_9_nodes_disjoint"]["train_n_nodes"]],
            ["panel", "panel_n_nodes", cond["canonical_panel_9_nodes_disjoint"]["panel_n_nodes"]]]
    add_sheet(wb, "recalibration_conditional_summary", ["group", "metric", "value"], rows)
    add_sheet(wb, "recalibration_conditional_panel",
              ["node", "mle", "true", "width", "covered"],
              _records(cond["canonical_panel_9_nodes_disjoint"]["rows"],
                       ["node", "mle", "true", "width", "covered"]))


def add_exploration(wb):
    exp = json.loads((EXPERIMENTS / "coverage_methods_exploration.json").read_text())
    rows = []
    for method, d in exp.items():
        rows.append([method, d.get("coverage"), d.get("n"), d.get("mean_width"), d.get("mean_shift"),
                     d.get("note") or d.get("comment") or d.get("verdict") or d.get("description")])
    add_sheet(wb, "coverage_exploration_summary",
              ["method", "coverage", "n", "mean_width", "mean_shift", "note"], rows)
    for method, d in exp.items():
        if "rows" not in d:
            continue
        k = list(d["rows"][0].keys()) if d["rows"] else []
        add_sheet(wb, f"coverage_row_{method[:19]}", k, _records(d["rows"], k))


def add_gestalt(wb):
    b = json.loads((RESULTS / "gestalt_real_consistency.json").read_text())
    audit = wb.create_sheet("gestalt_audit")
    for k, v in b.get("audit", {}).items():
        audit.append([k, v])

    ck = ["target", "lower_dp", "upper_dp", "mle_dp", "width_dp", "width_frac_of_T",
          "status", "depth", "boundary_clipped_artifact"]
    add_sheet(wb, "gestalt_certificates", ck, _records(b["certificates"], ck))
    add_sheet(wb, "gestalt_certificate_summary", ck, _records(b["certificate_summary"], ck))
    rk = ["target", "mle_dp", "conformal_lo_dp", "conformal_hi_dp", "in_domain", "width_dp"]
    add_sheet(wb, "gestalt_calibrated", rk, _records(b["recalibrated"], rk))
    add_sheet(wb, "gestalt_conclusions",
              ["conclusion", "value"],
              [[k, v] for k, v in b["conclusions"].items()])
    add_sheet(wb, "gestalt_effort",
              ["parameter", "value"],
              [[k, v] for k, v in b["effort"].items()])


def build_xlsx() -> Path:
    wb = Workbook()
    wb.remove(wb.active)
    add_manifest(wb)
    add_benchmark(wb)
    add_recalibration(wb)
    add_exploration(wb)
    add_gestalt(wb)
    wb.save(XLSX_OUT)
    return XLSX_OUT


def build_long_csv() -> Path:
    header = ["tree_id", "cell", "site", "state"]
    with gzip.open(LONG_OUT, "wt", newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        n = 0
        for csv_file in sorted((DATA / "convexml" / "extracted").rglob("*_character_matrix.csv")):
            tree_id = str(csv_file.relative_to(DATA / "convexml" / "extracted"))
            tree_id = tree_id.replace("_character_matrix.csv", "").replace("/trees/", "/")
            with csv_file.open() as fh:
                rd = csv.DictReader(fh)
                cell_col = "cell_id" if "cell_id" in rd.fieldnames else rd.fieldnames[0]
                for row in rd:
                    cell = row[cell_col]
                    for site_name, val in row.items():
                        if site_name == cell_col:
                            continue
                        site = int(site_name) if not site_name.startswith("s") else int(site_name[1:])
                        w.writerow([tree_id, cell, site, int(val)])
                        n += 1
        mf = DATA / "gestalt" / "GSE81713_fish_ADR1_PHYLIP_MIX_gte5_input.txt.gz"
        import b01_clock.data.adapters.gestalt as ga  # noqa
        m, _ = ga.load_mix_matrix(mf)
        for cell, vec in m.leaf_states.items():
            for site, val in enumerate(vec):
                w.writerow(["gestalt_ADR1_gte5", cell, site, int(val)])
                n += 1
    return LONG_OUT, n


def build_wide_csv() -> "tuple[Path, int]":
    """One plain (uncompressed) CSV with the ORIGINAL wide matrices + GESTALT.

    All ConvexML matrices are 39 sites wide, so that grid is preserved verbatim
    (tree_id, cell, site0..site38). GESTALT cells carry their original 1156-char
    MIX state string in a single `genome` column (as recorded in the source
    file) instead of being padded across 1156 site columns. Under 50 MB.
    """
    import b01_clock.data.adapters.gestalt as ga  # noqa

    metas = []
    for csv_file in sorted((DATA / "convexml" / "extracted").rglob("*_character_matrix.csv")):
        tree_id = str(csv_file.relative_to(DATA / "convexml" / "extracted"))
        tree_id = tree_id.replace("_character_matrix.csv", "").replace("/trees/", "/")
        with csv_file.open() as fh:
            fieldnames = next(csv.reader(fh))
        site_name = {}
        for f in fieldnames:
            if f in ("", "cell_id"):
                continue
            site_name[int(f[1:]) if f.startswith("s") else int(f)] = f
        metas.append({"tree_id": tree_id, "path": csv_file, "site_name": site_name,
                      "cell_col": "cell_id" if "cell_id" in fieldnames else None})

    m, _ = ga.load_mix_matrix(DATA / "gestalt" / "GSE81713_fish_ADR1_PHYLIP_MIX_gte5_input.txt.gz")
    n_sites = max((max(s["site_name"]) for s in metas), default=0) + 1
    genome_col = f"gestalt_mix_string_s{n_sites}+"

    header = ["tree_id", "cell"] + [f"site{i}" for i in range(n_sites)] + [genome_col]
    nrows = 0
    with WIDE_OUT.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        for meta in metas:
            with meta["path"].open() as fh:
                rd = csv.DictReader(fh)
                cell_col = meta["cell_col"] or rd.fieldnames[0]
                site_name = meta["site_name"]
                for row in rd:
                    vals = [""] * (n_sites + 1)
                    for s, name in site_name.items():
                        vals[s] = row[name]
                    w.writerow([meta["tree_id"], row[cell_col], *vals])
                    nrows += 1
        for cell, vec in m.leaf_states.items():
            row = [""] * (n_sites + 1)
            row[-1] = "".join(str(int(v)) for v in vec.tolist())
            w.writerow(["gestalt_ADR1_gte5", cell, *row])
            nrows += 1
    return WIDE_OUT, nrows


if __name__ == "__main__":
    xlsx = build_xlsx()
    long, n = build_long_csv()
    wide, wrows = build_wide_csv()
    print(f"xlsx: {xlsx} ({xlsx.stat().st_size:,} B)")
    print(f"long tidy csv.gz: {long} ({long.stat().st_size:,} B, {n:,} tidy rows)")
    print(f"wide original csv: {wide} ({wide.stat().st_size:,} B, {wrows:,} rows, "
          f"{(wide.stat().st_size < 50_000_000) and 'UNDER 50 MB' or 'OVER 50 MB'})")