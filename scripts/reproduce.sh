#!/usr/bin/env bash
# One-command reproduction of primary scout figure/table artifacts.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
# shellcheck disable=SC1091
source .venv/bin/activate
export PYTHONPATH="$ROOT${PYTHONPATH:+:$PYTHONPATH}"

python -m pytest -q
python scripts/run_full_research.py          # figures + 09_VALIDATION summaries + data audit
python scripts/run_convexml_reanalysis.py    # external-benchmark certificate on official Dryad tree
python scripts/validate_convexml_ground_truth.py   # ground-truth coverage vs known newick lengths
python -m b01_clock scout --config configs/scout.yaml --out results/scout
python -m b01_clock simulate --config configs/default.yaml --out results/adversarial.json
python -m b01_clock fit --config configs/default.yaml --out results/fit.json
python -m b01_clock audit --data-root data_store/convexml --out results/audit_convexml.json

echo "Reproduction complete. See B01_CRISPR_CLOCK/09_VALIDATION, 10_FIGURES, results/"
