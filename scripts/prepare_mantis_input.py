"""Prepare B01 lineage data as a Mantis space input (semantic metadata columns).

Writes to Downloads:
  B01_mantis_space.csv          - one row per cell; semantic columns for the
                                  Mantis 'metadata' labeling config type
                                  (Cell ID, Dataset, Sample Family, Analysis
                                  Subtree membership, Depth, embedding Label).
  B01_mantis_space_config.md    - Create-a-Space parameters in Mantis key/value
                                  style, ending in Generate Landmarks.
  B01_mantis_space_config.yaml  - same parameters (machine-readable) including a
                                  data_types map for the Mantis SDK
                                  (mantis.create_space(df, data_types=...)).

Reproduce:  python -m scripts.prepare_mantis_input
"""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data_store"
OUT = Path.home() / "Downloads"

CSV_OUT = OUT / "B01_mantis_space.csv"
MD_OUT = OUT / "B01_mantis_space_config.md"
YAML_OUT = OUT / "B01_mantis_space_config.yaml"

INDUCED_SEED = 13  # same 64-leaf induced-subtree protocol as the benchmark runs


def _depth(tree, node):
    d = 0
    while tree.nodes[node].parent is not None:
        node = tree.nodes[node].parent
        d += 1
    return d


def build_cells():
    from b01_clock.data.adapters.convexml import load_convexml_example
    from b01_clock.data.adapters.gestalt import GESTALT_TOTAL_TIME_DP, load_gestalt_newick

    bundle = load_convexml_example(DATA / "convexml", total_time=1.0, rep=0)
    ct = bundle.tree
    keep = sorted(np.random.default_rng(INDUCED_SEED).choice(
        sorted(ct.leaves()), size=64, replace=False).tolist())
    cti = ct.induced_subtree(keep)

    gt = load_gestalt_newick(
        str(DATA / "gestalt" / "GSE81713_fish_ADR1_PHYLIP_MIX_gte5_output.newick.txt.gz"),
        total_time=GESTALT_TOTAL_TIME_DP)
    gkeep = sorted(np.random.default_rng(INDUCED_SEED).choice(
        sorted(gt.leaves()), size=64, replace=False).tolist())
    gti = gt.induced_subtree(gkeep)

    rows = []
    for tree, subtree, dataset, sample, abbrev in [
        (ct, cti, "Cassiopeia-ConvexML", "convexml_default_tree_0", "convexml"),
        (gt, gti, "GESTALT zebrafish", "gestalt_ADR1_gte5", "gestalt"),
    ]:
        in_sub = set(subtree.leaves())
        for leaf in sorted(tree.leaves()):
            certified = leaf in in_sub
            rows.append({
                "cell_id": leaf,
                "tree_id": f"{abbrev}_tree",   # Mantis ingests one tree per space; keep constant per dataset
                "dataset": abbrev,
                "sample_family": sample if not certified else f"{sample}/certified_subtree",
                "analysis_subtree": "certified_64leaf" if certified else "full_tree",
                "depth_in_subtree": _depth(subtree, leaf) if certified else -1,
                "label": f"{dataset} cell {leaf} ({sample}{', certified subtree' if certified else ''})",
            })
    return rows


def write_csv(rows):
    headers = ["cell_id", "tree_id", "dataset", "sample_family",
               "analysis_subtree", "depth_in_subtree", "label"]
    with CSV_OUT.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=headers)
        w.writeheader()
        w.writerows(rows)
    return CSV_OUT


def write_configs():
    md = """Defaults use the Mantis vLLM embedding and instruct services.

Embedding
Method
llm
Service
vllm-embeddings
Model
auto
Batch Size
100

Reduction
Reduction Method
UMAP
End Dimension
2
N Neighbors
15
Min Dist
0.1

Clustering
Clustering Method
k_means
Metric
euclidean
Min Cluster Size
5

Optimal Depth Detection
Strategy
sse
Linkage Method
ward
N Clusters
auto

Labeling
Config Type
metadata
Service
vllm-instruct
Model
auto
Sampling Strategy
representative

Traversing Strategy
top_bottom

Data
Source
B01_mantis_space.csv
Label Column
label

Generate Landmarks
"""
    MD_OUT.write_text(md)

    yaml = f"""# Mantis Create-a-Space parameters for B01 lineage data
embedding:
  method: llm
  service: vllm-embeddings
  model: auto
  batch_size: 100
reduction:
  reduction_method: UMAP
  end_dimension: 2
  n_neighbors: 15
  min_dist: 0.1
clustering:
  clustering_method: k_means
  metric: euclidean
  min_cluster_size: 5
optimal_depth_detection:
  strategy: sse
  linkage_method: ward
  n_clusters: auto
labeling:
  config_type: metadata
  service: vllm-instruct
  model: auto
  sampling_strategy: representative
traversing_strategy: top_bottom
data:
  source: {CSV_OUT.name}
  label_column: label
data_types:
  cell_id: all_data
  label: labels
  dataset: category
  sample_family: category
  analysis_subtree: category
  depth_in_subtree: metadata
"""
    YAML_OUT.write_text(yaml)
    return MD_OUT, YAML_OUT


if __name__ == "__main__":
    rows = build_cells()
    write_csv(rows)
    write_configs()
    print(f"cells: {CSV_OUT} ({CSV_OUT.stat().st_size:,} B, {len(rows):,} rows)")
    print(f"config md: {MD_OUT}")
    print(f"config yaml: {YAML_OUT}")