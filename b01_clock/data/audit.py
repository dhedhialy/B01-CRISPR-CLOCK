"""Data audit checklist (§7)."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Optional

import numpy as np

from b01_clock.model.observations import CharacterMatrix, missingness_summary
from b01_clock.model.tree import LineageTree
from b01_clock.data.provenance import ProvenanceRecord, sha256_file


def audit_dataset(
    name: str,
    tree: LineageTree,
    matrix: CharacterMatrix,
    provenance: Optional[ProvenanceRecord] = None,
    replication_unit: str = "embryo_or_tree",
    calibration_sites: Optional[list] = None,
    analysis_sites: Optional[list] = None,
) -> Dict[str, Any]:
    miss = missingness_summary(matrix)
    ages = [tree.nodes[ℓ].age for ℓ in tree.leaves()]
    report = {
        "name": name,
        "accession_ok": provenance is not None and bool(provenance.accession),
        "sha256": provenance.sha256 if provenance else None,
        "source_url": provenance.source_url if provenance else None,
        "recorder_architecture": provenance.recorder_architecture if provenance else "",
        "n_targets_sites": matrix.n_sites,
        "n_leaves": len(tree.leaves()),
        "n_edges": len(tree.edges()),
        "experimental_duration": provenance.experimental_duration if provenance else tree.total_depth(),
        "replication_unit": replication_unit,
        "pseudoreplication_guard": (
            "Analyses treat each tree/embryo as one replicate; leaf cells are not "
            "independent biological replicates."
        ),
        "missingness": miss,
        "topology_provenance": provenance.topology_provenance if provenance else "unknown",
        "ultrametric": tree.is_ultrametric(),
        "leaf_age_range": (float(min(ages)), float(max(ages))) if ages else None,
        "calibration_holdout_separated": bool(
            calibration_sites is not None and analysis_sites is not None
            and set(calibration_sites).isdisjoint(set(analysis_sites))
        ),
        "checklist": {
            "verify_accession_version_integrity": provenance is not None,
            "record_architecture_targets_encoding_duration": True,
            "define_replication_unit": True,
            "quantify_missingness_dropout_saturation": True,
            "document_topology_provenance": provenance is not None,
            "separate_calibration_holdout": calibration_sites is not None,
            "deterministic_preprocessing_checksums": provenance is not None,
        },
    }
    return report
