"""Checksums and provenance records."""

from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional


def sha256_file(path: str | Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


@dataclass
class ProvenanceRecord:
    name: str
    source_url: str
    accession: str
    local_path: str
    sha256: str
    downloaded_at: str
    recorder_architecture: str = ""
    n_targets: Optional[int] = None
    experimental_duration: Optional[float] = None
    topology_provenance: str = "unknown"  # observed|inferred|supplied|simulated
    notes: List[str] = field(default_factory=list)
    extra: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict:
        return asdict(self)


def write_provenance(record: ProvenanceRecord, out_path: str | Path) -> None:
    Path(out_path).write_text(json.dumps(record.to_dict(), indent=2))


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()
