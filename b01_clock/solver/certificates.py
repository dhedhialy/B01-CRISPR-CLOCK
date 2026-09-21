"""Identification status and certificate packaging."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from enum import Enum
from typing import Any, Dict, List, Optional


class IdentificationStatus(str, Enum):
    UNIDENTIFIED = "unidentified"
    PARTIALLY_IDENTIFIED = "partially_identified"
    POINT_IDENTIFIED = "point_identified"
    INFEASIBLE = "infeasible"
    ABSTAIN = "abstain"


def classify_identification(
    lower: float,
    upper: float,
    total_time: float,
    point_tol: float = 1e-6,
    unidentified_frac: float = 0.99,
) -> IdentificationStatus:
    width = upper - lower
    if not (upper >= lower):
        return IdentificationStatus.INFEASIBLE
    if width <= point_tol * max(total_time, 1.0):
        return IdentificationStatus.POINT_IDENTIFIED
    if width >= unidentified_frac * total_time:
        return IdentificationStatus.UNIDENTIFIED
    return IdentificationStatus.PARTIALLY_IDENTIFIED


@dataclass
class CertificateReport:
    method: str
    total_time: float
    rate_bounds_summary: dict
    results: List[dict]
    topology_id: str = "default"
    seed: Optional[int] = None
    notes: List[str] = field(default_factory=list)

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)

    def narrowest(self) -> Optional[dict]:
        if not self.results:
            return None
        return min(self.results, key=lambda r: r["upper"] - r["lower"])
