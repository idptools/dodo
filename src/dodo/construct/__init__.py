"""Structure rebuilding: target dimensions, IDRs, loops, domain placement, backbone."""

from __future__ import annotations

from .backbone_refine import refine_backbone
from .ca_backbone import add_backbone_to_rebuilt, backbone_from_ca
from .dimensions import (
    DimensionTarget,
    albatross_available,
    predict_end_to_end,
    predict_radius_of_gyration,
    predict_scaling_exponent,
    target_dimensions,
)
from .pipeline import RebuildReport, RegionOutcome, build_from_sequence, rebuild
from .place import DomainPlacement, reposition_folded_domains

__all__ = [
    "DimensionTarget",
    "DomainPlacement",
    "RebuildReport",
    "RegionOutcome",
    "add_backbone_to_rebuilt",
    "albatross_available",
    "backbone_from_ca",
    "build_from_sequence",
    "predict_end_to_end",
    "predict_radius_of_gyration",
    "predict_scaling_exponent",
    "rebuild",
    "refine_backbone",
    "reposition_folded_domains",
    "target_dimensions",
]
