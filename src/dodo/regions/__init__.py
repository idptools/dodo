"""Identification of folded domains, IDRs and loops."""

from __future__ import annotations

from .contact import ContactProfile, contact_profile, is_loop_like, loop_contact_counts
from .identify import (
    RegionAssignment,
    Strategy,
    assign_regions,
    assign_regions_from_spec,
    find_runs,
    merge_blocks,
)

__all__ = [
    "ContactProfile",
    "RegionAssignment",
    "Strategy",
    "assign_regions",
    "assign_regions_from_spec",
    "contact_profile",
    "find_runs",
    "is_loop_like",
    "loop_contact_counts",
    "merge_blocks",
]
