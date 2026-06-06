"""M-3D Phase 4: schema accepts circular-inclined fidelity tags (plan §4)."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

SCHEMA: dict[str, Any] = json.loads(Path("data/catalogue.schema.json").read_text())


def _enum(field: str) -> list[Any]:
    enum: list[Any] = SCHEMA["items"]["properties"][field]["enum"]
    return enum


def test_orbit_and_vinf_fidelity_enums_include_circular_inclined() -> None:
    assert "circular-inclined" in _enum("orbit_fidelity")
    assert "circular-inclined" in _enum("vinf_fidelity")
    # existing tiers preserved
    for tier in ("circular-coplanar", "analytic-ephemeris", "real-de440", None):
        assert tier in _enum("orbit_fidelity")
