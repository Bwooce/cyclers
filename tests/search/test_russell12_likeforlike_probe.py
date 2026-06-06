"""#135 like-for-like seed-at-truth probe: pins the decisive diagnostic.

The Russell-12 closure campaign run on ``Ephemeris('circular')`` (coplanar vs
coplanar, against rows that are by-construction solutions OF the coplanar model)
still produces 0 CLOSE-AND-MATCH. The seed-at-truth probe explains why: the
row's OWN sourced ToF geometry is NOT a residual-zero point of our descriptor->
corrector genome -- the residual evaluated AT truth sits far above the 0.1 km/s
closure floor, so the corrector correctly walks away from it.

This test pins that diagnosis on the representative symmetric row
(``mcconaghy-2006-em-k2``, 153/153 d). It is NON-GOLDEN: the asserted numbers are
OUR residual at the SOURCED geometry, pinned to guard the genome-mapping finding
against silent regression -- NOT a rediscovery of any published value. See
``docs/notes/2026-06-06-russell12-likeforlike.md``.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path
from types import ModuleType
from typing import Any

import pytest
import yaml  # type: ignore[import-untyped]

REPO_ROOT = Path(__file__).resolve().parents[2]
CAMPAIGN = REPO_ROOT / "scripts" / "campaign_russell12.py"


def _load_campaign() -> ModuleType:
    spec = importlib.util.spec_from_file_location("campaign_russell12", CAMPAIGN)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _row(rid: str) -> dict[str, Any]:
    rows = yaml.safe_load((REPO_ROOT / "data" / "catalogue.yaml").read_text())
    return next(r for r in rows if r["id"] == rid)


@pytest.mark.slow
def test_sourced_geometry_is_not_a_circular_closure_point() -> None:
    """Seed-at-truth on circular: residual AT the sourced 153/153 d geometry is
    far above the closure floor -> the corrector WALKS AWAY (not a seeding gap)."""
    mod = _load_campaign()
    row = _row("mcconaghy-2006-em-k2")
    probe = mod.probe_at_truth(row, phase_epochs=64, model="circular")

    # The diagnosis: truth is not a closure point of our genome.
    assert probe["verdict"] == "WALKED-AWAY"
    assert not probe["truth_residual_below_floor"]
    # Residual at the sourced geometry is multi-km/s, far above the 0.1 floor
    # (pinned ~3.4 km/s; assert the regime, not a sourced value).
    assert probe["best_phase_truth_residual_kms"] > 1.0
    # The corrector does find a genuine closure -- just at a different geometry.
    assert probe["solved_converged"]
    assert probe["tof_drift_days"] > mod.TOL_TRANSIT_DAYS


@pytest.mark.slow
def test_genome_maps_truth_seed_to_sourced_tofs() -> None:
    """The truth seed IS the row's sourced geometry (segment + descriptor ToFs),
    so a WALKED-AWAY verdict implicates the genome mapping, not the seed."""
    mod = _load_campaign()
    row = _row("mcconaghy-2006-em-k2")
    genome = mod.build_genome(row)
    truth_free = mod._truth_seed(genome)
    # mcconaghy is symmetric 153/153 d with the E-E loop eliminated as slack.
    assert truth_free == [153.0, 153.0]
