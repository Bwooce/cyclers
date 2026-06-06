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


# ---------------------------------------------------------------------------
# #137 ACCEPTANCE GATE — the reworked free-return (radial-crossing) genome.
#
# The free-Lambert genome above makes truth NOT a closure point (WALKED-AWAY).
# The free-return genome (cyclerfinder.search.free_return) expresses the Mars
# transfer as a free-return arc on a single heliocentric ellipse; the per-body
# V_inf and leg ToFs EMERGE from the ellipse shape. Seeding (a, e) at the SOURCED
# S1L1 ellipse (a=1.30 AU, e=0.257; Rogers 2012, the SAME physical cycler per
# docs/notes/multi-arc-classification.md §12) at the best phase yields residual
# ~ 0 -- the sourced geometry IS now representable.
#
# GOLDEN DISCIPLINE: (a, e) is the SOURCED input (constraint); the asserted
# EVIDENCE is the EMERGED V_inf compared to the INDEPENDENTLY sourced V_inf
# anchors (Russell 4.99/5.10, McConaghy 4.7/5.0). V_inf is never imposed.
# ---------------------------------------------------------------------------

# Sourced S1L1 heliocentric ellipse (Rogers 2012 Table 1) — the constraint side.
_S1L1_A_AU = 1.30
_S1L1_E = 0.257


def _best_phase_t0(corr_mod: ModuleType, period_sec: float) -> float:
    """Scan t0 over one period; return the phase minimising the free-return
    residual at the SOURCED (a, e). Mirrors the probe's best-phase selection."""
    import numpy as np

    from cyclerfinder.core.ephemeris import Ephemeris

    ephem = Ephemeris("circular")
    best_t0, best_res = 0.0, float("inf")
    for frac in np.linspace(0.0, 1.0, 360, endpoint=False):
        t0 = float(frac) * period_sec
        res = corr_mod._residuals(
            np.array([_S1L1_A_AU, _S1L1_E, t0]),
            period_days=period_sec / 86400.0,
            ephem=ephem,
            bodies=("E", "M"),
            mu=132712440018.0,
        )
        m = max(abs(r) for r in res)
        if m < best_res:
            best_res, best_t0 = m, t0
    return best_t0


@pytest.mark.slow
@pytest.mark.parametrize(
    ("rid", "src_vinf_e", "src_vinf_m"),
    [
        ("mcconaghy-2006-em-k2", 4.7, 5.0),
        ("russell-ch4-4.991gG2", 4.99, 5.10),
    ],
)
def test_free_return_genome_makes_truth_a_closure_point(
    rid: str, src_vinf_e: float, src_vinf_m: float
) -> None:
    """#137 acceptance gate: the free-return genome reaches residual ~ 0 at the
    SOURCED geometry on the symmetric k=2 rows, and the EMERGED V_inf matches the
    INDEPENDENTLY sourced anchor (within the 0.5 km/s campaign tolerance)."""
    from cyclerfinder.core.ephemeris import Ephemeris
    from cyclerfinder.search import free_return as fr

    row = _row(rid)
    period_sec = float(row["period"]["years"]) * 365.25 * 86400.0
    t0 = _best_phase_t0(fr, period_sec)

    result = fr.free_return_correct(
        t0_seed_sec=t0,
        a_seed_au=_S1L1_A_AU,
        e_seed=_S1L1_E,
        period_sec=period_sec,
        ephem=Ephemeris("circular"),
        tol_kms=0.1,
    )

    # The decisive #137 gate: truth IS a residual-zero point now.
    assert result.converged
    assert result.max_residual_kms < 0.1
    # EVIDENCE (derived, not imposed): emerged V_inf matches the sourced anchor.
    assert abs(result.vinf_kms["E"] - src_vinf_e) <= 0.5
    assert abs(result.vinf_kms["M"] - src_vinf_m) <= 0.5
    # The emerged transfer ToF lands on the sourced ~150-153 d transit.
    assert abs(result.transfer_tof_days - 153.0) <= 5.0
