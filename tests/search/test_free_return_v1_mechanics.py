"""#137 Part 1 — §14 V1 mechanics evidence for the free-return genome (slow).

The free-return corrector (:mod:`cyclerfinder.search.free_return`) closes a single
heliocentric ellipse whose radial crossings reproduce a sourced Russell cycler's
DERIVED V_inf (the #137 breakthrough). This test applies the LITERAL §14 V1
mechanics to the closed geometry, like-for-like on the CIRCULAR ephemeris (a
circular-coplanar reproduction of a circular-coplanar source):

* path (a) — every leg re-solved with lamberthub izzo2015 + gooding1990, agreement
  < V1_TOLERANCE_MPS;
* path (c) — Kepler forward re-propagation residual < KEPLER_REPROP_TOL_KM;

both reused verbatim from :func:`cyclerfinder.verify.agreement.crosscheck_code_paths`.
A genome honesty gate additionally requires Mars-flyby V_inf continuity, so the
six rows whose single free-return ellipse does NOT close to Earth (their return
needs intermediate phasing loops — a ~24 km/s forced-Lambert discontinuity) are
correctly REFUSED V1.

GOLDEN DISCIPLINE: the (a, e) seed is the SOURCED input (constraint); the asserted
EVIDENCE is the §14 mechanics passing on the reconstructed arc (code-path
consistency + forward-propagation residual) — no published value is rediscovered
here. These three rows are the catalogue ``validation_level: V1`` evidence pointed
to by ``_LEVEL_EVIDENCE`` (mirrors the Aldrin V1 pattern). The other rows' refusal
is asserted too, so the gate keeps teeth against silent over-promotion.

See ``docs/notes/2026-06-07-russell12-freereturn-results.md`` (Part 1).
"""

from __future__ import annotations

import importlib.util
from pathlib import Path
from types import ModuleType
from typing import Any

import numpy as np
import pytest
import yaml  # type: ignore[import-untyped]

from cyclerfinder.core.constants import DAYS_PER_JULIAN_YEAR, MU_SUN_KM3_S2
from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.search.free_return import free_return_correct
from cyclerfinder.search.free_return_v1 import free_return_v1_mechanics
from cyclerfinder.verify.crosscheck import V1_TOLERANCE_MPS

REPO_ROOT = Path(__file__).resolve().parents[2]
CAMPAIGN = REPO_ROOT / "scripts" / "campaign_russell12.py"
DAY_S = 86400.0

# The rows whose single free-return ellipse forms a genuinely closed,
# V_inf-continuous E->M->E cycler. The f/F full-Mars-radius free returns plus the
# deep-aphelion 9.353Gg2 (promoted CLOSE-AND-MATCH by #137 Part 3's dense phase
# scan, then found to clear the V_inf-continuity gate too).
#
# 8.049gGf2 added here after the #200/#205 Lambert-accuracy fixes: its reconstructed
# free-return arc now closes V_inf-continuously (continuity 103 m/s — in-family with
# the accepted rows, which span 0.9-190.5 m/s; 5.75ggF3 is looser at 190.5), and the
# row is independently catalogue-validated to V3. The original #137 "multi-arc,
# refused V1" classification was superseded by the improved reconstruction; this is a
# census re-baseline to match ground truth, not a tolerance relaxation.
V1_ROWS = (
    "russell-ch4-5.30gGf3",
    "russell-ch4-9.94Gg3",
    "russell-ch4-5.75ggF3",
    "russell-ch4-9.353Gg2",
    "russell-ch4-8.049gGf2",
)
# Matched rows whose single free-return ellipse does NOT close to Earth (multi-arc;
# return needs phasing loops) — refused V1 on the V_inf-continuity gate — plus the
# Lambert-singular deep-aphelion row.
NON_V1_ROWS = (
    "mcconaghy-2006-em-k2",
    "russell-ch4-4.991gG2",
    "russell-ch4-3.64gGg3",
    "russell-ch4-3.78Gg3",
    "russell-ch4-3.66gfF3",
    "russell-ch4-5.30ggF3",
)


def _load_campaign() -> ModuleType:
    spec = importlib.util.spec_from_file_location("campaign_russell12", CAMPAIGN)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _row(rid: str) -> dict[str, Any]:
    rows = yaml.safe_load((REPO_ROOT / "data" / "catalogue.yaml").read_text())
    return next(r for r in rows if r["id"] == rid)


def _close(rid: str) -> tuple[Ephemeris, float, Any]:
    """Close the free-return geometry for *rid* on the circular ephemeris."""
    camp = _load_campaign()
    row = _row(rid)
    aphelion = row["orbit_elements"]["aphelion_au"]
    transit = row["invariants"]["transit_times_days"]
    period_sec = float(row["period"]["years"]) * DAYS_PER_JULIAN_YEAR * DAY_S
    a_seed, e_seed = camp._seed_ae_from_aphelion_transit(float(aphelion), float(transit[0]))
    ephem = Ephemeris("circular")
    best_t0, best_res = 0.0, float("inf")
    # Dense phase floor (campaign FR_PHASE_EPOCHS_FLOOR): the deep-aphelion high-e
    # rows have a narrow t0 residual basin a coarse grid steps over (#137 Part 3).
    for frac in np.linspace(0.0, 1.0, camp.FR_PHASE_EPOCHS_FLOOR, endpoint=False):
        t0 = float(frac) * period_sec
        res = camp._fr_residuals(
            np.array([a_seed, e_seed, t0]),
            period_days=period_sec / DAY_S,
            ephem=ephem,
            bodies=("E", "M"),
            mu=MU_SUN_KM3_S2,
        )
        m = max(abs(r) for r in res)
        if m < best_res:
            best_res, best_t0 = m, t0
    sol = free_return_correct(
        t0_seed_sec=best_t0,
        a_seed_au=a_seed,
        e_seed=e_seed,
        period_sec=period_sec,
        ephem=ephem,
        tol_kms=0.1,
    )
    return ephem, period_sec, sol


@pytest.mark.slow
@pytest.mark.parametrize("rid", V1_ROWS)
def test_free_return_rows_pass_section14_v1_mechanics(rid: str) -> None:
    """The three closing rows pass §14 V1: lamberthub agreement < V1_TOLERANCE_MPS
    AND Kepler forward re-propagation < 1 km, on a V_inf-continuous closed arc."""
    ephem, period_sec, sol = _close(rid)
    assert sol.converged, f"{rid}: free-return corrector did not close"
    v = free_return_v1_mechanics(sol, ephem, period_sec)
    assert v.built, v.detail
    # path (a): in-house Lambert vs lamberthub.
    assert v.lamberthub_passed, f"{rid}: lamberthub {v.lamberthub_max_diff_mps} m/s"
    assert v.lamberthub_max_diff_mps < V1_TOLERANCE_MPS
    # path (c): forward Kepler re-propagation residual.
    assert v.kepler_reprop_passed, f"{rid}: kepler reprop {v.kepler_reprop_max_residual_km} km"
    assert v.kepler_reprop_max_residual_km < 1.0
    # the genome honesty gate: the reconstructed arc is a genuinely closed,
    # V_inf-continuous cycler (not a vacuous forced-return).
    assert v.vinf_continuous, f"{rid}: Mars V_inf discontinuity {v.vinf_continuity_kms} km/s"
    assert v.v1_passed


@pytest.mark.slow
@pytest.mark.parametrize("rid", NON_V1_ROWS)
def test_multi_arc_rows_are_refused_v1(rid: str) -> None:
    """The matched rows whose single free-return ellipse does not close to Earth
    are correctly REFUSED V1 — either Lambert-singular on reconstruction, or the
    forced return leg breaks Mars V_inf continuity (multi-arc; needs phasing
    loops). Keeps the V1 gate honest against silent over-promotion."""
    ephem, period_sec, sol = _close(rid)
    v = free_return_v1_mechanics(sol, ephem, period_sec)
    assert not v.v1_passed, f"{rid}: unexpectedly passed V1 — {v.detail}"
    # Either the reconstruction is Lambert-singular, or V_inf continuity is broken
    # by far more than the tolerance (a genuinely multi-arc trajectory).
    assert (not v.built) or (not v.vinf_continuous)
