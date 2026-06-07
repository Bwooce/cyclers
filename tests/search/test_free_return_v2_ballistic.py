"""§14 V2-ballistic re-run — the four #137 free-return V1 rows (finding, slow).

The §14 V2 class-split amendment (2026-06-07) defines **V2-ballistic** as >=3
continuous laps with bounded rotating-frame drift, evaluated IN THE ROW'S DEFINING
MODEL (the circular ephemeris for these circular-coplanar rows). This module
RE-RUNS that gate on the four free-return rows the #137 work promoted to V1
(``russell-ch4-{5.30gGf3, 9.94Gg3, 5.75ggF3, 9.353Gg2}``) and records the
outcome.

THE FINDING (no promotion): the four rows do NOT pass V2-ballistic.
----------------------------------------------------------------------
These rows are ``cycler_class: multi-arc``. The #137 V1 evidence closes a single
heliocentric E->M->E free-return **ellipse slice** whose radial crossings
reproduce the sourced V_inf. But that slice spans only ~0.3-0.4 of the catalogue
cycler period (the E->M->E arc is ~620-708 d; the catalogue period is 4.27 yr
[k=2] or 6.41 yr [k=3], i.e. 2-3 E-M synodics). The FULL cycler includes the
Earth-to-Earth resonant phasing intervals between free-return passes (e.g. the
3:2 full-rev return descriptor ``f(3:2,...)``) that the single ellipse does NOT
represent. There is therefore no continuous >=3-lap trajectory to propagate for
these objects — the V2-ballistic gate is structurally inapplicable to a
single-arc slice of a multi-arc cycler.

Measured (NOT assumed): propagating the reconstructed single-ellipse arc over 3
laps (each lap = the catalogue cycler period) under the rotating-frame-repeat
metric gives a max drift of ~9.4e7-1.2e8 km — far above the M6a idealised
50,000 km tolerance — precisely because each lap re-pins the leg starts to the
incommensurately-phased planets (the same mechanism that makes the cross-cycle
metric the wrong instrument for the powered Aldrin, #134). This is a finding
about the object's structure (a multi-arc slice), not a numerical near-miss. The
four rows stay V1.

Discipline: ``slow`` (free-return correction + multi-lap propagation on the
circular ephemeris). The assertion records the qualitative finding with teeth
(drift far above the M6a tolerance) without pinning a brittle magnitude.
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
from cyclerfinder.search.free_return_v1 import build_free_return_cycler
from cyclerfinder.verify.propagate import DRIFT_TOLERANCE_KM, verify_long_term_stability

REPO_ROOT = Path(__file__).resolve().parents[2]
CAMPAIGN = REPO_ROOT / "scripts" / "campaign_russell12.py"
DAY_S = 86400.0

# The four #137 free-return rows promoted to V1 (the V2-ballistic candidates).
V1_ROWS = (
    "russell-ch4-5.30gGf3",
    "russell-ch4-9.94Gg3",
    "russell-ch4-5.75ggF3",
    "russell-ch4-9.353Gg2",
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
def test_free_return_rows_do_not_pass_v2_ballistic(rid: str) -> None:
    """FINDING (no promotion): each free-return V1 row is a single-ellipse slice of
    a multi-arc cycler, so the >=3-lap continuous-periodicity gate does not close
    in its defining (circular) model — drift far above the M6a tolerance. The row
    stays V1."""
    ephem, period_sec, sol = _close(rid)
    assert sol.converged, f"{rid}: free-return corrector did not close"
    # The row is multi-arc: the single E->M->E ellipse spans only a fraction of
    # the catalogue cycler period (the rest is Earth-to-Earth resonant phasing
    # the single ellipse does not represent).
    assert _row(rid)["cycler_class"] == "multi-arc"

    cycler = build_free_return_cycler(sol, ephem, period_sec)
    report = verify_long_term_stability(
        cycler,
        n_laps=3,  # spec §14 V2: >= 3 continuous laps
        ephem=ephem,
        t_start=sol.t0_sec,
        cycler_id=rid,
    )
    assert report.n_laps_propagated == 3
    # NOT bounded: the single-arc slice has no continuous 3-lap trajectory; the
    # rotating-frame-repeat metric re-pins to incommensurately-phased planets.
    assert report.stable is False
    assert report.max_drift_km > 10.0 * DRIFT_TOLERANCE_KM, report.max_drift_km
