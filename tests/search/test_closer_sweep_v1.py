"""Closer-sweep V1 acceptance gate — the NEW (2026-06-08) circular-coplanar rows.

The #137 free-return (radial-crossing) closer was already gated on the original 12
descriptor rows (``test_russell12_likeforlike_probe.py``: 8 CLOSE-AND-MATCH, 4
already V1). This file is its sibling for the rows ingested on 2026-06-08 — the 16
Russell 2004 Table 3.4 ``russell-ocampo-*`` cyclers and the 15 Rall 1970
``rall-1970-*`` free-fall orbits — none of which had been run through the closer.

The sweep (``docs/notes/2026-06-08-closer-sweep-v1-candidates.md``) found:

* 6 of the 16 Russell T3.4 rows CLOSE-AND-MATCH **and** pass the §14 V1 mechanics
  (closed, V_inf-continuous single free-return ellipse) — these are the new
  V1-promotable rows, pinned here.
* 2 more (4.7.1-2, 4.8.1+2) CLOSE-AND-MATCH on the emerged V_inf but FAIL the V1
  V_inf-continuity gate (the single ellipse does not close to Earth — multi-arc),
  pinned as CLOSE-AND-MATCH-but-not-V1 so the honesty gate cannot silently flip.
* 8 NO-CLOSE (long-transit / low aphelion-ratio rows whose aphelion+transit seed
  collapses to the e-floor — a single radial-crossing ellipse cannot reach Mars).
* All 15 Rall rows are NOT-REACHABLE: Rall Model I.B does not tabulate the per-arc
  ``aphelion_au``, so the closer cannot seed ``(a, e)``.

GOLDEN DISCIPLINE (project memory feedback_golden_tests_sourced_only): the EXPECTED
side of every match is the row's SOURCED V_inf anchor (Russell 2004 Table 3.4).
The closer's seed ``(a, e)`` is derived from the SOURCED aphelion + transit
(constraints); the per-body V_inf EMERGES from the converged ellipse and is the
EVIDENCE — never imposed. See the closer's own module docstring
(``cyclerfinder.search.free_return``) for the constraint-vs-evidence separation.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pytest
import yaml  # type: ignore[import-untyped]

from cyclerfinder.core.constants import DAYS_PER_JULIAN_YEAR, MU_SUN_KM3_S2
from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.search.free_return import _residuals as _fr_residuals
from cyclerfinder.search.free_return import (
    free_return_correct,
    seed_ae_from_aphelion_transit,
)
from cyclerfinder.search.free_return_v1 import free_return_v1_mechanics

REPO_ROOT = Path(__file__).resolve().parents[2]
DAY_S = 86400.0
MU = MU_SUN_KM3_S2

# Campaign defaults (mirrors scripts/campaign_russell12.py / data/discover.py) —
# never loosened per row.
TOL_VINF_KMS = 0.5
TOL_KMS = 0.1
# Dense t0 phase floor: deep-aphelion high-e rows have a narrow residual basin a
# coarse grid steps over (#137 Part 3). Pure residual evals (no Lambert), cheap.
PHASE_EPOCHS = 4096


def _row(rid: str) -> dict[str, Any]:
    rows = yaml.safe_load((REPO_ROOT / "data" / "catalogue.yaml").read_text())
    return next(r for r in rows if r["id"] == rid)


def _best_phase_t0(a_seed: float, e_seed: float, period_sec: float, ephem: Ephemeris) -> float:
    """Scan t0 over one period; return the phase minimising the free-return
    residual at the SOURCED ``(a, e)`` (mirrors the campaign / discover path)."""
    best_t0, best_res = 0.0, float("inf")
    for frac in np.linspace(0.0, 1.0, PHASE_EPOCHS, endpoint=False):
        t0 = float(frac) * period_sec
        try:
            res = _fr_residuals(
                np.array([a_seed, e_seed, t0]),
                period_days=period_sec / DAY_S,
                ephem=ephem,
                bodies=("E", "M"),
                mu=MU,
            )
        except Exception:
            continue
        m = max(abs(r) for r in res)
        if m < best_res:
            best_res, best_t0 = m, t0
    return best_t0


def _close_row(rid: str) -> tuple[Any, dict[str, float], dict[str, float]]:
    """Seed the closer from the row's SOURCED aphelion+transit, close the single
    free-return ellipse, return ``(sol, derived_vinf, sourced_vinf)``."""
    row = _row(rid)
    aphelion = float(row["orbit_elements"]["aphelion_au"])
    transit = row["invariants"]["transit_times_days"]
    sourced_vinf = {e["body"]: float(e["vinf_kms"]) for e in row["vinf_kms_at_encounters"]}
    period_sec = float(row["period"]["years"]) * DAYS_PER_JULIAN_YEAR * DAY_S
    ephem = Ephemeris("circular")
    a_seed, e_seed = seed_ae_from_aphelion_transit(aphelion, float(transit[0]), mu=MU)
    t0 = _best_phase_t0(a_seed, e_seed, period_sec, ephem)
    sol = free_return_correct(
        t0_seed_sec=t0,
        a_seed_au=a_seed,
        e_seed=e_seed,
        period_sec=period_sec,
        ephem=ephem,
        mu=MU,
        tol_kms=TOL_KMS,
    )
    derived = {k: float(v) for k, v in sol.vinf_kms.items()}
    return sol, derived, sourced_vinf


# ---------------------------------------------------------------------------
# V1-PROMOTABLE rows: CLOSE-AND-MATCH on the emerged V_inf AND pass the §14 V1
# mechanics (closed, V_inf-continuous single free-return ellipse). These are the
# acceptance tests a later V1 writeback will cite.
#
# EXPECTED side = the SOURCED Russell 2004 Table 3.4 (E, M) V_inf anchors.
# ---------------------------------------------------------------------------

V1_PROMOTABLE = [
    ("russell-ocampo-3.1.1+2", 5.4, 9.2),
    ("russell-ocampo-3.1.3+0", 5.1, 9.1),
    ("russell-ocampo-4.1.1-4", 5.5, 9.3),
    ("russell-ocampo-4.1.2-2", 5.2, 9.2),
    ("russell-ocampo-4.1.4-1", 5.1, 9.2),
    ("russell-ocampo-4.6.3+0", 6.4, 9.5),
]


@pytest.mark.slow
@pytest.mark.parametrize(("rid", "src_vinf_e", "src_vinf_m"), V1_PROMOTABLE)
def test_new_row_closes_matches_and_passes_v1(
    rid: str, src_vinf_e: float, src_vinf_m: float
) -> None:
    """The closer reaches residual ~ 0 at the sourced geometry, the EMERGED V_inf
    matches the INDEPENDENTLY sourced anchor (<= 0.5 km/s), and the §14 V1
    mechanics pass on the closed, V_inf-continuous reconstructed E->M->E arc."""
    sol, derived, sourced = _close_row(rid)

    # Closed at the corrector floor.
    assert sol.converged
    assert sol.max_residual_kms < TOL_KMS
    # EVIDENCE (derived, never imposed): emerged V_inf matches the sourced anchor.
    assert abs(derived["E"] - src_vinf_e) <= TOL_VINF_KMS
    assert abs(derived["M"] - src_vinf_m) <= TOL_VINF_KMS
    # Sanity: the parametrised anchor IS the row's sourced anchor (no drift).
    assert sourced["E"] == pytest.approx(src_vinf_e)
    assert sourced["M"] == pytest.approx(src_vinf_m)
    # §14 V1: closed, V_inf-continuous single free-return ellipse.
    period_sec = float(_row(rid)["period"]["years"]) * DAYS_PER_JULIAN_YEAR * DAY_S
    v1 = free_return_v1_mechanics(sol, Ephemeris("circular"), period_sec, mu=MU)
    assert v1.v1_passed, v1.detail


# ---------------------------------------------------------------------------
# CLOSE-AND-MATCH but NOT V1: the emerged V_inf matches the sourced anchor, yet a
# single free-return ellipse does NOT close to Earth (the reconstructed return
# leg breaks Mars V_inf continuity — multi-arc, needs phasing loops). Pinned so
# the honesty gate cannot silently flip these into the promotable set.
# ---------------------------------------------------------------------------

CLOSE_NOT_V1 = [
    ("russell-ocampo-4.7.1-2", 6.6, 11.4),
    ("russell-ocampo-4.8.1+2", 12.5, 10.7),
]


@pytest.mark.slow
@pytest.mark.parametrize(("rid", "src_vinf_e", "src_vinf_m"), CLOSE_NOT_V1)
def test_new_row_matches_but_fails_v1_continuity(
    rid: str, src_vinf_e: float, src_vinf_m: float
) -> None:
    """Emerged V_inf matches the sourced anchor (CLOSE-AND-MATCH), but the single
    ellipse breaks Mars V_inf continuity -> V1 mechanics REJECT it (multi-arc)."""
    sol, derived, _ = _close_row(rid)

    assert sol.converged
    assert abs(derived["E"] - src_vinf_e) <= TOL_VINF_KMS
    assert abs(derived["M"] - src_vinf_m) <= TOL_VINF_KMS

    period_sec = float(_row(rid)["period"]["years"]) * DAYS_PER_JULIAN_YEAR * DAY_S
    v1 = free_return_v1_mechanics(sol, Ephemeris("circular"), period_sec, mu=MU)
    # The honesty gate: V1 NOT awarded; the failure is the continuity break.
    assert not v1.v1_passed
    assert not v1.vinf_continuous


# ---------------------------------------------------------------------------
# NO-CLOSE rows: the aphelion+transit seed collapses to the e-floor (a single
# radial-crossing ellipse cannot reach Mars at this geometry). Pinned as a
# negative gate — the continuation campaign cannot seed these from a single
# ellipse.
# ---------------------------------------------------------------------------

NO_CLOSE = [
    "russell-ocampo-3.5.1+2",
    "russell-ocampo-4.1.1-6",
    "russell-ocampo-4.1.2-3",
    "russell-ocampo-4.6.1-4",
    "russell-ocampo-4.8.1+3",
    "russell-ocampo-4.9.1-3",
    "russell-ocampo-4.10.1-3",
    "russell-ocampo-4.12.1-2",
]


@pytest.mark.slow
@pytest.mark.parametrize("rid", NO_CLOSE)
def test_new_row_does_not_close_as_single_ellipse(rid: str) -> None:
    """A single free-return ellipse cannot represent these long-transit rows: the
    closer does not reach the residual floor (NOT a seeding gap — the geometry)."""
    sol, _, _ = _close_row(rid)
    assert not sol.converged


# ---------------------------------------------------------------------------
# NOT-REACHABLE rows: the Rall 1970 free-fall orbits (Model I.B) do not tabulate
# the per-arc aphelion, so the closer has nothing to seed (a, e) from. This is a
# DATA fact, not a closer run — pinned so the NOT-REACHABLE set is explicit and
# the continuation campaign knows which rows need Appendix C state.
# ---------------------------------------------------------------------------

RALL_NOT_REACHABLE = [
    "rall-1970-m4-1",
    "rall-1970-m6-1",
    "rall-1970-m6-2",
    "rall-1970-m6-3",
    "rall-1970-m4-1a",
    "rall-1970-m5-1a",
    "rall-1970-m5-1b",
    "rall-1970-m5-1c",
    "rall-1970-m5-1d",
    "rall-1970-m5-1e",
    "rall-1970-m5-2a",
    "rall-1970-m5-2b",
    "rall-1970-m5-2c",
    "rall-1970-m5-2d",
    "rall-1970-m5-2e",
]


@pytest.mark.parametrize("rid", RALL_NOT_REACHABLE)
def test_rall_row_is_not_closer_reachable(rid: str) -> None:
    """Rall Model I.B rows carry no per-arc aphelion -> the free-return closer
    cannot seed (a, e). The row is NOT-REACHABLE by the circular-coplanar closer."""
    row = _row(rid)
    aphelion = (row.get("orbit_elements") or {}).get("aphelion_au")
    assert aphelion is None
