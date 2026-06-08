"""Per-row V3 evidence gate over the continuation V1->V3 lift set (slow, DE440).

This mirrors ``tests/verify/test_aldrin_continuation_v3.py`` (#161) but generalised
across the whole lift set the #158 continuation can reach: the four Russell V1
free-return rows + Aldrin (``cyclerfinder.search.continuation_batch.LIFT_SET_RIDS``).
For each row the spec §14 **V3** evidence is the conjunction of two halves:

* **(a) phase-match** — the continued solution closes BALLISTICALLY at the true
  ephemeris (DE440), its emerged per-body V_inf matches the INDEPENDENTLY sourced
  anchor (the row's physical fingerprint), and its ``t0`` lands inside the sourced
  launch window (§14 gate). This is the batch driver's ``reaches_v3_closure``.
* **(b) bounded horizon TCM** — from the continued ``t0`` (the phase-matched real
  window), the ephemeris-mode horizon TCM over 3 laps (re-phasing one E-M synodic
  per lap, re-solving the in-family maintenance ΔV on DE440) is BOUNDED and every
  per-cycle ΔV stays within the engineering plausibility bar
  (``MAINTENANCE_DV_CONVENTION_KMS`` = 3.0 km/s, ``verify/plausibility.py``).

Sourced anchors are EXPECTED; our continued / maintenance values are EVIDENCE. The
gate is qualitative-with-teeth (sourced-floor discipline), NEVER an invented
absolute golden — the plausibility bar is the project's own CONVENTION, not a
sourced per-row maintenance budget (none is published; catalogue ``data_gaps``).

Lap cap (logged, no silent cap)
-------------------------------
Capped at :data:`N_LAPS` = 3 laps (the documented floor the task permits) for the
per-row wall budget — the per-lap DE440 maintenance re-solve is ~10-30 s and each
row also reproduces the cheap ``ladder=(3,)`` continuation. 3 laps is the §14 V3
minimum ("3-5 laps"); the un-run laps 4-5 are NOT silently dropped — this docstring
and the results note record the cap and the laps-4-5 behaviour observed in the
diagnostic (the 2009-2012 phasing window pushes several rows above the bar at laps
3-4, recorded in docs/notes/2026-06-08-continuation-batch-results.md).

Model-fidelity honesty
----------------------
Ephemeris-positions two-body + DE440 propagation = V3-class fidelity, NOT n-body
(the same class as #161). The maintenance chain re-solves Lambert/free-return legs
on DE440 planet states; it is the precursor to, not a substitute for, the planned
n-body harness (V4).

Break points are STRICT-XFAIL: a row that fails either half records its per-step
residuals / per-lap ΔV and is marked ``xfail(strict=True)`` — NEVER softened. The
day a break point passes, the marker must be removed (the row gained the missing
representability). No catalogue writeback is performed here (the main session
consolidates writebacks).
"""

from __future__ import annotations

from dataclasses import dataclass
from datetime import UTC, datetime, timedelta

import pytest
from _pytest.mark.structures import ParameterSet

from cyclerfinder.core.constants import SECONDS_PER_DAY
from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.search import continuation_batch as cb
from cyclerfinder.search.maintain import (
    MaintenanceOptimResult,
    optimise_aldrin_maintenance_dv,
    optimise_maintenance_dv,
)
from cyclerfinder.search.phase_match import PhaseSignature
from cyclerfinder.search.resonance import synodic_period_days
from cyclerfinder.verify.plausibility import MAINTENANCE_DV_CONVENTION_KMS
from cyclerfinder.verify.real_closure import _resolve_real_t_start

_J2000 = datetime(2000, 1, 1, 12, 0, 0, tzinfo=UTC)

# E-M synodic = the per-lap re-phasing cadence (the #161 Aldrin convention).
_SYNODIC_SEC = synodic_period_days("E", "M") * SECONDS_PER_DAY

# Lap cap (logged): the §14 V3 floor. See the module docstring for the laps-4-5
# note (un-run, not silently dropped).
N_LAPS = 3

# Per-row outbound (E->M) leg ToF guess (days), sourced from the catalogue's
# outbound transit (Russell tables) / the Aldrin 146 d. Drives the maintenance
# re-solve seed and the phase-signature leg duration.
_EM_TOF_DAYS: dict[str, float] = {
    "russell-ch4-5.30gGf3": 118.0,
    "russell-ch4-9.94Gg3": 82.0,
    "russell-ch4-5.75ggF3": 111.0,
    "russell-ch4-9.353Gg2": 85.0,
    "aldrin-classic-em-k1-outbound": 146.0,
}
_ME_TOF_DAYS_GUESS = 634.0
_TOF1_BOUNDS = (100.0, 250.0)
_TOF2_BOUNDS = (400.0, 900.0)
_EARTH_FLYBY_ALT_KM = 200.0


@pytest.fixture(scope="module")
def astropy_ephem() -> Ephemeris:
    return Ephemeris(model="astropy")


def _continued(seed: cb.BatchRowSeed, ephem: Ephemeris) -> cb.BatchRowResult:
    """The batch driver's continued result for one row (winning rung, cheap path).

    Uses ``ladder=(3,)`` (the documented winning rung — rng-free, deterministic,
    reproduces the full-ladder solution far cheaper) and the basin sweep that
    selects the V_inf-matched in-window basin.
    """
    return cb.run_continuation_for_seed(seed, ladder=(3,), final_ephemeris=ephem)


def _maintenance_at(
    rid: str, seed: cb.BatchRowSeed, pdate: datetime, ephem: Ephemeris
) -> MaintenanceOptimResult:
    """One lap's in-family maintenance re-solve on DE440 near ``pdate``.

    Aldrin uses its established resolver wrapper (the #161 path). The Russell rows
    use the body-agnostic optimiser seeded from the per-lap V_inf-matched DE440
    basin (resolved with the row's sourced V_inf fingerprint over a narrow window),
    closure_body="E" — the closest generic analog of the Aldrin path.
    """
    if rid == "aldrin-classic-em-k1-outbound":
        return optimise_aldrin_maintenance_dv(ephem, real_window_priority_date=pdate)
    tof1 = _EM_TOF_DAYS[rid]
    sig = PhaseSignature(
        bodies=("E", "M"),
        leg_durations_s=(tof1 * SECONDS_PER_DAY,),
        vinf_target_kms=(seed.sourced_vinf["E"], seed.sourced_vinf["M"]),
    )
    t0_guess = _resolve_real_t_start(sig, ephem, pdate, window_years=1.5, mismatch_cap_kms=20.0)
    if t0_guess is None:
        # No V_inf-matched basin near this lap: surface as a non-converged result
        # so the gate records the break honestly (never a silent substitution).
        t0_guess = (pdate - _J2000).total_seconds()
    return optimise_maintenance_dv(
        ["E", "M", "E"],
        ephem,
        t0_guess_sec=t0_guess,
        tof_days_guesses=(tof1, _ME_TOF_DAYS_GUESS),
        tof_bounds_days=(_TOF1_BOUNDS, _TOF2_BOUNDS),
        synodic_pair=("E", "M"),
        closure_body="E",
        closure_flyby_alt_km=_EARTH_FLYBY_ALT_KM,
        tof_jitter_half_days=(20.0, 60.0),
        n_starts=5,
        seed=0,
    )


def _horizon_tcm(rid: str, seed: cb.BatchRowSeed, start: datetime, ephem: Ephemeris) -> list[float]:
    """The per-lap maintenance ΔV (km/s) over :data:`N_LAPS` laps from ``start``."""
    dvs: list[float] = []
    for lap in range(N_LAPS):
        pdate = start + timedelta(seconds=lap * _SYNODIC_SEC)
        res = _maintenance_at(rid, seed, pdate, ephem)
        dvs.append(res.maintenance_dv_kms)
    return dvs


# ---------------------------------------------------------------------------
# The lift-set parametrisation: V3 rows pass both halves; break points strict-xfail
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class _Expect:
    """Per-row expected V3 outcome (the recorded diagnostic, not a fabricated
    golden): which half binds, for the strict-xfail markers."""

    rid: str
    reaches_v3: bool
    break_half: str | None  # None | "phase-match (a)" | "horizon-TCM (b)"


# Recorded outcomes from the 3-lap DE440 diagnostic (see the results note). These
# drive the strict-xfail markers; the ASSERTIONS below still recompute everything
# from sourced inputs — the table only decides which rows are gated as PASS vs
# strict-xfail break points (so a break point flipping to PASS fails loudly).
_EXPECT: dict[str, _Expect] = {
    # V3 (both halves): Aldrin (corroborates #161) + the 5.75ggF3 free-return row.
    "aldrin-classic-em-k1-outbound": _Expect(
        "aldrin-classic-em-k1-outbound", reaches_v3=True, break_half=None
    ),
    "russell-ch4-5.75ggF3": _Expect("russell-ch4-5.75ggF3", reaches_v3=True, break_half=None),
    # Break at half (b): phase-match closes & V_inf-matches, but the 3-lap horizon
    # TCM exceeds the 3.0 bar (lap 2 lands on a high-e ~0.44 member, dv ~13 km/s).
    "russell-ch4-5.30gGf3": _Expect(
        "russell-ch4-5.30gGf3", reaches_v3=False, break_half="horizon-TCM (b)"
    ),
    "russell-ch4-9.94Gg3": _Expect(
        "russell-ch4-9.94Gg3", reaches_v3=False, break_half="horizon-TCM (b)"
    ),
    # Break at half (a): no converged, V_inf-matched, in-window DE440 basin within
    # the seed sweep (the deep-aphelion high-e row continues off-fingerprint).
    "russell-ch4-9.353Gg2": _Expect(
        "russell-ch4-9.353Gg2", reaches_v3=False, break_half="phase-match (a)"
    ),
}


def _params() -> list[ParameterSet]:
    params: list[ParameterSet] = []
    for rid in cb.LIFT_SET_RIDS:
        exp = _EXPECT[rid]
        marks = []
        if not exp.reaches_v3:
            marks.append(
                pytest.mark.xfail(
                    strict=True,
                    reason=(
                        f"{rid}: continuation V3 break at {exp.break_half} — recorded "
                        f"break point (see docs/notes/2026-06-08-continuation-batch-"
                        f"results.md). NEVER soften: removing the marker requires the "
                        f"row to genuinely clear both V3 halves."
                    ),
                )
            )
        params.append(pytest.param(rid, id=rid, marks=marks))
    return params


@pytest.mark.slow
@pytest.mark.parametrize("rid", _params())
def test_continuation_row_reaches_v3(rid: str, astropy_ephem: Ephemeris) -> None:
    """V3 evidence for one lift-set row: phase-match (a) AND bounded horizon TCM (b).

    Rows recorded as reaching V3 must pass BOTH halves; rows recorded as break
    points are strict-xfail and MUST still fail one half (the day they pass, the
    marker comes off). Sourced anchors EXPECTED; continued / maintenance values
    EVIDENCE; the plausibility bar is the project CONVENTION (not a sourced budget).
    """
    seed = cb.seeds_for_ids((rid,))[0]
    continued = _continued(seed, astropy_ephem)

    # --- Half (a): phase-match (ballistic close + V_inf fingerprint + in-window).
    # The batch driver's reaches_v3_closure is the conjunction; assert its parts so
    # a failure pinpoints which sub-gate broke (recorded evidence on failure).
    half_a = continued.reaches_v3_closure
    assert continued.selected or not half_a, (rid, "no V_inf-matched in-window basin selected")
    if half_a:
        assert continued.ballistic_within_tol, (rid, continued.best_final_residual_kms)
        assert continued.vinf_matched, (rid, continued.vinf_offset_kms)
        assert continued.phase_matched, (rid, continued.window_offset_days)

    # If half (a) already broke, the horizon TCM is undefined — the row fails here
    # (strict-xfail captures it; the residual/offset are in the assert payloads).
    assert half_a, (
        rid,
        "phase-match (a) break",
        {
            "selected": continued.selected,
            "ballistic": continued.ballistic_within_tol,
            "best_final_residual_kms": continued.best_final_residual_kms,
            "vinf_offset_kms": continued.vinf_offset_kms,
            "window_offset_days": continued.window_offset_days,
        },
    )

    # --- Half (b): bounded horizon TCM over N_LAPS, every lap within the bar.
    per_lap_dv = _horizon_tcm(rid, seed, continued.emerged_t0, astropy_ephem)
    # BOUNDED: finite, non-negative, accumulating under n_laps * the per-cycle bar.
    horizon_tcm = sum(d for d in per_lap_dv if d == d)
    assert all(d == d for d in per_lap_dv), (rid, "NaN lap (no basin)", per_lap_dv)
    assert horizon_tcm >= 0.0, (rid, per_lap_dv)
    # WITHIN BAR: every per-cycle ΔV under the engineering plausibility convention.
    assert max(per_lap_dv) < MAINTENANCE_DV_CONVENTION_KMS, (rid, per_lap_dv)
    assert horizon_tcm < N_LAPS * MAINTENANCE_DV_CONVENTION_KMS, (rid, horizon_tcm, per_lap_dv)
