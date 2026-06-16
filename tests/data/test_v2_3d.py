"""Tests for the V2 long-span bounded-drift 3D gauntlet (#306 Phase 1 Part B).

Sourced golden discipline (per ``feedback_golden_tests_sourced_only``)
---------------------------------------------------------------------
V2's input is a CR3BP IC + period; the IC is taken from sourced inputs:

  * The JPL Periodic Orbits API L1 halo fixture
    (``tests/verify/fixtures/jpl_earth_moon_l1_halo_sample.json``) — sourced
    from NASA/JPL, mass-ratio 1.215058560962404e-02; used for the V2 POSITIVE
    test (a 3D orbit whose stability tolerates 3 same-model cycles).
  * The planar Braik-Ross C11a IC (catalogue row ``braik-ross-c11a-cycler-2026``)
    — sourced, but EXPLICITLY catalogued as V1 with the note "NOT V2: an
    unstable orbit cannot satisfy V2-ballistic's bounded-drift-over->=3-laps
    requirement". This is used as a V2 NEGATIVE / hyperbolic-amplification
    test — V2 correctly REJECTS it, mirroring the catalogue's own claim.

Test cases
----------
  1. JPL L1 halo (3D Earth-Moon) — V2 passes; drift across 3 cycles stays
     well below the 50,000 km same-model floor (it's a mildly unstable
     orbit but corrector-clean closure absorbs gracefully).
  2. Planar Braik-Ross C11a — V2 fails; the orbit is hyperbolic and the
     catalogue acknowledges it cannot earn V2-ballistic.
  3. 2x-period non-periodic IC — V2 fails (negative control).
  4. Argument validation: ``n_cycles=2`` is rejected (spec §14 minimum 3).
  5. Audit-trail fields and malformed-system guard.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.data.validation.v1_3d import run_v1_3d
from cyclerfinder.data.validation.v2_3d import (
    V2_DRIFT_FLOOR_KMS,
    V2_N_CYCLES_MIN,
    V2Verdict3D,
    run_v2_3d,
)
from cyclerfinder.search.cr3bp_general_periodic_3d import (
    FREE_VARS_SYMMETRIC_TULIP,
    RESIDUAL_PERPENDICULAR_HALF_PERIOD,
    correct_general_periodic_3d,
)

_REPO_ROOT = Path(__file__).resolve().parent.parent.parent

# Sourced planar golden (Braik-Ross 2026 C11a, in catalogue.yaml).
C11A_X0 = -0.8116406668238195
C11A_YDOT0 = -0.11859055759763637
C11A_PERIOD_TU = 9.69107744379376
EM_MU = 1.2150584270572e-2  # Braik-Ross 2026 Table 1
EM_L_KM = 384400.0
EM_T_S = 375699.8


def _em_system() -> cr3bp.CR3BPSystem:
    """Earth-Moon CR3BP at the sourced Braik-Ross mu."""
    return cr3bp.CR3BPSystem(
        mu=EM_MU,
        primary="earth",
        secondary="moon",
        l_km=EM_L_KM,
        t_s=EM_T_S,
    )


def _jpl_l1_halo_state_and_system() -> tuple[np.ndarray, float, cr3bp.CR3BPSystem]:
    """Read the lowest-stability JPL L1 halo from the test fixture.

    The fixture lives at ``tests/verify/fixtures/jpl_earth_moon_l1_halo_sample.json``.
    Returns ``(state0, period_tu, system)`` at the JPL fixture's mass ratio
    (1.215058560962404e-02 — slightly different from the Braik-Ross value).
    """
    path = _REPO_ROOT / "tests" / "verify" / "fixtures" / "jpl_earth_moon_l1_halo_sample.json"
    if not path.exists():
        pytest.skip(f"JPL L1 halo fixture missing at {path}")
    data = json.loads(path.read_text())
    rows = data["data"]
    # Lowest stability indicator is the closest to "neutral" — the most
    # forgiving member for V2's bounded-drift gate.
    best = min(rows, key=lambda r: float(r[8]))
    state0 = np.array([float(x) for x in best[:6]], dtype=np.float64)
    period_tu = float(best[7])
    system_meta = data["system"]
    system = cr3bp.CR3BPSystem(
        mu=float(system_meta["mass_ratio"]),
        primary="earth",
        secondary="moon",
        l_km=float(system_meta["lunit"]),
        t_s=float(system_meta["tunit"]),
    )
    return state0, period_tu, system


def test_v2_n_cycles_min_is_spec_value() -> None:
    """Fabrication guard: spec §14 says ≥3 continuous laps."""
    assert V2_N_CYCLES_MIN == 3


def test_v2_drift_floor_is_module_value() -> None:
    """The default V2 drift floor is the documented same-model 50,000 km.

    Sourced from ``propagate.DRIFT_TOLERANCE_KM``; this is the V2-ballistic
    same-model bar this project already uses for E-M class cyclers and
    absorbs the natural hyperbolic-instability amplification of a clean IC
    over 3 cycles.
    """
    assert V2_DRIFT_FLOOR_KMS == 50_000.0


def test_v2_jpl_l1_halo_passes() -> None:
    """V2 on a sourced JPL L1 halo (3D Earth-Moon).

    The fixture row is the lowest-stability member of the JPL L1 halo
    family — mildly hyperbolic but corrector-clean enough that 3 same-model
    cycles stay well below the 50,000 km floor. The IC, period and system
    parameters are read from the sourced fixture.

    This is the V2 POSITIVE test: a known-stable 3D periodic orbit must
    pass the bounded-drift gate.
    """
    state0_raw, period_tu, system = _jpl_l1_halo_state_and_system()
    # Re-close under the corrector (V1 best-practice) before V2.
    orbit = correct_general_periodic_3d(
        system,
        state0_raw,
        period_tu,
        tol=1e-12,
        independent_tol=1e-8,
    )
    verdict = run_v2_3d(
        "jpl-l1-halo-em",
        orbit.state0,
        orbit.T_TU,
        system,
        n_cycles=V2_N_CYCLES_MIN,
        notes="JPL Periodic Orbits API L1 halo fixture; mildly unstable 3D EM orbit",
    )
    assert isinstance(verdict, V2Verdict3D)
    assert verdict.converged_at_each_return
    assert verdict.n_cycles_propagated == V2_N_CYCLES_MIN
    assert verdict.passes_v2, (
        f"V2 FAIL on JPL L1 halo: max_drift={verdict.max_drift_kms:.3e} km > "
        f"floor={verdict.drift_floor_kms} km"
    )
    # Headline: drift stays well below the same-model floor.
    assert verdict.max_drift_kms < V2_DRIFT_FLOOR_KMS
    assert len(verdict.per_cycle_drift_kms) == V2_N_CYCLES_MIN


def test_v2_planar_braik_ross_c11a_correctly_rejected() -> None:
    """V2 correctly REJECTS the planar Braik-Ross C11a.

    Catalogue row ``braik-ross-c11a-cycler-2026`` carries the explicit note:
    "NOT V2: an unstable orbit cannot satisfy V2-ballistic's
    bounded-drift-over->=3-laps requirement".

    This test verifies that our V2 implementation enforces exactly that
    rejection — even a corrector-clean IC of an unstable orbit must FAIL
    V2 because hyperbolic amplification over 3 cycles drives same-model
    drift above the 50,000 km floor (we observe ~4e5 km).

    A V2 PASS here would be a false positive: an unstable orbit slipping
    past the gauntlet's ballistic-bounded-drift gate.
    """
    system = _em_system()
    state0 = np.array([C11A_X0, 0.0, 0.0, 0.0, C11A_YDOT0, 0.0], dtype=np.float64)
    # Re-close (V1) before V2 — V2 takes a corrector-clean IC.
    v1 = run_v1_3d(
        "c11a-pre-v2",
        state0,
        C11A_PERIOD_TU,
        system,
        free_vars=FREE_VARS_SYMMETRIC_TULIP,
        residual_indices=RESIDUAL_PERPENDICULAR_HALF_PERIOD,
        is_half_period_residual=True,
    )
    assert v1.passes_v1  # V1 must still pass — the IC is genuinely periodic.
    # Now V2 must reject — the orbit's hyperbolic instability drives
    # cumulative drift over the floor.
    state_clean = np.array([C11A_X0, 0.0, 0.0, 0.0, C11A_YDOT0, 0.0], dtype=np.float64)
    verdict = run_v2_3d(
        "braik-ross-c11a-planar",
        state_clean,
        C11A_PERIOD_TU,
        system,
        n_cycles=V2_N_CYCLES_MIN,
        notes="catalogue row braik-ross-c11a-cycler-2026 is V1 but explicitly NOT V2",
    )
    assert verdict.converged_at_each_return, "propagator must not fail outright"
    assert not verdict.passes_v2, (
        f"V2 false-positive on unstable Braik-Ross C11a: "
        f"max_drift={verdict.max_drift_kms:.3e} km — catalogue says this orbit "
        f"is NOT V2"
    )
    # The drift must be MUCH larger than the floor (hyperbolic amplification).
    assert verdict.max_drift_kms > V2_DRIFT_FLOOR_KMS


def test_v2_rejects_non_periodic_ic() -> None:
    """A non-periodic IC accumulates drift above the floor → V2 fails.

    Start from the planar C11a IC but assert a WRONG period (2x the real
    period). The propagator returns to a different state every cycle, so
    cumulative drift easily exceeds the floor.
    """
    system = _em_system()
    state0 = np.array([C11A_X0, 0.0, 0.0, 0.0, C11A_YDOT0, 0.0], dtype=np.float64)
    verdict = run_v2_3d(
        "wrong-period-non-periodic",
        state0,
        C11A_PERIOD_TU * 2.0,  # wrong: 2x real period
        system,
        n_cycles=V2_N_CYCLES_MIN,
        notes="negative control: 2x period, drift should breach floor",
    )
    if verdict.converged_at_each_return:
        assert verdict.max_drift_kms > V2_DRIFT_FLOOR_KMS, (
            f"V2 false-positive: 2x-period IC drift={verdict.max_drift_kms:.3e} km "
            f"is below floor — a same-model V2 must catch this"
        )
        assert not verdict.passes_v2
    else:
        # Propagator gave up — also a valid (non-pass) outcome.
        assert not verdict.passes_v2


def test_v2_rejects_n_cycles_below_spec_min() -> None:
    """``n_cycles=2`` is below spec §14's 3 — must raise."""
    system = _em_system()
    state0 = np.array([C11A_X0, 0.0, 0.0, 0.0, C11A_YDOT0, 0.0], dtype=np.float64)
    with pytest.raises(ValueError, match="n_cycles >= 3"):
        run_v2_3d(
            "below-min-cycles",
            state0,
            C11A_PERIOD_TU,
            system,
            n_cycles=2,
        )


def test_v2_rejects_malformed_system() -> None:
    """A CR3BPSystem with l_km=0 must fail loudly."""
    bad_system = cr3bp.CR3BPSystem(
        mu=EM_MU,
        primary="earth",
        secondary="moon",
        l_km=0.0,
        t_s=EM_T_S,
    )
    state0 = np.array([C11A_X0, 0.0, 0.0, 0.0, C11A_YDOT0, 0.0], dtype=np.float64)
    with pytest.raises(ValueError, match="invalid CR3BP system"):
        run_v2_3d(
            "bad-system",
            state0,
            C11A_PERIOD_TU,
            bad_system,
        )


def test_v2_verdict_carries_audit_fields() -> None:
    """The V2 verdict must carry the floors it was held against."""
    system = _em_system()
    state0 = np.array([C11A_X0, 0.0, 0.0, 0.0, C11A_YDOT0, 0.0], dtype=np.float64)
    verdict = run_v2_3d(
        "audit-trail-check",
        state0,
        C11A_PERIOD_TU,
        system,
        n_cycles=V2_N_CYCLES_MIN,
        notes="audit-trail check",
    )
    assert verdict.drift_floor_kms == V2_DRIFT_FLOOR_KMS
    assert verdict.n_cycles_min == V2_N_CYCLES_MIN
    assert verdict.candidate_id == "audit-trail-check"
    assert verdict.notes == "audit-trail check"
    assert len(verdict.per_cycle_drift_kms) == verdict.n_cycles_propagated
    assert len(verdict.per_cycle_velocity_drift_kms) == verdict.n_cycles_propagated
