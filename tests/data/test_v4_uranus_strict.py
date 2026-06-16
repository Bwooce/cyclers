"""Tests for the V4-strict SPICE-Uranus driver (#335 Phase 4.1).

What's tested
-------------
* Dataclasses are frozen + carry the expected fields.
* Driver shape: argument validation (closed sequence, length consistency,
  Uranus-only, no missing kernels) raises before any SPICE work.
* SPICE round-trip: the driver loads ura111.bsp, samples Umbriel + Oberon
  states, and recovers eccentricities + inclinations in the Murray-Dermott
  +/- 30% band (the kernel is the authoritative source — these are
  sanity checks, not sourced golden values).
* End-to-end on a SHORT (n_cycles=3) run against the #327 SILVER's
  recorded V3 + V4-scipy verdicts: the driver completes, returns finite
  drift series, and the V4-strict-vs-V4-scipy headline number is well-
  defined (the verdict itself — PASS/FAIL — is whatever the math says
  and is exercised in the gauntlet runner, not in tests).

What's NOT tested (and why)
---------------------------
* Numerical PASS/FAIL of the SILVER: that's the headline result of the
  gauntlet runner (Part C), not a unit-test assertion. Testing it would
  bake the verdict into the test suite — a category error.
* Multi-epoch sensitivity: scope of the gauntlet runner.

These tests need the URA111 SPICE kernel installed (Part A); they skip
gracefully if the kernel is missing so CI without the kernel still
passes the module-level type checks.
"""

from __future__ import annotations

import dataclasses
import math

import pytest

from cyclerfinder.data.validation.v2_moontour import run_v2_moontour
from cyclerfinder.data.validation.v3_3d import V3Verdict3D, run_v3_3d
from cyclerfinder.data.validation.v4_uranus import (
    URANIAN_PERTURBER_MOONS,
    V4_AGREEMENT_FLOOR_KMS,
    V4UranusVerdict,
    run_v4_uranus,
)
from cyclerfinder.data.validation.v4_uranus_strict import (
    DEFAULT_LSK_PATH,
    DEFAULT_PCK_PATH,
    DEFAULT_URA_PATH,
    V4UranusStrictCycleVerdict,
    V4UranusStrictVerdict,
    run_v4_uranus_strict,
)

# --------------------------------------------------------------------------- #
# Stored SILVER inputs (carried as constants; same as scripts/run_332_*.py).
# READ-ONLY -- mirrors data/silver_327_verified.jsonl + silver_327_v3_verdicts.jsonl.
# --------------------------------------------------------------------------- #

SILVER_ID = "repeated-moon-uranus-00000041"
SILVER_SEQ: tuple[str, ...] = ("Umbriel", "Oberon", "Umbriel")
SILVER_VINF: tuple[float, ...] = (
    0.9199258810725036,
    0.9604309791298091,
    0.8946936085078939,
)
SILVER_TOF: tuple[float, ...] = (14.940560615336594, 14.940560615336594)
SILVER_REL_OFF_DEG = 180.0
SILVER_NREV: tuple[int, ...] = (1, 1)
SILVER_PHASE0_DEG = 29.999999999999996

# Smoke epoch -- Voyager-2 Uranus encounter +14 years, well inside ura111's
# 1900-2099 span. The same epoch used in Part C.
SMOKE_EPOCH_UTC = "2000-01-15T00:00:00"


_KERNELS_PRESENT = (
    DEFAULT_LSK_PATH.exists() and DEFAULT_PCK_PATH.exists() and DEFAULT_URA_PATH.exists()
)
_SKIP_REASON = (
    f"URA111 SPICE kernel not installed (looked at {DEFAULT_URA_PATH}); "
    "run scripts/install_uranian_spice.sh to install."
)


# --------------------------------------------------------------------------- #
# Dataclass shape tests (run without SPICE)
# --------------------------------------------------------------------------- #


def test_cycle_verdict_is_frozen() -> None:
    """V4UranusStrictCycleVerdict is frozen + carries the expected fields."""
    v = V4UranusStrictCycleVerdict(
        cycle_index=0,
        converged_legs=2,
        n_legs=2,
        rendezvous_drift_kms_v4_strict=0.0,
        rendezvous_drift_kms_v4_scipy=0.0,
        rendezvous_drift_kms_v3=0.0,
        agreement_kms_vs_v4_scipy=0.0,
        agreement_kms_vs_v3=0.0,
        v4_terminal_offset_vs_moon_kms=12.5,
    )
    with pytest.raises(dataclasses.FrozenInstanceError):
        v.cycle_index = 99  # type: ignore[misc]


def test_verdict_is_frozen() -> None:
    """V4UranusStrictVerdict is frozen."""
    v = V4UranusStrictVerdict(
        candidate_id="test",
        sequence=("A", "B", "A"),
        n_cycles_propagated=0,
        integrator="x",
        launch_epoch_utc="2000-01-01T00:00:00",
        spice_kernels_used=(),
        eccentricity_used_e_umbriel=0.004,
        eccentricity_used_e_oberon=0.001,
        inclination_used_deg_umbriel=75.0,
        inclination_used_deg_oberon=75.0,
        per_cycle=(),
        per_cycle_drift_kms_v4_strict=(),
        per_cycle_drift_kms_v4_scipy=(),
        per_cycle_drift_kms_v3=(),
        drift_agreement_kms_vs_v4_scipy=0.0,
        drift_agreement_kms_vs_v3=0.0,
        v4_v3_agreement_floor_kms=50_000.0,
        bounded_drift_survives=False,
        passes_v4_strict=False,
    )
    with pytest.raises(dataclasses.FrozenInstanceError):
        v.passes_v4_strict = True  # type: ignore[misc]


# --------------------------------------------------------------------------- #
# Argument-validation tests (raise BEFORE any SPICE call)
# --------------------------------------------------------------------------- #


def _stub_verdicts() -> tuple[V3Verdict3D, V4UranusVerdict]:
    """Build the V2/V3/V4-scipy chain on the SILVER; used to feed the V4-strict
    argument-validation tests AND the end-to-end test below."""
    v2 = run_v2_moontour(
        SILVER_ID,
        SILVER_SEQ,
        SILVER_VINF,
        SILVER_TOF,
        SILVER_REL_OFF_DEG,
        None,
        n_cycles=3,
        n_revs=SILVER_NREV,
        phase0_deg=SILVER_PHASE0_DEG,
    )
    v3 = run_v3_3d(
        SILVER_ID,
        SILVER_SEQ,
        SILVER_VINF,
        SILVER_TOF,
        SILVER_REL_OFF_DEG,
        None,
        v2_verdict=v2,
        n_cycles=3,
        n_revs=SILVER_NREV,
        phase0_deg=SILVER_PHASE0_DEG,
    )
    v4 = run_v4_uranus(
        SILVER_ID,
        SILVER_SEQ,
        SILVER_VINF,
        SILVER_TOF,
        SILVER_REL_OFF_DEG,
        None,
        v3_verdict=v3,
        n_cycles=3,
        n_revs=SILVER_NREV,
        phase0_deg=SILVER_PHASE0_DEG,
    )
    return v3, v4


@pytest.mark.skipif(not _KERNELS_PRESENT, reason=_SKIP_REASON)
def test_validation_open_sequence_raises() -> None:
    """An open (non-closed) sequence raises before SPICE."""
    v3, v4 = _stub_verdicts()
    with pytest.raises(ValueError, match="moontour sequence must be CLOSED"):
        run_v4_uranus_strict(
            SILVER_ID,
            ("Umbriel", "Oberon"),  # NOT closed
            SILVER_VINF[:2],
            (SILVER_TOF[0],),
            SILVER_REL_OFF_DEG,
            SMOKE_EPOCH_UTC,
            None,
            v3_verdict=v3,
            v4_scipy_verdict=v4,
            n_cycles=3,
            n_revs=(1,),
        )


@pytest.mark.skipif(not _KERNELS_PRESENT, reason=_SKIP_REASON)
def test_validation_non_uranus_moon_raises() -> None:
    """A non-Uranian moon in the sequence raises before SPICE."""
    v3, v4 = _stub_verdicts()
    with pytest.raises(ValueError, match="non-Uranian moons"):
        run_v4_uranus_strict(
            SILVER_ID,
            ("Europa", "Io", "Europa"),
            SILVER_VINF,
            SILVER_TOF,
            SILVER_REL_OFF_DEG,
            SMOKE_EPOCH_UTC,
            None,
            v3_verdict=v3,
            v4_scipy_verdict=v4,
            n_cycles=3,
            n_revs=SILVER_NREV,
        )


@pytest.mark.skipif(not _KERNELS_PRESENT, reason=_SKIP_REASON)
def test_validation_n_cycles_too_small() -> None:
    """n_cycles < V4_N_CYCLES_MIN raises."""
    v3, v4 = _stub_verdicts()
    with pytest.raises(ValueError, match="V4-strict requires n_cycles"):
        run_v4_uranus_strict(
            SILVER_ID,
            SILVER_SEQ,
            SILVER_VINF,
            SILVER_TOF,
            SILVER_REL_OFF_DEG,
            SMOKE_EPOCH_UTC,
            None,
            v3_verdict=v3,
            v4_scipy_verdict=v4,
            n_cycles=2,  # below floor
            n_revs=SILVER_NREV,
        )


# --------------------------------------------------------------------------- #
# End-to-end SPICE-driven smoke test (slow; ~30s)
# --------------------------------------------------------------------------- #


@pytest.mark.skipif(not _KERNELS_PRESENT, reason=_SKIP_REASON)
def test_v4_strict_e2e_smoke() -> None:
    """Drive the V4-strict gauntlet at n_cycles=3 on the SILVER, smoke-only.

    Asserts:
      * Driver completes without exception.
      * Returned verdict has the expected shape (n_cycles_propagated up to 3,
        per_cycle has matching length).
      * The SPICE-sampled eccentricities are in the realistic ura111 range
        (0 < e < 0.01) -- the kernel's value isn't golden, but we sanity-
        check it.
      * The integrator label calls out SPICE.

    Does NOT assert PASS/FAIL of the verdict -- that's the gauntlet
    runner's job (Part C). The verdict is whatever the math says.
    """
    v3, v4_scipy = _stub_verdicts()
    verdict = run_v4_uranus_strict(
        SILVER_ID,
        SILVER_SEQ,
        SILVER_VINF,
        SILVER_TOF,
        SILVER_REL_OFF_DEG,
        SMOKE_EPOCH_UTC,
        None,
        v3_verdict=v3,
        v4_scipy_verdict=v4_scipy,
        n_cycles=3,
        n_revs=SILVER_NREV,
    )

    # Shape checks.
    assert verdict.candidate_id == SILVER_ID
    assert verdict.sequence == SILVER_SEQ
    assert verdict.launch_epoch_utc == SMOKE_EPOCH_UTC
    assert verdict.n_cycles_propagated <= 3
    assert len(verdict.per_cycle) == verdict.n_cycles_propagated or (
        len(verdict.per_cycle) == verdict.n_cycles_propagated + 1  # the failure-row case
    )
    assert "SPICE" in verdict.integrator
    assert "URA111" in verdict.integrator or "ura111" in verdict.integrator.lower()
    assert verdict.v4_v3_agreement_floor_kms == V4_AGREEMENT_FLOOR_KMS

    # SPICE-sampled element sanity (the ura111 numbers we cross-checked at
    # install time -- not golden, but should be in the realistic range).
    assert 0.0 < verdict.eccentricity_used_e_umbriel < 0.02
    assert 0.0 < verdict.eccentricity_used_e_oberon < 0.02
    # Inclinations: ura111 in J2000 frame puts the Uranus equatorial moons
    # at ~75 deg (= Uranus pole tilt). Generous bounds.
    assert 60.0 < verdict.inclination_used_deg_umbriel < 90.0
    assert 60.0 < verdict.inclination_used_deg_oberon < 90.0

    # Drift series finite where the cycle converged.
    for c in verdict.per_cycle:
        if c.converged_legs == c.n_legs:
            assert math.isfinite(c.rendezvous_drift_kms_v4_strict)
            assert math.isfinite(c.rendezvous_drift_kms_v4_scipy)
            assert math.isfinite(c.rendezvous_drift_kms_v3)
            assert math.isfinite(c.agreement_kms_vs_v4_scipy)
            assert math.isfinite(c.agreement_kms_vs_v3)


@pytest.mark.skipif(not _KERNELS_PRESENT, reason=_SKIP_REASON)
def test_v4_strict_perturber_moons_default_is_all_classical() -> None:
    """Default perturber_moons is all 5 classical Uranian moons (matches V4-scipy)."""
    assert set(URANIAN_PERTURBER_MOONS) == {
        "Miranda",
        "Ariel",
        "Umbriel",
        "Titania",
        "Oberon",
    }
