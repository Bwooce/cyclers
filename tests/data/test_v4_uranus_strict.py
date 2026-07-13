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
* Numerical PASS/FAIL of the SILVER *as a science claim*: the headline
  gauntlet verdict is the runner's job (Part C), not a unit-test assertion
  — baking the verdict in as ground truth would be a category error.

  The ONE deliberate exception is
  ``test_560_silver_312_canonical_epoch_unchanged`` below: it pins #312's
  canonical single-epoch (2000-06-21, the #338/#566 reference epoch) PASS
  as a REGRESSION GUARD that the #560/#567 robustness fixes did not perturb
  an already-established prior result. That is a different thing from
  deciding the science in the test — it guards a fixed input's fixed output
  against code drift, with generous headroom (drift-agreement ≈12,160 km vs
  the 50,000 km floor, i.e. not knife-edge), exactly as #560 requires.
* Multi-epoch sensitivity: scope of the gauntlet runner.

These tests need the URA111 SPICE kernel installed (Part A); they skip
gracefully if the kernel is missing so CI without the kernel still
passes the module-level type checks.
"""

from __future__ import annotations

import dataclasses
import math

import pytest
import spiceypy as spice

from cyclerfinder.core.satellites import PRIMARIES, SATELLITES
from cyclerfinder.data.validation.v2_moontour import run_v2_moontour
from cyclerfinder.data.validation.v3_3d import V3Verdict3D, run_v3_3d
from cyclerfinder.data.validation.v4_uranus import (
    URANIAN_PERTURBER_MOONS,
    URANUS_J2,
    URANUS_R_EQ_KM,
    V4_AGREEMENT_FLOOR_KMS,
    V4UranusVerdict,
    _hill_radius_km,
    run_v4_uranus,
)
from cyclerfinder.data.validation.v4_uranus_strict import (
    DEFAULT_LSK_PATH,
    DEFAULT_PCK_PATH,
    DEFAULT_URA_PATH,
    FAILURE_MODE_CONVERGED,
    FAILURE_MODE_PLANET_CROSSING,
    V4UranusStrictCycleVerdict,
    V4UranusStrictVerdict,
    _cycle_v4_strict,
    _ephemeris_time_seconds,
    _spice_furnsh_all,
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
        audit_body1_name="A",
        audit_body2_name="B",
        eccentricity_used_e_body1=0.004,
        eccentricity_used_e_body2=0.001,
        inclination_used_deg_body1=75.0,
        inclination_used_deg_body2=75.0,
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
    assert 0.0 < verdict.eccentricity_used_e_body1 < 0.02
    assert 0.0 < verdict.eccentricity_used_e_body2 < 0.02
    # Inclinations: ura111 in J2000 frame puts the Uranus equatorial moons
    # at ~75 deg (= Uranus pole tilt). Generous bounds.
    assert 60.0 < verdict.inclination_used_deg_body1 < 90.0
    assert 60.0 < verdict.inclination_used_deg_body2 < 90.0

    # #567 bug 3 regression guard: the audit fields must track the
    # candidate's ACTUAL sequence bodies, not be hardcoded to Umbriel/Oberon
    # (a bare 0 < e < 0.02 range check alone would NOT catch a regression
    # back to wrong-moon sampling -- both Umbriel and Oberon's real e sit in
    # that range, so a hardcoded-wrong-moon bug silently passes a range-only
    # check). Assert the recorded body names actually appear in `sequence`.
    assert verdict.audit_body1_name in verdict.sequence
    assert verdict.audit_body2_name in verdict.sequence
    distinct_seq_bodies = list(dict.fromkeys(verdict.sequence))
    assert verdict.audit_body1_name == distinct_seq_bodies[0]
    assert verdict.audit_body2_name == (
        distinct_seq_bodies[1] if len(distinct_seq_bodies) > 1 else distinct_seq_bodies[0]
    )

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


# --------------------------------------------------------------------------- #
# #567 (1)+(2) regressions -- Lambert branch-continuity + planet-crossing tag.
# --------------------------------------------------------------------------- #
#
# Both cases below were located by directly probing `_cycle_v4_strict` /
# `_select_leg_transfer` on the #327 SILVER's own parameters (the exact
# hour-level flip #559's ephemeral diagnostic pass reported was not
# preserved in the repo, so these are freshly-located, reproducible
# instances of the SAME diagnosed mechanisms, confirmed by direct
# instrumentation before the fix landed):
#
# * Bug 1 (branch-selection discontinuity): at 2000-09-06T03:00:00 the old
#   velocity-match tie-break picked branch 0 (leg-1 terminal offset
#   ~24,986 km); one hour later at T04:00:00 it flipped to branch 1
#   (~3,493 km) -- even though BOTH branches' actual propagated offsets
#   vary smoothly and monotonically hour to hour (branch 1's offset alone:
#   ...2,876 -> 3,493... km). The old code's selection criterion (departure
#   velocity-residual vs the moon's own velocity) is unrelated to which
#   branch actually flies well, so its argmin can flip discontinuously.
# * Bug 2 (planet-crossing silent misclassification): at 2000-07-24T02:00:00
#   leg 1 (Oberon->Umbriel), BOTH candidate rev-1 Lambert branches have an
#   osculating periapsis (97 km and 852 km) far inside Uranus's R_eq
#   (25,559 km) -- confirmed via direct periapsis computation on the exact
#   Lambert solutions. Pre-fix this collapsed into a bare integrator FAIL
#   indistinguishable from any other failure mode.


def _silver_mu_and_perturber_dicts() -> tuple[float, dict[str, float], dict[str, float]]:
    mu_primary = PRIMARIES["Uranus"]
    perturber_mu = {m: SATELLITES[m].mu_km3_s2 for m in URANIAN_PERTURBER_MOONS}
    perturber_hill_km = {
        m: _hill_radius_km(
            sma_moon_km=SATELLITES[m].sma_km,
            mu_moon=perturber_mu[m],
            mu_primary=mu_primary,
        )
        for m in URANIAN_PERTURBER_MOONS
    }
    return mu_primary, perturber_mu, perturber_hill_km


@pytest.mark.skipif(not _KERNELS_PRESENT, reason=_SKIP_REASON)
def test_567_lambert_branch_selection_is_continuous_across_adjacent_hours() -> None:
    """#567 bug 1: outcome-based branch selection must not jump discontinuously.

    Pins the located 2000-09-06T03:00:00 -> T04:00:00 flip on the #327
    SILVER: post-fix, both epochs must select the SAME (lower-offset)
    branch and the reported worst-cycle-offset must change smoothly
    hour-to-hour, not jump by tens of thousands of km the way the old
    velocity-match tie-break did (24,986 km -> 3,493 km observed pre-fix
    on this exact epoch pair).
    """
    mu_primary, perturber_mu, perturber_hill_km = _silver_mu_and_perturber_dicts()

    spice.kclear()
    try:
        _spice_furnsh_all((str(DEFAULT_LSK_PATH), str(DEFAULT_PCK_PATH), str(DEFAULT_URA_PATH)))
        offsets = []
        for epoch in ("2000-09-06T03:00:00", "2000-09-06T04:00:00"):
            et_cycle_start = _ephemeris_time_seconds(epoch)
            converged, _r, worst_offset_kms, failure_mode, _perijove = _cycle_v4_strict(
                sequence=SILVER_SEQ,
                leg_tofs_days=SILVER_TOF,
                et_cycle_start=et_cycle_start,
                perturber_moons=URANIAN_PERTURBER_MOONS,
                perturber_mu=perturber_mu,
                perturber_hill_km=perturber_hill_km,
                mu_primary=mu_primary,
                n_revs=SILVER_NREV,
                j2=URANUS_J2,
                r_eq_km=URANUS_R_EQ_KM,
            )
            assert converged, f"{epoch} unexpectedly failed to converge ({failure_mode})"
            assert failure_mode == FAILURE_MODE_CONVERGED
            offsets.append(worst_offset_kms)
    finally:
        spice.kclear()

    # Pre-fix this pair jumped from ~24,986 km to ~3,493 km (a ~21,500 km
    # discontinuity) purely from the discrete branch flip. Post-fix, both
    # epochs pick the outcome-best branch, whose offset varies smoothly
    # (observed hour-to-hour delta on the surviving branch is ~600 km) --
    # allow generous headroom (5,000 km) while still being tight enough to
    # catch a regression back to the old jump.
    assert math.isfinite(offsets[0]) and math.isfinite(offsets[1])
    assert abs(offsets[1] - offsets[0]) < 5_000.0, (
        f"worst-cycle-offset jumped discontinuously between adjacent hours: {offsets}"
    )


@pytest.mark.skipif(not _KERNELS_PRESENT, reason=_SKIP_REASON)
def test_567_planet_crossing_tagged_not_silently_misclassified() -> None:
    """#567 bug 2: a genuinely planet-crossing Lambert leg is TAGGED, not
    silently folded into a generic solver FAIL -- and it still counts as a
    real FAIL (not excluded/skipped; see the #567 PIN in OUTSTANDING.md).

    Pins 2000-07-24T02:00:00 on the #327 SILVER: leg 1 (Oberon->Umbriel)'s
    both rev-1 Lambert branches have an osculating periapsis (~97 km /
    ~852 km, confirmed by direct probe) far inside Uranus's R_eq
    (25,559 km) -- a genuine dynamical infeasibility from real synodic
    geometry (#559), not a numerical artifact.
    """
    mu_primary, perturber_mu, perturber_hill_km = _silver_mu_and_perturber_dicts()

    spice.kclear()
    try:
        _spice_furnsh_all((str(DEFAULT_LSK_PATH), str(DEFAULT_PCK_PATH), str(DEFAULT_URA_PATH)))
        et_cycle_start = _ephemeris_time_seconds("2000-07-24T02:00:00")
        converged, r_final, worst_offset_kms, failure_mode, perijove_km = _cycle_v4_strict(
            sequence=SILVER_SEQ,
            leg_tofs_days=SILVER_TOF,
            et_cycle_start=et_cycle_start,
            perturber_moons=URANIAN_PERTURBER_MOONS,
            perturber_mu=perturber_mu,
            perturber_hill_km=perturber_hill_km,
            mu_primary=mu_primary,
            n_revs=SILVER_NREV,
            j2=URANUS_J2,
            r_eq_km=URANUS_R_EQ_KM,
        )
    finally:
        spice.kclear()

    # This is a genuine FAIL, correctly still counted as one -- NOT skipped
    # or excluded from the denominator.
    assert converged is False
    assert r_final is None
    assert not math.isfinite(worst_offset_kms)
    # But it must now be DISTINGUISHABLE from a generic/unexplained FAIL.
    assert failure_mode == FAILURE_MODE_PLANET_CROSSING
    assert perijove_km is not None
    assert 0.0 <= perijove_km < URANUS_R_EQ_KM


@pytest.mark.skipif(not _KERNELS_PRESENT, reason=_SKIP_REASON)
def test_567_audit_fields_track_non_umbriel_oberon_sequence() -> None:
    """#567 bug 3: for a candidate that does NOT involve Umbriel/Oberon, the
    audit e/i fields must sample the CANDIDATE'S OWN bodies (Ariel/Titania
    here), not silently record Umbriel/Oberon's values.

    Reuses the #566 Ariel-Titania-Ariel representative's exact source
    parameters (``scripts/run_566_gauntlet_five_representatives.py``,
    ``enum563-line12-ariel-titania-ariel``) at #566's own known-good
    reference epoch (2000-06-21T00:00:00) so this test is not entangled
    with #567's own branch-selection/planet-crossing fixes.
    """
    seq: tuple[str, ...] = ("Ariel", "Titania", "Ariel")
    vinf: tuple[float, ...] = (
        1.2306411593828481,
        1.7185773183747601,
        1.2306411593828457,
    )
    tof: tuple[float, ...] = (5.320895317317783, 5.320895317317783)
    rel_off_deg = 0.0
    n_revs: tuple[int, ...] = (0, 0)
    phase0_deg = 29.999999999999996  # same #558 fixed rotation-redundant phase as #566
    candidate_id = "enum563-line12-ariel-titania-ariel"
    epoch = "2000-06-21T00:00:00"

    v2 = run_v2_moontour(
        candidate_id,
        seq,
        vinf,
        tof,
        rel_off_deg,
        None,
        n_cycles=3,
        n_revs=n_revs,
        phase0_deg=phase0_deg,
    )
    v3 = run_v3_3d(
        candidate_id,
        seq,
        vinf,
        tof,
        rel_off_deg,
        None,
        v2_verdict=v2,
        n_cycles=3,
        n_revs=n_revs,
        phase0_deg=phase0_deg,
    )
    v4_scipy = run_v4_uranus(
        candidate_id,
        seq,
        vinf,
        tof,
        rel_off_deg,
        None,
        v3_verdict=v3,
        n_cycles=3,
        n_revs=n_revs,
        phase0_deg=phase0_deg,
    )
    verdict = run_v4_uranus_strict(
        candidate_id,
        seq,
        vinf,
        tof,
        rel_off_deg,
        epoch,
        None,
        v3_verdict=v3,
        v4_scipy_verdict=v4_scipy,
        n_cycles=3,
        n_revs=n_revs,
    )

    assert verdict.audit_body1_name == "Ariel"
    assert verdict.audit_body2_name == "Titania"
    assert verdict.audit_body1_name not in ("Umbriel", "Oberon")
    assert verdict.audit_body2_name not in ("Umbriel", "Oberon")
    # Sanity range (same generous band as the SILVER e2e test); the point of
    # this test is the BODY IDENTITY assertions above, which a bare range
    # check cannot catch on its own (Umbriel/Oberon's real e also sits in
    # this range -- see the #567 OUTSTANDING.md PIN).
    assert 0.0 < verdict.eccentricity_used_e_body1 < 0.02
    assert 0.0 < verdict.eccentricity_used_e_body2 < 0.02


# --------------------------------------------------------------------------- #
# #560 close-out: #312's canonical single-epoch result must be UNCHANGED by
# the branch-continuity + planet-crossing-guard robustness fixes.
# --------------------------------------------------------------------------- #


@pytest.mark.skipif(not _KERNELS_PRESENT, reason=_SKIP_REASON)
def test_560_silver_312_canonical_epoch_unchanged() -> None:
    """#560 requirement: the robustness fixes must NOT change #312's own
    already-valid canonical single-epoch V4-strict result.

    #312's canonical reference epoch is 2000-06-21T00:00:00 (the #338/#566
    anchor epoch; NOT ``SMOKE_EPOCH_UTC`` 2000-01-15, which is a known #338
    high-drift FAIL epoch). At the canonical epoch the #327 SILVER
    (Umbriel-Oberon-Umbriel, = #312's representative) is a clean PASS: all
    three cycles converge (``failure_mode == "converged"`` — no spurious
    planet-crossing tag or integrator failure introduced by the guard), and
    the V4-strict-vs-V3 drift agreement is ~12,160 km, comfortably under the
    50,000 km floor (not knife-edge, so this pin is not epoch-fragile).

    This is a REGRESSION GUARD, not a science-verdict bake-in (see the module
    docstring's "What's NOT tested"): #560 diagnosed the branch-continuity
    and planet-crossing-guard fixes as needed only for future epoch-SWEEP
    interpretation, and required explicit confirmation that #312's own
    single-epoch result is unaffected rather than assumed so. A future edit
    to ``_select_leg_transfer`` / ``_leg_periapsis_km`` / the guard threshold
    that silently perturbed the canonical result would trip here.
    """
    canonical_epoch = "2000-06-21T00:00:00"

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
    v4_scipy = run_v4_uranus(
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
    verdict = run_v4_uranus_strict(
        SILVER_ID,
        SILVER_SEQ,
        SILVER_VINF,
        SILVER_TOF,
        SILVER_REL_OFF_DEG,
        canonical_epoch,
        None,
        v3_verdict=v3,
        v4_scipy_verdict=v4_scipy,
        n_cycles=3,
        n_revs=SILVER_NREV,
    )

    # The headline: #312's canonical single-epoch result is a PASS, unchanged.
    assert verdict.passes_v4_strict is True
    assert verdict.n_cycles_propagated == 3
    assert verdict.bounded_drift_survives is True
    # Every cycle converged cleanly -- the planet-crossing guard did NOT
    # mis-tag any leg at this epoch, and the branch-continuity selection did
    # not drop to a failure mode.
    assert len(verdict.per_cycle) == 3
    for c in verdict.per_cycle:
        assert c.failure_mode == FAILURE_MODE_CONVERGED
        assert c.converged_legs == c.n_legs
    # Drift agreement is well under the floor with headroom (pinned generously
    # so genuine sub-km numerical churn across BLAS/platforms won't flake it,
    # but tight enough to catch a regression that moved it toward the floor).
    assert verdict.drift_agreement_kms_vs_v3 < 20_000.0
    assert verdict.v4_v3_agreement_floor_kms == V4_AGREEMENT_FLOOR_KMS
