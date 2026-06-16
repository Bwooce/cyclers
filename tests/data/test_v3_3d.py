"""Tests for the V3 6D n-body independent-integrator gauntlet (#331 / #306 Phase 3).

Sourced golden discipline (per ``feedback_golden_tests_sourced_only``)
---------------------------------------------------------------------
The V3 verdict asserts INTEGRATOR AGREEMENT, never a specific drift
number. The SILVER's V_inf tuple + ToFs come from #327's JSONL; the V2
verdict comes from #330. Both are OUR computation — V3 asserts that
REBOUND IAS15 reproduces the V2 driver's per-cycle terminal positions
within :data:`V3_AGREEMENT_FLOOR_KMS`.

The CR3BP regression test members come from ``data/spike_287.jsonl``
and ``data/family_296_3d_subfamilies_299.jsonl`` — their stored
``independent_closure_residual`` / ``independent_closure_L2`` numbers
are an independent-Radau closure produced by the corrector run, used
here as a regression sentinel to lock in that the independent integrator
still reproduces it.

Test cases
----------
  1. #327 SILVER V3 (the load-bearing test): the SILVER's stored
     description (#327 + #330) under V3 IAS15 produces a verdict that
     records the actual integrator-vs-integrator agreement. Asserts
     the verdict is structurally valid; the ``passes_v3`` boolean is
     captured honestly (per ``feedback_orbit_closure_discipline`` —
     don't tune to pass).
  2. #287 spike 3D Braik-Ross (1,1) regression: a 3D member from the
     spike's JSONL re-propagated under an independent integrator must
     reproduce the stored ``independent_closure_L2`` to within a tight
     floor.
  3. #301 k=4 doubly-hyperbolic-pair regression: a member from the
     subfamily JSONL — same pattern.
  4. Negative: a deliberately-broken IC (random perturbation) should
     NOT pass V3 (closure residual exceeds the nondim floor).
  5. Argument validation: ``n_cycles=2`` rejected; non-closed sequence
     rejected; mismatched ToF length rejected; V2 verdict with too few
     cycles rejected.
  6. Audit-trail fields preserved.
  7. Spec constants: :data:`V3_AGREEMENT_FLOOR_KMS` + :data:`V3_N_CYCLES_MIN`
     are fabrication-guarded.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pytest

import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.data.validation.v2_moontour import run_v2_moontour
from cyclerfinder.data.validation.v3_3d import (
    V3_AGREEMENT_FLOOR_KMS,
    V3_N_CYCLES_MIN,
    V3CycleVerdict3D,
    V3PeriodicRegressionVerdict,
    V3Verdict3D,
    run_v3_3d,
    run_v3_periodic_regression,
)

# Repo root.
_REPO_ROOT = Path(__file__).resolve().parent.parent.parent

# #327 SILVER stored fields (data/silver_327_verified.jsonl rows 1 + 4).
SILVER_ID = "repeated-moon-uranus-00000041"
SILVER_SEQ = ("Umbriel", "Oberon", "Umbriel")
SILVER_VINF = (0.9199258810725036, 0.9604309791298091, 0.8946936085078939)
SILVER_TOF = (14.940560615336594, 14.940560615336594)
SILVER_REL_OFF = 180.0
SILVER_NREV = (1, 1)
SILVER_PHASE0 = 29.999999999999996

# Earth-Moon CR3BP for the #287 / #301 regressions.
EM_MU = 1.2150584270572e-2
EM_L_KM = 384400.0
EM_T_S = 375699.8


def _em_system() -> cr3bp.CR3BPSystem:
    """Earth-Moon CR3BP at the catalogued (sourced) Braik-Ross mu."""
    return cr3bp.CR3BPSystem(
        mu=EM_MU,
        primary="earth",
        secondary="moon",
        l_km=EM_L_KM,
        t_s=EM_T_S,
    )


# ---------------------------------------------------------------------------
# Test 1: SILVER V3 verdict (the load-bearing test).
# ---------------------------------------------------------------------------


def test_silver_v3_runs_end_to_end_and_produces_verdict() -> None:
    """SILVER row -> V3 driver produces a verdict; numbers are physically sensible.

    The load-bearing test of #331. V2 (DOP853+Lambert) is run first to
    produce the comparison series, then V3 (REBOUND IAS15) re-propagates
    each cycle. The verdict captures whether IAS15 agrees with the V2
    chain on the per-cycle terminal positions.

    PASS interpretation: V2's bounded-drift quasi_cycler signature is
    real (not integrator artefact). FAIL interpretation: it was numeric.
    """
    v2_verdict = run_v2_moontour(
        SILVER_ID,
        SILVER_SEQ,
        SILVER_VINF,
        SILVER_TOF,
        SILVER_REL_OFF,
        None,
        n_cycles=3,
        n_revs=SILVER_NREV,
        phase0_deg=SILVER_PHASE0,
    )
    v3_verdict = run_v3_3d(
        SILVER_ID,
        SILVER_SEQ,
        SILVER_VINF,
        SILVER_TOF,
        SILVER_REL_OFF,
        None,
        v2_verdict=v2_verdict,
        n_cycles=3,
        n_revs=SILVER_NREV,
        phase0_deg=SILVER_PHASE0,
    )
    assert isinstance(v3_verdict, V3Verdict3D)
    assert v3_verdict.candidate_id == SILVER_ID
    assert v3_verdict.sequence == SILVER_SEQ
    assert v3_verdict.n_cycles_propagated == 3
    # Integrator label set; either IAS15 (env has rebound) or LSODA fallback.
    assert v3_verdict.integrator in {"REBOUND IAS15", "scipy LSODA fallback"}
    # Per-cycle series populated and length-matched.
    assert len(v3_verdict.per_cycle) == 3
    assert len(v3_verdict.per_cycle_drift_kms_v3) == 3
    assert len(v3_verdict.per_cycle_drift_kms_v2) == 3
    # All cycle entries are V3CycleVerdict3D instances.
    for c in v3_verdict.per_cycle:
        assert isinstance(c, V3CycleVerdict3D)
        assert c.converged_legs == c.n_legs == 2
        assert math.isfinite(c.rendezvous_drift_kms_v3)
        assert math.isfinite(c.rendezvous_drift_kms_v2)
        assert math.isfinite(c.agreement_kms)
    # Cycle 0 drift is zero by construction (both V2 and V3).
    cycle0 = v3_verdict.per_cycle[0]
    assert cycle0.rendezvous_drift_kms_v3 == 0.0
    assert cycle0.rendezvous_drift_kms_v2 == 0.0
    assert cycle0.agreement_kms == 0.0
    # The drift_agreement_kms is finite and recorded honestly.
    assert math.isfinite(v3_verdict.drift_agreement_kms)


def test_silver_v3_ias15_passes_agreement_floor() -> None:
    """SILVER V3: REBOUND IAS15 reproduces V2 per-cycle drift to < 100 km.

    With ``ias15_epsilon=1e-12`` the per-cycle agreement should be sub-km;
    the 100 km agreement floor is the wide bar. This is the test that
    surfaces the #331 verdict: V3 PASS means the V2 bounded-drift
    signature is a REAL property of the shared model, not integrator
    artefact. The test asserts the PASS at the spec floor.

    If REBOUND is unavailable on the runner, this test still expects
    agreement at the wider LSODA floor; both architectures are
    independent of V2's DOP853+Lambert chain.
    """
    v2_verdict = run_v2_moontour(
        SILVER_ID,
        SILVER_SEQ,
        SILVER_VINF,
        SILVER_TOF,
        SILVER_REL_OFF,
        None,
        n_cycles=3,
        n_revs=SILVER_NREV,
        phase0_deg=SILVER_PHASE0,
    )
    v3_verdict = run_v3_3d(
        SILVER_ID,
        SILVER_SEQ,
        SILVER_VINF,
        SILVER_TOF,
        SILVER_REL_OFF,
        None,
        v2_verdict=v2_verdict,
        n_cycles=3,
        n_revs=SILVER_NREV,
        phase0_deg=SILVER_PHASE0,
    )
    # Agreement is well under the 100 km floor (IAS15 at 1e-12 epsilon
    # is integrator-noise level vs the analytic Lambert + closed-form
    # Kepler chain).
    assert v3_verdict.drift_agreement_kms < V3_AGREEMENT_FLOOR_KMS, (
        f"V3 disagreement = {v3_verdict.drift_agreement_kms:.3e} km "
        f"exceeds floor {V3_AGREEMENT_FLOOR_KMS} km — V2 signature integrator-noise"
    )
    assert v3_verdict.passes_v3, (
        f"V3 expected PASS at agreement={v3_verdict.drift_agreement_kms:.3e} "
        f"< floor={V3_AGREEMENT_FLOOR_KMS}"
    )
    # IAS15-vs-analytic-Kepler offset at each leg endpoint is tiny — the
    # per-leg numeric agreement of the two Kepler integrators.
    for c in v3_verdict.per_cycle:
        assert c.ias15_vs_analytic_kepler_kms < 1.0, (
            f"cycle {c.cycle_index}: IAS15 vs analytic Kepler "
            f"= {c.ias15_vs_analytic_kepler_kms:.3e} km (> 1 km bar)"
        )


# ---------------------------------------------------------------------------
# Test 2: #287 spike 3D Braik-Ross (1,1) periodic-orbit regression.
# ---------------------------------------------------------------------------


def _read_spike_287_3d_seed() -> tuple[np.ndarray, float, float]:
    """Read a genuinely-3D member from spike_287.jsonl + its stored closure."""
    path = _REPO_ROOT / "data" / "spike_287.jsonl"
    if not path.exists():
        pytest.skip(f"spike_287.jsonl missing at {path}")
    for raw in path.read_text().splitlines():
        if not raw.strip():
            continue
        row = json.loads(raw)
        if not row.get("converged"):
            continue
        z0 = row.get("z0", 0.0)
        if abs(float(z0)) < 1e-3:
            continue
        x0 = float(row["x0"])
        ydot0 = float(row["ydot0"])
        period_tu = float(row["T_TU"])
        closure = float(row.get("independent_closure_L2", row.get("corrector_residual", 1e-9)))
        state0 = np.array([x0, 0.0, float(z0), 0.0, ydot0, 0.0], dtype=np.float64)
        return state0, period_tu, closure
    pytest.skip("no genuinely-3D converged row found in spike_287.jsonl")


def test_v3_spike_287_periodic_regression_agrees_with_stored_closure() -> None:
    """V3 periodic-orbit regression on a #287 3D Braik-Ross member.

    The spike's stored ``independent_closure_L2`` came from an
    independent Radau closure in the corrector loop. V3 re-runs the
    independent propagation at tight tolerances; the resulting closure
    must reproduce the stored value to within numerical noise (the
    integrator is deterministic; this regression locks in that the
    stored 1e-8-ish closure stays 1e-8 across refactors).
    """
    state0, period_tu, stored_closure = _read_spike_287_3d_seed()
    system = _em_system()
    verdict = run_v3_periodic_regression(
        "spike-287-3d-member",
        state0,
        period_tu,
        system,
        stored_closure,
        closure_floor_nondim=1.0e-7,
    )
    assert isinstance(verdict, V3PeriodicRegressionVerdict)
    assert verdict.passes_v3, (
        f"V3 spike regression FAIL: closure_v3={verdict.closure_residual_nondim_v3:.3e} "
        f"> floor={verdict.closure_floor_nondim:.0e}"
    )
    # Tight regression: the V3 closure agrees with the spike's stored value
    # to better than 1e-7 nondim (both are independent re-propagations at
    # 1e-12-ish tolerance; differences should be integrator-noise level).
    assert verdict.agreement_nondim < 1.0e-7, (
        f"V3 closure {verdict.closure_residual_nondim_v3:.3e} vs "
        f"stored {verdict.closure_residual_nondim_stored:.3e}: "
        f"disagreement {verdict.agreement_nondim:.3e}"
    )


# ---------------------------------------------------------------------------
# Test 3: #301 k=4 doubly-hyperbolic-pair sub-family regression.
# ---------------------------------------------------------------------------


def _read_301_k4_doubly_hyperbolic_member() -> tuple[np.ndarray, float, float]:
    """Read a k=4 doubly-hyperbolic member from family_296 subfamilies."""
    path = _REPO_ROOT / "data" / "family_296_3d_subfamilies_299.jsonl"
    if not path.exists():
        pytest.skip(f"family_296 subfamilies JSONL missing at {path}")
    for raw in path.read_text().splitlines():
        if not raw.strip():
            continue
        row = json.loads(raw)
        bracket = row.get("bracket") or {}
        if bracket.get("k") != 4:
            continue
        members = row.get("subfamily_members") or []
        for m in members:
            tag = (m.get("stability_tag") or "").lower()
            # The subfamily here uses hyperbolic-pair tags; the doubly-hyperbolic
            # k=4 evidence is the bracket label (k=4 first-order resonance).
            if "hyper" not in tag:
                continue
            state_nd = m.get("state_nd")
            if not state_nd or len(state_nd) != 6:
                continue
            period_tu = float(m.get("T_TU", 0.0))
            if period_tu <= 0.0:
                continue
            closure = float(m.get("independent_closure_residual", 1e-9))
            return np.array(state_nd, dtype=np.float64), period_tu, closure
    pytest.skip("no k=4 doubly-hyperbolic member found in family_296 subfamilies JSONL")


def test_v3_301_k4_doubly_hyperbolic_regression() -> None:
    """V3 periodic-orbit regression on a #301 k=4 hyperbolic-pair member.

    The k=4 brackets in ``data/family_296_3d_subfamilies_299.jsonl``
    flag the doubly-hyperbolic sub-families. A member of one of those
    sub-families re-propagated under V3 must reproduce its stored
    independent closure.

    Hyperbolic-pair members have very large Floquet abs (~250 for the
    test row) — a strict floor on the absolute closure number would
    be over-tuned. The relevant check is INTEGRATOR AGREEMENT against
    the stored value, which is integrator-noise even for hyperbolic
    members at 1e-12 tolerance over a single period.
    """
    state0, period_tu, stored_closure = _read_301_k4_doubly_hyperbolic_member()
    system = _em_system()
    verdict = run_v3_periodic_regression(
        "family-296-301-k4-doubly-hyperbolic",
        state0,
        period_tu,
        system,
        stored_closure,
        closure_floor_nondim=1.0e-6,
    )
    # The V3 closure agrees with the stored independent_closure_residual
    # to integrator-noise level.
    assert verdict.agreement_nondim < 1.0e-6, (
        f"V3 closure {verdict.closure_residual_nondim_v3:.3e} vs "
        f"stored {verdict.closure_residual_nondim_stored:.3e}: "
        f"disagreement {verdict.agreement_nondim:.3e}"
    )


# ---------------------------------------------------------------------------
# Test 4: Negative — a broken IC must NOT pass V3.
# ---------------------------------------------------------------------------


def test_v3_periodic_negative_broken_ic_does_not_pass() -> None:
    """A deliberately broken IC fails V3: closure exceeds the nondim floor.

    Per ``feedback_orbit_closure_discipline``, a corrupted IC must NOT
    produce a PASS. We take the SILVER spike's IC and add a 1e-2 nondim
    perturbation; the closure residual blows up well past the floor.
    """
    state0, period_tu, _stored = _read_spike_287_3d_seed()
    # Deliberate perturbation: shift x0 by 1e-2 nondim (~4000 km Earth-Moon).
    state0_bad = state0.copy()
    state0_bad[0] += 1.0e-2
    system = _em_system()
    verdict = run_v3_periodic_regression(
        "spike-287-broken-ic",
        state0_bad,
        period_tu,
        system,
        closure_residual_nondim_stored=1.0e-9,  # what an unperturbed close would give
        closure_floor_nondim=1.0e-7,
    )
    assert not verdict.passes_v3, (
        f"broken IC passed V3 at closure {verdict.closure_residual_nondim_v3:.3e}"
    )


def test_v3_silver_lambert_failure_capped() -> None:
    """If the Lambert leg fails (silly ToF), V3 reports ``n_cycles_propagated < n_cycles``."""
    # First produce a valid V2 verdict at n_cycles=3 — needed as input.
    v2_ok = run_v2_moontour(
        SILVER_ID,
        SILVER_SEQ,
        SILVER_VINF,
        SILVER_TOF,
        SILVER_REL_OFF,
        None,
        n_cycles=3,
        n_revs=SILVER_NREV,
        phase0_deg=SILVER_PHASE0,
    )
    # Now run V3 with a tiny ToF that can't host an n_rev=1 Lambert.
    v3_bad = run_v3_3d(
        "lambert-fail",
        SILVER_SEQ,
        SILVER_VINF,
        (0.1, 0.1),
        SILVER_REL_OFF,
        None,
        v2_verdict=v2_ok,
        n_cycles=3,
        n_revs=SILVER_NREV,
        phase0_deg=SILVER_PHASE0,
    )
    assert v3_bad.n_cycles_propagated < 3
    assert not v3_bad.passes_v3


# ---------------------------------------------------------------------------
# Test 5: Argument validation.
# ---------------------------------------------------------------------------


def _silver_v2_for_validation() -> object:
    """Build a V2 verdict for the argument-validation tests."""
    return run_v2_moontour(
        SILVER_ID,
        SILVER_SEQ,
        SILVER_VINF,
        SILVER_TOF,
        SILVER_REL_OFF,
        None,
        n_cycles=3,
        n_revs=SILVER_NREV,
        phase0_deg=SILVER_PHASE0,
    )


def test_v3_rejects_n_cycles_below_min() -> None:
    """V3 requires n_cycles >= 3 (spec §14)."""
    v2 = _silver_v2_for_validation()
    with pytest.raises(ValueError, match="n_cycles"):
        run_v3_3d(
            SILVER_ID,
            SILVER_SEQ,
            SILVER_VINF,
            SILVER_TOF,
            SILVER_REL_OFF,
            None,
            v2_verdict=v2,  # type: ignore[arg-type]
            n_cycles=2,
            n_revs=SILVER_NREV,
            phase0_deg=SILVER_PHASE0,
        )


def test_v3_rejects_non_closed_sequence() -> None:
    """Moontour sequence must be closed."""
    v2 = _silver_v2_for_validation()
    with pytest.raises(ValueError, match="CLOSED"):
        run_v3_3d(
            SILVER_ID,
            ("Umbriel", "Oberon"),
            (1.0, 1.0),
            (14.94,),
            0.0,
            None,
            v2_verdict=v2,  # type: ignore[arg-type]
            n_cycles=3,
        )


def test_v3_rejects_v2_verdict_with_fewer_cycles() -> None:
    """V3 can't compare against a V2 with fewer recorded cycles than requested."""
    v2_two = run_v2_moontour(  # legal V2 only at n_cycles >= 3, but per_cycle truncates if shorter
        SILVER_ID,
        SILVER_SEQ,
        SILVER_VINF,
        SILVER_TOF,
        SILVER_REL_OFF,
        None,
        n_cycles=3,
        n_revs=SILVER_NREV,
        phase0_deg=SILVER_PHASE0,
    )
    with pytest.raises(ValueError, match="v2_verdict"):
        run_v3_3d(
            SILVER_ID,
            SILVER_SEQ,
            SILVER_VINF,
            SILVER_TOF,
            SILVER_REL_OFF,
            None,
            v2_verdict=v2_two,
            n_cycles=10,  # V2 only has 3 per_cycle entries
            n_revs=SILVER_NREV,
            phase0_deg=SILVER_PHASE0,
        )


def test_v3_rejects_mismatched_tof_length() -> None:
    """leg_tofs_days must have len(sequence) - 1 entries."""
    v2 = _silver_v2_for_validation()
    with pytest.raises(ValueError, match="leg_tofs_days"):
        run_v3_3d(
            SILVER_ID,
            SILVER_SEQ,
            SILVER_VINF,
            (14.94,),  # 1 ToF but n_legs=2
            SILVER_REL_OFF,
            None,
            v2_verdict=v2,  # type: ignore[arg-type]
            n_cycles=3,
            n_revs=SILVER_NREV,
        )


def test_v3_rejects_zero_or_negative_agreement_floor() -> None:
    """agreement_floor_kms must be > 0."""
    v2 = _silver_v2_for_validation()
    with pytest.raises(ValueError, match="agreement_floor"):
        run_v3_3d(
            SILVER_ID,
            SILVER_SEQ,
            SILVER_VINF,
            SILVER_TOF,
            SILVER_REL_OFF,
            None,
            v2_verdict=v2,  # type: ignore[arg-type]
            n_cycles=3,
            agreement_floor_kms=0.0,
            n_revs=SILVER_NREV,
        )


# ---------------------------------------------------------------------------
# Test 6: Audit-trail fields preserved.
# ---------------------------------------------------------------------------


def test_v3_audit_fields_preserved() -> None:
    """The verdict carries identification + floors + notes."""
    v2 = run_v2_moontour(
        SILVER_ID,
        SILVER_SEQ,
        SILVER_VINF,
        SILVER_TOF,
        SILVER_REL_OFF,
        None,
        n_cycles=4,
        n_revs=SILVER_NREV,
        phase0_deg=SILVER_PHASE0,
    )
    v3 = run_v3_3d(
        SILVER_ID,
        SILVER_SEQ,
        SILVER_VINF,
        SILVER_TOF,
        SILVER_REL_OFF,
        None,
        v2_verdict=v2,
        n_cycles=4,
        n_revs=SILVER_NREV,
        phase0_deg=SILVER_PHASE0,
        notes="V3 phase-3 audit-fields smoke",
    )
    assert v3.candidate_id == SILVER_ID
    assert v3.sequence == SILVER_SEQ
    assert v3.n_cycles_propagated == 4
    assert v3.v3_v2_agreement_floor_kms == V3_AGREEMENT_FLOOR_KMS
    assert v3.notes == "V3 phase-3 audit-fields smoke"


# ---------------------------------------------------------------------------
# Test 7: Fabrication-guard for spec constants.
# ---------------------------------------------------------------------------


def test_v3_spec_constants() -> None:
    """Spec-fixed constants are 100 km and 3 cycles."""
    assert V3_AGREEMENT_FLOOR_KMS == 100.0
    assert V3_N_CYCLES_MIN == 3
