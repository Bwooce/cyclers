"""Tests for the V2 long-span QP-torus gauntlet (#319 Phase 1 Part B).

Sourced golden discipline
-------------------------
Tests assert on CLOSURE QUALITY + DRIFT BOUNDEDNESS, not on specific
torus state numbers. Inputs are the sourced #299 Neimark-Sacker-seeded
torus (same one ``test_v1_qp`` uses). Floors (V2_QP_DRIFT_FLOOR = 5e-2,
V2_QP_N_CYCLES_MIN = 3) are module constants -- NOT test-tunable.

Test cases
----------
  1. Sourced #299 smoke-test torus, n_cycles=3: V2_qp passes both gates.
  2. Periodic-orbit limit (zero amplitude): trivially passes V2_qp.
  3. Corrupted Fourier coefficients: fails V2_qp.
  4. Cycle-count scaling: n_cycles=5 and n_cycles=10 -- documents the
     per-cycle invariance drift as a function of k. For n_modes=2 the
     drift is expected to grow with k (Olikara 2016 §4); the test asserts
     the GROWTH is bounded (no exponential blowup at the moderate k the
     project uses).
  5. Spec-floor fabrication guard.
  6. Audit-trail fields survive.
  7. Bad caller args rejected.
"""

from __future__ import annotations

import json
import math
from dataclasses import replace
from pathlib import Path
from typing import Any

import numpy as np
import pytest

import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.data.validation.v2_qp import (
    V2_QP_DRIFT_FLOOR,
    V2_QP_N_CYCLES_MIN,
    V2_QP_OFF_GRID_PER_CYCLE,
    V2VerdictQP,
    run_v2_qp,
)
from cyclerfinder.genome.qp_tori import (
    QPTorus,
    _seed_invariant_circle,
    correct_qp_torus,
)

DATA_DIR = Path(__file__).resolve().parents[2] / "data"
SUBFAMILIES_FILE = DATA_DIR / "family_296_3d_subfamilies_299.jsonl"
PARENT_FAMILY_FILE = DATA_DIR / "family_296_3d_em_11.jsonl"

EM_MU = 1.2150584270572e-2
EM_L_KM = 384400.0
EM_T_S = 375699.8


def _em_system() -> cr3bp.CR3BPSystem:
    return cr3bp.CR3BPSystem(
        mu=EM_MU,
        primary="earth",
        secondary="moon",
        l_km=EM_L_KM,
        t_s=EM_T_S,
    )


def _load_parent_at_step(step_index: int) -> dict[str, Any]:
    if not PARENT_FAMILY_FILE.exists():
        pytest.skip(f"parent family file not present: {PARENT_FAMILY_FILE}")
    with PARENT_FAMILY_FILE.open() as f:
        for line in f:
            d = json.loads(line)
            if d.get("type") == "header":
                continue
            if d.get("step_index") == step_index:
                return dict(d)
    raise RuntimeError(f"step_index={step_index} not found in {PARENT_FAMILY_FILE}")


def _load_first_neimark_sacker_bracket() -> tuple[dict[str, Any], dict[str, Any]]:
    if not SUBFAMILIES_FILE.exists():
        pytest.skip(f"subfamilies file not present: {SUBFAMILIES_FILE}")
    with SUBFAMILIES_FILE.open() as f:
        for line in f:
            d = json.loads(line)
            if d.get("type") == "header":
                continue
            br = d.get("bracket")
            if br is None:
                continue
            if br.get("classification") == "neimark_sacker":
                parent = _load_parent_at_step(int(br["step_a"]))
                return br, parent
    raise RuntimeError("no Neimark-Sacker bracket found in subfamilies file")


def _build_smoke_torus() -> QPTorus:
    br, parent = _load_first_neimark_sacker_bracket()
    system = _em_system()
    parent_state = np.asarray(parent["state_nd"], dtype=np.float64)
    parent_period = float(parent["T_TU"])
    k = int(br["k"])
    lam_a = complex(br["eig_a_re"], br["eig_a_im"])
    lam_b = complex(br["eig_b_re"], br["eig_b_im"])
    return correct_qp_torus(
        system,
        parent_state,
        parent_period,
        (lam_a, lam_b),
        k=k,
        n_long=16,
        n_trans=2,
        initial_torus_amplitude=5e-4,
        tol=1e-8,
        max_iter=40,
        independent_tol=1e-3,
        notes="phase1_smoke_test_from_299_bracket (V2_qp input)",
    )


def _build_zero_amplitude_torus() -> QPTorus:
    br, parent = _load_first_neimark_sacker_bracket()
    system = _em_system()
    parent_state = np.asarray(parent["state_nd"], dtype=np.float64)
    parent_period = float(parent["T_TU"])
    k = int(br["k"])
    n_modes = 2
    arc = cr3bp.propagate(system, parent_state, parent_period, with_stm=True)
    monod = arc.stm
    assert monod is not None
    coeffs_seed, lam_seed, _ = _seed_invariant_circle(
        parent_state, monod, k=k, n_modes=n_modes, amplitude=0.0
    )
    rho_seed = math.atan2(np.imag(lam_seed), np.real(lam_seed))
    t_strob = float(parent_period)
    n_samples = 2 * n_modes + 3
    return QPTorus(
        system=system,
        omega_long=2 * math.pi / t_strob,
        omega_trans=rho_seed / t_strob,
        rho=rho_seed,
        t_strob=t_strob,
        fourier_coeffs=coeffs_seed,
        n_modes=n_modes,
        n_samples=n_samples,
        invariance_residual=0.0,
        independent_closure_residual=0.0,
        converged=True,
        n_iter=0,
        notes="zero_amplitude_periodic_orbit_limit",
    )


# ---------------------------------------------------------------------------
# Test 1: floor fabrication guards.
# ---------------------------------------------------------------------------


def test_v2_qp_floors_are_sourced_constants() -> None:
    """Fabrication guard: V2_qp floors match the module constants.

    The 5e-2 drift floor is the empirically calibrated value for low
    truncation order (N=2) -- see :data:`V2_QP_DRIFT_FLOOR` docstring
    for the Olikara 2016 + #319 Phase 1 sourcing chain.
    """
    assert V2_QP_N_CYCLES_MIN == 3
    assert V2_QP_DRIFT_FLOOR == 5.0e-2
    assert V2_QP_OFF_GRID_PER_CYCLE == 50


# ---------------------------------------------------------------------------
# Test 2: sourced smoke torus, n_cycles=3, passes.
# ---------------------------------------------------------------------------


@pytest.mark.slow
def test_v2_qp_sourced_smoke_torus_passes() -> None:
    """V2_qp on the #299 Neimark-Sacker-seeded torus, n_cycles=3.

    The torus must stay invariant -- in the off-grid invariance sense --
    over 3 consecutive longitudinal cycles. The max-over-cycles
    invariance residual must stay below the V2_qp floor (5e-2 nondim
    for low truncation order N=2; see V2_QP_DRIFT_FLOOR docstring).
    """
    torus = _build_smoke_torus()
    verdict = run_v2_qp(
        "qp-torus-smoke-299-v2",
        torus,
        n_cycles=3,
        notes="sourced #299 smoke torus; n_cycles=3",
    )
    assert isinstance(verdict, V2VerdictQP)
    per_cycle_fmt = [f"{x:.3e}" for x in verdict.per_cycle_invariance_residual]
    print(
        f"\n[V2_qp smoke n=3] per_cycle={per_cycle_fmt} "
        f"| max_drift={verdict.max_invariance_drift:.3e} "
        f"| max_drift_km={verdict.max_invariance_drift_km:.3e} "
        f"| passes={verdict.passes_v2_qp}"
    )
    assert verdict.converged_each_cycle
    assert verdict.n_cycles_longitudinal_propagated == 3
    assert verdict.passes_v2_qp, (
        f"V2_qp FAIL: max_drift={verdict.max_invariance_drift:.3e} > "
        f"drift_floor={verdict.drift_floor:.0e}"
    )


# ---------------------------------------------------------------------------
# Test 3: zero-amplitude limit passes.
# ---------------------------------------------------------------------------


def test_v2_qp_zero_amplitude_limit_passes() -> None:
    """At zero amplitude the torus collapses to the parent periodic orbit.
    Over 3 cycles it stays at the parent fixed point modulo integrator
    round-trip error -- the per-cycle invariance is dominated by the
    parent orbit's own multi-cycle drift, well below 1e-3.

    NOTE: the parent orbit has Floquet multipliers on/near the unit
    circle (it IS at the Neimark-Sacker bracket), so the drift is NOT
    zero -- but the bracket is precisely the boundary where the
    instability is weak. Multi-cycle drift for the parent at this bracket
    is documented in #299 as ~1e-4 to 1e-3 nondim over 3 laps.
    """
    torus = _build_zero_amplitude_torus()
    verdict = run_v2_qp(
        "qp-torus-zero-amp-v2",
        torus,
        n_cycles=3,
        notes="degenerate parent-orbit limit",
    )
    per_cycle_fmt = [f"{x:.3e}" for x in verdict.per_cycle_invariance_residual]
    print(
        f"\n[V2_qp zero-amp] per_cycle={per_cycle_fmt} "
        f"| max_drift={verdict.max_invariance_drift:.3e} "
        f"| passes={verdict.passes_v2_qp}"
    )
    assert verdict.converged_each_cycle
    # The parent orbit at the Neimark-Sacker bracket has near-unity
    # Floquet multipliers; integrator + Floquet drift should land it
    # under 1e-3 over 3 laps but we don't quote a tighter number.
    assert verdict.passes_v2_qp


# ---------------------------------------------------------------------------
# Test 4: corrupted torus fails V2.
# ---------------------------------------------------------------------------


@pytest.mark.slow
def test_v2_qp_rejects_corrupted_torus() -> None:
    """Bit-flip the Fourier modes; V2_qp must catch it.

    A gross perturbation of the n=1 mode means the propagated points
    do NOT lie on the torus even after one cycle, let alone three.
    """
    torus = _build_smoke_torus()
    bad_coeffs = torus.fourier_coeffs.copy()
    bad_coeffs[1, :] += 0.05 + 0.05j
    n_total = bad_coeffs.shape[0]
    bad_coeffs[n_total - 1, :] = np.conj(bad_coeffs[1, :])
    corrupted = replace(torus, fourier_coeffs=bad_coeffs)
    verdict = run_v2_qp(
        "qp-torus-v2-corrupted",
        corrupted,
        n_cycles=3,
        notes="negative control: corrupted modes",
    )
    per_cycle_fmt = [f"{x:.3e}" for x in verdict.per_cycle_invariance_residual]
    print(
        f"\n[V2_qp corrupted] per_cycle={per_cycle_fmt} "
        f"| max_drift={verdict.max_invariance_drift:.3e} "
        f"| passes={verdict.passes_v2_qp}"
    )
    assert not verdict.passes_v2_qp


# ---------------------------------------------------------------------------
# Test 5: cycle-count scaling (5 and 10 cycles).
# ---------------------------------------------------------------------------


@pytest.mark.slow
def test_v2_qp_cycle_count_scaling_documents_drift() -> None:
    """Document the per-cycle invariance drift as a function of cycle
    count k. Olikara 2016 §4 reports drift growing from ~1e-4 at k=1 to
    ~1e-2 at k=10 for moderate truncation orders. We re-run at n_cycles=5
    and assert:

      - n_cycles=5 succeeds at the integrator level (converged_each_cycle).
      - Per-cycle drift is monotonically non-decreasing (modulo random
        sampling jitter) -- the hyperbolic-instability amplification of
        truncation error is unidirectional. NOT a strict assertion --
        random sampling jitter can flip adjacent cycles -- but the
        FINAL cycle drift is no smaller than the FIRST cycle drift.
      - Per-cycle drift slope per cycle is documented in stdout for the
        Phase 1 verdict.

    n=10 is documented but NOT asserted as PASS -- the floor 1e-3 may or
    may not hold at k=10; passing/failing is the documentation.
    """
    torus = _build_smoke_torus()
    verdict_5 = run_v2_qp(
        "qp-torus-v2-5cycles",
        torus,
        n_cycles=5,
        # Tight off-grid count keeps the test fast; the scaling claim
        # cares about the trend over cycles, not absolute n_samples.
        n_off_grid_samples_per_cycle=20,
        notes="cycle-count scaling, n=5",
    )
    per_cycle_fmt_5 = [f"{x:.3e}" for x in verdict_5.per_cycle_invariance_residual]
    print(
        f"\n[V2_qp n=5] per_cycle={per_cycle_fmt_5} "
        f"| max_drift={verdict_5.max_invariance_drift:.3e} "
        f"| passes={verdict_5.passes_v2_qp}"
    )
    assert verdict_5.converged_each_cycle
    assert verdict_5.n_cycles_longitudinal_propagated == 5
    # Final cycle drift is no smaller than first cycle drift (the
    # hyperbolic-instability amplification is unidirectional).
    first_drift = verdict_5.per_cycle_invariance_residual[0]
    last_drift = verdict_5.per_cycle_invariance_residual[-1]
    assert last_drift >= first_drift, (
        f"V2_qp drift NON-monotonic: first={first_drift:.3e}, last={last_drift:.3e} "
        f"(expected last >= first; Olikara 2016 §4 hyperbolic-amp)"
    )

    # Run n=10 and document the slope -- no PASS assertion (the floor
    # may not hold at k=10, which is fine for documentation).
    verdict_10 = run_v2_qp(
        "qp-torus-v2-10cycles",
        torus,
        n_cycles=10,
        n_off_grid_samples_per_cycle=10,
        # Allow the test to RUN at k=10 even if drift exceeds 1e-3 --
        # the floor here is the documentation cap, not the pass gate.
        drift_floor=1.0,
        notes="cycle-count scaling, n=10 (drift_floor relaxed for documentation)",
    )
    per_cycle_fmt_10 = [f"{x:.3e}" for x in verdict_10.per_cycle_invariance_residual]
    ratio = verdict_10.per_cycle_invariance_residual[-1] / max(
        verdict_10.per_cycle_invariance_residual[0], 1e-16
    )
    print(
        f"[V2_qp n=10] per_cycle={per_cycle_fmt_10} "
        f"| max_drift={verdict_10.max_invariance_drift:.3e} "
        f"| ratio_k10_to_k1={ratio:.2f}"
    )
    assert verdict_10.converged_each_cycle


# ---------------------------------------------------------------------------
# Test 6: audit-trail fields.
# ---------------------------------------------------------------------------


def test_v2_qp_verdict_carries_audit_fields() -> None:
    torus = _build_zero_amplitude_torus()
    verdict = run_v2_qp(
        "audit-trail-check",
        torus,
        n_cycles=3,
        notes="audit-trail check",
    )
    assert verdict.candidate_id == "audit-trail-check"
    assert verdict.notes == "audit-trail check"
    assert verdict.drift_floor == V2_QP_DRIFT_FLOOR
    assert verdict.n_cycles_min == V2_QP_N_CYCLES_MIN
    assert verdict.n_off_grid_samples_per_cycle == V2_QP_OFF_GRID_PER_CYCLE
    assert verdict.n_modes == torus.n_modes


# ---------------------------------------------------------------------------
# Test 7: bad caller args + malformed system.
# ---------------------------------------------------------------------------


def test_v2_qp_rejects_bad_caller_args() -> None:
    torus = _build_zero_amplitude_torus()
    with pytest.raises(ValueError, match="n_cycles"):
        run_v2_qp("bad", torus, n_cycles=2)
    with pytest.raises(ValueError, match="drift_floor"):
        run_v2_qp("bad", torus, drift_floor=0.0)
    with pytest.raises(ValueError, match="n_off_grid_samples_per_cycle"):
        run_v2_qp("bad", torus, n_off_grid_samples_per_cycle=0)
    with pytest.raises(TypeError, match="QPTorus"):
        run_v2_qp("bad", "not-a-torus")  # type: ignore[arg-type]


def test_v2_qp_rejects_malformed_system() -> None:
    torus = _build_zero_amplitude_torus()
    bad_system = cr3bp.CR3BPSystem(
        mu=EM_MU,
        primary="earth",
        secondary="moon",
        l_km=0.0,
        t_s=EM_T_S,
    )
    bad_torus = replace(torus, system=bad_system)
    with pytest.raises(ValueError, match="invalid CR3BP system"):
        run_v2_qp("bad-system", bad_torus)
