"""Tests for the V3-3D periodic-orbit independent-integrator gauntlet (#306).

The PARITY TEST gates everything: if REBOUND IAS15 cannot match the project's
scipy DOP853 on a known CR3BP orbit to a tight tolerance over one period, the
V3 verdict is meaningless. That test runs first and is the load-bearing check
for the velocity-dependent rotating-frame force callback (design Risk R1).

Candidate under test: the #444 C21 spatial 2:1 Earth-Moon orbit
(braik_ross_system), which PASSED V1 + V2 in the 3D gauntlet.
"""

from __future__ import annotations

import numpy as np
import pytest

import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.data.validation.v2_3d import run_v2_3d
from cyclerfinder.data.validation.v3_3d_periodic import (
    V3_PERIODIC_CLOSURE_FLOOR_NONDIM,
    _ias15_propagate_cr3bp_rotating,
    run_v3_3d_periodic,
)
from cyclerfinder.search.reachable_representatives import braik_ross_system

# C21 candidate: spatial 2:1 EM orbit, Floquet-stable (#444).
C21_STATE0 = np.array(
    [0.7440212218499672, 0.0, -0.2057098355650995, 0.0, 0.35368280201143637, 0.0],
    dtype=np.float64,
)
C21_PERIOD = 18.167169790651315


def _dop853_one_period(state0: np.ndarray, period: float, system: cr3bp.CR3BPSystem) -> np.ndarray:
    arc = cr3bp.propagate(system, state0, period, with_stm=False, rtol=1e-13, atol=1e-13)
    return np.asarray(arc.state_f, dtype=np.float64)


# --------------------------------------------------------------------------
# 1. PARITY GATE — IAS15 must match DOP853 on a known orbit over one period.
#    This is the mandatory mitigation for the velocity-dependent Coriolis
#    callback (design §4 R1). If it fails, every verdict below is meaningless.
# --------------------------------------------------------------------------


def test_parity_ias15_rotating_callback_matches_dop853_one_period() -> None:
    """Mode (A) rotating-frame callback: IAS15 vs DOP853 terminal state agree."""
    system = braik_ross_system()
    ref = _dop853_one_period(C21_STATE0, C21_PERIOD, system)
    ias15, label = _ias15_propagate_cr3bp_rotating(
        C21_STATE0, C21_PERIOD, system.mu, mode="rotating_callback"
    )
    residual = float(np.linalg.norm(ias15 - ref))
    # Tight parity bar: the two integrator architectures must agree to far
    # below the V3 closure floor (1e-7 nondim) on this known orbit.
    assert residual < 1.0e-9, (
        f"IAS15-vs-DOP853 parity FAILED ({residual:.3e} nondim) with {label}; "
        "the rotating-frame force callback is not faithful — V3 verdict is meaningless"
    )
    assert "REBOUND IAS15" in label or "LSODA" in label


def test_parity_residual_well_under_closure_floor() -> None:
    """The parity residual is comfortably under the V3 closure floor."""
    system = braik_ross_system()
    ref = _dop853_one_period(C21_STATE0, C21_PERIOD, system)
    ias15, _ = _ias15_propagate_cr3bp_rotating(
        C21_STATE0, C21_PERIOD, system.mu, mode="rotating_callback"
    )
    residual = float(np.linalg.norm(ias15 - ref))
    assert residual < V3_PERIODIC_CLOSURE_FLOOR_NONDIM


# --------------------------------------------------------------------------
# 2. C21 V3 VERDICT — wraps the IAS15 cross-check into a frozen verdict.
# --------------------------------------------------------------------------


def test_c21_passes_v3() -> None:
    """C21 PASSES V3: IAS15 reproduces the V1 closure + the V2 drift signature."""
    system = braik_ross_system()
    v2 = run_v2_3d("c21", C21_STATE0, C21_PERIOD, system, n_cycles=3)
    assert v2.passes_v2  # precondition: V2 must have produced a drift series

    v3 = run_v3_3d_periodic("c21", C21_STATE0, C21_PERIOD, system, v2_verdict=v2, n_cycles=3)

    assert v3.passes_v3, (
        f"C21 V3 FAILED: closure={v3.closure_residual_nondim_ias15:.3e} nondim, "
        f"drift_agreement={v3.drift_agreement_kms:.3e} km, integrator={v3.integrator}"
    )
    assert v3.closure_residual_nondim_ias15 < v3.closure_floor_nondim
    assert v3.drift_agreement_kms < v3.drift_agreement_floor_kms
    assert v3.n_cycles_propagated == 3
    assert v3.converged_at_each_return
    assert "REBOUND IAS15" in v3.integrator


def test_c21_drift_series_agrees_with_v2() -> None:
    """The IAS15 per-cycle drift series tracks the V2 (DOP853) series tightly."""
    system = braik_ross_system()
    v2 = run_v2_3d("c21", C21_STATE0, C21_PERIOD, system, n_cycles=3)
    v3 = run_v3_3d_periodic("c21", C21_STATE0, C21_PERIOD, system, v2_verdict=v2, n_cycles=3)
    assert len(v3.per_cycle_drift_kms_ias15) == 3
    assert len(v3.per_cycle_drift_kms_dop853) == 3
    # Every per-cycle pair agrees to well under 1 km.
    for ias15, dop853 in zip(
        v3.per_cycle_drift_kms_ias15, v3.per_cycle_drift_kms_dop853, strict=True
    ):
        assert abs(ias15 - dop853) < 1.0


# --------------------------------------------------------------------------
# 3. NEGATIVE GUARD — a non-periodic IC must FAIL V3 (no vacuous pass).
# --------------------------------------------------------------------------


def test_nonperiodic_ic_fails_v3() -> None:
    """A deliberately non-periodic IC blows the closure floor → FAIL."""
    system = braik_ross_system()
    # Perturb the C21 IC hard so it is no longer periodic at C21_PERIOD.
    bad_state0 = C21_STATE0 + np.array([0.05, 0.0, 0.05, 0.0, 0.05, 0.0])
    # Build the V2 verdict on the same bad IC so the comparison is honest.
    v2 = run_v2_3d("bad", bad_state0, C21_PERIOD, system, n_cycles=3)
    v3 = run_v3_3d_periodic("bad", bad_state0, C21_PERIOD, system, v2_verdict=v2, n_cycles=3)
    # The one-period closure must blow the floor (the IC does not return).
    assert not v3.passes_v3
    assert v3.closure_residual_nondim_ias15 > v3.closure_floor_nondim


# --------------------------------------------------------------------------
# 4. INPUT VALIDATION
# --------------------------------------------------------------------------


def test_rejects_sub_minimum_cycles() -> None:
    system = braik_ross_system()
    v2 = run_v2_3d("c21", C21_STATE0, C21_PERIOD, system, n_cycles=3)
    with pytest.raises(ValueError, match="n_cycles >= 3"):
        run_v3_3d_periodic("c21", C21_STATE0, C21_PERIOD, system, v2_verdict=v2, n_cycles=2)


def test_rejects_bad_state_shape() -> None:
    system = braik_ross_system()
    v2 = run_v2_3d("c21", C21_STATE0, C21_PERIOD, system, n_cycles=3)
    with pytest.raises(ValueError, match="6D"):
        run_v3_3d_periodic("c21", np.zeros(5), C21_PERIOD, system, v2_verdict=v2)


def test_rejects_v2_with_too_few_cycles() -> None:
    system = braik_ross_system()
    v2 = run_v2_3d("c21", C21_STATE0, C21_PERIOD, system, n_cycles=3)
    with pytest.raises(ValueError, match="cycles"):
        run_v3_3d_periodic("c21", C21_STATE0, C21_PERIOD, system, v2_verdict=v2, n_cycles=4)
