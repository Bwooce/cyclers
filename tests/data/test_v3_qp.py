"""#319 Phase 2 / #320 follow-on — V3 independent-integrator QP-torus gauntlet.

V3_qp re-checks a QP-torus's Fourier-mode invariance under REBOUND IAS15 (an
independent integrator from the V1/V2 scipy DOP853 flow) and cross-checks the IAS15
endpoint against DOP853. A PASS means the V1/V2 torus signature is a real CR3BP
property, not a DOP853 artifact. Mirrors `test_v2_qp` / `v3_3d_periodic`'s IAS15
cross-check. The sourced smoke torus is the #299 Neimark-Sacker bracket.
"""

from __future__ import annotations

import pytest

from cyclerfinder.data.validation.v2_qp import V2_QP_DRIFT_FLOOR, V2_QP_N_CYCLES_MIN
from cyclerfinder.data.validation.v3_qp import (
    V3_QP_AGREEMENT_FLOOR,
    V3VerdictQP,
    run_v3_qp,
)

# Reuse the V2_qp fixture builders (same #299 bracket source, no duplication).
from tests.data.test_v2_qp import _build_smoke_torus, _build_zero_amplitude_torus


def test_v3_qp_floors_are_sourced_constants() -> None:
    """The V3_qp floors are the documented constants (not silently fabricated)."""
    # V3 reuses V2's empirically-calibrated invariance floor (V3 not more permissive).
    assert pytest.approx(5e-2) == V2_QP_DRIFT_FLOOR
    # The integrator-agreement floor is one order under the invariance floor.
    assert pytest.approx(1e-3) == V3_QP_AGREEMENT_FLOOR
    assert 0.0 < V3_QP_AGREEMENT_FLOOR < V2_QP_DRIFT_FLOOR


def test_v3_qp_rejects_bad_caller_args() -> None:
    """Argument validation mirrors V2_qp (fail loud, never silently mis-gate)."""
    torus = _build_zero_amplitude_torus()
    with pytest.raises(TypeError):
        run_v3_qp("x", object())  # type: ignore[arg-type]
    with pytest.raises(ValueError):
        run_v3_qp("x", torus, n_cycles=V2_QP_N_CYCLES_MIN - 1)
    with pytest.raises(ValueError):
        run_v3_qp("x", torus, drift_floor=0.0)
    with pytest.raises(ValueError):
        run_v3_qp("x", torus, agreement_floor=0.0)
    with pytest.raises(ValueError):
        run_v3_qp("x", torus, n_off_grid_samples_per_cycle=0)


def test_v3_qp_zero_amplitude_limit_passes() -> None:
    """At the periodic-orbit limit (amplitude 0) IAS15 reproduces DOP853 → V3_qp PASS.

    The invariant circle collapses to the parent periodic orbit; IAS15 and DOP853
    integrate the same CR3BP EOM, so both the invariance residual and the
    integrator disagreement are tiny. Fast (no torus regeneration).
    """
    torus = _build_zero_amplitude_torus()
    v = run_v3_qp("qp-zero-amp", torus, n_off_grid_samples_per_cycle=8)
    assert isinstance(v, V3VerdictQP)
    assert v.n_cycles_longitudinal_propagated >= V2_QP_N_CYCLES_MIN
    assert v.max_invariance_drift_ias15 <= v.drift_floor, v.per_cycle_invariance_residual_ias15
    assert v.max_integrator_disagreement <= v.agreement_floor, v.per_cycle_integrator_disagreement
    assert v.passes_v3_qp
    # The IAS15 label records which propagator actually ran.
    assert "IAS15" in v.integrator or "LSODA" in v.integrator


def test_v3_qp_verdict_carries_audit_fields() -> None:
    """The verdict serialises the integrator + both metrics for the audit trail."""
    torus = _build_zero_amplitude_torus()
    v = run_v3_qp("qp-audit", torus, n_off_grid_samples_per_cycle=6)
    assert v.candidate_id == "qp-audit"
    assert len(v.per_cycle_invariance_residual_ias15) == v.n_cycles_longitudinal_propagated
    assert len(v.per_cycle_integrator_disagreement) == v.n_cycles_longitudinal_propagated
    assert v.max_invariance_drift_ias15_km >= 0.0
    assert v.max_integrator_disagreement_km >= 0.0
    assert v.n_modes == torus.n_modes


@pytest.mark.slow
def test_v3_qp_sourced_smoke_torus_passes() -> None:
    """The #299 sourced smoke torus passes V3_qp: IAS15 reproduces the V2 invariance.

    Load-bearing evidence (slow: ~20s torus regeneration). IAS15 invariance tracks
    DOP853's V2 invariance and the integrators agree to << the floor — so the torus
    signature is integrator-independent.
    """
    torus = _build_smoke_torus()
    v = run_v3_qp("qp-torus-smoke-299-v3", torus, n_off_grid_samples_per_cycle=12)
    assert v.passes_v3_qp, (
        f"inv={v.max_invariance_drift_ias15:.3e}/{v.drift_floor} "
        f"disagree={v.max_integrator_disagreement:.3e}/{v.agreement_floor}"
    )
    # IAS15 and DOP853 agree to well under the floor (signature is real, not noise).
    assert v.max_integrator_disagreement < 1e-4, v.per_cycle_integrator_disagreement
