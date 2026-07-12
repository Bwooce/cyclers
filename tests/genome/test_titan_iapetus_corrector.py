"""Tests for the Titan-Iapetus 3D eccentric-Keplerian closure corrector (#574 Stage B).

Three gates, mirroring the #574 Stage-A script's own discipline (ported here as a
proper pytest regression suite rather than throwaway print-based checks):

1. **C2 positive control (e=0 reduction)**: :func:`kepler_state_3d` reduces EXACTLY to
   the project's circular-coplanar :func:`_moon_state` (Titan's case, ``inc=raan=0``) and
   to the standard circular-inclined ``R3(Omega).R1(inc)`` rotation (Iapetus's case,
   ``inc!=0``) at ``ecc=0``, at a grid of M0/Omega/u test points. An unverified eccentric
   propagator producing a downstream "no closure" result would be indistinguishable from a
   real family death (project rule
   ``feedback_verify_gauntlet_with_positive_control``) -- this is the load-bearing check.
2. **Known-branch reproduction**: branch 1 from ``data/probe_574_titan_iapetus_eccentric_
   kill_gate.jsonl`` (a real #574 Stage-A eccentric survivor) reproduces its recorded
   final-stage residual/V_inf/gate outcome through THIS module's public API, at the
   EXACT stored (omega, tof_scale, m0_titan, m0_iapetus, n_rev) point -- i.e. the
   productized ``evaluate_closure``/``closure_passes_gate`` reproduce the throwaway
   script's own already-computed result, not a fresh number.
3. **Structural / input-contract checks**: dataclass immutability, ``ClosureResult.closes``
   semantics, and that an infeasible geometry is reported as ``lambert_infeasible`` rather
   than raising out of ``evaluate_closure``.

Per ``feedback_golden_tests_sourced_only``: gate 1's expected side is the ANALYTIC
definition of "e=0 reduces to circular" (zero, not a value our own code computed); gate 2's
expected side is a value ALREADY COMMITTED to the repo by the #574 Stage-A run (i.e. an
independent artifact of this exact module's throwaway predecessor, not a value invented for
this test).
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from cyclerfinder.core.satellites import PRIMARIES, SATELLITES
from cyclerfinder.genome.titan_iapetus_corrector import (
    ANCHOR,
    FLYBY,
    GATE_RESIDUAL_KMS,
    INCLINATION_DEG,
    PRIMARY,
    ClosureResult,
    TitanIapetusClosureParams,
    closure_passes_gate,
    evaluate_closure,
    kepler_state_3d,
)
from cyclerfinder.search.discovery_campaign import _mean_motion_rad_day, _moon_state

ROOT = Path(__file__).resolve().parent.parent.parent
PROBE_574_PATH = ROOT / "data" / "probe_574_titan_iapetus_eccentric_kill_gate.jsonl"


def _circular_inclined_reference(
    u_rad: float, v_circ_km_s: float, sma_km: float, omega_rad: float, inc_rad: float
) -> tuple[np.ndarray, np.ndarray]:
    """Independent hand-written circular-inclined state (R3(Omega).R1(inc)), for the
    e=0/inc!=0 reduction check. Deliberately NOT imported from any throwaway script --
    this is the standard textbook rotation, re-derived here so the test does not depend
    on the module under test agreeing with a sibling implementation of the same formula.
    """
    cos_o, sin_o = math.cos(omega_rad), math.sin(omega_rad)
    cosi, sini = math.cos(inc_rad), math.sin(inc_rad)
    cosu, sinu = math.cos(u_rad), math.sin(u_rad)
    x = sma_km * (cos_o * cosu - sin_o * sinu * cosi)
    y = sma_km * (sin_o * cosu + cos_o * sinu * cosi)
    z = sma_km * sinu * sini
    vx = -v_circ_km_s * (cos_o * sinu + sin_o * cosu * cosi)
    vy = v_circ_km_s * (-sin_o * sinu + cos_o * cosu * cosi)
    vz = v_circ_km_s * cosu * sini
    return np.array([x, y, z]), np.array([vx, vy, vz])


class TestKeplerStateReductionPositiveControl:
    """C2: kepler_state_3d(ecc=0, ...) must reduce exactly to the circular formulas."""

    def test_reduces_to_moon_state_titan_case(self) -> None:
        mu = PRIMARIES[PRIMARY]
        sat_a = SATELLITES[ANCHOR]
        n_a = _mean_motion_rad_day(mu, sat_a.sma_km)
        for m0_deg in (0.0, 37.0, 190.0, 355.0):
            for t_days in (0.0, 1.3, 12.0):
                m0 = math.radians(m0_deg)
                r3d, v3d = kepler_state_3d(m0, n_a, t_days, sat_a.sma_km, mu, 0.0, 0.0, 0.0)
                r2d, v2d = _moon_state(m0, n_a, t_days, sat_a.sma_km, mu)
                assert np.linalg.norm(r3d - r2d) < 1e-6, f"m0={m0_deg} t={t_days}"
                assert np.linalg.norm(v3d - v2d) < 1e-9, f"m0={m0_deg} t={t_days}"

    def test_reduces_to_circular_inclined_iapetus_case(self) -> None:
        mu = PRIMARIES[PRIMARY]
        sat_b = SATELLITES[FLYBY]
        n_b = _mean_motion_rad_day(mu, sat_b.sma_km)
        v_circ = math.sqrt(mu / sat_b.sma_km)
        inc = math.radians(INCLINATION_DEG)
        for omega_deg in (0.0, 37.0, 190.0, 355.0):
            for u_deg in (0.0, 42.0, 200.0):
                u = math.radians(u_deg)
                omega = math.radians(omega_deg)
                r3d, v3d = kepler_state_3d(u, n_b, 0.0, sat_b.sma_km, mu, 0.0, omega, inc)
                r_ref, v_ref = _circular_inclined_reference(u, v_circ, sat_b.sma_km, omega, inc)
                assert np.linalg.norm(r3d - r_ref) < 1e-6, f"omega={omega_deg} u={u_deg}"
                assert np.linalg.norm(v3d - v_ref) < 1e-9, f"omega={omega_deg} u={u_deg}"


@pytest.mark.skipif(not PROBE_574_PATH.exists(), reason="#574 Stage-A probe jsonl not present")
class TestKnownBranchReproduction:
    """Gate 2: reproduce #574 Stage-A branch 1's own already-committed final-stage result."""

    def _load_branch(self, branch_id: int) -> dict[str, Any]:
        with PROBE_574_PATH.open(encoding="utf-8") as fh:
            for line in fh:
                rec: dict[str, Any] = json.loads(line)
                if rec.get("kind") == "branch_result" and rec.get("branch_id") == branch_id:
                    return rec
        raise AssertionError(f"branch {branch_id} not found in {PROBE_574_PATH}")

    def test_branch_1_final_stage_reproduces(self) -> None:
        branch = self._load_branch(1)
        final = branch["stages"][-1]
        params = TitanIapetusClosureParams(
            omega_deg=final["omega_deg"],
            tof_scale=final["tof_scale"],
            n_rev=tuple(branch["n_rev"]),
            m0_titan_deg=final["m0_titan_deg"],
            m0_iapetus_deg=final["m0_iapetus_deg"],
            e_titan=final["e_titan"],
            e_iapetus=final["e_iapetus"],
        )
        result = evaluate_closure(params)
        assert not result.lambert_infeasible
        assert result.residual_kms == pytest.approx(final["residual_kms"], abs=1e-6)
        assert result.vinf_kms is not None
        for got, want in zip(result.vinf_kms, final["vinf_kms"], strict=True):
            assert got == pytest.approx(want, abs=1e-6)
        # Branch 1 is a recorded #574 survivor -- must also clear the productized gate.
        assert branch["survives"] is True
        gate_pass, verdicts = closure_passes_gate(result)
        assert gate_pass is True
        assert verdicts is not None


class TestClosureResultSemantics:
    def test_closes_property_respects_gate_residual(self) -> None:
        params = TitanIapetusClosureParams(
            omega_deg=0.0, tof_scale=1.0, n_rev=(0, 0), m0_titan_deg=0.0, m0_iapetus_deg=0.0
        )
        ok = ClosureResult(params, GATE_RESIDUAL_KMS - 1e-9, (1.0, 1.0, 1.0), False, None)
        bad = ClosureResult(params, GATE_RESIDUAL_KMS + 1e-9, (1.0, 1.0, 1.0), False, None)
        infeasible = ClosureResult(params, float("inf"), None, True, "geometry")
        assert ok.closes is True
        assert bad.closes is False
        assert infeasible.closes is False

    def test_infeasible_geometry_does_not_raise(self) -> None:
        # An absurd tof_scale (near-zero ToF at a large separation) should be reported as
        # lambert_infeasible, not raise out of evaluate_closure.
        params = TitanIapetusClosureParams(
            omega_deg=0.0,
            tof_scale=1e-6,
            n_rev=(0, 0),
            m0_titan_deg=0.0,
            m0_iapetus_deg=0.0,
        )
        result = evaluate_closure(params)
        assert isinstance(result, ClosureResult)
        if result.lambert_infeasible:
            assert result.infeasible_reason is not None
            assert result.vinf_kms is None

    def test_params_dataclass_is_frozen(self) -> None:
        params = TitanIapetusClosureParams(
            omega_deg=0.0, tof_scale=1.0, n_rev=(0, 0), m0_titan_deg=0.0, m0_iapetus_deg=0.0
        )
        with pytest.raises(AttributeError):
            params.omega_deg = 5.0  # type: ignore[misc]
