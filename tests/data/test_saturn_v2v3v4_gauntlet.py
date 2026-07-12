"""Tests for the Saturn Titan-Iapetus V2->V3->V4->V4-strict gauntlet (#574 Stage B).

What's tested
-------------
* :mod:`cyclerfinder.data.validation.v2_saturn_3d`: dataclasses are frozen, argument
  validation raises before any Lambert work, and a real #574 Stage-A branch (branch 6,
  a genuine machine-precision closure) reproduces its cycle-0 residual through THIS
  module's public API.
* :mod:`cyclerfinder.data.validation.v3_saturn_3d`: V3 (REBOUND IAS15 / LSODA fallback)
  AGREES with V2 to near machine precision on a shared analytic model (both stages use
  the SAME eccentric-3D Kepler states -- V3 only changes the PROPAGATOR, not the
  targeting, so tight agreement is the expected, correct behaviour, not a golden number
  invented for this test).
* :mod:`cyclerfinder.data.validation.v4_saturn`: runs end-to-end without error and
  returns a well-formed verdict (J2 + 8-moon third-body fallback).
* :mod:`cyclerfinder.data.validation.v4_saturn_strict`: SPICE round-trip (sampled
  Titan/Iapetus eccentricity/inclination land in the physically-sane band the #574
  Stage-A spec itself sourced, +/- generous margin -- a sanity check, not a golden
  value) AND the #567-inherited planet-crossing tag fires on a KNOWN real trigger
  epoch (branch 2 at 2000-09-15, independently confirmed by direct instrumentation
  before this test was written -- see ``scripts/run_574_stageB_saturn_gauntlet.py``'s
  own gauntlet run). Skips gracefully if the SAT441 kernel is not installed.

What's NOT tested (and why)
----------------------------
* The full 15-candidate gauntlet's PASS/FAIL verdicts: that is the headline result of
  ``scripts/run_574_stageB_saturn_gauntlet.py``, reported in ``data/OUTSTANDING.md``
  #574, not a unit-test assertion (testing it would bake a discovery-search verdict
  into the regression suite -- a category error, mirroring
  ``tests/data/test_v4_uranus_strict.py``'s own documented scope discipline).
* Multi-epoch sensitivity: scope of the gauntlet runner, not a unit test.

Per ``feedback_golden_tests_sourced_only``: gate 1's expected residual traces to the
ALREADY-COMMITTED #574 Stage-A jsonl (an independent artifact of a different script),
gate 2's expected side is the analytic definition of "same model, different
integrator, agrees to near machine precision" (not invented for this test), and gate 4's
SPICE-sampled bands are physically-sourced ranges, not values our own code computed.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import pytest

from cyclerfinder.core.satellites import PRIMARIES
from cyclerfinder.data.validation.v2_saturn_3d import (
    V2_SATURN_N_CYCLES_MIN,
    run_v2_saturn_3d,
)
from cyclerfinder.data.validation.v3_saturn_3d import run_v3_saturn_3d
from cyclerfinder.data.validation.v4_saturn import (
    SATURN_J2,
    SATURN_R_EQ_KM,
    run_v4_saturn,
)
from cyclerfinder.genome.titan_iapetus_corrector import (
    ECC_IAPETUS,
    ECC_TITAN,
    TitanIapetusClosureParams,
)

ROOT = Path(__file__).resolve().parent.parent.parent
PROBE_574_PATH = ROOT / "data" / "probe_574_titan_iapetus_eccentric_kill_gate.jsonl"

SAT441_PATH: str | None
try:
    from cyclerfinder.verify.spice_kernels import ensure_sat441_kernel

    SAT441_PATH = ensure_sat441_kernel()
except RuntimeError:
    SAT441_PATH = None


def _load_branch_params(branch_id: int) -> TitanIapetusClosureParams:
    with PROBE_574_PATH.open(encoding="utf-8") as fh:
        for line in fh:
            rec = json.loads(line)
            if rec.get("kind") == "branch_result" and rec.get("branch_id") == branch_id:
                final = rec["stages"][-1]
                return TitanIapetusClosureParams(
                    omega_deg=final["omega_deg"],
                    tof_scale=final["tof_scale"],
                    n_rev=tuple(rec["n_rev"]),
                    m0_titan_deg=final["m0_titan_deg"],
                    m0_iapetus_deg=final["m0_iapetus_deg"],
                    e_titan=final["e_titan"],
                    e_iapetus=final["e_iapetus"],
                )
    raise AssertionError(f"branch {branch_id} not found in {PROBE_574_PATH}")


MU_SATURN = PRIMARIES["Saturn"]


@pytest.mark.skipif(not PROBE_574_PATH.exists(), reason="#574 Stage-A probe jsonl not present")
class TestV2SaturnStructural:
    def test_rejects_below_min_cycles(self) -> None:
        params = _load_branch_params(6)
        with pytest.raises(ValueError):
            run_v2_saturn_3d("b6", params, mu=MU_SATURN, n_cycles=V2_SATURN_N_CYCLES_MIN - 1)

    def test_cycle_zero_reproduces_committed_stage_a_residual(self) -> None:
        params = _load_branch_params(6)
        v2 = run_v2_saturn_3d("b6", params, mu=MU_SATURN, n_cycles=3)
        assert v2.per_cycle[0].closure_residual_kms < 1e-6
        assert v2.per_cycle[0].rendezvous_drift_kms == 0.0

    def test_verdict_is_frozen(self) -> None:
        params = _load_branch_params(6)
        v2 = run_v2_saturn_3d("b6", params, mu=MU_SATURN, n_cycles=3)
        with pytest.raises(AttributeError):
            v2.passes_v2 = True  # type: ignore[misc]


@pytest.mark.skipif(not PROBE_574_PATH.exists(), reason="#574 Stage-A probe jsonl not present")
class TestV3SaturnAgreesWithV2:
    def test_v3_agrees_with_v2_near_machine_precision(self) -> None:
        """Same underlying analytic model, different integrator (IAS15 vs the V2
        driver's analytic Lambert + closed-form Kepler chain) -- agreement should be
        tight (well under the 100 km floor), not a specific golden number."""
        params = _load_branch_params(6)
        v2 = run_v2_saturn_3d("b6", params, mu=MU_SATURN, n_cycles=3)
        v3 = run_v3_saturn_3d("b6", params, mu=MU_SATURN, v2_verdict=v2, n_cycles=3)
        assert v3.n_cycles_propagated == 3
        assert v3.drift_agreement_kms < 1.0  # km -- near machine precision in practice


@pytest.mark.skipif(not PROBE_574_PATH.exists(), reason="#574 Stage-A probe jsonl not present")
class TestV4SaturnRunsEndToEnd:
    def test_v4_returns_well_formed_verdict(self) -> None:
        params = _load_branch_params(6)
        v2 = run_v2_saturn_3d("b6", params, mu=MU_SATURN, n_cycles=3)
        v3 = run_v3_saturn_3d("b6", params, mu=MU_SATURN, v2_verdict=v2, n_cycles=3)
        v4 = run_v4_saturn("b6", params, mu_primary=MU_SATURN, v3_verdict=v3, n_cycles=3)
        assert v4.n_cycles_propagated <= 3
        assert len(v4.per_cycle_drift_kms_v4) == v4.n_cycles_propagated
        assert "Saturn J2" in v4.integrator

    def test_saturn_constants_are_physically_sane(self) -> None:
        # Sourced-value sanity band, not a golden equality (see module docstring).
        assert 0.01 < SATURN_J2 < 0.02
        assert 55_000.0 < SATURN_R_EQ_KM < 65_000.0


@pytest.mark.skipif(SAT441_PATH is None, reason="sat441.bsp SPICE kernel not installed")
@pytest.mark.skipif(not PROBE_574_PATH.exists(), reason="#574 Stage-A probe jsonl not present")
class TestV4SaturnStrictSpice:
    def _chain_to_v4(self, branch_id: int) -> tuple[TitanIapetusClosureParams, Any, Any, Any]:
        from cyclerfinder.data.validation.v4_saturn_strict import run_v4_saturn_strict

        params = _load_branch_params(branch_id)
        v2 = run_v2_saturn_3d(f"b{branch_id}", params, mu=MU_SATURN, n_cycles=3)
        v3 = run_v3_saturn_3d(f"b{branch_id}", params, mu=MU_SATURN, v2_verdict=v2, n_cycles=3)
        v4 = run_v4_saturn(f"b{branch_id}", params, mu_primary=MU_SATURN, v3_verdict=v3, n_cycles=3)
        return params, v3, v4, run_v4_saturn_strict

    def test_spice_sampled_eccentricity_inclination_are_physically_sane(self) -> None:
        params, v3, v4, run_v4_saturn_strict = self._chain_to_v4(6)
        v4s = run_v4_saturn_strict(
            "b6",
            params,
            "2000-06-21T00:00:00",
            mu_primary=MU_SATURN,
            v3_verdict=v3,
            v4_scipy_verdict=v4,
            n_cycles=3,
        )
        # Physically-sane band around the #574 Stage-A spec's own sourced JPL SSD
        # mean-element values (ECC_TITAN/ECC_IAPETUS) -- SPICE osculating elements at
        # one epoch will differ from the mean elements, so this is a wide sanity band,
        # not an equality.
        assert 0.5 * ECC_TITAN < v4s.eccentricity_used_e_titan < 2.0 * ECC_TITAN
        assert 0.5 * ECC_IAPETUS < v4s.eccentricity_used_e_iapetus < 2.0 * ECC_IAPETUS
        assert 5.0 < v4s.inclination_used_deg_iapetus < 25.0
        assert v4s.launch_epoch_utc == "2000-06-21T00:00:00"

    def test_planet_crossing_tag_fires_on_known_trigger_epoch(self) -> None:
        """Branch 2 at 2000-09-15 is a directly-instrumented real
        FAILURE_MODE_PLANET_CROSSING trigger (confirmed by
        ``scripts/run_574_stageB_saturn_gauntlet.py``'s own epoch-sensitivity spot
        check before this test was written) -- pins the #567-inherited fix's tagging
        behaviour against a genuine case, not a synthetic one."""
        from cyclerfinder.data.validation.v4_saturn_strict import (
            FAILURE_MODE_PLANET_CROSSING,
            run_v4_saturn_strict,
        )

        params = _load_branch_params(2)
        v2 = run_v2_saturn_3d("b2", params, mu=MU_SATURN, n_cycles=3)
        v3 = run_v3_saturn_3d("b2", params, mu=MU_SATURN, v2_verdict=v2, n_cycles=3)
        v4 = run_v4_saturn("b2", params, mu_primary=MU_SATURN, v3_verdict=v3, n_cycles=3)
        v4s = run_v4_saturn_strict(
            "b2",
            params,
            "2000-09-15T00:00:00",
            mu_primary=MU_SATURN,
            v3_verdict=v3,
            v4_scipy_verdict=v4,
            n_cycles=3,
        )
        assert v4s.passes_v4_strict is False
        failure_modes = [c.failure_mode for c in v4s.per_cycle]
        assert FAILURE_MODE_PLANET_CROSSING in failure_modes
        crossing_cycle = next(
            c for c in v4s.per_cycle if c.failure_mode == FAILURE_MODE_PLANET_CROSSING
        )
        assert crossing_cycle.perijove_km is not None
        assert crossing_cycle.perijove_km < SATURN_R_EQ_KM
        # Never silently excluded -- the FAIL is real and finite periapsis is recorded.
        assert math.isfinite(crossing_cycle.perijove_km)
