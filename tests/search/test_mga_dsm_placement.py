"""#307 Task 2: Vasile-Conway DSM placement tests.

Core correctness claim: arc 2 is a Lambert to the arrival body, so the arrival
POSITION closes by construction (the Lambert/propagate consistency residual is ~0
for any eta / vinf). The optimiser then drives the *arrival* V∞ to a target while
minimising the DSM Δv.
"""

from __future__ import annotations

import numpy as np

from cyclerfinder.core.constants import AU_KM, MU_SUN_KM3_S2
from cyclerfinder.search.mga_dsm_placement import evaluate_dsm_leg, optimize_dsm_leg

_DAY_S = 86400.0


def _circular_geometry() -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, float]:
    """Earth-like (1 AU) -> Mars-like (1.52 AU at 100 deg) heliocentric legs with
    circular body velocities (perpendicular, prograde)."""
    r_dep = np.array([AU_KM, 0.0, 0.0])
    ang = np.radians(100.0)
    r_arr = 1.52 * AU_KM * np.array([np.cos(ang), np.sin(ang), 0.0])
    v_dep_mag = float(np.sqrt(MU_SUN_KM3_S2 / np.linalg.norm(r_dep)))
    v_arr_mag = float(np.sqrt(MU_SUN_KM3_S2 / np.linalg.norm(r_arr)))
    v_body_dep = v_dep_mag * np.array([0.0, 1.0, 0.0])  # prograde at +x position
    v_body_arr = v_arr_mag * np.array([-np.sin(ang), np.cos(ang), 0.0])  # prograde
    tof_s = 250.0 * _DAY_S
    return r_dep, v_body_dep, r_arr, v_body_arr, tof_s


def test_dsm_leg_closes_position_by_construction() -> None:
    """For ANY (eta, vinf), the Lambert arc 2 lands on the arrival body."""
    r_dep, v_body_dep, r_arr, v_body_arr, tof_s = _circular_geometry()
    for eta in (0.2, 0.4, 0.6):
        res = evaluate_dsm_leg(
            r_dep,
            v_body_dep,
            r_arr,
            v_body_arr,
            tof_s,
            eta=eta,
            vinf_dep_vec_kms=np.array([1.0, 0.5, 0.3]),
        )
        # Lambert targets r_arr; propagate-forward must agree to integrator precision.
        assert res.arrival_pos_residual_km < 1.0, f"eta={eta} pos_res={res.arrival_pos_residual_km}"
        assert np.isfinite(res.dsm_dv_kms) and res.dsm_dv_kms >= 0.0
        assert np.isfinite(res.arrival_vinf_kms) and res.arrival_vinf_kms >= 0.0


def test_optimize_dsm_leg_recovers_ballistic_target_cheaply() -> None:
    """Targeting the BALLISTIC arrival V∞ is reachable at ~zero DSM (the optimiser
    should find the near-ballistic solution: it can set vinf_dep to the ballistic
    departure and ride the ballistic arc, so DSM Δv ≈ 0)."""
    from cyclerfinder.core.lambert import lambert

    r_dep, v_body_dep, r_arr, v_body_arr, tof_s = _circular_geometry()
    ball = lambert(r_dep, r_arr, tof_s, prograde=True, max_revs=0)[0]
    ball_arr_vinf = float(np.linalg.norm(ball.v2 - v_body_arr))
    ball_dep_vinf = float(np.linalg.norm(ball.v1 - v_body_dep))

    res = optimize_dsm_leg(
        r_dep,
        v_body_dep,
        r_arr,
        v_body_arr,
        tof_s,
        target_arrival_vinf_kms=ball_arr_vinf,
        vinf_dep_max_kms=ball_dep_vinf + 3.0,
        seed=0,
        maxiter=60,
        popsize=15,
    )
    assert res.arrival_pos_residual_km < 1.0
    assert abs(res.arrival_vinf_kms - ball_arr_vinf) < 0.3, (
        f"arrival_vinf={res.arrival_vinf_kms} vs ballistic target {ball_arr_vinf}"
    )
    # Reaching the ballistic arrival V∞ needs essentially no DSM.
    assert res.dsm_dv_kms < 1.0, f"expected near-ballistic DSM, got {res.dsm_dv_kms}"
