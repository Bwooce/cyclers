"""Coplanar two-arc construction of the S1L1 two-synodic-period Earth-Mars cycler.

Sourced structure (Russell 2004 descriptor 4.991gG, schema v4.1 free_return_arcs):
  arc 1 = 1.4612 yr = 533 d  -> the MARS arc: E->M (154 d) + M->E (~379 d)
  arc 2 = 2.8096 yr = 1026 d -> the PHASING arc: a multi-rev E->E loop
Total 4.27 yr = 2 synodic periods. Published V_inf: 5.65 (Earth), 3.05 (Mars).

The single-ellipse construction gave V_inf ~ 4.90/4.98 and never closed on
5.65/3.05. This script builds the cycler as the patched sequence of three Lambert
legs, scans the Earth-Mars phase, and reports V_inf at every encounter plus the
flyby V_inf-continuity residuals (a ballistic cycler needs equal V_inf across
each flyby). Coplanar circular model first (Mars eccentricity, expected to matter
for the low 3.05 Mars value, is the next step).

Golden discipline: leg TOFs are sourced/derived inputs (154 d sourced; 379 d =
arc1 533 d - 154 d; 1026 d = arc2). The published 5.65/3.05 is the independent
cross-check target; every V_inf below is OUR computation.
"""

from __future__ import annotations

import numpy as np

from cyclerfinder.core.constants import AU_KM, MU_SUN_KM3_S2, PLANETS
from cyclerfinder.core.lambert import (
    LambertConvergenceError,
    LambertGeometryError,
    lambert,
)

DAY_S = 86400.0
VINF_E_PUB = 5.65
VINF_M_PUB = 3.05

# Sourced / derived leg structure (days).
TOF_EM = 154.0  # E->M outbound (sourced)
TOF_ME = 379.0  # M->E return  (arc1 533 d - 154 d)
TOF_EE = 1026.0  # E->E phasing (arc2 2.8096 yr)

R_EARTH = PLANETS["E"].sma_au * AU_KM
R_MARS = PLANETS["M"].sma_au * AU_KM
N_EARTH = np.radians(PLANETS["E"].mean_motion_deg_day)  # rad/day
N_MARS = np.radians(PLANETS["M"].mean_motion_deg_day)


def _state(radius_km: float, n_rad_day: float, lon0_rad: float, t_day: float) -> tuple:
    """Coplanar circular position+velocity at time t for a body of given radius."""
    ang = lon0_rad + n_rad_day * t_day
    r = np.array([radius_km * np.cos(ang), radius_km * np.sin(ang), 0.0])
    v_circ = float(np.sqrt(MU_SUN_KM3_S2 / radius_km))
    v = np.array([-v_circ * np.sin(ang), v_circ * np.cos(ang), 0.0])
    return r, v


def _best_leg(r1, r2, tof_day: float, max_revs: int, v_ref_mag: float | None):
    """Lambert leg; pick the branch whose departure V_inf is closest to v_ref_mag
    (or the lowest-energy solution when no reference is given). Returns the
    LambertSolution or None."""
    try:
        sols = lambert(r1, r2, tof_day * DAY_S, max_revs=max_revs)
    except (LambertConvergenceError, LambertGeometryError):
        return None
    return sols


def construct(phi_deg: float) -> dict | None:
    """Build the 3-leg cycler for Earth-Mars phase phi (Mars longitude at t0)."""
    phi = np.radians(phi_deg)
    # Encounter epochs (days from t0).
    t_e0, t_m, t_e1, t_e2 = 0.0, TOF_EM, TOF_EM + TOF_ME, TOF_EM + TOF_ME + TOF_EE

    r_e0, v_e0 = _state(R_EARTH, N_EARTH, 0.0, t_e0)
    r_m, v_m = _state(R_MARS, N_MARS, phi, t_m)
    r_e1, v_e1 = _state(R_EARTH, N_EARTH, 0.0, t_e1)
    r_e2, v_e2 = _state(R_EARTH, N_EARTH, 0.0, t_e2)

    leg_em = _best_leg(r_e0, r_m, TOF_EM, 0, None)
    leg_me = _best_leg(r_m, r_e1, TOF_ME, 1, None)
    leg_ee = _best_leg(r_e1, r_e2, TOF_EE, 2, None)
    if not (leg_em and leg_me and leg_ee):
        return None

    def vinf(v_sc, v_pl):
        return float(np.linalg.norm(np.asarray(v_sc) - np.asarray(v_pl)))

    # Outbound: single-rev.
    s_em = leg_em[0]
    vinf_e0 = vinf(s_em.v1, v_e0)
    vinf_m_in = vinf(s_em.v2, v_m)

    # Return + phasing: choose branch minimising Mars / Earth V_inf continuity.
    best = None
    for s_me in leg_me:
        vinf_m_out = vinf(s_me.v1, v_m)
        vinf_e1_in = vinf(s_me.v2, v_e1)
        for s_ee in leg_ee:
            vinf_e1_out = vinf(s_ee.v1, v_e1)
            vinf_e2_in = vinf(s_ee.v2, v_e2)
            res_m = abs(vinf_m_in - vinf_m_out)
            res_e1 = abs(vinf_e1_in - vinf_e1_out)
            res_cycle = abs(vinf_e2_in - vinf_e0)
            total = res_m + res_e1 + res_cycle
            cand = {
                "phi": phi_deg,
                "vinf_e0": vinf_e0,
                "vinf_m_in": vinf_m_in,
                "vinf_m_out": vinf_m_out,
                "vinf_e1_in": vinf_e1_in,
                "vinf_e1_out": vinf_e1_out,
                "vinf_e2_in": vinf_e2_in,
                "res_m": res_m,
                "res_e1": res_e1,
                "res_cycle": res_cycle,
                "res_total": total,
                "me_branch": (s_me.n_revs, s_me.branch),
                "ee_branch": (s_ee.n_revs, s_ee.branch),
            }
            if best is None or total < best["res_total"]:
                best = cand
    return best


def main() -> None:
    print("Coplanar two-arc S1L1 construction (Russell 4.991gG free-return arcs)")
    print(f"Legs: E->M {TOF_EM}d, M->E {TOF_ME}d, E->E {TOF_EE}d (multi-rev)")
    print(f"Target V_inf: Earth {VINF_E_PUB}, Mars {VINF_M_PUB} km/s\n")
    best = None
    for phi in np.arange(0.0, 360.0, 1.0):
        c = construct(phi)
        if c is None:
            continue
        if best is None or c["res_total"] < best["res_total"]:
            best = c
    if best is None:
        print("No feasible construction found.")
        return
    b = best
    print(
        f"Best-closure phase phi = {b['phi']:.0f} deg "
        f"(continuity residual {b['res_total']:.3f} km/s)"
    )
    print(f"  M->E branch {b['me_branch']}, E->E branch {b['ee_branch']}")
    print(
        "  Earth V_inf at encounters: "
        f"E0={b['vinf_e0']:.3f}  E1_in={b['vinf_e1_in']:.3f}  "
        f"E1_out={b['vinf_e1_out']:.3f}  E2_in={b['vinf_e2_in']:.3f}  "
        f"(target {VINF_E_PUB})"
    )
    print(
        f"  Mars  V_inf: in={b['vinf_m_in']:.3f}  out={b['vinf_m_out']:.3f}  (target {VINF_M_PUB})"
    )
    print(f"  residuals: Mars={b['res_m']:.3f}  E1={b['res_e1']:.3f}  cycle={b['res_cycle']:.3f}")
    earth_mean = np.mean([b["vinf_e0"], b["vinf_e1_in"], b["vinf_e1_out"], b["vinf_e2_in"]])
    mars_mean = np.mean([b["vinf_m_in"], b["vinf_m_out"]])
    print(f"\n  mean Earth V_inf {earth_mean:.2f} (target 5.65; single-ellipse gave 4.90)")
    print(f"  mean Mars  V_inf {mars_mean:.2f} (target 3.05; single-ellipse gave 4.98)")


if __name__ == "__main__":
    main()
