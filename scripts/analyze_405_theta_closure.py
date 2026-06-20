"""#411 theta-closure FEASIBILITY map for the patched-CR3BP cross-system cycle.

Necessary-condition check (NOT a closed-orbit claim): can the ~2.2 rad single-revolution
theta-time-consistency gap (quantified in the 2026-06-20 #411 analysis, see
data/negative_results.yaml cross_system_se_em_L2_patched_cr3bp) be nulled (mod 2pi) by
adding integer revolutions on the EM-L2 and/or SE-L2 orbits before manifold departure?

Each extra revolution on an orbit adds T_orbit to the cycle, advancing the SE-EM relative
phase by omega_rel * T_orbit. We grid-search (n_em, n_se) for the minimal mod-2pi residual
and report the smallest feasible revolution counts + the resulting cycle duration. This
GATES whether the full multi-rev corrector is worth building; it asserts no orbit.

Run: uv run python scripts/analyze_405_theta_closure.py
"""

from __future__ import annotations

import math

from cyclerfinder.genome.cross_system_cycle import (
    em_moon_system,
    se_earth_system,
    theta_commensurability,
)
from cyclerfinder.genome.heteroclinic_cycle import LyapunovNode

# Single-revolution theta-time-consistency gap, from the committed #411 analysis
# (forward 52-day transit advances theta 4.88 rad; theta_f+omega_rel*t_fwd vs theta_r).
SINGLE_REV_GAP_RAD = 2.2
CANALIAS_C_SE = 3.000863625


def main() -> None:
    se = se_earth_system()
    em = em_moon_system()
    # omega_rel = relative SE-EM line rate (rad/s); CR3BP mean motion = 1 per nondim time
    # => omega = 1/t_s (rad/s).
    omega_rel = 1.0 / em.t_s - 1.0 / se.t_s  # rad/s (EM faster than SE)

    # Orbit periods (nondim -> seconds) for the two cycle nodes (cheap to build).
    em_l2 = LyapunovNode.from_libration(
        em, x0_guess=1.18, jacobi=3.15, period_guess=3.4, label="EM-L2", ydot0_sign=-1.0
    )
    se_l2 = LyapunovNode.from_libration(
        se, x0_guess=1.009, jacobi=CANALIAS_C_SE, period_guess=3.06, label="SE-L2", ydot0_sign=-1.0
    )
    t_em_s = em_l2.period * em.t_s
    t_se_s = se_l2.period * se.t_s
    dtheta_em = omega_rel * t_em_s  # rad of relative phase per extra EM-L2 revolution
    dtheta_se = omega_rel * t_se_s  # rad per extra SE-L2 revolution
    yr = 365.25 * 86400.0

    two_pi = 2 * math.pi
    de, ds = dtheta_em % two_pi, dtheta_se % two_pi
    print("=== #411 theta-closure feasibility (patched-CR3BP cross-system cycle) ===")
    print(f"omega_rel      = {omega_rel:.6e} rad/s")
    print(f"EM-L2 period   = {t_em_s / 86400:.2f} d  -> dtheta_em = {de:.4f} rad/rev")
    print(f"SE-L2 period   = {t_se_s / yr:.3f} yr -> dtheta_se = {ds:.4f} rad/rev")
    print(f"single-rev gap = {SINGLE_REV_GAP_RAD:.3f} rad (from committed #411 analysis)")

    for tol in (5e-2, 1e-2):
        n_em, n_se, res, feasible = theta_commensurability(
            SINGLE_REV_GAP_RAD, dtheta_em, dtheta_se, n_max=60, tol_rad=tol
        )
        period_yr = (n_em * t_em_s + n_se * t_se_s) / yr
        verdict = "FEASIBLE" if feasible else "not feasible"
        print(
            f"tol={tol:.0e} rad -> best (n_em={n_em}, n_se={n_se}) residual={res:.4f} rad "
            f"[{verdict}], cycle duration ~ {period_yr:.2f} yr"
        )


if __name__ == "__main__":
    main()
