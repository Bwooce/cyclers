"""Two-arc diagnostic for the S1L1 two-synodic-period Earth-Mars cycler.

Hypothesis (2026-06-03 research): S1L1 is NOT a single repeating ellipse but a
patched sequence of Lambert "generic-return" arcs (McConaghy 2006; Russell 2004
descriptor 4.991gG). The single-(a,e) construction gave V_inf ~ 4.90/4.98 km/s
and never the published 5.65/3.05. A generic 154-day E->M transfer has a free
heliocentric transfer angle theta; its V_inf falls out of the Lambert solution.

This script scans theta for the coplanar-circular 154-day Earth->Mars Lambert arc
and reports V_inf at Earth (departure) and Mars (arrival), to see whether a single
theta reproduces the published 5.65/3.05 anchors. Pure diagnostic; prints only.

Sourced inputs: 154-day E->M transit (spec.md section 9 / McConaghy 2006); V_inf
anchors 5.65/3.05 (spec.md section 9). The Lambert output is OUR computation
cross-checked against those independently-published anchors (golden cross-check,
not circular).
"""

from __future__ import annotations

import numpy as np

from cyclerfinder.core.constants import AU_KM, MU_SUN_KM3_S2, PLANETS
from cyclerfinder.core.lambert import LambertConvergenceError, LambertGeometryError, lambert

DAY_S = 86400.0
TOF_DAYS = 154.0
VINF_E_PUB = 5.65
VINF_M_PUB = 3.05

R_EARTH = PLANETS["E"].sma_au * AU_KM
R_MARS = PLANETS["M"].sma_au * AU_KM
V_EARTH = float(np.sqrt(MU_SUN_KM3_S2 / R_EARTH))  # Earth circular speed
V_MARS = float(np.sqrt(MU_SUN_KM3_S2 / R_MARS))  # Mars circular speed


def vinf_pair(theta_deg: float, tof_days: float) -> tuple[float, float] | None:
    """V_inf (Earth, Mars) for a coplanar 154-d E->M Lambert at transfer angle theta."""
    th = np.radians(theta_deg)
    r1 = np.array([R_EARTH, 0.0, 0.0])
    r2 = np.array([R_MARS * np.cos(th), R_MARS * np.sin(th), 0.0])
    vel_earth = np.array([0.0, V_EARTH, 0.0])
    vel_mars = np.array([-V_MARS * np.sin(th), V_MARS * np.cos(th), 0.0])
    try:
        sols = lambert(r1, r2, tof_days * DAY_S, max_revs=0)
    except (LambertConvergenceError, LambertGeometryError):
        return None
    sol = sols[0]
    vinf_e = float(np.linalg.norm(np.asarray(sol.v1) - vel_earth))
    vinf_m = float(np.linalg.norm(np.asarray(sol.v2) - vel_mars))
    return vinf_e, vinf_m


def main() -> None:
    print(
        f"Coplanar circular model: r_E={R_EARTH / AU_KM:.4f} AU (v={V_EARTH:.3f}), "
        f"r_M={R_MARS / AU_KM:.4f} AU (v={V_MARS:.3f})"
    )
    print(
        f"Published S1L1 anchors: V_inf,E={VINF_E_PUB}, V_inf,M={VINF_M_PUB} km/s; "
        f"E->M TOF={TOF_DAYS} d\n"
    )
    print(f"{'theta':>7} {'Vinf_E':>8} {'Vinf_M':>8} {'|dE|':>7} {'|dM|':>7} {'sum':>7}")
    best = None
    for theta in np.arange(60.0, 300.5, 2.0):
        pair = vinf_pair(theta, TOF_DAYS)
        if pair is None:
            continue
        ve, vm = pair
        d_e, d_m = abs(ve - VINF_E_PUB), abs(vm - VINF_M_PUB)
        s = d_e + d_m
        if best is None or s < best[0]:
            best = (s, theta, ve, vm)
        if theta % 10 < 2:
            print(f"{theta:7.0f} {ve:8.3f} {vm:8.3f} {d_e:7.3f} {d_m:7.3f} {s:7.3f}")
    if best is not None:
        s, theta, ve, vm = best
        print(
            f"\nBEST coplanar single-arc match to 5.65/3.05: theta={theta:.0f} deg -> "
            f"V_inf,E={ve:.3f}, V_inf,M={vm:.3f} (|sum-err|={s:.3f} km/s)"
        )
        m_best = min(
            (p for th in np.arange(60.0, 300.5, 1.0) if (p := vinf_pair(th, TOF_DAYS)) is not None),
            key=lambda p: abs(p[1] - VINF_M_PUB),
        )
        print(
            f"theta best-matching Mars=3.05 alone -> V_inf,E={m_best[0]:.3f}, "
            f"V_inf,M={m_best[1]:.3f}"
        )


if __name__ == "__main__":
    main()
