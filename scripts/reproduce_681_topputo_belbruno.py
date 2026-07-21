#!/usr/bin/env python3
"""Positive-control reproduction for task #681 (Topputo & Belbruno 2015).

Reproduces the checkable numbers from Topputo & Belbruno, "Earth-Mars transfers
with ballistic capture," CeMDA 121:329 (2015), doi:10.1007/s10569-015-9605-8,
that validate the Sun-Mars physical constants, heliocentric geometry, and the
elliptic-restricted-three-body ballistic-capture machinery in
`cyclerfinder.core.sunmars_wsb` BEFORE any repeating-capture search is trusted:

  (a) Table 5 -- four bitangential Hohmann baselines (analytic patched conics).
  (b) Table 3 -- the H3-equivalent Hohmann orbit-insertion cost dV2 vs r_p
      (analytic; grows with r_p).
  (c) Table 3 -- the ballistic-capture cost dV_c vs r_p (the 3-body machinery;
      the paper's headline "flat ~2.03-2.04 km/s" signature).

Run: uv run python scripts/reproduce_681_topputo_belbruno.py
"""

from __future__ import annotations

import math

import numpy as np
from scipy.integrate import solve_ivp

import cyclerfinder.core.sunmars_wsb as sw
from cyclerfinder.core.lambert import LambertError, lambert

F0 = math.pi / 4.0  # paper's first-quadrant f0 for Table 3
ECC = 0.99

# Paper targets (digest / Table 5, Table 3).
_TABLE5 = {
    ("peri", "peri"): (2.179, 3.388, 234),
    ("peri", "apo"): (3.398, 2.090, 278),
    ("apo", "peri"): (2.414, 3.163, 239),
    ("apo", "apo"): (3.629, 1.881, 283),
}
_TABLE3_RP = [49896.0, 73896.0, 91897.0, 113897.0]
_TABLE3_DV2 = [2.116, 2.267, 2.344, 2.414]  # H3-equivalent Hohmann insertion
_TABLE3_DVC = [2.033, 2.036, 2.039, 2.041]  # ballistic-capture cost
_VINF_H3 = 3.163  # Table 5 H3 arrival v-infinity


def reproduce_hohmann() -> None:
    print("== (a) Table 5 bitangential Hohmann baselines ==")
    print(f"{'case':<11}{'dV1':>18}{'dV2,inf':>18}{'tof(d)':>16}")
    for (ea, ma), (p1, p2, ptof) in _TABLE5.items():
        h = sw.hohmann_baseline(earth_apsis=ea, mars_apsis=ma)  # type: ignore[arg-type]
        print(
            f"{ea + '->' + ma:<11}"
            f"{h.dv1_kms:8.3f} (pap {p1:<5})"
            f"{h.dv2_inf_kms:8.3f} (pap {p2:<5})"
            f"{h.tof_days:7.1f} (pap {ptof})"
        )


def hohmann_insertion_dv2(rp_km: float, vinf: float = _VINF_H3, ecc: float = ECC) -> float:
    """H3-equivalent Hohmann orbit-insertion cost into a (r_p, e) orbit (Table 3)."""
    v_hyp = math.sqrt(vinf**2 + 2.0 * sw.MU_MARS_KM3_S2 / rp_km)
    v_tar = math.sqrt(sw.MU_MARS_KM3_S2 * (1.0 + ecc) / rp_km)
    return v_hyp - v_tar


def reproduce_dv2() -> None:
    print("\n== (b) Table 3 H3-equivalent Hohmann insertion dV2 vs r_p ==")
    for rp, pap in zip(_TABLE3_RP, _TABLE3_DV2, strict=True):
        print(f"  r_p={rp:7.0f} km  dV2={hohmann_insertion_dv2(rp):.3f}  (paper {pap})")


def _find_xc(
    theta: float, branch: str, rp: float, target_km: float = 1.0e6, max_days: float = 600.0
) -> tuple[float, np.ndarray] | None:
    """Backward-integrate a capture periapsis to its post-escape x_c (~target_km)."""
    st = sw.capture_periapsis_state(r_p_km=rp, ecc=ECC, theta=theta, f0=F0, branch=branch)  # type: ignore[arg-type]

    def rhs(tau: float, y: np.ndarray) -> np.ndarray:
        return -sw.sunmars_eom(-tau, y, F0)

    sol = solve_ivp(
        rhs,
        (0.0, max_days * 86400.0),
        st,
        method="DOP853",
        rtol=1e-11,
        atol=1e-9,
        max_step=0.01 * sw.MARS_PERIOD_S,
        dense_output=True,
    )
    ts = np.linspace(0.0, sol.t[-1], 6000)
    y_hist = sol.sol(ts)
    e2 = np.array([sw.mars_kepler_energy(y_hist[:, k], -ts[k], F0) for k in range(len(ts))])
    dist = np.array([sw.mars_distance(y_hist[:, k], -ts[k], F0) for k in range(len(ts))])
    esc = np.where(e2 > 0.0)[0]
    if not len(esc):
        return None
    after = np.where(dist[esc[0] :] >= target_km)[0]
    if not len(after):
        return None
    k = esc[0] + after[0]
    return -ts[k], y_hist[:, k]


def dvc_min(rp: float, n_theta: int = 16) -> float:
    """Minimum ballistic-capture cost dV_c over a coarse capture/transfer grid.

    dV_c = min |v_capture(x_c) - v_arrival(x_c)| over capture family (theta,
    branch) and Earth->x_c Lambert transfers (departure phase, TOF).
    """
    best = math.inf
    for branch in ("prograde", "retrograde"):
        for theta in np.linspace(0.0, 2.0 * math.pi, n_theta, endpoint=False):
            res = _find_xc(theta, branch, rp)
            if res is None:
                continue
            t_c, y = res
            r_c = np.array([y[0], y[1], 0.0])
            v_cap = np.array([y[2], y[3], 0.0])
            for eph in np.linspace(0.0, 2.0 * math.pi, 48, endpoint=False):
                for tof_d in np.linspace(150.0, 320.0, 28):
                    re, _ = sw.body_state(
                        t_c - tof_d * 86400.0, a_km=sw.EARTH_A_KM, ecc=sw.EARTH_E, f0=eph
                    )
                    r1 = np.array([re[0], re[1], 0.0])
                    try:
                        sols = lambert(r1, r_c, tof_d * 86400.0, mu=sw.MU_SUN_KM3_S2, prograde=True)
                    except LambertError:
                        continue
                    for s in sols:
                        dv = float(np.linalg.norm(v_cap - np.array(s.v2)))
                        best = min(best, dv)
    return best


def reproduce_dvc() -> None:
    print("\n== (c) Table 3 ballistic-capture cost dV_c vs r_p (3-body machinery) ==")
    print("     (min over a COARSE capture-family x transfer grid; the paper's")
    print("      exact 2.03-2.04 needs its full optimized capture-set pipeline)")
    vals = []
    for rp, pap in zip(_TABLE3_RP, _TABLE3_DVC, strict=True):
        d = dvc_min(rp)
        vals.append(d)
        print(f"  r_p={rp:7.0f} km  dV_c(min)={d:.3f}  (paper {pap})")
    print(
        f"  --> range [{min(vals):.3f}, {max(vals):.3f}] km/s, roughly FLAT vs r_p "
        f"(paper 2.033-2.041); contrast dV2 grows 2.12->2.41."
    )


def main() -> None:
    reproduce_hohmann()
    reproduce_dv2()
    reproduce_dvc()


if __name__ == "__main__":
    main()
