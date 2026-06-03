"""Real-ephemeris (DE440) two-arc construction of the S1L1 cycler.

Follow-up to scripts/construct_s1l1_twoarc.py (coplanar), which showed the two-arc
STRUCTURE is right but the circular-coplanar FIDELITY cannot reach the published
Mars V_inf of 3.05 km/s (it floors ~4.1). This script repeats the patched
free-return construction on the real DE440 ephemeris (eccentric, inclined Earth
and Mars) to test whether Mars's real eccentricity delivers the low Mars V_inf,
the same effect that made Aldrin's 146-day leg cheap only on the real ephemeris.

Sourced structure (Russell descriptor 4.991gG, schema v4.1 free_return_arcs):
  arc1 1.4612 yr = 533 d  -> Mars arc: E->M 154 d + M->E 379 d
  arc2 2.8096 yr = 1026 d -> phasing arc: multi-rev E->E
Injection window 2030-2034 (Sanchez Net 2022, SSPE). Published V_inf 5.65 / 3.05.

The leg TOFs are sourced/derived construction inputs; 5.65/3.05 is the independent
published cross-check. Every V_inf below is OUR computation (golden cross-check).
"""

from __future__ import annotations

from datetime import UTC, datetime, timedelta

import numpy as np

from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.core.lambert import (
    LambertConvergenceError,
    LambertGeometryError,
    lambert,
)

DAY_S = 86400.0
J2000 = datetime(2000, 1, 1, 12, 0, 0, tzinfo=UTC)
VINF_E_PUB = 5.65
VINF_M_PUB = 3.05

TOF_EM, TOF_ME, TOF_EE = 154.0, 379.0, 1026.0  # days (sourced/derived)

EPHEM = Ephemeris("astropy")


def _t_sec(dt: datetime) -> float:
    return (dt - J2000).total_seconds()


def _vinf(v_sc, v_pl) -> float:
    return float(np.linalg.norm(np.asarray(v_sc) - np.asarray(v_pl)))


def construct(t0: datetime) -> dict | None:
    """Patched 3-leg real-ephemeris construction launched at t0."""
    t_e0 = _t_sec(t0)
    t_m = _t_sec(t0 + timedelta(days=TOF_EM))
    t_e1 = _t_sec(t0 + timedelta(days=TOF_EM + TOF_ME))
    t_e2 = _t_sec(t0 + timedelta(days=TOF_EM + TOF_ME + TOF_EE))

    r_e0, v_e0 = EPHEM.state("E", t_e0)
    r_m, v_m = EPHEM.state("M", t_m)
    r_e1, v_e1 = EPHEM.state("E", t_e1)
    r_e2, v_e2 = EPHEM.state("E", t_e2)

    try:
        leg_em = lambert(r_e0, r_m, TOF_EM * DAY_S, max_revs=0)
        leg_me = lambert(r_m, r_e1, TOF_ME * DAY_S, max_revs=1)
        leg_ee = lambert(r_e1, r_e2, TOF_EE * DAY_S, max_revs=2)
    except (LambertConvergenceError, LambertGeometryError):
        return None

    s_em = leg_em[0]
    vinf_e0 = _vinf(s_em.v1, v_e0)
    vinf_m_in = _vinf(s_em.v2, v_m)

    best = None
    for s_me in leg_me:
        vinf_m_out = _vinf(s_me.v1, v_m)
        vinf_e1_in = _vinf(s_me.v2, v_e1)
        for s_ee in leg_ee:
            vinf_e1_out = _vinf(s_ee.v1, v_e1)
            vinf_e2_in = _vinf(s_ee.v2, v_e2)
            res = (
                abs(vinf_m_in - vinf_m_out)
                + abs(vinf_e1_in - vinf_e1_out)
                + abs(vinf_e2_in - vinf_e0)
            )
            cand = {
                "t0": t0,
                "res_total": res,
                "vinf_e0": vinf_e0,
                "vinf_m_in": vinf_m_in,
                "vinf_m_out": vinf_m_out,
                "vinf_e1_in": vinf_e1_in,
                "vinf_e1_out": vinf_e1_out,
                "vinf_e2_in": vinf_e2_in,
            }
            if best is None or res < best["res_total"]:
                best = cand
    return best


def main() -> None:
    print("Real-ephemeris (DE440) two-arc S1L1 construction")
    print(f"Legs E->M {TOF_EM}d, M->E {TOF_ME}d, E->E {TOF_EE}d; window 2030-2034")
    print(f"Target V_inf: Earth {VINF_E_PUB}, Mars {VINF_M_PUB} km/s\n")
    start = datetime(2030, 1, 1, tzinfo=UTC)
    best_res = None
    best_mars = None  # epoch giving the lowest Mars V_inf (eccentricity test)
    for day in range(0, int(4 * 365.25), 20):
        c = construct(start + timedelta(days=day))
        if c is None:
            continue
        if best_res is None or c["res_total"] < best_res["res_total"]:
            best_res = c
        mars_v = min(c["vinf_m_in"], c["vinf_m_out"])
        if best_mars is None or mars_v < best_mars[0]:
            best_mars = (mars_v, c)
    if best_res is None:
        print("No feasible construction in window.")
        return
    b = best_res
    print(f"Best V_inf-continuity epoch: {b['t0'].date()} (residual {b['res_total']:.3f} km/s)")
    print(
        f"  Earth V_inf: E0={b['vinf_e0']:.2f} E1_in={b['vinf_e1_in']:.2f} "
        f"E1_out={b['vinf_e1_out']:.2f} E2_in={b['vinf_e2_in']:.2f} (target {VINF_E_PUB})"
    )
    print(f"  Mars  V_inf: in={b['vinf_m_in']:.2f} out={b['vinf_m_out']:.2f} (target {VINF_M_PUB})")
    mv, mc = best_mars
    print(
        f"\nLowest Mars V_inf in window: {mv:.2f} km/s at {mc['t0'].date()} "
        f"(coplanar floored ~4.1; target 3.05)"
    )
    print(
        f"  (that epoch's Earth V_inf E0={mc['vinf_e0']:.2f}, "
        f"Mars in/out={mc['vinf_m_in']:.2f}/{mc['vinf_m_out']:.2f})"
    )


if __name__ == "__main__":
    main()
