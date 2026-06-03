"""Two-arc differential corrector for the S1L1 cycler (real ephemeris).

The constrained optimiser (close_s1l1_realeph.py) collapsed into a degenerate
basin and never reached the Mars 2.83 km/s that the direct construction
(construct_s1l1_twoarc_realeph.py) finds. This builds the corrector directly on
that working construction: Newton / least-squares on the leg ToFs + launch epoch
to drive the three flyby V_inf-continuity residuals to zero (ballistic closure),
seeded from the direct two-arc solution.

Variables: x = [t0_offset_days, T_EM, T_ME]; the phasing leg T_EE is pinned by the
fixed 2-synodic period (arcs 1.4612 + 2.8096 yr = 4.2708 yr). Three residuals:
  r1 = |Vinf_M_in| - |Vinf_M_out|        (Mars flyby conserves energy)
  r2 = |Vinf_E1_in| - |Vinf_E1_out|      (intermediate Earth flyby)
  r3 = |Vinf_E2_in| - |Vinf_E0|          (cycle closure at Earth)
Square 3x3 system. A flyby only ROTATES V_inf, so magnitude continuity is the
closure condition; the required turn angle is checked post-hoc for feasibility.

Sourced: 5.65/3.05 anchors (spec §9) are the independent cross-check (NOT fitted);
leg-ToF seed + period from Russell 4.991gG arcs. Every V_inf is OUR computation.
"""

from __future__ import annotations

from datetime import UTC, datetime, timedelta

import numpy as np
from scipy.optimize import least_squares

from cyclerfinder.core.constants import PLANETS
from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.core.lambert import (
    LambertConvergenceError,
    LambertGeometryError,
    lambert,
)

DAY_S = 86400.0
J2000 = datetime(2000, 1, 1, 12, 0, 0, tzinfo=UTC)
VINF_E_PUB, VINF_M_PUB = 5.65, 3.05
PERIOD_DAYS = (1.4612 + 2.8096) * 365.25  # 2-synodic, sourced arcs
SEED = {"t0": datetime(2030, 3, 22, tzinfo=UTC), "t_em": 154.0, "t_me": 379.0}


def _t_sec(dt: datetime) -> float:
    return (dt - J2000).total_seconds()


def _pick(sols: list, n_revs: int, branch: str):
    for s in sols:
        if s.n_revs == n_revs and s.branch == branch:
            return s
    return sols[0]


def _legs(t0_day: float, t_em: float, t_me: float, t_ee: float, ephem: Ephemeris):
    t0 = SEED["t0"] + timedelta(days=t0_day)
    te0 = _t_sec(t0)
    tm = _t_sec(t0 + timedelta(days=t_em))
    te1 = _t_sec(t0 + timedelta(days=t_em + t_me))
    te2 = _t_sec(t0 + timedelta(days=t_em + t_me + t_ee))
    r_e0, v_e0 = ephem.state("E", te0)
    r_m, v_m = ephem.state("M", tm)
    r_e1, v_e1 = ephem.state("E", te1)
    r_e2, v_e2 = ephem.state("E", te2)
    em = lambert(r_e0, r_m, t_em * DAY_S, max_revs=0)
    me = lambert(r_m, r_e1, t_me * DAY_S, max_revs=1)
    ee = lambert(r_e1, r_e2, t_ee * DAY_S, max_revs=2)
    return (v_e0, v_m, v_e1, v_e2), em, me, ee


# Branch labels fixed from the seed (set in main()).
_ME_BR = {"n": 0, "b": "single"}
_EE_BR = {"n": 2, "b": "low"}


def _state_vinf(t0_day, t_em, t_me, ephem):
    t_ee = PERIOD_DAYS - t_em - t_me
    (v_e0, v_m, v_e1, v_e2), em, me, ee = _legs(t0_day, t_em, t_me, t_ee, ephem)
    s_em = em[0]
    s_me = _pick(me, _ME_BR["n"], _ME_BR["b"])
    s_ee = _pick(ee, _EE_BR["n"], _EE_BR["b"])

    def vi(vsc, vpl):
        return np.asarray(vsc) - np.asarray(vpl)

    return {
        "m_in": vi(s_em.v2, v_m),
        "m_out": vi(s_me.v1, v_m),
        "e0": vi(s_em.v1, v_e0),
        "e1_in": vi(s_me.v2, v_e1),
        "e1_out": vi(s_ee.v1, v_e1),
        "e2_in": vi(s_ee.v2, v_e2),
    }


def _residuals(x, ephem):
    try:
        d = _state_vinf(x[0], x[1], x[2], ephem)
    except (LambertConvergenceError, LambertGeometryError, ValueError):
        return [1e3, 1e3, 1e3]
    n = np.linalg.norm
    return [
        n(d["m_in"]) - n(d["m_out"]),
        n(d["e1_in"]) - n(d["e1_out"]),
        n(d["e2_in"]) - n(d["e0"]),
    ]


def _max_bend_deg(vinf_kms: float, body: str) -> float:
    pl = PLANETS[body]
    mu = pl.mu_km3_s2
    r_p = pl.radius_eq_km + pl.safe_alt_km
    e = 1.0 + r_p * vinf_kms * vinf_kms / mu
    return float(np.degrees(2.0 * np.arcsin(1.0 / e)))


def _bend_deg(v_in, v_out) -> float:
    a, b = np.asarray(v_in), np.asarray(v_out)
    c = float(np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b)))
    return float(np.degrees(np.arccos(max(-1.0, min(1.0, c)))))


def _solve(t0_off: float, ephem: Ephemeris) -> dict | None:
    """Run the corrector from one epoch offset; return the closed solution data."""
    sol = least_squares(
        _residuals,
        [t0_off, SEED["t_em"], SEED["t_me"]],
        args=(ephem,),
        method="lm",
        max_nfev=80,
        xtol=1e-9,
        ftol=1e-9,
    )
    x = sol.x
    res = _residuals(x, ephem)
    max_res = max(abs(r) for r in res)
    d = _state_vinf(x[0], x[1], x[2], ephem)
    n = np.linalg.norm
    ve = float(np.mean([n(d["e0"]), n(d["e1_in"]), n(d["e1_out"]), n(d["e2_in"])]))
    vm = float(np.mean([n(d["m_in"]), n(d["m_out"])]))
    bend_m, max_m = _bend_deg(d["m_in"], d["m_out"]), _max_bend_deg(n(d["m_in"]), "M")
    bend_e = _bend_deg(d["e1_in"], d["e1_out"]), _max_bend_deg(n(d["e1_in"]), "E")
    return {
        "t0": SEED["t0"] + timedelta(days=float(x[0])),
        "t_em": float(x[1]),
        "t_me": float(x[2]),
        "t_ee": float(PERIOD_DAYS - x[1] - x[2]),
        "max_res": max_res,
        "vinf_e": ve,
        "vinf_m": vm,
        "feas": bend_m <= max_m and bend_e[0] <= bend_e[1],
    }


def main() -> None:
    ephem = Ephemeris("astropy")
    print("Two-arc differential corrector for S1L1 (real DE440) — multi-start")
    print(f"Period {PERIOD_DAYS:.1f} d (2-synodic). Objective: flyby V_inf continuity ONLY.")
    print(f"Anchors E={VINF_E_PUB} M={VINF_M_PUB} km/s are the CROSS-CHECK (not fitted).\n")
    closed: list[dict] = []
    # Scan launch epochs across ~3 synodic windows, over the converging EE branches.
    for branch in [(2, "low"), (1, "low")]:
        _EE_BR["n"], _EE_BR["b"] = branch
        for off in range(-180, 1100, 80):
            try:
                s = _solve(float(off), ephem)
            except (LambertConvergenceError, LambertGeometryError, ValueError):
                continue
            if s is None or s["max_res"] >= 0.1:
                continue
            s["branch"] = branch
            closed.append(s)
            tag = (
                "  <-- near anchors"
                if (abs(s["vinf_e"] - VINF_E_PUB) < 0.5 and abs(s["vinf_m"] - VINF_M_PUB) < 0.5)
                else ""
            )
            print(
                f"  closed: EE{branch} t0={s['t0'].date()} "
                f"E={s['vinf_e']:.2f} M={s['vinf_m']:.2f} feas={s['feas']}{tag}"
            )
    if not closed:
        print("No ballistic closure found across the scan.")
        return
    nearest = min(
        closed,
        key=lambda s: abs(s["vinf_e"] - VINF_E_PUB) + abs(s["vinf_m"] - VINF_M_PUB),
    )
    print(f"\n{len(closed)} closed ballistic cyclers found.")
    print(f"Distinct Earth V_inf families: {sorted({round(s['vinf_e'], 1) for s in closed})}")
    print(f"Distinct Mars  V_inf families: {sorted({round(s['vinf_m'], 1) for s in closed})}")
    print(
        f"\nNEAREST to published 5.65/3.05: EE{nearest['branch']} t0={nearest['t0'].date()} "
        f"T_EM={nearest['t_em']:.0f} T_ME={nearest['t_me']:.0f} T_EE={nearest['t_ee']:.0f} d"
    )
    print(
        f"  V_inf E={nearest['vinf_e']:.2f} (t {VINF_E_PUB})  "
        f"M={nearest['vinf_m']:.2f} (t {VINF_M_PUB})  feasible_flybys={nearest['feas']}"
    )
    within = abs(nearest["vinf_e"] - VINF_E_PUB) < 0.5 and abs(nearest["vinf_m"] - VINF_M_PUB) < 0.5
    print(f"Reproduces S1L1 anchors within 0.5 km/s: {within}")


if __name__ == "__main__":
    main()
