"""#968 amendment 6, sec. 5.2: GanCal#5 in a patched conic WITH the unscheduled Ganymede pass.

Kepler legs about Jupiter in the R-S circular model (Callisto massless); the two scheduled Ganymede
flybys stay zero-SOI as in R-S; on the B->C leg the distant Ganymede pass is an impulsive turn of
the Ganymede-relative velocity at the Kepler closest approach, delta = 2 atan(mu / (d v^2)), toward
Ganymede. Unknowns: the three leg start dates and the B->C departure velocity. Residuals: Callisto
hit at its date, V_inf magnitude match at B and A, vector match at Callisto. Writes
data/968_control/kick.json.

    uv run python scripts/run_968_kick.py
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import numpy as np
from numpy.typing import NDArray
from scipy.optimize import least_squares, minimize_scalar

from cyclerfinder.core.lambert import lambert
from cyclerfinder.data.method_capability import MethodCapability
from cyclerfinder.data.preflight import preflight_search
from cyclerfinder.search.two_working_body import cycle_flybys, kepler_step

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "968_control"
DAY = 86400.0
KEY = "k3|LGanymede>Ganymede/1l|LGanymede>Callisto/1h|LCallisto>Ganymede/0s"
R_GAN_RS_KM = 2634.0
Arr = NDArray[np.float64]


def _enum() -> Any:
    import importlib.util

    spec = importlib.util.spec_from_file_location(
        "run_942_enumerate", ROOT / "scripts" / "run_942_enumerate.py"
    )
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def lambert_leg(
    circ: Any, frm: str, to: str, t0: float, t1: float, nrev: int, branch: str
) -> tuple[Arr, Arr]:
    """Departure and arrival V_inf vectors of one Lambert leg."""
    r1, w1 = circ.state(frm, t0)
    r2, w2 = circ.state(to, t1)
    sols = lambert(r1, r2, t1 - t0, mu=circ.mu, max_revs=nrev)
    sol = next(s for s in sols if s.n_revs == nrev and s.branch == branch)
    return sol.v1 - w1, sol.v2 - w2


def kicked_leg(
    circ: Any, u: Arr, t0: float, t1: float, kick: bool
) -> tuple[Arr, Arr, dict[str, float]]:
    """B->C leg from Ganymede at t0 with inertial velocity u; arrival state at t1."""
    mu_g = circ.body("Ganymede").mu_km3_s2
    r0, _ = circ.state("Ganymede", t0)

    def dist(t: float) -> float:
        r, _ = kepler_step(r0, u, t - t0, circ.mu)
        return float(np.linalg.norm(r - circ.state("Ganymede", t)[0]))

    ts = np.arange(t0 + 2 * DAY, t1 - 2 * DAY, 0.02 * DAY)
    ds = np.array([dist(float(t)) for t in ts])
    i = int(np.argmin(ds))
    res = minimize_scalar(
        dist,
        bounds=(ts[max(i - 1, 0)], ts[min(i + 1, len(ts) - 1)]),
        method="bounded",
        options={"xatol": 1.0},
    )
    t_ca = float(res.x)
    r, v = kepler_step(r0, u, t_ca - t0, circ.mu)
    rg, vg = circ.state("Ganymede", t_ca)
    rr, vr = r - rg, v - vg
    d = float(np.linalg.norm(rr))
    vmag = float(np.linalg.norm(vr))
    delta = 2.0 * math.atan(mu_g / (d * vmag * vmag)) if kick else 0.0
    n = np.cross(rr, vr)
    n = n / float(np.linalg.norm(n))
    vr_new = (
        vr * math.cos(delta)
        + np.cross(n, vr) * math.sin(delta)
        + n * float(n @ vr) * (1 - math.cos(delta))
    )
    r1, v1 = kepler_step(r, vg + vr_new, t1 - t_ca, circ.mu)
    return (
        r1,
        v1,
        {"t_ca_days": t_ca / DAY, "d_km": d, "v_rel_kms": vmag, "delta_deg": math.degrees(delta)},
    )


def residual(circ: Any, y: Arr, period: float, kick: bool) -> tuple[Arr, dict[str, Any]]:
    ta, tb, tc = y[0] * DAY, y[1] * DAY, y[2] * DAY
    u = y[3:6]
    _, vin_b = lambert_leg(circ, "Ganymede", "Ganymede", ta, tb, 1, "low")  # A -> B
    vout_a_next, _ = lambert_leg(circ, "Ganymede", "Ganymede", ta, tb, 1, "low")
    vdep_c, varr_a = lambert_leg(
        circ, "Callisto", "Ganymede", tc, ta + period, 0, "single"
    )  # C -> A'
    rc, vc_body = circ.state("Callisto", tc)
    r1, v1, ca = kicked_leg(circ, u, tb, tc, kick)
    _, vg_b = circ.state("Ganymede", tb)
    vout_b = u - vg_b
    varr_c = v1 - vc_body
    # V_inf at A' (= A one period later) compared with A's departure, rotated by Ganymede's advance.
    res = np.concatenate(
        [
            (r1 - rc) / 1e3,  # Callisto hit (1000 km units)
            [float(np.linalg.norm(vin_b) - np.linalg.norm(vout_b))],
            varr_c - vdep_c,
            [float(np.linalg.norm(varr_a) - np.linalg.norm(vout_a_next))],
        ]
    )
    info = {
        "vinf_B_in": float(np.linalg.norm(vin_b)),
        "vinf_B_out": float(np.linalg.norm(vout_b)),
        "vinf_A_in": float(np.linalg.norm(varr_a)),
        "vinf_A_out": float(np.linalg.norm(vout_a_next)),
        "callisto_speed": float(np.linalg.norm(varr_c)),
        "turn_B_deg": math.degrees(
            math.acos(
                max(
                    -1.0,
                    min(
                        1.0,
                        float(vin_b @ vout_b) / (np.linalg.norm(vin_b) * np.linalg.norm(vout_b)),
                    ),
                )
            )
        ),
        "pass": ca,
    }
    return res, info


def alt_km(mu: float, vinf: float, turn_deg: float) -> float:
    s = math.sin(math.radians(turn_deg) / 2)
    return mu / vinf**2 * (1 / s - 1) - R_GAN_RS_KM


def main() -> None:
    preflight_search(
        task_no=968,
        region_id="gancal5-patched-conic-with-pass-kick",
        method=MethodCapability(
            genome="R-S 2009 GanCal#5 (one structure) with an impulsive distant-pass kick",
            corrector="least squares on dates + B->C departure velocity, Kepler legs",
            capability_tags=frozenset({"ballistic", "patched-conic"}),
            git_sha="working-tree",
        ),
        script_path=Path(__file__),
        n_points=1,
    )
    enum = _enum()
    circ, a, b = enum.cell_system("gc1")
    cands = {
        c["key"]: c
        for c in json.loads((ROOT / "data/943_cell_gc_gauntlet.json").read_text())["candidates"]
    }
    _, cyc = enum.parse_cycle_key(KEY, circ, a, b)
    x = np.asarray(cands[KEY]["x_days"]) * DAY
    fl = cycle_flybys(circ, cyc, x)
    assert fl is not None
    fb = fl[0]
    period = cyc.period_s
    mu_g = circ.body("Ganymede").mu_km3_s2
    _, vg_b = circ.state("Ganymede", fb.t_s)
    y0 = np.concatenate([np.asarray(cands[KEY]["x_days"]), vg_b + fb.vinf_out])
    out: dict[str, Any] = {
        "note": "docs/notes/2026-10-07-968-jovian-nbody-positive-control.md sec. 5.2"
    }
    for kick in (False, True):
        sol = least_squares(
            lambda y, k=kick: residual(circ, y, period, k)[0],
            y0,
            method="lm",
            xtol=1e-15,
            ftol=1e-15,
            gtol=1e-15,
            max_nfev=4000,
        )
        r, info = residual(circ, sol.x, period, kick)
        vin = 0.5 * (info["vinf_B_in"] + info["vinf_B_out"])
        row = {
            "kick": kick,
            "max_abs_residual": float(np.max(np.abs(r))),
            "dates_days": sol.x[:3].tolist(),
            **info,
            "alt_B_km": alt_km(mu_g, vin, info["turn_B_deg"]),
        }
        # A's turn: A' arrival vs A departure, the departure rotated by Ganymede's advance over T.
        ta, tb = sol.x[0] * DAY, sol.x[1] * DAY
        vout_a, _ = lambert_leg(circ, "Ganymede", "Ganymede", ta, tb, 1, "low")
        _, varr_a = lambert_leg(
            circ, "Callisto", "Ganymede", sol.x[2] * DAY, ta + period, 0, "single"
        )
        adv = 2 * math.pi * period / circ.period_s("Ganymede")
        c, s_ = math.cos(adv), math.sin(adv)
        rot = np.array([[c, -s_, 0], [s_, c, 0], [0, 0, 1]])
        vo = rot @ vout_a
        turn_a = math.degrees(
            math.acos(
                max(
                    -1.0,
                    min(1.0, float(varr_a @ vo) / (np.linalg.norm(varr_a) * np.linalg.norm(vo))),
                )
            )
        )
        row["turn_A_deg"] = turn_a
        row["alt_A_km"] = alt_km(mu_g, info["vinf_A_in"], turn_a)
        out["kick" if kick else "no_kick"] = row
        print(json.dumps(row), flush=True)
        y0 = sol.x
    (OUT / "kick.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
