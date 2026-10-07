"""#1041: the continuous-gravity EGGIE in the paper's ideal model, quasi-periodic criterion
(pre-registration: docs/notes/2026-10-08-1041-eggie-quasiperiodic-continuous.md).

Reuses the #968 rung-(b) open-chain machinery (``scripts/run_968_rungb.py``: moon-relative nodes,
forward-backward legs, periapsis gauges, pinned end asymptotes, damped Newton with a step cap,
sigma continuation), loaded read-only, with this file's model: Io, Europa and Ganymede on the
paper's circular ideal orbits (``resonant_conic.ideal_t_syn``, Laplace angle 180 deg), all three
massive and scaled by sigma.

    uv run python scripts/run_1041_eggie.py --stage pc --n-cycles 3
    uv run python scripts/run_1041_eggie.py --stage sigma --n-cycles 3
    uv run python scripts/run_1041_eggie.py --stage down --n-cycles 3
    uv run python scripts/run_1041_eggie.py --stage ias15 --n-cycles 3
    uv run python scripts/run_1041_eggie.py --stage drift --n-cycles 3

Checkpoints in data/1041_eggie/; live log in data/1041_eggie/live/ (gitignored).
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import math
import sys
from pathlib import Path
from typing import Any

import numpy as np
from numpy.typing import NDArray

from cyclerfinder.core.lambert import lambert
from cyclerfinder.core.satellites import SATELLITES
from cyclerfinder.data.method_capability import MethodCapability
from cyclerfinder.data.preflight import preflight_search
from cyclerfinder.nbody.jovian import MU_JUPITER_KM3_S2, JovianRestrictedNBody, periapsis_node
from cyclerfinder.search.two_working_body import (
    Cycle,
    Flyby,
    LambertLeg,
    correct_dates,
    cycle_flybys,
    encounter_self_consistency,
    eval_lambert_legs,
    gate_cycle,
)

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "1041_eggie"
DAY = 86400.0
Arr = NDArray[np.float64]


def _load(name: str, path: Path) -> Any:
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


R = _load("run_968_rungb", ROOT / "scripts" / "run_968_rungb.py")
E23 = _load("run_1023_eggie", ROOT / "scripts" / "run_1023_eggie.py")
IO, EUR, GAN = "Io", "Europa", "Ganymede"
MOONS = (IO, EUR, GAN)
R.OUT = OUT
R.ALL4 = MOONS
R.FLOOR_KM = {IO: 25.0, EUR: 25.0, GAN: 25.0}  # the paper's floor (pre-registration sec. 3)


class CircEphem:
    """The paper's circular ideal ephemeris with the ``state``/``position`` interface."""

    def __init__(self) -> None:
        self.circ, self.period = E23.system("coded", 180.0)

    def state(self, moon: str, t_sec: float) -> tuple[Arr, Arr]:
        r, v = self.circ.state(moon, t_sec)
        return np.asarray(r, dtype=np.float64), np.asarray(v, dtype=np.float64)

    def position(self, moon: str, t_sec: float) -> Arr:
        return self.state(moon, t_sec)[0]


class Chain3(R.Chain):  # type: ignore[misc, name-defined]
    def set_sigma(self, s: float) -> None:
        self.mus = {k: s * SATELLITES[k].mu_km3_s2 for k in self.force}
        self.surf = {k: s * SATELLITES[k].radius_eq_km for k in self.force}


def advance(eph: CircEphem, moon: str, dt: float) -> Arr:
    """Rotation by the moon's own advance over dt (about the z axis)."""
    a = 2.0 * math.pi * dt / eph.circ.period_s(moon)
    c, s = math.cos(a), math.sin(a)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])


def stage_pc(n: int) -> None:
    eph = CircEphem()
    root = json.loads((ROOT / "data/1023_eggie/pc_roots_coded.json").read_text())[0]
    assert root["gate25"] == "pass" and abs(root["dist_table4"] - 0.173) < 1e-3
    legs_one = (
        LambertLeg(EUR, GAN, 0, "single"),
        LambertLeg(GAN, GAN, 1, "high"),
        LambertLeg(GAN, IO, 1, "low"),
        LambertLeg(IO, EUR, 1, "high"),
    )
    cyc = Cycle(legs_one * n, n * eph.period)
    x1 = np.asarray(root["x_days"]) * DAY
    x0 = np.concatenate([x1 + k * eph.period for k in range(n)])
    sol = correct_dates(eph.circ, cyc, x0, tol_kms=1e-9)
    fl = cycle_flybys(eph.circ, cyc, sol.x) if sol.converged else None
    miss = encounter_self_consistency(eph.circ, cyc, sol.x) if sol.converged else float("inf")
    gate = gate_cycle(eph.circ, fl, alt_floor_km=25.0).status if fl else "none"
    out: dict[str, Any] = {
        "n_cycles": n,
        "converged": bool(sol.converged),
        "max_residual_kms": sol.max_abs_residual_kms,
        "miss_km": miss,
        "gate25": gate,
        "x_days": (sol.x / DAY).tolist(),
    }
    if fl:
        legs_ev = eval_lambert_legs(eph.circ, cyc, sol.x)
        assert legs_ev is not None
        out["first_departure_vinf"] = np.asarray(legs_ev[0].vinf_dep).tolist()
        out["flybys"] = [
            {
                "body": f.body,
                "t_s": float(f.t_s),
                "vinf_in": np.asarray(f.vinf_in).tolist(),
                "vinf_out": np.asarray(f.vinf_out).tolist(),
                "vinf_kms": float(f.vinf_kms),
                "turn_deg": float(f.turn_deg),
            }
            for f in sorted(fl, key=lambda f: f.t_s)
        ]
    (OUT / f"pc_chain_n{n}.json").write_text(json.dumps(out, indent=1))
    R._log(
        f"pc n{n}: converged={out['converged']} res {out['max_residual_kms']:.2e} miss {miss:.2e} "
        f"gate25 {gate} vinf {[round(f['vinf_kms'], 3) for f in out.get('flybys', [])]}"
    )


def build(n: int) -> Any:
    """Open chain E0, (G, G, I, E) x n from the patched-conic n-cycle chain (pre-registration 2)."""
    eph = CircEphem()
    pc = json.loads((OUT / f"pc_chain_n{n}.json").read_text())
    fl = pc["flybys"]  # time order: G, G, I, ..., E (the closing Europa flyby at t0 + nT)
    t0 = pc["x_days"][0] * DAY
    e_last = fl[-1]
    rot = advance(eph, EUR, n * eph.period)
    # E0: inbound = the closing flyby's inbound rotated back by Europa's own advance (per moon).
    e0 = {
        "body": EUR,
        "t_s": t0,
        "vinf_in": rot.T @ np.asarray(e_last["vinf_in"]),
        "vinf_out": np.asarray(pc["first_departure_vinf"]),  # the chain's first departure
        "vinf_kms": e_last["vinf_kms"],
    }
    e_end = {
        "body": EUR,
        "t_s": e_last["t_s"],
        "vinf_in": np.asarray(e_last["vinf_in"]),
        "vinf_out": rot @ np.asarray(pc["first_departure_vinf"]),
        "vinf_kms": e_last["vinf_kms"],
    }
    seq = [e0, *fl[:-1], e_end]
    z = []
    for f in seq:
        # No SOI cap (#1043): a small turn needs a periapsis beyond 0.6 SOI, and the lane's default
        # clamp there would give the seed a wrong turn (0.7 deg at Io, about 0.1 km/s).
        r, v, _ = periapsis_node(
            f["body"],
            f["t_s"],
            np.asarray(f["vinf_in"]),
            np.asarray(f["vinf_out"]),
            eph,  # type: ignore[arg-type]
            max_offset_km=1.0e12,
        )
        rm, vm = eph.state(f["body"], f["t_s"])
        z.extend([*(r - rm), *(v - vm), f["t_s"]])
    c = Chain3(
        eph,
        np.zeros(6),
        t0 + DAY,
        [f["body"] for f in seq],
        np.asarray(z, dtype=np.float64),
        [float(f["vinf_kms"]) for f in seq],
        np.asarray(seq[0]["vinf_in"], dtype=np.float64),
        np.asarray(seq[-1]["vinf_out"], dtype=np.float64),
        force=MOONS,
    )
    return c


R.build = build


def unscheduled(c: Any, z: Arr, sigma: float) -> list[dict[str, Any]]:
    from scipy.integrate import solve_ivp

    nd = R.nodes(c, z)
    found = []
    for k in range(c.m - 1):
        x0, _, t0, _ = nd[k]
        t1 = nd[k + 1][2]
        sol = solve_ivp(
            lambda t, y: R.accel(c, y, t),
            (t0, t1),
            x0,
            method="DOP853",
            rtol=R.RTOL,
            atol=R.ATOL,
            dense_output=True,
        )
        ts = np.linspace(t0, t1, int((t1 - t0) / (0.002 * DAY)) + 2)
        ys = sol.sol(ts)
        for moon in MOONS:
            mu = SATELLITES[moon].mu_km3_s2 * sigma
            hill = SATELLITES[moon].sma_km * (mu / (3 * MU_JUPITER_KM3_S2)) ** (1 / 3)
            dist = np.array(
                [
                    np.linalg.norm(ys[:3, i] - c.eph.position(moon, float(t)))
                    for i, t in enumerate(ts)
                ]
            )
            for i in range(1, len(ts) - 1):
                if dist[i] <= dist[i - 1] and dist[i] <= dist[i + 1] and dist[i] < hill:
                    found.append(
                        {
                            "leg": k,
                            "moon": moon,
                            "t_days": float(ts[i] / DAY),
                            "dist_km": float(dist[i]),
                        }
                    )
    return found


R.unscheduled = unscheduled


def stage_ias15(n: int) -> None:
    c = build(n)
    path = OUT / f"sigma_n{n}.json"
    rec = json.loads(path.read_text())
    q = rec["points"][-1]
    assert q["sigma"] >= 1.0
    z = np.asarray(q["z"])
    c.set_sigma(1.0)
    nd = R.nodes(c, z)
    prop = JovianRestrictedNBody()
    rows = []
    for k in range(c.m - 1):
        xp, _, tp, _ = nd[k]
        xn, _, tn, _ = nd[k + 1]
        tm = 0.5 * (tp + tn)
        for lab, x0, t0 in (("fwd", xp, tp), ("bwd", xn, tn)):
            xd, _ = R.arc(c, x0, t0, tm)
            a = prop.propagate(
                x0[:3], x0[3:], t0, tm, moons=MOONS, cache=c.eph, max_wall_sec=3000.0
            )  # type: ignore[arg-type]
            if not a.converged or abs(a.t1_sec - tm) > 1e-6:
                rows.append({"leg": k, "dir": lab, "status": "timeout"})
                continue
            rows.append(
                {
                    "leg": k,
                    "dir": lab,
                    "status": "ok",
                    "dr_km": float(np.linalg.norm(a.r_km - xd[:3])),
                    "dv_kms": float(np.linalg.norm(a.v_km_s - xd[3:])),
                }
            )
    ok = all(r["status"] == "ok" and r["dr_km"] < 1e-2 and r["dv_kms"] < 1e-7 for r in rows)
    rec["ias15"] = {"rows": rows, "pass": ok}
    path.write_text(json.dumps(rec))
    R._log(f"ias15 n{n}: pass={ok} max dr {max(r.get('dr_km', float('nan')) for r in rows):.2e}")


def stage_drift(n: int) -> None:
    """Quasi-periodicity, measured per moon at sigma = 1 (pre-registration sec. 3)."""
    c = build(n)
    rec = json.loads((OUT / f"sigma_n{n}.json").read_text())
    q = rec["points"][-1]
    assert q["sigma"] >= 1.0
    z = np.asarray(q["z"])
    c.set_sigma(1.0)
    rows = []
    for k in range(c.m - 4):
        moon = c.moons[k]
        assert c.moons[k + 4] == moon
        rot = advance(c.eph, moon, c.eph.period)
        a, b = z[7 * k : 7 * k + 6], z[7 * (k + 4) : 7 * (k + 4) + 6]
        da = np.concatenate([rot @ a[:3], rot @ a[3:]])
        nk, nk4 = q["nodes"][k], q["nodes"][k + 4]
        rows.append(
            {
                "node": k,
                "moon": moon,
                "d_vinf_kms": nk4["vinf_kms"] - nk["vinf_kms"],
                "d_alt_km": nk4["alt_km"] - nk["alt_km"],
                "d_rel_pos_km": float(np.linalg.norm(b[:3] - da[:3])),
                "d_rel_vel_kms": float(np.linalg.norm(b[3:] - da[3:])),
            }
        )
    rec["drift"] = rows
    (OUT / f"sigma_n{n}.json").write_text(json.dumps(rec))
    R._log(f"drift n{n}: {json.dumps(rows)}")


MARCH_LEGS = (
    (EUR, GAN, 0, "single"),
    (GAN, GAN, 1, "high"),
    (GAN, IO, 1, "low"),
    (IO, EUR, 1, "high"),
)


def _leg(
    eph: CircEphem, frm: str, to: str, nrev: int, br: str, t0: float, t1: float
) -> tuple[Arr, Arr] | None:
    r1, w1 = eph.state(frm, t0)
    r2, w2 = eph.state(to, t1)
    if t1 <= t0:
        return None
    try:
        sols = lambert(r1, r2, t1 - t0, mu=eph.circ.mu, max_revs=nrev)
    except Exception:
        return None
    sol = next((q for q in sols if q.n_revs == nrev and q.branch == br), None)
    return None if sol is None else (np.asarray(sol.v1 - w1), np.asarray(sol.v2 - w2))


def _cycle_legs(eph: CircEphem, dates: list[float]) -> list[tuple[Arr, Arr]] | None:
    """dates = [E_k, G1, G2, I, E_{k+1}] (s); the four legs' (V_inf dep, V_inf arr)."""
    out = []
    for (frm, to, nr, br), t0, t1 in zip(MARCH_LEGS, dates[:-1], dates[1:], strict=True):
        lg = _leg(eph, frm, to, nr, br, t0, t1)
        if lg is None:
            return None
        out.append(lg)
    return out


def march(start: str, n_max: int = 10) -> None:
    from scipy.optimize import least_squares

    eph = CircEphem()
    root = json.loads((ROOT / "data/1023_eggie/pc_roots_coded.json").read_text())[0]
    x = np.asarray(root["x_days"]) * DAY  # [E0, G1, G2, I] (Lambert-leg starts)
    t_cyc = eph.period
    cycles = []
    if start == "root":
        dates = [*x.tolist(), x[0] + t_cyc]
        legs = _cycle_legs(eph, dates)
        assert legs is not None
        # E_0's inbound: root (iii)'s closing I>E arrival, rotated back by Europa's own advance.
        v_in_e0 = advance(eph, EUR, t_cyc).T @ legs[3][1]
    else:  # Table 4 seed, E_1 fixed at E_0 + 28.22 d, no inbound at E_0
        tofs = np.array([1.59, 8.60, 7.34]) * DAY
        e0 = float(x[0])
        e1 = e0 + 28.22 * DAY

        def res1(y: Arr) -> Arr:
            d = [e0, *y.tolist(), e1]
            lg = _cycle_legs(eph, d)
            if lg is None:
                return np.full(3, 1e3)
            return np.array(
                [np.linalg.norm(lg[k][1]) - np.linalg.norm(lg[k + 1][0]) for k in range(3)]
            )

        sol = least_squares(
            res1, e0 + np.cumsum(tofs), x_scale=DAY, xtol=1e-15, ftol=1e-15, gtol=1e-15
        )
        dates = [e0, *sol.x.tolist(), e1]
        legs = _cycle_legs(eph, dates)
        assert legs is not None and float(np.max(np.abs(res1(sol.x)))) < 1e-9, (
            "Table-4 cycle 1 did not solve"
        )
        v_in_e0 = None
    for k in range(n_max):
        if k > 0:
            e_k = dates[-1]
            v_in = legs[3][1]  # previous I>E arrival
            prev = np.asarray(dates[1:]) + t_cyc

            def res(y: Arr, e_k: float = e_k, v_in: Arr = v_in) -> Arr:
                d = [e_k, *y.tolist()]
                lg = _cycle_legs(eph, d)
                if lg is None:
                    return np.full(4, 1e3)
                r = [np.linalg.norm(lg[0][0]) - np.linalg.norm(v_in)]
                r += [np.linalg.norm(lg[j][1]) - np.linalg.norm(lg[j + 1][0]) for j in range(3)]
                return np.asarray(r)

            sol = least_squares(
                res, prev, x_scale=DAY, xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=2000
            )
            rmax = float(np.max(np.abs(res(sol.x))))
            if rmax > 1e-9:
                R._log(
                    f"march {start}: cycle {k + 1} corrector failed (max residual {rmax:.2e}); stop"
                )
                cycles.append({"cycle": k + 1, "converged": False, "max_residual": rmax})
                break
            dates = [e_k, *sol.x.tolist()]
            legs = _cycle_legs(eph, dates)
            assert legs is not None
            v_in_e = v_in
        else:
            v_in_e = v_in_e0
        fl = []
        if v_in_e is not None:
            fl.append(Flyby(EUR, dates[0], np.asarray(v_in_e), legs[0][0]))
        for j, body in ((1, GAN), (2, GAN), (3, IO)):
            fl.append(Flyby(body, dates[j], legs[j - 1][1], legs[j][0]))
        row: dict[str, Any] = {
            "cycle": k + 1,
            "converged": True,
            "dates_days": [d / DAY for d in dates],
        }
        for floor_lab, floor in (("25km", 25.0), ("project", None)):
            rep = gate_cycle(eph.circ, fl, alt_floor_km=floor)
            row[f"gate_{floor_lab}"] = rep.status
            row[f"encounters_{floor_lab}"] = [
                {
                    "body": e.body,
                    "vinf": 0.5 * (e.vinf_in_kms + e.vinf_out_kms),
                    "turn": e.demanded_turn_deg,
                    "ratio": e.ratio,
                    "status": e.status,
                }
                for e in rep.gate.encounters
            ]
        cycles.append(row)
        enc = row["encounters_25km"]
        R._log(
            f"march {start} cycle {k + 1}: gate25 {row['gate_25km']} project {row['gate_project']} "
            + " ".join(
                f"{e['body'][0]} {e['vinf']:.3f}/{e['turn']:.1f}/{e['ratio']:.2f}" for e in enc
            )
        )
    first = {
        lab: next((c["cycle"] for c in cycles if c.get(f"gate_{lab}") not in (None, "pass")), None)
        for lab in ("25km", "project")
    }
    (OUT / f"march_{start}.json").write_text(
        json.dumps({"start": start, "cycles": cycles, "first_gate_fail": first}, indent=1)
    )
    R._log(f"march {start}: first gate fail {first}")


def main() -> None:
    preflight_search(
        task_no=1041,
        region_id="eggie-paper-model-quasiperiodic-continuous",
        method=MethodCapability(
            genome="EGGIE n-cycle open chain (one structure), the paper's ideal Galilean model",
            corrector="forward-backward shooting, DOP853 + STM; sigma continuation; pinned ends",
            capability_tags=frozenset({"ballistic", "n-body"}),
            git_sha="working-tree",
        ),
        script_path=Path(__file__),
        n_points=1,
    )
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--stage", choices=["pc", "sigma", "down", "ias15", "drift", "march"], required=True
    )
    ap.add_argument("--n-cycles", type=int, default=3)
    ap.add_argument("--start", choices=["root", "table4"], default="root")
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    n = args.n_cycles
    if args.stage == "pc":
        stage_pc(n)
    elif args.stage == "sigma":
        R.stage_sigma(n)
    elif args.stage == "down":
        R.stage_down(n)
    elif args.stage == "march":
        march(args.start)
    elif args.stage == "ias15":
        stage_ias15(n)
    else:
        stage_drift(n)


if __name__ == "__main__":
    main()
