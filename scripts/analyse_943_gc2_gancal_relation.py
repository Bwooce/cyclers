"""#943: gc-2 and the Russell-Strange GanCal family (note sec. 6.20, pre-registered).

H1a: gate verdicts of gc-2 and GanCal#5 against a Callisto GM scale s. The date residuals do not
contain any GM, so the solutions themselves do not move; only the turn capacity does.
H1b: along GanCal#5's G->C 1-rev leg, the closest approach to Callisto (where gc-2 inserts its extra
Callisto encounter), and the conic elements of the legs involved.

Ideal circular gc model (both moons massive). Writes data/943_gc2_gancal_relation.json.
"""

from __future__ import annotations

import dataclasses
import importlib.util
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from cyclerfinder.search.two_working_body import (
    CircularSystem,
    Cycle,
    cycle_flybys,
    eval_lambert_legs,
    gate_cycle,
    kepler_step,
    sphere_of_influence_km,
)

DAY = 86400.0
ROOT = Path(__file__).resolve().parents[1]
GAUNTLET = ROOT / "data" / "943_cell_gc_gauntlet.json"
OUT = ROOT / "data" / "943_gc2_gancal_relation.json"
GC2_KEY = (
    "k3|LGanymede>Ganymede/1l|LGanymede>Callisto/0s|LCallisto>Callisto/1h|LCallisto>Ganymede/0s"
)
GANCAL5_KEY = "k3|LGanymede>Ganymede/1l|LGanymede>Callisto/1h|LCallisto>Ganymede/0s"


def _load(name: str, path: Path) -> Any:
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


ENUM = _load("run_942_enumerate", ROOT / "scripts" / "run_942_enumerate.py")


def scaled(circ: CircularSystem, s: float) -> CircularSystem:
    """The same system with Callisto's GM multiplied by ``s``."""
    cal = circ.body("Callisto")
    over = dict(circ.flyby_overrides)
    over["Callisto"] = dataclasses.replace(cal, mu_km3_s2=cal.mu_km3_s2 * s)
    return dataclasses.replace(circ, flyby_overrides=over)


def verdict(circ: CircularSystem, cycle: Cycle, x: np.ndarray, s: float) -> dict[str, Any]:
    sysm = scaled(circ, s)
    fl = cycle_flybys(sysm, cycle, x)
    assert fl is not None
    rep = gate_cycle(sysm, fl)
    cal = [e for f, e in zip(fl, rep.gate.encounters, strict=True) if f.body == "Callisto"]
    return {
        "s": s,
        "status": rep.status,
        "worst_ratio": rep.gate.worst_ratio,
        "callisto_turn_deg": [f.turn_deg for f in fl if f.body == "Callisto"],
        "callisto_ratio": max((e.ratio for e in cal), default=0.0),
        "callisto_status": [e.status for e in cal],
    }


def elements(r: np.ndarray, v: np.ndarray, mu: float) -> tuple[float, float]:
    """(a km, e) of the conic through state (r, v)."""
    rn = float(np.linalg.norm(r))
    a = 1.0 / (2.0 / rn - float(v @ v) / mu)
    e_vec = ((float(v @ v) - mu / rn) * r - float(r @ v) * v) / mu
    return a, float(np.linalg.norm(e_vec))


def main() -> None:
    cands = {c["key"]: c for c in json.loads(GAUNTLET.read_text())["candidates"]}
    circ, a, b = ENUM.cell_system("gc")
    note = "docs/notes/2026-10-05-942-943-two-working-body-generator.md 6.20"
    out: dict[str, Any] = {"note": note}
    solved = {}
    for name, key in (("gc-2", GC2_KEY), ("GanCal#5", GANCAL5_KEY)):
        _, cyc = ENUM.parse_cycle_key(key, circ, a, b)
        x = np.asarray(cands[key]["x_days"]) * DAY
        solved[name] = (cyc, x)

    # H1a: verdicts against the Callisto GM scale.
    grid = [float(v) for v in np.logspace(-3, 0, 13)]
    for name, (cyc, x) in solved.items():
        rows = [verdict(circ, cyc, x, s) for s in grid]
        out[f"H1a_{name}"] = rows
        for r in rows:
            print(
                f"{name:9s} s={r['s']:.4f} status={r['status']:13s} "
                f"Callisto ratio {r['callisto_ratio']:.3f} turns "
                f"{[round(t, 2) for t in r['callisto_turn_deg']]} {r['callisto_status']}",
                flush=True,
            )
    cyc, x = solved["gc-2"]

    def feasible(s: float) -> bool:
        return bool(verdict(circ, cyc, x, s)["callisto_ratio"] <= 1.0)

    lo, hi = 1e-3, 1.0
    assert feasible(hi) and not feasible(lo)
    while hi / lo - 1.0 > 1e-4:
        mid = math.sqrt(lo * hi)
        if feasible(mid):
            hi = mid
        else:
            lo = mid
    out["H1a_gc2_s_star"] = hi
    out["H1a_gc2_at_s_star"] = verdict(circ, cyc, x, hi)
    print(f"gc-2: Callisto encounters feasible down to s* = {hi:.5f}", flush=True)

    # H1b: closest approach of GanCal#5's G->C 1-rev leg to Callisto.
    cyc5, x5 = solved["GanCal#5"]
    legs5 = eval_lambert_legs(circ, cyc5, x5)
    assert legs5 is not None
    lam5 = [cyc5.legs[i] for i in cyc5.lambert_index]
    j = next(k for k, lg in enumerate(lam5) if lg.frm == "Ganymede" and lg.to == "Callisto")
    leg = legs5[j]
    r0, w0 = circ.state("Ganymede", leg.t_dep)
    v0 = w0 + leg.vinf_dep
    ts = np.arange(leg.t_dep + 0.01 * DAY, leg.t_arr - 0.5 * DAY, 0.01 * DAY)
    best = (math.inf, 0.0)
    for t in ts:
        r, _ = kepler_step(r0, v0, float(t - leg.t_dep), circ.mu)
        rc, _ = circ.state("Callisto", float(t))
        d = float(np.linalg.norm(r - rc))
        if d < best[0]:
            best = (d, float(t))
    soi = sphere_of_influence_km(circ, "Callisto")
    a5, e5 = elements(r0, v0, circ.mu)
    legs2 = eval_lambert_legs(circ, cyc, x)
    assert legs2 is not None
    lam2 = [cyc.legs[i] for i in cyc.lambert_index]
    el2 = {}
    for k, lg in enumerate(lam2):
        rr, ww = circ.state(lg.frm, legs2[k].t_dep)
        el2[f"{lg.frm}>{lg.to}"] = {
            "days": (legs2[k].t_arr - legs2[k].t_dep) / DAY,
            "a_km_e": elements(rr, ww + legs2[k].vinf_dep, circ.mu),
        }
    out["H1b"] = {
        "gancal5_GC_leg_days": (leg.t_arr - leg.t_dep) / DAY,
        "gancal5_GC_leg_a_km_e": (a5, e5),
        "min_distance_to_callisto_km": best[0],
        "at_days_after_departure": (best[1] - leg.t_dep) / DAY,
        "callisto_soi_km": soi,
        "gc2_insertion_days_after_departure": el2["Ganymede>Callisto"]["days"],
        "gc2_legs": el2,
    }
    print(json.dumps(out["H1b"], indent=1), flush=True)
    OUT.write_text(json.dumps(out, indent=1, default=float))
    print(f"wrote {OUT}", flush=True)


if __name__ == "__main__":
    main()
