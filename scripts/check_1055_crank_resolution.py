"""#1055 sec. 13.6 instrument gate (gc-1 e2, jup365): the crank RESOLUTION of each down point.

For each point of data/1055_gc1e2/down_forced.json (read-only; the running downj stage owns that
file), at its sigma: the residual's projection on the weakest column-scaled left singular vector,
the Newton move it implies along the weakest right singular vector (`weak_direction_move`), and
that move's effect on every node's V_inf in / out angle against the chain tool's (deg).
Amendment 14 is judged only on points whose largest such effect is below 0.03 deg.

Diagnostic (no verdict weight): the signed pump and crank of the Callisto -> Callisto leg's V_inf
(C1 outbound) about Callisto's jup365 velocity (crank 0 = radially outward in Callisto's orbit
plane, 90 = along Callisto's orbit normal), for the chain and each point.

    uv run python scripts/check_1055_crank_resolution.py [sigma ...]

Writes data/1055_gc1e2/crank_resolution.json (one entry per sigma, merged on re-runs).
"""

from __future__ import annotations

import importlib.util
import json
import math
import sys
from pathlib import Path
from typing import Any

import numpy as np

from cyclerfinder.nbody.degenerate_family import weak_direction_move

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "1055_gc1e2"


def _load(name: str, path: Path) -> Any:
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


DRV = _load("run_1055_seedfix", ROOT / "scripts" / "run_1055_seedfix.py")


def pump_crank(t: Any, k: int, z: Any, which: str) -> tuple[float, float]:
    """Pump and crank of node k's V_inf (in or out) about the moon's velocity."""
    moon = t.moons[k]
    tt = float(z[7 * k + 6])
    rm, vm = t.c.eph.state(moon, tt)
    r, v = np.asarray(z[7 * k : 7 * k + 3]), np.asarray(z[7 * k + 3 : 7 * k + 6])
    vec = t.rb.asymptote(t.c.mus[moon], r, v, which)
    vh = np.asarray(vm) / float(np.linalg.norm(vm))
    nh = np.cross(rm, vm)
    nh = nh / float(np.linalg.norm(nh))
    rh = np.cross(vh, nh)
    pump = math.degrees(
        math.acos(max(-1.0, min(1.0, float(vec @ vh) / float(np.linalg.norm(vec)))))
    )
    crank = math.degrees(math.atan2(float(vec @ nh), float(vec @ rh)))
    return pump, crank


def main() -> None:
    rec = json.loads((OUT / "down_forced.json").read_text())
    want = [float(x) for x in sys.argv[1:]]
    t = DRV.Target("gc1e2", OUT)
    path = OUT / "crank_resolution.json"
    res: dict[str, Any] = json.loads(path.read_text()) if path.exists() else {}
    t.set_sigma(1.0)
    seed = np.asarray(t.seed)
    res["chain_C1_out_pump_crank"] = pump_crank(t, 1, seed, "out")
    res["chain_C2_in_pump_crank"] = pump_crank(t, 2, seed, "in")
    for p in rec["points"]:
        s = float(p["sigma"])
        if want and not any(abs(s - w) < 1e-12 for w in want):
            continue
        z = np.asarray(p["z"])
        t.set_sigma(s)
        r, j = t.rj(z)
        wm = weak_direction_move(r, j, z, len(t.moons))
        base = DRV._member_rows(t, z, s)
        moved = DRV._member_rows(t, z + wm.dz, s)
        eff = [
            max(
                abs(b["vinf_in_angle_deg"] - a["vinf_in_angle_deg"]),
                abs(b["vinf_out_angle_deg"] - a["vinf_out_angle_deg"]),
            )
            for a, b in zip(base, moved, strict=True)
        ]
        t.set_sigma(s)
        row = {
            "sigma": s,
            "accept": p["accept"],
            "sv_min": wm.sv_min,
            "u_dot_r": wm.u_dot_r,
            "weak_move_scaled": wm.move_scaled,
            "angle_effect_deg_per_node": eff,
            "resolution_deg": max(eff),
            "resolved": max(eff) < 0.03,
            "C1_out_pump_crank": pump_crank(t, 1, z, "out"),
            "C2_in_pump_crank": pump_crank(t, 2, z, "in"),
            "angles_in_out_deg": [(n["vinf_in_angle_deg"], n["vinf_out_angle_deg"]) for n in base],
        }
        res[str(s)] = row
        path.write_text(json.dumps(res, indent=1))
        print(
            f"sigma {s}: sv_min {wm.sv_min:.3e} u.r {wm.u_dot_r:.3e} resolution "
            f"{max(eff):.4f} deg (per node {[round(e, 4) for e in eff]}); "
            f"C1 out pump/crank {row['C1_out_pump_crank']} "
            f"(chain {res['chain_C1_out_pump_crank']})",
            flush=True,
        )


if __name__ == "__main__":
    main()
