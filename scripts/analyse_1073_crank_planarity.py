"""#1073 pre-registration analysis: are the 171 strong_1065 Pluto members planar, and do their
planar-crank siblings pass the #1065 strong filters?

Read-only on the committed #1065 gauntlets (no search, no solve). For each strong_1065 member:
1. the largest out-of-plane fraction |v_z| / |V_inf| over its flyby V_inf vectors, with the
   resonant-leg cranks chosen by the generator's own convention (``cycle_flybys``, minimax);
2. the planar-crank sibling: every resonant (full-revolution) leg's crank set to one of its two
   in-plane values (the combination with the lowest worst demanded/available turn ratio per block),
   then the #1065 strong filters on that sibling: the #888/#937 gate (near-180 rule included), the
   largest required Charon periapsis inside the run's Charon SOI, the moving-Pluto surface, and the
   #973 unscheduled-pass screen.

The dates are the member's own (the date residual does not depend on the cranks).

    uv run python scripts/analyse_1073_crank_planarity.py \
        --out data/1073_pluto_cg/crank_planarity.json
"""

from __future__ import annotations

import argparse
import importlib.util
import itertools
import json
import math
import sys
from pathlib import Path
from typing import Any

import numpy as np

import cyclerfinder.search.two_working_body as twb
from cyclerfinder.search.two_working_body_binary import primary_min_distance_km

REPO = Path(__file__).resolve().parents[1]
DAY = 86400.0
PLUTO_RADIUS_KM = 1188.3  # as scripts/gauntlet_998.py (Nimmo et al. 2017)


def _load(name: str, path: Path) -> Any:
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


R1065 = _load(
    "run_1065_barycentric_enumerate", REPO / "scripts" / "run_1065_barycentric_enumerate.py"
)


def planar_flybys(system: Any, cyc: Any, x: np.ndarray, **_: Any) -> list[twb.Flyby]:
    """Flybys with every resonant crank in-plane (best in-plane combination per block)."""
    legs = twb.eval_lambert_legs(system, cyc, x)
    assert legs is not None
    out: list[twb.Flyby] = []
    for blk in twb._blocks(system, cyc, legs):
        body = system.body(blk.body)
        if body.massless and not blk.fixed:
            out.append(twb.Flyby(blk.body, blk.t_arr, blk.v_in, blk.v_out))
            continue
        vinf = float(np.linalg.norm(blk.v_in))
        avail = twb._available(body, vinf)
        opts = []
        for leg, t in blk.fixed:
            circ = twb.resonant_circle(system, leg, t, vinf)
            assert circ is not None
            _, _, _, e1, e2 = circ
            ph = math.atan2(-e1[2], e2[2])  # u_z(phi) = 0 (the resonant circle meets the plane)
            opts.append([ph, ph + math.pi])
        best: tuple[float, list[Any]] | None = None
        for combo in itertools.product(*opts) if opts else [()]:
            seq = twb._block_vectors(system, blk, list(combo), [])
            assert seq is not None
            rat = max(
                twb.bend_angle(seq[2 * i], seq[2 * i + 1]) / avail for i in range(len(seq) // 2)
            )
            if best is None or rat < best[0]:
                best = (rat, seq)
        assert best is not None
        out.extend(twb._flybys_from_seq(system, blk, best[1]))
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    # the #973 screen recomputes the flybys through its module-level cycle_flybys
    minimax = R1065.R973.cycle_flybys
    rows = []
    for cell, moon in R1065.CELLS.items():
        system = R1065.bary_system(moon)
        soi = twb.sphere_of_influence_km(system, "Charon")
        mu_c = system.body("Charon").mu_km3_s2
        gaunt = json.loads(
            (REPO / "data" / "1065_barycentric" / f"{cell}_gauntlet.json").read_text()
        )
        for c in gaunt["candidates"]:
            if not c.get("strong_1065"):
                continue
            _, cyc = R1065.ENUM.parse_cycle_key(c["key"], system, "Charon", moon)
            x = np.asarray(c["x_days"]) * DAY
            fl_mm = minimax(system, cyc, x)
            oop = max(
                max(abs(f.vinf_in[2]), abs(f.vinf_out[2])) / float(np.linalg.norm(f.vinf_in))
                for f in fl_mm
            )
            fl = planar_flybys(system, cyc, x)
            rep = twb.gate_cycle(system, fl)
            rps = []
            for f in fl:
                if f.body != "Charon":
                    continue
                d = math.radians(f.turn_deg)
                rps.append(
                    mu_c / f.vinf_kms**2 * (1 / math.sin(d / 2) - 1) if d > 1e-9 else math.inf
                )
            pa = primary_min_distance_km(system, cyc, x, fl).min_km
            R1065.R973.cycle_flybys = planar_flybys
            try:
                un = R1065._unscheduled(system, cyc, x)["unscheduled_status"]
            finally:
                R1065.R973.cycle_flybys = minimax
            strong = bool(
                rep.status == "pass"
                and max(rps) <= soi
                and pa > PLUTO_RADIUS_KM
                and un.startswith("pass")
            )
            rows.append(
                {
                    "cell": cell,
                    "index": c["index"],
                    "key": c["key"],
                    "n_resonant_legs": c["key"].count("|R"),
                    "max_out_of_plane_fraction": oop,
                    "planar_member": bool(oop < 1e-9),
                    "minimax_worst_ratio": c["worst_ratio"],
                    "sibling_worst_ratio": rep.gate.worst_ratio,
                    "sibling_gate": rep.status,
                    "sibling_max_rp_charon_km": max(rps),
                    "sibling_pluto_min_km": pa,
                    "sibling_unscheduled": un,
                    "sibling_strong": strong,
                }
            )
    oop = np.array([r["max_out_of_plane_fraction"] for r in rows])
    summary = {
        "n_members": len(rows),
        "n_planar": sum(r["planar_member"] for r in rows),
        "n_with_resonant_leg": sum(r["n_resonant_legs"] > 0 for r in rows),
        "out_of_plane_fraction_percentiles_0_10_50_90_100": np.percentile(
            oop, [0, 10, 50, 90, 100]
        ).tolist(),
        "n_sibling_ratio_below_1": sum(r["sibling_worst_ratio"] < 1 for r in rows),
        "n_sibling_gate_pass": sum(r["sibling_gate"] == "pass" for r in rows),
        "n_sibling_strong": sum(r["sibling_strong"] for r in rows),
        "n_spatial_with_strong_sibling": sum(
            r["sibling_strong"] and not r["planar_member"] for r in rows
        ),
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps({"summary": summary, "members": rows}, indent=1, default=float))
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
