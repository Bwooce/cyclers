"""#942 note 6.57: deterministic dV of the ev-B / ev-A / ev-C real-ephemeris chains under the #415
bands (pre-registered in docs/notes/2026-10-05-942-943-two-working-body-generator.md 6.57).

At every interior flyby of each EXISTING converged Standish chain (no re-solve) the full-rev
directions are the chain tool's minimax ones (``interior_gate``) and the #888 gate prices the turn
deficit two ways (``impulse_beyond_bend_kms``, Oberth-credited ``impulse_periapsis_kms``). The
chain's own mid-course closure dV is added where the chain tool reports one (ev-A on DE440,
note 6.25).
Totals are per chain, per cycle and per 7 cycles (pro-rata, Russell's basis), with the band.

Usage: uv run python scripts/dv_band_942.py
"""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path
from typing import Any

import numpy as np

from cyclerfinder.search.two_working_body import (
    Cycle,
    _blocks,
    eval_lambert_legs,
    gate_cycle,
    optimise_block,
)
from cyclerfinder.verify.dv_band_acceptance import RUSSELL_BASIS_CYCLES

REPO = Path(__file__).resolve().parents[1]
DAY = 86400.0
OUT = REPO / "data" / "942_dv_band.json"

KEYS = {
    "ev-B": "k3|RE/1:1|LE>V/0s|LV>V/1l|RV/3:2|LV>E/0s",
    "ev-A": "k2|LE>V/0s|RV/1:1|LV>V/1h|LV>E/0s",
    "ev-C": "k2|LE>V/0s|LV>V/1h|LV>E/0s",
}


def _load(name: str, path: Path) -> Any:
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod  # dataclasses in the loaded script need their module registered
    spec.loader.exec_module(mod)
    return mod


CHAIN = _load("run_942_realeph_chain", REPO / "scripts" / "run_942_realeph_chain.py")


def band(per7_mps: float) -> str:
    """The #415 band of a total over 7 cycles (m/s)."""
    if per7_mps < 1.0:
        return "strictly-ballistic"
    if per7_mps < 10.0:
        return "essentially-ballistic"
    if per7_mps < 300.0:
        return "low-maintenance"
    return "above-low-maintenance"


def price_chain(real: Any, key: str, n: int, x0_days: float, y: list[float]) -> dict[str, Any]:
    """Turn-deficit prices at every interior flyby of one Standish chain (minimax directions)."""
    circ, a, b = CHAIN.ENUM.cell_system("ev")
    _, one = CHAIN.ENUM.parse_cycle_key(key, circ, a, b)
    yy = np.asarray(y)
    cyc = Cycle(one.legs * n, float(yy[-1]) * DAY)
    x = np.concatenate([[x0_days], yy[:-1]]) * DAY
    legs = eval_lambert_legs(real, cyc, x)
    assert legs is not None
    fl = []
    for blk in _blocks(real, cyc, legs)[:-1]:
        res = optimise_block(real, blk)
        assert res is not None
        fl.extend(res[0])
    rep = gate_cycle(real, fl)
    beyond = sum(e.impulse_beyond_bend_kms for e in rep.gate.encounters) * 1000.0
    peri_vals = [e.impulse_periapsis_kms for e in rep.gate.encounters]
    peri = sum(v for v in peri_vals if np.isfinite(v)) * 1000.0
    return {
        "worst_ratio": rep.gate.worst_ratio,
        "status": rep.status,
        "n_flybys": len(fl),
        "deficit_beyond_bend_mps": beyond,
        "deficit_periapsis_mps": peri,
        "periapsis_nonfinite": int(sum(not np.isfinite(v) for v in peri_vals)),
    }


def row(
    obj: str,
    model: str,
    ep: dict[str, Any],
    n: int,
    priced: dict[str, Any] | None,
    closure_mps: float,
    note: str = "",
) -> dict[str, Any]:
    beyond = priced["deficit_beyond_bend_mps"] if priced else 0.0
    peri = priced["deficit_periapsis_mps"] if priced else 0.0
    tot_b, tot_p = beyond + closure_mps, peri + closure_mps
    s = RUSSELL_BASIS_CYCLES / n
    out = {
        "object": obj,
        "model": model,
        "epoch_jd": ep["epoch_jd"],
        "cycles": n,
        "recorded_worst_ratio": ep.get("recorded_worst"),
        "priced": priced,
        "closure_dv_mps": closure_mps,
        "per_chain_mps": {"beyond_bend": tot_b, "periapsis": tot_p},
        "per_cycle_mps": {"beyond_bend": tot_b / n, "periapsis": tot_p / n},
        "per_7_cycles_mps": {"beyond_bend": tot_b * s, "periapsis": tot_p * s},
        "band_beyond_bend": band(tot_b * s),
        "band_periapsis": band(tot_p * s),
        "note": note,
    }
    print(
        f"{obj:5s} {model:8s} JD {ep['epoch_jd']:.1f} n={n} worst(rec)={ep.get('recorded_worst')} "
        f"worst(now)={priced['worst_ratio'] if priced else None} "
        f"per7: beyond {tot_b * s:.3f} m/s ({out['band_beyond_bend']}), "
        f"periapsis {tot_p * s:.3f} m/s ({out['band_periapsis']}) {note}",
        flush=True,
    )
    return out


def _last_worst(ep: dict[str, Any]) -> float | None:
    conv = [s for s in ep.get("steps", []) if s.get("converged")]
    return conv[-1].get("worst_ratio") if conv else None


def main() -> None:
    real = CHAIN.real_ephemeris("ev", "mean")
    direct = json.loads((REPO / "data" / "942_direct_route_and_gm_fix.json").read_text())
    rows: list[dict[str, Any]] = []

    # ev-B, Standish: epochs 0 and 3 from 6.23 (4-cycle passes), 2 and 4 from 6.55b, 1 = 3 cycles.
    evb_623 = direct["evB_direct"]["epochs"]
    for ie in (0, 3):
        ep = dict(evb_623[ie], recorded_worst=_last_worst(evb_623[ie]))
        pr = price_chain(real, KEYS["ev-B"], 4, ep["x0_days"], ep["final_y_dates"])
        rows.append(row("ev-B", "standish", ep, 4, pr, 0.0))
    for ie in (2, 4):
        rec = json.loads(
            (REPO / "data" / "942_evB_ramp" / f"e{ie}_step1" / "realeph_chain.json").read_text()
        )[0]
        ep = dict(rec, recorded_worst=_last_worst(rec))
        pr = price_chain(real, KEYS["ev-B"], 4, ep["x0_days"], ep["final_y_dates"])
        rows.append(row("ev-B", "standish", ep, 4, pr, 0.0))
    rec = json.loads(
        (REPO / "data" / "942_evB_ramp" / "e1_step2" / "realeph_chain.json").read_text()
    )[0]
    k3 = next(g for g in rec["grow"] if g["k"] == 3)
    # the last converged length is k = 3; its dates are not stored separately when k = 4
    # fails, so re-price only if the stored dates are the 3-cycle chain
    y = rec["final_y_dates"]
    ncl = sum(isinstance(lg, CHAIN.LambertLeg) for lg in _one_legs(KEYS["ev-B"]))
    if len(y) - 1 == 3 * ncl - 1:
        ep = dict(rec, recorded_worst=k3.get("worst_ratio"))
        pr = price_chain(real, KEYS["ev-B"], 3, ep["x0_days"], y)
        rows.append(
            row(
                "ev-B",
                "standish",
                ep,
                3,
                pr,
                0.0,
                "FLAG: 3-cycle chain (no 4-cycle chain at this epoch), scaled x 7/3",
            )
        )
    else:
        print("ev-B epoch 1: stored dates are not the 3-cycle chain; not priced", flush=True)

    # ev-A and ev-C, Standish (6.23), 5 cycles
    for obj, key in (("ev-A", "evA_direct"), ("ev-C", "evC_standish_rerun")):
        for ep0 in direct[key]["epochs"]:
            ep = dict(ep0, recorded_worst=_last_worst(ep0))
            pr = price_chain(real, KEYS[obj], 5, ep["x0_days"], ep["final_y_dates"])
            rows.append(row(obj, "standish", ep, 5, pr, 0.0))

    # DE440: ev-A closure dV (6.25; the minimax gate at lambda = 1 passes, so no deficit),
    # ev-C (6.23 rerun; ballistic, gate pass at every epoch, no closure dV)
    de = json.loads((REPO / "data" / "942_evAB_de440_closure_dv.json").read_text())
    for c in de["evA_de440"]["closure_dv"]:
        g = c["gate_at_lambda1_minimax"]
        assert g["status"] == "pass", g
        cl = float(sum(leg["dv_mid_ms"] for leg in c["legs"]))
        ep = {"epoch_jd": c["epoch_jd"], "recorded_worst": g["worst_ratio"]}
        rows.append(
            row(
                "ev-A",
                "de440",
                ep,
                5,
                None,
                cl,
                "closure dV of the full-rev legs (6.25); minimax gate pass -> deficit 0",
            )
        )
    for ep0 in direct["evC_de440_rerun"]["epochs"]:
        assert ep0["rung_pass"], ep0["epoch_jd"]
        ep = dict(ep0, recorded_worst=_last_worst(ep0))
        rows.append(
            row(
                "ev-C",
                "de440",
                ep,
                5,
                None,
                0.0,
                "ballistic shoot, gate pass -> deficit 0, no closure dV",
            )
        )

    OUT.write_text(json.dumps({"note": "6.57", "rows": rows}, indent=1, default=float))
    print("wrote", OUT.relative_to(REPO), flush=True)


def _one_legs(key: str) -> tuple:
    circ, a, b = CHAIN.ENUM.cell_system("ev")
    _, one = CHAIN.ENUM.parse_cycle_key(key, circ, a, b)
    return one.legs


if __name__ == "__main__":
    main()
