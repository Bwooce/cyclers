"""#968 rung (b), step 1 (note sec. 9, amendment 8): reconstruct the GanEur#316@2019 real-ephemeris
patched-conic chain from the stored #943 result, with the chain tool PINNED at commit 1e4b7aeb (the
version that produced data/943_ganeur316_realeph/n10_rs2019_rel/), loaded from ``git show`` into a
temporary module (the live file is being edited by another agent and is not imported).

Positive control of the reconstruction: max residual < 1e-6, gate worst 0.804 (to 1e-3), minimum
required altitude 982 km (to 1 km). Writes data/968_rungb/seed_chain.json (flyby list, chain period,
tool hash).

    uv run python scripts/run_968_rungb_seed.py
"""

from __future__ import annotations

import importlib.util
import json
import math
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any

import numpy as np

from cyclerfinder.data.method_capability import MethodCapability
from cyclerfinder.data.preflight import preflight_search

ROOT = Path(__file__).resolve().parents[1]
TOOL_HASH = "1e4b7aeb"
SRC = ROOT / "data" / "943_ganeur316_realeph" / "n10_rs2019_rel"
OUT = ROOT / "data" / "968_rungb"
KEY = (
    "k7|LGanymede>Ganymede/1h|HGanymede/3,1,p|LGanymede>Ganymede/1h"
    "|LGanymede>Europa/1l|LEuropa>Ganymede/2l"
)
FIRST_EPOCH_JD = 2458597.5
N_CYCLES = 10
DAY = 86400.0


def pinned_tool() -> Any:
    text = subprocess.run(
        ["git", "show", f"{TOOL_HASH}:scripts/run_942_realeph_chain.py"],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    ).stdout
    tmp = Path(tempfile.mkdtemp()) / "chain_pinned.py"
    tmp.write_text(text)
    # The pinned tool resolves its sibling scripts relative to its own file; point it at the repo.
    text = text.replace("Path(__file__).resolve().parents[1]", f"Path({str(ROOT)!r})")
    tmp.write_text(text)
    spec = importlib.util.spec_from_file_location("chain_pinned", tmp)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    sys.modules["chain_pinned"] = mod
    spec.loader.exec_module(mod)
    return mod


def main() -> None:
    preflight_search(
        task_no=968,
        region_id="ganeur316-rs2019-chain-reconstruction",
        method=MethodCapability(
            genome="GanEur#316 10-cycle real-ephemeris patched-conic chain (stored result)",
            corrector="none (re-evaluation of the stored solution, pinned chain tool)",
            capability_tags=frozenset({"ballistic", "patched-conic", "real-ephemeris"}),
            git_sha=TOOL_HASH,
        ),
        script_path=Path(__file__),
        n_points=1,
    )
    ch = pinned_tool()
    enum = ch.ENUM
    stored = json.loads((SRC / "realeph_chain.json").read_text())
    rec = stored[0] if isinstance(stored, list) else stored
    circ, a, b = enum.cell_system("ge")
    _, one = enum.parse_cycle_key(KEY, circ, a, b)
    recall = json.loads((ROOT / "data/943_ganeur316_recall.json").read_text())
    seed = {c["key"]: c for c in recall}[KEY]
    x1 = np.asarray(seed["x_days"]) * DAY
    real = ch.real_ephemeris("ge", "auto")
    near_s = (FIRST_EPOCH_JD - 2440000.0) * DAY
    basis = ch.plane_basis(real, (a, b), near_s)
    th = {c: 2 * math.pi * x1[0] / circ.period_s(c) for c in (a, b)}
    target = (th[b] - th[a]) % (2 * math.pi)

    def g(t_s: float) -> float:
        d = ch.in_plane_longitude(real, basis, b, t_s) - ch.in_plane_longitude(real, basis, a, t_s)
        return ((d - target) % (2 * math.pi) + math.pi) % (2 * math.pi) - math.pi

    syn = circ.synodic_s(a, b)
    grid = np.linspace(near_s, near_s + 1.1 * syn, 440)
    vals = [g(float(t)) for t in grid]
    te = None
    for i in range(len(grid) - 1):
        if vals[i] * vals[i + 1] < 0 and abs(vals[i] - vals[i + 1]) < math.pi:
            lo, hi = float(grid[i]), float(grid[i + 1])
            for _ in range(60):
                mid = 0.5 * (lo + hi)
                if g(lo) * g(mid) <= 0:
                    hi = mid
                else:
                    lo = mid
            te = 0.5 * (lo + hi)
            break
    assert te is not None
    rot = ch.in_plane_longitude(real, basis, a, te) - th[a]
    t_shift = te - x1[0]
    sysm = ch.Blend(circ, real, basis, t_shift, rot, 1.0)
    legs = one.legs * N_CYCLES
    start = json.loads((SRC / "shoot_start_epoch0.json").read_text())
    x0 = float(start["x0"])
    y = np.asarray(rec["final_y"], dtype=float)
    ev = ch.chain_eval(sysm, legs, x0, y)
    assert ev is not None
    gate = ch.gate_eval(sysm, ev)
    max_res = float(np.max(np.abs(ev.residual)))
    flybys = []
    for blk in ev.block_flybys:
        for f in blk:
            flybys.append(
                {
                    "body": f.body,
                    "t_s_past_jd2440000": float(f.t_s),
                    "vinf_in": np.asarray(f.vinf_in).tolist(),
                    "vinf_out": np.asarray(f.vinf_out).tolist(),
                    "vinf_kms": float(f.vinf_kms),
                    "turn_deg": float(f.turn_deg),
                }
            )
    # The first leg's departure (x0) is a Ganymede encounter with no inbound leg in the open chain.
    seg0 = ev.segments[0]
    _, w0 = sysm.state(seg0[0], seg0[1])
    out = {
        "tool_hash": TOOL_HASH,
        "epoch_te_s_past_jd2440000": te,
        "x0_s_past_jd2440000": x0,
        "first_departure": {
            "body": seg0[0],
            "t_s_past_jd2440000": float(seg0[1]),
            "vinf_out": (np.asarray(seg0[2]) - w0).tolist(),
        },
        "max_residual": max_res,
        "gate": gate,
        "control_pass": bool(
            max_res < 1e-6
            and abs(float(gate["worst_ratio"]) - 0.804) < 1e-3
            and abs(float(gate["min_required_alt_km"]) - 982.0) < 1.0
        ),
        "n_flybys": len(flybys),
        "flybys": flybys,
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "seed_chain.json").write_text(json.dumps(out, indent=1, default=float))
    print(json.dumps({k: v for k, v in out.items() if k != "flybys"}, default=float))


if __name__ == "__main__":
    main()
