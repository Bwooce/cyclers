"""#1033: descending continuation (sigma 1 -> 0.005) from the full-mass G-C-C-G orbit that gc-2's
seed converges to at sigma = 1.
Pre-registration: docs/notes/2026-10-08-1033-gc2-fullmass-descent.md.
Reuses the #1025 driver's corrector, monitor step and IAS15 unchanged.

    uv run python scripts/run_1033_gc2_descent.py

Checkpoint data/1033_gc2_descent/descent.json; live log in .../live/ (gitignored).
"""

from __future__ import annotations

import importlib.util
import itertools
import json
import math
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np

from cyclerfinder.data.method_capability import MethodCapability
from cyclerfinder.data.preflight import preflight_search

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "1033_gc2_descent"
PC = (3.617, 3.039, 3.039, 3.617)  # gc-2 patched conic (G, C, C, G)


def _load(name: str, path: Path) -> Any:
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


B = _load("run_1025_sigma_batch", ROOT / "scripts" / "run_1025_sigma_batch.py")
B.OUT = OUT


def _gc2_branches() -> dict[str, list[dict[str, Any]]]:
    r = json.loads((ROOT / "data" / "1025_sigma" / "gc-2.json").read_text())
    return {"lower": r["points"], "upper": r["upper"][2:]}


def _annotate(q: dict[str, Any], br: dict[str, list[dict[str, Any]]]) -> None:
    vs = [n["vinf_kms"] for n in q["nodes"]]
    q["max_dvinf_vs_gc2_pc"] = max(abs(v - p) for v, p in zip(vs, PC, strict=True))
    for lab, pts in br.items():
        near = [t for t in pts if abs(t["sigma"] - q["sigma"]) <= 0.02 * q["sigma"]]
        if near:
            t = min(near, key=lambda t: abs(t["sigma"] - q["sigma"]))
            q[f"max_dvinf_vs_gc2_{lab}"] = max(
                abs(v - n["vinf_kms"]) for v, n in zip(vs, t["nodes"], strict=True)
            )


def main() -> None:
    preflight_search(
        task_no=1033,
        region_id="gc2-fullmass-gccg-descent",
        method=MethodCapability(
            genome="full-mass symmetric G-C-C-G orbit from gc-2's seed, R-S circular, joint sigma",
            corrector="forward-backward shooting, DOP853 + STM; descending sigma; monitor; IAS15",
            capability_tags=frozenset({"ballistic", "n-body"}),
            git_sha="working-tree",
        ),
        script_path=Path(__file__),
        n_points=1,
    )
    OUT.mkdir(parents=True, exist_ok=True)
    deadline = time.monotonic() + 440.0
    path = OUT / "descent.json"
    br = _gc2_branches()
    p = B.build("gc-2")
    rec: dict[str, Any] = json.loads(path.read_text()) if path.exists() else {"member": "gc2-desc"}

    def save() -> None:
        path.write_text(json.dumps(rec))

    if "points" not in rec:
        p.set_sigma(1.0)
        z, info = B.solve(p, p.seed, 25, 0, False)
        if not B.converged(info):
            rec["outcome"] = "start failed at sigma = 1"
            save()
            return
        q = B._point(p, z, 1.0, info)
        _annotate(q, br)
        rec.update(points=[q], fac=1.2, stage="down")
        rec["ias15_sigma1"] = {}
        B.ias15(p, z, rec["ias15_sigma1"], deadline)
        B._log(f"gc2-desc: sigma=1 vinf {[round(n['vinf_kms'], 4) for n in q['nodes']]}")
        save()
    while time.monotonic() < deadline and rec["stage"] in ("down", "fold"):
        pts = rec["points"]
        if rec["stage"] == "fold":
            up = rec["upper"]
            if up[-1]["sigma"] > rec["sigma_stop"] * (1 + 1e-3):
                rec.update(
                    stage="done", outcome="(iii) fold on the descent", sigma_f=rec["sigma_stop"]
                )
            elif up[-1]["sigma"] < rec["sigma_stop"] * (1 - 1e-3):
                pts.append(up[-1])
                rec.pop("upper")
                rec.pop("mon_ratio", None)
                rec["stage"], rec["fac"] = "down", 1.01
                B._log("gc2-desc: monitor passed the stop with sigma falling; resume")
            elif not B.monitor_step(p, rec):
                rec.update(stage="done", outcome="numerical stop (monitor)")
            elif len(up) > 32:
                rec.update(stage="done", outcome="numerical stop (30 monitor points)")
            save()
            continue
        s1, z1 = pts[-1]["sigma"], np.asarray(pts[-1]["z"])
        if s1 <= 0.005 + 1e-12:
            rec["stage"] = "end"
            break
        fac = rec["fac"]
        sg = max(0.005, s1 / fac)
        if (s1 - sg) / s1 < 1e-5:
            rec.update(stage="fold", sigma_stop=s1, upper=[pts[-2], pts[-1]])
            B._log(f"gc2-desc: step below 1e-5 at sigma {s1:.6f}; monitor")
            save()
            continue
        if len(pts) >= 2:
            s0, z0 = pts[-2]["sigma"], np.asarray(pts[-2]["z"])
            zs = z1 + (z1 - z0) * math.log(sg / s1) / math.log(s1 / s0)
        else:
            zs = B.scale_offsets(p, z1, sg / s1)
        p.set_sigma(sg)
        z, info = B.solve(p, zs, 25, 15, True)
        if B.converged(info) or B.noise_ok(info):
            q = B._point(p, z, sg, info)
            if any(n["impact"] for n in q["nodes"]):
                rec.update(
                    stage="done", outcome="(iii) impact on the descent", sigma_i=sg, impact_point=q
                )
            else:
                _annotate(q, br)
                pts.append(q)
                rec["fac"] = min(1.2, 1.0 + 2.0 * (fac - 1.0))
                vs = [round(n["vinf_kms"], 4) for n in q["nodes"]]
                rps = [round(n["rp_over_sigma"]) for n in q["nodes"]]
                B._log(f"gc2-desc: sigma={sg:.6g} vinf {vs} rp/s {rps}")
        else:
            rec["fac"] = 1.0 + 0.5 * (fac - 1.0)
            B._log(f"gc2-desc: sigma={sg:.6g} failed; factor {rec['fac']:.6g}")
        save()
    if rec.get("stage") == "end":
        pts = rec["points"]
        sh = [q["max_dvinf_vs_gc2_pc"] for q in pts[-3:]]
        mono = all(b <= a + 1e-4 for a, b in itertools.pairwise(sh))
        q = pts[-1]
        p.set_sigma(q["sigma"])
        rec["ias15_end"] = {}
        if B.ias15(p, np.asarray(q["z"]), rec["ias15_end"], deadline + 120.0):
            ok = q["max_dvinf_vs_gc2_pc"] < 0.005 and mono
            rec["outcome"] = "(i) reaches gc-2's generating orbit" if ok else "(ii) different limit"
            rec["stage"] = "done"
    save()
    B._log(f"gc2-desc: stage {rec.get('stage')} outcome {rec.get('outcome')}")


if __name__ == "__main__":
    main()
