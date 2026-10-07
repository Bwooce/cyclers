"""#1043: one-cycle EGGIE in continuous gravity (pre-registration:
docs/notes/2026-10-08-1043-eggie-one-cycle-continuous.md). A thin wrapper of
``scripts/run_1041_eggie.py`` with n = 1 and its own data directory.

    uv run python scripts/run_1043_eggie.py --stage pc|sigma|down|ias15

Checkpoints in data/1043_eggie/; live log in data/1043_eggie/live/ (gitignored).
"""

from __future__ import annotations

import argparse
import importlib.util
import sys
from pathlib import Path
from typing import Any

from cyclerfinder.data.method_capability import MethodCapability
from cyclerfinder.data.preflight import preflight_search

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "1043_eggie"


def _load(name: str, path: Path) -> Any:
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


Q = _load("run_1041_eggie", ROOT / "scripts" / "run_1041_eggie.py")
Q.OUT = OUT
Q.R.OUT = OUT
# #1043's criteria (note sec. 3) have no SOI bound: a 0.9-deg Io turn needs r_p beyond the SOI
# even in the patched conic (amendment 2). Floors and the unscheduled-pass check still apply.
Q.R.SOI_LIMIT = float("inf")


def _js(c: Any, z: Any, sg: float) -> Any:
    out = []
    for f in (1.0 + 1e-6, 1.0 - 1e-6):
        c.set_sigma(sg * f)
        out.append(Q.R.residual_and_jac(c, z, False)[0])
    c.set_sigma(sg)
    return (out[0] - out[1]) / (2e-6 * sg)


def stage_monitor(n_points: int) -> None:
    """Fold check (pre-registration sec. 4): pass the turning point with sigma as an unknown and
    the fastest-changing scaled state component stepped and fixed by an extra row."""
    import json

    import numpy as np

    R = Q.R  # noqa: N806
    c = Q.build(1)
    path = OUT / "monitor_n1.json"
    if path.exists():
        rec = json.loads(path.read_text())
    else:
        pts = json.loads((OUT / "sigma_n1.json").read_text())["points"]
        rec = {"points": [pts[-2], pts[-1]]}
    for _ in range(n_points):
        a, b = rec["points"][-2], rec["points"][-1]
        za, zb = np.asarray(a["z"]), np.asarray(b["z"])
        c.set_sigma(b["sigma"])
        _, jz, _ = R.residual_and_jac(c, zb)
        d = np.linalg.norm(jz, axis=0)
        d[d == 0.0] = 1.0
        j = int(np.argmax(np.abs((zb - za) * d)))
        h = float(zb[j] - za[j])
        ok = False
        for _half in range(6):
            target = zb[j] + h
            y = np.concatenate(
                [
                    zb + (zb - za) * (h / float(zb[j] - za[j])),
                    [b["sigma"] + (b["sigma"] - a["sigma"]) * (h / float(zb[j] - za[j]))],
                ]
            )
            for _it in range(20):
                sg = float(y[-1])
                c.set_sigma(sg)
                try:
                    r, jzz, info = R.residual_and_jac(c, y[:-1])
                except RuntimeError:
                    break
                if R.converged(info) and abs(y[j] - target) < 1e-9 * max(1.0, abs(target)):
                    ok = True
                    break
                js = _js(c, y[:-1], sg)
                row = np.zeros(len(y))
                row[j] = 1.0
                jac = np.vstack([np.hstack([jzz, js[:, None]]), row[None, :]])
                f = np.concatenate([r, [y[j] - target]])
                dc = np.linalg.norm(jac, axis=0)
                dc[dc == 0.0] = 1.0
                step = np.linalg.lstsq(jac / dc, -f, rcond=None)[0] / dc
                f0 = float(np.linalg.norm(f))
                for alpha in [0.5**i for i in range(8)]:
                    yt = y + alpha * step
                    c.set_sigma(float(yt[-1]))
                    rt = R.residual_and_jac(c, yt[:-1], False)[0]
                    if float(np.linalg.norm(np.concatenate([rt, [yt[j] - target]]))) < f0:
                        y = yt
                        break
                else:
                    break
            if ok:
                break
            h *= 0.5
        if not ok:
            R._log("monitor: six halvings without convergence: NUMERICAL stop")
            break
        sg = float(y[-1])
        c.set_sigma(sg)
        desc = R.describe(c, y[:-1], sg)
        rec["points"].append({"sigma": sg, "z": y[:-1].tolist(), **desc, "monitor_index": j})
        path.write_text(json.dumps(rec))
        R._log(f"monitor sigma={sg:.6f} vinf {[round(n['vinf_kms'], 4) for n in desc['nodes']]}")


def stage_two(sigmas: list[float]) -> None:
    import json

    import numpy as np

    R = Q.R  # noqa: N806
    c = Q.build(1)
    lower = json.loads((OUT / "sigma_n1.json").read_text())["points"]
    upper = json.loads((OUT / "monitor_n1.json").read_text())["points"][2:]
    out = {}
    for sg in sigmas:
        res = {}
        for lab, pts in (("lower", lower), ("upper", upper)):
            q = min(pts, key=lambda t: abs(t["sigma"] - sg))
            c.set_sigma(sg)
            z, info = R.newton(c, np.asarray(q["z"]), 30, f"two {lab} {sg}")
            res[lab] = R.describe(c, z, sg) if R.converged(info) else None
        row: dict[str, Any] = dict(res)
        if res["lower"] and res["upper"]:
            vl = [n["vinf_kms"] for n in res["lower"]["nodes"]]
            vu = [n["vinf_kms"] for n in res["upper"]["nodes"]]
            rl = [n["rp_km"] for n in res["lower"]["nodes"]]
            ru = [n["rp_km"] for n in res["upper"]["nodes"]]
            row["max_d_vinf"] = max(abs(x - y) for x, y in zip(vl, vu, strict=True))
            row["max_d_rp_km"] = max(abs(x - y) for x, y in zip(rl, ru, strict=True))
            row["fold_supported"] = row["max_d_vinf"] > 0.02 and row["max_d_rp_km"] > 100.0 * sg
        out[str(sg)] = row
        R._log(
            f"two sigma={sg}: {row.get('max_d_vinf')} {row.get('max_d_rp_km')} "
            f"{row.get('fold_supported')}"
        )
    (OUT / "two_n1.json").write_text(json.dumps(out, indent=1))


def main() -> None:
    preflight_search(
        task_no=1043,
        region_id="eggie-one-cycle-continuous-paper-model",
        method=MethodCapability(
            genome="EGGIE one-cycle open chain (one structure), the paper's ideal Galilean model",
            corrector="forward-backward shooting, DOP853 + STM; sigma continuation; pinned ends",
            capability_tags=frozenset({"ballistic", "n-body"}),
            git_sha="working-tree",
        ),
        script_path=Path(__file__),
        n_points=1,
    )
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--stage", choices=["pc", "sigma", "down", "ias15", "monitor", "two"], required=True
    )
    ap.add_argument("--points", type=int, default=10)
    ap.add_argument("--at", default="0.5,0.49")
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    if args.stage == "pc":
        Q.stage_pc(1)
    elif args.stage == "sigma":
        Q.R.stage_sigma(1)
    elif args.stage == "down":
        Q.R.stage_down(1)
    elif args.stage == "monitor":
        stage_monitor(args.points)
    elif args.stage == "two":
        stage_two([float(v) for v in args.at.split(",")])
    else:
        Q.stage_ias15(1)


if __name__ == "__main__":
    main()
