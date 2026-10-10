"""#1025: rung (d) and #415 dV bands for the 22 clean #973 Earth-Venus members.

Pre-registration: docs/notes/2026-10-11-1025-ev-realeph-rung-d.md (route ladder sec. 2, controls
sec. 3, dV band sec. 4). The shared tools are called as subprocesses, one epoch per call under a
timeout; nothing in them is changed:
  scripts/run_942_realeph_chain.py   (chain: Standish ramp / DE440 blend+shoot)
  scripts/check_942_realeph_chain.py (DOP853 re-fly of a passing chain)
  scripts/analyse_942_fullrev_closure_dv.py (DE440 full-rev closure dV, 6.24a)
  scripts/dv_band_942.py::price_chain, row (turn-deficit prices, 6.57)

Subcommands:
  control C1|C2|C3            the route controls (sec. 3)
  run MEMBER std|de [--epochs 0,1,..]  the ladder for one member and model (resumable)
  summary                     verdicts and bands for every member run so far

Output: data/1025_ev_realeph/<member>/<route>/e<i>/ (realeph_chain.json, tool.log, check.log) and
data/1025_ev_realeph/runlog.txt (static, appended).

Usage: uv run python scripts/ev_realeph_1025.py run ev4-0 std
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import re
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

REPO = Path(__file__).resolve().parents[1]
OUT = REPO / "data" / "1025_ev_realeph"
CHAIN_TOOL = REPO / "scripts" / "run_942_realeph_chain.py"
CHECK_TOOL = REPO / "scripts" / "check_942_realeph_chain.py"
CLOSURE_TOOL = REPO / "scripts" / "analyse_942_fullrev_closure_dv.py"
TIMEOUT_S = 420
EPOCHS = range(5)

#: Route ladder (sec. 2): name -> extra chain-tool flags.
ROUTES = {
    "S1": ["--real", "mean", "--direct"],
    "S2": ["--real", "mean"],
    "S3": ["--real", "mean", "--grow-chain", "--grow-start", "2"],
    "D1": ["--real", "de440", "--shoot-jac", "sparse"],
    "D2": ["--real", "de440", "--direct", "--shoot-jac", "sparse"],
}

#: Route controls (sec. 3): name -> (key, x_days, cycles, route, epochs, recorded passes).
CONTROLS = {
    "C1": (
        "k2|RE/1:1|LE>V/0s|RV/1:1|RV/1:1|LV>E/0s",
        "474.59719334719347,1100.6183800415804",
        5,
        "S1",
        list(EPOCHS),
        [1, 2],
    ),
    "C2": (
        "k2|LE>V/0s|RV/1:1|LV>V/1h|LV>E/0s",
        "583.944386694386,984.2215547169349,1576.2577920606109",
        5,
        "S1",
        list(EPOCHS),
        [0, 1, 2, 3, 4],
    ),
    "C3": (
        "k3|RE/1:1|LE>V/0s|LV>V/1l|RV/3:2|LV>E/0s",
        "547.5969022070199,597.6931769759217,1778.7965570094584",
        4,
        "S2",
        [0],
        [0],
    ),
}


def log(msg: str) -> None:
    line = f"{time.strftime('%Y-%m-%dT%H:%M:%S%z')} {msg}"
    print(line, flush=True)
    OUT.mkdir(parents=True, exist_ok=True)
    with (OUT / "runlog.txt").open("a") as f:
        f.write(line + "\n")


def members() -> dict[str, tuple[str, str, int]]:
    """name -> (key, x_days csv, cycles)."""
    out = {}
    for k, idx, n in ((4, [0, 1, 3, 4, 5, 6, 7, 8, 9, 10, 11], 3), (5, list(range(11)), 2)):
        g = json.loads((REPO / "data" / "973_ev" / f"k{k}_gauntlet.json").read_text())
        for i in idx:
            c = g["candidates"][i]
            assert c["screen_status"] == "pass", (k, i)
            out[f"ev{k}-{i}"] = (c["key"], ",".join(repr(float(v)) for v in c["x_days"]), n)
    return out


def run_route(
    tag: str, key: str, xd: str, n: int, route: str, ie: int, extra: list[str] | None = None
) -> dict[str, Any]:
    """One chain-tool call (one epoch); then the checker if the tool reports a rung pass."""
    d = OUT / tag / route / f"e{ie}"
    res_path = d / "realeph_chain.json"
    if res_path.exists():
        return classify(json.loads(res_path.read_text())[0], route, n, d)
    d.mkdir(parents=True, exist_ok=True)
    cmd = [
        "uv", "run", "python", str(CHAIN_TOOL), "--cell", "ev", "--key", key, "--x-days", xd,
        "--n-cycles", str(n), "--only-epoch", str(ie), "--out", str(d), *ROUTES[route],
        *(extra or []),
    ]  # fmt: skip
    t0 = time.time()
    try:
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=TIMEOUT_S, cwd=REPO)
        (d / "tool.log").write_text(p.stdout + p.stderr)
    except subprocess.TimeoutExpired as exc:
        (d / "tool.log").write_text(f"TIMEOUT after {TIMEOUT_S} s\n{exc.stdout or ''}")
        log(f"{tag} {route} e{ie}: NOT RUN (timeout {TIMEOUT_S} s)")
        return {"route": route, "epoch": ie, "outcome": "NOT RUN", "reason": "timeout"}
    if not res_path.exists():
        log(f"{tag} {route} e{ie}: NOT RUN (no output; exit {p.returncode})")
        return {"route": route, "epoch": ie, "outcome": "NOT RUN", "reason": f"exit {p.returncode}"}
    r = classify(json.loads(res_path.read_text())[0], route, n, d)
    log(
        f"{tag} {route} e{ie} JD {r.get('epoch_jd', 0):.1f}: {r['outcome']} "
        f"worst={r.get('worst_ratio')} refly={r.get('refly_miss_km')} ({time.time() - t0:.0f}s)"
    )
    return r


def _check(d: Path, key: str, n: int, real: str) -> tuple[float | None, float | None]:
    cl = d / "check.log"
    if not cl.exists():
        p = subprocess.run(
            ["uv", "run", "python", str(CHECK_TOOL), "--cell", "ev", "--key", key,
             "--chain", str(d / "realeph_chain.json"), "--n-cycles", str(n), "--real", real],
            capture_output=True, text=True, timeout=TIMEOUT_S, cwd=REPO,
        )  # fmt: skip
        cl.write_text(p.stdout + p.stderr)
    txt = cl.read_text()
    m = re.search(r"max arrival miss ([0-9.eE+-]+) km", txt)
    v = re.search(r"max V_inf vector error ([0-9.eE+-]+) km/s", txt)
    return (float(m.group(1)) if m else None, float(v.group(1)) if v else None)


def classify(ep: dict[str, Any], route: str, n: int, d: Path) -> dict[str, Any]:
    """PASS / FAIL / NOT REACHED for one route at one epoch (sec. 2)."""
    steps = [s for s in ep.get("steps", []) if s.get("converged")]
    out: dict[str, Any] = {"route": route, "epoch": int(d.name[1:]), "epoch_jd": ep["epoch_jd"]}
    full_len = False
    worst = status = None
    if route == "S3":
        grow = ep.get("grow") or []
        if grow and grow[-1]["k"] == n and grow[-1].get("converged"):
            full_len, worst, status = True, grow[-1].get("worst_ratio"), grow[-1].get("status")
    elif route.startswith("S"):
        if steps and steps[-1]["lambda"] >= 1.0:
            full_len, worst, status = True, steps[-1].get("worst_ratio"), steps[-1].get("status")
    else:  # DE440: the shoot closure decides
        sb = ep.get("shoot_best")
        if sb is not None:
            full_len, worst, status = True, sb.get("worst_ratio"), sb.get("status")
        out["blend_reached_lambda1"] = bool(steps) and steps[-1]["lambda"] >= 1.0
        if out["blend_reached_lambda1"]:
            out["minimax_gate_lambda1"] = {
                "status": steps[-1].get("status"),
                "worst_ratio": steps[-1].get("worst_ratio"),
            }
    out["last_converged_lambda"] = ep.get("last_converged_lambda")
    out["worst_ratio"] = worst
    if not full_len:
        out["outcome"] = "NOT REACHED"
        return out
    if status != "pass":
        out["outcome"] = "FAIL"
        out["gate_status"] = status
        return out
    miss, verr = _check(d, _KEY[0], n, "mean" if route.startswith("S") else "de440")
    out["refly_miss_km"], out["refly_vinf_err_kms"] = miss, verr
    ok = miss is not None and miss < 1.0
    if route.startswith("S"):
        ok = ok and verr is not None and verr < 1e-6
    out["outcome"] = "PASS" if ok else "FAIL"
    if not ok:
        out["gate_status"] = "pass, re-fly failed"
    return out


_KEY = [""]  # the current member's key, for the checker


def ladder(tag: str, key: str, xd: str, n: int, model: str, epochs: list[int]) -> list[dict]:
    _KEY[0] = key
    res = []
    for ie in epochs:
        if model == "std":
            r1 = run_route(tag, key, xd, n, "S1", ie)
            tried = [r1]
            if r1["outcome"] != "PASS":
                r2 = run_route(tag, key, xd, n, "S2", ie)
                tried.append(r2)
                if r2["outcome"] == "NOT REACHED" and n >= 3:
                    tried.append(run_route(tag, key, xd, n, "S3", ie))
        else:
            r1 = run_route(tag, key, xd, n, "D1", ie)
            tried = [r1]
            if r1["outcome"] == "NOT REACHED" and not r1.get("blend_reached_lambda1"):
                tried.append(run_route(tag, key, xd, n, "D2", ie))
        res.append({"epoch": ie, "routes": tried, "outcome": epoch_outcome(tried)})
    return res


def epoch_outcome(tried: list[dict[str, Any]]) -> str:
    outs = [t["outcome"] for t in tried]
    for o in ("PASS", "FAIL", "NOT REACHED"):
        if o in outs:
            return o
    return "NOT RUN"


def cmd_control(argv: list[str]) -> None:
    name = argv[0]
    key, xd, n, route, epochs, recorded = CONTROLS[name]
    _KEY[0] = key
    got = []
    for ie in epochs:
        r = run_route(f"control-{name}", key, xd, n, route, ie)
        got.append(r)
    passes = [r["epoch"] for r in got if r["outcome"] == "PASS"]
    same = passes == recorded
    log(f"CONTROL {name}: passes at epochs {passes}; recorded {recorded}; reproduced={same}")


def cmd_run(argv: list[str]) -> None:
    ap = argparse.ArgumentParser(prog="run")
    ap.add_argument("member")
    ap.add_argument("model", choices=["std", "de"])
    ap.add_argument("--epochs", default="0,1,2,3,4")
    a = ap.parse_args(argv)
    key, xd, n = members()[a.member]
    res = ladder(a.member, key, xd, n, a.model, [int(v) for v in a.epochs.split(",")])
    for r in res:
        log(f"{a.member} {a.model} e{r['epoch']}: {r['outcome']}")


def main() -> None:
    cmds = {"control": cmd_control, "run": cmd_run}
    if len(sys.argv) < 2 or sys.argv[1] not in cmds:
        raise SystemExit(f"usage: ev_realeph_1025.py {{{','.join(cmds)}}} ...")
    cmds[sys.argv[1]](sys.argv[2:])


def _load(name: str, path: Path) -> Any:
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


if __name__ == "__main__":
    main()
