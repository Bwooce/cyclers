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
#: Overridable by `run --root/--timeout` (the lead-launched DE440 shoot, note sec. 7).
ROOT = [OUT]
TIMEOUT = [TIMEOUT_S]
EPOCHS = range(5)

#: Route ladder (sec. 2): name -> extra chain-tool flags.
ROUTES = {
    "S1": ["--real", "mean", "--direct"],
    "S2": ["--real", "mean"],
    "S3": ["--real", "mean", "--grow-chain", "--grow-start", "2"],
    "D1": ["--real", "de440", "--shoot-jac", "sparse"],
    "D2": ["--real", "de440", "--direct", "--shoot-jac", "sparse"],
    # amendment A1 (note sec. 5a): blend to lambda 1, no shoot (the 6.18/6.24 reading)
    "DB1": ["--real", "de440", "--shoot-restarts", "0"],
    "DB2": ["--real", "de440", "--direct", "--shoot-restarts", "0"],
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
    # #1057 step 1 (note sec. 9): ev-A, 5 cycles, the C2 key and dates, as a DE440 shoot control
    out["ev-A"] = (CONTROLS["C2"][0], CONTROLS["C2"][1], CONTROLS["C2"][2])
    return out


def run_route(
    tag: str, key: str, xd: str, n: int, route: str, ie: int, extra: list[str] | None = None
) -> dict[str, Any]:
    """One chain-tool call (one epoch); then the checker if the tool reports a rung pass."""
    d = ROOT[0] / tag / route / f"e{ie}"
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
        # live log (gitignored under de440_shoot/), copied to the static tool.log at the end
        with (d / "tool.live.log").open("w") as live:
            p = subprocess.run(
                cmd, stdout=live, stderr=subprocess.STDOUT, timeout=TIMEOUT[0], cwd=REPO
            )
        (d / "tool.log").write_text((d / "tool.live.log").read_text())
    except subprocess.TimeoutExpired:
        (d / "tool.log").write_text(
            f"TIMEOUT after {TIMEOUT[0]} s\n" + (d / "tool.live.log").read_text()
        )
        log(f"{tag} {route} e{ie}: NOT RUN (timeout {TIMEOUT[0]} s)")
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


# ---------------------------------------------------------------------------
# DE440 under amendment A1 (note sec. 5a): blend to lambda 1, no shoot
# ---------------------------------------------------------------------------


def _deb_eval(tag: str, key: str, n: int, route: str, ie: int) -> dict[str, Any] | None:
    """Reached / minimax gate / closure dV (Earth, Venus) / deficit for one DB route and epoch."""
    d = OUT / tag / route / f"e{ie}"
    p = d / "realeph_chain.json"
    if not p.exists():
        return None
    ep = json.loads(p.read_text())[0]
    steps = [s for s in ep.get("steps", []) if s.get("converged")]
    reached = bool(steps) and steps[-1]["lambda"] >= 1.0
    out: dict[str, Any] = {
        "route": route,
        "epoch": ie,
        "epoch_jd": ep["epoch_jd"],
        "reached": reached,
    }
    out["last_converged_lambda"] = ep.get("last_converged_lambda")
    if not reached:
        return out
    out["minimax_gate"] = {
        "status": steps[-1].get("status"),
        "worst_ratio": steps[-1].get("worst_ratio"),
    }
    cj = d / "closure_dv.json"
    if not cj.exists():
        subprocess.run(
            ["uv", "run", "python", str(CLOSURE_TOOL), "--cell", "ev", "--key", key,
             "--n-cycles", str(n), "--real", "de440", "--chain", str(p), "--out", str(cj)],
            capture_output=True, text=True, timeout=TIMEOUT_S, cwd=REPO, check=True,
        )  # fmt: skip
    rows = json.loads(cj.read_text())[0]["legs"]
    out["closure_dv_earth_mps"] = sum(r["dv_mid_ms"] for r in rows if r["body"] == "E")
    out["closure_dv_venus_mps"] = sum(r["dv_mid_ms"] for r in rows if r["body"] == "V")
    out["max_miss_earth_km"] = max((r["miss_km"] for r in rows if r["body"] == "E"), default=0.0)
    out["max_miss_venus_km"] = max((r["miss_km"] for r in rows if r["body"] == "V"), default=0.0)
    dvb = _load("dv_band_942", REPO / "scripts" / "dv_band_942.py")
    # lambda = 1 blend (real states, model periods for the fixed-leg timing), as in the
    # closure-dV script; DE440Planets alone has no period
    import numpy as np

    circ = dvb.CHAIN.ENUM.cell_system("ev")[0]
    sys1 = dvb.CHAIN.Blend(circ, dvb.CHAIN.real_ephemeris("ev", "de440"), np.eye(3), 0.0, 0.0, 1.0)
    pr = dvb.price_chain(sys1, key, n, ep["x0_days"], ep["final_y_dates"])
    out["deficit_beyond_bend_mps"] = pr["deficit_beyond_bend_mps"]
    out["n_priced_flybys"] = pr["n_flybys"]
    s7 = 7.0 / n
    tot = out["closure_dv_earth_mps"] + out["closure_dv_venus_mps"] + pr["deficit_beyond_bend_mps"]
    out["per7_total_mps"] = tot * s7
    out["band_total_upper_bound"] = dvb.band(tot * s7)
    ven = out["closure_dv_venus_mps"] + pr["deficit_beyond_bend_mps"]
    out["per7_venus_plus_deficit_mps"] = ven * s7
    out["band_venus_plus_deficit"] = dvb.band(ven * s7)
    return out


def ladder_deb(tag: str, key: str, xd: str, n: int, epochs: list[int]) -> list[dict]:
    _KEY[0] = key
    res = []
    for ie in epochs:
        run_route(tag, key, xd, n, "DB1", ie)
        ev = _deb_eval(tag, key, n, "DB1", ie)
        if ev is None or not ev["reached"]:
            run_route(tag, key, xd, n, "DB2", ie)
            ev2 = _deb_eval(tag, key, n, "DB2", ie)
            ev = ev2 if ev2 is not None else ev
        res.append(ev or {"epoch": ie, "reached": False, "route": None, "not_run": True})
        e = res[-1]
        log(
            f"{tag} deb e{ie}: reached={e.get('reached')} route={e.get('route')} "
            f"gate={(e.get('minimax_gate') or {}).get('status')} "
            f"worst={(e.get('minimax_gate') or {}).get('worst_ratio')} "
            f"closure E/V={e.get('closure_dv_earth_mps')}/{e.get('closure_dv_venus_mps')} m/s "
            f"per7={e.get('per7_total_mps')}"
        )
    (OUT / tag / "de440_blend.json").write_text(json.dumps(res, indent=1, default=float))
    return res


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
    ap.add_argument("model", choices=["std", "de", "deb"])
    ap.add_argument("--epochs", default="0,1,2,3,4")
    ap.add_argument("--root", default=None, help="subdirectory of data/1025_ev_realeph")
    ap.add_argument("--timeout", type=int, default=TIMEOUT_S, help="seconds per epoch and route")
    a = ap.parse_args(argv)
    if a.root:
        ROOT[0] = OUT / a.root
    TIMEOUT[0] = a.timeout
    key, xd, n = members()[a.member]
    if a.model == "deb":
        ladder_deb(a.member, key, xd, n, [int(v) for v in a.epochs.split(",")])
        return
    res = ladder(a.member, key, xd, n, a.model, [int(v) for v in a.epochs.split(",")])
    for r in res:
        log(f"{a.member} {a.model} e{r['epoch']}: {r['outcome']}")


def _chain_dates(route: str, ep: dict[str, Any], n: int) -> tuple[float, list[float]] | None:
    """(x0_days, y) of the full-length lambda = 1 chain a route produced, else None."""
    steps = [s for s in ep.get("steps", []) if s.get("converged")]
    if route == "S3":
        grow = ep.get("grow") or []
        if grow and grow[-1]["k"] == n and grow[-1].get("converged"):
            return ep["x0_days"], ep["final_y_dates"]
        return None
    if steps and steps[-1]["lambda"] >= 1.0:
        return ep["x0_days"], ep["final_y_dates"]
    return None


def standish_band(tag: str) -> dict[str, Any]:
    """Sec. 4: per epoch, price the passing chain (else the lowest-worst converged chain)."""
    dvb = _load("dv_band_942", REPO / "scripts" / "dv_band_942.py")
    real = dvb.CHAIN.real_ephemeris("ev", "mean")
    key, _, n = members()[tag]
    per_epoch = []
    for ie in EPOCHS:
        cands = []
        for route in ("S1", "S2", "S3"):
            d = OUT / tag / route / f"e{ie}"
            p = d / "realeph_chain.json"
            if not p.exists():
                continue
            ep = json.loads(p.read_text())[0]
            cl = classify(ep, route, n, d)
            ch = _chain_dates(route, ep, n)
            if ch is not None:
                cands.append((cl["outcome"] != "PASS", cl.get("worst_ratio") or 9.0, route, ep, ch))
        if not cands:
            per_epoch.append({"epoch": ie, "per7_beyond_bend_mps": None, "band": "n/a"})
            continue
        _, _, route, ep, (x0, y) = min(cands, key=lambda c: (c[0], c[1]))
        pr = dvb.price_chain(real, key, n, x0, y)
        r = dvb.row(tag, "standish", {"epoch_jd": ep["epoch_jd"]}, n, pr, 0.0, route)
        per_epoch.append(
            {
                "epoch": ie,
                "route": route,
                "epoch_jd": ep["epoch_jd"],
                "worst_ratio": pr["worst_ratio"],
                "n_priced_flybys": pr["n_flybys"],
                "per7_beyond_bend_mps": r["per_7_cycles_mps"]["beyond_bend"],
                "per7_periapsis_mps": r["per_7_cycles_mps"]["periapsis"],
                "band": r["band_beyond_bend"],
            }
        )
    vals = [e["per7_beyond_bend_mps"] for e in per_epoch if e["per7_beyond_bend_mps"] is not None]
    best = min(
        (e for e in per_epoch if e["per7_beyond_bend_mps"] is not None),
        key=lambda e: e["per7_beyond_bend_mps"],
        default=None,
    )
    return {
        "per_epoch": per_epoch,
        "best": best,
        "range_per7_mps": [min(vals), max(vals)] if vals else None,
    }


def cmd_summary(argv: list[str]) -> None:
    out: dict[str, Any] = {}
    for tag, (key, _, n) in members().items():
        if not (OUT / tag).exists():
            continue
        rec: dict[str, Any] = {"key": key, "cycles": n}
        outs = []
        for ie in EPOCHS:
            tried = []
            for route in ("S1", "S2", "S3"):
                d = OUT / tag / route / f"e{ie}"
                if (d / "realeph_chain.json").exists():
                    _KEY[0] = key
                    tried.append(
                        classify(json.loads((d / "realeph_chain.json").read_text())[0], route, n, d)
                    )
            outs.append(
                {
                    "epoch": ie,
                    "outcome": epoch_outcome(tried) if tried else "NOT RUN",
                    "routes": tried,
                }
            )
        rec["standish"] = outs
        rec["standish_counts"] = {
            o: sum(e["outcome"] == o for e in outs)
            for o in ("PASS", "FAIL", "NOT REACHED", "NOT RUN")
        }
        rec["standish_band"] = standish_band(tag)
        dbj = OUT / tag / "de440_blend.json"
        if dbj.exists():
            db = json.loads(dbj.read_text())
            reached = [e for e in db if e.get("reached")]
            best_db = min(reached, key=lambda e: e["per7_total_mps"], default=None)
            rec["de440_blend"] = {
                "per_epoch": db,
                "reached": len(reached),
                "minimax_gate_pass": sum(
                    (e.get("minimax_gate") or {}).get("status") == "pass" for e in reached
                ),
                "best": best_db,
                "range_per7_total_mps": [
                    min(e["per7_total_mps"] for e in reached),
                    max(e["per7_total_mps"] for e in reached),
                ]
                if reached
                else None,
            }
        out[tag] = rec
        c = rec["standish_counts"]
        b = rec["standish_band"]["best"]
        best_s = (
            "n/a"
            if b is None
            else "{:.2f} m/s ({}, e{})".format(b["per7_beyond_bend_mps"], b["band"], b["epoch"])
        )
        print(
            f"{tag}: Standish {c['PASS']} pass / {c['FAIL']} fail / {c['NOT REACHED']} not reached"
            f" / {c['NOT RUN']} not run; best per-7 {best_s}"
            f"; range {rec['standish_band']['range_per7_mps']}",
            flush=True,
        )
    (OUT / "summary.json").write_text(json.dumps(out, indent=1, default=float))


#: The lead-launched DE440 shoot (note sec. 7): the Standish >= 4/5 members.
SHOOT_MEMBERS = [
    "ev4-1", "ev4-4", "ev4-6", "ev4-7", "ev5-4", "ev5-5", "ev5-6", "ev5-7", "ev5-8", "ev5-9",
]  # fmt: skip


def cmd_shootcheck(argv: list[str]) -> None:
    """Outcome per member and epoch of the DE440 shoot under data/1025_ev_realeph/de440_shoot/;
    an epoch with no result file is "NOT RUN" (resume by re-running the launch command)."""
    ap = argparse.ArgumentParser(prog="shootcheck")
    ap.add_argument("--root", default="de440_shoot")
    ap.add_argument("--members", default=",".join(SHOOT_MEMBERS))
    a = ap.parse_args(argv)
    ROOT[0] = OUT / a.root
    allm = members()
    res: dict[str, Any] = {}
    for tag in a.members.split(","):
        key, _, n = allm[tag]
        _KEY[0] = key
        outs = []
        for ie in EPOCHS:
            tried = []
            for route in ("D1", "D2"):
                d = ROOT[0] / tag / route / f"e{ie}"
                if (d / "realeph_chain.json").exists():
                    ep = json.loads((d / "realeph_chain.json").read_text())[0]
                    tried.append(classify(ep, route, n, d))
            outs.append({"epoch": ie, "outcome": epoch_outcome(tried) if tried else "NOT RUN",
                         "routes": tried})  # fmt: skip
        res[tag] = outs
        print(tag, " ".join(f"e{o['epoch']}:{o['outcome']}" for o in outs), flush=True)
    ROOT[0].mkdir(parents=True, exist_ok=True)
    (ROOT[0] / "summary.json").write_text(json.dumps(res, indent=1, default=float))


def main() -> None:
    cmds = {
        "control": cmd_control,
        "run": cmd_run,
        "summary": cmd_summary,
        "shootcheck": cmd_shootcheck,
    }
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
