"""#1025 rung-d note sec. 11: V2-ballistic judgement of the 10 Earth-Venus #973 rows under the
final spec sec. 14 wording (#1007 secs. 8.4b, 8.17).

The criteria code is imported from scripts/v2_samephase_1007.py (as of 5721cc8b) and not edited:
``final_record`` (criteria 1-3 from v2_chain_1007.judge, criterion 4 from final_c4) and
``final_controls`` (H&M orbit 1 positive, broken gc-1 negative). L is computed by
cyclerfinder.verify.geometric_repeat at 0.02 rev. The 10-cycle Standish chains come from
``scripts/ev_realeph_1025.py run MEMBER std --n-cycles 10 --root v2_10cyc`` (rung (d) ladder).

Usage:
  uv run python scripts/v2_ev10_1025.py controls
  uv run python scripts/v2_ev10_1025.py judge [MEMBER ...]
  uv run python scripts/v2_ev10_1025.py summary

Output: data/1025_ev_realeph/v2_10cyc/{controls,judge_<member>,summary}.json
"""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path
from typing import Any

from cyclerfinder.search.two_working_body import LambertLeg
from cyclerfinder.verify.geometric_repeat import GEOMETRIC_REPEAT_TOL_REV, repeat_for_pair

REPO = Path(__file__).resolve().parents[1]
OUT = REPO / "data" / "1025_ev_realeph" / "v2_10cyc"
N_CYCLES = 10
ROWS = ["ev4-1", "ev4-4", "ev4-6", "ev4-7", "ev5-4", "ev5-5", "ev5-6", "ev5-7", "ev5-8", "ev5-9"]
ROUTES = ("S1", "S2", "S3")


def _load(name: str, path: Path) -> Any:
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


SP = _load("v2_samephase_1007", REPO / "scripts" / "v2_samephase_1007.py")
EVR = _load("ev_realeph_1025", REPO / "scripts" / "ev_realeph_1025.py")
CHAIN = SP.CHAIN


def _write(name: str, obj: Any) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / name).write_text(json.dumps(obj, indent=1, default=float))
    print("wrote", (OUT / name).relative_to(REPO), flush=True)


def cmd_controls() -> None:
    _write("controls.json", SP.final_controls())


def judge_member(member: str, real: Any) -> dict[str, Any]:
    key, xd, _ = EVR.members()[member]
    k = int(key.split("|")[0][1:])
    big_l = repeat_for_pair(k, "E", "V")
    out: dict[str, Any] = {"member": member, "key": key, "k": k, "L": big_l,
                           "tol_rev": GEOMETRIC_REPEAT_TOL_REV, "cycles": N_CYCLES}  # fmt: skip
    judgeable = 2 * big_l <= N_CYCLES
    circ, a, b = CHAIN.ENUM.cell_system("ev")
    _, one = CHAIN.ENUM.parse_cycle_key(key, circ, a, b)
    nl = sum(isinstance(lg, LambertLeg) for lg in one.legs)
    x_tpl = [float(v) for v in xd.split(",")]
    epochs = []
    for ie in range(5):
        rec: dict[str, Any] = {"epoch": ie, "routes": []}
        for route in ROUTES:
            p = OUT / member / route / f"e{ie}" / "realeph_chain.json"
            if not p.exists():
                continue
            ep = json.loads(p.read_text())[0]
            rec["epoch_jd"] = ep["epoch_jd"]
            dates = EVR._chain_dates(route, ep, N_CYCLES)
            if dates is None:
                rec["routes"].append({"route": route, "reached": False})
                continue
            x0, y = dates
            lab = f"{member} {route} JD {ep['epoch_jd']:.1f}"
            r = SP.final_record(lab, "ev", key, x_tpl, N_CYCLES, "ramp", x0, list(y), real, nl,
                                big_l)  # fmt: skip
            rec["routes"].append({"route": route, "reached": True, "judged": r})
        judged = [t["judged"] for t in rec["routes"] if t.get("reached")]
        rec["pass"] = any(j["pass"] for j in judged)
        rec["outcome"] = "PASS" if rec["pass"] else ("FAIL" if judged else "NO 10-CYCLE CHAIN")
        epochs.append(rec)
    n_pass = sum(e["pass"] for e in epochs)
    if not judgeable:
        # sec. 11.7: the 10-cycle chains are a diagnostic of criteria 1-3 only
        for e in epochs:
            e["c123"] = any(
                t["judged"]["c1_continuity"] and t["judged"]["c2_gate"] and t["judged"]["c3_refly"]
                for t in e["routes"]
                if t.get("reached") and "c1_continuity" in t["judged"]
            )
            if e["outcome"] == "FAIL":
                e["outcome"] = "C1-3 PASS" if e["c123"] else "C1-3 FAIL"
        reason = f"n = {N_CYCLES} < 2L = {2 * big_l} (sec. 11.1); chains are a c1-3 diagnostic"
        out |= {"judgeable": False, "epochs": epochs, "passes": 0, "verdict": "NOT JUDGEABLE",
                "reason": reason, "c123_epochs": sum(e["c123"] for e in epochs)}  # fmt: skip
        print(f"{member}: NOT JUDGEABLE (L = {big_l}); c1-3 at {out['c123_epochs']}/5", flush=True)
        return out
    out |= {"judgeable": True, "epochs": epochs, "passes": n_pass,
            "verdict": "PASS" if n_pass >= 3 else "FAIL"}  # fmt: skip
    print(f"{member}: V2 (final wording) {out['verdict']} ({n_pass}/5)", flush=True)
    return out


def cmd_judge(members: list[str]) -> None:
    real = CHAIN.real_ephemeris("ev", "mean")
    for m in members or ROWS:
        _write(f"judge_{m}.json", judge_member(m, real))


def cmd_summary() -> None:
    rows = {}
    for m in ROWS:
        p = OUT / f"judge_{m}.json"
        if p.exists():
            j = json.loads(p.read_text())
            rows[m] = {"L": j["L"], "verdict": j["verdict"], "passes": j["passes"],
                       "c123_epochs": j.get("c123_epochs"),
                       "epochs": [e["outcome"] for e in j.get("epochs", [])]}  # fmt: skip
            print(m, rows[m], flush=True)
    _write("summary.json", rows)


def main() -> None:
    if len(sys.argv) < 2 or sys.argv[1] not in ("controls", "judge", "summary"):
        raise SystemExit("usage: v2_ev10_1025.py {controls,judge,summary} ...")
    if sys.argv[1] == "controls":
        cmd_controls()
    elif sys.argv[1] == "judge":
        cmd_judge(sys.argv[2:])
    else:
        cmd_summary()


if __name__ == "__main__":
    main()
