"""#1065: full re-enumeration of the #998 Pluto cells in the #1024 barycentric binary model.

Pre-registration: docs/notes/2026-10-11-1024-barycentric-binary-model.md sec. 9. The cells are the
#998 ones (scripts/run_998_enumerate.py::pluto_system) with the #1024 option on
(two_working_body_binary.barycentric_binary, Charon mass fraction 106.1 / 975.5), named psb, pnb,
pkb, phb. Charon's half-rev returns are left out of the catalogue: half_rev_vectors refuses a
non-Keplerian body (#1024 amendment A1). Everything else is the #998 setting.

Subcommands:
  enumerate ARGS    the run_942 driver (same flags and outputs), these cells, #1065's preflight
  check DIR [DIR ...]
                    coverage validator: every structure of the cell done once, no duplicated zeros
  gauntlet --cell pnb DIR [DIR ...] --out F
                    every gate-passing merged group (scripts/analyse_942_enumeration.py rules),
                    re-assessed with the #1024 reassess() (gate, re-fly, moving-Pluto surface,
                    Charon SOI, deep-dip diagnostics) plus the #973 unscheduled-pass screen
  controls --cell pnb DIR [DIR ...] --out F
                    sec. 9 controls: the #1024-continued members must be re-found with the same
                    verdict; the #998 default-model zeros must not appear

Usage:
  uv run python scripts/run_1065_barycentric_enumerate.py enumerate --cell pnb --k 1,2,3,4,5,6 \
      --out DIR --max-returns 2,0 --transfer-revs 0,1 --generic-revs 1,2 --n-phase 24 \
      --n-split 10 --n-refine 60 --shard i/n --timing-pilot-s 1.3
"""

from __future__ import annotations

import argparse
import functools
import importlib.util
import json
import sys
from pathlib import Path
from typing import Any

import numpy as np

from cyclerfinder.data.candidate_ids import candidate_cid
from cyclerfinder.data.preflight import preflight_search
from cyclerfinder.search.two_working_body import CircularSystem, sphere_of_influence_km
from cyclerfinder.search.two_working_body_binary import barycentric_binary
from cyclerfinder.search.two_working_body_enum import CatalogueSpec

REPO = Path(__file__).resolve().parents[1]
DAY = 86400.0
GM_CHARON = 106.1  # registry
VINF_TOL_KMS = 1e-4  # sec. 9.4: re-found = same key, Charon V_inf within 0.1 m/s
CELLS = {"psb": "Styx", "pnb": "Nix", "pkb": "Kerberos", "phb": "Hydra"}


def _load(name: str, path: Path) -> Any:
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


ENUM = _load("run_942_enumerate", REPO / "scripts" / "run_942_enumerate.py")
DRV998 = _load("run_998_enumerate", REPO / "scripts" / "run_998_enumerate.py")
R1024 = _load("run_1024_barycentric", REPO / "scripts" / "run_1024_barycentric.py")
R973 = _load("run_973_enumerate", REPO / "scripts" / "run_973_enumerate.py")
ANALYSE = _load("analyse_942_enumeration", REPO / "scripts" / "analyse_942_enumeration.py")
_RUN_942_CELL_SYSTEM = ENUM.cell_system


def bary_system(moon: str) -> CircularSystem:
    off = DRV998.pluto_system(moon)
    return barycentric_binary(off, "Charon", GM_CHARON / off.mu)


def cell_system(cell: str) -> tuple[CircularSystem, str, str]:
    if cell in CELLS:
        return bary_system(CELLS[cell]), "Charon", CELLS[cell]
    result: tuple[CircularSystem, str, str] = _RUN_942_CELL_SYSTEM(cell)
    return result


def _preflight_1065(**kwargs: Any) -> None:
    preflight_search(task_no=1065, script_path=Path(__file__), **kwargs)


def _patch_enum() -> None:
    """run_942's main() and check look up the cell and the catalogue through module globals."""
    ENUM.cell_system = cell_system
    # no half-rev returns at the non-Keplerian Charon (#1024 amendment A1)
    ENUM.CatalogueSpec = functools.partial(CatalogueSpec, half_revs=())


def cmd_enumerate(argv: list[str]) -> None:
    pre = argparse.ArgumentParser(add_help=False)
    pre.add_argument("--cell", required=True)
    cell = pre.parse_known_args(argv)[0].cell
    if cell not in CELLS:
        raise SystemExit(f"#1065 cells are {sorted(CELLS)}")
    _patch_enum()
    sys.argv = [sys.argv[0], *argv]
    ENUM.main(preflight=_preflight_1065, x1=False)


def cmd_check(argv: list[str]) -> None:
    """Coverage validator (run_942 check) over the shard directories of one cell."""
    _patch_enum()
    ENUM.cmd_check(argv)


def _unscheduled(system: Any, cyc: Any, x: np.ndarray) -> dict[str, Any]:
    """The #973 unscheduled-pass screen (amendments 2.6/2.7) for the two cell bodies."""
    cp = R973.closest_unscheduled(system, cyc, x)
    radius = {c: system.body(c).radius_km for c in system.bodies}
    fixed = [p for p in cp if not p["leg"].startswith("R")]
    free = [p for p in cp if p["leg"].startswith("R")]
    impact = [p for p in fixed if p["min_km"] < radius[p["body"]]]
    inside = [p for p in fixed if p["min_km"] < p["soi_km"]]
    flagged = [p for p in free if p["min_km"] < p["soi_km"]]
    status = (
        "reject: body impact"
        if impact
        else "model-invalid"
        if inside
        else "pass, direction-dependent"
        if flagged
        else "pass"
    )
    return {"unscheduled_status": status, "unscheduled_inside_soi": inside + flagged}


def cmd_gauntlet(argv: list[str]) -> None:
    ap = argparse.ArgumentParser(prog="gauntlet")
    ap.add_argument("--cell", required=True, choices=sorted(CELLS))
    ap.add_argument("dirs", nargs="+", type=Path)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args(argv)
    system, a, b = cell_system(args.cell)
    n_struct, recs, groups, passing = ANALYSE.collect(args.dirs)
    out = []
    for idx, ((k, _seq), g) in enumerate(
        sorted(passing.items(), key=lambda kv: (kv[0][0], kv[0][1]))
    ):
        best = min(
            (m for m in g["members"] if ANALYSE.is_pass(m)),
            key=lambda m: m["worst_ratio"] if m["worst_ratio"] is not None else 9,
        )
        _, cyc = ENUM.parse_cycle_key(best["key"], system, a, b)
        x = np.asarray(best["x_days"]) * DAY
        rec = {"index": idx, "k": k, "key": best["key"], "x_days": best["x_days"]}
        # #1071: stable id (new writes only; committed files are addressed by index)
        rec["cid"] = candidate_cid(args.cell, idx, best["x_days"])
        rec |= R1024.reassess(system, cyc, x)
        rec |= _unscheduled(system, cyc, x)
        rec["strong_1065"] = bool(
            rec.get("strong") and rec["unscheduled_status"].startswith("pass")
        )
        out.append(rec)
    summary = {
        "cell": args.cell,
        "n_structures": n_struct,
        "n_zeros": len(recs),
        "n_physical": len(groups),
        "n_gate_passing": len(passing),
        "charon_soi_km": sphere_of_influence_km(system, "Charon"),
        "candidates": out,
    }
    args.out.write_text(json.dumps(summary, indent=1, default=float))
    print(
        f"{args.cell}: structures {n_struct}, zeros {len(recs)}, gate-passing {len(passing)}, "
        f"strong {sum(r['strong_1065'] for r in out)}",
        flush=True,
    )


def cmd_controls(argv: list[str]) -> None:
    """Sec. 9 controls on a finished cell (after the enumeration)."""
    ap = argparse.ArgumentParser(prog="controls")
    ap.add_argument("--cell", required=True, choices=sorted(CELLS))
    ap.add_argument("dirs", nargs="+", type=Path)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args(argv)
    old_cell = args.cell[:2]
    zeros = [z for d in args.dirs for z in map(json.loads, (d / "zeros.jsonl").open())]
    by_key: dict[str, list[dict[str, Any]]] = {}
    for z in zeros:
        by_key.setdefault(z["key"], []).append(z)
    cont = json.loads(
        (REPO / "data" / "1024_barycentric" / f"{old_cell}_continue.json").read_text()
    )

    def near(key: str, v_charon: float) -> list[dict[str, Any]]:
        return [
            z
            for z in by_key.get(key, [])
            if "Charon" in (z.get("vinf_kms") or {})
            and abs(float(z["vinf_kms"]["Charon"]) - v_charon) < VINF_TOL_KMS
        ]

    refound = []
    stale = []
    for r in cont:
        old = near(r["key"], float(r["vinf_charon_998"]))
        if old:
            stale.append(
                {"index": r["index"], "key": r["key"], "vinf_charon_998": r["vinf_charon_998"]}
            )
        if r["outcome"] != "CONTINUES":
            continue
        hits = near(r["key"], float(r["vinf_kms"]["Charon"]))
        same = [z for z in hits if ANALYSE.is_pass(z) == bool(r["pass"])]
        refound.append(
            {"index": r["index"], "key": r["key"], "found": bool(hits), "same_verdict": bool(same)}
        )
    res = {
        "cell": args.cell,
        "continued_members": len(refound),
        "refound": sum(r["found"] for r in refound),
        "refound_same_verdict": sum(r["same_verdict"] for r in refound),
        "not_refound": [r for r in refound if not r["found"]],
        "old_model_zeros_present": stale,
    }
    args.out.write_text(json.dumps(res, indent=1, default=float))
    print(
        f"{args.cell}: re-found {res['refound']} of {res['continued_members']} "
        f"(same verdict {res['refound_same_verdict']}); old-model zeros present: {len(stale)}",
        flush=True,
    )


def main() -> None:
    cmds = {
        "enumerate": cmd_enumerate,
        "check": cmd_check,
        "gauntlet": cmd_gauntlet,
        "controls": cmd_controls,
    }
    if len(sys.argv) < 2 or sys.argv[1] not in cmds:
        raise SystemExit(f"usage: run_1065_barycentric_enumerate.py {{{','.join(cmds)}}} ...")
    cmds[sys.argv[1]](sys.argv[2:])


if __name__ == "__main__":
    main()
