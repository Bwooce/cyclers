"""#827: digit-grade reproduction of Kumar et al. (2026) Table-5 3:1 -> 2:1 states.

Runs :func:`cyclerfinder.search.kumar_em_resonant_heteroclinics.
reproduce_table5_intersection` at each of the four registered target Jacobi
constants (``C in {3.00, 3.05, 3.10, 3.15}`` -- the paper's own Table-1 set) and
records, per C: the converged connection parameters, `#822`'s full verification
battery, and the digit-grade perigee-state comparison against the PRINTED
Table-5 intersection (achieved distance, per-component diffs, runner-up
separation). Results to
``data/found/827_kumar_table5_reproduction/results.json``, checkpointed per C.

Named ``screen_*`` (not ``run_*``): a fixed published-table reproduction
against two already-reproduced orbit families, with no region_id/n_points
sweep-region concept to preflight (same category and precedent as
``scripts/screen_822_vaquero_em_free_transfer.py``), not a catalogue-region
discovery sweep.

Foreground, single-process. Runtime: ~4-5 min per C point.
"""

from __future__ import annotations

import argparse
import json
import time
from dataclasses import asdict
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import cyclerfinder.search.kumar_em_resonant_heteroclinics as keh

OUT_DIR = (
    Path(__file__).resolve().parent.parent / "data" / "found" / "827_kumar_table5_reproduction"
)

_COMPONENTS = ("x", "y", "xdot", "ydot")
_PLANAR_IDX = (0, 1, 3, 4)


def _row(res: keh.Table5Reproduction, elapsed_s: float, settings: dict[str, Any]) -> dict[str, Any]:
    row: dict[str, Any] = {
        "jacobi": res.jacobi,
        "transfer_type": res.transfer_type,
        "printed_state": list(keh.KUMAR_TABLE5_31_TO_21[res.jacobi][1]),
        "matched": res.matched,
        "match_distance": res.match_distance,
        "runner_up_distance": res.runner_up_distance,
        "n_candidates": res.n_candidates,
        "n_refined": res.n_refined,
        "n_converged": res.n_converged,
        "elapsed_s": round(elapsed_s, 1),
        "settings": settings,
        "notes": res.notes,
    }
    if res.connection is not None and res.evidence is not None and res.matched_state is not None:
        conn = asdict(res.connection)
        conn["crossing_xv"] = [float(v) for v in res.connection.crossing_xv]
        ev = asdict(res.evidence)
        ev["crossing_state"] = [float(v) for v in res.evidence.crossing_state]
        row["connection"] = conn
        row["evidence"] = ev
        row["matched_perigee_state"] = [float(v) for v in res.matched_state]
        row["matched_t_from_crossing"] = res.matched_t
        printed = keh.kumar_table5_state6(res.jacobi)
        row["component_diffs"] = {
            name: float(res.matched_state[i] - printed[i])
            for name, i in zip(_COMPONENTS, _PLANAR_IDX, strict=True)
        }
    return row


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--c-values",
        nargs="*",
        type=float,
        default=None,
        help="subset of target C values to (re)run (default: all not-yet-recorded)",
    )
    ap.add_argument("--n-tau", type=int, default=48)
    ap.add_argument("--n-periods", type=float, default=5.0)
    ap.add_argument("--max-refine", type=int, default=8)
    ap.add_argument("--epsilon", type=float, default=keh.KUMAR_EPSILON)
    args = ap.parse_args()

    system = keh.kumar_system()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    path = OUT_DIR / "results.json"
    if path.exists():
        out: dict[str, Any] = json.loads(path.read_text(encoding="utf-8"))
    else:
        out = {
            "task": "#827",
            "generated_utc": datetime.now(UTC).isoformat(),
            "source": (
                "Kumar, Rawat, Rosengren & Ross (2026), Adv. Space Res. 77(3):3815, "
                "DOI 10.1016/j.asr.2025.12.005, Table 5 (3:1 to 2:1 block) / Table 6"
            ),
            "mu": system.mu,
            "match_tol": keh.KUMAR_MATCH_TOL,
            "rows": [],
        }
    done = {row["jacobi"] for row in out["rows"]}
    todo = (
        [round(float(c), 4) for c in args.c_values]
        if args.c_values is not None
        else [c for c in keh.KUMAR_REPRODUCTION_CS if c not in done]
    )
    settings = {
        "n_tau": args.n_tau,
        "n_periods": args.n_periods,
        "max_refine": args.max_refine,
        "epsilon": args.epsilon,
    }

    print(f"#827 reproduction: C points {todo}", flush=True)
    for c in todo:
        print(f"[{datetime.now(UTC).isoformat(timespec='seconds')}] C={c} ...", flush=True)
        t0 = time.time()
        res = keh.reproduce_table5_intersection(
            system,
            c,
            n_tau=args.n_tau,
            n_periods=args.n_periods,
            max_refine=args.max_refine,
            epsilon=args.epsilon,
        )
        row = _row(res, time.time() - t0, settings)
        if res.matched:
            print(
                f"  C={c}: MATCHED dist={res.match_distance:.3e} "
                f"runner_up={res.runner_up_distance:.3e} "
                f"(type {res.transfer_type}; candidates={res.n_candidates} "
                f"refined={res.n_refined} converged={res.n_converged}) "
                f"[{row['elapsed_s']}s]",
                flush=True,
            )
        else:
            print(
                f"  C={c}: NO MATCH best={res.match_distance:.3e} "
                f"(candidates={res.n_candidates} refined={res.n_refined} "
                f"converged={res.n_converged}) [{row['elapsed_s']}s]",
                flush=True,
            )
        out["rows"] = [r for r in out["rows"] if r["jacobi"] != c] + [row]
        out["rows"].sort(key=lambda r: r["jacobi"])
        path.write_text(json.dumps(out, indent=2), encoding="utf-8")  # checkpoint per C

    n_ok = sum(1 for r in out["rows"] if r["matched"])
    print(f"\nmatched {n_ok}/{len(out['rows'])} recorded C points; results at {path}", flush=True)


if __name__ == "__main__":
    main()
