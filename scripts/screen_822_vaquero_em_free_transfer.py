"""#822: Vaquero 2:1 <-> 3:1 Earth-Moon free-transfer sweep across the overlap band.

Runs :func:`cyclerfinder.search.vaquero_em_cycler_connections.find_free_transfer`
(Wu(2:1) -> Ws(3:1), same-C heteroclinic manifold connection) at every shared
0.01-grid Jacobi constant in the two families' printed overlap band
``C in [2.54, 2.66]`` -- the band where Vaquero (2013, pp. 171-172) asserts a
free transfer exists but prints no trajectory. Additionally demonstrates the
REVERSE direction (Wu(3:1) -> Ws(2:1)) at the mid-band C=2.60 via the module's
own primitives (the CR3BP time-reversal symmetry guarantees it given the
forward hit; demonstrated directly rather than merely asserted).

Writes the full per-C record (connection parameters + complete verification
evidence) to ``data/found/822_vaquero_em_free_transfer/results.json``,
checkpointing after every C grid point (incremental, restart-safe to read).

Named ``screen_*`` (not ``run_*``): a fixed published-claim reproduction
against two already-reproduced families, with no region_id/n_points
sweep-region concept to preflight (same category and precedent as
``scripts/screen_799_vaquero_em_cycler_families.py``), not a catalogue-region
discovery sweep.

Foreground, single-process. Runtime: ~2-4 min per C point (13 points + the
reverse demo), ~40 min total.
"""

from __future__ import annotations

import argparse
import json
import math
import time
from dataclasses import asdict
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import cyclerfinder.core.cr3bp as cr3bp
import cyclerfinder.search.vaquero_em_cycler_connections as vcc

OUT_DIR = Path(__file__).resolve().parent.parent / "data" / "found" / "822_vaquero_em_free_transfer"

#: The 13 shared 0.01-grid points of the overlap band (both endpoints inclusive).
C_GRID = [round(2.54 + 0.01 * i, 4) for i in range(13)]

#: Mid-band C for the reverse-direction (Wu(3:1) -> Ws(2:1)) demonstration.
REVERSE_C = 2.60


def _conn_dict(conn: Any) -> dict[str, Any]:
    d = asdict(conn)
    d["crossing_xv"] = [float(v) for v in conn.crossing_xv]
    return d


def _ev_dict(ev: Any) -> dict[str, Any]:
    d = asdict(ev)
    d["crossing_state"] = [float(v) for v in ev.crossing_state]
    return d


def _sweep_one(
    system: cr3bp.CR3BPSystem,
    c: float,
    *,
    n_tau: int = 36,
    n_periods: float = 9.0,
    max_refine: int = 8,
    epsilon: float = vcc.EPSILON,
) -> dict[str, Any]:
    t0 = time.time()
    res = vcc.find_free_transfer(
        system, c, n_tau=n_tau, n_periods=n_periods, max_refine=max_refine, epsilon=epsilon
    )
    dt = time.time() - t0
    row: dict[str, Any] = {
        "jacobi": c,
        "direction": "Wu(2:1)->Ws(3:1)",
        "n_tau": n_tau,
        "n_periods": n_periods,
        "max_refine": max_refine,
        "epsilon": epsilon,
        "elapsed_s": round(dt, 1),
        "n_seeds": res.n_seeds,
        "n_refined": res.n_refined,
        "n_converged": res.n_converged,
        "verified": res.connection is not None,
        "notes": res.notes,
    }
    if res.connection is not None and res.evidence is not None:
        row["connection"] = _conn_dict(res.connection)
        row["evidence"] = _ev_dict(res.evidence)
        ev = res.evidence
        days = system.t_s / 86400.0
        cxv = res.connection.crossing_xv
        print(
            f"  C={c}: VERIFIED residual={res.connection.residual:.2e} "
            f"crossing=({cxv[0]:+.6f}, {cxv[1]:+.6f}) "
            f"k=({res.connection.k_u},{res.connection.k_s}) "
            f"branches=({res.connection.branch_u:+d},{res.connection.branch_s:+d}) "
            f"t_u={ev.t_u * days:.1f}d t_s={abs(ev.t_s) * days:.1f}d "
            f"full_gap={ev.full_state_gap:.1e} radau={ev.radau_gap:.1e} "
            f"ghost=({ev.ghost_distance_from:.3f},{ev.ghost_distance_to:.3f}) [{dt:.0f}s]",
            flush=True,
        )
    else:
        print(
            f"  C={c}: no verified connection (seeds={res.n_seeds} refined={res.n_refined} "
            f"converged={res.n_converged}) [{dt:.0f}s]",
            flush=True,
        )
    return row


def _reverse_demo(
    system: cr3bp.CR3BPSystem,
    c: float,
    *,
    n_tau: int = 36,
    n_periods: float = 9.0,
    max_refine: int = 12,
    epsilon: float = vcc.EPSILON,
) -> dict[str, Any]:
    """Wu(3:1) -> Ws(2:1) at ``c`` via the module primitives (nodes swapped),
    with the same seed-diversity discipline as ``find_free_transfer``."""
    t0 = time.time()
    node21 = vcc.build_vaquero_overlap_node(system, 2, c)
    node31 = vcc.build_vaquero_overlap_node(system, 3, c)
    cross_u = {
        b: vcc.manifold_section_crossings(
            system,
            node31,
            direction="unstable",
            branch=b,
            n_tau=n_tau,
            n_periods=n_periods,
            epsilon=epsilon,
        )
        for b in (+1, -1)
    }
    cross_s = {
        b: vcc.manifold_section_crossings(
            system,
            node21,
            direction="stable",
            branch=b,
            n_tau=n_tau,
            n_periods=n_periods,
            epsilon=epsilon,
        )
        for b in (+1, -1)
    }
    seeds: list[vcc.ConnectionSeed] = []
    for bu in (+1, -1):
        for bs in (+1, -1):
            seeds.extend(
                vcc.find_connection_seeds(cross_u[bu], cross_s[bs], branch_u=bu, branch_s=bs)
            )
    seeds.sort(key=lambda s: s.distance)
    row: dict[str, Any] = {
        "jacobi": c,
        "direction": "Wu(3:1)->Ws(2:1)",
        "n_tau": n_tau,
        "n_periods": n_periods,
        "max_refine": max_refine,
        "epsilon": epsilon,
        "n_seeds": len(seeds),
        "verified": False,
    }
    n_ref = n_conv = 0
    tried: list[tuple[float, float]] = []
    for seed in seeds:
        if n_ref >= max_refine:
            break
        if any(math.hypot(seed.x - x, seed.xdot - xd) < 0.05 for x, xd in tried):
            continue
        tried.append((seed.x, seed.xdot))
        n_ref += 1
        conn = vcc.refine_connection(system, node31, node21, seed, epsilon=epsilon)
        if not conn.converged:
            continue
        n_conv += 1
        ev = vcc.verify_connection(system, node31, node21, conn, epsilon=epsilon)
        if ev.passed:
            row["verified"] = True
            row["connection"] = _conn_dict(conn)
            row["evidence"] = _ev_dict(ev)
            print(
                f"  REVERSE C={c}: VERIFIED residual={conn.residual:.2e} "
                f"crossing=({conn.crossing_xv[0]:+.6f}, {conn.crossing_xv[1]:+.6f}) "
                f"k=({conn.k_u},{conn.k_s}) [{time.time() - t0:.0f}s]",
                flush=True,
            )
            break
    row["n_refined"] = n_ref
    row["n_converged"] = n_conv
    row["elapsed_s"] = round(time.time() - t0, 1)
    if not row["verified"]:
        print(f"  REVERSE C={c}: no verified connection [{time.time() - t0:.0f}s]", flush=True)
    return row


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--c-values",
        nargs="*",
        type=float,
        default=None,
        help="subset of C grid points to (re)run this invocation (default: all "
        "not-yet-recorded points); the results file is merged/checkpointed per C, "
        "so the sweep can be driven in bounded foreground chunks",
    )
    ap.add_argument(
        "--reverse-demo",
        action="store_true",
        help=f"run the reverse-direction demonstration at C={REVERSE_C}",
    )
    ap.add_argument("--n-tau", type=int, default=36, help="manifold phase-grid density")
    ap.add_argument(
        "--max-refine",
        type=int,
        default=8,
        help="Newton-refinement budget per C (diversity-filtered, closest seeds first)",
    )
    ap.add_argument(
        "--epsilon",
        type=float,
        default=vcc.EPSILON,
        help="manifold-offset magnitude; the high-C end of the band needs a LARGER "
        "offset (e.g. 1e-4) so the legs span fewer Floquet amplification periods and "
        "the re-derivation evidence floor stays under the verify_connection gates "
        "(the gates themselves are never loosened) -- recorded per row",
    )
    ap.add_argument(
        "--n-periods",
        type=float,
        default=9.0,
        help="manifold horizon in orbit periods; the low-C end of the band needs "
        "more (the 2:1 saddle weakens to |lambda|=3.15 at C=2.54, so epsilon-scale "
        "offsets need ~13 periods to reach O(1) separation)",
    )
    args = ap.parse_args()

    system = cr3bp.cr3bp_system("Earth", "Moon")
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    path = OUT_DIR / "results.json"
    if path.exists():
        out: dict[str, Any] = json.loads(path.read_text(encoding="utf-8"))
    else:
        out = {
            "task": "#822",
            "generated_utc": datetime.now(UTC).isoformat(),
            "mu": system.mu,
            "l_km": system.l_km,
            "t_s": system.t_s,
            "epsilon": vcc.EPSILON,
            "settings": {"n_tau": 36, "n_periods": 9.0, "max_refine": 8, "tol": 1e-9},
            "sweep": [],
            "reverse_demo": None,
        }
    done = {row["jacobi"] for row in out["sweep"]}
    todo = args.c_values if args.c_values is not None else [c for c in C_GRID if c not in done]

    print(f"#822 sweep: C points {todo} of band {vcc.OVERLAP_C_RANGE}", flush=True)
    for c in todo:
        c = round(float(c), 4)
        print(f"[{datetime.now(UTC).isoformat(timespec='seconds')}] C={c} ...", flush=True)
        row = _sweep_one(
            system,
            c,
            n_tau=args.n_tau,
            n_periods=args.n_periods,
            max_refine=args.max_refine,
            epsilon=args.epsilon,
        )
        out["sweep"] = [r for r in out["sweep"] if r["jacobi"] != c] + [row]
        out["sweep"].sort(key=lambda r: r["jacobi"])
        path.write_text(json.dumps(out, indent=2), encoding="utf-8")  # checkpoint per C

    if args.reverse_demo:
        print(
            f"[{datetime.now(UTC).isoformat(timespec='seconds')}] reverse demo C={REVERSE_C} ...",
            flush=True,
        )
        out["reverse_demo"] = _reverse_demo(
            system,
            REVERSE_C,
            n_tau=args.n_tau,
            n_periods=args.n_periods,
            max_refine=args.max_refine,
            epsilon=args.epsilon,
        )
        path.write_text(json.dumps(out, indent=2), encoding="utf-8")

    n_ok = sum(1 for r in out["sweep"] if r["verified"])
    print(f"\nverified {n_ok}/{len(out['sweep'])} recorded C points; results at {path}", flush=True)


if __name__ == "__main__":
    main()
