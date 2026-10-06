"""#942/#943 ideal-model enumeration driver (two-working-body generator).

Cells (circular-coplanar ideal model):
  ev  Earth-Venus, both massive              (#942 cell b; Venus 0.61520 yr)
  em  Earth-Mars, both massive               (#942 cell a; Mars 1.875 yr, Russell)
  vm  Venus-Mars, Mars massless              (#942 cell c; R-S 2007 Table 2 constants)
  vm2 Venus-Mars, both massive               (#942 cell c, two working bodies; same constants)
  gc  Ganymede-Callisto, both massive        (#943 X1; R-S 2009 Table 2 constants)
  ge  Ganymede-Europa, both massive          (#943 X1)
  gc1 / ge1  the one-body limits (Callisto / Europa massless), recall controls

Each structure (cycle template) is solved from a seed grid and every distinct
zero is assessed (gate at the project floor and at H&M's 1.1 radii, near-180
rule, encounter self-consistency). Output: one JSON line per structure in
``<out>/structures.jsonl`` (resumable: finished structure keys are skipped)
and one per zero in ``<out>/zeros.jsonl``. Progress with ETA every structure.

Usage:
  uv run python scripts/run_942_enumerate.py --cell ev --k 2 --out DIR [--count-only]
      [--max-returns 2,2] [--transfer-revs 0] [--generic-revs 1] [--n-phase 24]
      [--n-split 10] [--n-refine 60] [--sample N --seed S] [--shard i/n]
"""

from __future__ import annotations

import argparse
import json
import math
import random
import time
from collections.abc import Callable
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

from cyclerfinder.data.method_capability import MethodCapability
from cyclerfinder.data.preflight import preflight_search
from cyclerfinder.search.two_working_body import (
    CircularSystem,
    Cycle,
    FlybyBody,
    HalfRevLeg,
    LambertLeg,
    ResonantLeg,
    heliocentric_circular,
)
from cyclerfinder.search.two_working_body_enum import (
    CatalogueSpec,
    assess,
    flyby_table,
    solve_structure,
    structures,
)
from cyclerfinder.verify.turn_gate_closures import RS_BODIES, RS_PRIMARY_GM

DAY = 86400.0


def log(msg: str) -> None:
    print(f"{datetime.now(UTC).isoformat(timespec='seconds')} {msg}", flush=True)


def rs_moon_system(moons: list[str], massless: list[str]) -> CircularSystem:
    """Jupiter moons with Russell & Strange 2009 Table 2 constants (p.148)."""
    mu = RS_PRIMARY_GM["Jupiter"]
    bodies = {}
    over = {}
    for m in moons:
        b = RS_BODIES[m]
        a = (mu * (b.ideal_period_s / (2.0 * math.pi)) ** 2) ** (1.0 / 3.0)
        bodies[m] = (a, b.ideal_period_s, 0.0)
        # the project floor (registry safe_alt_km) with the paper's GM and radius
        from cyclerfinder.verify.turn_gate import body_constants

        over[m] = FlybyBody(m, b.gm_km3_s2, b.radius_km, body_constants(m).alt_floor_km)
    return CircularSystem(mu, bodies, frozenset(massless), flyby_overrides=over)


def rs2007_venus_mars(*, mars_massless: bool = True) -> CircularSystem:
    """Venus-Mars ideal model with Russell & Strange 2007 (AAS 07-118) Table 2
    constants (p.8): Sun mu 1.3271244e11, Venus period 19,414,153 s, Mars
    59,354,429 s, Venus mu 324,860 and radius 6,052 km; Mars massless. The
    same model as their VenMar#45, so that row is the cell's recall control."""
    mu = 1.3271244e11
    bodies = {}
    for c, per in (("V", 19_414_153.0), ("M", 59_354_429.0)):
        bodies[c] = ((mu * (per / (2.0 * math.pi)) ** 2) ** (1.0 / 3.0), per, 0.0)
    from cyclerfinder.verify.turn_gate import body_constants

    over = {
        "V": FlybyBody("V", 324_860.0, 6052.0, body_constants("V").alt_floor_km),
        "M": FlybyBody("M", 42_828.3, 3399.0, body_constants("M").alt_floor_km),
    }
    massless = frozenset({"M"}) if mars_massless else frozenset()
    return CircularSystem(mu, bodies, massless, flyby_overrides=over)


def cell_system(cell: str) -> tuple[CircularSystem, str, str]:
    if cell == "ev":
        return heliocentric_circular({"E": 1.0, "V": 0.61520}), "E", "V"
    if cell == "em":
        return heliocentric_circular({"E": 1.0, "M": 1.875}), "E", "M"
    if cell == "vm":
        return rs2007_venus_mars(), "V", "M"
    if cell == "vm2":
        return rs2007_venus_mars(mars_massless=False), "V", "M"
    if cell == "gc":
        return rs_moon_system(["Ganymede", "Callisto"], []), "Ganymede", "Callisto"
    if cell == "ge":
        return rs_moon_system(["Ganymede", "Europa"], []), "Ganymede", "Europa"
    if cell == "gc1":
        return rs_moon_system(["Ganymede", "Callisto"], ["Callisto"]), "Ganymede", "Callisto"
    if cell == "ge1":
        return rs_moon_system(["Ganymede", "Europa"], ["Europa"]), "Ganymede", "Europa"
    raise ValueError(cell)


def leg_key(leg: object) -> str:
    if isinstance(leg, LambertLeg):
        return f"L{leg.frm}>{leg.to}/{leg.nrev}{leg.branch[0]}"
    if isinstance(leg, ResonantLeg):
        return f"R{leg.body}/{leg.body_revs}:{leg.sc_revs}"
    if isinstance(leg, HalfRevLeg):
        return f"H{leg.body}/{leg.half_periods},{leg.k_sc},{'p' if leg.via_peri else 'a'}"
    raise TypeError(leg)


def parse_leg(token: str) -> LambertLeg | ResonantLeg | HalfRevLeg:
    """Inverse of :func:`leg_key`."""
    kind, rest = token[0], token[1:]
    if kind == "L":
        pair, rev = rest.split("/")
        frm, to = pair.split(">")
        n = int(rev[:-1])
        br = {"s": "single", "l": "low", "h": "high"}[rev[-1]]
        return LambertLeg(frm, to, n, br)
    if kind == "R":
        body, nm = rest.split("/")
        n, m = nm.split(":")
        return ResonantLeg(body, int(n), int(m))
    if kind == "H":
        body, args = rest.split("/")
        h, k, pa = args.split(",")
        return HalfRevLeg(body, int(h), int(k), pa == "p", 0)
    raise ValueError(token)


def parse_cycle_key(key: str, system: CircularSystem, a: str, b: str) -> tuple[int, Cycle]:
    """Inverse of :func:`cycle_key`: ``(k, Cycle)``."""
    parts = key.split("|")
    k = int(parts[0][1:])
    legs = tuple(parse_leg(t) for t in parts[1:])
    return k, Cycle(legs, k * system.synodic_s(a, b))


def cycle_key(c: Cycle, k: int) -> str:
    return f"k{k}|" + "|".join(leg_key(lg) for lg in c.legs)


#: Cells that belong to #943 (X1, Jovian); they run through scripts/run_943_enumerate.py.
X1_CELLS = frozenset({"gc", "ge", "gc1", "ge1"})


def _preflight_942(**kwargs: Any) -> None:
    preflight_search(task_no=942, script_path=Path(__file__), **kwargs)


def main(preflight: Callable[..., None] = _preflight_942, x1: bool = False) -> None:
    """Run one cell. ``preflight`` is the task-specific gate (#942 here, #943 in
    run_943_enumerate.py); ``x1`` selects which cells this entry point accepts."""
    ap = argparse.ArgumentParser()
    ap.add_argument("--cell", required=True)
    ap.add_argument("--k", type=str, required=True, help="comma list of synodic multiples")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--max-returns", type=str, default="2,2")
    ap.add_argument("--transfer-revs", type=str, default="0")
    ap.add_argument("--generic-revs", type=str, default="1")
    ap.add_argument("--visits", type=int, default=1)
    ap.add_argument(
        "--resonant-only",
        type=str,
        default="",
        help="comma list of bodies whose returns are full-revolution (n:m resonant) only",
    )
    ap.add_argument("--n-phase", type=int, default=24)
    ap.add_argument("--n-split", type=int, default=10)
    ap.add_argument("--n-refine", type=int, default=60)
    ap.add_argument("--count-only", action="store_true")
    ap.add_argument("--sample", type=int, default=0)
    ap.add_argument("--seed", type=int, default=942)
    ap.add_argument("--shard", type=str, default="0/1")
    ap.add_argument(
        "--timing-pilot-s",
        type=float,
        default=None,
        help="measured seconds per structure from a --sample pilot (preflight gate)",
    )
    args = ap.parse_args()
    if (args.cell in X1_CELLS) != x1:
        owner = (
            "#943 (run_943_enumerate.py)"
            if args.cell in X1_CELLS
            else "#942 (run_942_enumerate.py)"
        )
        raise SystemExit(f"cell {args.cell!r} belongs to {owner}")

    system, a, b = cell_system(args.cell)
    ma, mb = (int(v) for v in args.max_returns.split(","))
    if system.body(b).massless:
        mb = 0
    generic = tuple(int(v) for v in args.generic_revs.split(","))
    res_only = {c for c in args.resonant_only.split(",") if c}
    spec = {
        c: (
            CatalogueSpec(half_revs=(), generic_revs=())
            if c in res_only
            else CatalogueSpec(generic_revs=generic)
        )
        for c in (a, b)
    }
    t_revs = tuple(int(v) for v in args.transfer_revs.split(","))
    ks = [int(v) for v in args.k.split(",")]
    shard_i, shard_n = (int(v) for v in args.shard.split("/"))
    syn = system.synodic_s(a, b)

    todo: list[tuple[int, Cycle]] = []
    for k in ks:
        n_k = 0
        for c in structures(
            system,
            a,
            b,
            k,
            max_returns={a: ma, b: mb},
            spec=spec,
            transfer_revs=t_revs,
            visits=args.visits,
        ):
            n_k += 1
            todo.append((k, c))
        log(f"cell {args.cell} k={k}: {n_k} structures (T = {k * syn / DAY:.2f} d)")
    if args.sample:
        random.Random(args.seed).shuffle(todo)
        todo = todo[: args.sample]
    todo = [t for i, t in enumerate(todo) if i % shard_n == shard_i]
    log(f"total to run: {len(todo)} (shard {shard_i}/{shard_n})")
    if args.count_only:
        return
    preflight(
        region_id=(
            f"{args.cell}-two-working-body-ideal-k{args.k.replace(',', '-')}"
            f"-ret{args.max_returns.replace(',', '-')}-tr{args.transfer_revs.replace(',', '-')}"
            f"-gen{args.generic_revs.replace(',', '-')}-v{args.visits}"
        ),
        method=MethodCapability(
            genome=(
                "two-working-body cycle templates [A-block, A->B, B-block, B->A] x visits: "
                "full-rev n:m, half-rev n-pi and generic returns; Lambert transfers"
            ),
            corrector="two_working_body.correct_dates (date residual) + minimax turn gate",
            capability_tags=frozenset({"ballistic", "coplanar", "patched-conic", "circular"}),
            git_sha="working-tree",
        ),
        n_points=len(todo),
        timing_pilot_seconds_per_point=args.timing_pilot_s,
    )

    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "settings.json").write_text(
        json.dumps(vars(args) | {"synodic_days": syn / DAY}, default=str, indent=1)
    )
    done_path = args.out / "structures.jsonl"
    done = set()
    if done_path.exists():
        for line in done_path.read_text().splitlines():
            d = json.loads(line)
            if not d.get("error"):  # errored structures are retried on resume
                done.add(d["key"])
    t0 = time.time()
    n_new = 0
    n_zero = n_pass = n_err = n_zero_err = 0
    for i, (k, c) in enumerate(todo):
        key = cycle_key(c, k)
        if key in done:
            continue
        ts = time.time()
        try:
            zeros = solve_structure(
                system,
                c,
                phase_period_s=syn,
                n_phase=args.n_phase,
                n_split=args.n_split,
                n_refine=args.n_refine,
            )
        except Exception as exc:  # one structure must never kill the cell
            n_err += 1
            with done_path.open("a") as fd:
                fd.write(json.dumps({"key": key, "k": k, "error": repr(exc)}) + "\n")
            log(f"[{i + 1}/{len(todo)}] {key}: ERROR in solve: {exc!r} (errors={n_err})")
            continue
        statuses = []
        with (args.out / "zeros.jsonl").open("a") as fz:
            for z in zeros:
                try:
                    ass = assess(system, z)
                except Exception as exc:  # recorded, never counted as a fail
                    n_zero_err += 1
                    fz.write(
                        json.dumps(
                            {
                                "key": key,
                                "k": k,
                                "x_days": (z.x / DAY).tolist(),
                                "residual_kms": z.residual_kms,
                                "status": "error",
                                "error": repr(exc),
                            }
                        )
                        + "\n"
                    )
                    statuses.append("error")
                    log(f"    zero assessment ERROR in {key}: {exc!r}")
                    continue
                rep = ass.report
                rec = {
                    "key": key,
                    "k": k,
                    "x_days": (z.x / DAY).tolist(),
                    "residual_kms": z.residual_kms,
                    "status": ass.status,
                    "status_hm_floor": ass.report_hm_floor.status if ass.report_hm_floor else None,
                    "worst_ratio": rep.gate.worst_ratio if rep else None,
                    "max_turn_deg": rep.max_turn_deg if rep else None,
                    "min_required_alt_km": rep.gate.min_required_alt_km if rep else None,
                    "max_encounter_miss_km": ass.max_encounter_miss_km,
                    "min_soi_km": ass.min_soi_km,
                    "r_min_km": ass.r_min_km,
                    "r_max_km": ass.r_max_km,
                    "vinf_kms": ass.vinf_kms,
                    "flybys": flyby_table(system, ass),
                }
                fz.write(json.dumps(rec, default=float) + "\n")
                statuses.append(ass.status)
        n_zero += len(zeros)
        n_pass += statuses.count("pass")
        with done_path.open("a") as fd:
            fd.write(
                json.dumps(
                    {
                        "key": key,
                        "k": k,
                        "n_zeros": len(zeros),
                        "statuses": statuses,
                        "secs": time.time() - ts,
                    }
                )
                + "\n"
            )
        n_new += 1
        el = time.time() - t0
        eta = el / n_new * (len(todo) - i - 1)
        log(
            f"[{i + 1}/{len(todo)}] {key}: {len(zeros)} zeros {statuses.count('pass')} pass "
            f"({time.time() - ts:.1f}s) totals zeros={n_zero} pass={n_pass} eta {eta / 60:.1f} min"
        )
    log(
        f"DONE structures={len(todo)} zeros={n_zero} pass={n_pass} "
        f"structure_errors={n_err} zero_assessment_errors={n_zero_err}"
    )


if __name__ == "__main__":
    main()
