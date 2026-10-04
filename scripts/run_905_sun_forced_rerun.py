"""#905 -- rerun of the #884 Sun-forced Earth-Moon periodic-orbit search in the corrected models.

Driver: :mod:`cyclerfinder.search.sun_forced_905` (multiple shooting over N Sun periods, the
Rhouma-Chicone screens, the Melnikov zeros in the Sun phase, continuation in eps from 0 to 1 with
fold / branch-point / eps = 0 crossing records, bifurcation detection in the family walk, refined
passes with surface exclusion, one orbit per symmetry class). Published positive controls are in
``tests/search/test_sun_forced_905.py`` and must pass before this is run.

Stages (each writes one JSON per task under ``data/found/905_sun_forced_rerun/<stage>/``;
re-running skips finished tasks, so an interrupted run resumes):

  --stage controls   Leiva & Briozzo (2008) 5/2 members 180A_1, 180A_2, 357 (C32, C32, C31) and
                     the 4T orbit 013 in the coherent model at the paper's mass ratio, every
                     Melnikov phase, full diagnostics (the members they obtained only as arcs).
  --stage families   walk every catalogued Earth-Moon family in the three-body problem with
                     bifurcation detection; members at the commensurate periods.
  --stage forced     per member: screens, Melnikov zeros, and per model (bcr4bp, qbcp) and per
                     simple zero the continuation in eps and the diagnostics.
  --stage summary    summary.json: per member and model, branches reached, symmetry classes,
                     largest multipliers, periselene ranges, flags.

Runlog: ``runlog.jsonl`` in the output directory, one flushed line per finished task with the
UTC time, stage, task, status, seconds, running counts and an ETA from the mean task time.

Launch (the coordinator owns this):
  timeout 36h uv run python scripts/run_905_sun_forced_rerun.py --stage all --workers 4 \\
      > data/found/905_sun_forced_rerun/run.log 2>&1
"""

from __future__ import annotations

import argparse
import dataclasses
import datetime
import json
import math
import sys
import time
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any

import numpy as np
import yaml

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from cyclerfinder.core import qbcp  # noqa: E402
from cyclerfinder.data.method_capability import MethodCapability  # noqa: E402
from cyclerfinder.data.preflight import preflight_search  # noqa: E402
from cyclerfinder.search import sun_forced_905 as sf  # noqa: E402

OUT = ROOT / "data" / "found" / "905_sun_forced_rerun"
EM_MU_ROWS = (0.0121505, 0.0121506)

_REGION_ID = "earth-moon-sun-forced-commensurate-po-bcr4bp-qbcp-corrected-sense-905-2026-10-05"
_METHOD = MethodCapability(
    genome=(
        "catalogued Earth-Moon CR3BP periodic families at Sun-commensurate periods "
        "(k Tg / 2, 8 Tg / 3, 5 Tg / 6), Melnikov-zero Sun phases, continued in the Sun's "
        "strength to the corrected bicircular and coherent models"
    ),
    corrector="sun_forced_905.continue_in_eps (multiple shooting over N Sun periods)",
    capability_tags=frozenset(
        {
            "bcr4bp",
            "qbcp",
            "multiple-shooting",
            "family-continuation",
            "melnikov",
            "planar",
            "spatial",
        }
    ),
    git_sha="working-tree",
)

# One seed per catalogued family (same list as #884).
FAMILY_SEEDS: dict[str, str] = {
    "C11 (Braik-Ross C11a)": "braik-ross-c11a-cycler-2026",
    "C11 (Braik-Ross C11b)": "braik-ross-c11b-cycler-2026",
    "C11 (Ross-RT 1:1)": "ross-rt-em-cycler-11-2025",
    "C21 planar (Ross-RT 2:1)": "ross-rt-em-cycler-21-2025",
    "C21 spatial (#438)": "em-cycler-21-3d-spatial-2026",
    "C21 3D corridor (#682)": "braik-ross-c21-3d-corridor-05-2026",
    "C31 (Ross-RT 3:1)": "ross-rt-em-cycler-31-2025",
    "C32 (Ross-RT 3:2)": "ross-rt-em-cycler-32-2025",
    "C32 (Braik-Ross)": "braik-ross-c32-cycler-2026",
    "C33 (Ross-RT 3:3)": "ross-rt-em-cycler-33-2025",
    "C32 3D corridor (#682)": "braik-ross-c32-3d-corridor-05-2026",
    "L1 Lyapunov 3D corridor (#682)": "lyapunov3d-l1-corridor-04-2026",
    "R21-S planar resonant (Braik-Ross)": "braik-ross-planar-r21-s-corridor-2026",
    "R31-S planar resonant (Braik-Ross)": "braik-ross-planar-r31-s-corridor-2026",
    "R52-S planar resonant (Braik-Ross)": "braik-ross-planar-r52-s-corridor-2026",
    "Casoliva 1:2 (c)": "casoliva-1-2c-em-resonant-po-2010",
    "Casoliva 1:2 (d)": "casoliva-1-2d-em-resonant-po-2010",
    "Casoliva 2:1 (a)": "casoliva-2-1a-em-resonant-po-2010",
    "Casoliva 2:1 (b)": "casoliva-2-1b-em-resonant-po-2010",
    "Casoliva 3:2 (c)": "casoliva-3-2c-em-resonant-po-2010",
    "Casoliva 7:3 (a)": "casoliva-7-3a-em-cycler-2010",
    "Casoliva 7:3 (b)": "casoliva-7-3b-em-cycler-2010",
    "Casoliva 7:3 (c)": "casoliva-7-3c-em-cycler-2010",
    "Vaquero 2:1": "vaquero-21-c246-em-cycler-2013",
    "Vaquero 3:1": "vaquero-31-c254-em-cycler-2013",
}

# Leiva & Briozzo (2008) Table 1 (p234): h, y, ydot, p, q in their frame (Earth at +mu).
LB_MU = 0.0121505482
LB_SECTION_X = 0.836915310
LB_TABLE_1 = {
    "180A_1": (-1.59005198, 0.00520342002, 0.0479318298, 5, 2),
    "180A_2": (-1.57583831, -0.0283283340, 0.117872065, 5, 2),
    "357": (-1.52791268, 0.129037155, 0.115396082, 5, 2),
    "013": (-1.58740571, -0.0399746624, -0.0441622383, 4, 1),
}
# Their Table 3 arcs at the same members (t_i, printed d_M km, |s1|) for comparison only.
LB_ARCS = {
    "180A_1": [(1.51986327, 7371, 49.9), (4.91546020, 7365, 49.5)],
    "180A_2": [(1.17849132, 23966, 5153.0)],
    "357": [(2.82792318, 14713, 5365.7), (6.22352012, 14716, 5377.8)],
}


def utc() -> str:
    return datetime.datetime.now(datetime.UTC).strftime("%Y-%m-%dT%H:%M:%SZ")


def log(msg: str) -> None:
    print(f"[{utc()}] {msg}", flush=True)


def _write(path: Path, obj: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".tmp")
    tmp.write_text(json.dumps(obj, indent=1, default=_json_default))
    tmp.replace(path)


def _json_default(o: Any) -> Any:
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, np.generic):
        return o.item()
    if isinstance(o, sf.SymmetricMember):
        return dataclasses.asdict(o)
    raise TypeError(type(o))


class Runlog:
    """Append-and-flush runlog with running counts and an ETA from the mean task time."""

    def __init__(self, path: Path, stage: str, total: int) -> None:
        self.path = path
        self.stage = stage
        self.total = total
        self.done = 0
        self.counts: dict[str, int] = {}
        self.t0 = time.monotonic()
        path.parent.mkdir(parents=True, exist_ok=True)
        self._emit({"event": "stage_start", "tasks": total})

    def _emit(self, rec: dict[str, Any]) -> None:
        rec = {"utc": utc(), "stage": self.stage, **rec}
        with self.path.open("a", encoding="utf-8") as fh:
            fh.write(json.dumps(rec, default=_json_default) + "\n")
            fh.flush()

    def task(self, name: str, status: str, seconds: float) -> None:
        self.done += 1
        self.counts[status] = self.counts.get(status, 0) + 1
        elapsed = time.monotonic() - self.t0
        eta_s = elapsed / self.done * (self.total - self.done)
        eta = (datetime.datetime.now(datetime.UTC) + datetime.timedelta(seconds=eta_s)).strftime(
            "%Y-%m-%dT%H:%M:%SZ"
        )
        self._emit(
            {
                "task": name,
                "status": status,
                "seconds": round(seconds, 1),
                "done": self.done,
                "total": self.total,
                "counts": self.counts,
                "eta_utc": eta,
            }
        )
        log(
            f"[{self.stage}] {self.done}/{self.total} {name}: {status} ({seconds:.0f} s) "
            f"counts={self.counts} ETA {eta}"
        )


# ---------------------------------------------------------------------------
# Seeds and targets.
# ---------------------------------------------------------------------------


def _perpendicular_crossing(mu: float, state: list[float], period: float) -> list[float] | None:
    from scipy.integrate import solve_ivp

    from cyclerfinder.core import cr3bp

    x = np.asarray(state, dtype=float)
    if abs(x[1]) < 1e-9 and abs(x[3]) < 1e-6 and abs(x[5]) < 1e-6:
        return [float(x[0]), float(x[2]), float(x[4])]

    def ev(t: float, y: Any) -> float:
        return float(y[1])

    sol = solve_ivp(
        lambda t, y: cr3bp.cr3bp_eom(t, y, mu),
        (0, period),
        x,
        events=ev,
        rtol=1e-12,
        atol=1e-12,
        method="DOP853",
    )
    for y in sol.y_events[0]:
        if abs(y[3]) < 1e-6 and abs(y[5]) < 1e-6:
            return [float(y[0]), float(y[2]), float(y[4])]
    return None


def catalogue_seed(row_id: str) -> dict[str, Any] | None:
    rows = yaml.safe_load((ROOT / "data" / "catalogue.yaml").read_text())
    for r in rows:
        if r.get("id") != row_id:
            continue
        c = ((r.get("orbit_elements") or {}).get("cr3bp")) or {}
        if c.get("period_nd") is None or not c.get("state_nd"):
            return None
        mu = float(c["mass_ratio"])
        if not (EM_MU_ROWS[0] < mu < EM_MU_ROWS[1]):
            return None
        return {
            "row": row_id,
            "period": float(c["period_nd"]),
            "seed": _perpendicular_crossing(mu, c["state_nd"], float(c["period_nd"])),
            "orbit_class": r.get("orbit_class"),
        }
    return None


def targets(tg: float) -> list[float]:
    return [k * tg / 2 for k in range(1, 13)] + [8 * tg / 3, 5 * tg / 6]


# ---------------------------------------------------------------------------
# Tasks.
# ---------------------------------------------------------------------------


def _branch_record(
    model: sf.SunModel, parent: sf.Parent, zz: dict[str, Any], wall_s: float
) -> tuple[dict[str, Any], Any]:
    prob, xs0 = sf.forced_problem(model, parent, zz["tau"])
    br = sf.continue_in_eps(prob, xs0, wall_s=wall_s)
    rec: dict[str, Any] = {"zero": zz, "n_seg": prob.n_seg, "branch": br.summary()}
    samples = None
    if br.stop_reason == "reached_target" and br.final_nodes is not None:
        d = sf.orbit_diagnostics(prob, br.final_nodes, 1.0)
        samples = d.pop("_clock_samples")
        rec["diagnostics"] = d
    return rec, samples


def _forced_for_parent(
    model: sf.SunModel, parent: sf.Parent, dtdc: float | None, dtdc_src: str, wall_s: float
) -> dict[str, Any]:
    out: dict[str, Any] = {
        "model": model.kind,
        "parent": {
            "label": parent.label,
            "state": parent.state,
            "period": parent.period,
            "laps_M": parent.laps,
            "sun_periods_N": parent.n_sun,
            "jacobi": sf.jacobi_pv(parent.state, model.mu),
        },
    }
    out["screens"] = sf.screens(model, parent, dtdc=dtdc, dtdc_source=dtdc_src)
    mz = sf.melnikov_zeros(sf.melnikov(model, parent))
    mz.pop("grid")
    out["melnikov"] = mz
    if mz["identically_zero"]:
        out["status"] = "melnikov identically zero: theorem silent, not continued"
        return out
    branches = []
    samples = []
    for zz in mz["zeros"]:
        if not zz["simple"]:
            branches.append({"zero": zz, "skipped": "not a simple zero"})
            continue
        rec, smp = _branch_record(model, parent, zz, wall_s)
        branches.append(rec)
        if smp is not None:
            samples.append((len(branches) - 1, smp))
    if samples:
        cls = sf.symmetry_classes([s for _, s in samples], 128)
        for (i, _), c in zip(samples, cls, strict=True):
            branches[i]["symmetry_class"] = c
    out["branches"] = branches
    out["n_reached"] = sum(1 for b in branches if "diagnostics" in b)
    out["n_classes"] = len(set(b.get("symmetry_class", -1) for b in branches) - {-1})
    out["status"] = "done"
    return out


def control_task(name: str, wall_s: float) -> dict[str, Any]:
    sf.set_threads(2)
    model = sf.qbcp_model(dataclasses.replace(qbcp.qbcp_default(), mu=LB_MU))
    h, y, ydot, p, q = LB_TABLE_1[name]
    x, y, ydot = LB_SECTION_X, -y, -ydot
    r1 = math.hypot(x + LB_MU, y)
    r2 = math.hypot(x - 1.0 + LB_MU, y)
    xdot = math.sqrt(2 * h + x * x + y * y + 2 * (1 - LB_MU) / r1 + 2 * LB_MU / r2 - ydot**2)
    pv, res = sf.correct_parent(model, np.array([x, y, 0.0, xdot, ydot, 0.0]), p / q * model.tg)
    parent = sf.Parent(f"leiva-briozzo-2008-{name}", pv, p / q * model.tg, q, p, residual=res)
    out = _forced_for_parent(model, parent, None, "", wall_s)
    out["parent_residual"] = res
    out["leiva_briozzo_arcs"] = LB_ARCS.get(name, [])
    return out


def family_task(name: str, row_id: str, wall_s: float) -> dict[str, Any]:
    sf.set_threads(2)
    model = sf.bcr4bp_model()
    seed = catalogue_seed(row_id)
    rec: dict[str, Any] = {"family": name, "row": row_id}
    if seed is None or seed["seed"] is None:
        rec["status"] = "no catalogue state or no perpendicular crossing; skipped"
        return rec
    rec["orbit_class"] = seed["orbit_class"]
    members: list[dict[str, Any]] = []
    walks = []
    for direction in (1.0, -1.0):
        found, info = sf.walk_family(
            model,
            np.array(seed["seed"]),
            seed["period"],
            targets(model.tg),
            direction,
            max_steps=600,
            wall_s=wall_s,
            t_window=(max(0.5, seed["period"] - 7), seed["period"] + 7),
        )
        walks.append({"direction": direction, **info})
        for f in found:
            m = f["member"]
            if any(
                abs(m.period - o["period"]) < 1e-9 and abs(m.x0 - o["x0"]) < 1e-7 for o in members
            ):
                continue
            members.append(
                {
                    "x0": m.x0,
                    "z0": m.z0,
                    "vy0": m.vy0,
                    "period": m.period,
                    "residual": m.residual,
                    "dT_dC": f["dT_dC"],
                    "events_before": f["events_before"],
                    "after_bifurcation": bool(f["events_before"]),
                }
            )
    rec["walks"] = walks
    rec["members"] = members
    rec["status"] = "done"
    return rec


def forced_task(family: str, idx: int, member: dict[str, Any], wall_s: float) -> dict[str, Any]:
    sf.set_threads(2)
    out: dict[str, Any] = {"family": family, "member_index": idx, "member": member}
    state = np.array([member["x0"], 0.0, member["z0"], 0.0, member["vy0"], 0.0])
    for model in (sf.bcr4bp_model(), sf.qbcp_model()):
        n, m, _ = sf.commensurability(member["period"], model.tg)
        pv, res = sf.correct_parent(model, state, member["period"])
        parent = sf.Parent(f"{family}#{idx}", pv, member["period"], m, n, residual=res)
        if res > 1e-10:
            out[model.kind] = {"status": f"parent did not correct in this model (res {res:.1e})"}
            continue
        out[model.kind] = _forced_for_parent(model, parent, member["dT_dC"], "family walk", wall_s)
    out["status"] = "done"
    return out


# ---------------------------------------------------------------------------
# Stages.
# ---------------------------------------------------------------------------


def _run_pool(stage: str, jobs: list[tuple[str, Path, Any, tuple[Any, ...]]], workers: int) -> None:
    todo = [j for j in jobs if not j[1].exists()]
    rl = Runlog(OUT / "runlog.jsonl", stage, len(todo))
    log(f"[{stage}] {len(jobs)} tasks, {len(jobs) - len(todo)} already done")
    if not todo:
        return
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = {}
        for name, path, fn, args in todo:
            futs[ex.submit(_timed, fn, args)] = (name, path)
        for fut in as_completed(futs):
            name, path = futs[fut]
            secs, rec, err = fut.result()
            if err is not None:
                rec = {"status": "error", "error": err}
            _write(path, rec)
            rl.task(name, str(rec.get("status", "?"))[:40], secs)


def _timed(fn: Any, args: tuple[Any, ...]) -> tuple[float, Any, str | None]:
    t = time.monotonic()
    try:
        rec = fn(*args)
        return time.monotonic() - t, rec, None
    except Exception:  # the error is recorded in the task file
        return time.monotonic() - t, None, traceback.format_exc()


def stage_controls(workers: int, wall_s: float, only: str | None) -> None:
    jobs = [
        (n, OUT / "controls" / f"{n}.json", control_task, (n, wall_s))
        for n in LB_TABLE_1
        if not only or only in n
    ]
    _run_pool("controls", jobs, workers)


def stage_families(workers: int, wall_s: float, only: str | None) -> None:
    jobs = [
        (f, OUT / "families" / f"{row}.json", family_task, (f, row, wall_s))
        for f, row in FAMILY_SEEDS.items()
        if not only or only in f
    ]
    _run_pool("families", jobs, workers)


def stage_forced(workers: int, wall_s: float, only: str | None, max_members: int | None) -> None:
    jobs = []
    for path in sorted((OUT / "families").glob("*.json")):
        fam = json.loads(path.read_text())
        if fam.get("status") != "done":
            continue
        if only and only not in fam["family"]:
            continue
        for i, m in enumerate(fam["members"]):
            name = f"{fam['row']}#{i}"
            jobs.append(
                (
                    name,
                    OUT / "forced" / f"{fam['row']}__{i}.json",
                    forced_task,
                    (fam["family"], i, m, wall_s),
                )
            )
    if max_members is not None:
        jobs = jobs[:max_members]
    _run_pool("forced", jobs, workers)


def stage_summary() -> None:
    rows = []
    for path in sorted((OUT / "forced").glob("*.json")):
        rec = json.loads(path.read_text())
        if rec.get("status") != "done":
            rows.append({"file": path.name, "status": rec.get("status")})
            continue
        for kind in ("bcr4bp", "qbcp"):
            r = rec.get(kind, {})
            if r.get("status") != "done":
                rows.append({"family": rec["family"], "model": kind, "status": r.get("status")})
                continue
            reached = [b for b in r["branches"] if "diagnostics" in b]
            rows.append(
                {
                    "family": rec["family"],
                    "member": rec["member_index"],
                    "model": kind,
                    "N_over_M": f"{r['parent']['sun_periods_N']}/{r['parent']['laps_M']}",
                    "jacobi": r["parent"]["jacobi"],
                    "screens_pass": r["screens"]["pass"],
                    "after_bifurcation": rec["member"]["after_bifurcation"],
                    "melnikov_amplitude": r["melnikov"]["amplitude"],
                    "zeros": len(r["melnikov"]["zeros"]),
                    "reached": len(reached),
                    "symmetry_classes": r["n_classes"],
                    "stops": [b["branch"]["stop_reason"] for b in r["branches"] if "branch" in b],
                    "folds": [b["branch"]["folds"] for b in r["branches"] if "branch" in b],
                    "max_abs_floquet": [b["diagnostics"]["floquet"]["max_abs"] for b in reached],
                    "periselene_km": [b["diagnostics"]["periselene_km"] for b in reached],
                    "below_moon_surface": any(
                        b["diagnostics"]["below_moon_surface"] for b in reached
                    ),
                    "cycler_class": [b["diagnostics"]["cycler_class"] for b in reached],
                }
            )
    _write(OUT / "summary.json", rows)
    log(f"[summary] {len(rows)} rows -> {OUT / 'summary.json'}")


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument(
        "--stage", required=True, choices=["controls", "families", "forced", "summary", "all"]
    )
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument(
        "--branch-wall-s",
        type=float,
        default=900.0,
        help="wall-clock cap per continuation branch / family walk direction",
    )
    ap.add_argument(
        "--only", default=None, help="families / forced / controls stages: name substring filter"
    )
    ap.add_argument(
        "--max-members",
        type=int,
        default=None,
        help="forced stage: at most this many members (smoke runs)",
    )
    ap.add_argument("--out", default=None, help="output directory (default data/found/905_...)")
    args = ap.parse_args()
    global OUT
    if args.out is not None:
        OUT = Path(args.out)
    preflight_search(
        task_no=905,
        region_id=_REGION_ID,
        method=_METHOD,
        script_path=Path(__file__),
        n_points=len(FAMILY_SEEDS) + len(LB_TABLE_1),
    )
    OUT.mkdir(parents=True, exist_ok=True)
    stages = ["controls", "families", "forced", "summary"] if args.stage == "all" else [args.stage]
    for st in stages:
        t = time.monotonic()
        log(f"stage {st} start")
        if st == "controls":
            stage_controls(args.workers, args.branch_wall_s, args.only)
        elif st == "families":
            stage_families(args.workers, args.branch_wall_s, args.only)
        elif st == "forced":
            stage_forced(args.workers, args.branch_wall_s, args.only, args.max_members)
        else:
            stage_summary()
        log(f"stage {st} done in {time.monotonic() - t:.0f} s")


if __name__ == "__main__":
    main()
