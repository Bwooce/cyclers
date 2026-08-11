"""Task #810: Pluto-Charon (5,1) fixed-hc C-sweep re-attempt.

`#656`'s higher-(k1,k2) grid sweep found a genuine prograde (5,1) seed at
x0=-0.6685146994, C=3.05 (its own grid FLOOR), T=28.1427 TU — but the
downstream `c_sweep_find_nu_zero` was called with `hc=None`, auto-redetecting
the half-crossing index at every C step, and walked OFF the genuine (5,1)
branch onto the unrelated retrograde (4,0) family before finding a nu=0
crossing there (`#656`'s own hand diagnosis, replayed deterministically by
`#808`). `#504`'s (3,2) positive control avoids exactly this failure mode by
sweeping UPWARD ONLY with the half-crossing count held FIXED through the
whole sweep. `#810` re-runs the (5,1) topology under that same convention.

This is the last asterisk on the #504/#549/#656 15-topology Pluto-Charon
census: (5,1) is the one k2<=k1<=5 topology whose negative rests on a known
search-method gap (the hc=None branch loss) rather than an
exhausted-within-budget seed search.

Design (all phases run FOREGROUND, checkpointed per step, resumable):

1. ``--phase control``: mandatory positive control — `#504`'s own
   ``sweep_32_positive_control()``, UNMODIFIED. If this fails, STOP.
2. ``--phase seed``: reconverge the recorded (5,1) seed and verify its
   topology INDEPENDENTLY (``winding_topology``: k1=5, k2=1, prograde,
   reaches_secondary) — do not trust the recorded claim — and MEASURE the
   seed's own half-crossing index (the perpendicular crossing nearest T/2),
   which the sweep then holds fixed.
3. ``--phase up``: instrumented fixed-hc C-sweep upward from the seed's C to
   C_L1(PC)-0.002 (the #504 convention). Every step records (C, x0, T,
   Barden nu, winding topology) to a JSONL checkpoint; stability is judged
   per-step by |nu|<1 directly (strictly stronger than the shared
   machinery's nu=0-sign-change detection, which misses a stable window
   whose nu never changes sign). Branch identity is CHECKED at every step —
   the exact thing the hc=None path could not do.
4. ``--phase down``: bounded downward diagnostic below the seed's C
   (BEYOND the #504 upward-only convention, labelled as such): the seed sat
   AT `#656`'s grid floor C=3.05, so the family plausibly extends lower.
5. ``--phase verdict``: read both checkpoints, print the verdict; with
   ``--stamp``, append the empty-region record iff the verdict is a clean
   negative (refuses to stamp if any gate-passing stable member was found).

Usage
-----
  uv run python scripts/run_810_pc_51_fixed_hc_sweep.py --phase control
  uv run python scripts/run_810_pc_51_fixed_hc_sweep.py --phase seed
  uv run python scripts/run_810_pc_51_fixed_hc_sweep.py --phase up
  uv run python scripts/run_810_pc_51_fixed_hc_sweep.py --phase down
  uv run python scripts/run_810_pc_51_fixed_hc_sweep.py --phase verdict [--stamp]

``--phase up``/``down`` resume from their checkpoint file automatically, so
they can be run in bounded foreground chunks.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

_SRC = Path(__file__).parent.parent / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

import numpy as np  # noqa: E402
from scipy.optimize import brentq  # noqa: E402

import cyclerfinder.search.cr3bp_periodic as cp  # noqa: E402
from cyclerfinder.data.empty_regions import (  # noqa: E402
    DEFAULT_EMPTY_REGIONS_PATH,
    EmptyRegionReport,
    append_empty_region,
)
from cyclerfinder.data.method_capability import MethodCapability  # noqa: E402
from cyclerfinder.data.preflight import preflight_search  # noqa: E402
from cyclerfinder.search.binary_star_search import winding_topology  # noqa: E402
from cyclerfinder.search.pluto_charon_kk_sweep import (  # noqa: E402
    PC_MU,
    _build_result,
    _c_l1,
    _run_with_timeout,
    make_pluto_charon_system,
    sweep_32_positive_control,
)

# #656's recorded genuine prograde (5,1) grid seed (data/OUTSTANDING.md, #656
# bullet; independently re-verified by --phase seed before any sweep trusts it).
SEED_X0 = -0.6685146994
SEED_C = 3.05
SEED_T = 28.1427

#: Sweep step. #504's own convention is n_coarse=60 over the family's C band
#: (dC ~ 0.01 for a 0.6-wide band); dC=0.005 over the ~0.57-wide [3.05,
#: C_L1-0.002] band is 2x finer than that convention, not coarser.
DC = 0.005

#: Downward-diagnostic floor (beyond-convention; the seed sat AT #656's grid
#: floor C=3.05, so the family plausibly extends below it).
C_DOWN_FLOOR = 2.90

#: Per-corrector-call SIGALRM bound (s). The T~28 TU (5,1) orbit is ~2.4x
#: longer than the (3,2) positive control's; #656's grid used 3-4 s calls.
PER_CALL_TIMEOUT = 30

#: Consecutive failed/timed-out correction steps before declaring the branch
#: ended (the shared machinery just silently `continue`s forever; we bound it).
MAX_CONSEC_FAIL = 5

OUT_DIR = Path(__file__).parent.parent / "data" / "found" / "810_pc51_fixed_hc_sweep"
RAW_PATH = Path(__file__).parent.parent / "docs" / "notes" / "scratch" / "810_pc51_fixedhc_raw.txt"

_REGION_ID = "pluto-charon-51-fixed-hc-cycler-sweep-2026-08-11"
_METHOD = MethodCapability(
    genome=(
        "Pluto-Charon (5,1) CR3BP cycler family, #504/#549/#656's validated "
        "fixed-Jacobi symmetric corrector + Barden stability + winding-topology "
        "classifier + independent Radau crosscheck, re-swept from #656's genuine "
        "prograde (5,1) grid seed with the half-crossing index held FIXED "
        "through the whole C-sweep (#504's (3,2) positive-control convention) "
        "instead of the hc=None auto-redetection that lost the branch in #656; "
        "per-step branch-identity (winding) verification + per-step |nu|<1 "
        "stability judgement (strictly stronger than nu=0 sign-change "
        "detection); bounded downward extension below the seed's grid-floor C"
    ),
    corrector=(
        "correct_symmetric_fixed_jacobi with half_crossings fixed at the "
        "seed's own measured perpendicular-crossing index, SIGALRM-bounded "
        "per call, instrumented per-step checkpointing"
    ),
    capability_tags=frozenset(
        {
            "cr3bp",
            "binary-cycler",
            "k1k2-genome",
            "real-binary",
            "grid-search",
            "pluto-charon",
            "fixed-hc-branch-tracking",
        }
    ),
    git_sha="working-tree",
)


def _log(lines: list[str]) -> None:
    RAW_PATH.parent.mkdir(parents=True, exist_ok=True)
    with RAW_PATH.open("a") as f:
        f.write("\n".join(lines) + "\n")
    for line in lines:
        print(line, flush=True)


def _stamp() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%S%z")


# ---------------------------------------------------------------------------
# Phase: control
# ---------------------------------------------------------------------------


def phase_control() -> bool:
    """Mandatory positive control: #504's own (3,2) sweep, unmodified."""
    t0 = time.time()
    r = sweep_32_positive_control()
    ok = bool(
        r.stable_found and r.topology_ok and r.prograde and r.reaches_secondary and r.crosscheck_ok
    )
    _log(
        [
            f"[{_stamp()}] POSITIVE CONTROL (PC 3,2): stable={r.stable_found} "
            f"C={r.jacobi_mid} x0={r.x0_mid} T={r.period_mid} nu={r.nu_mid} "
            f"topo_ok={r.topology_ok} xcheck={r.crosscheck_ok} "
            f"[{time.time() - t0:.1f}s] PASS={ok}"
        ]
    )
    if not ok:
        _log(["POSITIVE CONTROL FAILED -- STOP, do not trust anything else."])
    return ok


# ---------------------------------------------------------------------------
# Phase: seed  (reconverge + independent topology verification + measure hc)
# ---------------------------------------------------------------------------


def measure_half_crossing_index(
    system: Any, orbit: cp.SymmetricOrbit
) -> tuple[int, int, np.ndarray]:
    """Measure the 1-based index of the perpendicular crossing nearest T/2.

    Returns (hc_index, n_crossings_within_1p25T, crossing_times). This is the
    index ``correct_symmetric_fixed_jacobi`` fixes internally when called with
    ``half_crossings=None`` on this orbit — measured here explicitly so the
    sweep can hold it fixed (the #504 convention #656's path lacked).
    """
    state0 = np.array([orbit.x0, 0.0, 0.0, 0.0, orbit.ydot0, 0.0])
    times, _states = cp._xaxis_crossings(
        system, state0, 1.25 * orbit.period, with_stm=False, rtol=1e-12, atol=1e-12
    )
    if len(times) == 0:
        raise RuntimeError("no x-axis crossings found on the reconverged seed orbit")
    hc = int(np.argmin(np.abs(times - 0.5 * orbit.period))) + 1
    return hc, len(times), times


def reconverge_seed() -> tuple[cp.SymmetricOrbit, int, Any]:
    """Reconverge #656's recorded (5,1) seed and verify it independently."""
    sys_pc = make_pluto_charon_system()
    orbit = cp.correct_symmetric_fixed_jacobi(
        sys_pc,
        SEED_X0,
        SEED_C,
        SEED_T,
        ydot0_sign=-1.0,
        half_crossings=None,
        tol=1e-10,
    )
    if not orbit.converged:
        raise RuntimeError(
            f"recorded (5,1) seed failed to reconverge (residual={orbit.crossing_residual:.2e})"
        )
    state0 = np.array([orbit.x0, 0.0, 0.0, 0.0, orbit.ydot0, 0.0])
    topo = winding_topology(sys_pc.mu, state0, orbit.period)
    hc, n_cross, _times = measure_half_crossing_index(sys_pc, orbit)
    _log(
        [
            f"[{_stamp()}] SEED reconverged: x0={orbit.x0:.10f} C={orbit.jacobi:.10f} "
            f"T={orbit.period:.5f} TU  residual={orbit.crossing_residual:.2e}",
            f"[{_stamp()}] SEED topology: (k1,k2)=({topo.k1},{topo.k2}) "
            f"w1={topo.w1:+.4f} w2={topo.w2:+.4f} prograde={topo.prograde} "
            f"reaches_secondary={topo.reaches_secondary}",
            f"[{_stamp()}] SEED half-crossing index (nearest T/2): hc={hc} "
            f"(of {n_cross} crossings within 1.25T)",
        ]
    )
    if not (topo.k1 == 5 and topo.k2 == 1 and topo.prograde and topo.reaches_secondary):
        raise RuntimeError(
            "reconverged seed is NOT a prograde (5,1) orbit -- recorded claim "
            f"failed independent verification: ({topo.k1},{topo.k2}), "
            f"prograde={topo.prograde}, reaches_secondary={topo.reaches_secondary}"
        )
    return orbit, hc, topo


# ---------------------------------------------------------------------------
# Instrumented fixed-hc sweep (checkpointed, resumable)
# ---------------------------------------------------------------------------


def _correct_bounded(
    system: Any, x0: float, c: float, t_guess: float, hc: int
) -> cp.SymmetricOrbit | None:
    """One SIGALRM-bounded fixed-hc correction; None on timeout/failure."""

    def _fn() -> cp.SymmetricOrbit:
        return cp.correct_symmetric_fixed_jacobi(
            system,
            x0,
            c,
            t_guess,
            ydot0_sign=-1.0,
            half_crossings=hc,
            tol=1e-10,
        )

    try:
        o = _run_with_timeout(_fn, seconds=PER_CALL_TIMEOUT)
    except (ValueError, RuntimeError):
        return None
    if o is None or not o.converged:
        return None
    return o


def _step_record(system: Any, c: float, orbit: cp.SymmetricOrbit, wall_s: float) -> dict[str, Any]:
    nu, _lam = cp.barden_stability(system, orbit, rtol=1e-13, atol=1e-13)
    state0 = np.array([orbit.x0, 0.0, 0.0, 0.0, orbit.ydot0, 0.0])
    topo = winding_topology(system.mu, state0, orbit.period)
    return {
        "c": float(c),
        "converged": True,
        "x0": float(orbit.x0),
        "ydot0": float(orbit.ydot0),
        "period": float(orbit.period),
        "nu": float(nu),
        "stable": bool(abs(nu) < 1.0),
        "k1": int(topo.k1),
        "k2": int(topo.k2),
        "w1": float(topo.w1),
        "w2": float(topo.w2),
        "prograde": bool(topo.prograde),
        "reaches_secondary": bool(topo.reaches_secondary),
        "topo_ok": bool(topo.k1 == 5 and topo.k2 == 1 and topo.prograde),
        "wall_s": round(wall_s, 2),
    }


def _load_ckpt(path: Path) -> list[dict[str, Any]]:
    if not path.exists():
        return []
    return [json.loads(line) for line in path.read_text().splitlines() if line.strip()]


def _append_ckpt(path: Path, rec: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a") as f:
        f.write(json.dumps(rec) + "\n")


def _refine_nu_zero(
    system: Any,
    hc: int,
    prev: dict[str, Any],
    cur: dict[str, Any],
) -> dict[str, Any] | None:
    """Brentq-refine a nu sign change between two topology-verified steps.

    hc stays FIXED through the refinement (the whole point of #810). Returns a
    full gate-check record at the nu=0 midpoint, or None if brentq fails.
    """
    x0_b, t_b = prev["x0"], prev["period"]

    def _nu_at(c_val: float) -> float:
        o = _correct_bounded(system, x0_b, c_val, t_b, hc)
        if o is None:
            return float("nan")
        nu, _ = cp.barden_stability(system, o, rtol=1e-13, atol=1e-13)
        return float(nu)

    try:
        c_mid = brentq(_nu_at, prev["c"], cur["c"], xtol=1e-10, rtol=1e-10, maxiter=60)
    except ValueError:
        return None
    o_mid = _correct_bounded(system, x0_b, float(c_mid), t_b, hc)
    if o_mid is None:
        return None
    res = _build_result(5, 1, system, o_mid, "fixed_hc_c_sweep_brentq_midpoint")
    return {
        "c_mid": float(c_mid),
        "x0": res.x0_mid,
        "ydot0": res.ydot0_mid,
        "period": res.period_mid,
        "period_days": res.period_days,
        "nu": res.nu_mid,
        "stable_found": res.stable_found,
        "topology_ok": res.topology_ok,
        "prograde": res.prograde,
        "reaches_secondary": res.reaches_secondary,
        "crosscheck_ok": res.crosscheck_ok,
        "crosscheck_dj": res.crosscheck_dj,
    }


def run_sweep(direction: str, budget_s: float) -> None:
    """Instrumented fixed-hc sweep, resumable from its JSONL checkpoint.

    ``direction`` is "up" (seed C -> C_L1-0.002, the #504 convention) or
    "down" (seed C -> C_DOWN_FLOOR, beyond-convention diagnostic). Stops
    cleanly when ``budget_s`` elapses so callers can re-invoke to resume.
    """
    sys_pc = make_pluto_charon_system()
    c_l1 = _c_l1(PC_MU)

    seed_orbit, hc, _topo = reconverge_seed()

    if direction == "up":
        c_values = np.arange(SEED_C, c_l1 - 0.002 + 1e-12, DC)
        ckpt = OUT_DIR / "sweep_up.jsonl"
    else:
        c_values = np.arange(SEED_C, C_DOWN_FLOOR - 1e-12, -DC)
        ckpt = OUT_DIR / "sweep_down.jsonl"

    done = _load_ckpt(ckpt)
    done_steps = [r for r in done if "c" in r]
    n_done = len(done_steps)
    # Carriers: last converged record, else the seed.
    x0_cur, t_cur = seed_orbit.x0, seed_orbit.period
    consec_fail = 0
    for rec in done_steps:
        if rec.get("converged"):
            x0_cur, t_cur = rec["x0"], rec["period"]
            consec_fail = 0
        else:
            consec_fail += 1
    if any(rec.get("event") == "end" for rec in done):
        _log([f"[{_stamp()}] sweep {direction}: already complete ({n_done} steps), nothing to do"])
        return

    _log(
        [
            f"[{_stamp()}] sweep {direction}: hc FIXED at {hc}, "
            f"{len(c_values)} C steps total, resuming at step {n_done}, "
            f"budget {budget_s:.0f}s"
        ]
    )
    t_start = time.time()
    prev_ok: dict[str, Any] | None = None
    for rec in reversed(done_steps):
        if rec.get("converged") and rec.get("topo_ok"):
            prev_ok = rec
            break

    for i in range(n_done, len(c_values)):
        if time.time() - t_start > budget_s:
            _log([f"[{_stamp()}] sweep {direction}: budget reached at step {i}, resume later"])
            return
        c = float(c_values[i])
        t0 = time.time()
        o = _correct_bounded(sys_pc, x0_cur, c, t_cur, hc)
        if o is None:
            consec_fail += 1
            rec = {
                "c": c,
                "converged": False,
                "consec_fail": consec_fail,
                "wall_s": round(time.time() - t0, 2),
            }
            _append_ckpt(ckpt, rec)
            _log([f"[{_stamp()}]   step {i} C={c:.5f}: FAILED ({consec_fail} consecutive)"])
            if consec_fail >= MAX_CONSEC_FAIL:
                _append_ckpt(ckpt, {"event": "end", "reason": f"{MAX_CONSEC_FAIL} consecutive"})
                _log([f"[{_stamp()}] sweep {direction}: branch ended (consecutive failures)"])
                return
            continue
        consec_fail = 0
        rec = _step_record(sys_pc, c, o, time.time() - t0)
        _append_ckpt(ckpt, rec)
        _log(
            [
                f"[{_stamp()}]   step {i} C={c:.5f}: x0={rec['x0']:.8f} "
                f"T={rec['period']:.4f} nu={rec['nu']:+.4e} stable={rec['stable']} "
                f"({rec['k1']},{rec['k2']}) prograde={rec['prograde']} "
                f"topo_ok={rec['topo_ok']} [{rec['wall_s']}s]"
            ]
        )
        if not rec["topo_ok"]:
            _append_ckpt(
                ckpt,
                {
                    "event": "end",
                    "reason": (
                        f"branch identity lost DESPITE fixed hc={hc}: recovered "
                        f"({rec['k1']},{rec['k2']}) w1={rec['w1']:+.3f} "
                        f"w2={rec['w2']:+.3f} at C={c:.5f}"
                    ),
                },
            )
            _log([f"[{_stamp()}] sweep {direction}: branch identity lost at C={c:.5f} -- STOP"])
            return
        if prev_ok is not None and prev_ok["nu"] * rec["nu"] < 0.0:
            _log([f"[{_stamp()}]   nu SIGN CHANGE in [{prev_ok['c']:.5f},{c:.5f}] -- refining"])
            mid = _refine_nu_zero(sys_pc, hc, prev_ok, rec)
            _append_ckpt(ckpt, {"event": "nu_zero", "bracket": [prev_ok["c"], c], "mid": mid})
            _log([f"[{_stamp()}]   nu=0 midpoint: {json.dumps(mid)}"])
        prev_ok = rec
        x0_cur, t_cur = rec["x0"], rec["period"]

    _append_ckpt(ckpt, {"event": "end", "reason": "C range exhausted"})
    _log([f"[{_stamp()}] sweep {direction}: COMPLETE ({len(c_values)} steps)"])


# ---------------------------------------------------------------------------
# Phase: verdict (+ optional empty-region stamp)
# ---------------------------------------------------------------------------


def summarize() -> dict[str, Any]:
    """Aggregate both checkpoints into the verdict summary."""
    out: dict[str, Any] = {"up": {}, "down": {}}
    for direction in ("up", "down"):
        recs = _load_ckpt(OUT_DIR / f"sweep_{direction}.jsonl")
        steps = [r for r in recs if "c" in r]
        conv = [r for r in steps if r.get("converged")]
        topo_ok = [r for r in conv if r.get("topo_ok")]
        stable = [r for r in topo_ok if r.get("stable")]
        nu_zero = [r for r in recs if r.get("event") == "nu_zero"]
        ends = [r for r in recs if r.get("event") == "end"]
        out[direction] = {
            "n_steps": len(steps),
            "n_converged": len(conv),
            "n_topo_ok": len(topo_ok),
            "n_stable_topo_ok": len(stable),
            "stable_records": stable,
            "nu_zero_events": nu_zero,
            "end_reason": ends[0]["reason"] if ends else "INCOMPLETE",
            "c_range_converged": (
                [min(r["c"] for r in conv), max(r["c"] for r in conv)] if conv else None
            ),
            "nu_range_topo_ok": (
                [min(r["nu"] for r in topo_ok), max(r["nu"] for r in topo_ok)] if topo_ok else None
            ),
            "min_abs_nu_topo_ok": (min(abs(r["nu"]) for r in topo_ok) if topo_ok else None),
        }
    return out


def phase_verdict(stamp: bool) -> None:
    s = summarize()
    _log([f"[{_stamp()}] VERDICT SUMMARY:", json.dumps(s, indent=2)])
    n_stable = s["up"]["n_stable_topo_ok"] + s["down"]["n_stable_topo_ok"]
    incomplete = "INCOMPLETE" in (s["up"]["end_reason"], s["down"]["end_reason"])
    if incomplete:
        _log(["VERDICT: sweeps incomplete -- run --phase up/down to completion first"])
        return
    if n_stable > 0:
        _log(
            [
                f"VERDICT: {n_stable} stable (|nu|<1) topology-verified (5,1) member(s) "
                "found -- POSITIVE CANDIDATE. Mandatory next steps: independent "
                "reproduction, literature_check.py novelty gate, adjudication. "
                "Do NOT stamp an empty region; do NOT write the catalogue from here."
            ]
        )
        return
    _log(["VERDICT: clean negative -- family exists but no stable member in the swept C range"])
    if not stamp:
        return
    git_sha = subprocess.run(
        ["git", "rev-parse", "--short", "HEAD"],
        capture_output=True,
        text=True,
        cwd=Path(__file__).parent.parent,
        check=False,
    ).stdout.strip()
    total_steps = s["up"]["n_steps"] + s["down"]["n_steps"]
    wall = 0.0
    for direction in ("up", "down"):
        for r in _load_ckpt(OUT_DIR / f"sweep_{direction}.jsonl"):
            wall += r.get("wall_s", 0.0)
    report = EmptyRegionReport(
        region_id=_REGION_ID,
        family="pluto-charon-(5,1)-cr3bp-cycler",
        centre="Pluto-Charon barycentre (CR3BP, mu=0.10876473603280369)",
        topologies=({"k1": 5, "k2": 1, "prograde": True},),
        method_capability=_METHOD,
        search_extent={
            "points_total": total_steps,
            "c_up": s["up"]["c_range_converged"],
            "c_down": s["down"]["c_range_converged"],
            "dc": DC,
            "hc_fixed": True,
            "seed": {"x0": SEED_X0, "c": SEED_C, "period_tu": SEED_T},
            "not_covered": (
                "C below the down-sweep's own convergence floor; retrograde "
                "(5,1); non-symmetric (non-perpendicular-crossing) members of "
                "the family; 3D members"
            ),
        },
        prune_gates=(
            "converged corrector (fixed hc)",
            "winding topology (5,1) prograde per step",
            "Barden |nu|<1 per step",
            "independent Radau crosscheck (on any nu=0 midpoint)",
        ),
        result={
            "up": {k: v for k, v in s["up"].items() if k != "stable_records"},
            "down": {k: v for k, v in s["down"].items() if k != "stable_records"},
        },
        verdict=(
            "clean negative: the genuine prograde (5,1) family exists at PC mu "
            "(seed independently re-verified) and the fixed-hc sweep held its "
            "branch identity, but no member with Barden |nu|<1 exists in the "
            "swept C range -- closing the #656 (5,1) search-method-gap asterisk "
            "on the #504/#549/#656 15-topology census"
        ),
        interpretation=(
            "#810: #656's hc=None auto-redetection branch loss was the ONLY "
            "reason (5,1) stayed UNSETTLED; re-sweeping with the seed's own "
            "measured half-crossing index held fixed (per #504's (3,2) "
            "positive-control convention), verifying winding topology at every "
            "step, and judging stability per-step by |nu|<1 directly (stronger "
            "than nu=0 sign-change detection) finds the family has no stable "
            "member. Conditional on the planar CR3BP, symmetric "
            "(perpendicular-crossing) family members, and the swept C range. "
            "RESWEEP CONDITION: a method covering non-symmetric or 3D family "
            "members, or extending convergence below the down-sweep floor, "
            "could subsume this."
        ),
        source_anchors=(
            "scripts/run_810_pc_51_fixed_hc_sweep.py + "
            "data/found/810_pc51_fixed_hc_sweep/sweep_up.jsonl + sweep_down.jsonl "
            "(full per-step record, reproducible); "
            "docs/notes/2026-08-11-810-pc-51-fixed-hc-sweep.md; #656/#808 "
            "diagnosis in data/OUTSTANDING.md + "
            "docs/notes/2026-08-09-808-real-binary-grid-topology-gate.md; prior "
            "region pluto-charon-kk-45-cycler-sweep-2026-07-19 (the (5,1) row "
            "this supersedes)"
        ),
        run={
            "date": time.strftime("%Y-%m-%d"),
            "git_sha": git_sha,
            "task": 810,
            "wall_s": round(wall, 1),
            "note": (
                "foreground chunked resumable run; #504 (3,2) positive control "
                "passed before any sweep; no catalogue write"
            ),
        },
    )
    append_empty_region(DEFAULT_EMPTY_REGIONS_PATH, report)
    _log([f"[{_stamp()}] stamped empty region {_REGION_ID} to {DEFAULT_EMPTY_REGIONS_PATH}"])


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--phase", choices=["control", "seed", "up", "down", "verdict"], required=True)
    ap.add_argument("--budget-s", type=float, default=420.0, help="per-invocation sweep budget")
    ap.add_argument("--stamp", action="store_true", help="verdict phase: stamp empty region")
    args = ap.parse_args()

    n_up = int(np.ceil((_c_l1(PC_MU) - 0.002 - SEED_C) / DC)) + 1
    n_down = int(np.ceil((SEED_C - C_DOWN_FLOOR) / DC)) + 1
    preflight_search(
        task_no=810,
        region_id=_REGION_ID,
        method=_METHOD,
        script_path=Path(__file__),
        n_points=n_up + n_down,
        override_reason=(
            "#810 IS registered in data/OUTSTANDING.md (its bullet uses the "
            "backtick list style, not the bold style the allocation regex "
            "matches -- the same false positive #656 recorded); the region "
            "intentionally overlaps pluto-charon-kk-45-cycler-sweep-2026-07-19 "
            "because THIS method (fixed-hc branch tracking) is strictly more "
            "capable on the (5,1) row that region itself flags as UNSETTLED "
            "due to the hc=None branch loss; grid is small and bounded "
            "(~145 C steps), reusing the #504/#549/#656-validated harness with "
            "a mandatory positive control."
        ),
    )

    if args.phase == "control":
        if not phase_control():
            raise SystemExit(1)
    elif args.phase == "seed":
        reconverge_seed()
    elif args.phase in ("up", "down"):
        run_sweep(args.phase, args.budget_s)
    elif args.phase == "verdict":
        phase_verdict(args.stamp)


if __name__ == "__main__":
    main()
