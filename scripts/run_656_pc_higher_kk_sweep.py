"""Task #656: Pluto-Charon higher-(k1,k2) topology extension — runnable script.

`#504`/`#549` only ever swept `(k1,k2)` in {(1,1),(2,1),(2,2),(3,1),(3,2),(3,3)}
at Pluto-Charon's own mass ratio mu=0.10876473603280369
(`cyclerfinder.search.pluto_charon_kk_sweep.PC_MU`, sourced from
`core/cr3bp.cr3bp_system("Pluto","Charon")` — GM_Charon/GM_PlutoSystem via
DE440 planetary constants in `core/satellites.py`). Ross-RT 2026's own
Table I stops at (3,3) too, so nothing beyond that has ever been swept at
PC's mass ratio — the most direct extension of the ONE method that has
already produced a real hit ((3,2), catalogue row `ross-rt-pc-cycler-32-2026`).

This script extends the topology list to (4,1) through (5,5), matching
#549's own enumeration convention exactly: k1 = winding number about the
PRIMARY (Pluto), k2 = winding number about the SECONDARY (Charon), with
k2 <= k1 always (every anchor AND every prior grid target in #504/#549 has
k2 <= k1 — e.g. (3,1),(3,2),(3,3) but never (1,3) or (2,3)). "(4,1) through
(5,5)" is therefore the 9-topology set {(4,1),(4,2),(4,3),(4,4),(5,1),(5,2),
(5,3),(5,4),(5,5)}.

None of these 9 topologies has a Ross-RT 2026 Table-I anchor (the table
stops at (3,3)), so — exactly like #504's own (2,1)/(2,2) and #549's own
(2,1)/(2,2) at other mass ratios — every one of them runs via
`real_binary_kk_sweep.sweep_family_grid`: a bounded (x0, C, hc) brute-force
seed search (SIGALRM-bounded per call) followed by the same C-sweep +
Barden-stability + winding-topology + independent-Radau machinery #504/#549
already validated. Grid box and hc_list are scaled per-topology (see
`_grid_for_topology` below) since higher (k1,k2) orbits wind more and take
longer periods — a single one-size-fits-all grid (as #504/#549 used for
their lower topologies) would systematically miss the higher ones.

Step 1 (mandatory positive control): re-find the admitted Pluto-Charon (3,2)
cycler via #504's OWN `sweep_32_positive_control()`, UNMODIFIED. If this
fails, STOP — do not trust anything else in this script's output.

Usage
-----
  uv run python scripts/run_656_pc_higher_kk_sweep.py --phase init
  uv run python scripts/run_656_pc_higher_kk_sweep.py --phase sweep
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

_SRC = Path(__file__).parent.parent / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

import numpy as np  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402

from cyclerfinder.data.method_capability import MethodCapability  # noqa: E402
from cyclerfinder.data.preflight import preflight_search  # noqa: E402
from cyclerfinder.search.pluto_charon_kk_sweep import (  # noqa: E402
    PC_MU,
    SweepResult,
    _c_l1,
    make_pluto_charon_system,
    sweep_32_positive_control,
)
from cyclerfinder.search.real_binary_kk_sweep import sweep_family_grid  # noqa: E402

OUT_PATH = Path(__file__).parent.parent / "docs" / "notes" / "scratch" / "656_pc_kk45_sweep_raw.txt"

_REGION_ID = "pluto-charon-kk-45-cycler-sweep-2026-07-19"
_METHOD = MethodCapability(
    genome=(
        "Pluto-Charon (k1,k2) CR3BP cycler genome, #504/#549's already-validated "
        "fixed-Jacobi symmetric corrector + Barden stability + winding-topology "
        "classifier + independent Radau crosscheck, extended to the 9 higher "
        "topologies (4,1)-(5,5) never swept before (Ross-RT 2026 Table I and "
        "#504/#549 both stop at (3,3))"
    ),
    corrector=(
        "bounded (x0,C,hc) grid search (sweep_family_grid, same machinery as "
        "#504's own (2,1)/(2,2) grid search) — no Table-I anchor exists for "
        "any of these 9 topologies"
    ),
    capability_tags=frozenset(
        {"cr3bp", "binary-cycler", "k1k2-genome", "real-binary", "grid-search", "pluto-charon"}
    ),
    git_sha="working-tree",
)

# (k1,k2) with k2<=k1, matching #549's own enumeration convention exactly
# (every anchor/grid target #504/#549 ever swept has k2<=k1: (1,1),(2,1),
# (2,2),(3,1),(3,2),(3,3) — never a k2>k1 pair).
TOPOLOGIES: list[tuple[int, int]] = [
    (4, 1),
    (4, 2),
    (4, 3),
    (4, 4),
    (5, 1),
    (5, 2),
    (5, 3),
    (5, 4),
    (5, 5),
]


def _grid_for_topology(k1: int, k2: int) -> tuple[np.ndarray, np.ndarray, tuple[int, ...], float]:
    """Per-topology (x0_grid, c_grid, hc_list, period_guess).

    x0/C ranges follow #504's own (2,1)/(2,2) grid box (x0 in roughly
    [-0.85,-0.30], C in [3.05, C_L1(PC)-0.01]) widened slightly to cover
    both anchors' spans. hc_list and period_guess are SCALED with
    (k1+k2) since higher-winding orbits take more half-crossings and a
    longer period: #504's own sourced correlation is hc ~= (k1+k2) or
    (k1+k2)+1 ((3,2): hc=6, k1+k2=5; (2,1): hc in (3,4), k1+k2=3; (2,2):
    hc in (4,5), k1+k2=4) and period(TU) ~= 2.4*(k1+k2) ((3,2) mid
    T=11.83 TU at k1+k2=5, giving the ratio 11.83/5=2.366 used here).
    """
    c_l1 = _c_l1(PC_MU)
    x0_grid = np.linspace(-0.85, -0.30, 8)
    c_grid = np.linspace(3.05, c_l1 - 0.01, 8)
    s = k1 + k2
    hc_list = tuple(h for h in range(s - 1, s + 3) if h >= 1)
    period_guess = 2.4 * s
    return x0_grid, c_grid, hc_list, period_guess


def _fmt(r: SweepResult) -> str:
    if r.stable_found:
        return (
            f"STABLE  C={r.jacobi_mid:.7f}  x0={r.x0_mid:.9f}  "
            f"T={r.period_mid:.5f} TU ({r.period_days:.3f} d)  "
            f"nu={r.nu_mid:.2e}  topo_ok={r.topology_ok}  xcheck={r.crosscheck_ok}  "
            f"method={r.method!r}"
        )
    return f"negative  method={r.method!r}  note={r.note!r}"


def _run_positive_control() -> tuple[str, bool]:
    t0 = time.time()
    r = sweep_32_positive_control()
    elapsed = time.time() - t0
    ok = r.stable_found and r.topology_ok and r.prograde and r.reaches_secondary and r.crosscheck_ok
    line = f"POSITIVE CONTROL (PC 3,2 @ mu={PC_MU})  {_fmt(r)}  [{elapsed:.1f}s]  PASS={ok}"
    return line, ok


JobResult = tuple[int, int, SweepResult, float]


def _run_job(k1: int, k2: int, *, per_call_timeout: int) -> JobResult:
    pc = make_pluto_charon_system()
    x0_grid, c_grid, hc_list, period_guess = _grid_for_topology(k1, k2)
    t0 = time.time()
    r = sweep_family_grid(
        pc,
        k1,
        k2,
        x0_grid=x0_grid,
        c_grid=c_grid,
        hc_list=hc_list,
        period_guess=period_guess,
        per_call_timeout=per_call_timeout,
    )
    return (k1, k2, r, time.time() - t0)


def _append(lines: list[str]) -> None:
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    with OUT_PATH.open("a") as f:
        f.write("\n".join(lines) + "\n")


def _print_and_append(lines: list[str]) -> None:
    for line in lines:
        print(line, flush=True)
    _append(lines)


def phase_init() -> bool:
    """Fresh header + mandatory positive control. Overwrites OUT_PATH. Returns PASS."""
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    stamp = time.strftime("%Y-%m-%dT%H:%M:%SZ")
    OUT_PATH.write_text(f"Task #656 Pluto-Charon higher-(k1,k2) sweep — run {stamp}\n")

    print("Step 1: MANDATORY positive control (PC 3,2)...", flush=True)
    pc_line, ok = _run_positive_control()
    _print_and_append([pc_line])
    if not ok:
        _print_and_append(["POSITIVE CONTROL FAILED -- STOPPING, do not trust anything else."])
    return ok


def phase_sweep(*, per_call_timeout: int = 3, n_jobs: int = -1) -> None:
    jobs = [delayed(_run_job)(k1, k2, per_call_timeout=per_call_timeout) for k1, k2 in TOPOLOGIES]
    print(f"\nStep 2: dispatching {len(jobs)} topology jobs ({TOPOLOGIES})...", flush=True)
    t0 = time.time()
    for k1, k2 in TOPOLOGIES:
        stamp = time.strftime("%Y-%m-%dT%H:%M:%SZ")
        print(f"  [{stamp}] queued ({k1},{k2})", flush=True)
    results: list[JobResult] = Parallel(n_jobs=n_jobs, verbose=10)(jobs)
    elapsed = time.time() - t0
    print(f"all {len(jobs)} topology jobs done in {elapsed:.1f}s", flush=True)
    lines = [f"\nResults (elapsed {elapsed:.1f}s)"]
    for k1, k2, r, job_elapsed in sorted(results, key=lambda t: (t[0], t[1])):
        lines.append(f"({k1},{k2})  {_fmt(r)}  [{job_elapsed:.1f}s]")
    _print_and_append(lines)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--phase", choices=["init", "sweep"], required=True)
    ap.add_argument("--per-call-timeout", type=int, default=3)
    ap.add_argument("--n-jobs", type=int, default=-1)
    args = ap.parse_args()

    preflight_search(
        task_no=656,
        region_id=_REGION_ID,
        method=_METHOD,
        script_path=Path(__file__),
        n_points=len(TOPOLOGIES) * 8 * 8 * 4,
        override_reason=(
            "reuses #504/#549's already-validated binary-cycler harness verbatim "
            "(mandatory positive control re-finds the committed PC (3,2) row before "
            "trusting anything else); the (topology) job count is small and bounded "
            "by construction (9 topologies x an 8x8x<=4 grid), not an open-ended "
            "discovery grid needing a timing pilot."
        ),
    )

    if args.phase == "init":
        ok = phase_init()
        if not ok:
            raise SystemExit(1)
    elif args.phase == "sweep":
        phase_sweep(per_call_timeout=args.per_call_timeout, n_jobs=args.n_jobs)

    print(f"\n[phase={args.phase}] appended to {OUT_PATH}")


if __name__ == "__main__":
    main()
