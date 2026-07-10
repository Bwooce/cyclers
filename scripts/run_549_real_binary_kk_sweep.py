"""Task #549: real-binary (k1,k2)-cycler genome sweep — runnable script.

Reuses the #494/#504 binary-cycler CR3BP genome (fixed-Jacobi symmetric
corrector + Barden stability + winding-topology classifier + independent
Radau crosscheck) verbatim via
:mod:`cyclerfinder.search.real_binary_kk_sweep`, which only generalizes
#504's mu-continuation driver to finish in an arbitrary target system
instead of a hardcoded Pluto-Charon system.

Step 1 (mandatory positive control): re-find the admitted Pluto-Charon
(3,2) cycler via #504's OWN `sweep_32_positive_control()`, unmodified.

Step 2: sweep (k1,k2) in {(1,1),(2,1),(2,2),(3,1),(3,2),(3,3)} — the same
six topologies #504 swept — at four sourced real-binary mass ratios:
Patroclus-Menoetius, Didymos-Dimorphos, Orcus-Vanth, Eris-Dysnomia (see
`REAL_BINARY_SYSTEMS` for citations). Anchor-seeded families run via
mu-continuation from the Ross-RT 2026 Table-I anchors; (2,1)/(2,2) (not in
Table-I) run via a bounded (x0,C,hc) grid search with a per-call SIGALRM
timeout, exactly as #504 did.

Runs the (system, topology) jobs in parallel (joblib) since they are fully
independent. Output: a per-system/topology result table to stdout AND
incrementally appended to docs/notes/scratch/549_kk_sweep_raw.txt (so a
partial run is never silently lost).

Usage
-----
  uv run python scripts/run_549_real_binary_kk_sweep.py
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
from cyclerfinder.search.pluto_charon_kk_sweep import PC_MU, sweep_32_positive_control  # noqa: E402
from cyclerfinder.search.real_binary_kk_sweep import (  # noqa: E402
    REAL_BINARY_SYSTEMS,
    SweepResult,
    sweep_family,
    sweep_family_grid,
)

OUT_PATH = Path(__file__).parent.parent / "docs" / "notes" / "scratch" / "549_kk_sweep_raw.txt"

_REGION_ID = "real-binary-kk-cycler-sweep-2026-07-10"
_METHOD = MethodCapability(
    genome=(
        "Real-binary (k1,k2) CR3BP cycler genome (#494/#504's fixed-Jacobi symmetric "
        "corrector + Barden stability + winding-topology classifier + independent Radau "
        "crosscheck), swept at four sourced real-binary mass ratios (Patroclus-Menoetius, "
        "Didymos-Dimorphos, Orcus-Vanth, Eris-Dysnomia) across 6 (k1,k2) topologies"
    ),
    corrector=(
        "anchor-seeded mu-continuation from Ross-RT 2026 Table-I anchors for "
        "(1,1)/(3,1)/(3,2)/(3,3); bounded (x0,C,hc) grid search for (2,1)/(2,2)"
    ),
    capability_tags=frozenset(
        {"cr3bp", "binary-cycler", "k1k2-genome", "real-binary", "mu-continuation"}
    ),
    git_sha="working-tree",
)

# anchor_key per (k1,k2) for the anchor-seeded families; (2,1)/(2,2) have no
# Table-I anchor (per #504) and run via grid search instead.
ANCHOR_TOPOLOGIES: list[tuple[int, int, str]] = [
    (1, 1, "mu001_11"),
    (1, 1, "mu01215_11"),
    (1, 1, "mu05_11"),
    (3, 2, "mu01_32"),
    (3, 1, "mu03_31"),
    (3, 3, "mu01215_33"),
]
GRID_TOPOLOGIES: list[tuple[int, int]] = [(2, 1), (2, 2)]


def _fmt(r: SweepResult) -> str:
    if r.stable_found:
        return (
            f"STABLE  C={r.jacobi_mid:.7f}  x0={r.x0_mid:.9f}  "
            f"T={r.period_mid:.5f} TU ({r.period_days:.3f} d)  "
            f"nu={r.nu_mid:.2e}  topo_ok={r.topology_ok}  xcheck={r.crosscheck_ok}  "
            f"method={r.method!r}"
        )
    return f"negative  method={r.method!r}  note={r.note!r}"


def _run_positive_control() -> str:
    t0 = time.time()
    r = sweep_32_positive_control()
    elapsed = time.time() - t0
    line = f"POSITIVE CONTROL (PC 3,2 @ mu={PC_MU})  {_fmt(r)}  [{elapsed:.1f}s]"
    return line


JobResult = tuple[str, int, int, str, SweepResult, float]


def _run_anchor_job(sys_key: str, k1: int, k2: int, anchor_key: str) -> JobResult:
    target = REAL_BINARY_SYSTEMS[sys_key].to_cr3bp_system()
    t0 = time.time()
    r = sweep_family(target, anchor_key)
    return (sys_key, k1, k2, anchor_key, r, time.time() - t0)


def _run_grid_job(sys_key: str, k1: int, k2: int) -> JobResult:
    target = REAL_BINARY_SYSTEMS[sys_key].to_cr3bp_system()
    mu = target.mu
    # Grid scaled to bracket the whole region between beyond-P1 and
    # near-P2 for ANY mu in [0.001, 0.5] (broader than #504's PC-specific
    # box, since these target mu values span a much wider range).
    x0_grid = np.linspace(-1.05, (1.0 - mu) - 0.05, 8)
    from cyclerfinder.search.pluto_charon_kk_sweep import _c_l1

    c_l1 = _c_l1(mu)
    c_grid = np.linspace(max(2.6, c_l1 - 1.0), c_l1 - 0.01, 6)
    hc_list = (2, 3, 4)
    t0 = time.time()
    r = sweep_family_grid(
        target,
        k1,
        k2,
        x0_grid=x0_grid,
        c_grid=c_grid,
        hc_list=hc_list,
        period_guess=14.0,
        per_call_timeout=2,
    )
    return (sys_key, k1, k2, "grid_search", r, time.time() - t0)


def _append(lines: list[str]) -> None:
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    with OUT_PATH.open("a") as f:
        f.write("\n".join(lines) + "\n")


def _print_and_append(lines: list[str]) -> None:
    for line in lines:
        print(line, flush=True)
    _append(lines)


def phase_init() -> None:
    """Fresh header + positive control + sourced-system table. Overwrites OUT_PATH."""
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    stamp = time.strftime("%Y-%m-%dT%H:%M:%SZ")
    OUT_PATH.write_text(f"Task #549 real-binary (k1,k2) sweep — run {stamp}\n")

    print("Step 1: positive control (PC 3,2)...", flush=True)
    pc_line = _run_positive_control()
    _print_and_append([pc_line])

    print("\nSourced systems:")
    for key, s in REAL_BINARY_SYSTEMS.items():
        lines = [
            f"SYSTEM {key}: mu={s.mu:.10f} l_km={s.l_km} t_s={s.t_s:.4f}",
            f"  mu_source: {s.mu_source}",
            f"  l_source: {s.l_source}",
            f"  t_source: {s.t_source}",
        ]
        if s.caveat:
            lines.append(f"  CAVEAT: {s.caveat}")
        _print_and_append(lines)


def phase_anchors(system_keys: list[str]) -> None:
    jobs = [
        delayed(_run_anchor_job)(sys_key, k1, k2, anchor_key)
        for sys_key in system_keys
        for k1, k2, anchor_key in ANCHOR_TOPOLOGIES
    ]
    print(
        f"\nStep 2a (anchor-seeded): dispatching {len(jobs)} jobs in parallel ({system_keys})...",
        flush=True,
    )
    t0 = time.time()
    results = Parallel(n_jobs=-1, verbose=10)(jobs)
    print(f"anchor jobs done in {time.time() - t0:.1f}s")
    _emit_results(results)


def phase_grid(system_keys: list[str]) -> None:
    jobs = [
        delayed(_run_grid_job)(sys_key, k1, k2)
        for sys_key in system_keys
        for k1, k2 in GRID_TOPOLOGIES
    ]
    print(
        f"\nStep 2b (grid-seeded): dispatching {len(jobs)} jobs in parallel ({system_keys})...",
        flush=True,
    )
    t0 = time.time()
    results = Parallel(n_jobs=-1, verbose=10)(jobs)
    print(f"grid jobs done in {time.time() - t0:.1f}s")
    _emit_results(results)


def _emit_results(results: list[JobResult]) -> None:
    lines = ["\nResults"]
    for sys_key, k1, k2, method_key, r, elapsed in results:
        lines.append(f"[{sys_key}] ({k1},{k2}) [{method_key}]  {_fmt(r)}  [{elapsed:.1f}s]")
    _print_and_append(lines)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--phase", choices=["init", "anchors", "grid"], required=True)
    ap.add_argument(
        "--systems",
        nargs="*",
        default=list(REAL_BINARY_SYSTEMS),
        help="system keys to sweep",
    )
    args = ap.parse_args()

    preflight_search(
        task_no=549,
        region_id=_REGION_ID,
        method=_METHOD,
        script_path=Path(__file__),
        n_points=len(REAL_BINARY_SYSTEMS) * len(ANCHOR_TOPOLOGIES)
        + len(REAL_BINARY_SYSTEMS) * len(GRID_TOPOLOGIES) * 8 * 6 * 3,
        override_reason=(
            "reuses #494/#504's already-validated binary-cycler harness verbatim "
            "(positive control re-finds the committed PC (3,2) row to 9 sig figs); "
            "the (system, topology) job count is small and bounded by construction, "
            "not an open-ended discovery grid needing a timing pilot."
        ),
    )

    if args.phase == "init":
        phase_init()
    elif args.phase == "anchors":
        phase_anchors(args.systems)
    elif args.phase == "grid":
        phase_grid(args.systems)

    print(f"\n[phase={args.phase}] appended to {OUT_PATH}")


if __name__ == "__main__":
    main()
