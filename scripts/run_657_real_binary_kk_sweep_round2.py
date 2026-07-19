"""Task #657: real-binary (k1,k2)-cycler genome, round 2 -- runnable script.

`#654` shortlist item 3. Two independent pieces, reusing #549's/#504's
already-validated CR3BP binary-cycler harness
(:mod:`cyclerfinder.search.real_binary_kk_sweep`) VERBATIM -- no changes to
`real_binary_kk_sweep.py`'s corrector/sweep machinery, only the new sourced
`REAL_BINARY_SYSTEMS` entries (sila-nunam/antiope/lempo-hiisi) and this
thin driver.

(a) UNCAP #549's 4 compute-capped INCONCLUSIVE probes. #549's own ad-hoc
    finer-resolution robustness pass (10 of the 32 original anchor-seeded
    negatives re-run at <=0.001-mu-step resolution) hit an unprincipled
    240-second-per-job wall-clock cap for 4 of them -- Patroclus-Menoetius
    (1,1) [all 3 Table-I anchors], (3,1), (3,2), and Eris-Dysnomia (3,2) --
    without converging OR failing outright. Per
    `[[feedback_long_runs_acceptable]]`, re-run these EXACT probes (same
    `sweep_family` call, same fine mu-step resolution) with NO artificial
    per-job cap, letting `mu_step_to_system`'s own Newton corrector run to
    either convergence or a genuine failure (ValueError / non-convergence
    returned as `None`) at every step.

(b) Sweep 3 new near-equal-mu real binary systems #549 never tried --
    Sila-Nunam, Antiope, Lempo-Hiisi (as an isolated binary core, Paha's
    perturbation explicitly unmodeled) -- across the SAME six (k1,k2)
    topologies #549 swept ((1,1) from all 3 Table-I anchors, (3,1), (3,2),
    (3,3) anchor-seeded; (2,1),(2,2) grid-seeded), using the identical
    `ANCHOR_TOPOLOGIES`/`GRID_TOPOLOGIES` lists #549's own
    `scripts/run_549_real_binary_kk_sweep.py` defined.

Mandatory positive control (both pieces gated on this): re-find Pluto-
Charon (3,2) via #504's own `sweep_32_positive_control()`, UNMODIFIED,
before trusting anything below.

Usage
-----
  uv run python scripts/run_657_real_binary_kk_sweep_round2.py --phase positive_control
  uv run python scripts/run_657_real_binary_kk_sweep_round2.py --phase uncapped_recheck
  uv run python scripts/run_657_real_binary_kk_sweep_round2.py --phase new_systems_anchors
  uv run python scripts/run_657_real_binary_kk_sweep_round2.py --phase new_systems_grid
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
    ANCHORS,
    REAL_BINARY_SYSTEMS,
    SweepResult,
    sweep_family,
    sweep_family_grid,
)

OUT_PATH = (
    Path(__file__).parent.parent / "docs" / "notes" / "scratch" / "657_kk_sweep_round2_raw.txt"
)

_REGION_ID = "real-binary-kk-cycler-sweep-round2-2026-07-19"
_METHOD = MethodCapability(
    genome=(
        "Real-binary (k1,k2) CR3BP cycler genome (#494/#504/#549's fixed-Jacobi symmetric "
        "corrector + Barden stability + winding-topology classifier + independent Radau "
        "crosscheck), reused verbatim: (a) uncapped fine-resolution mu-continuation re-check "
        "of #549's 4 compute-capped inconclusive probes; (b) swept at three NEW sourced "
        "near-equal-mu real-binary mass ratios (Sila-Nunam, Antiope, Lempo-Hiisi-as-binary-core) "
        "across the same 6 (k1,k2) topologies #549 used"
    ),
    corrector=(
        "anchor-seeded mu-continuation from Ross-RT 2026 Table-I anchors for "
        "(1,1)/(3,1)/(3,2)/(3,3), fine (<=0.001 mu-step) resolution for the uncapped recheck; "
        "bounded (x0,C,hc) grid search for (2,1)/(2,2)"
    ),
    capability_tags=frozenset(
        {"cr3bp", "binary-cycler", "k1k2-genome", "real-binary", "mu-continuation"}
    ),
    git_sha="working-tree",
)

# Identical (k1,k2,anchor_key) topology set #549's own script defined --
# reused verbatim for the 3 new systems (task scope: same topologies, new
# mass ratios).
ANCHOR_TOPOLOGIES: list[tuple[int, int, str]] = [
    (1, 1, "mu001_11"),
    (1, 1, "mu01215_11"),
    (1, 1, "mu05_11"),
    (3, 2, "mu01_32"),
    (3, 1, "mu03_31"),
    (3, 3, "mu01215_33"),
]
GRID_TOPOLOGIES: list[tuple[int, int]] = [(2, 1), (2, 2)]

NEW_SYSTEM_KEYS = ["sila-nunam", "antiope", "lempo-hiisi"]

# The 4 probes #549's own ad-hoc fine-resolution robustness pass left
# INCONCLUSIVE at a 240s-per-job cap (see data/OUTSTANDING.md's #549 entry).
# Each tuple is (system_key, anchor_key) -- the anchor_key encodes (k1,k2)
# via ANCHORS[anchor_key].
UNCAPPED_RECHECK_PROBES: list[tuple[str, str]] = [
    ("patroclus-menoetius", "mu001_11"),  # (1,1)
    ("patroclus-menoetius", "mu01215_11"),  # (1,1)
    ("patroclus-menoetius", "mu05_11"),  # (1,1)
    ("patroclus-menoetius", "mu03_31"),  # (3,1)
    ("patroclus-menoetius", "mu01_32"),  # (3,2)
    ("eris-dysnomia", "mu01_32"),  # (3,2)
]


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


def _fine_n_steps(anchor_mu: float, target_mu: float, max_step: float = 0.001) -> int:
    """Steps needed to keep every mu-continuation increment <= max_step.

    Matches #549's own ad-hoc "each step <=0.001 in mu" fine-resolution
    robustness-pass convention (its own bullet: "up to 426 steps").
    """
    return max(1, int(np.ceil(abs(target_mu - anchor_mu) / max_step)))


def _append(lines: list[str]) -> None:
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    with OUT_PATH.open("a") as f:
        f.write("\n".join(lines) + "\n")


def _print_and_append(lines: list[str]) -> None:
    for line in lines:
        print(line, flush=True)
    _append(lines)


def phase_positive_control() -> None:
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    stamp = time.strftime("%Y-%m-%dT%H:%M:%SZ")
    OUT_PATH.write_text(f"Task #657 real-binary (k1,k2) sweep round 2 -- run {stamp}\n")
    print("Step 1: positive control (PC 3,2)...", flush=True)
    pc_line = _run_positive_control()
    _print_and_append([pc_line])

    print("\nSourced NEW systems (task #657):")
    for key in NEW_SYSTEM_KEYS:
        s = REAL_BINARY_SYSTEMS[key]
        lines = [
            f"SYSTEM {key}: mu={s.mu:.10f} l_km={s.l_km} t_s={s.t_s:.4f}",
            f"  mu_source: {s.mu_source}",
            f"  l_source: {s.l_source}",
            f"  t_source: {s.t_source}",
        ]
        if s.caveat:
            lines.append(f"  CAVEAT: {s.caveat}")
        _print_and_append(lines)


def _run_uncapped_job(sys_key: str, anchor_key: str) -> tuple[str, str, SweepResult, int, float]:
    target = REAL_BINARY_SYSTEMS[sys_key].to_cr3bp_system()
    anchor_mu = ANCHORS[anchor_key]["anchor_mu"]
    n_steps = _fine_n_steps(anchor_mu, target.mu)
    t0 = time.time()
    r = sweep_family(target, anchor_key, n_steps=n_steps)
    elapsed = time.time() - t0
    return (sys_key, anchor_key, r, n_steps, elapsed)


def phase_uncapped_recheck() -> None:
    """(a): re-run #549's 4 compute-capped probes, uncapped, fine resolution."""
    print(
        f"\nUncapped fine-resolution recheck: {len(UNCAPPED_RECHECK_PROBES)} probes "
        "(NO artificial per-job cap)...",
        flush=True,
    )
    jobs = [
        delayed(_run_uncapped_job)(sys_key, anchor_key)
        for sys_key, anchor_key in UNCAPPED_RECHECK_PROBES
    ]
    t0 = time.time()
    results = Parallel(n_jobs=-1, verbose=10)(jobs)
    total = time.time() - t0
    lines = [f"\nUncapped recheck results (total wall {total:.1f}s):"]
    for sys_key, anchor_key, r, n_steps, elapsed in results:
        a = ANCHORS[anchor_key]
        lines.append(
            f"[{sys_key}] ({a['k1']},{a['k2']}) [{anchor_key}, n_steps={n_steps}]  "
            f"{_fmt(r)}  [{elapsed:.1f}s]"
        )
    _print_and_append(lines)


def _run_anchor_job(
    sys_key: str, k1: int, k2: int, anchor_key: str
) -> tuple[str, int, int, str, SweepResult, float]:
    target = REAL_BINARY_SYSTEMS[sys_key].to_cr3bp_system()
    t0 = time.time()
    r = sweep_family(target, anchor_key)
    return (sys_key, k1, k2, anchor_key, r, time.time() - t0)


def _run_grid_job(sys_key: str, k1: int, k2: int) -> tuple[str, int, int, str, SweepResult, float]:
    target = REAL_BINARY_SYSTEMS[sys_key].to_cr3bp_system()
    mu = target.mu
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


def _emit_results(results: list[tuple[str, int, int, str, SweepResult, float]]) -> None:
    lines = ["\nResults"]
    for sys_key, k1, k2, method_key, r, elapsed in results:
        lines.append(f"[{sys_key}] ({k1},{k2}) [{method_key}]  {_fmt(r)}  [{elapsed:.1f}s]")
    _print_and_append(lines)


def phase_new_systems_anchors() -> None:
    jobs = [
        delayed(_run_anchor_job)(sys_key, k1, k2, anchor_key)
        for sys_key in NEW_SYSTEM_KEYS
        for k1, k2, anchor_key in ANCHOR_TOPOLOGIES
    ]
    print(
        f"\nNew-systems (anchor-seeded): dispatching {len(jobs)} jobs in parallel "
        f"({NEW_SYSTEM_KEYS})...",
        flush=True,
    )
    t0 = time.time()
    results = Parallel(n_jobs=-1, verbose=10)(jobs)
    print(f"anchor jobs done in {time.time() - t0:.1f}s")
    _emit_results(results)


def phase_new_systems_grid() -> None:
    jobs = [
        delayed(_run_grid_job)(sys_key, k1, k2)
        for sys_key in NEW_SYSTEM_KEYS
        for k1, k2 in GRID_TOPOLOGIES
    ]
    print(
        f"\nNew-systems (grid-seeded): dispatching {len(jobs)} jobs in parallel "
        f"({NEW_SYSTEM_KEYS})...",
        flush=True,
    )
    t0 = time.time()
    results = Parallel(n_jobs=-1, verbose=10)(jobs)
    print(f"grid jobs done in {time.time() - t0:.1f}s")
    _emit_results(results)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--phase",
        choices=["positive_control", "uncapped_recheck", "new_systems_anchors", "new_systems_grid"],
        required=True,
    )
    args = ap.parse_args()

    n_points = (
        1  # positive control
        + len(UNCAPPED_RECHECK_PROBES)
        + len(NEW_SYSTEM_KEYS) * len(ANCHOR_TOPOLOGIES)
        + len(NEW_SYSTEM_KEYS) * len(GRID_TOPOLOGIES) * 8 * 6 * 3
    )
    preflight_search(
        task_no=657,
        region_id=_REGION_ID,
        method=_METHOD,
        script_path=Path(__file__),
        n_points=n_points,
        override_reason=(
            "reuses #494/#504/#549's already-validated binary-cycler harness verbatim "
            "(positive control re-finds the committed PC (3,2) row); the uncapped recheck "
            "targets 6 SPECIFIC previously-run probes (no new region), and the 3 new systems' "
            "(system, topology) job count is small and bounded by construction, not an "
            "open-ended discovery grid needing a timing pilot."
        ),
    )

    if args.phase == "positive_control":
        phase_positive_control()
    elif args.phase == "uncapped_recheck":
        phase_uncapped_recheck()
    elif args.phase == "new_systems_anchors":
        phase_new_systems_anchors()
    elif args.phase == "new_systems_grid":
        phase_new_systems_grid()

    print(f"\n[phase={args.phase}] appended to {OUT_PATH}")


if __name__ == "__main__":
    main()
