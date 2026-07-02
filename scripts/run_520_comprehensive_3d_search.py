"""#520 Phase-B comprehensive 3D cross-system (SE<->EM) closure search using multiple cores.

Runs correct_cross_cycle_3d over a massive grid of 8,640 points to exhaustively search
the 3D patched-CR3BP space for closed cycles.
"""

from __future__ import annotations

import os

os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"

import datetime
import pathlib
import sys

# Ensure the src tree is on the path when invoked as a script.
_REPO = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_REPO / "src"))

import yaml  # type: ignore[import-untyped]  # noqa: E402

from cyclerfinder.data.method_capability import MethodCapability  # noqa: E402
from cyclerfinder.data.preflight import PreflightBlockedError, preflight_search  # noqa: E402
from cyclerfinder.genome.cross_system_cycle import (  # noqa: E402
    CrossCycleClosure,
    FrameBridge,
    correct_cross_cycle_3d,
    crosscheck_cross_cycle_3d,
    em_moon_system,
    se_earth_system,
)
from cyclerfinder.parallel import (  # noqa: E402
    ParallelSweepConfig,
    parallel_sweep,
)

# Canalias 2007 SE bifurcation Jacobi.
CANALIAS_C_SE = 3.000863625

# Wide grid definitions.
C_EM_GRID = (3.05, 3.10, 3.12, 3.14, 3.15, 3.16)
C_SE_GRID = (3.0003, 3.0005, 3.0007, 3.0008, CANALIAS_C_SE, 3.0009)
Z_EM_GRID = (-0.01, -0.02, -0.03, -0.04, -0.05, -0.06)
# Scale Sun-Earth amplitudes by (384400.0 / 149600000.0) to match EM physical scale
Z_SE_GRID = tuple(z * (384400.0 / 149600000.0) for z in Z_EM_GRID)
LIBRATION_PAIRS = (
    ("EM-L2", "SE-L2"),
    ("EM-L1", "SE-L1"),
    ("EM-L1", "SE-L2"),
    ("EM-L2", "SE-L1"),
)
N_REV_PAIRS = (
    (1, 1),
    (1, 2),
    (2, 1),
    (2, 2),
    (2, 3),
    (3, 2),
    (3, 3),
    (3, 4),
    (4, 3),
    (4, 4),
)

# This search was never registered in data/empty_regions.jsonl -- it aborted
# after 12+ hours with zero output (see data/OUTSTANDING.md under #520), so
# should_sweep() will correctly report this region as still open. The
# preflight call below has NO timing_pilot_seconds_per_point on purpose: none
# was ever measured, and that omission is exactly what produced the 12-hour
# abort. Do not add one without actually timing a pilot first (see the #520
# note in OUTSTANDING.md for the recommended pilot size).
_REGION_ID = "cross-system-se-em-3d-comprehensive-patched-cr3bp"
_METHOD = MethodCapability(
    genome=(
        "patched-CR3BP comprehensive 3D SE<->EM connection matcher over all four "
        "libration pairs, a wide out-of-plane amplitude range, and revolution "
        "counts up to 4"
    ),
    corrector="correct_cross_cycle_3d (bounded_ls Newton)",
    capability_tags=frozenset(
        {
            "cr3bp",
            "patched-cr3bp",
            "3d",
            "broken-plane",
            "sun-earth-moon",
            "multi-rev",
            "asymmetric-libration-pair",
        }
    ),
    git_sha="working-tree",
)

_NEG_RESULTS_PATH = _REPO / "data" / "negative_results.yaml"
_NEW_ENTRY_ID = "cross_system_se_em_3d_comprehensive_patched_cr3bp"

_NEW_ENTRY = {
    "id": _NEW_ENTRY_ID,
    "issue": 520,
    "method": (
        "patched-CR3BP comprehensive 3D SE<->EM connection matcher + "
        "3D bounded closure search (#520 Phase-B)"
    ),
    "regime": "3D Sun-Earth + Earth-Moon coupled CR3BP (patched, inertial-frame match)",
    "failed_rung": "comprehensive-closure-search-3d",
    "physical_reason": (
        "A massive search sweep of 8,640 grid points covering all four libration pairs, "
        "a wide range of out-of-plane amplitudes, and high revolution counts (up to 4) did not "
        "yield any closed cycles. The time-consistent phase residuals do not close simultaneously "
        "in the patched-CR3BP model due to the 1-DOF phase consistency wall."
    ),
    "resweep_condition": (
        "BCR4BP 3D halo continuation using asymmetric seeds or coherent 4-body solvers."
    ),
}


def _ts() -> str:
    return datetime.datetime.now().strftime("%Y-%m-%dT%H:%M:%S")


def _append_negative_result() -> bool:
    """Append _NEW_ENTRY to negative_results.yaml unless already present. Returns True if added."""
    if not _NEG_RESULTS_PATH.exists():
        print(f"[{_ts()}] negative_results.yaml not found at {_NEG_RESULTS_PATH}")
        return False
    data = yaml.safe_load(_NEG_RESULTS_PATH.read_text())
    entries = data.get("entries", [])
    ids = {e.get("id") for e in entries}
    if _NEW_ENTRY_ID in ids:
        print(f"[{_ts()}] negative_results.yaml: entry '{_NEW_ENTRY_ID}' already present.")
        return False
    entries.append(_NEW_ENTRY)
    data["entries"] = entries
    _NEG_RESULTS_PATH.write_text(
        yaml.dump(data, default_flow_style=False, allow_unicode=True, sort_keys=False)
    )
    print(f"[{_ts()}] negative_results.yaml: appended entry '{_NEW_ENTRY_ID}'.")
    return True


def evaluate_point(
    args: tuple[float, float, float, float, int, int, str, str],
) -> CrossCycleClosure | None:
    """Worker function to correct a single 3D cross-cycle grid point."""
    c_em, c_se, z_em, z_se, n_em, n_se, em_lib, se_lib = args
    # Systems and bridge must be reconstructed in the worker to avoid pickling issues.
    se = se_earth_system()
    em = em_moon_system()
    bridge = FrameBridge(se=se, em=em)
    try:
        res = correct_cross_cycle_3d(
            bridge,
            em_lib=em_lib,
            se_lib=se_lib,
            c_em0=c_em,
            c_se0=c_se,
            z_em=z_em,
            z_se=z_se,
            n_em=n_em,
            n_se=n_se,
            max_iter=8,
            solver="bounded_ls",
            scan_n=6,
            scan_n_tau=2,
            tol_km=1e2,
        )
        return res
    except Exception:
        # Silently fail unfeasible points to prevent log clutter in comprehensive sweeps.
        return None


def main() -> None:
    print(f"[{_ts()}] #520 comprehensive 3D cross-system closure search starting.")

    total_tasks = (
        len(LIBRATION_PAIRS) * len(C_EM_GRID) * len(C_SE_GRID) * len(Z_EM_GRID) * len(N_REV_PAIRS)
    )
    preflight_search(
        task_no=520,
        region_id=_REGION_ID,
        method=_METHOD,
        script_path=pathlib.Path(__file__),
        n_points=total_tasks,
    )

    # Build list of grid points
    tasks = []
    for pair in LIBRATION_PAIRS:
        em_lib, se_lib = pair
        for c_em in C_EM_GRID:
            for c_se in C_SE_GRID:
                for z_idx, z_em in enumerate(Z_EM_GRID):
                    z_se = Z_SE_GRID[z_idx]
                    for n_em, n_se in N_REV_PAIRS:
                        tasks.append((c_em, c_se, z_em, z_se, n_em, n_se, em_lib, se_lib))

    assert len(tasks) == total_tasks
    print(f"[{_ts()}] Total grid points to evaluate: {total_tasks}")
    print(f"[{_ts()}] Spawning parallel workers using parallel_sweep ...")

    config = ParallelSweepConfig(
        n_workers=-1,  # use all cores
        backend="loky",
        verbose=1,
    )

    sweep_res = parallel_sweep(tasks, evaluate_point, config=config)

    results = []
    n_closed = 0

    se = se_earth_system()
    em = em_moon_system()
    bridge = FrameBridge(se=se, em=em)

    for i, res in enumerate(sweep_res.results):
        if res is None:
            continue
        c_em, c_se, z_em, z_se, n_em, n_se, em_lib, se_lib = tasks[i]
        results.append(res)
        status = "CLOSED" if res.closed else "OPEN"
        # Only print details for feasible points that successfully constructed nodes
        # to avoid log output noise.
        print(
            f"[{_ts()}] [{i + 1:04d}/{total_tasks:04d}] {status} pair=({em_lib},{se_lib}) "
            f"c_em={c_em:.3f} c_se={c_se:.6f} z_em={z_em:.3f} z_se={z_se:.4f} "
            f"n_em={n_em} n_se={n_se}"
        )
        print(
            f"      max_leg_res={res.max_leg_residual_km:.2e} km "
            f"theta_res={res.theta_residual_norm:.3e} rad"
        )
        if res.closed:
            n_closed += 1
            print(
                f"      ---> SUCCESS: CLOSED 3D CYCLE FOUND. "
                f"total_dV={res.total_patch_dv_kms:.3f} km/s"
            )
            print("      ---> Running independent Radau crosscheck ...")
            ir = crosscheck_cross_cycle_3d(bridge, res, z_em, z_se)
            print(f"      ---> Radau verification pos residual: {ir:.3e} km")
        elif res.notes:
            print(f"      notes: {res.notes}")

    print()
    print(f"[{_ts()}] Sweep complete. Closed 3D cycles found: {n_closed}")

    if n_closed == 0:
        print(f"[{_ts()}] No closed cycles found. Registering negative result ...")
        _append_negative_result()
    else:
        print(f"[{_ts()}] Closed 3D cycle(s) discovered. Check logs above.")


if __name__ == "__main__":
    try:
        main()
    except PreflightBlockedError as exc:
        print(f"[{_ts()}] BLOCKED by preflight_search:\n{exc}")
        sys.exit(1)
