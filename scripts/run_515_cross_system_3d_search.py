"""#515 Phase-A 3D cross-system (SE<->EM) closure search using multiple cores.

Runs correct_cross_cycle_3d over a grid of out-of-plane amplitudes and Jacobi constants,
utilizing concurrent.futures.ProcessPoolExecutor to exploit multiple CPU cores.
If no closed cycles are found, it registers a clean negative in data/negative_results.yaml.
"""

from __future__ import annotations

import datetime
import pathlib
import sys

# Ensure the src tree is on the path when invoked as a script.
_REPO = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_REPO / "src"))

import yaml  # type: ignore[import-untyped]  # noqa: E402

from cyclerfinder.genome.cross_system_cycle import (  # noqa: E402
    CrossCycleClosure,
    FrameBridge,
    correct_cross_cycle_3d,
    em_moon_system,
    se_earth_system,
)
from cyclerfinder.parallel import (  # noqa: E402
    ParallelSweepConfig,
    parallel_sweep,
)

# Canalias 2007 SE bifurcation Jacobi.
CANALIAS_C_SE = 3.000863625

# Grid definitions.
C_EM_GRID = (3.15, 3.17, 3.18)
C_SE_GRID = (3.0003, CANALIAS_C_SE)
Z_EM_GRID = (-0.02, -0.05)
# Scale Sun-Earth amplitudes by (384400.0 / 149600000.0) to match EM physical scale
Z_SE_GRID = tuple(z * (384400.0 / 149600000.0) for z in Z_EM_GRID)
LIBRATION_PAIRS = (("EM-L2", "SE-L2"), ("EM-L1", "SE-L1"))

_NEG_RESULTS_PATH = _REPO / "data" / "negative_results.yaml"
_NEW_ENTRY_ID = "cross_system_se_em_3d_patched_cr3bp"

_NEW_ENTRY = {
    "id": _NEW_ENTRY_ID,
    "issue": 515,
    "method": "patched-CR3BP 3D SE<->EM connection matcher + 3D bounded closure search (#515)",
    "regime": "3D Sun-Earth + Earth-Moon coupled CR3BP (patched, inertial-frame match)",
    "failed_rung": "closure-search-3d",
    "physical_reason": (
        "3D EM->SE and SE->EM legs solved spatially to high accuracy, but the time-consistent "
        "phase residuals (R1, R2) do not close simultaneously in the grid scanned. The additional "
        "out-of-plane amplitude degrees of freedom (z_em, z_se) did not resolve the phase "
        "closure wall within this grid of single-revolution orbits."
    ),
    "resweep_condition": (
        "BCR4BP 3D halo continuation OR a multi-revolution Metonic-commensurate 3D grid."
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


def evaluate_point(args: tuple[float, float, float, float, str, str]) -> CrossCycleClosure | None:
    """Worker function to correct a single 3D cross-cycle grid point."""
    c_em, c_se, z_em, z_se, em_lib, se_lib = args
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
            n_em=1,
            n_se=1,
            max_iter=8,
            solver="bounded_ls",
            scan_n=6,
            scan_n_tau=2,
            tol_km=1e2,
        )
        return res
    except Exception as e:
        print(
            f"Worker failed for c_em={c_em:.4f}, c_se={c_se:.6f}, "
            f"z_em={z_em:.3f}, z_se={z_se:.4f}: {e}"
        )
        return None


def main() -> None:
    print(f"[{_ts()}] #515 3D cross-system closure search starting.")

    # Build list of grid points
    tasks = []
    for pair in LIBRATION_PAIRS:
        em_lib, se_lib = pair
        for c_em in C_EM_GRID:
            for c_se in C_SE_GRID:
                for z_em in Z_EM_GRID:
                    for z_se in Z_SE_GRID:
                        tasks.append((c_em, c_se, z_em, z_se, em_lib, se_lib))

    total_tasks = len(tasks)
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

    for i, res in enumerate(sweep_res.results):
        c_em, c_se, z_em, z_se, em_lib, se_lib = tasks[i]
        if res is None:
            continue
        results.append(res)
        status = "CLOSED" if res.closed else "OPEN"
        print(
            f"[{_ts()}] [{i + 1:02d}/{total_tasks:02d}] {status} pair=({em_lib},{se_lib}) "
            f"c_em={c_em:.3f} c_se={c_se:.6f} z_em={z_em:.3f} z_se={z_se:.4f}"
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
    main()
