"""#480 M1 — second STM corrector chunk, reseeded from chunk-1 corrected states.

Rebuilds the periapsis seed, runs chunk-1 (max_nfev=60) to get corrected states,
saves them, builds a new ShootingSeed from the corrected states (same epochs,
sequence, V∞ vectors carried from the original periapsis seed), then runs a
second chunk (max_nfev=60). Reports final defect + maintenance ΔV.

Usage::

    timeout 540 uv run python scripts/_v2_stm_chunk2_480.py [runlog_path]
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import spiceypy

from cyclerfinder.nbody.jovian import (
    JovianEphemeris,
    JovianRailsCache,
    jovian_defect_residual,
    jovian_shoot,
)
from cyclerfinder.nbody.shooter import ShootingSeed
from cyclerfinder.search.ieg_seed import paper_departure_et
from cyclerfinder.verify.spice_kernels import ensure_jup365_kernel, ensure_leapseconds_kernel

_SPD = 86400.0
_TABLE4_VINF_KMS = {"Europa": 9.12, "Ganymede": 7.07, "Io": 8.38}


def main() -> None:
    runlog_path = (
        Path(sys.argv[1])
        if len(sys.argv) > 1
        else Path(
            "/tmp/claude-1000/-home-bruce-dev-cyclers/"
            "b2ef8d69-0098-4763-909e-37ee4edbdf86/scratchpad/l3stm_chunk2.log"
        )
    )
    runlog_path.parent.mkdir(parents=True, exist_ok=True)
    log = open(runlog_path, "w", buffering=1)  # noqa: SIM115

    def write(msg: str) -> None:
        ts = time.strftime("%Y-%m-%dT%H:%M:%S")
        line = f"{ts} {msg}"
        print(line, flush=True)
        log.write(line + "\n")
        log.flush()

    write("START _v2_stm_chunk2_480.py")
    t_global = time.monotonic()

    kernel_path = ensure_jup365_kernel()
    spiceypy.furnsh(ensure_leapseconds_kernel())
    jeph = JovianEphemeris(kernel_path)
    moons = ("Io", "Europa", "Ganymede")

    # Import seed builder from the chunk-1 script (same module, avoids duplication).
    sys.path.insert(0, str(Path(__file__).parent))
    from _v2_stm_480 import build_realeph_periapsis_seed  # type: ignore[import]

    write("Building periapsis seed at paper_et + 3.78 d ...")
    paper_et = paper_departure_et()
    departure_et = paper_et + 3.78 * _SPD
    seed = build_realeph_periapsis_seed(departure_et, jeph)
    write(f"  n_nodes={len(seed.sequence)} sequence={seed.sequence}")

    # --- Chunk 1 ---
    write("=== CHUNK 1: jovian_shoot(stm, max_nfev=60) ===")
    t0 = time.monotonic()
    result1 = jovian_shoot(
        seed,
        kernel_path=kernel_path,
        moons=moons,
        accuracy=1e-11,
        max_nfev=60,
        max_wall_sec=500.0,
        jacobian="stm",
    )
    write(f"CHUNK 1: done in {time.monotonic() - t0:.1f}s, nfev={result1.n_iterations}")
    write(f"  converged={result1.converged}")
    write(f"  defect_norm={result1.defect_norm:.6e}")
    write(f"  correction_dv_ms={result1.correction_dv_kms * 1000:.2f}")

    # Save corrected states for chunk 2.
    cs1 = result1.corrected_states
    np.save(
        runlog_path.parent / "l3stm_chunk1_states.npy",
        np.array([np.asarray(s, dtype=np.float64) for s in cs1]),
    )
    write(f"Corrected states saved ({len(cs1)} nodes)")

    # --- Build reseeded ShootingSeed from chunk-1 corrected states ---
    # Carry epochs, sequence, tofs, vinf_in/out from original seed.
    # The corrected states are the new node_states (better starting point).
    seed2 = ShootingSeed(
        node_states=[np.asarray(s, dtype=np.float64) for s in cs1],
        epochs=list(seed.epochs),
        tofs=list(seed.tofs),
        sequence=seed.sequence,
        slack_leg=seed.slack_leg,
        period_days=seed.period_days,
        vinf_in=[np.asarray(v, dtype=np.float64) for v in seed.vinf_in],
        vinf_out=[np.asarray(v, dtype=np.float64) for v in seed.vinf_out],
    )

    # Verify seed2 defect == chunk-1 final defect.
    cache2 = JovianRailsCache(moons, jeph, min(seed2.epochs), max(seed2.epochs))
    res2_seed = jovian_defect_residual(
        seed2, ephem=jeph, cache=cache2, moons=moons, accuracy=1e-11, max_wall_sec=60.0
    )
    write(
        f"Chunk-2 seed defect_norm = {float(np.linalg.norm(res2_seed)):.6e} (expect ~chunk-1 final)"
    )

    # --- Chunk 2 ---
    write("=== CHUNK 2: jovian_shoot(stm, max_nfev=60) from chunk-1 states ===")
    t0 = time.monotonic()
    result2 = jovian_shoot(
        seed2,
        kernel_path=kernel_path,
        moons=moons,
        accuracy=1e-11,
        max_nfev=60,
        max_wall_sec=500.0,
        jacobian="stm",
    )
    write(f"CHUNK 2: done in {time.monotonic() - t0:.1f}s, nfev={result2.n_iterations}")
    write(f"  converged={result2.converged}")
    write(f"  seed_defect_norm={result2.seed_defect_norm:.6e}")
    write(f"  final_defect_norm={result2.defect_norm:.6e}")
    write(f"  correction_dv_ms={result2.correction_dv_kms * 1000:.2f}")
    write(f"  bend_feasible={result2.bend_feasible}")

    write("CHUNK 2 Corrected V∞ per encounter:")
    for i, (body, vinf_km) in enumerate(
        zip(result2.sequence, result2.vinf_per_encounter_kms, strict=True)
    ):
        tgt = _TABLE4_VINF_KMS.get(body, float("nan"))
        write(f"  node {i} {body}: V∞={vinf_km:.3f} km/s  (Table-4 {tgt:.2f})")

    # Save chunk-2 corrected states.
    np.save(
        runlog_path.parent / "l3stm_chunk2_states.npy",
        np.array([np.asarray(s, dtype=np.float64) for s in result2.corrected_states]),
    )
    write("Chunk-2 corrected states saved.")

    # --- Summary ---
    write("=== SUMMARY ===")
    write(f"chunk1_defect: {result1.defect_norm:.6e}")
    write(f"chunk2_defect: {result2.defect_norm:.6e}")
    write(f"chunk2_converged: {result2.converged}")
    write(f"chunk2_correction_dv_ms: {result2.correction_dv_kms * 1000:.2f}")
    write(f"TOTAL wall: {time.monotonic() - t_global:.1f}s")

    summary = {
        "chunk1": {
            "defect_norm": result1.defect_norm,
            "correction_dv_ms": result1.correction_dv_kms * 1000,
            "converged": result1.converged,
        },
        "chunk2": {
            "defect_norm": result2.defect_norm,
            "seed_defect_norm": result2.seed_defect_norm,
            "correction_dv_ms": result2.correction_dv_kms * 1000,
            "converged": result2.converged,
            "vinf_kms": list(result2.vinf_per_encounter_kms),
            "sequence": list(result2.sequence),
        },
        "table4_vinf_kms": _TABLE4_VINF_KMS,
    }
    summary_path = runlog_path.with_suffix(".json")
    summary_path.write_text(json.dumps(summary, indent=2))
    write(f"Summary written to {summary_path}")
    log.close()


if __name__ == "__main__":
    main()
