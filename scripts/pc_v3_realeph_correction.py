"""#511: Pluto-Charon (3,2) real-ephemeris differential-correction lever.

The PROPER V2->V3 real-eph lever for ``ross-rt-pc-cycler-32-2026``, as scoped
by #506 (``docs/notes/2026-07-01-506-pluto-charon-v3-scope.md``): NOT a naive
propagation test (rejected -- it only measures a ~756 m model-mismatch
oscillation from Charon's tiny eccentricity, not stability), but a
DIFFERENTIAL CORRECTION of the CR3BP (3,2) periodic orbit into a real-eph
analog.

Pipeline:
  1. Resolve ``plu060.bsp`` (#510) and read Charon's REAL osculating orbital
     elements relative to Pluto at a chosen epoch via SPICE
     (:func:`cyclerfinder.verify.pluto_charon_realeph.charon_osculating_elements`).
  2. Re-converge the CR3BP (3,2) seed exactly at e=0 in the ER3BP pulsating
     frame (exact bridge, #441 Sec. 1).
  3. Continue the family from e=0 to the REAL eccentricity via the existing,
     independently-validated (#293) ER3BP e-continuator
     (:func:`cyclerfinder.genome.er3bp_continuation.continue_er3bp_family_in_e`).
  4. Gate on BOTH the corrector residual (symmetric half-period condition
     only) and the INDEPENDENT full-orbit-closure residual (the #441
     period_f-trap gate) -- a converged corrector residual with a large
     independent residual is a false-positive-closure signature, not a
     real-eph analog orbit.

Usage:
    uv run python scripts/pc_v3_realeph_correction.py
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import cyclerfinder.core.cr3bp as cr3bp
import cyclerfinder.search.cr3bp_periodic as cp
from cyclerfinder.verify.pluto_charon_realeph import (
    charon_osculating_elements,
    differential_correct_pc32_to_eccentricity,
)
from cyclerfinder.verify.spice_kernels import ensure_pluto_kernel

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from scripts.pc_v2_longspan import (
    _HALF_CROSSINGS,
    _YDOT0_SIGN,
    PC_C,
    PC_T_GUESS,
    PC_X0,
    ROW_ID,
)

METHOD_TAG = "pc-v3-realeph-correction-v1"

# Representative epoch for the SPICE osculating-element read. Charon's
# osculating eccentricity is essentially epoch-invariant (2.17e-4 - 2.33e-4
# sampled 2000-2050, see the #511 verdict note), so any epoch in the kernel's
# 1800-2199 coverage is representative; "today" per project convention.
DEFAULT_EPOCH_ISO = "2026-07-01T00:00:00"

# Independent-closure gate (#441): production tolerance for "actually closes".
INDEPENDENT_TOL = 1e-8


def run_pc_v3_realeph(epoch_iso: str = DEFAULT_EPOCH_ISO, n_steps: int = 20):
    t0 = time.monotonic()
    kernel_path = ensure_pluto_kernel()

    charon = charon_osculating_elements(epoch_iso, kernel_path)

    system = cr3bp.cr3bp_system("Pluto", "Charon")
    orbit = cp.correct_symmetric_fixed_jacobi(
        system,
        PC_X0,
        PC_C,
        PC_T_GUESS,
        ydot0_sign=_YDOT0_SIGN,
        half_crossings=_HALF_CROSSINGS,
        tol=1e-10,
    )
    if not orbit.converged:
        raise RuntimeError(
            f"PC (3,2) CR3BP seed did not converge (res={orbit.crossing_residual:.2e})"
        )

    result = differential_correct_pc32_to_eccentricity(
        charon.eccentricity,
        x0_seed=orbit.x0,
        ydot0_seed=orbit.ydot0,
        period_seed=orbit.period,
        mu=system.mu,
        n_steps=n_steps,
        independent_tol=INDEPENDENT_TOL,
    )

    elapsed = time.monotonic() - t0
    return charon, orbit, result, elapsed


def main() -> None:
    print(f"# {METHOD_TAG}: PC (3,2) real-ephemeris differential correction (#511)")
    print(f"# row: {ROW_ID}")
    print(f"[{time.strftime('%H:%M:%S')}] resolving kernel + running ...", flush=True)

    charon, orbit, result, elapsed = run_pc_v3_realeph()

    print(f"  -> Charon real osculating orbit ({DEFAULT_EPOCH_ISO}):")
    print(
        f"     a={charon.a_km:.5f} km, e={charon.eccentricity:.6e}, "
        f"T={charon.period_days:.6f} d, r={charon.r_km:.3f} km"
    )
    print(
        f"  -> CR3BP seed: T_nd={orbit.period:.10f}, "
        f"T/(2*pi)={result.period_ratio_sc_to_charon:.6f}"
    )
    print(
        "     (non-integer => structurally incommensurate with Charon's own "
        "orbital period under the ER3BP pulsating frame -- #441 period_f trap)"
    )
    print(
        f"  -> e=0 bridge: corrector_res={result.seed_corrector_residual:.3e}, "
        f"independent_res={result.seed_independent_residual:.3e}"
    )
    print(
        f"  -> e={charon.eccentricity:.3e} target: corrector_res="
        f"{result.target_corrector_residual:.3e}, independent_res="
        f"{result.target_independent_residual:.3e} (gate {INDEPENDENT_TOL:.1e})"
    )
    print(f"  -> CONVERGED (real-eph analog found): {result.converged}")
    print(f"  -> elapsed: {elapsed:.1f}s")


if __name__ == "__main__":
    main()
