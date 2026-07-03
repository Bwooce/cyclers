"""#532 fine-grid follow-up: the coarse pilot (`run_532_resonance_overlap_pilot.py`,
0.03 C-step) found a NEAR-MISS, not a clean non-overlap -- the 3:2 Hilda
family's unstable orbits sampled up to C=3.1413, and the 2:1 interior
family's one sampled unstable orbit sits at C=3.1453, just 0.004 beyond the
coarse grid's own edge. That gap is smaller than the coarse step itself, so
the coarse "DO NOT OVERLAP" verdict may simply be a sampling-resolution
artifact, not a real non-overlap.

This script re-scans BOTH families at a 5x finer C-step (0.005 instead of
0.03) restricted to the C=[3.10, 3.18] window bracketing the near-miss,
reusing the exact same seed derivation, enumeration, correction, and
monodromy-stability machinery as the coarse pilot (imported directly, not
reimplemented) -- only the C-band resolution/window changes.
"""

from __future__ import annotations

import dataclasses
import datetime
import pathlib
import sys
import time

_REPO = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_REPO / "src"))
sys.path.insert(0, str(_REPO))

import cyclerfinder.core.cr3bp as cr3bp  # noqa: E402
from cyclerfinder.data.preflight import PreflightBlockedError, preflight_search  # noqa: E402
from scripts.run_532_resonance_overlap_pilot import (  # noqa: E402
    _METHOD,
    GRID,
    HILDA,
    INTERIOR_2_1,
    N_RANGE,
    _positive_control,
    _scan_family,
    _seed_from_period_ratio,
)

# Fine C-band: 3.10 to 3.18 in steps of 0.005 (17 points), bracketing the
# coarse pilot's near-miss (Hilda unstable up to 3.1413, interior unstable at
# 3.1453) with margin on both sides.
FINE_C_WINDOW = tuple(round(3.10 + 0.005 * i, 4) for i in range(17))


def _ts() -> str:
    return datetime.datetime.now().strftime("%Y-%m-%dT%H:%M:%S")


def _fine_spec(base, mu: float):
    _, _, c_seed = _seed_from_period_ratio(mu, base.period_ratio)
    offsets = tuple(round(c - c_seed, 6) for c in FINE_C_WINDOW)
    return dataclasses.replace(base, c_offsets=offsets)


def main() -> None:
    print(f"[{_ts()}] #532 fine-grid overlap follow-up starting (window {FINE_C_WINDOW}).")
    system = cr3bp.cr3bp_system("Sun", "Jupiter")
    mu = float(system.mu)

    hilda_fine = _fine_spec(HILDA, mu)
    interior_fine = _fine_spec(INTERIOR_2_1, mu)

    n_points = 2 * len(FINE_C_WINDOW) * len(N_RANGE) * GRID[0] * GRID[1]
    preflight_search(
        task_no=532,
        region_id="sun-jupiter-21-32-mmr-unstable-c-overlap-finegrid",
        method=_METHOD,
        script_path=pathlib.Path(__file__),
        n_points=n_points,
        # Coarse pilot measured ~451.6s / 10080 points = 0.0448 s/point on the
        # identical enumerator/grid/correction pipeline; this run is the same
        # machinery, just a finer/narrower C-band.
        timing_pilot_seconds_per_point=0.05,
    )

    _positive_control(system)

    t0 = time.time()
    hilda_results = _scan_family(system, hilda_fine)
    interior_results = _scan_family(system, interior_fine)
    dt = time.time() - t0

    hilda_unstable_c = sorted(
        c for c, orbits in hilda_results.items() if any(r["unstable"] for r in orbits)
    )
    interior_unstable_c = sorted(
        c for c, orbits in interior_results.items() if any(r["unstable"] for r in orbits)
    )

    print()
    print(f"[{_ts()}] Fine-grid scan complete in {dt:.1f}s.")
    print(f"[{_ts()}] Hilda (3:2) unstable C values in window: {hilda_unstable_c}")
    print(f"[{_ts()}] Interior (2:1) unstable C values in window: {interior_unstable_c}")

    shared = sorted(
        set(round(c, 3) for c in hilda_unstable_c) & set(round(c, 3) for c in interior_unstable_c)
    )
    print(f"[{_ts()}] Exact-C matches (rounded to 0.001): {shared}")
    if shared:
        print(
            f"[{_ts()}] GO: at least one Jacobi constant hosts an unstable orbit in BOTH "
            f"families -- #532's full cross-family connection search is now worth building."
        )
    elif hilda_unstable_c and interior_unstable_c:
        lo = max(min(hilda_unstable_c), min(interior_unstable_c))
        hi = min(max(hilda_unstable_c), max(interior_unstable_c))
        if lo <= hi:
            print(
                f"[{_ts()}] CLOSE: both families show unstable orbits within the same "
                f"[{lo:.4f}, {hi:.4f}] sub-band at this resolution, but no exact-C match -- "
                f"a still-finer scan or a direct continuation bridging the two nearest "
                f"unstable orbits is the next step, not yet a confirmed shared-Jacobi pair."
            )
        else:
            print(
                f"[{_ts()}] NO-GO at this finer resolution: unstable sub-ranges "
                f"[{min(hilda_unstable_c):.4f},{max(hilda_unstable_c):.4f}] (Hilda) and "
                f"[{min(interior_unstable_c):.4f},{max(interior_unstable_c):.4f}] (interior) "
                f"still do not overlap."
            )
    else:
        print(
            f"[{_ts()}] NO-GO: at least one family shows no unstable member in this "
            f"finer window at all."
        )


if __name__ == "__main__":
    try:
        main()
    except PreflightBlockedError as exc:
        print(f"[{_ts()}] BLOCKED by preflight_search:\n{exc}")
        sys.exit(1)
