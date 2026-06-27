"""#480 M1 SCIENCE VERDICT — IEG EGGIE fine epoch scan (multi-rev seed).

Both infra fixes are now in (Task 4a ``jovian_shoot`` + Task 4b multi-rev
``ieg_eggie_seed``).  This scan sweeps ``departure_et`` finely over +/- 1
synodic period (~7 d) around the paper's 2020-Oct-02 EGGIE departure, building
the multi-rev Lambert-real seed at each step and recording its
``seed_defect_norm`` via ``jovian_defect_residual`` — the SAME residual the
corrector minimises (so the chosen best epoch is the corrector's warm-start).

A flushed JSONL runlog is written one line per epoch to
``/tmp/ieg_rerun_scan_480.jsonl``.
"""

from __future__ import annotations

import datetime as _dt
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
)
from cyclerfinder.search.ieg_seed import ieg_eggie_seed
from cyclerfinder.verify.spice_kernels import ensure_jup365_kernel, ensure_leapseconds_kernel

_SPD = 86400.0
_T_SYN_DAYS = 7.05  # EGGIE synodic period (digest Table 4).


def paper_departure_et() -> float:
    spiceypy.furnsh(ensure_leapseconds_kernel())
    return float(spiceypy.str2et("2020-OCT-02 12:00 TDB"))


def main() -> None:
    runlog = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("/tmp/ieg_rerun_scan_480.jsonl")
    step_days = float(sys.argv[2]) if len(sys.argv) > 2 else 0.1
    half_window_days = float(sys.argv[3]) if len(sys.argv) > 3 else _T_SYN_DAYS

    et0 = paper_departure_et()
    print(f"MARK paper_et {et0}", flush=True)

    moons = ("Io", "Europa", "Ganymede")
    jeph = JovianEphemeris(ensure_jup365_kernel())

    offsets = np.arange(-half_window_days, half_window_days + 0.5 * step_days, step_days)
    best = (float("inf"), float("nan"), float("nan"))

    with runlog.open("w") as fh:
        for k, off_d in enumerate(offsets):
            et = et0 + float(off_d) * _SPD
            seed = ieg_eggie_seed(departure_et=et)
            cache = JovianRailsCache(moons, jeph, min(seed.epochs), max(seed.epochs))
            res = jovian_defect_residual(
                seed, ephem=jeph, cache=cache, moons=moons, max_wall_sec=30.0
            )
            sdn = float(np.linalg.norm(res))
            ts = _dt.datetime.now(_dt.UTC).isoformat()
            record = {"et": et, "off_days": float(off_d), "seed_defect_norm": sdn, "ts": ts}
            fh.write(json.dumps(record) + "\n")
            fh.flush()
            if sdn < best[0]:
                best = (sdn, et, float(off_d))
            print(
                f"MARK scan {k + 1}/{len(offsets)} off={off_d:+.3f}d sdn={sdn:.6e}",
                flush=True,
            )

    print(
        f"MARK best seed_defect_norm={best[0]:.6e} at et={best[1]} off={best[2]:+.3f}d",
        flush=True,
    )


if __name__ == "__main__":
    t0 = time.time()
    main()
    print(f"MARK done wall={time.time() - t0:.1f}s", flush=True)
