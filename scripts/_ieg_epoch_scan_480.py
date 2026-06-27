"""#480 Task 4 — IEG EGGIE departure-epoch scan against real Galilean ephemeris.

Sweeps ``departure_et`` across roughly +/- 2 synodic periods around the paper's
2020-Oct-02 departure (Hernandez-Jones-Jesick 2017, AAS 17-608, Table 4 EGGIE),
building the Lambert-real seed at each step and recording its multiple-shooting
``seed_defect_norm`` (the L2 norm of the per-node continuity residual under the
n-body propagator).  The lowest-residual epoch is the warm-start for the
corrector.

Observability
-------------
A flushed JSONL runlog is written one line per epoch::

    {"et": <float>, "seed_defect_norm": <float>, "ts": <iso8601>}

so a long scan is observable line-by-line (per feedback_incremental_progress_reports).

SCIENCE FINDING (#480 Task 4 — the load-bearing negative)
---------------------------------------------------------
``cyclerfinder.nbody.shooter.shoot`` / ``defect_residual`` integrate about a
**heliocentric** central mass (``MU_SUN_KM3_S2``, ~1047x Jupiter's GM) with planet
perturbers looked up in ``core.constants.PLANETS``.  The IEG seed is **Jupiter-centred**
with Galilean moon nodes, which are NOT in ``PLANETS`` — so feeding it to ``shoot()``
raises ``KeyError`` inside the REBOUND ctypes callback (swallowed, leaving the moon
third-body terms dropped and the wrong central GM applied).  The Jupiter-correct
propagator (``nbody.jovian.JovianRestrictedNBody``, central ``MU_JUPITER``, moons on
rails) is NOT wired into ``shoot()``.  **The #480 corrector therefore cannot run on
the Jovian seed at all.**

This scan consequently uses the Jovian-correct ``JovianRestrictedNBody`` directly to
measure the per-epoch ``seed_defect_norm`` (the only physically meaningful residual
available).  Even with the correct propagator, the single-revolution Lambert EGGIE
seed at the paper epoch sits FAR from the basin (residual ~2.3e6, V_inf 3-27 km/s vs
the paper's 9.12/7.07/7.07/8.38) — a characterised negative: the construction does not
land in the EGGIE basin, independent of the corrector gap.
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
    JovianRestrictedNBody,
)
from cyclerfinder.search.ieg_seed import ieg_eggie_seed
from cyclerfinder.verify.spice_kernels import ensure_jup365_kernel, ensure_leapseconds_kernel

_SPD = 86400.0
_T_SYN_DAYS = 7.05  # EGGIE synodic period (digest Table 4).


def paper_departure_et() -> float:
    """SPICE ET (TDB s past J2000) for 2020-OCT-02 12:00 TDB."""
    spiceypy.furnsh(ensure_leapseconds_kernel())
    return float(spiceypy.str2et("2020-OCT-02 12:00 TDB"))


def main() -> None:
    runlog = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("/tmp/ieg_epoch_scan_480.jsonl")
    step_days = 0.25
    half_window_days = 2.0 * _T_SYN_DAYS  # +/- 2 synodic periods (~14 d).

    et0 = paper_departure_et()
    print(f"MARK paper_et {et0}", flush=True)

    moons = ("Io", "Europa", "Ganymede")
    jeph = JovianEphemeris(ensure_jup365_kernel())
    prop = JovianRestrictedNBody()

    offsets = np.arange(-half_window_days, half_window_days + 0.5 * step_days, step_days)
    best = (float("inf"), float("nan"))

    with runlog.open("w") as fh:
        for k, off_d in enumerate(offsets):
            et = et0 + float(off_d) * _SPD
            seed = ieg_eggie_seed(departure_et=et)
            cache = JovianRailsCache(moons, jeph, min(seed.epochs), max(seed.epochs))
            res: list[float] = []
            n = len(seed.sequence)
            for i in range(n - 1):
                s_i = seed.node_states[i]
                arc = prop.propagate(
                    np.asarray(s_i[:3]),
                    np.asarray(s_i[3:]),
                    t0_sec=seed.epochs[i],
                    t1_sec=seed.epochs[i + 1],
                    moons=moons,
                    cache=cache,
                    max_wall_sec=30.0,
                )
                s_next = seed.node_states[i + 1]
                res.extend(float(x) for x in (arc.r_km - np.asarray(s_next[:3])))
                res.extend(float(x) for x in (arc.v_km_s - np.asarray(s_next[3:])))
            sdn = float(np.linalg.norm(res))
            ts = _dt.datetime.now(_dt.UTC).isoformat()
            record = {"et": et, "off_days": float(off_d), "seed_defect_norm": sdn, "ts": ts}
            fh.write(json.dumps(record) + "\n")
            fh.flush()
            if sdn < best[0]:
                best = (sdn, et)
            print(f"MARK scan {k + 1}/{len(offsets)} off={off_d:+.2f}d sdn={sdn:.6e}", flush=True)

    print(f"MARK best seed_defect_norm={best[0]:.6e} at et={best[1]}", flush=True)


if __name__ == "__main__":
    t0 = time.time()
    main()
    print(f"MARK done wall={time.time() - t0:.1f}s", flush=True)
