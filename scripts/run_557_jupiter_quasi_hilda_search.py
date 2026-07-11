"""#557 Sun-Jupiter quasi-Hilda transient-capture ``quasi_cycler`` search.

Points #535's validated Hill-sphere-return detector
(``search/hill_sphere_return_detector.py``) at the Sun-Jupiter quasi-Hilda
population. Admission criterion settled in writing FIRST -- see
``docs/notes/2026-07-11-557-jupiter-quasi-hilda-admission-criterion.md`` (the
single source of truth this script implements exactly).

CRITICAL DESIGN POINT (plan Sec 5 / Fable review, the latent bug this avoids):
the criterion is DIMENSIONLESS, in units of one rotating-frame period (= 2*pi
nondim CR3BP time). #535's Earth script hardcodes ``t = years * 2*pi`` -- true
ONLY because Sun-Earth's period is 1 year. Cloning that literally for Jupiter
would integrate ~500 revolutions instead of ~55. Here EVERYTHING is in PERIODS;
years appear only in human-readable reporting via T_JUPITER_YEARS. One period =
ONE_PERIOD = 2*pi nondim, always.

Positive control: 82P/Gehrels 3, real JPL SBDB elements converted to a planar
rotating-frame IC (see the criterion note Sec 6). The anchor enters Jupiter's
Hill sphere (0.12 R_hill) -- pipeline validated -- but a single deep capture is
not >=3 distinct returns, so the anchor motivates rather than is the candidate.
The scan looks for recurrent-capture ICs in the neck-open band C in [3.00, 3.038]
(strictly below C_L1 = 3.0388, so the L1/L2 necks are open -- NOT the closed-neck
C=3.14 family the #527/#530/#531 entries already showed produces zero Hill
encounters by construction).
"""

from __future__ import annotations

import datetime
import pathlib
import sys
import time

_REPO = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_REPO / "src"))

import numpy as np  # noqa: E402
from scipy.integrate import solve_ivp  # noqa: E402

import cyclerfinder.core.cr3bp as cr3bp  # noqa: E402
from cyclerfinder.data.method_capability import MethodCapability  # noqa: E402
from cyclerfinder.data.preflight import PreflightBlockedError, preflight_search  # noqa: E402
from cyclerfinder.search.hill_sphere_return_detector import (  # noqa: E402
    find_admission_windows,
    find_returns,
)

# One rotating-frame period = one full revolution = 2*pi nondim CR3BP time.
# ALL timescales in this script are expressed in these period units.
ONE_PERIOD = 2.0 * np.pi

# Verified in-repo this session (see criterion note Sec 1): Jupiter's heliocentric
# period, used ONLY to convert period-units to years for human-readable reporting.
T_JUPITER_YEARS = 11.868

# Criterion (dimensionless, in periods) -- criterion note Secs 2-5.
MIN_SEPARATION_PERIODS = 1.0  # floor = 1 Jupiter period (11.868 yr)
WINDOW_LO_PERIODS, WINDOW_HI_PERIODS = 10.0, 15.0  # window = 10-15 periods (119-178 yr)
N_RETURNS_LO, N_RETURNS_HI = 3, 15
GEOMETRY_FACTOR = 3.0
HORIZON_PERIODS = 55.0  # ~4x the window midpoint, matching #535's horizon/window ratio
SAMPLES_PER_PERIOD = 500  # per-REVOLUTION density (never per calendar year)

# Neck-open capture band: C strictly below C_L1 = 3.0388 (criterion note Sec 1).
C_VALUES = (3.000, 3.005, 3.010, 3.015, 3.020, 3.025, 3.030, 3.035)
# x0 spanning the interior 3:2 resonance region down through the capture bubble
# (the Gehrels 3 anchor sits at x0 = 0.696).
X0_GRID = np.round(np.linspace(0.60, 0.92, 33), 5)

# 82P/Gehrels 3 anchor IC (criterion note Sec 6): sourced JPL SBDB elements
# -> vis-viva-at-perihelion rotating-frame IC. Jacobi C = 3.02943 matches the
# SBDB Tisserand 3.027 to 0.002 (independent cross-check).
ANCHOR_X0 = 0.69591
ANCHOR_XDOT0 = 0.0
ANCHOR_YDOT0 = 0.57308

_REGION_ID = "sun-jupiter-quasi-hilda-transient-capture-quasi-cycler"
_CRITERION_VERSION = "557-optionA-system-period-relative-2026-07-11"
_METHOD = MethodCapability(
    genome="Direct long-duration CR3BP propagation + repeated-Hill-sphere-return detection "
    "(search/hill_sphere_return_detector.py) over the Sun-Jupiter neck-open (C<C_L1) "
    "quasi-Hilda energy band -- targets aperiodic transient-capture trajectories, "
    "structurally distinct from the #527/#530/#531 periodic-orbit/manifold Hilda searches",
    corrector="none (no periodicity requirement by design -- see #557 admission criterion note)",
    capability_tags=frozenset(
        {"cr3bp", "ballistic", "coplanar", "single-arc", "hill-sphere-return-detection"}
    ),
    git_sha="working-tree",
)

_RUNLOG = _REPO / "data" / "scan_557_runlog.jsonl"


def _ts() -> str:
    return datetime.datetime.now().strftime("%Y-%m-%dT%H:%M:%S")


def _ydot0_from_jacobi(
    x0: float, xdot0: float, c_target: float, mu: float, sign: float
) -> float | None:
    r1 = abs(x0 + mu)
    r2 = abs(x0 - 1.0 + mu)
    rad = x0 * x0 + 2.0 * (1.0 - mu) / r1 + 2.0 * mu / r2 - c_target - xdot0 * xdot0
    if rad < 0.0:
        return None
    return float(sign) * float(np.sqrt(rad))


def _propagate_and_detect(
    system: cr3bp.CR3BPSystem, x0: float, xdot0: float, ydot0: float, r_hill: float
) -> tuple[int, list] | None:
    """Propagate ``HORIZON_PERIODS`` revolutions and apply the admission criterion.

    ``x0/xdot0/ydot0`` are a rotating-frame planar IC; the caller supplies
    ``ydot0`` directly (anchor) or via :func:`_ydot0_from_jacobi` (grid).
    """
    mu = float(system.mu)
    state0 = np.array([x0, 0.0, 0.0, xdot0, ydot0, 0.0])
    t_max = HORIZON_PERIODS * ONE_PERIOD
    t_eval = np.linspace(0.0, t_max, int(HORIZON_PERIODS * SAMPLES_PER_PERIOD))
    sol = solve_ivp(
        cr3bp.cr3bp_eom,
        (0.0, t_max),
        state0,
        args=(mu,),
        method="DOP853",
        rtol=1e-11,
        atol=1e-11,
        t_eval=t_eval,
        max_step=0.5,
    )
    if not sol.success:
        return None
    pos_rel_jupiter = np.stack([sol.y[0] - (1.0 - mu), sol.y[1]], axis=1)
    returns = find_returns(
        sol.t,
        pos_rel_jupiter,
        r_hill=r_hill,
        min_separation=MIN_SEPARATION_PERIODS * ONE_PERIOD,
    )
    windows = find_admission_windows(
        returns,
        float(sol.t[0]),
        float(sol.t[-1]),
        window_lo=WINDOW_LO_PERIODS * ONE_PERIOD,
        window_hi=WINDOW_HI_PERIODS * ONE_PERIOD,
        n_returns_lo=N_RETURNS_LO,
        n_returns_hi=N_RETURNS_HI,
        geometry_factor=GEOMETRY_FACTOR,
    )
    admissible = [w for w in windows if w.geometry_ok]
    return len(returns), admissible


def main() -> None:
    print(f"[{_ts()}] #557 Sun-Jupiter quasi-Hilda transient-capture search starting.", flush=True)

    system = cr3bp.cr3bp_system("Sun", "Jupiter")
    mu = float(system.mu)
    r_hill = (mu / 3.0) ** (1.0 / 3.0)
    print(
        f"[{_ts()}] mu={mu:.7e} r_hill={r_hill:.6f} nondim  "
        f"horizon={HORIZON_PERIODS:.0f} periods (~{HORIZON_PERIODS * T_JUPITER_YEARS:.0f} yr)  "
        f"window={WINDOW_LO_PERIODS:.0f}-{WINDOW_HI_PERIODS:.0f} periods "
        f"(~{WINDOW_LO_PERIODS * T_JUPITER_YEARS:.0f}-"
        f"{WINDOW_HI_PERIODS * T_JUPITER_YEARS:.0f} yr)  "
        f"floor={MIN_SEPARATION_PERIODS:.0f} period (~{T_JUPITER_YEARS:.1f} yr)",
        flush=True,
    )

    n_points = len(X0_GRID) * len(C_VALUES)
    preflight_search(
        task_no=557,
        region_id=_REGION_ID,
        method=_METHOD,
        script_path=pathlib.Path(__file__),
        n_points=n_points,
        # Measured this session by direct timing of one propagate+detect at the
        # 55-period horizon (see the preflight pilot in the #557 run log): ~1 s/point,
        # the same per-point cost as #535's Earth run (same nondim-revolution count).
        timing_pilot_seconds_per_point=1.5,
    )

    # Positive control: the sourced Gehrels 3 anchor must reproduce a Hill-sphere
    # capture (criterion note Sec 6). It need NOT be admissible (a single deep
    # capture is not >=3 distinct returns) -- but zero Hill entries would mean the
    # pipeline is broken; stop in that case.
    print(f"[{_ts()}] Positive control: Gehrels 3 anchor (x0={ANCHOR_X0}) ...", flush=True)
    anchor = _propagate_and_detect(system, ANCHOR_X0, ANCHOR_XDOT0, ANCHOR_YDOT0, r_hill)
    if anchor is None or anchor[0] == 0:
        raise RuntimeError(
            "POSITIVE CONTROL FAILED: the Gehrels 3 anchor produced zero Hill-sphere "
            "returns. The elements->IC->propagation pipeline is broken; do not trust the scan."
        )
    n_anchor, adm_anchor = anchor
    print(
        f"[{_ts()}] Positive control PASSED: anchor has {n_anchor} Hill-sphere return(s), "
        f"{len(adm_anchor)} admissible window(s) (a single deep capture is expected, "
        f"admissibility is not required of the anchor).",
        flush=True,
    )

    # Fresh run log (checkpoint per point: append + flush, running counts, ETA).
    _RUNLOG.write_text("", encoding="utf-8")
    t0 = time.time()
    findings: list[dict] = []
    done = 0
    for c_target in C_VALUES:
        for x0 in X0_GRID:
            ydot0 = _ydot0_from_jacobi(float(x0), 0.0, float(c_target), mu, 1.0)
            done += 1
            rec: dict = {
                "ts": _ts(),
                "i": done,
                "n": n_points,
                "x0": float(x0),
                "c_target": float(c_target),
            }
            if ydot0 is None:
                rec["status"] = "no-real-ydot0"
            else:
                result = _propagate_and_detect(system, float(x0), 0.0, ydot0, r_hill)
                if result is None:
                    rec["status"] = "integration-failed"
                else:
                    n_returns, admissible = result
                    rec["status"] = "ok"
                    rec["n_returns"] = n_returns
                    rec["n_admissible"] = len(admissible)
                    if admissible:
                        best = min(w.geometry_ratio for w in admissible)
                        rec["best_geometry_ratio"] = float(best)
                        finding = {
                            "x0": float(x0),
                            "c_target": float(c_target),
                            "ydot0": float(ydot0),
                            "n_returns_total": n_returns,
                            "n_admissible_windows": len(admissible),
                            "best_geometry_ratio": float(best),
                        }
                        findings.append(finding)
                        print(
                            f"[{_ts()}]   ADMISSIBLE: x0={x0:.5f} C={c_target:.4f} "
                            f"n_returns={n_returns} n_windows={len(admissible)} "
                            f"geom={best:.3f}",
                            flush=True,
                        )
            elapsed = time.time() - t0
            eta = (elapsed / done) * (n_points - done)
            rec["elapsed_s"] = round(elapsed, 1)
            rec["eta_s"] = round(eta, 1)
            rec["admissible_so_far"] = len(findings)
            with _RUNLOG.open("a", encoding="utf-8") as fh:
                import json

                fh.write(json.dumps(rec) + "\n")
                fh.flush()
            if done % 20 == 0 or done == n_points:
                print(
                    f"[{_ts()}] progress {done}/{n_points}  "
                    f"admissible={len(findings)}  elapsed={elapsed:.0f}s  eta={eta:.0f}s",
                    flush=True,
                )

    dt = time.time() - t0
    print(flush=True)
    print(
        f"[{_ts()}] Scan complete in {dt:.1f}s over {n_points} points "
        f"(criterion version {_CRITERION_VERSION}). Admissible candidates: {len(findings)}.",
        flush=True,
    )
    for f in findings:
        print(f"    {f}", flush=True)
    if not findings:
        print(
            f"[{_ts()}] CLEAN NULL: no admissible transient-capture quasi_cycler in the "
            f"neck-open band under the Option-A system-period-relative criterion. "
            f"Register in data/empty_regions.jsonl with this criterion version.",
            flush=True,
        )


if __name__ == "__main__":
    try:
        main()
    except PreflightBlockedError as exc:
        print(f"[{_ts()}] BLOCKED by preflight_search:\n{exc}", flush=True)
        sys.exit(1)
