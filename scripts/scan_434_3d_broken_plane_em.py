"""#434 — 3D broken-plane cycler discovery sweep (Earth-Moon).

Campaign goal: discover novel out-of-plane (z != 0) cycler families by
continuing genuinely-3D orbits and screening for novel structure.

SPIKE-CONFIRMED CAVEAT (Task 0, 2026-06-24): small-z0 lifts of a planar root
COLLAPSE back to the planar manifold (|z0| -> < 1e-14); the 3D branch only locks
in from a seed already well out-of-plane (the #287 spike's C11a 3D branch locks
near z0 ~ -0.24). So this sweep uses TWO seed routes that genuinely lock onto 3D
structure:

  (i)  vertical-Lyapunov / halo generator -- ``lyapunov_seed_3d`` produces
       genuinely out-of-plane ICs at the collinear point; continue each.
  (ii) z0-amplitude lock per planar cycler -- recover ALL FOUR Braik-Ross
       Earth-Moon cyclers (C11a, C11b, C21, C32) via
       ``recover_all_cyclers_braik_ross`` and, for each CONFIRMED planar root,
       step |z0| up from a MODERATE start in {0.05,0.10,0.15,0.20,0.24} and let
       ``correct_general_periodic_3d`` find the 3D branch (try symmetric-tulip
       free vars first; fall back to full-asymmetric); record which cyclers HAVE
       a 3D extension (lock) vs collapse-only. C11b/C21/C32 carry DIFFERENT
       planar winding topology than C11a, so their 3D extensions are the genuine
       novel-family candidates (#438 unlock of the #434 frontier).

For each converged 3D seed the sweep runs the pseudo-arclength family tracer and
logs every closure-verified member.

Report-only -- NO catalogue writeback. Per the literature_check discipline a
"not-found" literature status is NECESSARY-NOT-SUFFICIENT for novelty; the
V0-V5 gauntlet still governs (Task 5, controller).

Usage::

    uv run python scripts/scan_434_3d_broken_plane_em.py            # full sweep
    SCAN_434_SMOKE=1 uv run python scripts/scan_434_3d_broken_plane_em.py  # tiny

The smoke mode (env var ``SCAN_434_SMOKE=1`` or ``--smoke``) recovers the four
cyclers + lifts only the FIRST confirmed one at one z0 (plus the L1 vertical-
Lyapunov seed) with ``n_steps_max=20`` to confirm the pipeline recovers, lifts,
converges, continues, and writes records. The full sweep is launched by the
controller (Task 5), NOT by the smoke.
"""

from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path

import numpy as np
from joblib import Parallel, delayed

from cyclerfinder.core.cr3bp import CR3BPSystem
from cyclerfinder.search.binary_star_search import topology_3d
from cyclerfinder.search.cr3bp_3d_family_tracer import Family3D, continue_general_3d_family
from cyclerfinder.search.cr3bp_general_periodic_3d import (
    FREE_VARS_FULL_ASYMMETRIC,
    FREE_VARS_SYMMETRIC_TULIP,
    RESIDUAL_PERPENDICULAR_HALF_PERIOD,
    Periodic3DOrbit,
    correct_general_periodic_3d,
)
from cyclerfinder.search.cr3bp_seed_generator import lyapunov_seed_3d
from cyclerfinder.search.reachable_representatives import (
    Representative,
    braik_ross_system,
    recover_all_cyclers_braik_ross,
)

_DATA_DIR = Path(__file__).resolve().parents[1] / "data"
_OUT_PATH = _DATA_DIR / "scan_434_3d_broken_plane_em.jsonl"

# Moderate-start |z0| ladder for route (ii). Small-z0 lifts collapse to the
# planar manifold (spike caveat), so we start at 0.05 and the spike's known
# 3D-branch amplitude (~0.24) is the last rung.
_Z0_LADDER: tuple[float, ...] = (0.05, 0.10, 0.15, 0.20, 0.24)

# Collinear points whose vertical (out-of-plane) mode is oscillatory.
_LYAPUNOV_POINTS: tuple[str, ...] = ("L1", "L2")

# Parallel workers for the independent per-family lift+continue units. Each unit
# (one route-(i) libration point, or one route-(ii) (cycler, z0) pair) is fully
# independent, so they fan out across cores. loky's backend auto-caps inner
# (BLAS) threads per worker, so no manual thread-env capping is needed.
_N_JOBS = int(os.environ.get("CYCLERS_434_NJOBS", "12"))


# Planar Braik-Ross Earth-Moon cycler roots are recovered at run time via
# ``recover_all_cyclers_braik_ross`` (C11a, C11b, C21, C32). Each Representative
# carries the perpendicular-crossing planar IC (x0, 0, 0, 0, ydot0, 0), the
# nondim full period, the per-family Jacobi, and a ``.confirmed`` flag (period
# matched the Braik-Ross Table-2 sourced value within tolerance). ONLY confirmed
# cyclers are lifted; the rest are logged. The four cyclers do NOT share planar
# winding topology, so C11b/C21/C32's 3D extensions are the genuine novel-family
# candidates the #438 broadening is after.
def _print_progress(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def _topology_counts(system: CR3BPSystem, orbit: Periodic3DOrbit) -> tuple[int, int, int]:
    try:
        topo = topology_3d(system.mu, orbit.state0, orbit.T_TU)
    except Exception:  # topology is diagnostic-only; never abort the sweep
        return (-1, -1, -1)
    return (topo.k1, topo.k2, topo.k_z)


def _records_for_family(
    system: CR3BPSystem,
    family: Family3D,
    seed_route: str,
    seed_label: str,
) -> list[dict[str, object]]:
    records: list[dict[str, object]] = []
    for m in family.members:
        orbit = m.orbit
        k1, k2, k_z = _topology_counts(system, orbit)
        records.append(
            {
                "seed_route": seed_route,
                "seed_label": seed_label,
                "step_index": m.step_index,
                "x0": float(orbit.state0[0]),
                "z0": float(orbit.state0[2]),
                "ydot0": float(orbit.state0[4]),
                "T_TU": float(orbit.T_TU),
                "jacobi": float(orbit.jacobi),
                "k1": k1,
                "k2": k2,
                "k_z": k_z,
                "corrector_residual": float(orbit.corrector_residual),
                "independent_closure_residual": float(orbit.independent_closure_residual),
                "floquet_tag": m.stability_tag,
            }
        )
    return records


def _continue_seed(
    system: CR3BPSystem,
    seed_state: np.ndarray,
    seed_period: float,
    *,
    n_steps_max: int,
    monodromy_eval: bool,
) -> Family3D:
    return continue_general_3d_family(
        system,
        seed_state,
        seed_period,
        free_vars=FREE_VARS_SYMMETRIC_TULIP,
        residual_indices=RESIDUAL_PERPENDICULAR_HALF_PERIOD,
        continuation="pseudo_arclength",
        step=0.01,
        n_steps_max=n_steps_max,
        direction="both",
        monodromy_eval=monodromy_eval,
    )


def _lock_planar_root(
    system: CR3BPSystem,
    root: Representative,
    z0: float,
) -> Periodic3DOrbit | None:
    """Attempt the z0-amplitude lock for one (cycler, z0); return a non-degenerate
    3D orbit, or None on collapse / non-convergence.

    Tries the symmetric-tulip free vars first; if that collapses to the planar
    manifold, retries with the full-asymmetric mask.
    """
    x0 = float(root.state0[0])
    ydot0 = float(root.state0[4])
    seed_state = np.array([x0, 0.0, z0, 0.0, ydot0, 0.0], dtype=np.float64)
    # Symmetric-tulip first (perpendicular half-period residual); fall back to
    # the full-asymmetric mask (full-state closure at T) if it collapses.
    tulip = correct_general_periodic_3d(
        system,
        seed_state,
        root.period,
        free_vars=FREE_VARS_SYMMETRIC_TULIP,
        residual_indices=RESIDUAL_PERPENDICULAR_HALF_PERIOD,
        is_half_period_residual=True,
    )
    if tulip.converged and not tulip.degenerate_planar:
        return tulip
    asym = correct_general_periodic_3d(
        system,
        seed_state,
        root.period,
        free_vars=FREE_VARS_FULL_ASYMMETRIC,
        is_half_period_residual=False,
    )
    if asym.converged and not asym.degenerate_planar:
        return asym
    return None


def _worker_lyapunov(
    system: CR3BPSystem,
    point: str,
    *,
    n_steps_max: int,
    monodromy_eval: bool,
) -> dict[str, object]:
    """Route-(i) work unit: one libration point -> lyapunov_seed_3d -> continue.

    Runs in a joblib worker process. Returns a small picklable result dict; the
    parent aggregates counts/records (no shared mutable state).
    """
    label = f"lyapunov3d-{point}"
    try:
        state0, period = lyapunov_seed_3d(system, point=point)
    except ValueError as exc:
        return {"kind": "lyapunov", "label": label, "status": "seed_failed", "detail": str(exc)}
    family = _continue_seed(
        system, state0, period, n_steps_max=n_steps_max, monodromy_eval=monodromy_eval
    )
    recs = _records_for_family(system, family, "lyapunov3d", label)
    return {
        "kind": "lyapunov",
        "label": label,
        "status": "ok",
        "z0": float(state0[2]),
        "period": float(period),
        "records": recs,
    }


def _worker_z0_lift(
    system: CR3BPSystem,
    root: Representative,
    z0: float,
    *,
    n_steps_max: int,
    monodromy_eval: bool,
) -> dict[str, object]:
    """Route-(ii) work unit: one (cycler, z0) -> lock -> if locked, continue.

    Runs in a joblib worker process. Returns a small picklable result dict; the
    parent aggregates the lock/collapse tally and records (no shared mutable
    state).
    """
    cycler_label = f"braik-ross-{root.label}-em"
    orbit = _lock_planar_root(system, root, -z0)  # negative z0 per spike branch
    if orbit is None:
        return {
            "kind": "z0_lift",
            "cycler": root.label,
            "cycler_label": cycler_label,
            "z0": z0,
            "status": "collapse",
        }
    seed_label = f"{cycler_label}-z0_{z0:.2f}"
    family = _continue_seed(
        system,
        orbit.state0,
        orbit.T_TU,
        n_steps_max=n_steps_max,
        monodromy_eval=monodromy_eval,
    )
    recs = _records_for_family(system, family, "z0_lift", seed_label)
    return {
        "kind": "z0_lift",
        "cycler": root.label,
        "cycler_label": cycler_label,
        "z0": z0,
        "status": "lock",
        "seed_label": seed_label,
        "locked_z0": float(orbit.state0[2]),
        "locked_T": float(orbit.T_TU),
        "records": recs,
    }


def _closure_summary(records: list[dict[str, object]]) -> dict[str, float]:
    if not records:
        return {"n": 0, "max": float("nan"), "median": float("nan")}
    vals = np.array(
        [float(r["independent_closure_residual"]) for r in records],  # type: ignore[arg-type]
        dtype=np.float64,
    )
    return {
        "n": int(vals.size),
        "max": float(vals.max()),
        "median": float(np.median(vals)),
    }


def main(*, smoke: bool = False) -> None:
    t0 = time.time()
    n_steps_max = 20 if smoke else 100
    points = (_LYAPUNOV_POINTS[0],) if smoke else _LYAPUNOV_POINTS
    z0_ladder = (_Z0_LADDER[0],) if smoke else _Z0_LADDER

    _print_progress(
        f"#434 3D broken-plane discovery (Earth-Moon){' [SMOKE]' if smoke else ''} "
        f"n_steps_max={n_steps_max}"
    )
    # Use the Braik-Ross / Ross-RT nondimensional Earth-Moon system (the mu the
    # cycler recovery is pinned to). The (2,1) C21 family has Jacobi extent
    # ~4e-12 (a single point in C) and only converges at the Braik-Ross mu; the
    # generic cr3bp_system("Earth","Moon") mu is ~1.2e-10 off and drops C21,
    # which defeats the #438 broadening. The mu difference is negligible for the
    # L1 vertical-Lyapunov route-(i) seed.
    system = braik_ross_system()

    # --- Recover the four planar Braik-Ross cyclers (route (ii) seed source) -
    recovered = recover_all_cyclers_braik_ross(system)
    confirmed_roots: list[Representative] = []
    for rep in recovered:
        status = (
            "CONFIRMED"
            if rep.confirmed
            else ("converged-UNCONFIRMED" if rep.converged else "NOT-CONVERGED")
        )
        _print_progress(
            f"recover {rep.label}: {status} "
            f"period={rep.period_days:.3f}d (sourced {rep.sourced_period_days:.3f}d) "
            f"jacobi={rep.jacobi:.6f} x0={rep.state0[0]:.6f} ydot0={rep.state0[4]:.6f}"
        )
        if rep.confirmed:
            confirmed_roots.append(rep)
    if smoke:
        confirmed_roots = confirmed_roots[:1]

    # --- Build the flat list of INDEPENDENT lift+continue work units --------
    # route (i): one libration point each; route (ii): one (cycler, z0) each.
    monodromy_eval = not smoke
    n_route_i = len(points)
    n_route_ii = len(confirmed_roots) * len(z0_ladder)
    _print_progress(
        f"dispatching {n_route_i + n_route_ii} independent units "
        f"(route(i)={n_route_i}, route(ii)={n_route_ii}) across n_jobs={_N_JOBS}"
    )

    jobs = [
        delayed(_worker_lyapunov)(
            system, point, n_steps_max=n_steps_max, monodromy_eval=monodromy_eval
        )
        for point in points
    ]
    jobs += [
        delayed(_worker_z0_lift)(
            system, root, z0, n_steps_max=n_steps_max, monodromy_eval=monodromy_eval
        )
        for root in confirmed_roots
        for z0 in z0_ladder
    ]
    results: list[dict[str, object]] = Parallel(n_jobs=_N_JOBS, backend="loky")(jobs)

    # --- Aggregate in the parent (order may differ; that's fine) ------------
    all_records: list[dict[str, object]] = []
    per_seed_counts: list[tuple[str, int]] = []
    lock_tally: dict[str, tuple[int, int]] = {  # cycler -> (locks, collapses)
        root.label: (0, 0) for root in confirmed_roots
    }

    for res in results:
        if res["kind"] == "lyapunov":
            label = str(res["label"])
            if res["status"] == "seed_failed":
                _print_progress(f"route(i) {label}: seed FAILED ({res['detail']}); skipping")
                continue
            recs = res["records"]
            assert isinstance(recs, list)
            _print_progress(
                f"route(i) {label}: seed locked z0={float(res['z0']):.4f} "  # type: ignore[arg-type]
                f"T={float(res['period']):.4f}; {len(recs)} converged member(s)"  # type: ignore[arg-type]
            )
            all_records.extend(recs)
            per_seed_counts.append((label, len(recs)))
        else:  # z0_lift
            cycler = str(res["cycler"])
            locks, collapses = lock_tally[cycler]
            if res["status"] == "collapse":
                lock_tally[cycler] = (locks, collapses + 1)
                _print_progress(
                    f"route(ii) {res['cycler_label']} |z0|={float(res['z0']):.2f}: "  # type: ignore[arg-type]
                    f"COLLAPSE / no 3D lock"
                )
                continue
            lock_tally[cycler] = (locks + 1, collapses)
            seed_label = str(res["seed_label"])
            recs = res["records"]
            assert isinstance(recs, list)
            _print_progress(
                f"route(ii) {seed_label}: LOCK z0={float(res['locked_z0']):.4f} "  # type: ignore[arg-type]
                f"T={float(res['locked_T']):.4f}; {len(recs)} converged member(s)"  # type: ignore[arg-type]
            )
            all_records.extend(recs)
            per_seed_counts.append((seed_label, len(recs)))

    # --- Write JSONL --------------------------------------------------------
    _OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    with _OUT_PATH.open("w") as f:
        for rec in all_records:
            f.write(json.dumps(rec) + "\n")
    _print_progress(f"Wrote {_OUT_PATH.relative_to(_DATA_DIR.parent)} ({len(all_records)} records)")

    # --- Report -------------------------------------------------------------
    _print_progress("per-seed family member counts:")
    for label, n in per_seed_counts:
        _print_progress(f"  {label}: {n}")

    _print_progress("route-(ii) lock-vs-collapse tally:")
    for label, (locks, collapses) in lock_tally.items():
        _print_progress(f"  {label}: locks={locks} collapses={collapses}")

    summary = _closure_summary(all_records)
    _print_progress(
        f"closure distribution: n={summary['n']} "
        f"max={summary['max']:.3e} median={summary['median']:.3e}"
    )
    _print_progress(f"Campaign complete in {time.time() - t0:.1f}s")


if __name__ == "__main__":
    smoke_mode = (
        "--smoke" in sys.argv
        or os.environ.get("SCAN_434_SMOKE") == "1"
        or os.environ.get("CYCLERS_434_SMOKE") == "1"
    )
    main(smoke=smoke_mode)
