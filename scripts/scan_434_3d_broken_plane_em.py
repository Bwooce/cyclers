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
  (ii) z0-amplitude lock per planar root -- for each planar Aldrin/Braik-Ross
       root, step |z0| up from a MODERATE start in {0.05,0.10,0.15,0.20,0.24}
       and let ``correct_general_periodic_3d`` find the 3D branch (try
       symmetric-tulip free vars first; fall back to full-asymmetric); record
       which roots HAVE a 3D extension (lock) vs collapse-only.

For each converged 3D seed the sweep runs the pseudo-arclength family tracer and
logs every closure-verified member.

Report-only -- NO catalogue writeback. Per the literature_check discipline a
"not-found" literature status is NECESSARY-NOT-SUFFICIENT for novelty; the
V0-V5 gauntlet still governs (Task 5, controller).

Usage::

    uv run python scripts/scan_434_3d_broken_plane_em.py            # full sweep
    SCAN_434_SMOKE=1 uv run python scripts/scan_434_3d_broken_plane_em.py  # tiny

The smoke mode (env var ``SCAN_434_SMOKE=1`` or ``--smoke``) runs only the L1
vertical-Lyapunov seed + one planar root with ``n_steps_max=20`` to confirm the
pipeline lifts, converges, continues, and writes records. The full sweep is
launched by the controller (Task 5), NOT by the smoke.
"""

from __future__ import annotations

import json
import sys
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from cyclerfinder.core.cr3bp import CR3BPSystem, cr3bp_system
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

_DATA_DIR = Path(__file__).resolve().parents[1] / "data"
_OUT_PATH = _DATA_DIR / "scan_434_3d_broken_plane_em.jsonl"

# Moderate-start |z0| ladder for route (ii). Small-z0 lifts collapse to the
# planar manifold (spike caveat), so we start at 0.05 and the spike's known
# 3D-branch amplitude (~0.24) is the last rung.
_Z0_LADDER: tuple[float, ...] = (0.05, 0.10, 0.15, 0.20, 0.24)

# Collinear points whose vertical (out-of-plane) mode is oscillatory.
_LYAPUNOV_POINTS: tuple[str, ...] = ("L1", "L2")


# Planar Aldrin/Braik-Ross Earth-Moon roots (rotating-frame CR3BP ICs).
#
# SOURCING: the catalogue rows (data/catalogue.yaml) carry HELIOCENTRIC mission
# elements, not rotating-frame CR3BP ICs, so they are not machine-loadable as
# corrector seeds. The canonical planar Earth-Moon Braik-Ross C11a (1,1) root is
# taken from the #287 spike's published-on-disk planar baseline
# (data/spike_287.jsonl, case "planar_baseline_z0eq0"): the perpendicular-
# crossing IC (x0, 0, 0, 0, ydot0, 0) at jacobi 3.1294, the same baseline used by
# the #287 spike and the v3_3d / v2_3d Braik-Ross (1,1) validation closures.
@dataclass(frozen=True)
class PlanarRoot:
    label: str
    x0: float
    ydot0: float
    period: float
    source: str


_PLANAR_ROOTS: tuple[PlanarRoot, ...] = (
    PlanarRoot(
        label="braik-ross-C11a-em-k1",
        x0=-0.8116406668238195,
        ydot0=-0.11859055759763637,
        period=9.69107744379376,
        source=(
            "#287 spike planar baseline (data/spike_287.jsonl case "
            "planar_baseline_z0eq0); Braik-Ross (1,1) Earth-Moon C11a at jacobi 3.1294"
        ),
    ),
)


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
    root: PlanarRoot,
    z0: float,
) -> Periodic3DOrbit | None:
    """Attempt the z0-amplitude lock for one (root, z0); return a non-degenerate
    3D orbit, or None on collapse / non-convergence.

    Tries the symmetric-tulip free vars first; if that collapses to the planar
    manifold, retries with the full-asymmetric mask.
    """
    seed_state = np.array([root.x0, 0.0, z0, 0.0, root.ydot0, 0.0], dtype=np.float64)
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
    roots = (_PLANAR_ROOTS[0],) if smoke else _PLANAR_ROOTS
    z0_ladder = (_Z0_LADDER[0],) if smoke else _Z0_LADDER

    _print_progress(
        f"#434 3D broken-plane discovery (Earth-Moon){' [SMOKE]' if smoke else ''} "
        f"n_steps_max={n_steps_max}"
    )
    system = cr3bp_system("Earth", "Moon")

    all_records: list[dict[str, object]] = []
    per_seed_counts: list[tuple[str, int]] = []

    # --- Route (i): vertical-Lyapunov / halo seeds --------------------------
    for point in points:
        label = f"lyapunov3d-{point}"
        try:
            state0, period = lyapunov_seed_3d(system, point=point)
        except ValueError as exc:
            _print_progress(f"route(i) {label}: seed FAILED ({exc}); skipping")
            continue
        _print_progress(
            f"route(i) {label}: seed locked z0={state0[2]:.4f} T={period:.4f}; continuing"
        )
        family = _continue_seed(
            system, state0, period, n_steps_max=n_steps_max, monodromy_eval=not smoke
        )
        recs = _records_for_family(system, family, "lyapunov3d", label)
        all_records.extend(recs)
        per_seed_counts.append((label, len(recs)))
        _print_progress(f"route(i) {label}: {len(recs)} converged member(s)")

    # --- Route (ii): z0-amplitude lock per planar root ----------------------
    lock_tally: dict[str, tuple[int, int]] = {}  # label -> (locks, collapses)
    for root in roots:
        locks = 0
        collapses = 0
        for z0 in z0_ladder:
            orbit = _lock_planar_root(system, root, -z0)  # negative z0 per spike branch
            if orbit is None:
                collapses += 1
                _print_progress(f"route(ii) {root.label} |z0|={z0:.2f}: COLLAPSE / no 3D lock")
                continue
            locks += 1
            seed_label = f"{root.label}-z0_{z0:.2f}"
            _print_progress(
                f"route(ii) {seed_label}: LOCK z0={orbit.state0[2]:.4f} "
                f"T={orbit.T_TU:.4f}; continuing"
            )
            family = _continue_seed(
                system,
                orbit.state0,
                orbit.T_TU,
                n_steps_max=n_steps_max,
                monodromy_eval=not smoke,
            )
            recs = _records_for_family(system, family, "z0_lift", seed_label)
            all_records.extend(recs)
            per_seed_counts.append((seed_label, len(recs)))
            _print_progress(f"route(ii) {seed_label}: {len(recs)} converged member(s)")
        lock_tally[root.label] = (locks, collapses)

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
    import os

    smoke_mode = "--smoke" in sys.argv or os.environ.get("SCAN_434_SMOKE") == "1"
    main(smoke=smoke_mode)
