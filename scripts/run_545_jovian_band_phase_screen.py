#!/usr/bin/env python3
"""Search runner for Task #545: Jupiter-Europa L1/L2 Jacobi-band + phase-offset screen.

Task #536 tested the Jupiter-Europa L1<->L2 quasi-halo torus heteroclinic-
connection question at exactly ONE Jacobi constant (C=3.0015) and found 0
connections -- per this project's "no X found is conditional on the search
formulation" discipline, a single-point probe does not certify that region
empty. This script generalizes #536/#534's single-point pipeline into a
genuine BAND sweep:

  * a Jacobi-constant BAND (multiple matched-C L1/L2 torus pairs, built via
    ``search.nrho_continuation.continue_nrho_family`` natural-parameter
    continuation -- cheap, since NRHO correction itself is ~0.1-0.3s; the
    EXPENSIVE step, ~60-90s, is ``genome.qp_tori.correct_qp_torus`` per
    member, so torus corrections are checkpointed to disk (see
    ``_CHECKPOINT_DIR``, gitignored under ``data/runlogs/``) and this script
    is safely re-invocable/resumable across many short sessions without
    re-paying an already-converged torus's cost);
  * crossed with a SYNODIC-PHASE-OFFSET sweep: the CR3BP is autonomous, so
    the only meaningful extra phase freedom between an independently-built
    departure (L1) torus and arrival (L2) torus is WHERE along the L1
    torus's own longitudinal family angle the manifold is launched, relative
    to the L2 reference grid's own zero. This script adds a
    ``theta0_offset`` shift to the departure grid's longitudinal sampling
    (0 in the original #534/#536 scripts) and sweeps it explicitly;
  * with candidate connections refined via #524's deflated-Newton root
    enumerator (``search.deflated_newton.enumerate_roots``) on a genuine
    scalar residual -- ``genome.qp_torus_heteroclinic.closest_curve_distance``
    (new, added for #545) -- instead of just reporting the coarse linear
    scan's sign-change midpoints as candidates.

REQUIRED POSITIVE CONTROL (run first, always): the identical band+phase
pipeline is pointed at the Earth-Moon L1/L2 pair task #534 already
established and committed (``scripts/run_534_torus_connection.py``,
mu=0.012153643, the same seeds) -- if the generalized pipeline cannot
reproduce #534's own crossing/linking-number behavior on a system already
characterized, a Jovian negative from the same pipeline would be
meaningless.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import json
import math
import pathlib
import pickle
import time
from dataclasses import dataclass
from typing import Any

import numpy as np
import scipy.integrate
from numpy.typing import NDArray

import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.data.method_capability import MethodCapability
from cyclerfinder.data.preflight import preflight_search
from cyclerfinder.genome.qp_tori import QPTorus, correct_qp_torus
from cyclerfinder.genome.qp_torus_heteroclinic import (
    closest_curve_distance,
    scan_linking_number,
)
from cyclerfinder.genome.qp_torus_manifold import ManifoldGrid, local_stability, torus_point_stm
from cyclerfinder.search.bifurcation_detector import floquet_multipliers, monodromy
from cyclerfinder.search.deflated_newton import enumerate_roots
from cyclerfinder.search.nrho_continuation import (
    SymmetricNRHO,
    continue_nrho_family,
    correct_symmetric_nrho,
)

_REPO_ROOT = pathlib.Path(__file__).resolve().parent.parent
_CHECKPOINT_DIR = _REPO_ROOT / "data" / "runlogs" / "545_torus_checkpoint"
_REPORT_PATH = _REPO_ROOT / "data" / "runlogs" / "545_band_screen_report.json"

_REGION_ID = "jupiter-europa-l1-l2-torus-jacobi-band-phase-offset-2026-07-10"
_METHOD = MethodCapability(
    genome=(
        "Jupiter-Europa L1/L2 quasi-halo torus stable/unstable manifold crossing grids "
        "(10x10 long/lat) swept over a Jacobi-constant BAND (3 matched-C pairs via "
        "continue_nrho_family natural-parameter continuation) crossed with a "
        "synodic longitudinal phase-offset sweep (4 offsets), candidate connections "
        "refined by #524 deflated-Newton root enumeration on the closest-curve-distance "
        "residual -- strictly supersedes #536's single-point (C=3.0015, no phase sweep, "
        "no root refinement) negative"
    ),
    corrector=(
        "continue_nrho_family + build_manifold_grids-style custom grid + "
        "scan_linking_number + deflated_newton.enumerate_roots on closest_curve_distance"
    ),
    capability_tags=frozenset(
        {
            "cr3bp",
            "qp-torus",
            "heteroclinic",
            "linking-number",
            "jovian",
            "jacobi-band",
            "phase-offset",
            "deflated-newton",
        }
    ),
    git_sha="working-tree",
)

_N_LONG, _N_LAT = 10, 10
_T_MAX_UNSTABLE, _T_MAX_STABLE = 15.0, 25.0
_D_SCAN_N = 30
_PHASE_OFFSETS = (0.0, math.pi / 2.0, math.pi, 3.0 * math.pi / 2.0)


@dataclass(frozen=True)
class SeedSpec:
    x0: float
    z0: float
    ydot0: float
    period: float
    k: int


EM_MU = 0.012153643
EM_SEED_L1 = SeedSpec(0.836314, -0.145689, 0.231349, 2.751255, k=5)
EM_SEED_L2 = SeedSpec(1.023731, 0.183250, -0.106950, 1.533637, k=7)
JE_SEED_L1 = SeedSpec(0.980935577267, 0.019677883305, 0.030609455324, 2.87928189, k=5)
JE_SEED_L2 = SeedSpec(1.018932698495, 0.020127929170, -0.030157907163, 2.91762625, k=7)


def _em_system() -> cr3bp.CR3BPSystem:
    return cr3bp.CR3BPSystem(
        mu=EM_MU, primary="Earth", secondary="Moon", l_km=384400.0, t_s=382981.0
    )


def _je_system() -> cr3bp.CR3BPSystem:
    return cr3bp.cr3bp_system("Jupiter", "Europa")


# ---------------------------------------------------------------------------
# Step 1: Jacobi-band construction (cheap -- NRHO correction only).
# ---------------------------------------------------------------------------


def _build_family_members(
    system: cr3bp.CR3BPSystem, seed: SeedSpec, *, d_x0: float, n_steps: int
) -> list[SymmetricNRHO]:
    """Bidirectional natural-parameter x0-continuation, seed included once."""
    seed_orbit = correct_symmetric_nrho(system, seed.x0, seed.z0, seed.ydot0, seed.period)
    members = [seed_orbit]
    for direction in (1, -1):
        branch = continue_nrho_family(
            seed_orbit, system, direction=direction, d_x0=d_x0, n_steps_max=n_steps
        )
        members.extend(m for m in branch.members[1:] if m.converged)
    return members


def _match_band(
    l1_members: list[SymmetricNRHO], l2_members: list[SymmetricNRHO], n_c: int
) -> list[tuple[SymmetricNRHO, SymmetricNRHO, float]]:
    """Nearest-neighbor-match ``n_c`` evenly spaced target Jacobi constants
    within the overlap of the two achieved C ranges. Returns
    ``(l1_member, l2_member, delta_C)`` triples, ``delta_C`` reported honestly
    (this is nearest-neighbor matching on a coarse continuation grid, NOT
    #534's tight ~1e-5 isoenergetic refinement -- adequate for a first-pass
    band SCREEN, not a final converged connection).
    """
    c1 = np.array([m.jacobi for m in l1_members])
    c2 = np.array([m.jacobi for m in l2_members])
    lo = max(c1.min(), c2.min())
    hi = min(c1.max(), c2.max())
    if lo >= hi:
        raise ValueError(
            f"no Jacobi-band overlap: L1 range [{c1.min()},{c1.max()}], "
            f"L2 range [{c2.min()},{c2.max()}]"
        )
    targets = np.linspace(lo, hi, n_c)
    out = []
    seen: set[tuple[int, int]] = set()
    for target in targets:
        i1 = int(np.argmin(np.abs(c1 - target)))
        i2 = int(np.argmin(np.abs(c2 - target)))
        if (i1, i2) in seen:
            continue
        seen.add((i1, i2))
        out.append((l1_members[i1], l2_members[i2], float(abs(c1[i1] - c2[i2]))))
    return out


# ---------------------------------------------------------------------------
# Step 2: torus correction (expensive -- checkpointed).
# ---------------------------------------------------------------------------


def _floquet_pair(system: cr3bp.CR3BPSystem, orbit: SymmetricNRHO) -> tuple[complex, complex]:
    state0 = np.array([orbit.x0, 0.0, orbit.z0, 0.0, orbit.ydot0, 0.0])
    mono = monodromy(system, state0, orbit.T_TU)
    eigs = floquet_multipliers(mono)
    cands = [e for e in eigs if abs(e - 1.0) > 1e-3 and abs(e.imag) > 1e-4]
    if len(cands) < 2:
        raise ValueError("fewer than 2 candidate Neimark-Sacker eigenvalues found")
    return (cands[0], cands[1])


def _checkpoint_path(system_key: str, side: str, c_index: int) -> pathlib.Path:
    return _CHECKPOINT_DIR / f"{system_key}_{side}_c{c_index}.pkl"


def _load_or_correct_torus(
    system: cr3bp.CR3BPSystem,
    orbit: SymmetricNRHO,
    seed: SeedSpec,
    *,
    system_key: str,
    side: str,
    c_index: int,
    budget: _CorrectionBudget,
) -> QPTorus | None:
    """Load a cached torus correction, or spend one unit of ``budget`` to
    compute + cache it. Returns ``None`` if the budget is exhausted (caller
    must re-invoke the script to resume)."""
    path = _checkpoint_path(system_key, side, c_index)
    if path.exists():
        with path.open("rb") as fh:
            cached: QPTorus = pickle.load(fh)
            return cached
    if not budget.spend():
        return None
    pair = _floquet_pair(system, orbit)
    t0 = time.time()
    torus = correct_qp_torus(
        system,
        np.array([orbit.x0, 0.0, orbit.z0, 0.0, orbit.ydot0, 0.0]),
        orbit.T_TU,
        pair,
        k=seed.k,
        n_trans=8,
        initial_torus_amplitude=5e-4,
        tol=1e-5,
    )
    print(
        f"  [{system_key}/{side}/c{c_index}] torus corrected in {time.time() - t0:.1f}s, "
        f"invariance_residual={torus.invariance_residual:.2e}"
    )
    _CHECKPOINT_DIR.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as fh:
        pickle.dump(torus, fh)
    return torus


class _CorrectionBudget:
    """Caps how many NEW (uncached) torus corrections one invocation performs
    -- each costs ~60-90s, so this is how the script stays inside a bounded
    per-invocation wall-clock window across repeated resumed invocations.
    ``None`` means unlimited (a from-scratch, single full run)."""

    def __init__(self, max_new: int | None) -> None:
        self.remaining = max_new

    def spend(self) -> bool:
        if self.remaining is None:
            return True
        if self.remaining <= 0:
            return False
        self.remaining -= 1
        return True


# ---------------------------------------------------------------------------
# Step 3: manifold-grid construction with a longitudinal phase offset.
# ---------------------------------------------------------------------------


def _compute_custom_row(
    args: tuple[int, float, NDArray[np.float64], QPTorus, str, float, float, float, float, float],
) -> tuple[int, NDArray[np.float64], NDArray[np.float64], NDArray[np.bool_]]:
    (i, theta_long, thetas_lat, torus, branch, eps, sign, t_max, surface_x, mu) = args
    n_lat = len(thetas_lat)
    row_origins = np.zeros((n_lat, 2))
    row_endpoints = np.full((n_lat, 6), np.nan)
    row_hyperbolic = np.zeros(n_lat, dtype=bool)

    prev_vec = None
    for j, theta_trans in enumerate(thetas_lat):
        row_origins[j, 0] = theta_long
        row_origins[j, 1] = theta_trans
        state, stm = torus_point_stm(torus, float(theta_long), float(theta_trans))
        stab = local_stability(
            state,
            stm,
            hyperbolicity_tol=1e-4,
            prev_vec_u=prev_vec if branch == "unstable" else None,
            prev_vec_s=prev_vec if branch == "stable" else None,
        )
        vec = stab.vec_u if branch == "unstable" else stab.vec_s
        if vec is None:
            continue
        if prev_vec is None and vec[0] * sign < 0.0:
            vec = -vec
        row_hyperbolic[j] = True
        prev_vec = vec
        perturbed = state + eps * vec
        direction = 1.0 if branch == "unstable" else -1.0
        t_span = (0.0, math.copysign(t_max, direction))

        def x_event(t: float, y: NDArray[np.float64], mu: float) -> float:
            return float(y[0] - surface_x)

        x_event.terminal = False  # type: ignore[attr-defined]
        x_event.direction = 0.0  # type: ignore[attr-defined]
        sol = scipy.integrate.solve_ivp(  # type: ignore[call-overload]
            cr3bp.cr3bp_eom,
            t_span,
            perturbed,
            args=(mu,),
            method="DOP853",
            rtol=1e-12,
            atol=1e-12,
            events=x_event,
            max_step=abs(t_max) / 20.0,
        )
        crossing = None
        if sol.t_events is not None and len(sol.t_events[0]) > 0:
            for y_event in sol.y_events[0]:
                if y_event[4] < 0.0:
                    crossing = np.asarray(y_event, dtype=np.float64)
                    break
        if crossing is not None:
            row_endpoints[j, :] = crossing

    return i, row_origins, row_endpoints, row_hyperbolic


def build_phase_shifted_grid(
    torus: QPTorus, branch: str, sign: float, t_max: float, surface_x: float, theta0_offset: float
) -> ManifoldGrid:
    """The #534/#536 ``build_custom_grid`` (ydot<0 first-crossing filter),
    extended with a longitudinal ``theta0_offset`` (#545's synodic
    phase-offset sweep -- the only extra phase freedom an autonomous CR3BP
    single-system pair admits, see module docstring)."""
    thetas_long = (2.0 * math.pi * np.arange(_N_LONG) / _N_LONG + theta0_offset) % (2.0 * math.pi)
    thetas_lat = 2.0 * math.pi * np.arange(_N_LAT) / _N_LAT

    origins = np.zeros((_N_LONG, _N_LAT, 2), dtype=np.float64)
    endpoints = np.full((_N_LONG, _N_LAT, 6), np.nan, dtype=np.float64)
    hyperbolic = np.zeros((_N_LONG, _N_LAT), dtype=np.bool_)

    tasks = [
        (i, theta_long, thetas_lat, torus, branch, 1e-5, sign, t_max, surface_x, torus.system.mu)
        for i, theta_long in enumerate(thetas_long)
    ]
    with concurrent.futures.ProcessPoolExecutor(max_workers=12) as executor:
        results = list(executor.map(_compute_custom_row, tasks))

    for i, row_origins, row_endpoints, row_hyperbolic in results:
        origins[i] = row_origins
        endpoints[i] = row_endpoints
        hyperbolic[i] = row_hyperbolic

    return ManifoldGrid(origins=origins, endpoints=endpoints, hyperbolic=hyperbolic)


# ---------------------------------------------------------------------------
# Step 4: linking-number scan + deflated-Newton candidate refinement.
# ---------------------------------------------------------------------------


def _refine_candidates(
    stable_grid: ManifoldGrid, unstable_grid: ManifoldGrid, d_lo: float, d_hi: float
) -> list[float]:
    """Enumerate DISTINCT roots of ``closest_curve_distance(D) = 0`` via
    #524's deflated Newton, seeded across the overlap range (not just the
    coarse scan's sign-change midpoints) -- "enumerate, don't hope one seed
    lands in the right basin"."""

    def residual_fn(z: np.ndarray) -> np.ndarray | None:
        d = float(z[0])
        if d < d_lo or d > d_hi:
            return None
        dist = closest_curve_distance(
            stable_grid,
            unstable_grid,
            scanning_component="z",
            curve_components=("y", "ydot", "zdot"),
            d=d,
        )
        if dist is None:
            return None
        return np.array([dist])

    span = d_hi - d_lo
    seeds = [np.array([d_lo + span * f]) for f in np.linspace(0.05, 0.95, 8)]
    roots = enumerate_roots(
        residual_fn,
        seeds,
        tol=1e-6,
        max_iter=30,
        dedup_tol=max(span * 0.01, 1e-8),
        step_cap=span * 0.25,
    )
    return sorted(float(r[0]) for r in roots)


def _run_c_phase_combo(
    stable_torus: QPTorus, unstable_torus: QPTorus, theta0_offset: float, label: str
) -> dict[str, Any]:
    surface_x = 1.0 - stable_torus.system.mu
    unstable_grid = build_phase_shifted_grid(
        unstable_torus, "unstable", 1.0, _T_MAX_UNSTABLE, surface_x, theta0_offset
    )
    stable_grid = build_phase_shifted_grid(
        stable_torus, "stable", -1.0, _T_MAX_STABLE, surface_x, 0.0
    )
    n_cross_u = int(np.sum(np.isfinite(unstable_grid.endpoints[:, :, 0])))
    n_cross_s = int(np.sum(np.isfinite(stable_grid.endpoints[:, :, 0])))

    z_u = unstable_grid.endpoints[:, :, 2]
    z_s = stable_grid.endpoints[:, :, 2]
    z_u_f, z_s_f = z_u[np.isfinite(z_u)], z_s[np.isfinite(z_s)]

    result = {
        "label": label,
        "theta0_offset": theta0_offset,
        "n_cross_unstable": n_cross_u,
        "n_cross_stable": n_cross_s,
        "overlap": None,
        "sign_change_locations": [],
        "refined_candidates": [],
    }
    if len(z_u_f) == 0 or len(z_s_f) == 0:
        print(f"    {label}: no finite crossings (u={n_cross_u}, s={n_cross_s})")
        return result

    overlap_min = max(z_s_f.min(), z_u_f.min())
    overlap_max = min(z_s_f.max(), z_u_f.max())
    if overlap_min >= overlap_max:
        print(f"    {label}: no z-overlap (u={n_cross_u}, s={n_cross_s})")
        return result

    result["overlap"] = [float(overlap_min), float(overlap_max)]
    d_values = np.linspace(overlap_min, overlap_max, _D_SCAN_N)
    scan = scan_linking_number(
        stable_grid,
        unstable_grid,
        scanning_component="z",
        curve_components=("y", "ydot", "zdot"),
        d_values=d_values,
    )
    sign_changes = scan.sign_change_locations()
    result["sign_change_locations"] = sign_changes

    refined = _refine_candidates(stable_grid, unstable_grid, overlap_min, overlap_max)
    result["refined_candidates"] = refined

    print(
        f"    {label}: u={n_cross_u}/{_N_LONG * _N_LAT} s={n_cross_s}/{_N_LONG * _N_LAT} "
        f"overlap=[{overlap_min:.5f},{overlap_max:.5f}] sign_changes={len(sign_changes)} "
        f"refined_roots={len(refined)}"
    )
    return result


# ---------------------------------------------------------------------------
# Orchestration.
# ---------------------------------------------------------------------------


def sweep_system(
    system_key: str,
    system: cr3bp.CR3BPSystem,
    seed_l1: SeedSpec,
    seed_l2: SeedSpec,
    *,
    d_x0: float,
    n_steps_each_dir: int,
    n_c: int,
    budget: _CorrectionBudget,
) -> tuple[list[dict[str, Any]], bool]:
    """Returns ``(per-combo results, all_corrections_ready)``. If
    ``all_corrections_ready`` is False the caller should re-invoke this
    script (budget exhausted mid-way) before trusting the results list."""
    print(f"\n=== {system_key}: building Jacobi band (natural-parameter continuation) ===")
    l1_members = _build_family_members(system, seed_l1, d_x0=d_x0, n_steps=n_steps_each_dir)
    l2_members = _build_family_members(system, seed_l2, d_x0=d_x0, n_steps=n_steps_each_dir)
    band = _match_band(l1_members, l2_members, n_c)
    for i, (m1, m2, dc) in enumerate(band):
        print(
            f"  C-pair {i}: L1 x0={m1.x0:.8f} C={m1.jacobi:.6f}  "
            f"L2 x0={m2.x0:.8f} C={m2.jacobi:.6f}  deltaC={dc:.3e}"
        )

    all_ready = True
    results: list[dict[str, Any]] = []
    for c_index, (l1_orbit, l2_orbit, delta_c) in enumerate(band):
        torus_l1 = _load_or_correct_torus(
            system,
            l1_orbit,
            seed_l1,
            system_key=system_key,
            side="l1",
            c_index=c_index,
            budget=budget,
        )
        torus_l2 = _load_or_correct_torus(
            system,
            l2_orbit,
            seed_l2,
            system_key=system_key,
            side="l2",
            c_index=c_index,
            budget=budget,
        )
        if torus_l1 is None or torus_l2 is None:
            all_ready = False
            print(
                f"  C-pair {c_index}: correction budget exhausted, resume with another invocation"
            )
            continue
        print(
            f"  C-pair {c_index} (deltaC={delta_c:.3e}): "
            f"sweeping {len(_PHASE_OFFSETS)} phase offsets"
        )
        for offset in _PHASE_OFFSETS:
            label = f"{system_key} C-pair{c_index} phase={offset:.4f}"
            combo = _run_c_phase_combo(torus_l2, torus_l1, offset, label)
            combo["system"] = system_key
            combo["c_index"] = c_index
            combo["jacobi_l1"] = l1_orbit.jacobi
            combo["jacobi_l2"] = l2_orbit.jacobi
            combo["delta_c"] = delta_c
            results.append(combo)
    return results, all_ready


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--max-new-corrections",
        type=int,
        default=None,
        help=(
            "Cap on NEW (uncached) torus corrections this invocation performs "
            "(~60-90s each). Omit for an unlimited from-scratch run; pass a small "
            "number (e.g. 1-2) to keep one invocation's wall-clock bounded and "
            "resume via a subsequent invocation (checkpoints persist under "
            f"{_CHECKPOINT_DIR})."
        ),
    )
    args = parser.parse_args()
    budget = _CorrectionBudget(args.max_new_corrections)

    preflight_search(
        task_no=545,
        region_id=_REGION_ID,
        method=_METHOD,
        script_path=pathlib.Path(__file__),
        n_points=3 * len(_PHASE_OFFSETS),  # 3 Jupiter-Europa C-pairs x 4 phase offsets = 12
    )

    print("=" * 60)
    print("POSITIVE CONTROL: Earth-Moon L1/L2 (task #534's own system/seeds)")
    print("=" * 60)
    em_results, em_ready = sweep_system(
        "earth_moon",
        _em_system(),
        EM_SEED_L1,
        EM_SEED_L2,
        d_x0=2e-4,
        n_steps_each_dir=6,
        n_c=3,
        budget=budget,
    )

    print("\n" + "=" * 60)
    print("TARGET SWEEP: Jupiter-Europa L1/L2 Jacobi band + phase offsets")
    print("=" * 60)
    je_results, je_ready = sweep_system(
        "jupiter_europa",
        _je_system(),
        JE_SEED_L1,
        JE_SEED_L2,
        d_x0=1e-4,
        n_steps_each_dir=6,
        n_c=3,
        budget=budget,
    )

    report = {
        "earth_moon_positive_control": em_results,
        "earth_moon_all_corrections_ready": em_ready,
        "jupiter_europa": je_results,
        "jupiter_europa_all_corrections_ready": je_ready,
    }
    _REPORT_PATH.parent.mkdir(parents=True, exist_ok=True)
    with _REPORT_PATH.open("w", encoding="utf-8") as fh:
        json.dump(report, fh, indent=2)
    print(f"\nWrote {_REPORT_PATH}")

    if not (em_ready and je_ready):
        print(
            "\nNOT ALL torus corrections are cached yet -- re-invoke this script "
            "(same or smaller --max-new-corrections) to resume; this run's "
            "results above are PARTIAL."
        )
        return

    n_refined_total = sum(len(c["refined_candidates"]) for c in em_results + je_results)
    print(f"\nTotal refined-candidate connections across both systems: {n_refined_total}")


if __name__ == "__main__":
    main()
