"""Pickle-safe per-cell closure for the #284 asymmetric-corrector scan (#343).

This module hosts the top-level closure that :func:`parallel_sweep` dispatches
across workers. Keeping it at module level (not nested in the driver script)
satisfies the loky pickle-safety contract per #321 -- lambdas / nested closures
fail at submission with :class:`PicklingError`.

The closure rebuilds heavyweight state (system, flagger) inside each worker
lazily, caching it in a process-global dict. This avoids re-importing /
re-fitting per cell while still pickling cleanly: only the cell tuple and the
top-level function reference cross the IPC boundary.

The full scan pipeline (corrector -> closure gate -> topology gate ->
independent DOP853 cross-check -> known-family check -> ML flagger ->
prioritizer score -> literature-corpus footprint) lives here. The driver script
just builds the grid, calls :func:`parallel_sweep`, and writes the JSONL.

NO catalogue writeback. NO novelty claims. The closure emits a row per
CONVERGED + dedup-survivor cell; the driver applies the final aggregation.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from datetime import UTC, datetime
from typing import Any

import numpy as np
from scipy.integrate import solve_ivp

import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.ml.falsepos_flagger import FalsePosFlagger
from cyclerfinder.ml.falsepos_labels import build_training_set
from cyclerfinder.search.binary_star_search import winding_topology
from cyclerfinder.search.cr3bp_general_periodic import (
    GeneralPeriodicOrbit,
    _ydot0_general,
    correct_general_periodic,
)
from cyclerfinder.search.literature_check import (
    CandidateSignature,
    _candidate_anchors,
)
from cyclerfinder.search.reachable_representatives import TU_DAYS, braik_ross_system
from cyclerfinder.search.single_orbit_prioritizer import score_single_orbit

# Per-process lazy cache. Each worker fork sees its own copy after fork() so
# this is a cheap-once-per-worker fit, not a per-cell fit.
_WORKER_STATE: dict[str, Any] = {}

# Scan-wide constants (mirror #284 phase 1 exactly so the comparison is honest).
PERIOD_GUESS_TU = 14.0
CORRECTOR_TOL = 1e-11
CORRECTOR_MAX_ITER = 20
INDEPENDENT_CLOSURE_TOL = 1e-6


@dataclass(frozen=True)
class AsymCell:
    """One grid cell -- pickle-safe (frozen, primitives only)."""

    k1: int
    k2: int
    jacobi_request: float
    seed_x0: float
    seed_xdot0: float
    ydot0_sign: float
    half_crossings: int


def _get_system() -> cr3bp.CR3BPSystem:
    sys_ = _WORKER_STATE.get("system")
    if sys_ is None:
        sys_ = braik_ross_system()
        _WORKER_STATE["system"] = sys_
    return sys_


def _get_flagger() -> FalsePosFlagger | None:
    if "flagger_fitted" not in _WORKER_STATE:
        try:
            clf = FalsePosFlagger()
            x_train, y_train, _ = build_training_set()
            clf.fit(x_train, y_train)
            _WORKER_STATE["flagger"] = clf
        except Exception:
            _WORKER_STATE["flagger"] = None
        _WORKER_STATE["flagger_fitted"] = True
    return _WORKER_STATE["flagger"]  # type: ignore[no-any-return]


def _independent_closure(mu: float, orbit: GeneralPeriodicOrbit) -> tuple[float, bool]:
    """DOP853 rtol=atol=1e-12 re-propagation closure (independent cross-check)."""
    if not math.isfinite(orbit.ydot0):
        return float("inf"), False
    state0 = np.array([orbit.x0, 0.0, 0.0, orbit.xdot0, orbit.ydot0, 0.0])
    sol = solve_ivp(
        cr3bp.cr3bp_eom,
        (0.0, orbit.period),
        state0,
        args=(mu,),
        method="DOP853",
        rtol=1e-12,
        atol=1e-12,
    )
    if not sol.success:
        return float("inf"), False
    return float(np.linalg.norm(sol.y[:, -1] - state0)), True


def _signature_for_em_cycler(k1: int, k2: int) -> CandidateSignature:
    """Structural fingerprint identical to #284 (Earth-Moon CR3BP cycler)."""
    return CandidateSignature(
        primary="Earth",
        sequence=("Moon",),
        period_k=None,
        period_years=None,
        vinf_per_encounter_kms=(),
        resonances=(f"{k1}:{k2}",),
        n_rev=(k1, k2),
    )


def _offline_em_corpus_match(sig: CandidateSignature) -> list[dict[str, str]]:
    """Offline corpus anchors (necessary-not-sufficient, same as #284)."""
    anchors = _candidate_anchors(sig)
    return [{"name": a.name, "citation": a.citation, "doi": a.doi or ""} for a in anchors]


def _fp_flagger_score(clf: FalsePosFlagger | None, orbit: GeneralPeriodicOrbit) -> float:
    """Identical wiring to #284 -- ML flagger fed CR3BP-analogue fields."""
    if clf is None:
        return float("nan")
    rec: dict[str, Any] = {
        "max_residual_kms": float(orbit.closure_residual),
        "bend_feasible": True,
        "topology_match": True,
        "period_days": float(orbit.period * TU_DAYS),
        "closure_method_version": "23b980e",
        "closure_date": "2026-06-17",
        "model_assumption": "cr3bp_asymmetric_general_periodic",
        "cross_check_shared_with_primary": False,
    }
    try:
        return float(clf.score(rec))
    except Exception:
        return float("nan")


def _is_known_em_family(k1: int, k2: int, period_days: float) -> bool:
    """Same sourced (k1, k2, period) band as #284."""
    em_known: list[tuple[int, int, float]] = [
        (1, 1, 42.14),  # C11a
        (1, 1, 55.96),  # C11b
        (2, 1, 84.53),  # C21
        (3, 2, 78.61),  # C32
    ]
    for kk1, kk2, p in em_known:
        if k1 == kk1 and k2 == kk2 and abs(period_days - p) / p < 0.02:
            return True
    return False


def process_cell(cell: AsymCell) -> dict[str, Any] | None:
    """Run one cell through the full guard chain.

    Returns ``None`` for cells that fail the corrector or the early gates; the
    driver aggregates over non-``None`` results. Pickle-safe: the cell is a
    frozen dataclass over primitives, the return is a plain dict.
    """
    system = _get_system()
    mu = system.mu
    flagger = _get_flagger()

    # Jacobi-feasibility short-circuit.
    try:
        _ydot0_general(
            cell.seed_x0,
            cell.seed_xdot0,
            cell.jacobi_request,
            mu,
            cell.ydot0_sign,
        )
    except Exception:
        return None

    # Corrector.
    try:
        orbit = correct_general_periodic(
            system,
            cell.seed_x0,
            cell.seed_xdot0,
            cell.jacobi_request,
            PERIOD_GUESS_TU,
            half_crossings=cell.half_crossings,
            ydot0_sign=cell.ydot0_sign,
            tol=CORRECTOR_TOL,
            max_iter=CORRECTOR_MAX_ITER,
            rtol=1e-12,
            atol=1e-12,
        )
    except (ValueError, RuntimeError):
        return None
    except Exception:
        return None

    if not orbit.converged:
        return None
    if not math.isfinite(orbit.closure_residual):
        return None
    if orbit.closure_residual > INDEPENDENT_CLOSURE_TOL:
        return None

    # Topology classification.
    state0 = np.array([orbit.x0, 0.0, 0.0, orbit.xdot0, orbit.ydot0, 0.0])
    try:
        topo = winding_topology(mu, state0, orbit.period)
    except Exception:
        return None

    period_days = orbit.period * TU_DAYS
    topology_match = topo.k1 == cell.k1 and topo.k2 == cell.k2

    # Independent DOP853 closure cross-check.
    ind_closure, ind_ok = _independent_closure(mu, orbit)
    independent_pass = ind_ok and ind_closure < INDEPENDENT_CLOSURE_TOL

    # Literature signature (offline; necessary-not-sufficient).
    sig = _signature_for_em_cycler(topo.k1, topo.k2)
    anchors = _offline_em_corpus_match(sig)
    known_family = _is_known_em_family(topo.k1, topo.k2, period_days)
    literature_fresh = (not known_family) and independent_pass

    # ML false-positive flagger.
    p_fp = _fp_flagger_score(flagger, orbit)
    ml_low_fp = math.isfinite(p_fp) and p_fp < 0.5

    # Prioritizer score (#310 single-orbit adapter, default FiveTierPrioritizer
    # = Tier 0 only). The adapter returns NaN for skipped tiers; we surface the
    # Tier 0 ΔV-to-surrogate (lower-is-better) as the rank signal.
    prioritizer_tier0: float | None = None
    prioritizer_neighbor: str | None = None
    prioritizer_notes = ""
    try:
        score = score_single_orbit(
            state0,
            float(orbit.period),
            system,
            candidate_id=f"k{cell.k1}{cell.k2}_C{orbit.jacobi:.4f}_T{period_days:.2f}d",
            candidate_jacobi=float(orbit.jacobi),
            use_surrogate_neighbor=True,
        )
        prioritizer_tier0 = score.tier_0_score
        prioritizer_neighbor = score.surrogate_pair_neighbor_id
        prioritizer_notes = score.notes
    except Exception as exc:
        prioritizer_notes = f"prioritizer raised: {exc!r}"

    novelty_claimable = topology_match and independent_pass and literature_fresh and ml_low_fp

    return {
        "k_target": [cell.k1, cell.k2],
        "k_classified": [topo.k1, topo.k2],
        "topology_match": bool(topology_match),
        "c_request": round(cell.jacobi_request, 6),
        "seed_x0": round(cell.seed_x0, 6),
        "seed_xdot0": round(cell.seed_xdot0, 6),
        "ydot0_sign": cell.ydot0_sign,
        "half_crossings": cell.half_crossings,
        "x0": round(orbit.x0, 9),
        "xdot0": round(orbit.xdot0, 9),
        "ydot0": round(orbit.ydot0, 9),
        "jacobi": round(orbit.jacobi, 10),
        "period_TU": round(orbit.period, 9),
        "period_days": round(period_days, 6),
        "asymmetry": round(orbit.asymmetry, 9),
        "x_min": round(topo.x_min, 6),
        "x_max": round(topo.x_max, 6),
        "prograde": bool(topo.prograde),
        "reaches_secondary": bool(topo.reaches_secondary),
        "corrector_residual": float(orbit.residual),
        "corrector_closure_radau": float(orbit.closure_residual),
        "independent_closure_dop853": (float(ind_closure) if ind_ok else None),
        "independent_closure_pass": bool(independent_pass),
        "n_newton_iter": int(orbit.n_iter),
        "literature_offline_anchors": anchors,
        "known_em_family_collision": bool(known_family),
        "literature_fresh_offline": bool(literature_fresh),
        "p_false_positive": (round(p_fp, 4) if math.isfinite(p_fp) else None),
        "ml_low_fp": bool(ml_low_fp),
        "prioritizer_tier0_dv_kms": (
            round(prioritizer_tier0, 6) if prioritizer_tier0 is not None else None
        ),
        "prioritizer_neighbor_id": prioritizer_neighbor,
        "prioritizer_notes": prioritizer_notes,
        "novelty_claimable": bool(novelty_claimable),
        "scan_id": "343",
        "scan_method_version": "asymmetric_corrector(23b980e)+prioritizer(310)+parallel(321)",
        "scan_timestamp_utc": datetime.now(UTC).isoformat(),
    }


__all__ = [
    "INDEPENDENT_CLOSURE_TOL",
    "AsymCell",
    "process_cell",
]
