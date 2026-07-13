"""#582: mandatory downstream pipeline for ``isolated_3d_asymmetric_fitness``
GA survivors.

Per #582's spec (``data/OUTSTANDING.md``) a GA cluster is a BASIN INDICATOR,
NOT a converged orbit -- before any "found an orbit" claim a survivor must be
routed, IN ORDER, through:

  1. :func:`refine_ga_candidate` -- the EXISTING asymmetric 3D corrector
     (``search/cr3bp_general_periodic_3d.py::correct_general_periodic_3d``,
     #291) with independent closure verification. The GA's new capability is
     asymmetric/3D SEED GENERATION, not asymmetric correction; this module
     does NOT add a new corrector.
  2. :func:`classify_symmetry` -- an explicit mirror-image / perpendicular-
     crossing test. A converged point that happens to sit ON a known
     SYMMETRIC orbit is trivially periodic too; without this test a symmetric
     orbit could be misreported as a novel asymmetric one.
  3. :func:`build_candidate_signature` -- populates a
     ``literature_check.CandidateSignature`` so the existing 3D literature
     matcher can engage (it only fires on anchors sharing the same
     ``primary``; see the function's own docstring for how the mu=0.001
     system reaches its one matching anchor).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.search.binary_star_search import topology_3d as compute_topology_3d
from cyclerfinder.search.cr3bp_general_periodic_3d import (
    FREE_VARS_FULL_ASYMMETRIC,
    RESIDUAL_FULL_STATE_AT_T,
    Periodic3DOrbit,
    correct_general_periodic_3d,
)
from cyclerfinder.search.literature_check import CandidateSignature

# ---------------------------------------------------------------------------
# Step 1: refine through the EXISTING asymmetric 3D corrector.
# ---------------------------------------------------------------------------


def refine_ga_candidate(
    system: cr3bp.CR3BPSystem,
    genome: NDArray[np.float64],
    *,
    tol: float = 1e-10,
    independent_tol: float = 1e-6,
    max_iter: int = 60,
    require_monotone_decrease: bool = True,
) -> Periodic3DOrbit:
    """Refine one ``isolated_3d_asymmetric_fitness`` genome via the existing corrector.

    ``genome`` is ``(x0, z0, xdot0, ydot0, zdot0, T)`` (same layout the GA
    uses). Routes through
    :func:`~cyclerfinder.search.cr3bp_general_periodic_3d.correct_general_periodic_3d`
    with the full-asymmetric free-var/residual bundle -- NOT a new corrector.

    ``require_monotone_decrease`` defaults True per that function's OWN
    guidance ("Set True for ill-conditioned long-arc closures (asymmetric
    full-6D-at-T mode)"): a GA seed already close to periodic (small blind
    residual) can still overshoot into an unrelated basin under the default
    blind-Newton step in this 6-residual/7-unknown under-determined mode --
    confirmed empirically while building #582's positive control (blind
    Newton: residual 0.30/independent 2.85, landed on a totally different
    point; monotone-decrease: residual 1.3e-12/independent 1.2e-12, landed
    near the known #440 member). Callers with a different failure profile can
    still pass ``False`` explicitly.
    """
    g = np.asarray(genome, dtype=np.float64)
    x0, z0, xdot0, ydot0, zdot0, t = (float(v) for v in g)
    state0_guess = np.array([x0, 0.0, z0, xdot0, ydot0, zdot0])
    return correct_general_periodic_3d(
        system,
        state0_guess,
        t,
        free_vars=FREE_VARS_FULL_ASYMMETRIC,
        residual_indices=RESIDUAL_FULL_STATE_AT_T,
        is_half_period_residual=False,
        tol=tol,
        independent_tol=independent_tol,
        max_iter=max_iter,
        require_monotone_decrease=require_monotone_decrease,
    )


# ---------------------------------------------------------------------------
# Step 2: explicit symmetry classification.
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class SymmetryClassification:
    """Verdict: does a converged orbit ALSO sit on a symmetric periodic orbit?

    The CR3BP has the standard mirror symmetry
    ``(x,y,z,xdot,ydot,zdot,t) -> (x,-y,z,-xdot,ydot,-zdot,-t)``. A periodic
    orbit is invariant under this map iff it crosses the x-z plane (``y=0``)
    PERPENDICULARLY (``xdot=0`` AND ``zdot=0`` at the crossing) -- exactly the
    IC convention ``correct_symmetric_fixed_jacobi`` /
    ``correct_symmetric_nrho`` use. A candidate converged by the general
    (unconstrained) asymmetric corrector can land on such an orbit purely by
    coincidence, parameterized from an arbitrary phase where the perpendicular
    condition is not obviously visible at ``t=0``. This function densely
    samples the FULL period looking for ANY such crossing, so a symmetric
    orbit is never misreported as novel-asymmetric just because its
    corrector-returned IC phase does not itself look symmetric.
    """

    is_symmetric: bool
    best_crossing_residual: float
    """``sqrt(xdot^2 + zdot^2)`` (nondim speed) at the most-perpendicular y=0
    crossing found; ``inf`` if no crossing exists (propagation never returns
    to ``y=0``, e.g. a genuinely 3D orbit confined to one side)."""
    best_crossing_time: float | None
    n_crossings_checked: int
    tol: float


def classify_symmetry(
    system: cr3bp.CR3BPSystem,
    state0: NDArray[np.float64],
    period: float,
    *,
    tol: float = 1e-6,
    rtol: float = 1e-12,
    atol: float = 1e-12,
    n_samples: int = 4000,
) -> SymmetryClassification:
    """Perpendicular x-z-plane-crossing test -- the real symmetry classifier.

    Propagates one full period, finds every ``y=0`` sign-change (bracketed),
    root-refines each crossing time with Brent's method on ``y(t)``, and
    evaluates ``(xdot, zdot)`` there via dense output. If ANY crossing has
    ``sqrt(xdot^2 + zdot^2) < tol`` the orbit is classified symmetric: it
    sits on the mirror-invariant manifold regardless of which phase the
    general corrector happened to converge at. ``tol`` default 1e-6 nondim
    speed is the same order as ``correct_general_periodic_3d``'s own
    ``independent_tol`` floor.
    """
    sol = solve_ivp(
        cr3bp.cr3bp_eom,
        (0.0, period),
        np.asarray(state0, dtype=np.float64),
        args=(system.mu,),
        method="DOP853",
        rtol=rtol,
        atol=atol,
        max_step=max(period / n_samples, 1e-9),
        dense_output=True,
    )
    if not sol.success or sol.sol is None:
        return SymmetryClassification(
            is_symmetric=False,
            best_crossing_residual=float("inf"),
            best_crossing_time=None,
            n_crossings_checked=0,
            tol=tol,
        )
    dense = sol.sol

    y = sol.y[1]
    sign = np.sign(y)
    crossing_idx = np.where(np.diff(sign) != 0.0)[0]

    def y_of(t: float) -> float:
        return float(dense(t)[1])

    best_res = float("inf")
    best_t: float | None = None
    n_checked = 0
    for i in crossing_idx:
        t_lo, t_hi = float(sol.t[i]), float(sol.t[i + 1])
        y_lo, y_hi = y_of(t_lo), y_of(t_hi)
        if y_lo == 0.0 or y_hi == 0.0 or (y_lo > 0.0) == (y_hi > 0.0):
            continue
        try:
            t_cross = brentq(y_of, t_lo, t_hi, xtol=1e-13)
        except ValueError:
            continue
        state_c = dense(t_cross)
        xdot_c, zdot_c = float(state_c[3]), float(state_c[5])
        res = math.hypot(xdot_c, zdot_c)
        n_checked += 1
        if res < best_res:
            best_res = res
            best_t = t_cross

    return SymmetryClassification(
        is_symmetric=best_res < tol,
        best_crossing_residual=best_res,
        best_crossing_time=best_t,
        n_crossings_checked=n_checked,
        tol=tol,
    )


# ---------------------------------------------------------------------------
# Step 3: populate the literature-matcher's CandidateSignature.
# ---------------------------------------------------------------------------


def build_candidate_signature(
    system: cr3bp.CR3BPSystem,
    orbit: Periodic3DOrbit,
    *,
    p: int,
    q: int,
    rtol: float = 1e-11,
    atol: float = 1e-11,
) -> CandidateSignature:
    """Populate the literature-matcher's structural fingerprint for a survivor.

    ``primary="Earth"``, ``sequence=("Moon",)`` is a MATCHER-ENGAGEMENT LABEL,
    not a physical claim about this candidate. The only spatial-CR3BP anchor
    at mu=0.001 (the #579-corrected Antoniadou & Libert 2019 record in
    ``genome/known_corpus_3d.py``) is itself pinned under
    ``primary="Earth"``, ``body_set={"Moon"}`` even though its own docstring
    explicitly annotates it ``mu=0.001`` (a generic Jupiter-mass-planet MMR
    system, NOT physical Earth-Moon mu=0.0121506683) -- the corpus schema
    has no separate "generic system" primary slot, so the anchor reuses the
    Earth/Moon labels purely to have a matchable ``(primary, body_set)`` key.
    ``literature_check._candidate_anchors``' own filter is a literal
    ``anchor.primary != sig.primary`` check (see its docstring: "it only
    fires on anchors sharing the same primary"), so a signature for THIS
    mu=0.001 system must reuse the same label to reach the anchor that
    actually covers it; anything else would silently never match.
    """
    topo = compute_topology_3d(system.mu, orbit.state0, orbit.T_TU, rtol=rtol, atol=atol)
    return CandidateSignature(
        primary="Earth",
        sequence=("Moon",),
        resonances=(f"{p}:{q}",),
        topology_label=frozenset({"resonant"}),
        topology_3d={
            "k1": topo.k1,
            "k2": topo.k2,
            "k_z": topo.k_z,
            "jacobi": float(orbit.jacobi),
        },
    )


def literature_anchors_engaged(sig: CandidateSignature) -> list[str]:
    """Names of ``literature_check`` corpus anchors this signature would search.

    Proves the matcher ACTUALLY ENGAGES (non-empty return) for a given
    ``CandidateSignature`` -- the check #582 requires before any novelty
    claim. This is NOT the literature check itself (``check_literature``
    needs a live, injected web-search callable); it verifies the structural
    fingerprint reaches the right anchor pool, which is the load-bearing
    prerequisite a wrongly-labeled signature would silently fail.
    """
    from cyclerfinder.search.literature_check import _candidate_anchors

    return [a.name for a in _candidate_anchors(sig)]
