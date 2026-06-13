"""One-DSM-per-leg genome (the Takao eta-coordinate) — multi-arc-return primitive.

This is the architectural sibling of :mod:`cyclerfinder.search.free_return` (the
#137 radial-crossing genome). Where the free-return module models an
Earth->Mars->Earth transfer as a *single* heliocentric ellipse (one shape ``(a, e)``
shared by both legs), this module supplies the primitive the single-ellipse genome
provably lacks: an *interior impulse* on each leg. That interior impulse is what
lets a leg follow a different ballistic arc on its front fraction than on its back
fraction — i.e. a leg that is NOT a piece of one repeating ellipse. The S1L1 /
Jones multi-arc closure blocker (``docs/notes/multi-arc-classification.md`` §7/§12;
the #137 / MBH Gate-3 negative result) is exactly the case where the sourced
geometry is *two generic-return arcs*, not one ellipse, so a single-ellipse genome
has no sourced-anchor basin. The one-DSM-per-leg genome can represent that.

The transcription
-----------------
Y. Takao, "Mission Analysis for the First-Ever Saturn Trojan 2019 UO14,"
arXiv:2501.06586 (astro-ph.EP), 2025. Per-leg genome
``Y_i = [V_inf, alpha, beta, tau, eta]_i`` (Eqs.1-2). The DSM leg evaluator
(Eqs.6-7, transcribed in ``docs/notes/2026-06-07-takao-2025-mpga-1dsm-mining.md``):

* propagate the outgoing state ballistically (heliocentric 2-body) for ``eta*tau``,
  giving the DSM position ``r_12`` and velocity ``v_12`` (the velocity arriving at
  the DSM, BEFORE the impulse);
* solve Lambert from the DSM position to the next body over the remaining
  ``(1-eta)*tau``, giving the post-impulse departure ``v_21`` and the arrival
  ``v_22``;
* the DSM impulse magnitude is the velocity mismatch at the DSM point,
  ``dV_DSM = ||v_21 - v_12||`` (Eq.6);
* the incoming hyperbolic velocity at the next body is ``v_22 - v_planet`` (Eq.7).

``eta in [0, 1]`` is the only genuinely new genome coordinate vs the free-return /
ballistic correctors; everything else (2-body propagate, Lambert) is reused from
:mod:`cyclerfinder.core`.

Conventions mirrored from ``free_return.py``
--------------------------------------------
* A frozen dataclass result carrying both the CONSTRAINED quantity (the residual /
  total dV the corrector drives) and the FREE / EMERGED evidence (per-leg dV
  breakdown, per-body emerged V_inf). The sourced V_inf is NEVER imposed; it emerges
  from the converged genome and is comparable as evidence (the golden-rule
  separation).
* ``converged`` decided by residual MAGNITUDE alone (``< tol_kms``), by design; the
  ``least_squares`` ``success`` flag is kept only for the audit trail.
* A full audit trail (per-leg dV, DSM states, solver nfev) on the result.
* RNG-free / deterministic: no randomness anywhere in this module (MBH supplies the
  hops; this module only evaluates and corrects).

Pure: depends only on core/constants, core/kepler (propagate), core/lambert.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass, field, replace
from typing import Any

import numpy as np
from numpy.typing import NDArray
from scipy.optimize import least_squares

from cyclerfinder.core.constants import (
    AU_KM,
    MU_SUN_KM3_S2,
    PLANETS,
    SAFE_PERIHELION_KM,
    SECONDS_PER_DAY,
)
from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.core.fbs_match_point import (
    FbsLeg,
    match_point_defect,
    match_point_defect_jacobian,
)
from cyclerfinder.core.flyby import flyby_dv
from cyclerfinder.core.kepler import KeplerError, propagate
from cyclerfinder.core.lambert import LambertError, lambert
from cyclerfinder.search.mbh import MBHStep

Vec3 = NDArray[np.float64]  # shape (3,), dtype float64
DAY_S = SECONDS_PER_DAY


# ---------------------------------------------------------------------------
# The leg primitive (Takao Eq.6-7)
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class DsmLegResult:
    """Outcome of a single propagate-then-Lambert (one-DSM) leg.

    Attributes
    ----------
    v_arrive:
        Heliocentric velocity at the target endpoint (``v_22`` in Takao Eq.6),
        km/s. The post-flyby caller forms the incoming V_inf as
        ``v_arrive - v_planet_target`` (Eq.7).
    v_depart_post_dsm:
        Heliocentric departure velocity from the DSM point AFTER the impulse
        (``v_21``), km/s. The Lambert start velocity.
    v_arrive_pre_dsm:
        Heliocentric velocity arriving AT the DSM point BEFORE the impulse
        (``v_12``), km/s. The ballistic propagation endpoint velocity.
    dv_dsm_kms:
        The DSM impulse magnitude ``||v_21 - v_12||`` (Eq.6), km/s.
    r_dsm:
        Heliocentric position of the DSM point (``r_12 = r_21``), km. Audit/viz.
    t_dsm_sec:
        Time of the DSM along the leg, ``eta * tof`` seconds (relative to leg
        start). Audit/viz.
    n_revs_chosen:
        The revolution count of the back-arc Lambert branch actually selected,
        ``0`` for the single-revolution (direct) branch. When ``max_revs == 0``
        this is always ``0`` (the legacy single-rev path); when ``max_revs > 0``
        it records which multi-rev branch minimised :attr:`dv_dsm_kms`. Audit.
    branch_chosen:
        The branch label of the selected back-arc Lambert solution
        (``"single"`` / ``"low"`` / ``"high"``). Audit; pairs with
        :attr:`n_revs_chosen`.
    """

    v_arrive: Vec3
    v_depart_post_dsm: Vec3
    v_arrive_pre_dsm: Vec3
    dv_dsm_kms: float
    r_dsm: Vec3
    t_dsm_sec: float
    n_revs_chosen: int = 0
    branch_chosen: str = "single"


# eta is clamped this far from the singular endpoints {0, 1}: at exactly eta=0 the
# ballistic front arc has zero duration (DSM == leg start) and at exactly eta=1 the
# Lambert back arc has zero duration (singular t->0 Lambert). The degeneracy gate
# probes eta APPROACHING these endpoints (e.g. 1e-4), which is well-posed; this
# guard only rejects the exact-endpoint singular calls.
_ETA_EPS: float = 1.0e-9


def dsm_leg(
    r0: Vec3,
    v0: Vec3,
    tof: float,
    eta: float,
    target_r: Vec3,
    *,
    mu: float = MU_SUN_KM3_S2,
    prograde: bool = True,
    max_revs: int = 0,
    rev_branch: tuple[int, str] | None = None,
) -> DsmLegResult:
    """One interior-impulse (DSM) leg: propagate ``eta*tof`` then Lambert the rest.

    Takao Eq.6-7. Starting from the outgoing heliocentric state ``(r0, v0)`` at the
    leg's first body, the trajectory is propagated ballistically (heliocentric
    2-body) for ``eta*tof`` to the DSM point; Lambert is then solved from the DSM
    position to ``target_r`` over the remaining ``(1-eta)*tof``. The DSM impulse is
    the velocity mismatch at the DSM point.

    Parameters
    ----------
    r0, v0:
        Heliocentric inertial state at the leg start, ``(3,)`` float64, km and km/s.
    tof:
        Total leg time of flight, seconds. Must be strictly positive.
    eta:
        DSM timing fraction in ``[0, 1]`` (Takao). ``eta`` of the ToF is flown
        ballistically before the impulse; ``(1-eta)`` is the Lambert arc. Must be
        strictly inside ``(0, 1)`` up to :data:`_ETA_EPS` (the exact endpoints are
        singular -- the degeneracy gate probes eta APPROACHING them).
    target_r:
        Heliocentric position of the target body at leg arrival, ``(3,)`` km.
    mu:
        Central-body gravitational parameter, km^3/s^2 (heliocentric default).
    prograde:
        Lambert transfer sense for the back arc (default prograde / short-way).
    max_revs:
        Maximum number of full revolutions to consider for the back-arc Lambert
        (passed straight through to :func:`cyclerfinder.core.lambert.lambert`).
        The default ``0`` keeps the historical single-revolution-only behaviour
        (bit-identical results). When ``max_revs > 0`` the solver enumerates the
        single-rev branch plus every feasible multi-rev ``low``/``high`` branch
        up to ``max_revs`` and selects the one that MINIMISES the DSM impulse
        ``dV_DSM`` (the leg objective). The resonant loop arcs of a multi-arc
        cycler are >1-period transfers, for which the single-rev branch is forced
        onto a degenerate near-radial high-energy solution; allowing multi-rev
        branches is what makes those arcs representable (the #153 diagnosis).
    rev_branch:
        Optional ``(n_revs, branch)`` selector. When given, the back arc uses
        exactly that branch (e.g. ``(2, "low")``) instead of minimising over the
        enumerated set; ``max_revs`` is widened to ``n_revs`` if needed so the
        requested branch is produced. A :class:`LambertError` is raised if the
        requested branch is infeasible for this geometry/ToF. ``None`` (default)
        selects the dV-minimising branch as described above.

    Returns
    -------
    DsmLegResult
        Arrival velocity, post/pre-DSM velocities, the DSM impulse magnitude, the
        DSM state, and the chosen back-arc revolution/branch for the audit trail.

    Raises
    ------
    ValueError
        On non-positive ``tof`` or ``eta`` at/outside the singular endpoints, or a
        negative ``max_revs``.
    LambertError
        If the back-arc Lambert solve fails (degenerate geometry / no
        convergence), or the explicitly requested ``rev_branch`` is infeasible.
    """
    if tof <= 0.0:
        raise ValueError(f"tof must be positive, got {tof}")
    if not (_ETA_EPS <= eta <= 1.0 - _ETA_EPS):
        raise ValueError(
            f"eta must lie in [{_ETA_EPS}, {1.0 - _ETA_EPS}] (the exact endpoints "
            f"0/1 are singular ballistic/Lambert degeneracies), got {eta}"
        )
    if max_revs < 0:
        raise ValueError(f"max_revs must be non-negative, got {max_revs}")

    r0_arr = np.asarray(r0, dtype=np.float64)
    v0_arr = np.asarray(v0, dtype=np.float64)
    target_arr = np.asarray(target_r, dtype=np.float64)

    t_front = eta * tof
    t_back = (1.0 - eta) * tof

    # Front arc: ballistic 2-body propagation to the DSM point (Eq.6, v_12).
    r_dsm, v12 = propagate(r0_arr, v0_arr, t_front, mu)

    # Back arc: Lambert from the DSM position to the target over (1-eta)*tof.
    # max_revs=0 returns ONLY the single-rev (direct) branch -> the historical
    # path is bit-identical (sols[0] is the single-rev solution). max_revs>0 also
    # returns the feasible multi-rev low/high branches, among which the chosen
    # branch is the one minimising the DSM impulse (the resonant loop arcs need
    # this; the single-rev branch on a >1-period leg is degenerate -- #153).
    enumerate_revs = max_revs
    if rev_branch is not None:
        # Widen the enumeration so the explicitly requested branch is produced.
        enumerate_revs = max(max_revs, int(rev_branch[0]))
    sols = lambert(r_dsm, target_arr, t_back, mu=mu, prograde=prograde, max_revs=enumerate_revs)

    if rev_branch is not None:
        want_n, want_b = int(rev_branch[0]), str(rev_branch[1])
        chosen = next(
            (s for s in sols if s.n_revs == want_n and s.branch == want_b),
            None,
        )
        if chosen is None:
            raise LambertError(
                f"requested rev_branch ({want_n}, {want_b!r}) is infeasible for this "
                f"geometry/ToF (available: "
                f"{[(s.n_revs, s.branch) for s in sols]})"
            )
    elif max_revs == 0:
        # Legacy path: the single-rev (direct) branch, identical to before.
        chosen = sols[0]
    else:
        # Pick the branch minimising the DSM impulse ||v21 - v12||.
        chosen = min(sols, key=lambda s: float(np.linalg.norm(s.v1 - v12)))

    v21 = chosen.v1  # post-impulse departure from the DSM (Eq.6, v_21)
    v22 = chosen.v2  # arrival at the target (Eq.6, v_22)

    dv_dsm = float(np.linalg.norm(v21 - v12))

    return DsmLegResult(
        v_arrive=np.asarray(v22, dtype=np.float64),
        v_depart_post_dsm=np.asarray(v21, dtype=np.float64),
        v_arrive_pre_dsm=np.asarray(v12, dtype=np.float64),
        dv_dsm_kms=dv_dsm,
        r_dsm=np.asarray(r_dsm, dtype=np.float64),
        t_dsm_sec=float(t_front),
        n_revs_chosen=int(chosen.n_revs),
        branch_chosen=str(chosen.branch),
    )


# ---------------------------------------------------------------------------
# OPT-IN: Lambert-free single-leg corrector using the FBS analytic Jacobian
# (Ellison 2018 Path-B; #226). Additive — the Lambert dsm_leg above and the
# chain corrector below are untouched and remain the default everywhere.
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class FbsLegCorrectionResult:
    """Outcome of a Lambert-free FBS match-point single-leg correction (#226).

    Attributes
    ----------
    dv_dsm:
        The converged interior impulse vector ``Δv`` (km/s) at the DSM point. The
        magnitude ``‖Δv‖`` is the same quantity the Lambert :func:`dsm_leg`
        reports as :attr:`DsmLegResult.dv_dsm_kms`.
    v_arrive:
        The converged heliocentric arrival velocity ``vf`` at the target (km/s) —
        the FBS analogue of :attr:`DsmLegResult.v_arrive` (Lambert's ``v_22``).
    r_dsm:
        Heliocentric position of the DSM / match point (km).
    max_residual:
        ``max(|match-point defect|)`` at the solution (the position rows in km,
        velocity rows in km/s; the convergence quantity).
    converged:
        ``max_residual < tol`` (residual-magnitude acceptance, by design, mirroring
        the Lambert correctors).
    solver_success, solver_nfev, solver_njev:
        ``least_squares`` diagnostics (audit trail only).
    """

    dv_dsm: Vec3
    v_arrive: Vec3
    r_dsm: Vec3
    max_residual: float
    converged: bool
    solver_success: bool = True
    solver_nfev: int = 0
    solver_njev: int = 0


def dsm_leg_correct_fbs(
    r0: Vec3,
    v0: Vec3,
    target_r: Vec3,
    tof: float,
    eta: float,
    dv0: Vec3,
    vf0: Vec3,
    *,
    mu: float = MU_SUN_KM3_S2,
    use_analytic_jac: bool = True,
    tol: float = 1.0e-9,
    max_nfev: int = 200,
) -> FbsLegCorrectionResult:
    r"""Correct a single DSM leg to match-point closure WITHOUT Lambert (#226).

    The Ellison-2018 forward-backward-shooting (FBS) alternative to the Lambert
    :func:`dsm_leg`: instead of solving a Lambert back arc, the interior impulse
    ``Δv`` is an explicit decision variable and the leg's 6-vector match-point
    defect (:func:`cyclerfinder.core.fbs_match_point.match_point_defect`) is driven
    to zero. The free variables are ``x = [Δv (3); vf (3)]`` (the impulse and the
    heliocentric arrival velocity) — 6 unknowns for the 6 defect rows — with the
    boundary positions ``r0``/``target_r``, the departure velocity ``v0``, the ToF
    and the burn fraction ``alpha = eta`` held fixed.

    At the Lambert leg's own solution the defect is exactly zero with
    ``Δv = v21 - v12`` and ``vf = v22``, so this corrector converges to the SAME
    leg geometry (the parity cross-check; see the Phase 6 test) — no Lambert solve
    anywhere on this path.

    Parameters
    ----------
    r0, v0:
        Heliocentric departure state, km / km/s.
    target_r:
        Heliocentric target position at arrival, km (the right-boundary position).
    tof, eta:
        Total leg ToF (s) and DSM burn fraction (``alpha``); same meaning as
        :func:`dsm_leg`.
    dv0, vf0:
        Seeds for the impulse and the arrival velocity (km/s). A perturbed seed
        still converges to the leg solution (demonstrated in the test).
    use_analytic_jac:
        When ``True`` (default for this opt-in entry point) the analytic FBS
        Jacobian (:func:`...match_point_defect_jacobian`, columns ``∂c/∂Δv`` and
        ``∂c/∂vf``) is supplied to ``least_squares`` as ``jac=``. When ``False`` the
        solver finite-differences the residual instead (the comparison baseline);
        both reach the same solution.
    mu, tol, max_nfev:
        Central body, convergence tolerance on the NON-DIMENSIONALISED
        ``max|defect|`` (position rows scaled by AU, velocity rows by the local
        circular speed), and solver evaluation cap.

    Returns
    -------
    FbsLegCorrectionResult
        The converged impulse, arrival velocity, DSM position, residual, and the
        solver diagnostics.
    """
    r0_a = np.asarray(r0, dtype=np.float64)
    v0_a = np.asarray(v0, dtype=np.float64)
    rf_a = np.asarray(target_r, dtype=np.float64)

    # Non-dimensionalise the defect so least_squares sees a well-scaled 6-vector:
    # position rows (km, ~1e8) by AU, velocity rows (km/s) by the local circular
    # speed. The match-point defect is multi-rooted (like Lambert — the backward
    # arc can reach the DSM via alternate conics), so good scaling keeps the solver
    # in the seed's basin and makes tol a meaningful convergence threshold.
    v_scale = float(np.sqrt(mu / float(np.linalg.norm(r0_a))))
    scale = np.array([AU_KM, AU_KM, AU_KM, v_scale, v_scale, v_scale], dtype=np.float64)

    def _leg(vf: Vec3) -> FbsLeg:
        return FbsLeg(
            r0=r0_a,
            v0=v0_a,
            rf=rf_a,
            vf=np.asarray(vf, dtype=np.float64),
            tof_s=float(tof),
            alpha=float(eta),
            mu=mu,
        )

    def _res(x: NDArray[np.float64]) -> NDArray[np.float64]:
        dv = x[0:3]
        vf = x[3:6]
        return match_point_defect(_leg(vf), dv) / scale

    def _jac(x: NDArray[np.float64]) -> NDArray[np.float64]:
        dv = x[0:3]
        vf = x[3:6]
        full = match_point_defect_jacobian(_leg(vf), dv)  # 6x9 [Δv | v0 | vf]
        # Free variables are [Δv; vf] -> stack the Δv (cols 0:3) and vf (cols 6:9),
        # row-scaled to match the non-dimensionalised residual.
        return np.column_stack([full[:, 0:3], full[:, 6:9]]) / scale[:, None]

    x0 = np.concatenate([np.asarray(dv0, dtype=np.float64), np.asarray(vf0, dtype=np.float64)])
    jac_arg: Any = _jac if use_analytic_jac else "2-point"
    sol = least_squares(
        _res, x0, jac=jac_arg, method="lm", max_nfev=max_nfev, xtol=1e-14, ftol=1e-14
    )

    dv_sol = np.asarray(sol.x[0:3], dtype=np.float64)
    vf_sol = np.asarray(sol.x[3:6], dtype=np.float64)
    leg = _leg(vf_sol)
    defect = match_point_defect(leg, dv_sol)
    # Convergence on the SCALED defect (position/AU, velocity/v_circ) so the
    # km-scale position rows do not dominate; tol is the scaled threshold.
    max_res = float(np.max(np.abs(defect / scale)))
    r_dsm, _ = propagate(r0_a, v0_a, float(eta) * float(tof), mu)
    return FbsLegCorrectionResult(
        dv_dsm=dv_sol,
        v_arrive=vf_sol,
        r_dsm=np.asarray(r_dsm, dtype=np.float64),
        max_residual=max_res,
        converged=bool(max_res < tol),
        solver_success=bool(sol.success),
        solver_nfev=int(sol.nfev),
        solver_njev=int(getattr(sol, "njev", 0) or 0),
    )


# ---------------------------------------------------------------------------
# Chained multi-leg evaluator + sequence-keyed bounds
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class DsmChainResult:
    """Outcome + audit of a chained one-DSM-per-leg evaluation/correction.

    Mirrors :class:`cyclerfinder.search.free_return.FreeReturnClosureResult`'s
    constraint-vs-evidence separation.

    Attributes
    ----------
    total_dv_kms:
        The CONSTRAINED objective: ``sum(dV_DSM) + dV_arrive`` (Takao Eq.15 without
        the powered-flyby term, which the catalogue's flyby surrogate owns), km/s.
        This is the residual the corrector / MBH drives.
    max_residual_kms:
        ``max(|residual_vector|)`` — the SAME quantity the ``converged`` decision
        is made on (the name the ``free_return`` / ``ballistic`` correctors use
        for their residual). In the default path ``residual_vector`` is the
        length-1 ``[total_dv_kms]`` so this still equals ``total_dv_kms``; in the
        charged path it is the worst per-component residual, NOT the sum
        (which is ``total_dv_kms``).
    dv_dsm_per_leg_kms:
        EMERGED per-leg DSM impulse magnitudes (the audit breakdown), km/s.
    dv_arrive_kms:
        EMERGED terminal-arrival V_inf magnitude (rendezvous) or 0 (flyby), km/s.
    vinf_in_kms, vinf_out_kms:
        EMERGED per-body incoming / outgoing V_inf magnitudes, keyed by leg index.
        Evidence -- never imposed (Takao Eq.7).
    eta_per_leg:
        The converged DSM fractions (the new genome coordinate), per leg.
    tof_days_per_leg:
        The converged per-leg ToFs (days).
    n_revs_per_leg:
        EMERGED back-arc revolution count selected on each leg (audit). All zeros
        when the chain was evaluated single-rev (``max_revs == 0``); otherwise the
        per-leg multi-rev branch the dV-minimising selection landed on.
    branch_per_leg:
        EMERGED back-arc Lambert branch label per leg (audit), pairs with
        :attr:`n_revs_per_leg`.
    t0_sec:
        Converged departure epoch (seconds).
    vinf_out0_kms, alpha0, beta0:
        The converged departure-V_inf genome (magnitude + azimuth + elevation,
        Takao Eq.5). These ENTER the dV objective through the departure state
        ``v_planet,0 + V_inf_out0`` and are moved by the corrector; the adapter
        echoes them so the landed decision vector is in the seed's coordinates.
    dsm_states:
        Per-leg ``(r_dsm_km, t_dsm_sec)`` for the audit trail / viz.
    converged:
        ``total_dv_kms < tol_kms`` -- residual-magnitude acceptance, BY DESIGN
        (mirrors ``free_return_correct``; the solver flag is secondary).
    solver_success, solver_nfev:
        DIAGNOSTIC ``least_squares`` outcome (audit trail only, never gates).
    residual_vector:
        The vector residual the corrector drives (#162). In the default
        (``charge_flyby_continuity=False``) path this is the length-1 scalar
        ``[total_dv_kms]`` -- bit-identical to the historical single-scalar
        objective. In the ``charge_flyby_continuity=True`` path it is the multi-term
        root-find vector ``[*dV_DSM_per_leg, *flyby_dv_per_intermediate_flyby,
        arrival]`` (design §3): every DSM impulse, every intermediate-flyby bend
        cost (:func:`_flyby_continuity_residual`), and the terminal arrival V_inf.
        ``converged`` in that mode is ``max(|residual_vector|) < tol_kms`` (every
        DSM AND every flyby_dv below tol), still residual-magnitude only.
    """

    total_dv_kms: float
    max_residual_kms: float
    dv_dsm_per_leg_kms: tuple[float, ...]
    dv_arrive_kms: float
    vinf_in_kms: dict[int, float]
    vinf_out_kms: dict[int, float]
    eta_per_leg: tuple[float, ...]
    tof_days_per_leg: tuple[float, ...]
    t0_sec: float
    vinf_out0_kms: float = 0.0
    alpha0: float = 0.0
    beta0: float = 0.0
    n_revs_per_leg: tuple[int, ...] = field(default_factory=tuple)
    branch_per_leg: tuple[str, ...] = field(default_factory=tuple)
    dsm_states: tuple[tuple[Vec3, float], ...] = field(default_factory=tuple)
    converged: bool = False
    solver_success: bool = True
    solver_nfev: int = 0
    residual_vector: NDArray[np.float64] = field(
        default_factory=lambda: np.zeros(1, dtype=np.float64)
    )
    alpha_int_per_leg: tuple[float, ...] = field(default_factory=tuple)
    beta_int_per_leg: tuple[float, ...] = field(default_factory=tuple)
    flyby_dv_per_flyby_kms: tuple[float, ...] = field(default_factory=tuple)


def _vinf_out_dir(v_inf: float, alpha: float, beta: float) -> Vec3:
    """Outgoing V_inf vector from magnitude + azimuth/elevation (Takao Eq.5)."""
    return np.array(
        [
            v_inf * np.cos(alpha) * np.cos(beta),
            v_inf * np.sin(alpha) * np.cos(beta),
            v_inf * np.sin(beta),
        ],
        dtype=np.float64,
    )


def _flyby_continuity_residual(
    vinf_in_vec: Vec3,
    vinf_out_vec: Vec3,
    target_code: str,
) -> float:
    """Per-flyby V_inf-continuity + bend-feasibility residual (#162, design §3).

    Returns the powered-flyby surrogate Delta V (:func:`core.flyby.flyby_dv`,
    Russell Eq.5.5) needed to convert the incoming hyperbolic excess ``vinf_in_vec``
    into the outgoing ``vinf_out_vec`` at the body ``target_code`` -- exactly ``0.0``
    when the pair is ballistic-feasible (equal magnitude AND turn within the bend
    cone), strictly positive otherwise. The planet ``(mu, rp_min)`` are resolved from
    the constants registry the same way :func:`core.flyby.flyby_dv_for` does (no
    hardcoded planet numbers). This is the bend-feasibility term the default scalar
    objective omits; in the ``charge_flyby_continuity`` path the chain genome frees
    the departure-V_inf DIRECTION so the magnitude is inherited by construction and
    the only flyby residual is this bend cost.
    """
    mu_planet = PLANETS[target_code].mu_km3_s2
    rp_min = SAFE_PERIHELION_KM[target_code]
    return flyby_dv(vinf_in_vec, vinf_out_vec, mu_planet, rp_min)


# Allowed values of the ``gradient`` opt-in (additive — default is "lambert").
_GRADIENT_LAMBERT = "lambert"
_GRADIENT_FBS_ANALYTIC = "fbs-analytic"
_GRADIENT_CHOICES = (_GRADIENT_LAMBERT, _GRADIENT_FBS_ANALYTIC)


def _eval_leg_fbs_analytic(
    r_curr: Vec3,
    v_depart: Vec3,
    tof_s: float,
    eta: float,
    r_target: Vec3,
    *,
    mu: float,
    max_revs: int = 0,
    rev_branch: tuple[int, str] | None = None,
) -> DsmLegResult:
    """Evaluate ONE leg with the FBS analytic-gradient corrector (#244 opt-in).

    The opt-in backbone: the interior impulse and the arrival velocity are produced
    by :func:`dsm_leg_correct_fbs` — the Ellison-2018 forward-backward-shooting
    match-point corrector driven by the #226 ANALYTIC Jacobian. The match-point
    defect is multi-rooted (the backward arc can reach the DSM via alternate
    conics; documented in :func:`dsm_leg_correct_fbs`), so a cold seed can converge
    to a DIFFERENT conic than the Lambert leg. To land the SAME leg the corrector
    is seeded from the Lambert ballistic solution (exactly the #243 trial design —
    "each problem's ballistic seed is built by Lambert; the FBS analytic gradient
    is the refinement / optimisation engine"), then refined with the analytic
    Jacobian. This is the honest scope of the opt-in: FBS-analytic GRADIENTS as the
    refinement engine, not a Lambert-free single-leg root-finder (the #242 finding:
    cold-seed FBS is no better as a feasibility solver).

    Raises :class:`LambertError` (caught by the chain loop as an infeasible hop) if
    EITHER the Lambert seed solve OR the FBS refinement fails — mirroring the
    Lambert lane's failure handling.
    """
    r_curr_a = np.asarray(r_curr, dtype=np.float64)
    v_dep_a = np.asarray(v_depart, dtype=np.float64)
    r_tgt_a = np.asarray(r_target, dtype=np.float64)
    # Ballistic Lambert seed for the leg's basin (the #243 trial seeding design).
    lam = dsm_leg(
        r_curr_a, v_dep_a, tof_s, eta, r_tgt_a, mu=mu, max_revs=max_revs, rev_branch=rev_branch
    )
    dv_seed = lam.v_depart_post_dsm - lam.v_arrive_pre_dsm
    res = dsm_leg_correct_fbs(
        r_curr_a,
        v_dep_a,
        r_tgt_a,
        tof_s,
        eta,
        dv_seed,
        lam.v_arrive,
        mu=mu,
        use_analytic_jac=True,
        tol=1.0e-7,
        max_nfev=200,
    )
    if not res.converged:
        raise LambertError(
            f"FBS analytic leg corrector did not converge "
            f"(max_residual={res.max_residual:.3e}); treated as an infeasible hop"
        )
    return DsmLegResult(
        v_arrive=res.v_arrive,
        v_depart_post_dsm=res.dv_dsm,  # not used downstream; kept for shape parity
        v_arrive_pre_dsm=np.zeros(3, dtype=np.float64),
        dv_dsm_kms=float(np.linalg.norm(res.dv_dsm)),
        r_dsm=res.r_dsm,
        t_dsm_sec=float(eta * tof_s),
        n_revs_chosen=int(lam.n_revs_chosen),
        branch_chosen="fbs",
    )


def evaluate_dsm_chain(
    *,
    sequence: tuple[str, ...],
    t0_sec: float,
    vinf_out0_kms: float,
    alpha0: float,
    beta0: float,
    tof_days_per_leg: tuple[float, ...],
    eta_per_leg: tuple[float, ...],
    ephem: Ephemeris,
    mu: float = MU_SUN_KM3_S2,
    rendezvous: bool = False,
    tol_kms: float = 0.1,
    max_revs: int = 0,
    rev_branch_per_leg: tuple[tuple[int, str] | None, ...] | None = None,
    charge_flyby_continuity: bool = False,
    alpha_int_per_leg: tuple[float, ...] = (),
    beta_int_per_leg: tuple[float, ...] = (),
    gradient: str = _GRADIENT_LAMBERT,
) -> DsmChainResult:
    """Evaluate a chained one-DSM-per-leg trajectory (Takao Eqs.3-7, 14-15).

    The genome strings ``len(sequence) - 1`` legs. The departure state at body 0 is
    ``v_planet,0 + V_inf_out0`` with ``V_inf_out0`` from (``vinf_out0_kms``,
    ``alpha0``, ``beta0``) (Eq.4-5). Each leg runs :func:`dsm_leg`; the heliocentric
    arrival velocity is the spacecraft's incoming velocity at the next body.

    Default path (``charge_flyby_continuity=False``)
    ------------------------------------------------
    The departing velocity for the NEXT leg is that same incoming heliocentric
    velocity (a ballistic flyby preserves heliocentric speed continuity; the
    powered-flyby surrogate, owned by ``core/flyby.py``, is applied separately at
    scoring time and is NOT charged here -- this evaluator's objective is the sum of
    DSM impulses plus the terminal arrival, Eq.15 without the P-FB term). The
    ``residual_vector`` is the length-1 scalar ``[total_dv_kms]``, BIT-IDENTICAL to
    the historical behaviour.

    Charged path (``charge_flyby_continuity=True``, #162 / design §3)
    ----------------------------------------------------------------
    The next leg's departure V_inf DIRECTION at each intermediate flyby is a free
    genome coordinate (``alpha_int_per_leg[k]`` / ``beta_int_per_leg[k]`` for the
    ``k``-th intermediate body, ``k = 0 .. n_legs-2``); its MAGNITUDE is inherited
    from the arrival V_inf (so V_inf magnitude continuity holds by construction).
    Per intermediate flyby the bend-feasibility cost
    :func:`_flyby_continuity_residual` (``flyby_dv``, 0 when ballistic-feasible) is
    charged, and the result's ``residual_vector`` becomes the multi-term root-find
    vector ``[*dV_DSM_per_leg, *flyby_dv_per_flyby, arrival]``. ``converged`` is then
    ``max(|residual_vector|) < tol_kms`` (every DSM AND every flyby_dv below tol).
    This is the only path that rewards the sourced bend-feasible low-V_inf basin; it
    is default-off and adds zero cost to every existing caller.

    Body epochs chain as ``t_j = t0 + sum(tau_i)`` (Eq.3); body states pulled from
    ``ephem`` at those epochs.

    Returns
    -------
    DsmChainResult
        ``total_dv_kms`` is the scalar objective (back-compat / audit);
        ``residual_vector`` is what the corrector drives. The per-leg DSM impulses,
        emerged V_inf, and DSM states are the evidence/audit trail.
    """
    if gradient not in _GRADIENT_CHOICES:
        raise ValueError(f"gradient must be one of {_GRADIENT_CHOICES}, got {gradient!r}")
    n_legs = len(sequence) - 1
    if n_legs < 1:
        raise ValueError("sequence must name at least two bodies (one leg)")
    if len(tof_days_per_leg) != n_legs or len(eta_per_leg) != n_legs:
        raise ValueError(
            f"tof_days_per_leg and eta_per_leg must each have {n_legs} entries "
            f"(one per leg), got {len(tof_days_per_leg)} / {len(eta_per_leg)}"
        )
    n_flybys = n_legs - 1  # intermediate flybys (bodies 1 .. n_legs-1)
    if charge_flyby_continuity and (
        len(alpha_int_per_leg) != n_flybys or len(beta_int_per_leg) != n_flybys
    ):
        raise ValueError(
            f"charge_flyby_continuity requires alpha_int_per_leg / beta_int_per_leg "
            f"of length {n_flybys} (one per intermediate flyby), got "
            f"{len(alpha_int_per_leg)} / {len(beta_int_per_leg)}"
        )

    r0, v_planet0 = ephem.state(sequence[0], t0_sec)
    v_out0 = _vinf_out_dir(vinf_out0_kms, alpha0, beta0)
    v_depart = np.asarray(v_planet0, dtype=np.float64) + v_out0

    dv_dsm: list[float] = []
    dsm_states: list[tuple[Vec3, float]] = []
    n_revs_legs: list[int] = []
    branch_legs: list[str] = []
    flyby_dv_per_flyby: list[float] = []
    vinf_in: dict[int, float] = {}
    vinf_out: dict[int, float] = {0: float(vinf_out0_kms)}

    t_cursor = t0_sec
    r_curr = np.asarray(r0, dtype=np.float64)
    feasible = True
    for i in range(n_legs):
        tof_s = tof_days_per_leg[i] * DAY_S
        t_arrive = t_cursor + tof_s
        target_body = sequence[i + 1]
        r_target, v_planet_target = ephem.state(target_body, t_arrive)
        rb = rev_branch_per_leg[i] if rev_branch_per_leg is not None else None
        try:
            if gradient == _GRADIENT_FBS_ANALYTIC:
                # FBS match point needs alpha strictly inside (0,1); clamp the eta
                # seed the same way the Lambert dsm_leg does (the exact endpoints are
                # singular). The Lambert seed inside _eval_leg_fbs_analytic honours
                # max_revs / rev_branch so the FBS refinement lands the same basin.
                eta_clamped = min(max(eta_per_leg[i], _ETA_EPS), 1.0 - _ETA_EPS)
                leg = _eval_leg_fbs_analytic(
                    r_curr,
                    v_depart,
                    tof_s,
                    eta_clamped,
                    np.asarray(r_target),
                    mu=mu,
                    max_revs=max_revs,
                    rev_branch=rb,
                )
            else:
                leg = dsm_leg(
                    r_curr,
                    v_depart,
                    tof_s,
                    eta_per_leg[i],
                    np.asarray(r_target),
                    mu=mu,
                    max_revs=max_revs,
                    rev_branch=rb,
                )
        except (LambertError, KeplerError, ValueError):
            # A hop into a too-energetic departure can drive the ballistic
            # propagation hyperbolic past Newton convergence, or the back-arc
            # Lambert degenerate; report the chain INFEASIBLE rather than crashing
            # the whole MBH search (the optimiser simply rejects this hop).
            feasible = False
            break
        dv_dsm.append(leg.dv_dsm_kms)
        dsm_states.append((leg.r_dsm, leg.t_dsm_sec))
        n_revs_legs.append(leg.n_revs_chosen)
        branch_legs.append(leg.branch_chosen)
        v_planet_target_arr = np.asarray(v_planet_target, dtype=np.float64)
        v_inf_in_vec = leg.v_arrive - v_planet_target_arr
        vinf_in_mag = float(np.linalg.norm(v_inf_in_vec))
        vinf_in[i + 1] = vinf_in_mag

        is_intermediate = i < n_legs - 1
        if charge_flyby_continuity and is_intermediate:
            # The next leg departs on a FREE-direction V_inf at the inherited
            # magnitude (V_inf-magnitude continuity holds by construction); the bend
            # cost of turning v_inf_in -> v_inf_out at this body is the per-flyby
            # residual term (0 iff ballistic-feasible).
            v_inf_out_vec = _vinf_out_dir(vinf_in_mag, alpha_int_per_leg[i], beta_int_per_leg[i])
            flyby_dv_per_flyby.append(
                _flyby_continuity_residual(v_inf_in_vec, v_inf_out_vec, target_body)
            )
            v_depart = v_planet_target_arr + v_inf_out_vec
            vinf_out[i + 1] = float(np.linalg.norm(v_inf_out_vec))
        else:
            # Ballistic-flyby heliocentric continuity for the next leg's departure.
            # (The powered-flyby bend cost is the catalogue scorer's job in this
            # default path; here the next leg simply departs on the arrival velocity.)
            v_depart = np.asarray(leg.v_arrive, dtype=np.float64)
            vinf_out[i + 1] = vinf_in_mag
        r_curr = np.asarray(r_target, dtype=np.float64)
        t_cursor = t_arrive

    if not feasible:
        return DsmChainResult(
            total_dv_kms=float("inf"),
            max_residual_kms=float("inf"),
            dv_dsm_per_leg_kms=tuple(dv_dsm),
            dv_arrive_kms=float("inf"),
            vinf_in_kms=vinf_in,
            vinf_out_kms=vinf_out,
            eta_per_leg=tuple(eta_per_leg),
            tof_days_per_leg=tuple(tof_days_per_leg),
            t0_sec=float(t0_sec),
            vinf_out0_kms=float(vinf_out0_kms),
            alpha0=float(alpha0),
            beta0=float(beta0),
            n_revs_per_leg=tuple(n_revs_legs),
            branch_per_leg=tuple(branch_legs),
            dsm_states=tuple(dsm_states),
            converged=False,
            residual_vector=np.array([float("inf")], dtype=np.float64),
            alpha_int_per_leg=tuple(alpha_int_per_leg),
            beta_int_per_leg=tuple(beta_int_per_leg),
            flyby_dv_per_flyby_kms=tuple(flyby_dv_per_flyby),
        )

    dv_arrive = vinf_in[n_legs] if rendezvous else 0.0
    total_dv = float(sum(dv_dsm) + dv_arrive)
    if charge_flyby_continuity:
        residual_vector = np.array([*dv_dsm, *flyby_dv_per_flyby, dv_arrive], dtype=np.float64)
        converged = bool(np.max(np.abs(residual_vector)) < tol_kms)
    else:
        residual_vector = np.array([total_dv], dtype=np.float64)
        converged = total_dv < tol_kms
    return DsmChainResult(
        total_dv_kms=total_dv,
        # The quantity `converged` was decided on (NOT the total_dv sum; in the
        # charged path those differ and callers were getting the wrong one).
        max_residual_kms=float(np.max(np.abs(residual_vector))),
        dv_dsm_per_leg_kms=tuple(dv_dsm),
        dv_arrive_kms=float(dv_arrive),
        vinf_in_kms=vinf_in,
        vinf_out_kms=vinf_out,
        eta_per_leg=tuple(eta_per_leg),
        tof_days_per_leg=tuple(tof_days_per_leg),
        t0_sec=float(t0_sec),
        vinf_out0_kms=float(vinf_out0_kms),
        alpha0=float(alpha0),
        beta0=float(beta0),
        n_revs_per_leg=tuple(n_revs_legs),
        branch_per_leg=tuple(branch_legs),
        dsm_states=tuple(dsm_states),
        converged=converged,
        residual_vector=residual_vector,
        alpha_int_per_leg=tuple(alpha_int_per_leg),
        beta_int_per_leg=tuple(beta_int_per_leg),
        flyby_dv_per_flyby_kms=tuple(flyby_dv_per_flyby),
    )


# ---------------------------------------------------------------------------
# Sequence-keyed automatic bounds (Takao Appendix A.1-A.3)
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class DsmBounds:
    """Box bounds for the DSM-chain decision vector (Takao Appendix A.1-A.3).

    Layout matches :func:`dsm_chain_decision_vector`:
    ``[t0_sec, vinf_out0, alpha0, beta0, *tof_days_per_leg, *eta_per_leg]`` in the
    default path, with ``2*(n_legs-1)`` intermediate-flyby direction coords
    ``[*alpha_int, *beta_int]`` appended when ``charge_flyby_continuity`` is on
    (#162).
    """

    lower: NDArray[np.float64]
    upper: NDArray[np.float64]


def _synodic_period_days(inner: str, outer: str) -> float:
    """Synodic period of two bodies (days), from their sidereal periods.

    A same-body pair (e.g. the E-E leg of a resonant return, as produced by the
    descriptor-seed lane) has equal periods: the pair never laps, and the synodic
    period's mathematically-correct limit is +inf. Return that limit explicitly
    rather than evaluating ``1/0`` (which yields the same +inf but trips numpy's
    divide-by-zero RuntimeWarning — #217). Callers that need finite box bounds
    must cap the resulting infinite A.2 ToF upper bound themselves (the
    descriptor-seed lane already does).
    """
    p_in = 2.0 * np.pi * np.sqrt((PLANETS[inner].sma_au * AU_KM) ** 3 / MU_SUN_KM3_S2) / DAY_S
    p_out = 2.0 * np.pi * np.sqrt((PLANETS[outer].sma_au * AU_KM) ** 3 / MU_SUN_KM3_S2) / DAY_S
    if p_in == p_out:
        return float("inf")
    return float(abs(1.0 / (1.0 / p_in - 1.0 / p_out)))


def _hohmann_tof_days(inner: str, outer: str) -> float:
    """Hohmann (half-period of the transfer ellipse) ToF inner<->outer (days)."""
    a_t_km = 0.5 * (PLANETS[inner].sma_au + PLANETS[outer].sma_au) * AU_KM
    return float(np.pi * np.sqrt(a_t_km**3 / MU_SUN_KM3_S2) / DAY_S)


def sequence_keyed_bounds(
    *,
    sequence: tuple[str, ...],
    t0_window_sec: tuple[float, float],
    vinf_out0_bounds_kms: tuple[float, float] = (1.0, 5.1),
    charge_flyby_continuity: bool = False,
) -> DsmBounds:
    """Automatic box bounds from the body sequence (Takao Appendix A.1-A.3).

    Only the departure epoch window and the departure V_inf magnitude bound need to
    be supplied (Takao: "users need to specify only the epoch and the V_inf
    magnitude at launch"); the rest are sequence-keyed:

    * ``alpha in [-pi, pi]``, ``beta in [-pi/2, pi/2]``, ``eta in [0, 1]`` (A.1).
    * Per-leg ToF: an inner pair (both bodies inner planets, here the E<->M class)
      gets ``[30 d, P_s + P_H]`` (A.2); any leg touching an outer body gets
      ``[0.3 P_H, 1.3 P_H]`` (A.3). ``P_s`` synodic, ``P_H`` Hohmann for that leg's
      pair. A same-body inner leg (e.g. an E-E resonant return) has ``P_s = +inf``
      — the A.2 upper bound is then ``+inf`` and callers needing a finite box must
      cap it (#217; the descriptor-seed lane does).

    The departure V_inf default ``[1, 5.1]`` km/s is Takao's Earth-departure window
    (5.1 km/s = the 1:2 Earth-resonant cap; 1 km/s prevents low-velocity flybys).

    When ``charge_flyby_continuity`` is on (#162) the box is extended by the
    ``2*(n_legs-1)`` intermediate-flyby departure-direction coords, each boxed like
    the leg-0 direction (``alpha_int in [-pi, pi]``, ``beta_int in [-pi/2, pi/2]``).
    """
    n_legs = len(sequence) - 1
    if n_legs < 1:
        raise ValueError("sequence must name at least two bodies (one leg)")

    t0_lo, t0_hi = t0_window_sec
    vinf_lo, vinf_hi = vinf_out0_bounds_kms

    tof_lo: list[float] = []
    tof_hi: list[float] = []
    # An "inner" pair: both endpoints are at/inside Mars (the E<->M resonant class
    # Takao's A.2 covers). Anything else (Jupiter/Saturn etc.) uses the A.3 window.
    inner_codes = {"Me", "V", "E", "M"}
    for i in range(n_legs):
        a, b = sequence[i], sequence[i + 1]
        # Order the pair by semi-major axis for the synodic/Hohmann helpers.
        inner, outer = (a, b) if PLANETS[a].sma_au <= PLANETS[b].sma_au else (b, a)
        p_h = _hohmann_tof_days(inner, outer)
        if a in inner_codes and b in inner_codes:
            p_s = _synodic_period_days(inner, outer)
            tof_lo.append(30.0)
            tof_hi.append(p_s + p_h)
        else:
            tof_lo.append(0.3 * p_h)
            tof_hi.append(1.3 * p_h)

    n_flybys = n_legs - 1 if charge_flyby_continuity else 0
    dir_lo = [-np.pi] * n_flybys + [-0.5 * np.pi] * n_flybys
    dir_hi = [np.pi] * n_flybys + [0.5 * np.pi] * n_flybys

    lower = np.array(
        [t0_lo, vinf_lo, -np.pi, -0.5 * np.pi, *tof_lo, *([0.0] * n_legs), *dir_lo],
        dtype=np.float64,
    )
    upper = np.array(
        [t0_hi, vinf_hi, np.pi, 0.5 * np.pi, *tof_hi, *([1.0] * n_legs), *dir_hi],
        dtype=np.float64,
    )
    return DsmBounds(lower=lower, upper=upper)


def dsm_chain_decision_vector(
    *,
    t0_sec: float,
    vinf_out0_kms: float,
    alpha0: float,
    beta0: float,
    tof_days_per_leg: tuple[float, ...],
    eta_per_leg: tuple[float, ...],
    alpha_int_per_leg: tuple[float, ...] = (),
    beta_int_per_leg: tuple[float, ...] = (),
) -> NDArray[np.float64]:
    """Pack the chain genome into the flat decision vector the corrector consumes.

    Layout ``[t0_sec, vinf_out0, alpha0, beta0, *tof_days_per_leg, *eta_per_leg]``
    (matches :class:`DsmBounds`), with ``[*alpha_int_per_leg, *beta_int_per_leg]``
    appended when the intermediate-flyby direction coords are supplied (#162; the
    ``charge_flyby_continuity`` path). When both are empty the vector is
    bit-identical to the historical layout.
    """
    return np.array(
        [
            t0_sec,
            vinf_out0_kms,
            alpha0,
            beta0,
            *tof_days_per_leg,
            *eta_per_leg,
            *alpha_int_per_leg,
            *beta_int_per_leg,
        ],
        dtype=np.float64,
    )


def _unpack(
    x: NDArray[np.float64], n_legs: int, *, charge_flyby_continuity: bool = False
) -> dict[str, Any]:
    """Inverse of :func:`dsm_chain_decision_vector`.

    With ``charge_flyby_continuity`` the trailing ``2*(n_legs-1)`` coords are read
    back as the intermediate-flyby departure-direction genome
    (``alpha_int_per_leg`` / ``beta_int_per_leg``).
    """
    base: dict[str, Any] = {
        "t0_sec": float(x[0]),
        "vinf_out0_kms": float(x[1]),
        "alpha0": float(x[2]),
        "beta0": float(x[3]),
        "tof_days_per_leg": tuple(float(v) for v in x[4 : 4 + n_legs]),
        "eta_per_leg": tuple(float(v) for v in x[4 + n_legs : 4 + 2 * n_legs]),
    }
    if charge_flyby_continuity:
        n_flybys = n_legs - 1
        off = 4 + 2 * n_legs
        base["alpha_int_per_leg"] = tuple(float(v) for v in x[off : off + n_flybys])
        base["beta_int_per_leg"] = tuple(float(v) for v in x[off + n_flybys : off + 2 * n_flybys])
    return base


def dsm_chain_correct(
    x0: NDArray[np.float64],
    *,
    sequence: tuple[str, ...],
    ephem: Ephemeris,
    bounds: DsmBounds | None = None,
    mu: float = MU_SUN_KM3_S2,
    rendezvous: bool = False,
    tol_kms: float = 0.1,
    max_nfev: int = 200,
    max_revs: int = 0,
    rev_branch_per_leg: tuple[tuple[int, str] | None, ...] | None = None,
    charge_flyby_continuity: bool = False,
    gradient: str = _GRADIENT_LAMBERT,
) -> DsmChainResult:
    """Drive the chained-DSM residual to a minimum with bounded least-squares.

    Free variables = the full decision vector
    ``[t0_sec, vinf_out0, alpha0, beta0, *tof, *eta]`` (see
    :func:`dsm_chain_decision_vector`), with the ``2*(n_legs-1)``
    intermediate-flyby direction coords appended when ``charge_flyby_continuity``.

    Default path (``charge_flyby_continuity=False``): the residual is the length-1
    scalar ``total_dv_kms`` (Takao Eq.15 minus the P-FB term) -- bit-identical to
    the historical corrector. ``least_squares`` minimises it inside ``bounds``;
    converged iff ``total_dv_kms < tol_kms``.

    Charged path (``charge_flyby_continuity=True``, #162): the residual is the
    VECTOR ``result.residual_vector`` = ``[*dV_DSM, *flyby_dv_per_flyby, arrival]``,
    a true multi-term root-find (``least_squares`` natively handles vector
    residuals). Converged iff ``max(|residual_vector|) < tol_kms`` -- every DSM AND
    every flyby bend cost below tol. This is the only mode that rewards the sourced
    bend-feasible low-V_inf basin.

    This is the local-solve primitive the MBH wrapper drives via
    :func:`make_dsm_chain_step`; MBH supplies the basin hops, this supplies the
    refinement.

    Gradient backbone (``gradient``, #244 opt-in, default ``"lambert"``)
    -------------------------------------------------------------------
    ``"lambert"`` (default): each leg is evaluated by the Lambert ``dsm_leg`` and
    the outer ``least_squares`` finite-differences over the genome — BIT-IDENTICAL
    to the historical corrector. ``"fbs-analytic"``: each leg is evaluated by the
    Lambert-free Ellison-2018 forward-backward-shooting match-point corrector
    (:func:`_eval_leg_fbs_analytic` → :func:`dsm_leg_correct_fbs`), whose inner
    root-find is driven by the #226 ANALYTIC Jacobian — no Lambert anywhere on the
    leg. At a leg's own solution the FBS defect is zero with the same impulse /
    arrival velocity Lambert produces (#226 parity), so the two backbones reach the
    same chain geometry; the opt-in only changes HOW each leg is solved. Default-off
    so every existing caller's result is unchanged.
    """
    n_legs = len(sequence) - 1

    def _eval(x: NDArray[np.float64]) -> DsmChainResult:
        params = _unpack(
            np.asarray(x, dtype=np.float64),
            n_legs,
            charge_flyby_continuity=charge_flyby_continuity,
        )
        return evaluate_dsm_chain(
            sequence=sequence,
            ephem=ephem,
            mu=mu,
            rendezvous=rendezvous,
            tol_kms=tol_kms,
            max_revs=max_revs,
            rev_branch_per_leg=rev_branch_per_leg,
            charge_flyby_continuity=charge_flyby_continuity,
            gradient=gradient,
            **params,
        )

    # Residual length the optimiser expects (so an INFEASIBLE hop returns a finite
    # penalty vector of the right shape rather than crashing least_squares).
    res_len = (n_legs + (n_legs - 1) + 1) if charge_flyby_continuity else 1

    def _res(x: NDArray[np.float64]) -> NDArray[np.float64]:
        r = _eval(x)
        vec = np.asarray(r.residual_vector, dtype=np.float64)
        if vec.shape != (res_len,) or not np.all(np.isfinite(vec)):
            # Infeasible chain (broke early) -> uniform large penalty of the
            # expected length; the optimiser rejects this direction.
            return np.full(res_len, 1.0e3, dtype=np.float64)
        return vec

    x0_arr = np.asarray(x0, dtype=np.float64)
    if bounds is not None:
        # An MBH hop can perturb the seed just outside the box; ``trf`` rejects an
        # out-of-bounds x0 outright. Clip into the box (strictly interior, to avoid
        # the ``lb >= ub`` / on-edge degeneracies trf also rejects) so the hop's
        # intent is preserved rather than crashing the whole search.
        span = bounds.upper - bounds.lower
        eps = np.where(span > 0.0, 1.0e-9 * span, 0.0)
        x0_arr = np.clip(x0_arr, bounds.lower + eps, bounds.upper - eps)
        sol = least_squares(
            _res,
            x0_arr,
            bounds=(bounds.lower, bounds.upper),
            method="trf",
            max_nfev=max_nfev,
            xtol=1e-12,
            ftol=1e-12,
        )
    else:
        sol = least_squares(_res, x0_arr, method="trf", max_nfev=max_nfev, xtol=1e-12, ftol=1e-12)

    result = _eval(np.asarray(sol.x, dtype=np.float64))
    # The evaluator already set ``converged`` from the residual; carry the solver
    # diagnostics onto the otherwise-complete result (replace keeps every field,
    # including the converged departure-V_inf genome, in sync).
    return replace(result, solver_success=bool(sol.success), solver_nfev=int(sol.nfev))


# ---------------------------------------------------------------------------
# MBH adapter: wrap dsm_chain_correct as an objective_and_solve closure.
# This imports and CALLS the MBH machinery; it never edits mbh.py (the wrapper
# is generic and takes this closure -- see docs/notes/2026-06-07-mbh-wrapper.md).
# ---------------------------------------------------------------------------


def make_dsm_chain_step(
    *,
    sequence: tuple[str, ...],
    ephem: Ephemeris,
    bounds: DsmBounds | None = None,
    mu: float = MU_SUN_KM3_S2,
    rendezvous: bool = False,
    tol_kms: float = 0.1,
    max_nfev: int = 200,
    max_revs: int = 0,
    rev_branch_per_leg: tuple[tuple[int, str] | None, ...] | None = None,
    charge_flyby_continuity: bool = False,
) -> Callable[[np.ndarray, np.random.Generator], MBHStep]:
    """Adapter: :func:`dsm_chain_correct` as an MBH local-solve closure.

    Genome ``x = [t0_sec, vinf_out0, alpha0, beta0, *tof_days, *eta]`` (see
    :func:`dsm_chain_decision_vector`), with the ``2*(n_legs-1)`` intermediate-flyby
    direction coords appended when ``charge_flyby_continuity``. Objective =
    ``total_dv_kms`` (the scalar audit sum); feasible = ``converged`` (which is
    ``max(|residual_vector|) < tol_kms`` in the charged path). The emerged per-body
    V_inf is carried in ``info`` as the EVIDENCE the sourced-anchor gate compares
    against (it is never the objective, which would impose it). Deterministic:
    ``rng`` accepted for the generic signature and ignored.
    """
    n_legs = len(sequence) - 1

    def step(x: np.ndarray, rng: np.random.Generator) -> MBHStep:
        r = dsm_chain_correct(
            np.asarray(x, dtype=np.float64),
            sequence=sequence,
            ephem=ephem,
            bounds=bounds,
            mu=mu,
            rendezvous=rendezvous,
            tol_kms=tol_kms,
            max_nfev=max_nfev,
            max_revs=max_revs,
            rev_branch_per_leg=rev_branch_per_leg,
            charge_flyby_continuity=charge_flyby_continuity,
        )
        # The corrector moves the WHOLE genome -- including the departure V_inf
        # (vinf_out0/alpha0/beta0) and (charged path) the intermediate-flyby
        # departure directions -- so the landed vector takes the corrector's
        # converged values, not the seed's. This is what lets MBH hop in the
        # departure-energy direction too (the 6.44Gg3 probe needs that lever).
        landed_x = dsm_chain_decision_vector(
            t0_sec=r.t0_sec,
            vinf_out0_kms=r.vinf_out0_kms,
            alpha0=r.alpha0,
            beta0=r.beta0,
            tof_days_per_leg=r.tof_days_per_leg,
            eta_per_leg=r.eta_per_leg,
            alpha_int_per_leg=r.alpha_int_per_leg,
            beta_int_per_leg=r.beta_int_per_leg,
        )
        return MBHStep(
            x=landed_x,
            objective=float(r.total_dv_kms),
            feasible=bool(r.converged),
            info={
                "total_dv_kms": float(r.total_dv_kms),
                "max_residual_kms": float(np.max(np.abs(r.residual_vector))),
                "dv_dsm_per_leg_kms": tuple(r.dv_dsm_per_leg_kms),
                "flyby_dv_per_flyby_kms": tuple(r.flyby_dv_per_flyby_kms),
                "dv_arrive_kms": float(r.dv_arrive_kms),
                "vinf_in_kms": dict(r.vinf_in_kms),
                "vinf_out_kms": dict(r.vinf_out_kms),
                "eta_per_leg": tuple(r.eta_per_leg),
                "tof_days_per_leg": tuple(r.tof_days_per_leg),
                "n_revs_per_leg": tuple(r.n_revs_per_leg),
                "branch_per_leg": tuple(r.branch_per_leg),
                "alpha_int_per_leg": tuple(r.alpha_int_per_leg),
                "beta_int_per_leg": tuple(r.beta_int_per_leg),
                "solver_nfev": int(r.solver_nfev),
            },
        )

    _ = n_legs  # documented genome arity; kept for adapter symmetry/readability
    return step
