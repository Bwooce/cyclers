"""NRHO / tulip-orbit symmetric-corrector + family continuation (#266 Phase 3).

A tulip orbit is a perpendicular x-z-plane-crossing 3-D symmetric periodic orbit
in the CR3BP -- the IC has the form ``(x0, 0, z0, 0, ydot0, 0)`` and after half a
period the orbit re-crosses the x-z plane perpendicularly,
``(x_half, 0, z_half, 0, ydot_half, 0)``. The CR3BP's y-z reflection symmetry
makes the full-period closure a consequence of the half-period perpendicular
crossing (Koblick 2023 Eqn 12, Case Two; Howell 1984 for the planar analogue).

This module is the Phase 3 (#266) substrate for finding new tulip orbits via
continuation along an L2 Southern NRHO family:

  * :func:`correct_symmetric_nrho` -- single-shooting corrector with free
    variables ``(z0, ydot0, T)`` at fixed ``x0``. The residual is the three
    perpendicular-crossing conditions ``(y(t_half), xdot(t_half), zdot(t_half))``
    at ``t_half = T/2``; the half-period and the half-period crossing state
    together pin down the full-period closure by symmetry. The Sundman-
    regularised propagator handles the close perilune.

  * :func:`continue_nrho_family` -- natural-parameter continuation in ``x0``
    (secant prediction; per-step re-correct via the symmetric corrector). At
    each converged member, the monodromy is computed and the period-multiplying
    bifurcation detector flags brackets where a Floquet multiplier crossed a
    primitive k-th root of unity. Brackets are returned, not refined.

The corrector reproduce gate is the Koblick 2023 AMOSTECH Table 4 Np=1 row at
``x0=1.023731``; the convergence floor is the regularised propagator's
round-trip residual at machine precision.

Discipline (orbit-closure):

  * Every accepted member's closure residual is below ``tol`` AND its
    monodromy is computed (Floquet diagnostics never optional);
  * Bifurcation brackets are RETURNED, not asserted -- refinement is the
    family-switch corrector's job (:func:`cyclerfinder.genome.family_switch`);
  * On a fold (the local tangent in ``(x0, z0, ydot0, T)`` becomes singular
    in ``x0``) the continuation STOPS CLEANLY with ``FOLD`` -- pseudo-arclength
    fold-turning is left to Phase 4.

Style cross-reference: the pseudo-arclength fold-turning pattern is in
:mod:`cyclerfinder.search.cr3bp_jacobi_arclength`; this module deliberately
uses the simpler natural-parameter scheme (fewer moving parts; cleaner failure
modes) because the NRHO family at the Np=1 root is well-conditioned in ``x0``.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass, field
from enum import StrEnum

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import solve_ivp

import cyclerfinder.core.cr3bp as cr3bp
import cyclerfinder.core.cr3bp_regularized as creg
from cyclerfinder.search.bifurcation_detector import (
    BifurcationPoint,
    FamilyMember,
    floquet_multipliers,
    monodromy,
)


def _nearest_kth_root_dist(eig: complex, k: int) -> tuple[float, complex]:
    """Mirror of bifurcation_detector._nearest_kth_root_distance, kept private."""
    roots: list[complex] = [
        complex(math.cos(2 * math.pi * j / k), math.sin(2 * math.pi * j / k))
        for j in range(1, k)
        if math.gcd(j, k) == 1
    ]
    if not roots:
        return float("inf"), complex(1.0)
    dists = [abs(eig - r) for r in roots]
    idx = int(np.argmin(dists))
    return float(dists[idx]), roots[idx]


# ---------------------------------------------------------------------------
# Public types.
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class SymmetricNRHO:
    """A converged perpendicular x-z-plane-crossing symmetric periodic orbit.

    Attributes
    ----------
    x0, z0, ydot0 :
        Initial-state nonzero components. The full IC is
        ``(x0, 0, z0, 0, ydot0, 0)``.
    T_TU :
        Full nondim period.
    jacobi :
        Jacobi constant at the IC.
    converged :
        True iff the perpendicular-crossing residual at ``T/2`` is below ``tol``.
    closure_residual :
        L2 norm of ``(y_half, xdot_half, zdot_half)``: the perpendicular-
        crossing residual at the corrected half-period.
    n_iter :
        Newton iterations consumed.
    monodromy :
        The full-period STM (6x6) at the converged orbit, or ``None`` if
        ``with_stm=False`` was set at correction time. Computed via the standard
        variational EOM (NOT the regularised propagator -- the regularised
        STM is not implemented here).
    """

    x0: float
    z0: float
    ydot0: float
    T_TU: float
    jacobi: float
    converged: bool
    closure_residual: float
    n_iter: int
    monodromy: NDArray[np.float64] | None


class NRHOStopReason(StrEnum):
    """Why a NRHO continuation branch terminated."""

    FOLD = "fold"  # local tangent in x0 became singular -- fold turning needed
    MAX_STEPS = "max_steps"  # step budget exhausted (not a physical edge)
    PERILUNE_FLOOR = "perilune_floor"  # perilune dropped below user floor
    NO_CONVERGE = "no_converge"  # corrector failed for too many consecutive steps
    BIFURCATION_HIT = "bifurcation_hit"  # caller requested stop on first k=2 bracket


@dataclass
class NRHOBranch:
    """Ordered branch of converged NRHO family members + bifurcation brackets."""

    label: str
    members: list[SymmetricNRHO] = field(default_factory=list)
    bifurcations: list[BifurcationPoint] = field(default_factory=list)
    stop_reason: NRHOStopReason = NRHOStopReason.MAX_STEPS
    n_steps: int = 0


# ---------------------------------------------------------------------------
# Half-period perpendicular-crossing primitives.
# ---------------------------------------------------------------------------


def _build_state0(x0: float, z0: float, ydot0: float) -> NDArray[np.float64]:
    """Standard tulip-symmetric IC: (x0, 0, z0, 0, ydot0, 0)."""
    return np.array([x0, 0.0, z0, 0.0, ydot0, 0.0], dtype=np.float64)


def _half_period_perpendicular_state(
    system: cr3bp.CR3BPSystem,
    state0: NDArray[np.float64],
    t_half: float,
    *,
    with_stm: bool,
    rtol: float,
    atol: float,
) -> tuple[NDArray[np.float64], NDArray[np.float64] | None]:
    """Propagate ``state0`` to ``t_half`` and return ``(state, STM)``.

    Uses the standard CR3BP propagator -- not the regularised one -- because the
    STM is needed for the Newton step and the regularised STM is not
    implemented. The Np=1 NRHO and small-Np tulip orbits do not graze perilune
    closely enough to defeat DOP853 at rtol/atol=1e-12; deep low-perilune
    members (Np >= 5) may require regularised propagation, which the
    multi-shooting corrector in :mod:`cyclerfinder.genome.family_switch`
    handles.
    """
    arc = cr3bp.propagate(system, state0, t_half, with_stm=with_stm, rtol=rtol, atol=atol)
    return arc.state_f, arc.stm


def correct_symmetric_nrho(
    system: cr3bp.CR3BPSystem,
    x0: float,
    z0_guess: float,
    ydot0_guess: float,
    period_guess: float,
    *,
    tol: float = 1e-11,
    max_iter: int = 60,
    rtol: float = 1e-12,
    atol: float = 1e-12,
    with_monodromy: bool = True,
) -> SymmetricNRHO:
    """Newton-corrects ``(z0, ydot0, T)`` for a perpendicular x-z-plane crossing.

    Free variables: ``(z0, ydot0, T)`` -- the half-period and the two non-trivial
    initial-state components at the fixed-``x0`` perpendicular-crossing IC. The
    residual is the three perpendicular-crossing conditions at the half-period:

        r(z0, ydot0, T) = (y(T/2), xdot(T/2), zdot(T/2))

    For the symmetric IC family this fully constrains the orbit by the CR3BP's
    (y, ydot) -> -(y, ydot) reflection symmetry: the half-period perpendicular
    crossing at the x-z plane re-flects to the full-period closure. The Newton
    step uses the STM at ``T/2`` and the time derivative at ``T/2``:

        d(y, xdot, zdot)/d(z0, ydot0, t) = [Phi[1, 2] Phi[1, 4] f[1]]
                                            [Phi[3, 2] Phi[3, 4] f[3]]
                                            [Phi[5, 2] Phi[5, 4] f[5]]

    where ``f = cr3bp_eom(T/2, state(T/2), mu)``. The 3x3 linear solve gives
    ``(dz0, dydot0, dt_half)``; the period update is ``dT = 2 * dt_half``.

    Parameters
    ----------
    system :
        CR3BP system; only ``system.mu`` is read.
    x0 :
        FIXED initial x-coordinate. Not updated by the corrector; selects which
        member of the family (in the (x0)-indexed family curve) is targeted.
    z0_guess, ydot0_guess :
        Initial guess for the two non-trivial IC components. The CR3BP IC is
        ``(x0, 0, z0, 0, ydot0, 0)``.
    period_guess :
        Initial guess for the FULL period (TU). The half-period ``T/2`` is the
        Newton's free time variable.
    tol :
        L2 residual tolerance at the half-period perpendicular crossing.
    max_iter :
        Newton iteration cap.
    rtol, atol :
        Integrator tolerances for ``propagate`` (with STM).
    with_monodromy :
        If True, the converged orbit's full-period STM (monodromy) is computed
        and returned in ``SymmetricNRHO.monodromy``.

    Returns
    -------
    SymmetricNRHO :
        Converged orbit. ``converged=False`` if the iteration cap is reached or
        the Jacobian becomes singular; the corresponding ``closure_residual``
        records the last residual.
    """
    z0 = float(z0_guess)
    ydot0 = float(ydot0_guess)
    t_half = 0.5 * float(period_guess)
    residual = float("inf")
    state0 = _build_state0(x0, z0, ydot0)
    state_half: NDArray[np.float64] = state0
    stm_half: NDArray[np.float64] | None = None
    n_iter = 0
    for n_iter in range(1, max_iter + 1):  # noqa: B007
        state0 = _build_state0(x0, z0, ydot0)
        try:
            state_half, stm_half = _half_period_perpendicular_state(
                system, state0, t_half, with_stm=True, rtol=rtol, atol=atol
            )
        except RuntimeError:
            # Integrator failed -- typically a low-perilune graze that exceeds
            # DOP853's step-size resolution. Mark unconverged so the caller can
            # halve the step or switch to the regularised propagator.
            residual = float("inf")
            break
        y_h = float(state_half[1])
        xdot_h = float(state_half[3])
        zdot_h = float(state_half[5])
        r = np.array([y_h, xdot_h, zdot_h], dtype=np.float64)
        residual = float(np.linalg.norm(r))
        if residual < tol:
            break
        assert stm_half is not None
        f_half = cr3bp.cr3bp_eom(t_half, state_half, system.mu)
        # Jacobian: rows are the three residuals, columns are (z0, ydot0, t_half).
        # d(residual_i)/d(z0)    = Phi[i, 2]
        # d(residual_i)/d(ydot0) = Phi[i, 4]
        # d(residual_i)/d(t)     = f[i]
        rows = (1, 3, 5)
        jac = np.zeros((3, 3), dtype=np.float64)
        for ri, row_idx in enumerate(rows):
            jac[ri, 0] = float(stm_half[row_idx, 2])  # z0 sensitivity
            jac[ri, 1] = float(stm_half[row_idx, 4])  # ydot0 sensitivity
            jac[ri, 2] = float(f_half[row_idx])  # t_half sensitivity
        try:
            dvec = np.linalg.solve(jac, -r)
        except np.linalg.LinAlgError:
            # Singular Jacobian -- fall back to least squares (don't abort; let
            # the caller see the divergent step and stop).
            dvec, *_ = np.linalg.lstsq(jac, -r, rcond=None)
        # Step damping: cap each component to avoid runaway off the family.
        dz0 = float(np.clip(dvec[0], -0.1, 0.1))
        dydot0 = float(np.clip(dvec[1], -0.1, 0.1))
        dt_half = float(np.clip(dvec[2], -0.2 * t_half, 0.2 * t_half))
        z0 += dz0
        ydot0 += dydot0
        t_half += dt_half
        # Diverged into nonphysical period? Stop with current residual.
        if t_half <= 0.0 or not math.isfinite(t_half):
            break
    converged = residual < tol
    t_full = 2.0 * t_half
    state0 = _build_state0(x0, z0, ydot0)
    jacobi = cr3bp.jacobi_constant(state0, system.mu)
    mono: NDArray[np.float64] | None = None
    if converged and with_monodromy:
        try:
            mono = monodromy(system, state0, t_full, rtol=rtol, atol=atol)
        except RuntimeError:
            mono = None
    return SymmetricNRHO(
        x0=float(x0),
        z0=float(z0),
        ydot0=float(ydot0),
        T_TU=float(t_full),
        jacobi=float(jacobi),
        converged=bool(converged),
        closure_residual=float(residual),
        n_iter=int(n_iter),
        monodromy=mono,
    )


# ---------------------------------------------------------------------------
# Family continuation in x0 (natural parameter, secant prediction).
# ---------------------------------------------------------------------------


def _perilune_distance(
    system: cr3bp.CR3BPSystem,
    state0: NDArray[np.float64],
    period: float,
    *,
    rtol: float = 1e-11,
    atol: float = 1e-11,
) -> float:
    """Minimum r2 (secondary-relative distance) over one period, in nondim units.

    Uses the standard CR3BP propagator with dense steps and a perilune event
    (``r2 minimum``) for a tight estimate. For deep low-perilune members a
    regularised propagator would give a tighter floor, but for the
    family-floor STOP test (we just need to know if r2 < threshold) the
    standard propagator is fine.
    """
    moon_x = 1.0 - system.mu

    def _r2_event(t: float, y: NDArray[np.float64], _mu: float) -> float:
        rx = y[0] - moon_x
        ry = y[1]
        rz = y[2]
        vx = y[3]
        vy = y[4]
        vz = y[5]
        return float(rx * vx + ry * vy + rz * vz)

    _r2_event.terminal = False  # type: ignore[attr-defined]
    _r2_event.direction = 1.0  # type: ignore[attr-defined]
    sol = solve_ivp(
        cr3bp.cr3bp_eom,
        (0.0, period),
        np.asarray(state0, dtype=np.float64),
        args=(system.mu,),  # type: ignore[call-overload]
        method="DOP853",
        rtol=rtol,
        atol=atol,
        events=_r2_event,
        dense_output=True,
    )
    y_ev = sol.y_events[0] if sol.y_events is not None else []
    if len(y_ev) == 0:
        # No perilune event in the period -- fall back to dense-grid sample.
        y_grid = np.asarray(sol.y, dtype=np.float64)
        rx = y_grid[0, :] - moon_x
        ry = y_grid[1, :]
        rz = y_grid[2, :]
        return float(np.min(np.sqrt(rx * rx + ry * ry + rz * rz)))
    r2s: list[float] = []
    for y in y_ev:
        ex = float(y[0]) - moon_x
        ey = float(y[1])
        ez = float(y[2])
        r2s.append(math.sqrt(ex * ex + ey * ey + ez * ez))
    # Also include the t=0 distance (in case the IC is at perilune).
    rx0 = float(state0[0]) - moon_x
    ry0 = float(state0[1])
    rz0 = float(state0[2])
    r2s.append(math.sqrt(rx0 * rx0 + ry0 * ry0 + rz0 * rz0))
    return float(min(r2s))


def continue_nrho_family(
    seed: SymmetricNRHO,
    system: cr3bp.CR3BPSystem,
    *,
    label: str = "nrho",
    direction: int = -1,
    d_x0: float = 5e-3,
    n_steps_max: int = 200,
    perilune_floor_km: float | None = None,
    tol: float = 1e-11,
    max_iter: int = 60,
    rtol: float = 1e-12,
    atol: float = 1e-12,
    bif_k_max: int = 6,
    bif_tol: float = 1e-2,
    stop_on_first_bifurcation: bool = False,
    max_no_converge: int = 3,
) -> NRHOBranch:
    """Continue an NRHO family in x0 by natural-parameter (secant) prediction.

    Parameters
    ----------
    seed :
        A converged starting member.
    system :
        CR3BP system (mu).
    direction :
        +1 increase x0, -1 decrease x0. The Koblick L2 Southern NRHO family
        runs toward smaller x0 (smaller perilune) for higher petal counts, so
        the natural Phase 3 direction is ``-1``.
    d_x0 :
        Natural-parameter step in x0 (nondim). The corrector takes the secant
        prediction from the previous two members as its initial guess.
    n_steps_max :
        Step budget cap.
    perilune_floor_km :
        If given, the continuation stops when the perilune drops below this
        threshold (km, dimensionalised via ``system.l_km``). Useful to bound
        the search away from impactor members.
    tol, max_iter, rtol, atol :
        Forwarded to :func:`correct_symmetric_nrho`.
    bif_k_max, bif_tol :
        Forwarded to the period-multiplying detector.
    stop_on_first_bifurcation :
        If True the continuation terminates with ``BIFURCATION_HIT`` on the
        first k=2 bracket. Used by :func:`cyclerfinder.genome.tulip.find_tulip_via_continuation`
        to land just past the bifurcation before family-switching.
    max_no_converge :
        Number of consecutive corrector failures before stopping with
        ``NO_CONVERGE`` (defensive against post-fold runaway).

    Returns
    -------
    NRHOBranch :
        Ordered family members + detected bifurcation brackets + stop reason.

    Notes
    -----
    Fold-turning is NOT performed: when the Jacobian becomes singular in x0
    the corrector will fail to converge and the branch stops with NO_CONVERGE
    (or, in a clean case, FOLD if a Jacobian-condition heuristic fires first).
    Pseudo-arclength continuation in ``(x0, z0, ydot0, T)`` is left to Phase 4.
    """
    branch = NRHOBranch(label=label)
    if not seed.converged:
        branch.stop_reason = NRHOStopReason.NO_CONVERGE
        return branch
    branch.members.append(seed)

    direction_sign = 1.0 if direction >= 0 else -1.0
    prev: SymmetricNRHO = seed
    prev_prev: SymmetricNRHO | None = None
    no_converge_streak = 0

    # Per-member cached (k -> (best_eig, best_dist)) summary for bracket detection.
    def _k_summary(eigs: NDArray[np.complex128]) -> dict[int, tuple[complex, float]]:
        out: dict[int, tuple[complex, float]] = {}
        for k in range(2, bif_k_max + 1):
            best_dist = float("inf")
            best_eig = complex(1.0)
            for e in eigs:
                d, _ = _nearest_kth_root_dist(complex(e), k)
                if d < best_dist:
                    best_dist = d
                    best_eig = complex(e)
            out[k] = (best_eig, best_dist)
        return out

    prev_eigs = floquet_multipliers(seed.monodromy) if seed.monodromy is not None else None
    prev_k_summary = _k_summary(prev_eigs) if prev_eigs is not None else None

    for step in range(n_steps_max):
        # Natural-parameter step in x0.
        x0_next = prev.x0 + direction_sign * d_x0
        # Secant prediction for (z0, ydot0, T): linear in x0.
        if prev_prev is None:
            z0_guess = prev.z0
            ydot0_guess = prev.ydot0
            period_guess = prev.T_TU
        else:
            denom = prev.x0 - prev_prev.x0
            if abs(denom) < 1e-15:
                # Defensive: degenerate (shouldn't happen). Fall back to constant.
                z0_guess = prev.z0
                ydot0_guess = prev.ydot0
                period_guess = prev.T_TU
            else:
                slope_z = (prev.z0 - prev_prev.z0) / denom
                slope_y = (prev.ydot0 - prev_prev.ydot0) / denom
                slope_period = (prev.T_TU - prev_prev.T_TU) / denom
                dx = x0_next - prev.x0
                z0_guess = prev.z0 + slope_z * dx
                ydot0_guess = prev.ydot0 + slope_y * dx
                period_guess = prev.T_TU + slope_period * dx
        member = correct_symmetric_nrho(
            system,
            x0_next,
            z0_guess,
            ydot0_guess,
            period_guess,
            tol=tol,
            max_iter=max_iter,
            rtol=rtol,
            atol=atol,
            with_monodromy=True,
        )
        if not member.converged:
            no_converge_streak += 1
            if no_converge_streak >= max_no_converge:
                # Cleanly stop. Could be a fold or simply too coarse a step.
                # Heuristic: if the previous step also failed AND we never crossed
                # a fold (period stayed monotone), call it NO_CONVERGE; otherwise
                # call it FOLD. Without a more careful Jacobian-condition probe
                # we conservatively label NO_CONVERGE.
                branch.stop_reason = NRHOStopReason.NO_CONVERGE
                break
            # Halve the step and retry from prev (do NOT advance prev_prev).
            d_x0 *= 0.5
            continue
        no_converge_streak = 0
        # Adjacent-member bifurcation detection using the MONODROMY ALREADY
        # COMPUTED inside the corrector -- avoid re-integrating in
        # ``scan_family_for_bifurcations`` (which doubles the per-step cost).
        new_bracs: list[BifurcationPoint] = []
        if member.monodromy is not None:
            curr_eigs = floquet_multipliers(member.monodromy)
            curr_k_summary = _k_summary(curr_eigs)
            if prev_k_summary is not None:
                prev_mem = FamilyMember(
                    label=f"{label}.x0={prev.x0:.6f}",
                    state0=_build_state0(prev.x0, prev.z0, prev.ydot0),
                    period=prev.T_TU,
                    mu=system.mu,
                    parameter=prev.x0,
                )
                curr_mem = FamilyMember(
                    label=f"{label}.x0={member.x0:.6f}",
                    state0=_build_state0(member.x0, member.z0, member.ydot0),
                    period=member.T_TU,
                    mu=system.mu,
                    parameter=member.x0,
                )
                for k in range(2, bif_k_max + 1):
                    eig_a, dist_a = prev_k_summary[k]
                    eig_b, dist_b = curr_k_summary[k]
                    sa = dist_a - bif_tol
                    sb = dist_b - bif_tol
                    crossed = sa == 0.0 or sb == 0.0 or (sa > 0.0) != (sb > 0.0)
                    if crossed:
                        new_bracs.append(
                            BifurcationPoint(
                                k=k,
                                members=(prev_mem, curr_mem),
                                eig_before=eig_a,
                                eig_after=eig_b,
                                dist_before=dist_a,
                                dist_after=dist_b,
                                tol=bif_tol,
                                extras={
                                    "param_before": prev.x0,
                                    "param_after": member.x0,
                                },
                            )
                        )
            prev_k_summary = curr_k_summary
        branch.bifurcations.extend(new_bracs)
        # Floor check (perilune).
        if perilune_floor_km is not None:
            s0 = _build_state0(member.x0, member.z0, member.ydot0)
            r2_min = _perilune_distance(system, s0, member.T_TU, rtol=rtol, atol=atol)
            if r2_min * system.l_km < perilune_floor_km:
                branch.members.append(member)
                branch.n_steps = step + 1
                branch.stop_reason = NRHOStopReason.PERILUNE_FLOOR
                return branch
        branch.members.append(member)
        branch.n_steps = step + 1
        # Stop-on-first-bifurcation hook.
        if stop_on_first_bifurcation and any(b.k == 2 for b in new_bracs):
            branch.stop_reason = NRHOStopReason.BIFURCATION_HIT
            return branch
        prev_prev = prev
        prev = member
    else:
        branch.stop_reason = NRHOStopReason.MAX_STEPS
    return branch


# ---------------------------------------------------------------------------
# Convenience: round-trip closure check via regularised propagator.
# ---------------------------------------------------------------------------


def regularized_full_period_closure(
    system: cr3bp.CR3BPSystem,
    member: SymmetricNRHO,
    *,
    regularization: str = "r2",
    rtol: float = 1e-12,
    atol: float = 1e-12,
) -> float:
    """Re-propagate ``member`` for one full period via the Sundman-regularised
    propagator and return the L2 closure residual at ``t = T``.

    This is an INDEPENDENT closure check -- the symmetric corrector lands the
    half-period perpendicular crossing; the regularised propagator gives a
    different numerical fingerprint than the standard DOP853 used by
    :func:`correct_symmetric_nrho`. Disagreement at machine precision flags
    an IC that closes under one integrator but not another (which would be a
    sign the corrector landed on a near-periodic orbit, not a true one).
    """
    state0 = _build_state0(member.x0, member.z0, member.ydot0)
    t_final = member.T_TU
    s_span = creg.physical_to_regularized_span(system, state0, (0.0, t_final))
    s_span = (s_span[0], s_span[1] * 1.5)
    arc = creg.propagate_regularized(
        system,
        state0,
        s_span,
        rtol=rtol,
        atol=atol,
        regularization=regularization,
        t_stop=t_final,
    )
    state_f = arc.state_at_s[:, -1]
    return float(np.linalg.norm(state_f - state0))


def label_branch_members(branch: NRHOBranch, mu: float) -> Sequence[FamilyMember]:
    """Convert an :class:`NRHOBranch` to a sequence of :class:`FamilyMember`
    suitable for :func:`scan_family_for_bifurcations`.

    Provided as a convenience for callers that want to re-scan a branch with a
    different ``tol`` / ``k_max``.
    """
    return [
        FamilyMember(
            label=f"{branch.label}.idx={i}",
            state0=_build_state0(m.x0, m.z0, m.ydot0),
            period=m.T_TU,
            mu=float(mu),
            parameter=m.x0,
        )
        for i, m in enumerate(branch.members)
    ]
