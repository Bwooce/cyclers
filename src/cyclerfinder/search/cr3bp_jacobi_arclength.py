"""Pseudo-arclength continuation of CR3BP symmetric (k1,k2) cyclers in the
Jacobi constant C at FIXED mass parameter mu (fold-turning).

Motivation (#249, Braik-Ross reachable-set scorer gate)
-------------------------------------------------------
At a common Jacobi constant the network's symmetric cyclers split across a
saddle-center FOLD in C: the recoverable stable branch (e.g. C11b ~ 56 d) and
the unrecovered members (C11a ~ 42 d, C21 ~ 85 d) sit on opposite sides of a
turning point ``dx0/dC -> inf``. Two existing tools cannot cross it:

  - the 1-DOF perpendicular-x-crossing symmetric corrector
    (:func:`cyclerfinder.search.cr3bp_periodic.correct_symmetric_fixed_jacobi`)
    lands only on the stable branch;
  - natural-parameter Jacobi continuation (``continue_family``: fix C, re-solve
    x0) DIVERGES at the fold.

This module follows the family curve by PSEUDO-ARCLENGTH continuation -- the same
fold-turning scheme as :func:`cyclerfinder.search.mu_continuation.continue_in_mu`
but with mu FROZEN, so the free vector is ``z = (x0, C)`` (one fewer unknown).

What is continued
-----------------
A perpendicular-x-axis-crossing symmetric cycler has IC ``(x0, 0, 0, 0, ydot0,
0)`` with ``ydot0`` fixed algebraically from the Jacobi constant ``C`` (Ross
Eq. 9). For a given ``k1,k2`` (fixed by the crossing index + velocity sign) the
single periodicity condition is ``xdot(t_half) = 0`` at the family's half-period
crossing (Ross Eq. 11). With mu held fixed the solution set of

    r(x0, C) = xdot(t_half; x0, C, mu) = 0

is a 1-D *family curve* in the ``(x0, C)`` plane. We follow it by predicting
along the local tangent (null vector of ``dr/dz``), then correcting back onto
``{r = 0, tangent . (z - z_pred) = 0}``. The arclength constraint replaces the
fixed-C condition that natural-parameter continuation uses, so the walk turns
the fold (where ``dx0/dC`` blows up but ``ds`` stays finite). To LAND on a
requested ``C_target`` after turning the fold, :func:`land_at_jacobi` then
Newton-projects ``x0`` onto ``r(x0; C_target, mu) = 0`` at the held C.

This is the 2-variable specialisation of ``mu_continuation``'s 3-variable
``z = (x0, C, mu)`` scheme. The small residual / tangent / corrector helpers are
replicated here in their 2-var form (rather than importing the 3-var ones) so
this module is self-contained and the frozen-mu math is explicit; the
member-building reuses ``MuMember`` / ``crosscheck_periodic`` directly.

Discipline (orbit-closure):
  - every kept member is a genuine periodic orbit: ``|xdot(t_half)| < tol`` AND
    an independent-Radau full-period re-closure (different integrator than the
    DOP853 corrector) with bounded Jacobi drift;
  - stability ``nu = 1/2 (lambda + 1/lambda)`` from the Barden half-period
    monodromy (Ross Eqs. 13-15); ``|nu| < 1`` is linearly stable;
  - the recovered ``(C, T, IC, nu)`` are OUR OWN computed values -- DISCOVERIES,
    not sourced rows; NO catalogue writeback.

Model: pure planar CR3BP (PCR3BP). Pure (math/numpy/scipy + the CR3BP core).
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass, field
from enum import StrEnum

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import solve_ivp

import cyclerfinder.core.cr3bp as cr3bp
import cyclerfinder.search.cr3bp_periodic as cp
from cyclerfinder.search.mu_continuation import MuMember

# Earth-Moon scales used to build the CR3BPSystem for the Radau cross-check; the
# pure dynamics depend only on mu (these set the dimensional reporting frame).
_EM_L_KM = 384400.0
_EM_T_S = 375699.8


class JacobiStopReason(StrEnum):
    """Why a Jacobi (fixed-mu) arclength branch terminated."""

    TARGET_REACHED = "target_reached"  # landed on the requested C_target
    MAX_STEPS = "max_steps"  # step budget exhausted (not a physical edge)
    STEP_UNDERFLOW = "step_underflow"  # arclength step halved below floor (edge)
    NO_MEMBER = "no_member"  # the symmetric IC lost its k-th crossing (topology loss)
    TOPOLOGY_JUMP = "topology_jump"  # period jumped vs previous member (branch switch)


@dataclass
class JacobiBranch:
    """The ordered branch produced by continuing one seed in C at fixed mu."""

    label: str
    mu: float
    half_crossings: int
    ydot0_sign: float
    c_start: float
    c_target: float
    members: list[MuMember] = field(default_factory=list)
    stop_reason: JacobiStopReason = JacobiStopReason.MAX_STEPS
    n_steps: int = 0


# --------------------------------------------------------------------------- #
# 2-variable (x0, C) helpers at fixed mu.
# These are the frozen-mu specialisations of mu_continuation's 3-var ``_ydot0``,
# ``_half_crossing``, ``_residual_jac``, ``_nu_at``, ``_correct``, ``_tangent``.
# --------------------------------------------------------------------------- #


def _ydot0(x0: float, c: float, mu: float, sign: float) -> float | None:
    rad = cp._ubar_x_at_axis(x0, mu) - c
    if rad < 0.0:
        return None
    return float(sign) * float(np.sqrt(rad))


def _half_crossing(
    x0: float,
    c: float,
    mu: float,
    sign: float,
    hc: int,
    t_hi: float,
    *,
    with_stm: bool,
    rtol: float,
    atol: float,
) -> tuple[float, NDArray[np.float64]] | None:
    """Integrate the symmetric IC and return (t, state[+stm]) at the hc-th y=0 crossing."""
    ydot0 = _ydot0(x0, c, mu, sign)
    if ydot0 is None:
        return None
    st = np.array([x0, 0.0, 0.0, 0.0, ydot0, 0.0])

    def _ev(t: float, y: NDArray[np.float64], _mu: float) -> float:
        return float(y[1])

    _ev.direction = 0.0  # type: ignore[attr-defined]
    if with_stm:
        y0 = np.concatenate([st, np.eye(6).reshape(36)])
        rhs: Callable[[float, NDArray[np.float64], float], NDArray[np.float64]] = (
            cr3bp.cr3bp_stm_eom
        )
    else:
        y0 = st
        rhs = cr3bp.cr3bp_eom
    sol = solve_ivp(
        rhs,
        (0.0, t_hi),
        y0,
        args=(mu,),  # type: ignore[call-overload]
        method="DOP853",
        rtol=rtol,
        atol=atol,
        events=_ev,
    )
    pairs = [
        (t, y) for t, y in zip(sol.t_events[0], sol.y_events[0], strict=True) if t > 1e-6 * t_hi
    ]
    if len(pairs) < hc:
        return None
    t_half, yf = pairs[hc - 1]
    return float(t_half), yf


def _residual_jac(
    z: NDArray[np.float64],
    mu: float,
    sign: float,
    hc: int,
    t_hi: float,
    *,
    rtol: float,
    atol: float,
) -> tuple[float, float, NDArray[np.float64]] | None:
    """Return (r0=xdot(t_half), t_half, dr/dz) for z=(x0,C) by finite differences (no STM)."""
    base = _half_crossing(z[0], z[1], mu, sign, hc, t_hi, with_stm=False, rtol=rtol, atol=atol)
    if base is None:
        return None
    t_half, yf = base
    r0 = float(yf[3])
    grad = np.zeros(2)
    h = (1e-7, 1e-8)
    for i in range(2):
        zp = z.copy()
        zp[i] += h[i]
        rp = _half_crossing(zp[0], zp[1], mu, sign, hc, t_hi, with_stm=False, rtol=rtol, atol=atol)
        if rp is None:
            return None
        grad[i] = (float(rp[1][3]) - r0) / h[i]
    return r0, t_half, grad


def _nu_at(
    z: NDArray[np.float64],
    mu: float,
    sign: float,
    hc: int,
    t_hi: float,
    *,
    rtol: float,
    atol: float,
) -> tuple[float, complex, float] | None:
    """Barden half-period monodromy nu at the member; returns (nu, lambda, period)."""
    m = _half_crossing(z[0], z[1], mu, sign, hc, t_hi, with_stm=True, rtol=rtol, atol=atol)
    if m is None:
        return None
    t_half, yf = m
    phi = yf[6:].reshape(6, 6)
    idx = [0, 1, 3, 4]
    phi4 = phi[np.ix_(idx, idx)]
    g4 = np.diag([1.0, -1.0, -1.0, 1.0])
    monodromy = g4 @ np.linalg.inv(phi4) @ g4 @ phi4
    eigs = np.linalg.eigvals(monodromy)
    order = np.argsort(np.abs(eigs - 1.0))
    nontrivial = eigs[order[2:]]
    lam = complex(nontrivial[np.argmax(np.abs(nontrivial))])
    nu = float((0.5 * (lam + 1.0 / lam)).real)
    return nu, lam, 2.0 * t_half


def tangent(
    z: NDArray[np.float64],
    mu: float,
    sign: float,
    hc: int,
    t_hi: float,
    *,
    prev: NDArray[np.float64] | None,
    rtol: float = 1e-12,
    atol: float = 1e-12,
) -> NDArray[np.float64] | None:
    """Unit tangent = null vector of dr/dz (z=(x0,C)), oriented for continuity.

    The family curve in the ``(x0, C)`` plane has tangent direction in the null
    space of the 1x2 Jacobian ``dr/dz``. Returned as a unit vector; oriented to
    keep a positive dot product with ``prev`` (or, at the start, toward
    increasing C).
    """
    rj = _residual_jac(z, mu, sign, hc, t_hi, rtol=rtol, atol=atol)
    if rj is None:
        return None
    _r0, _t_half, grad = rj
    _, _, vt = np.linalg.svd(grad.reshape(1, 2))
    tan = np.asarray(vt[-1], dtype=np.float64)
    if prev is None:
        if tan[1] < 0:  # orient toward increasing C at the start
            tan = -tan
    elif tan @ prev < 0:
        tan = -tan
    return tan


def _correct(
    z_pred: NDArray[np.float64],
    tan: NDArray[np.float64],
    mu: float,
    sign: float,
    hc: int,
    t_hi: float,
    *,
    tol: float,
    max_iter: int,
    rtol: float,
    atol: float,
    dx0_cap: float = 0.05,
    dc_cap: float = 0.05,
) -> NDArray[np.float64] | None:
    """Newton onto {r(z)=0, tangent.(z-z_pred)=0}; z=(x0,C) (2 eqns, 2 unknowns)."""
    z = z_pred.copy()
    for _ in range(max_iter):
        rj = _residual_jac(z, mu, sign, hc, t_hi, rtol=rtol, atol=atol)
        if rj is None:
            return None
        r0, _t_half, grad = rj
        arc = float(tan @ (z - z_pred))
        if abs(r0) < tol and abs(arc) < 1e-12:
            return z
        jmat = np.vstack([grad, tan])  # 2 x 2
        try:
            dz = np.linalg.solve(jmat, -np.array([r0, arc]))
        except np.linalg.LinAlgError:
            dz = np.linalg.lstsq(jmat, -np.array([r0, arc]), rcond=None)[0]
        dz = np.asarray(dz, dtype=np.float64)
        dz[0] = float(np.clip(dz[0], -dx0_cap, dx0_cap))
        dz[1] = float(np.clip(dz[1], -dc_cap, dc_cap))
        z = z + dz
    return None


def land_at_jacobi(
    z_pred: NDArray[np.float64],
    c_target: float,
    mu: float,
    sign: float,
    hc: int,
    t_hi: float,
    *,
    tol: float,
    max_iter: int,
    rtol: float = 1e-12,
    atol: float = 1e-12,
    dx0_cap: float = 0.05,
) -> NDArray[np.float64] | None:
    """Land exactly on ``C = c_target`` by projecting ``x0`` onto ``r=0`` at the held C.

    Holds ``C = c_target`` and Newton-projects ``x0`` onto the periodicity curve
    ``r(x0, c_target) = 0`` (1 unknown, 1 equation). This is the fixed-C analogue
    of ``mu_continuation._land_at_mu`` (which holds mu) -- it gives a unique
    nearby member exactly at the requested Jacobi constant after the arclength
    walk has carried the curve (possibly through a fold) to ``c_target``.

    Returns the converged ``z = (x0, c_target)`` or ``None`` on failure (the
    crossing topology was lost, or ``dr/dx0`` vanished -- i.e. ``c_target`` sits
    at a fold where x0 is not single-valued in C).
    """
    z = z_pred.copy()
    z[1] = float(c_target)
    for _ in range(max_iter):
        rj = _residual_jac(z, mu, sign, hc, t_hi, rtol=rtol, atol=atol)
        if rj is None:
            return None
        r0, _t_half, grad = rj
        if abs(r0) < tol:
            return z
        gx = float(grad[0])  # only x0 moves
        if abs(gx) < 1e-30:
            return None
        d0 = -r0 / gx
        d0 = float(np.clip(d0, -dx0_cap, dx0_cap))
        z = z + np.array([d0, 0.0])
    return None


def _make_member(
    z: NDArray[np.float64],
    mu: float,
    sign: float,
    hc: int,
    t_hi: float,
    *,
    rtol: float,
    atol: float,
    radau_closure_tol: float,
    radau_jacobi_tol: float,
) -> MuMember | None:
    """Build a fully-characterised member (nu + independent-Radau cross-check)."""
    nr = _nu_at(z, mu, sign, hc, t_hi, rtol=rtol, atol=atol)
    if nr is None:
        return None
    nu, lam, period = nr
    ydot0 = _ydot0(z[0], z[1], mu, sign)
    if ydot0 is None:
        return None
    state0 = np.array([z[0], 0.0, 0.0, 0.0, ydot0, 0.0])
    system = cr3bp.CR3BPSystem(
        mu=float(mu), primary="Earth", secondary="Moon", l_km=_EM_L_KM, t_s=_EM_T_S
    )
    jac = cr3bp.jacobi_constant(state0, system.mu)
    rj = _half_crossing(z[0], z[1], mu, sign, hc, t_hi, with_stm=False, rtol=rtol, atol=atol)
    cross_res = abs(float(rj[1][3])) if rj is not None else float("nan")
    po = cp.PeriodicOrbit(
        state0=state0, period=period, jacobi=jac, converged=True, closure_residual=cross_res
    )
    _radau_ok, radau_dj = cp.crosscheck_periodic(
        system, po, closure_tol=radau_closure_tol, jacobi_tol=radau_jacobi_tol
    )
    return MuMember(
        mu=float(mu),
        state0=state0,
        x0=float(z[0]),
        ydot0=float(ydot0),
        jacobi=float(jac),
        period=float(period),
        nu=float(nu),
        abs_lambda=float(abs(lam)),
        crossing_residual=float(cross_res),
        radau_djacobi=float(radau_dj),
        stable=abs(float(nu)) < 1.0,
    )


def continue_in_jacobi(
    seed: cp.SymmetricOrbit,
    *,
    mu: float,
    half_crossings: int,
    ydot0_sign: float,
    c_target: float,
    label: str = "",
    ds0: float = 8e-3,
    ds_max: float = 4e-2,
    ds_min: float = 1e-6,
    max_steps: int = 5000,
    corrector_tol: float = 1e-11,
    corrector_max_iter: int = 80,
    t_hi_frac: float = 1.8,
    period_jump_frac: float = 0.0,
    radau_closure_tol: float = 1e-3,
    radau_jacobi_tol: float = 1e-8,
    rtol: float = 1e-12,
    atol: float = 1e-12,
    record_every: int = 1,
) -> JacobiBranch:
    """Pseudo-arclength continuation of a symmetric cycler in C at FIXED mu.

    Walks the ``(x0, C)`` family curve from ``seed`` toward ``c_target``,
    turning folds (where natural-parameter Jacobi continuation diverges). When a
    full arclength step would carry C past ``c_target`` the walk lands EXACTLY on
    ``c_target`` via :func:`land_at_jacobi`.

    Parameters
    ----------
    seed:
        A converged :class:`~cyclerfinder.search.cr3bp_periodic.SymmetricOrbit`
        produced at ``mu`` (e.g. from ``correct_symmetric_fixed_jacobi``).
    mu:
        The (frozen) mass parameter. Held fixed along the entire continuation.
    half_crossings, ydot0_sign:
        Identify the (k1,k2) branch (the half-period crossing index and the sign
        of ``ydot0``); held fixed along the continuation.
    c_target:
        The Jacobi constant to continue toward (possibly across a fold).
    ds0, ds_max, ds_min:
        Initial / maximum / minimum arclength step (adaptive: grows on success,
        halves on corrector failure; underflow below ``ds_min`` stops the walk).
    period_jump_frac:
        If > 0, reject (and shrink the step on) any corrected member whose period
        differs from the previous member by more than this fraction -- a
        branch-switch / topology-jump guard that keeps the same crossing-index
        family. ``0.0`` (default) disables it (the bare arclength follows the
        branch through folds, which can switch the crossing topology).

    Returns
    -------
    JacobiBranch
        Ordered kept members (seed first), stop reason and step count.

    Notes
    -----
    Because the curve may FOLD in C, ``c_target`` is not monotone in arclength:
    the walk tracks the SIGNED progress of C toward the target and only attempts
    the fixed-C landing when the local tangent points toward the target. After
    turning a fold the second (e.g. unstable) branch is reached, and the landing
    projection then converges on the target-C member of THAT branch.
    """
    mu = float(mu)
    z = np.array([seed.x0, seed.jacobi])
    per = seed.period

    branch = JacobiBranch(
        label=label,
        mu=mu,
        half_crossings=half_crossings,
        ydot0_sign=ydot0_sign,
        c_start=float(seed.jacobi),
        c_target=float(c_target),
    )

    t_hi = t_hi_frac * per
    first = _make_member(
        z,
        mu,
        ydot0_sign,
        half_crossings,
        t_hi,
        rtol=rtol,
        atol=atol,
        radau_closure_tol=radau_closure_tol,
        radau_jacobi_tol=radau_jacobi_tol,
    )
    if first is None:
        branch.stop_reason = JacobiStopReason.NO_MEMBER
        return branch
    branch.members.append(first)

    tan = tangent(z, mu, ydot0_sign, half_crossings, t_hi, prev=None, rtol=rtol, atol=atol)
    if tan is None:
        branch.stop_reason = JacobiStopReason.NO_MEMBER
        return branch

    ds = ds0
    for _step in range(max_steps):
        if abs(z[1] - c_target) < 1e-12:
            branch.stop_reason = JacobiStopReason.TARGET_REACHED
            break
        # If a full step (in the current tangent direction) would carry C across
        # c_target, land EXACTLY on c_target by the fixed-C projection.
        if tan[1] != 0.0:
            ds_to_target = (c_target - z[1]) / tan[1]
            if 0.0 < ds_to_target <= ds:
                z_pred = z + ds_to_target * tan
                zland = land_at_jacobi(
                    z_pred,
                    c_target,
                    mu,
                    ydot0_sign,
                    half_crossings,
                    t_hi_frac * max(per, seed.period),
                    tol=corrector_tol,
                    max_iter=corrector_max_iter,
                    rtol=rtol,
                    atol=atol,
                )
                if zland is not None:
                    z = zland
                    branch.n_steps += 1
                    branch.stop_reason = JacobiStopReason.TARGET_REACHED
                    break
        z_pred = z + ds * tan
        t_hi = t_hi_frac * max(per, seed.period)
        zc = _correct(
            z_pred,
            tan,
            mu,
            ydot0_sign,
            half_crossings,
            t_hi,
            tol=corrector_tol,
            max_iter=corrector_max_iter,
            rtol=rtol,
            atol=atol,
        )
        if zc is None:
            ds *= 0.5
            if ds < ds_min:
                branch.stop_reason = JacobiStopReason.STEP_UNDERFLOW
                break
            continue
        ntan = tangent(zc, mu, ydot0_sign, half_crossings, t_hi, prev=tan, rtol=rtol, atol=atol)
        if ntan is None:
            ds *= 0.5
            if ds < ds_min:
                branch.stop_reason = JacobiStopReason.STEP_UNDERFLOW
                break
            continue
        # Topology / branch-switch guard: a corrected member whose period jumps
        # vs the previous member is the corrector landing on a different crossing
        # branch. Reject and shrink the step rather than silently switch families.
        cheap = _half_crossing(
            zc[0], zc[1], mu, ydot0_sign, half_crossings, t_hi, with_stm=False, rtol=rtol, atol=atol
        )
        new_per = 2.0 * cheap[0] if cheap is not None else per
        if period_jump_frac > 0.0 and abs(new_per - per) > period_jump_frac * per:
            ds *= 0.5
            if ds < ds_min:
                branch.stop_reason = JacobiStopReason.TOPOLOGY_JUMP
                break
            continue
        z, tan, per = zc, ntan, new_per
        branch.n_steps += 1
        ds = min(ds * 1.3, ds_max)
        if branch.n_steps % record_every == 0:
            mem = _make_member(
                z,
                mu,
                ydot0_sign,
                half_crossings,
                t_hi_frac * per,
                rtol=rtol,
                atol=atol,
                radau_closure_tol=radau_closure_tol,
                radau_jacobi_tol=radau_jacobi_tol,
            )
            if mem is not None:
                branch.members.append(mem)
    else:
        branch.stop_reason = JacobiStopReason.MAX_STEPS

    # Always record the final landed member (even if not on a record_every step).
    final = _make_member(
        z,
        mu,
        ydot0_sign,
        half_crossings,
        t_hi_frac * per,
        rtol=rtol,
        atol=atol,
        radau_closure_tol=radau_closure_tol,
        radau_jacobi_tol=radau_jacobi_tol,
    )
    if final is not None and (
        not branch.members or abs(branch.members[-1].jacobi - final.jacobi) > 1e-12
    ):
        branch.members.append(final)
    return branch
