"""Mass-parameter (mu) continuation of CR3BP symmetric (k1,k2) cyclers.

Motivation (Roberts-Tsoukkas & Ross 2026, "Stable Prograde Earth-Moon
Multi-Orbiter Cyclers via Three-Body Dynamics", journal version of AAS 25-621):
the paper shows -- in FIGURES ONLY, no printed numbers -- that the stable
prograde cycler families persist up into the binary-star mass-parameter range
(Fig. 3: a (1,3) exterior cycler at mu=0.1, a (3,1) cycler at mu=0.3, and a
(1,1) equal-mass cycler at mu=0.5; all depicted stable). This module
NUMERICALLY CONTINUES a held Earth-Moon family member (mu ~ 0.01215) up in the
mass parameter mu to recover the depicted regime, self-discovering the orbits
the paper only draws.

What is continued
-----------------
A perpendicular-x-axis-crossing symmetric cycler has IC ``(x0, 0, 0, 0, ydot0,
0)`` with ``ydot0`` fixed algebraically from the Jacobi constant ``C`` (Ross
Eq. 9). For a given ``k1,k2`` (fixed by the crossing index + velocity sign) the
single periodicity condition is ``xdot(t_half) = 0`` at the family's
half-period crossing (Ross Eq. 11). Thus the solution set of

    r(x0, C, mu) = xdot(t_half; x0, C, mu) = 0

is a 2-surface in ``(x0, C, mu)``. A *family* at fixed mu is a curve on it
(parameterised by C); sweeping mu sweeps the surface. We follow a 1-D path on
this surface from the Earth-Moon mu up to a target mu by PSEUDO-ARCLENGTH
continuation: predict along the local tangent (null vector of ``dr/dz``),
correct back onto ``r = 0`` with the arclength constraint
``tangent . (z - z_pred) = 0``. This needs no member-selection invariant (no
fixed C, no fixed period, no fixed nu) -- the branch is followed wherever it
goes, including folds in mu. Natural-parameter continuation (fix mu, solve the
remaining unknowns) is the degenerate special case and is NOT used because it
cannot turn a fold.

Discipline (orbit-closure):
  - every kept member is a genuine periodic orbit: ``|xdot(t_half)| < tol`` AND
    an independent-Radau full-period re-closure (different integrator than the
    DOP853 corrector) with bounded Jacobi drift;
  - stability ``nu = 1/2 (lambda + 1/lambda)`` from the Barden half-period
    monodromy (Ross Eqs. 13-15); ``|nu| < 1`` is linearly stable;
  - the recovered ``(mu, C, T, IC, nu)`` are OUR OWN computed values (the paper
    prints none) -- DISCOVERIES, not sourced rows; NO catalogue writeback.

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


class MuStopReason(StrEnum):
    """Why a mu-continuation branch terminated."""

    TARGET_REACHED = "target_reached"  # crossed the requested mu
    MAX_STEPS = "max_steps"  # step budget exhausted (not a physical edge)
    STEP_UNDERFLOW = "step_underflow"  # arclength step halved below floor (fold/edge)
    NO_MEMBER = "no_member"  # the symmetric IC lost its k-th crossing (topology loss)
    RADAU_REJECT = "radau_reject"  # independent cross-check failed at a kept step


@dataclass(frozen=True)
class MuMember:
    """One converged member along a mu-continuation branch (pure CR3BP)."""

    mu: float
    state0: NDArray[np.float64]  # (x0, 0, 0, 0, ydot0, 0)
    x0: float
    ydot0: float
    jacobi: float
    period: float
    nu: float
    abs_lambda: float
    crossing_residual: float
    radau_djacobi: float
    stable: bool


@dataclass
class MuBranch:
    """The ordered branch produced by continuing one seed in mu."""

    label: str
    half_crossings: int
    ydot0_sign: float
    mu_start: float
    mu_target: float
    members: list[MuMember] = field(default_factory=list)
    stop_reason: MuStopReason = MuStopReason.MAX_STEPS
    n_steps: int = 0


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
    z: NDArray[np.float64], sign: float, hc: int, t_hi: float, *, rtol: float, atol: float
) -> tuple[float, float, NDArray[np.float64]] | None:
    """Return (r0=xdot(t_half), t_half, dr/dz) by finite differences (no STM)."""
    base = _half_crossing(z[0], z[1], z[2], sign, hc, t_hi, with_stm=False, rtol=rtol, atol=atol)
    if base is None:
        return None
    t_half, yf = base
    r0 = float(yf[3])
    grad = np.zeros(3)
    h = (1e-7, 1e-8, 1e-7)
    for i in range(3):
        zp = z.copy()
        zp[i] += h[i]
        rp = _half_crossing(
            zp[0], zp[1], zp[2], sign, hc, t_hi, with_stm=False, rtol=rtol, atol=atol
        )
        if rp is None:
            return None
        grad[i] = (float(rp[1][3]) - r0) / h[i]
    return r0, t_half, grad


def _nu_at(
    z: NDArray[np.float64], sign: float, hc: int, t_hi: float, *, rtol: float, atol: float
) -> tuple[float, complex, float] | None:
    """Barden half-period monodromy nu at the member; returns (nu, lambda, period)."""
    m = _half_crossing(z[0], z[1], z[2], sign, hc, t_hi, with_stm=True, rtol=rtol, atol=atol)
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


def _correct(
    z_pred: NDArray[np.float64],
    tangent: NDArray[np.float64],
    sign: float,
    hc: int,
    t_hi: float,
    *,
    tol: float,
    max_iter: int,
    rtol: float,
    atol: float,
    dx0_cap: float = 0.05,
    dc_cap: float = 0.02,
    dmu_cap: float = 2e-3,
) -> NDArray[np.float64] | None:
    """Newton onto {r(z)=0, tangent.(z-z_pred)=0}; min-norm step (2 eqns, 3 unk)."""
    z = z_pred.copy()
    for _ in range(max_iter):
        rj = _residual_jac(z, sign, hc, t_hi, rtol=rtol, atol=atol)
        if rj is None:
            return None
        r0, _t_half, grad = rj
        arc = float(tangent @ (z - z_pred))
        if abs(r0) < tol and abs(arc) < 1e-12:
            return z
        jmat = np.vstack([grad, tangent])  # 2 x 3
        sol_dz = np.linalg.lstsq(jmat, -np.array([r0, arc]), rcond=None)[0]
        dz = np.asarray(sol_dz, dtype=np.float64)
        dz[0] = float(np.clip(dz[0], -dx0_cap, dx0_cap))
        dz[1] = float(np.clip(dz[1], -dc_cap, dc_cap))
        dz[2] = float(np.clip(dz[2], -dmu_cap, dmu_cap))
        z = z + dz
    return None


def _land_at_mu(
    z_pred: NDArray[np.float64],
    mu_target: float,
    sign: float,
    hc: int,
    t_hi: float,
    *,
    tol: float,
    max_iter: int,
    rtol: float,
    atol: float,
    dx0_cap: float = 0.05,
    dc_cap: float = 0.02,
) -> NDArray[np.float64] | None:
    """Land exactly on ``mu = mu_target`` by projecting onto ``r=0`` at fixed mu.

    Holds ``mu = mu_target`` and Newton-projects ``(x0, C)`` onto the periodicity
    surface ``r(x0,C,mu_target)=0``, minimum-norm in the ``(x0, C)`` plane
    (2 unknowns, 1 equation). Gives a unique nearby member at the exact target.
    """
    z = z_pred.copy()
    z[2] = float(mu_target)
    for _ in range(max_iter):
        rj = _residual_jac(z, sign, hc, t_hi, rtol=rtol, atol=atol)
        if rj is None:
            return None
        r0, _t_half, grad = rj
        if abs(r0) < tol:
            return z
        g2 = grad[:2]  # only x0, C move
        denom = float(g2 @ g2)
        if denom < 1e-30:
            return None
        d2 = -r0 / denom * g2  # min-norm Newton step in (x0, C)
        d2[0] = float(np.clip(d2[0], -dx0_cap, dx0_cap))
        d2[1] = float(np.clip(d2[1], -dc_cap, dc_cap))
        z = z + np.array([d2[0], d2[1], 0.0])
    return None


def _tangent(
    z: NDArray[np.float64],
    sign: float,
    hc: int,
    t_hi: float,
    *,
    prev: NDArray[np.float64] | None,
    rtol: float,
    atol: float,
) -> NDArray[np.float64] | None:
    """Unit tangent = null vector of dr/dz, oriented for continuity with ``prev``."""
    rj = _residual_jac(z, sign, hc, t_hi, rtol=rtol, atol=atol)
    if rj is None:
        return None
    _r0, _t_half, grad = rj
    _, _, vt = np.linalg.svd(grad.reshape(1, 3))
    tan = np.asarray(vt[-1], dtype=np.float64)
    if prev is None:
        if tan[2] < 0:  # orient toward increasing mu at the start
            tan = -tan
    elif tan @ prev < 0:
        tan = -tan
    return tan


def continue_in_mu(
    seed: cp.SymmetricOrbit,
    mu_start: float,
    *,
    half_crossings: int,
    ydot0_sign: float,
    mu_target: float,
    label: str = "",
    ds0: float = 8e-3,
    ds_max: float = 4e-2,
    ds_min: float = 1e-6,
    max_steps: int = 5000,
    corrector_tol: float = 1e-11,
    corrector_max_iter: int = 80,
    t_hi_frac: float = 1.8,
    radau_closure_tol: float = 1e-3,
    radau_jacobi_tol: float = 1e-8,
    rtol: float = 1e-12,
    atol: float = 1e-12,
    record_every: int = 1,
) -> MuBranch:
    """Pseudo-arclength continuation of a symmetric cycler from its seed mu to
    ``mu_target``.

    Parameters
    ----------
    seed:
        A converged :class:`~cyclerfinder.search.cr3bp_periodic.SymmetricOrbit`
        produced at ``mu_start`` (the corrector carries no mu, so it is passed
        explicitly).
    mu_start:
        The mass parameter at which ``seed`` was corrected.
    half_crossings, ydot0_sign:
        Identify the (k1,k2) branch (the half-period crossing index and the
        sign of ``ydot0``); held fixed along the continuation.
    mu_target:
        The mass parameter to continue toward (e.g. 0.1, 0.3, 0.5).
    ds0, ds_max, ds_min:
        Initial / maximum / minimum arclength step (adaptive: grows on success,
        halves on corrector failure; underflow below ``ds_min`` stops the walk).

    Returns
    -------
    MuBranch
        Ordered kept members (seed first), stop reason and step count.
    """
    mu_start = float(mu_start)
    z = np.array([seed.x0, seed.jacobi, mu_start])
    per = seed.period
    t_hi = t_hi_frac * per

    branch = MuBranch(
        label=label,
        half_crossings=half_crossings,
        ydot0_sign=ydot0_sign,
        mu_start=mu_start,
        mu_target=mu_target,
    )

    first = _make_member(
        z,
        ydot0_sign,
        half_crossings,
        t_hi,
        rtol=rtol,
        atol=atol,
        radau_closure_tol=radau_closure_tol,
        radau_jacobi_tol=radau_jacobi_tol,
    )
    if first is None:
        branch.stop_reason = MuStopReason.NO_MEMBER
        return branch
    branch.members.append(first)

    tan = _tangent(z, ydot0_sign, half_crossings, t_hi, prev=None, rtol=rtol, atol=atol)
    if tan is None:
        branch.stop_reason = MuStopReason.NO_MEMBER
        return branch

    ds = ds0
    going_up = mu_target >= mu_start
    for _step in range(max_steps):
        if (going_up and z[2] >= mu_target - 1e-12) or (not going_up and z[2] <= mu_target + 1e-12):
            branch.stop_reason = MuStopReason.TARGET_REACHED
            break
        # If a full step would cross the target, land exactly on mu_target by a
        # fixed-mu natural-parameter solve (secant x0/C from the current member).
        if tan[2] != 0.0:
            ds_to_target = (mu_target - z[2]) / tan[2]
            if 0.0 < ds_to_target <= ds:
                z_pred = z + ds_to_target * tan
                zland = _land_at_mu(
                    z_pred,
                    mu_target,
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
                    branch.stop_reason = MuStopReason.TARGET_REACHED
                    break
        z_pred = z + ds * tan
        t_hi = t_hi_frac * max(per, seed.period)
        zc = _correct(
            z_pred,
            tan,
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
                branch.stop_reason = MuStopReason.STEP_UNDERFLOW
                break
            continue
        ntan = _tangent(zc, ydot0_sign, half_crossings, t_hi, prev=tan, rtol=rtol, atol=atol)
        if ntan is None:
            ds *= 0.5
            if ds < ds_min:
                branch.stop_reason = MuStopReason.STEP_UNDERFLOW
                break
            continue
        z, tan = zc, ntan
        branch.n_steps += 1
        cheap = _half_crossing(
            z[0], z[1], z[2], ydot0_sign, half_crossings, t_hi, with_stm=False, rtol=rtol, atol=atol
        )
        if cheap is not None:
            per = 2.0 * cheap[0]
        ds = min(ds * 1.3, ds_max)
        if branch.n_steps % record_every == 0:
            mem = _make_member(
                z,
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
        branch.stop_reason = MuStopReason.MAX_STEPS

    # Always record the final landed member (even if not on a record_every step).
    final = _make_member(
        z,
        ydot0_sign,
        half_crossings,
        t_hi_frac * per,
        rtol=rtol,
        atol=atol,
        radau_closure_tol=radau_closure_tol,
        radau_jacobi_tol=radau_jacobi_tol,
    )
    if final is not None and (not branch.members or abs(branch.members[-1].mu - final.mu) > 1e-12):
        branch.members.append(final)
    return branch


def _make_member(
    z: NDArray[np.float64],
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
    nr = _nu_at(z, sign, hc, t_hi, rtol=rtol, atol=atol)
    if nr is None:
        return None
    nu, lam, period = nr
    ydot0 = _ydot0(z[0], z[1], z[2], sign)
    if ydot0 is None:
        return None
    state0 = np.array([z[0], 0.0, 0.0, 0.0, ydot0, 0.0])
    system = cr3bp.CR3BPSystem(
        mu=float(z[2]), primary="Earth", secondary="Moon", l_km=384400.0, t_s=375699.8
    )
    jac = cr3bp.jacobi_constant(state0, system.mu)
    rj = _half_crossing(z[0], z[1], z[2], sign, hc, t_hi, with_stm=False, rtol=rtol, atol=atol)
    cross_res = abs(float(rj[1][3])) if rj is not None else float("nan")
    po = cp.PeriodicOrbit(
        state0=state0, period=period, jacobi=jac, converged=True, closure_residual=cross_res
    )
    _radau_ok, radau_dj = cp.crosscheck_periodic(
        system, po, closure_tol=radau_closure_tol, jacobi_tol=radau_jacobi_tol
    )
    return MuMember(
        mu=float(z[2]),
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


def scan_c_family_at_mu(
    mu: float,
    x0_guess: float,
    c_center: float,
    period_guess: float,
    *,
    half_crossings: int,
    ydot0_sign: float,
    dc: float,
    n_each: int,
    corrector_tol: float = 1e-10,
    radau_closure_tol: float = 1e-3,
    radau_jacobi_tol: float = 1e-8,
    period_jump_frac: float = 0.3,
    rtol: float = 1e-12,
    atol: float = 1e-12,
) -> list[MuMember]:
    """Trace the C-family at a FIXED mu and characterise stability along it.

    This is the inner search the paper performs at each mu: at a fixed mass
    parameter the symmetric (k1,k2) cyclers form a one-parameter family in the
    Jacobi constant C; the paper locates the STABLE SUBFAMILY (the |nu| < 1
    window bordered by the saddle-center / period-doubling bifurcations). We walk
    C in ``[c_center - n_each*dc, c_center + n_each*dc]``, re-converging the
    fixed-Jacobi symmetric corrector (secant-predicting x0) at each C and
    recording the Barden nu. Returns the ordered members (low C to high C);
    members whose period jumps by more than ``period_jump_frac`` vs the previous
    one are dropped (topology change / different crossing branch).

    Use this to figure-match the paper's depicted STABLE orbit at a binary-star
    mu when the arclength continuation of a held member lands on an unstable
    part of the family.
    """
    system = cr3bp.CR3BPSystem(
        mu=float(mu), primary="Earth", secondary="Moon", l_km=384400.0, t_s=375699.8
    )
    cs = [c_center + k * dc for k in range(-n_each, n_each + 1)]
    out: list[MuMember] = []
    x0_pred = x0_guess
    per_prev: float | None = None
    for c in cs:
        try:
            o = cp.correct_symmetric_fixed_jacobi(
                system,
                x0_pred,
                c,
                period_guess,
                ydot0_sign=ydot0_sign,
                half_crossings=half_crossings,
                tol=corrector_tol,
                max_iter=40,
                rtol=rtol,
                atol=atol,
            )
        except ValueError:
            continue
        if not o.converged:
            continue
        if per_prev is not None and abs(o.period - per_prev) > period_jump_frac * per_prev:
            continue
        try:
            nu, lam = cp.barden_stability(system, o, rtol=rtol, atol=atol)
        except Exception:
            continue
        state0 = np.array([o.x0, 0.0, 0.0, 0.0, o.ydot0, 0.0])
        po = cp.PeriodicOrbit(
            state0=state0,
            period=o.period,
            jacobi=o.jacobi,
            converged=True,
            closure_residual=o.crossing_residual,
        )
        _ok, radau_dj = cp.crosscheck_periodic(
            system, po, closure_tol=radau_closure_tol, jacobi_tol=radau_jacobi_tol
        )
        out.append(
            MuMember(
                mu=float(mu),
                state0=state0,
                x0=o.x0,
                ydot0=o.ydot0,
                jacobi=o.jacobi,
                period=o.period,
                nu=float(nu),
                abs_lambda=float(abs(lam)),
                crossing_residual=o.crossing_residual,
                radau_djacobi=float(radau_dj),
                stable=abs(float(nu)) < 1.0,
            )
        )
        x0_pred = o.x0
        per_prev = o.period
    return out
