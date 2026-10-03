"""Sun-forced continuation of Earth-Moon CR3BP periodic orbits into the BCR4BP (#884).

Question
--------
Which catalogued Earth-Moon cycler-type periodic orbits survive when the Sun's
periodic forcing is added? Brown, Peterson, Henry & Scheeres (SIAM J. Appl. Dyn.
Syst. 24(1):346-375, 2025; filed in the private paper corpus as
``brown-peterson-henry-scheeres-2025-periodic-orbit-families-hill-restricted-4-body-problem-siads-24-1-346-doi-10.1137-24M1637301-published.pdf``)
state that a CR3BP periodic orbit of period ``T*`` continues into a periodically
forced model only if ``a T* = b Tg`` (``Tg`` the forcing period) and only at zeros
of a Melnikov-type function of the orbit phase. They test libration-point families
in the Hill restricted four-body problem (HR4BP). This module applies the same
programme to the project's bicircular model (:mod:`cyclerfinder.core.bcr4bp`).

Model and homotopy
------------------
The forced model is the standard incoherent BCR4BP with the Andreu / Rosales-Jorba
constants. A homotopy parameter ``eps`` multiplies the WHOLE solar term (direct and
indirect), i.e. ``mu_sun = eps * mu_sun_phys``; ``eps = 0`` is the CR3BP and
``eps = 1`` the physical BCR4BP. The Sun's synodic rate ``omega_S`` is held at its
physical value throughout, so the forcing period is ``Tg = 2 pi / omega_S`` (the
synodic month, 6.79119 TU in this model) at every ``eps``. (Brown et al. instead
hold the forcing period fixed in a time unit that itself depends on their parameter
``m``; their starting CR3BP members are therefore not ours. See the #884 note.)

Stroboscopic map
----------------
For a forced period ``P = n Tg`` a periodic orbit is a fixed point of the map
``X -> phi_eps(X; t0 = 0 -> P, theta_S(0) = theta0)``. No phase condition is
needed (the system is non-autonomous); ``theta0`` is a parameter. The map is
evaluated by multiple shooting (``N`` equal segments) for the long unstable
cyclers; the fixed-point equations are ``phi(X_k; t_k -> t_{k+1}) - X_{k+1} = 0``
cyclically.

Phase selection (Melnikov)
--------------------------
At ``eps = 0`` every point of a commensurate CR3BP orbit is a fixed point. The
Lyapunov-Schmidt solvability condition on the cokernel of ``M^a - I`` (spanned by
the Jacobi-constant gradient, since ``grad C^T M = grad C^T`` on the orbit) gives
the bifurcation function

    Mel(theta0) = d/d eps [ C(phi_eps^P(x(0))) - C(x(0)) ] at eps = 0
                = -2 * integral_0^P v(t) . a_sun(x(t), theta0 + omega_S t) dt

(per unit ``eps``; ``a_sun`` at ``mu_sun = mu_sun_phys``), the work of the solar
force along the unperturbed orbit, i.e. Brown et al.'s Eq. 2.7 up to the factor
``-2``. A shift of the start point along the orbit is equivalent to a shift of
``theta0`` (their Proposition 1), so it suffices to scan ``theta0`` at a fixed
start point. With ``T* = (n/a) Tg`` the function is ``2 pi / a``-periodic in
``theta0``; each zero in ``[0, 2 pi / a)`` is a candidate dynamical equivalent.

Discipline
----------
Pure numerics; no catalogue writeback; reuses only the published BCR4BP constants
from :mod:`cyclerfinder.core.bcr4bp` (read-only).
"""

from __future__ import annotations

import math
import time
from collections.abc import Callable
from dataclasses import dataclass, field
from typing import Any

import numpy as np
from numba import njit
from numpy.typing import NDArray
from scipy.integrate import solve_ivp

import cyclerfinder.core.bcr4bp as bcr4bp

FloatArr = NDArray[np.float64]

#: Earth-Moon distance used to convert lengths to km (the catalogue's EM unit).
EM_LENGTH_KM: float = 384_400.0
#: Lunar Laplace sphere-of-influence radius used by ``data/validate.py`` to tell a
#: cycler (periselene inside) from a resonant periodic orbit (outside), km.
LUNAR_SOI_KM: float = 66_182.9


@dataclass(frozen=True)
class ForcedModel:
    """BCR4BP with a homotopy parameter on the solar term."""

    mu: float
    mu_sun_phys: float
    a_sun: float
    omega_sun: float

    @property
    def tg(self) -> float:
        """Forcing period: the Sun's synodic period ``2 pi / omega_S`` (TU)."""
        return 2.0 * math.pi / self.omega_sun


def default_model() -> ForcedModel:
    """The Andreu / Rosales-Jorba BCR4BP constants of :func:`bcr4bp.andreu_default`."""
    s = bcr4bp.andreu_default()
    return ForcedModel(
        mu=s.mu, mu_sun_phys=s.mu_sun, a_sun=s.a_sun_nondim, omega_sun=s.omega_sun_nondim
    )


# ---------------------------------------------------------------------------
# Equations of motion (numba).
# ---------------------------------------------------------------------------


@njit(cache=True)  # type: ignore[untyped-decorator]
def _rhs(
    t: float, y: FloatArr, mu: float, mus: float, a_s: float, w_s: float, th0: float
) -> FloatArr:
    x, yy, z, vx, vy, vz = y[0], y[1], y[2], y[3], y[4], y[5]
    r1 = math.sqrt((x + mu) ** 2 + yy * yy + z * z)
    r2 = math.sqrt((x - 1.0 + mu) ** 2 + yy * yy + z * z)
    r13 = r1**3
    r23 = r2**3
    ax = 2.0 * vy + x - (1.0 - mu) * (x + mu) / r13 - mu * (x - 1.0 + mu) / r23
    ay = -2.0 * vx + yy - (1.0 - mu) * yy / r13 - mu * yy / r23
    az = -(1.0 - mu) * z / r13 - mu * z / r23
    if mus != 0.0:
        th = th0 + w_s * t
        sx = a_s * math.cos(th)
        sy = a_s * math.sin(th)
        dx = x - sx
        dy = yy - sy
        d2 = dx * dx + dy * dy + z * z
        d3 = d2 * math.sqrt(d2)
        a3 = a_s**3
        ax += -mus * dx / d3 - mus * sx / a3
        ay += -mus * dy / d3 - mus * sy / a3
        az += -mus * z / d3
    out = np.empty(6)
    out[0] = vx
    out[1] = vy
    out[2] = vz
    out[3] = ax
    out[4] = ay
    out[5] = az
    return out


@njit(cache=True)  # type: ignore[untyped-decorator]
def _rhs_var(
    t: float,
    y: FloatArr,
    mu: float,
    mus_phys: float,
    eps: float,
    a_s: float,
    w_s: float,
    th0: float,
) -> FloatArr:
    """State (6) + STM (36, row-major) + d state / d eps (6)."""
    mus = eps * mus_phys
    f = _rhs(t, y[:6], mu, mus, a_s, w_s, th0)
    x, yy, z = y[0], y[1], y[2]
    r1 = math.sqrt((x + mu) ** 2 + yy * yy + z * z)
    r2 = math.sqrt((x - 1.0 + mu) ** 2 + yy * yy + z * z)
    r13, r23 = r1**3, r2**3
    r15, r25 = r1**5, r2**5
    om1 = 1.0 - mu
    uxx = (
        1 - om1 / r13 - mu / r23 + 3 * om1 * (x + mu) ** 2 / r15 + 3 * mu * (x - 1 + mu) ** 2 / r25
    )
    uyy = 1 - om1 / r13 - mu / r23 + 3 * om1 * yy * yy / r15 + 3 * mu * yy * yy / r25
    uzz = -om1 / r13 - mu / r23 + 3 * om1 * z * z / r15 + 3 * mu * z * z / r25
    uxy = 3 * om1 * (x + mu) * yy / r15 + 3 * mu * (x - 1 + mu) * yy / r25
    uxz = 3 * om1 * (x + mu) * z / r15 + 3 * mu * (x - 1 + mu) * z / r25
    uyz = 3 * om1 * yy * z / r15 + 3 * mu * yy * z / r25
    th = th0 + w_s * t
    sx = a_s * math.cos(th)
    sy = a_s * math.sin(th)
    dx = x - sx
    dy = yy - sy
    d2 = dx * dx + dy * dy + z * z
    d3 = d2 * math.sqrt(d2)
    d5 = d3 * d2
    a3 = a_s**3
    # Unit-eps solar acceleration (mu_sun = mus_phys) for the parameter sensitivity.
    gx = -mus_phys * dx / d3 - mus_phys * sx / a3
    gy = -mus_phys * dy / d3 - mus_phys * sy / a3
    gz = -mus_phys * z / d3
    if mus != 0.0:
        uxx += -mus * (1.0 / d3 - 3.0 * dx * dx / d5)
        uyy += -mus * (1.0 / d3 - 3.0 * dy * dy / d5)
        uzz += -mus * (1.0 / d3 - 3.0 * z * z / d5)
        uxy += 3.0 * mus * dx * dy / d5
        uxz += 3.0 * mus * dx * z / d5
        uyz += 3.0 * mus * dy * z / d5
    a = np.zeros((6, 6))
    a[0, 3] = 1.0
    a[1, 4] = 1.0
    a[2, 5] = 1.0
    a[3, 0] = uxx
    a[3, 1] = uxy
    a[3, 2] = uxz
    a[4, 0] = uxy
    a[4, 1] = uyy
    a[4, 2] = uyz
    a[5, 0] = uxz
    a[5, 1] = uyz
    a[5, 2] = uzz
    a[3, 4] = 2.0
    a[4, 3] = -2.0
    phi = y[6:42].reshape((6, 6))
    dphi = a @ phi
    s = y[42:48]
    ds = a @ s
    ds[3] += gx
    ds[4] += gy
    ds[5] += gz
    out = np.empty(48)
    out[:6] = f
    out[6:42] = dphi.reshape(36)
    out[42:48] = ds
    return out


def sun_acc_unit(model: ForcedModel, r: FloatArr, theta: FloatArr) -> FloatArr:
    """Vectorised solar acceleration at ``mu_sun = mu_sun_phys`` (unit ``eps``).

    ``r`` has shape (..., 3) and ``theta`` broadcasts against ``r[..., 0]``.
    """
    sx = model.a_sun * np.cos(theta)
    sy = model.a_sun * np.sin(theta)
    dx = r[..., 0] - sx
    dy = r[..., 1] - sy
    dz = r[..., 2]
    d3 = (dx * dx + dy * dy + dz * dz) ** 1.5
    a3 = model.a_sun**3
    m = model.mu_sun_phys
    return np.stack([-m * dx / d3 - m * sx / a3, -m * dy / d3 - m * sy / a3, -m * dz / d3], axis=-1)


# ---------------------------------------------------------------------------
# Propagation.
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Arc:
    """End state of a propagation; STM and eps-sensitivity when requested."""

    state_f: FloatArr
    stm: FloatArr | None
    dxdeps: FloatArr | None


def propagate(
    model: ForcedModel,
    eps: float,
    theta0: float,
    x0: FloatArr,
    t0: float,
    t1: float,
    *,
    variational: bool = False,
    method: str = "DOP853",
    rtol: float = 1e-12,
    atol: float = 1e-13,
    dense: bool = False,
) -> Any:
    """Propagate the eps-scaled BCR4BP from ``t0`` to ``t1`` (Sun angle ``theta0 + w t``).

    Returns an :class:`Arc`, or the raw ``solve_ivp`` solution when ``dense``.
    ``method`` may be any ``solve_ivp`` method; ``Radau`` / ``LSODA`` provide the
    independent-integrator cross-check.
    """
    mu, ms, a_s, w_s = model.mu, model.mu_sun_phys, model.a_sun, model.omega_sun
    if variational:
        y0 = np.concatenate([np.asarray(x0, dtype=np.float64), np.eye(6).reshape(36), np.zeros(6)])

        def fun(t: float, y: FloatArr) -> FloatArr:
            return _rhs_var(t, y, mu, ms, eps, a_s, w_s, theta0)  # type: ignore[no-any-return]

    else:
        y0 = np.asarray(x0, dtype=np.float64).copy()
        mus = eps * ms

        def fun(t: float, y: FloatArr) -> FloatArr:
            return _rhs(t, y, mu, mus, a_s, w_s, theta0)  # type: ignore[no-any-return]

    sol = solve_ivp(  # type: ignore[call-overload]
        fun, (t0, t1), y0, method=method, rtol=rtol, atol=atol, dense_output=dense
    )
    if not sol.success:
        raise RuntimeError(f"propagation failed at t={sol.t[-1]:.6f}: {sol.message}")
    if dense:
        return sol
    yf = sol.y[:, -1]
    if variational:
        return Arc(
            state_f=yf[:6].copy(), stm=yf[6:42].reshape(6, 6).copy(), dxdeps=yf[42:48].copy()
        )
    return Arc(state_f=yf[:6].copy(), stm=None, dxdeps=None)


def jacobi(state: FloatArr, mu: float) -> float:
    """Jacobi constant ``C = x^2 + y^2 + 2(1-mu)/r1 + 2 mu/r2 - v^2`` (repo convention)."""
    x, y, z, vx, vy, vz = (float(v) for v in state)
    r1 = math.sqrt((x + mu) ** 2 + y * y + z * z)
    r2 = math.sqrt((x - 1.0 + mu) ** 2 + y * y + z * z)
    return x * x + y * y + 2 * (1 - mu) / r1 + 2 * mu / r2 - (vx * vx + vy * vy + vz * vz)


def grad_jacobi(state: FloatArr, mu: float) -> FloatArr:
    """Gradient of :func:`jacobi` with respect to the 6-state."""
    x, y, z, vx, vy, vz = (float(v) for v in state)
    r1 = math.sqrt((x + mu) ** 2 + y * y + z * z)
    r2 = math.sqrt((x - 1.0 + mu) ** 2 + y * y + z * z)
    c1 = 2 * (1 - mu) / r1**3
    c2 = 2 * mu / r2**3
    return np.array(
        [
            2 * x - c1 * (x + mu) - c2 * (x - 1 + mu),
            2 * y - c1 * y - c2 * y,
            -c1 * z - c2 * z,
            -2 * vx,
            -2 * vy,
            -2 * vz,
        ]
    )


def to_canonical_matrix() -> FloatArr:
    """Linear map (x, v) -> (q, p) with ``p = v + omega x r`` (``px = vx - y``, ``py = vy + x``)."""
    m = np.eye(6)
    m[3, 1] = -1.0
    m[4, 0] = 1.0
    return m


def symplectic_defect(stm: FloatArr) -> float:
    """``max |Mc^T J Mc - J|`` with ``Mc`` the STM in canonical coordinates."""
    c = to_canonical_matrix()
    mc = c @ stm @ np.linalg.inv(c)
    j = np.zeros((6, 6))
    j[:3, 3:] = np.eye(3)
    j[3:, :3] = -np.eye(3)
    return float(np.max(np.abs(mc.T @ j @ mc - j)))


# ---------------------------------------------------------------------------
# CR3BP family member at a prescribed period (symmetric half-period shooting).
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class SymmetricOrbit:
    """x-axis / xz-plane symmetric CR3BP periodic orbit ``(x0, 0, z0, 0, vy0, 0)``."""

    x0: float
    z0: float
    vy0: float
    period: float
    residual: float

    @property
    def state(self) -> FloatArr:
        return np.array([self.x0, 0.0, self.z0, 0.0, self.vy0, 0.0])


def _half_residual(
    model: ForcedModel, u: FloatArr, period: float, spatial: bool
) -> tuple[FloatArr, FloatArr, FloatArr]:
    """Residual at T/2 and its Jacobians w.r.t. the unknowns and the period."""
    x0 = np.array([u[0], 0.0, u[1] if spatial else 0.0, 0.0, u[-1], 0.0])
    arc = propagate(model, 0.0, 0.0, x0, 0.0, 0.5 * period, variational=True)
    assert arc.stm is not None
    xf = arc.state_f
    f = _rhs(0.5 * period, xf, model.mu, 0.0, model.a_sun, model.omega_sun, 0.0)
    rows = [1, 3, 5] if spatial else [1, 3]
    cols = [0, 2, 4] if spatial else [0, 4]
    res = xf[rows]
    jac = arc.stm[np.ix_(rows, cols)]
    dt = 0.5 * f[rows]
    return res, jac, dt


def correct_symmetric_fixed_period(
    model: ForcedModel,
    guess: FloatArr,
    period: float,
    *,
    tol: float = 1e-12,
    max_iter: int = 30,
) -> SymmetricOrbit:
    """Newton on ``(x0, [z0], vy0)`` at fixed period. ``guess`` = (x0, z0, vy0)."""
    spatial = abs(float(guess[1])) > 0.0
    u = np.array([guess[0], guess[1], guess[2]]) if spatial else np.array([guess[0], guess[2]])
    res = np.array([np.inf])
    for _ in range(max_iter):
        res, jac, _ = _half_residual(model, u, period, spatial)
        if float(np.max(np.abs(res))) < tol:
            break
        u = u - np.linalg.solve(jac, res)
    nres = float(np.max(np.abs(res)))
    return SymmetricOrbit(
        x0=float(u[0]),
        z0=float(u[1]) if spatial else 0.0,
        vy0=float(u[-1]),
        period=float(period),
        residual=nres,
    )


def continue_to_period(
    model: ForcedModel,
    seed: FloatArr,
    seed_period: float,
    target_period: float,
    **kw: Any,
) -> tuple[SymmetricOrbit | None, dict[str, Any]]:
    """Walk a symmetric family from ``seed`` toward ``target_period``; member there or None."""
    direction = 1.0 if target_period > seed_period else -1.0
    members, info = walk_family(model, seed, seed_period, [target_period], direction, **kw)
    return (members[0] if members else None), info


def walk_family(
    model: ForcedModel,
    seed: FloatArr,
    seed_period: float,
    targets: list[float],
    direction: float,
    *,
    ds: float = 0.005,
    ds_max: float = 0.05,
    max_steps: int = 400,
    tol: float = 1e-11,
    t_window: tuple[float, float] = (0.5, 60.0),
    wall_s: float = math.inf,
    log: Callable[[str], None] | None = None,
) -> tuple[list[SymmetricOrbit], dict[str, Any]]:
    """Pseudo-arclength along a symmetric family, collecting members at target periods.

    ``seed`` = (x0, z0, vy0) at period ``seed_period``, corrected first at that
    period. ``direction`` (+1/-1) sets the initial sense of the period. Every
    crossing of a period in ``targets`` is corrected at exactly that period and
    returned (a family folding in period can cross a target more than once).
    A step is accepted only if the corrector lands within half a step of the
    predictor (guards against jumping to another family). The info dict records
    the period range visited, a subsampled trace ``(T, C, x0, vy0)`` and the stop
    reason.
    """
    spatial = abs(float(seed[1])) > 0.0
    first = correct_symmetric_fixed_period(model, seed, seed_period)
    found: list[SymmetricOrbit] = []
    if first.residual > 1e-9:
        return found, {"reason": f"seed did not correct (res {first.residual:.2e})"}
    u = np.array([first.x0, first.z0, first.vy0]) if spatial else np.array([first.x0, first.vy0])
    z = np.concatenate([u, [seed_period]])
    _, jac, dt = _half_residual(model, u, seed_period, spatial)
    tan = np.linalg.svd(np.column_stack([jac, dt]))[2][-1]
    if tan[-1] * direction < 0:
        tan = -tan
    t_min = t_max = seed_period
    trace: list[list[float]] = [[seed_period, jacobi(first.state, model.mu), first.x0, first.vy0]]
    info: dict[str, Any] = {"reason": "max_steps", "seed_corrected": first.state.tolist()}
    step = 0
    t_start = time.monotonic()
    for step in range(max_steps):
        if time.monotonic() - t_start > wall_s:
            info["reason"] = f"wall-clock limit at T={z[-1]:.6f}"
            break
        pred = z + ds * tan
        zz = pred.copy()
        ok = False
        try:
            for _ in range(8):
                res, jac, dt = _half_residual(model, zz[:-1], zz[-1], spatial)
                if float(np.max(np.abs(res))) < tol:
                    ok = True
                    break
                big = np.vstack([np.column_stack([jac, dt]), tan])
                rhs = np.concatenate([res, [tan @ (zz - pred)]])
                zz = zz - np.linalg.solve(big, rhs)
        except (np.linalg.LinAlgError, RuntimeError, ValueError):
            ok = False
        if ok and float(np.max(np.abs(zz - pred))) > 0.5 * ds:
            ok = False  # corrector jumped: refuse and shorten the step
        if not ok:
            ds *= 0.5
            if ds < 1e-7:
                info["reason"] = f"step collapse at T={z[-1]:.6f}"
                break
            continue
        newtan = np.linalg.svd(np.column_stack([jac, dt]))[2][-1]
        if newtan @ tan < 0:
            newtan = -newtan
        t_old, t_new = float(z[-1]), float(zz[-1])
        z_prev = z
        z, tan = zz, newtan
        t_min, t_max = min(t_min, t_new), max(t_max, t_new)
        if step % 5 == 0:
            st = np.array([z[0], 0.0, z[1] if spatial else 0.0, 0.0, z[-2], 0.0])
            trace.append([t_new, jacobi(st, model.mu), float(z[0]), float(z[-2])])
        if log is not None and step % 50 == 0:
            log(f"  family step {step}: T={t_new:.6f} x0={z[0]:.6f}")
        for tp in targets:
            if (t_old - tp) * (t_new - tp) <= 0.0 and t_old != tp:
                w = (tp - t_old) / (t_new - t_old)
                zi = (1 - w) * z_prev + w * z  # chord interpolation
                guess = np.array([zi[0], zi[1] if spatial else 0.0, zi[-2]])
                try:
                    member = correct_symmetric_fixed_period(model, guess, tp)
                except (np.linalg.LinAlgError, RuntimeError):
                    continue
                near = abs(member.x0 - guess[0]) + abs(member.vy0 - guess[2]) < 1e-3
                if member.residual < 1e-9 and near:
                    found.append(member)
        if len(found) and len(targets) == 1:
            info["reason"] = "reached"
            break
        if abs(z[0]) > 5 or abs(z[-2]) > 10 or not (t_window[0] < t_new < t_window[1]):
            info["reason"] = f"left domain at T={t_new:.6f}"
            break
        r1 = math.hypot(z[0] + model.mu, z[1] if spatial else 0.0)
        r2 = math.hypot(z[0] - 1 + model.mu, z[1] if spatial else 0.0)
        if min(r1, r2) < 0.005:
            info["reason"] = f"start point near a primary at T={t_new:.6f}"
            break
        ds = min(ds * 1.3, ds_max)
    info.update({"T_min": t_min, "T_max": t_max, "steps": step, "trace": trace})
    return found, info


# ---------------------------------------------------------------------------
# Melnikov function.
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class MelnikovSamples:
    """Unperturbed orbit sampled on a uniform grid over the forced period ``P``."""

    t: FloatArr
    r: FloatArr
    v: FloatArr
    weights: FloatArr  # composite Simpson weights including h / 3


def melnikov_samples(
    model: ForcedModel,
    x0: FloatArr,
    period_forced: float,
    *,
    period_orbit: float | None = None,
    samples_per_tu: int = 2000,
) -> MelnikovSamples:
    """Sample the CR3BP orbit over ``[0, P]``.

    With ``period_orbit = T*`` only ONE lap is integrated and the orbit is
    extended periodically (``x(t) = x(t mod T*)``). For unstable orbits with
    ``a = P / T* > 1`` this is required: integrating ``a`` laps in one shot
    departs from the periodic orbit by the growth of round-off along the
    unstable direction.
    """
    span = period_orbit if period_orbit is not None else period_forced
    sol = propagate(model, 0.0, 0.0, x0, 0.0, span, dense=True)
    n = int(samples_per_tu * period_forced)
    n = n + 1 if n % 2 == 0 else n
    t = np.linspace(0.0, period_forced, n)
    tm = np.mod(t, span) if period_orbit is not None else t
    if period_orbit is not None:
        tm[-1] = span if abs(t[-1] - round(t[-1] / span) * span) < 1e-9 else tm[-1]
    st = sol.sol(tm)
    h = t[1] - t[0]
    w = np.ones(n)
    w[1:-1:2] = 4.0
    w[2:-1:2] = 2.0
    return MelnikovSamples(t=t, r=st[:3].T.copy(), v=st[3:].T.copy(), weights=w * h / 3.0)


def melnikov_eval(model: ForcedModel, smp: MelnikovSamples, thetas: FloatArr) -> FloatArr:
    """``Mel(theta0) = -2 int_0^P v . a_sun(x(t), theta0 + w t) dt`` (Simpson)."""
    out = np.empty(len(thetas))
    for i, th in enumerate(np.atleast_1d(thetas)):
        acc = sun_acc_unit(model, smp.r, th + model.omega_sun * smp.t)
        out[i] = float(smp.weights @ (-2.0 * np.sum(smp.v * acc, axis=1)))
    return out


def melnikov_scan(
    model: ForcedModel,
    x0: FloatArr,
    period_forced: float,
    thetas: FloatArr,
    *,
    period_orbit: float | None = None,
    samples_per_tu: int = 2000,
) -> FloatArr:
    """Convenience wrapper: sample the orbit, then evaluate the Melnikov function."""
    smp = melnikov_samples(
        model, x0, period_forced, period_orbit=period_orbit, samples_per_tu=samples_per_tu
    )
    return melnikov_eval(model, smp, thetas)


def melnikov_variational(
    model: ForcedModel, x0: FloatArr, period_forced: float, theta0: float
) -> float:
    """Melnikov value from the eps-sensitivity ODE: ``grad C(x(P)) . dX/deps(P)``."""
    arc = propagate(model, 0.0, theta0, x0, 0.0, period_forced, variational=True)
    assert arc.dxdeps is not None
    return float(grad_jacobi(arc.state_f, model.mu) @ arc.dxdeps)


def find_zeros(
    model: ForcedModel,
    x0: FloatArr,
    period_forced: float,
    a: int,
    *,
    period_orbit: float | None = None,
    n_grid: int = 240,
    samples_per_tu: int = 2000,
) -> tuple[FloatArr, FloatArr, list[float]]:
    """Scan ``theta0`` over one Melnikov period (length ``2 pi / a``); refine sign changes.

    The grid is offset by half a cell so that zeros at the symmetric phases
    (multiples of ``pi / a``) fall strictly inside a cell.
    """
    span = 2.0 * math.pi / a
    smp = melnikov_samples(
        model, x0, period_forced, period_orbit=period_orbit, samples_per_tu=samples_per_tu
    )
    cell = span / n_grid
    th = np.asarray(-0.5 * cell + cell * np.arange(n_grid + 1), dtype=np.float64)
    mel = melnikov_eval(model, smp, th)
    zeros: list[float] = []
    for i in range(n_grid):
        f0, f1 = mel[i], mel[i + 1]
        if f0 * f1 < 0:
            lo, hi, flo = th[i], th[i + 1], f0
            for _ in range(45):
                mid = 0.5 * (lo + hi)
                fm = float(melnikov_eval(model, smp, np.array([mid]))[0])
                if fm * flo <= 0:
                    hi = mid
                else:
                    lo, flo = mid, fm
            zeros.append(float(0.5 * (lo + hi)) % span)
    return th, mel, sorted(zeros)


# ---------------------------------------------------------------------------
# Multiple-shooting fixed point of the stroboscopic map; continuation in eps.
# ---------------------------------------------------------------------------


@dataclass
class ShootingProblem:
    """Fixed point of the ``P``-stroboscopic map with ``N`` segments at Sun phase ``theta0``."""

    model: ForcedModel
    period: float
    theta0: float
    n_seg: int

    def times(self) -> FloatArr:
        return np.linspace(0.0, self.period, self.n_seg + 1)

    def evaluate(
        self, xs: FloatArr, eps: float
    ) -> tuple[FloatArr, FloatArr, FloatArr, list[FloatArr]]:
        """Residual (6N), Jacobian (6N x 6N), d residual / d eps (6N), segment STMs."""
        n = self.n_seg
        tt = self.times()
        res = np.empty(6 * n)
        jac = np.zeros((6 * n, 6 * n))
        jeps = np.empty(6 * n)
        stms: list[FloatArr] = []
        for k in range(n):
            arc = propagate(self.model, eps, self.theta0, xs[k], tt[k], tt[k + 1], variational=True)
            assert arc.stm is not None and arc.dxdeps is not None
            kn = (k + 1) % n
            res[6 * k : 6 * k + 6] = arc.state_f - xs[kn]
            jac[6 * k : 6 * k + 6, 6 * k : 6 * k + 6] += arc.stm
            jac[6 * k : 6 * k + 6, 6 * kn : 6 * kn + 6] -= np.eye(6)
            jeps[6 * k : 6 * k + 6] = arc.dxdeps
            stms.append(arc.stm)
        return res, jac, jeps, stms

    def monodromy(self, stms: list[FloatArr]) -> FloatArr:
        m = np.eye(6)
        for s in stms:
            m = s @ m
        return m

    def nodes_from_orbit(
        self, x0: FloatArr, eps: float = 0.0, period_orbit: float | None = None
    ) -> FloatArr:
        """Nodes at the segment times of the trajectory from ``x0``.

        With ``period_orbit = T*`` (unperturbed periodic orbit, ``eps = 0``) one lap
        is integrated and extended periodically, which keeps the nodes on the
        orbit for unstable orbits spanning several laps.
        """
        tt = self.times()[:-1]
        if period_orbit is not None:
            sol = propagate(self.model, 0.0, 0.0, x0, 0.0, period_orbit, dense=True)
            return np.asarray(sol.sol(np.mod(tt, period_orbit)).T, dtype=np.float64)
        sol = propagate(self.model, eps, self.theta0, x0, 0.0, self.period, dense=True)
        return np.asarray(sol.sol(tt).T, dtype=np.float64)


def newton_fixed_eps(
    prob: ShootingProblem, xs: FloatArr, eps: float, *, tol: float = 1e-11, max_iter: int = 15
) -> tuple[FloatArr, float, bool]:
    """Newton (least squares) on the shooting residual at fixed ``eps``."""
    nrm = math.inf
    for _ in range(max_iter):
        res, jac, _, _ = prob.evaluate(xs, eps)
        nrm = float(np.max(np.abs(res)))
        if nrm < tol:
            return xs, nrm, True
        dx = np.linalg.lstsq(jac, -res, rcond=None)[0]
        step = float(np.max(np.abs(dx)))
        if step > 0.1:
            dx *= 0.1 / step
        xs = xs + dx.reshape(xs.shape)
    res, _, _, _ = prob.evaluate(xs, eps)
    nrm = float(np.max(np.abs(res)))
    return xs, nrm, nrm < tol


@dataclass
class Branch:
    """One continuation branch in ``eps`` at fixed Sun phase."""

    theta0: float
    eps: list[float] = field(default_factory=list)
    nodes: list[FloatArr] = field(default_factory=list)
    max_abs_eig: list[float] = field(default_factory=list)
    stop_reason: str = ""
    folds: list[float] = field(default_factory=list)


def _tangent(jac: FloatArr, jeps: FloatArr, prev: FloatArr | None, sign_eps: float) -> FloatArr:
    full = np.column_stack([jac, jeps])
    t = np.linalg.svd(full)[2][-1]
    if prev is not None:
        if t @ prev < 0:
            t = -t
    elif t[-1] * sign_eps < 0:
        t = -t
    return np.asarray(t, dtype=np.float64)


def continue_in_eps(
    prob: ShootingProblem,
    xs0: FloatArr,
    *,
    eps_start: float = 1e-4,
    ds0: float = 0.02,
    ds_max: float = 0.1,
    ds_min: float = 1e-6,
    max_steps: int = 400,
    tol: float = 1e-10,
    log: Callable[[str], None] | None = None,
    eps_target: float = 1.0,
    reverse_from: tuple[FloatArr, float] | None = None,
    wall_s: float = math.inf,
) -> Branch:
    """Pseudo-arclength in ``(nodes, eps)`` from the CR3BP orbit (or a given point).

    Default start: nodes ``xs0`` on the unperturbed orbit at a Melnikov zero; a
    fixed-``eps`` Newton solve at ``eps_start`` from the first-order predictor
    ``x + eps Y`` (``Y`` the minimum-norm solution of ``J Y = -J_eps``) leaves the
    degenerate ``eps = 0`` plane before arclength stepping begins. Stops on
    reaching ``eps_target`` (corrected exactly there), on returning to ``eps <= 0``,
    or on step-size collapse.
    ``reverse_from=(nodes, eps)`` instead continues from a converged point toward
    ``eps = 0`` (reversibility gate).
    """
    n6 = 6 * prob.n_seg
    br = Branch(theta0=prob.theta0)
    if reverse_from is None:
        _, jac0, jeps0, _ = prob.evaluate(xs0, 0.0)
        y = np.linalg.lstsq(jac0, -jeps0, rcond=None)[0]
        ok = False
        nrm = math.inf
        xs = xs0
        e_try = eps_start
        for e_try in (eps_start, 0.1 * eps_start, 0.01 * eps_start):
            # The first-order predictor is valid for eps << 1 / max|lambda|; very
            # unstable orbits need a smaller first step.
            xs, nrm, ok = newton_fixed_eps(prob, xs0 + e_try * y.reshape(xs0.shape), e_try, tol=tol)
            if ok:
                break
        if not ok:
            br.stop_reason = f"start Newton failed down to eps={e_try:.1e} (res {nrm:.2e})"
            return br
        eps = e_try
        sign = 1.0
    else:
        xs, eps = reverse_from
        sign = -1.0
        eps_target = 0.0
    res, jac, jeps, stms = prob.evaluate(xs, eps)
    br.eps.append(eps)
    br.nodes.append(xs.copy())
    br.max_abs_eig.append(float(np.max(np.abs(np.linalg.eigvals(prob.monodromy(stms))))))
    z = np.concatenate([xs.reshape(n6), [eps]])
    tan = _tangent(jac, jeps, None, sign)
    ds = ds0
    t_start = time.monotonic()
    for step in range(max_steps):
        if time.monotonic() - t_start > wall_s:
            br.stop_reason = f"wall-clock limit at eps={z[-1]:.6g}"
            return br
        pred = z + ds * tan
        zz = pred.copy()
        ok = False
        nrm = math.inf
        for _ in range(8):
            res, jac, jeps, stms = prob.evaluate(zz[:n6].reshape(-1, 6), float(zz[-1]))
            nrm = float(np.max(np.abs(res)))
            arc_res = float(tan @ (zz - pred))
            if nrm < tol and abs(arc_res) < 1e-10:
                ok = True
                break
            big = np.vstack([np.column_stack([jac, jeps]), tan])
            try:
                dz = np.linalg.solve(big, -np.concatenate([res, [arc_res]]))
            except np.linalg.LinAlgError:
                break
            if not np.all(np.isfinite(dz)) or float(np.max(np.abs(dz))) > 0.5:
                break
            zz = zz + dz
        if not ok:
            ds *= 0.5
            if ds < ds_min:
                br.stop_reason = f"step collapse at eps={z[-1]:.6g}"
                return br
            continue
        newtan = _tangent(jac, jeps, tan, sign)
        if newtan[-1] * tan[-1] < 0:
            br.folds.append(float(zz[-1]))
        e_old, e_new = float(z[-1]), float(zz[-1])
        z, tan = zz, newtan
        br.eps.append(e_new)
        br.nodes.append(z[:n6].reshape(-1, 6).copy())
        br.max_abs_eig.append(float(np.max(np.abs(np.linalg.eigvals(prob.monodromy(stms))))))
        if log is not None:
            log(
                f"    step {step}: eps={e_new:.6f} ds={ds:.3g} "
                f"max|lam|={br.max_abs_eig[-1]:.3e} folds={len(br.folds)}"
            )
        hit = (e_old - eps_target) * (e_new - eps_target) <= 0.0
        if hit:
            # Correct exactly at the target from a linear interpolation.
            w = (eps_target - e_old) / (e_new - e_old) if e_new != e_old else 1.0
            guess = br.nodes[-2] + w * (br.nodes[-1] - br.nodes[-2])
            xs_t, nrm_t, ok_t = newton_fixed_eps(prob, guess, eps_target, tol=tol)
            if ok_t:
                _, _, _, stms_t = prob.evaluate(xs_t, eps_target)
                br.eps.append(eps_target)
                br.nodes.append(xs_t)
                br.max_abs_eig.append(
                    float(np.max(np.abs(np.linalg.eigvals(prob.monodromy(stms_t)))))
                )
                br.stop_reason = "reached_target"
            else:
                br.stop_reason = f"target correction failed (res {nrm_t:.2e})"
            return br
        if reverse_from is None and e_new <= 0.0:
            br.stop_reason = "returned_to_eps0"
            return br
        if float(np.max(np.abs(z[:n6]))) > 20:
            br.stop_reason = f"left domain at eps={e_new:.6g}"
            return br
        ds = min(ds * 1.5, ds_max)
    br.stop_reason = f"max_steps at eps={z[-1]:.6g}"
    return br


# ---------------------------------------------------------------------------
# Diagnostics of a converged forced orbit.
# ---------------------------------------------------------------------------


def orbit_diagnostics(
    prob: ShootingProblem, xs: FloatArr, eps: float, *, samples_per_tu: int = 400
) -> dict[str, Any]:
    """Closure (DOP853 nodes, Radau nodes), monodromy, Floquet, periselene/perigee."""
    res, _, _, stms = prob.evaluate(xs, eps)
    mono = prob.monodromy(stms)
    eig = np.linalg.eigvals(mono)
    tt = prob.times()
    radau = 0.0
    r_moon = math.inf
    r_earth = math.inf
    mu = prob.model.mu
    for k in range(prob.n_seg):
        arc = propagate(prob.model, eps, prob.theta0, xs[k], tt[k], tt[k + 1], method="Radau")
        radau = max(radau, float(np.max(np.abs(arc.state_f - xs[(k + 1) % prob.n_seg]))))
        sol = propagate(prob.model, eps, prob.theta0, xs[k], tt[k], tt[k + 1], dense=True)
        ts = np.linspace(tt[k], tt[k + 1], max(50, int(samples_per_tu * (tt[k + 1] - tt[k]))))
        st = sol.sol(ts)
        dm = np.sqrt((st[0] - 1 + mu) ** 2 + st[1] ** 2 + st[2] ** 2)
        de = np.sqrt((st[0] + mu) ** 2 + st[1] ** 2 + st[2] ** 2)
        r_moon = min(r_moon, float(dm.min()))
        r_earth = min(r_earth, float(de.min()))
    order = np.argsort(-np.abs(eig))
    eig = eig[order]
    return {
        "closure_dop853": float(np.max(np.abs(res))),
        "closure_radau": radau,
        "floquet": [[float(e.real), float(e.imag)] for e in eig],
        "max_abs_floquet": float(np.max(np.abs(eig))),
        "symplectic_defect": symplectic_defect(mono),
        "periselene_km": r_moon * EM_LENGTH_KM,
        "perigee_km": r_earth * EM_LENGTH_KM,
        "cycler_class": bool(r_moon * EM_LENGTH_KM <= LUNAR_SOI_KM),
        "stability": classify(eig),
    }


def classify(eig: NDArray[Any]) -> str:
    """Linear stability of the forced orbit from its monodromy eigenvalues."""
    mags = np.abs(eig)
    if float(np.max(mags)) < 1.0 + 1e-6:
        return "linearly_stable"
    n_unst = int(np.sum(mags > 1.0 + 1e-6))
    return f"unstable({n_unst})"
