"""Stroboscopic-map invariant circles, whiskers and connections for the planar CCR4BP.

#882. Replaces the flow-time formulation of :mod:`ccr4bp_heteroclinic_search`
(whose "connections" matched ``(x, y, vx, vy)`` at two DIFFERENT forcing phases,
see ``docs/notes/2026-10-03-882-umbriel-torus-row-adversarial-review.md``) with a
formulation in which every object lives at ONE forcing phase.

Conventions
-----------
* Model: :mod:`cyclerfinder.core.ccr4bp`, planar. States are 4-vectors
  ``(x, y, vx, vy)`` in the base (planet-moon) synodic frame; the 6-vector
  ``(x, y, 0, vx, vy, 0)`` is the embedding used by ``ccr4bp_eom``.
* The equations are time-periodic with period ``P = system.ganymede_synodic_period``
  (the perturber's synodic period). ``F`` is the time-``P`` flow map started at
  the ABSOLUTE time ``t0``; it maps the plane ``{t = t0 mod P}`` to itself.
  ``F^n`` with ``n < 0`` integrates backwards. ``t0`` is part of every object
  below and two objects can only be joined when their ``t0`` agree mod ``P``.
* Invariant circle: ``u(theta)``, ``theta`` in ``[0, 2*pi)``, stored as an ODD
  number ``N`` of nodes ``u_j = u(2*pi*j/N)`` and evaluated by trigonometric
  interpolation; ``F(u(theta)) = u(theta + rho)``.
* Bundles: ``DF(u(theta)) v(theta) = lam * v(theta + rho)`` with a real ``lam``.
  ``v`` is normalised so that the RMS over nodes of ``|v_j|`` is 1 (a single
  global scale, NOT per node, so that ``lam`` is a constant).
* Manifold points: ``Wu(theta, s, n) = F^n(u(theta) + sign*s*eps*v_u(theta))``
  and ``Ws(theta, s, n) = F^-n(u(theta) + sign*s*eps*v_s(theta))``. The
  fundamental domain of ``s`` is ``[1, |lam_u|)`` on both branches (for the
  stable branch ``1/|lam_s| = |lam_u|`` by symplecticity). The consistency
  identity is ``Wu(theta, s, n+1) = Wu(theta + rho, lam_u*s, n) + O(eps^2)``.
* A connection for ``(n_u, n_s)`` is a zero of ``Wu(theta_u, s_u, n_u) -
  Ws(theta_s, s_s, n_s)``, four equations in four unknowns, both sides at the
  same forcing phase by construction.

Integration
-----------
Many trajectories are integrated together with a vectorised planar RHS (and
optionally the planar 4x4 STM) in ONE ``solve_ivp`` call. ``solve_ivp``'s error
norm is an RMS over all components, so a batch controls each member's error
less tightly than a single integration at the same tolerance; batched results
are used for the circle corrector (tight tolerance) and for coarse clouds. All
verification numbers come from single-trajectory integrations.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field, replace
from typing import Any, Literal

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import solve_ivp
from scipy.optimize import minimize_scalar
from scipy.spatial import cKDTree

from cyclerfinder.core.ccr4bp import CCR4BPSystem

FloatArray = NDArray[np.float64]
Branch = Literal["unstable", "stable"]

_TWO_PI = 2.0 * math.pi
_PLANAR = (0, 1, 3, 4)


# ---------------------------------------------------------------------------
# Vectorised planar equations of motion.
# ---------------------------------------------------------------------------


def forcing_period(system: CCR4BPSystem) -> float:
    """The forcing period ``P`` (perturber synodic period, TU)."""
    return system.ganymede_synodic_period


def _accel_and_hessian(
    t: float,
    x: FloatArray,
    y: FloatArray,
    system: CCR4BPSystem,
    with_hessian: bool,
) -> tuple[FloatArray, FloatArray, FloatArray | None, FloatArray | None, FloatArray | None]:
    """Gravity + centrifugal part of the planar acceleration and its Hessian.

    Returns ``(gx, gy, uxx, uxy, uyy)`` where the full acceleration is
    ``ax = gx + 2*vy``, ``ay = gy - 2*vx``. Arrays are broadcast over trajectories.
    """
    mu = system.mu
    om1 = 1.0 - mu
    dx1 = x + mu
    dx2 = x - 1.0 + mu
    r1sq = dx1 * dx1 + y * y
    r2sq = dx2 * dx2 + y * y
    r1 = np.sqrt(r1sq)
    r2 = np.sqrt(r2sq)
    inv_r13 = 1.0 / (r1sq * r1)
    inv_r23 = 1.0 / (r2sq * r2)
    gx = x - om1 * dx1 * inv_r13 - mu * dx2 * inv_r23
    gy = y - om1 * y * inv_r13 - mu * y * inv_r23
    uxx = uxy = uyy = None
    if with_hessian:
        inv_r15 = inv_r13 / r1sq
        inv_r25 = inv_r23 / r2sq
        common = 1.0 - om1 * inv_r13 - mu * inv_r23
        uxx = common + 3.0 * om1 * dx1 * dx1 * inv_r15 + 3.0 * mu * dx2 * dx2 * inv_r25
        uyy = common + 3.0 * om1 * y * y * inv_r15 + 3.0 * mu * y * y * inv_r25
        uxy = 3.0 * om1 * dx1 * y * inv_r15 + 3.0 * mu * dx2 * y * inv_r25
    if system.mu_gan != 0.0:
        th = system.theta_gan0 + system.omega_gan * t
        px = system.a_gan * math.cos(th)
        py = system.a_gan * math.sin(th)
        dgx = x - px
        dgy = y - py
        d2 = dgx * dgx + dgy * dgy
        inv_d3 = 1.0 / (d2 * np.sqrt(d2))
        a3 = system.a_gan**3
        gx = gx - system.mu_gan * dgx * inv_d3 - system.mu_gan * px / a3
        gy = gy - system.mu_gan * dgy * inv_d3 - system.mu_gan * py / a3
        if with_hessian:
            assert uxx is not None and uxy is not None and uyy is not None
            inv_d5 = inv_d3 / d2
            uxx = uxx - system.mu_gan * (inv_d3 - 3.0 * dgx * dgx * inv_d5)
            uyy = uyy - system.mu_gan * (inv_d3 - 3.0 * dgy * dgy * inv_d5)
            uxy = uxy + system.mu_gan * 3.0 * dgx * dgy * inv_d5
    return gx, gy, uxx, uxy, uyy


def planar_rhs_batch(
    t: float, yflat: FloatArray, system: CCR4BPSystem, m: int, with_stm: bool
) -> FloatArray:
    """Vectorised planar CCR4BP RHS for ``m`` trajectories.

    Layout of ``yflat``: ``(4, m)`` row-major without STM; ``(20, m)`` with STM,
    rows 4..19 holding each trajectory's 4x4 STM in row-major order. The planar
    block is exact because ``z = vz = 0`` is invariant and the 6x6 variational
    matrix decouples the ``(z, vz)`` block on that plane.
    """
    rows = 20 if with_stm else 4
    yy = yflat.reshape(rows, m)
    x, y, vx, vy = yy[0], yy[1], yy[2], yy[3]
    gx, gy, uxx, uxy, uyy = _accel_and_hessian(t, x, y, system, with_stm)
    out = np.empty_like(yy)
    out[0] = vx
    out[1] = vy
    out[2] = gx + 2.0 * vy
    out[3] = gy - 2.0 * vx
    if with_stm:
        phi = yy[4:].reshape(4, 4, m)
        dphi = out[4:].reshape(4, 4, m)
        dphi[0] = phi[2]
        dphi[1] = phi[3]
        dphi[2] = uxx * phi[0] + uxy * phi[1] + 2.0 * phi[3]
        dphi[3] = uxy * phi[0] + uyy * phi[1] - 2.0 * phi[2]
    return out.reshape(-1)


def planar_jacobian(t: float, state4: FloatArray, system: CCR4BPSystem) -> FloatArray:
    """4x4 variational matrix ``A(t, x)`` of the planar vector field."""
    x = np.asarray([state4[0]], dtype=np.float64)
    y = np.asarray([state4[1]], dtype=np.float64)
    _, _, uxx, uxy, uyy = _accel_and_hessian(t, x, y, system, True)
    assert uxx is not None and uxy is not None and uyy is not None
    a = np.zeros((4, 4))
    a[0, 2] = a[1, 3] = 1.0
    a[2, 0], a[2, 1], a[2, 3] = float(uxx[0]), float(uxy[0]), 2.0
    a[3, 0], a[3, 1], a[3, 2] = float(uxy[0]), float(uyy[0]), -2.0
    return a


def to_state6(state4: FloatArray) -> FloatArray:
    """Embed a planar 4-state into the 6-state ``(x, y, 0, vx, vy, 0)``."""
    s = np.asarray(state4, dtype=np.float64)
    return np.array([s[0], s[1], 0.0, s[2], s[3], 0.0])


def to_state4(state6: FloatArray) -> FloatArray:
    """Planar part ``(x, y, vx, vy)`` of a 6-state."""
    s = np.asarray(state6, dtype=np.float64)
    return np.array([s[0], s[1], s[3], s[4]])


@dataclass(frozen=True)
class CollisionRadii:
    """Collision radii (nondimensional) of the planet, base moon and perturber.

    A trajectory that comes closer than the radius to a body is marked as a
    collision. Zero disables a body's check.
    """

    planet: float = 0.0
    base_moon: float = 0.0
    perturber: float = 0.0

    def as_array(self) -> FloatArray:
        return np.array([self.planet, self.base_moon, self.perturber])


def body_distances(system: CCR4BPSystem, t: FloatArray, x: FloatArray, y: FloatArray) -> FloatArray:
    """Distances ``(3, ...)`` to planet, base moon and perturber at times ``t``."""
    mu = system.mu
    d_p = np.hypot(x + mu, y)
    d_m = np.hypot(x - 1.0 + mu, y)
    th = system.theta_gan0 + system.omega_gan * t
    d_g = np.hypot(x - system.a_gan * np.cos(th), y - system.a_gan * np.sin(th))
    return np.stack([d_p, d_m, d_g])


@dataclass(frozen=True)
class StrobOrbit:
    """Per-period iterates of a batch of trajectories.

    ``states[k, i]`` is trajectory ``i`` after ``k`` periods (``k = 0..n``; NaN
    once the trajectory is dropped). ``stms[k, i]`` is the cumulative 4x4 STM
    from ``k = 0`` (None unless requested). ``min_dist[i]`` is the minimum
    distance to (planet, base moon, perturber) over the accepted integrator
    steps. ``dropped_at[i]`` is the period index in which a collision was
    detected (-1 if never).
    """

    states: FloatArray
    stms: FloatArray | None
    min_dist: FloatArray
    dropped_at: NDArray[np.int64]
    t0: float
    direction: int

    @property
    def alive(self) -> NDArray[np.bool_]:
        return np.asarray(self.dropped_at < 0)


def strob_iterates(
    system: CCR4BPSystem,
    states: FloatArray,
    *,
    n: int,
    t0: float = 0.0,
    with_stm: bool = False,
    radii: CollisionRadii | None = None,
    rtol: float = 1e-13,
    atol: float = 1e-13,
    batch_size: int = 512,
) -> StrobOrbit:
    """Iterate ``F^sign(n)`` ``|n|`` times for a batch of 4-states.

    ``states`` is ``(4,)`` or ``(m, 4)``. Integration restarts exactly at every
    period boundary (no dense-output interpolation in the iterates). With
    ``radii``, a trajectory that approaches a body inside its radius during a
    period is dropped from that period onwards.
    """
    s = np.atleast_2d(np.asarray(states, dtype=np.float64))
    m = s.shape[0]
    nabs = abs(int(n))
    direction = 1 if n >= 0 else -1
    period = forcing_period(system)
    out = np.full((nabs + 1, m, 4), np.nan)
    out[0] = s
    stms = np.full((nabs + 1, m, 4, 4), np.nan) if with_stm else None
    if stms is not None:
        stms[0] = np.eye(4)
    min_dist = np.full((m, 3), np.inf)
    dropped = np.full(m, -1, dtype=np.int64)
    rad = radii.as_array() if radii is not None else np.zeros(3)
    for start in range(0, m, batch_size):
        idx_all = np.arange(start, min(start + batch_size, m))
        cur = s[idx_all].copy()
        cur_phi = np.repeat(np.eye(4)[None], idx_all.size, axis=0) if with_stm else None
        active = np.ones(idx_all.size, dtype=bool)
        for k in range(1, nabs + 1):
            ids = np.nonzero(active)[0]
            if ids.size == 0:
                break
            mm = ids.size
            ta = t0 + direction * (k - 1) * period
            tb = t0 + direction * k * period
            if with_stm:
                assert cur_phi is not None
                y0 = np.concatenate([cur[ids].T, cur_phi[ids].transpose(1, 2, 0).reshape(16, mm)])
            else:
                y0 = cur[ids].T.copy()
            sol = solve_ivp(
                planar_rhs_batch,
                (ta, tb),
                y0.reshape(-1),
                args=(system, mm, with_stm),
                method="DOP853",
                rtol=rtol,
                atol=atol,
            )
            if not sol.success:
                raise RuntimeError(f"batch integration failed: {sol.message}")
            rows = 20 if with_stm else 4
            ys = sol.y.reshape(rows, mm, -1)
            d = body_distances(system, sol.t[None, :], ys[0], ys[1])  # (3, mm, nt)
            dmin = d.min(axis=2).T  # (mm, 3)
            gi = idx_all[ids]
            min_dist[gi] = np.minimum(min_dist[gi], dmin)
            end = ys[:, :, -1]
            cur[ids] = end[:4].T
            if with_stm:
                assert cur_phi is not None
                cur_phi[ids] = end[4:].reshape(4, 4, mm).transpose(2, 0, 1)
            hit = np.any(dmin < rad[None, :], axis=1)
            out[k, gi] = cur[ids]
            if stms is not None and cur_phi is not None:
                stms[k, gi] = cur_phi[ids]
            if np.any(hit):
                dropped[gi[hit]] = k
                out[k, gi[hit]] = np.nan
                if stms is not None:
                    stms[k, gi[hit]] = np.nan
                active[ids[hit]] = False
    return StrobOrbit(
        states=out, stms=stms, min_dist=min_dist, dropped_at=dropped, t0=t0, direction=direction
    )


def strob_map(
    system: CCR4BPSystem,
    state4: FloatArray,
    *,
    n: int = 1,
    t0: float = 0.0,
    with_stm: bool = False,
    rtol: float = 1e-13,
    atol: float = 1e-13,
) -> tuple[FloatArray, FloatArray | None]:
    """``F^n(state4)`` started at absolute time ``t0`` (single trajectory).

    Returns ``(state4_out, stm_4x4 or None)``. ``n < 0`` maps backwards.
    """
    orb = strob_iterates(system, state4, n=n, t0=t0, with_stm=with_stm, rtol=rtol, atol=atol)
    st = orb.states[-1, 0].copy()
    phi = orb.stms[-1, 0].copy() if orb.stms is not None else None
    return st, phi


# ---------------------------------------------------------------------------
# Trigonometric interpolation on an odd node set.
# ---------------------------------------------------------------------------


def _wavenumbers(n_nodes: int) -> NDArray[np.int64]:
    if n_nodes % 2 != 1:
        raise ValueError(f"node count must be odd, got {n_nodes}")
    k: NDArray[np.int64] = np.rint(np.fft.fftfreq(n_nodes, 1.0 / n_nodes)).astype(np.int64)
    return k


def node_angles(n_nodes: int) -> FloatArray:
    """``theta_j = 2*pi*j/N``."""
    return _TWO_PI * np.arange(n_nodes) / n_nodes


def fourier_eval(nodes: FloatArray, theta: float | FloatArray, deriv: int = 0) -> FloatArray:
    """Trigonometric interpolant of ``nodes`` (``(N, d)``) or its derivative.

    Scalar ``theta`` -> ``(d,)``; array -> ``(len, d)``.
    """
    n_nodes = nodes.shape[0]
    k = _wavenumbers(n_nodes)
    c = np.fft.fft(nodes, axis=0) / n_nodes
    th = np.atleast_1d(np.asarray(theta, dtype=np.float64))
    e = np.exp(1j * np.outer(th, k)) * (1j * k[None, :]) ** deriv
    val = np.real(e @ c)
    if np.ndim(theta) == 0:
        return np.asarray(val[0], dtype=np.float64)
    return np.asarray(val, dtype=np.float64)


def shift_matrix(n_nodes: int, delta: float) -> FloatArray:
    """Real ``N x N`` matrix ``S`` with ``(S f)_j = f_interp(theta_j + delta)``."""
    k = _wavenumbers(n_nodes)
    th = node_angles(n_nodes)
    e = np.exp(1j * np.outer(th + delta, k))
    f = np.exp(-1j * np.outer(k, th)) / n_nodes
    return np.asarray(np.real(e @ f), dtype=np.float64)


def fourier_amplitudes(nodes: FloatArray) -> FloatArray:
    """Per-harmonic amplitude ``a_k = max_c (|c_k| + |c_-k|)`` for ``k = 0..K``."""
    n_nodes = nodes.shape[0]
    kk = (n_nodes - 1) // 2
    c = np.fft.fft(nodes, axis=0) / n_nodes
    amp = np.abs(c)
    a = np.empty(kk + 1)
    a[0] = amp[0].max()
    for k in range(1, kk + 1):
        a[k] = (amp[k] + amp[n_nodes - k]).max()
    return a


def fourier_tail(nodes: FloatArray) -> float:
    """Ratio of the largest amplitude in the top quarter of harmonics to the
    largest non-constant amplitude (a resolution measure: small = resolved)."""
    a = fourier_amplitudes(nodes)
    kk = a.size - 1
    k_tail = max(1, math.ceil(0.75 * kk))
    return float(a[k_tail:].max() / max(a[1:].max(), 1e-300))


# ---------------------------------------------------------------------------
# Invariant circles of the stroboscopic map.
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class InvariantCircle:
    """Invariant circle ``F(u(theta)) = u(theta + rho)`` of the time-``P`` map at ``t0``.

    ``residual`` is the max-abs node residual ``|R(-rho) F(u_j) - u_j|`` at the
    returned nodes (NaN if never evaluated); ``converged`` means ``residual <=
    tol`` of the corrector that produced it.
    """

    system: CCR4BPSystem
    t0: float
    rho: float
    nodes: FloatArray
    residual: float = float("nan")
    converged: bool = False
    n_iter: int = 0
    residual_history: tuple[float, ...] = ()
    notes: str = ""

    @property
    def n_nodes(self) -> int:
        return int(self.nodes.shape[0])

    @property
    def thetas(self) -> FloatArray:
        return node_angles(self.n_nodes)

    def state(self, theta: float | FloatArray) -> FloatArray:
        return fourier_eval(self.nodes, theta)

    def dstate(self, theta: float | FloatArray) -> FloatArray:
        return fourier_eval(self.nodes, theta, deriv=1)

    def fourier_tail(self) -> float:
        return fourier_tail(self.nodes)

    def fourier_amplitudes(self) -> FloatArray:
        return fourier_amplitudes(self.nodes)


def circle_residual(
    circle: InvariantCircle, *, rtol: float = 1e-13, atol: float = 1e-13
) -> tuple[float, FloatArray]:
    """Max-abs invariance residual and the ``(N, 4)`` residual array."""
    orb = strob_iterates(circle.system, circle.nodes, n=1, t0=circle.t0, rtol=rtol, atol=atol)
    fu = orb.states[1]
    g = shift_matrix(circle.n_nodes, -circle.rho) @ fu - circle.nodes
    return float(np.max(np.abs(g))), g


def correct_invariant_circle(
    system: CCR4BPSystem,
    seed_nodes: FloatArray,
    rho: float,
    *,
    t0: float = 0.0,
    fix_rho: bool = True,
    tol: float = 1e-10,
    max_iter: int = 15,
    rtol: float = 1e-13,
    atol: float = 1e-13,
    phase_ref: FloatArray | None = None,
    verbose: bool = False,
) -> InvariantCircle:
    """GMOS-type Gauss-Newton corrector for an invariant circle of the time-``P`` map.

    Unknowns: the ``4N`` node values (and ``rho`` when ``fix_rho=False``).
    Equations: ``R(-rho) F(u_j) - u_j = 0`` (``4N``), the phase condition
    ``sum_j <u_j - ref_j, ref'(theta_j)> = 0`` (the forced map has no time-shift
    symmetry, only the shift in ``theta``), and, when ``rho`` is free, the
    family-selection condition ``mean_j x_j = mean_j x_j(seed)``. Each step is
    a least-squares (minimum-norm) solve, so rank deficiencies are tolerated.
    A step that raises the residual by more than a factor of 2 is halved (at
    most 6 times).
    """
    nodes = np.array(seed_nodes, dtype=np.float64)
    n_nodes = nodes.shape[0]
    ref = nodes.copy() if phase_ref is None else np.asarray(phase_ref, dtype=np.float64)
    dref = fourier_eval(ref, node_angles(n_nodes), deriv=1)
    prow = dref.reshape(-1) / np.linalg.norm(dref)
    xmean0 = float(nodes[:, 0].mean())
    rho_c = float(rho)
    history: list[float] = []

    def evaluate(u: FloatArray, r: float) -> tuple[FloatArray, FloatArray, FloatArray]:
        orb = strob_iterates(system, u, n=1, t0=t0, with_stm=True, rtol=rtol, atol=atol)
        assert orb.stms is not None
        fu = orb.states[1]
        dfs = orb.stms[1]
        g = shift_matrix(n_nodes, -r) @ fu - u
        return g, fu, dfs

    g, fu, dfs = evaluate(nodes, rho_c)
    res = float(np.max(np.abs(g)))
    history.append(res)
    it = 0
    while res > tol and it < max_iter:
        it += 1
        smat = shift_matrix(n_nodes, -rho_c)
        jac = np.einsum("jl,lab->jalb", smat, dfs).reshape(4 * n_nodes, 4 * n_nodes)
        jac -= np.eye(4 * n_nodes)
        phase_val = float(np.dot(prow, (nodes - ref).reshape(-1)))
        if fix_rho:
            a = np.vstack([jac, prow[None, :]])
            b = -np.concatenate([g.reshape(-1), [phase_val]])
        else:
            # d/drho of R(-rho) F(u) evaluated at the nodes = -(d/dtheta of it).
            drho = -fourier_eval(smat @ fu, node_angles(n_nodes), deriv=1).reshape(-1)
            mrow = np.zeros(4 * n_nodes + 1)
            mrow[0 : 4 * n_nodes : 4] = 1.0 / n_nodes
            a = np.zeros((4 * n_nodes + 2, 4 * n_nodes + 1))
            a[: 4 * n_nodes, : 4 * n_nodes] = jac
            a[: 4 * n_nodes, -1] = drho
            a[4 * n_nodes, : 4 * n_nodes] = prow
            a[4 * n_nodes + 1] = mrow
            b = -np.concatenate([g.reshape(-1), [phase_val, float(nodes[:, 0].mean()) - xmean0]])
        step = np.linalg.lstsq(a, b, rcond=None)[0]
        alpha = 1.0
        for _ in range(7):
            trial = nodes + alpha * step[: 4 * n_nodes].reshape(n_nodes, 4)
            trial_rho = rho_c if fix_rho else rho_c + alpha * float(step[-1])
            try:
                g_t, fu_t, dfs_t = evaluate(trial, trial_rho)
                res_t = float(np.max(np.abs(g_t)))
            except RuntimeError:
                res_t = float("inf")
            if res_t <= 2.0 * res or alpha < 0.02:
                break
            alpha *= 0.5
        if not np.isfinite(res_t):
            break
        nodes, rho_c, g, fu, dfs, res = trial, trial_rho, g_t, fu_t, dfs_t, res_t
        history.append(res)
        if verbose:
            print(f"  circle GN iter {it}: residual {res:.3e} (alpha {alpha:g})", flush=True)
        if len(history) >= 4 and res > 0.5 * history[-4]:
            break  # stagnation
    return InvariantCircle(
        system=system,
        t0=t0,
        rho=rho_c,
        nodes=nodes,
        residual=res,
        converged=res <= tol,
        n_iter=it,
        residual_history=tuple(history),
    )


def symmetric_periodic_orbit(
    mu: float,
    x0: float,
    vy0: float,
    half_period: float,
    *,
    tol: float = 1e-12,
    max_iter: int = 50,
    max_step: float = 0.05,
) -> tuple[FloatArray, float, float]:
    """x-axis-symmetric planar CR3BP periodic orbit at FIXED crossing point ``x0``.

    Unknowns ``(vy0, tau)``; conditions ``y(tau) = vx(tau) = 0`` (perpendicular
    re-crossing). Returns ``(state4, period = 2*tau, residual)``. Starting from
    a Keplerian guess on either apse selects the stable or the unstable member
    of a resonant pair.
    """
    base = CCR4BPSystem(mu=mu, mu_gan=0.0, a_gan=2.0, omega_gan=-0.5)
    vy = float(vy0)
    tau = float(half_period)
    res = float("inf")
    for _ in range(max_iter):
        y0 = np.concatenate([[x0, 0.0, 0.0, vy], np.eye(4).reshape(-1)])
        sol = solve_ivp(
            planar_rhs_batch,
            (0.0, tau),
            y0,
            args=(base, 1, True),
            method="DOP853",
            rtol=1e-13,
            atol=1e-13,
        )
        f = sol.y[:, -1]
        phi = f[4:].reshape(4, 4)
        d = planar_rhs_batch(tau, f[:4].copy(), base, 1, False)
        g = np.array([f[1], f[2]])
        res = float(np.linalg.norm(g))
        if res < tol:
            break
        jac = np.array([[phi[1, 3], d[1]], [phi[2, 3], d[2]]])
        dz = np.linalg.solve(jac, -g)
        nrm = float(np.linalg.norm(dz))
        if nrm > max_step:
            dz *= max_step / nrm
        vy += float(dz[0])
        tau += float(dz[1])
    return np.array([x0, 0.0, 0.0, vy]), 2.0 * tau, res


def seed_circle_from_periodic_orbit(
    system: CCR4BPSystem,
    orbit_state: FloatArray,
    orbit_period: float,
    n_nodes: int,
    *,
    t0: float = 0.0,
    rtol: float = 1e-13,
    atol: float = 1e-13,
) -> InvariantCircle:
    """Invariant circle of the UNPERTURBED (``mu_gan = 0``) map from a periodic orbit.

    ``u(theta) = orbit(theta * T / (2*pi))`` and ``rho = 2*pi*P/T mod 2*pi``.
    Exact (up to interpolation) when ``system.mu_gan == 0``; a seed otherwise.
    ``orbit_state`` may be a 4- or 6-state.
    """
    s = np.asarray(orbit_state, dtype=np.float64)
    s4 = to_state4(s) if s.size == 6 else s
    base = CCR4BPSystem(mu=system.mu, mu_gan=0.0, a_gan=system.a_gan, omega_gan=system.omega_gan)
    times = node_angles(n_nodes) * orbit_period / _TWO_PI
    nodes = np.empty((n_nodes, 4))
    nodes[0] = s4
    cur = s4.copy()
    for j in range(1, n_nodes):
        seg = solve_ivp(
            planar_rhs_batch,
            (float(times[j - 1]), float(times[j])),
            cur,
            args=(base, 1, False),
            method="DOP853",
            rtol=rtol,
            atol=atol,
        )
        cur = seg.y[:, -1].copy()
        nodes[j] = cur
    rho = (_TWO_PI * forcing_period(system) / orbit_period) % _TWO_PI
    return InvariantCircle(system=system, t0=t0, rho=rho, nodes=nodes, notes="periodic_orbit_seed")


def continue_circle_in_mass(
    circle: InvariantCircle,
    target_mu_gan: float,
    fractions: tuple[float, ...] = (0.03, 0.1, 0.25, 0.5, 0.75, 1.0),
    *,
    tol: float = 1e-10,
    max_iter: int = 12,
    verbose: bool = False,
) -> list[InvariantCircle]:
    """Natural-parameter continuation of a circle in the perturber mass at fixed ``rho``.

    Steps ``mu_gan`` through ``fractions * target_mu_gan``, warm-starting each
    correction from the previous circle and anchoring its phase condition to
    it. Stops at the first non-converged step (which is returned as the last
    element, so the caller sees the failure).
    """
    out: list[InvariantCircle] = []
    cur = circle
    for frac in fractions:
        sys_k = replace(circle.system, mu_gan=float(frac * target_mu_gan))
        nxt = correct_invariant_circle(
            sys_k,
            cur.nodes,
            cur.rho,
            t0=cur.t0,
            tol=tol,
            max_iter=max_iter,
            phase_ref=cur.nodes,
            verbose=verbose,
        )
        out.append(nxt)
        if not nxt.converged:
            break
        cur = nxt
    return out


def seed_circle_from_pseudospectral(
    torus: object, n_nodes: int, theta1: float = 0.0
) -> InvariantCircle:
    """Seed circle from a :class:`CCR4BPTorusVariationalResult`.

    Nodes ``u_j = evaluate_torus_state(torus, theta1, theta_j)``, ``rho =
    torus.rho_strob`` and ``t0 = theta1 / torus.omega1`` (the module ties
    ``theta1 = omega1 * t``). Requires ``torus.period_multiple == 1`` so that the
    torus period equals the forcing period.
    """
    from cyclerfinder.search.variational_ccr4bp_torus import (
        CCR4BPTorusVariationalResult,
        evaluate_torus_state,
    )

    if not isinstance(torus, CCR4BPTorusVariationalResult):
        raise TypeError("torus must be a CCR4BPTorusVariationalResult")
    if torus.period_multiple != 1:
        raise ValueError("only period_multiple == 1 tori map to the forcing-period map")
    th = node_angles(n_nodes)
    nodes = np.asarray(evaluate_torus_state(torus, np.full(n_nodes, theta1), th), dtype=np.float64)
    return InvariantCircle(
        system=torus.system,
        t0=theta1 / torus.omega1,
        rho=float(torus.rho_strob) % _TWO_PI,
        nodes=nodes,
        notes="pseudospectral_seed",
    )


def shift_circle_phase(circle: InvariantCircle, dt: float) -> InvariantCircle:
    """The same invariant torus seen at forcing phase ``t0 + dt``: nodes flowed by ``dt``.

    Exact (no correction needed): if ``u`` is invariant for the map at ``t0``,
    ``phi_{t0 -> t0+dt}(u)`` is invariant for the map at ``t0 + dt`` with the
    same ``rho``.
    """
    sol = solve_ivp(
        planar_rhs_batch,
        (circle.t0, circle.t0 + dt),
        circle.nodes.T.reshape(-1).copy(),
        args=(circle.system, circle.n_nodes, False),
        method="DOP853",
        rtol=1e-13,
        atol=1e-13,
    )
    nodes = sol.y[:, -1].reshape(4, circle.n_nodes).T.copy()
    return InvariantCircle(
        system=circle.system,
        t0=circle.t0 + dt,
        rho=circle.rho,
        nodes=nodes,
        residual=circle.residual,
        converged=circle.converged,
        notes=f"phase_shifted_by_{dt:.6g}",
    )


def distance_to_circle(circle: InvariantCircle, state4: FloatArray) -> tuple[float, float]:
    """``min_theta |state4 - u(theta)|`` and the minimising ``theta``."""
    s = np.asarray(state4, dtype=np.float64)
    m = 8 * circle.n_nodes
    th = _TWO_PI * np.arange(m) / m
    pts = circle.state(th)
    d = np.linalg.norm(pts - s[None, :], axis=1)
    i = int(np.argmin(d))
    h = _TWO_PI / m

    def f(t: float) -> float:
        return float(np.linalg.norm(circle.state(t) - s))

    r = minimize_scalar(
        f, bounds=(th[i] - h, th[i] + h), method="bounded", options={"xatol": 1e-12}
    )
    best = min(float(r.fun), float(d[i]))
    return best, float(r.x) % _TWO_PI


# ---------------------------------------------------------------------------
# Hyperbolic bundles.
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class HyperbolicBundles:
    """Stable / unstable bundles of an invariant circle.

    ``v_u``, ``v_s`` are ``(N, 4)`` node values. ``residual_u/s`` are the
    OFF-GRID invariance residuals ``max |DF(u(theta)) v(theta) - lam v(theta+rho)|``
    at the midpoints ``theta_j + pi/N`` (fresh STM integrations, not the
    eigenproblem's own residual). ``tail_u/s`` are the eigenfunctions' Fourier
    tails. ``spectrum`` holds the full eigenvalue list of the circle operator.
    """

    lam_u: float
    lam_s: float
    v_u: FloatArray
    v_s: FloatArray
    residual_u: float
    residual_s: float
    tail_u: float
    tail_s: float
    spectrum: NDArray[np.complex128] = field(
        repr=False, default_factory=lambda: np.zeros(0, dtype=np.complex128)
    )

    def vector(self, branch: Branch, theta: float | FloatArray, deriv: int = 0) -> FloatArray:
        v = self.v_u if branch == "unstable" else self.v_s
        return fourier_eval(v, theta, deriv=deriv)

    def lam(self, branch: Branch) -> float:
        return self.lam_u if branch == "unstable" else self.lam_s


def _pick_real_eigen(
    vals: NDArray[np.complex128],
    vecs: NDArray[np.complex128],
    n_nodes: int,
    want_unstable: bool,
    imag_tol: float,
) -> tuple[float, FloatArray, float]:
    mod = np.abs(vals)
    real_mask = np.abs(vals.imag) <= imag_tol * np.maximum(mod, 1e-300)
    side = mod > 1.0 + 1e-6 if want_unstable else mod < 1.0 - 1e-6
    cand = np.nonzero(real_mask & side)[0]
    if cand.size == 0:
        raise ValueError("no real eigenvalue off the unit circle: the circle is not hyperbolic")
    best: tuple[float, float, int] | None = None
    for ic in cand:
        v = vecs[:, ic]
        # rotate the global complex phase out, take the real eigenfunction
        j = int(np.argmax(np.abs(v)))
        vr = np.real(v * np.exp(-1j * np.angle(v[j]))).reshape(n_nodes, 4)
        tail = fourier_tail(vr)
        key = (tail, float(-mod[ic] if want_unstable else mod[ic]), int(ic))
        if best is None or key[:2] < best[:2]:
            best = key
    assert best is not None
    i = best[2]
    v = vecs[:, i]
    j = int(np.argmax(np.abs(v)))
    vr = np.real(v * np.exp(-1j * np.angle(v[j]))).reshape(n_nodes, 4)
    vr = vr / math.sqrt(float(np.mean(np.sum(vr * vr, axis=1))))
    return float(np.real(vals[i])), vr, best[0]


def hyperbolic_bundles(
    circle: InvariantCircle,
    *,
    imag_tol: float = 1e-9,
    ref_u: FloatArray | None = None,
    ref_s: FloatArray | None = None,
    check_offgrid: bool = True,
    rtol: float = 1e-13,
    atol: float = 1e-13,
) -> HyperbolicBundles:
    """Stable / unstable bundles from the eigenproblem of ``R(-rho) blockdiag(DF(u_j))``.

    The spectrum of that ``4N x 4N`` operator comes in families ``lam *
    exp(i*k*rho)`` (eigenfunctions ``v(theta) exp(i*k*theta)``). The selected
    ``lam_u`` is the real eigenvalue with ``|lam| > 1`` whose eigenfunction has
    the smallest Fourier tail (tie: largest modulus); ``lam_s`` likewise with
    ``|lam| < 1``. Signs: aligned with ``ref_u``/``ref_s`` by node-wise dot
    product when given (for continuation), else the largest-magnitude component
    of node 0 is made positive.
    """
    n_nodes = circle.n_nodes
    orb = strob_iterates(
        circle.system, circle.nodes, n=1, t0=circle.t0, with_stm=True, rtol=rtol, atol=atol
    )
    assert orb.stms is not None
    dfs = orb.stms[1]
    smat = shift_matrix(n_nodes, -circle.rho)
    op = np.einsum("jl,lab->jalb", smat, dfs).reshape(4 * n_nodes, 4 * n_nodes)
    vals, vecs = np.linalg.eig(op)
    lam_u, v_u, tail_u = _pick_real_eigen(vals, vecs, n_nodes, True, imag_tol)
    lam_s, v_s, tail_s = _pick_real_eigen(vals, vecs, n_nodes, False, imag_tol)

    def orient(v: FloatArray, ref: FloatArray | None) -> FloatArray:
        if ref is not None:
            return v if float(np.sum(v * ref)) >= 0.0 else -v
        k = int(np.argmax(np.abs(v[0])))
        return v if v[0, k] >= 0.0 else -v

    v_u = orient(v_u, ref_u)
    v_s = orient(v_s, ref_s)
    res_u = res_s = float("nan")
    if check_offgrid:
        th = node_angles(n_nodes) + math.pi / n_nodes
        pts = circle.state(th)
        orb2 = strob_iterates(
            circle.system, pts, n=1, t0=circle.t0, with_stm=True, rtol=rtol, atol=atol
        )
        assert orb2.stms is not None
        df_mid = orb2.stms[1]
        for name, lam, v in (("u", lam_u, v_u), ("s", lam_s, v_s)):
            lhs = np.einsum("jab,jb->ja", df_mid, fourier_eval(v, th))
            rhs = lam * fourier_eval(v, th + circle.rho)
            r = float(np.max(np.abs(lhs - rhs)))
            if name == "u":
                res_u = r
            else:
                res_s = r
    return HyperbolicBundles(
        lam_u=lam_u,
        lam_s=lam_s,
        v_u=v_u,
        v_s=v_s,
        residual_u=res_u,
        residual_s=res_s,
        tail_u=tail_u,
        tail_s=tail_s,
        spectrum=vals,
    )


# ---------------------------------------------------------------------------
# Manifold points and clouds.
# ---------------------------------------------------------------------------


def manifold_departure(
    circle: InvariantCircle,
    bundles: HyperbolicBundles,
    branch: Branch,
    theta: float,
    s: float,
    *,
    eps: float,
    sign: float = 1.0,
) -> tuple[FloatArray, FloatArray]:
    """Departure state ``u(theta) + sign*s*eps*v(theta)`` and its ``(4, 2)``
    derivative with respect to ``(theta, s)``."""
    v = bundles.vector(branch, theta)
    dv = bundles.vector(branch, theta, deriv=1)
    x0 = circle.state(theta) + sign * s * eps * v
    d = np.empty((4, 2))
    d[:, 0] = circle.dstate(theta) + sign * s * eps * dv
    d[:, 1] = sign * eps * v
    return x0, d


def manifold_point(
    circle: InvariantCircle,
    bundles: HyperbolicBundles,
    branch: Branch,
    theta: float,
    s: float,
    n: int,
    *,
    eps: float,
    sign: float = 1.0,
    with_jac: bool = False,
    rtol: float = 1e-13,
    atol: float = 1e-13,
) -> tuple[FloatArray, FloatArray | None]:
    """``Wu(theta, s, n)`` (forward ``n`` periods) or ``Ws(theta, s, n)``
    (backward ``n`` periods), with the optional ``(4, 2)`` Jacobian in
    ``(theta, s)``."""
    x0, d0 = manifold_departure(circle, bundles, branch, theta, s, eps=eps, sign=sign)
    if n == 0:
        return x0, (d0 if with_jac else None)
    nn = n if branch == "unstable" else -n
    xf, phi = strob_map(
        circle.system, x0, n=nn, t0=circle.t0, with_stm=with_jac, rtol=rtol, atol=atol
    )
    if with_jac:
        assert phi is not None
        return xf, phi @ d0
    return xf, None


@dataclass(frozen=True)
class ManifoldCloud:
    """Sampled manifold: ``points[n, i]`` = ``W(theta_i, s_i, n)``, NaN if dropped."""

    branch: Branch
    sign: float
    eps: float
    t0: float
    thetas: FloatArray
    s_values: FloatArray
    points: FloatArray
    dropped_at: NDArray[np.int64]
    min_dist: FloatArray

    @property
    def n_dropped(self) -> int:
        return int(np.sum(self.dropped_at >= 0))


def manifold_cloud(
    circle: InvariantCircle,
    bundles: HyperbolicBundles,
    branch: Branch,
    *,
    n_theta: int,
    n_s: int,
    n_max: int,
    eps: float,
    sign: float = 1.0,
    radii: CollisionRadii | None = None,
    rtol: float = 1e-11,
    atol: float = 1e-11,
    batch_size: int = 512,
) -> ManifoldCloud:
    """Grid ``theta`` uniformly and ``s`` geometrically across ``[1, |lam_u|)``,
    and iterate every departure ``n_max`` periods in one batched integration
    (forward for ``unstable``, backward for ``stable``)."""
    lam = abs(bundles.lam_u) if branch == "unstable" else 1.0 / abs(bundles.lam_s)
    th = _TWO_PI * np.arange(n_theta) / n_theta
    sv = lam ** (np.arange(n_s) / n_s)
    tgrid, sgrid = np.meshgrid(th, sv, indexing="ij")
    tt: FloatArray = tgrid.reshape(-1)
    ss: FloatArray = sgrid.reshape(-1)
    v = bundles.vector(branch, tt)
    x0 = circle.state(tt) + sign * eps * ss[:, None] * v
    nn = n_max if branch == "unstable" else -n_max
    orb = strob_iterates(
        circle.system,
        x0,
        n=nn,
        t0=circle.t0,
        radii=radii,
        rtol=rtol,
        atol=atol,
        batch_size=batch_size,
    )
    return ManifoldCloud(
        branch=branch,
        sign=sign,
        eps=eps,
        t0=circle.t0,
        thetas=tt,
        s_values=ss,
        points=orb.states,
        dropped_at=orb.dropped_at,
        min_dist=orb.min_dist,
    )


# ---------------------------------------------------------------------------
# Coarse intersections and refinement.
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class ConnectionCandidate:
    theta_u: float
    s_u: float
    theta_s: float
    s_s: float
    n_u: int
    n_s: int
    sign_u: float
    sign_s: float
    distance: float


def _same_phase(system: CCR4BPSystem, t_a: float, t_b: float, tol: float = 1e-9) -> bool:
    if system.mu_gan == 0.0:
        return True  # autonomous: every phase is equivalent
    p = forcing_period(system)
    d = (t_a - t_b) / p
    return abs(d - round(d)) * p <= tol


def coarse_intersections(
    cloud_u: ManifoldCloud,
    cloud_s: ManifoldCloud,
    *,
    pairs: list[tuple[int, int]] | None = None,
    n_best: int = 10,
    vel_weight: float = 1.0,
    min_excursion: float = 0.0,
    circle: InvariantCircle | None = None,
) -> list[ConnectionCandidate]:
    """Nearest unstable/stable cloud points for each ``(n_u, n_s)`` pair.

    Distance metric: ``sqrt(dx^2 + dy^2 + w^2 (dvx^2 + dvy^2))`` with ``w =
    vel_weight``. With ``circle`` and ``min_excursion > 0`` points closer than
    ``min_excursion`` to the circle's node cloud are excluded (removes the
    trivial match near the torus). Returns the ``n_best`` closest overall.
    """
    w = np.array([1.0, 1.0, vel_weight, vel_weight])
    n_max_u = cloud_u.points.shape[0] - 1
    n_max_s = cloud_s.points.shape[0] - 1
    if pairs is None:
        pairs = [(a, b) for a in range(1, n_max_u + 1) for b in range(1, n_max_s + 1)]
    near_tree = None
    if circle is not None and min_excursion > 0.0:
        dense = circle.state(_TWO_PI * np.arange(16 * circle.n_nodes) / (16 * circle.n_nodes))
        near_tree = cKDTree(dense * w)
    out: list[ConnectionCandidate] = []
    for n_u, n_s in pairs:
        pu = cloud_u.points[n_u] * w
        ps = cloud_s.points[n_s] * w
        iu = np.nonzero(np.all(np.isfinite(pu), axis=1))[0]
        is_ = np.nonzero(np.all(np.isfinite(ps), axis=1))[0]
        if near_tree is not None:
            if iu.size:
                iu = iu[near_tree.query(pu[iu])[0] >= min_excursion]
            if is_.size:
                is_ = is_[near_tree.query(ps[is_])[0] >= min_excursion]
        if iu.size == 0 or is_.size == 0:
            continue
        tree = cKDTree(ps[is_])
        dq, jq = tree.query(pu[iu])
        d = np.asarray(dq, dtype=np.float64)
        j = np.asarray(jq, dtype=np.int64)
        order = np.argsort(d)[:n_best]
        for o in order:
            a = int(iu[o])
            b = int(is_[j[o]])
            out.append(
                ConnectionCandidate(
                    theta_u=float(cloud_u.thetas[a]),
                    s_u=float(cloud_u.s_values[a]),
                    theta_s=float(cloud_s.thetas[b]),
                    s_s=float(cloud_s.s_values[b]),
                    n_u=int(n_u),
                    n_s=int(n_s),
                    sign_u=cloud_u.sign,
                    sign_s=cloud_s.sign,
                    distance=float(d[o]),
                )
            )
    out.sort(key=lambda c: c.distance)
    return out[:n_best]


@dataclass(frozen=True)
class Connection:
    """A refined zero of ``Wu(theta_u, s_u, n_u) - Ws(theta_s, s_s, n_s)``."""

    theta_u: float
    s_u: float
    theta_s: float
    s_s: float
    n_u: int
    n_s: int
    sign_u: float
    sign_s: float
    eps: float
    residual: float
    singular_values: FloatArray
    junction_state: FloatArray
    converged: bool
    n_iter: int
    residual_history: tuple[float, ...] = ()


def refine_connection(
    circle_u: InvariantCircle,
    bundles_u: HyperbolicBundles,
    circle_s: InvariantCircle,
    bundles_s: HyperbolicBundles,
    candidate: ConnectionCandidate,
    *,
    eps: float,
    tol: float = 1e-11,
    max_iter: int = 40,
    max_dtheta: float = 0.05,
    s_window: tuple[float, float] = (-1.0, 2.0),
    rtol: float = 1e-13,
    atol: float = 1e-13,
) -> Connection:
    """Gauss-Newton (least-squares, minimum-norm steps) on ``(theta_u, ln s_u,
    theta_s, ln s_s)`` with the STM Jacobian.

    ``s`` is solved for in log form and each ``ln s`` is clipped to
    ``[s_window[0], s_window[1]] * ln|lam_u|`` (by default one fundamental
    domain below to two above ``[1, |lam_u|)``). Without this the iteration
    can slide to the TRIVIAL zero ``s_u = s_s = 0`` (both points on the circle
    itself, ``u(theta_u + n_u*rho) = u(theta_s - n_s*rho)``), which it did in
    the first R1 run. Both circles must be at the same forcing phase (raises
    ``ValueError`` otherwise: the defect this module exists to prevent).
    ``singular_values`` are those of the column-scaled Jacobian in ``(theta,
    ln s)`` at the final iterate (columns scaled to unit norm; a near-zero
    value means a non-isolated or tangential intersection)."""
    if not _same_phase(circle_u.system, circle_u.t0, circle_s.t0):
        raise ValueError("unstable and stable circles are at different forcing phases")
    if candidate.s_u <= 0.0 or candidate.s_s <= 0.0:
        raise ValueError("s must be positive (the offset sign is carried by sign_u/sign_s)")
    z = np.array(
        [candidate.theta_u, math.log(candidate.s_u), candidate.theta_s, math.log(candidate.s_s)]
    )
    llu = math.log(abs(bundles_u.lam_u))
    lls = math.log(1.0 / abs(bundles_s.lam_s))
    lo = np.array([-np.inf, s_window[0] * llu, -np.inf, s_window[0] * lls])
    hi = np.array([np.inf, s_window[1] * llu, np.inf, s_window[1] * lls])
    hist: list[float] = []
    sv = np.full(4, np.nan)
    g = np.full(4, np.nan)
    xu = np.full(4, np.nan)
    n_done = 0
    for it in range(1, max_iter + 1):
        n_done = it
        xu, ju = manifold_point(
            circle_u,
            bundles_u,
            "unstable",
            z[0],
            math.exp(z[1]),
            candidate.n_u,
            eps=eps,
            sign=candidate.sign_u,
            with_jac=True,
            rtol=rtol,
            atol=atol,
        )
        xs, js = manifold_point(
            circle_s,
            bundles_s,
            "stable",
            z[2],
            math.exp(z[3]),
            candidate.n_s,
            eps=eps,
            sign=candidate.sign_s,
            with_jac=True,
            rtol=rtol,
            atol=atol,
        )
        assert ju is not None and js is not None
        g = xu - xs
        res = float(np.linalg.norm(g))
        hist.append(res)
        jac = np.hstack([ju, -js])
        jac[:, 1] *= math.exp(z[1])
        jac[:, 3] *= math.exp(z[3])
        cn = np.linalg.norm(jac, axis=0)
        cn[cn == 0.0] = 1.0
        sv = np.linalg.svd(jac / cn, compute_uv=False)
        if res <= tol:
            break
        dz = np.linalg.lstsq(jac / cn, -g, rcond=1e-12)[0] / cn
        scale = max(abs(dz[0]) / max_dtheta, abs(dz[2]) / max_dtheta, 1.0)
        z = np.clip(z + dz / scale, lo, hi)
    return Connection(
        theta_u=float(z[0]) % _TWO_PI,
        s_u=math.exp(float(z[1])),
        theta_s=float(z[2]) % _TWO_PI,
        s_s=math.exp(float(z[3])),
        n_u=candidate.n_u,
        n_s=candidate.n_s,
        sign_u=candidate.sign_u,
        sign_s=candidate.sign_s,
        eps=eps,
        residual=hist[-1] if hist else float("nan"),
        singular_values=sv,
        junction_state=xu,
        converged=bool(hist and hist[-1] <= tol),
        n_iter=n_done,
        residual_history=tuple(hist),
    )


# ---------------------------------------------------------------------------
# Verification: one continuous trajectory.
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class VerificationThresholds:
    junction_residual: float = 1e-9
    final_over_excursion: float = 1e-3
    excursion_over_offset: float = 100.0
    integrator_agreement: float = 1e-6


@dataclass(frozen=True)
class ConnectionVerification:
    """Outcome of :func:`verify_trajectory` / :func:`verify_connection`.

    ``dist_to_circle[k]`` is the distance of the DOP853 state after ``k``
    periods to the stable-side circle (minimised over ``theta``).
    ``integrator_disagreement`` is ``|x_DOP853 - x_Radau|`` at the end.
    ``min_dist`` = (planet, base moon, perturber) over the whole trajectory.
    ``jacobi_drift`` = max |C - C(start)| of the base CR3BP Jacobi constant
    along the trajectory (a conserved quantity only when ``mu_gan == 0``).
    """

    genuine: bool
    reasons: tuple[str, ...]
    phase_consistent: bool
    junction_residual: float
    dist_to_circle: FloatArray
    max_excursion: float
    final_distance: float
    offset_size: float
    integrator_disagreement: float
    min_dist: FloatArray
    jacobi_drift: float
    jacobi_start: float
    states: FloatArray


def _jacobi4(states: FloatArray, mu: float) -> FloatArray:
    x, y, vx, vy = states[..., 0], states[..., 1], states[..., 2], states[..., 3]
    r1 = np.hypot(x + mu, y)
    r2 = np.hypot(x - 1.0 + mu, y)
    return np.asarray(x * x + y * y + 2.0 * (1.0 - mu) / r1 + 2.0 * mu / r2 - vx * vx - vy * vy)


def _integrate_single(
    system: CCR4BPSystem,
    x0: FloatArray,
    t_start: float,
    n_periods: int,
    method: Literal["DOP853", "Radau"],
    rtol: float,
    atol: float,
) -> tuple[FloatArray, FloatArray, float, FloatArray]:
    """One continuous trajectory over ``n_periods``.

    The integrator is restarted at every period boundary from the exact end
    state of the previous period (no state is patched; only the step-size
    controller restarts), so the per-period states are integrator step
    endpoints rather than dense-output interpolants. Returns per-period states
    ``(n+1, 4)``, body min distances ``(3,)`` (accepted steps plus 3 dense
    samples per step), the base Jacobi-constant drift, and the dense states.
    """
    period = forcing_period(system)

    def rhs(t: float, y: FloatArray) -> FloatArray:
        return planar_rhs_batch(t, y, system, 1, False)

    def jac(t: float, y: FloatArray) -> FloatArray:
        return planar_jacobian(t, y, system)

    solver: Any = solve_ivp
    per = np.empty((n_periods + 1, 4))
    per[0] = np.asarray(x0, dtype=np.float64)
    dense: list[FloatArray] = []
    dense_times: list[FloatArray] = []
    sub = np.linspace(0.0, 1.0, 4, endpoint=False)
    for k in range(n_periods):
        span = (t_start + k * period, t_start + (k + 1) * period)
        y0 = per[k].copy()
        if method == "Radau":
            sol = solver(
                rhs, span, y0, method="Radau", rtol=rtol, atol=atol, dense_output=True, jac=jac
            )
        else:
            sol = solver(rhs, span, y0, method=method, rtol=rtol, atol=atol, dense_output=True)
        if not sol.success:
            raise RuntimeError(f"{method} integration failed: {sol.message}")
        per[k + 1] = sol.y[:, -1]
        tt = sol.t
        dt_k = (tt[:-1, None] + np.diff(tt)[:, None] * sub[None, :]).reshape(-1)
        dense_times.append(dt_k)
        dense.append(np.asarray(sol.sol(dt_k)))
    if dense:
        dense_t = np.concatenate([*dense_times, np.array([t_start + n_periods * period])])
        ys = np.concatenate([*dense, per[-1][:, None]], axis=1)
    else:
        dense_t = np.array([t_start])
        ys = per[0][:, None].copy()
    d = body_distances(system, dense_t, ys[0], ys[1]).min(axis=1)
    cj = _jacobi4(ys.T, system.mu)
    return per, d, float(np.max(np.abs(cj - cj[0]))), ys.T


def verify_trajectory(
    system: CCR4BPSystem,
    x0: FloatArray,
    t_start: float,
    n_u: int,
    n_s: int,
    circle_s: InvariantCircle,
    *,
    junction_target: FloatArray,
    offset_size: float,
    thresholds: VerificationThresholds | None = None,
    rtol: float = 1e-13,
    atol: float = 1e-13,
    radau_rtol: float = 1e-12,
    radau_atol: float = 1e-13,
) -> ConnectionVerification:
    """Integrate ``x0`` (at absolute time ``t_start``) through ``n_u + n_s``
    periods in ONE integration and test that it is a single connecting trajectory.

    ``junction_target`` is the independently computed stable-side point
    ``Ws(theta_s, s_s, n_s)``; the junction residual is its distance to the
    continuous trajectory's state after ``n_u`` periods. ``genuine`` requires:
    phase consistency of ``t_start`` and ``circle_s.t0``; junction residual,
    final-distance/excursion, excursion/offset and DOP853-vs-Radau agreement
    within ``thresholds``.
    """
    th = thresholds or VerificationThresholds()
    n_tot = n_u + n_s
    per, dmin, jdrift, _ = _integrate_single(system, x0, t_start, n_tot, "DOP853", rtol, atol)
    per_r, dmin_r, _, _ = _integrate_single(
        system, x0, t_start, n_tot, "Radau", radau_rtol, radau_atol
    )
    phase_ok = _same_phase(system, t_start, circle_s.t0)
    dists = np.array([distance_to_circle(circle_s, p)[0] for p in per])
    jres = float(np.linalg.norm(per[n_u] - np.asarray(junction_target)))
    max_exc = float(dists.max())
    final = float(dists[-1])
    disagree = float(np.linalg.norm(per[-1] - per_r[-1]))
    reasons: list[str] = []
    if not phase_ok:
        reasons.append("forcing phase of departure differs from the stable circle's phase")
    if not jres <= th.junction_residual:
        reasons.append(f"junction residual {jres:.3e} > {th.junction_residual:.1e}")
    if not final <= th.final_over_excursion * max_exc:
        reasons.append(
            f"final distance {final:.3e} > {th.final_over_excursion:g} x excursion {max_exc:.3e}"
        )
    if not max_exc >= th.excursion_over_offset * offset_size:
        reasons.append(
            f"excursion {max_exc:.3e} < {th.excursion_over_offset:g} x offset {offset_size:.3e}"
        )
    if not disagree <= th.integrator_agreement:
        reasons.append(f"DOP853/Radau disagreement {disagree:.3e} > {th.integrator_agreement:.1e}")
    c0 = float(_jacobi4(np.asarray(x0), system.mu))
    return ConnectionVerification(
        genuine=not reasons,
        reasons=tuple(reasons),
        phase_consistent=phase_ok,
        junction_residual=jres,
        dist_to_circle=dists,
        max_excursion=max_exc,
        final_distance=final,
        offset_size=offset_size,
        integrator_disagreement=disagree,
        min_dist=np.minimum(dmin, dmin_r),
        jacobi_drift=jdrift,
        jacobi_start=c0,
        states=per,
    )


def verify_connection(
    circle_u: InvariantCircle,
    bundles_u: HyperbolicBundles,
    circle_s: InvariantCircle,
    bundles_s: HyperbolicBundles,
    connection: Connection,
    *,
    thresholds: VerificationThresholds | None = None,
    rtol: float = 1e-13,
    atol: float = 1e-13,
) -> ConnectionVerification:
    """Verify a refined connection as ONE trajectory (see :func:`verify_trajectory`).

    The departure is ``u_u(theta_u) + sign_u*s_u*eps*v_u(theta_u)`` at
    ``circle_u.t0``; the junction target is an independent evaluation of
    ``Ws(theta_s, s_s, n_s)``. ``offset_size = eps * s_u``.
    """
    c = connection
    x0, _ = manifold_departure(
        circle_u, bundles_u, "unstable", c.theta_u, c.s_u, eps=c.eps, sign=c.sign_u
    )
    target, _ = manifold_point(
        circle_s,
        bundles_s,
        "stable",
        c.theta_s,
        c.s_s,
        c.n_s,
        eps=c.eps,
        sign=c.sign_s,
        rtol=rtol,
        atol=atol,
    )
    return verify_trajectory(
        circle_u.system,
        x0,
        circle_u.t0,
        c.n_u,
        c.n_s,
        circle_s,
        junction_target=target,
        offset_size=c.eps * abs(c.s_u),
        thresholds=thresholds,
        rtol=rtol,
        atol=atol,
    )
