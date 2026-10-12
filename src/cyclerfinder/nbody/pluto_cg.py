"""#1073: continuous-gravity lane for a binary with one massive secondary (Pluto-Charon).

Pre-registration: docs/notes/2026-10-12-1073-pluto-continuous-gravity-lane-prereg.md.

Model (sec. 2-3 of the note). Barycentric inertial frame, km, s, km^3/s^2. The primary and the
secondary ride circles about the barycentre at the binary period (secondary at (1 - mu_S) a_rel,
angle ``theta0 + n t``; primary opposite at mu_S a_rel); massless targets ride their own circles
(kinematic waypoints). With ``n^2 a_rel^3 = GM_total`` the rails are the exact two-body solution, so
at sigma = 1 the field is the circular restricted three-body problem.

The sigma homotopy is three point masses that co-rotate with the secondary:
``(1 - sigma) GM_total`` at the barycentre, ``sigma GM_primary`` on the primary's circle and
``sigma GM_secondary`` on the secondary's circle. At sigma = 0 the field is the monopole of the
#1024 patched conic; the field is time-independent in the rotating frame at every sigma, so a Jacobi
integral exists at every sigma (:func:`jacobi`).

Softening (a divergence guard, never a physical model): the secondary at ``sigma R_S``, the primary
at ``R_P`` (unscaled: the primary is not a flyby node), the barycentre point at 1 km.

The corrector is the #1025 closed chain generalised to three node kinds:
- ``flyby``: secondary-relative position and velocity offsets plus the epoch (7 unknowns), with the
  periapsis gauge row;
- ``pin``: the massless target's position at the node epoch (plus an optional out-of-plane miss
  ``b_z``, the note's M2 diagnostic), with the velocity relative to the target and the epoch free
  (4 or 5 unknowns);
- ``free``: a full inertial state at a fixed epoch (6 unknowns; the published-orbit control).
Forward-backward half-arcs are matched at each leg's mid-time; the last leg closes on the first node
rotated by the secondary's advance over the period (the Q wrap). Optionally the period is an
unknown and components of the first node's state are fixed by extra rows (the PC1 control).
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass, field
from typing import Any

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares, minimize_scalar

Arr = NDArray[np.float64]

W = np.array([1.0, 1.0, 1.0, 1e3, 1e3, 1e3])
W_GAUGE = 1e4
W_EXTRA = 1.0
RTOL, ATOL = 1e-13, 1e-12
BARY_CLAMP_KM = 1.0
FLOOR_DR_KM, FLOOR_DV_KMS, FLOOR_GAUGE = 1e-3, 1e-6, 1e-9
NOISE_DR_KM, NOISE_GAUGE = 1e-2, 1e-7


def _rot6(theta: float) -> Arr:
    c, s = math.cos(theta), math.sin(theta)
    q = np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])
    q6 = np.zeros((6, 6))
    q6[:3, :3] = q
    q6[3:, 3:] = q
    return q6


def _circle(radius: float, period: float, theta0: float, t: float) -> tuple[Arr, Arr, Arr]:
    n = 2.0 * math.pi / period
    th = theta0 + n * t
    c, s = math.cos(th), math.sin(th)
    r = np.array([radius * c, radius * s, 0.0])
    v = np.array([-radius * n * s, radius * n * c, 0.0])
    return r, v, -(n**2) * r


@dataclass(frozen=True)
class BinaryCG:
    """A binary (primary + one massive secondary) on circular rails, plus massless targets.

    ``targets`` maps a name to ``(radius_km, period_s, theta0_rad)`` about the barycentre.
    """

    gm_total: float
    gm_secondary: float
    a_rel_km: float
    period_s: float
    secondary_radius_km: float
    primary_radius_km: float
    targets: dict[str, tuple[float, float, float]] = field(default_factory=dict)
    theta0: float = 0.0

    @property
    def mass_fraction(self) -> float:
        return self.gm_secondary / self.gm_total

    @property
    def gm_primary(self) -> float:
        return self.gm_total - self.gm_secondary

    @property
    def n(self) -> float:
        return 2.0 * math.pi / self.period_s

    def secondary_state(self, t: float) -> tuple[Arr, Arr, Arr]:
        """Position, velocity, acceleration of the secondary at ``t``."""
        return _circle((1.0 - self.mass_fraction) * self.a_rel_km, self.period_s, self.theta0, t)

    def primary_state(self, t: float) -> tuple[Arr, Arr, Arr]:
        return _circle(self.mass_fraction * self.a_rel_km, self.period_s, self.theta0 + math.pi, t)

    def target_state(self, name: str, t: float) -> tuple[Arr, Arr, Arr]:
        a, per, th0 = self.targets[name]
        return _circle(a, per, th0, t)

    def masses(self, sigma: float, t: float) -> list[tuple[float, Arr, float]]:
        """``(GM, position, softening radius)`` of the three point masses at ``t``."""
        return [
            ((1.0 - sigma) * self.gm_total, np.zeros(3), BARY_CLAMP_KM),
            (sigma * self.gm_primary, self.primary_state(t)[0], self.primary_radius_km),
            (
                sigma * self.gm_secondary,
                self.secondary_state(t)[0],
                sigma * self.secondary_radius_km,
            ),
        ]


def accel_and_gradient(m: BinaryCG, sigma: float, r: Arr, t: float) -> tuple[Arr, Arr]:
    """Acceleration and its 3x3 gradient at ``r`` (km) and ``t`` (s)."""
    a = np.zeros(3)
    g = np.zeros((3, 3))
    eye = np.eye(3)
    for gm, p, soft in m.masses(sigma, t):
        if gm == 0.0:
            continue
        d = r - p
        dn = math.sqrt(float(d @ d))
        if dn < soft:
            a -= gm * d / soft**3
            g -= gm / soft**3 * eye
        else:
            a -= gm * d / dn**3
            g -= gm / dn**3 * (eye - 3.0 * np.outer(d, d) / dn**2)
    return a, g


def potential(m: BinaryCG, sigma: float, r: Arr, t: float) -> float:
    """U with a = grad U (unsoftened; for the Jacobi integral away from the bodies)."""
    u = 0.0
    for gm, p, _ in m.masses(sigma, t):
        if gm != 0.0:
            u += gm / float(np.linalg.norm(r - p))
    return u


def jacobi(m: BinaryCG, sigma: float, x: Arr, t: float) -> float:
    """Jacobi integral in inertial coordinates: 2U + 2n (x v_y - y v_x) - v^2 (km^2/s^2)."""
    r, v = x[:3], x[3:]
    return float(2.0 * potential(m, sigma, r, t) + 2.0 * m.n * (r[0] * v[1] - r[1] * v[0]) - v @ v)


def wrap_q6(m: BinaryCG, dt: float) -> Arr:
    """Rotation by the secondary's advance over ``dt`` (6x6, acting on a state)."""
    return _rot6(m.n * dt)


def _rhs(m: BinaryCG, sigma: float, stm: bool) -> Any:
    def f(t: float, y: Arr) -> Arr:
        a, g = accel_and_gradient(m, sigma, y[:3], t)
        if not stm:
            return np.concatenate([y[3:6], a])
        phi = y[6:].reshape(6, 6)
        dphi = np.empty((6, 6))
        dphi[:3] = phi[3:]
        dphi[3:] = g @ phi[:3]
        return np.concatenate([y[3:6], a, dphi.ravel()])

    return f


def _solve(
    m: BinaryCG,
    sigma: float,
    y0: Arr,
    t0: float,
    t1: float,
    stm: bool,
    rtol: float,
    atol: float,
    dense: bool,
) -> Any:
    sol = solve_ivp(
        _rhs(m, sigma, stm),
        (t0, t1),
        y0,
        method="DOP853",
        rtol=rtol,
        atol=atol,
        dense_output=dense,
    )
    if not sol.success:
        raise RuntimeError(f"propagation failed: {sol.message}")
    return sol


def propagate(
    m: BinaryCG,
    sigma: float,
    x0: Arr,
    t0: float,
    t1: float,
    *,
    stm: bool = True,
    rtol: float = RTOL,
    atol: float = ATOL,
) -> tuple[Arr, Arr]:
    """DOP853 state (+ STM) from ``t0`` to ``t1``. Returns ``(x1, Phi)`` (Phi = I without stm).

    Raises RuntimeError when the integrator fails.
    """
    if t1 == t0:
        return np.asarray(x0, dtype=np.float64).copy(), np.eye(6)
    y0 = np.concatenate([x0, np.eye(6).ravel()]) if stm else np.asarray(x0, dtype=np.float64)
    y1 = np.asarray(_solve(m, sigma, y0, t0, t1, stm, rtol, atol, False).y[:, -1])
    return y1[:6], (y1[6:].reshape(6, 6) if stm else np.eye(6))


def propagate_dense(
    m: BinaryCG,
    sigma: float,
    x0: Arr,
    t0: float,
    t1: float,
    *,
    rtol: float = RTOL,
    atol: float = ATOL,
) -> Any:
    """The dense-output interpolant (callable on t, returns the 6-state) of one arc."""
    return _solve(m, sigma, np.asarray(x0, dtype=np.float64), t0, t1, False, rtol, atol, True).sol


def flyby_seed_offsets(gm: float, vinf_in: Arr, vinf_out: Arr) -> tuple[Arr, Arr]:
    """Periapsis offsets (position, velocity) of the hyperbola about a body of ``gm`` that turns
    ``vinf_in`` into ``vinf_out`` (the ``jovian.periapsis_node`` construction, unclamped)."""
    vi, vo = float(np.linalg.norm(vinf_in)), float(np.linalg.norm(vinf_out))
    vmag = 0.5 * (vi + vo)
    e_in, e_out = vinf_in / vi, vinf_out / vo
    delta = math.acos(max(-1.0, min(1.0, float(e_in @ e_out))))
    sh = math.sin(0.5 * delta)
    if sh < 1e-12:
        raise ValueError("zero turn: no periapsis node")
    rp = gm * (1.0 / sh - 1.0) / vmag**2
    e_p = (e_in - e_out) / float(np.linalg.norm(e_in - e_out))
    v_dir = e_in + e_out
    v_dir = v_dir - float(v_dir @ e_p) * e_p
    e_v = v_dir / float(np.linalg.norm(v_dir))
    return rp * e_p, math.sqrt(vmag**2 + 2.0 * gm / rp) * e_v


@dataclass(frozen=True)
class NodeSpec:
    kind: str  # "flyby" | "pin" | "free"
    target: str = ""  # pin: the target name
    t_fixed: float = 0.0  # free: the node epoch


@dataclass
class Chain:
    """Closed forward-backward chain over one period (sec. 4.1 of the note)."""

    model: BinaryCG
    specs: Sequence[NodeSpec]
    period_s: float
    sigma: float
    free_bz: bool = False
    free_period: bool = False
    extra: Sequence[tuple[int, float]] = ()  # (component of node 0's state, value): fixed rows

    def __post_init__(self) -> None:
        self.offsets: list[int] = []
        k = 0
        for s in self.specs:
            self.offsets.append(k)
            k += {"flyby": 7, "pin": 5 if self.free_bz else 4, "free": 6}[s.kind]
        self._n = k + (1 if self.free_period else 0)

    @property
    def n_unknowns(self) -> int:
        return self._n

    def period(self, z: Arr) -> float:
        return float(z[-1]) if self.free_period else self.period_s

    def fd_steps(self, z: Arr) -> Arr:
        st = np.empty(self._n)
        for s, o in zip(self.specs, self.offsets, strict=True):
            if s.kind == "flyby":
                st[o : o + 3], st[o + 3 : o + 6], st[o + 6] = 1e-2, 1e-6, 10.0
            elif s.kind == "pin":
                st[o : o + 3], st[o + 3] = 1e-6, 10.0
                if self.free_bz:
                    st[o + 4] = 1e-2
            else:
                st[o : o + 3], st[o + 3 : o + 6] = 1e-2, 1e-6
        if self.free_period:
            st[-1] = 10.0
        return st

    def nodes(self, z: Arr) -> list[tuple[Arr, Arr, float, Arr]]:
        """Per node ``(x, dx/dz, t, dt/dz)``, plus the wrapped first node at the end."""
        n = self._n
        m = self.model
        out = []
        for s, o in zip(self.specs, self.offsets, strict=True):
            dx = np.zeros((6, n))
            dt = np.zeros(n)
            if s.kind == "flyby":
                t = float(z[o + 6])
                rs, vs, as_ = m.secondary_state(t)
                x = np.concatenate([rs + z[o : o + 3], vs + z[o + 3 : o + 6]])
                dx[:, o : o + 6] = np.eye(6)
                dx[:3, o + 6], dx[3:, o + 6] = vs, as_
                dt[o + 6] = 1.0
            elif s.kind == "pin":
                t = float(z[o + 3])
                rm, vm, am = m.target_state(s.target, t)
                bz = float(z[o + 4]) if self.free_bz else 0.0
                x = np.concatenate([rm + np.array([0.0, 0.0, bz]), vm + z[o : o + 3]])
                dx[3:, o : o + 3] = np.eye(3)
                dx[:3, o + 3], dx[3:, o + 3] = vm, am
                if self.free_bz:
                    dx[2, o + 4] = 1.0
                dt[o + 3] = 1.0
            else:
                t = s.t_fixed
                x = np.asarray(z[o : o + 6], dtype=np.float64).copy()
                dx[:, o : o + 6] = np.eye(6)
            out.append((x, dx, t, dt))
        per = self.period(z)
        q6 = wrap_q6(m, per)
        x0, dx0, t0, dt0 = out[0]
        dxw = q6 @ dx0
        dtw = dt0.copy()
        if self.free_period:
            dq = m.n * (wrap_q6(m, per + 1e-3 / m.n) - wrap_q6(m, per - 1e-3 / m.n)) / (2e-3)
            dxw[:, -1] += dq @ x0
            dtw[-1] += 1.0
        out.append((q6 @ x0, dxw, t0 + per, dtw))
        return out

    def f(self, x: Arr, t: float) -> Arr:
        return np.concatenate([x[3:], accel_and_gradient(self.model, self.sigma, x[:3], t)[0]])

    def _gauge(self, x: Arr, t: float) -> tuple[float, Arr, float]:
        rs, vs, as_ = self.model.secondary_state(t)
        dr, dv = x[:3] - rs, x[3:] - vs
        a, b = float(np.linalg.norm(dr)), float(np.linalg.norm(dv))
        g = float(dr @ dv) / (a * b)
        g_dr = dv / (a * b) - g * dr / a**2
        g_dv = dr / (a * b) - g * dv / b**2
        return g, np.concatenate([g_dr, g_dv]), float(-(g_dr @ vs) - (g_dv @ as_))

    def residual_and_jac(self, z: Arr, want_jac: bool = True) -> tuple[Arr, Arr, dict[str, Any]]:
        nd = self.nodes(z)
        res: list[float] = []
        rows: list[Arr] = []
        info: dict[str, Any] = {"dr": [], "dv": [], "gauge": [], "extra": []}
        n_leg = len(self.specs)
        for k in range(n_leg):
            xp, dxp, tp, dtp = nd[k]
            xn, dxn, tn, dtn = nd[k + 1]
            tm = 0.5 * (tp + tn)
            try:
                xf, pf = propagate(self.model, self.sigma, xp, tp, tm, stm=want_jac)
                xb, pb = propagate(self.model, self.sigma, xn, tn, tm, stm=want_jac)
            except RuntimeError:
                if want_jac:
                    raise
                bad = {"dr": [np.inf], "dv": [np.inf], "gauge": [np.inf], "extra": [np.inf]}
                return np.full(self.n_rows, np.inf), np.zeros((0, self._n)), bad
            d = xf - xb
            info["dr"].append(float(np.linalg.norm(d[:3])))
            info["dv"].append(float(np.linalg.norm(d[3:])))
            res.extend(W * d)
            if want_jac:
                ff, fb = self.f(xf, tm), self.f(xb, tm)
                dxf = (
                    pf @ dxp
                    + np.outer(-pf @ self.f(xp, tp) + 0.5 * ff, dtp)
                    + np.outer(0.5 * ff, dtn)
                )
                dxb = (
                    pb @ dxn
                    + np.outer(-pb @ self.f(xn, tn) + 0.5 * fb, dtn)
                    + np.outer(0.5 * fb, dtp)
                )
                rows.append(W[:, None] * (dxf - dxb))
        for k, s in enumerate(self.specs):
            if s.kind != "flyby":
                continue
            x, dx, t, dt = nd[k]
            g, gx, gt = self._gauge(x, t)
            info["gauge"].append(g)
            res.append(W_GAUGE * g)
            if want_jac:
                rows.append(W_GAUGE * (gx @ dx + gt * dt)[None, :])
        x0, dx0, _, _ = nd[0]
        for comp, val in self.extra:
            info["extra"].append(float(x0[comp] - val))
            res.append(W_EXTRA * float(x0[comp] - val) * (W[comp]))
            if want_jac:
                rows.append(W_EXTRA * W[comp] * dx0[comp][None, :])
        jac = np.vstack(rows) if want_jac else np.zeros((0, self._n))
        return np.asarray(res), jac, info

    @property
    def n_rows(self) -> int:
        return 6 * len(self.specs) + sum(s.kind == "flyby" for s in self.specs) + len(self.extra)

    def scale_offsets(self, z: Arr, f: float) -> Arr:
        """Scale every flyby node's position offset by ``f`` (the sigma scaling of #1025)."""
        z = z.copy()
        for s, o in zip(self.specs, self.offsets, strict=True):
            if s.kind == "flyby":
                z[o : o + 3] *= f
        return z

    def describe(self, z: Arr) -> list[dict[str, Any]]:
        """Per node: epoch, V_inf, periapsis, r_p / sigma, turn (flyby); speed and b_z (pin)."""
        out: list[dict[str, Any]] = []
        nd = self.nodes(z)
        m = self.model
        for k, s in enumerate(self.specs):
            x, _, t, _ = nd[k]
            rec: dict[str, Any] = {"kind": s.kind, "t_days": t / 86400.0}
            if s.kind == "flyby":
                rs, vs, _ = m.secondary_state(t)
                d = float(np.linalg.norm(x[:3] - rs))
                gm = self.sigma * m.gm_secondary
                vinf2 = float(np.linalg.norm(x[3:] - vs)) ** 2 - 2.0 * gm / d
                e = 1.0 + d * vinf2 / gm if gm > 0 else math.inf
                rec |= {
                    "vinf_kms": math.sqrt(max(vinf2, 0.0)),
                    "rp_km": d,
                    "rp_over_sigma": d / self.sigma if self.sigma > 0 else math.inf,
                    "turn_deg": math.degrees(2.0 * math.asin(1.0 / e)) if e > 1.0 else float("nan"),
                    "impact": d <= self.sigma * m.secondary_radius_km,
                    "floor_ok": d >= self.sigma * (m.secondary_radius_km + 100.0),
                }
            elif s.kind == "pin":
                rm, vm, _ = m.target_state(s.target, t)
                rec |= {
                    "target": s.target,
                    "vinf_kms": float(np.linalg.norm(x[3:] - vm)),
                    "b_z_km": float(x[2] - rm[2]),
                }
            out.append(rec)
        return out

    def arcs(self, z: Arr) -> list[tuple[Arr, float, float]]:
        """``(x_start, t_start, t_end)`` of every half-arc (forward then backward per leg)."""
        nd = self.nodes(z)
        out = []
        for k in range(len(self.specs)):
            xp, _, tp, _ = nd[k]
            xn, _, tn, _ = nd[k + 1]
            tm = 0.5 * (tp + tn)
            out += [(xp, tp, tm), (xn, tn, tm)]
        return out


def converged(info: dict[str, Any]) -> bool:
    return bool(
        max(info["dr"]) < FLOOR_DR_KM
        and max(info["dv"]) < FLOOR_DV_KMS
        and max((abs(g) for g in info["gauge"]), default=0.0) < FLOOR_GAUGE
        and max((abs(e) for e in info["extra"]), default=0.0) < FLOOR_DR_KM
    )


def noise_ok(info: dict[str, Any]) -> bool:
    """The #968 noise-floor rule for intermediate points (sigma < 1)."""
    return bool(
        max(info["dr"]) < NOISE_DR_KM
        and max(info["dv"]) < FLOOR_DV_KMS
        and max((abs(g) for g in info["gauge"]), default=0.0) < NOISE_GAUGE
        and max((abs(e) for e in info["extra"]), default=0.0) < NOISE_DR_KM
    )


class _DoneError(Exception):
    pass


BAD_INFO: dict[str, Any] = {"dr": [np.inf], "dv": [np.inf], "gauge": [np.inf], "extra": [np.inf]}


def solve(
    ch: Chain, z0: Arr, newton_iters: int = 40, trf_nfev: int = 200, accept_noise: bool = False
) -> tuple[Arr, dict[str, Any], dict[str, Any]]:
    """The #1004/#1025 solve: damped Newton (column-scaled lstsq, flyby offset step capped at half
    the offset, backtracking), then trf with the analytic Jacobian.

    Returns ``(z, info, stats)``; ``stats`` holds the iteration counts and the residual norms.
    """

    def ok(info: dict[str, Any]) -> bool:
        return converged(info) or (accept_noise and noise_ok(info))

    stats: dict[str, Any] = {"newton": 0, "trf_nfev": 0, "norms": []}
    z = np.asarray(z0, dtype=np.float64).copy()
    info = BAD_INFO
    for _ in range(newton_iters):
        try:
            r, j, info = ch.residual_and_jac(z)
        except RuntimeError:
            return z, BAD_INFO, stats
        stats["newton"] += 1
        stats["norms"].append(float(np.linalg.norm(r)))
        if ok(info):
            return z, info, stats
        d = np.linalg.norm(j, axis=0)
        d[d == 0.0] = 1.0
        step = np.linalg.lstsq(j / d, -r, rcond=None)[0] / d
        cap = 1.0
        for s, o in zip(ch.specs, ch.offsets, strict=True):
            if s.kind != "flyby":
                continue
            rr = float(np.linalg.norm(z[o : o + 3]))
            dd = float(np.linalg.norm(step[o : o + 3]))
            if dd > 0.5 * rr:
                cap = min(cap, 0.5 * rr / dd)
        step = cap * step
        f0 = float(np.linalg.norm(r))
        for alpha in [0.5**i for i in range(8)]:
            rn = ch.residual_and_jac(z + alpha * step, False)[0]
            if np.all(np.isfinite(rn)) and float(np.linalg.norm(rn)) < f0:
                z = z + alpha * step
                break
        else:
            break
    if trf_nfev > 0:
        cache: dict[str, Any] = {}

        def fun(y: Arr) -> Arr:
            try:
                r, j, inf = ch.residual_and_jac(y)
            except RuntimeError:
                return np.full(ch.n_rows, 1e12)
            stats["trf_nfev"] += 1
            cache["z"], cache["j"] = y.copy(), j
            if ok(inf):
                cache["done"] = y.copy()
                raise _DoneError
            return r

        def jac(y: Arr) -> Arr:
            if "z" in cache and np.array_equal(cache["z"], y):
                return np.asarray(cache["j"])
            return ch.residual_and_jac(y)[1]

        try:
            sol = least_squares(
                fun,
                z,
                jac=jac,
                method="trf",
                x_scale="jac",
                xtol=1e-15,
                ftol=1e-15,
                gtol=1e-15,
                max_nfev=trf_nfev,
            )
            z = np.asarray(sol.x)
        except _DoneError:
            z = np.asarray(cache["done"])
        except RuntimeError:
            return z, BAD_INFO, stats
    try:
        r, _, info = ch.residual_and_jac(z, False)
        stats["norms"].append(float(np.linalg.norm(r)))
    except RuntimeError:
        info = BAD_INFO
    return z, info, stats


def closest_approaches(ch: Chain, z: Arr, step_days: float = 0.002) -> dict[str, float]:
    """Smallest distance (km) to the primary, the barycentre and the secondary over every half-arc:
    sampled every ``step_days`` on the dense output, every interior sample minimum refined."""
    m = ch.model
    keys = ("primary", "barycentre", "secondary")
    best = dict.fromkeys(keys, math.inf)

    def dist(sol: Any, t: float, key: str) -> float:
        r = sol(t)[:3]
        if key == "primary":
            return float(np.linalg.norm(r - m.primary_state(t)[0]))
        if key == "secondary":
            return float(np.linalg.norm(r - m.secondary_state(t)[0]))
        return float(np.linalg.norm(r))

    for x0, t0, t1 in ch.arcs(z):
        if t1 == t0:
            continue
        sol = propagate_dense(m, ch.sigma, x0, t0, t1)
        ns = max(3, int(abs(t1 - t0) / (step_days * 86400.0)) + 1)
        ts = np.linspace(t0, t1, ns)
        r = np.asarray(sol(ts))[:3]
        th = m.theta0 + m.n * ts
        cs, sn = np.cos(th), np.sin(th)
        a_s = (1.0 - m.mass_fraction) * m.a_rel_km
        a_p = m.mass_fraction * m.a_rel_km
        d = {
            "secondary": np.sqrt((r[0] - a_s * cs) ** 2 + (r[1] - a_s * sn) ** 2 + r[2] ** 2),
            "primary": np.sqrt((r[0] + a_p * cs) ** 2 + (r[1] + a_p * sn) ** 2 + r[2] ** 2),
            "barycentre": np.sqrt(r[0] ** 2 + r[1] ** 2 + r[2] ** 2),
        }
        for key in keys:
            dk = d[key]
            best[key] = min(best[key], float(dk.min()))
            for i in range(1, ns - 1):
                if dk[i] <= dk[i - 1] and dk[i] <= dk[i + 1]:
                    a, b = sorted((float(ts[i - 1]), float(ts[i + 1])))
                    res = minimize_scalar(
                        lambda t, k=key, so=sol: dist(so, t, k),
                        bounds=(a, b),
                        method="bounded",
                        options={"xatol": 0.1},
                    )
                    best[key] = min(best[key], float(res.fun))
    return best


@dataclass(frozen=True)
class IdentityTolerances:
    """Sec. 6.2 of the note: amendment 10 scaled to the binary's speeds and time scales."""

    vinf_rel: float = 0.004
    rp_rel: float = 0.05
    turn_abs_deg: float = 0.2
    turn_rel: float = 0.05
    epoch_frac_of_flyby_time: float = 0.025
    max_last_sigma: float = 0.0025


def identity_check(
    ladder: Sequence[dict[str, Any]],
    ref: Sequence[dict[str, Any]],
    tol: IdentityTolerances | None = None,
) -> dict[str, Any]:
    """Member identity by the first-order (Richardson) limit of the two smallest converged sigma.

    ``ladder``: converged points ``{"sigma", "nodes"}`` (nodes as :meth:`Chain.describe`), any
    order. ``ref``: the patched-conic member per node, ``{"kind", "vinf_kms", "t_days"}`` plus
    ``"rp_km"`` (sigma = 1) and ``"turn_deg"`` for flyby nodes. Verdict ``PASS``, ``SIBLING``
    (V_inf and epochs pass, a periapsis or turn fails) or ``FAIL``.
    """
    tol = tol or IdentityTolerances()
    pts = sorted(ladder, key=lambda p: p["sigma"])
    out: dict[str, Any] = {"verdict": "FAIL", "checks": []}
    if len(pts) < 2:
        out["reason"] = "fewer than two converged ladder points"
        return out
    lo, hi = pts[0], pts[1]
    if abs(hi["sigma"] / lo["sigma"] - 2.0) > 1e-6 or lo["sigma"] > tol.max_last_sigma:
        out["reason"] = (
            f"ladder ends at sigma {lo['sigma']} (need a halving pair, <= {tol.max_last_sigma})"
        )
        return out
    fly = [r for r in ref if r["kind"] == "flyby"]
    t_fly = float(np.median([r["rp_km"] / r["vinf_kms"] for r in fly])) if fly else 3600.0
    tol_t_days = tol.epoch_frac_of_flyby_time * t_fly / 86400.0
    ok = {"vinf": True, "epoch": True, "rp": True, "turn": True, "ratio": True}
    for k, r in enumerate(ref):
        qs: list[tuple[str, float, float]] = [
            ("vinf", r["vinf_kms"], tol.vinf_rel * r["vinf_kms"]),
            ("epoch", r["t_days"], tol_t_days),
        ]
        if r["kind"] == "flyby":
            qs += [
                ("rp", r["rp_km"], tol.rp_rel * r["rp_km"]),
                ("turn", r["turn_deg"], max(tol.turn_abs_deg, tol.turn_rel * r["turn_deg"])),
            ]
        for name, refv, tq in qs:
            field_name = {"vinf": "vinf_kms", "epoch": "t_days", "rp": "rp_over_sigma"}.get(
                name, "turn_deg"
            )
            q_lo, q_hi = float(lo["nodes"][k][field_name]), float(hi["nodes"][k][field_name])
            lim = 2.0 * q_lo - q_hi
            row: dict[str, Any] = {
                "node": k,
                "q": name,
                "ref": refv,
                "limit": lim,
                "last": q_lo,
                "tol": tq,
                "limit_ok": abs(lim - refv) <= tq,
                "last_ok": abs(q_lo - refv) <= 2.0 * tq,
            }
            if name == "vinf":
                d_hi, d_lo = abs(q_hi - refv), abs(q_lo - refv)
                small = d_hi < 0.1 * tq and d_lo < 0.1 * tq
                ratio = d_hi / d_lo if d_lo > 0 else math.inf
                row["ratio"] = ratio
                row["ratio_ok"] = small or 1.4 <= ratio <= 2.8
                ok["ratio"] &= bool(row["ratio_ok"])
            ok[name] &= bool(row["limit_ok"] and row["last_ok"])
            out["checks"].append(row)
    out["ok"] = ok
    out["sigma_pair"] = [lo["sigma"], hi["sigma"]]
    out["epoch_tol_s"] = tol_t_days * 86400.0
    if ok["vinf"] and ok["epoch"] and ok["ratio"]:
        out["verdict"] = "PASS" if ok["rp"] and ok["turn"] else "SIBLING"
    return out
