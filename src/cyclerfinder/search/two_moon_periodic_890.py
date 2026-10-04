"""Two-moon symmetric periodic orbits of the planar CCR4BP (#890).

Question: is the one symmetric Titania-Oberon-Titania closure that passes the
#888 demanded-turn gate (ideal circular-coplanar patched conic) the shadow of a
real periodic orbit of the restricted four-body problem in which BOTH moons are
massive? This module supplies the model wrapper, the frame conversions, a
symmetric multiple-shooting corrector, continuation in a moon-mass scale, and
encounter diagnostics. It reuses :mod:`cyclerfinder.core.ccr4bp` and the
vectorised planar right-hand side of :mod:`cyclerfinder.search.ccr4bp_strob_connection`
unmodified.

Model and units
---------------
* Planar CCR4BP, base moon ``B`` (Titania) as the second primary of the rotating
  frame, perturber ``P`` (Oberon) on a concentric circle (no moon-moon force).
* Mass scale ``lam`` multiplies both moon GMs. The planet GM is
  ``GM_sys - lam (GM_B + GM_P)``, so the total stays the registry system GM
  (``core.satellites.PRIMARIES``, which already contains the inner moons, which
  orbit inside the spacecraft and act on it approximately as central mass). At
  ``lam = 0`` the model is the two-body problem with GM_sys and the moons move at
  exactly the patched-conic mean motions used by the #888 enumeration.
* Nondimensional: length ``L = a_B``, time ``1/n_B`` with
  ``n_B = sqrt((GM_planet + GM_B)/L^3)``, frame rate 1. Planet at ``(-mu, 0)``,
  base moon at ``(1 - mu, 0)``, perturber at angle ``theta0 + omega tau`` on a
  circle of radius ``a_P/a_B`` about the origin (``ccr4bp`` convention,
  ``omega = two_body_synodic_rate``). The forcing period is ``2 pi/|omega|``.

Inertial frame
--------------
"Inertial" below means the non-rotating frame centred on the planet-base-moon
barycentre whose x-axis is the base moon's direction at ``tau = 0``. Rotating
and inertial axes coincide at ``tau = 0``. Planet-centred vectors (the patched
conic's frame) are shifted by the planet's barycentric position, which is
zero at ``lam = 0``.

Symmetry
--------
With ``theta0 in {0, pi}`` the planar equations are invariant under
``(tau, x, y, vx, vy) -> (-tau, x, -y, -vx, vy)`` (both moons lie on the x-axis
at ``tau = 0``). If additionally the perturber is on the x-axis at ``tau = T``
(here ``T = 2.5`` forcing periods, where its angle is ``-5 pi = pi`` mod
``2 pi``), the same reflection about ``tau = T`` is a symmetry. An orbit that
crosses the x-axis perpendicularly at ``tau = 0`` and at ``tau = T`` is then
periodic with period ``2 T`` (five forcing periods): a fixed point of the
five-period stroboscopic map.
"""

from __future__ import annotations

import math
from collections.abc import Callable, Sequence
from dataclasses import dataclass, field
from typing import Any

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.integrate import solve_ivp
from scipy.optimize import minimize_scalar

from cyclerfinder.core.ccr4bp import CCR4BPSystem, two_body_synodic_rate
from cyclerfinder.core.satellites import PRIMARIES, SATELLITES
from cyclerfinder.search.ccr4bp_strob_connection import planar_jacobian, planar_rhs_batch

FloatArray = NDArray[np.float64]
DAY_S = 86400.0


# ---------------------------------------------------------------------------
# Model
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class TwoMoonModel:
    """Planar CCR4BP for a planet, a base moon and a perturbing moon, mass scale ``lam``."""

    primary: str = "Uranus"
    base: str = "Titania"
    pert: str = "Oberon"
    lam: float = 1.0
    theta0: float = 0.0

    @property
    def gm_sys(self) -> float:
        return float(PRIMARIES[self.primary])

    @property
    def gm_base(self) -> float:
        return self.lam * SATELLITES[self.base].mu_km3_s2

    @property
    def gm_pert(self) -> float:
        return self.lam * SATELLITES[self.pert].mu_km3_s2

    @property
    def gm_planet(self) -> float:
        return self.gm_sys - self.gm_base - self.gm_pert

    @property
    def length_km(self) -> float:
        return float(SATELLITES[self.base].sma_km)

    @property
    def n_base(self) -> float:
        """Base-moon mean motion, rad/s (the frame rate)."""
        return math.sqrt((self.gm_planet + self.gm_base) / self.length_km**3)

    @property
    def time_unit_s(self) -> float:
        return 1.0 / self.n_base

    @property
    def vel_unit_kms(self) -> float:
        return self.length_km * self.n_base

    @property
    def system(self) -> CCR4BPSystem:
        denom = self.gm_planet + self.gm_base
        mu = self.gm_base / denom
        mu_gan = self.gm_pert / denom
        a_gan = SATELLITES[self.pert].sma_km / self.length_km
        return CCR4BPSystem(
            mu=mu,
            mu_gan=mu_gan,
            a_gan=a_gan,
            omega_gan=two_body_synodic_rate(mu, mu_gan, a_gan),
            theta_gan0=self.theta0,
        )

    @property
    def forcing_period(self) -> float:
        """Perturber synodic period in the rotating frame (nondim)."""
        return self.system.ganymede_synodic_period

    def days(self, tau: float) -> float:
        return tau * self.time_unit_s / DAY_S

    def tau(self, days: float) -> float:
        return days * DAY_S / self.time_unit_s

    # -- body positions -----------------------------------------------------

    def moon_rot(self, body: str, tau: float) -> tuple[FloatArray, FloatArray]:
        """Position and velocity (nondim) of a moon in the ROTATING frame."""
        sy = self.system
        if body == self.base:
            return np.array([1.0 - sy.mu, 0.0]), np.zeros(2)
        if body == self.pert:
            th = sy.theta_gan0 + sy.omega_gan * tau
            r = sy.a_gan * np.array([math.cos(th), math.sin(th)])
            return r, sy.omega_gan * np.array([-r[1], r[0]])
        if body == self.primary:
            return np.array([-sy.mu, 0.0]), np.zeros(2)
        raise KeyError(body)

    def moon_inertial(self, body: str, tau: float) -> tuple[FloatArray, FloatArray]:
        """Barycentric inertial position (km) and velocity (km/s) of a body."""
        r, v = self.moon_rot(body, tau)
        return inertial_from_rot(self, tau, np.concatenate([r, v]))


def _rot(theta: float) -> FloatArray:
    c, s = math.cos(theta), math.sin(theta)
    return np.array([[c, -s], [s, c]])


def rot_from_inertial(
    model: TwoMoonModel, tau: float, r_km: ArrayLike, v_kms: ArrayLike
) -> FloatArray:
    """Barycentric inertial (km, km/s, planar) -> rotating nondim 4-state at time ``tau``."""
    r = np.asarray(r_km, dtype=np.float64)[:2] / model.length_km
    v = np.asarray(v_kms, dtype=np.float64)[:2] / model.vel_unit_kms
    rt = _rot(-tau)
    x = rt @ r
    # v_rot = R(-tau) (v - omega x r), omega = 1 (nondim)
    vr = rt @ (v - np.array([-r[1], r[0]]))
    return np.array([x[0], x[1], vr[0], vr[1]])


def inertial_from_rot(
    model: TwoMoonModel, tau: float, s4: ArrayLike
) -> tuple[FloatArray, FloatArray]:
    """Rotating nondim 4-state -> barycentric inertial position (km), velocity (km/s)."""
    s = np.asarray(s4, dtype=np.float64)
    rr = _rot(tau)
    r = rr @ s[:2]
    v = rr @ (s[2:4] + np.array([-s[1], s[0]]))
    return r * model.length_km, v * model.vel_unit_kms


def planet_offset_inertial(model: TwoMoonModel, tau: float) -> tuple[FloatArray, FloatArray]:
    """Planet's barycentric inertial position and velocity (km, km/s)."""
    return model.moon_inertial(model.primary, tau)


def mirror(s4: ArrayLike) -> FloatArray:
    """The reflection ``(x, y, vx, vy) -> (x, -y, -vx, vy)`` (time reversal implied)."""
    s = np.asarray(s4, dtype=np.float64)
    return np.array([s[0], -s[1], -s[2], s[3]])


# ---------------------------------------------------------------------------
# Propagation
# ---------------------------------------------------------------------------


@dataclass
class Arc:
    state: FloatArray
    stm: FloatArray | None
    sol: Any


def propagate(
    model: TwoMoonModel,
    s4: ArrayLike,
    tau0: float,
    tau1: float,
    *,
    with_stm: bool = False,
    method: str = "DOP853",
    rtol: float = 1e-13,
    atol: float = 1e-13,
    dense: bool = False,
) -> Arc:
    """Integrate the planar CCR4BP from ``tau0`` to ``tau1`` (either direction).

    ``method="Radau"`` (no STM) uses the analytic planar Jacobian.
    """
    sy = model.system
    s = np.asarray(s4, dtype=np.float64)
    if with_stm:
        if method != "DOP853":
            raise ValueError("STM propagation is DOP853 only")
        y0 = np.concatenate([s, np.eye(4).reshape(-1)])
        sol = solve_ivp(
            planar_rhs_batch,
            (tau0, tau1),
            y0,
            args=(sy, 1, True),
            method="DOP853",
            rtol=rtol,
            atol=atol,
            dense_output=dense,
        )
    else:

        def rhs(t: float, y: FloatArray) -> FloatArray:
            return planar_rhs_batch(t, y, sy, 1, False)

        if method == "Radau":

            def jac(t: float, y: FloatArray) -> FloatArray:
                return planar_jacobian(t, y, sy)

            sol = solve_ivp(
                rhs, (tau0, tau1), s, method="Radau", rtol=rtol, atol=atol,
                dense_output=dense, jac=jac,
            )  # fmt: skip
        else:
            sol = solve_ivp(
                rhs, (tau0, tau1), s, method="DOP853", rtol=rtol, atol=atol, dense_output=dense
            )
    if not sol.success:
        raise RuntimeError(f"propagation failed at tau={sol.t[-1]}: {sol.message}")
    yf = sol.y[:, -1]
    return Arc(state=yf[:4].copy(), stm=yf[4:].reshape(4, 4).copy() if with_stm else None, sol=sol)


# ---------------------------------------------------------------------------
# Kepler helper (patched-conic guesses and the two-body reduction)
# ---------------------------------------------------------------------------


def kepler_propagate(r0: FloatArray, v0: FloatArray, dt_s: float, gm: float) -> FloatArray:
    """Two-body propagation (planar, km, km/s) by numerical integration at tight tolerance."""

    def rhs(_t: float, y: FloatArray) -> FloatArray:
        r = y[:2]
        return np.concatenate([y[2:], -gm * r / float(np.linalg.norm(r)) ** 3])

    if dt_s == 0.0:
        return np.concatenate([r0[:2], v0[:2]])
    sol = solve_ivp(
        rhs,
        (0.0, dt_s),
        np.concatenate([r0[:2], v0[:2]]),
        method="DOP853",
        rtol=1e-13,
        atol=1e-9,
    )
    return np.asarray(sol.y[:, -1], dtype=np.float64)


# ---------------------------------------------------------------------------
# Patched-conic flyby hyperbola
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Hyperbola:
    """Moon-centric hyperbola that turns ``v_in`` onto ``v_out`` (planar, km, km/s)."""

    gm: float
    vinf: float
    turn_rad: float
    ecc: float
    rp_km: float
    vp_kms: float
    rp_vec: FloatArray
    vp_vec: FloatArray


def flyby_hyperbola(v_in: ArrayLike, v_out: ArrayLike, gm: float) -> Hyperbola:
    """Hyperbola of an unpowered flyby turning ``v_in`` onto ``v_out``.

    Periapsis direction ``(v_in_hat - v_out_hat)`` (the impulse points at the
    moon), periapsis velocity along ``v_in_hat + v_out_hat``. Uses the mean of
    the two magnitudes.
    """
    a = np.asarray(v_in, dtype=np.float64)[:2]
    b = np.asarray(v_out, dtype=np.float64)[:2]
    va, vb = float(np.linalg.norm(a)), float(np.linalg.norm(b))
    vinf = 0.5 * (va + vb)
    ah, bh = a / va, b / vb
    turn = math.acos(max(-1.0, min(1.0, float(ah @ bh))))
    ecc = 1.0 / math.sin(0.5 * turn)
    rp = gm * (ecc - 1.0) / vinf**2
    vp = math.sqrt(vinf**2 + 2.0 * gm / rp)
    d = ah - bh
    w = ah + bh
    return Hyperbola(
        gm=gm,
        vinf=vinf,
        turn_rad=turn,
        ecc=ecc,
        rp_km=rp,
        vp_kms=vp,
        rp_vec=rp * d / float(np.linalg.norm(d)),
        vp_vec=vp * w / float(np.linalg.norm(w)),
    )


# ---------------------------------------------------------------------------
# Symmetric multiple shooting
# ---------------------------------------------------------------------------


@dataclass
class ShootingResult:
    z: FloatArray
    converged: bool
    residual_history: list[float]
    final_residual: float
    jac_logdet_sign: float
    node_taus: FloatArray
    message: str = ""
    extra: dict[str, Any] = field(default_factory=dict)


class SymmetricShooter:
    """Multiple shooting for an orbit perpendicular to the x-axis at ``tau = 0`` and ``T``.

    Unknowns ``z = (x0, vy0, s_1, ..., s_M, xT, vyT)`` with interior nodes
    ``s_k`` (4-states) at ``tau_k = k T/(M+1)``. Segments: node 0 forward to
    ``tau_1``; node ``k`` forward to ``tau_{k+1}``; the end node backward from
    ``T`` to ``tau_M``. Residuals are the 4-vector mismatches at
    ``tau_1..tau_M`` (4M + 4 equations, 4M + 4 unknowns).
    """

    def __init__(self, model: TwoMoonModel, half_period: float, n_interior: int) -> None:
        self.model = model
        self.T = float(half_period)
        self.M = int(n_interior)
        self.taus = np.array([k * self.T / (self.M + 1) for k in range(self.M + 2)])

    @property
    def n(self) -> int:
        return 4 * self.M + 4

    def unpack(self, z: FloatArray) -> tuple[FloatArray, list[FloatArray], FloatArray]:
        s0 = np.array([z[0], 0.0, 0.0, z[1]])
        mids = [z[2 + 4 * k : 6 + 4 * k].copy() for k in range(self.M)]
        s_end = np.array([z[-2], 0.0, 0.0, z[-1]])
        return s0, mids, s_end

    @staticmethod
    def pack(s0: FloatArray, mids: Sequence[FloatArray], s_end: FloatArray) -> FloatArray:
        return np.concatenate([[s0[0], s0[3]], *mids, [s_end[0], s_end[3]]])

    def evaluate(
        self, z: FloatArray, *, with_jac: bool = True
    ) -> tuple[FloatArray, FloatArray | None]:
        s0, mids, s_end = self.unpack(z)
        n = self.n
        res = np.zeros(n)
        jac = np.zeros((n, n)) if with_jac else None
        tk = self.taus
        # segment 0: node 0 -> tau_1, equation block 0
        a = propagate(self.model, s0, 0.0, tk[1], with_stm=with_jac)
        res[0:4] = a.state - mids[0]
        if jac is not None:
            assert a.stm is not None
            jac[0:4, 0] = a.stm[:, 0]
            jac[0:4, 1] = a.stm[:, 3]
            jac[0:4, 2:6] = -np.eye(4)
        for k in range(1, self.M):
            a = propagate(self.model, mids[k - 1], tk[k], tk[k + 1], with_stm=with_jac)
            rows = slice(4 * k, 4 * k + 4)
            res[rows] = a.state - mids[k]
            if jac is not None:
                assert a.stm is not None
                jac[rows, 2 + 4 * (k - 1) : 6 + 4 * (k - 1)] = a.stm
                jac[rows, 2 + 4 * k : 6 + 4 * k] = -np.eye(4)
        # last: end node backward to tau_M
        a = propagate(self.model, s_end, self.T, tk[self.M], with_stm=with_jac)
        rows = slice(4 * self.M, 4 * self.M + 4)
        res[rows] = a.state - mids[self.M - 1]
        if jac is not None:
            assert a.stm is not None
            jac[rows, n - 2] = a.stm[:, 0]
            jac[rows, n - 1] = a.stm[:, 3]
            jac[rows, 2 + 4 * (self.M - 1) : 6 + 4 * (self.M - 1)] = -np.eye(4)
        return res, jac

    def newton(
        self,
        z0: FloatArray,
        *,
        tol: float = 1e-11,
        max_iter: int = 25,
        max_step: float = 0.05,
        log: Callable[[str], None] | None = None,
    ) -> ShootingResult:
        z = np.asarray(z0, dtype=np.float64).copy()
        hist: list[float] = []
        sign = float("nan")
        msg = "max_iter"
        converged = False
        for it in range(max_iter):
            try:
                res, jac = self.evaluate(z)
            except RuntimeError as exc:
                msg = f"integration failure: {exc}"
                break
            assert jac is not None
            r = float(np.max(np.abs(res)))
            hist.append(r)
            sign, _ = np.linalg.slogdet(jac)
            if log is not None:
                log(f"    newton it {it}: max|res| {r:.3e}")
            if r < tol:
                converged = True
                msg = "converged"
                break
            dz = np.linalg.solve(jac, -res)
            step = float(np.max(np.abs(dz)))
            scale = 1.0 if step <= max_step else max_step / step
            # backtracking on the residual norm
            accepted = False
            for _ in range(8):
                ztry = z + scale * dz
                try:
                    rtry, _j = self.evaluate(ztry, with_jac=False)
                except RuntimeError:
                    scale *= 0.5
                    continue
                if float(np.linalg.norm(rtry)) < float(np.linalg.norm(res)) or scale < 1e-3:
                    accepted = True
                    break
                scale *= 0.5
            if not accepted:
                msg = "line search failed"
                z = z + scale * dz
                break
            z = ztry
        return ShootingResult(
            z=z,
            converged=converged,
            residual_history=hist,
            final_residual=hist[-1] if hist else float("inf"),
            jac_logdet_sign=float(sign),
            node_taus=self.taus.copy(),
            message=msg,
        )


# ---------------------------------------------------------------------------
# Encounter diagnostics
# ---------------------------------------------------------------------------


def laplace_soi_km(primary: str, moon: str) -> float:
    """Laplace sphere of influence ``a (m/M)^(2/5)``."""
    sat = SATELLITES[moon]
    return float(sat.sma_km * (sat.mu_km3_s2 / PRIMARIES[primary]) ** 0.4)


def hill_radius_km(primary: str, moon: str) -> float:
    sat = SATELLITES[moon]
    return float(sat.sma_km * (sat.mu_km3_s2 / (3.0 * PRIMARIES[primary])) ** (1.0 / 3.0))


def relative_inertial(
    model: TwoMoonModel, body: str, tau: float, s4: FloatArray
) -> tuple[FloatArray, FloatArray]:
    """Spacecraft position/velocity relative to ``body`` in inertial axes (km, km/s)."""
    r, v = inertial_from_rot(model, tau, s4)
    rb, vb = model.moon_inertial(body, tau)
    return r - rb, v - vb


def osculating_flyby(gm: float, dr: FloatArray, dv: FloatArray) -> dict[str, float]:
    """Moon-centric osculating two-body elements (planar)."""
    r = float(np.linalg.norm(dr))
    v2 = float(dv @ dv)
    energy = 0.5 * v2 - gm / r
    h = float(dr[0] * dv[1] - dr[1] * dv[0])
    ecc = math.sqrt(max(0.0, 1.0 + 2.0 * energy * h * h / gm**2)) if gm > 0 else float("inf")
    out = {"energy": energy, "ecc": ecc, "h": h}
    if energy > 0 and gm > 0:
        out["vinf_kms"] = math.sqrt(2.0 * energy)
        out["turn_deg"] = math.degrees(2.0 * math.asin(min(1.0, 1.0 / ecc)))
    else:
        out["vinf_kms"] = float("nan")
        out["turn_deg"] = float("nan")
    return out


def osculating_asymptotes(
    gm: float, dr: FloatArray, dv: FloatArray
) -> tuple[FloatArray, FloatArray]:
    """Incoming and outgoing V-infinity vectors of the osculating hyperbola (planar)."""
    r = float(np.linalg.norm(dr))
    v2 = float(dv @ dv)
    energy = 0.5 * v2 - gm / r
    if energy <= 0:
        raise ValueError("not hyperbolic")
    vinf = math.sqrt(2.0 * energy)
    h = float(dr[0] * dv[1] - dr[1] * dv[0])
    evec = (v2 - gm / r) * dr / gm - float(dr @ dv) * dv / gm
    ecc = float(np.linalg.norm(evec))
    p_hat = evec / ecc
    q_hat = np.sign(h) * np.array([-p_hat[1], p_hat[0]])
    # true anomaly of the asymptotes: cos f = -1/e
    f = math.acos(-1.0 / ecc)

    def vel_dir(fa: float) -> FloatArray:
        # velocity direction at true anomaly fa (unit), perifocal: (-sin f, e + cos f)
        u = -math.sin(fa) * p_hat + (ecc + math.cos(fa)) * q_hat
        return np.asarray(u / float(np.linalg.norm(u)), dtype=np.float64)

    return vinf * vel_dir(-f), vinf * vel_dir(f)


def find_encounters(
    model: TwoMoonModel,
    sol: Any,
    *,
    bodies: Sequence[str],
    within_km: float,
    n_samples: int = 40000,
) -> list[dict[str, Any]]:
    """Every local minimum of the distance to each body below ``within_km``.

    ``sol`` must carry dense output. Each minimum is refined with a bounded
    scalar minimisation of the dense interpolant. Reports moon-centric
    osculating elements at the minimum (periapsis).
    """
    t0, t1 = float(sol.t[0]), float(sol.t[-1])
    ts = np.linspace(t0, t1, n_samples)
    ys = sol.sol(ts)
    out: list[dict[str, Any]] = []
    for body in bodies:
        d = np.array(
            [
                float(np.linalg.norm(relative_inertial(model, body, float(t), ys[:4, i])[0]))
                for i, t in enumerate(ts)
            ]
        )
        idx = np.nonzero((d[1:-1] <= d[:-2]) & (d[1:-1] <= d[2:]))[0] + 1
        cand: list[int] = [int(i) for i in idx]
        # endpoints count as minima too (a periapsis at the start of the arc)
        if d[0] < d[1]:
            cand.append(0)
        if d[-1] < d[-2]:
            cand.append(len(d) - 1)
        for i in cand:
            if d[i] > within_km:
                continue
            lo = float(ts[max(i - 1, 0)])
            hi = float(ts[min(i + 1, len(ts) - 1)])

            def dist(t: float, body: str = body) -> float:
                return float(np.linalg.norm(relative_inertial(model, body, t, sol.sol(t)[:4])[0]))

            if hi > lo:
                rr = minimize_scalar(
                    dist,
                    bounds=(min(lo, hi), max(lo, hi)),
                    method="bounded",
                    options={"xatol": 1e-12},
                )
                tm = float(rr.x)
            else:
                tm = float(ts[i])
            s = sol.sol(tm)[:4]
            dr, dv = relative_inertial(model, body, tm, s)
            gm = model.gm_base if body == model.base else model.gm_pert
            osc = osculating_flyby(gm, dr, dv)
            sat = SATELLITES[body]
            out.append(
                {
                    "body": body,
                    "tau": tm,
                    "t_days": model.days(tm),
                    "dist_km": float(np.linalg.norm(dr)),
                    "alt_km": float(np.linalg.norm(dr)) - sat.radius_eq_km,
                    "rel_speed_kms": float(np.linalg.norm(dv)),
                    "hill_radii": float(np.linalg.norm(dr)) / hill_radius_km(model.primary, body),
                    **{f"osc_{k}": v for k, v in osc.items()},
                }
            )
    out.sort(key=lambda e: e["tau"])
    return out


def soi_crossing_velocities(
    model: TwoMoonModel, sol: Any, body: str, tau_peri: float, radius_km: float
) -> dict[str, Any]:
    """Moon-relative inertial velocity where the arc crosses ``radius_km`` before/after periapsis.

    Searches outward from ``tau_peri`` on the dense output in steps of 0.01 TU
    (about 2 h at Titania) and refines each crossing by bisection.
    """

    def dist(t: float) -> float:
        return float(np.linalg.norm(relative_inertial(model, body, t, sol.sol(t)[:4])[0]))

    t_lo, t_hi = float(sol.t[0]), float(sol.t[-1])
    if t_lo > t_hi:
        t_lo, t_hi = t_hi, t_lo
    found: dict[str, Any] = {}
    for label, sgn in (("in", -1.0), ("out", 1.0)):
        t_prev = tau_peri
        t = tau_peri
        hit = None
        while True:
            t = t_prev + sgn * 0.01
            if t < t_lo or t > t_hi:
                break
            if dist(t) >= radius_km:
                hit = (t_prev, t)
                break
            t_prev = t
        if hit is None:
            found[label] = None
            continue
        a, b = hit
        for _ in range(60):
            m = 0.5 * (a + b)
            if dist(m) < radius_km:
                a = m
            else:
                b = m
        tm = 0.5 * (a + b)
        dr, dv = relative_inertial(model, body, tm, sol.sol(tm)[:4])
        found[label] = {"tau": tm, "t_days": model.days(tm), "dr_km": dr, "dv_kms": dv}
    return found


# ---------------------------------------------------------------------------
# Guess from the patched-conic closure
# ---------------------------------------------------------------------------


def _local_frame(r: FloatArray) -> tuple[FloatArray, FloatArray]:
    """(radial, along-track) unit vectors of a prograde circular orbit at ``r`` (planar)."""
    rh = r[:2] / float(np.linalg.norm(r[:2]))
    return rh, np.array([-rh[1], rh[0]])


def to_local(v: FloatArray, r_body: FloatArray) -> FloatArray:
    rh, th = _local_frame(r_body)
    return np.array([float(v[:2] @ rh), float(v[:2] @ th)])


def from_local(vl: FloatArray, r_body: FloatArray) -> FloatArray:
    rh, th = _local_frame(r_body)
    return np.asarray(vl[0] * rh + vl[1] * th, dtype=np.float64)


@dataclass
class ClosureGeometry:
    """Patched-conic data needed for guesses (planet-centred, km, km/s)."""

    tof_days: float
    gm_sys: float
    r_dep: FloatArray  # base moon at t = 0
    v_dep: FloatArray  # spacecraft at t = 0 (leg 0 departure)
    base_in_local: FloatArray  # closing (wrap) V-inf at the base moon, base-local axes
    base_out_local: FloatArray  # leg-0 departure V-inf, base-local axes
    pert_in_local: FloatArray  # middle encounter, perturber-local axes
    pert_out_local: FloatArray


def closure_geometry(closure: Any, primary: str = "Uranus") -> ClosureGeometry:
    """Extract the guess data from a ``turn_gate_closures.SymmetricClosure``."""
    from cyclerfinder.verify.turn_gate_closures import _moon_state

    mu = float(PRIMARIES[primary])
    tof = float(closure.tof_days)
    r0, w0 = _moon_state(closure.anchor, 0.0, 0.0, mu)
    r1, _w1 = _moon_state(closure.flyby, math.radians(closure.rel_offset_deg), tof, mu)
    r2, _w2 = _moon_state(closure.anchor, 0.0, 2.0 * tof, mu)
    vec = {k: np.asarray(v, dtype=np.float64) for k, v in closure.vectors.items()}
    return ClosureGeometry(
        tof_days=tof,
        gm_sys=mu,
        r_dep=np.asarray(r0, dtype=np.float64)[:2],
        v_dep=(np.asarray(w0, dtype=np.float64) + vec["out0"])[:2],
        base_in_local=to_local(vec["in2"], r2),
        base_out_local=to_local(vec["out0"], r0),
        pert_in_local=to_local(vec["in1"], r1),
        pert_out_local=to_local(vec["out1"], r1),
    )


def flyby_node_state(
    model: TwoMoonModel, body: str, tau: float, vin_local: FloatArray, vout_local: FloatArray
) -> tuple[FloatArray, Hyperbola]:
    """Rotating 4-state at the periapsis of the patched-conic hyperbola at ``body``."""
    rb, vb = model.moon_inertial(body, tau)
    gm = model.gm_base if body == model.base else model.gm_pert
    hyp = flyby_hyperbola(from_local(vin_local, rb), from_local(vout_local, rb), gm)
    s = rot_from_inertial(model, tau, rb + hyp.rp_vec, vb + hyp.vp_vec)
    return s, hyp


def symmetric_guess(
    model: TwoMoonModel, geom: ClosureGeometry, n_interior: int
) -> tuple[SymmetricShooter, FloatArray, dict[str, Any]]:
    """Multiple-shooting guess: hyperbola periapses at both flybys, leg-0 conic between.

    Interior node ``k`` is the leg-0 Kepler conic at the same FRACTION of the
    leg as the node is of the model's half period (the model's forcing period
    differs from the patched conic's by O(moon mass)).
    """
    half = 2.5 * model.forcing_period
    sh = SymmetricShooter(model, half, n_interior)
    s0, h0 = flyby_node_state(model, model.base, 0.0, geom.base_in_local, geom.base_out_local)
    s_end, h_end = flyby_node_state(
        model, model.pert, half, geom.pert_in_local, geom.pert_out_local
    )
    mids = []
    for k in range(1, n_interior + 1):
        frac = k / (n_interior + 1)
        st = kepler_propagate(geom.r_dep, geom.v_dep, frac * geom.tof_days * DAY_S, geom.gm_sys)
        tau_k = sh.taus[k]
        rp, vp = planet_offset_inertial(model, tau_k)
        mids.append(rot_from_inertial(model, tau_k, st[:2] + rp, st[2:] + vp))
    z = SymmetricShooter.pack(s0, mids, s_end)
    info = {
        "s0_dropped_y_vx": [float(s0[1]), float(s0[2])],
        "sT_dropped_y_vx": [float(s_end[1]), float(s_end[2])],
        "hyp_base": {"rp_km": h0.rp_km, "turn_deg": math.degrees(h0.turn_rad), "vinf": h0.vinf},
        "hyp_pert": {
            "rp_km": h_end.rp_km,
            "turn_deg": math.degrees(h_end.turn_rad),
            "vinf": h_end.vinf,
        },
    }
    return sh, z, info
