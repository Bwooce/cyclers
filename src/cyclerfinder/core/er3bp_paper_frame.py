"""ER3BP dynamics in the Antoniadou & Libert (2018) NON-pulsating rotating frame.

Reference
---------
Antoniadou, K.I. & Libert, A.-S. (2018), "Spatial resonant periodic orbits in
the restricted three-body problem", arXiv:1805.00288 (MNRAS 2018). The dynamics
of the massless body P1 are governed by their Lagrangian **Eq. 1**, integrated
in **TIME t** (not true anomaly) in a rotating-but-NON-pulsating frame Oxy. Both
primaries physically slide on the Ox axis along their mutual Keplerian ellipse
(separation ``r(theta)`` oscillates with eccentricity ``e2``); the star P0 sits
at ``x = -mu*r`` and the planet P2 at ``x = (1-mu)*r``.

This is a DIFFERENT frame from :mod:`cyclerfinder.core.er3bp` (the Szebehely
Ch. 10 pulsating / Nechvile frame with the independent variable being the true
anomaly ``f`` and both primaries fixed). Here the published ``(e1, e2)``
configurations of the paper are the natural seed coordinates, with ``e1`` the
osculating heliocentric eccentricity of P1 and ``e2`` the primaries' mutual
eccentricity.

Provenance of the paper-frame EOM (all from ``d/dt(dL/dqdot) - dL/dq = 0`` on
Eq. 1, planar slice):

* ``+thetaddot*y / -thetaddot*x`` : Euler term (frame angular ACCELERATION;
  ZERO in CR3BP since ``thetadot`` is constant) -- the genuinely new ER3BP term.
* ``+2*thetadot*vy / -2*thetadot*vx`` : Coriolis (reduces to +/-2 at e2=0).
* ``+thetadot^2 * x / +thetadot^2 * y`` : centrifugal (reduces to +x, +y).
* ``-g_x / -g_y`` : gravity from the two MOVING primaries.

At ``e2 = 0`` (``r=1``, ``thetadot=1``, ``thetaddot=0``) the EOM reduces EXACTLY
to the planar slice of :func:`cyclerfinder.core.cr3bp.cr3bp_eom` (verified to
machine precision by :func:`tests/core/test_er3bp_paper_frame`).

This module promotes the #442 capability work: the validated Eq. 1 EOM plus the
4-vector full-period least_squares corrector that tracks the CONNECTED 3/1 (pi,0)
resonant family from ``e2=0`` to ``e2=0.90``.
"""

from __future__ import annotations

import math

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares

__all__ = [
    "correct_resonant_member",
    "osculating_a1",
    "paper_frame_eom",
    "primary_kepler",
]


def primary_kepler(theta: float, e2: float) -> tuple[float, float, float, float]:
    """Return ``(r, rdot, thetadot, thetaddot)`` of the primaries' mutual orbit.

    ``theta`` is the primaries' true anomaly. The mutual orbit has ``a2 = 1`` and
    ``G(m0+m2) = 1`` so the specific angular momentum is
    ``h = sqrt(a2*(1-e2^2)) = sqrt(1-e2^2)``.
    """
    one_m_e2sq = 1.0 - e2 * e2
    r = one_m_e2sq / (1.0 + e2 * math.cos(theta))
    h = math.sqrt(one_m_e2sq)  # = r^2 * thetadot
    thetadot = h / (r * r)
    # rdot = dr/dtheta * thetadot ; dr/dtheta = r^2*e2*sin/(1-e2^2)
    rdot = (r * r) * e2 * math.sin(theta) / one_m_e2sq * thetadot
    # thetaddot from d/dt(r^2 thetadot) = 0 -> 2 r rdot thetadot + r^2 thetaddot = 0
    thetaddot = -2.0 * rdot * thetadot / r
    return r, rdot, thetadot, thetaddot


def paper_frame_eom(
    t: float, state5: NDArray[np.float64], mu: float, e2: float
) -> NDArray[np.float64]:
    """Antoniadou & Libert (2018) Eq. 1 EOM in the non-pulsating rotating frame.

    Independent variable is TIME ``t``. The state is
    ``[x, y, vx, vy, theta]`` where ``theta`` is the primaries' true anomaly
    (carried as a 5th coordinate so the augmented system is autonomous and any
    ``t``-grid integrates cleanly). Returns ``d/dt`` of the state.

    The primaries slide on the Ox axis: the star P0 (mass ``1-mu``) at
    ``x = -mu*r`` and the planet P2 (mass ``mu``) at ``x = (1-mu)*r``, with
    ``r = r(theta)`` the current mutual separation.

    Reference: arXiv:1805.00288, Eq. 1. See the module docstring for term
    provenance and the ``e2=0`` CR3BP reduction.
    """
    x, y, vx, vy, theta = (float(v) for v in state5)
    r, _rdot, thetadot, thetaddot = primary_kepler(theta, e2)

    # Moving primaries on the Ox axis:
    x_star = -mu * r  # P0 (mass 1-mu) at  -mu*r
    x_plan = (1.0 - mu) * r  # P2 (mass mu)   at (1-mu)*r
    r1 = math.sqrt((x - x_star) ** 2 + y * y)  # = sqrt((x+mu r)^2 + y^2)
    r2 = math.sqrt((x - x_plan) ** 2 + y * y)  # = sqrt((x-(1-mu) r)^2 + y^2)
    r1_3 = r1**3
    r2_3 = r2**3

    g_x = (1.0 - mu) * (x - x_star) / r1_3 + mu * (x - x_plan) / r2_3
    g_y = (1.0 - mu) * y / r1_3 + mu * y / r2_3

    ax = thetaddot * y + 2.0 * thetadot * vy + thetadot * thetadot * x - g_x
    ay = -thetaddot * x - 2.0 * thetadot * vx + thetadot * thetadot * y - g_y

    return np.array([vx, vy, ax, ay, thetadot], dtype=np.float64)


def osculating_a1(x: float, vy: float, theta0: float, mu: float, e2: float) -> float:
    """Osculating heliocentric (about star P0) semi-major axis at a crossing IC.

    For a perpendicular x-axis crossing ``[x, 0, 0, vy, theta0]`` the inertial
    tangential speed about the star is ``vy + thetadot*x`` (de-rotation), and the
    heliocentric radius is ``|x - x_star|`` with ``x_star = -mu*r``. Vis-viva
    about the star (GM = ``1-mu``) then yields ``a1``. This is the quantity the
    paper plots as ``a1/a2`` in its DS-maps.
    """
    r_prim, _rdot, thetadot, _thetaddot = primary_kepler(theta0, e2)
    x_star = -mu * r_prim
    r1 = abs(x - x_star)
    v_inert = vy + thetadot * x
    inv_a = 2.0 / r1 - v_inert * v_inert / (1.0 - mu)
    if abs(inv_a) < 1e-14:
        return float("inf")
    return 1.0 / inv_a


def _period_residual(
    free: NDArray[np.float64],
    theta0: float,
    mu: float,
    e2: float,
    period: float,
    rtol: float,
    atol: float,
) -> NDArray[np.float64]:
    """4-vector full-period periodicity residual for the (x, vy) free variables.

    Residual ``[x(T)-x, y(T), vx(T), vy(T)-vy]`` for the doubly-symmetric ansatz
    ``[x, y=0, vx=0, vy, theta0]`` integrated over one period ``T``. Returns a
    large constant vector on integrator blow-up so the optimiser steers away.
    """
    x, vy = float(free[0]), float(free[1])
    if not (abs(x) < 50.0 and abs(vy) < 50.0):
        return np.array([1e3, 1e3, 1e3, 1e3])
    s0 = np.array([x, 0.0, 0.0, vy, theta0], dtype=np.float64)
    sol = solve_ivp(
        paper_frame_eom,
        (0.0, period),
        s0,
        args=(mu, e2),
        method="DOP853",
        rtol=rtol,
        atol=atol,
    )
    if not sol.success:
        return np.array([1e3, 1e3, 1e3, 1e3])
    sf = sol.y[:, -1]
    return np.array([sf[0] - x, sf[1], sf[2], sf[3] - vy])


def correct_resonant_member(
    x0: float,
    vy0: float,
    *,
    theta0: float,
    mu: float,
    e2: float,
    period: float = 2.0 * math.pi,
    rtol: float = 1e-12,
    atol: float = 1e-12,
    max_nfev: int = 400,
) -> dict[str, float | str]:
    """Converge a doubly-symmetric resonant member in the paper frame.

    Drives the 4-vector full-period periodicity objective
    ``[x(T)-x, y(T), vx(T), vy(T)-vy]`` to zero with ``scipy.least_squares``,
    trying the ``lm`` / ``trf`` / ``dogbox`` methods and keeping the best. The IC
    is the perpendicular crossing ``[x, y=0, vx=0, vy, theta0]`` with period
    ``T = period`` (``2*pi`` for the 3/1 (pi,0) family).

    This is the corrector that (per #442) replaced the stalling 2-var half-period
    damped-Newton and tracked the CONNECTED 3/1 (pi,0) family cleanly from
    ``e2=0`` to ``e2=0.90``.

    Returns a dict with the converged ``x``, ``vy``, the ``method`` that won, and
    ``residual`` -- the INDEPENDENT full-period closure residual
    ``max|state(T) - state(0)|`` recomputed from the converged IC (not the
    optimiser's internal cost), so it is a faithful closure gauge.
    """
    best: dict[str, float | str] | None = None
    for method in ("lm", "trf", "dogbox"):
        try:
            res = least_squares(
                _period_residual,
                x0=np.array([x0, vy0], dtype=np.float64),
                args=(theta0, mu, e2, period, rtol, atol),
                method=method,
                xtol=1e-15,
                ftol=1e-15,
                gtol=1e-15,
                max_nfev=max_nfev,
            )
        except Exception:
            continue
        xc, vyc = float(res.x[0]), float(res.x[1])
        # Independent re-integration of the converged IC (faithful closure gauge).
        indep = _period_residual(np.array([xc, vyc]), theta0, mu, e2, period, rtol, atol)
        residual = float(np.max(np.abs(indep)))
        cand: dict[str, float | str] = {
            "x": xc,
            "vy": vyc,
            "residual": residual,
            "method": method,
        }
        if best is None or residual < float(best["residual"]):
            best = cand
    if best is None:
        raise RuntimeError(
            f"correct_resonant_member: all least_squares methods failed (e2={e2}, theta0={theta0})"
        )
    return best
