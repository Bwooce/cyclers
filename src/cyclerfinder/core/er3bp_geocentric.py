"""Earth-centered (geocentric) pulsating ER3BP frame — Gurfil & Kasdin (2002).

Reference: Gurfil, P. & Kasdin, N.J. (2002), "Niching genetic algorithms-based
characterization of geocentric orbits in the 3D elliptic restricted three-body
problem", Comput. Methods Appl. Mech. Engrg. 191, 5683-5706,
DOI 10.1016/S0045-7825(02)00481-4. Digest:
``docs/notes/2026-07-12-digest-gurfil-kasdin-2002-er3bp-geocentric-orbits.md``.

Frame (paper Fig. 1 + Eq. 9-11): origin at the *secondary* (Earth) center,
x radially outward from the Sun, y along Earth's velocity, z normal to the
ecliptic; rotating-pulsating with true anomaly ``theta`` as the independent
variable. The Sun (mass fraction ``1 - mu``) sits at ``(-1, 0, 0)``;
``mu = mu_E / (mu_E + mu_S)`` is the *Earth* (secondary) mass fraction.

This is NOT a new frame implementation: the Gurfil-Kasdin frame is exactly the
Szebehely/Nechvile barycentric pulsating frame of :mod:`cyclerfinder.core.er3bp`
translated by ``+(1 - mu)`` along x (origin moved from the barycenter to the
secondary). Substituting ``x_bary = x_geo + (1 - mu)`` into the barycentric
EOM reproduces the paper's Eq. 9-11 term-for-term, including the ``+ (1 - mu)``
indirect term. Verified numerically to machine precision (2026-07-12, and again
in ``tests/core/test_er3bp_geocentric.py`` against an independently-coded
Eq. 9-11 right-hand side).

State convention here is ``[x, y, z, x', y', z']`` (prime = d/dtheta), matching
``er3bp.py`` and the paper's Eq. 13. NOTE: the paper's own Tables 2 and 3 print
their vectors in INTERLEAVED order ``[x, x', y, y', z, z']`` despite Eq. 13 —
established by cross-checking every Table 3 IC against Table 4's r0/v0 values
(e.g. family J: r0 = 7,527,807 km matches exactly under interleaved ordering
and is off by 32% under the Eq. 13 ordering). Use
:func:`table_interleaved_to_state` when consuming the paper's printed vectors.
"""

from __future__ import annotations

import math
from typing import Literal

import numpy as np
from numpy.typing import NDArray

from cyclerfinder.core.cr3bp import StmMode
from cyclerfinder.core.er3bp import ER3BPSystem, er3bp_eom, propagate_er3bp

# Paper p. 5686: mu = mu_E / (mu_E + mu_S) for the Sun-Earth system.
MU_SUN_EARTH_GURFIL_KASDIN: float = 3.0034495182e-6

# Paper p. 5685: Earth heliocentric eccentricity.
E_SUN_EARTH_GURFIL_KASDIN: float = 0.0167

# Paper p. 5685: a = 1 AU = 1.496e8 km (the paper's own rounded value; kept
# verbatim so Table 3/4 cross-checks use the source's unit convention).
A_AU_KM_GURFIL_KASDIN: float = 1.496e8

# Mean Earth radius for the collision constraint (paper Eq. 17), normalized.
R_EARTH_KM: float = 6378.0
R_EARTH_NORM: float = R_EARTH_KM / A_AU_KM_GURFIL_KASDIN

SUN_EARTH_ER3BP = ER3BPSystem(
    mu=MU_SUN_EARTH_GURFIL_KASDIN,
    e=E_SUN_EARTH_GURFIL_KASDIN,
    primary_name="Sun",
    secondary_name="Earth",
)


def geocentric_to_barycentric(state6: NDArray[np.float64], mu: float) -> NDArray[np.float64]:
    """Translate a geocentric pulsating-frame state to the barycentric frame.

    Both primaries are fixed in the pulsating frame, so the transform is the
    constant translation ``x_bary = x_geo + (1 - mu)`` (velocities unchanged).
    """
    out = np.asarray(state6, dtype=float).copy()
    out[0] += 1.0 - mu
    return out


def barycentric_to_geocentric(state6: NDArray[np.float64], mu: float) -> NDArray[np.float64]:
    """Inverse of :func:`geocentric_to_barycentric`."""
    out = np.asarray(state6, dtype=float).copy()
    out[0] -= 1.0 - mu
    return out


def er3bp_geocentric_eom(
    f: float, state6: NDArray[np.float64], mu: float, e: float
) -> NDArray[np.float64]:
    """Gurfil-Kasdin Eq. 9-11 via the barycentric ``er3bp_eom`` + offset.

    A constant translation leaves all derivative components unchanged, so the
    geocentric EOM is the barycentric EOM evaluated at the shifted state.
    """
    return er3bp_eom(f, geocentric_to_barycentric(state6, mu), mu, e)


def propagate_er3bp_geocentric(
    state6: NDArray[np.float64],
    f_span: tuple[float, float],
    sys: ER3BPSystem = SUN_EARTH_ER3BP,
    rtol: float = 1e-12,
    atol: float = 1e-12,
    with_stm: bool = False,
    stm_mode: StmMode = "variable",
    method: Literal["RK23", "RK45", "DOP853", "Radau", "BDF", "LSODA"] = "DOP853",
) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]]:
    """Propagate a geocentric-frame state via the existing barycentric machinery.

    Same return contract as :func:`cyclerfinder.core.er3bp.propagate_er3bp`,
    with the state history translated back to the geocentric frame. (The STM
    is translation-invariant, so it needs no frame correction.)
    """
    bary0 = geocentric_to_barycentric(np.asarray(state6, dtype=float), sys.mu)
    f_eval, states, stm = propagate_er3bp(
        bary0,
        f_span,
        sys,
        rtol=rtol,
        atol=atol,
        with_stm=with_stm,
        stm_mode=stm_mode,
        method=method,
    )
    states = states.copy()
    states[0, :] -= 1.0 - sys.mu
    return f_eval, states, stm


def table_interleaved_to_state(vec6: NDArray[np.float64]) -> NDArray[np.float64]:
    """Convert a paper Table 2/3 vector ``[x, x', y, y', z, z']`` to state order.

    See the module docstring for why the paper's printed tables are interleaved.
    """
    x, xp, y, yp, z, zp = (float(v) for v in vec6)
    return np.array([x, y, z, xp, yp, zp])


def gurfil_kasdin_fitness(
    state6_geo: NDArray[np.float64],
    theta0: float,
    sys: ER3BPSystem = SUN_EARTH_ER3BP,
    n_rev: float = 1.0,
    rtol: float = 1e-9,
    atol: float = 1e-9,
    n_samples: int = 4000,
    escape_radius: float = 0.5,
) -> float:
    """Paper Eq. 15 objective: ``1 / ((rmax - rmin)^2 + 1)``, maximized.

    ``rmin``/``rmax`` are min/max normalized geocentric distance over a
    ``n_rev``-revolution (1 rev = 1 year) integration from ``theta0``.
    Constraint Eq. 17 (``rmin > R_Earth``) and escape beyond ``escape_radius``
    (normalized; guards against runaway heliocentric drift dominating the
    integrator budget) are enforced as death penalties (fitness 0), matching
    the paper's exclusion of collision/escape trajectories.
    """
    from scipy.integrate import solve_ivp

    mu, e = sys.mu, sys.e
    bary0 = geocentric_to_barycentric(np.asarray(state6_geo, dtype=float), mu)
    offset = 1.0 - mu

    def geo_r2(f: float, s: NDArray[np.float64]) -> float:
        dx = s[0] - offset
        return float(dx * dx + s[1] * s[1] + s[2] * s[2])

    def collision(f: float, s: NDArray[np.float64], *_a: float) -> float:
        return geo_r2(f, s) - R_EARTH_NORM**2

    def escape(f: float, s: NDArray[np.float64], *_a: float) -> float:
        return geo_r2(f, s) - escape_radius**2

    collision.terminal = True  # type: ignore[attr-defined]
    escape.terminal = True  # type: ignore[attr-defined]

    f_end = theta0 + n_rev * 2.0 * math.pi
    sol = solve_ivp(
        er3bp_eom,
        (theta0, f_end),
        bary0,
        args=(mu, e),
        method="DOP853",
        rtol=rtol,
        atol=atol,
        t_eval=np.linspace(theta0, f_end, n_samples),
        events=(collision, escape),
    )
    if not sol.success and sol.status != 1:
        return 0.0
    if sol.status == 1:  # terminated on collision or escape
        return 0.0
    if sol.y.shape[1] < 2:
        return 0.0

    dx = sol.y[0] - offset
    r = np.sqrt(dx * dx + sol.y[1] ** 2 + sol.y[2] ** 2)
    rmin = float(r.min())
    rmax = float(r.max())
    if rmin <= R_EARTH_NORM:
        return 0.0
    return 1.0 / ((rmax - rmin) ** 2 + 1.0)
