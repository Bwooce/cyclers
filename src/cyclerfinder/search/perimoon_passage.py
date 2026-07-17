"""Perimoon-passage geometry for a full-period planar CR3BP symmetric orbit.

Task #627 pilot need: after mu-continuing a Ross-Roberts-Tsoukkas 2026 (k1,k2)
family member down to a real planet-moon mass ratio, "the family survives" is
NOT by itself evidence the resulting orbit is a *useful* cycler encounter — a
periodic orbit can pass arbitrarily close to the secondary (a near-collision,
unflyable) or stay so far from it that the "encounter" is dynamically no
different from a generic heliocentric pass. This module answers that
quantitative question directly from the corrected orbit's IC: propagate one
full period, find the closest approach to the secondary (the moon), and report
the periapsis altitude + local relative speed there so the existing
:mod:`cyclerfinder.search.physical_sanity` bend gate can judge whether the
passage is a genuine gravity-assist-scale encounter.

The "local V_inf" used here is the ROTATING-FRAME relative speed at closest
approach: since the secondary is stationary in the CR3BP rotating frame, the
spacecraft's rotating-frame velocity AT the secondary's location *is* its
velocity relative to the secondary (exact, not an approximation) -- this is
the same frame convention the rest of the (k1,k2) genome uses (Jacobi
constant, Barden stability, winding topology all live in the rotating frame).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import solve_ivp
from scipy.optimize import minimize_scalar

import cyclerfinder.core.cr3bp as cr3bp


@dataclass(frozen=True)
class PerimoonPassage:
    """Closest approach to the secondary body over one full orbital period."""

    t_periapsis: float  # nondim time of closest approach (within [0, period])
    r2_nd: float  # nondim distance to secondary at closest approach
    r2_km: float  # dimensional distance to secondary, km
    altitude_km: float  # r2_km - secondary_radius_km (negative == below surface)
    speed_rel_nd: float  # nondim rotating-frame speed at closest approach
    speed_rel_kms: float  # dimensional rotating-frame relative speed, km/s
    below_surface: bool  # altitude_km < 0 (a literal collision, unphysical)


def find_perimoon_passage(
    system: cr3bp.CR3BPSystem,
    state0: NDArray[np.float64],
    period: float,
    secondary_radius_km: float,
    *,
    n_coarse: int = 4000,
    rtol: float = 1e-12,
    atol: float = 1e-12,
) -> PerimoonPassage:
    """Find the closest approach to the secondary over one full period.

    Propagates ``state0`` for one full ``period`` (DOP853, dense output),
    coarse-samples ``n_coarse`` points to bracket the global minimum distance
    to the secondary (fixed at ``(1-mu, 0)`` in the rotating frame), then
    refines with a bounded local minimization. Returns both the nondim and
    dimensional (km, km/s) closest-approach geometry.
    """
    mu = system.mu
    sec = np.array([1.0 - mu, 0.0])

    sol = solve_ivp(
        cr3bp.cr3bp_eom,
        (0.0, float(period)),
        np.asarray(state0, dtype=np.float64),
        args=(mu,),
        method="DOP853",
        rtol=rtol,
        atol=atol,
        dense_output=True,
    )
    if not sol.success:
        raise RuntimeError(f"perimoon propagation failed: {sol.message}")
    # scipy's stubs type ``OdeResult.sol`` as ``Optional[OdeSolution]`` because
    # it depends on the runtime ``dense_output`` flag; we always pass
    # ``dense_output=True`` above (and just confirmed the solve succeeded),
    # so a dense interpolant is guaranteed to exist here. Bind it to a
    # concretely-typed local so mypy can follow the non-None narrowing into
    # the ``_r2_at`` closure below (narrowing of a captured outer variable
    # does not propagate into nested functions).
    assert sol.sol is not None
    dense_sol = sol.sol

    ts = np.linspace(0.0, float(period), n_coarse)
    ys = dense_sol(ts)
    dx = ys[0] - sec[0]
    dy = ys[1] - sec[1]
    r2 = np.sqrt(dx * dx + dy * dy)
    i0 = int(np.argmin(r2))

    # Bounded local refinement around the coarse minimum.
    lo = ts[max(i0 - 2, 0)]
    hi = ts[min(i0 + 2, n_coarse - 1)]

    def _r2_at(t: float) -> float:
        y = dense_sol(t)
        ddx = y[0] - sec[0]
        ddy = y[1] - sec[1]
        return float(math.sqrt(ddx * ddx + ddy * ddy))

    if hi > lo:
        res = minimize_scalar(_r2_at, bounds=(lo, hi), method="bounded")
        t_star = float(res.x)
    else:
        t_star = float(ts[i0])

    y_star = dense_sol(t_star)
    ddx = float(y_star[0] - sec[0])
    ddy = float(y_star[1] - sec[1])
    r2_star = math.sqrt(ddx * ddx + ddy * ddy)
    speed_nd = float(math.hypot(y_star[3], y_star[4]))

    r2_km = r2_star * system.l_km
    speed_kms = speed_nd * system.l_km / system.t_s
    altitude_km = r2_km - secondary_radius_km

    return PerimoonPassage(
        t_periapsis=t_star,
        r2_nd=r2_star,
        r2_km=r2_km,
        altitude_km=altitude_km,
        speed_rel_nd=speed_nd,
        speed_rel_kms=speed_kms,
        below_surface=altitude_km < 0.0,
    )
