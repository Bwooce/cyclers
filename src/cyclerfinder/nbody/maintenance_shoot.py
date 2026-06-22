"""M7 leg targeting — the n-body fixed-time position-targeting Newton solver.

This is the inner kernel of the M7 horizon-TCM measurement
(``docs/notes/2026-06-23-m7-implementation-plan.md``). Given a departure position
``r0`` at epoch ``t0`` and a required arrival position ``r_target`` at ``t1``, it
solves for the departure velocity ``v_dep`` whose **n-body** (typically
Mars-perturbed) propagation arrives at ``r_target``. This is the n-body analogue of
a fixed-time Lambert solve: it is exactly the "B-plane targeting" of an M7 flyby
expressed as position targeting (aiming the flyby = choosing the post-flyby velocity
= hitting the next encounter planet).

It is a Newton iteration on the propagator's co-integrated state-transition matrix:
the residual is ``r_propagated(t1) - r_target`` and the Jacobian is
``d r_f / d v_0 = STM[0:3, 3:6]`` (the upper-right 3x3 block of the 6x6 STM returned
by :meth:`cyclerfinder.nbody.propagator.RestrictedNBody.propagate` with
``with_stm=True``). The STM is co-integrated (REBOUND variational particles with the
Mars gravity-gradient tensor applied — see ``reference_rebound_variation_custom_force_gotcha``),
so the Jacobian is correct under the perturber, not just Sun-only.

Honesty discipline: if the propagation diverges (``converged=False``) or the Newton
step stalls, the solver returns ``converged=False`` and the caller leaves that leg /
row unmeasured (off-anchor / high-V∞ families stay V0, never forced).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.nbody.forces import RailsEphemerisCache
from cyclerfinder.nbody.propagator import RestrictedNBody

Vec3 = NDArray[np.float64]

# Default convergence tolerance on the arrival miss. 1 km is far inside any flyby
# SOI band (Mars 3-SOI ~ 1.7e6 km) and well below the patched-conic miss the chain
# already records as evidence; tightening past ~1 km buys nothing physically and
# costs Newton iterations against integrator noise.
DEFAULT_TOL_KM = 1.0
DEFAULT_MAX_ITER = 30


@dataclass(frozen=True)
class TargetLegResult:
    """Outcome of one :func:`target_leg` solve.

    ``v_dep_km_s`` is the departure velocity whose n-body propagation arrives at
    ``r_target`` (only meaningful when ``converged``). ``miss_km`` is the final
    arrival miss; ``n_iter`` the Newton steps taken; ``converged`` the honesty
    flag (False => propagation diverged or Newton stalled; do not use the velocity).
    """

    v_dep_km_s: Vec3
    miss_km: float
    converged: bool
    n_iter: int


def target_leg(
    prop: RestrictedNBody,
    r0_km: Vec3,
    t0_sec: float,
    t1_sec: float,
    r_target_km: Vec3,
    v_guess_km_s: Vec3,
    *,
    bodies: tuple[str, ...] = ("M",),
    ephem: Ephemeris | None = None,
    cache: RailsEphemerisCache | None = None,
    tol_km: float = DEFAULT_TOL_KM,
    max_iter: int = DEFAULT_MAX_ITER,
    accuracy: float = 1e-10,
    max_wall_sec: float = 90.0,
) -> TargetLegResult:
    """Solve for the departure velocity whose n-body arc hits ``r_target_km`` at ``t1``.

    Newton on the propagator STM: at each step propagate ``(r0, v)`` to ``t1`` with
    ``with_stm=True``, form the residual ``r_f - r_target`` and the Jacobian
    ``J = STM[0:3, 3:6]``, and take ``v += solve(J, -residual)``. Converges when the
    arrival miss falls below ``tol_km``.

    Parameters
    ----------
    prop:
        The :class:`RestrictedNBody` propagator instance.
    r0_km, t0_sec, t1_sec, r_target_km:
        Leg boundary conditions (departure position + epoch, arrival epoch + target
        position), km and seconds.
    v_guess_km_s:
        Initial departure-velocity guess (e.g. the Lambert / sourced cycler value).
    bodies:
        Perturber body codes for the propagation (``("M",)`` = Mars-perturbed;
        ``()`` = Sun-only two-body, used by the golden cross-check vs ``lambert``).
    ephem:
        Ephemeris backend (required when ``bodies`` is non-empty and ``cache`` is
        None, so the rails cache can be built).
    cache:
        Optional pre-built :class:`RailsEphemerisCache` spanning ``[t0, t1]`` — the
        chain driver builds one per leg / itinerary and passes it in to avoid the
        ~0.5 s per-call spline build.
    tol_km, max_iter:
        Newton convergence tolerance (arrival miss, km) and iteration cap.
    accuracy, max_wall_sec:
        Forwarded to :meth:`RestrictedNBody.propagate`.

    Returns
    -------
    TargetLegResult
    """
    v = np.asarray(v_guess_km_s, dtype=np.float64).copy()
    r0 = np.asarray(r0_km, dtype=np.float64)
    r_target = np.asarray(r_target_km, dtype=np.float64)

    miss_km = float("inf")
    for n_iter in range(1, max_iter + 1):
        arc = prop.propagate(
            r0,
            v,
            t0_sec,
            t1_sec,
            bodies=bodies,
            ephem=ephem,
            cache=cache,
            accuracy=accuracy,
            max_wall_sec=max_wall_sec,
            with_stm=True,
        )
        if not arc.converged or arc.stm is None:
            # Propagation budgeted out / diverged — honest failure, no velocity.
            return TargetLegResult(
                v_dep_km_s=v, miss_km=float("inf"), converged=False, n_iter=n_iter
            )

        residual = np.asarray(arc.r_km, dtype=np.float64) - r_target
        miss_km = float(np.linalg.norm(residual))
        if miss_km < tol_km:
            return TargetLegResult(v_dep_km_s=v, miss_km=miss_km, converged=True, n_iter=n_iter)

        jac = np.asarray(arc.stm, dtype=np.float64)[0:3, 3:6]  # d r_f / d v_0
        try:
            step = np.linalg.solve(jac, -residual)
        except np.linalg.LinAlgError:
            # Singular sensitivity — cannot target this leg.
            return TargetLegResult(v_dep_km_s=v, miss_km=miss_km, converged=False, n_iter=n_iter)
        v = v + step

    # Exhausted iterations without reaching tol_km.
    return TargetLegResult(v_dep_km_s=v, miss_km=miss_km, converged=False, n_iter=max_iter)
