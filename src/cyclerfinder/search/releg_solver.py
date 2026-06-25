"""Releg leg-solver seam — swap a ballistic Lambert leg for a powered one (#449).

A moon-tour cycle is a chain ``M0 -> M1 -> ... -> M0`` closed today by re-solving
a ballistic Lambert between consecutive moons and requiring the hyperbolic excess
speed ``V_inf`` to be *continuous* across each flyby (the unbridged ``|ΔV_inf|`` is
the closure residual, gated at 0.05 km/s — ``discovery_campaign._close_one_phasing``
lines 494-528). For the ice/gas-giant systems that residual is structurally large
(high-V_inf basin, or disjoint Tisserand contours), so the ballistic tour cannot
close. A **powered releg** introduces the missing knob: a budgeted ΔV inside the
leg that *changes* V_inf — the V_inf-leveraging maneuver (VILM) physics — turning
the feasibility wall into a budget question.

This module is the **swap target**: a single ``Releg`` protocol with backends.

* :class:`BallisticReleg` — the zero-ΔV backend. Calls ``core.lambert.lambert``
  and picks the lowest-energy branch, exactly reproducing today's ballistic leg
  (the regression lock; the ballistic leg is the ``dv_kms == 0`` limit a powered
  leg strictly subsumes).
* :class:`DsmReleg` — the primary powered backend. Wraps ``search.dsm_leg.dsm_leg``
  (one deep-space maneuver per leg, the Takao eta-coordinate / Vasile-Campagnola
  BS3 model) and optimises the DSM timing fraction ``eta`` to retarget the arrival
  V_inf to the value the next leg needs, spending ``dv_dsm_kms``. The delivered ΔV
  is golden-anchored to the analytic VILM leveraging floor
  (``search.vilm.vilm_dv_min``), which is already golden-validated against
  Campagnola-Russell Endgame Part-1 Tables 1/2.

Reuse over rebuild: no optimiser is written here and no corrector is changed; the
DSM leg solver (#307) and the VILM cost model are imported, not re-derived.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol, runtime_checkable

import numpy as np
from numpy.typing import NDArray

from cyclerfinder.core.lambert import LambertError, lambert

Vec3 = NDArray[np.float64]


@dataclass(frozen=True)
class RelegResult:
    """Outcome of one releg (powered or ballistic) leg solve.

    Attributes
    ----------
    vinf_out:
        Departure hyperbolic excess speed at moon A (km/s) — ``||v_depart -
        v_a_moon||`` on the leg's first sub-arc. ``inf`` on an infeasible leg.
    vinf_in:
        Arrival hyperbolic excess speed at moon B (km/s) — ``||v_arrive -
        v_b_moon||``. For a powered leg this is the *retargeted* arrival V_inf
        (the value the maneuver delivered). ``inf`` on an infeasible leg.
    dv_kms:
        DELIVERED ΔV inside the leg (km/s). Exactly ``0.0`` for the ballistic
        leg; the DSM impulse magnitude for the powered leg.
    feasible:
        ``True`` iff the leg solver produced a usable arc (a Lambert/DSM
        solution at the requested rev existed and, for the powered leg, the
        retarget converged within the leg's reach).
    """

    vinf_out: float
    vinf_in: float
    dv_kms: float
    feasible: bool


@runtime_checkable
class Releg(Protocol):
    """The swappable moon-tour leg-solver contract (design §2).

    A backend takes the planet-frame departure/arrival moon states, the leg
    time-of-flight, the primary ``mu`` and a revolution count, and returns the
    departure/arrival V_inf magnitudes plus the delivered in-leg ΔV. The
    ballistic backend delivers ``dv_kms == 0``; the powered backends additionally
    accept a ``vinf_target_in`` they retarget the arrival V_inf to.
    """

    def solve(
        self,
        r_a: Vec3,
        v_a: Vec3,
        r_b: Vec3,
        v_b: Vec3,
        tof_s: float,
        mu: float,
        *,
        n_rev: int = 0,
        vinf_target_in: float | None = None,
    ) -> RelegResult:
        """Solve one moon-to-moon leg; return its V_inf chain + delivered ΔV."""
        ...


def _infeasible() -> RelegResult:
    return RelegResult(
        vinf_out=float("inf"), vinf_in=float("inf"), dv_kms=float("inf"), feasible=False
    )


class BallisticReleg:
    """Zero-ΔV ballistic leg — reproduces today's ``_close_one_phasing`` leg.

    Exactly the lines 494-501 logic: solve Lambert at the requested ``n_rev``,
    keep that rev, pick the lowest-energy (smallest departure V_inf) branch, and
    read the departure/arrival V_inf off it. ``dv_kms`` is always ``0.0``;
    ``vinf_target_in`` is ignored (a ballistic leg has no knob to retarget).
    """

    def solve(
        self,
        r_a: Vec3,
        v_a: Vec3,
        r_b: Vec3,
        v_b: Vec3,
        tof_s: float,
        mu: float,
        *,
        n_rev: int = 0,
        vinf_target_in: float | None = None,
    ) -> RelegResult:
        r_a = np.asarray(r_a, dtype=np.float64)
        v_a = np.asarray(v_a, dtype=np.float64)
        r_b = np.asarray(r_b, dtype=np.float64)
        v_b = np.asarray(v_b, dtype=np.float64)
        try:
            sols = lambert(r_a, r_b, tof_s, mu=mu, max_revs=max(0, n_rev))
        except (LambertError, ValueError):
            return _infeasible()
        wanted = [s for s in sols if s.n_revs == n_rev]
        if not wanted:
            return _infeasible()
        best = min(wanted, key=lambda s: float(np.linalg.norm(s.v1 - v_a)))
        vinf_out = float(np.linalg.norm(best.v1 - v_a))
        vinf_in = float(np.linalg.norm(best.v2 - v_b))
        return RelegResult(vinf_out=vinf_out, vinf_in=vinf_in, dv_kms=0.0, feasible=True)
