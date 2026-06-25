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

import math
from collections.abc import Sequence
from dataclasses import dataclass
from typing import Protocol, runtime_checkable

import numpy as np
from numpy.typing import NDArray
from scipy.optimize import minimize

from cyclerfinder.core.kepler import KeplerError
from cyclerfinder.core.lambert import LambertError, lambert
from cyclerfinder.search.dsm_leg import dsm_leg

Vec3 = NDArray[np.float64]

# DSM timing fraction is clamped this far from the singular endpoints {0, 1}
# (see ``dsm_leg._ETA_EPS``): at eta=0 the ballistic front arc is degenerate, at
# eta=1 the Lambert back arc is singular. The retarget sweep stays strictly
# inside this open interval.
_ETA_LO: float = 0.02
_ETA_HI: float = 0.98


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
        vinf_depart_mag: float | None = None,
    ) -> RelegResult:
        """Solve one moon-to-moon leg; return its V_inf chain + delivered ΔV.

        ``vinf_target_in`` retargets the arrival V_inf (powered backends only);
        ``vinf_depart_mag`` pins the departure V_inf magnitude (powered backends
        only, used by the moon-tour driver for by-construction flyby continuity).
        The ballistic backend ignores both.
        """
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
        vinf_depart_mag: float | None = None,
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


class DsmReleg:
    """Powered one-DSM leg — retargets arrival V_inf, spends a budgeted ΔV (#449).

    The primary powered backend. Wraps :func:`cyclerfinder.search.dsm_leg.dsm_leg`
    (one deep-space maneuver splitting the leg into two Lambert sub-arcs joined at
    a free interior point — the Takao eta-coordinate / Vasile-Campagnola BS3
    "one-DSM-per-leg" model). The departure velocity is seeded from the ballistic
    Lambert leg connecting the two moons; the DSM timing fraction ``eta`` is then
    optimised so the *arrival* V_inf hits the next leg's required departure V_inf
    (``vinf_target_in``), satisfying flyby continuity AFTER the maneuver. The
    delivered ΔV is the DSM impulse magnitude ``dv_dsm_kms``.

    When ``vinf_target_in`` is ``None`` the leg minimises the delivered ΔV over
    ``eta`` (the cheapest powered close), which in the limit recovers a near-zero
    impulse where the ballistic leg already closes.

    The delivered ΔV is golden-anchored from below by the analytic VILM leveraging
    floor (:func:`cyclerfinder.search.vilm.vilm_dv_min`): a single in-leg DSM that
    performs the inter-moon V_inf change cannot beat the theoretical-minimum VILM
    ΔV for that transfer (the floor is the published Campagnola-Russell Endgame
    Part-1 Table 1/2 ΔV_min, already golden-validated in ``vilm.py``).
    """

    def __init__(self, *, max_revs: int = 0, n_eta: int = 25) -> None:
        self._max_revs = max_revs
        self._n_eta = max(3, n_eta)

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
        vinf_depart_mag: float | None = None,
    ) -> RelegResult:
        r_a = np.asarray(r_a, dtype=np.float64)
        v_a = np.asarray(v_a, dtype=np.float64)
        r_b = np.asarray(r_b, dtype=np.float64)
        v_b = np.asarray(v_b, dtype=np.float64)

        # The ballistic Lambert leg fixes the *seed* departure V_inf direction.
        # A pure-eta DSM with that exact departure reconstructs the same arc
        # (DSM impulse == 0): the leg geometry is pinned by (r_a, r_b, tof). The
        # retarget lever is the DEPARTURE V_inf — its magnitude factor + an
        # in-plane direction perturbation (the Takao vinf_out0/alpha0 genome) —
        # which moves the DSM point so the back-arc Lambert delivers a different
        # arrival V_inf, the DSM impulse paying for the mismatch.
        try:
            sols = lambert(r_a, r_b, tof_s, mu=mu, max_revs=max(0, n_rev))
        except (LambertError, ValueError):
            return _infeasible()
        wanted = [s for s in sols if s.n_revs == n_rev]
        if not wanted:
            return _infeasible()
        seed = min(wanted, key=lambda s: float(np.linalg.norm(s.v1 - v_a)))
        vinf_seed = np.asarray(seed.v1, dtype=np.float64) - v_a
        vinf_seed_mag = float(np.linalg.norm(vinf_seed))
        if vinf_seed_mag < 1.0e-9:
            return _infeasible()
        u_hat = vinf_seed / vinf_seed_mag
        # Orthonormal in-plane perturbation axis (the moon orbits are coplanar in
        # the +z plane, so the +z normal gives a well-defined in-plane rotation).
        z_hat = np.array([0.0, 0.0, 1.0])
        e1 = np.cross(u_hat, z_hat)
        n_e1 = float(np.linalg.norm(e1))
        e1 = e1 / n_e1 if n_e1 > 1.0e-12 else np.array([1.0, 0.0, 0.0])

        def _depart(mag_factor: float, ang: float) -> Vec3:
            v_dir = math.cos(ang) * u_hat + math.sin(ang) * e1
            return v_a + (vinf_seed_mag * mag_factor) * v_dir

        pinned_mag_factor = vinf_depart_mag / vinf_seed_mag if vinf_depart_mag is not None else None

        def _eval(params: Sequence[float]) -> tuple[float, float, float] | None:
            """``(mag_factor, ang, eta) -> (vinf_out, vinf_in, dv_dsm)`` or None."""
            mag_factor, ang, eta = float(params[0]), float(params[1]), float(params[2])
            if pinned_mag_factor is not None:
                # Departure V_inf magnitude is pinned (flyby continuity by
                # construction); the optimiser only moves the direction + eta.
                mag_factor = pinned_mag_factor
            eta = min(max(eta, _ETA_LO), _ETA_HI)
            v_depart = _depart(mag_factor, ang)
            try:
                leg = dsm_leg(r_a, v_depart, tof_s, float(eta), r_b, mu=mu, max_revs=self._max_revs)
            except (LambertError, KeplerError, ValueError):
                # A too-energetic departure can drive the ballistic front arc
                # hyperbolic past Newton convergence, or the back-arc Lambert
                # degenerate; treat that sample as infeasible (the search skips
                # it) rather than crashing the leg solve.
                return None
            vinf_out = float(np.linalg.norm(v_depart - v_a))
            vinf_in = float(np.linalg.norm(leg.v_arrive - v_b))
            return vinf_out, vinf_in, float(leg.dv_dsm_kms)

        def _cost(params: Sequence[float]) -> float:
            ev = _eval(params)
            if ev is None:
                return 1.0e9
            _, vinf_in, dv = ev
            if vinf_target_in is None:
                # Cheapest powered close: minimise the delivered ΔV (recovers the
                # ballistic dv≈0 leg where it already closes).
                return dv
            # Retarget the arrival V_inf to the requested value (flyby continuity
            # after the maneuver); the ΔV is whatever the DSM must spend to do so.
            return abs(vinf_in - vinf_target_in)

        # Coarse grid over (magnitude factor, direction angle, eta) — the cost is
        # multi-modal in all three — then a Nelder-Mead refine from the best node.
        # When the caller PINS the departure V_inf magnitude (``vinf_depart_mag``,
        # used by the moon-tour driver to make every interior flyby V_inf-continuous
        # by construction), the magnitude factor is fixed so departure V_inf equals
        # the pinned value exactly; only the direction + eta are free.
        if pinned_mag_factor is not None:
            # _eval overrides the magnitude factor to the pinned value; the grid
            # node value is irrelevant, so a single placeholder suffices.
            mag_grid = np.array([pinned_mag_factor])
        else:
            mag_grid = np.linspace(0.6, 1.4, 5)
        ang_grid = np.linspace(-0.5, 0.5, 7)
        eta_grid = np.linspace(_ETA_LO, _ETA_HI, max(3, self._n_eta // 4))
        best_params: tuple[float, float, float] | None = None
        best_cost = 1.0e9
        for mf in mag_grid:
            for ag in ang_grid:
                for et in eta_grid:
                    c = _cost((float(mf), float(ag), float(et)))
                    if c < best_cost:
                        best_cost = c
                        best_params = (float(mf), float(ag), float(et))
        if best_params is None or best_cost >= 1.0e9:
            return _infeasible()

        def _cost_arr(x: NDArray[np.float64]) -> float:
            return _cost([float(x[0]), float(x[1]), float(x[2])])

        refined = minimize(
            _cost_arr,
            np.array(best_params, dtype=np.float64),
            method="Nelder-Mead",
            options={"xatol": 1e-4, "fatol": 1e-6, "maxiter": 400},
        )
        refined_x = [float(v) for v in refined.x]
        final_params = refined_x if _cost(refined_x) <= best_cost else list(best_params)
        ev = _eval(final_params)
        if ev is None:
            return _infeasible()
        vinf_out, vinf_in, dv_dsm = ev
        return RelegResult(vinf_out=vinf_out, vinf_in=vinf_in, dv_kms=dv_dsm, feasible=True)
