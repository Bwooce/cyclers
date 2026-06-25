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
* :class:`LowThrustReleg` — the secondary powered backend (#449 Task 7). Wraps the
  Sims-Flanagan N-segment low-thrust leg (``core.sims_flanagan.SimsFlanaganLeg`` +
  ``search.lowthrust.solve_leg_min_dv``, #309): the boundary states are pinned to
  the moon encounter states retargeted to the requested departure/arrival V_inf,
  and the deliverable ΔV is distributed across the thrust train (so a low-thrust
  leg is strictly more ΔV-efficient than the single-impulse DSM for the same V_inf
  change). Gated behind the DSM branch (design §7): SF is slower to converge and
  its leg-model golden is bracket-only — the SF delivered ΔV is asserted to BRACKET
  the DSM/VILM-floor result (a different transcription, so bracket not equal), with
  the same sourced VILM floor as the non-circular lower bound.

Reuse over rebuild: no optimiser is written here and no corrector is changed; the
DSM leg solver (#307), the SF low-thrust leg solver (#309) and the VILM cost model
are imported, not re-derived.
"""

from __future__ import annotations

import math
from collections.abc import Callable, Sequence
from dataclasses import dataclass
from typing import Protocol, runtime_checkable

import numpy as np
from numpy.typing import NDArray
from scipy.optimize import minimize

from cyclerfinder.core.constants import MU_SUN_KM3_S2
from cyclerfinder.core.kepler import KeplerError
from cyclerfinder.core.lambert import LambertError, LambertSolution, lambert
from cyclerfinder.core.sims_flanagan import SimsFlanaganError, SimsFlanaganLeg, segment_dv_bounds
from cyclerfinder.search.dsm_leg import dsm_leg
from cyclerfinder.search.leveraging_chain import walk_vinf_down
from cyclerfinder.search.leveraging_leg import LeveragingLegResult
from cyclerfinder.search.lowthrust import solve_leg_min_dv

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
    chain_hops: tuple[LeveragingLegResult, ...] = ()


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
        arrival_moon: str | None = None,
    ) -> RelegResult:
        """Solve one moon-to-moon leg; return its V_inf chain + delivered ΔV.

        ``vinf_target_in`` retargets the arrival V_inf (powered backends only);
        ``vinf_depart_mag`` pins the departure V_inf magnitude (powered backends
        only, used by the moon-tour driver for by-construction flyby continuity).
        The ballistic backend ignores both. ``arrival_moon`` names the arrival
        flyby body (the leverage body); only the multi-rev leveraging backend
        uses it — the ballistic/DSM/SF backends read V_inf off the Lambert seed
        and ignore it (a backwards-compatible protocol add, #465).
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
        arrival_moon: str | None = None,
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
        arrival_moon: str | None = None,
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


class LowThrustReleg:
    r"""Powered Sims-Flanagan low-thrust leg — distributed-ΔV V_inf retarget (#449).

    The secondary powered backend (Task 7), wrapping the #309 Sims-Flanagan leg
    solver as a third :class:`Releg` behind the same protocol as :class:`DsmReleg`
    (swappable into the moon-tour driver with no driver rewrite). Where the DSM
    backend delivers the V_inf retarget with a SINGLE in-leg impulse, the
    low-thrust backend distributes the deliverable ΔV across an ``n_segments``
    thrust train, which is strictly more ΔV-efficient for the same V_inf change
    (the physics the whole #449 low-thrust bet rests on).

    Boundary states (the moon-tour driver's continuity-by-construction contract):

    * Departure ``v0 = v_a_moon + vinf_depart_mag * û_depart`` — the departure
      V_inf magnitude is PINNED to ``vinf_depart_mag`` (the driver's common flyby
      target ``T``) along the ballistic-Lambert departure direction.
    * Arrival ``vf = v_b_moon + vinf_target_in * û_arrive`` — the arrival V_inf
      magnitude is RETARGETED to ``vinf_target_in`` (again ``T``) along the
      ballistic-Lambert arrival direction.

    The SF leg then solves the minimum-ΔV thrust schedule connecting those fixed
    boundary states over the leg ToF (``search.lowthrust.solve_leg_min_dv``); the
    delivered ΔV is ``total_dv_kms``. When neither V_inf is retargeted (both
    ``None``) the boundary states are the ballistic Lambert endpoints and the leg's
    all-zero schedule closes at ΔV ≈ 0 — the coplanar/zero-retarget regression
    limit that reproduces :class:`BallisticReleg` (the SF leg's zero-thrust limit
    IS the ballistic leg).

    Thrust capability is auto-scaled (``thrust_reach``) so the per-leg
    segment-ΔV budget comfortably exceeds the requested V_inf change: a leg with
    too little thrust capability simply cannot reach the retargeted boundary state
    (the defect never closes), so the backend sizes ``tmax_kn`` from the V_inf
    delta. This is a *capability* knob, not a spacecraft-specific number: the
    backend reports the ΔV the distributed train must spend to deliver the
    retarget, which is what the driver scores against the powered dv-band.

    Golden discipline (design §6/§8). The SF leg model has no clean state-level
    literature anchor (its own ``lowthrust.py`` docstring: "there is no usable
    literature anchor for the leg model"; the Vasile-Campagnola 2009 DFET
    transcription "DOES NOT MAP" to our SF leg, digest
    ``2026-06-07-vasile-campagnola-dfet-method-mining.md`` §2.6). So the SF golden
    is BRACKET-only: the delivered ΔV must be ≥ the same sourced VILM floor
    (``vilm.vilm_dv_min``) the DSM branch uses (a powered leg cannot beat the
    theoretical-minimum VILM ΔV for the transfer) and ≤ the single-impulse DSM
    cost (distributed thrust is more efficient than one impulse). The non-circular
    lower bound is the published Campagnola-Russell floor, never a value SF
    computed.
    """

    def __init__(
        self,
        *,
        n_segments: int = 10,
        m0_kg: float = 1000.0,
        isp_s: float = 3000.0,
        n_starts: int = 2,
        thrust_reach: float = 4.0,
        min_tmax_kn: float = 1.0e-3,
    ) -> None:
        self._n_segments = max(2, n_segments)
        self._m0_kg = m0_kg
        self._isp_s = isp_s
        self._n_starts = max(1, n_starts)
        # The per-leg segment-ΔV budget is sized to ``thrust_reach`` times the
        # requested V_inf delta so the leg has the reach to deliver the retarget
        # (a too-weak train cannot close the boundary defect). >= 1 is required for
        # any reach; 4x leaves comfortable margin without making the box wild.
        self._thrust_reach = max(1.5, thrust_reach)
        self._min_tmax_kn = min_tmax_kn

    def _tmax_kn_for_reach(self, dv_reach_kms: float, tof_s: float) -> float:
        r"""Choose ``tmax_kn`` so the leg's total segment-ΔV budget ≈ ``dv_reach``.

        ``segment_dv_bounds`` returns ``Δv_max,i = (T_max / m_i) * dt_seg``; summed
        over ``N`` segments the (constant-mass) total budget is
        ``N * (T_max / m0) * (tof / N) = (T_max / m0) * tof``. Inverting for the
        thrust that yields a total budget of ``dv_reach`` km/s gives
        ``T_max = m0 * dv_reach / tof`` (kN, since the unit system is km/km/s/s).
        """
        tmax = self._m0_kg * dv_reach_kms / tof_s
        return max(self._min_tmax_kn, tmax)

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
        arrival_moon: str | None = None,
    ) -> RelegResult:
        r_a = np.asarray(r_a, dtype=np.float64)
        v_a = np.asarray(v_a, dtype=np.float64)
        r_b = np.asarray(r_b, dtype=np.float64)
        v_b = np.asarray(v_b, dtype=np.float64)

        # The ballistic Lambert leg fixes the departure/arrival V_inf DIRECTIONS;
        # the requested magnitudes (pin/retarget) set the boundary states. With no
        # retarget the boundary states are exactly the ballistic Lambert endpoints
        # (the zero-thrust SF leg then closes at ΔV ~ 0 — the regression limit).
        try:
            sols = lambert(r_a, r_b, tof_s, mu=mu, max_revs=max(0, n_rev))
        except (LambertError, ValueError):
            return _infeasible()
        wanted = [s for s in sols if s.n_revs == n_rev]
        if not wanted:
            return _infeasible()
        seed = min(wanted, key=lambda s: float(np.linalg.norm(s.v1 - v_a)))

        vinf_dep_seed = np.asarray(seed.v1, dtype=np.float64) - v_a
        vinf_arr_seed = np.asarray(seed.v2, dtype=np.float64) - v_b
        dep_mag_seed = float(np.linalg.norm(vinf_dep_seed))
        arr_mag_seed = float(np.linalg.norm(vinf_arr_seed))
        if dep_mag_seed < 1.0e-9 or arr_mag_seed < 1.0e-9:
            return _infeasible()
        u_dep = vinf_dep_seed / dep_mag_seed
        u_arr = vinf_arr_seed / arr_mag_seed

        dep_mag = vinf_depart_mag if vinf_depart_mag is not None else dep_mag_seed
        arr_mag = vinf_target_in if vinf_target_in is not None else arr_mag_seed
        v0 = v_a + dep_mag * u_dep
        vf = v_b + arr_mag * u_arr

        # The boundary V_inf change the thrust train must bridge: the magnitude
        # change at departure + arrival vs the ballistic endpoints (the leg geometry
        # is otherwise the same arc). Size the thrust budget to cover it with reach.
        dv_reach = self._thrust_reach * (
            abs(dep_mag - dep_mag_seed) + abs(arr_mag - arr_mag_seed) + 0.5 * (dep_mag + arr_mag)
        )
        # Floor the reach so even the zero-retarget leg carries a small live thrust
        # budget (so the optimiser has room to drive the defect to zero, recovering
        # the ~0-ΔV ballistic close rather than being clamped to coast-only).
        dv_reach = max(dv_reach, 0.5)
        tmax_kn = self._tmax_kn_for_reach(dv_reach, tof_s)

        try:
            leg = SimsFlanaganLeg(
                r0=r_a,
                v0=v0,
                rf=r_b,
                vf=vf,
                tof_s=tof_s,
                n_segments=self._n_segments,
                m0_kg=self._m0_kg,
                isp_s=self._isp_s,
                tmax_kn=tmax_kn,
                mu=mu,
            )
        except SimsFlanaganError:
            return _infeasible()

        # Guard: if the sized budget still cannot span the boundary V_inf change the
        # leg can never close; report infeasible rather than a non-converged ΔV.
        cap_total = float(np.sum(segment_dv_bounds(leg, np.zeros((self._n_segments, 3)))))
        if cap_total < (abs(dep_mag - dep_mag_seed) + abs(arr_mag - arr_mag_seed)):
            return _infeasible()

        try:
            res = solve_leg_min_dv(leg, n_starts=self._n_starts, use_de=False)
        except (SimsFlanaganError, KeplerError, ValueError):
            return _infeasible()
        if not res.converged:
            return _infeasible()

        vinf_out = float(np.linalg.norm(v0 - v_a))
        vinf_in = float(np.linalg.norm(vf - v_b))
        return RelegResult(
            vinf_out=vinf_out, vinf_in=vinf_in, dv_kms=float(res.total_dv_kms), feasible=True
        )


class MultiRevLeveragingReleg:
    """Powered multi-rev leveraging leg — a CHAIN of resonant hops (#465).

    The third powered :class:`Releg`, behind the same protocol as
    :class:`DsmReleg` / :class:`LowThrustReleg` (swappable into the moon-tour
    driver with no driver rewrite). Where the DSM backend sheds the leg's whole
    V_inf defect in ONE impulse — paying the single-VILM *maximum* (``vilm`` Eq.14)
    — this backend internally CHAINS N resonant-hop legs
    (:func:`cyclerfinder.search.leveraging_chain.walk_vinf_down`, each one a #179
    apse VILM) to walk the arrival V_inf from its natural high value DOWN to the
    common flyby target step by step. The total delivered ΔV is the sum of the
    per-hop apse burns, which approaches the Eq.(13) multi-VILM *minimum* — roughly
    an order of magnitude cheaper than the single-impulse shed (design §1.2).

    Contract (preserving the driver's by-construction continuity):

    * Seed: the ballistic Lambert leg fixes the geometry + the NATURAL arrival
      V_inf magnitude ``V_inf_H`` (lowest-energy branch — the same seed
      :class:`DsmReleg` / :class:`LowThrustReleg` use).
    * ``vinf_target_in`` RETARGETS the arrival V_inf to the common flyby target
      ``T``: the chain walks ``V_inf_H`` down to ``T`` at the ARRIVAL moon.
    * ``vinf_depart_mag`` PINS the departure V_inf magnitude to ``T`` (continuity
      by construction — every leg departs/arrives at the same ``T``).
    * Zero-retarget limit: ``vinf_target_in is None`` ⇒ no hops ⇒ ``dv_kms = 0``,
      reproducing :class:`BallisticReleg` (the regression limit every backend
      honours, ``feedback_orbit_closure_discipline``).

    The arrival flyby body (the leverage body) is named by ``arrival_moon`` (the
    driver threads ``sequence[k+1]``); a constructor default ``moon`` is the
    fallback for direct callers. The chain ΔV is golden-anchored from below by the
    Eq.(13) leverage quadrature (a finite integer-resonance chain cannot beat the
    continuous minimum) and from above by the published finite-chain penalty
    (``leveraging_chain`` golden) — the same sourced floor the DSM branch uses.

    No new optimiser, no new cost model: the per-hop primitive
    (``leveraging_leg``, #179) and the floor (``vilm``, golden) both exist; this
    backend is the chain ORCHESTRATION (choose hop resonances, sum ΔV).
    """

    def __init__(
        self,
        *,
        moon: str | None = None,
        exterior: bool = True,
        max_hops: int = 500,
        max_revs: int = 5000,
    ) -> None:
        self._moon = moon
        self._exterior = exterior
        self._max_hops = max_hops
        self._max_revs = max_revs

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
        arrival_moon: str | None = None,
    ) -> RelegResult:
        r_a = np.asarray(r_a, dtype=np.float64)
        v_a = np.asarray(v_a, dtype=np.float64)
        r_b = np.asarray(r_b, dtype=np.float64)
        v_b = np.asarray(v_b, dtype=np.float64)

        # The ballistic Lambert leg fixes the geometry + the NATURAL arrival V_inf.
        try:
            sols = lambert(r_a, r_b, tof_s, mu=mu, max_revs=max(0, n_rev))
        except (LambertError, ValueError):
            return _infeasible()
        wanted = [s for s in sols if s.n_revs == n_rev]
        if not wanted:
            return _infeasible()
        seed = min(wanted, key=lambda s: float(np.linalg.norm(s.v1 - v_a)))
        vinf_out_seed = float(np.linalg.norm(seed.v1 - v_a))
        vinf_in_seed = float(np.linalg.norm(seed.v2 - v_b))

        # No retarget ⇒ ballistic-equivalent (no leveraging walk needed).
        if vinf_target_in is None:
            return RelegResult(
                vinf_out=vinf_out_seed, vinf_in=vinf_in_seed, dv_kms=0.0, feasible=True
            )

        moon = arrival_moon if arrival_moon is not None else self._moon
        if moon is None:
            # The leverage body must be named (the driver threads sequence[k+1]).
            return _infeasible()

        # CHAIN: walk the arrival V_inf from its natural high value down to the
        # common flyby target via N resonant hops at the arrival moon.
        chain = walk_vinf_down(
            moon,
            vinf_in_seed,
            vinf_target_in,
            exterior=self._exterior,
            max_hops=self._max_hops,
            max_revs=self._max_revs,
        )
        if not chain.converged:
            return _infeasible()

        # Departure V_inf is pinned to T by construction (the cycle is continuous
        # because every leg departs/arrives at the same T); arrival V_inf is the
        # retargeted T the chain reached.
        vinf_out = vinf_depart_mag if vinf_depart_mag is not None else vinf_out_seed
        return RelegResult(
            vinf_out=vinf_out,
            vinf_in=chain.vinf_end_kms,
            dv_kms=chain.total_dv_kms,
            feasible=True,
            chain_hops=chain.hops,
        )


# ---------------------------------------------------------------------------
# Swap-seam adapter: expose a Releg backend at the injected-`lambert` call site
# (discovery_campaign._close_one_phasing line 463 / v2_moontour._cycle_residual)
# WITHOUT changing the default. The ballistic adapter is the regression lock —
# it must reproduce the baseline lambert path bit-for-bit (design §3 / §7 risk
# "the two swap sites drift out of sync").
# ---------------------------------------------------------------------------


def ballistic_lambert_adapter(releg: BallisticReleg) -> Callable[..., list[LambertSolution]]:
    """A ``lambert``-signature callable backed by ``BallisticReleg`` (the seam).

    ``_close_one_phasing`` / ``_cycle_residual`` call the injected solver as
    ``lambert(r_a, r_b, tof_s, mu=mu, max_revs=...) -> list[LambertSolution]``
    and select the lowest-energy branch themselves. The ballistic releg is the
    ``dv_kms == 0`` zero-impulse limit that wraps exactly ``core.lambert.lambert``
    with the same lowest-energy branch pick, so the adapter returns the SAME
    ``LambertSolution`` list — proving the swap is a drop-in for the ballistic
    backend (bit-for-bit). The ``releg`` argument is the explicit backend the
    adapter is bound to (so a future powered adapter is a parallel constructor,
    not a hidden default).
    """
    if not isinstance(releg, BallisticReleg):  # pragma: no cover - guard
        raise TypeError(
            "ballistic_lambert_adapter requires a BallisticReleg (the zero-ΔV "
            "regression-locked backend); powered backends use the releg_moontour "
            "driver, not the lambert-compatible seam"
        )

    def _adapter(
        r1: Vec3,
        r2: Vec3,
        tof: float,
        *,
        mu: float = MU_SUN_KM3_S2,
        prograde: bool = True,
        max_revs: int = 0,
    ) -> list[LambertSolution]:
        return lambert(r1, r2, tof, mu=mu, prograde=prograde, max_revs=max_revs)

    return _adapter
