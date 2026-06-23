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

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from cyclerfinder.core.constants import MU_SUN_KM3_S2, PLANETS
from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.core.flyby import flyby_dv
from cyclerfinder.core.lambert import LambertError, lambert
from cyclerfinder.nbody.forces import RailsEphemerisCache
from cyclerfinder.nbody.propagator import RestrictedNBody

Vec3 = NDArray[np.float64]

# Per-body minimum unpowered-flyby altitude now lives in the single source of truth
# PLANETS[body].safe_alt_km (sourced per body: Earth/Mars 200 km, Russell 2004 p.165;
# giants radiation/ring-limited). The chain uses that per-body floor by default
# (min_flyby_alt_km=None); pass a scalar override only to reproduce a source designed
# under a different floor. NOTE: a too-high floor spuriously charges marginal flybys —
# the old unsourced 300 km charged ~40 m/s to S1L1's one near-limit Earth flyby.

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
    ``r_target`` (only meaningful when ``converged``). ``v_arrive_km_s`` is the
    heliocentric velocity at arrival (``t1``) on that converged arc — the chain
    driver subtracts the arrival-planet velocity from it to get the incoming
    ``V_inf`` at the next flyby node. ``miss_km`` is the final arrival miss;
    ``n_iter`` the Newton steps taken; ``converged`` the honesty flag (False =>
    propagation diverged or Newton stalled; do not use the velocities).
    """

    v_dep_km_s: Vec3
    v_arrive_km_s: Vec3
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
    zero3 = np.zeros(3, dtype=np.float64)

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
                v_dep_km_s=v,
                v_arrive_km_s=zero3,
                miss_km=float("inf"),
                converged=False,
                n_iter=n_iter,
            )

        residual = np.asarray(arc.r_km, dtype=np.float64) - r_target
        miss_km = float(np.linalg.norm(residual))
        if miss_km < tol_km:
            return TargetLegResult(
                v_dep_km_s=v,
                v_arrive_km_s=np.asarray(arc.v_km_s, dtype=np.float64),
                miss_km=miss_km,
                converged=True,
                n_iter=n_iter,
            )

        jac = np.asarray(arc.stm, dtype=np.float64)[0:3, 3:6]  # d r_f / d v_0
        try:
            step = np.linalg.solve(jac, -residual)
        except np.linalg.LinAlgError:
            # Singular sensitivity — cannot target this leg.
            return TargetLegResult(
                v_dep_km_s=v,
                v_arrive_km_s=zero3,
                miss_km=miss_km,
                converged=False,
                n_iter=n_iter,
            )
        v = v + step

    # Exhausted iterations without reaching tol_km.
    return TargetLegResult(
        v_dep_km_s=v, v_arrive_km_s=zero3, miss_km=miss_km, converged=False, n_iter=max_iter
    )


# ---------------------------------------------------------------------------
# The continuous maintenance chain (M7 Task 2)
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class MaintenanceNode:
    """One encounter of the continuous maintenance chain.

    Interior (``is_flyby``) nodes carry a maintenance ``dv_kms`` — the part of the
    arrival-to-departure ``V_inf`` transition a free ballistic flyby cannot supply
    (|Δ|V_inf|| magnitude term + un-bendable excess-bend impulse, via
    :func:`cyclerfinder.core.flyby.flyby_dv_for`). The first node (initial departure)
    and last node (final arrival) have no flyby and ``dv_kms == 0.0``. ``converged``
    is False if either adjoining leg's targeting solve failed — then ``dv_kms`` is
    ``inf`` (unmeasurable) and the whole chain is flagged ``diverged``.
    """

    index: int
    body: str
    is_flyby: bool
    miss_km: float
    vinf_in_kms: float
    vinf_out_kms: float
    dv_kms: float
    converged: bool


@dataclass(frozen=True)
class MaintenanceChain:
    """Result of :func:`continuous_maintenance_chain` over a node sequence.

    ``horizon_tcm_mps`` is the summed maintenance ΔV across all flyby nodes (m/s),
    or ``inf`` if any leg's targeting diverged (the row is then unmeasurable, left
    V0 — never forced). ``per_cycle_tcm_mps`` divides by ``n_cycles`` when given.
    """

    nodes: tuple[MaintenanceNode, ...]
    horizon_tcm_mps: float
    per_cycle_tcm_mps: float
    diverged: bool
    n_legs_converged: int
    n_legs_total: int


def continuous_maintenance_chain(
    node_epochs_sec: Sequence[float],
    node_bodies: Sequence[str],
    ephem: Ephemeris,
    prop: RestrictedNBody,
    *,
    cruise_perturbers: tuple[str, ...] = (),
    exclude_endpoint_bodies: bool = True,
    leg_v_guess: Sequence[Vec3] | None = None,
    min_flyby_alt_km: float | None = None,
    n_cycles: int = 1,
    tol_km: float = DEFAULT_TOL_KM,
    max_iter: int = DEFAULT_MAX_ITER,
    accuracy: float = 1e-10,
    max_wall_sec: float = 90.0,
) -> MaintenanceChain:
    """Walk a cycler's encounter sequence as ONE continuous n-body trajectory.

    The patched-conic-node model with REAL flybys: each node sits at the real planet
    position (position continuity is exact); each leg is propagated in restricted
    n-body and **targeted** to hit the next node's planet via :func:`target_leg`
    (Lambert-seeded). At each interior node the arrival ``V_inf`` (incoming leg) and the
    departure ``V_inf`` (targeting the outgoing leg) differ; the maintenance ΔV is what
    a free ballistic flyby cannot bridge (:func:`flyby_dv_for`). Summing those over the
    node span is the horizon TCM.

    Generalises the Sun-only proxy
    (:func:`cyclerfinder.search.appc_corrected.continuous_chain`) by (a) propagating
    legs with the configured perturbers instead of Kepler and (b) targeting the next
    planet position (n-body Lambert) instead of using the sourced ``V_inf`` direction —
    which is what makes the perturbed chain self-consistent rather than divergent.

    **Perturber rule (the general, artifact-free model).** ``cruise_perturbers`` lists
    the *system's* significant gravitating bodies (beyond the heliocentric central
    body) that perturb the cruise — e.g. ``()`` for Earth-Mars (no third body
    materially perturbs the 1-2 AU heliocentric cruise), or other planets for a
    multi-flyby heliocentric cycler (a Venus-Earth-Mars chain, Jupiter, ...). With
    ``exclude_endpoint_bodies=True`` (default) each leg's perturber set has the leg's
    OWN two endpoint flyby bodies removed: those flybys are modelled patched-conic AT
    the nodes, so integrating an endpoint body continuously from its own node centre
    would start the spacecraft inside that body's softened core — the divergent
    patched-conic handoff artifact, not a fuel cost. Excluding only the endpoints keeps
    every genuine third-body perturbation (a body that is never this leg's endpoint)
    while avoiding the artifact. For Earth-Mars this reduces to Sun-cruise; for a moon
    tour or a Venus-flyby chain it captures the real cross-cruise perturbation.

    Note: this is the heliocentric-central-body lane (the Lambert seed uses ``MU_SUN``
    and the propagator is Sun-central). Planet-central moon-tour M7 (Saturn/Uranus as
    the central body) is a separate generalisation; those rows currently validate via
    their own V3/V4 lanes, not :func:`real_closure.verify_real_closure`.

    Parameters
    ----------
    node_epochs_sec, node_bodies:
        Parallel sequences of encounter epochs (J2000 seconds) and one-letter body
        codes (e.g. ``("E", "M", "E", ...)``). Length ``n`` => ``n-1`` legs =>
        ``n-2`` interior flybys.
    ephem, prop:
        Ephemeris backend and :class:`RestrictedNBody` propagator.
    cruise_perturbers:
        System significant perturber body codes (see the perturber rule above).
        ``()`` => Sun-only cruise.
    exclude_endpoint_bodies:
        When True (default, the correct model) each leg drops its own endpoint flyby
        bodies from the perturber set. Set False ONLY to reproduce the naive-handoff
        divergence artifact in tests.
    leg_v_guess:
        Optional per-leg departure-velocity seeds (length ``n-1``, km/s). REQUIRED in
        practice for multi-rev resonant cyclers: pass the rev-correct sourced /
        constructed departure velocities so the Newton stays in the intended transfer
        basin. Without it the chain seeds from a single-rev Lambert, which lands a
        high-energy single-rev solution on multi-rev legs (wildly inflated TCM).
    min_flyby_alt_km:
        Scalar override (km) for the per-flyby bend feasibility / ΔV floor. Default
        ``None`` uses each body's sourced ``PLANETS[body].safe_alt_km`` (Earth/Mars
        200 km, Russell 2004 p.165). Pass a scalar only to reproduce a source designed
        under a different floor — a too-high floor spuriously charges marginal flybys.
    n_cycles:
        Number of cycles the node span covers (for ``per_cycle_tcm_mps``).
    tol_km, max_iter, accuracy, max_wall_sec:
        Forwarded to :func:`target_leg`.

    Returns
    -------
    MaintenanceChain
    """
    n = len(node_epochs_sec)
    if n != len(node_bodies):
        raise ValueError("node_epochs_sec and node_bodies must have equal length")
    if n < 2:
        raise ValueError(f"need >= 2 nodes, got {n}")
    if leg_v_guess is not None and len(leg_v_guess) != n - 1:
        raise ValueError(f"leg_v_guess must have length n-1={n - 1}, got {len(leg_v_guess)}")

    states = [
        tuple(
            np.asarray(x, dtype=np.float64) for x in ephem.state(node_bodies[i], node_epochs_sec[i])
        )
        for i in range(n)
    ]

    # Target each leg i -> i+1: record the converged departure / arrival velocities.
    leg_vdep: list[Vec3] = []
    leg_varr: list[Vec3] = []
    leg_miss: list[float] = []
    leg_conv: list[bool] = []
    for i in range(n - 1):
        r_i, v_i = states[i]
        r_j, _v_j = states[i + 1]
        tof = node_epochs_sec[i + 1] - node_epochs_sec[i]
        if leg_v_guess is not None:
            # Caller-supplied departure-velocity seed (the rev-correct sourced /
            # constructed value). ESSENTIAL for multi-rev resonant legs: a single-rev
            # Lambert guess lands the Newton on a high-energy single-rev transfer (wildly
            # wrong V_inf), so the sourced seed keeps it in the intended basin.
            v_guess = np.asarray(leg_v_guess[i], dtype=np.float64)
        else:
            try:
                v_guess = np.asarray(
                    lambert(r_i, r_j, tof, mu=MU_SUN_KM3_S2)[0].v1, dtype=np.float64
                )
            except (LambertError, ValueError):
                v_guess = v_i  # degenerate geometry — Newton works from the planet velocity
        # Per-leg perturbers: the system set minus this leg's endpoint flyby bodies
        # (those are patched-conic at the nodes — see the perturber rule in the docstring).
        if exclude_endpoint_bodies:
            endpoints = {node_bodies[i], node_bodies[i + 1]}
            leg_bodies = tuple(b for b in cruise_perturbers if b not in endpoints)
        else:
            leg_bodies = tuple(cruise_perturbers)
        res = target_leg(
            prop,
            r_i,
            node_epochs_sec[i],
            node_epochs_sec[i + 1],
            r_j,
            v_guess,
            bodies=leg_bodies,
            ephem=ephem,
            cache=None,
            tol_km=tol_km,
            max_iter=max_iter,
            accuracy=accuracy,
            max_wall_sec=max_wall_sec,
        )
        leg_vdep.append(res.v_dep_km_s)
        leg_varr.append(res.v_arrive_km_s)
        leg_miss.append(res.miss_km)
        leg_conv.append(res.converged)

    nodes: list[MaintenanceNode] = []
    horizon_kms = 0.0
    for i in range(n):
        _r_i, v_pl = states[i]
        is_flyby = 0 < i < n - 1
        vin = (leg_varr[i - 1] - v_pl) if i > 0 else np.zeros(3)
        vout = (leg_vdep[i] - v_pl) if i < n - 1 else np.zeros(3)
        miss_km = leg_miss[i - 1] if i > 0 else 0.0
        dv = 0.0
        conv = True
        if is_flyby:
            conv = leg_conv[i - 1] and leg_conv[i]
            if conv:
                pdata = PLANETS[node_bodies[i]]
                # Per-body sourced floor (single source of truth) unless overridden.
                alt = pdata.safe_alt_km if min_flyby_alt_km is None else min_flyby_alt_km
                rp_min = pdata.radius_eq_km + alt
                dv = flyby_dv(vin, vout, pdata.mu_km3_s2, rp_min)
                horizon_kms += dv
            else:
                dv = float("inf")
        nodes.append(
            MaintenanceNode(
                index=i,
                body=node_bodies[i],
                is_flyby=is_flyby,
                miss_km=miss_km,
                vinf_in_kms=float(np.linalg.norm(vin)) if i > 0 else 0.0,
                vinf_out_kms=float(np.linalg.norm(vout)) if i < n - 1 else 0.0,
                dv_kms=dv,
                converged=conv,
            )
        )

    diverged = not all(leg_conv)
    horizon_mps = float("inf") if diverged else horizon_kms * 1000.0
    per_cycle = horizon_mps / n_cycles if (n_cycles > 0 and not diverged) else float("inf")
    return MaintenanceChain(
        nodes=tuple(nodes),
        horizon_tcm_mps=horizon_mps,
        per_cycle_tcm_mps=per_cycle,
        diverged=diverged,
        n_legs_converged=sum(leg_conv),
        n_legs_total=n - 1,
    )
