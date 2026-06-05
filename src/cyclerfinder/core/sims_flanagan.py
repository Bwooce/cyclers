r"""Sims-Flanagan low-thrust leg model (Phase 1 of the v2 low-thrust scope).

A low-thrust leg connects a start state to an end state over a fixed time of
flight. Following the Sims-Flanagan transcription, the continuous thrust arc is
discretised into ``N`` segments; the thrust over each segment is modelled as a
single impulsive ``Delta V`` applied at the segment midpoint, with ballistic
(two-body Kepler) coast arcs in between. The leg is propagated **forward** from
its start state and **backward** from its end state to a common **match point**;
the position/velocity/mass mismatch there (the *defect*) must vanish for the leg
to be dynamically consistent.

This module implements the leg model only — pure, typed, no optimiser coupling.
The per-segment ``Delta V`` vectors are the manoeuvre decision variables; an
all-zero schedule degenerates to a pure ballistic Kepler arc (the key
zero-thrust regression anchor in the tests).

References
----------
Yam, C. H., Di Lorenzo, D., and Izzo, D., "Constrained global optimization of
low-thrust interplanetary trajectories," *IEEE Congress on Evolutionary
Computation*, 2010, DOI 10.1109/cec.2010.5586019. The Sims-Flanagan framing,
the per-segment thrust bound ``Delta V_max = (T_max/m)(t_f - t_0)/N`` (Eq. 1),
the rocket-equation mass update ``m_{i+1} = m_i * exp(-Delta V_i/(g0*Isp))``
(Eq. 5), and the ``neq = 7`` (3-D position + velocity + mass) match-point
constraint count are taken from the in-repo extraction
``docs/v2-future-references.md`` §1 (the PDF is not mirrored locally).

Sims, J. A., and Flanagan, S. N., "Preliminary design of low-thrust
interplanetary missions," AAS/AIAA Astrodynamics Specialist Conference, 1999
(the originating transcription cited by Yam 2010).

Conventions
-----------
* Units: km, km/s, s, kg throughout. ``Isp`` in seconds; ``T_max`` in kN so
  that ``T_max / m`` has units km/s^2 (see :func:`SimsFlanaganLeg.__post_init__`
  notes). Acceleration over a segment of duration ``dt_seg`` seconds yields a
  ``Delta V`` capability in km/s.
* ``mu`` defaults to the heliocentric :data:`MU_SUN_KM3_S2`.
* The match point defaults to the temporal midpoint of the leg.

Plan: ``docs/superpowers/plans/2026-06-05-sims-flanagan-lowthrust.md`` (Phase 1).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from math import exp

import numpy as np
from numpy.typing import NDArray

from cyclerfinder.core.constants import MU_SUN_KM3_S2, STANDARD_GRAVITY_KM_S2
from cyclerfinder.core.kepler import propagate

Vec3 = NDArray[np.float64]  # shape (3,), dtype float64 — matches kepler.py/flyby.py
VecN3 = NDArray[np.float64]  # shape (N, 3), per-segment Delta V schedule


class SimsFlanaganError(Exception):
    """Base class for errors raised by the Sims-Flanagan leg model."""


@dataclass(frozen=True)
class SimsFlanaganLeg:
    r"""Configuration of a single low-thrust leg (Sims-Flanagan transcription).

    The leg carries everything fixed about the transcription — the boundary
    states, the segmentation, the spacecraft mass and propulsion parameters, and
    the match-point location. The per-segment ``Delta V`` schedule is *not* part
    of the config; it is passed separately to the propagation functions as the
    manoeuvre decision variable.

    Attributes
    ----------
    r0, v0:
        Inertial start position/velocity, ``(3,)`` float64, km and km/s.
    rf, vf:
        Inertial end position/velocity, ``(3,)`` float64, km and km/s.
    tof_s:
        Total time of flight across the leg, seconds. Must be positive.
    n_segments:
        Number of Sims-Flanagan segments ``N``. Must be positive. Each segment
        spans ``tof_s / n_segments`` seconds and carries one impulsive
        ``Delta V`` at its midpoint.
    m0_kg:
        Spacecraft mass at the leg start, kg. Must be positive.
    isp_s:
        Specific impulse ``Isp``, seconds. Must be positive.
    tmax_kn:
        Maximum thrust ``T_max``, kN (kg*km/s^2). With mass in kg this gives an
        acceleration ``T_max / m`` in km/s^2 consistent with the km/km/s/s unit
        system used everywhere else. Must be non-negative; ``0.0`` models a pure
        ballistic (coast-only) leg whose only feasible schedule is all-zero.
    mu:
        Central-body gravitational parameter, km^3/s^2. Defaults to the
        heliocentric :data:`MU_SUN_KM3_S2`.
    match_index:
        Segment boundary index ``k`` at which the forward and backward
        propagations meet, ``0 <= k <= n_segments``. The forward propagation
        covers segments ``[0, k)`` and the backward propagation covers segments
        ``[k, n_segments)``. Defaults to the temporal midpoint
        (``n_segments // 2``), per the usual Sims-Flanagan choice.

    Notes
    -----
    The per-segment time step is ``dt_seg = tof_s / n_segments``. Within a
    segment the trajectory coasts a half step, receives the segment's impulsive
    ``Delta V``, then coasts the remaining half step (the midpoint-impulse
    convention). The per-segment thrust capability (Yam Eq. 1, applied with the
    instantaneous segment mass rather than a constant mass) is
    ``Delta V_max,i = (T_max / m_i) * dt_seg``.
    """

    r0: Vec3
    v0: Vec3
    rf: Vec3
    vf: Vec3
    tof_s: float
    n_segments: int
    m0_kg: float
    isp_s: float
    tmax_kn: float
    mu: float = MU_SUN_KM3_S2
    match_index: int = field(default=-1)

    def __post_init__(self) -> None:
        if self.tof_s <= 0.0:
            raise SimsFlanaganError(f"tof_s must be positive, got {self.tof_s}")
        if self.n_segments <= 0:
            raise SimsFlanaganError(f"n_segments must be positive, got {self.n_segments}")
        if self.m0_kg <= 0.0:
            raise SimsFlanaganError(f"m0_kg must be positive, got {self.m0_kg}")
        if self.isp_s <= 0.0:
            raise SimsFlanaganError(f"isp_s must be positive, got {self.isp_s}")
        if self.tmax_kn < 0.0:
            raise SimsFlanaganError(f"tmax_kn must be non-negative, got {self.tmax_kn}")
        # Resolve the default match point (temporal midpoint) lazily; the frozen
        # dataclass forces object.__setattr__ for the in-place fixup.
        resolved = self.n_segments // 2 if self.match_index < 0 else self.match_index
        if not (0 <= resolved <= self.n_segments):
            raise SimsFlanaganError(
                f"match_index must be in [0, {self.n_segments}], got {resolved}"
            )
        object.__setattr__(self, "match_index", resolved)

    @property
    def dt_seg_s(self) -> float:
        """Per-segment time step, seconds (``tof_s / n_segments``)."""
        return self.tof_s / self.n_segments


def _validate_schedule(leg: SimsFlanaganLeg, dvs: VecN3) -> Vec3:
    """Coerce and shape-check the per-segment ``Delta V`` schedule.

    Returns the schedule as a contiguous ``(N, 3)`` float64 array.
    """
    arr = np.asarray(dvs, dtype=np.float64)
    if arr.shape != (leg.n_segments, 3):
        raise SimsFlanaganError(f"dvs must have shape ({leg.n_segments}, 3), got {arr.shape}")
    return arr


def segment_dv_bounds(leg: SimsFlanaganLeg, dvs: VecN3) -> NDArray[np.float64]:
    r"""Per-segment ``Delta V`` capability bound along a forward mass profile.

    Returns the length-``N`` array of ``Delta V_max,i = (T_max / m_i) *
    dt_seg`` where ``m_i`` is the spacecraft mass *entering* segment ``i`` under
    the supplied schedule (mass decreases segment to segment via the rocket
    equation). This is the inequality-constraint bound the optimiser will apply
    in Phase 3, and the invariant the tests assert against
    (``|Delta V_i| <= Delta V_max,i``).

    Yam Eq. 1 states the bound with a constant mass; applying it with the
    *instantaneous* segment mass is a tighter, physically honest variant (a
    spacecraft that has already burned propellant accelerates faster, so its
    per-segment capability rises — the bound grows monotonically along the leg).

    Parameters
    ----------
    leg:
        Leg configuration.
    dvs:
        Per-segment ``Delta V`` schedule, ``(N, 3)`` km/s. Determines the mass
        profile via the rocket equation.

    Returns
    -------
    NDArray[np.float64]
        Length-``N`` array of per-segment ``Delta V`` capability bounds, km/s.
    """
    arr = _validate_schedule(leg, dvs)
    dt = leg.dt_seg_s
    g0_isp = STANDARD_GRAVITY_KM_S2 * leg.isp_s
    bounds = np.empty(leg.n_segments, dtype=np.float64)
    mass = leg.m0_kg
    for i in range(leg.n_segments):
        bounds[i] = (leg.tmax_kn / mass) * dt
        dv_mag = float(np.linalg.norm(arr[i]))
        mass *= exp(-dv_mag / g0_isp)
    return bounds


@dataclass(frozen=True)
class MatchPointState:
    """State at the match point reached by one half of a leg propagation.

    Attributes
    ----------
    r, v:
        Inertial position/velocity at the match point, ``(3,)`` float64, km and
        km/s.
    mass_kg:
        Spacecraft mass at the match point, kg. For the backward propagation
        this is the mass back-propagated from the (caller-supplied) end mass,
        i.e. the mass the spacecraft *must* have had at the match point to arrive
        at the end mass after the remaining burns.
    """

    r: Vec3
    v: Vec3
    mass_kg: float


def propagate_forward(leg: SimsFlanaganLeg, dvs: VecN3) -> MatchPointState:
    r"""Propagate the leg forward from ``(r0, v0, m0)`` to the match point.

    For each segment ``i`` in ``[0, match_index)``: coast a half segment
    (``dt_seg / 2``) under two-body dynamics, apply the impulsive ``Delta V_i``
    (subtracting mass via the rocket equation), then coast the remaining half
    segment. Returns the state and mass at the match point.

    Parameters
    ----------
    leg:
        Leg configuration.
    dvs:
        Per-segment ``Delta V`` schedule, ``(N, 3)`` km/s. Only the first
        ``match_index`` rows are consumed here.

    Returns
    -------
    MatchPointState
        Forward state ``(r, v, mass)`` at the match point.
    """
    arr = _validate_schedule(leg, dvs)
    half = 0.5 * leg.dt_seg_s
    g0_isp = STANDARD_GRAVITY_KM_S2 * leg.isp_s

    r = np.asarray(leg.r0, dtype=np.float64).copy()
    v = np.asarray(leg.v0, dtype=np.float64).copy()
    mass = leg.m0_kg

    for i in range(leg.match_index):
        r, v = propagate(r, v, half, leg.mu)
        dv = arr[i]
        v = v + dv
        dv_mag = float(np.linalg.norm(dv))
        mass *= exp(-dv_mag / g0_isp)
        r, v = propagate(r, v, half, leg.mu)

    return MatchPointState(r=r, v=v, mass_kg=mass)


def propagate_backward(leg: SimsFlanaganLeg, dvs: VecN3, mf_kg: float) -> MatchPointState:
    r"""Propagate the leg backward from ``(rf, vf, mf)`` to the match point.

    The time-reverse of :func:`propagate_forward`. For each segment ``i`` from
    ``n_segments - 1`` down to ``match_index``: coast backward a half segment
    (``-dt_seg / 2``), undo the impulsive ``Delta V_i`` (subtracting it from the
    velocity, since forward applied ``+Delta V_i``), then coast backward the
    remaining half segment. Mass is propagated *backward*: a forward burn
    *reduces* mass, so going backward the mass *increases* by the inverse
    rocket-equation factor ``m_before = m_after * exp(+|Delta V_i| / (g0*Isp))``.

    Parameters
    ----------
    leg:
        Leg configuration.
    dvs:
        Per-segment ``Delta V`` schedule, ``(N, 3)`` km/s. Only rows
        ``[match_index, n_segments)`` are consumed here.
    mf_kg:
        Spacecraft mass at the leg end, kg. Must be positive. For a feasible leg
        this equals the forward-propagated mass after all ``N`` burns; supplying
        it lets the backward pass produce a match-point mass directly comparable
        to the forward one.

    Returns
    -------
    MatchPointState
        Backward state ``(r, v, mass)`` at the match point.
    """
    if mf_kg <= 0.0:
        raise SimsFlanaganError(f"mf_kg must be positive, got {mf_kg}")
    arr = _validate_schedule(leg, dvs)
    half = 0.5 * leg.dt_seg_s
    g0_isp = STANDARD_GRAVITY_KM_S2 * leg.isp_s

    r = np.asarray(leg.rf, dtype=np.float64).copy()
    v = np.asarray(leg.vf, dtype=np.float64).copy()
    mass = mf_kg

    for i in range(leg.n_segments - 1, leg.match_index - 1, -1):
        r, v = propagate(r, v, -half, leg.mu)
        dv = arr[i]
        v = v - dv
        dv_mag = float(np.linalg.norm(dv))
        mass *= exp(dv_mag / g0_isp)
        r, v = propagate(r, v, -half, leg.mu)

    return MatchPointState(r=r, v=v, mass_kg=mass)


def final_mass(leg: SimsFlanaganLeg, dvs: VecN3) -> float:
    r"""Spacecraft mass at the leg end after the full ``Delta V`` schedule.

    Applies the rocket equation across all ``N`` segments forward from
    ``m0_kg``. This is the value to pass as ``mf_kg`` to
    :func:`propagate_backward` for a self-consistent (zero mass-defect) leg.

    Parameters
    ----------
    leg:
        Leg configuration.
    dvs:
        Per-segment ``Delta V`` schedule, ``(N, 3)`` km/s.

    Returns
    -------
    float
        Final spacecraft mass, kg.
    """
    arr = _validate_schedule(leg, dvs)
    g0_isp = STANDARD_GRAVITY_KM_S2 * leg.isp_s
    mass = leg.m0_kg
    for i in range(leg.n_segments):
        mass *= exp(-float(np.linalg.norm(arr[i])) / g0_isp)
    return mass


def match_point_defect(leg: SimsFlanaganLeg, dvs: VecN3, mf_kg: float) -> NDArray[np.float64]:
    r"""Match-point defect 7-vector ``[Delta r (3), Delta v (3), Delta m (1)]``.

    The forward propagation (from the start) and the backward propagation (from
    the end) must agree at the match point for the leg to be dynamically
    consistent. This returns ``S_mf - S_mb`` as a length-7 vector — the
    ``neq = 7`` constraint count of the 3-D + mass problem recorded from Yam §1.
    Components: position mismatch (km, 3), velocity mismatch (km/s, 3), and mass
    mismatch (kg, 1). At convergence this vector is ~0.

    Parameters
    ----------
    leg:
        Leg configuration.
    dvs:
        Per-segment ``Delta V`` schedule, ``(N, 3)`` km/s.
    mf_kg:
        Spacecraft mass at the leg end, kg (see :func:`propagate_backward`).

    Returns
    -------
    NDArray[np.float64]
        Length-7 defect vector ``[dr_x, dr_y, dr_z, dv_x, dv_y, dv_z, dm]``.
    """
    fwd = propagate_forward(leg, dvs)
    bwd = propagate_backward(leg, dvs, mf_kg)
    defect = np.empty(7, dtype=np.float64)
    defect[0:3] = fwd.r - bwd.r
    defect[3:6] = fwd.v - bwd.v
    defect[6] = fwd.mass_kg - bwd.mass_kg
    return defect
