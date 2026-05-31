"""Cycler / Leg / Encounter dataclasses and per-cycler derived metrics.

Three ``frozen=True`` dataclasses represent a cycler as a tree:

* :class:`Encounter` — a flyby of a single body at a single instant, carrying
  the planet's heliocentric state and the spacecraft :math:`V_\\infty`
  in/out vectors.
* :class:`Leg` — the heliocentric arc between two consecutive encounters,
  carrying the spacecraft departure and arrival velocities and the Lambert
  branch selection (``n_revs``, ``branch``).
* :class:`Cycler` — the ordered sequence of encounters and legs plus the
  total period, with read-only metric methods.

The dataclasses are immutable so a derived metric (closure residual,
maintenance ΔV, etc.) cannot silently drift relative to the geometry it
was computed from. Construction is delegated to
:mod:`cyclerfinder.search.construct`; ``__post_init__`` validation is kept
deliberately thin (numeric shape checks, not physics) because the source
of truth for "this is a valid cycler" is the constructor, not the
dataclass.

Closure convention
------------------
``Cycler.closure_residual(omega)`` rotates the final-encounter departure
velocity by ``omega * period`` (the synodic frame's advance over one
period) and compares it to the first-encounter departure velocity in
km/s. The residual is therefore a **velocity** residual; position closure
is automatic by construction of the encounter times against the
ephemeris (the last encounter happens at body B at ``t = period``; one
period later, body A sits at the same rotating-frame angular position as
the first encounter). What can fail and what we measure is the velocity
match — i.e. the maintenance ΔV in the idealised limit.

The default ``omega`` is :func:`~cyclerfinder.core.frames.synodic_omega`
of ``"E"`` (Earth's mean motion), which is the natural frame for E-M
cyclers per spec §3. Callers analysing a Venus-anchored frame must pass
``omega`` explicitly.

References
----------

* Spec §3 (closure definition), §6 (interfaces), §10 (closure-frame risk).
* Plan: ``docs/phases/m3-model-construct/plan.md`` §3.1.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import sqrt

import numpy as np
from numpy.typing import NDArray

from cyclerfinder.core.constants import AU_KM, MU_SUN_KM3_S2
from cyclerfinder.core.frames import synodic_omega, to_rotating

Vec3 = NDArray[np.float64]  # shape (3,), dtype float64


@dataclass(frozen=True)
class Encounter:
    """A single flyby in the cycler.

    Attributes
    ----------
    body:
        One-letter body code matching
        :data:`cyclerfinder.core.constants.PLANETS` (``"V"``, ``"E"``, ``"M"``).
    t:
        Encounter epoch in seconds from the cycler's ``t = 0`` reference.
    r:
        Heliocentric inertial position of the planet at ``t``, km.
    v_planet:
        Heliocentric inertial velocity of the planet at ``t``, km/s.
    vinf_in:
        Spacecraft :math:`V_\\infty` approach vector: ``v_sc_arrive - v_planet``.
    vinf_out:
        Spacecraft :math:`V_\\infty` departure vector: ``v_sc_depart - v_planet``.

    Notes
    -----
    For the *first* and *last* encounters of an open sequence,
    :mod:`cyclerfinder.search.construct` sets ``vinf_in == vinf_out`` by
    convention (no preceding/succeeding leg exists). This makes
    :meth:`Cycler.maintenance_dv` consistent: a ballistic flyby contributes
    zero to the sum even at the boundary encounters of an open sequence.
    """

    body: str
    t: float
    r: Vec3
    v_planet: Vec3
    vinf_in: Vec3
    vinf_out: Vec3


@dataclass(frozen=True)
class Leg:
    """A heliocentric arc between two consecutive encounters.

    Attributes
    ----------
    from_body, to_body:
        One-letter body codes of the departure and arrival planets.
    t_depart, t_arrive:
        Departure/arrival epochs in seconds (same reference as
        :attr:`Encounter.t`); enforced ``t_arrive > t_depart`` in
        :func:`cyclerfinder.search.construct.construct_cycler`.
    v_depart, v_arrive:
        Spacecraft heliocentric velocities at departure/arrival, km/s.
        (Not :math:`V_\\infty` vectors — those live on the encounter.)
    n_revs:
        Number of full heliocentric revolutions on this leg. ``0`` is the
        single-revolution / direct case used in M3.
    branch:
        Lambert branch label per spec §12(b): ``"single"`` when ``n_revs == 0``
        (only one solution exists), ``"low"`` or ``"high"`` for ``n_revs >= 1``.
        Defaults to ``"single"`` because the M3 constructor only requests
        ``n_revs = 0`` solutions; M4 will exercise the multi-rev branches.
    """

    from_body: str
    to_body: str
    t_depart: float
    t_arrive: float
    v_depart: Vec3
    v_arrive: Vec3
    n_revs: int = 0
    branch: str = "single"


@dataclass(frozen=True)
class Cycler:
    """A closed-periodic flyby trajectory.

    Attributes
    ----------
    bodies:
        Ordered flyby sequence (e.g. ``["E", "M", "E"]``).
    period:
        Total period in seconds. Typically ``k * synodic_period(pair)`` from
        :func:`cyclerfinder.search.resonance.synodic_period_days`.
    encounters:
        ``len(encounters) == len(bodies)``; ordered by time.
    legs:
        ``len(legs) == len(bodies) - 1`` for an open sequence.

    No ``__eq__`` override
    ---------------------
    Frozen dataclasses generate ``__eq__`` from field-wise ``==``, which on
    ``numpy.ndarray`` returns an array and triggers an ambiguous-truth-value
    error. Two cyclers with identical state therefore compare *unequal*
    under ``==``. This is documented (not fixed) at the M3 layer; the
    catalogue's M7 signature-based identity is the correct identity for
    cyclers and will be implemented there.
    """

    bodies: list[str]
    period: float
    encounters: list[Encounter]
    legs: list[Leg]

    # ------------------------------------------------------------------
    # Metric methods
    # ------------------------------------------------------------------

    def maintenance_dv(self) -> float:
        """Sum of per-encounter velocity discontinuities, km/s.

        ``sum_i ||vinf_out_i - vinf_in_i||``. M3 returns the raw velocity
        mismatch sum without invoking
        :func:`cyclerfinder.core.flyby.flyby_dv` — the bend-feasibility
        version (which discounts a bend a real flyby can deliver
        ballistically) arrives in M4. For an idealised cycler this is
        exactly zero by construction; the constructor enforces
        ``vinf_in == vinf_out`` at the open-sequence boundaries so a
        cyclical-but-open representation does not spuriously charge those
        endpoints.
        """
        total = 0.0
        for enc in self.encounters:
            total += float(np.linalg.norm(enc.vinf_out - enc.vinf_in))
        return total

    def closure_residual(self, omega_rad_per_s: float | None = None) -> float:
        """Velocity-only closure residual in the rotating frame, km/s.

        After one period the final encounter's spacecraft departure velocity,
        viewed in the synodic rotating frame at angular rate ``omega``,
        must equal the first encounter's spacecraft departure velocity.
        The residual is ``|| v_first_rot - v_last_rot ||`` where
        ``v_*_rot`` is obtained by transforming the inertial
        ``(r, v_sc_depart)`` to the rotating frame at the encounter's epoch.

        Parameters
        ----------
        omega_rad_per_s:
            Frame angular rate (rad/s). Defaults to
            :func:`~cyclerfinder.core.frames.synodic_omega` of ``"E"`` —
            Earth's mean motion — which is the natural frame for E-M
            cyclers. Callers analysing a Venus-anchored or VEM cycler must
            pass ``omega`` explicitly (the default is appropriate for the
            M3 Earth-Mars case only).

        Returns
        -------
        float
            Velocity residual in km/s, non-negative. Zero indicates exact
            geometric closure in the supplied rotating frame; small
            non-zero values arise from Lambert numerical error in the
            circular-coplanar regime.
        """
        if omega_rad_per_s is None:
            omega_rad_per_s = synodic_omega("E")
        first_enc = self.encounters[0]
        last_enc = self.encounters[-1]
        first_leg = self.legs[0]
        last_leg = self.legs[-1]
        # Spacecraft heliocentric velocity at first departure and last arrival.
        # We compare *departure* velocities at matching geometric phase: the
        # spacecraft leaving the first encounter on lap N+1 must equal the
        # spacecraft leaving the cycler's first encounter on lap N. Under
        # idealised cycling, "departing the last encounter at t=period" is
        # equivalent to "departing the first encounter at t=0" up to the
        # rotating-frame advance.
        _r0_rot, v0_rot = to_rotating(first_enc.r, first_leg.v_depart, first_enc.t, omega_rad_per_s)
        # Use last leg's arrival state (spacecraft state at the last encounter,
        # which closes back to the first encounter one frame-period later).
        _r1_rot, v1_rot = to_rotating(last_enc.r, last_leg.v_arrive, last_enc.t, omega_rad_per_s)
        return float(np.linalg.norm(v0_rot - v1_rot))

    def radial_span(self) -> tuple[float, float]:
        """``(min_perihelion_AU, max_aphelion_AU)`` over all legs.

        Each leg's analytic perihelion/aphelion is computed from its
        departure state ``(r, v_depart)`` via vis-viva (``a = 1/(2/r -
        v^2/mu)``) and the eccentricity-vector magnitude. The output is in
        AU to match the spec §9 anchors.

        Returns
        -------
        (min_peri_AU, max_apo_AU):
            Minimum perihelion and maximum aphelion radii encountered on
            any leg of the cycler, in AU.

        Raises
        ------
        ValueError
            If the cycler has no legs (a single-encounter "cycler" has no
            heliocentric arc to analyse).
        """
        if not self.legs:
            raise ValueError("radial_span requires at least one leg")
        peris_au: list[float] = []
        apos_au: list[float] = []
        for leg, enc in zip(self.legs, self.encounters[:-1], strict=True):
            a_km, e = _orbit_from_rv(enc.r, leg.v_depart, MU_SUN_KM3_S2)
            peris_au.append(a_km * (1.0 - e) / AU_KM)
            apos_au.append(a_km * (1.0 + e) / AU_KM)
        return (min(peris_au), max(apos_au))

    def max_vinf(self) -> float:
        """Largest :math:`|V_\\infty|` across all encounters, km/s.

        Uses ``vinf_in`` only — in a steady cycler
        ``|vinf_in| == |vinf_out|`` at each encounter (ballistic flyby
        invariant), so taking only one avoids double-counting.
        """
        return max(float(np.linalg.norm(enc.vinf_in)) for enc in self.encounters)


# ---------------------------------------------------------------------------
# Local helpers
# ---------------------------------------------------------------------------


def _orbit_from_rv(r: Vec3, v: Vec3, mu: float) -> tuple[float, float]:
    """Return ``(a_km, e)`` for a heliocentric state via vis-viva + e-vector.

    A small algebra helper used by :meth:`Cycler.radial_span`. Not promoted to
    the public ``core`` surface in M3 because a proper "orbital elements"
    module belongs to M5 once we need inclination / arg-of-periapsis as well.

    Parameters
    ----------
    r, v:
        Heliocentric state, km and km/s.
    mu:
        Gravitational parameter of the central body, km^3/s^2.

    Returns
    -------
    (a_km, e):
        Semi-major axis (km, signed: positive for elliptic) and eccentricity
        magnitude. Eccentricity returned as ``[0, 1)`` for elliptic, ``>= 1``
        for hyperbolic.
    """
    r_n = float(np.linalg.norm(r))
    v_n = float(np.linalg.norm(v))
    # Vis-viva: 1/a = 2/r - v^2/mu.
    inv_a = 2.0 / r_n - v_n * v_n / mu
    a_km = 1.0 / inv_a
    # Eccentricity vector: e = ((v^2 - mu/r) * r - (r . v) * v) / mu.
    rdotv = float(np.dot(r, v))
    e_vec = ((v_n * v_n - mu / r_n) * r - rdotv * v) / mu
    e = float(np.linalg.norm(e_vec))
    # Catch near-parabolic / numerical-overflow cases sanely; in M3 we don't
    # build hyperbolic cyclers, so just return e and let radial_span absorb it.
    return a_km, e


def orbit_elements_au(r: Vec3, v: Vec3, mu: float = MU_SUN_KM3_S2) -> tuple[float, float]:
    """Return ``(a_AU, e)`` for a heliocentric state.

    Thin AU-scaled wrapper around :func:`_orbit_from_rv`. Exposed at module
    scope (rather than buried in the test file) because both
    :meth:`Cycler.radial_span` and the Aldrin gate test consume the same
    formula and the leak would otherwise be a copy of the algebra.

    Parameters
    ----------
    r, v:
        Heliocentric state, km and km/s.
    mu:
        Gravitational parameter, km^3/s^2. Defaults to the solar value.

    Returns
    -------
    (a_AU, e):
        Semi-major axis in astronomical units; eccentricity dimensionless.
    """
    a_km, e = _orbit_from_rv(r, v, mu)
    return a_km / AU_KM, e


__all__ = ["Cycler", "Encounter", "Leg", "orbit_elements_au"]


# Defensive: import-time sanity that helper math is well-defined for a
# canonical Earth circular orbit (a = 1 AU, e = 0). Cheap and catches a
# refactor that flips a sign.
def _self_check() -> None:
    a_km = 1.0 * AU_KM
    v_circ = sqrt(MU_SUN_KM3_S2 / a_km)
    r = np.array([a_km, 0.0, 0.0], dtype=np.float64)
    v = np.array([0.0, v_circ, 0.0], dtype=np.float64)
    a_au, e = orbit_elements_au(r, v)
    assert abs(a_au - 1.0) < 1.0e-12, a_au
    assert e < 1.0e-12, e


_self_check()
