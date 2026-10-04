"""Demanded-turn gate for patched-conic flyby chains (#888).

A patched-conic chain is flyable without propellant only if, at every flyby,
the spacecraft can be turned from the incoming V-infinity vector onto the
outgoing one by the body's gravity alone. That needs two things:

1. equal magnitudes (an unpowered hyperbola about the body keeps |V_inf|), and
2. a DEMANDED turn (the angle between the incoming and the outgoing vector)
   no larger than the largest bend the body can supply at the altitude floor.

This module checks the second condition explicitly, and reports the first
separately. It exists because the project's earlier gates did not: the #324
physical gate (:mod:`cyclerfinder.search.physical_sanity`) tests only that the
bend CAPACITY ``max_bend`` at the floor exceeds 5 degrees, and the moon-tour
validation lane matched V-infinity MAGNITUDES and restarted every leg from its
own Lambert solution, so no incoming vector was ever compared with an outgoing
one. Six catalogued Uranian rows passed that lane while demanding 1.8 to 28
times the available bend (#888; ``docs/notes/2026-10-04-888-demanded-turn-gate.md``).

Definitions (per encounter)
---------------------------
``demanded_turn``
    ``arccos(v_in . v_out / (|v_in| |v_out|))`` in ``[0, pi]``: frame-invariant,
    any common frame works (inertial, rotating, body-local), provided both
    vectors are expressed in it. :func:`cyclerfinder.core.flyby.bend_angle`.
``available_bend``
    ``2 asin(1 / (1 + r_p v^2 / GM))`` at ``r_p = R + alt_floor`` and
    ``v = min(|v_in|, |v_out|)`` (the smaller speed gives the larger, i.e. the
    more generous, bend; the magnitude mismatch is reported separately).
    :func:`cyclerfinder.core.flyby.max_bend`.
``ratio``
    ``demanded / available``; > 1 means the flyby cannot be flown unpowered.
``required_alt_km``
    The periapsis altitude at which an unpowered flyby at ``v`` turns by
    exactly ``demanded``: ``r_p = GM/v^2 (1/sin(demanded/2) - 1)``, minus R.
    Below the floor (possibly negative, i.e. inside the body) when infeasible;
    ``+inf`` when the demanded turn is zero.
``impulse_beyond_bend_kms``
    Minimum-impulse estimate when the flyby is powered. MODEL: a single
    impulse applied outside the sphere of influence after (equivalently
    before) an unpowered flyby that rotates ``v_in`` by up to
    ``available_bend`` towards ``v_out`` in the plane of the two vectors;
    the impulse is ``|v_out - R(min(demanded, available)) v_in|``. With equal
    magnitudes this is ``2 v sin((demanded - available)/2)`` (Strange and
    Longuski 2002, eq. 9; :func:`cyclerfinder.core.flyby.dv_from_turn_deficit`),
    and it also pays any magnitude mismatch. It is the generalisation of
    ``uranus_bridge_885.entry_correction_kms`` to 3-D.
``impulse_periapsis_kms``
    A cheaper alternative estimate, a tangential impulse pair at periapsis
    (Oberth credit, :func:`cyclerfinder.core.flyby.dv_powered_flyby_periapsis`);
    for equal-magnitude encounters only (``nan`` otherwise).

Overall verdict: ``turn_feasible`` iff every encounter's demanded turn is no
larger than its available bend (to ``turn_tol_rad``); ``ballistic`` iff in
addition every magnitude mismatch is within ``mag_tol_kms``.

Three-way ``status`` (#937, #906): ``"pass"`` / ``"fail"`` follow
``turn_feasible`` except where the patched conic itself cannot decide, which
is ``"indeterminate"`` (never a rejection):

* the TIDAL TURN SCALE ``(n r_soi / v)^2`` (n the body's mean motion about its
  primary, ``r_soi`` its Laplace sphere of influence, v the smaller
  V-infinity) is a radian or more: at such a V-infinity the primary's tide
  over the transit of the sphere turns the velocity as much as the body does,
  and the patched-conic V-infinity has no meaning (the slow-pass regime of
  Henon, Brjuno and the low-energy transfers);
* ``|demanded - available|`` is within ``TIDAL_BAND_FACTOR`` tidal turn scales:
  the verdict could be flipped by the part of the turn the patched conic
  leaves out. The factor 3 is a convention set from the full-model controls
  of #937 (the non-body part of the integrated turn of every hyperbolic lunar
  pass of Casoliva's 7-3a/b/c and Leiva-Briozzo 2005 is 0.3 to 1.9 tidal
  scales; ``verify/turn_gate_fullmodel.py``).

``turn_feasible`` itself is unchanged, so every existing caller keeps its
verdict. Full-model orbits (CR3BP, bicircular) must be measured with
:mod:`cyclerfinder.verify.turn_gate_fullmodel`, which converts the rotating-
frame velocities to inertial axes before differencing them; differencing
rotating-frame components taken a day apart adds the frame's own rotation
(13.2 degrees per day in the Earth-Moon system), the error of the #899 review.

Positive controls and regression: ``tests/verify/test_turn_gate.py``
(McConaghy-Longuski-Byrnes 2002 Table 4; Russell-Strange 2009 Tables 3-6;
the six withdrawn Uranian rows).
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from cyclerfinder.core.constants import AU_KM, MU_SUN_KM3_S2, PLANETS
from cyclerfinder.core.flyby import (
    bend_angle,
    dv_powered_flyby_periapsis,
    max_bend,
)
from cyclerfinder.core.satellites import SATELLITES

Vec = NDArray[np.float64]

#: Default magnitude-mismatch tolerance for the ``ballistic`` verdict, km/s.
#: The #558/#563 closure gate is 0.05 km/s; symmetric closures match to 1e-12.
DEFAULT_MAG_TOL_KMS: float = 1.0e-3
#: Slack on ``demanded <= available`` (rad), absorbing round-off only.
DEFAULT_TURN_TOL_RAD: float = 1.0e-9
#: ``status`` is ``"indeterminate"`` when ``|demanded - available|`` is within this many
#: tidal turn scales (convention from the #937 full-model controls; module docstring).
TIDAL_BAND_FACTOR: float = 3.0
#: ``status`` is ``"indeterminate"`` when the tidal turn scale reaches this (rad).
TIDAL_SCALE_LIMIT_RAD: float = 1.0


def tidal_speed_kms(body: str) -> float:
    """``n r_soi`` (km/s) of a flyby body about its primary: the speed scale of the
    tidal turn ``(n r_soi / v)^2``. ``nan`` when the body is not in the registries."""
    pl = PLANETS.get(body)
    if pl is None:
        pl = next((p for p in PLANETS.values() if p.name == body), None)
    if pl is not None:
        a_km = pl.sma_au * AU_KM
        gm_ratio = pl.mu_km3_s2 / MU_SUN_KM3_S2
        n = math.radians(pl.mean_motion_deg_day) / 86400.0
        return float(n * a_km * gm_ratio**0.4)
    s = SATELLITES.get(body)
    if s is None:
        return math.nan
    primary = next((p for p in PLANETS.values() if p.name == s.primary), None)
    if primary is None:
        return math.nan
    n = math.radians(s.mean_motion_deg_day) / 86400.0
    return float(n * s.sma_km * (s.mu_km3_s2 / primary.mu_km3_s2) ** 0.4)


@dataclass(frozen=True)
class BodyConstants:
    """GM, equatorial radius and project altitude floor of a flyby body."""

    name: str
    mu_km3_s2: float
    radius_km: float
    alt_floor_km: float


def body_constants(body: str) -> BodyConstants:
    """Resolve a body from the project registries.

    Heliocentric planets by one-letter code (``"E"``, ``"M"``, ...) or name
    (``"Earth"``) from :data:`cyclerfinder.core.constants.PLANETS`; moons by
    full name from :data:`cyclerfinder.core.satellites.SATELLITES`. The floor is
    the registry ``safe_alt_km`` (sourced per body there, e.g. Earth 200 km from
    Russell 2004 p.165; Uranian tour moons 50 km; Titan 1500 km).
    """
    if body in PLANETS:
        pl = PLANETS[body]
        return BodyConstants(pl.name, pl.mu_km3_s2, pl.radius_eq_km, pl.safe_alt_km)
    for pl in PLANETS.values():
        if pl.name == body:
            return BodyConstants(pl.name, pl.mu_km3_s2, pl.radius_eq_km, pl.safe_alt_km)
    if body in SATELLITES:
        s = SATELLITES[body]
        return BodyConstants(s.name, s.mu_km3_s2, s.radius_eq_km, s.safe_alt_km)
    raise KeyError(f"unknown flyby body {body!r} (not in PLANETS or SATELLITES)")


@dataclass(frozen=True)
class Encounter:
    """One flyby of a patched-conic chain: incoming and outgoing V-infinity.

    ``vinf_in`` / ``vinf_out`` are body-relative velocities (km/s) in any
    common frame, 2-D or 3-D.
    """

    body: str
    mu_km3_s2: float
    radius_km: float
    alt_floor_km: float
    vinf_in: Vec
    vinf_out: Vec
    label: str = ""
    #: ``n r_soi`` of the body (km/s), for the tidal turn scale; ``nan`` disables
    #: the ``"indeterminate"`` band (``status`` then follows ``turn_feasible``).
    tidal_speed_kms: float = math.nan

    @classmethod
    def for_body(
        cls,
        body: str,
        vinf_in: Sequence[float] | Vec,
        vinf_out: Sequence[float] | Vec,
        *,
        alt_floor_km: float | None = None,
        label: str = "",
    ) -> Encounter:
        """Build from the registry constants; ``alt_floor_km=None`` uses the project floor."""
        bc = body_constants(body)
        return cls(
            body=body,
            mu_km3_s2=bc.mu_km3_s2,
            radius_km=bc.radius_km,
            alt_floor_km=bc.alt_floor_km if alt_floor_km is None else float(alt_floor_km),
            vinf_in=_vec3(vinf_in),
            vinf_out=_vec3(vinf_out),
            label=label,
            tidal_speed_kms=tidal_speed_kms(body),
        )

    def with_floor(self, alt_floor_km: float) -> Encounter:
        """The same encounter judged at a different altitude floor."""
        return Encounter(
            self.body,
            self.mu_km3_s2,
            self.radius_km,
            float(alt_floor_km),
            self.vinf_in,
            self.vinf_out,
            self.label,
            self.tidal_speed_kms,
        )


@dataclass(frozen=True)
class EncounterTurn:
    """Gate result for one encounter (angles in degrees, speeds in km/s)."""

    body: str
    label: str
    vinf_in_kms: float
    vinf_out_kms: float
    magnitude_mismatch_kms: float
    demanded_turn_deg: float
    available_bend_deg: float
    alt_floor_km: float
    ratio: float
    required_alt_km: float
    impulse_beyond_bend_kms: float
    impulse_periapsis_kms: float
    turn_feasible: bool
    tidal_turn_scale_deg: float = math.nan
    status: str = ""

    def as_dict(self) -> dict[str, object]:
        return {
            "body": self.body,
            "label": self.label,
            "vinf_in_kms": self.vinf_in_kms,
            "vinf_out_kms": self.vinf_out_kms,
            "magnitude_mismatch_kms": self.magnitude_mismatch_kms,
            "demanded_turn_deg": self.demanded_turn_deg,
            "available_bend_deg": self.available_bend_deg,
            "alt_floor_km": self.alt_floor_km,
            "ratio": _json_float(self.ratio),
            "required_alt_km": _json_float(self.required_alt_km),
            "impulse_beyond_bend_kms": self.impulse_beyond_bend_kms,
            "impulse_periapsis_kms": _json_float(self.impulse_periapsis_kms),
            "turn_feasible": self.turn_feasible,
            "tidal_turn_scale_deg": _json_float(self.tidal_turn_scale_deg),
            "status": self.status,
        }


@dataclass(frozen=True)
class TurnGateReport:
    """Gate result for a whole chain."""

    encounters: tuple[EncounterTurn, ...]
    mag_tol_kms: float

    @property
    def turn_feasible(self) -> bool:
        """Every demanded turn fits inside the available bend."""
        return all(e.turn_feasible for e in self.encounters)

    @property
    def status(self) -> str:
        """``"fail"`` if any encounter fails, else ``"indeterminate"`` if any is, else
        ``"pass"`` (module docstring)."""
        st = {e.status for e in self.encounters}
        if "fail" in st:
            return "fail"
        if "indeterminate" in st:
            return "indeterminate"
        return "pass"

    @property
    def max_magnitude_mismatch_kms(self) -> float:
        return max((e.magnitude_mismatch_kms for e in self.encounters), default=0.0)

    @property
    def ballistic(self) -> bool:
        """Turn-feasible AND every magnitude matches within ``mag_tol_kms``."""
        return self.turn_feasible and self.max_magnitude_mismatch_kms <= self.mag_tol_kms

    @property
    def worst_ratio(self) -> float:
        return max((e.ratio for e in self.encounters), default=0.0)

    @property
    def binding(self) -> EncounterTurn | None:
        """The encounter with the largest demanded/available ratio."""
        return max(self.encounters, key=lambda e: e.ratio, default=None)

    @property
    def min_required_alt_km(self) -> float:
        return min((e.required_alt_km for e in self.encounters), default=math.inf)

    @property
    def total_impulse_beyond_bend_kms(self) -> float:
        return float(sum(e.impulse_beyond_bend_kms for e in self.encounters))

    def as_dict(self) -> dict[str, object]:
        return {
            "turn_feasible": self.turn_feasible,
            "status": self.status,
            "ballistic": self.ballistic,
            "worst_ratio": _json_float(self.worst_ratio),
            "max_magnitude_mismatch_kms": self.max_magnitude_mismatch_kms,
            "min_required_alt_km": _json_float(self.min_required_alt_km),
            "total_impulse_beyond_bend_kms": self.total_impulse_beyond_bend_kms,
            "encounters": [e.as_dict() for e in self.encounters],
        }


# ---------------------------------------------------------------------------
# Primitives
# ---------------------------------------------------------------------------


def _vec3(v: Sequence[float] | Vec) -> Vec:
    a = np.asarray(v, dtype=np.float64).reshape(-1)
    if a.size == 2:
        return np.array([a[0], a[1], 0.0])
    if a.size != 3:
        raise ValueError(f"V-infinity must have 2 or 3 components, got {a.size}")
    return a


def _json_float(x: float) -> float | str:
    if math.isinf(x):
        return "inf" if x > 0 else "-inf"
    if math.isnan(x):
        return "nan"
    return float(x)


def available_bend_rad(mu_km3_s2: float, radius_km: float, alt_km: float, vinf_kms: float) -> float:
    """Largest unpowered bend at periapsis altitude ``alt_km`` (wraps ``core.flyby.max_bend``)."""
    return max_bend(mu_km3_s2, radius_km + alt_km, vinf_kms)


def required_periapsis_alt_km(
    mu_km3_s2: float, radius_km: float, vinf_kms: float, turn_rad: float
) -> float:
    """Periapsis altitude at which an unpowered flyby turns ``vinf`` by ``turn_rad``.

    Inverse of the bend formula: ``r_p = GM/v^2 (1/sin(turn/2) - 1)``. Returns
    ``+inf`` for a zero turn, and ``-radius_km`` (periapsis at the centre) for a
    180 degree turn; at ``vinf == 0`` any turn is available at any altitude
    (returns ``+inf``).
    """
    if turn_rad <= 0.0 or vinf_kms <= 0.0:
        return math.inf
    s = math.sin(0.5 * min(turn_rad, math.pi))
    r_p = mu_km3_s2 / (vinf_kms * vinf_kms) * (1.0 / s - 1.0)
    return r_p - radius_km


def rotate_towards(v_from: Vec, v_to: Vec, angle_rad: float) -> Vec:
    """Rotate ``v_from`` by ``angle_rad`` towards ``v_to`` in their common plane.

    For anti-parallel vectors the plane is undefined; any perpendicular axis
    is used (every choice gives the same distance to ``v_to``).
    """
    a = np.asarray(v_from, dtype=np.float64)
    b = np.asarray(v_to, dtype=np.float64)
    axis = np.cross(a, b)
    n = float(np.linalg.norm(axis))
    if n < 1e-300:
        # Parallel or anti-parallel: pick any axis perpendicular to a.
        trial = (
            np.array([1.0, 0.0, 0.0])
            if abs(a[0]) < 0.9 * np.linalg.norm(a)
            else np.array([0.0, 1.0, 0.0])
        )
        axis = np.cross(a, trial)
        n = float(np.linalg.norm(axis))
        if n == 0.0:
            return a.copy()
    k = axis / n
    c, s = math.cos(angle_rad), math.sin(angle_rad)
    # Rodrigues' rotation formula.
    return a * c + np.cross(k, a) * s + k * float(np.dot(k, a)) * (1.0 - c)


def impulse_beyond_bend_kms(v_in: Vec, v_out: Vec, bend_rad: float) -> float:
    """Single-impulse estimate after a maximal unpowered bend (see module docstring)."""
    a = np.asarray(v_in, dtype=np.float64)
    b = np.asarray(v_out, dtype=np.float64)
    if float(np.linalg.norm(a)) == 0.0 or float(np.linalg.norm(b)) == 0.0:
        return float(np.linalg.norm(b - a))
    turn = bend_angle(a, b)
    rotated = rotate_towards(a, b, min(turn, bend_rad))
    return float(np.linalg.norm(b - rotated))


def evaluate_encounter(
    enc: Encounter,
    *,
    turn_tol_rad: float = DEFAULT_TURN_TOL_RAD,
    mag_tol_kms: float = DEFAULT_MAG_TOL_KMS,
) -> EncounterTurn:
    """Apply the gate to one encounter."""
    v_in = _vec3(enc.vinf_in)
    v_out = _vec3(enc.vinf_out)
    m_in = float(np.linalg.norm(v_in))
    m_out = float(np.linalg.norm(v_out))
    if m_in == 0.0 or m_out == 0.0:
        raise ValueError(f"zero V-infinity at {enc.body} {enc.label!r}: no flyby is defined")
    v = min(m_in, m_out)
    demanded = bend_angle(v_in, v_out)
    available = available_bend_rad(enc.mu_km3_s2, enc.radius_km, enc.alt_floor_km, v)
    if demanded <= turn_tol_rad:
        ratio = 0.0
    elif available <= 0.0:
        ratio = math.inf
    else:
        ratio = demanded / available
    feasible = demanded <= available + turn_tol_rad
    if abs(m_in - m_out) <= mag_tol_kms:
        periapsis_dv = dv_powered_flyby_periapsis(
            v, demanded, available, enc.mu_km3_s2, enc.radius_km + enc.alt_floor_km
        )
    else:
        periapsis_dv = math.nan
    tidal = (enc.tidal_speed_kms / v) ** 2
    if not math.isnan(tidal) and (
        tidal >= TIDAL_SCALE_LIMIT_RAD or abs(demanded - available) <= TIDAL_BAND_FACTOR * tidal
    ):
        status = "indeterminate"
    else:
        status = "pass" if feasible else "fail"
    return EncounterTurn(
        body=enc.body,
        label=enc.label,
        vinf_in_kms=m_in,
        vinf_out_kms=m_out,
        magnitude_mismatch_kms=abs(m_in - m_out),
        demanded_turn_deg=math.degrees(demanded),
        available_bend_deg=math.degrees(available),
        alt_floor_km=enc.alt_floor_km,
        ratio=ratio,
        required_alt_km=required_periapsis_alt_km(enc.mu_km3_s2, enc.radius_km, v, demanded),
        impulse_beyond_bend_kms=impulse_beyond_bend_kms(v_in, v_out, available),
        impulse_periapsis_kms=periapsis_dv,
        turn_feasible=bool(feasible),
        tidal_turn_scale_deg=math.degrees(tidal),
        status=status,
    )


def demanded_turn_gate(
    encounters: Sequence[Encounter],
    *,
    turn_tol_rad: float = DEFAULT_TURN_TOL_RAD,
    mag_tol_kms: float = DEFAULT_MAG_TOL_KMS,
    alt_floor_km: float | None = None,
) -> TurnGateReport:
    """Apply the demanded-turn gate to every encounter of a chain.

    ``alt_floor_km`` overrides every encounter's own floor (e.g. a uniform
    50 km re-judgement); ``None`` keeps each encounter's floor.
    """
    encs = [e if alt_floor_km is None else e.with_floor(alt_floor_km) for e in encounters]
    return TurnGateReport(
        encounters=tuple(
            evaluate_encounter(e, turn_tol_rad=turn_tol_rad, mag_tol_kms=mag_tol_kms) for e in encs
        ),
        mag_tol_kms=mag_tol_kms,
    )


# ---------------------------------------------------------------------------
# Frames and adapters
# ---------------------------------------------------------------------------


def to_body_local(v: Sequence[float] | Vec, r_body: Vec, v_body: Vec) -> Vec:
    """Components of ``v`` along (radial, along-track, orbit-normal) of a body.

    Two V-infinity vectors measured at different epochs of a periodic chain
    (e.g. the arrival that closes one cycle and the departure that opens the
    next) can be compared in this frame when the chain repeats in the body's
    rotating frame, which is the case for circular-coplanar ideal models.
    """
    r = np.asarray(r_body, dtype=np.float64)
    w = np.asarray(v_body, dtype=np.float64)
    r_hat = r / np.linalg.norm(r)
    h = np.cross(r, w)
    n_hat = h / np.linalg.norm(h)
    t_hat = np.cross(n_hat, r_hat)
    vv = _vec3(v)
    return np.array([float(vv @ r_hat), float(vv @ t_hat), float(vv @ n_hat)])


def encounters_from_vinf_nodes(
    nodes: Mapping[str, Sequence[float] | Vec],
    sequence: Sequence[str],
    *,
    wrap_rotation_rad: float | None = None,
    alt_floor_km: float | None = None,
) -> list[Encounter]:
    """Encounters from the ``b{i}_in`` / ``b{i}_out`` node dict of
    :func:`cyclerfinder.search.correct._vinf_nodes`.

    The intermediate encounters ``1 .. n-2`` are flybys. With
    ``wrap_rotation_rad`` given, the periodicity wrap is appended: the arrival
    ``b{n-1}_in`` is compared with the departure ``b0_out`` rotated about +z by
    that angle (the home body's advance over one period, as in
    :func:`cyclerfinder.search.turn_ratio_check.wrap_node_turn`).
    """
    out: list[Encounter] = []
    last = len(sequence) - 1
    for i in range(1, last):
        out.append(
            Encounter.for_body(
                sequence[i],
                nodes[f"b{i}_in"],
                nodes[f"b{i}_out"],
                alt_floor_km=alt_floor_km,
                label=f"b{i}",
            )
        )
    if wrap_rotation_rad is not None:
        c, s = math.cos(wrap_rotation_rad), math.sin(wrap_rotation_rad)
        rot = np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])
        v_out = rot @ _vec3(nodes["b0_out"])
        out.append(
            Encounter.for_body(
                sequence[last],
                nodes[f"b{last}_in"],
                v_out,
                alt_floor_km=alt_floor_km,
                label=f"b{last}-wrap",
            )
        )
    return out
