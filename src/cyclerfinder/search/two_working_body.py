"""Two-working-body cycler generator and date-residual corrector (#942 R1, #943 X1).

A cycler here is a closed chain of patched-conic legs between the encounters
of TWO massive ("working") bodies, either of which may host returns and take
part of the turn. It generalises Russell's Earth-only generic-return model
(:mod:`cyclerfinder.search.generic_return`), where the second body is a
massless target, to the class Hollister & Menning 1970 solved for Earth and
Venus ("Periodic Swing-By Orbits between Earth and Venus", J. Spacecraft 7(10),
pp. 1193-1198).

Method (Hollister & Menning 1970, pp. 1194-1195)
-----------------------------------------------
The trajectory is fixed by the encounter dates. Legs are of two kinds:

* LAMBERT legs (interplanetary transfers and same-body "symmetric"/generic
  returns), whose start dates are the unknowns;
* FIXED legs of known duration. A full-revolution (``n:m`` resonant) return
  lasts ``n`` body periods; its V-infinity has the block's magnitude and a
  direction that is free on a circle (H&M p.1194, "a double infinity of such
  orbits"). A half-revolution (n-pi) return lasts an odd number of body
  half-periods and has discrete directions at a given V-infinity.

"The number of independent dates N equals the total number of interplanetary
transfers and symmetric returns" (p.1194); the N residuals are the V-infinity
magnitude differences at the N junctions between consecutive Lambert legs
(H&M "differences in relative velocity magnitudes at each flyby"). A massless
encounter (Russell's target, R1 cell (c)'s Mars) contributes the full vector
difference instead (the spacecraft passes it on one conic). The square
system is solved by least squares on the dates. The direction of each fixed
full-revolution leg is then chosen to minimise the largest turn ratio of the
flybys that bracket it (a minimax over the free circle angles), and every
flyby is judged by the #888/#937 demanded-turn gate
(:mod:`cyclerfinder.verify.turn_gate`) on inertial-axis vectors.

The ephemeris is swappable: :class:`CircularSystem` is the circular-coplanar
ideal model; :class:`MeanElementSystem` is the fixed-element inclined-elliptic
model used for the Hollister-Menning positive control. Controls and candidates
run through the same residual and gate code.
"""

from __future__ import annotations

import math
from collections.abc import Callable, Sequence
from dataclasses import dataclass, field
from typing import Protocol

import numpy as np
from numpy.typing import NDArray
from scipy.optimize import least_squares, minimize

from cyclerfinder.core.constants import AU_KM, MU_SUN_KM3_S2, PLANETS, SECONDS_PER_DAY
from cyclerfinder.core.ephemeris import inclined_planets
from cyclerfinder.core.flyby import bend_angle
from cyclerfinder.core.kepler import propagate
from cyclerfinder.core.lambert import LambertError, lambert
from cyclerfinder.core.satellites import SATELLITES
from cyclerfinder.verify.turn_gate import (
    Encounter,
    TurnGateReport,
    available_bend_rad,
    body_constants,
    demanded_turn_gate,
    tidal_speed_kms,
)

Vec = NDArray[np.float64]
#: Encounter start dates (s): a sequence or an array.
Dates = Sequence[float] | NDArray[np.float64]

#: Earth Mean Orbital Speed, km/s: H&M's V-infinity unit (Table 3 "EMOS");
#: the value used by the project's transcription
#: ``data/sources/hollister-menning-1970-table3.yaml``.
EMOS_KMS: float = 29.785
#: A demanded turn within this of 180 degrees is a rejection (owner ruling
#: 2026-10-05, recorded under #937 / #906).
NEAR_180_DEG: float = 1.0


# ---------------------------------------------------------------------------
# Ephemerides
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class FlybyBody:
    """Physical constants of an encounter body (km, km^3/s^2)."""

    code: str
    mu_km3_s2: float
    radius_km: float
    alt_floor_km: float
    massless: bool = False


class System(Protocol):
    """A central body plus encounter bodies with a state function."""

    mu: float

    def state(self, code: str, t_s: float) -> tuple[Vec, Vec]: ...

    def period_s(self, code: str) -> float: ...

    def body(self, code: str) -> FlybyBody: ...

    def wrap_rotation(self, code: str, dt_s: float) -> Vec: ...


def _flyby_body(code: str, massless: bool) -> FlybyBody:
    bc = body_constants(code)
    return FlybyBody(code, bc.mu_km3_s2, bc.radius_km, bc.alt_floor_km, massless)


@dataclass
class CircularSystem:
    """Circular-coplanar ideal model.

    ``bodies`` maps a code to ``(a_km, period_s, theta0_rad)``: each body rides
    a prograde circle in the xy-plane at angle ``theta0 + 2 pi t / period``.
    ``massless`` names target bodies that take no part in the turn;
    ``flyby_overrides`` replaces registry flyby constants (e.g. a paper's own).
    """

    mu: float
    bodies: dict[str, tuple[float, float, float]]
    massless: frozenset[str] = frozenset()
    flyby_overrides: dict[str, FlybyBody] = field(default_factory=dict)
    _fb: dict[str, FlybyBody] = field(default_factory=dict, repr=False)

    def state(self, code: str, t_s: float) -> tuple[Vec, Vec]:
        a, per, th0 = self.bodies[code]
        th = th0 + 2.0 * math.pi * t_s / per
        v = math.sqrt(self.mu / a)
        return (
            np.array([a * math.cos(th), a * math.sin(th), 0.0]),
            np.array([-v * math.sin(th), v * math.cos(th), 0.0]),
        )

    def period_s(self, code: str) -> float:
        return self.bodies[code][1]

    def body(self, code: str) -> FlybyBody:
        if code not in self._fb:
            if code in self.flyby_overrides:
                fb = self.flyby_overrides[code]
                self._fb[code] = FlybyBody(
                    fb.code, fb.mu_km3_s2, fb.radius_km, fb.alt_floor_km, code in self.massless
                )
            else:
                self._fb[code] = _flyby_body(code, code in self.massless)
        return self._fb[code]

    def wrap_rotation(self, code: str, dt_s: float) -> Vec:
        """The model is invariant under rotation: after ``dt_s`` every body has
        advanced, and when ``dt_s`` is a multiple of their synodic period the
        whole configuration is the start rotated by the advance of ``code``."""
        return _rot_z(2.0 * math.pi * dt_s / self.period_s(code))

    def synodic_s(self, a: str, b: str) -> float:
        return 1.0 / abs(1.0 / self.period_s(a) - 1.0 / self.period_s(b))


def heliocentric_circular(
    periods_yr: dict[str, float],
    *,
    massless: Sequence[str] = (),
    theta0_rad: dict[str, float] | None = None,
) -> CircularSystem:
    """Russell-style heliocentric ideal model in km/s units.

    Earth is 1 AU with a 1-yr period (365.25 d, :data:`generic_return.YEAR_DAYS`);
    every other body has ``a = P_yr^(2/3)`` AU (Kepler III with Earth as the
    unit), the convention of :meth:`RussellModel.sma_au`. ``mu`` is set so that
    the 1-AU circle has exactly the 1-yr period, which keeps the period ratios
    exact (Russell's 1.875-yr Mars is load-bearing for his tables).
    """
    year_s = 365.25 * SECONDS_PER_DAY
    mu = 4.0 * math.pi**2 * AU_KM**3 / year_s**2
    th = theta0_rad or {}
    bodies = {
        c: (AU_KM * p ** (2.0 / 3.0), p * year_s, th.get(c, 0.0)) for c, p in periods_yr.items()
    }
    return CircularSystem(mu, bodies, frozenset(massless))


def moon_circular(
    primary_mu: float,
    moons: Sequence[str],
    *,
    massless: Sequence[str] = (),
    theta0_rad: dict[str, float] | None = None,
) -> CircularSystem:
    """Planet-centred ideal model of moons on circles (registry sma, Kepler-III period)."""
    th = theta0_rad or {}
    bodies = {}
    for m in moons:
        a = SATELLITES[m].sma_km
        bodies[m] = (a, 2.0 * math.pi * math.sqrt(a**3 / primary_mu), th.get(m, 0.0))
    return CircularSystem(primary_mu, bodies, frozenset(massless))


@dataclass
class MeanElementSystem:
    """Fixed mean-element Keplerian ephemeris (inclined, elliptic), heliocentric.

    Elements at J2000 from Standish & Williams, "Approximate Positions of the
    Planets" Table 1 (a, e, varpi, L from :data:`PLANETS`; i, Omega from
    :func:`cyclerfinder.core.ephemeris.inclined_planets`). The mean motion is
    Kepler III from ``a`` with the project ``MU_SUN``. No secular rates: over
    the 1970-1986 span of the H&M control they move the orbits by about 0.1
    degree, below the paper's own one-day date rounding.
    ``t_s`` is seconds past JD 2440000.0 (H&M's "Julian - 2440000" day zero).

    ``periods_days`` optionally overrides the period of a body (the mean
    longitude is then exact at ``anchor_jd`` and advances at the overridden
    rate). H&M "assume exact periodicity of the solar system" over 16 yr
    (p.1194), i.e. Earth 5844/16 d and Venus 5844/26 d.
    """

    mu: float = MU_SUN_KM3_S2
    massless: frozenset[str] = frozenset()
    _fb: dict[str, FlybyBody] = field(default_factory=dict, repr=False)

    #: JD of ``t_s = 0``.
    jd0: float = 2440000.0
    periods_days: dict[str, float] | None = None
    anchor_jd: float = 2451545.0

    def _elements(self, code: str) -> tuple[float, float, float, float, float, float]:
        p = inclined_planets()[code]
        return (
            p.sma_au * AU_KM,
            p.ecc,
            math.radians(p.inc_deg),
            math.radians(p.lan_deg),
            math.radians(p.varpi_deg),
            math.radians(p.L0_deg),
        )

    def _kepler_period_s(self, code: str) -> float:
        a = self._elements(code)[0]
        return 2.0 * math.pi * math.sqrt(a**3 / self.mu)

    def period_s(self, code: str) -> float:
        if self.periods_days and code in self.periods_days:
            return self.periods_days[code] * SECONDS_PER_DAY
        return self._kepler_period_s(code)

    def state(self, code: str, t_s: float) -> tuple[Vec, Vec]:
        a, e, inc, lan, varpi, l0 = self._elements(code)
        n0 = 2.0 * math.pi / self._kepler_period_s(code)
        n = 2.0 * math.pi / self.period_s(code)
        anchor_s = (self.anchor_jd - 2451545.0) * SECONDS_PER_DAY
        t_j2000 = t_s + (self.jd0 - 2451545.0) * SECONDS_PER_DAY
        m_anom = (l0 - varpi + n0 * anchor_s + n * (t_j2000 - anchor_s)) % (2.0 * math.pi)
        ecc_anom = m_anom
        for _ in range(50):
            d = (ecc_anom - e * math.sin(ecc_anom) - m_anom) / (1.0 - e * math.cos(ecc_anom))
            ecc_anom -= d
            if abs(d) < 1e-15:
                break
        cos_e, sin_e = math.cos(ecc_anom), math.sin(ecc_anom)
        b = math.sqrt(1.0 - e * e)
        x, y = a * (cos_e - e), a * b * sin_e
        rdot = n * a / (1.0 - e * cos_e)
        vx, vy = -rdot * sin_e, rdot * b * cos_e
        argp = varpi - lan
        rot = _rot_z(lan) @ _rot_x(inc) @ _rot_z(argp)
        return rot @ np.array([x, y, 0.0]), rot @ np.array([vx, vy, 0.0])

    def body(self, code: str) -> FlybyBody:
        if code not in self._fb:
            self._fb[code] = _flyby_body(code, code in self.massless)
        return self._fb[code]

    def wrap_rotation(self, code: str, dt_s: float) -> Vec:
        """Identity: the cycle is compared in inertial axes (H&M's exact-periodicity
        assumption; with ``periods_days`` set to commensurate values it is exact)."""
        return np.eye(3)


def _rot_z(a: float) -> Vec:
    c, s = math.cos(a), math.sin(a)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])


def _rot_x(a: float) -> Vec:
    c, s = math.cos(a), math.sin(a)
    return np.array([[1.0, 0.0, 0.0], [0.0, c, -s], [0.0, s, c]])


def kepler_step(r0: Vec, v0: Vec, dt: float, mu: float) -> tuple[Vec, Vec]:
    """Two-body propagation by ``dt``.

    Elliptic orbits use the classical eccentric-anomaly Lagrange f and g
    solution (Kepler's equation by Newton, ``dt`` reduced modulo the period);
    :func:`cyclerfinder.core.kepler.propagate` fails to converge on some
    single-revolution heliocentric arcs here (reported as a papercut).
    Non-elliptic orbits fall back to :func:`propagate`.
    """
    r0 = np.asarray(r0, dtype=np.float64)
    v0 = np.asarray(v0, dtype=np.float64)
    r0n = float(np.linalg.norm(r0))
    alpha = 2.0 / r0n - float(v0 @ v0) / mu
    if alpha <= 0.0:
        r, v = propagate(r0, v0, dt, mu)
        return np.asarray(r, dtype=np.float64), np.asarray(v, dtype=np.float64)
    a = 1.0 / alpha
    n = math.sqrt(mu / a**3)
    sigma0 = float(r0 @ v0) / math.sqrt(mu * a)  # e sin E0
    ecos0 = 1.0 - r0n / a  # e cos E0
    e0 = math.atan2(sigma0, ecos0)
    ecc = math.hypot(sigma0, ecos0)
    m0 = e0 - sigma0
    dm = math.fmod(n * dt, 2.0 * math.pi)
    m1 = m0 + dm
    big_e = m1 if ecc < 0.8 else math.pi * (1.0 if math.sin(m1) >= 0 else -1.0) + m1 - math.pi
    for _ in range(100):
        f = big_e - ecc * math.sin(big_e) - m1
        d = f / (1.0 - ecc * math.cos(big_e))
        big_e -= d
        if abs(d) < 1e-15:
            break
    de = big_e - e0
    r1n = a + (r0n - a) * math.cos(de) + sigma0 * a * math.sin(de)
    f_c = 1.0 - a / r0n * (1.0 - math.cos(de))
    g_c = (dm - (de - math.sin(de))) / n
    fdot = -math.sqrt(mu * a) / (r1n * r0n) * math.sin(de)
    gdot = 1.0 - a / r1n * (1.0 - math.cos(de))
    return f_c * r0 + g_c * v0, fdot * r0 + gdot * v0


# ---------------------------------------------------------------------------
# Legs
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class LambertLeg:
    """A conic from ``frm`` to ``to`` with ``nrev`` complete revolutions.

    ``branch`` is ``"single"`` for ``nrev == 0`` and ``"low"``/``"high"`` for the
    two multi-revolution solutions of :func:`cyclerfinder.core.lambert.lambert`.
    A same-body leg with ``nrev >= 1`` is a generic (H&M "symmetric") return.
    """

    frm: str
    to: str
    nrev: int = 0
    branch: str = "single"

    @property
    def fixed_tof(self) -> bool:
        return False


@dataclass(frozen=True)
class ResonantLeg:
    """Full-revolution return: ``body_revs`` body periods, ``sc_revs`` spacecraft revolutions.

    Duration ``body_revs * P``; the spacecraft period is ``P body_revs / sc_revs``.
    Its V-infinity leaves and returns unchanged, on the circle
    ``|v_body + u| = v_F`` (Russell Eq. 2.13/2.17).
    """

    body: str
    body_revs: int = 1
    sc_revs: int = 1

    @property
    def fixed_tof(self) -> bool:
        return True


@dataclass(frozen=True)
class HalfRevLeg:
    """Half-revolution (n-pi) return after ``half_periods`` (odd) body half-periods.

    The spacecraft crosses the body orbit at the antipode with semi-latus
    rectum equal to the body radius. ``k_sc`` is the number of complete
    spacecraft revolutions in the flight and ``via_peri`` whether the final
    half-revolution passes perihelion. ``sign`` picks one of the two mirror
    solutions about the body orbit plane (circular model only).
    """

    body: str
    half_periods: int = 1
    k_sc: int = 0
    via_peri: bool = True
    sign: int = 1

    @property
    def fixed_tof(self) -> bool:
        return True


Leg = LambertLeg | ResonantLeg | HalfRevLeg


def leg_body_from(leg: Leg) -> str:
    return leg.frm if isinstance(leg, LambertLeg) else leg.body


def leg_body_to(leg: Leg) -> str:
    return leg.to if isinstance(leg, LambertLeg) else leg.body


def fixed_duration_s(system: System, leg: Leg) -> float:
    if isinstance(leg, ResonantLeg):
        return leg.body_revs * system.period_s(leg.body)
    if isinstance(leg, HalfRevLeg):
        return 0.5 * leg.half_periods * system.period_s(leg.body)
    raise TypeError("Lambert legs have no fixed duration")


@dataclass(frozen=True)
class Cycle:
    """A closed chain of legs. Encounter ``i`` is the start of leg ``i``; the
    end of the last leg is encounter 0 one period later."""

    legs: tuple[Leg, ...]
    period_s: float

    def __post_init__(self) -> None:
        n = len(self.legs)
        for i in range(n):
            if leg_body_to(self.legs[i]) != leg_body_from(self.legs[(i + 1) % n]):
                raise ValueError(
                    f"leg {i} ends at {leg_body_to(self.legs[i])} but leg {i + 1} starts elsewhere"
                )
        if not any(isinstance(leg, LambertLeg) for leg in self.legs):
            raise ValueError("a cycle needs at least one Lambert leg")

    @property
    def lambert_index(self) -> list[int]:
        return [i for i, leg in enumerate(self.legs) if isinstance(leg, LambertLeg)]


# ---------------------------------------------------------------------------
# Fixed-leg V-infinity geometry
# ---------------------------------------------------------------------------


def _perp_basis(w: Vec) -> tuple[Vec, Vec, Vec]:
    """Unit ``w_hat`` and two unit vectors completing a right-handed basis.

    ``e1`` is chosen in the plane of ``w`` and +z where possible, so for a
    coplanar body velocity ``e1 = z_hat`` (out of plane) and ``e2`` in plane.
    """
    w_hat = w / np.linalg.norm(w)
    ref = np.array([0.0, 0.0, 1.0])
    if abs(float(w_hat @ ref)) > 0.9:
        ref = np.array([1.0, 0.0, 0.0])
    e1 = ref - (ref @ w_hat) * w_hat
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(w_hat, e1)
    return w_hat, e1, e2


def resonant_circle(
    system: System, leg: ResonantLeg, t_s: float, vinf: float
) -> tuple[float, float, Vec, Vec, Vec] | None:
    """Circle of feasible full-revolution V-infinity vectors at epoch ``t_s``.

    Returns ``(z, rho, w_hat, e1, e2)`` so that ``u(phi) = z w_hat + rho (cos phi
    e1 + sin phi e2)``, or ``None`` when the circle does not exist at this
    ``|v_inf|`` (``|z| > vinf``) or the resonance is unbound.
    """
    r, w = system.state(leg.body, t_s)
    rn = float(np.linalg.norm(r))
    a_body = (system.mu * (system.period_s(leg.body) / (2.0 * math.pi)) ** 2) ** (1.0 / 3.0)
    a_sc = a_body * (leg.body_revs / leg.sc_revs) ** (2.0 / 3.0)
    vf2 = system.mu * (2.0 / rn - 1.0 / a_sc)
    if vf2 <= 0.0:
        return None
    wn = float(np.linalg.norm(w))
    z = (vf2 - vinf * vinf - wn * wn) / (2.0 * wn)
    if abs(z) > vinf:
        return None
    w_hat, e1, e2 = _perp_basis(w)
    return z, math.sqrt(max(vinf * vinf - z * z, 0.0)), w_hat, e1, e2


def resonant_vec(circle: tuple[float, float, Vec, Vec, Vec], phi: float) -> Vec:
    z, rho, w_hat, e1, e2 = circle
    return z * w_hat + rho * (math.cos(phi) * e1 + math.sin(phi) * e2)


def half_rev_tof_s(mu: float, r: float, e: float, k_sc: int, via_peri: bool) -> float:
    """Flight time of a half-revolution conic with ``p = r`` (true anomaly -90 -> +90
    through perihelion, or +90 -> 270 through aphelion), plus ``k_sc`` periods."""
    if not 0.0 <= e < 1.0:
        raise ValueError("half-rev conic must be elliptic")
    a = r / (1.0 - e * e)
    n = math.sqrt(mu / a**3)
    # eccentric anomaly at true anomaly 90 deg: cos E = e
    big_e = math.acos(e)
    m90 = big_e - e * math.sin(big_e)
    t_peri = 2.0 * m90 / n
    period = 2.0 * math.pi / n
    return (t_peri if via_peri else period - t_peri) + k_sc * period


def half_rev_vectors(system: System, leg: HalfRevLeg, t_s: float, vinf: float) -> list[Vec]:
    """Departure V-infinity vectors of a half-revolution return (circular body only).

    With ``p = r`` the spacecraft speed has transverse part ``sqrt(mu/r)`` (equal
    to the circular body speed ``W``) and radial part ``+/- W e``. Tilting the
    conic's plane by ``alpha`` about the radius gives
    ``|v_inf|^2 = W^2 (2 - 2 cos alpha + e^2)``; ``e`` is fixed by the flight
    time, so ``alpha`` follows from ``|v_inf|``. The radial sign is ``-`` when the
    first half-revolution heads to perihelion. Returns the (up to two) vectors
    ``+alpha`` and ``-alpha`` filtered by ``leg.sign`` (``0`` returns both).
    The arrival V-infinity is the mirror image (radial part reversed) at the
    antipode; see :func:`half_rev_arrival`.
    """
    r, w = system.state(leg.body, t_s)
    rn = float(np.linalg.norm(r))
    big_w = float(np.linalg.norm(w))
    target = 0.5 * leg.half_periods * system.period_s(leg.body)
    e = _solve_half_rev_e(system.mu, rn, target, leg.k_sc, leg.via_peri)
    if e is None:
        return []
    cos_alpha = 1.0 - 0.5 * ((vinf / big_w) ** 2 - e * e)
    if abs(cos_alpha) > 1.0:
        return []
    alpha = math.acos(cos_alpha)
    r_hat = r / rn
    w_hat = w / big_w
    h_hat = np.cross(r_hat, w_hat)
    radial = (-1.0 if leg.via_peri else 1.0) * big_w * e
    out = []
    for s in (1, -1):
        if leg.sign not in (0, s):
            continue
        t_hat = math.cos(alpha) * w_hat + s * math.sin(alpha) * h_hat
        v_sc = big_w * t_hat + radial * r_hat
        out.append(v_sc - w)
    return out


def half_rev_arrival(system: System, leg: HalfRevLeg, t_dep: float, u_dep: Vec) -> Vec:
    """Arrival V-infinity of a half-revolution return that departed with ``u_dep``.

    Propagates the conic by the leg duration (an independent two-body step) and
    subtracts the body velocity at arrival.
    """
    r0, w0 = system.state(leg.body, t_dep)
    dt = fixed_duration_s(system, leg)
    _, v1 = kepler_step(r0, w0 + u_dep, dt, system.mu)
    _, w1 = system.state(leg.body, t_dep + dt)
    return np.asarray(v1 - w1, dtype=np.float64)


def _solve_half_rev_e(
    mu: float, r: float, target: float, k_sc: int, via_peri: bool
) -> float | None:
    def f(e: float) -> float:
        return half_rev_tof_s(mu, r, e, k_sc, via_peri) - target

    lo, hi = 0.0, 0.999
    grid = np.linspace(lo, hi, 400)
    vals = [f(float(g)) for g in grid]
    for i in range(len(grid) - 1):
        if vals[i] == 0.0:
            return float(grid[i])
        if vals[i] * vals[i + 1] < 0.0:
            a, b = float(grid[i]), float(grid[i + 1])
            for _ in range(80):
                m = 0.5 * (a + b)
                if f(a) * f(m) <= 0.0:
                    b = m
                else:
                    a = m
            return 0.5 * (a + b)
    return None


# ---------------------------------------------------------------------------
# Date-residual corrector
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class LegEval:
    """One evaluated Lambert leg: epochs (s) and V-infinity vectors (km/s)."""

    t_dep: float
    t_arr: float
    vinf_dep: Vec
    vinf_arr: Vec


def lambert_epochs(system: System, cycle: Cycle, x: Dates) -> list[tuple[float, float]]:
    """``(t_dep, t_arr)`` of each Lambert leg from its start dates ``x``.

    ``x[j]`` is the start of the ``j``-th Lambert leg. The arrival of Lambert leg
    ``j`` is the next Lambert start minus the fixed legs between them (plus the
    period when the chain wraps).
    """
    li = cycle.lambert_index
    n = len(li)
    out = []
    for j in range(n):
        i0, i1 = li[j], li[(j + 1) % n]
        between = 0.0
        k = (i0 + 1) % len(cycle.legs)
        while k != i1:
            between += fixed_duration_s(system, cycle.legs[k])
            k = (k + 1) % len(cycle.legs)
        t_next = x[(j + 1) % n] + (cycle.period_s if j == n - 1 else 0.0)
        out.append((float(x[j]), float(t_next - between)))
    return out


def eval_lambert_legs(system: System, cycle: Cycle, x: Dates) -> list[LegEval] | None:
    """Solve every Lambert leg; ``None`` if any is infeasible."""
    res = []
    for leg_i, (t0, t1) in zip(cycle.lambert_index, lambert_epochs(system, cycle, x), strict=True):
        leg = cycle.legs[leg_i]
        assert isinstance(leg, LambertLeg)
        if t1 <= t0:
            return None
        r1, w1 = system.state(leg.frm, t0)
        r2, w2 = system.state(leg.to, t1)
        try:
            sols = lambert(r1, r2, t1 - t0, mu=system.mu, max_revs=leg.nrev)
        except (LambertError, ValueError):
            return None
        sol = next((s for s in sols if s.n_revs == leg.nrev and s.branch == leg.branch), None)
        if sol is None:
            return None
        res.append(LegEval(t0, t1, sol.v1 - w1, sol.v2 - w2))
    return res


def date_residual(system: System, cycle: Cycle, x: Dates) -> Vec | None:
    """H&M residual: at each junction between consecutive Lambert legs, the
    arrival-minus-departure V-infinity magnitude (km/s); at a massless body the
    full vector difference (3 components)."""
    legs = eval_lambert_legs(system, cycle, x)
    if legs is None:
        return None
    n = len(legs)
    out: list[float] = []
    for j in range(n):
        arr = legs[j].vinf_arr
        dep = legs[(j + 1) % n].vinf_dep
        body = leg_body_to(cycle.legs[cycle.lambert_index[j]])
        if system.body(body).massless:
            if j == n - 1:
                raise ValueError("the wrap junction must be at a massive body")
            out.extend((arr - dep).tolist())
        else:
            out.append(float(np.linalg.norm(arr) - np.linalg.norm(dep)))
    return np.asarray(out)


@dataclass(frozen=True)
class Solution:
    """A converged date vector and its residual."""

    x: Vec
    residual: Vec
    converged: bool

    @property
    def max_abs_residual_kms(self) -> float:
        return float(np.max(np.abs(self.residual))) if self.residual.size else 0.0


def correct_dates(
    system: System,
    cycle: Cycle,
    x0: Dates,
    *,
    tol_kms: float = 1e-9,
    day_scale: float = SECONDS_PER_DAY,
    max_nfev: int = 400,
) -> Solution:
    """Newton/least-squares on the start dates (H&M's Newton-Raphson stage).

    Infeasible Lambert evaluations return a large residual so the solver backs
    off. ``converged`` requires every residual component below ``tol_kms``.
    """
    x0a = np.asarray(x0, dtype=np.float64)
    r0 = date_residual(system, cycle, x0a)
    if r0 is None:
        return Solution(x0a, np.array([np.inf]), False)
    n_res = len(r0)
    if n_res == 0:
        return Solution(x0a, np.array([np.inf]), False)

    def fun(y: Vec) -> Vec:
        r = date_residual(system, cycle, y * day_scale)
        if r is None:
            return np.full(n_res, 1e3)
        return r

    sol = least_squares(
        fun,
        x0a / day_scale,
        method="lm",
        xtol=1e-14,
        ftol=1e-14,
        gtol=1e-14,
        max_nfev=max_nfev * len(x0a),
    )
    x = sol.x * day_scale
    r = date_residual(system, cycle, x)
    if r is None:
        return Solution(x, np.array([np.inf]), False)
    return Solution(x, r, bool(np.max(np.abs(r)) < tol_kms))


# ---------------------------------------------------------------------------
# Flybys: free directions, minimax and the turn gate
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Flyby:
    """One flyby of a solved cycle (inertial-axis V-infinity vectors, km/s)."""

    body: str
    t_s: float
    vinf_in: Vec
    vinf_out: Vec

    @property
    def turn_deg(self) -> float:
        return math.degrees(bend_angle(self.vinf_in, self.vinf_out))

    @property
    def vinf_kms(self) -> float:
        return 0.5 * float(np.linalg.norm(self.vinf_in) + np.linalg.norm(self.vinf_out))


def rmin_radii(body: FlybyBody, vinf_kms: float, turn_rad: float) -> float:
    """Periapsis distance from the body centre, in body radii, for an unpowered
    turn ``turn_rad`` at ``vinf_kms``: ``r_p = mu/v^2 (1/sin(turn/2) - 1)``."""
    if turn_rad <= 0.0:
        return math.inf
    rp = body.mu_km3_s2 / vinf_kms**2 * (1.0 / math.sin(0.5 * turn_rad) - 1.0)
    return rp / body.radius_km


@dataclass(frozen=True)
class _Block:
    """Fixed legs between Lambert leg ``j`` (arriving) and ``j+1`` (departing)."""

    body: str
    t_arr: float
    v_in: Vec
    v_out: Vec
    fixed: tuple[tuple[Leg, float], ...]  # (leg, departure epoch)


def _blocks(system: System, cycle: Cycle, legs: list[LegEval]) -> list[_Block]:
    out = []
    li = cycle.lambert_index
    n = len(li)
    for j in range(n):
        i0, i1 = li[j], li[(j + 1) % n]
        t = legs[j].t_arr
        fixed = []
        k = (i0 + 1) % len(cycle.legs)
        while k != i1:
            fixed.append((cycle.legs[k], t))
            t += fixed_duration_s(system, cycle.legs[k])
            k = (k + 1) % len(cycle.legs)
        v_out = legs[(j + 1) % n].vinf_dep
        if j == n - 1:
            # the departure that opens the next cycle, in the arrival's axes
            v_out = system.wrap_rotation(leg_body_to(cycle.legs[i0]), cycle.period_s) @ v_out
        out.append(
            _Block(
                leg_body_to(cycle.legs[i0]),
                legs[j].t_arr,
                legs[j].vinf_arr,
                v_out,
                tuple(fixed),
            )
        )
    return out


def _block_vectors(
    system: System, blk: _Block, phis: Sequence[float], halfrev_choice: Sequence[int]
) -> list[Vec] | None:
    """V-infinity sequence ``[v_in, u_1dep, u_1arr, ..., v_out]`` through a block."""
    vinf = float(np.linalg.norm(blk.v_in))
    seq: list[Vec] = [blk.v_in]
    ip = ih = 0
    for leg, t in blk.fixed:
        if isinstance(leg, ResonantLeg):
            circ = resonant_circle(system, leg, t, vinf)
            if circ is None:
                return None
            u = resonant_vec(circ, phis[ip])
            ip += 1
            seq += [u, u]
        elif isinstance(leg, HalfRevLeg):
            cands = half_rev_vectors(system, leg, t, vinf)
            if halfrev_choice[ih] >= len(cands):
                return None
            u = cands[halfrev_choice[ih]]
            ih += 1
            seq += [u, half_rev_arrival(system, leg, t, u)]
        else:  # pragma: no cover - blocks hold only fixed legs
            raise TypeError(leg)
    seq.append(blk.v_out)
    return seq


def _available(body: FlybyBody, vinf: float) -> float:
    return available_bend_rad(body.mu_km3_s2, body.radius_km, body.alt_floor_km, vinf)


def optimise_block(
    system: System, blk: _Block, *, n_grid: int = 72
) -> tuple[list[Flyby], float] | None:
    """Choose the free directions of a block's fixed legs to minimise the largest
    demanded/available turn ratio over the block's flybys (minimax).

    Returns the flybys and the minimax ratio, or ``None`` if no choice exists
    (a circle or half-rev solution missing at this V-infinity).
    """
    body = system.body(blk.body)
    n_phi = sum(isinstance(leg, ResonantLeg) for leg, _ in blk.fixed)
    n_half = sum(isinstance(leg, HalfRevLeg) for leg, _ in blk.fixed)
    vinf = float(np.linalg.norm(blk.v_in))
    avail = _available(body, vinf) if not body.massless else math.pi

    def epochs() -> list[float]:
        ts = [blk.t_arr]
        for leg, t in blk.fixed:
            ts.append(t + fixed_duration_s(system, leg))
        return ts

    def ratios(seq: list[Vec]) -> list[float]:
        return [bend_angle(seq[2 * i], seq[2 * i + 1]) / avail for i in range(len(seq) // 2)]

    best: tuple[float, list[Vec], list[float], tuple[int, ...]] | None = None
    for hc in np.ndindex(*([2] * n_half)) if n_half else [()]:
        if n_phi == 0:
            seq = _block_vectors(system, blk, [], list(hc))
            if seq is None:
                continue
            val = max(ratios(seq))
            if best is None or val < best[0]:
                best = (val, seq, [], hc)
            continue
        grids = np.linspace(
            0.0, 2.0 * math.pi, n_grid if n_phi == 1 else max(12, 36 // n_phi), endpoint=False
        )
        starts = []
        for combo in np.ndindex(*([len(grids)] * n_phi)):
            phis = [float(grids[c]) for c in combo]
            seq = _block_vectors(system, blk, phis, list(hc))
            if seq is None:
                break
            starts.append((max(ratios(seq)), phis))
        if not starts:
            continue
        starts.sort(key=lambda s: s[0])
        for _, p0 in starts[:4]:
            # epigraph form: minimise s subject to ratio_i(phi) <= s
            def obj(z: Vec) -> float:
                return float(z[-1])

            def cons(z: Vec, hc: tuple[int, ...] = hc) -> Vec:
                seq = _block_vectors(system, blk, list(z[:-1]), list(hc))
                if seq is None:
                    return np.full(len(blk.fixed) + 1, -1.0)
                return np.asarray(z[-1] - np.asarray(ratios(seq)), dtype=np.float64)

            seq0 = _block_vectors(system, blk, p0, list(hc))
            assert seq0 is not None
            z0 = np.array([*p0, max(ratios(seq0))])
            r = minimize(  # type: ignore[call-overload]
                obj,
                z0,
                constraints=[{"type": "ineq", "fun": cons}],
                method="SLSQP",
                options={"ftol": 1e-12, "maxiter": 300},
            )
            seq = _block_vectors(system, blk, list(r.x[:-1]), list(hc))
            if seq is None:
                continue
            val = max(ratios(seq))
            if best is None or val < best[0]:
                best = (val, seq, list(r.x[:-1]), hc)
    if best is None:
        return None
    seq = best[1]
    if n_phi >= 1 and len(seq) // 2 >= 3:
        seq = _balance(system, blk, best[0], best[2], best[3], ratios) or seq
    ts = epochs()
    flybys = [Flyby(blk.body, ts[i], seq[2 * i], seq[2 * i + 1]) for i in range(len(seq) // 2)]
    return flybys, best[0]


def _balance(
    system: System,
    blk: _Block,
    s_star: float,
    phis: list[float],
    hc: tuple[int, ...],
    ratios: Callable[[list[Vec]], list[float]],
) -> list[Vec] | None:
    """Tie-break of the minimax: with the largest ratio held at ``s_star``, raise the
    smallest (spread the turn as evenly as the circles allow).

    The minimax leaves the inner flybys of a block with two or more free
    directions undetermined; this fixes them. It cannot change the largest
    ratio, so it never changes a gate verdict.
    """
    cap = s_star * (1.0 + 1e-9) + 1e-12

    def obj(z: Vec) -> float:
        return -float(z[-1])

    def cons(z: Vec) -> Vec:
        seq = _block_vectors(system, blk, list(z[:-1]), list(hc))
        if seq is None:
            return np.full(2 * (len(blk.fixed) + 1), -1.0)
        rr = np.asarray(ratios(seq))
        return np.asarray(np.concatenate([rr - z[-1], cap - rr]), dtype=np.float64)

    seq0 = _block_vectors(system, blk, phis, list(hc))
    if seq0 is None:
        return None
    z0 = np.array([*phis, min(ratios(seq0))])
    r = minimize(  # type: ignore[call-overload]
        obj,
        z0,
        constraints=[{"type": "ineq", "fun": cons}],
        method="SLSQP",
        options={"ftol": 1e-12, "maxiter": 300},
    )
    seq = _block_vectors(system, blk, list(r.x[:-1]), list(hc))
    if seq is None or max(ratios(seq)) > cap * (1.0 + 1e-6):
        return None
    return seq


@dataclass(frozen=True)
class CycleReport:
    """Every flyby of a solved cycle with the demanded-turn gate applied."""

    flybys: tuple[Flyby, ...]
    gate: TurnGateReport
    near_180: bool

    @property
    def status(self) -> str:
        """``fail`` on any near-180-degree demand (owner ruling), else the gate's."""
        return "fail" if self.near_180 else self.gate.status


def cycle_flybys(system: System, cycle: Cycle, x: Dates, *, n_grid: int = 72) -> list[Flyby] | None:
    """All flybys of a solved cycle, with minimax-chosen free directions."""
    legs = eval_lambert_legs(system, cycle, x)
    if legs is None:
        return None
    out: list[Flyby] = []
    for blk in _blocks(system, cycle, legs):
        if system.body(blk.body).massless and not blk.fixed:
            out.append(Flyby(blk.body, blk.t_arr, blk.v_in, blk.v_out))
            continue
        res = optimise_block(system, blk, n_grid=n_grid)
        if res is None:
            return None
        out.extend(res[0])
    return out


def gate_cycle(
    system: System, flybys: Sequence[Flyby], *, alt_floor_km: float | None = None
) -> CycleReport:
    """#888/#937 demanded-turn gate on every massive flyby (inertial axes)."""
    encs = []
    near = False
    for i, f in enumerate(flybys):
        b = system.body(f.body)
        if b.massless:
            continue
        if f.turn_deg >= 180.0 - NEAR_180_DEG:
            near = True
        encs.append(
            Encounter(
                body=f.body,
                mu_km3_s2=b.mu_km3_s2,
                radius_km=b.radius_km,
                alt_floor_km=b.alt_floor_km if alt_floor_km is None else alt_floor_km,
                vinf_in=np.asarray(f.vinf_in),
                vinf_out=np.asarray(f.vinf_out),
                label=f"f{i}@{f.t_s / SECONDS_PER_DAY:.2f}d",
                tidal_speed_kms=tidal_speed_kms(f.body),
            )
        )
    return CycleReport(tuple(flybys), demanded_turn_gate(encs), near)


def encounter_self_consistency(system: System, cycle: Cycle, x: Dates) -> float:
    """Largest miss distance (km) when each Lambert leg's departure state is
    propagated independently (two-body universal variables, not the Lambert
    solver) to its arrival epoch and compared with the arrival body's position.

    Patched conics put each encounter at the body centre, so this must be far
    inside the body's sphere of influence (#480 constructed-tour rule).
    """
    legs = eval_lambert_legs(system, cycle, x)
    if legs is None:
        return math.inf
    worst = 0.0
    for leg_i, ev in zip(cycle.lambert_index, legs, strict=True):
        leg = cycle.legs[leg_i]
        assert isinstance(leg, LambertLeg)
        r1, w1 = system.state(leg.frm, ev.t_dep)
        r_end, _ = kepler_step(r1, w1 + ev.vinf_dep, ev.t_arr - ev.t_dep, system.mu)
        r2, _ = system.state(leg.to, ev.t_arr)
        worst = max(worst, float(np.linalg.norm(r_end - r2)))
    return worst


def sphere_of_influence_km(system: System, code: str) -> float:
    """Laplace sphere of influence ``a (m/M)^(2/5)`` of an encounter body."""
    if code in PLANETS:
        a = PLANETS[code].sma_au * AU_KM
    elif code in SATELLITES:
        a = SATELLITES[code].sma_km
    else:
        raise KeyError(code)
    return float(a * (system.body(code).mu_km3_s2 / system.mu) ** 0.4)
