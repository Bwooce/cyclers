"""Bridge from a published Uranus tour end-state into a catalogued quasi-cycler (#885).

Planar (circular-coplanar) patched-conic primitives about Uranus, used by
``scripts/screen_885_uranus_tour_to_cycler_bridge.py``:

* moon states on circular coplanar orbits, with longitudes optionally anchored
  to URA111 at an epoch (projected onto the Uranian equatorial plane);
* planar Kepler propagation, next-apoapsis timing, arc minimum radius;
* the flyby (rotation of V-infinity by at most the project's ``max_bend`` at a
  stated altitude);
* :func:`landau_leg`, the tour step of Landau, Davis & Karimi (AAS 23-460,
  "Moon Tour Trade Space", steps 3-4 and 7): propagate from the flyby to the
  next apoapsis, burn impulsively there, Lambert (at most one revolution) to the
  next moon; and
* :func:`tisserand_ladder`, the ideal-model, phase-free flyby-count search on
  the Tisserand graph (a flyby keeps V-infinity at its own moon and rotates the
  pump angle by at most the maximum bend).

The Campagnola-Russell phase-free VILM layer (:mod:`cyclerfinder.search.vilm`)
is deliberately NOT used: its Eq. 25 Gamma has a pole near adimensional
V-infinity 0.41 and a zero near 1.38 (exterior), and every V-infinity of this
problem lies in or across that band. Only Kepler arithmetic and Lambert arcs
are used here.

Frames: a planar inertial frame in the moons' common orbit plane, x along an
arbitrary fixed direction, prograde moon motion counter-clockwise. Units: km,
km/s, seconds.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass, field

import numpy as np

from cyclerfinder.core.flyby import max_bend
from cyclerfinder.core.lambert import LambertConvergenceError, LambertGeometryError, lambert
from cyclerfinder.core.satellites import PRIMARIES, SATELLITES

MU_URANUS: float = PRIMARIES["Uranus"]
TOUR_MOONS: tuple[str, ...] = ("Ariel", "Umbriel", "Titania", "Oberon")
#: Landau, Davis & Karimi 2023 (AAS 23-460) Table 2, "min. range to Uranus".
RP_FLOOR_KM: float = 103000.0
DAY_S: float = 86400.0


def moon_sma(moon: str) -> float:
    return SATELLITES[moon].sma_km


def moon_speed(moon: str) -> float:
    """Circular speed of ``moon`` about Uranus, km/s."""
    return math.sqrt(MU_URANUS / moon_sma(moon))


def moon_mean_motion(moon: str) -> float:
    """Mean motion of ``moon`` about Uranus, rad/s (Kepler III, registry sma)."""
    return math.sqrt(MU_URANUS / moon_sma(moon) ** 3)


def moon_period_days(moon: str) -> float:
    return 2.0 * math.pi / moon_mean_motion(moon) / DAY_S


@dataclass(frozen=True)
class MoonPhases:
    """Moon longitudes (rad) at t = 0; circular coplanar motion afterwards."""

    theta0: dict[str, float]

    def state(self, moon: str, t_s: float) -> tuple[np.ndarray, np.ndarray]:
        th = self.theta0[moon] + moon_mean_motion(moon) * t_s
        a = moon_sma(moon)
        vc = moon_speed(moon)
        c, s = math.cos(th), math.sin(th)
        return np.array([a * c, a * s]), np.array([-vc * s, vc * c])

    def longitude(self, moon: str, t_s: float) -> float:
        return self.theta0[moon] + moon_mean_motion(moon) * t_s


def phases_from_ura111(epoch_utc: str, kernel_paths: Sequence[str]) -> MoonPhases:
    """Moon longitudes at ``epoch_utc`` from URA111, in the Uranian equatorial plane.

    The plane normal is Ariel's orbital angular momentum at the epoch (the regular
    moons are all within ~0.4 deg of the equator). Longitudes are measured from
    the projection of Ariel's position, so Ariel sits at 0 rad.
    """
    import spiceypy as spice

    spice.kclear()
    try:
        for p in kernel_paths:
            spice.furnsh(p)
        et = float(spice.str2et(epoch_utc))
        states: dict[str, np.ndarray] = {}
        for m in TOUR_MOONS:
            st, _ = spice.spkezr(m.upper(), et, "J2000", "NONE", "URANUS")
            states[m] = np.asarray(st, dtype=np.float64)
    finally:
        spice.kclear()
    ra = states["Ariel"]
    zhat = np.cross(ra[:3], ra[3:])
    zhat /= np.linalg.norm(zhat)
    xhat = ra[:3] - np.dot(ra[:3], zhat) * zhat
    xhat /= np.linalg.norm(xhat)
    yhat = np.cross(zhat, xhat)
    theta0 = {}
    for m, st in states.items():
        theta0[m] = math.atan2(float(np.dot(st[:3], yhat)), float(np.dot(st[:3], xhat)))
    return MoonPhases(theta0=theta0)


# ---------------------------------------------------------------------------
# Planar Kepler
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Elements:
    a: float
    e: float
    omega: float  # longitude of periapsis, rad
    h: float  # specific angular momentum (z), km^2/s
    mean_anomaly: float  # rad, in [0, 2 pi)

    @property
    def rp(self) -> float:
        return self.a * (1.0 - self.e)

    @property
    def ra(self) -> float:
        return self.a * (1.0 + self.e)

    @property
    def n(self) -> float:
        return math.sqrt(MU_URANUS / self.a**3)

    @property
    def period_s(self) -> float:
        return 2.0 * math.pi / self.n


def elements(r: np.ndarray, v: np.ndarray) -> Elements | None:
    """Planar elements of a bound prograde orbit; None if unbound or retrograde."""
    rn = float(np.hypot(r[0], r[1]))
    v2 = float(v[0] ** 2 + v[1] ** 2)
    energy = 0.5 * v2 - MU_URANUS / rn
    h = float(r[0] * v[1] - r[1] * v[0])
    if energy >= 0.0 or h <= 0.0:
        return None
    a = -MU_URANUS / (2.0 * energy)
    rv = float(r[0] * v[0] + r[1] * v[1])
    ev = ((v2 - MU_URANUS / rn) * r - rv * v) / MU_URANUS
    e = float(np.hypot(ev[0], ev[1]))
    if e < 1e-12:
        omega = 0.0
        nu = math.atan2(r[1], r[0])
    else:
        omega = math.atan2(ev[1], ev[0])
        nu = math.atan2(r[1], r[0]) - omega
    ecc_anom = 2.0 * math.atan2(
        math.sqrt(1.0 - e) * math.sin(nu / 2), math.sqrt(1.0 + e) * math.cos(nu / 2)
    )
    m_anom = (ecc_anom - e * math.sin(ecc_anom)) % (2.0 * math.pi)
    return Elements(a=a, e=e, omega=omega, h=h, mean_anomaly=m_anom)


def _solve_kepler(m: float, e: float) -> float:
    ea = m if e < 0.8 else math.pi
    for _ in range(60):
        f = ea - e * math.sin(ea) - m
        d = f / (1.0 - e * math.cos(ea))
        ea -= d
        if abs(d) < 1e-14:
            break
    return ea


def state_at_mean_anomaly(el: Elements, m_anom: float) -> tuple[np.ndarray, np.ndarray]:
    ea = _solve_kepler(m_anom % (2.0 * math.pi), el.e)
    nu = 2.0 * math.atan2(
        math.sqrt(1.0 + el.e) * math.sin(ea / 2), math.sqrt(1.0 - el.e) * math.cos(ea / 2)
    )
    p = el.a * (1.0 - el.e**2)
    rn = p / (1.0 + el.e * math.cos(nu))
    th = el.omega + nu
    k = math.sqrt(MU_URANUS / p)
    vr = k * el.e * math.sin(nu)
    vt = k * (1.0 + el.e * math.cos(nu))
    c, s = math.cos(th), math.sin(th)
    return np.array([rn * c, rn * s]), np.array([vr * c - vt * s, vr * s + vt * c])


def propagate(r: np.ndarray, v: np.ndarray, dt_s: float) -> tuple[np.ndarray, np.ndarray]:
    el = elements(r, v)
    if el is None:
        raise ValueError("propagate: orbit is unbound or retrograde")
    return state_at_mean_anomaly(el, el.mean_anomaly + el.n * dt_s)


def time_to_next_apoapsis(el: Elements) -> float:
    return ((math.pi - el.mean_anomaly) % (2.0 * math.pi)) / el.n


def arc_min_radius(r1: np.ndarray, v1: np.ndarray, tof_s: float) -> float:
    """Minimum radius along the conic from (r1, v1) over ``tof_s`` seconds."""
    el = elements(r1, v1)
    if el is None:
        return 0.0
    t_peri = ((2.0 * math.pi - el.mean_anomaly) % (2.0 * math.pi)) / el.n
    r1n = float(np.hypot(*r1))
    if t_peri <= tof_s:
        return el.rp
    r2, _ = propagate(r1, v1, tof_s)
    return min(r1n, float(np.hypot(*r2)))


# ---------------------------------------------------------------------------
# Flyby
# ---------------------------------------------------------------------------


def max_bend_rad(moon: str, vinf_kms: float, alt_km: float) -> float:
    """Maximum ballistic turn at ``moon`` (project gate :func:`core.flyby.max_bend`)."""
    sat = SATELLITES[moon]
    return max_bend(sat.mu_km3_s2, sat.radius_eq_km + alt_km, vinf_kms)


def rotate(vec: np.ndarray, ang: float) -> np.ndarray:
    c, s = math.cos(ang), math.sin(ang)
    return np.array([c * vec[0] - s * vec[1], s * vec[0] + c * vec[1]])


def signed_angle(u: np.ndarray, w: np.ndarray) -> float:
    return math.atan2(float(u[0] * w[1] - u[1] * w[0]), float(u[0] * w[0] + u[1] * w[1]))


def entry_correction_kms(vinf_in: np.ndarray, vinf_need: np.ndarray, bend_max: float) -> float:
    """Smallest impulse that turns ``vinf_in`` into ``vinf_need`` with one flyby.

    The flyby rotates ``vinf_in`` by at most ``bend_max`` keeping its magnitude;
    the impulse (applied at the sphere of influence, as McAdams et al. 2011 do)
    covers what remains. Returns the impulse magnitude, km/s.
    """
    ang = signed_angle(vinf_in, vinf_need)
    turn = max(-bend_max, min(bend_max, ang))
    return float(np.linalg.norm(vinf_need - rotate(vinf_in, turn)))


# ---------------------------------------------------------------------------
# The Landau-Davis-Karimi tour step
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Leg:
    dep_moon: str
    arr_moon: str
    t_dep: float
    t_burn: float
    t_arr: float
    r_dep: np.ndarray
    v_dep: np.ndarray  # spacecraft velocity just after the departure flyby
    r_burn: np.ndarray
    v_burn_before: np.ndarray
    v_burn_after: np.ndarray
    r_arr: np.ndarray
    v_arr: np.ndarray  # spacecraft velocity at arrival (before the arrival flyby)
    dv_kms: float
    vinf_in: np.ndarray  # arrival V-infinity vector at arr_moon
    min_radius_km: float
    n_revs: int
    branch: str = "single"
    extra: dict[str, float] = field(default_factory=dict)


def landau_leg(
    *,
    dep_moon: str,
    t_dep: float,
    vinf_out: np.ndarray,
    arr_moon: str,
    t_arr: float,
    phases: MoonPhases,
    max_revs: int = 1,
    rp_floor_km: float = RP_FLOOR_KM,
) -> list[Leg]:
    """All feasible Landau tour steps for one (departure, arrival time) pair.

    Departure state: ``dep_moon`` at ``t_dep`` with V-infinity ``vinf_out`` (after
    the flyby). Propagate to the next apoapsis, burn there, Lambert to
    ``arr_moon`` at ``t_arr`` with at most ``max_revs`` revolutions. Every branch
    whose pre-burn and post-burn arcs stay above ``rp_floor_km`` is returned,
    cheapest first. An empty list means infeasible.
    """
    r0, vm0 = phases.state(dep_moon, t_dep)
    v0 = vm0 + vinf_out
    el = elements(r0, v0)
    if el is None:
        return []
    dt_apo = time_to_next_apoapsis(el)
    if el.rp < rp_floor_km and arc_min_radius(r0, v0, dt_apo) < rp_floor_km:
        return []
    t_burn = t_dep + dt_apo
    tof = t_arr - t_burn
    if tof <= 3600.0:
        return []
    rb, vb = state_at_mean_anomaly(el, math.pi)
    r1, vm1 = phases.state(arr_moon, t_arr)
    try:
        sols = lambert(
            np.array([rb[0], rb[1], 0.0]),
            np.array([r1[0], r1[1], 0.0]),
            tof,
            mu=MU_URANUS,
            prograde=True,
            max_revs=max_revs,
        )
    except (LambertGeometryError, LambertConvergenceError, ValueError):
        return []
    out: list[Leg] = []
    for s in sols:
        va = np.asarray(s.v1[:2], dtype=np.float64)
        v2 = np.asarray(s.v2[:2], dtype=np.float64)
        rmin = arc_min_radius(rb, va, tof)
        if rmin < rp_floor_km:
            continue
        out.append(
            Leg(
                dep_moon=dep_moon,
                arr_moon=arr_moon,
                t_dep=t_dep,
                t_burn=t_burn,
                t_arr=t_arr,
                r_dep=r0,
                v_dep=v0,
                r_burn=rb,
                v_burn_before=vb,
                v_burn_after=va,
                r_arr=r1,
                v_arr=v2,
                dv_kms=float(np.linalg.norm(va - vb)),
                vinf_in=v2 - vm1,
                min_radius_km=rmin,
                n_revs=int(s.n_revs),
                branch=str(s.branch),
            )
        )
    out.sort(key=lambda leg: leg.dv_kms)
    return out


# ---------------------------------------------------------------------------
# Ideal model: Tisserand points and the phase-free flyby ladder
# ---------------------------------------------------------------------------


def orbit_from_vinf_pair(m1: str, v1: float, m2: str, v2: float) -> tuple[float, float]:
    """(a, e) of the coplanar orbit with V-infinity ``v1`` at ``m1`` and ``v2`` at ``m2``.

    Solves the two Tisserand relations T = r_M / a + 2 sqrt(p / r_M) for 1/a and
    sqrt(p).
    """
    r1, r2 = moon_sma(m1), moon_sma(m2)
    t1 = 3.0 - (v1 / moon_speed(m1)) ** 2
    t2 = 3.0 - (v2 / moon_speed(m2)) ** 2
    mat = np.array([[r1, 2.0 / math.sqrt(r1)], [r2, 2.0 / math.sqrt(r2)]])
    inv_a, sqrt_p = np.linalg.solve(mat, [t1, t2])
    a = 1.0 / float(inv_a)
    p = float(sqrt_p) ** 2
    return a, math.sqrt(max(0.0, 1.0 - p / a))


def vinf_at(moon: str, a: float, e: float) -> float:
    """V-infinity at ``moon`` for the coplanar orbit (a, e); nan if it does not cross."""
    r = moon_sma(moon)
    if not (a * (1.0 - e) <= r <= a * (1.0 + e)):
        return float("nan")
    t = r / a + 2.0 * math.sqrt(a * (1.0 - e * e) / r)
    return moon_speed(moon) * math.sqrt(max(0.0, 3.0 - t))


def arrival_pump(moon: str, a: float, e: float) -> tuple[float, float] | None:
    """(V-infinity, pump angle in [0, pi]) on arrival at ``moon`` for orbit (a, e)."""
    r = moon_sma(moon)
    if not (a * (1.0 - e) <= r <= a * (1.0 + e)):
        return None
    h = math.sqrt(MU_URANUS * a * (1.0 - e * e))
    vt = h / r
    vr = math.sqrt(max(0.0, MU_URANUS * (2.0 / r - 1.0 / a) - vt * vt))
    vx = vt - moon_speed(moon)
    return math.hypot(vx, vr), math.atan2(vr, vx)


def orbit_from_pump(moon: str, vinf: float, pump: float) -> tuple[float, float] | None:
    """(a, e) of the orbit leaving ``moon`` with V-infinity ``vinf`` at pump angle ``pump``."""
    vt = moon_speed(moon) + vinf * math.cos(pump)
    vr = vinf * math.sin(pump)
    r = moon_sma(moon)
    energy = 0.5 * (vt * vt + vr * vr) - MU_URANUS / r
    if energy >= 0.0 or vt <= 0.0:
        return None
    a = -MU_URANUS / (2.0 * energy)
    p = (r * vt) ** 2 / MU_URANUS
    return a, math.sqrt(max(0.0, 1.0 - p / a))


@dataclass(frozen=True)
class LadderResult:
    n_flybys: int | None  # flybys INCLUDING the final one onto the target orbit
    path: tuple[tuple[str, float, float], ...]  # (moon, vinf, pump) arrival states
    n_states: int


def tisserand_ladder(
    start: tuple[str, float, float],
    targets: dict[str, tuple[float, float]],
    *,
    moons: Sequence[str] = TOUR_MOONS,
    alt_km: float = 50.0,
    rp_floor_km: float = RP_FLOOR_KM,
    bend_fractions: Sequence[float] = (1.0, 0.5, 0.0, -0.5, -1.0),
    vinf_tol_kms: float = 0.01,
    vinf_bucket_kms: float = 0.002,
    pump_bucket_rad: float = math.radians(0.05),
    max_depth: int = 120,
    max_front: int = 600_000,
) -> LadderResult:
    """Breadth-first minimum flyby count on the coplanar Tisserand graph.

    A node is an ARRIVAL state (moon, V-infinity, pump angle). An edge is one
    ballistic flyby at that moon (pump angle rotated by a sampled fraction of the
    maximum bend at ``alt_km``) followed by the next encounter with any moon the
    new orbit crosses. Orbits with periapsis under ``rp_floor_km`` are discarded.
    Phase-free: it assumes every encounter can be phased, so the count is a
    LOWER BOUND on ballistic flybys at these bend samples, not a trajectory.
    ``targets[moon] = (vinf, pump)``: reached when a flyby at that moon can turn
    the arrival state onto the target (same V-infinity within ``vinf_tol_kms``,
    pump angle within the maximum bend).
    """

    def key(s: tuple[str, float, float]) -> tuple[str, int, int]:
        return (s[0], round(s[1] / vinf_bucket_kms), round(s[2] / pump_bucket_rad))

    parent: dict[tuple[str, int, int], tuple[str, int, int] | None] = {key(start): None}
    states: dict[tuple[str, int, int], tuple[str, float, float]] = {key(start): start}
    front = [start]
    depth = 0
    while front and depth < max_depth:
        nxt: list[tuple[str, float, float]] = []
        for s in front:
            m, v, al = s
            dmax = max_bend_rad(m, v, alt_km)
            if m in targets:
                vt, at = targets[m]
                if abs(v - vt) < vinf_tol_kms and abs(al - at) <= dmax:
                    path = [s]
                    k = parent[key(s)]
                    while k is not None:
                        path.append(states[k])
                        k = parent[k]
                    return LadderResult(depth + 1, tuple(reversed(path)), len(parent))
            for f in bend_fractions:
                al2 = min(max(al + f * dmax, 0.0), math.pi)
                orb = orbit_from_pump(m, v, al2)
                if orb is None:
                    continue
                a, e = orb
                if a * (1.0 - e) < rp_floor_km:
                    continue
                for j in moons:
                    arr = arrival_pump(j, a, e)
                    if arr is None:
                        continue
                    ns = (j, arr[0], arr[1])
                    k2 = key(ns)
                    if k2 in parent:
                        continue
                    parent[k2] = key(s)
                    states[k2] = ns
                    nxt.append(ns)
        depth += 1
        front = nxt
        if len(front) > max_front:
            break
    return LadderResult(None, (), len(parent))


# ---------------------------------------------------------------------------
# Cycler entry geometry and the realisation beam search
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class CyclerGeometry:
    """Two-moon (1,1) symmetric quasi-cycler in the circular-coplanar model.

    ``rel_offset_rad`` is the second moon's longitude minus the first moon's at
    the first moon's departure (the catalogue / #566 convention, ``theta0[moon_b]
    = phase0 + rel_offset``); ``leg_s`` the equal leg time; ``n_revs`` per leg.
    """

    moon_a: str
    moon_b: str
    leg_s: float
    rel_offset_rad: float
    n_revs: int = 0

    def departure_times(self, phases: MoonPhases, t_lo: float, t_hi: float) -> list[float]:
        """Epochs in [t_lo, t_hi] at which moon_a departures start a cycle."""
        na, nb = moon_mean_motion(self.moon_a), moon_mean_motion(self.moon_b)
        rel0 = phases.theta0[self.moon_b] - phases.theta0[self.moon_a] - self.rel_offset_rad
        t_syn = 2.0 * math.pi / abs(na - nb)
        t0 = (
            (rel0 % (2.0 * math.pi)) / (na - nb)
            if na > nb
            else (-rel0 % (2.0 * math.pi)) / (nb - na)
        )
        k0 = math.ceil((t_lo - t0) / t_syn)
        out = []
        t = t0 + k0 * t_syn
        while t <= t_hi:
            out.append(t)
            t += t_syn
        return out

    def leg_vinf(
        self, phases: MoonPhases, t_dep: float, dep: str, arr: str
    ) -> tuple[np.ndarray, np.ndarray]:
        """(V-infinity out at ``dep``, V-infinity in at ``arr``) of one cycler leg."""
        r0, v0 = phases.state(dep, t_dep)
        r1, v1 = phases.state(arr, t_dep + self.leg_s)
        sols = lambert(
            np.array([r0[0], r0[1], 0.0]),
            np.array([r1[0], r1[1], 0.0]),
            self.leg_s,
            mu=MU_URANUS,
            prograde=True,
            max_revs=self.n_revs,
        )
        want = [s for s in sols if s.n_revs == self.n_revs]
        best = min(want, key=lambda s: float(np.linalg.norm(np.asarray(s.v1[:2]) - v0)))
        return np.asarray(best.v1[:2]) - v0, np.asarray(best.v2[:2]) - v1

    def entry_targets(
        self, phases: MoonPhases, t_lo: float, t_hi: float
    ) -> list[tuple[str, float, np.ndarray]]:
        """(moon, epoch, V-infinity the spacecraft must LEAVE with) for every entry.

        Entry at moon_a at a cycle start (leave on the a->b leg), or at moon_b one
        leg later (leave on the b->a leg).
        """
        out: list[tuple[str, float, np.ndarray]] = []
        for tc in self.departure_times(phases, t_lo - self.leg_s, t_hi):
            if t_lo <= tc <= t_hi:
                out.append(
                    (self.moon_a, tc, self.leg_vinf(phases, tc, self.moon_a, self.moon_b)[0])
                )
            tb = tc + self.leg_s
            if t_lo <= tb <= t_hi:
                out.append(
                    (self.moon_b, tb, self.leg_vinf(phases, tb, self.moon_b, self.moon_a)[0])
                )
        return out
