"""Jovian-moon n-body lane — VALIDATION INFRASTRUCTURE ONLY (#223).

Jupiter-centered analogue of the heliocentric restricted n-body harness
(:mod:`cyclerfinder.nbody.propagator` / :mod:`cyclerfinder.nbody.forces`),
built for the Liang et al. 2024 Member D same-family re-propagation lane
(`liang-2024-cgcec-ephemeris-2033`). Nothing here is consumed by the
production construct/score/verify pipeline (seeds-not-tracks intact, mirror
``nbody/__init__.py``).

Model (recorded precisely — the lane's verdict is conditional on it):

* **Ephemeris** — Galilean moon states from the JPL NAIF JUP365 SPK kernel
  (generic kernels; the SAME kernel Liang et al. cite for Member D), read via
  spiceypy as Jupiter-centered (observer ``JUPITER``, 599) J2000-equatorial
  km / km/s. The time axis is TDB seconds since J2000 — identical to SPICE ET
  and to the harness axis (``nbody/convert.py``), so ``t_sec`` is passed to
  ``spkezr`` unmodified. The kernel path is NOT committed; callers pass it
  (default from the ``CYCLERFINDER_JUP365`` environment variable).
* **Frame** — J2000 equatorial throughout (Lambert's prograde convention is
  satisfied: the Galilean orbit normals have +z components ~0.90 in J2000).
  The paper works in ecliptic J2000; both are inertial, and every quantity
  compared (ToF, V∞ magnitude, altitude, defect Δv) is frame-invariant.
* **GM values** — registry conventions (``core/satellites.py``): moon GMs
  from :data:`cyclerfinder.core.satellites.SATELLITES` (JPL SSD), Jupiter
  central GM from :data:`cyclerfinder.core.satellites.PRIMARIES`
  (1.26686534e8 km^3/s^2). NOTE: that registry value matches the JUP365
  comment-area *Jupiter-alone* GM (1.266865319003704e8, fractional diff
  1.6e-8) — the correct central mass when the moons enter as separate point
  masses (no double count); the JUP365 *system* GM is 1.267127618414429e8.
  The JUP365 comment-area moon GMs differ from the registry SSD values by
  <= 8e-6 fractional (Europa worst: 3202.712 vs 3202.739) — immaterial at
  this lane's tolerances; the registry values are used per convention.
* **n-body dynamics** — REBOUND / IAS15: massless spacecraft about a central
  Jupiter point mass, with Io+Europa+Ganymede+Callisto as point masses ON
  JUP365 RAILS (cubic-spline cache, 0.02-day grid; spline error ~1e-2 km,
  pinned by the smoke test) including the Jupiter-frame indirect term.
  Deliberately EXCLUDED: Jupiter J2, solar tide, smaller moons — the paper's
  own model has none of them (Kepler legs + impulsive flybys), so including
  them would measure a different gap than the patched-conic -> continuous-
  moon-gravity gap this lane exists to measure. The rails force is softened
  at the MOON SURFACE (``radius_eq_km``), not at the safe altitude — flybys
  down to 100 km must be integrated for real; a below-surface dive is a
  divergent-seed signal, not a modelling device.

Construction route (mirrors Liang's own ephemeris-model method, Sec. III):
patched-conic chain on real JUP365 geometry (multi-rev Lambert legs between
moon encounter positions, per-cycle local optimization of the flyby epochs to
minimize the powered-flyby defect Δv, chained cycle-by-cycle with the inbound
V∞ inherited) -> periapsis-node conversion -> REBOUND re-propagation +
per-cycle multiple-shooting correction. See the results note
``docs/notes/2026-06-13-liang-member-d-nbody.md`` for the ceiling assessment
(same-family V1-class evidence, NOT V3) and the run record.
"""

from __future__ import annotations

import math
import os
from collections.abc import Sequence
from dataclasses import dataclass
from typing import TYPE_CHECKING, cast

import numpy as np
from numpy.typing import NDArray
from scipy.interpolate import CubicSpline

from cyclerfinder.core.constants import SECONDS_PER_DAY
from cyclerfinder.core.lambert import lambert
from cyclerfinder.core.satellites import PRIMARIES, SATELLITES

if TYPE_CHECKING:
    from cyclerfinder.nbody.shooter import ShootingSeed, ShootResult

Vec3 = NDArray[np.float64]

MU_JUPITER_KM3_S2: float = PRIMARIES["Jupiter"]
"""Central Jupiter GM (km^3/s^2) — registry value; equals the JUP365
comment-area Jupiter-alone GM to 1.6e-8 fractional (module docstring)."""

GALILEAN: tuple[str, ...] = ("Io", "Europa", "Ganymede", "Callisto")

JUP365_ENV_VAR = "CYCLERFINDER_JUP365"

#: CGCEC flyby sequence per ~100 d cycle (Liang et al. 2024, Fig. 4).
CGCEC: tuple[str, ...] = ("Callisto", "Ganymede", "Callisto", "Europa", "Callisto")


def jup365_kernel_path() -> str | None:
    """Resolve the JUP365 kernel path from the environment (never committed)."""
    path = os.environ.get(JUP365_ENV_VAR)
    if path and os.path.exists(path):
        return path
    return None


def tdb_sec_from_iso(iso: str) -> float:
    """ISO calendar string (TDB scale assumed) -> TDB seconds since J2000.

    Routed through astropy (a core dependency) rather than SPICE ``str2et``
    so no leapseconds kernel is needed: the axis is TDB end-to-end. The
    paper does not state the time scale of its printed epochs; TDB is
    assumed (a UTC reading would shift by ~69 s — immaterial for a
    same-family candidate of our own construction, recorded in the note).
    """
    from astropy.time import Time

    return float((Time(iso, scale="tdb") - Time(2451545.0, format="jd", scale="tdb")).sec)


# --- Kernel-backed moon ephemeris ------------------------------------------------

_FURNISHED: set[str] = set()


class JovianEphemeris:
    """Jupiter-centered J2000 moon states from a furnished JUP365 kernel."""

    def __init__(self, kernel_path: str) -> None:
        import spiceypy

        if not os.path.exists(kernel_path):
            raise FileNotFoundError(kernel_path)
        if kernel_path not in _FURNISHED:
            spiceypy.furnsh(kernel_path)
            _FURNISHED.add(kernel_path)
        self.kernel_path = kernel_path

    def state(self, moon: str, t_sec: float) -> tuple[Vec3, Vec3]:
        """Moon ``(r_km, v_km_s)`` relative to Jupiter (599), J2000 frame.

        ``t_sec`` is TDB seconds since J2000 == SPICE ET, passed unmodified.
        """
        import spiceypy

        st, _ = spiceypy.spkezr(moon.upper(), float(t_sec), "J2000", "NONE", "JUPITER")
        arr = np.asarray(st, dtype=np.float64)
        return arr[:3].copy(), arr[3:].copy()


class JovianRailsCache:
    """Spline-interpolated moons-on-rails positions (mirror ``RailsEphemerisCache``).

    Grid step 0.02 day: Europa (3.55 d period) is the fastest perturber; the
    cubic-spline position error at that sampling is ~1e-2 km (4th-derivative
    bound; pinned numerically by the smoke test), far below flyby tolerances.
    """

    def __init__(
        self,
        moons: Sequence[str],
        ephem: JovianEphemeris,
        t0_sec: float,
        t1_sec: float,
        *,
        step_days: float = 0.02,
        pad_days: float = 2.0,
    ) -> None:
        lo = min(t0_sec, t1_sec) - pad_days * SECONDS_PER_DAY
        hi = max(t0_sec, t1_sec) + pad_days * SECONDS_PER_DAY
        n = max(4, int((hi - lo) / (step_days * SECONDS_PER_DAY)) + 1)
        grid = np.linspace(lo, hi, n)
        self._splines: dict[str, CubicSpline] = {}
        for moon in moons:
            samples = np.empty((n, 3), dtype=np.float64)
            for i, t in enumerate(grid):
                samples[i], _ = ephem.state(moon, float(t))
            self._splines[moon] = CubicSpline(grid, samples, axis=0)

    def position(self, moon: str, t_sec: float) -> Vec3:
        return np.asarray(self._splines[moon](t_sec), dtype=np.float64)


# --- Jupiter-centered restricted n-body propagator --------------------------------


@dataclass(frozen=True)
class JovianArc:
    """Frozen result of one Jupiter-centered propagation."""

    r_km: Vec3
    v_km_s: Vec3
    t1_sec: float
    moons: tuple[str, ...]
    converged: bool


class JovianRestrictedNBody:
    """REBOUND/IAS15: massless spacecraft, central Jupiter, moons on rails.

    Mirror of :class:`cyclerfinder.nbody.propagator.RestrictedNBody` with the
    Sun/DE440 swapped for Jupiter/JUP365 and the softening clamp moved to the
    moon SURFACE (module docstring — flybys must be integrated for real).
    """

    def propagate(
        self,
        r0_km: Vec3,
        v0_km_s: Vec3,
        t0_sec: float,
        t1_sec: float,
        *,
        moons: Sequence[str] = GALILEAN,
        cache: JovianRailsCache,
        accuracy: float = 1e-11,
        max_wall_sec: float = 60.0,
    ) -> JovianArc:
        import time

        import rebound

        moons = tuple(moons)
        r0 = np.asarray(r0_km, dtype=np.float64)
        v0 = np.asarray(v0_km_s, dtype=np.float64)

        sim = rebound.Simulation()
        sim.G = MU_JUPITER_KM3_S2  # central Jupiter of mass 1.0 => G*M = mu_J
        sim.integrator = "ias15"
        sim.integrator.epsilon = accuracy
        sim.add(m=1.0, x=0.0, y=0.0, z=0.0, vx=0.0, vy=0.0, vz=0.0)
        sim.add(
            m=0.0,
            x=float(r0[0]),
            y=float(r0[1]),
            z=float(r0[2]),
            vx=float(v0[0]),
            vy=float(v0[1]),
            vz=float(v0[2]),
        )
        sim.t = float(t0_sec)

        mus = {m: SATELLITES[m].mu_km3_s2 for m in moons}
        surf = {m: SATELLITES[m].radius_eq_km for m in moons}

        def additional_forces(reb_sim_pointer: object) -> None:
            reb_sim = reb_sim_pointer.contents  # type: ignore[attr-defined]
            t_sec = float(reb_sim.t)
            sc = reb_sim.particles[1]
            r = np.array([sc.x, sc.y, sc.z], dtype=np.float64)
            ax = ay = az = 0.0
            for moon in moons:
                r_m = cache.position(moon, t_sec)
                d = r_m - r
                d_norm = float(np.linalg.norm(d))
                # Soften at the moon SURFACE only: a real >=100 km flyby is
                # integrated exactly; a below-surface dive stays finite and is
                # surfaced as a (huge) defect, never a NaN crash.
                d_eff = max(d_norm, surf[moon])
                rm3 = float(np.linalg.norm(r_m)) ** 3
                acc = mus[moon] * (d / d_eff**3 - r_m / rm3)
                ax += float(acc[0])
                ay += float(acc[1])
                az += float(acc[2])
            sc.ax += ax
            sc.ay += ay
            sc.az += az

        sim.additional_forces = additional_forces
        sim.force_is_velocity_dependent = 0

        t_target = float(t1_sec)
        arc_days = abs(t_target - float(sim.t)) / SECONDS_PER_DAY
        n_chunks = int(min(64, max(2, arc_days / 5.0)))
        chunk = (t_target - float(sim.t)) / n_chunks
        wall_start = time.monotonic()
        converged = True
        try:
            for i in range(n_chunks):
                target = t_target if i == n_chunks - 1 else float(sim.t) + chunk
                sim.integrate(target)
                sc_p = sim.particles[1]
                if not (np.isfinite(sc_p.x) and np.isfinite(sc_p.y) and np.isfinite(sc_p.z)):
                    converged = False
                    break
                if (time.monotonic() - wall_start) > max_wall_sec:
                    converged = False
                    break
        except Exception:
            converged = False

        sc = sim.particles[1]
        return JovianArc(
            r_km=np.array([sc.x, sc.y, sc.z], dtype=np.float64),
            v_km_s=np.array([sc.vx, sc.vy, sc.vz], dtype=np.float64),
            t1_sec=float(sim.t),
            moons=moons,
            converged=converged,
        )


# --- Patched-conic chain on real JUP365 geometry (Liang's own model) ---------------


def flyby_min_dv(
    vinf_in: Vec3,
    vinf_out: Vec3,
    moon: str,
    *,
    min_alt_km: float = 50.0,
) -> tuple[float, float, float]:
    """Powered-flyby defect: ``(dv_kms, turn_needed_rad, turn_max_rad)``.

    Standard turn-bounded defect (the optimizer-facing analogue of Liang
    Eqs. 5-13): rotate ``vinf_in`` toward ``vinf_out`` in their common plane
    by ``min(angle, turn_max)`` at constant magnitude; the defect is the
    remaining vector difference. ``turn_max`` from
    ``sin(d/2) = mu/(r_p_min v^2 + mu)`` at ``r_p_min = R_moon + min_alt``
    (Liang's 50 km optimizer constraint, p. 11). Zero iff an unpowered flyby
    above ``min_alt`` can deliver the turn AND the magnitudes match — i.e.
    agrees with Liang's defect exactly at ballistic convergence.
    """
    sat = SATELLITES[moon]
    vi = float(np.linalg.norm(vinf_in))
    vo = float(np.linalg.norm(vinf_out))
    if vi <= 0.0 or vo <= 0.0:
        return float(np.linalg.norm(vinf_out - vinf_in)), 0.0, 0.0
    cos_d = float(np.dot(vinf_in, vinf_out) / (vi * vo))
    delta = math.acos(max(-1.0, min(1.0, cos_d)))
    r_p_min = sat.radius_eq_km + min_alt_km
    sin_half_max = sat.mu_km3_s2 / (r_p_min * vi * vi + sat.mu_km3_s2)
    delta_max = 2.0 * math.asin(min(1.0, sin_half_max))
    turn = min(delta, delta_max)
    if delta < 1e-15:
        rotated = vinf_in
    else:
        # Rotate vinf_in by `turn` toward vinf_out in their common plane.
        e1 = vinf_in / vi
        perp = vinf_out / vo - cos_d * e1
        p_norm = float(np.linalg.norm(perp))
        if p_norm < 1e-15:
            rotated = vinf_in
        else:
            e2 = perp / p_norm
            rotated = vi * (math.cos(turn) * e1 + math.sin(turn) * e2)
    return float(np.linalg.norm(vinf_out - rotated)), delta, delta_max


def flyby_altitude_km(vinf_in: Vec3, vinf_out: Vec3, moon: str) -> float:
    """Patched-conic flyby periapsis altitude for the needed turn (km).

    ``r_p = mu (1/sin(d/2) - 1) / vinf^2`` — on REAL geometry this is the
    physical patched-conic periapsis (unlike the Tables 3/5/7 fiction it is
    still the same formula, but here the in/out asymptotes come from the real
    ephemeris legs, which is exactly the paper's own ephemeris-member
    convention for "all flybys above 100 km"). Near-zero turns map to huge
    altitudes (no real encounter needed) — reported as-is, flagged by callers.
    """
    sat = SATELLITES[moon]
    vi = float(np.linalg.norm(vinf_in))
    vo = float(np.linalg.norm(vinf_out))
    cos_d = float(np.dot(vinf_in, vinf_out) / (vi * vo))
    delta = math.acos(max(-1.0, min(1.0, cos_d)))
    sin_half = math.sin(0.5 * delta)
    if sin_half <= 1e-12:
        return float("inf")
    v2 = vi * vo  # symmetric in the (tiny) in/out magnitude mismatch
    return sat.mu_km3_s2 * (1.0 / sin_half - 1.0) / v2 - sat.radius_eq_km


def flyby_maneuver_dv(
    vinf_in: Vec3,
    vinf_out: Vec3,
    moon: str,
    *,
    alt_min: float = 25.0,
    alt_max: float = 70000.0,
) -> tuple[float, float, bool]:
    """Paper's zero-radius-SOI powered-flyby maneuver ΔV (Hernandez 2017 Eqs 3-5).

    This is the EXACT flyby model in which the paper defines EGGIE's 0.70 m/s total
    ΔV (AAS 17-608, pp.6-7), distinct from the turn-bounded leftover-vector defect
    :func:`flyby_min_dv`. Here the *bend* is provided ballistically by the gravity
    assist, and ΔV is charged ONLY for the V∞-magnitude mismatch via a tangential
    periapsis maneuver:

    * ``δ = angle(vinf_in, vinf_out)`` — the required bend (Eq 3).
    * Solve ``r_p`` from ``asin(μ/(μ + r_p·|vinf_in|²)) +
      asin(μ/(μ + r_p·|vinf_out|²)) = δ`` (Eq 4). The left side is strictly
      decreasing in ``r_p`` (from ``π`` at ``r_p→0`` to ``0`` at ``r_p→∞``), so for
      any ``δ ∈ (0, π)`` there is a unique root, found by bisection.
    * ``v_p- = sqrt(|vinf-|² + 2μ/r_p)`` (Eq 5); ``ΔV = |v_p_out - v_p_in|`` - the
      tangential burn that fixes the energy (magnitude) mismatch at periapsis.

    Feasibility: altitude ``= r_p - R_moon`` must lie in ``[alt_min, alt_max]`` km
    (paper window 25-70000 km). Equal-magnitude flybys give ``ΔV -> 0`` at a feasible
    ``r_p`` (the ballistic-cycler property: adjacent same-body encounters at equal
    V∞). A near-zero bend pushes ``r_p → ∞`` (altitude above ``alt_max`` →
    infeasible: no real encounter); a bend exceeding what ``alt_min`` can supply
    pushes ``r_p`` below ``R_moon + alt_min`` (infeasible: cannot bend enough).

    Returns ``(dv_ms, alt_km, feasible)``: ΔV in m/s, periapsis altitude in km, and
    whether the altitude lies in the window.
    """
    sat = SATELLITES[moon]
    mu = sat.mu_km3_s2
    vi = float(np.linalg.norm(vinf_in))
    vo = float(np.linalg.norm(vinf_out))
    if vi <= 0.0 or vo <= 0.0:
        return 0.0, float("inf"), False
    cos_d = float(np.dot(vinf_in, vinf_out) / (vi * vo))
    delta = math.acos(max(-1.0, min(1.0, cos_d)))

    def bend_of_rp(r_p: float) -> float:
        a_in = math.asin(min(1.0, mu / (mu + r_p * vi * vi)))
        a_out = math.asin(min(1.0, mu / (mu + r_p * vo * vo)))
        return a_in + a_out

    if delta <= 0.0:
        # No bend required → no encounter (r_p → ∞); ΔV is the pure magnitude burn,
        # but with no real flyby it is infeasible (altitude above the window).
        return abs(vo - vi) * 1.0e3, float("inf"), False

    # Bisect the strictly-decreasing bend_of_rp(r_p) - delta for the unique r_p.
    lo, hi = 1.0e-6, 1.0e15
    # Guard: if even r_p→0 cannot reach delta (delta ≥ π) it is unachievable.
    if bend_of_rp(lo) < delta:
        # delta ~ π: deepest possible pass still under-bends → infeasible.
        r_p = lo
    else:
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            if bend_of_rp(mid) > delta:
                lo = mid
            else:
                hi = mid
            if hi - lo < 1.0e-9 * max(1.0, hi):
                break
        r_p = 0.5 * (lo + hi)

    alt_km = r_p - sat.radius_eq_km
    v_p_in = math.sqrt(vi * vi + 2.0 * mu / r_p)
    v_p_out = math.sqrt(vo * vo + 2.0 * mu / r_p)
    dv_ms = abs(v_p_out - v_p_in) * 1.0e3
    feasible = alt_min <= alt_km <= alt_max
    return dv_ms, alt_km, feasible


@dataclass(frozen=True)
class ConicCycle:
    """One converged patched-conic CGCEC cycle on real JUP365 geometry."""

    index: int  # 1-based cycle number
    epochs_sec: tuple[float, ...]  # 5 flyby epochs (C, G, C, E, C)
    tofs_days: tuple[float, float, float, float]
    cycle_tof_days: float
    revs: tuple[int, ...]
    branches: tuple[str, ...]
    vinf_kms: tuple[float, ...]  # |V_inf| inbound at each flyby (out at flyby 0)
    defects_ms: tuple[float, ...]  # defect dv (m/s) at the flybys charged to this cycle
    altitudes_km: tuple[float, ...]  # patched-conic periapsis altitude at interior flybys
    sum_defect_ms: float
    converged: bool


def _solve_cycle_legs(
    epochs: Sequence[float],
    ephem: JovianEphemeris,
    branch_plan: Sequence[tuple[int, str]],
) -> tuple[list[Vec3], list[Vec3]] | None:
    """Lambert-solve the 4 legs at the given epochs; returns (vinf_out, vinf_in) lists.

    ``vinf_out[k]`` is the departure V∞ vector at flyby k (k=0..3);
    ``vinf_in[k]`` the arrival V∞ at flyby k (k=1..4, list index k-1).
    Returns None when any leg has no Lambert solution on its planned branch.
    """
    states = [ephem.state(m, t) for m, t in zip(CGCEC, epochs, strict=True)]
    vinf_out: list[Vec3] = []
    vinf_in: list[Vec3] = []
    for k in range(4):
        tof = epochs[k + 1] - epochs[k]
        if tof <= 0.0:
            return None
        n_revs, branch = branch_plan[k]
        try:
            sols = lambert(
                states[k][0],
                states[k + 1][0],
                tof,
                mu=MU_JUPITER_KM3_S2,
                max_revs=max(n_revs, 1),
            )
        except Exception:
            return None
        match = [s for s in sols if s.n_revs == n_revs and s.branch == branch]
        if not match:
            return None
        sol = match[0]
        vinf_out.append(np.asarray(sol.v1, dtype=np.float64) - states[k][1])
        vinf_in.append(np.asarray(sol.v2, dtype=np.float64) - states[k + 1][1])
    return vinf_out, vinf_in


#: Per-leg (n_revs, branch) plan validated for the idealized members (#222):
#: every Liang leg is a 1-rev Lambert solution, branch pattern high/low/high/low.
BRANCH_PLAN: tuple[tuple[int, str], ...] = ((1, "high"), (1, "low"), (1, "high"), (1, "low"))

_DEFECT_SENTINEL_MS = 1.0e6  # honest large-finite defect for failed Lambert legs


def optimize_cycle(
    t_start_sec: float,
    tof_seed_days: Sequence[float],
    ephem: JovianEphemeris,
    *,
    cycle_index: int,
    vinf_in_prev: Vec3 | None,
    bound_days: float = 3.0,
    min_alt_km: float = 50.0,
) -> tuple[ConicCycle, Vec3]:
    """Locally optimize one cycle's interior epochs to minimize defect Δv.

    Mirrors Liang Sec. III.C (ephemeris variant, Eq. 20): the start epoch and
    inherited inbound V∞ are FIXED; the four downstream flyby epochs move
    within ``+/- bound_days`` of the ToF seed; objective = the defect Δv at
    the flybys charged to this cycle (start Callisto when an inbound V∞ is
    inherited, plus Ganymede / Callisto / Europa interior flybys). Solver:
    Nelder-Mead on the summed defect (the defect surface is piecewise-smooth
    with flat zero plateaus — exactly what their DE solver tolerated; NM is
    the lightweight local analogue), polished by a second NM restart.

    Returns the converged cycle record and the outbound V∞ vector at the
    cycle-end Callisto flyby's INBOUND side (i.e. the arrival V∞ the next
    cycle inherits).
    """
    from scipy.optimize import minimize

    tof_seed = np.asarray(tof_seed_days, dtype=np.float64)

    def defects_of(x_tofs: NDArray[np.float64]) -> tuple[list[float], object]:
        epochs = [t_start_sec]
        for tof in x_tofs:
            epochs.append(epochs[-1] + float(tof) * SECONDS_PER_DAY)
        legs = _solve_cycle_legs(epochs, ephem, BRANCH_PLAN)
        if legs is None:
            return [_DEFECT_SENTINEL_MS] * 4, None
        vinf_out, vinf_in = legs
        ds: list[float] = []
        if vinf_in_prev is not None:
            dv, _, _ = flyby_min_dv(vinf_in_prev, vinf_out[0], CGCEC[0], min_alt_km=min_alt_km)
            ds.append(dv * 1.0e3)
        for k in (1, 2, 3):  # G, C, E interior flybys
            dv, _, _ = flyby_min_dv(vinf_in[k - 1], vinf_out[k], CGCEC[k], min_alt_km=min_alt_km)
            ds.append(dv * 1.0e3)
        return ds, (epochs, vinf_out, vinf_in)

    def objective(x: NDArray[np.float64]) -> float:
        if np.any(np.abs(x - tof_seed) > bound_days):
            return _DEFECT_SENTINEL_MS * 10.0
        ds, _ = defects_of(x)
        return float(sum(ds))

    x = tof_seed.copy()
    for _ in range(2):
        res = minimize(
            objective,
            x,
            method="Nelder-Mead",
            options={"xatol": 1e-10, "fatol": 1e-12, "maxiter": 4000, "adaptive": True},
        )
        x = np.asarray(res.x, dtype=np.float64)

    ds, detail = defects_of(x)
    if detail is None:
        empty: Vec3 = np.zeros(3)
        return (
            ConicCycle(
                index=cycle_index,
                epochs_sec=(t_start_sec,) * 5,
                tofs_days=(0.0, 0.0, 0.0, 0.0),
                cycle_tof_days=0.0,
                revs=tuple(r for r, _ in BRANCH_PLAN),
                branches=tuple(b for _, b in BRANCH_PLAN),
                vinf_kms=(0.0,) * 5,
                defects_ms=tuple(ds),
                altitudes_km=(0.0,) * 3,
                sum_defect_ms=float(sum(ds)),
                converged=False,
            ),
            empty,
        )
    epochs, vinf_out, vinf_in = cast("tuple[list[float], list[Vec3], list[Vec3]]", detail)
    vinfs = [float(np.linalg.norm(vinf_out[0]))]
    vinfs += [float(np.linalg.norm(vinf_in[k])) for k in range(4)]
    alts = tuple(flyby_altitude_km(vinf_in[k - 1], vinf_out[k], CGCEC[k]) for k in (1, 2, 3))
    cycle = ConicCycle(
        index=cycle_index,
        epochs_sec=tuple(float(t) for t in epochs),
        tofs_days=tuple(float(t) for t in x),  # type: ignore[arg-type]
        cycle_tof_days=float((epochs[-1] - epochs[0]) / SECONDS_PER_DAY),
        revs=tuple(r for r, _ in BRANCH_PLAN),
        branches=tuple(b for _, b in BRANCH_PLAN),
        vinf_kms=tuple(vinfs),
        defects_ms=tuple(ds),
        altitudes_km=alts,
        sum_defect_ms=float(sum(ds)),
        converged=bool(sum(ds) < 1.0),  # < 1 m/s total — reported, not asserted
    )
    return cycle, vinf_in[3]


def chain_cycles(
    t0_sec: float,
    ephem: JovianEphemeris,
    *,
    n_cycles: int = 10,
    tof_seed_days: Sequence[float] = (31.8973, 18.1697, 29.9343, 19.9747),
    bound_days: float = 3.0,
    min_alt_km: float = 50.0,
    progress: bool = False,
) -> list[ConicCycle]:
    """Chain ``n_cycles`` CGCEC cycles from ``t0_sec`` (Liang's chained scheme).

    The default ToF seed is idealized Member A's printed first-cycle leg ToFs
    (Table 3) — the validated scaffold (#222) phased to the real epoch; each
    later cycle re-seeds from the previous cycle's converged ToFs.
    """
    cycles: list[ConicCycle] = []
    t_start = t0_sec
    vinf_prev: Vec3 | None = None
    seed = list(tof_seed_days)
    for i in range(1, n_cycles + 1):
        cycle, vinf_end = optimize_cycle(
            t_start,
            seed,
            ephem,
            cycle_index=i,
            vinf_in_prev=vinf_prev,
            bound_days=bound_days,
            min_alt_km=min_alt_km,
        )
        cycles.append(cycle)
        if progress:
            print(
                f"  cycle {i}: tof={cycle.cycle_tof_days:.4f} d, "
                f"sum_defect={cycle.sum_defect_ms:.3e} m/s, "
                f"vinf={[f'{v:.4f}' for v in cycle.vinf_kms]}",
                flush=True,
            )
        if not cycle.converged and cycle.cycle_tof_days <= 0.0:
            break
        t_start = cycle.epochs_sec[-1]
        vinf_prev = vinf_end
        seed = list(cycle.tofs_days)
    return cycles


# --- Patched-conic -> n-body node conversion ---------------------------------------


def periapsis_node(
    moon: str,
    t_sec: float,
    vinf_in: Vec3,
    vinf_out: Vec3,
    ephem: JovianEphemeris,
    *,
    max_offset_km: float | None = None,
) -> tuple[Vec3, Vec3, float]:
    """Jupiter-frame spacecraft state at the flyby periapsis; returns (r, v, d_km).

    Builds the moon-centered hyperbola that turns ``vinf_in`` into
    ``vinf_out`` (patched-conic convention): periapsis distance from the turn
    angle, periapsis direction along the asymptote-bisector
    ``(vin_hat - vout_hat)``, periapsis speed ``sqrt(vinf^2 + 2 mu/r_p)``
    perpendicular to it in the flyby plane. For near-zero turns the implied
    ``r_p`` exceeds the moon's gravitational reach; the offset is clamped to
    ``max_offset_km`` (default 0.6x the moon's Laplace SOI) — there the moon
    term is negligible, so the clamp only keeps the node geometrically sane.
    """
    sat = SATELLITES[moon]
    mu = sat.mu_km3_s2
    r_m, v_m = ephem.state(moon, t_sec)
    vi = float(np.linalg.norm(vinf_in))
    vo = float(np.linalg.norm(vinf_out))
    vmag = 0.5 * (vi + vo)
    cos_d = float(np.dot(vinf_in, vinf_out) / (vi * vo))
    delta = math.acos(max(-1.0, min(1.0, cos_d)))
    sin_half = math.sin(0.5 * delta)
    soi_km = sat.sma_km * (mu / (3.0 * MU_JUPITER_KM3_S2)) ** (1.0 / 3.0)
    cap = max_offset_km if max_offset_km is not None else 0.6 * soi_km
    r_p = mu * (1.0 / sin_half - 1.0) / (vmag * vmag) if sin_half > 1e-12 else float("inf")
    d = min(r_p, cap)
    e_in = vinf_in / vi
    e_out = vinf_out / vo
    bisector = e_in - e_out
    b_norm = float(np.linalg.norm(bisector))
    if b_norm < 1e-12:  # zero turn: offset perpendicular to the asymptote
        ref = np.array([0.0, 0.0, 1.0])
        perp = np.cross(e_in, ref)
        bisector = perp
        b_norm = float(np.linalg.norm(perp))
    e_p = bisector / b_norm
    # Periapsis velocity direction: perpendicular to e_p in the (e_in, e_out)
    # plane, pointing from the inbound toward the outbound asymptote.
    v_dir = e_in + e_out
    v_dir = v_dir - float(np.dot(v_dir, e_p)) * e_p
    v_norm = float(np.linalg.norm(v_dir))
    if v_norm < 1e-12:
        v_dir = e_in - float(np.dot(e_in, e_p)) * e_p
        v_norm = float(np.linalg.norm(v_dir))
    e_v = v_dir / v_norm
    v_p = math.sqrt(vmag * vmag + 2.0 * mu / d)
    return r_m + d * e_p, v_m + v_p * e_v, d


# --- Per-cycle multiple shooting in the Jovian n-body model ------------------------

_W_VEL = 1.0e3  # velocity residual weight: 1 m/s ~ 1 km position defect


@dataclass(frozen=True)
class JovianShootResult:
    """Per-cycle multiple-shooting record (honest, never raised)."""

    cycle_index: int
    converged: bool
    seed_leg_defects_km: tuple[float, ...]  # |dr| per leg at the seed
    seed_leg_defects_kms: tuple[float, ...]  # |dv| per leg at the seed
    final_leg_defects_km: tuple[float, ...]
    final_leg_defects_kms: tuple[float, ...]
    moon_distances_km: tuple[float, ...]  # |r_node - r_moon| at nodes 1..4 (final)
    boundary_dv_ms: float  # V∞-magnitude mismatch at the cycle-start flyby (m/s)
    node_states: tuple[Vec3, ...]
    epochs_sec: tuple[float, ...]
    nfev: int


def _cycle_residual(
    x: NDArray[np.float64],
    node0: Vec3,
    epochs_seed: Sequence[float],
    moons_at_nodes: Sequence[str],
    d_caps_km: Sequence[float],
    prop: JovianRestrictedNBody,
    cache: JovianRailsCache,
    ephem: JovianEphemeris,
    vinf_in_mag_prev: float | None,
    accuracy: float,
) -> NDArray[np.float64]:
    """Residual: leg continuity + encounter-proximity hinges + boundary V∞ pin.

    Free vector ``x`` = [node0 velocity (3), nodes 1..4 full states (24),
    epoch offsets for nodes 1..4 in days (4)] — 31 vars. Node 0's POSITION is
    pinned (the patched-conic periapsis at the published departure epoch /
    the previous cycle's corrected end node). Residual blocks:

    1. 4 legs x 6 continuity defects (km; velocity x ``_W_VEL``).
    2. nodes 1..4: hinge ``max(0, |r - r_moon| - d_cap)`` (km) — the
       encounter-band constraint that makes "it closed" mean "it closed
       THROUGH the moons" (orbit-closure discipline: all binding constraints
       in the residual).
    3. when a previous cycle hands over an inbound V∞ magnitude: the
       cycle-start V∞-magnitude mismatch (m/s ~ km weight) — an unpowered
       flyby preserves |V∞|; the residual mismatch IS the per-cycle Δv.
    """
    v0 = x[0:3]
    nodes = [np.concatenate([node0[:3], v0])]
    for k in range(4):
        nodes.append(x[3 + 6 * k : 3 + 6 * (k + 1)])
    epochs = [epochs_seed[0]]
    for k in range(4):
        epochs.append(epochs_seed[k + 1] + x[27 + k] * SECONDS_PER_DAY)

    res: list[float] = []
    for k in range(4):
        arc = prop.propagate(
            nodes[k][:3],
            nodes[k][3:],
            epochs[k],
            epochs[k + 1],
            cache=cache,
            accuracy=accuracy,
        )
        if arc.converged:
            dr = arc.r_km - nodes[k + 1][:3]
            dv = arc.v_km_s - nodes[k + 1][3:]
            res.extend(float(c) for c in dr)
            res.extend(float(c) * _W_VEL for c in dv)
        else:
            res.extend([1e7] * 6)

    for k in range(1, 5):
        r_m, _ = ephem.state(moons_at_nodes[k], epochs[k])
        dist = float(np.linalg.norm(nodes[k][:3] - r_m))
        res.append(max(0.0, dist - d_caps_km[k]))

    if vinf_in_mag_prev is not None:
        _, v_m = ephem.state(moons_at_nodes[0], epochs[0])
        vinf0 = float(np.linalg.norm(nodes[0][3:] - v_m))
        res.append((vinf0 - vinf_in_mag_prev) * _W_VEL)

    return np.asarray(res, dtype=np.float64)


def shoot_cycle(
    cycle_nodes: Sequence[tuple[Vec3, Vec3]],  # 5 (r, v) periapsis states
    epochs_sec: Sequence[float],
    ephem: JovianEphemeris,
    cache: JovianRailsCache,
    *,
    cycle_index: int,
    d_caps_km: Sequence[float],
    vinf_in_mag_prev: float | None,
    accuracy: float = 1e-11,
    max_nfev: int = 60,
) -> JovianShootResult:
    """Multiple-shoot one CGCEC cycle to n-body continuity (see _cycle_residual)."""
    from scipy.optimize import least_squares

    prop = JovianRestrictedNBody()
    node0 = np.concatenate([cycle_nodes[0][0], cycle_nodes[0][1]])
    x0 = np.concatenate(
        [cycle_nodes[0][1]] + [np.concatenate([r, v]) for r, v in cycle_nodes[1:]] + [np.zeros(4)]
    )

    def residual(x: NDArray[np.float64]) -> NDArray[np.float64]:
        return _cycle_residual(
            x,
            node0,
            epochs_sec,
            CGCEC,
            d_caps_km,
            prop,
            cache,
            ephem,
            vinf_in_mag_prev,
            accuracy,
        )

    f0 = residual(x0)
    seed_dr = tuple(float(np.linalg.norm(f0[6 * k : 6 * k + 3])) for k in range(4))
    seed_dv = tuple(float(np.linalg.norm(f0[6 * k + 3 : 6 * k + 6]) / _W_VEL) for k in range(4))

    sol = least_squares(
        residual,
        x0,
        method="trf",
        x_scale="jac",
        max_nfev=max_nfev,
        diff_step=1e-7,  # type: ignore[arg-type]
    )
    xf = np.asarray(sol.x, dtype=np.float64)
    ff = residual(xf)
    fin_dr = tuple(float(np.linalg.norm(ff[6 * k : 6 * k + 3])) for k in range(4))
    fin_dv = tuple(float(np.linalg.norm(ff[6 * k + 3 : 6 * k + 6]) / _W_VEL) for k in range(4))

    nodes = [np.concatenate([node0[:3], xf[0:3]])]
    for k in range(4):
        nodes.append(xf[3 + 6 * k : 3 + 6 * (k + 1)])
    epochs = [epochs_sec[0]]
    for k in range(4):
        epochs.append(epochs_sec[k + 1] + xf[27 + k] * SECONDS_PER_DAY)
    dists = []
    for k in range(1, 5):
        r_m, _ = ephem.state(CGCEC[k], epochs[k])
        dists.append(float(np.linalg.norm(nodes[k][:3] - r_m)))
    boundary_dv = 0.0
    if vinf_in_mag_prev is not None:
        _, v_m = ephem.state(CGCEC[0], epochs[0])
        boundary_dv = abs(float(np.linalg.norm(nodes[0][3:] - v_m)) - vinf_in_mag_prev) * 1e3

    pos_ok = all(d < 1.0 for d in fin_dr)  # km
    vel_ok = all(d < 1e-5 for d in fin_dv)  # km/s
    return JovianShootResult(
        cycle_index=cycle_index,
        converged=bool(pos_ok and vel_ok),
        seed_leg_defects_km=seed_dr,
        seed_leg_defects_kms=seed_dv,
        final_leg_defects_km=fin_dr,
        final_leg_defects_kms=fin_dv,
        moon_distances_km=tuple(dists),
        boundary_dv_ms=boundary_dv,
        node_states=tuple(nodes),
        epochs_sec=tuple(epochs),
        nfev=int(sol.nfev),
    )


# --- Jupiter-central multiple-shooting corrector (#480 M1) -------------------------
#
# A Jupiter-centred analogue of :func:`cyclerfinder.nbody.shooter.shoot`. The shared
# heliocentric ``shoot()`` is hardwired to ``RestrictedNBody`` (central MU_SUN, the
# ``PLANETS`` perturber registry) and so cannot integrate a Galilean-moon tour
# (KeyError 'Io' + IAS15 step-collapse). Rather than parameterise that core module
# — which would ripple through ``defect_residual`` / ``_stm_jacobian`` / the worker
# pool and every frozen heliocentric golden — this is a CONTAINED corrector that
# reuses the same multiple-shooting structure (``ShootingSeed`` / ``ShootResult`` /
# ``least_squares``) but propagates each leg with :class:`JovianRestrictedNBody`.
# The heliocentric ``shoot()`` is untouched.

_STATE_DIM = 6

# SNOPT-analogue continuity floors (mirror shooter.py: Jones AAS 17-577 §2.5,
# 1e-3 km position / 1e-6 km/s velocity); the Jovian leg defects are compared to
# the same per-component floor scaled by sqrt(n_legs).
_POS_CONTINUITY_KM = 1.0e-3
_VEL_CONTINUITY_KMS = 1.0e-6

# Velocity-residual weight so a 1 m/s velocity defect carries ~1 km of position
# defect in the least-squares norm (mirror shoot_cycle's _W_VEL above).
_W_VEL = 1.0e3


def _jovian_flyby_hinge_km(vinf_in: Vec3, vinf_out: Vec3, moon: str) -> float:
    """Bend-feasibility shortfall (km) for a Galilean flyby (mirror shooter hinge).

    An unpowered flyby turns ``vinf_in`` into ``vinf_out`` at constant magnitude;
    the periapsis needed is ``r_p = (mu/v∞^2)(1/sin(δ/2) - 1)``. Returns
    ``max(0, r_safe - r_p)`` where ``r_safe = radius_eq_km + min_alt`` (100 km, the
    paper's flyby floor). Zero when the bend is feasible above the surface.
    """
    sat = SATELLITES[moon]
    r_safe = sat.radius_eq_km + 100.0
    vi = float(np.linalg.norm(vinf_in))
    vo = float(np.linalg.norm(vinf_out))
    vmag = 0.5 * (vi + vo)
    if vmag <= 0.0:
        return 0.0
    cos_d = float(np.dot(vinf_in, vinf_out) / (vi * vo)) if vi > 0 and vo > 0 else 1.0
    cos_d = max(-1.0, min(1.0, cos_d))
    delta = math.acos(cos_d)
    sin_half = math.sin(0.5 * delta)
    if sin_half <= 1e-12:
        return 0.0
    r_p = (sat.mu_km3_s2 / (vmag * vmag)) * (1.0 / sin_half - 1.0)
    return max(0.0, r_safe - r_p)


def jovian_defect_residual(
    seed: ShootingSeed,
    *,
    ephem: JovianEphemeris,
    cache: JovianRailsCache,
    moons: Sequence[str] = GALILEAN,
    accuracy: float = 1e-11,
    max_wall_sec: float = 30.0,
) -> NDArray[np.float64]:
    """Full-state multiple-shooting residual in the Jupiter-central model.

    Jupiter-centred analogue of
    :func:`cyclerfinder.nbody.shooter.defect_residual`: same residual layout
    (leg continuity, flyby hinges, periodicity wrap) but each interior leg is
    propagated with :class:`JovianRestrictedNBody` against moons-on-jup365-rails.
    Velocity defects are weighted by ``_W_VEL`` so position and velocity enter the
    least-squares norm on comparable scales. Divergence is an honest large-finite
    sentinel, never a NaN/raise (mirror the heliocentric residual).
    """
    prop = JovianRestrictedNBody()
    n = len(seed.sequence)
    res: list[float] = []

    # 1. Leg continuity defects.
    for i in range(n - 1):
        s_i = seed.node_states[i]
        r0 = np.asarray(s_i[:3], dtype=np.float64)
        v0 = np.asarray(s_i[3:], dtype=np.float64)
        arc = prop.propagate(
            r0,
            v0,
            seed.epochs[i],
            seed.epochs[i + 1],
            moons=moons,
            cache=cache,
            accuracy=accuracy,
            max_wall_sec=max_wall_sec,
        )
        s_next = seed.node_states[i + 1]
        if arc.converged and np.all(np.isfinite(arc.r_km)) and np.all(np.isfinite(arc.v_km_s)):
            dr = arc.r_km - np.asarray(s_next[:3], dtype=np.float64)
            dv = arc.v_km_s - np.asarray(s_next[3:], dtype=np.float64)
            res.extend(float(x) for x in dr)
            res.extend(float(x) * _W_VEL for x in dv)
        else:
            res.extend([1e9] * _STATE_DIM)

    # 2. Flyby bend-feasibility hinges (interior encounters only).
    for i in range(1, n - 1):
        res.append(_jovian_flyby_hinge_km(seed.vinf_in[i], seed.vinf_out[i], seed.sequence[i]))

    # 3. Periodicity wrap in the home-moon-relative frame (a pure moon-ephemeris
    #    shift of the home moon is not charged as a defect).
    r_home, v_home = ephem.state(seed.sequence[0], seed.epochs[0])
    r_wrap_pl, v_wrap_pl = ephem.state(seed.sequence[-1], seed.epochs[-1])
    s0 = seed.node_states[0]
    sn = seed.node_states[-1]
    rel0_r = np.asarray(s0[:3], dtype=np.float64) - np.asarray(r_home, dtype=np.float64)
    rel0_v = np.asarray(s0[3:], dtype=np.float64) - np.asarray(v_home, dtype=np.float64)
    reln_r = np.asarray(sn[:3], dtype=np.float64) - np.asarray(r_wrap_pl, dtype=np.float64)
    reln_v = np.asarray(sn[3:], dtype=np.float64) - np.asarray(v_wrap_pl, dtype=np.float64)
    res.extend(float(x) for x in (reln_r - rel0_r))
    res.extend(float(x) * _W_VEL for x in (reln_v - rel0_v))

    return np.asarray(res, dtype=np.float64)


def jovian_shoot(
    seed: ShootingSeed,
    *,
    kernel_path: str | None = None,
    moons: Sequence[str] = GALILEAN,
    accuracy: float = 1e-11,
    max_nfev: int = 40,
    max_wall_sec: float = 30.0,
) -> ShootResult:
    """Multiple-shooting corrector in the Jupiter-central model (#480 M1).

    Drop-in Jupiter analogue of :func:`cyclerfinder.nbody.shooter.shoot` for a
    Galilean-moon tour (e.g. the Hernandez 2017 EGGIE cycler). Free variables are
    the per-node full Cartesian states (6 each); node epochs are held at the seed
    epochs (same documented choice as the heliocentric shooter). Each leg is
    propagated with :class:`JovianRestrictedNBody` (central MU_JUPITER, moons on
    jup365 rails) — so the corrector RUNS where the heliocentric ``shoot()`` raises
    KeyError on 'Io'. Returns a :class:`~cyclerfinder.nbody.shooter.ShootResult`;
    divergence is recorded honestly (mirror ``shoot()``), never raised.

    ``kernel_path`` resolves the JUP365 SPK (default: the ensured generic kernel,
    the same one ``Ephemeris(center='Jupiter', model='spice')`` furnishes), so the
    propagator's rails match the seed's moon geometry.
    """
    from scipy.optimize import least_squares

    from cyclerfinder.nbody.correction_dv import node_impulse_correction_dv
    from cyclerfinder.nbody.shooter import (
        ShootResult,
        _seed_with_states,
        _states_to_x,
        _x_to_states,
    )

    if kernel_path is None:
        from cyclerfinder.verify.spice_kernels import ensure_jup365_kernel

        kernel_path = ensure_jup365_kernel()
    jeph = JovianEphemeris(kernel_path)

    moons = tuple(moons)
    n = len(seed.sequence)
    cache = JovianRailsCache(moons, jeph, min(seed.epochs), max(seed.epochs))

    def residual_of_x(x: NDArray[np.float64]) -> NDArray[np.float64]:
        trial = _seed_with_states(seed, _x_to_states(x, n))
        return jovian_defect_residual(
            trial,
            ephem=jeph,
            cache=cache,
            moons=moons,
            accuracy=accuracy,
            max_wall_sec=max_wall_sec,
        )

    x0 = _states_to_x(seed.node_states)
    seed_res = residual_of_x(x0)
    seed_defect_norm = float(np.linalg.norm(seed_res))

    sol = least_squares(residual_of_x, x0, method="trf", x_scale="jac", max_nfev=max_nfev)
    corrected_states = _x_to_states(np.asarray(sol.x, dtype=np.float64), n)
    final_res = residual_of_x(np.asarray(sol.x, dtype=np.float64))
    defect_norm = float(np.linalg.norm(final_res))

    # Continuity acceptance reads only the leg-defect block (velocity weighted).
    n_leg_defects = (n - 1) * _STATE_DIM
    leg_block = final_res[:n_leg_defects]
    floor_vec = ([_POS_CONTINUITY_KM] * 3 + [_VEL_CONTINUITY_KMS * _W_VEL] * 3) * (n - 1)
    continuity_floor = float(np.linalg.norm(floor_vec))
    converged = bool(sol.success) and float(np.linalg.norm(leg_block)) < continuity_floor

    vinf_mag: list[float] = []
    for state, body, epoch in zip(corrected_states, seed.sequence, seed.epochs, strict=True):
        _, v_m = jeph.state(body, epoch)
        vinf_mag.append(
            float(np.linalg.norm(np.asarray(state[3:], dtype=np.float64) - np.asarray(v_m)))
        )

    seed_nodes = {f"b{i}": np.asarray(seed.node_states[i][3:]) for i in range(n)}
    corr_nodes = {f"b{i}": np.asarray(corrected_states[i][3:]) for i in range(n)}
    cdv = node_impulse_correction_dv(seed_nodes, corr_nodes)

    bend_feasible = all(
        _jovian_flyby_hinge_km(seed.vinf_in[i], seed.vinf_out[i], seed.sequence[i]) <= 0.0
        for i in range(1, n - 1)
    )

    return ShootResult(
        converged=converged,
        defect_norm=defect_norm,
        seed_defect_norm=seed_defect_norm,
        corrected_states=corrected_states,
        vinf_per_encounter_kms=vinf_mag,
        correction_dv_kms=cdv.total_kms,
        bend_feasible=bend_feasible,
        sequence=seed.sequence,
        integrator_accuracy=accuracy,
        n_iterations=int(sol.nfev),
    )


__all__ = [
    "BRANCH_PLAN",
    "CGCEC",
    "GALILEAN",
    "MU_JUPITER_KM3_S2",
    "ConicCycle",
    "JovianArc",
    "JovianEphemeris",
    "JovianRailsCache",
    "JovianRestrictedNBody",
    "JovianShootResult",
    "chain_cycles",
    "flyby_altitude_km",
    "flyby_min_dv",
    "jovian_defect_residual",
    "jovian_shoot",
    "jup365_kernel_path",
    "optimize_cycle",
    "periapsis_node",
    "shoot_cycle",
    "tdb_sec_from_iso",
]
