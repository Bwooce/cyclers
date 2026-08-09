"""Vaquero 2013 Earth-Moon 2:1 / 3:1 "Periodic Cycler" family reproduction (#799).

Direct CR3BP family-continuation reproduction of the two planar resonant
"Earth-Moon Periodic Cycler" families of Vaquero (2013), Purdue Ph.D.
dissertation, Sec. 4.4.7 (pp. 169-172, Figs. 4.43-4.44). Full section
digest: ``docs/notes/2026-08-08-787-vaquero-2013-earth-moon-chapter-digest.md``.

Why continuation is the ONLY route here (`#787` digest, headline finding 3):
neither Vaquero nor Casoliva print a digit-grade IC table for these two
specific families anywhere -- Fig. 4.44's ``x0`` axis (the only per-member
numeric handle she ever plots) carries no printed value labels. The sourced
ground truth is therefore her PROSE: the two families' own Jacobi-constant
ranges, the Earth->Moon time-of-flight values at the four range endpoints,
and the four qualitative selection criteria (below). Reproduction means:
converge a family whose members land inside her printed ``C`` range AND
satisfy her own printed criteria AND hit her printed endpoint TOFs.

Sourced prose values (Sec. 4.4.7, p. 171; digest section 4)
-----------------------------------------------------------
========  ======================  =============================
family    Jacobi-constant range   Earth->Moon TOF at endpoints
========  ======================  =============================
2:1       ``C in [1.98, 2.66]``   6.39 d at C=1.98, 4.91 d at C=2.66
3:1       ``C in [2.54, 3.13]``   4.90 d at C=2.54, 5.04 d at C=3.13
========  ======================  =============================

Vaquero's own four selection criteria (p. 169-170, quoted in the digest):

1. **Close Earth approach**: insertion from a circular LEO with
   ``180 km <= r - r_E <= 35,786 km`` (GEO altitude) -- Earth-center
   perigee in ``[6558.137, 42164.137]`` km (the sums are DERIVED from her
   printed altitudes + the conventional ``R_E = 6378.137`` km, not
   themselves printed -- same derivation `#798`'s note documents).
2. **Close Moon approach**: surface contact or a small-amplitude L1/L2
   LPO connection (3:1 is her cislunar/L1 family, 2:1 her circumlunar/L2
   family).
3. **Period**: Earth->Moon time-of-flight ``<= 7`` days.
4. **Stability**: "only resonant cyclers that are stable or possess small
   unstable modes are considered" (direct quote, p. 170 -- NOT a strict
   stability requirement; Fig. 4.44's own stability coloring spans the
   instability threshold within both plotted families).

Seeding strategy (two-body, per Vaquero's own Ch. 3.5.1 method lineage)
-----------------------------------------------------------------------
Her ``p:q`` convention (Eq. 3.1, p. 83-84; digest section 6): the
spacecraft completes ``p`` revolutions while the Moon completes ``q``;
``p > q`` = interior. Both cycler families are ``p:1`` interior resonances,
so the two-body ellipse has ``a = p**(-2/3)`` (in Earth-Moon distance
units) and full nondimensional period ``2*pi*q = 2*pi``. Rather than the
project's plain :func:`~cyclerfinder.search.jovian_resonant_families.two_body_resonant_seed`
(periapsis pinned AT the secondary's radius -- an exterior-flavored
construction that the four-system lineage note,
:data:`~cyclerfinder.search.earth_moon_resonant_families.TWO_BODY_SEED_LINEAGE_NOTE`,
shows never lands on the intended topology), the seed here pins the
ellipse's APOAPSIS perpendicular crossing on the ``+x`` axis (the Moon
side) and fixes the eccentricity from the target Jacobi constant via the
Tisserand relation ``C ~= 1/a + 2*sqrt(a*(1-e**2))`` (two-body
approximation of the Jacobi constant, prograde planar) -- exactly the
geometry of Vaquero's own families: perigee deep in the LEO-GEO band,
apogee at the Moon's vicinity (3:1 apogee ~0.89 inside the Moon's orbit =
cislunar; 2:1 apogee ~1.19 beyond it = circumlunar), which is also
precisely Fig. 4.44's plotted ``x0 in [~0.6, ~1.2+]`` band. Vaquero's
frame convention (Purdue/Howell, Earth at ``(-mu, 0)``, Moon at
``(1-mu, 0)``) matches this project's own :mod:`cyclerfinder.core.cr3bp`,
so ``x0`` values are directly comparable with no flip (unlike the
Casoliva-convention flip `#780` documents).

The seed is corrected at FIXED Jacobi constant by the existing
:func:`~cyclerfinder.search.cr3bp_periodic.correct_symmetric_fixed_jacobi`
and continued across Vaquero's own printed ``C`` range by the existing
:func:`~cyclerfinder.search.cr3bp_continuation.continue_family` (which
runs every member through the full convergence gauntlet including the
independent-Radau cross-check) -- no new continuation machinery.

Scope / novelty
---------------
This module REPRODUCES a described-but-untabulated PUBLISHED family.
Nothing here is novel and nothing is claimed novel; the
``literature_check.py`` novelty gate is therefore not in play (same
precedent as `#780`'s own module: reproduction of a named published
family, scoped away from any novelty claim). No catalogue writeback is
performed here (writeback is a separate, later task).

Pure: math/numpy/scipy + :mod:`cyclerfinder.core.cr3bp`,
:mod:`cyclerfinder.core.crnbp` (collinear-point helper),
:mod:`cyclerfinder.search.cr3bp_periodic`,
:mod:`cyclerfinder.search.cr3bp_continuation`.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import solve_ivp

import cyclerfinder.core.cr3bp as cr3bp
import cyclerfinder.core.crnbp as crnbp
import cyclerfinder.search.cr3bp_continuation as cc
import cyclerfinder.search.cr3bp_periodic as cp

TWO_PI = 2.0 * math.pi

# ---------------------------------------------------------------------------
# Sourced constants (Vaquero 2013, Sec. 4.4.7, p. 171 prose -- via the #787
# digest; there is NO printed table for these families, see module docstring).
# ---------------------------------------------------------------------------

#: 2:1 family Jacobi-constant range (SOURCED, p. 171).
VAQUERO_C_RANGE_21: tuple[float, float] = (1.98, 2.66)

#: 3:1 family Jacobi-constant range (SOURCED, p. 171).
VAQUERO_C_RANGE_31: tuple[float, float] = (2.54, 3.13)

#: Earth->Moon TOF at the 2:1 family's own printed range endpoints,
#: ``{C: days}`` (SOURCED, p. 171: "6.39 days (C=1.98)" / "4.91 days (C=2.66)").
VAQUERO_TOF_DAYS_21: dict[float, float] = {1.98: 6.39, 2.66: 4.91}

#: Earth->Moon TOF at the 3:1 family's own printed range endpoints (SOURCED,
#: p. 171: "4.90 days (C=2.54)" / "5.04 days (C=3.13)").
VAQUERO_TOF_DAYS_31: dict[float, float] = {2.54: 4.90, 3.13: 5.04}

#: Earth->Moon TOF ceiling, days (SOURCED criterion 3, p. 169-170).
VAQUERO_TOF_CEILING_DAYS = 7.0

#: Conventional Earth equatorial radius, km (convention, not a Vaquero print).
EARTH_RADIUS_KM = 6378.137

#: Earth-center perigee band implied by criterion 1 (DERIVED:
#: ``R_E + 180`` .. ``R_E + 35786`` km -- the same sums `#798`'s note uses).
PERIGEE_BAND_KM: tuple[float, float] = (EARTH_RADIUS_KM + 180.0, EARTH_RADIUS_KM + 35786.0)

#: Half-crossing index for both families' apoapsis-seeded symmetric orbits,
#: measured (not assumed) on the converged seeds: each family's orbit crosses
#: ``y = 0`` six times per period, so the half-period perpendicular
#: re-crossing is the 3rd crossing. Verified by :func:`half_crossing_index`
#: at seed time (the seeding helper asserts agreement rather than trusting
#: this constant).
VAQUERO_HALF_CROSSINGS = 3


def earth_moon_system() -> cr3bp.CR3BPSystem:
    """This project's canonical Earth-Moon CR3BP system (registry mu/l*/t*)."""
    return cr3bp.cr3bp_system("Earth", "Moon")


# ---------------------------------------------------------------------------
# Two-body apoapsis seed (Tisserand-matched eccentricity at the target C).
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class ApoapsisSeed:
    """Two-body ``p:1`` interior resonant ellipse pinned at its ``+x``-axis
    apoapsis, eccentricity fixed by the Tisserand relation at Jacobi ``c``.

    All quantities nondimensional (Earth-Moon distance / synodic units).
    ``ydot_rot`` is the rotating-frame ``ydot`` at apoapsis (negative: the
    spacecraft moves slower than the rotating frame there).
    """

    p: int
    c: float
    a: float
    e: float
    r_apo: float
    r_peri: float
    ydot_rot: float


def tisserand_apoapsis_seed(p: int, c: float) -> ApoapsisSeed:
    """Two-body ``p:1`` resonant-ellipse apoapsis seed at Tisserand-Jacobi ``c``.

    ``a = p**(-2/3)`` (Kepler, moon-period units); ``e`` from
    ``C ~= 1/a + 2*sqrt(a*(1-e**2))`` (planar prograde Tisserand
    approximation of the Jacobi constant, mu -> 0 limit).

    Raises
    ------
    ValueError
        If no prograde ellipse of this ``a`` attains the requested ``c``
        (Tisserand radicand out of ``(0, 1)``).
    """
    if p <= 1:
        raise ValueError(f"interior p:1 resonance needs p >= 2, got p={p}")
    a = float(p) ** (-2.0 / 3.0)
    s = (c - 1.0 / a) / 2.0  # sqrt(a*(1-e^2))
    one_minus_e2 = s * s / a
    if s <= 0.0 or not (0.0 < one_minus_e2 < 1.0):
        raise ValueError(f"no prograde p={p} two-body ellipse at Tisserand C={c}")
    e = math.sqrt(1.0 - one_minus_e2)
    r_apo = a * (1.0 + e)
    r_peri = a * (1.0 - e)
    v_apo = math.sqrt(2.0 / r_apo - 1.0 / a)  # vis-viva, GM=1
    return ApoapsisSeed(p=p, c=c, a=a, e=e, r_apo=r_apo, r_peri=r_peri, ydot_rot=v_apo - r_apo)


def half_crossing_index(
    system: cr3bp.CR3BPSystem, orbit: cp.SymmetricOrbit, *, rtol: float = 1e-12, atol: float = 1e-12
) -> int:
    """1-based index of the half-period crossing among the orbit's ``y = 0``
    crossings -- the ``half_crossings`` value that pins
    :func:`~cyclerfinder.search.cr3bp_continuation.continue_family` to this
    family's branch."""
    state0 = np.array([orbit.x0, 0.0, 0.0, 0.0, orbit.ydot0, 0.0])
    times, _states = cp._xaxis_crossings(
        system, state0, 1.05 * orbit.t_half, with_stm=False, rtol=rtol, atol=atol
    )
    if len(times) == 0:
        raise RuntimeError("no x-axis crossings found on converged orbit")
    return int(np.argmin(np.abs(times - orbit.t_half))) + 1


def seed_vaquero_family(
    system: cr3bp.CR3BPSystem,
    p: int,
    c_seed: float,
    *,
    tol: float = 1e-11,
    rtol: float = 1e-13,
    atol: float = 1e-13,
) -> tuple[cp.SymmetricOrbit, ApoapsisSeed]:
    """Converge one family member at ``c_seed`` from the two-body apoapsis seed.

    Returns the converged symmetric orbit and the two-body geometry it was
    seeded from. Raises ``RuntimeError`` if the corrector fails or lands on
    an unexpected crossing topology (the half-crossing index is MEASURED and
    checked against :data:`VAQUERO_HALF_CROSSINGS`, not trusted).
    """
    geom = tisserand_apoapsis_seed(p, c_seed)
    orbit = cp.correct_symmetric_fixed_jacobi(
        system,
        geom.r_apo,
        c_seed,
        TWO_PI,
        ydot0_sign=-1.0,
        half_crossings=None,
        tol=tol,
        rtol=rtol,
        atol=atol,
    )
    if not orbit.converged:
        raise RuntimeError(
            f"seed corrector failed for p={p} at C={c_seed}: residual {orbit.crossing_residual:.3e}"
        )
    idx = half_crossing_index(system, orbit)
    if idx != VAQUERO_HALF_CROSSINGS:
        raise RuntimeError(
            f"unexpected half-crossing topology for p={p} at C={c_seed}: "
            f"measured index {idx}, expected {VAQUERO_HALF_CROSSINGS}"
        )
    return orbit, geom


# ---------------------------------------------------------------------------
# Per-member criteria report (Vaquero's own four printed criteria).
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class MemberReport:
    """One family member evaluated against Vaquero's own printed criteria.

    Geometry from a dense full-period propagation (100k samples + parabolic
    refinement of each extremum); stability from the Barden half-period
    monodromy (planar ``nu``) plus the decoupled out-of-plane 2x2 monodromy
    block (``nu_z = (M_z[0,0] + M_z[1,1]) / 2`` -- the planar-CR3BP analogue
    of Fig. 4.44's own ``nu_3D`` coloring).
    """

    jacobi: float
    x0: float
    ydot0: float
    period: float
    period_over_2pi: float
    perigee_km: float
    apogee_km: float
    moon_min_km: float
    a_two_body_lu: float  # (r_peri + r_apo)/2 in LU -- resonance check
    tof_earth_moon_days: float
    nu_planar: float
    abs_lambda: float
    nu_out_of_plane: float
    stable_planar: bool
    perigee_in_band: bool  # criterion 1
    tof_within_ceiling: bool  # criterion 3
    moon_min_vs_l1_dist: float  # moon_min / |L1 - Moon| (criterion 2 scale)
    moon_min_vs_l2_dist: float  # moon_min / |L2 - Moon| (criterion 2 scale)


def _refine_extremum_circular(
    t: NDArray[np.float64], r: NDArray[np.float64], i: int, period: float
) -> tuple[float, float]:
    """Parabolic refinement of a sampled extremum at circular index ``i``.

    The samples cover one period on ``[0, T)`` (endpoint excluded), so the
    neighbors wrap; the refined time is reported modulo the period.
    """
    n = len(r)
    im, ip = (i - 1) % n, (i + 1) % n
    r0, r1, r2 = float(r[im]), float(r[i]), float(r[ip])
    denom = r0 - 2.0 * r1 + r2
    if denom == 0.0:
        return float(t[i]), r1
    h = period / n  # uniform circular spacing
    delta = 0.5 * (r0 - r2) / denom
    t_v = (float(t[i]) + delta * h) % period
    r_v = r1 - 0.125 * (r0 - r2) ** 2 / denom
    return t_v, r_v


def _local_minima_circular(r: NDArray[np.float64]) -> NDArray[np.intp]:
    """Indices of strict local minima of a one-period circular sample array.

    The array covers ``[0, T)`` with the ``t = T`` sample excluded (it
    duplicates ``t = 0`` for a periodic orbit), so minima at the seed point
    itself (e.g. the apoapsis/Moon-encounter crossing at ``t = 0``) are
    detected rather than lost as array endpoints.
    """
    return np.where((r < np.roll(r, 1)) & (r < np.roll(r, -1)))[0]


def member_report(
    system: cr3bp.CR3BPSystem,
    state0: NDArray[np.float64],
    period: float,
    *,
    n_samples: int = 100_000,
    rtol: float = 1e-13,
    atol: float = 1e-13,
) -> MemberReport:
    """Evaluate one converged member against Vaquero's own printed criteria."""
    s0 = np.asarray(state0, float)
    # Circular sampling: [0, T) with the duplicate t=T sample excluded, so
    # extrema at the seed crossing itself (t=0) are not lost as endpoints.
    sol = solve_ivp(
        cr3bp.cr3bp_eom,
        (0.0, period),
        s0,
        args=(system.mu,),
        method="DOP853",
        rtol=rtol,
        atol=atol,
        t_eval=np.linspace(0.0, period, n_samples, endpoint=False),
    )
    if not sol.success:
        raise RuntimeError(f"member_report propagation failed: {sol.message}")
    t = sol.t
    x, y = sol.y[0], sol.y[1]
    r1 = np.hypot(x + system.mu, y)
    r2 = np.hypot(x - 1.0 + system.mu, y)

    # Perigee / apogee (refined) -- criterion 1 + resonance bookkeeping.
    i_peri = int(np.argmin(r1))
    i_apo = int(np.argmax(r1))
    _t_peri_g, r_peri = _refine_extremum_circular(t, r1, i_peri, period)
    _t_apo, r_apo_neg = _refine_extremum_circular(t, -r1, i_apo, period)
    r_apo = -r_apo_neg

    # Closest lunar approach (refined) -- criterion 2.
    i_moon = int(np.argmin(r2))
    _t_moon, r_moon = _refine_extremum_circular(t, r2, i_moon, period)

    # Earth->Moon TOF -- criterion 3: shortest periodic-wraparound time from
    # a USABLE perigee (any local r1 minimum inside the LEO-GEO insertion
    # band of criterion 1, or within 5% of the global perigee when none is
    # in-band) to a genuine close lunar approach (any local r2 minimum
    # within 5% of the global). Vaquero's own printed endpoint TOFs pair an
    # in-band perigee with the Moon encounter -- e.g. her 4.90 d at C=2.54
    # (3:1) is the in-band 8,300-km-class perigee leg, NOT the deeper
    # out-of-phase global perigee half a period from the encounter.
    band_lo = PERIGEE_BAND_KM[0] / system.l_km
    band_hi = PERIGEE_BAND_KM[1] / system.l_km
    all_peri = list(_local_minima_circular(r1))
    peri_idx = [i for i in all_peri if band_lo <= r1[i] <= band_hi]
    if not peri_idx:
        peri_idx = [i for i in all_peri if r1[i] <= 1.05 * r1[i_peri]] or [i_peri]
    moon_idx = [i for i in _local_minima_circular(r2) if r2[i] <= 1.05 * r2[i_moon]]
    if not moon_idx:
        moon_idx = [i_moon]

    def _wrap_dt(ta: float, tb: float) -> float:
        d = abs(ta - tb) % period
        return min(d, period - d)

    tof_nd = min(
        _wrap_dt(
            _refine_extremum_circular(t, r2, j, period)[0],
            _refine_extremum_circular(t, r1, i, period)[0],
        )
        for i in peri_idx
        for j in moon_idx
    )
    tof_days = tof_nd * system.t_s / 86400.0

    # Stability: planar Barden nu via the half-period monodromy; out-of-plane
    # nu_z from the decoupled (z, vz) block of the FULL-period STM.
    x0 = float(s0[0])
    ydot0 = float(s0[4])
    jac = cr3bp.jacobi_constant(s0, system.mu)
    orbit = cp.SymmetricOrbit(
        x0=x0,
        ydot0=ydot0,
        jacobi=jac,
        t_half=0.5 * period,
        period=period,
        converged=True,
        crossing_residual=0.0,
        n_iter=0,
    )
    nu, lam = cp.barden_stability(system, orbit, rtol=rtol, atol=atol)
    arc = cr3bp.propagate(system, s0, period, with_stm=True, rtol=rtol, atol=atol)
    assert arc.stm is not None
    mz = arc.stm[np.ix_([2, 5], [2, 5])]
    nu_z = 0.5 * float(mz[0, 0] + mz[1, 1])

    # Moon-relative libration-point distances (criterion 2 scale).
    l1x = crnbp.cr3bp_collinear_point(system.mu, "L1")
    l2x = crnbp.cr3bp_collinear_point(system.mu, "L2")
    moon_x = 1.0 - system.mu
    d_l1 = abs(moon_x - l1x)
    d_l2 = abs(l2x - moon_x)

    perigee_km = r_peri * system.l_km
    moon_min_km = r_moon * system.l_km
    return MemberReport(
        jacobi=jac,
        x0=x0,
        ydot0=ydot0,
        period=period,
        period_over_2pi=period / TWO_PI,
        perigee_km=perigee_km,
        apogee_km=r_apo * system.l_km,
        moon_min_km=moon_min_km,
        a_two_body_lu=0.5 * (r_peri + r_apo),
        tof_earth_moon_days=tof_days,
        nu_planar=float(nu),
        abs_lambda=float(abs(lam)),
        nu_out_of_plane=nu_z,
        stable_planar=abs(float(nu)) < 1.0,
        perigee_in_band=PERIGEE_BAND_KM[0] <= perigee_km <= PERIGEE_BAND_KM[1],
        tof_within_ceiling=tof_days <= VAQUERO_TOF_CEILING_DAYS,
        moon_min_vs_l1_dist=r_moon / d_l1,
        moon_min_vs_l2_dist=r_moon / d_l2,
    )


# ---------------------------------------------------------------------------
# Family reproduction: seed + two-directional continuation across Vaquero's
# own printed C range, via the existing continue_family machinery.
# ---------------------------------------------------------------------------


@dataclass
class FamilyReproduction:
    """Both continuation branches from one seed plus per-member reports."""

    label: str
    p: int
    c_seed: float
    c_range: tuple[float, float]
    seed_geom: ApoapsisSeed
    branch_down: cc.FamilyBranch
    branch_up: cc.FamilyBranch
    members: list[cc.BranchMember]
    reports: list[MemberReport]


#: Continuation endpoint slack: ``continue_family`` accumulates ``C`` by
#: repeated addition, so the exact printed endpoint can round to one ULP
#: outside the bound; a 1e-6 margin keeps the printed endpoints reachable
#: without admitting any member meaningfully outside Vaquero's own range.
_C_BOUND_EPS = 1e-6


def reproduce_vaquero_family(
    system: cr3bp.CR3BPSystem,
    p: int,
    *,
    c_seed: float,
    c_range: tuple[float, float],
    d_jacobi: float = 0.01,
    label: str = "",
    rtol: float = 1e-13,
    atol: float = 1e-13,
) -> FamilyReproduction:
    """Seed at ``c_seed`` and continue to both of Vaquero's printed endpoints.

    Members are merged (down-branch reversed + seed + up-branch, deduplicated
    on Jacobi) and each is evaluated with :func:`member_report`.

    Tolerances default to 1e-13 (tighter than ``continue_family``'s own
    1e-12 default): the low-C end of both families dives to ~7,000-8,500 km
    perigees, where a 1e-12 DOP853 propagation leaks just past the
    gauntlet's 1e-10 Jacobi-conservation gate (measured: the C=2.54 3:1
    member converges at residual 7e-13 but fails that one gate at 1e-12).
    """
    lo, hi = min(c_range), max(c_range)
    seed, geom = seed_vaquero_family(system, p, c_seed, rtol=rtol, atol=atol)
    n_down = math.ceil((c_seed - lo) / d_jacobi) + 2
    n_up = math.ceil((hi - c_seed) / d_jacobi) + 2
    branch_down = cc.continue_family(
        system,
        seed,
        direction=-1,
        d_jacobi=d_jacobi,
        n_steps=n_down,
        min_jacobi=lo - _C_BOUND_EPS,
        max_jacobi=hi + _C_BOUND_EPS,
        half_crossings=VAQUERO_HALF_CROSSINGS,
        ydot0_sign=-1.0,
        seed_label=label or f"vaquero_{p}to1",
        rtol=rtol,
        atol=atol,
    )
    branch_up = cc.continue_family(
        system,
        seed,
        direction=+1,
        d_jacobi=d_jacobi,
        n_steps=n_up,
        min_jacobi=lo - _C_BOUND_EPS,
        max_jacobi=hi + _C_BOUND_EPS,
        half_crossings=VAQUERO_HALF_CROSSINGS,
        ydot0_sign=-1.0,
        seed_label=label or f"vaquero_{p}to1",
        rtol=rtol,
        atol=atol,
    )
    members: list[cc.BranchMember] = []
    seen: set[float] = set()
    ordered = list(reversed(branch_down.members)) + branch_up.members
    for m in ordered:
        key = round(m.jacobi, 9)
        if key in seen:
            continue
        seen.add(key)
        members.append(m)
    members.sort(key=lambda m: m.jacobi)
    reports = [member_report(system, m.state0, m.period, rtol=rtol, atol=atol) for m in members]
    return FamilyReproduction(
        label=label or f"vaquero_{p}to1",
        p=p,
        c_seed=c_seed,
        c_range=(lo, hi),
        seed_geom=geom,
        branch_down=branch_down,
        branch_up=branch_up,
        members=members,
        reports=reports,
    )


__all__ = [
    "EARTH_RADIUS_KM",
    "PERIGEE_BAND_KM",
    "TWO_PI",
    "VAQUERO_C_RANGE_21",
    "VAQUERO_C_RANGE_31",
    "VAQUERO_HALF_CROSSINGS",
    "VAQUERO_TOF_CEILING_DAYS",
    "VAQUERO_TOF_DAYS_21",
    "VAQUERO_TOF_DAYS_31",
    "ApoapsisSeed",
    "FamilyReproduction",
    "MemberReport",
    "earth_moon_system",
    "half_crossing_index",
    "member_report",
    "reproduce_vaquero_family",
    "seed_vaquero_family",
    "tisserand_apoapsis_seed",
]
