"""Second-species seeds continued in the mass ratio to the Earth-Moon value (#899 step 2).

Reproduction of the method of Casoliva, Mondelo, Villac, Mease, Barrabes & Olle (JGCD 33(5),
2010, section IV.C; AIAA 2008-6434, section IV.C): periodic orbits of the planar circular
restricted problem at a small mass ratio (1e-6), seeded from second-species solutions, are
continued in the mass ratio to the Earth-Moon value by a three-step strategy (mass at fixed
period; Jacobi constant at fixed mass to raise the periselene when the family heads for the
Moon; mass again at fixed Jacobi constant), and the resonant members (period 2 pi q) found at
the Earth-Moon mass are compared with the printed Table 3 rows of the 2010 paper. Digests:
``docs/notes/2026-10-04-digest-casoliva-2008-aiaa-families-cycler-trajectories-seeds.md`` and
``docs/notes/2026-07-27-725-casoliva-earth-moon-cycler-families-digest.md``.

Seeds: the ten corrected second-species orbits of the 2008 paper's Table 2 at mu = 1e-6
(:data:`TABLE2_SEEDS`, transcribed from the 2008 digest), which are the seeds of the method as
published. The p-q labels follow Barrabes & Gomez (p spacecraft revolutions per q lunar
revolutions, Kepler a = (q/p)^(2/3), period near 2 pi q), as the digest establishes.

Frames: both Casoliva tables put the Earth at (+mu, 0) and the Moon at (mu - 1, 0); the project
frame (:mod:`cyclerfinder.core.cr3bp`) has the Earth at (-mu, 0) and the Moon at (1 - mu, 0).
The two are related by the rotation by pi, (x, y, u, v) -> (-x, -y, -u, -v) (validated in
:mod:`cyclerfinder.search.earth_moon_resonant_families`). The Jacobi constant convention is the
same in both (2 Omega - V^2, no mu(1 - mu) term).

Integration: every segment is propagated with the planar moon-centred Levi-Civita propagator of
:mod:`cyclerfinder.search.second_species_lc` (compiled; checked against the #928 KS propagator
of :mod:`cyclerfinder.core.cr3bp_ks` to 1e-13 in state over a third of a period and to the
noise of the unstable orbits over a full one), regular at the Moon, so passes of 1e-4 lunar
distances at mu = 1e-6 need no special handling; the Earth is not regularised (the perigees of
the catalogued rows are at least 0.017). The planar 4 x 4 transition matrix gives the in-plane
index and the vertical (z, vz) block the out-of-plane index (Casoliva prints one or the other;
see :class:`cyclerfinder.search.earth_moon_resonant_families.StabilityIndex`).

Corrector: general (asymmetric-capable) multiple shooting on N segments of equal duration
T / N, unknowns the N planar node states and T, equations the N continuity conditions, the
section condition y = 0 at the first node, and one or two parameter conditions (fixed period,
fixed Jacobi constant). The one redundant continuity condition (the Jacobi integral) makes the
system consistent but overdetermined, so each Newton step is a least-squares (minimum-norm)
step. Continuation in the mass is natural-parameter in log(mu) with a secant predictor and step
halving; continuation in the Jacobi constant at fixed mass is pseudo-arclength (tangent from the
SVD of the Jacobian).

Stability: two scales are logged and must not be mixed. ``k_par`` and ``k_perp`` are
Casoliva's Eq. 8 indices lambda + 1/lambda (critical at |k| = 2) from the full-period planar
and vertical blocks; the :mod:`cyclerfinder.core.floquet_classes` index is the normalised
k = -b/2 (critical at |k| = 1).

Validity numbers logged at every member (Gomez & Olle, Guillaume, Perko digests): the
periselene r_p, nu_eff = ln(r_p)/ln(mu), mu |ln mu| / V^3 and r_p / sqrt(mu), with
V = sqrt(3 - C) the encounter speed of the mu -> 0 problem.
"""

from __future__ import annotations

import math
from collections.abc import Callable
from dataclasses import dataclass, field

import numpy as np
from numpy.typing import NDArray

from cyclerfinder.core.floquet_classes import classify_planar_monodromy
from cyclerfinder.core.satellites import SATELLITES
from cyclerfinder.search.second_species_lc import (
    LCPropagationError,
    jacobi_planar,
    planar_accel,
    propagate_lc,
)

FloatArray = NDArray[np.float64]

PLANAR = (0, 1, 3, 4)

#: The mass ratio printed in Casoliva et al. 2010 (p.1629), the target of their continuation.
CASOLIVA_MU_2010 = 0.0121529529

#: The mass ratio of the 2008 Table 2 seeds (the printed states reproduce the printed Jacobi
#: constants to 1e-16 at this value only; 2008 digest section 6.1).
SEED_MU = 1.0e-6

#: Moon mean radius in lunar distances (registry radius over registry semi-major axis).
MOON_RADIUS_ND = SATELLITES["Moon"].radius_eq_km / SATELLITES["Moon"].sma_km

#: Earth equatorial radius (WGS84) in lunar distances.
EARTH_RADIUS_ND = 6378.137 / SATELLITES["Moon"].sma_km


# ---------------------------------------------------------------------------------------------
# Seeds: Casoliva et al. 2008 Table 2 (mu = 1e-6), Casoliva frame


@dataclass(frozen=True)
class Table2Seed:
    """One printed row of the 2008 Table 2 (Casoliva frame: Moon at (mu - 1, 0))."""

    designation: str
    p: int  # spacecraft revolutions (Barrabes & Gomez labelling)
    q: int  # lunar revolutions; period near 2 pi q
    c_j: float
    period: float
    x_i: float
    u_i: float
    v_i: float
    k: float  # printed stability index, lambda + 1/lambda (planar block)

    def state_project(self) -> FloatArray:
        """Planar project-frame state (x, y, vx, vy) at the printed y = 0 crossing."""
        return np.array([-self.x_i, 0.0, -self.u_i, -self.v_i])


#: 2008 Table 2, transcribed digit by digit in the 2008 digest (section 4); y_i = 0 for all.
TABLE2_SEEDS: tuple[Table2Seed, ...] = (
    Table2Seed("12a", 1, 2, -0.4048949508278787, 12.5729004699819580,
               -0.9997842236429277, -1.0475686168407430, 1.5221048236500363, 994.8214),
    Table2Seed("21a", 2, 1, 0.3044238301466371, 6.2807379868905375,
               -0.9994542367695188, 0.0000000000159674, 1.6429377286480178, 2.0214),
    Table2Seed("23a", 2, 3, -1.4624706218555543, 18.8645742008117736,
               -0.9997040086087932, -0.0000000000029493, 2.1140593042383320, 1.8405),
    Table2Seed("23b", 2, 3, -0.3856270265962789, 18.8497995179817757,
               -0.9953809710199844, -0.9554688344071062, 1.5726409617865593, -0.7753),
    Table2Seed("32a", 3, 2, -0.3403221450450835, 12.5363376944500722,
               -0.9999035023356472, 0.0000000000131068, 1.8333742370247279, 6.3818),
    Table2Seed("32b", 3, 2, 2.0635340336761394, 12.5660196280208911,
               -1.0188148478462549, -0.7183851038146349, 0.6492612209948345, 1.5681),
    Table2Seed("52a", 5, 2, 1.0461882704974470, 12.5651492405106922,
               -0.9994423797258251, 0.0000000000459148, 1.3990717541095201, 2.0366),
    Table2Seed("54a", 5, 4, -0.5902501452788234, 25.1304852528305673,
               -0.9988849982363450, -0.2709574260065351, 1.8758004354491455, -1.3192),
    Table2Seed("54b", 5, 4, -0.6598717597930506, 25.1321450585584110,
               -0.9905419593393083, -0.1198961389475177, 1.9094434196183308, 1.9943),
    Table2Seed("73a", 7, 3, 0.8957501590757784, 18.8492803402344329,
               -0.9954265899784440, -0.2486030886355384, 1.4293154529931373, 1.8799),
)  # fmt: skip


def table2_seed(designation: str) -> Table2Seed:
    for s in TABLE2_SEEDS:
        if s.designation == designation:
            return s
    raise ValueError(f"unknown 2008 Table 2 designation {designation!r}")


# ---------------------------------------------------------------------------------------------
# Propagation (planar, moon-centred Levi-Civita; see second_species_lc)


@dataclass(frozen=True)
class Segment:
    end: FloatArray  # planar (x, y, vx, vy)
    stm4: FloatArray | None  # planar 4 x 4 fixed-time transition matrix
    stm_z: FloatArray | None  # vertical (z, vz) 2 x 2 block
    r_min: float  # closest approach to the Moon on the segment
    earth_min: float  # closest approach to the Earth on the segment


def planar_eom(state4: FloatArray, mu: float) -> FloatArray:
    return planar_accel(state4, mu)


def jacobi4(state4: FloatArray, mu: float) -> float:
    return jacobi_planar(state4, mu)


def propagate_segment(
    mu: float,
    state4: FloatArray,
    dt: float,
    *,
    with_stm: bool = True,
    rtol: float = 1e-12,
    atol: float = 1e-14,
) -> Segment:
    """Propagate a planar state for ``dt`` (regular at the Moon)."""
    arc = propagate_lc(mu, state4, dt, with_stm=with_stm, rtol=rtol, atol=atol)
    return Segment(
        end=arc.end, stm4=arc.stm4, stm_z=arc.stm_z, r_min=arc.r_min, earth_min=arc.earth_min
    )


# ---------------------------------------------------------------------------------------------
# Multiple-shooting corrector


@dataclass
class MSOrbit:
    """A periodic orbit as N multiple-shooting nodes (planar states) and a period."""

    mu: float
    nodes: FloatArray  # (N, 4); node 0 lies on y = 0
    period: float
    residual: float = math.inf
    iterations: int = 0
    converged: bool = False
    segments: list[Segment] = field(default_factory=list, repr=False)

    @property
    def n(self) -> int:
        return int(self.nodes.shape[0])

    @property
    def jacobi(self) -> float:
        return jacobi4(self.nodes[0], self.mu)

    def unknowns(self) -> FloatArray:
        return np.concatenate([self.nodes.ravel(), [self.period]])

    def monodromy4(self) -> FloatArray:
        """Planar one-period monodromy from node 0."""
        m = np.eye(4)
        for seg in self.segments:
            assert seg.stm4 is not None
            m = seg.stm4 @ m
        return m

    def monodromy_z(self) -> FloatArray:
        """Vertical (z, vz) one-period monodromy from node 0."""
        m = np.eye(2)
        for seg in self.segments:
            assert seg.stm_z is not None
            m = seg.stm_z @ m
        return m


class CorrectionError(RuntimeError):
    pass


_STEP_FAILURES = (
    CorrectionError,
    LCPropagationError,
    np.linalg.LinAlgError,
    ZeroDivisionError,
    ValueError,
)


def nodes_from_state(mu: float, state4: FloatArray, period: float, n: int) -> FloatArray:
    """N equally spaced nodes along the orbit of ``state4`` (node 0 = ``state4``)."""
    nodes = np.empty((n, 4))
    nodes[0] = state4
    for i in range(1, n):
        nodes[i] = propagate_segment(mu, nodes[i - 1], period / n, with_stm=False).end
    return nodes


def renode(orbit: MSOrbit, n: int) -> MSOrbit:
    """The same orbit with ``n`` equally spaced nodes from node 0 (uncorrected)."""
    return MSOrbit(
        orbit.mu, nodes_from_state(orbit.mu, orbit.nodes[0], orbit.period, n), orbit.period
    )


def _evaluate(
    mu: float, nodes: FloatArray, period: float
) -> tuple[FloatArray, FloatArray, list[Segment]]:
    """Continuity residuals (4N) and their Jacobian with respect to (nodes, T)."""
    n = nodes.shape[0]
    dt = period / n
    res = np.empty(4 * n)
    jac = np.zeros((4 * n, 4 * n + 1))
    segs = []
    for i in range(n):
        seg = propagate_segment(mu, nodes[i], dt, with_stm=True)
        segs.append(seg)
        j = (i + 1) % n
        res[4 * i : 4 * i + 4] = seg.end - nodes[j]
        assert seg.stm4 is not None
        jac[4 * i : 4 * i + 4, 4 * i : 4 * i + 4] = seg.stm4
        jac[4 * i : 4 * i + 4, 4 * j : 4 * j + 4] -= np.eye(4)
        jac[4 * i : 4 * i + 4, 4 * n] = planar_eom(seg.end, mu) / n
    return res, jac, segs


def _jacobi_grad(state4: FloatArray, mu: float) -> FloatArray:
    x, y, vx, vy = state4
    r1 = math.hypot(x + mu, y)
    r2 = math.hypot(x - 1.0 + mu, y)
    ox = x - (1.0 - mu) * (x + mu) / r1**3 - mu * (x - 1.0 + mu) / r2**3
    oy = y - (1.0 - mu) * y / r1**3 - mu * y / r2**3
    return np.array([2.0 * ox, 2.0 * oy, -2.0 * vx, -2.0 * vy])


def correct(
    orbit: MSOrbit,
    *,
    fix_period: float | None = None,
    fix_jacobi: float | None = None,
    arclength: tuple[FloatArray, FloatArray, float] | None = None,
    tol: float = 1e-10,
    accept_tol: float = 1e-8,
    max_iter: int = 12,
) -> MSOrbit:
    """Gauss-Newton multiple-shooting correction at fixed ``orbit.mu``.

    Conditions: continuity, y = 0 at node 0, and any of ``fix_period``, ``fix_jacobi`` and
    ``arclength = (u_prev, tangent, ds)`` (pseudo-arclength: tangent . (U - u_prev) = ds,
    with U = (nodes, T)). With neither fix nor arclength the system is underdetermined (the
    family at fixed mu) and the minimum-norm step is taken. Raises :class:`CorrectionError`
    when the residual does not fall below ``tol``; a residual below ``accept_tol`` at which
    Newton stalls (the integration-noise floor) is accepted and reported in ``residual``.
    """
    mu = orbit.mu
    nodes = orbit.nodes.copy()
    period = float(orbit.period)
    n = nodes.shape[0]
    nu = 4 * n + 1
    last_norm = math.inf
    for it in range(max_iter + 1):
        res, jac, segs = _evaluate(mu, nodes, period)
        rows_r = [res, np.array([nodes[0, 1]])]
        rows_j = [jac, np.eye(nu)[1:2]]
        if fix_period is not None:
            rows_r.append(np.array([period - fix_period]))
            rows_j.append(np.eye(nu)[nu - 1 : nu])
        if fix_jacobi is not None:
            g = np.zeros((1, nu))
            g[0, :4] = _jacobi_grad(nodes[0], mu)
            rows_r.append(np.array([jacobi4(nodes[0], mu) - fix_jacobi]))
            rows_j.append(g)
        if arclength is not None:
            u_prev, tan, ds = arclength
            u = np.concatenate([nodes.ravel(), [period]])
            rows_r.append(np.array([float(tan @ (u - u_prev)) - ds]))
            rows_j.append(tan.reshape(1, nu))
        f = np.concatenate(rows_r)
        jm = np.vstack(rows_j)
        norm = float(np.max(np.abs(f)))
        if norm < tol:
            return MSOrbit(
                mu=mu,
                nodes=nodes,
                period=period,
                residual=norm,
                iterations=it,
                converged=True,
                segments=segs,
            )
        stalled = it >= 1 and norm > 0.5 * last_norm
        if stalled and norm < accept_tol:
            # integration-noise floor reached (the floor grows with the orbit's instability)
            return MSOrbit(
                mu=mu,
                nodes=nodes,
                period=period,
                residual=norm,
                iterations=it,
                converged=True,
                segments=segs,
            )
        if it == max_iter or not math.isfinite(norm) or (it >= 4 and norm > 0.9 * last_norm):
            raise CorrectionError(f"multiple shooting did not converge (residual {norm:.3e})")
        last_norm = norm
        step = np.linalg.lstsq(jm, -f, rcond=None)[0]
        nodes = nodes + step[:-1].reshape(n, 4)
        period = period + float(step[-1])
        if not period > 0.0:
            raise CorrectionError("multiple shooting: period became non-positive")
    raise CorrectionError("unreachable")


def family_tangent(orbit: MSOrbit) -> FloatArray:
    """Unit tangent of the family at fixed mu (null vector of continuity + section)."""
    _res, jac, _ = _evaluate(orbit.mu, orbit.nodes, orbit.period)
    nu = jac.shape[1]
    jm = np.vstack([jac, np.eye(nu)[1:2]])
    # one continuity row is redundant (Jacobi integral): the null space is the
    # right singular vector of the smallest singular value
    _, _, vt = np.linalg.svd(jm)
    return np.asarray(vt[-1], dtype=np.float64)


# ---------------------------------------------------------------------------------------------
# Diagnostics


@dataclass(frozen=True)
class Diagnostics:
    mu: float
    jacobi: float
    period: float
    periselene: float  # lunar distances, from the Moon centre
    perigee: float  # lunar distances, from the Earth centre
    k_par: float  # Casoliva Eq. 8, planar block: trace - 2 (critical |k| = 2)
    k_perp: float  # Casoliva Eq. 8, vertical (z, vz) block trace (critical |k| = 2)
    floquet_regime: str  # floquet_classes on reduced_monodromy (normalised k, critical |k| = 1)
    floquet_k: tuple[float, float]  # normalised indices
    floquet_delta: float
    nu_eff: float  # ln(r_p) / ln(mu)
    validity_speed: float  # mu |ln mu| / V^3, V = sqrt(3 - C)
    validity_rp: float  # r_p / sqrt(mu)
    impact_moon: bool  # periselene below the Moon's radius
    impact_earth: bool  # perigee below the Earth's radius


def reduced_monodromy(k_par: float, k_perp: float) -> FloatArray:
    """A 4 x 4 matrix whose reciprocal characteristic polynomial carries the two nontrivial
    pairs of a planar periodic orbit: the in-plane pair (Casoliva index ``k_par``) and the
    vertical pair (``k_perp``), as companion blocks of lambda^2 - k lambda + 1.

    The planar 4 x 4 block itself always carries the trivial pair (lambda = 1, 1), which
    :mod:`cyclerfinder.core.floquet_classes` would report as "boundary" for every orbit; the
    6 x 6 monodromy reduced by that pair is the four-dimensional problem those classes
    describe. Normalised indices there are k_par / 2 and k_perp / 2.
    """
    m = np.zeros((4, 4))
    m[0, 1], m[1, 0], m[1, 1] = -1.0, 1.0, k_par
    m[2, 3], m[3, 2], m[3, 3] = -1.0, 1.0, k_perp
    return m


def diagnose(orbit: MSOrbit) -> Diagnostics:
    """Periselene, perigee, stability indices and validity numbers of a corrected orbit."""
    mu = orbit.mu
    c = orbit.jacobi
    rp = min(seg.r_min for seg in orbit.segments)
    perigee = min(seg.earth_min for seg in orbit.segments)
    m4 = orbit.monodromy4()
    mz = orbit.monodromy_z()
    k_par = float(np.trace(m4)) - 2.0
    k_perp = float(np.trace(mz))
    fc = classify_planar_monodromy(reduced_monodromy(k_par, k_perp))
    v = math.sqrt(max(3.0 - c, 1e-300))
    return Diagnostics(
        mu=mu,
        jacobi=c,
        period=orbit.period,
        periselene=rp,
        perigee=perigee,
        k_par=k_par,
        k_perp=k_perp,
        floquet_regime=fc.regime,
        floquet_k=(float(fc.k1.real), float(fc.k2.real)),
        floquet_delta=float(fc.delta),
        nu_eff=math.log(rp) / math.log(mu),
        validity_speed=mu * abs(math.log(mu)) / v**3,
        validity_rp=rp / math.sqrt(mu),
        impact_moon=rp < MOON_RADIUS_ND,
        impact_earth=perigee < EARTH_RADIUS_ND,
    )


# ---------------------------------------------------------------------------------------------
# Reference orbits: the printed 2010 Table 3 states corrected in place


def table3_reference(designation: str, mu: float = CASOLIVA_MU_2010, n: int = 4) -> MSOrbit:
    """The printed 2010 Table 3 crossing of ``designation`` corrected at fixed T = 2 pi q.

    At the paper's mass ratio the printed states close with corrections of about 1e-10 in
    every component (the printed digits), while at the registry mass ratio the corrections are
    about 1e-6; so the paper's mu is the target of the continuation (#899 step 2).
    """
    from cyclerfinder.search.earth_moon_resonant_families import table3_row

    row = table3_row(designation)
    st = np.array([-row.x_i, 0.0, -row.u_i, -row.v_i])
    period = 2.0 * math.pi * row.q
    nodes = nodes_from_state(mu, st, period, n)
    return correct(MSOrbit(mu, nodes, period), fix_period=period)


def casoliva_k(designation: str, diag: Diagnostics) -> float:
    """The index Casoliva prints for a 2010 Table 3 row: the larger in modulus of k_par and
    k_perp, except the rows of the #801 override set, where the printed value is k_perp
    (:class:`cyclerfinder.search.earth_moon_resonant_families.StabilityIndex`)."""
    from cyclerfinder.search.earth_moon_resonant_families import StabilityIndex

    return StabilityIndex(
        k_par=diag.k_par,
        k_perp=diag.k_perp,
        k_eig=math.nan,
        lam=complex(math.nan),
        agree=True,
        designation=designation,
    ).k_signed


# ---------------------------------------------------------------------------------------------
# Continuation


@dataclass(frozen=True)
class Member:
    """One accepted continuation member with its diagnostics."""

    leg: str  # "mu@T", "C@mu", "mu@C", "walk@mu"
    orbit: MSOrbit
    diag: Diagnostics
    step: float  # step used to reach this member (log mu, or arclength)


StopFn = Callable[[Member], str | None]
LogFn = Callable[[Member], None]


def _max_node_change(a: MSOrbit, b: MSOrbit) -> float:
    return float(np.max(np.abs(a.nodes - b.nodes)))


@dataclass
class LegResult:
    members: list[Member]
    reason: str  # "target", "stop:<why>", "min_step", "max_steps"


def continue_mu(
    orbit: MSOrbit,
    mu_target: float,
    *,
    fix: str,
    leg: str = "",
    dlog0: float = 0.1,
    dlog_min: float = 1e-4,
    dlog_max: float = 0.4,
    max_jump: float = 0.05,
    max_steps: int = 2000,
    stop: StopFn | None = None,
    log: LogFn | None = None,
) -> LegResult:
    """Natural-parameter continuation in log(mu) at fixed period (``fix="period"``) or fixed
    Jacobi constant (``fix="jacobi"``) from ``orbit`` to ``mu_target`` (secant predictor, step
    halving on a failed correction or on a node jump above ``max_jump``)."""
    if fix not in ("period", "jacobi"):
        raise ValueError("fix must be 'period' or 'jacobi'")
    value = orbit.period if fix == "period" else orbit.jacobi
    leg = leg or ("mu@T" if fix == "period" else "mu@C")
    members: list[Member] = []
    cur = orbit
    prev: MSOrbit | None = None
    lt = math.log(mu_target)
    dlog = math.copysign(dlog0, lt - math.log(cur.mu))
    for _ in range(max_steps):
        lc = math.log(cur.mu)
        if abs(lt - lc) < 1e-14:
            return LegResult(members, "target")
        step = dlog if abs(dlog) < abs(lt - lc) else lt - lc
        while True:
            ln = lc + step
            mu_new = mu_target if ln == lt else math.exp(ln)
            nodes = cur.nodes.copy()
            period = cur.period
            if prev is not None:
                lp = math.log(prev.mu)
                f = step / (lc - lp)
                nodes = nodes + f * (cur.nodes - prev.nodes)
                period = period + f * (cur.period - prev.period)
            if fix == "period":
                period = value
            trial = MSOrbit(mu_new, nodes, period)
            try:
                new = correct(
                    trial,
                    fix_period=value if fix == "period" else None,
                    fix_jacobi=value if fix == "jacobi" else None,
                    max_iter=8,
                )
                ok = _max_node_change(new, cur) <= max_jump
            except _STEP_FAILURES:
                ok = False
            if ok:
                break
            step *= 0.5
            if abs(step) < dlog_min:
                return LegResult(members, "min_step")
        mem = Member(leg, new, diagnose(new), step)
        members.append(mem)
        if log is not None:
            log(mem)
        prev, cur = cur, new
        if new.iterations <= 4 and abs(step) == abs(dlog):
            dlog = math.copysign(min(abs(dlog) * 1.5, dlog_max), dlog)
        elif abs(step) < abs(dlog):
            dlog = step
        if stop is not None:
            why = stop(mem)
            if why:
                return LegResult(members, "stop:" + why)
    return LegResult(members, "max_steps")


def _jacobi_rate(orbit: MSOrbit, tan: FloatArray) -> float:
    return float(_jacobi_grad(orbit.nodes[0], orbit.mu) @ tan[:4])


def continue_jacobi(
    orbit: MSOrbit,
    direction: float,
    *,
    leg: str = "C@mu",
    ds0: float = 0.01,
    ds_min: float = 1e-6,
    ds_max: float = 0.1,
    max_jump: float = 0.05,
    max_steps: int = 500,
    stop: StopFn | None = None,
    log: LogFn | None = None,
) -> LegResult:
    """Pseudo-arclength continuation of the family through ``orbit`` at fixed mu, starting in
    the direction of increasing (``direction > 0``) or decreasing Jacobi constant; the
    direction is then kept by tangent orientation (folds in C are passed)."""
    members: list[Member] = []
    cur = orbit
    tan = family_tangent(cur)
    if _jacobi_rate(cur, tan) * direction < 0.0:
        tan = -tan
    ds = ds0
    for _ in range(max_steps):
        while True:
            u0 = cur.unknowns()
            up = u0 + ds * tan
            n = cur.n
            trial = MSOrbit(cur.mu, up[:-1].reshape(n, 4), float(up[-1]))
            try:
                new = correct(trial, arclength=(u0, tan, ds), max_iter=8)
                ok = _max_node_change(new, cur) <= max_jump
            except _STEP_FAILURES:
                ok = False
            if ok:
                break
            ds *= 0.5
            if ds < ds_min:
                return LegResult(members, "min_step")
        mem = Member(leg, new, diagnose(new), ds)
        members.append(mem)
        if log is not None:
            log(mem)
        new_tan = family_tangent(new)
        if float(new_tan @ tan) < 0.0:
            new_tan = -new_tan
        cur, tan = new, new_tan
        if new.iterations <= 4:
            ds = min(ds * 1.5, ds_max)
        if stop is not None:
            why = stop(mem)
            if why:
                return LegResult(members, "stop:" + why)
    return LegResult(members, "max_steps")


# ---------------------------------------------------------------------------------------------
# Section crossings and comparison with printed rows


def y_crossings(orbit: MSOrbit, *, steps_per_segment: int = 60) -> list[tuple[float, FloatArray]]:
    """All crossings of y = 0 over one period, as (time from node 0, planar state)."""
    mu = orbit.mu
    out: list[tuple[float, FloatArray]] = []
    dt = orbit.period / orbit.n
    h = dt / steps_per_segment
    for i in range(orbit.n):
        s = orbit.nodes[i]
        t = i * dt
        for _ in range(steps_per_segment):
            s2 = propagate_segment(mu, s, h, with_stm=False).end
            if s[1] == 0.0 and not (i == 0 and t == 0.0 and out):
                if not out or abs(out[-1][0] - t) > 1e-9:
                    out.append((t, s.copy()))
            elif s[1] * s2[1] < 0.0:
                # Newton on y(t) = 0 from s with ydot = vy
                tau = -s[1] / s[3] if s[3] != 0.0 else 0.5 * h
                tau = min(max(tau, 0.0), h)
                sc = s
                for _ in range(30):
                    sc = propagate_segment(mu, s, tau, with_stm=False).end if tau > 0 else s
                    if abs(sc[1]) < 1e-14:
                        break
                    tau -= sc[1] / sc[3]
                    tau = min(max(tau, 0.0), h)
                out.append((t + tau, sc.copy()))
            s, t = s2, t + h
    return out


def mirror(state4: FloatArray) -> FloatArray:
    """The time-reversal mirror of a planar state about the x axis."""
    return np.array([state4[0], -state4[1], -state4[2], state4[3]])


@dataclass(frozen=True)
class RowMatch:
    designation: str
    state_distance: float  # min over crossings (and mirrors) of max-abs state difference
    mirrored: bool
    crossing_time: float
    jacobi_rel_err: float  # against the printed C_J
    period_rel_err: float  # against the printed T
    x_rel_err: float
    u_abs_err: float
    v_rel_err: float
    k_printed: float
    k_par: float
    k_perp: float


def match_row(orbit: MSOrbit, designation: str) -> RowMatch:
    """Compare a continued orbit with a printed 2010 Table 3 row at its printed crossing."""
    from cyclerfinder.search.earth_moon_resonant_families import table3_row

    row = table3_row(designation)
    target = np.array([-row.x_i, 0.0, -row.u_i, -row.v_i])
    best: tuple[float, bool, float, FloatArray] = (math.inf, False, math.nan, target)
    for t, st in y_crossings(orbit):
        for mir, cand in ((False, st), (True, mirror(st))):
            dist = float(np.max(np.abs(cand - target)))
            if dist < best[0]:
                best = (dist, mir, t, cand)
    dg = diagnose(orbit) if orbit.segments else None
    _, _, _, st = best
    return RowMatch(
        designation=designation,
        state_distance=best[0],
        mirrored=best[1],
        crossing_time=best[2],
        jacobi_rel_err=abs(orbit.jacobi - row.c_j) / abs(row.c_j),
        period_rel_err=abs(orbit.period - row.period) / abs(row.period),
        x_rel_err=abs(st[0] - target[0]) / abs(target[0]),
        u_abs_err=abs(st[2] - target[2]),
        v_rel_err=abs(st[3] - target[3]) / abs(target[3]),
        k_printed=row.k,
        k_par=dg.k_par if dg else math.nan,
        k_perp=dg.k_perp if dg else math.nan,
    )
