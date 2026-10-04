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
from collections.abc import Callable, Sequence
from dataclasses import dataclass, field
from itertools import pairwise

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
    nodes: FloatArray  # (N, 4)
    period: float
    residual: float = math.inf
    iterations: int = 0
    converged: bool = False
    segments: list[Segment] = field(default_factory=list, repr=False)
    fractions: FloatArray | None = None  # segment durations / T (None: equal)
    phase_mode: str = "y0"  # "y0": node 0 on y = 0; "flow": node 0 on the plane through the
    # guess's node 0 orthogonal to the flow there (Poincare phase condition)

    @property
    def n(self) -> int:
        return int(self.nodes.shape[0])

    @property
    def fracs(self) -> FloatArray:
        if self.fractions is None:
            return np.full(self.n, 1.0 / self.n)
        return self.fractions

    def like(self, mu: float, nodes: FloatArray, period: float) -> MSOrbit:
        """An uncorrected orbit with this orbit's segment fractions and phase mode."""
        return MSOrbit(mu, nodes, period, fractions=self.fractions, phase_mode=self.phase_mode)

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
    mu: float, nodes: FloatArray, period: float, fracs: FloatArray | None = None
) -> tuple[FloatArray, FloatArray, list[Segment]]:
    """Continuity residuals (4N) and their Jacobian with respect to (nodes, T)."""
    n = nodes.shape[0]
    fr = np.full(n, 1.0 / n) if fracs is None else fracs
    res = np.empty(4 * n)
    jac = np.zeros((4 * n, 4 * n + 1))
    segs = []
    for i in range(n):
        seg = propagate_segment(mu, nodes[i], period * fr[i], with_stm=True)
        segs.append(seg)
        j = (i + 1) % n
        res[4 * i : 4 * i + 4] = seg.end - nodes[j]
        assert seg.stm4 is not None
        jac[4 * i : 4 * i + 4, 4 * i : 4 * i + 4] = seg.stm4
        jac[4 * i : 4 * i + 4, 4 * j : 4 * j + 4] -= np.eye(4)
        jac[4 * i : 4 * i + 4, 4 * n] = planar_eom(seg.end, mu) * fr[i]
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
    phase: tuple[FloatArray, FloatArray] | None = None,
    damped: bool = False,
    tol: float = 1e-10,
    accept_tol: float = 1e-8,
    max_iter: int = 12,
) -> MSOrbit:
    """Gauss-Newton multiple-shooting correction at fixed ``orbit.mu``.

    Conditions: continuity, y = 0 at node 0 (or, with ``phase = (point, direction)``,
    direction . (node 0 - point) = 0), and any of ``fix_period``, ``fix_jacobi`` and
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
    fr = orbit.fracs
    if phase is None and orbit.phase_mode == "flow":
        f0 = planar_eom(orbit.nodes[0], mu)
        phase = (orbit.nodes[0].copy(), f0 / np.linalg.norm(f0))

    def _done(norm: float, it: int, segs: list[Segment]) -> MSOrbit:
        return MSOrbit(
            mu=mu,
            nodes=nodes,
            period=period,
            residual=norm,
            iterations=it,
            converged=True,
            segments=segs,
            fractions=orbit.fractions,
            phase_mode=orbit.phase_mode,
        )

    def _residual_norm(tn: FloatArray, tp: float) -> float:
        worst = 0.0
        for i in range(n):
            end = propagate_segment(mu, tn[i], tp * fr[i], with_stm=False).end
            worst = max(worst, float(np.max(np.abs(end - tn[(i + 1) % n]))))
        if phase is None:
            worst = max(worst, abs(float(tn[0, 1])))
        else:
            worst = max(worst, abs(float(phase[1] @ (tn[0] - phase[0]))))
        if fix_period is not None:
            worst = max(worst, abs(tp - fix_period))
        if fix_jacobi is not None:
            worst = max(worst, abs(jacobi4(tn[0], mu) - fix_jacobi))
        if arclength is not None:
            u = np.concatenate([tn.ravel(), [tp]])
            worst = max(worst, abs(float(arclength[1] @ (u - arclength[0])) - arclength[2]))
        return worst

    for it in range(max_iter + 1):
        res, jac, segs = _evaluate(mu, nodes, period, fr)
        if phase is None:
            rows_r = [res, np.array([nodes[0, 1]])]
            rows_j = [jac, np.eye(nu)[1:2]]
        else:
            ref_pt, ref_dir = phase
            ph = np.zeros((1, nu))
            ph[0, :4] = ref_dir
            rows_r = [res, np.array([float(ref_dir @ (nodes[0] - ref_pt))])]
            rows_j = [jac, ph]
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
            return _done(norm, it, segs)
        stalled = it >= 1 and norm > 0.5 * last_norm
        if stalled and norm < accept_tol:
            # integration-noise floor reached (the floor grows with the orbit's instability)
            return _done(norm, it, segs)
        if (
            it == max_iter
            or not math.isfinite(norm)
            or (not damped and it >= 4 and norm > 0.9 * last_norm)
        ):
            raise CorrectionError(f"multiple shooting did not converge (residual {norm:.3e})")
        last_norm = norm
        step = np.linalg.lstsq(jm, -f, rcond=None)[0]
        lam = 1.0
        if damped and norm > 1e-6:
            # backtracking on the max-norm residual (state-only propagation)
            for _ in range(8):
                tn = nodes + lam * step[:-1].reshape(n, 4)
                tp = period + lam * float(step[-1])
                try:
                    if tp > 0.0 and _residual_norm(tn, tp) < norm:
                        break
                except _STEP_FAILURES:
                    pass
                lam *= 0.5
        nodes = nodes + lam * step[:-1].reshape(n, 4)
        period = period + lam * float(step[-1])
        if not period > 0.0:
            raise CorrectionError("multiple shooting: period became non-positive")
    raise CorrectionError("unreachable")


def family_tangent(
    orbit: MSOrbit, previous: FloatArray | None = None, *, gap: float = 1e-4
) -> FloatArray:
    """Unit tangent of the family at fixed mu (null vector of continuity + section).

    One continuity row is redundant (the Jacobi integral), so the null vector is the right
    singular vector of the smallest singular value. At a branch point (a pair at k = +2,
    Casoliva scale) the null space is two-dimensional; when the two smallest singular values
    are both below ``gap`` times the third, the ``previous`` tangent projected onto that plane
    is returned, which keeps the continuation on the branch it came along.
    """
    _res, jac, _ = _evaluate(orbit.mu, orbit.nodes, orbit.period, orbit.fracs)
    nu = jac.shape[1]
    if orbit.phase_mode == "flow":
        ph = np.zeros((1, nu))
        f0 = planar_eom(orbit.nodes[0], orbit.mu)
        ph[0, :4] = f0 / np.linalg.norm(f0)
    else:
        ph = np.eye(nu)[1:2]
    jm = np.vstack([jac, ph])
    _, sv, vt = np.linalg.svd(jm)
    if previous is not None and sv[-2] < gap * sv[-3]:
        basis = vt[-2:]
        proj = basis.T @ (basis @ previous)
        norm = float(np.linalg.norm(proj))
        if norm > 0.0:
            return np.asarray(proj / norm, dtype=np.float64)
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
    reason: str  # "target", "stop:<why>", "min_step", "max_steps", "earth_impact", "stalled"
    rejected_newton: int = 0  # trial steps rejected because the corrector failed
    rejected_jump: int = 0  # trial steps rejected by the node-jump guard


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
    n_newton = n_jump = 0
    for _ in range(max_steps):
        lc = math.log(cur.mu)
        if abs(lt - lc) < 1e-14:
            return LegResult(members, "target", n_newton, n_jump)
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
            trial = cur.like(mu_new, nodes, period)
            try:
                new = correct(
                    trial,
                    fix_period=value if fix == "period" else None,
                    fix_jacobi=value if fix == "jacobi" else None,
                    max_iter=8,
                )
                ok = _max_node_change(new, cur) <= max_jump
                n_jump += 0 if ok else 1
            except _STEP_FAILURES:
                ok = False
                n_newton += 1
            if ok:
                break
            step *= 0.5
            if abs(step) < dlog_min:
                return LegResult(members, "min_step", n_newton, n_jump)
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
                return LegResult(members, "stop:" + why, n_newton, n_jump)
    return LegResult(members, "max_steps", n_newton, n_jump)


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
    stall_window: int = 10,
    stall_dc: float = 1e-7,
    earth_stop: float = EARTH_RADIUS_ND,
    stop: StopFn | None = None,
    log: LogFn | None = None,
) -> LegResult:
    """Pseudo-arclength continuation of the family through ``orbit`` at fixed mu, starting in
    the direction of increasing (``direction > 0``) or decreasing Jacobi constant; the
    direction is then kept by tangent orientation (folds in C are passed). Stops at an Earth
    impact (perigee below ``earth_stop``, default the Earth's radius; the primary is not
    regularised) and when C moves by less than ``stall_dc`` over
    ``stall_window`` members (a walk along a degenerate direction)."""
    members: list[Member] = []
    cur = orbit
    tan = family_tangent(cur)
    if _jacobi_rate(cur, tan) * direction < 0.0:
        tan = -tan
    ds = ds0
    n_newton = n_jump = 0
    for _ in range(max_steps):
        while True:
            u0 = cur.unknowns()
            up = u0 + ds * tan
            n = cur.n
            trial = cur.like(cur.mu, up[:-1].reshape(n, 4), float(up[-1]))
            try:
                new = correct(trial, arclength=(u0, tan, ds), max_iter=8)
                ok = _max_node_change(new, cur) <= max_jump
                n_jump += 0 if ok else 1
            except _STEP_FAILURES:
                ok = False
                n_newton += 1
            if ok:
                break
            ds *= 0.5
            if ds < ds_min:
                return LegResult(members, "min_step", n_newton, n_jump)
        mem = Member(leg, new, diagnose(new), ds)
        members.append(mem)
        if log is not None:
            log(mem)
        if mem.diag.perigee < earth_stop:
            return LegResult(members, "earth_impact", n_newton, n_jump)
        if len(members) > stall_window and (
            abs(members[-1].diag.jacobi - members[-1 - stall_window].diag.jacobi) < stall_dc
        ):
            return LegResult(members, "stalled", n_newton, n_jump)
        new_tan = family_tangent(new, tan)
        if float(new_tan @ tan) < 0.0:
            new_tan = -new_tan
        cur, tan = new, new_tan
        if new.iterations <= 4:
            ds = min(ds * 1.5, ds_max)
        if stop is not None:
            why = stop(mem)
            if why:
                return LegResult(members, "stop:" + why, n_newton, n_jump)
    return LegResult(members, "max_steps", n_newton, n_jump)


# ---------------------------------------------------------------------------------------------
# Section crossings and comparison with printed rows


def y_crossings(orbit: MSOrbit, *, steps_per_segment: int = 60) -> list[tuple[float, FloatArray]]:
    """All crossings of y = 0 over one period, as (time from node 0, planar state)."""
    mu = orbit.mu
    out: list[tuple[float, FloatArray]] = []
    fr = orbit.fracs
    starts = np.concatenate([[0.0], np.cumsum(fr)[:-1]]) * orbit.period
    for i in range(orbit.n):
        k = max(4, round(steps_per_segment * orbit.n * float(fr[i])))
        h = orbit.period * float(fr[i]) / k
        s = orbit.nodes[i]
        t = float(starts[i])
        for _ in range(k):
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


def orbit_distance(orbit: MSOrbit, state4: FloatArray) -> tuple[float, bool, float]:
    """Smallest max-abs difference between ``state4`` (on y = 0) and a y = 0 crossing of
    ``orbit`` or of its mirror image: (distance, mirrored, crossing time)."""
    best = (math.inf, False, math.nan)
    for t, st in y_crossings(orbit):
        for mir, cand in ((False, st), (True, mirror(st))):
            dist = float(np.max(np.abs(cand - state4)))
            if dist < best[0]:
                best = (dist, mir, t)
    return best


def same_orbit(a: MSOrbit, b: MSOrbit, tol: float = 1e-7) -> bool:
    """True when a y = 0 crossing of ``b`` lies on ``a`` (or on its mirror image) within
    ``tol``."""
    crossings = y_crossings(b)
    if not crossings:
        return orbit_distance(a, b.nodes[0])[0] < tol

    def clearance(st: FloatArray) -> float:
        return min(math.hypot(st[0] + b.mu, st[1]), math.hypot(st[0] - 1.0 + b.mu, st[1]))

    # the crossing farthest from both primaries: a crossing at a fast perigee or lunar pass
    # is the least well located by the sampled search
    probe = max((st for _, st in crossings), key=clearance)
    return orbit_distance(a, probe)[0] < tol


# ---------------------------------------------------------------------------------------------
# Resonant members at fixed mu


@dataclass(frozen=True)
class ResonantHit:
    orbit: MSOrbit
    diag: Diagnostics
    direction: float  # walk direction in C that found it
    member_index: int  # index of the walk member before the crossing


def resonant_crossings(
    orbit: MSOrbit,
    q: int,
    *,
    directions: tuple[float, ...] = (1.0, -1.0),
    max_steps: int = 300,
    ds_max: float = 0.05,
    period_window: float = 0.25,
    log: LogFn | None = None,
) -> tuple[list[ResonantHit], dict[float, LegResult]]:
    """Walk the family through ``orbit`` at fixed mu in C (both directions) and correct every
    crossing of T = 2 pi q at that period. The walk stops at an Earth impact, a stall, when
    |T - 2 pi q| exceeds ``period_window`` times 2 pi q, or after ``max_steps``."""
    tq = 2.0 * math.pi * q

    def stop(m: Member) -> str | None:
        if abs(m.diag.period - tq) > period_window * tq:
            return "period_window"
        return None

    hits: list[ResonantHit] = []
    legs: dict[float, LegResult] = {}
    for direction in directions:
        leg = continue_jacobi(
            orbit, direction, leg="walk@mu", ds_max=ds_max, max_steps=max_steps, stop=stop, log=log
        )
        legs[direction] = leg
        chain = [orbit, *[m.orbit for m in leg.members]]
        for i in range(len(chain) - 1):
            a, b = chain[i], chain[i + 1]
            fa, fb = a.period - tq, b.period - tq
            if fa == 0.0 or fa * fb < 0.0:
                w = fa / (fa - fb) if fa != fb else 0.0
                guess = a.like(a.mu, a.nodes + w * (b.nodes - a.nodes), tq)
                try:
                    hit = correct(guess, fix_period=tq, tol=1e-11)
                except _STEP_FAILURES:
                    continue
                if any(same_orbit(h.orbit, hit) for h in hits):
                    continue
                hits.append(ResonantHit(hit, diagnose(hit), direction, i))
    return hits, legs


# ---------------------------------------------------------------------------------------------
# The three-step strategy (Casoliva et al. 2010 section IV.C)


def moon_radius_scaled(mu: float) -> float:
    """The Moon's radius scaled at constant density to a secondary of mass ratio ``mu``:
    R_M (mu / mu_M)^(1/3) (the physical radius at the paper's mu). DECISION (#899 step 2):
    the surface exclusion used as the impact trigger during a continuation, since the
    physical radius has no meaning at mu = 1e-6, where the seeds pass at 1e-4 lunar
    distances."""
    return float(MOON_RADIUS_ND * (mu / CASOLIVA_MU_2010) ** (1.0 / 3.0))


def impact_trigger(window: int = 3) -> StopFn:
    """Stop a mass leg when the periselene is below :func:`moon_radius_scaled` and has fallen
    over the last ``window`` members (the family heading for the Moon)."""
    history: list[float] = []

    def stop(m: Member) -> str | None:
        history.append(m.diag.periselene)
        if m.diag.periselene < moon_radius_scaled(m.diag.mu) and len(history) > window:
            recent = history[-window - 1 :]
            if all(b < a for a, b in pairwise(recent)):
                return "lunar_impact_trend"
        return None

    return stop


@dataclass
class PathResult:
    legs: list[tuple[str, LegResult]]
    reached: bool
    final: MSOrbit | None
    reason: str = ""  # "target", "no_c_leg", "creep", "max_switches"


def three_step(
    orbit: MSOrbit,
    mu_target: float = CASOLIVA_MU_2010,
    *,
    first_fix: str = "period",
    max_switches: int = 4,
    c_walk_dc: float = 0.05,
    rp_gain: float = 1.5,
    c_walk_steps: int = 60,
    dlog_max: float = 0.2,
    min_mu_gain: float = 1.1,
    log: LogFn | None = None,
) -> PathResult:
    """Casoliva's strategy: continue in mu (first at ``first_fix``); when the leg ends short of
    ``mu_target`` (impact trend, fold, corrector failure), continue in C at fixed mu in the
    direction that raises the periselene, until it has risen by ``rp_gain`` or C has moved by
    ``c_walk_dc``; then resume in mu at fixed C. At most ``max_switches`` C legs; the path is
    abandoned ("creep") when a resumed mu leg gains less than a factor ``min_mu_gain`` in mu."""
    legs: list[tuple[str, LegResult]] = []
    cur = orbit
    fix = first_fix
    for switch in range(max_switches + 1):
        mu_start = cur.mu
        leg = continue_mu(
            cur, mu_target, fix=fix, dlog_max=dlog_max, stop=impact_trigger(), log=log
        )
        legs.append(("mu@T" if fix == "period" else "mu@C", leg))
        if leg.reason == "target":
            return PathResult(legs, True, leg.members[-1].orbit, "target")
        base = leg.members[-1].orbit if leg.members else cur
        if switch > 0 and base.mu < min_mu_gain * mu_start:
            return PathResult(legs, False, None, "creep")
        rp0 = diagnose(base).periselene
        c0 = base.jacobi
        best: tuple[float, LegResult | None] = (-math.inf, None)
        for direction in (1.0, -1.0):

            def stop(m: Member, c0: float = c0, rp0: float = rp0) -> str | None:
                if m.diag.periselene >= rp_gain * rp0:
                    return "periselene_raised"
                if abs(m.diag.jacobi - c0) >= c_walk_dc:
                    return "dc_limit"
                return None

            walk = continue_jacobi(
                base, direction, max_steps=c_walk_steps, ds_max=0.02, stop=stop, log=log
            )
            if walk.members:
                gain = walk.members[-1].diag.periselene / rp0
                if gain > best[0]:
                    best = (gain, walk)
        if best[1] is None or not best[1].members:
            return PathResult(legs, False, None, "no_c_leg")
        legs.append(("C@mu", best[1]))
        cur = best[1].members[-1].orbit
        fix = "jacobi"
    return PathResult(legs, False, None, "max_switches")


# ---------------------------------------------------------------------------------------------
# Second-species seeds from the matched in/out maps (Casoliva 2008 Eqs. 14-18, Barrabes &
# Gomez 2003 Eqs. 45-46, planar case)


def jacobi_interval(p: int, q: int) -> tuple[float, float]:
    """Eq. 15: the C_J interval of the p-q family (p spacecraft, q lunar revolutions)."""
    inv_a = (p / q) ** (2.0 / 3.0)
    half = 2.0 * math.sqrt(2.0 - inv_a)
    return inv_a - half, inv_a + half


def seed_directions(p: int, q: int, c_j: float) -> tuple[float, float]:
    """Eq. 16: the two polar angles psi of the velocity, (asin s, pi - asin s)."""
    inv_a = (p / q) ** (2.0 / 3.0)
    s = (2.0 - c_j + inv_a) / (2.0 * math.sqrt(3.0 - c_j))
    s = min(1.0, max(-1.0, s))
    psi = math.asin(s)
    return psi, math.pi - psi


def _eq17(theta: float, psi: float, c_j: float) -> float:
    d = theta - psi
    return (2.0 / math.sqrt(3.0 - c_j)) * (
        math.cos(theta) * math.sin(d) ** 2 - math.cos(psi) * math.cos(d)
    ) + math.sin(d) ** 3


def seed_angle(psi: float, c_j: float, samples: int = 3600) -> float:
    """Eq. 17: the angle theta on the circle about the Moon, with cos(theta - psi) > 0."""
    from scipy.optimize import brentq

    grid = np.linspace(0.0, 2.0 * math.pi, samples + 1)
    vals = [_eq17(t, psi, c_j) for t in grid]
    roots = []
    for i in range(samples):
        if vals[i] == 0.0 or vals[i] * vals[i + 1] < 0.0:
            r = float(brentq(_eq17, grid[i], grid[i + 1], args=(psi, c_j), xtol=1e-15))
            if math.cos(r - psi) > 0.0:
                roots.append(r)
    if not roots:
        raise ValueError("seed_angle: no root of Eq. 17 with cos(theta - psi) > 0")
    return roots[0]


def seed_state(
    p: int, q: int, c_j: float, branch: int, mu: float, alpha: float = 0.4
) -> FloatArray:
    """Eq. 14 initial condition on the circle of radius mu^alpha about the Moon, project frame.

    The speed comes from the exact Jacobi relation at ``mu`` (2008 digest recipe, step 4);
    ``branch`` 0 or 1 picks psi = asin(s) or pi - asin(s). The only printed alpha is 0.4
    (Barrabes & Gomez); the period guess is 2 pi q (Eq. 18)."""
    psi = seed_directions(p, q, c_j)[branch]
    theta = seed_angle(psi, c_j)
    rho = mu**alpha
    x = mu - 1.0 + rho * math.cos(theta)  # Casoliva frame: Moon at (mu - 1, 0)
    y = rho * math.sin(theta)
    r1 = math.hypot(x - mu, y)
    r2 = math.hypot(x - mu + 1.0, y)
    two_omega = x * x + y * y + 2.0 * (1.0 - mu) / r1 + 2.0 * mu / r2
    v = math.sqrt(two_omega - c_j)
    u, w = v * math.cos(psi), v * math.sin(psi)
    return np.array([-x, -y, -u, -w])  # rotation by pi to the project frame


def first_far_crossing(mu: float, state4: FloatArray, period: float, n: int) -> MSOrbit:
    """Nodes for a seed that is not on y = 0: the seed trajectory is followed for ``period``
    and node 0 is put at its y = 0 crossing farthest from the Moon."""
    tmp = MSOrbit(mu, nodes_from_state(mu, state4, period, n), period)
    crossings = y_crossings(tmp)
    if not crossings:
        raise ValueError("first_far_crossing: the seed trajectory does not cross y = 0")
    _, st = max(crossings, key=lambda c: math.hypot(c[1][0] - 1.0 + mu, c[1][1]))
    return MSOrbit(mu, nodes_from_state(mu, st, period, n), period)


def corrected_seed(
    p: int,
    q: int,
    c_j: float,
    branch: int,
    *,
    mu: float = SEED_MU,
    alpha: float = 0.4,
    n: int = 12,
    max_iter: int = 30,
) -> MSOrbit:
    """A seed of :func:`seed_state` corrected at fixed C_J and ``mu`` (2008 p.9: "the initial
    conditions had to be differentially corrected using a grid of C_J values"). Node 0 stays
    at the seed point on the circle about the Moon (phase condition orthogonal to the flow
    there), where the matching of the in and out maps is made; an undamped Gauss-Newton
    attempt is followed by a backtracking one if it fails; the corrected orbit is then
    re-noded from its y = 0 crossing farthest from the Moon for continuation."""
    st = seed_state(p, q, c_j, branch, mu, alpha)
    period = 2.0 * math.pi * q
    guess = MSOrbit(mu, nodes_from_state(mu, st, period, n), period)
    flow = planar_eom(st, mu)
    phase = (st, flow / np.linalg.norm(flow))
    try:
        orbit = correct(guess, fix_jacobi=c_j, phase=phase, max_iter=max_iter)
    except _STEP_FAILURES:
        orbit = correct(guess, fix_jacobi=c_j, phase=phase, damped=True, max_iter=max_iter)
    far = first_far_crossing(mu, orbit.nodes[0], orbit.period, n)
    return correct(far, fix_jacobi=c_j)


# ---------------------------------------------------------------------------------------------
# Chains of returning collision arcs (two or more lunar encounters per period)


@dataclass(frozen=True)
class ReturningArc:
    """A Kepler arc that leaves the Moon's position and returns to it after ``i`` lunar
    revolutions and ``j`` particle revolutions (Henon's same-point arcs, A = i/j; inverse
    semi-major axis (j/i)^(2/3)); ``sign`` is the sign of the radial relative velocity at the
    Moon (+1 outward). At mu = 0 the arc returns with the velocity it left with."""

    i: int
    j: int
    sign: int


def arc_relative_velocity(arc: ReturningArc, c_j: float) -> FloatArray:
    """Rotating-frame relative velocity at the Moon (project frame: radial +x, tangential +y)
    of a returning arc at Jacobi constant ``c_j`` (mu = 0): V^2 = 3 - C and the tangential
    component (C - 2 - 1/a)/2 from the vis-viva speed sqrt(2 - 1/a)."""
    inv_a = (arc.j / arc.i) ** (2.0 / 3.0)
    v2 = 3.0 - c_j
    vt = (c_j - 2.0 - inv_a) / 2.0
    vr2 = v2 - vt * vt
    if vr2 < 0.0:
        raise ValueError(f"arc {arc} does not reach the Moon at C = {c_j}")
    return np.array([arc.sign * math.sqrt(vr2), vt])


def flyby_states(
    v_in: FloatArray, v_out: FloatArray, mu: float, rho: float
) -> tuple[FloatArray, FloatArray, float]:
    """Moon-centred two-body hyperbola turning ``v_in`` into ``v_out`` (equal magnitudes): the
    states (relative position, velocity) where it crosses the circle of radius ``rho`` inbound
    and outbound in the rotating frame, and the time from that circle to periapsis. First-order
    matching (Gomez & Olle; Guillaume 1975): sin(delta/2) = 1/e, r_p = mu (e - 1)/V^2. The
    Earth's tide over the passage is neglected."""
    v = float(np.linalg.norm(v_in))
    u1 = v_in / v
    u2 = v_out / float(np.linalg.norm(v_out))
    cosd = float(np.clip(u1 @ u2, -1.0, 1.0))
    delta = math.acos(cosd)
    if delta < 1e-12:
        raise ValueError("flyby_states: zero turn (a single-arc orbit, not a chain)")
    e = 1.0 / math.sin(0.5 * delta)
    rp = mu * (e - 1.0) / (v * v)
    if rho <= rp:
        raise ValueError("flyby_states: circle inside the periapsis")
    p = rp * (1.0 + e)
    phat = (u1 - u2) / np.linalg.norm(u1 - u2)
    qhat = (u1 + u2) / np.linalg.norm(u1 + u2)
    nu0 = math.acos((p / rho - 1.0) / e)
    k = math.sqrt(mu / p)

    def state(nu: float) -> tuple[FloatArray, FloatArray]:
        r = p / (1.0 + e * math.cos(nu))
        pos = r * (math.cos(nu) * phat + math.sin(nu) * qhat)
        vel = k * (-math.sin(nu) * phat + (e + math.cos(nu)) * qhat)
        return pos, vel

    cosh_f = (e + math.cos(nu0)) / (1.0 + e * math.cos(nu0))
    big_f = math.acosh(cosh_f)
    a_abs = mu / (v * v)
    t_half = (e * math.sinh(big_f) - big_f) / math.sqrt(mu / a_abs**3)

    def rotating(nu: float, tau: float) -> FloatArray:
        # the hyperbola lives in the moon-centred non-rotating frame aligned with the rotating
        # frame at periapsis; at time tau from periapsis, r_rot = R(-tau) r and
        # v_rot = R(-tau) v - omega x r_rot (the frame rotation moves the impact parameter by
        # about V tau^2, larger than the impact parameter itself at mu = 1e-6)
        pos, vel = state(nu)
        c, s_ = math.cos(-tau), math.sin(-tau)
        rot = np.array([[c, -s_], [s_, c]])
        pr = rot @ pos
        vr = rot @ vel - np.array([-pr[1], pr[0]])
        return np.concatenate([pr, vr])

    return rotating(-nu0, -t_half), rotating(nu0, t_half), t_half


def chain_seed(
    arcs: Sequence[ReturningArc],
    c_j: float,
    mu: float,
    *,
    alpha: float = 0.4,
    nodes_per_arc: int = 6,
) -> MSOrbit:
    """Multiple-shooting guess for the periodic orbit generated by a chain of returning arcs
    joined by lunar flybys at Jacobi constant ``c_j``: for each junction, nodes on the circle
    of radius mu^alpha inbound and outbound on the matched hyperbola, the passage between them
    as one segment, and ``nodes_per_arc`` nodes along each arc (the first half propagated
    forward from the outbound node, the second half backward from the inbound one). Period
    2 pi sum(i). Node 0 is the first outbound node, with a flow phase condition.
    """
    m = len(arcs)
    if m < 2:
        raise ValueError("chain_seed: a chain needs at least two arcs")
    rho = mu**alpha
    moon = np.array([1.0 - mu, 0.0])
    vel = [arc_relative_velocity(a, c_j) for a in arcs]
    # junction k joins arc k-1 (inbound) to arc k (outbound)
    junctions = [flyby_states(vel[k - 1], vel[k], mu, rho) for k in range(m)]
    period = 2.0 * math.pi * sum(a.i for a in arcs)
    nodes: list[FloatArray] = []
    durations: list[float] = []
    for k, arc in enumerate(arcs):
        _, out_k, t_k = junctions[k]
        in_next, _, t_next = junctions[(k + 1) % m]
        start = np.array([moon[0] + out_k[0], moon[1] + out_k[1], out_k[2], out_k[3]])
        end = np.array([moon[0] + in_next[0], moon[1] + in_next[1], in_next[2], in_next[3]])
        arc_time = 2.0 * math.pi * arc.i - t_k - t_next
        h = arc_time / nodes_per_arc
        half = nodes_per_arc // 2
        # first half forward from the outbound node, second half backward from the inbound
        # one (time reversal: mirror, propagate, mirror), so the seed's mismatch sits mid-arc
        for jn in range(nodes_per_arc):
            if jn == 0:
                node = start
            elif jn < half:
                node = propagate_segment(mu, start, jn * h, with_stm=False).end
            else:
                back = propagate_segment(mu, mirror(end), (nodes_per_arc - jn) * h, with_stm=False)
                node = mirror(back.end)
            nodes.append(node)
            durations.append(h)
        nodes.append(end)
        durations.append(2.0 * t_next)
    fr = np.array(durations) / period
    return MSOrbit(mu, np.array(nodes), period, fractions=fr / fr.sum(), phase_mode="flow")


def corrected_chain(
    arcs: Sequence[ReturningArc],
    c_j: float,
    *,
    mu: float = SEED_MU,
    alpha: float = 0.4,
    nodes_per_arc: int = 6,
    max_iter: int = 30,
) -> MSOrbit:
    """:func:`chain_seed` corrected at fixed C_J (undamped, then with backtracking)."""
    guess = chain_seed(arcs, c_j, mu, alpha=alpha, nodes_per_arc=nodes_per_arc)
    try:
        return correct(guess, fix_jacobi=c_j, max_iter=max_iter)
    except _STEP_FAILURES:
        return correct(guess, fix_jacobi=c_j, damped=True, max_iter=max_iter)


def two_arc_chains(p: int, q: int) -> list[tuple[ReturningArc, ReturningArc]]:
    """Two-arc chains of returning arcs for the p-q resonance: (i1, j1) + (i2, j2) with
    i1 + i2 = q, j1 + j2 = p, each pair coprime (a non-primitive returning arc meets the Moon
    half-way) and each arc able to reach the Moon's orbit (inverse semi-major axis below 2);
    radial-sign pairs (+,+) and (+,-) only, since (-,-) and (-,+) are their mirror images."""
    out = []
    for i1 in range(1, q):
        i2 = q - i1
        for j1 in range(1, p):
            j2 = p - j1
            if (i1, j1) > (i2, j2):
                continue  # the cyclic order is irrelevant
            if math.gcd(i1, j1) != 1 or math.gcd(i2, j2) != 1:
                continue
            if (j1 / i1) ** (2.0 / 3.0) >= 2.0 or (j2 / i2) ** (2.0 / 3.0) >= 2.0:
                continue
            for s2 in (1, -1):
                out.append((ReturningArc(i1, j1, 1), ReturningArc(i2, j2, s2)))
    return out
