"""#899 step 2: second-species seeds continued in the mass ratio (Casoliva et al. 2008, 2010).

Expected values are the printed tables only: 2008 Table 2 (AIAA 2008-6434, p.10; seeds at
mu = 1e-6, transcribed in ``docs/notes/2026-10-04-digest-casoliva-2008-aiaa-families-cycler-
trajectories-seeds.md``) and 2010 Table 3 (JGCD 33(5) p.1630, vendored in
:mod:`cyclerfinder.search.earth_moon_resonant_families`). The planar Levi-Civita propagator is
checked against the independent #928 KS propagator of :mod:`cyclerfinder.core.cr3bp_ks`.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from cyclerfinder.core.cr3bp_ks import MoonCentredCR3BP, propagate_ks
from cyclerfinder.search import earth_moon_resonant_families as emrf
from cyclerfinder.search import second_species_continuation as ssc
from cyclerfinder.search.second_species_lc import propagate_lc

_PLANAR = [0, 1, 3, 4]


def _ks(mu: float, st4: np.ndarray, t: float) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    s6 = np.array([st4[0], st4[1], 0.0, st4[2], st4[3], 0.0])
    arc = propagate_ks(MoonCentredCR3BP(mu), s6, t, with_stm=True, rtol=1e-13, atol=1e-15)
    assert arc.stm is not None
    end = np.array([arc.state[0], arc.state[1], arc.state[3], arc.state[4]])
    return end, arc.stm[np.ix_(_PLANAR, _PLANAR)], arc.stm[np.ix_([2, 5], [2, 5])], arc.r_min


@pytest.mark.parametrize(
    ("mu", "state", "t"),
    [
        # 2008 seed 32a at mu = 1e-6: a lunar pass at 9.5e-5 lunar distances
        (ssc.SEED_MU, ssc.table2_seed("32a").state_project(), ssc.table2_seed("32a").period / 3),
        # 2008 seed 73a, a third of the period
        (ssc.SEED_MU, ssc.table2_seed("73a").state_project(), ssc.table2_seed("73a").period / 3),
        # 2010 row 2-1b at the paper's mu, full period (perigee 0.0177)
        (
            emrf.CASOLIVA_MU_2010,
            emrf.table3_seed_state(emrf.table3_row("2-1b"))[[0, 1, 3, 4]],
            2.0 * math.pi,
        ),
    ],
)
def test_lc_propagator_matches_ks(mu: float, state: np.ndarray, t: float) -> None:
    # both at rtol 1e-13: at 1e-12 each integrator alone is off by up to 1.5e-9 over the 73a
    # arc (measured against its own 1e-14 run), and the two then differ at that level
    arc = propagate_lc(mu, state, t, rtol=1e-13, atol=1e-15)
    end, m4, mz, r_min = _ks(mu, state, t)
    assert arc.stm4 is not None
    assert arc.stm_z is not None
    assert np.max(np.abs(arc.end - end)) < 5e-11
    assert np.max(np.abs(arc.stm4 - m4)) / np.max(np.abs(m4)) < 1e-7
    assert np.max(np.abs(arc.stm_z - mz)) / np.max(np.abs(mz)) < 1e-9
    assert abs(arc.r_min - r_min) / r_min < 1e-9
    assert abs(arc.jacobi_drift) < 5e-11


def test_lc_propagator_mu_zero_limit_kepler_period() -> None:
    """mu -> 0 (1e-12) with the particle far from the secondary: a rotating-frame Kepler
    circle of radius 0.5 about the primary returns after its synodic period
    2 pi / (n - 1), n = 0.5^-1.5 (closed form)."""
    mu = 1e-12
    r = 0.5
    n = r**-1.5
    st = np.array([r - mu, 0.0, 0.0, r * (n - 1.0)])
    t_syn = 2.0 * math.pi / (n - 1.0)
    arc = propagate_lc(mu, st, t_syn, with_stm=False)
    assert np.max(np.abs(arc.end - st)) < 1e-9


@pytest.mark.parametrize("seed", ssc.TABLE2_SEEDS, ids=lambda s: s.designation)
def test_table2_seed_closes_at_mu_1e6(seed: ssc.Table2Seed) -> None:
    """The ten printed 2008 seeds close at mu = 1e-6 (fixed printed C_J) with corrections
    below 1e-7 in state and 1e-9 in period, and the planar index lambda + 1/lambda matches the
    printed k to its four printed decimals (relative 3e-5 for 12a, k = 994.8)."""
    orbit = ssc.correct(
        ssc.renode(ssc.MSOrbit(ssc.SEED_MU, np.array([seed.state_project()]), seed.period), 4),
        fix_jacobi=seed.c_j,
    )
    assert orbit.converged
    assert np.max(np.abs(orbit.nodes[0] - seed.state_project())) < 1e-7
    assert abs(orbit.period - seed.period) < 1e-9
    k_par = ssc.diagnose(orbit).k_par
    assert abs(k_par - seed.k) <= max(1e-4, 3e-5 * abs(seed.k))


@pytest.mark.parametrize("designation", emrf.TABLE3_VALID_DESIGNATIONS)
def test_table3_row_closes_at_paper_mu(designation: str) -> None:
    """The nine catalogued 2010 rows, corrected at the paper's mu with T = 2 pi q fixed, move
    by at most 3e-10 from the printed crossing (the printed digits), keep the printed C_J to
    2e-10, and reproduce the printed k (the block chosen by Casoliva, per #801) to 1e-7."""
    row = emrf.table3_row(designation)
    ref = ssc.table3_reference(designation)
    printed = np.array([-row.x_i, 0.0, -row.u_i, -row.v_i])
    assert np.max(np.abs(ref.nodes[0] - printed)) < 3e-10
    assert abs(ref.jacobi - row.c_j) < 2e-10
    g = ssc.diagnose(ref)
    k = ssc.casoliva_k(designation, g)
    assert abs(k - row.k) / abs(row.k) < 1e-7


def test_registry_mu_misses_the_printed_digits() -> None:
    """At the registry mu the same correction moves the printed 7-3a crossing by 1e-6 or more:
    the printed digits belong to the paper's mu, which is the continuation target."""
    row = emrf.table3_row("7-3a")
    ref = ssc.table3_reference("7-3a", mu=emrf.earth_moon_system().mu)
    printed = np.array([-row.x_i, 0.0, -row.u_i, -row.v_i])
    assert np.max(np.abs(ref.nodes[0] - printed)) > 1e-6


def test_walk_at_mu_finds_the_7_3a_crossing() -> None:
    """Crossing finder control with a known answer: from a member of the 7-3a family a few
    steps away in C, the walk at the paper's mu finds T = 6 pi and the correction lands on
    the printed 7-3a crossing."""
    row = emrf.table3_row("7-3a")
    ref = ssc.table3_reference("7-3a")
    away = ssc.continue_jacobi(ref, 1.0, max_steps=5, ds_max=0.02).members[-1].orbit
    assert abs(away.period - 6.0 * math.pi) > 1e-3
    hits, _ = ssc.resonant_crossings(away, 3, directions=(-1.0,), max_steps=20)
    printed = np.array([-row.x_i, 0.0, -row.u_i, -row.v_i])
    dists = [ssc.orbit_distance(h.orbit, printed)[0] for h in hits]
    assert min(dists) < 1e-9


def test_7_3b_and_7_3c_are_one_orbit_and_its_mirror() -> None:
    """The printed 7-3c crossing is the mirror image of a crossing of the 7-3b orbit (same
    printed C_J and k): one asymmetric orbit, not two."""
    b = ssc.table3_reference("7-3b")
    c = ssc.table3_reference("7-3c")
    dist, mirrored, _ = ssc.orbit_distance(b, c.nodes[0])
    assert dist < 1e-9
    assert mirrored


def test_eq14_seed_generator_reproduces_table2_seeds() -> None:
    """Positive control of the Eq. 14-17 seed generator: at the printed C_J of the 2008 seeds
    54a (psi branch 0) and 73a (branch 1), the corrected seed is the printed Table 2 orbit."""
    for name, branch in (("54a", 0), ("73a", 1)):
        s = ssc.table2_seed(name)
        orbit = ssc.corrected_seed(s.p, s.q, s.c_j, branch)
        assert ssc.orbit_distance(orbit, s.state_project())[0] < 1e-8
        assert abs(orbit.period - s.period) < 1e-8


def test_jacobi_interval_matches_eq15() -> None:
    lo, hi = ssc.jacobi_interval(1, 2)
    assert abs(lo - (-1.711013183)) < 1e-9
    assert abs(hi - 2.970934233) < 1e-9


def test_two_arc_chain_converges_at_mu_1e6() -> None:
    """A (1,2) + (2,5) returning-arc chain at C = 1.1 closes at mu = 1e-6 with two lunar
    passes per period at a periselene of order mu (second species) and T near 6 pi; the sign
    pairs (+,+) and (-,-) give one orbit and its mirror."""
    arc = ssc.ReturningArc
    a = ssc.corrected_chain([arc(1, 2, 1), arc(2, 5, 1)], 1.1)
    b = ssc.corrected_chain([arc(1, 2, -1), arc(2, 5, -1)], 1.1)
    for o in (a, b):
        g = ssc.diagnose(o)
        assert 0.5 * ssc.SEED_MU < g.periselene < 10.0 * ssc.SEED_MU
        assert abs(o.period - 6.0 * math.pi) < 1e-4
    assert abs(ssc.diagnose(a).k_par - ssc.diagnose(b).k_par) / ssc.diagnose(a).k_par < 1e-4
    assert ssc.same_orbit(a, b, tol=1e-6)


def test_forward_reproduction_7_3b_from_a_second_species_chain() -> None:
    """#899 step 2 forward reproduction: a (1,2) + (2,5) chain seed at mu = 1e-6 and a grid
    value C = 1.05 (not the printed one), continued in mu at fixed C to the paper's mu, then
    walked in C to T = 6 pi, lands on the printed 7-3b crossing (and 7-3c is its mirror)."""
    arc = ssc.ReturningArc
    seed = ssc.corrected_chain([arc(1, 2, -1), arc(2, 5, -1)], 1.05)
    leg = ssc.continue_mu(seed, ssc.CASOLIVA_MU_2010, fix="jacobi", dlog_max=0.2)
    assert leg.reason == "target"
    hits, _ = ssc.resonant_crossings(leg.members[-1].orbit, 3, directions=(1.0,), max_steps=15)
    row_b = emrf.table3_row("7-3b")
    row_c = emrf.table3_row("7-3c")
    best = min(
        hits,
        key=lambda h: ssc.orbit_distance(
            h.orbit, np.array([-row_b.x_i, 0.0, -row_b.u_i, -row_b.v_i])
        )[0],
    )
    db = ssc.orbit_distance(best.orbit, np.array([-row_b.x_i, 0.0, -row_b.u_i, -row_b.v_i]))
    dc = ssc.orbit_distance(best.orbit, np.array([-row_c.x_i, 0.0, -row_c.u_i, -row_c.v_i]))
    assert db[0] < 1e-9
    assert dc[0] < 1e-9
    assert abs(best.diag.jacobi - row_b.c_j) < 1e-9
    assert abs(best.diag.k_par - row_b.k) / row_b.k < 1e-7


def test_forward_reproduction_2_1a_from_table2_seed_21a() -> None:
    """#899 step 2 forward reproduction: the 2008 seed 21a (mu = 1e-6) continued in mu at its
    fixed C_J to the paper's mu, then walked in C to T = 2 pi, lands on the printed 2-1a
    crossing (its mirror image)."""
    s = ssc.table2_seed("21a")
    start = ssc.MSOrbit(ssc.SEED_MU, np.array([s.state_project()]), s.period)
    seed = ssc.correct(ssc.renode(start, 12), fix_jacobi=s.c_j)
    leg = ssc.continue_mu(seed, ssc.CASOLIVA_MU_2010, fix="jacobi", dlog_max=0.4)
    assert leg.reason == "target"
    hits, _ = ssc.resonant_crossings(
        leg.members[-1].orbit, 1, directions=(-1.0,), max_steps=200, ds_max=0.2
    )
    row = emrf.table3_row("2-1a")
    printed = np.array([-row.x_i, 0.0, -row.u_i, -row.v_i])
    dists = [ssc.orbit_distance(h.orbit, printed)[0] for h in hits]
    assert min(dists) < 1e-9


@pytest.mark.slow
def test_forward_reproduction_3_2c_from_the_family_of_seed_32a() -> None:
    """#899 step 2 forward reproduction: the family of the 2008 seed 32a at mu = 1e-6, walked
    in C (decreasing first; it folds and turns up through a near-collision into resonant
    orbits with T near 4 pi), sampled at the first member past C = 0.27, continued in mu at
    that fixed C to the paper's mu, then walked in C to T = 4 pi, lands on the printed 3-2c
    crossing."""
    s = ssc.table2_seed("32a")
    start = ssc.MSOrbit(ssc.SEED_MU, np.array([s.state_project()]), s.period)
    seed = ssc.correct(ssc.renode(start, 12), fix_jacobi=s.c_j)
    walk = ssc.continue_jacobi(seed, -1.0, max_steps=300, ds_max=0.05)
    sample = next(
        m.orbit
        for k, m in enumerate(walk.members)
        if k > 0 and walk.members[k - 1].diag.jacobi < 0.27 <= m.diag.jacobi
    )
    leg = ssc.continue_mu(sample, ssc.CASOLIVA_MU_2010, fix="jacobi", dlog_max=0.2)
    assert leg.reason == "target"
    hits, _ = ssc.resonant_crossings(leg.members[-1].orbit, 2, max_steps=150)
    row = emrf.table3_row("3-2c")
    printed = np.array([-row.x_i, 0.0, -row.u_i, -row.v_i])
    best = min(hits, key=lambda h: ssc.orbit_distance(h.orbit, printed)[0])
    assert ssc.orbit_distance(best.orbit, printed)[0] < 1e-9
    assert abs(best.diag.jacobi - row.c_j) < 2e-10
    assert abs(ssc.casoliva_k("3-2c", best.diag) - row.k) / abs(row.k) < 1e-7


def test_matching_validity_speed_is_infinite_at_and_above_c_3() -> None:
    """C >= 3 means zero relative speed at collision; the 899 scan crashed dividing by V^3 = 0."""
    from cyclerfinder.search.second_species_continuation import matching_validity_speed

    mu = 0.0121529529
    assert matching_validity_speed(mu, 3.0) == math.inf
    assert matching_validity_speed(mu, 3.2) == math.inf
    assert matching_validity_speed(mu, 3.0 - 1e-300) == math.inf
    v = math.sqrt(3.0 - 1.0)
    assert matching_validity_speed(mu, 1.0) == pytest.approx(mu * abs(math.log(mu)) / v**3)
