"""Tests for the Jupiter-Europa resonant-orbit family module (#753).

Reproduce-before-trust: sourced constants (mu, C_flyby, Table 1 targets) are
checked verbatim against Anderson & Lo 2011's own stated values (verified
directly against the PDF text layer this task -- see the module docstring).
The Table 1 gate is HONEST: only ``5:6-LI`` is asserted to pass; the other
three rows are asserted to correctly FAIL (a documented, real negative
result -- not silently dropped). See
``docs/notes/2026-07-28-753-jupiter-europa-resonant-families-table1-gate.md``
for the full search history.
"""

from __future__ import annotations

import itertools
import math

import pytest

import cyclerfinder.search.cr3bp_periodic as cp
import cyclerfinder.search.jovian_resonant_families as jrf


@pytest.fixture(scope="module")
def system() -> object:
    return jrf.jupiter_europa_system()


# ---------------------------------------------------------------------------
# (1) Sourced-constant reproduce-before-trust gate.
# ---------------------------------------------------------------------------


def test_mu_matches_paper_p169() -> None:
    """Anderson & Lo 2011 p.169: mu = 2.5266448850435e-5, verbatim."""
    assert jrf.ANDERSON_LO_MU == 2.5266448850435e-5


def test_c_flyby_matches_paper_p179() -> None:
    """Anderson & Lo 2011 p.179: Cflyby = 2.99163956830415, verbatim."""
    assert jrf.ANDERSON_LO_C_FLYBY == 2.99163956830415


def test_table1_targets_match_paper_p184() -> None:
    """Anderson & Lo 2011 Table 1, p.184, verbatim (all four rows)."""
    assert jrf.TABLE1_TARGETS == {
        "3:4-LO": 1036.116088,
        "5:6-LI": 1.000008,
        "5:6-LO": 4445.387515,
        "5:6-NO": 28178.258323,
    }


def test_mu_differs_from_registry_by_small_known_amount(system: object) -> None:
    """The paper's mu is NOT this project's DE440-registry Jupiter-Europa mu
    (a documented ~0.034% GM-vintage difference, #745 digest) -- confirms
    :func:`jupiter_europa_system` really uses the paper's own value, not the
    registry's.
    """
    import cyclerfinder.core.cr3bp as cr3bp

    registry_sys = cr3bp.cr3bp_system("Jupiter", "Europa")
    rel = abs(registry_sys.mu - jrf.ANDERSON_LO_MU) / jrf.ANDERSON_LO_MU
    assert 1e-5 < rel < 1e-3, f"unexpected mu delta {rel:.2e} (expected ~3.4e-4)"
    assert system.mu == jrf.ANDERSON_LO_MU  # type: ignore[attr-defined]


# ---------------------------------------------------------------------------
# (2) Two-body p:q resonant-ellipse seed sanity (Anderson & Lo Fig. 1 / Eq. 5-6).
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(("p", "q"), [(3, 4), (5, 6)])
def test_two_body_seed_period_matches_pq_convention(p: int, q: int) -> None:
    """T_full = 2*pi*q (rotating-frame closure period) and the two-body
    spacecraft period ratio matches [Europa period]:[spacecraft period] = p:q
    (Anderson & Lo Eq. 5-6, verified directly against the paper's text)."""
    seed = jrf.two_body_resonant_seed(p, q)
    assert seed.period_full == pytest.approx(2.0 * math.pi * q)
    # Europa's own two-body period is 2*pi in these units; ratio must be p:q.
    europa_period = 2.0 * math.pi
    ratio = europa_period / seed.period_two_body
    assert ratio == pytest.approx(p / q, rel=1e-9)


def test_two_body_seed_periapse_at_secondary_radius() -> None:
    """Fig. 1 caption: periapse intersects Europa's orbit (r=1, GM=1 barycentric)."""
    seed = jrf.two_body_resonant_seed(3, 4, x0_sign=-1)
    assert abs(seed.x0) == pytest.approx(1.0)
    assert seed.eccentricity > 0.0  # a genuine ellipse, not a circle
    assert seed.semi_major_axis > 1.0  # exterior (spacecraft period > Europa's)


def test_two_body_seed_rejects_bad_pq() -> None:
    with pytest.raises(ValueError, match="positive integers"):
        jrf.two_body_resonant_seed(0, 4)
    with pytest.raises(ValueError, match="x0_sign"):
        jrf.two_body_resonant_seed(3, 4, x0_sign=0)


# ---------------------------------------------------------------------------
# (2b) Conjugate-apse ("encounter-phase") seed (#860 Sec. 4(c1), #861).
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(("p", "q"), [(3, 4), (4, 5), (5, 6), (4, 3), (5, 4), (6, 5)])
def test_conjugate_apse_seed_period_matches_pq_convention(p: int, q: int) -> None:
    """Same rotating-frame closure period convention as two_body_resonant_seed
    (T_full = 2*pi*q) -- the target resonance, not the seed point, fixes this."""
    seed = jrf.two_body_conjugate_apse_seed(p, q)
    assert seed.period_full == pytest.approx(2.0 * math.pi * q)
    assert seed.period_two_body == pytest.approx(2.0 * math.pi * q / p)


@pytest.mark.parametrize(("p", "q"), [(3, 4), (4, 5), (5, 6), (4, 3), (5, 4), (6, 5)])
def test_conjugate_apse_seed_at_far_apse_of_same_ellipse(p: int, q: int) -> None:
    """x0 = -(2a-1), the OTHER apse of the SAME physical ellipse
    two_body_resonant_seed(p, q, x0_sign=+1) uses (near apse at r=1)."""
    seed = jrf.two_body_conjugate_apse_seed(p, q)
    a = (q / p) ** (2.0 / 3.0)
    assert seed.semi_major_axis == pytest.approx(a)
    assert seed.r_near == pytest.approx(1.0)
    assert seed.r_far == pytest.approx(2.0 * a - 1.0)
    assert seed.x0 == pytest.approx(-seed.r_far)
    # Genuinely a different seed point than the opposition (x0_sign=-1) seed.
    opposition = jrf.two_body_resonant_seed(p, q, x0_sign=-1)
    assert seed.x0 != pytest.approx(opposition.x0, rel=1e-3)


@pytest.mark.parametrize(("p", "q"), [(3, 4), (4, 5), (5, 6), (4, 3), (5, 4), (6, 5)])
def test_conjugate_apse_seed_is_self_consistent_vis_viva_point(p: int, q: int) -> None:
    """The constructed IC's own two-body specific energy at r_far matches the
    vis-viva energy at the SAME semi-major-axis, GM=1 (self-consistency of the
    derivation, not a value asserted against itself): -1/(2a) = 0.5*v_far^2 -
    1/r_far, i.e. v_far^2 = 2/r_far - 1/a, with v_far the INERTIAL tangential
    speed |ydot0_rotating + omega x r| at the seed point (omega=1)."""
    seed = jrf.two_body_conjugate_apse_seed(p, q)
    v_far_inertial = seed.r_far - seed.ydot0  # invert ydot0 = r_far - v_far
    energy_vis_viva = -1.0 / (2.0 * seed.semi_major_axis)
    energy_from_state = 0.5 * v_far_inertial**2 - 1.0 / seed.r_far
    assert energy_from_state == pytest.approx(energy_vis_viva, rel=1e-9)


def test_conjugate_apse_seed_rejects_bad_pq() -> None:
    with pytest.raises(ValueError, match="positive integers"):
        jrf.two_body_conjugate_apse_seed(0, 4)
    # a too small -> r_far <= 0, same domain limit as two_body_resonant_seed
    with pytest.raises(ValueError):
        jrf.two_body_conjugate_apse_seed(5, 1)


def test_conjugate_apse_seed_converges_a_genuine_symmetric_orbit(system: object) -> None:
    """The seed is a valid perpendicular-crossing IC (xdot0=0 by construction)
    that the existing corrector converges from -- confirms the new seed slots
    into the SAME machinery two_body_resonant_seed feeds, unmodified."""
    import numpy as np

    import cyclerfinder.core.cr3bp as cr3bp

    seed = jrf.two_body_conjugate_apse_seed(4, 5)
    state0 = np.array([seed.x0, 0.0, 0.0, 0.0, seed.ydot0, 0.0])
    seed_jacobi = float(cr3bp.jacobi_constant(state0, system.mu))  # type: ignore[attr-defined]
    orbit = cp.correct_symmetric_fixed_jacobi(
        system,  # type: ignore[arg-type]
        seed.x0,
        seed_jacobi,
        seed.period_full,
        ydot0_sign=1.0 if seed.ydot0 >= 0.0 else -1.0,
        half_crossings=None,
    )
    assert orbit.converged
    assert orbit.crossing_residual < 1e-8  # default corrector tol


# ---------------------------------------------------------------------------
# (3) Hardcoded Table-1 candidate seeds: standing regression (still converge).
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("label", sorted(jrf.TABLE1_TARGETS))
def test_hardcoded_seed_converges(label: str, system: object) -> None:
    cand = jrf.recover_table1_candidate(label, system)  # type: ignore[arg-type]
    assert cand.crossing_residual < 1e-9
    assert cand.jacobi == pytest.approx(jrf.ANDERSON_LO_C_FLYBY, abs=1e-8)


def test_5_6_li_period_confirms_genuine_resonance_lineage(system: object) -> None:
    """5:6-LI's period/2pi must be ~6 (q=6) for the recovery to mean anything
    beyond a coincidental eigenvalue match -- this is the check that
    distinguishes the CONFIRMED candidate from the three NOT-confirmed ones
    (whose periods do not land on a clean 2*pi*q multiple, see the module
    docstring / results note)."""
    cand = jrf.recover_table1_candidate("5:6-LI", system)  # type: ignore[arg-type]
    assert cand.period_over_2pi == pytest.approx(6.0, abs=1e-3)


# ---------------------------------------------------------------------------
# (4) Independent cross-check: barden_stability vs _planar_floquet agree for
# every LARGE-eigenvalue candidate (feedback_orbit_closure_discipline).
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("label", ["3:4-LO", "5:6-LO", "5:6-NO"])
def test_large_eigenvalue_candidates_barden_matches_planar_floquet(
    label: str, system: object
) -> None:
    """For |lambda| >> 1, barden_stability and _planar_floquet are two
    independently-derived eigenvalue extractions (Barden 1994 half-period
    identity vs a direct full-period monodromy eigendecomposition) that must
    agree closely if the classification is trustworthy."""
    cand = jrf.recover_table1_candidate(label, system)  # type: ignore[arg-type]
    rel = abs(cand.max_eigenvalue - cand.planar_floquet_eigenvalue) / cand.max_eigenvalue
    assert rel < 1e-4, f"{label}: barden={cand.max_eigenvalue} pf={cand.planar_floquet_eigenvalue}"


# ---------------------------------------------------------------------------
# (5) Methodological regression: _planar_floquet's largest-magnitude
# heuristic is demonstrably unreliable near |lambda| ~ 1 (module docstring's
# "HONEST FINDING") -- barden_stability correctly reveals a genuinely
# COMPLEX (marginally/neutrally stable) nontrivial pair that
# ``_planar_floquet`` mis-reports as a trivial real ~1.0.
# ---------------------------------------------------------------------------


def test_planar_floquet_degenerate_eigenvalue_pitfall_is_real(system: object) -> None:
    orbit = cp.correct_symmetric_fixed_jacobi(
        system,  # type: ignore[arg-type]
        -0.970806,
        jrf.ANDERSON_LO_C_FLYBY,
        37.6975,
        ydot0_sign=1.0,
        half_crossings=6,
        tol=1e-11,
    )
    assert orbit.converged
    lam_barden, max_eig, is_real, lam_pf = jrf._classify(system, orbit)  # type: ignore[arg-type]
    # _planar_floquet reports an essentially-trivial real ~1.0 ...
    assert lam_pf == pytest.approx(1.0, abs=1e-3)
    # ... but barden_stability correctly reveals the TRUE nontrivial pair is
    # complex (a unit-modulus, marginally-stable pair) -- NOT the same value.
    assert not is_real
    assert abs(lam_barden.imag) > 0.05
    assert max_eig == pytest.approx(1.0, abs=1e-3)  # |lambda| for the complex pair


# ---------------------------------------------------------------------------
# (6) Seed sweep sanity: the general grid+bisection tool converges SOMETHING
# nontrivial on a small, fast grid (not a Table-1 claim -- just that the tool
# works, per its own docstring's disclosed limitations at C_flyby directly).
# ---------------------------------------------------------------------------


def test_survey_candidates_finds_nontrivial_unstable_orbit(system: object) -> None:
    found = jrf.survey_candidates(
        system,  # type: ignore[arg-type]
        jacobi=3.3,
        ydot0_sign=-1.0,
        half_crossings=6,
        t_hi=50.0,
        x0_lo=1.03,
        x0_hi=1.8,
        n_grid=16,
        label_prefix="test_",
    )
    assert len(found) >= 1
    assert any(c.max_eigenvalue > 1.01 for c in found)


# ---------------------------------------------------------------------------
# (7) The honest Table-1 gate itself: exactly 1/4 rows pass at
# TABLE1_GATE_REL_TOL. This is the mandatory gate from #753's dispatch spec,
# reported family-by-family, with no fudged tolerance.
# ---------------------------------------------------------------------------


def test_table1_gate_honest_report(system: object) -> None:
    rows = {r.label: r for r in jrf.gate_report(system)}  # type: ignore[arg-type]
    assert set(rows) == set(jrf.TABLE1_TARGETS)

    # CONFIRMED: 5:6-LI recovers to well within the 1e-3 gate tolerance AND
    # its period lands on the clean q=6 multiple (dual criterion, #753/#755).
    assert rows["5:6-LI"].passed, rows["5:6-LI"]
    assert rows["5:6-LI"].rel_err < 1e-3
    assert rows["5:6-LI"].eigenvalue_confirmed
    assert rows["5:6-LI"].period_confirmed

    # NOT CONFIRMED under the dual-criterion gate (documented negative --
    # see the #755 results note). Asserted FALSE deliberately, so this stays
    # an honest, tracked negative rather than a silently-dropped claim; a
    # future fix that finds the true families should update this test, not
    # work around it.
    for label in ("3:4-LO", "5:6-LO", "5:6-NO"):
        assert not rows[label].passed, (
            f"{label} unexpectedly passed ({rows[label]}) -- if a real fix "
            "was found, update this test's assertion and the results note, "
            "don't silently leave this assert in place"
        )

    n_passed = sum(r.passed for r in rows.values())
    assert n_passed == 1

    # #755's own striking finding, encoded as a standing regression: 3:4-LO's
    # RECOVERED EIGENVALUE matches to near-machine precision (a MUCH tighter
    # match than any of #753's original near-misses, which were 2-27% off)
    # -- but its period does NOT land on the clean q=4 multiple (a real,
    # reproducible ~2% offset, not tolerance noise). This is a qualitatively
    # different, stronger candidate than a plain near-miss; kept as its own
    # assertion so a future reviewer/fix does not conflate it with the
    # 5:6-LO/5:6-NO rows (which fail on EIGENVALUE, not just period).
    row_34lo = rows["3:4-LO"]
    assert row_34lo.eigenvalue_confirmed, row_34lo
    assert row_34lo.rel_err < 1e-6, row_34lo
    assert not row_34lo.period_confirmed, row_34lo
    assert row_34lo.period_rel_err > 1e-2, row_34lo
    assert not row_34lo.passed

    # #758: 5:6-LO's eigenvalue is now ALSO confirmed to extreme precision
    # (rel_err 3.4e-7, six orders of magnitude tighter than #753's original
    # 1.98% seed) via a genuinely new, Table-2-derived seed strategy -- see
    # the module docstring's #758 update. Same pattern as 3:4-LO: eigenvalue
    # matches, period does not (real ~2.83% offset), so the strict
    # dual-criterion gate still honestly reports FAIL pending reviewer
    # judgment (see test_758_table2_seeded_candidate_* below for the full
    # evidentiary regression).
    row_56lo = rows["5:6-LO"]
    assert row_56lo.eigenvalue_confirmed, row_56lo
    assert row_56lo.rel_err < 1e-5, row_56lo
    assert not row_56lo.period_confirmed, row_56lo
    assert row_56lo.period_rel_err > 1e-2, row_56lo
    assert not row_56lo.passed

    # 5:6-NO still fails on EIGENVALUE (unlike 3:4-LO and 5:6-LO, both of
    # which now match to extreme precision) -- no dedicated targeted search
    # has been run for this row since #753's original wide sweep.
    assert not rows["5:6-NO"].eigenvalue_confirmed


# ---------------------------------------------------------------------------
# (8) Continuation-tooling demonstration: the 5:6-LI family (the confirmed,
# well-behaved candidate) continues smoothly in Jacobi constant using the
# EXISTING natural-parameter continuation gauntlet.
# ---------------------------------------------------------------------------


def test_5_6_li_continues_smoothly_toward_c_flyby(system: object) -> None:
    # Small step count deliberately: this candidate's ~130-crossing trajectory
    # makes each continuation step (STM propagation + independent-Radau
    # cross-check) expensive (~5-6s/step observed) -- 6 steps is enough to
    # demonstrate smooth, gauntlet-passing continuation without a slow test.
    c_start = jrf.ANDERSON_LO_C_FLYBY + 0.006
    branch = jrf.continue_candidate_toward_c_flyby(
        system,  # type: ignore[arg-type]
        "5:6-LI",
        c_start=c_start,
        d_jacobi=0.001,
        n_steps=6,
    )
    assert len(branch.members) >= 4, (
        f"expected the well-behaved 5:6-LI family to continue smoothly; "
        f"got {len(branch.members)} members, stop_reason={branch.stop_reason}"
    )
    last = branch.members[-1]
    assert abs(last.jacobi - jrf.ANDERSON_LO_C_FLYBY) < 0.01


# ---------------------------------------------------------------------------
# (9) #755 strategy 2: the two-body flyby-VECTOR-ROTATION seed (Anderson &
# Lo "Designing Flybys Using the Two-Body Approximations", pp.172-174,
# Fig. 2) -- distinct from and more sophisticated than the plain resonant-
# ellipse seed in (2)/(3) above. Did not itself locate a Table-1 match in
# the time available (see the #755 results note) but is a genuine,
# geometrically-verified, reusable alternative seed strategy.
# ---------------------------------------------------------------------------


def test_flyby_rotation_seed_preserves_v_infinity_magnitude() -> None:
    """The hyperbolic flyby only ROTATES V-infinity; its magnitude (the
    two-body hyperbolic excess speed relative to the secondary) must be
    identical before and after, for any turn_sign/r_periapsis."""
    seed = jrf.two_body_flyby_rotation_seed(3, 4, 5, 6, r_periapsis=0.01, turn_sign=1)
    # Reconstruct v_infinity_after from the returned rotating-frame state:
    # v_inertial = v_rot + omega x r: (vx_rot - y, vy_rot + x0); secondary's
    # own inertial velocity there is (0, 1).
    vx_inertial = seed.xdot - seed.y0
    vy_inertial = seed.ydot + seed.x0
    vinf_after = math.hypot(vx_inertial - 0.0, vy_inertial - 1.0)
    assert vinf_after == pytest.approx(abs(seed.v_infinity), rel=1e-9)


def test_flyby_rotation_seed_turn_angle_decreases_with_periapsis_radius() -> None:
    """Standard hyperbolic-flyby geometry: a MORE DISTANT closest approach
    produces a SMALLER turn angle (weaker gravitational deflection)."""
    seed_close = jrf.two_body_flyby_rotation_seed(3, 4, 5, 6, r_periapsis=0.01)
    seed_far = jrf.two_body_flyby_rotation_seed(3, 4, 5, 6, r_periapsis=0.1)
    assert 0.0 < seed_far.turn_angle < seed_close.turn_angle < math.pi


def test_flyby_rotation_seed_rejects_bad_pq() -> None:
    with pytest.raises(ValueError, match="positive integers"):
        jrf.two_body_flyby_rotation_seed(0, 4, 5, 6)
    with pytest.raises(ValueError, match="turn_sign"):
        jrf.two_body_flyby_rotation_seed(3, 4, 5, 6, turn_sign=0)


def test_flyby_rotation_seed_backs_off_from_the_singularity() -> None:
    """The exact periapsis (x=1.0 barycentric) is only ``mu`` away from the
    CRTBP's own singularity at ``1-mu`` -- the paper's own text (p.174)
    documents this exact close-approach hazard for its "crudest method" and
    fixes it by backing the patchpoints off slightly; this module does the
    same (see the function's own docstring)."""
    seed = jrf.two_body_flyby_rotation_seed(3, 4, 5, 6, safety_margin=0.01)
    assert seed.x0 == pytest.approx(0.99)


# ---------------------------------------------------------------------------
# (10) #756: relaxed-period-criterion search for 5:6-LO. The coordinating
# session's #755 reviewer ruling (3:4-LO CONFIRMED despite a 2.1% period
# offset) raised the hypothesis that #755's own 5:6-LO search implicitly
# favored period_over_2pi near q=6 while scanning. #756 redid the search
# ranking candidates PURELY by eigenvalue closeness, with period treated
# as a secondary corroboration signal. HONEST RESULT: this did not find a
# better candidate than #753's pre-existing one, and NEITHER the
# pre-existing candidate NOR any of #756's new near-misses has a plausible
# period OR a close Europa approach (unlike 3:4-LO's genuine corroboration
# on both counts) -- a well-evidenced, continued negative, not a silently
# dropped claim (see the module docstring's #756 update and the results
# note for the full search log, 159 candidates, checkpointed at
# data/found/756_jupiter_europa_5_6_lo_relaxed_period/candidates.jsonl).
# ---------------------------------------------------------------------------


def test_europa_closest_approach_confirms_3_4_lo_close_flyby(system: object) -> None:
    """3:4-LO's confirmed candidate makes a genuine close Europa flyby
    (~1641 km, matching the paper's own attributed instability mechanism,
    p.177-178) -- the positive control for :func:`europa_closest_approach`
    itself."""
    cand = jrf.recover_table1_candidate("3:4-LO", system)  # type: ignore[arg-type]
    dist = jrf.europa_closest_approach(system, cand.x0, cand.ydot0, cand.period)  # type: ignore[arg-type]
    assert dist == pytest.approx(0.002445, abs=2e-5)


@pytest.mark.parametrize("label", sorted(jrf._756_RELAXED_SEARCH_NEAR_MISSES))
def test_756_near_miss_converges_and_reproduces_recorded_eigenvalue(
    label: str, system: object
) -> None:
    """Standing regression: each #756 near-miss seed still converges to a
    tightly-residual, real-unstable orbit at its own recorded eigenvalue --
    a genuine periodic orbit, just not a 5:6-LO match (see the module's own
    per-seed comments for the recorded target values)."""
    cand = jrf.recover_756_near_miss(label, system)  # type: ignore[arg-type]
    assert cand.crossing_residual < 1e-9
    assert cand.is_real_unstable
    assert cand.jacobi == pytest.approx(jrf.ANDERSON_LO_C_FLYBY, abs=1e-8)
    # Barden vs planar_floquet cross-check (feedback_orbit_closure_discipline).
    rel = abs(cand.max_eigenvalue - cand.planar_floquet_eigenvalue) / cand.max_eigenvalue
    assert rel < 1e-4, f"{label}: barden={cand.max_eigenvalue} pf={cand.planar_floquet_eigenvalue}"


def test_756_relaxed_search_near_misses_do_not_beat_753_original_candidate(
    system: object,
) -> None:
    """The core #756 finding (at the time it ran): relaxing the
    period-proximity search criterion and ranking purely by eigenvalue
    closeness did NOT surface anything better than #753's own ORIGINAL
    5:6-LO candidate (x0=0.81360506, rel_err=1.98%). That original
    candidate has since been SUPERSEDED in ``_TABLE1_CANDIDATE_SEEDS`` by
    #758's dramatically better Table-2-seeded find (rel_err=3.4e-7, see
    the module docstring's #758 update) -- this test therefore compares
    against the historical 1.98% baseline directly (not via
    ``recover_table1_candidate("5:6-LO")``, which now returns the #758
    candidate) so it keeps testing #756's own actual finding rather than
    silently becoming vacuous."""
    target = jrf.TABLE1_TARGETS["5:6-LO"]
    historical_753_rel_err = 0.0198  # #753's original x0=0.81360506 seed
    for label in sorted(jrf._756_RELAXED_SEARCH_NEAR_MISSES):
        cand = jrf.recover_756_near_miss(label, system)  # type: ignore[arg-type]
        rel_err = abs(cand.max_eigenvalue - target) / target
        assert rel_err > historical_753_rel_err, (
            f"{label} rel_err={rel_err} unexpectedly beat #753's original "
            f"candidate's {historical_753_rel_err} -- if a real improvement "
            "was found, update this test and the #756 results note, don't "
            "just loosen this assertion"
        )


def test_756_near_misses_lack_close_europa_approach_corroboration(system: object) -> None:
    """Unlike 3:4-LO (0.00245 nondim, a genuine close flyby), NEITHER
    #753's ORIGINAL 5:6-LO candidate (x0=0.81360506, since superseded in
    ``_TABLE1_CANDIDATE_SEEDS`` by #758's find -- see the module
    docstring's #758 update) NOR any #756 near-miss gets anywhere close to
    Europa -- all stay an order of magnitude farther away, with no
    qualitative close-flyby signature to corroborate family lineage
    (contrast :func:`test_758_table2_seeded_candidate_makes_a_close_europa_approach`,
    which shows the #758 candidate DOES)."""
    historical_753_cand = jrf.converge_candidate(
        system,  # type: ignore[arg-type]
        "753-historical-5:6-LO",
        0.81360506,
        jrf.ANDERSON_LO_C_FLYBY,
        101.2145,
        ydot0_sign=1.0,
        half_crossings=2,
    )
    assert historical_753_cand is not None
    historical_753_dist = jrf.europa_closest_approach(
        system,  # type: ignore[arg-type]
        historical_753_cand.x0,
        historical_753_cand.ydot0,
        historical_753_cand.period,
    )
    assert historical_753_dist > 0.01

    for label in sorted(jrf._756_RELAXED_SEARCH_NEAR_MISSES):
        cand = jrf.recover_756_near_miss(label, system)  # type: ignore[arg-type]
        dist = jrf.europa_closest_approach(system, cand.x0, cand.ydot0, cand.period)  # type: ignore[arg-type]
        assert dist > 0.01, f"{label}: unexpectedly close Europa approach {dist}"


def test_756_near_misses_lack_plausible_period(system: object) -> None:
    """None of #756's near-miss candidates has a period even loosely near
    the clean q=6 multiple (unlike 3:4-LO's 2.1% offset) -- all are 2x+
    the naive value, a qualitatively different (much weaker) situation."""
    for label in sorted(jrf._756_RELAXED_SEARCH_NEAR_MISSES):
        cand = jrf.recover_756_near_miss(label, system)  # type: ignore[arg-type]
        period_rel_err = abs(cand.period_over_2pi - 6.0) / 6.0
        assert period_rel_err > 0.5, f"{label}: period_over_2pi={cand.period_over_2pi}"


def test_flyby_rotation_symmetric_seed_runs_quickly(system: object) -> None:
    """Propagating the (backed-off) post-flyby state forward to its next
    perpendicular crossing must complete quickly (not grind through a
    near-singular close approach -- see the function's own ``max_step``
    docstring note) and return a sensible ``(x0, jacobi)`` pair or ``None``.
    """
    seed = jrf.two_body_flyby_rotation_seed(3, 4, 5, 6, r_periapsis=0.01)
    result = jrf.flyby_rotation_symmetric_seed(
        system,  # type: ignore[arg-type]
        seed,
        t_hi=2.0 * math.pi * 6.0 * 1.25,
    )
    assert result is not None, "expected a perpendicular crossing within t_hi"
    x0, jacobi = result
    assert -2.0 < x0 < 2.0
    assert math.isfinite(jacobi)


# ---------------------------------------------------------------------------
# (11) #758: Table-2-seeded 5:6-LO candidate. A genuinely new, sourced seed
# strategy (Anderson & Lo's own Table 2 homoclinic-intersection state, p.190,
# plus their own stated ~8.0e-5 x-offset to the 5:6 orbit, pp.184/190) found
# a STRONG candidate -- eigenvalue rel_err 3.4e-7 (six orders of magnitude
# tighter than #753's original 1.98% seed), an x0 offset from Table 2's own
# point matching the paper's stated value to ~4%, a closer Europa flyby than
# the CONFIRMED 3:4-LO orbit, basin robustness across the searched window,
# and a clean independent Radau cross-check. The one thing that does NOT
# confirm is the period (real ~2.83% offset from q=6, same qualitative
# pattern the #755 reviewer ruling already accepted for 3:4-LO's own
# comparable period offset). This module does NOT unilaterally flip the
# strict dual-criterion gate -- see
# ``docs/notes/2026-07-28-758-jupiter-europa-5-6-lo-table2-seeded-search.md``
# for the full evidentiary writeup and "candidate found, reviewer judgment
# invited" framing.
# ---------------------------------------------------------------------------


def test_table2_homoclinic_state_matches_paper_p190() -> None:
    """Anderson & Lo 2011 Table 2, p.190, verbatim (x, xdot, ydot; y=0)."""
    assert jrf.TABLE2_HOMOCLINIC_X == -1.28427733
    assert jrf.TABLE2_HOMOCLINIC_XDOT == 0.00000009
    assert jrf.TABLE2_HOMOCLINIC_YDOT == 0.46372205


def test_table2_5_6_lo_x_offset_sourced_matches_paper_p184_p190() -> None:
    """Anderson & Lo 2011 pp.184/190, verbatim: "a difference in x position
    from the 5:6 orbit of approximately 8.0 x 10^-5"."""
    assert jrf.TABLE2_5_6_LO_X_OFFSET_SOURCED == 8.0e-5


def test_758_table2_seeded_candidate_converges_and_matches_recorded_values(
    system: object,
) -> None:
    """Standing regression: the hardcoded #758 seed still converges to a
    tightly-residual, real-unstable orbit at its recorded state."""
    cand = jrf.recover_758_table2_seeded_candidate(system)  # type: ignore[arg-type]
    assert cand.crossing_residual < 1e-9
    assert cand.jacobi == pytest.approx(jrf.ANDERSON_LO_C_FLYBY, abs=1e-8)
    assert cand.is_real_unstable
    assert cand.x0 == pytest.approx(-1.2842003, abs=1e-6)
    assert cand.period_over_2pi == pytest.approx(6.169686, abs=1e-4)


def test_758_table2_seeded_candidate_eigenvalue_confirms_under_table1_gate(
    system: object,
) -> None:
    """The headline #758 finding: the recovered eigenvalue matches Table 1's
    5:6-LO target FAR inside TABLE1_GATE_REL_TOL, cross-checked between
    barden_stability and _planar_floquet (feedback_orbit_closure_discipline)
    -- a much tighter match than #753's original 1.98% seed by six orders of
    magnitude."""
    cand = jrf.recover_758_table2_seeded_candidate(system)  # type: ignore[arg-type]
    target = jrf.TABLE1_TARGETS["5:6-LO"]
    rel_err = abs(cand.max_eigenvalue - target) / target
    assert rel_err < jrf.TABLE1_GATE_REL_TOL
    assert rel_err == pytest.approx(3.4e-7, abs=2e-7)
    cross_rel = abs(cand.max_eigenvalue - cand.planar_floquet_eigenvalue) / cand.max_eigenvalue
    assert cross_rel < 1e-4, f"barden={cand.max_eigenvalue} pf={cand.planar_floquet_eigenvalue}"


def test_758_table2_seeded_candidate_x0_offset_matches_paper_stated_value(
    system: object,
) -> None:
    """The recovered candidate's x0 sits close to TABLE2_HOMOCLINIC_X, and
    the offset itself matches the paper's OWN stated ~8.0e-5 value to ~4%
    relative -- independent, paper-sourced numeric corroboration beyond the
    eigenvalue match alone (not just a qualitative shape/mechanism
    inference)."""
    cand = jrf.recover_758_table2_seeded_candidate(system)  # type: ignore[arg-type]
    offset = cand.x0 - jrf.TABLE2_HOMOCLINIC_X
    assert offset == pytest.approx(7.7e-5, abs=5e-6)
    rel = abs(abs(offset) - jrf.TABLE2_5_6_LO_X_OFFSET_SOURCED) / jrf.TABLE2_5_6_LO_X_OFFSET_SOURCED
    assert rel < 0.1, f"offset {offset} vs paper's stated {jrf.TABLE2_5_6_LO_X_OFFSET_SOURCED}"


def test_758_table2_seeded_candidate_makes_a_closer_europa_approach_than_3_4_lo(
    system: object,
) -> None:
    """Corroboration signal #2: the #758 candidate makes a genuine close
    Europa flyby (~668 km) -- CLOSER than the CONFIRMED 3:4-LO orbit's own
    1641 km approach, strongly matching the paper's attributed instability
    mechanism (p.177-178). Contrast #753's original seed and every #756
    near-miss, none of which get within an order of magnitude of Europa
    (see test_756_near_misses_lack_close_europa_approach_corroboration)."""
    cand = jrf.recover_758_table2_seeded_candidate(system)  # type: ignore[arg-type]
    dist = jrf.europa_closest_approach(system, cand.x0, cand.ydot0, cand.period)  # type: ignore[arg-type]
    assert dist == pytest.approx(0.000996, abs=2e-5)

    three_four_lo = jrf.recover_table1_candidate("3:4-LO", system)  # type: ignore[arg-type]
    three_four_lo_dist = jrf.europa_closest_approach(
        system,  # type: ignore[arg-type]
        three_four_lo.x0,
        three_four_lo.ydot0,
        three_four_lo.period,
    )
    assert dist < three_four_lo_dist


def test_758_table2_seeded_candidate_period_not_confirmed(system: object) -> None:
    """Corroboration signal that does NOT confirm (honestly reported, same
    pattern the #755 reviewer ruling already accepted for 3:4-LO's own
    comparable period offset): period_over_2pi is a real, tightly-converged
    ~2.83% offset from the clean q=6 multiple -- fails
    TABLE1_PERIOD_REL_TOL. This is real signal, not tolerance noise (the
    crossing_residual is ~1e-12, see
    test_758_table2_seeded_candidate_converges_and_matches_recorded_values)."""
    cand = jrf.recover_758_table2_seeded_candidate(system)  # type: ignore[arg-type]
    period_rel_err = abs(cand.period_over_2pi - 6.0) / 6.0
    assert period_rel_err > jrf.TABLE1_PERIOD_REL_TOL
    assert period_rel_err == pytest.approx(0.0283, abs=2e-3)


def test_758_table2_seeded_candidate_independent_radau_crosscheck(system: object) -> None:
    """Independent-integrator cross-check (feedback_orbit_closure_discipline):
    re-propagating the #758 candidate's own state with Radau (a DIFFERENT
    integrator from the DOP853 the corrector uses) must re-close and
    conserve the Jacobi constant tightly."""
    import numpy as np

    cand = jrf.recover_758_table2_seeded_candidate(system)  # type: ignore[arg-type]
    state0 = np.array([cand.x0, 0.0, 0.0, 0.0, cand.ydot0, 0.0])
    orbit = cp.PeriodicOrbit(
        state0=state0,
        period=cand.period,
        jacobi=cand.jacobi,
        converged=True,
        closure_residual=cand.crossing_residual,
    )
    ok, dj = cp.crosscheck_periodic(
        system,  # type: ignore[arg-type]
        orbit,
        method="Radau",
        rtol=1e-12,
        atol=1e-12,
        closure_tol=1e-6,
        jacobi_tol=1e-6,
    )
    assert ok, f"Radau cross-check failed: jacobi drift={dj}"
    assert dj < 1e-10


def test_758_table2_seeded_candidate_supersedes_753_seed_in_shared_table(
    system: object,
) -> None:
    """``_TABLE1_CANDIDATE_SEEDS["5:6-LO"]`` now points at the #758
    candidate, not #753's original weaker seed -- ``recover_table1_candidate``
    and ``recover_758_table2_seeded_candidate`` must agree."""
    shared = jrf.recover_table1_candidate("5:6-LO", system)  # type: ignore[arg-type]
    direct = jrf.recover_758_table2_seeded_candidate(system)  # type: ignore[arg-type]
    assert shared.x0 == pytest.approx(direct.x0, abs=1e-9)
    assert shared.max_eigenvalue == pytest.approx(direct.max_eigenvalue, rel=1e-9)


def test_basin_robustness_scan_shows_dominant_basin(system: object) -> None:
    """The #758 candidate is NOT an isolated numerical fluke: a majority of
    evenly-spaced seeds across the sourced +-2e-4 window converge Newton-
    directly to the exact same point (mirrors #756's own basin-robustness
    check on the pre-existing 5:6-LO seed, and the same
    survey_candidates-bracket-scan-misses-it tooling nuance #756
    documented)."""
    results = jrf.basin_robustness_scan(
        system,  # type: ignore[arg-type]
        x0_lo=jrf.TABLE2_HOMOCLINIC_X - 2e-4,
        x0_hi=jrf.TABLE2_HOMOCLINIC_X + 2e-4,
        n_seeds=11,
        jacobi=jrf.ANDERSON_LO_C_FLYBY,
        ydot0_sign=1.0,
        half_crossings=2,
        period_guess=40.0,
    )
    n_converged = sum(1 for _, c in results if c is not None)
    assert n_converged >= 8, f"expected most seeds to converge, got {n_converged}/11"
    target_x0 = jrf.recover_758_table2_seeded_candidate(system).x0  # type: ignore[arg-type]
    n_at_target = sum(1 for _, c in results if c is not None and abs(c.x0 - target_x0) < 1e-6)
    assert n_at_target >= 6, (
        f"expected the #758 candidate's basin to dominate the window, got "
        f"{n_at_target}/11 seeds landing there"
    )


def test_survey_candidates_bracket_scan_misses_758_root_at_coarse_resolution(
    system: object,
) -> None:
    """Honest tooling-nuance regression (same one #756 documented for the
    pre-existing 5:6-LO seed): survey_candidates' own bracket/sign-flip scan
    does NOT detect the #758 root at a modest grid resolution in this
    narrow window -- direct Newton convergence (basin_robustness_scan /
    recover_758_table2_seeded_candidate) is the ground truth, not this
    scan. Documented so a future reader does not conclude the region is
    empty just because this particular tool found nothing here."""
    found = jrf.survey_candidates(
        system,  # type: ignore[arg-type]
        jacobi=jrf.ANDERSON_LO_C_FLYBY,
        ydot0_sign=1.0,
        half_crossings=2,
        t_hi=40.0,
        x0_lo=jrf.TABLE2_HOMOCLINIC_X - 2e-4,
        x0_hi=jrf.TABLE2_HOMOCLINIC_X + 2e-4,
        n_grid=400,
        label_prefix="758_scan_test_",
    )
    assert len(found) == 0


# ---------------------------------------------------------------------------
# (13) #761: the confirmed 3:4-LO and Kumar et al. 2021's own "arbitrarily
# chosen" 3:4 seed (C=3.0041) are the SAME continuous family -- established
# by a clean, fold-free, gauntlet-validated continuation between the two
# paper-anchored endpoints, corroborated by Kumar 2021's own published
# closest-approach number. See the module docstring's #761 update and
# docs/notes/2026-07-29-761-torus-seed-continuation-tractability.md.
# ---------------------------------------------------------------------------


def test_761_kumar_constants_match_paper_p8() -> None:
    """Kumar et al. 2021 (AAS 21-651) p.8, verbatim: Jacobi constant 3.0041
    ("arbitrarily chosen") and a 22052 km Europa closest approach at mu3=0."""
    assert jrf.KUMAR_2021_C == 3.0041
    assert jrf.KUMAR_2021_CLOSEST_APPROACH_KM == 22052.0


def test_761_continuation_connects_3_4_lo_to_kumar_2021_seed(system: object) -> None:
    """The #761 finding as a standing regression: the confirmed 3:4-LO
    continues smoothly (no fold, no branch/topology jump, every member
    through the full continue_family gauntlet) from C_flyby up to Kumar
    2021's own C=3.0041, and the endpoint member reproduces the paper's own
    published 22,052 km Europa closest approach (the mu3=0 PCRTBP value,
    exactly this model) to <0.2% relative -- independent, published-number
    evidence that the two papers' points lie on one continuous family."""
    branch, endpoint = jrf.continue_34lo_to_kumar_c(system)  # type: ignore[arg-type]

    # The walk reached the Jacobi bound (not a fold/no-converge/jump stop),
    # with a genuinely resolved path, monotone in x0 (smooth family, no
    # back-tracking or basin hopping along the way).
    assert len(branch.members) >= 20
    x0s = [m.x0 for m in branch.members]
    assert all(b > a for a, b in itertools.pairwise(x0s))

    # Endpoint: a tightly-converged, genuinely saddle-unstable member at
    # exactly Kumar's stated Jacobi constant (values recorded #761;
    # |lambda|~54.59 -- far weaker than 3:4-LO's own 1036 at C_flyby, the
    # smooth |lambda| decay along the branch is part of the finding).
    assert endpoint.jacobi == pytest.approx(jrf.KUMAR_2021_C, abs=1e-9)
    assert endpoint.crossing_residual < 1e-10
    assert endpoint.is_real_unstable
    assert endpoint.max_eigenvalue == pytest.approx(54.5898, rel=1e-3)
    assert endpoint.x0 == pytest.approx(-1.3852484456, abs=1e-8)

    # Sourced golden corroboration (Kumar 2021 p.8's own number; expected
    # side is published, never our own computation -- reproduced #761 at
    # 22,035.8 km, 0.073% relative).
    import cyclerfinder.core.cr3bp as cr3bp

    assert isinstance(system, cr3bp.CR3BPSystem)
    ca_km = (
        jrf.europa_closest_approach(system, endpoint.x0, endpoint.ydot0, endpoint.period)
        * system.l_km
    )
    assert ca_km == pytest.approx(jrf.KUMAR_2021_CLOSEST_APPROACH_KM, rel=2e-3)
