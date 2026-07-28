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

    # 5:6-LO and 5:6-NO fail on EIGENVALUE too (unlike 3:4-LO) -- #755 spent
    # substantial additional targeted search on 5:6-LO specifically and
    # found nothing closer than #753's original candidate.
    assert not rows["5:6-LO"].eigenvalue_confirmed
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
