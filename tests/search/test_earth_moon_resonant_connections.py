"""Tests for `#783`: Earth-Moon He1 homoclinic-connection stage
(Barrabes-Mondelo-Olle 2009 continuation-of-homoclinic-connections method).

Honest headline (see the module's own docstring and
``docs/notes/2026-08-08-783-earth-moon-homoclinic-connection.md`` for the full account):
this is a REPRODUCTION ATTEMPT of Casoliva 2010's own published He1 connection that resulted
in a CLEAN NEGATIVE on the connection itself -- both correctors below make genuine, measured
residual progress but do not converge within this task's own compute budget (a "conditioning
wall" diagnosed quantitatively in the module docstring). What DOES reproduce cleanly: the
periodic-orbit anchor itself (already `#780`'s own result, re-confirmed here), the monodromy
eigenvector (Table 4's own printed ``V^u``, OCR-corrected this task), and Table 6 row 8's own
printed LEO-rendezvous delta-v self-consistency. Newton-corrector tests below use SMALL
iteration/segment budgets deliberately (this module's own full investigation used much larger
budgets -- see the results note) -- they assert genuine residual REDUCTION and the expected
``converged=False`` outcome, not a fabricated full convergence.
"""

from __future__ import annotations

import numpy as np
import pytest

import cyclerfinder.core.cr3bp as cr3bp
import cyclerfinder.search.earth_moon_resonant_connections as emc
import cyclerfinder.search.earth_moon_resonant_families as emf


@pytest.fixture(scope="module")
def system() -> cr3bp.CR3BPSystem:
    return emf.earth_moon_system()


@pytest.fixture(scope="module")
def seed(system: cr3bp.CR3BPSystem) -> emc.HE1ConnectionSeed:
    return emc.build_he1_connection_seed(system)


# ---------------------------------------------------------------------------
# (1) Table 4 continuation-variable constants, OCR-corrected this task.
# ---------------------------------------------------------------------------


def test_he1_vu_raw_pxpy_matches_image_read_signs() -> None:
    """Casoliva 2010 Table 4's own printed V^u, read directly off PDF page 13
    (module docstring's "OCR SIGN HAZARD") -- components 1 and 3 (0-indexed 0, 2)
    positive, components 2 and 4 (0-indexed 1, 3) negative."""
    dx, dy, dpx, dpy = emc.HE1_VU_RAW_PXPY
    assert dx > 0.0
    assert dy < 0.0
    assert dpx > 0.0
    assert dpy < 0.0
    assert abs(dx - 0.04936117474608325) < 1e-15
    assert abs(dy - (-0.5669396006141868)) < 1e-13
    assert abs(dpx - 0.8204920465924677) < 1e-13
    assert abs(dpy - (-0.05418270168269977)) < 1e-15


def test_he1_vs_raw_pxpy_footnote_sign_rule() -> None:
    """Table 4's own footnote: V^s obtained from V^u by flipping the signs of its own
    2nd and 3rd (x,y,px,py) components -- i.e. y and p_x, NOT x or p_y."""
    vu = emc.HE1_VU_RAW_PXPY
    vs = emc.HE1_VS_RAW_PXPY
    assert vs[0] == vu[0]
    assert vs[1] == -vu[1]
    assert vs[2] == -vu[2]
    assert vs[3] == vu[3]


def test_he1_theta_and_time_signs() -> None:
    """Barrabes-Mondelo-Olle's own stated convention: theta^s = -theta^u; T, T^u > 0 and
    T^s < 0 (module docstring -- their OWN words, not a guess this module made)."""
    assert emc.HE1_THETA_S == -emc.HE1_THETA_U
    assert emc.HE1_T_U > 0.0
    assert emc.HE1_T_S < 0.0


def test_he1_vu_raw_pxpy_is_unit_norm_in_casolivas_own_basis() -> None:
    """Table 4's own equation 'V^u_2 - 1 = 0' is stated in Casoliva's own (x, y, p_x, p_y)
    canonical-momentum basis (where the printed vector is defined) -- confirmed directly
    here. NOTE: the (x, y, vx, vy) basis change (``_pxpy_to_vxvy``) is NOT orthogonal, so
    :data:`HE1_VU_PROJECT_PLANAR`'s own norm need NOT be 1 (and empirically is not, ~0.63) --
    this is expected, not a bug."""
    norm = float(np.linalg.norm(emc.HE1_VU_RAW_PXPY))
    assert abs(norm - 1.0) < 1e-9


# ---------------------------------------------------------------------------
# (2) Table 5 / Table 6: vendored, image-read (module docstring), structural checks.
# ---------------------------------------------------------------------------


def test_table5_has_4_rows_spanning_half_period_symmetry() -> None:
    assert len(emc.TABLE5_ROWS) == 4
    r1, r2, _r3, r4 = emc.TABLE5_ROWS
    # Mirror symmetry about the half-period point (label 3): label 2 and label 4 share
    # (r, v, a, e), with opposite-sign omega.
    assert r2.r_km == r4.r_km
    assert r2.v_kms == r4.v_kms
    assert r2.a_km == r4.a_km
    assert r2.e == r4.e
    assert r2.omega_deg == -r4.omega_deg
    # Label 1 (periselene, t=0) has the smallest r among the 4 -- Casoliva's own text:
    # "periselene 1 for the Lyapunov orbit is 6331 km."
    assert r1.r_km == min(row.r_km for row in emc.TABLE5_ROWS)
    assert abs(r1.r_km - 6331.184) < 1e-9


def test_table6_has_19_rows_closing_on_itself() -> None:
    """Label 1 and label 19 (periselene 1 and periselene 19) share identical printed
    (r, v, a, e) -- Casoliva's own text: 'the time needed... from periselene 1 to
    periselene 19 is 113.6319 days,' the connection's own homoclinic closure."""
    assert len(emc.TABLE6_ROWS) == 19
    r1 = emc.TABLE6_ROWS[0]
    r19 = emc.TABLE6_ROWS[-1]
    assert r1.label == 1
    assert r19.label == 19
    assert r1.body == r19.body == "Moon"
    assert r1.r_km == r19.r_km
    assert r1.v_kms == r19.v_kms
    assert r1.a_km == r19.a_km
    assert r1.e == r19.e
    assert abs(r19.t_flight_days - 113.632) < 1e-6


def test_table6_body_labels_match_moon_then_earth_then_moon() -> None:
    bodies = [row.body for row in emc.TABLE6_ROWS]
    assert bodies[:5] == ["Moon"] * 5
    assert bodies[5:14] == ["Earth"] * 9
    assert bodies[14:] == ["Moon"] * 5


def test_table6_earth_rows_are_bound_ellipses_moon_rows_are_hyperbolic() -> None:
    """Earth-relative rows: a > 0, e < 1 (bound ellipse). Moon-relative rows: a < 0, e > 1
    (hyperbolic, two-body-relative-to-Moon) -- Casoliva's own text: 'the periselene of a
    hyperbolic orbit of eccentricity 1.126.'"""
    for row in emc.TABLE6_ROWS:
        if row.body == "Earth":
            assert row.a_km > 0.0
            assert row.e < 1.0
        else:
            assert row.a_km < 0.0
            assert row.e > 1.0


def test_table6_row8_leo_dv_matches_casolivas_own_printed_717_5_mps() -> None:
    """Text-stated (2010, p. 1635): 717.5 m/s at r=67,869 km. Computed DIRECTLY from
    Table 6 row 8's own printed (r, v), NOT from this module's own corrector (module
    docstring's honesty discipline)."""
    dv_mps, v_circ = emc.table6_row8_leo_dv_check()
    assert abs(dv_mps - emf.HE1_LEO_DV_MPS_2010) < 1.0
    assert v_circ > 0.0


# ---------------------------------------------------------------------------
# (3) He1 connection seed: the already-converged periodic orbit + eigenpair.
# ---------------------------------------------------------------------------


def test_build_he1_connection_seed_reproduces_780(seed: emc.HE1ConnectionSeed) -> None:
    """Reuses `#780`'s own recover_he1_lyapunov -- lam_u should match
    emf.HE1_TARGET_LAMBDA_U to the same tight tolerance `#780`'s own gate already
    established (5.3e-8 relative)."""
    rel_err = abs(seed.lam_u - emf.HE1_TARGET_LAMBDA_U) / emf.HE1_TARGET_LAMBDA_U
    assert rel_err < 1e-6
    # Near machine precision -- NOT the source of this module's own non-convergence
    # (module docstring's "THE CONDITIONING WALL").
    assert seed.closure_residual < 1e-9


def test_eigenvector_reproduction_check_passes(
    system: cr3bp.CR3BPSystem, seed: emc.HE1ConnectionSeed
) -> None:
    """The ONE piece of Table 4 this module DOES reproduce cleanly (module docstring) --
    this module's own independently-recovered monodromy eigenvector matches Table 4's own
    printed, OCR-corrected V^u to a tight tolerance."""
    check = emc.eigenvector_reproduction_check(system, seed)
    assert check.passed
    assert check.cos_angle > 1.0 - 1e-4
    assert check.rel_err_unit < 1e-3
    assert check.monodromy_eig_spread_rel < 1e-6


# ---------------------------------------------------------------------------
# (4) analytic_manifold_seed: Barrabes-Mondelo-Olle Eq. (4), algebraic sanity.
# ---------------------------------------------------------------------------


def test_manifold_seed_at_theta_zero_is_state_plus_xi_v(
    system: cr3bp.CR3BPSystem, seed: emc.HE1ConnectionSeed
) -> None:
    """At theta=0, Phi_0 = Identity and lam^0 = 1, so
    psi(0, xi) = z0 + xi * v0 exactly (module docstring's Eq. (4))."""
    state6, stm6 = emc.analytic_manifold_seed(system, seed, 0.0, direction="unstable")
    expected = seed.state0.copy()
    v6 = np.array([seed.v_u[0], seed.v_u[1], 0.0, seed.v_u[2], seed.v_u[3], 0.0])
    expected = expected + emc.HE1_XI0 * v6
    assert np.allclose(state6, expected, atol=1e-13)
    assert np.allclose(stm6, np.eye(6), atol=1e-9)


def test_manifold_seed_at_theta_2pi_returns_to_similar_offset(
    system: cr3bp.CR3BPSystem, seed: emc.HE1ConnectionSeed
) -> None:
    """Eq. (5)'s own periodic-renormalisation property: psi(2*pi, xi) should closely
    match psi(0, xi) again (both project onto the SAME point + comparable offset scale),
    since Phi_T(z0)=z0 and lam^(-1) exactly cancels the eigenvector's own natural growth
    over one full period."""
    s0, _ = emc.analytic_manifold_seed(system, seed, 0.0, direction="unstable")
    s2pi, _ = emc.analytic_manifold_seed(system, seed, 2.0 * np.pi, direction="unstable")
    assert np.linalg.norm(s2pi - s0) < 1e-6


def test_manifold_seed_direction_must_be_valid(
    system: cr3bp.CR3BPSystem, seed: emc.HE1ConnectionSeed
) -> None:
    with pytest.raises(ValueError):
        emc.analytic_manifold_seed(system, seed, 0.0, direction="sideways")


# ---------------------------------------------------------------------------
# (5) Newton correctors: genuine residual reduction, expected honest non-convergence.
# Deliberately SMALL budgets (this module's own full investigation used much larger
# ones -- see the results note); these assert the ALGORITHM works, not full closure.
# ---------------------------------------------------------------------------


def test_single_shoot_reduces_residual_from_casolivas_own_guess(
    system: cr3bp.CR3BPSystem, seed: emc.HE1ConnectionSeed
) -> None:
    result = emc.correct_connection_single_shoot(system, seed, max_iter=5)
    assert result.method == "single_shoot"
    assert abs(result.initial_residual - 0.393663) < 1e-3
    assert result.residual < result.initial_residual
    # NOT asserting convergence -- module docstring's own honest headline result.
    assert result.n_iter > 0


def test_multi_shoot_reduces_residual_and_reports_segment_diagnostics(
    system: cr3bp.CR3BPSystem, seed: emc.HE1ConnectionSeed
) -> None:
    result, stm_norms = emc.correct_connection_multi_shoot(system, seed, n_u=4, n_s=4, max_iter=5)
    assert result.method == "multi_shoot"
    assert result.residual < result.initial_residual
    assert len(stm_norms) == (4 - 1) + (4 - 1) + 2  # interior continuities + final match legs
    assert all(sn > 0.0 for sn in stm_norms)


def test_connection_residual_single_shoot_matches_corrector_initial_residual(
    system: cr3bp.CR3BPSystem, seed: emc.HE1ConnectionSeed
) -> None:
    """Standalone residual function agrees with the corrector's own internal evaluation
    at the same (theta_u, T_u, theta_s, T_s)."""
    resid = emc.connection_residual_single_shoot(
        system, seed, emc.HE1_THETA_U, emc.HE1_T_U, emc.HE1_THETA_S, emc.HE1_T_S
    )
    result = emc.correct_connection_single_shoot(system, seed, max_iter=0)
    assert abs(float(np.linalg.norm(resid)) - result.initial_residual) < 1e-6


def test_multi_shoot_rejects_invalid_segment_counts(
    system: cr3bp.CR3BPSystem, seed: emc.HE1ConnectionSeed
) -> None:
    with pytest.raises(ValueError):
        emc.correct_connection_multi_shoot(system, seed, n_u=0, n_s=4, max_iter=1)
