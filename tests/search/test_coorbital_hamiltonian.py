"""Tests for :mod:`cyclerfinder.search.coorbital_hamiltonian` (task #666).

The mandatory validation gate comes first (`test_l4_libration_frequency_matches_classical`):
this module's numerical-averaging machinery must reproduce a known, independently-derived
closed-form result (the classical CR3BP L4/L5 tadpole libration frequency) before any of its
QS/HS/tadpole classification is trusted for anything downstream (per
`[[feedback_verify_gauntlet_with_positive_control]]`).
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from cyclerfinder.search.coorbital_hamiltonian import (
    averaged_disturbing_function,
    build_phase_portrait,
    classify_point,
    ecliptic_longitude_from_elements,
    hamiltonian,
    hill_radius_delta,
    kepler_eccentric_anomaly,
    kepler_energy_rotating_frame,
    true_anomaly_from_mean,
)

MU_SUN_EARTH = 3.0034805950690393e-06  # cyclerfinder.core.cr3bp.cr3bp_system("Sun", "Earth").mu


# ---------------------------------------------------------------------------
# Kepler solver sanity
# ---------------------------------------------------------------------------


def test_kepler_circular_orbit_true_anomaly_equals_mean_anomaly() -> None:
    """At e=0, E=M and nu=M identically (no eccentricity, no anomaly offset)."""
    m = np.linspace(0.0, 2.0 * math.pi, 37)
    ecc_anom = kepler_eccentric_anomaly(m, 0.0)
    nu = true_anomaly_from_mean(m, 0.0)
    np.testing.assert_allclose(ecc_anom, np.mod(m, 2.0 * math.pi), atol=1e-12)
    np.testing.assert_allclose(nu, np.mod(m, 2.0 * math.pi), atol=1e-12)


@pytest.mark.parametrize("e", [0.05, 0.10269, 0.3, 0.6])
def test_kepler_equation_residual_is_solved(e: float) -> None:
    """``E - e*sin(E) - M == 0`` at the solver's output, for a spread of M."""
    m = np.linspace(0.1, 2.0 * math.pi - 0.1, 23)
    ecc_anom = kepler_eccentric_anomaly(m, e)
    residual = ecc_anom - e * np.sin(ecc_anom) - m
    np.testing.assert_allclose(residual, 0.0, atol=1e-10)


def test_true_anomaly_periapsis_and_apoapsis() -> None:
    """M=0 -> nu=0 (periapsis); M=pi -> nu=pi (apoapsis), any e."""
    for e in (0.1, 0.5, 0.8):
        nu0 = float(true_anomaly_from_mean(np.array([0.0]), e)[0])
        nu_pi = float(true_anomaly_from_mean(np.array([math.pi]), e)[0])
        assert abs(nu0) < 1e-9
        assert abs(nu_pi - math.pi) < 1e-9


# ---------------------------------------------------------------------------
# Kepler energy term
# ---------------------------------------------------------------------------


def test_kepler_energy_critical_at_resonance() -> None:
    """H0'(a=1) == 0 -- a = a_planet is always a critical point (the resonance itself)."""
    h = 1e-6
    d_plus = kepler_energy_rotating_frame(1.0 + h)
    d_minus = kepler_energy_rotating_frame(1.0 - h)
    derivative = (d_plus - d_minus) / (2 * h)
    assert abs(derivative) < 1e-6


def test_kepler_energy_quadratic_coefficient_is_minus_three_eighths() -> None:
    """``H0(1+delta) - H0(1) ~= -(3/8) delta**2`` for small delta (module docstring derivation)."""
    delta = 1e-3
    h0_1 = kepler_energy_rotating_frame(1.0)
    h0_delta = kepler_energy_rotating_frame(1.0 + delta)
    expected = -(3.0 / 8.0) * delta**2
    actual = h0_delta - h0_1
    # O(delta**3) truncation residual at delta=1e-3 is ~1e-9 relative to the
    # ~3.75e-7 leading term -- i.e. a real ~0.1% deviation from the pure
    # quadratic, not numerical noise; rel=2e-2 comfortably separates that
    # from an actual coefficient bug (e.g. a sign or factor-of-2 error).
    assert actual == pytest.approx(expected, rel=2e-2)


# ---------------------------------------------------------------------------
# MANDATORY validation gate: L4/L5 tadpole libration frequency
# ---------------------------------------------------------------------------


def test_l4_libration_frequency_matches_classical_closed_form() -> None:
    """omega_lib(L4, e->0) must match n_planet * sqrt(27/4 * mu * (1-mu)).

    This is THIS MODULE'S positive control on its own numerical machinery,
    independent of any real-object data: the classical linear-stability
    result for the CR3BP triangular Lagrange points is completely standard
    (e.g. Murray & Dermott Ch. 3) and is NOT itself derived via this
    module's averaging -- if the numerical averaged Hamiltonian's curvature
    at sigma=60deg, e~0, delta=0 does not reproduce this frequency, the
    averaging or the reduced-Hamiltonian kinetic term has a bug, and
    nothing built on top of it (the Kamo'oalewa positive control included)
    should be trusted.
    """
    mu = MU_SUN_EARTH
    e = 0.0
    sigma0 = math.radians(60.0)
    h = 1e-4
    vals = {
        ds: averaged_disturbing_function(1.0, e, sigma0 + ds, mu, n_theta=400, n_varpi=8)
        for ds in (-2 * h, -h, 0.0, h, 2 * h)
    }
    # 5-point central second derivative.
    second_deriv = (-vals[-2 * h] + 16 * vals[-h] - 30 * vals[0.0] + 16 * vals[h] - vals[2 * h]) / (
        12 * h * h
    )
    omega_model = math.sqrt(3.0 * abs(second_deriv))
    omega_classical = math.sqrt(27.0 / 4.0 * mu * (1.0 - mu))
    assert omega_model == pytest.approx(omega_classical, rel=2e-3)


# ---------------------------------------------------------------------------
# Symmetry (L4 <-> L5 mirror)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("sigma_deg", [10.0, 30.0, 60.0, 90.0, 150.0, 180.0])
def test_disturbing_function_mirror_symmetric(sigma_deg: float) -> None:
    """R_avg(sigma, e) == R_avg(360 - sigma, e): mirror symmetry of the coplanar,
    circular-planet restricted problem (L4 <-> L5)."""
    mu = MU_SUN_EARTH
    e = 0.10269
    v1 = averaged_disturbing_function(1.0, e, math.radians(sigma_deg), mu, n_theta=150, n_varpi=150)
    v2 = averaged_disturbing_function(
        1.0, e, math.radians(360.0 - sigma_deg), mu, n_theta=150, n_varpi=150
    )
    assert v1 == pytest.approx(v2, rel=1e-9, abs=1e-15)


def test_hamiltonian_zero_referenced_at_delta_zero() -> None:
    """``H(sigma, delta=0, e) == -R_avg(sigma, e, delta=0)`` exactly (H0 term
    cancels at delta=0 by construction; the MINUS sign is the standard
    disturbing-function convention -- see :func:`hamiltonian`'s own comment
    for the derivation and the L4-stability check that caught it)."""
    mu = MU_SUN_EARTH
    e = 0.05
    sigma = math.radians(45.0)
    h_val = hamiltonian(sigma, 0.0, e, mu, n_theta=100, n_varpi=100)
    r_avg = averaged_disturbing_function(1.0, e, sigma, mu, n_theta=100, n_varpi=100)
    assert h_val == pytest.approx(-r_avg, rel=1e-9)


# ---------------------------------------------------------------------------
# QS island existence (qualitative, known-theory check): QS requires
# "enough" eccentricity (e.g. Namouni 1999) -- at very low e, sigma=0 is NOT
# an isolated stable island; at Kamo'oalewa's real e it clearly is.
# ---------------------------------------------------------------------------


def test_hill_radius_delta_scale_sun_earth() -> None:
    """Sanity: the Sun-Earth Hill radius is ~0.01 (nondim), matching the well-known
    ~1.5e6 km figure."""
    hr = hill_radius_delta(MU_SUN_EARTH)
    assert hr == pytest.approx(0.01, rel=0.05)


def test_quasi_satellite_island_requires_eccentricity() -> None:
    """At very low e (0.001), sigma=0/delta=0 is NOT an isolated stable island
    (it sits on the same connected component as the wide non-librating
    region) -- but at Kamo'oalewa's real e=0.10269 it clearly is. This is the
    standard QS-needs-eccentricity result (e.g. Namouni 1999); the exact
    critical e is resolution/mu-dependent and not asserted here, only the
    qualitative low-e/high-e contrast."""
    mu = MU_SUN_EARTH
    hr = hill_radius_delta(mu)
    portrait_low_e = build_phase_portrait(
        0.001, mu, delta_max=6 * hr, n_sigma=91, n_delta=61, n_theta=80, n_varpi=80
    )
    regime_low_e, _ = classify_point(portrait_low_e, 0.0, 0.0)
    assert regime_low_e != "quasi_satellite"

    portrait_kamo_e = build_phase_portrait(
        0.10269, mu, delta_max=6 * hr, n_sigma=91, n_delta=61, n_theta=80, n_varpi=80
    )
    regime_kamo_e, _ = classify_point(portrait_kamo_e, 0.0, 0.0)
    assert regime_kamo_e == "quasi_satellite"


# ---------------------------------------------------------------------------
# Ecliptic-longitude projection helper
# ---------------------------------------------------------------------------


def test_ecliptic_longitude_reduces_to_planar_case_at_zero_inclination() -> None:
    """At i=0, the projected longitude is just Omega + omega + nu (no projection correction)."""
    e = 0.1
    raan = math.radians(40.0)
    argp = math.radians(70.0)
    m = math.radians(15.0)
    lam = ecliptic_longitude_from_elements(1.0, e, 0.0, raan, argp, m)
    nu = float(true_anomaly_from_mean(np.array([m]), e)[0])
    expected = math.atan2(math.sin(raan + argp + nu), math.cos(raan + argp + nu))
    assert math.atan2(math.sin(lam), math.cos(lam)) == pytest.approx(expected, abs=1e-9)


def test_ecliptic_longitude_at_node_crossings_matches_raan() -> None:
    """At u = argp+nu = 0 (ascending node crossing), projected longitude == Omega exactly,
    for any inclination (the projection formula's own defining property)."""
    raan = math.radians(65.79)
    inc = math.radians(20.0)
    e = 0.1
    # Choose mean anomaly so nu = -argp (u=0): argp=0 for simplicity here.
    lam = ecliptic_longitude_from_elements(1.0, e, inc, raan, 0.0, 0.0)
    assert math.atan2(math.sin(lam - raan), math.cos(lam - raan)) == pytest.approx(0.0, abs=1e-9)


# ---------------------------------------------------------------------------
# classify_point: known qualitative regimes
# ---------------------------------------------------------------------------


def test_classify_point_low_e_near_l4_is_tadpole() -> None:
    """A small-amplitude point right at L4 (sigma=60deg, delta=0), low e, is NOT the QS island
    and does not reach sigma=0 -- the classic tadpole regime."""
    mu = MU_SUN_EARTH
    hr = hill_radius_delta(mu)
    portrait = build_phase_portrait(
        0.02, mu, delta_max=6 * hr, n_sigma=121, n_delta=81, n_theta=90, n_varpi=90
    )
    regime, diag = classify_point(portrait, math.radians(60.0), 0.0)
    assert regime == "tadpole"
    assert diag["reaches_sigma0"] == 0.0


def test_classify_point_kamooalewa_real_point_is_quasi_satellite_and_narrow() -> None:
    """Kamo'oalewa's real (sigma~-3.3deg, delta~+9.4e-4, e=0.10269) point classifies as
    ``quasi_satellite`` in a NARROW (< 15 deg wide) libration island around sigma=0 --
    directly matching its published CURRENT QS state (de la Fuente Marcos & de la Fuente
    Marcos 2016, MNRAS 462:3441). The narrowness is itself notable: a small perturbation
    (from the physics this averaged, fixed-e model does not capture -- Earth's own orbital
    eccentricity, other planets) is plausibly enough to carry the real trajectory across
    this island's boundary, consistent with the real object's observed QS<->HS switching.
    """
    mu = MU_SUN_EARTH
    e = 0.10269
    sigma_rad = math.radians(-3.3245)
    delta = 0.00094
    hr = hill_radius_delta(mu)
    portrait = build_phase_portrait(
        e, mu, delta_max=6 * hr, n_sigma=181, n_delta=121, n_theta=120, n_varpi=120
    )
    regime, diag = classify_point(portrait, sigma_rad, delta)
    assert regime == "quasi_satellite"
    assert diag["sigma_span_deg"] < 15.0


def test_classify_point_2010so16_real_point_is_horseshoe() -> None:
    """(419624) 2010 SO16 is a well-known, long-lived Earth HORSESHOE co-orbital (Christou
    & Asher 2011, MNRAS 414:2965 -- "A long-lived horseshoe companion to the Earth"), a
    genuinely different regime from Kamo'oalewa's QS. Its real (sigma, delta, e) at epoch
    2017-Sep-04 (JD 2458000.5; a=1.0028 au, e=0.0754, i=14.52deg, RAAN=40.397deg,
    argp=108.99deg, M=173.30deg -- JPL SBDB via Wikipedia) should classify as a WIDE
    (``horseshoe``) region, not an isolated QS island -- an independent cross-check of this
    module's classifier on a second, previously-published, differently-labelled real object.
    """
    mu = MU_SUN_EARTH
    e = 0.0754
    sigma_rad = math.radians(-17.234574037658483)  # see run_666 script for the derivation
    delta = 0.00279738269883123
    hr = hill_radius_delta(mu)
    portrait = build_phase_portrait(
        e, mu, delta_max=6 * hr, n_sigma=121, n_delta=81, n_theta=100, n_varpi=100
    )
    regime, diag = classify_point(portrait, sigma_rad, delta)
    assert regime == "horseshoe"
    assert diag["sigma_span_deg"] > 180.0
