"""Tests for the Uranus-Umbriel-Titania CCR4BP system constructor (`#701`).

Mirrors `#695`'s own ``test_ccr4bp_io_europa.py`` / `#696`'s
``test_ccr4bp_io_ganymede.py`` sourcing/reduction discipline, applied to the
new ``uranus_umbriel_titania_default()`` constructor
(``core.ccr4bp_umbriel_titania``) instead of building anything new in
``core.ccr4bp`` itself (that module is reused UNMODIFIED, per `#701`'s scope).
No value computed by our own code sits on the EXPECTED side of a structural
assertion -- sourced constants (``core.satellites`` / `#693`'s own screening
table) or the already-proven ``core.ccr4bp``/``core.cr3bp`` reduction
structure only.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

import cyclerfinder.core.ccr4bp as ccr4bp
import cyclerfinder.core.ccr4bp_umbriel_titania as ut
import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.core.satellites import PRIMARIES, SATELLITES

_SAMPLE_STATE = np.array([0.6, 0.2, 0.0, 0.05, 0.35, 0.0], dtype=np.float64)


# ---------------------------------------------------------------------------
# Sourced-constant cross-checks.
# ---------------------------------------------------------------------------


def test_constants_cross_check_against_satellites_registry() -> None:
    """mu/mu_gan/a_gan trace exactly to core.satellites (no independent re-sourcing)."""
    sys_c = ut.uranus_umbriel_titania_default()
    gm_u = PRIMARIES["Uranus"]
    gm_umbriel = SATELLITES["Umbriel"].mu_km3_s2
    gm_titania = SATELLITES["Titania"].mu_km3_s2
    sma_umbriel = SATELLITES["Umbriel"].sma_km
    sma_titania = SATELLITES["Titania"].sma_km
    denom = gm_u + gm_umbriel
    assert sys_c.mu == pytest.approx(gm_umbriel / denom, rel=0.0, abs=1e-18)
    assert sys_c.mu_gan == pytest.approx(gm_titania / denom, rel=0.0, abs=1e-18)
    assert sys_c.a_gan == pytest.approx(sma_titania / sma_umbriel, rel=0.0, abs=1e-15)


def test_cross_check_against_693_screening_numbers() -> None:
    """Cross-check against #693's own independently-derived screening-pass
    table (mu_base=1.469e-5, mu_pert=3.916e-5, P ratio 2.101)."""
    sys_c = ut.uranus_umbriel_titania_default()
    assert sys_c.mu == pytest.approx(1.469e-5, abs=5e-8)
    assert sys_c.mu_gan == pytest.approx(3.916e-5, abs=5e-8)
    # Physical period ratio T_titania / T_umbriel == n_umbriel/n_titania == (1+omega_gan)^-1.
    period_ratio = (1.0 + sys_c.omega_gan) ** (-1.0)
    assert period_ratio == pytest.approx(2.101, abs=0.01)
    assert 2.09 < period_ratio < 2.12  # near, not exactly, 2:1


def test_a_gan_greater_than_one_titania_is_outer() -> None:
    """Titania (the perturber) really is farther from Uranus than Umbriel (the
    base moon) -- a_gan > 1, the model's own structural assumption."""
    sys_c = ut.uranus_umbriel_titania_default()
    assert sys_c.a_gan > 1.0


def test_mu_gan_comparable_order_to_jeg_within_documented_2x_gap() -> None:
    """`#693`'s headline claim for this pair: mu_gan is only ~2x below the
    already-validated JEG mu_gan (~7.8e-5), the best-conditioned non-Jovian
    candidate surveyed."""
    ut_sys = ut.uranus_umbriel_titania_default()
    jeg = ccr4bp.jupiter_europa_ganymede_default()
    ratio = jeg.mu_gan / ut_sys.mu_gan
    assert 1.5 < ratio < 2.5, ratio


# ---------------------------------------------------------------------------
# Structural reduction: mu_gan -> 0 reduces to the bare Uranus-Umbriel CR3BP,
# exactly as #689 proved for Jupiter-Europa (core.ccr4bp is UNMODIFIED --
# this just exercises it at a different mu).
# ---------------------------------------------------------------------------


def test_mugan0_reduces_to_uranus_umbriel_cr3bp() -> None:
    """At mu_gan=0 the CCR4BP EOM (already proven generic by #689) reduces to
    the plain Uranus-Umbriel CR3BP at THIS system's mu -- the structural check
    on the CONSTRUCTOR's mu value, reusing #689's already-validated reduction."""
    sys_c = ut.uranus_umbriel_titania_default()
    sys0 = ccr4bp.CCR4BPSystem(
        mu=sys_c.mu, mu_gan=0.0, a_gan=sys_c.a_gan, omega_gan=sys_c.omega_gan
    )
    rhs_c = ccr4bp.ccr4bp_eom(0.0, _SAMPLE_STATE, sys0)
    rhs_3 = cr3bp.cr3bp_eom(0.0, _SAMPLE_STATE, sys0.mu)
    assert np.allclose(rhs_c, rhs_3, rtol=0.0, atol=1e-14)


def test_mugan0_stm_matches_cr3bp() -> None:
    sys_c = ut.uranus_umbriel_titania_default()
    sys0 = ccr4bp.CCR4BPSystem(
        mu=sys_c.mu, mu_gan=0.0, a_gan=sys_c.a_gan, omega_gan=sys_c.omega_gan
    )
    y42 = np.concatenate([_SAMPLE_STATE, np.eye(6).reshape(36)])
    rhs_c = ccr4bp.ccr4bp_stm_eom(0.5, y42, sys0)
    rhs_3 = cr3bp.cr3bp_stm_eom(0.5, y42, sys0.mu)
    assert np.allclose(rhs_c, rhs_3, rtol=0.0, atol=1e-13)


def test_ganymede_role_synodic_position_and_period() -> None:
    sys_c = ut.uranus_umbriel_titania_default()
    for t in np.linspace(0.0, sys_c.ganymede_synodic_period, 11):
        gx, gy, gz = ccr4bp._ganymede_position(float(t), sys_c)
        assert gz == 0.0
        assert abs(float(np.hypot(gx, gy)) - sys_c.a_gan) < 1e-12
    p0 = ccr4bp._ganymede_position(0.0, sys_c)
    p1 = ccr4bp._ganymede_position(sys_c.ganymede_synodic_period, sys_c)
    assert np.allclose(p0, p1, atol=1e-10)


# ---------------------------------------------------------------------------
# L_KM / v_unit_km_s.
# ---------------------------------------------------------------------------


def test_v_unit_km_s_matches_manual_two_body_formula() -> None:
    """v_unit_km_s() == L_KM * sqrt((GM_Uranus + GM_Umbriel) / L_KM**3)
    (manual, independent recomputation -- guards against a copy-paste
    unit-formula bug like the one #695 documented in #694's
    Europa-hardcoded ``_v_unit_km_s``)."""
    gm_u = PRIMARIES["Uranus"]
    gm_umbriel = SATELLITES["Umbriel"].mu_km3_s2
    l_km = SATELLITES["Umbriel"].sma_km
    expected = l_km * math.sqrt((gm_u + gm_umbriel) / l_km**3)
    assert ut.v_unit_km_s() == pytest.approx(expected, rel=1e-15)
    assert pytest.approx(SATELLITES["Umbriel"].sma_km, rel=0.0, abs=0.0) == ut.L_KM


def test_l_km_differs_from_jeg_and_galilean_units() -> None:
    """A guardrail against silently reusing another system's length unit:
    Umbriel's SMA is neither Europa's (671,100 km) nor Io's (421,800 km)."""
    assert pytest.approx(265_986.0, abs=1.0) == ut.L_KM
    assert abs(ut.L_KM - 671_100.0) > 1000.0
    assert abs(ut.L_KM - 421_800.0) > 1000.0
