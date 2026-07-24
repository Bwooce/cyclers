"""Tests for the Jupiter-Europa-Callisto CCR4BP system constructor (`#703`).

Mirrors `#695`'s own ``test_ccr4bp_io_europa.py`` / `#696`'s
``test_ccr4bp_io_ganymede.py`` / `#701`'s ``test_ccr4bp_umbriel_titania.py``
sourcing/reduction discipline, applied to the new
``jupiter_europa_callisto_default()`` constructor
(``core.ccr4bp_europa_callisto``) instead of building anything new in
``core.ccr4bp`` itself (that module is reused UNMODIFIED, per `#703`'s scope).
No value computed by our own code sits on the EXPECTED side of a structural
assertion -- sourced constants (``core.satellites`` / `#693`'s and `#700`'s
own screening tables) or the already-proven ``core.ccr4bp``/``core.cr3bp``
reduction structure only.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

import cyclerfinder.core.ccr4bp as ccr4bp
import cyclerfinder.core.ccr4bp_europa_callisto as ec
import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.core.satellites import PRIMARIES, SATELLITES

_SAMPLE_STATE = np.array([0.6, 0.2, 0.0, 0.05, 0.35, 0.0], dtype=np.float64)


# ---------------------------------------------------------------------------
# Sourced-constant cross-checks.
# ---------------------------------------------------------------------------


def test_constants_cross_check_against_satellites_registry() -> None:
    """mu/mu_gan/a_gan trace exactly to core.satellites (no independent re-sourcing)."""
    sys_c = ec.jupiter_europa_callisto_default()
    gm_j = PRIMARIES["Jupiter"]
    gm_europa = SATELLITES["Europa"].mu_km3_s2
    gm_callisto = SATELLITES["Callisto"].mu_km3_s2
    sma_europa = SATELLITES["Europa"].sma_km
    sma_callisto = SATELLITES["Callisto"].sma_km
    denom = gm_j + gm_europa
    assert sys_c.mu == pytest.approx(gm_europa / denom, rel=0.0, abs=1e-18)
    assert sys_c.mu_gan == pytest.approx(gm_callisto / denom, rel=0.0, abs=1e-18)
    assert sys_c.a_gan == pytest.approx(sma_callisto / sma_europa, rel=0.0, abs=1e-15)


def test_mu_matches_jeg_own_mu_exactly() -> None:
    """Europa's base role is UNCHANGED from JEG: this system's mu must equal
    the already-validated jupiter_europa_ganymede_default()'s own mu exactly
    (same numerator, same denominator -- both keyed off Europa/Jupiter)."""
    sys_c = ec.jupiter_europa_callisto_default()
    jeg = ccr4bp.jupiter_europa_ganymede_default()
    assert sys_c.mu == jeg.mu


def test_cross_check_against_693_700_screening_numbers() -> None:
    """Cross-check against #693's/#700's own independently-derived screening-pass
    table (mu_base=2.528e-5, mu_pert=5.667e-5, a_pert=2.8054, P ratio 4.734)."""
    sys_c = ec.jupiter_europa_callisto_default()
    assert sys_c.mu == pytest.approx(2.528e-5, abs=5e-8)
    assert sys_c.mu_gan == pytest.approx(5.667e-5, abs=5e-8)
    assert sys_c.a_gan == pytest.approx(2.8054, abs=1e-3)
    # Physical period ratio T_callisto / T_europa == n_europa/n_callisto == (1+omega_gan)^-1.
    period_ratio = (1.0 + sys_c.omega_gan) ** (-1.0)
    assert period_ratio == pytest.approx(4.734, abs=0.05)
    assert 4.6 < period_ratio < 4.9  # loose, NOT a clean low-integer commensurability


def test_a_gan_greater_than_one_callisto_is_outer() -> None:
    """Callisto (the perturber) really is farther from Jupiter than Europa
    (the base moon) -- a_gan > 1, the model's own structural assumption."""
    sys_c = ec.jupiter_europa_callisto_default()
    assert sys_c.a_gan > 1.0


def test_mu_gan_comparable_order_to_jeg() -> None:
    """`#693`'s/`#700`'s headline claim for this pair: mu_gan is comparable
    order-of-magnitude forcing strength to the already-validated JEG mu_gan
    (~7.8e-5) -- Callisto's own mass is comparable to Ganymede's."""
    ec_sys = ec.jupiter_europa_callisto_default()
    jeg = ccr4bp.jupiter_europa_ganymede_default()
    ratio = jeg.mu_gan / ec_sys.mu_gan
    assert 1.0 < ratio < 2.0, ratio


# ---------------------------------------------------------------------------
# Structural reduction: mu_gan -> 0 reduces to the bare Jupiter-Europa CR3BP,
# exactly as #689 proved (core.ccr4bp is UNMODIFIED -- this just exercises it
# at Europa's own mu, already validated by #689's own test module too).
# ---------------------------------------------------------------------------


def test_mugan0_reduces_to_jupiter_europa_cr3bp() -> None:
    sys_c = ec.jupiter_europa_callisto_default()
    sys0 = ccr4bp.CCR4BPSystem(
        mu=sys_c.mu, mu_gan=0.0, a_gan=sys_c.a_gan, omega_gan=sys_c.omega_gan
    )
    rhs_c = ccr4bp.ccr4bp_eom(0.0, _SAMPLE_STATE, sys0)
    rhs_3 = cr3bp.cr3bp_eom(0.0, _SAMPLE_STATE, sys0.mu)
    assert np.allclose(rhs_c, rhs_3, rtol=0.0, atol=1e-14)


def test_mugan0_stm_matches_cr3bp() -> None:
    sys_c = ec.jupiter_europa_callisto_default()
    sys0 = ccr4bp.CCR4BPSystem(
        mu=sys_c.mu, mu_gan=0.0, a_gan=sys_c.a_gan, omega_gan=sys_c.omega_gan
    )
    y42 = np.concatenate([_SAMPLE_STATE, np.eye(6).reshape(36)])
    rhs_c = ccr4bp.ccr4bp_stm_eom(0.5, y42, sys0)
    rhs_3 = cr3bp.cr3bp_stm_eom(0.5, y42, sys0.mu)
    assert np.allclose(rhs_c, rhs_3, rtol=0.0, atol=1e-13)


def test_callisto_role_synodic_position_and_period() -> None:
    sys_c = ec.jupiter_europa_callisto_default()
    for t in np.linspace(0.0, sys_c.ganymede_synodic_period, 11):
        gx, gy, gz = ccr4bp._ganymede_position(float(t), sys_c)
        assert gz == 0.0
        assert abs(float(np.hypot(gx, gy)) - sys_c.a_gan) < 1e-12
    p0 = ccr4bp._ganymede_position(0.0, sys_c)
    p1 = ccr4bp._ganymede_position(sys_c.ganymede_synodic_period, sys_c)
    assert np.allclose(p0, p1, atol=1e-10)


# ---------------------------------------------------------------------------
# L_KM / v_unit_km_s -- the "coincidence, verified not assumed" claim.
# ---------------------------------------------------------------------------


def test_l_km_equals_europa_sma_and_module_native_jeg_constant() -> None:
    """Unlike `#695`/`#696`/`#701`'s own non-Jovian-or-non-Europa-base
    systems, THIS system's base moon IS Europa -- so L_KM must equal
    671,100 km exactly (Europa's own SMA), the SAME value
    ``search.ccr4bp_heteroclinic_search`` hardcodes as ``_L_KM``. Verified
    directly against the registry, not asserted by construction alone."""
    assert pytest.approx(671_100.0, abs=1.0) == ec.L_KM
    assert SATELLITES["Europa"].sma_km == ec.L_KM


def test_v_unit_km_s_matches_manual_two_body_formula_and_module_native() -> None:
    """v_unit_km_s() == L_KM * sqrt((GM_Jupiter + GM_Europa) / L_KM**3)
    (manual, independent recomputation -- guards against a copy-paste
    unit-formula bug like the one #695 documented in #694's own
    ``_v_unit_km_s``), AND agrees with that module's own hardcoded
    (Europa-based) velocity unit for THIS system specifically."""
    gm_j = PRIMARIES["Jupiter"]
    gm_europa = SATELLITES["Europa"].mu_km3_s2
    l_km = SATELLITES["Europa"].sma_km
    expected = l_km * math.sqrt((gm_j + gm_europa) / l_km**3)
    assert ec.v_unit_km_s() == pytest.approx(expected, rel=1e-15)

    # Module-native _v_unit_km_s (hardcoded to Europa/Jupiter) computed
    # independently here (not imported -- that helper is private) to confirm
    # the two formulae agree for THIS system, not merely by shared code path.
    n_europa_native = math.sqrt((gm_j + gm_europa) / 671_100.0**3)
    v_unit_native = 671_100.0 * n_europa_native
    assert ec.v_unit_km_s() == pytest.approx(v_unit_native, rel=1e-12)
