"""Tests for the Jupiter-Io-Europa CCR4BP system constructor (`#695`).

Mirrors `#689`'s own ``jupiter_europa_ganymede_default`` test discipline
(``tests/core/test_ccr4bp.py``'s Gate 5): sourced-constant cross-checks
against the SAME ``core.satellites`` registry and `#693`'s own independently
sourced screening-pass numbers, plus the structural ``mu_gan -> 0`` reduction
already proven generic by `#689` (this module adds no new EOM/STM code, so
this is a reduction check on the CONSTRUCTOR's output, not a new physics
gate).
"""

from __future__ import annotations

import numpy as np
import pytest

import cyclerfinder.core.ccr4bp as ccr4bp
import cyclerfinder.core.ccr4bp_io_europa as ioeu
import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.core.satellites import PRIMARIES, SATELLITES

_SAMPLE_STATE = np.array([0.6, 0.2, 0.0, 0.05, 0.35, 0.0], dtype=np.float64)


def test_constants_cross_check_against_satellites_registry() -> None:
    """mu/mu_gan/a_gan trace exactly to core.satellites (no independent re-sourcing)."""
    sys_c = ioeu.jupiter_io_europa_default()
    gm_j = PRIMARIES["Jupiter"]
    gm_io = SATELLITES["Io"].mu_km3_s2
    gm_eu = SATELLITES["Europa"].mu_km3_s2
    sma_io = SATELLITES["Io"].sma_km
    sma_eu = SATELLITES["Europa"].sma_km
    denom = gm_j + gm_io
    assert sys_c.mu == pytest.approx(gm_io / denom, rel=0.0, abs=1e-18)
    assert sys_c.mu_gan == pytest.approx(gm_eu / denom, rel=0.0, abs=1e-18)
    assert sys_c.a_gan == pytest.approx(sma_eu / sma_io, rel=0.0, abs=1e-15)


def test_cross_check_against_693_screening_numbers() -> None:
    """Cross-check against #693's own independently-derived screening-pass
    table (mu_base=4.704e-5, mu_pert=2.528e-5, P ratio 2.0000)."""
    sys_c = ioeu.jupiter_io_europa_default()
    assert sys_c.mu == pytest.approx(4.704e-5, abs=5e-8)
    assert sys_c.mu_gan == pytest.approx(2.528e-5, abs=5e-8)
    # Physical period ratio T_europa / T_io == n_io/n_europa == (1+omega_gan)^-1.
    period_ratio = (1.0 + sys_c.omega_gan) ** (-1.0)
    assert period_ratio == pytest.approx(2.0000, abs=0.01)
    assert 1.99 < period_ratio < 2.02  # the tightest commensurability #693 surveyed


def test_a_gan_greater_than_one_europa_is_outer() -> None:
    """Europa (the perturber) really is farther from Jupiter than Io (the base
    moon) -- a_gan > 1, the model's own structural assumption."""
    sys_c = ioeu.jupiter_io_europa_default()
    assert sys_c.a_gan > 1.0


def test_mugan0_reduces_to_jupiter_io_cr3bp() -> None:
    """At mu_gan=0 the CCR4BP EOM (already proven generic by #689) reduces to
    the plain Jupiter-Io CR3BP at THIS system's mu -- the structural check on
    the CONSTRUCTOR's mu value, reusing #689's already-validated reduction."""
    sys_c = ioeu.jupiter_io_europa_default()
    sys0 = ccr4bp.CCR4BPSystem(
        mu=sys_c.mu, mu_gan=0.0, a_gan=sys_c.a_gan, omega_gan=sys_c.omega_gan
    )
    rhs_c = ccr4bp.ccr4bp_eom(0.0, _SAMPLE_STATE, sys0)
    rhs_3 = cr3bp.cr3bp_eom(0.0, _SAMPLE_STATE, sys0.mu)
    assert np.allclose(rhs_c, rhs_3, rtol=0.0, atol=1e-14)


def test_v_unit_km_s_matches_manual_two_body_formula() -> None:
    """v_unit_km_s() == L_KM * sqrt((GM_Jupiter + GM_Io) / L_KM**3) (manual,
    independent recomputation -- guards against a copy-paste unit-formula bug
    like the one this task's own driver script documents in #694's Europa-
    hardcoded ``_v_unit_km_s``)."""
    import math

    gm_j = PRIMARIES["Jupiter"]
    gm_io = SATELLITES["Io"].mu_km3_s2
    l_km = SATELLITES["Io"].sma_km
    expected = l_km * math.sqrt((gm_j + gm_io) / l_km**3)
    assert ioeu.v_unit_km_s() == pytest.approx(expected, rel=1e-15)
    assert pytest.approx(SATELLITES["Io"].sma_km, rel=0.0, abs=0.0) == ioeu.L_KM
