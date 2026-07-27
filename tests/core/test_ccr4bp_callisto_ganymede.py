"""Tests for the Jupiter-Callisto-Ganymede CCR4BP system constructor (`#715`).

Mirrors `#696`'s ``tests/core/test_ccr4bp_io_ganymede.py`` sourcing/reduction
discipline, applied to the new ``jupiter_callisto_ganymede_default()``
constructor (``core.ccr4bp_callisto_ganymede``). No value computed by our
own code sits on the EXPECTED side of a structural assertion -- sourced
constants (``core.satellites`` / the `#710` digest's own quoted Aryan &
Fitzgerald Table 2 numbers) or the already-proven ``core.ccr4bp``/
``core.cr3bp`` reduction structure only.

The genuinely NEW thing this file checks (not present in any prior sibling
test) is the ``a_gan < 1`` / ``omega_gan > 0`` INTERIOR-perturber case --
every other CCR4BP system built so far has an OUTER perturber. See the
module docstring of ``core.ccr4bp_callisto_ganymede`` for why this is
expected to still work (the EOM/synodic-rate formulae are algebraically
indifferent to the perturber's radial placement).
"""

from __future__ import annotations

import numpy as np
import pytest

import cyclerfinder.core.ccr4bp as ccr4bp
import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.core.ccr4bp_callisto_ganymede import (
    jupiter_callisto_ganymede_default,
    v_unit_km_s,
)
from cyclerfinder.genome.composed_moon_map import moon_config

_SAMPLE_STATE_PLANAR = np.array([0.6, 0.2, 0.0, 0.05, 0.35, 0.0], dtype=np.float64)


# ---------------------------------------------------------------------------
# Sourced-constant cross-checks.
# ---------------------------------------------------------------------------


def test_constants_cross_check_against_composed_moon_map() -> None:
    """The Callisto-Ganymede default parameters trace to the SAME JPL SSD
    registry values `#688`'s composed_moon_map uses -- radii, mass ratios."""
    sys_c = jupiter_callisto_ganymede_default()
    callisto = moon_config("Callisto")
    ganymede = moon_config("Ganymede")

    assert sys_c.mu == callisto.mu
    assert sys_c.a_gan == pytest.approx(ganymede.sma_km / callisto.sma_km, rel=0.0, abs=1e-15)
    assert sys_c.mu3_reduction == pytest.approx(ganymede.mu, rel=1e-4)

    # Ganymede is INTERIOR to Callisto -- unlike every prior CCR4BP system
    # this project has built, a_gan < 1 here.
    assert sys_c.a_gan < 1.0
    assert sys_c.omega_gan > 0.0  # faster inner perturber ADVANCES in the base frame


def test_constants_cross_check_against_693_screening_table() -> None:
    """Cross-check against `#693`'s own sourced numeric table (docs/notes/
    2026-07-23-693-ccr4bp-moonpair-screening.md), Ganymede->Callisto row read
    in the OTHER direction (mu_pert there is Callisto's ratio relative to
    Jupiter+Ganymede, not directly this module's mu/mu_gan convention -- so
    this check instead targets the raw Callisto/Ganymede GM ratios, which are
    convention-independent)."""
    sys_c = jupiter_callisto_ganymede_default()
    # #693's table: mu_base(Ganymede)=7.804e-5, mu_pert(Callisto)=5.667e-5,
    # both relative to Jupiter+Ganymede -- NOT the same denominator this
    # module uses (Jupiter+Callisto). Only the period ratio (denominator-
    # independent) transfers directly:
    period_ratio_693 = 2.333  # Callisto:Ganymede period ratio, #693 Table
    period_ratio_ccr4bp = 1.0 + sys_c.omega_gan  # n_ganymede/n_callisto = T_callisto/T_ganymede
    assert period_ratio_ccr4bp == pytest.approx(period_ratio_693, rel=1e-2)


def test_constants_cross_check_against_710_paper_table2() -> None:
    """Cross-check against the `#710` digest's own quoted Aryan & Fitzgerald
    Table 2 numbers (``docs/notes/2026-07-26-710-digest-aryan-fitzgerald-2024-
    jovian-pccfbp.md``): ``mu1=5.6623e-05`` (Callisto), ``mu2=7.7890e-05``
    (Ganymede) -- both at ``Theta3,0=0``. This project's independently-
    sourced JPL DE440 registry reproduces both to well under 0.1%/0.3%
    relative (live-observed: mu1 gap ~0.008%, mu2 gap ~0.2%), the same class
    of small cross-source residual `#689`'s own docstring documents for the
    Europa-Ganymede pair against Kumar 2021."""
    sys_c = jupiter_callisto_ganymede_default()
    assert sys_c.mu == pytest.approx(5.6623e-05, rel=2e-3)
    assert sys_c.mu_gan == pytest.approx(7.7890e-05, rel=3e-3)


# ---------------------------------------------------------------------------
# Structural reduction: mu_gan -> 0 reduces to the bare Jupiter-Callisto
# CR3BP, exactly as #689 proved for Jupiter-Europa (core.ccr4bp is
# UNMODIFIED -- this just exercises it at a different mu AND, for the first
# time, an a_gan < 1 interior perturber).
# ---------------------------------------------------------------------------


def test_mugan0_limit_eom_matches_cr3bp() -> None:
    d = jupiter_callisto_ganymede_default()
    sys0 = ccr4bp.CCR4BPSystem(mu=d.mu, mu_gan=0.0, a_gan=d.a_gan, omega_gan=d.omega_gan)
    rhs_c = ccr4bp.ccr4bp_eom(0.0, _SAMPLE_STATE_PLANAR, sys0)
    rhs_3 = cr3bp.cr3bp_eom(0.0, _SAMPLE_STATE_PLANAR, sys0.mu)
    assert np.allclose(rhs_c, rhs_3, rtol=0.0, atol=1e-14)


def test_mugan0_limit_stm_matches_cr3bp() -> None:
    d = jupiter_callisto_ganymede_default()
    sys0 = ccr4bp.CCR4BPSystem(mu=d.mu, mu_gan=0.0, a_gan=d.a_gan, omega_gan=d.omega_gan)
    y42 = np.concatenate([_SAMPLE_STATE_PLANAR, np.eye(6).reshape(36)])
    rhs_c = ccr4bp.ccr4bp_stm_eom(0.5, y42, sys0)
    rhs_3 = cr3bp.cr3bp_stm_eom(0.5, y42, sys0.mu)
    assert np.allclose(rhs_c, rhs_3, rtol=0.0, atol=1e-13)


def test_ganymede_synodic_position_and_period_interior() -> None:
    """The `#696` Io-Ganymede sibling test, applied here -- confirms the
    ``_ganymede_position`` circle-parametrisation and periodicity structure
    hold unmodified for an a_gan < 1 (interior) perturber."""
    sys_c = jupiter_callisto_ganymede_default()
    for t in np.linspace(0.0, sys_c.ganymede_synodic_period, 11):
        gx, gy, gz = ccr4bp._ganymede_position(float(t), sys_c)
        assert gz == 0.0
        assert abs(float(np.hypot(gx, gy)) - sys_c.a_gan) < 1e-12
    p0 = ccr4bp._ganymede_position(0.0, sys_c)
    p1 = ccr4bp._ganymede_position(sys_c.ganymede_synodic_period, sys_c)
    assert np.allclose(p0, p1, atol=1e-10)


def test_v_unit_km_s_positive_and_consistent_with_l_km() -> None:
    """Sanity check on the independently-recomputed velocity unit (mirrors
    `#703`'s own cross-check discipline for its L_KM/v_unit_km_s pair)."""
    v = v_unit_km_s()
    assert v > 0.0
    # A Jupiter-Callisto circular orbital speed is a few km/s (Callisto's
    # real orbital speed ~8.2 km/s) -- sanity range, not an exact target.
    assert 5.0 < v < 12.0
