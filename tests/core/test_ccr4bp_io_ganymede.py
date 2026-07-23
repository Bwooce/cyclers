"""Tests for the Jupiter-Io-Ganymede CCR4BP system constructor (#696).

Mirrors `#689`'s own ``tests/core/test_ccr4bp.py`` sourcing/reduction
discipline, applied to the new ``jupiter_io_ganymede_default()`` constructor
(``core.ccr4bp_io_ganymede``) instead of building anything new in
``core.ccr4bp`` itself (that module is reused UNMODIFIED, per `#696`'s scope).
No value computed by our own code sits on the EXPECTED side of a structural
assertion -- sourced constants (``core.satellites`` / `#693`'s own screening
table) or the already-proven ``core.ccr4bp``/``core.cr3bp`` reduction
structure only.
"""

from __future__ import annotations

import numpy as np
import pytest

import cyclerfinder.core.ccr4bp as ccr4bp
import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.core.ccr4bp_io_ganymede import jupiter_io_ganymede_default
from cyclerfinder.genome.composed_moon_map import moon_config

_SAMPLE_STATE_PLANAR = np.array([0.6, 0.2, 0.0, 0.05, 0.35, 0.0], dtype=np.float64)


# ---------------------------------------------------------------------------
# Sourced-constant cross-checks.
# ---------------------------------------------------------------------------


def test_constants_cross_check_against_composed_moon_map() -> None:
    """The Io-Ganymede default parameters trace to the SAME JPL SSD registry
    values `#688`'s composed_moon_map uses -- radii, mass ratios, period ratio."""
    sys_c = jupiter_io_ganymede_default()
    io = moon_config("Io")
    ganymede = moon_config("Ganymede")

    assert sys_c.mu == io.mu
    assert sys_c.a_gan == pytest.approx(ganymede.sma_km / io.sma_km, rel=0.0, abs=1e-15)
    assert sys_c.mu3_reduction == pytest.approx(ganymede.mu, rel=1e-4)

    ratio_composed = io.mean_motion_rad_s / ganymede.mean_motion_rad_s
    period_ratio_ccr4bp = (1.0 + sys_c.omega_gan) ** (-1.0)  # n_Io/n_Ganymede = T_Ganymede/T_Io
    assert ratio_composed == pytest.approx(4.04, abs=0.05)
    assert period_ratio_ccr4bp == pytest.approx(ratio_composed, abs=1e-3)
    assert 3.9 < period_ratio_ccr4bp < 4.1  # near, not exactly, the 4:1 Laplace-chain ratio


def test_constants_cross_check_against_693_screening_table() -> None:
    """Cross-check against `#693`'s own sourced numeric table (docs/notes/
    2026-07-23-693-ccr4bp-moonpair-screening.md): mu_base(Io)=4.704e-5,
    mu_pert(Ganymede)=7.805e-5, period ratio 4.059."""
    sys_c = jupiter_io_ganymede_default()
    assert sys_c.mu == pytest.approx(4.704e-5, rel=1e-3)
    assert sys_c.mu_gan == pytest.approx(7.805e-5, rel=1e-3)
    period_ratio = (1.0 + sys_c.omega_gan) ** (-1.0)
    # #693's own value (4.059) is JPL-fitted; ours is two-body-Kepler-on-SMA
    # (same ~0.4% class of gap #689's own docstring documents for JEG's
    # 2.014-vs-2.030 Europa:Ganymede ratio) -- agree to 1%, not to the digit.
    assert period_ratio == pytest.approx(4.059, rel=1e-2)


def test_mu_gan_essentially_identical_to_jeg() -> None:
    """`#693`'s headline claim: Io-Ganymede's mu_gan is essentially IDENTICAL
    forcing strength to the already-validated JEG mu_gan (same Ganymede mass,
    both denominators dominated by GM_Jupiter)."""
    io_gan = jupiter_io_ganymede_default()
    jeg = ccr4bp.jupiter_europa_ganymede_default()
    assert io_gan.mu_gan == pytest.approx(jeg.mu_gan, rel=2e-2)


# ---------------------------------------------------------------------------
# Structural reduction: mu_gan -> 0 reduces to the bare Jupiter-Io CR3BP,
# exactly as #689 proved for Jupiter-Europa (core.ccr4bp is UNMODIFIED --
# this just exercises it at a different mu).
# ---------------------------------------------------------------------------


def test_mugan0_limit_eom_matches_cr3bp() -> None:
    d = jupiter_io_ganymede_default()
    sys0 = ccr4bp.CCR4BPSystem(mu=d.mu, mu_gan=0.0, a_gan=d.a_gan, omega_gan=d.omega_gan)
    rhs_c = ccr4bp.ccr4bp_eom(0.0, _SAMPLE_STATE_PLANAR, sys0)
    rhs_3 = cr3bp.cr3bp_eom(0.0, _SAMPLE_STATE_PLANAR, sys0.mu)
    assert np.allclose(rhs_c, rhs_3, rtol=0.0, atol=1e-14)


def test_mugan0_limit_stm_matches_cr3bp() -> None:
    d = jupiter_io_ganymede_default()
    sys0 = ccr4bp.CCR4BPSystem(mu=d.mu, mu_gan=0.0, a_gan=d.a_gan, omega_gan=d.omega_gan)
    y42 = np.concatenate([_SAMPLE_STATE_PLANAR, np.eye(6).reshape(36)])
    rhs_c = ccr4bp.ccr4bp_stm_eom(0.5, y42, sys0)
    rhs_3 = cr3bp.cr3bp_stm_eom(0.5, y42, sys0.mu)
    assert np.allclose(rhs_c, rhs_3, rtol=0.0, atol=1e-13)


def test_ganymede_synodic_position_and_period() -> None:
    sys_c = jupiter_io_ganymede_default()
    for t in np.linspace(0.0, sys_c.ganymede_synodic_period, 11):
        gx, gy, gz = ccr4bp._ganymede_position(float(t), sys_c)
        assert gz == 0.0
        assert abs(float(np.hypot(gx, gy)) - sys_c.a_gan) < 1e-12
    p0 = ccr4bp._ganymede_position(0.0, sys_c)
    p1 = ccr4bp._ganymede_position(sys_c.ganymede_synodic_period, sys_c)
    assert np.allclose(p0, p1, atol=1e-10)
