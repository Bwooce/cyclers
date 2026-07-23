"""CCR4BP torus positive control for Jupiter-Io-Europa (`#695`).

Applies `#690`'s pseudospectral corrector (``search.variational_ccr4bp_torus``,
reused UNMODIFIED) to `#695`'s new Io-Europa system
(``core.ccr4bp_io_europa.jupiter_io_europa_default``). Mirrors `#690`'s own
test discipline (``tests/search/test_variational_ccr4bp_torus.py``) exactly:
the ``mu_gan=0`` machine-precision-flat regression FIRST, then the
physical-mass thin-torus check, before any manifold/whisker work trusts this
torus.

Base resonant orbit choice
--------------------------
`#690`/`#691`'s own JEG tests seed from the Jupiter-Europa 3:4 SYMMETRIC
resonant orbit (spacecraft:moon = 3:4, i.e. ``T_sc/T_Europa = 4/3``). For
Io-Europa, this task picks the spacecraft:Io = 1:2 resonance (``T_sc/T_Io =
2``, i.e. the spacecraft's own period equals TWICE Io's period) -- the
task's own dispatch describes this as "the symmetric 2:1 resonant CR3BP
periodic orbit" (period ratio T_sc:T_Io = 2:1). Two candidate resonances were
tried (``resonance_semimajor``-seeded Newton continuation, identical
scaffolding to `#690`'s own test module): the INTERIOR reading
(spacecraft:Io = 2:1, a0~0.63) does NOT converge with this damped-Newton
scheme (residual floor ~3.5e-2 after 79 iterations -- a genuine, reported
non-convergence, not swept under the rug); the EXTERIOR reading
(spacecraft:Io = 1:2, a0~1.70) converges cleanly to ~1.9e-14. This EXTERIOR
orbit is used throughout. It has a physically elegant property germane to
this specific system: since Io:Europa are themselves ~2:1 commensurate, this
orbit's own two-body semimajor-axis GUESS (1.587, in Io-SMA units) lands
almost exactly on Europa's actual position (``a_gan=1.591``) -- a genuine
consequence of the resonance chosen matching the moon pair's own physical
commensurability, not a coincidence engineered into the test. The
Newton-corrected orbit's true IC (a0=1.698 after correction) ranges
1.477-1.698 in Io-SMA units (a modest eccentricity ~0.07, apoapsis
extending slightly beyond Europa's own SMA) -- directly analogous to JEG's
own 3:4 orbit's IC (a0=1.610) extending slightly beyond Ganymede's SMA
(1.595). Periapsis (0.698 Io-SMA units from Io itself) stays far outside
Io's Hill sphere, so no close-approach singularity risk.
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy.integrate import solve_ivp

import cyclerfinder.core.ccr4bp as ccr4bp
import cyclerfinder.core.ccr4bp_io_europa as ioeu
import cyclerfinder.core.cr3bp as cr3bp
import cyclerfinder.search.variational_ccr4bp_torus as vt
from cyclerfinder.genome.composed_moon_map import resonance_semimajor


def _resonant_symmetric_orbit(
    mu: float, p_sc: int, q_moon: int, *, max_iter: int = 80, tol: float = 1e-12, cap: float = 0.05
) -> tuple[np.ndarray, float, float]:
    """Symmetric p:q resonant CR3BP periodic orbit. Identical test-only
    scaffolding to #690/#691/#694's own test/driver-script copies (no
    production code under test)."""
    a = resonance_semimajor(p_sc, q_moon)
    period = 2.0 * np.pi * q_moon
    th = 0.5 * period
    x0 = a - mu
    vy0 = float(np.sqrt((1.0 - mu) / a)) - x0
    res = np.inf
    for k in range(max_iter):
        s0 = np.array([x0, 0.0, 0.0, 0.0, vy0, 0.0])
        y42 = np.concatenate([s0, np.eye(6).reshape(36)])
        sol = solve_ivp(
            cr3bp.cr3bp_stm_eom, (0.0, th), y42, args=(mu,), method="DOP853", rtol=1e-12, atol=1e-12
        )
        sf = sol.y[:, -1]
        phi = sf[6:].reshape(6, 6)
        g = np.array([sf[1], sf[3]])
        res = float(np.linalg.norm(g))
        if res < tol:
            break
        jac = np.array([[phi[1, 0], phi[1, 4]], [phi[3, 0], phi[3, 4]]])
        dz = np.linalg.solve(jac, -g)
        dz = dz * (0.3 if k < 8 else 1.0)
        norm = float(np.linalg.norm(dz))
        if norm > cap:
            dz = dz / norm * cap
        x0 += dz[0]
        vy0 += dz[1]
    return np.array([x0, 0.0, 0.0, 0.0, vy0, 0.0]), period, res


def _theta1_content(result: vt.CCR4BPTorusVariationalResult) -> float:
    return float(np.sqrt(np.sum(result.coeffs[:, 1:, :] ** 2)))


@pytest.fixture(scope="module")
def system() -> ccr4bp.CCR4BPSystem:
    return ioeu.jupiter_io_europa_default()


@pytest.fixture(scope="module")
def orbit_12(system: ccr4bp.CCR4BPSystem) -> tuple[np.ndarray, float]:
    """The Jupiter-Io spacecraft:Io = 1:2 symmetric resonant periodic orbit
    (base CR3BP) -- see module docstring for why this reading (not the
    interior 2:1 reading, which fails to converge with this scaffolding)."""
    s0, period, res = _resonant_symmetric_orbit(system.mu, 1, 2)
    assert res < 1e-10, f"resonant orbit did not converge: perp residual {res:.2e}"
    return s0, period


def test_interior_2to1_reading_does_not_converge(system: ccr4bp.CCR4BPSystem) -> None:
    """Documents the negative: the OTHER plausible reading of "2:1"
    (spacecraft:Io = 2:1, interior, a0~0.63) does not converge with this
    damped-Newton scaffolding -- an honest reported non-convergence, the
    reason the exterior (1:2) reading is used throughout this task."""
    _s0, _period, res = _resonant_symmetric_orbit(system.mu, 2, 1)
    assert res > 1e-6, f"expected non-convergence, got residual {res:.2e}"


def test_base_orbit_stays_clear_of_io_and_near_europa_sma(
    system: ccr4bp.CCR4BPSystem, orbit_12: tuple[np.ndarray, float]
) -> None:
    """Sanity: the seed orbit's periapsis is far from Io (no close-approach
    singularity) and its extent brackets Europa's own SMA (the elegant
    consequence of picking a spacecraft:Io=1:2 resonance in a system where
    Io:Europa are themselves ~2:1 commensurate)."""
    s0, period = orbit_12
    mu = system.mu
    sol = solve_ivp(
        cr3bp.cr3bp_eom,
        (0.0, period),
        s0,
        args=(mu,),
        method="DOP853",
        rtol=1e-13,
        atol=1e-13,
        dense_output=True,
    )
    ts = np.linspace(0.0, period, 2000)
    xs = np.array([sol.sol(t) for t in ts])
    r_io = np.sqrt((xs[:, 0] - 1.0 + mu) ** 2 + xs[:, 1] ** 2)
    r_jup = np.sqrt((xs[:, 0] + mu) ** 2 + xs[:, 1] ** 2)
    assert r_io.min() > 0.3  # far outside Io's Hill sphere (~0.036)
    assert r_jup.min() < system.a_gan < r_jup.max() * 1.1  # brackets Europa's SMA


def test_mu_gan_zero_flat_regression(
    system: ccr4bp.CCR4BPSystem, orbit_12: tuple[np.ndarray, float]
) -> None:
    """At mu_gan=0 the corrector recovers the base 1:2 orbit as a theta1-FLAT
    torus -- the corrector-level analogue of #689's mu_gan->0 structural
    reduction, on THIS new system (mirrors #690's own test exactly)."""
    s0, period = orbit_12
    sys0 = ccr4bp.CCR4BPSystem(
        mu=system.mu, mu_gan=0.0, a_gan=system.a_gan, omega_gan=system.omega_gan
    )
    res = vt.discover_ccr4bp_torus_from_resonant_orbit(
        sys0,
        s0,
        period,
        n1=1,
        n2=20,
        tr_solver="exact",
        max_nfev=600,
        gauge_weight=30.0,
        rho_weight=100.0,
    )
    assert _theta1_content(res) < 1e-8, _theta1_content(res)
    assert res.residual_rms < 1e-3, res.residual_rms
    expected_rot = (2.0 * np.pi / period) / abs(system.omega_gan)
    assert res.rotation_number == pytest.approx(expected_rot, rel=1e-3)


def test_physical_mass_thin_torus(
    system: ccr4bp.CCR4BPSystem, orbit_12: tuple[np.ndarray, float]
) -> None:
    """The 1:2 Io resonant orbit persists as a genuine thin 2D invariant torus
    under Europa's PHYSICAL-mass forcing -- the actual discovery-attempt
    object this task builds on. Reports invariance residual, an INDEPENDENT
    closure residual, the stroboscopic rotation number, and theta1-forcing
    content appearing only when mu_gan > 0 (mirrors #690's own
    test_physical_mass_thin_torus exactly, on the new system)."""
    s0, period = orbit_12
    sys0 = ccr4bp.CCR4BPSystem(
        mu=system.mu, mu_gan=0.0, a_gan=system.a_gan, omega_gan=system.omega_gan
    )
    flat = vt.discover_ccr4bp_torus_from_resonant_orbit(
        sys0,
        s0,
        period,
        n1=1,
        n2=20,
        tr_solver="exact",
        max_nfev=600,
        gauge_weight=30.0,
        rho_weight=100.0,
    )
    phys = vt.discover_ccr4bp_torus_from_resonant_orbit(
        system,
        s0,
        period,
        n1=1,
        n2=20,
        tr_solver="exact",
        max_nfev=600,
        gauge_weight=30.0,
        rho_weight=100.0,
    )
    assert phys.system.mu_gan > 0.0
    assert phys.residual_rms < 1e-3, phys.residual_rms
    assert phys.closure_residual < 5e-3, phys.closure_residual
    assert _theta1_content(phys) > 100.0 * _theta1_content(flat)
    assert _theta1_content(phys) > 1e-3
    expected = 2.0 * np.pi * (2.0 * np.pi / period) / abs(system.omega_gan)
    assert phys.rho_strob == pytest.approx(expected, rel=1e-3)
