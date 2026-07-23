"""CCR4BP torus positive control for Jupiter-Io-Ganymede (#696).

Applies `#690`'s pseudospectral quasi-periodic-torus corrector
(``search.variational_ccr4bp_torus``, reused UNMODIFIED) to the new
Io-Ganymede system (``core.ccr4bp_io_ganymede.jupiter_io_ganymede_default``,
#696) as a genuine discovery attempt, not another capability-proof exercise.
Mirrors `#690`'s own test discipline exactly: the ``mu_gan=0`` regression is a
machine-precision-flat structural check; the physical-mass torus is reported
at its achieved invariance residual.

Base resonant orbit: NOT the naive "moon-ratio-matching" choice
(spacecraft:Io = 1:4, i.e. ``_resonant_symmetric_orbit(mu, 1, 4)``). That
orbit's semimajor axis (``resonance_semimajor(1,4) = 4^(2/3) ~ 2.5198``) sits
almost exactly AT Ganymede's own orbital radius (``a_gan ~ 2.5377`` -- only a
0.018 Io-SMA gap, ~7,500 km, comparable to Ganymede's own point-mass
Hill-radius scale ``mu_gan^(1/3) * a_gan ~ 0.11``). Propagating this orbit
under Ganymede's PHYSICAL forcing (over one Ganymede-synodic period) throws a
genuine ``RuntimeError`` from :func:`cyclerfinder.core.ccr4bp.propagate_ccr4bp`
("Required step size is less than spacing between numbers") -- a real
near-collision singularity in the direct Ganymede-acceleration term, NOT a
numerical artifact: the spacecraft and Ganymede's moving point both start on
the +x axis at essentially the same radius, and a family scan (varying the
Newton seed radius at the SAME fixed period=8*pi) found no nearby symmetric
family member with different eccentricity that clears Ganymede's orbit at
this exact resonance (the semimajor axis is fixed by Kepler's third law at
this period regardless of shape). This is a genuine, reportable dynamical
wall SPECIFIC to Io-Ganymede's near-4:1 Laplace-chain period ratio (unlike
JEG's Europa 3:4 seed, whose ``a=1.211`` sits comfortably clear of
``a_gan=1.595``) -- see :func:`test_1to4_resonant_orbit_collides_with_ganymede`.

This module instead uses the OTHER natural "4:1"-labeled resonant orbit,
spacecraft:Io = 4:1 (``_resonant_symmetric_orbit(mu, 4, 1)``,
``a = (1/4)^(2/3) ~ 0.3968``, INTERIOR to Io, comfortably far from Ganymede's
orbit -- minimum spacecraft-Ganymede distance over one full Ganymede-synodic
torus period is ``~2.14`` Io-SMA units, two orders of magnitude larger than
the 1:4 orbit's near-collision gap). This orbit converges cleanly as a
symmetric periodic CR3BP orbit (residual ~3e-14) and its torus converges
CLEANER than JEG's own Europa 3:4 torus (physical-mass residual ~1.5e-6 here
vs JEG's own ~1.6e-4), because it is much less eccentric.
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy.integrate import solve_ivp

import cyclerfinder.core.ccr4bp as ccr4bp
import cyclerfinder.core.cr3bp as cr3bp
import cyclerfinder.search.variational_ccr4bp_torus as vt
from cyclerfinder.core.ccr4bp_io_ganymede import jupiter_io_ganymede_default
from cyclerfinder.genome.composed_moon_map import resonance_semimajor


def _resonant_symmetric_orbit(
    mu: float, p_sc: int, q_moon: int, *, max_iter: int = 60, tol: float = 1e-12, cap: float = 0.05
) -> tuple[np.ndarray, float, float]:
    """Identical test-only scaffolding to `#690`/`#691`'s own test modules
    (base-CR3BP resonant-orbit seed; no production code under test)."""
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
    return jupiter_io_ganymede_default()


@pytest.fixture(scope="module")
def orbit_4to1(system: ccr4bp.CCR4BPSystem) -> tuple[np.ndarray, float]:
    """The Jupiter-Io 4:1 (spacecraft:Io) symmetric resonant periodic orbit
    (base CR3BP, INTERIOR, a~0.397) -- the tractable "4:1" seed (see module
    docstring for why the naive 1:4 exterior choice is a dead end here)."""
    s0, period, res = _resonant_symmetric_orbit(system.mu, 4, 1)
    assert res < 1e-10, f"resonant orbit did not converge: perp residual {res:.2e}"
    return s0, period


# ---------------------------------------------------------------------------
# The 1:4 "moon-ratio-matching" orbit is a genuine dead end: documented here
# so the choice of 4:1 (not 1:4) is falsifiable, not asserted.
# ---------------------------------------------------------------------------


def test_1to4_resonant_orbit_collides_with_ganymede(system: ccr4bp.CCR4BPSystem) -> None:
    """The OTHER natural "4:1"-labeled orbit (spacecraft:Io = 1:4, exterior,
    a~2.52) converges fine as a bare CR3BP periodic orbit, but its semimajor
    axis sits within 0.018 Io-SMA units of Ganymede's own orbital radius --
    propagating it under Ganymede's PHYSICAL forcing over one Ganymede-synodic
    period must fail (a genuine near-collision singularity, not solver noise;
    RuntimeError from ``propagate_ccr4bp``, e.g. "step size ... below spacing
    between numbers")."""
    s0, _period, res = _resonant_symmetric_orbit(system.mu, 1, 4)
    assert res < 1e-10  # the bare periodic orbit itself IS a valid solution
    a_1to4 = resonance_semimajor(1, 4)
    assert abs(a_1to4 - system.a_gan) < 0.02  # near-coincident with Ganymede's radius
    with pytest.raises(RuntimeError):
        ccr4bp.propagate_ccr4bp(system, s0, 2.0 * np.pi / abs(system.omega_gan), with_stm=True)


# ---------------------------------------------------------------------------
# Positive control: mu_gan=0 -> theta1-flat torus (the tractable 4:1 orbit).
# ---------------------------------------------------------------------------


def test_mu_gan_zero_flat_regression(
    system: ccr4bp.CCR4BPSystem, orbit_4to1: tuple[np.ndarray, float]
) -> None:
    """At mu_gan=0 the corrector recovers the base 4:1 orbit as a theta1-FLAT
    torus -- the corrector-level analogue of `#689`'s ``mu_gan -> 0``
    structural reduction, on this new system."""
    s0, period = orbit_4to1
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
    system: ccr4bp.CCR4BPSystem, orbit_4to1: tuple[np.ndarray, float]
) -> None:
    """The Io 4:1 resonant orbit persists as a genuine thin 2D invariant torus
    under Ganymede's PHYSICAL-mass forcing. Live-observed: this torus converges
    CLEANER than JEG's own Europa 3:4 torus (residual ~1.5e-6 here vs JEG's
    ~1.6e-4), because the base orbit is much less eccentric."""
    s0, period = orbit_4to1
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
    assert _theta1_content(phys) > 1e-8
    expected = 2.0 * np.pi * (2.0 * np.pi / period) / abs(system.omega_gan)
    assert phys.rho_strob == pytest.approx(expected, rel=1e-3)


def test_continue_ccr4bp_torus_mass_reaches_physical(
    system: ccr4bp.CCR4BPSystem, orbit_4to1: tuple[np.ndarray, float]
) -> None:
    """Mass homotopy from mu_gan=0 walks to the physical mass without diverging."""
    s0, period = orbit_4to1
    sys0 = ccr4bp.CCR4BPSystem(
        mu=system.mu, mu_gan=0.0, a_gan=system.a_gan, omega_gan=system.omega_gan
    )
    seed = vt.discover_ccr4bp_torus_from_resonant_orbit(
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
    steps = vt.continue_ccr4bp_torus_mass(
        seed,
        system,
        n_steps=3,
        tr_solver="exact",
        max_nfev=600,
        gauge_weight=30.0,
        rho_weight=100.0,
    )
    assert len(steps) == 3
    final = steps[-1]
    assert final.system.mu_gan == pytest.approx(system.mu_gan)
    assert final.residual_rms < 5e-3, final.residual_rms
    assert _theta1_content(final) > _theta1_content(seed)
