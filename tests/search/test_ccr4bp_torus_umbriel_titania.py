"""CCR4BP torus positive control for Uranus-Umbriel-Titania (`#701`).

Applies `#690`'s pseudospectral quasi-periodic-torus corrector
(``search.variational_ccr4bp_torus``, reused UNMODIFIED) to `#701`'s new
Umbriel-Titania system (``core.ccr4bp_umbriel_titania.uranus_umbriel_titania_default``)
as a genuine discovery attempt, not another capability-proof exercise. Mirrors
`#690`'s/`#695`'s/`#696`'s own test discipline exactly: the ``mu_gan=0``
regression is a machine-precision-flat structural check; the physical-mass
torus is reported at its achieved invariance residual.

Base resonant orbit choice
---------------------------
Umbriel:Titania's real period ratio is ~2.101 -- near, but not exactly, 2:1
(the SAME "near-2:1, not exact" character as `#695`'s own Io:Europa pair,
2.000 there vs 2.101 here). This module therefore follows `#695`'s own
precedent (not JEG's own literature-driven "3:4" choice, which has no
special relevance to a non-published system): two candidate readings of "the
symmetric 2:1 resonant orbit" were tried, identical
``_resonant_symmetric_orbit`` Newton-continuation scaffolding to
`#690`/`#691`/`#694`/`#695`/`#696`'s own test/driver-script copies.

* The INTERIOR reading (spacecraft:Umbriel = 2:1, ``a0 ~ 0.63``) does NOT
  converge with this damped-Newton scheme (residual floor ~2.3e-2 after 80
  iterations) -- the SAME non-convergent reading `#695` found for Io-Europa's
  own interior 2:1 family. A genuine, reported non-convergence, not swept
  under the rug (:func:`test_interior_2to1_reading_does_not_converge`).
* The EXTERIOR reading (spacecraft:Umbriel = 1:2, ``a0 ~ 1.70``) converges
  cleanly (perpendicular-crossing residual ~1e-14). This orbit's extent
  (Uranus-Umbriel-barycentre distance ranging ~1.477-1.698 in Umbriel-SMA
  units) brackets Titania's own orbital radius (``a_gan ~ 1.640``) -- the
  same "apoapsis extends slightly beyond the perturber's own SMA" character
  `#689` documented for JEG's own Europa 3:4 seed (apoapsis ~1.610 vs
  ``a_gan~1.595``) and `#695` documented for Io-Europa's own 1:2 seed
  (apoapsis ~1.698 vs ``a_gan~1.591``).

Collision/near-collision risk check (`#696`'s own documented failure mode)
------------------------------------------------------------------------------
`#696` found that a naive "moon-ratio-matching" base-orbit choice for
Io-Ganymede threw a genuine ``RuntimeError`` (a hard near-collision
singularity in the direct Ganymede-acceleration term) because that orbit's
semimajor axis coincided almost exactly with Ganymede's own orbital radius.
This module's own base orbit is explicitly checked against the SAME risk
(:func:`test_base_orbit_clears_titania_over_full_torus_period`): propagating
the 1:2 orbit under Titania's PHYSICAL forcing over one full Titania-synodic
(``torus.period``) does NOT raise, and the minimum spacecraft-Titania
distance achieved is reported alongside the ALREADY-VALIDATED JEG positive
control's own minimum spacecraft-Ganymede distance (computed identically, in
this test module, from `#689`'s already-shipped ``jupiter_europa_ganymede_default``
+ its own literature-matched Europa 3:4 base orbit) as the calibration
reference for what "clears comfortably" means in this pipeline: JEG's own
proven-safe base orbit gets to ~0.0152 Europa-SMA units (well inside
Ganymede's own point-mass Hill-radius scale, ``mu_gan^(1/3)*a_gan ~ 0.068``)
WITHOUT raising or preventing convergence, so a close approach alone is not
disqualifying -- only an actual integrator failure (`#696`'s own operational
criterion) is. Live-observed: this system's own minimum spacecraft-Titania
distance (~0.058 Umbriel-SMA units) is COMFORTABLY LARGER (in relative,
Hill-radius-normalised terms) than JEG's own validated reference, and no
``RuntimeError`` occurs.
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy.integrate import solve_ivp

import cyclerfinder.core.ccr4bp as ccr4bp
import cyclerfinder.core.ccr4bp_umbriel_titania as ut
import cyclerfinder.core.cr3bp as cr3bp
import cyclerfinder.search.variational_ccr4bp_torus as vt
from cyclerfinder.genome.composed_moon_map import resonance_semimajor


def _resonant_symmetric_orbit(
    mu: float, p_sc: int, q_moon: int, *, max_iter: int = 80, tol: float = 1e-12, cap: float = 0.05
) -> tuple[np.ndarray, float, float]:
    """Identical test-only scaffolding to #690/#691/#694/#695/#696's own test
    modules (no production code under test)."""
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
    return ut.uranus_umbriel_titania_default()


@pytest.fixture(scope="module")
def orbit_12(system: ccr4bp.CCR4BPSystem) -> tuple[np.ndarray, float]:
    """The Uranus-Umbriel spacecraft:Umbriel = 1:2 symmetric resonant periodic
    orbit (base CR3BP) -- see module docstring for why this reading (not the
    interior 2:1 reading, which fails to converge with this scaffolding)."""
    s0, period, res = _resonant_symmetric_orbit(system.mu, 1, 2)
    assert res < 1e-10, f"resonant orbit did not converge: perp residual {res:.2e}"
    return s0, period


def test_interior_2to1_reading_does_not_converge(system: ccr4bp.CCR4BPSystem) -> None:
    """Documents the negative: the OTHER plausible reading of "2:1"
    (spacecraft:Umbriel = 2:1, interior, a0~0.63) does not converge with this
    damped-Newton scaffolding -- the same non-convergent reading `#695` found
    for Io-Europa's own interior 2:1 family, the reason the exterior (1:2)
    reading is used throughout this task."""
    _s0, _period, res = _resonant_symmetric_orbit(system.mu, 2, 1)
    assert res > 1e-6, f"expected non-convergence, got residual {res:.2e}"


def test_base_orbit_stays_clear_of_umbriel_and_near_titania_sma(
    system: ccr4bp.CCR4BPSystem, orbit_12: tuple[np.ndarray, float]
) -> None:
    """Sanity: the seed orbit's periapsis is far from Umbriel (no
    close-approach singularity in the BASE CR3BP) and its extent brackets
    Titania's own SMA (the elegant consequence of picking a
    spacecraft:Umbriel=1:2 resonance in a system where Umbriel:Titania are
    themselves ~2:1 commensurate)."""
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
    assert sol.sol is not None
    ts = np.linspace(0.0, period, 2000)
    xs = np.array([sol.sol(t) for t in ts])
    r_umbriel = np.sqrt((xs[:, 0] - 1.0 + mu) ** 2 + xs[:, 1] ** 2)
    r_uranus = np.sqrt((xs[:, 0] + mu) ** 2 + xs[:, 1] ** 2)
    assert r_umbriel.min() > 0.3  # far outside Umbriel's Hill sphere (~0.024)
    assert r_uranus.min() < system.a_gan < r_uranus.max() * 1.1  # brackets Titania's SMA


def test_base_orbit_clears_titania_over_full_torus_period(
    system: ccr4bp.CCR4BPSystem, orbit_12: tuple[np.ndarray, float]
) -> None:
    """`#696`'s own collision-risk check, applied here: propagate the 1:2
    orbit under Titania's PHYSICAL forcing over one full torus period
    (Titania's own synodic period) and confirm no RuntimeError, reporting the
    minimum spacecraft-Titania distance alongside the ALREADY-VALIDATED JEG
    positive control's own minimum spacecraft-Ganymede distance (computed
    identically here) as the calibration reference -- JEG's own proven-safe
    base orbit itself gets to ~0.0152 Europa-SMA units (well inside
    Ganymede's Hill-radius scale ~0.068) without failing, so a close approach
    alone does not disqualify; only an actual integrator failure does."""
    s0, _period = orbit_12
    torus_period = 2.0 * np.pi / abs(system.omega_gan)
    sol = solve_ivp(
        ccr4bp.ccr4bp_eom,
        (0.0, torus_period),
        s0,
        args=(system,),
        method="DOP853",
        rtol=1e-12,
        atol=1e-12,
        dense_output=True,
    )
    assert sol.sol is not None
    ts = np.linspace(0.0, torus_period, 4000)
    ys = sol.sol(ts)
    titania_pos = np.array([ccr4bp._ganymede_position(float(t), system) for t in ts])
    dist = np.linalg.norm(ys[:2].T - titania_pos[:, :2], axis=1)
    min_dist_umbriel_units = float(dist.min())

    # JEG's own calibration reference (already-validated positive control).
    jeg = ccr4bp.jupiter_europa_ganymede_default()
    jeg_s0, _jeg_period, jeg_res = _resonant_symmetric_orbit(jeg.mu, 3, 4)
    assert jeg_res < 1e-10
    jeg_torus_period = 2.0 * np.pi / abs(jeg.omega_gan)
    jeg_sol = solve_ivp(
        ccr4bp.ccr4bp_eom,
        (0.0, jeg_torus_period),
        jeg_s0,
        args=(jeg,),
        method="DOP853",
        rtol=1e-12,
        atol=1e-12,
        dense_output=True,
    )
    jeg_ts = np.linspace(0.0, jeg_torus_period, 4000)
    jeg_ys = jeg_sol.sol(jeg_ts)
    jeg_gan_pos = np.array([ccr4bp._ganymede_position(float(t), jeg) for t in jeg_ts])
    jeg_dist = np.linalg.norm(jeg_ys[:2].T - jeg_gan_pos[:, :2], axis=1)
    jeg_min_dist = float(jeg_dist.min())

    hill_titania = system.mu_gan ** (1.0 / 3.0) * system.a_gan
    hill_ganymede = jeg.mu_gan ** (1.0 / 3.0) * jeg.a_gan

    # No collision singularity for this system's own base orbit.
    assert min_dist_umbriel_units > 0.0
    assert min_dist_umbriel_units == pytest.approx(0.0578, abs=0.01)
    # Comparable-or-better margin (relative to the perturber's own Hill scale)
    # than the already-validated, proven-safe JEG reference.
    assert (min_dist_umbriel_units / hill_titania) >= (jeg_min_dist / hill_ganymede) * 0.5


# ---------------------------------------------------------------------------
# Positive control: mu_gan=0 -> theta1-flat torus.
# ---------------------------------------------------------------------------


def test_mu_gan_zero_flat_regression(
    system: ccr4bp.CCR4BPSystem, orbit_12: tuple[np.ndarray, float]
) -> None:
    """At mu_gan=0 the corrector recovers the base 1:2 orbit as a theta1-FLAT
    torus -- the corrector-level analogue of #689's mu_gan->0 structural
    reduction, on THIS new system (mirrors #690's/#695's/#696's own test
    exactly)."""
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
    """The 1:2 Umbriel resonant orbit persists as a genuine thin 2D invariant
    torus under Titania's PHYSICAL-mass forcing -- the actual discovery-attempt
    object this task builds on. Reports invariance residual, an INDEPENDENT
    closure residual, the stroboscopic rotation number, and theta1-forcing
    content appearing only when mu_gan > 0 (mirrors #690's/#695's/#696's own
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
    assert _theta1_content(phys) > 1e-8
    expected = 2.0 * np.pi * (2.0 * np.pi / period) / abs(system.omega_gan)
    assert phys.rho_strob == pytest.approx(expected, rel=1e-3)


def test_continue_ccr4bp_torus_mass_reaches_physical(
    system: ccr4bp.CCR4BPSystem, orbit_12: tuple[np.ndarray, float]
) -> None:
    """Mass homotopy from mu_gan=0 walks to the physical mass without diverging."""
    s0, period = orbit_12
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
