"""CCR4BP torus positive control for Jupiter-Europa-Callisto (`#703`).

Applies `#690`'s pseudospectral quasi-periodic-torus corrector
(``search.variational_ccr4bp_torus``, reused UNMODIFIED) to `#703`'s new
Europa-Callisto system
(``core.ccr4bp_europa_callisto.jupiter_europa_callisto_default``) as a genuine
discovery attempt, not another capability-proof exercise. Mirrors `#690`'s/
`#695`'s/`#696`'s/`#701`'s own test discipline exactly: the ``mu_gan=0``
regression is a machine-precision-flat structural check; the physical-mass
torus is reported at its achieved invariance residual.

Base resonant orbit choice
---------------------------
Europa:Callisto's real period ratio is ~4.734 -- unlike every other CCR4BP
pair this project has built, this is NOT close to any clean low-integer
commensurability (best continued-fraction convergents: 4/1 off by 0.70, 5/1
off by 0.30; `#700`'s own deep literature/tractability note,
``docs/notes/2026-07-24-700-europa-callisto-deep-litcheck.md``). `#700`'s own
quick scratch check (not committed) found every seed ``(p_sc, q_moon)`` tried
converges to essentially the SAME near-degenerate near-circular family; this
module independently re-derives that finding rather than trusting the citation
(this project's own "reproduce, don't just cite" discipline).

* The INTERIOR spacecraft:Europa = 4:1 reading (``a0 ~ 0.397``, following
  `#696`'s own Io-Ganymede precedent for a near-4:1-labeled system) converges
  cleanly (perpendicular-crossing residual ~2.2e-14) and clears Callisto's own
  orbit by a WIDE margin -- minimum spacecraft-Callisto distance over one full
  torus period is ~2.409 Europa-SMA units, ~22.4x Callisto's own point-mass
  Hill-radius scale (``mu_gan^(1/3)*a_gan ~ 0.108``), a far more comfortable
  ratio than JEG's own proven-safe reference (~0.22x its own Hill scale) or
  even Umbriel-Titania's own (`#701`, comfortably-safe but nowhere near this
  wide -- see that module's own test).
* The naive "moon-ratio-matching" EXTERIOR spacecraft:Europa = 1:5 reading
  (``a0 ~ 2.924``, exterior to Callisto's own SMA ``a_gan ~ 2.805`` by a
  static-SMA measure) does NOT throw a hard ``RuntimeError`` under Callisto's
  physical forcing (unlike Io-Ganymede's own naive 1:4 choice, `#696`'s own
  documented hard near-collision singularity) -- but this task's own
  numerical re-check finds its minimum spacecraft-Callisto approach distance
  is only ~1.10x Callisto's Hill-radius scale, a THIN, uncomfortable margin
  compared to the interior 4:1 choice's ~22x. This module therefore uses the
  interior 4:1 reading throughout, per this task's own dispatch instruction
  to avoid the naive exterior choice and per this project's collision-risk
  discipline (`#696`'s own documented failure mode), even though the naive
  choice here happens not to be a hard integrator failure the way Io-Ganymede's
  was (:func:`test_naive_exterior_1to5_reading_has_thin_but_nonfatal_margin`).
"""

from __future__ import annotations

import contextlib

import numpy as np
import pytest
from scipy.integrate import solve_ivp

import cyclerfinder.core.ccr4bp as ccr4bp
import cyclerfinder.core.ccr4bp_europa_callisto as ec
import cyclerfinder.core.cr3bp as cr3bp
import cyclerfinder.search.variational_ccr4bp_torus as vt
from cyclerfinder.genome.composed_moon_map import resonance_semimajor


def _resonant_symmetric_orbit(
    mu: float, p_sc: int, q_moon: int, *, max_iter: int = 80, tol: float = 1e-12, cap: float = 0.05
) -> tuple[np.ndarray, float, float]:
    """Identical test-only scaffolding to #690/#691/#694/#695/#696/#701's own
    test modules (no production code under test)."""
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
    return ec.jupiter_europa_callisto_default()


@pytest.fixture(scope="module")
def orbit_4to1(system: ccr4bp.CCR4BPSystem) -> tuple[np.ndarray, float]:
    """The Jupiter-Europa spacecraft:Europa = 4:1 symmetric resonant periodic
    orbit (base CR3BP, INTERIOR, a~0.397) -- see module docstring for why this
    reading, not the naive exterior 1:5 reading."""
    s0, period, res = _resonant_symmetric_orbit(system.mu, 4, 1)
    assert res < 1e-10, f"resonant orbit did not converge: perp residual {res:.2e}"
    return s0, period


def test_base_orbit_stays_clear_of_europa_and_interior_to_callisto(
    system: ccr4bp.CCR4BPSystem, orbit_4to1: tuple[np.ndarray, float]
) -> None:
    """Sanity: the seed orbit's periapsis is far from Europa (no close-approach
    singularity in the BASE CR3BP) and stays comfortably INTERIOR to Callisto's
    own SMA (unlike JEG's/Umbriel-Titania's own seed orbits, which bracket
    their perturber's SMA -- the interior 4:1 reading has no reason to)."""
    s0, period = orbit_4to1
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
    r_europa = np.sqrt((xs[:, 0] - 1.0 + mu) ** 2 + xs[:, 1] ** 2)
    r_jup = np.sqrt((xs[:, 0] + mu) ** 2 + xs[:, 1] ** 2)
    assert r_europa.min() > 0.3  # far outside Europa's own Hill sphere (~0.008)
    assert r_jup.max() < 0.5 * system.a_gan  # comfortably interior to Callisto's SMA


def test_base_orbit_clears_callisto_over_full_torus_period_with_wide_margin(
    system: ccr4bp.CCR4BPSystem, orbit_4to1: tuple[np.ndarray, float]
) -> None:
    """`#696`'s own collision-risk check, applied here: propagate the 4:1
    orbit under Callisto's PHYSICAL forcing over one full torus period
    (Callisto's own synodic period) and confirm no RuntimeError, reporting the
    minimum spacecraft-Callisto distance -- live-observed ~2.409 Europa-SMA
    units, ~22.4x Callisto's own point-mass Hill-radius scale, WIDER (in
    relative, Hill-radius-normalised terms) than every other system this
    project has built (JEG's own proven-safe reference ~0.22x its Hill scale,
    Umbriel-Titania's own comfortably-safe-but-narrower margin), consistent
    with `#700`'s own quick citation (~2.32-2.46 Europa-SMA gap for interior
    candidates, larger than Io-Ganymede's own 4:1 choice's 2.14 Io-SMA-unit
    margin) -- re-derived here directly, not merely trusted."""
    s0, _period = orbit_4to1
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
    callisto_pos = np.array([ccr4bp._ganymede_position(float(t), system) for t in ts])
    dist = np.linalg.norm(ys[:2].T - callisto_pos[:, :2], axis=1)
    min_dist_europa_units = float(dist.min())

    hill_callisto = system.mu_gan ** (1.0 / 3.0) * system.a_gan

    assert min_dist_europa_units > 0.0
    assert min_dist_europa_units == pytest.approx(2.4085, abs=0.05)
    assert (min_dist_europa_units / hill_callisto) > 10.0  # a wide, comfortable margin


def test_naive_exterior_1to5_reading_has_thin_but_nonfatal_margin(
    system: ccr4bp.CCR4BPSystem,
) -> None:
    """Documents the negative that motivates using the interior 4:1 reading
    instead: the naive "moon-ratio-matching" spacecraft:Europa=1:5 orbit
    (exterior, a0~2.924, sits static-SMA-outside Callisto's own a_gan~2.805)
    converges fine as a bare CR3BP periodic orbit AND does not throw a hard
    RuntimeError under Callisto's physical forcing (unlike Io-Ganymede's own
    naive 1:4 choice, `#696`'s own hard near-collision singularity) -- but its
    minimum spacecraft-Callisto approach distance is only ~1.10x Callisto's
    own Hill-radius scale, a thin, uncomfortable margin an order of magnitude
    narrower than the chosen interior 4:1 reading's ~22x. Falsifiable, not
    asserted: this project's own "reproduce, don't just cite" discipline
    applied to `#700`'s own static-SMA-only citation."""
    s0, _period, res = _resonant_symmetric_orbit(system.mu, 1, 5)
    assert res < 1e-10  # the bare periodic orbit itself IS a valid solution
    a_1to5 = resonance_semimajor(1, 5)
    assert a_1to5 > system.a_gan  # static-SMA "exterior" per #700's own citation

    torus_period = 2.0 * np.pi / abs(system.omega_gan)
    arc_or_none = None
    with contextlib.suppress(RuntimeError):
        arc_or_none = ccr4bp.propagate_ccr4bp(system, s0, torus_period, with_stm=False)
    assert arc_or_none is not None  # does NOT hard-fail, unlike Io-Ganymede's own 1:4 orbit

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
    ts = np.linspace(0.0, torus_period, 8000)
    ys = sol.sol(ts)
    callisto_pos = np.array([ccr4bp._ganymede_position(float(t), system) for t in ts])
    dist = np.linalg.norm(ys[:2].T - callisto_pos[:, :2], axis=1)
    min_dist_europa_units = float(dist.min())
    hill_callisto = system.mu_gan ** (1.0 / 3.0) * system.a_gan
    ratio = min_dist_europa_units / hill_callisto
    assert ratio == pytest.approx(1.10, abs=0.1)
    assert ratio < 2.0  # thin margin, motivating the interior 4:1 choice instead


# ---------------------------------------------------------------------------
# Positive control: mu_gan=0 -> theta1-flat torus.
# ---------------------------------------------------------------------------


def test_mu_gan_zero_flat_regression(
    system: ccr4bp.CCR4BPSystem, orbit_4to1: tuple[np.ndarray, float]
) -> None:
    """At mu_gan=0 the corrector recovers the base 4:1 orbit as a theta1-FLAT
    torus -- the corrector-level analogue of #689's mu_gan->0 structural
    reduction, on THIS new system (mirrors #690's/#695's/#696's/#701's own
    test exactly)."""
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
    """The 4:1 Europa resonant orbit persists as a genuine thin 2D invariant
    torus under Callisto's PHYSICAL-mass forcing -- the actual discovery-
    attempt object this task builds on. Live-observed: this torus converges
    VERY cleanly (physical-mass residual_rms ~7.9e-7), comparable to or
    cleaner than every other system this project has built, because the base
    orbit is near-circular (mirrors #696's own finding for Io-Ganymede's own
    interior 4:1 choice)."""
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
