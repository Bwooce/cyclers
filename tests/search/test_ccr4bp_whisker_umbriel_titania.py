"""CCR4BP whisker/manifold-direction extraction tests for Uranus-Umbriel-Titania (`#701`).

Applies `#691`'s segment-anchored discrete-QR/CLV extraction
(``search.ccr4bp_whisker``, reused UNMODIFIED) to the Umbriel-Titania 1:2
torus (see ``test_ccr4bp_torus_umbriel_titania.py`` for why the tractable base
orbit is spacecraft:Umbriel = 1:2, exterior). Mirrors `#691`'s own go/no-go
discipline: run the perturbation-robustness diagnostic on this NEW system
before trusting the cheap segmented-CLV direction for the downstream
manifold-globalization + heteroclinic search (`#694`'s modules, run in the
discovery-attempt driver script).

Live-observed headline result: this torus's one-period unstable eigenvalue
ranges ``|lam_u| ~ 6.04-13.45`` across the 12 checked phases (dominant real
hyperbolic pair, separation ratio == ``|lam_u|`` since no other eigenvalue
approaches its magnitude at any checked phase) -- comparable order to JEG's
own reference system. The perturbation-robustness diagnostic is CLEANLY
trustworthy here: worst one-shot swing ~0.20 deg, worst segmented swing
~0.0035 deg (both far below the 10-degree go/no-go bar this project uses, and
the segmented method is, as expected, ~60x tighter than the one-shot method).
The expensive full Kumar bundle-solve whisker construction is NOT indicated.
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy.integrate import solve_ivp

import cyclerfinder.core.ccr4bp as ccr4bp
import cyclerfinder.core.ccr4bp_umbriel_titania as ut
import cyclerfinder.core.cr3bp as cr3bp
import cyclerfinder.search.ccr4bp_whisker as wk
import cyclerfinder.search.variational_ccr4bp_torus as vt
from cyclerfinder.genome.composed_moon_map import resonance_semimajor


def _resonant_symmetric_orbit(
    mu: float, p_sc: int, q_moon: int, *, max_iter: int = 80, tol: float = 1e-12, cap: float = 0.05
) -> tuple[np.ndarray, float, float]:
    """Identical test-only scaffolding to #690/#691/#694/#695/#696's own test modules."""
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


@pytest.fixture(scope="module")
def system() -> ccr4bp.CCR4BPSystem:
    return ut.uranus_umbriel_titania_default()


@pytest.fixture(scope="module")
def orbit_12(system: ccr4bp.CCR4BPSystem) -> tuple[np.ndarray, float]:
    s0, period, res = _resonant_symmetric_orbit(system.mu, 1, 2)
    assert res < 1e-10, f"resonant orbit did not converge: perp residual {res:.2e}"
    return s0, period


@pytest.fixture(scope="module")
def phys_torus(
    system: ccr4bp.CCR4BPSystem, orbit_12: tuple[np.ndarray, float]
) -> vt.CCR4BPTorusVariationalResult:
    s0, period = orbit_12
    return vt.discover_ccr4bp_torus_from_resonant_orbit(
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


# ---------------------------------------------------------------------------
# Step 1: one-period amplification report on the physical-mass torus.
# ---------------------------------------------------------------------------


def test_amplification_diagnostics_physical_torus(
    phys_torus: vt.CCR4BPTorusVariationalResult,
) -> None:
    """A clear, well-separated unstable direction at every phase checked.
    Live-observed: |lam_u| ~6.04-13.45 across phases, separation ratio equal
    to |lam_u| (no competing eigenvalue at any checked phase)."""
    checked = 0
    for theta1 in np.linspace(0.0, 2 * np.pi, 4, endpoint=False):
        for theta2 in np.linspace(0.0, 2 * np.pi, 3, endpoint=False):
            state0, stm = wk.one_period_stm(phys_torus, float(theta1), float(theta2))
            diag = wk.amplification_diagnostics(state0, stm)
            assert diag.lam_u is not None and abs(diag.lam_u) > 1.0, (theta1, theta2)
            assert 5.0 < abs(diag.lam_u) < 15.0, (theta1, theta2, diag.lam_u)
            assert diag.lam_s is not None
            assert diag.lam_u * diag.lam_s == pytest.approx(1.0, abs=5e-3)
            assert diag.separation_ratio > 5.0, (theta1, theta2, diag.separation_ratio)
            checked += 1
    assert checked == 12


# ---------------------------------------------------------------------------
# Step 2: the actual go/no-go perturbation-robustness diagnostic.
# ---------------------------------------------------------------------------


def test_perturbation_robustness_physical_torus(
    phys_torus: vt.CCR4BPTorusVariationalResult,
) -> None:
    """THE GO/NO-GO DIAGNOSTIC, run on this NEW system before trusting the
    cheap segmented-CLV direction for the downstream heteroclinic search.
    Checks an exploratory-sweep worst phase plus two representative ones,
    mirroring #691's/#695's/#696's own discipline exactly. Live-observed
    worst-case swings (5 seeds, 3 phases x 2 branches): one-shot ~0.20 deg,
    segmented ~0.0035 deg -- BOTH methods trustworthy."""
    phases = [(0.0, 0.0), (1.0, 2.0), (2.0943951023931953, 4.1887902047863905)]
    for theta1, theta2 in phases:
        for branch in ("unstable", "stable"):
            worst_one_shot = 0.0
            worst_segmented = 0.0
            for seed in range(5):
                s1 = wk.perturbation_swing_deg(
                    phys_torus, theta1, theta2, branch, method="one_shot", seed=seed
                )
                s2 = wk.perturbation_swing_deg(
                    phys_torus, theta1, theta2, branch, method="segmented", n_segments=32, seed=seed
                )
                if s1 is not None:
                    worst_one_shot = max(worst_one_shot, s1)
                if s2 is not None:
                    worst_segmented = max(worst_segmented, s2)
            assert worst_one_shot < 5.0, (theta1, theta2, branch, worst_one_shot)
            assert worst_segmented < 0.1, (theta1, theta2, branch, worst_segmented)
            assert worst_segmented <= worst_one_shot + 1e-9, (theta1, theta2, branch)
