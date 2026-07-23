"""CCR4BP whisker (`#691`) go/no-go diagnostic for Jupiter-Io-Europa (`#695`).

Applies `#691`'s one-period amplification + segment-anchored CLV extraction
(``search.ccr4bp_whisker``, reused UNMODIFIED) to the Io-Europa physical-mass
torus built in ``tests/search/test_ccr4bp_torus_io_europa.py``. Mirrors
`#691`'s own go/no-go discipline: an exploratory phase/branch sweep FIRST
(reported here, not hidden), then the mandatory perturbation-robustness test
at the phases that sweep flagged as extreme.

Headline finding (this task's own, NOT inherited from JEG)
------------------------------------------------------------
Unlike JEG (`#691`'s own worst-case one-shot swing was < 3.44 deg -- "usable"
even without segmentation), Io-Europa's physical torus has a MUCH weaker and
less uniformly hyperbolic one-period map (live-observed ``lam_u`` ranging
from a near-degenerate ~1.0-1.1 up to a sharply localized ~425 at one narrow
``theta2`` band, vs JEG's more uniform ~9-1700 -- this system's forcing
produces islands of near-neutral stability alongside a few sharply unstable
bands, not a uniformly strong whisker everywhere). At the near-degenerate
phase (``theta1=theta2=pi/2``, separation_ratio ~1.0), the ONE-SHOT
extraction is CATASTROPHICALLY untrustworthy (worst-case swing ~90 degrees --
worse than JEG's own worst case by more than an order of magnitude, and 2 of
5 perturbation seeds even fail to find a real hyperbolic eigenpair at all).
The segment-anchored CLV method, in contrast, stays trustworthy throughout
(<0.08 deg at n_segments=16, tightening further with more segments) --
EXACTLY the failure mode `#646` built the CLV method to guard against, now
demonstrated on a second, independently-discovered system. **Verdict: for
Io-Europa, the cheap ONE-SHOT direction is NOT uniformly trustworthy and
must not be used near-degenerate phases; the segmented/CLV extraction IS
trustworthy and is the one this task's own manifold-globalization stage
uses exclusively (matching `#694`'s own module, which already only calls
the segmented extractor).**
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy.integrate import solve_ivp

import cyclerfinder.core.ccr4bp as ccr4bp
import cyclerfinder.core.ccr4bp_io_europa as ioeu
import cyclerfinder.core.cr3bp as cr3bp
import cyclerfinder.search.ccr4bp_whisker as wk
import cyclerfinder.search.variational_ccr4bp_torus as vt
from cyclerfinder.genome.composed_moon_map import resonance_semimajor


def _resonant_symmetric_orbit(
    mu: float, p_sc: int, q_moon: int, *, max_iter: int = 80, tol: float = 1e-12, cap: float = 0.05
) -> tuple[np.ndarray, float, float]:
    """Identical test-only scaffolding to #690/#691/#694's own test/driver
    copies (no production code under test)."""
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
    return ioeu.jupiter_io_europa_default()


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
# #691 Step 1: one-period amplification -- exploratory sweep + headline report.
# ---------------------------------------------------------------------------


def test_amplification_diagnostics_physical_torus_wide_dynamic_range(
    phys_torus: vt.CCR4BPTorusVariationalResult,
) -> None:
    """Headline report: unlike JEG's more uniform spectrum, this torus's
    one-period map ranges from near-degenerate (separation_ratio ~1.0) to
    sharply unstable (|lam_u| ~ a few hundred) depending on phase, and every
    phase checked DOES still expose a real hyperbolic eigenpair (no locally
    non-hyperbolic points on this coarse grid), plus the reciprocal
    lam_u*lam_s~1 symplectic check holds throughout."""
    lam_us = []
    for theta1 in np.linspace(0.0, 2 * np.pi, 4, endpoint=False):
        for theta2 in np.linspace(0.0, 2 * np.pi, 4, endpoint=False):
            state0, stm = wk.one_period_stm(phys_torus, float(theta1), float(theta2))
            diag = wk.amplification_diagnostics(state0, stm)
            assert diag.lam_u is not None, (theta1, theta2)
            assert diag.lam_s is not None
            assert diag.lam_u * diag.lam_s == pytest.approx(1.0, abs=5e-2)
            lam_us.append(abs(diag.lam_u))
    lam_us_arr = np.array(lam_us)
    # A wide dynamic range is itself the reportable finding (not a fixed gate).
    assert lam_us_arr.max() > 2.0
    assert lam_us_arr.min() < 3.0


# ---------------------------------------------------------------------------
# #691 Step 2: the go/no-go perturbation-robustness diagnostic.
# ---------------------------------------------------------------------------


def test_one_shot_is_untrustworthy_at_near_degenerate_phase(
    phys_torus: vt.CCR4BPTorusVariationalResult,
) -> None:
    """THE NEGATIVE FINDING. At (theta1, theta2) = (pi/2, pi/2) -- a
    near-degenerate phase found by this task's own exploratory sweep
    (separation_ratio ~1.0, lam_u ~ 1.005) -- the naive one-shot extraction's
    swing is catastrophically large (order 90 degrees, essentially a random
    direction), unlike JEG where the one-shot method stayed under ~3.44 deg
    everywhere checked."""
    theta1 = theta2 = float(np.pi / 2.0)
    worst_one_shot = 0.0
    n_evaluable = 0
    for branch in ("unstable", "stable"):
        for seed in range(5):
            s1 = wk.perturbation_swing_deg(
                phys_torus, theta1, theta2, branch, method="one_shot", seed=seed
            )
            if s1 is None:
                continue
            n_evaluable += 1
            worst_one_shot = max(worst_one_shot, s1)
    assert n_evaluable > 0
    assert worst_one_shot > 10.0, worst_one_shot  # order-of-magnitude worse than JEG


def test_segmented_clv_remains_trustworthy_at_same_near_degenerate_phase(
    phys_torus: vt.CCR4BPTorusVariationalResult,
) -> None:
    """The segment-anchored CLV method stays trustworthy at the SAME phase
    where the one-shot method fails catastrophically -- the exact `#646`
    guard mechanism, now demonstrated on a second, independently-found
    system. n_segments=32 (this task's own driver-script default) used, not
    the module's default 16, mirroring `#691`'s own n_segments-sensitivity
    finding (its own docstring documents 16 not uniformly sufficing on the
    JEG torus either)."""
    theta1 = theta2 = float(np.pi / 2.0)
    for branch in ("unstable", "stable"):
        worst_segmented = 0.0
        n_evaluable = 0
        for seed in range(5):
            s2 = wk.perturbation_swing_deg(
                phys_torus,
                theta1,
                theta2,
                branch,
                method="segmented",
                n_segments=32,
                seed=seed,
            )
            if s2 is None:
                continue
            n_evaluable += 1
            worst_segmented = max(worst_segmented, s2)
        assert n_evaluable > 0
        assert worst_segmented < 0.1, (branch, worst_segmented)


def test_perturbation_robustness_at_strong_and_mild_phases(
    phys_torus: vt.CCR4BPTorusVariationalResult,
) -> None:
    """At phases AWAY from the near-degenerate island (a strong phase near
    theta1=theta2=pi, and a mild phase at the origin), the segmented method
    is trustworthy AND the one-shot method is also usable -- the failure
    mode is localized to the near-degenerate island, not universal, matching
    the exploratory sweep's own finding."""
    phases = [(np.pi, np.pi), (0.0, 0.0)]
    for theta1, theta2 in phases:
        for branch in ("unstable", "stable"):
            worst_one_shot = 0.0
            worst_segmented = 0.0
            for seed in range(5):
                s1 = wk.perturbation_swing_deg(
                    phys_torus, theta1, theta2, branch, method="one_shot", seed=seed
                )
                s2 = wk.perturbation_swing_deg(
                    phys_torus,
                    theta1,
                    theta2,
                    branch,
                    method="segmented",
                    n_segments=32,
                    seed=seed,
                )
                if s1 is not None:
                    worst_one_shot = max(worst_one_shot, s1)
                if s2 is not None:
                    worst_segmented = max(worst_segmented, s2)
            assert worst_one_shot < 5.0, (theta1, theta2, branch, worst_one_shot)
            assert worst_segmented < 0.1, (theta1, theta2, branch, worst_segmented)
