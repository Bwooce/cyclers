"""CCR4BP heteroclinic/homoclinic search regression test for Umbriel-Titania (`#701`).

`#701` ran the full pipeline (``scripts/screen_701_ccr4bp_umbriel_titania_search.py``)
as a genuine discovery attempt on the Umbriel-Titania physical torus. This test does
NOT re-run the expensive coarse KD-tree scan (that lives only in the driver script,
mirroring `#694`'s/`#695`'s/`#696`'s own screen scripts, likewise untested by a
dedicated pytest module -- their "test" IS the committed result.json). Instead it
regression-locks the driver script's own headline findings by re-deriving them
directly from the SAME seed phases the coarse search found, using `#694`'s
``refine_candidate``/``ghost_guard`` (as fixed by `#702`).

**The headline finding, AFTER `#702`'s ghost_guard ref_vec-anchoring fix: TWO
independently-confirmed genuine connections, not one near-miss plus one
suspected ghost artifact.**

1. A ROBUST, comfortably ghost-guard-verified near-miss family exists,
   reproducing near-identically across all 4 independent lobe-branch combos:
   best ``corrected_pos_gap_km ~43.1`` km / ``corrected_vel_gap_km_s ~2.58``
   m/s, with ``corrected_integrator_delta_km ~1.2e-7`` km -- SIX orders of
   magnitude below the 1.0 km ghost-guard gate (the SAME comfortable margin
   JEG's own positive control achieved, ~6.5e-7 km). This is analogous to
   (about 2.6x the magnitude of) `#695`'s own Io-Europa near-miss (16.7 km /
   0.54 m/s).
2. Separately, several candidates refine to NEAR-MACHINE-PRECISION residuals
   (``residual_norm ~1e-14``, ``pos_gap_corrected ~1e-9`` km). `#701`'s own
   FIRST-PASS investigation (before `#702`) flagged these as a suspected
   ghost artifact: the module-native ``ghost_guard`` Radau/DOP853 independent-
   integrator cross-check disagreed by up to ~0.98-21.5 km depending on lobe
   choice, right at the edge of (or over) the 1.0 km rejection gate.
   `#702` found and fixed the ROOT CAUSE: ``ghost_guard`` was re-deriving its
   Radau-side CLV sign-anchor (``ref_vec``) at the FINAL converged ``theta2``
   instead of reusing the SEED-anchored ``ref_vec`` that
   ``refine_candidate``'s own DOP853 states were actually computed with --
   the raw CLV sign is discontinuous between the seed and converged phase
   here, so the old code's Radau re-check was silently stepping the WRONG
   manifold lobe, comparing two DIFFERENT trajectories and reporting a large,
   physically-meaningless disagreement that had nothing to do with genuine
   chaos-amplified integrator sensitivity. With the fix (``ghost_guard`` now
   reuses ``RefinedConnection.ref_vec_u``/``ref_vec_s``, the SAME seed-
   anchored vectors ``refine_candidate`` itself used), this family's
   independent-integrator agreement is UNIFORMLY tiny (``~1e-7`` to ``~1e-10``
   km across every lobe combo and every seed that lands in this basin) --
   the SAME comfortable, six-plus-order-of-magnitude margin the ROBUST family
   and JEG's own positive control show, NOT a fragile coincidence. This is
   confirmed further by `#701`'s own re-run's seed-perturbation stability
   check (two different starting seeds for the same nominal candidate now
   agree to within ~5x, not ~1e7x) and mesh-refinement re-check (the
   candidate persists under a 4x-denser globalization grid). This family is
   the first apparently-novel discovery-grade CCR4BP homoclinic connection in
   this project's arc -- see
   ``data/found/701_ccr4bp_umbriel_titania_search/result.json``'s
   ``best_genuine_connection_corrected``/``best_robust_genuine_connection_corrected``
   (identical entries: this candidate clears BOTH the raw and the strict
   robust gate).

Physical-unit correction: `#694`'s own ``_L_KM``/``_v_unit_km_s`` are
hardcoded to the JEG system (Europa's SMA); this task's driver script
back-corrects to Umbriel's own SMA-based unit
(``core.ccr4bp_umbriel_titania.L_KM``) via pure post-hoc rescaling (see the
driver script's own module docstring for the full derivation) -- reproduced
here identically so this test stays consistent with the committed
``data/found/701_ccr4bp_umbriel_titania_search/result.json``.
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy.integrate import solve_ivp

import cyclerfinder.core.ccr4bp_umbriel_titania as ut
import cyclerfinder.core.cr3bp as cr3bp
import cyclerfinder.search.ccr4bp_heteroclinic_search as hs
import cyclerfinder.search.variational_ccr4bp_torus as vt
from cyclerfinder.genome.composed_moon_map import resonance_semimajor

# #694's own hardcoded (Europa-SMA-based) physical-unit constant -- needed
# ONLY as a sanity cross-check against core.ccr4bp_umbriel_titania.L_KM
# (see the driver script's own module docstring).
_JEG_L_KM = 671_100.0


def _resonant_symmetric_orbit(
    mu: float, p_sc: int, q_moon: int, *, max_iter: int = 80, tol: float = 1e-12, cap: float = 0.05
) -> tuple[np.ndarray, float, float]:
    """Identical test-only scaffolding to #690/#691/#694/#695/#696's own copies."""
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
def phys_torus() -> vt.CCR4BPTorusVariationalResult:
    system = ut.uranus_umbriel_titania_default()
    s0, period, res = _resonant_symmetric_orbit(system.mu, 1, 2)
    assert res < 1e-10
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


# #737 (2026-07-27): confirmed CI-resource-budget mismatch, not a code
# regression -- this test times out (>600s pytest-timeout) on CI (run
# 30261678515) but reproduces cleanly in ~70s locally on this 8-core Mac
# (well under budget, ~8.5x margin), consistent with #631's own documented
# precedent (CI runners have only 2 cores per pyproject.toml's own comment;
# full-suite -n auto parallel contention on a shared 2-core runner inflates
# wall-clock for CPU-heavy scipy/numpy tests -- here, refining across all 4
# independent lobe-branch combos plus the Radau/DOP853 ghost_guard cross-
# check -- far more than core-count alone would suggest). This file's other
# test, test_near_exact_candidate_confirmed_genuine_after_702_ref_vec_fix,
# runs in ~7s locally and did NOT time out on this same CI run, so it is
# left unmarked. Not previously classified -- first CI run of this exact
# code since #702's own ghost_guard fix.
@pytest.mark.slow
def test_robust_near_miss_is_the_trustworthy_headline_result(
    phys_torus: vt.CCR4BPTorusVariationalResult,
) -> None:
    """Regression-locks #701's own driver-script ROBUST finding: the seed
    (theta2_u=5.759586531581288, t_u=22.69457455316767,
    theta2_s=5.654866776461628, t_s=2.253362012371258, u_lobe=+1, s_lobe=-1
    -- the exact coarse-search hit from
    ``data/found/701_ccr4bp_umbriel_titania_search/result.json``'s
    ``best_robust_genuine_connection_corrected``) refines to a genuine LOCAL
    MINIMUM (``residual_norm`` NOT machine precision) with a comfortable
    ghost-guard pass (integrator_delta_km many orders of magnitude below the
    1.0 km gate -- the SAME comfortable-margin pattern JEG's own positive
    control showed)."""
    cand = hs.ManifoldCandidate(
        theta2_u=5.759586531581288,
        t_u=22.69457455316767,
        theta2_s=5.654866776461628,
        t_s=2.253362012371258,
        gap_planar=0.0010789032877761132,
    )
    refined = hs.refine_candidate(
        phys_torus, phys_torus, cand, lobe_sign_u=1.0, lobe_sign_s=-1.0, n_segments_dir=32
    )
    assert refined is not None
    # NOT machine precision (unlike JEG's own best ~6.8e-15): a genuine nonzero floor.
    assert refined.residual_norm > 1e-5, refined.residual_norm
    assert refined.residual_norm < 1e-2, refined.residual_norm

    guard = hs.ghost_guard(
        phys_torus, phys_torus, refined, lobe_sign_u=1.0, lobe_sign_s=-1.0, n_segments_dir=32
    )
    length_scale = ut.L_KM / _JEG_L_KM
    corrected_pos_gap_km = refined.pos_gap_km * length_scale
    corrected_off_torus_km = guard.off_torus_km * length_scale
    corrected_integrator_delta_km = guard.integrator_delta_km * length_scale

    assert corrected_pos_gap_km == pytest.approx(43.137, abs=0.05)
    assert corrected_integrator_delta_km < 1e-5  # comfortably clears the 1.0 km gate
    assert corrected_off_torus_km >= 1000.0  # genuinely far from the torus, not trivial
    corrected_genuine = corrected_integrator_delta_km < 1.0 and corrected_off_torus_km >= 1000.0
    assert corrected_genuine is True


def test_near_exact_candidate_confirmed_genuine_after_702_ref_vec_fix(
    phys_torus: vt.CCR4BPTorusVariationalResult,
) -> None:
    """Regression for `#702`: the near-machine-precision "exact" family
    (``residual_norm~1e-14``) is CONFIRMED genuine, not a suspected ghost
    artifact -- `#701`'s own first-pass ``ghost_sensitivity_check`` finding
    (two nominally-identical candidates showing wildly different Radau
    agreement, ~0 km vs ~0.98 km) was itself an artifact of `#694`'s
    ``ghost_guard`` re-deriving its Radau-side CLV ``ref_vec`` anchor at the
    FINAL converged ``theta2`` instead of reusing the SEED-anchored vector
    ``refine_candidate`` actually used -- exactly the case where the raw CLV
    sign flips between seed and converged phase, so the old code's Radau
    re-check silently stepped the WRONG manifold lobe on some (but not all)
    lobe/seed combinations. With the fix (``ghost_guard`` now reuses
    ``RefinedConnection.ref_vec_u``/``ref_vec_s``), this exact candidate
    (`#701`'s own ``best_genuine_connection_corrected`` /
    ``best_robust_genuine_connection_corrected`` seed) refines via the
    NORMAL, unmodified ``refine_candidate`` path and passes ``ghost_guard``
    with a comfortable, JEG-positive-control-like margin -- confirming the
    apparent "fragility" was the bug, not real chaos-amplified sensitivity."""
    seed = hs.ManifoldCandidate(
        theta2_u=3.665191429188092,
        t_u=19.314531534610783,
        theta2_s=3.5604716740684323,
        t_s=18.187850528425155,
        gap_planar=0.0011381326404778097,
    )
    refined = hs.refine_candidate(
        phys_torus, phys_torus, seed, lobe_sign_u=-1.0, lobe_sign_s=-1.0, n_segments_dir=32
    )
    assert refined is not None
    assert refined.converged
    # Near-machine-precision closure on the DOP853 side (loose bound, pinning
    # behaviour not an exact float -- see this project's own testing convention).
    assert refined.residual_norm < 1e-10, refined.residual_norm

    guard = hs.ghost_guard(
        phys_torus, phys_torus, refined, lobe_sign_u=-1.0, lobe_sign_s=-1.0, n_segments_dir=32
    )
    length_scale = ut.L_KM / _JEG_L_KM
    corrected_pos_gap_km = refined.pos_gap_km * length_scale
    corrected_off_torus_km = guard.off_torus_km * length_scale
    corrected_integrator_delta_km = guard.integrator_delta_km * length_scale

    assert corrected_pos_gap_km < 1e-6, corrected_pos_gap_km
    assert corrected_off_torus_km >= 1000.0, corrected_off_torus_km
    # Comfortable, JEG-like margin -- NOT the ~0.98 km marginal value the
    # pre-#702 buggy final-theta2-anchored Radau re-check reported for this
    # same candidate. Loose bound: pins "small and robust", not an exact float.
    assert corrected_integrator_delta_km < 1e-4, corrected_integrator_delta_km
    assert guard.genuine, guard.notes
