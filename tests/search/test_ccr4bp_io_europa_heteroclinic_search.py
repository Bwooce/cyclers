"""CCR4BP heteroclinic/homoclinic search regression test for Io-Europa (`#695`).

`#695` ran the full pipeline (``scripts/screen_695_ccr4bp_io_europa.py``) as a
genuine discovery attempt on the Io-Europa physical torus. This test does NOT
re-run the expensive coarse KD-tree scan (that lives only in the driver
script, mirroring `#694`'s own screen_694 script, which is likewise untested
by a dedicated pytest module -- its "test" IS the committed result.json).
Instead it regression-locks the driver script's own headline finding by
re-deriving it directly from the SAME seed phase the coarse search found,
using `#694`'s ``refine_candidate``/``ghost_guard`` reused UNMODIFIED:

**The headline finding**: unlike JEG's positive control (`#694`, best
candidate residual_norm ~6.8e-15, an essentially EXACT homoclinic
intersection), Io-Europa's best ghost-guard-passing candidate refines to a
genuine LOCAL MINIMUM with a NONZERO floor (residual_norm ~5.0e-5, a ~16.7 km
position / ~0.54 m/s velocity CLOSEST APPROACH, not an exact closure). Every
near-zero-residual candidate the coarse search turned up FAILED the
ghost-guard instead (either as a trivial near-departure artifact within 1000
km of the torus, or Radau/DOP853-inconsistent -- a chaos-amplified numerical
artifact) -- exactly the false-positive-catching behavior `#620`/`#626`'s
guard was built for. This is reported as an honest NEAR-MISS, not a genuine
connection: the two branches get close but do not intersect.

Physical-unit correction: `#694`'s own ``_L_KM``/``_v_unit_km_s`` are
hardcoded to the JEG system (Europa's SMA); this task's driver script
back-corrects to Io's own SMA-based unit (``core.ccr4bp_io_europa.L_KM``) via
pure post-hoc rescaling (see the driver script's own module docstring for
the full derivation) -- reproduced here identically so this test stays
consistent with the committed ``data/found/695_ccr4bp_io_europa_search/
result.json``.
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy.integrate import solve_ivp

import cyclerfinder.core.ccr4bp_io_europa as ioeu
import cyclerfinder.core.cr3bp as cr3bp
import cyclerfinder.search.ccr4bp_heteroclinic_search as hs
import cyclerfinder.search.variational_ccr4bp_torus as vt
from cyclerfinder.genome.composed_moon_map import resonance_semimajor

# #694's own hardcoded (Europa-SMA-based) physical-unit constants -- needed
# ONLY to back out the nondimensional value from its "_km"-suffixed outputs,
# identically to the driver script (see its own module docstring).
_JEG_L_KM = 671_100.0


def _resonant_symmetric_orbit(
    mu: float, p_sc: int, q_moon: int, *, max_iter: int = 80, tol: float = 1e-12, cap: float = 0.05
) -> tuple[np.ndarray, float, float]:
    """Identical test-only scaffolding to #690/#691/#694's own copies."""
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
    system = ioeu.jupiter_io_europa_default()
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


def test_best_candidate_is_a_robust_near_miss_not_an_exact_connection(
    phys_torus: vt.CCR4BPTorusVariationalResult,
) -> None:
    """Regression-locks #695's own driver-script finding: the seed
    (theta2_u=3.5605, t_u=12.775, theta2_s=3.4558, t_s=24.710, u_lobe=+1,
    s_lobe=-1 -- the exact coarse-search hit from
    ``data/found/695_ccr4bp_io_europa_search/result.json``) refines to a
    genuine LOCAL MINIMUM (``converged=True``) with residual_norm several
    orders of magnitude ABOVE machine precision -- i.e. NOT an exact
    homoclinic closure, unlike JEG's own ~6.8e-15 positive-control result."""
    cand = hs.ManifoldCandidate(
        theta2_u=3.5604716740684323,
        t_u=12.77546133712574,
        theta2_s=3.4557519189487724,
        t_s=24.710431796809,
        gap_planar=0.0004078508862963092,
    )
    refined = hs.refine_candidate(
        phys_torus, phys_torus, cand, lobe_sign_u=1.0, lobe_sign_s=-1.0, n_segments_dir=32
    )
    assert refined is not None
    assert refined.converged is True
    # NOT machine precision (JEG's own best was ~6.8e-15): a genuine nonzero floor.
    assert refined.residual_norm > 1e-6, refined.residual_norm
    assert refined.residual_norm < 1e-3, refined.residual_norm

    guard = hs.ghost_guard(
        phys_torus, phys_torus, refined, lobe_sign_u=1.0, lobe_sign_s=-1.0, n_segments_dir=32
    )
    # #694's own JEG-hardcoded km fields are the WRONG scale for this system
    # (see driver-script docstring) -- back-correct to Io's own SMA unit.
    length_scale = ioeu.L_KM / _JEG_L_KM
    corrected_pos_gap_km = refined.pos_gap_km * length_scale
    corrected_off_torus_km = guard.off_torus_km * length_scale
    corrected_integrator_delta_km = guard.integrator_delta_km * length_scale

    assert corrected_pos_gap_km == pytest.approx(16.689, abs=0.05)
    assert corrected_integrator_delta_km < 1.0  # Radau/DOP853 agree -- not a chaos artifact
    assert corrected_off_torus_km >= 1000.0  # genuinely far from the torus, not trivial
    corrected_genuine = corrected_integrator_delta_km < 1.0 and corrected_off_torus_km >= 1000.0
    assert corrected_genuine is True


def test_near_zero_residual_candidates_fail_the_ghost_guard(
    phys_torus: vt.CCR4BPTorusVariationalResult,
) -> None:
    """The OTHER seeds the coarse search found (which DO refine to
    near-machine-precision residuals) fail the ghost-guard instead -- either
    as trivial near-departure matches (< 1000 km off the torus) or
    Radau/DOP853-inconsistent chaos-amplified artifacts. This is the guard
    doing exactly the job `#620`/`#626` built it for, demonstrated on a
    second system: every "it closed!" candidate here is a false positive."""
    length_scale = ioeu.L_KM / _JEG_L_KM
    # Seed from result.json's (u_lobe=+1, s_lobe=+1) combo, entry 1 -- refines
    # to residual_norm ~6.8e-15 (near-exact) but fails the guard.
    cand = hs.ManifoldCandidate(
        theta2_u=3.4557519189487724,
        t_u=7.39631972149385,
        theta2_s=3.5604716740684323,
        t_s=5.379141615631891,
        gap_planar=0.00032932888632118673,
    )
    refined = hs.refine_candidate(
        phys_torus, phys_torus, cand, lobe_sign_u=1.0, lobe_sign_s=1.0, n_segments_dir=32
    )
    assert refined is not None
    assert refined.residual_norm < 1e-10  # near-exact closure numerically

    guard = hs.ghost_guard(
        phys_torus, phys_torus, refined, lobe_sign_u=1.0, lobe_sign_s=1.0, n_segments_dir=32
    )
    corrected_off_torus_km = guard.off_torus_km * length_scale
    corrected_integrator_delta_km = guard.integrator_delta_km * length_scale
    corrected_genuine = corrected_integrator_delta_km < 1.0 and corrected_off_torus_km >= 1000.0
    # Fails: trivially close to the torus departure point, not a real excursion.
    assert corrected_off_torus_km < 1000.0, corrected_off_torus_km
    assert corrected_genuine is False
