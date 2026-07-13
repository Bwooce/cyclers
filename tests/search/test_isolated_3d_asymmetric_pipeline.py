"""Tests for #582's mandatory downstream pipeline: refine -> classify symmetry
-> populate CandidateSignature -> confirm the literature matcher engages.
"""

from __future__ import annotations

import numpy as np
import pytest

import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.search.er3bp_isolated_seeds import ResonantSeed, all_mmr_seeds
from cyclerfinder.search.isolated_3d_asymmetric_pipeline import (
    build_candidate_signature,
    classify_symmetry,
    literature_anchors_engaged,
    refine_ga_candidate,
)

MU = 0.001


def _system() -> cr3bp.CR3BPSystem:
    return cr3bp.CR3BPSystem(mu=MU, primary="Sun", secondary="planet", l_km=1.0, t_s=1.0)


def _known_32() -> ResonantSeed:
    seeds = all_mmr_seeds(mu=MU)
    seed = next(s for s in seeds if s.p == 3 and s.q == 2)
    assert seed.converged
    return seed


# ---------------------------------------------------------------------------
# refine_ga_candidate
# ---------------------------------------------------------------------------


def test_refine_ga_candidate_recovers_known_seed_from_its_own_ic() -> None:
    """Feeding the corrector the #440 3:2 seed's OWN IC should converge cleanly.

    Not the GA path (that's the driver script's positive-control harness),
    just a regression check that ``refine_ga_candidate``'s
    ``require_monotone_decrease=True`` default converges tightly from a
    near-periodic seed, matching (x0, ydot0, T, C) to a tight tolerance since
    the seed is ALREADY essentially periodic.
    """
    seed = _known_32()
    system = _system()
    genome = np.array(
        [
            seed.state0[0],
            seed.state0[2],
            seed.state0[3],
            seed.state0[4],
            seed.state0[5],
            seed.period,
        ]
    )
    orbit = refine_ga_candidate(system, genome)
    assert orbit.converged
    assert orbit.corrector_residual < 1e-8
    assert orbit.independent_closure_residual < 1e-6
    assert orbit.state0[0] == pytest.approx(seed.state0[0], rel=1e-3)
    assert orbit.state0[4] == pytest.approx(seed.state0[4], rel=1e-3)
    assert pytest.approx(seed.period, rel=1e-3) == orbit.T_TU
    assert orbit.jacobi == pytest.approx(seed.jacobi, abs=1e-3)


def test_refine_ga_candidate_blind_newton_can_diverge_documented_reason() -> None:
    """Regression guard for the exact failure #582's build hit: blind Newton
    (require_monotone_decrease=False) on this under-determined 6-residual/
    7-unknown mode can overshoot a good seed into an unrelated point, which is
    why ``refine_ga_candidate`` defaults True. Not asserting divergence is
    guaranteed in general -- just documenting/pinning a case that demonstrates
    the default, so a future change to the default is a conscious decision,
    not silent regression.

    CI FIX (2026-07-13): the original genome offset (found by hand) converged
    on the monotone path but used 51 of the default 60 max_iter -- a narrow
    margin that flipped to non-convergence on CI's Linux runner (this repo's
    dev machine is a Mac/Accelerate-BLAS build; CI runs Linux) even though
    ``tol=1e-10`` never changed. This genome was instead found by a small
    local sweep over random offsets, selected for a LARGE convergence margin
    (converges in 7 of 60 iterations, residual 8.1e-14 -- roughly 4 orders of
    magnitude inside tol) while still reproducing a clean blind-Newton
    divergence (hits the iteration cap with residual ~1.0, not a near-miss).
    The offsets are arbitrary (this test pins a robustness PROPERTY, not a
    sourced physical value), so picking a more robust arbitrary seed is not a
    tolerance weakening -- see [[feedback_isolated_sweep_flips_suspect_artifact]].
    """
    seed = _known_32()
    system = _system()
    assert seed.state0[0] == pytest.approx(0.7310974, abs=1e-6)  # sanity: same known seed
    genome = np.array([0.745697, -0.0288198, -0.02400734, 0.41758785, 0.02295395, 11.93970051])
    blind = refine_ga_candidate(system, genome, require_monotone_decrease=False)
    monotone = refine_ga_candidate(system, genome, require_monotone_decrease=True)
    assert monotone.converged
    assert monotone.n_iter <= 30  # wide margin vs max_iter=60, not a borderline pass
    # The blind-Newton path is not required to fail forever (it is a
    # documented empirical risk, not a proof), but on THIS seed it does --
    # pin that so a silent behavior change is visible.
    assert not blind.converged


# ---------------------------------------------------------------------------
# classify_symmetry
# ---------------------------------------------------------------------------


def test_classify_symmetry_true_for_known_symmetric_seed() -> None:
    seed = _known_32()
    system = _system()
    result = classify_symmetry(system, seed.state0, seed.period)
    assert result.is_symmetric
    assert result.best_crossing_residual < 1e-9
    assert result.n_crossings_checked >= 1


def test_classify_symmetry_false_when_no_perpendicular_crossing() -> None:
    """A short partial arc broken away from y=0 finds no crossing at all."""
    seed = _known_32()
    system = _system()
    state = np.asarray(seed.state0, dtype=np.float64).copy()
    state[3] += 0.05  # xdot0
    state[5] += 0.05  # zdot0
    result = classify_symmetry(system, state, seed.period / 3.0)
    assert not result.is_symmetric
    assert result.n_crossings_checked == 0
    assert result.best_crossing_residual == float("inf")


def test_classify_symmetry_converged_asymmetric_corrector_output_is_symmetric() -> None:
    """The #582 positive-control finding: the general (unconstrained) corrector,
    fed a GA-realistic asymmetric seed for the 3:2 MMR, converges onto a
    perpendicular-crossing (symmetric) orbit -- exactly the case the
    classifier exists to catch (a converged general-corrector point can sit
    ON a known symmetric orbit; without this check it would be misreported
    novel-asymmetric).
    """
    seed = _known_32()
    system = _system()
    genome = np.array(
        [0.742418, 0.037257, 0.008324, 0.420981, 0.028177, 11.431951]
    )  # #582's own recorded GA-best genome for 3:2 (data/found/582_niching_ga)
    orbit = refine_ga_candidate(system, genome)
    assert orbit.converged
    result = classify_symmetry(system, orbit.state0, orbit.T_TU)
    assert result.is_symmetric
    assert result.best_crossing_residual < 1e-6
    del seed


# ---------------------------------------------------------------------------
# build_candidate_signature / literature_anchors_engaged
# ---------------------------------------------------------------------------


def test_build_candidate_signature_and_matcher_engagement() -> None:
    seed = _known_32()
    system = _system()
    genome = np.array(
        [
            seed.state0[0],
            seed.state0[2],
            seed.state0[3],
            seed.state0[4],
            seed.state0[5],
            seed.period,
        ]
    )
    orbit = refine_ga_candidate(system, genome)
    assert orbit.converged
    sig = build_candidate_signature(system, orbit, p=3, q=2)

    assert sig.primary == "Earth"
    assert sig.sequence == ("Moon",)
    assert sig.resonances == ("3:2",)
    assert sig.topology_3d is not None
    assert set(sig.topology_3d) >= {"k1", "k2", "k_z", "jacobi"}
    assert sig.topology_3d["jacobi"] == pytest.approx(orbit.jacobi)

    anchors = literature_anchors_engaged(sig)
    assert anchors, "the mu=0.001 candidate signature must reach at least one corpus anchor"
    assert any("Antoniadou" in name for name in anchors), (
        f"expected the #579-corrected Antoniadou & Libert 2019 anchor to engage; got {anchors}"
    )


def test_literature_anchors_engaged_empty_for_wrong_primary() -> None:
    """Sanity: a signature that does NOT reuse the Earth/Moon matcher label
    reaches nothing (proving the label is load-bearing, not a no-op)."""
    from cyclerfinder.search.literature_check import CandidateSignature

    sig = CandidateSignature(
        primary="Nonexistent-Primary",
        sequence=("Nonexistent-Body",),
        topology_3d={"k1": 1, "k2": 0, "k_z": 0, "jacobi": 3.0},
    )
    assert literature_anchors_engaged(sig) == []
