"""Tests for the five-tier accessibility prioritizer (#282).

Smoke-level: feed known input, return sane structure. Does NOT re-test each
underlying scorer's quality (the source modules already have those tests).
What we test here is the COMPOSITION:

* rank-product on heterogeneous numeric ranges produces a sensible ordering
  (the rank-product helper is independently testable);
* patched-conic leg mode runs Tier 0 and documents the architectural gap for
  Tiers 1-5 in the result dict;
* multi-leg aggregation produces sane per-candidate stats from per-leg dicts;
* the leg-reconstruction helper produces N-1 PatchedConicLeg records for an
  N-moon sequence (matching the #264 closure machinery).
"""

from __future__ import annotations

import math

import numpy as np
import pytest

import cyclerfinder.search.five_tier_prioritizer as fp


def test_rank_product_basic() -> None:
    """Two tiers, three candidates -- the best across the stack wins."""
    # candidate A: best in tier 0 (low ΔV), best in tier 1 (low cost) -> rank 1, rank 1
    # candidate B: middle in both
    # candidate C: worst in both
    scores = [
        [1.0, 2.0, 3.0],  # tier 0 (lower-is-better)
        [0.1, 0.2, 0.3],  # tier 1 (lower-is-better)
    ]
    composite = fp.rank_product_score(scores)
    # A < B < C in composite rank product.
    assert composite[0] < composite[1] < composite[2]
    assert math.isclose(composite[0], math.sqrt(1.0 * 1.0), rel_tol=1e-9)
    assert math.isclose(composite[1], math.sqrt(2.0 * 2.0), rel_tol=1e-9)
    assert math.isclose(composite[2], math.sqrt(3.0 * 3.0), rel_tol=1e-9)


def test_rank_product_descending_mixed() -> None:
    """One ascending tier (ΔV), one descending tier (corridor strength)."""
    # tier 0: ΔV in km/s; lower is better. candidate A is the best.
    # tier 1: corridor strength in [0,1]; higher is better. candidate A is best here too.
    scores = [
        [0.5, 1.0, 2.0],  # ascending (lower-is-better)
        [0.9, 0.5, 0.1],  # descending (higher-is-better)
    ]
    composite = fp.rank_product_score(scores, descending=[False, True])
    assert composite[0] < composite[1] < composite[2]


def test_rank_product_nan_sinks_to_worst() -> None:
    """NaN entries get the worst rank within their tier."""
    scores = [
        [1.0, float("nan"), 3.0],  # tier 0: NaN should be ranked worst (3rd)
        [0.1, 0.2, 0.3],
    ]
    composite = fp.rank_product_score(scores)
    # candidate B has NaN in tier 0 (rank 3) and rank 2 in tier 1 -> sqrt(6)
    # candidate A has rank 1 in tier 0 and rank 1 in tier 1 -> 1
    # candidate C has rank 2 in tier 0 (3.0 < nan in lower-is-better) and rank 3
    # in tier 1 -> sqrt(6)
    assert composite[0] < composite[1]
    assert composite[0] < composite[2]


def test_rank_product_input_mismatch_raises() -> None:
    """Mismatched per-tier lengths raise ValueError."""
    with pytest.raises(ValueError, match="same length"):
        fp.rank_product_score([[1.0, 2.0], [0.1, 0.2, 0.3]])


def test_score_leg_runs_tier0_and_documents_skipped_tiers() -> None:
    """One Lambert leg in: Tier 0 verdict + tiers_skipped audit out."""
    prioritizer = fp.FiveTierPrioritizer()
    # Synthetic heliocentric-scale leg (Earth-to-Mars ballpark). The NN's
    # admit_threshold defaults to 5 km/s; we are not asserting the value,
    # only the STRUCTURE.
    leg = fp.PatchedConicLeg(
        label_from="A",
        label_to="B",
        r1_m=np.array([1.5e11, 0.0, 0.0]),
        v1_m_s=np.array([0.0, 30000.0, 0.0]),
        r2_m=np.array([2.3e11, 0.0, 0.0]),
        v2_m_s=np.array([0.0, 24000.0, 0.0]),
        dt_s=240.0 * 86400.0,
        mu_m3_s2=1.32712440018e20,
    )
    result = prioritizer.score_leg(leg)
    assert "tier0_predicted_dv_kms" in result
    assert "tier0_admitted" in result
    assert "tier0_model_available" in result
    assert result["label_from"] == "A"
    assert result["label_to"] == "B"
    # The architectural gap MUST be documented.
    assert "tiers_skipped" in result
    skipped_ids = {s["tier"] for s in result["tiers_skipped"]}
    assert skipped_ids == {1, 2, 3, 4, 5}
    for s in result["tiers_skipped"]:
        assert "CR3BP" in s["reason"]


def test_score_candidate_legs_aggregates_per_leg_stats() -> None:
    """Multi-leg candidate aggregates max_dv, sum_dv, all_admitted flags."""
    prioritizer = fp.FiveTierPrioritizer()
    leg1 = fp.PatchedConicLeg(
        label_from="A",
        label_to="B",
        r1_m=np.array([1.5e11, 0.0, 0.0]),
        v1_m_s=np.array([0.0, 30000.0, 0.0]),
        r2_m=np.array([2.3e11, 0.0, 0.0]),
        v2_m_s=np.array([0.0, 24000.0, 0.0]),
        dt_s=240.0 * 86400.0,
        mu_m3_s2=1.32712440018e20,
    )
    leg2 = fp.PatchedConicLeg(
        label_from="B",
        label_to="C",
        r1_m=np.array([2.3e11, 0.0, 0.0]),
        v1_m_s=np.array([0.0, 24000.0, 0.0]),
        r2_m=np.array([1.5e11, 0.0, 0.0]),
        v2_m_s=np.array([0.0, 30000.0, 0.0]),
        dt_s=240.0 * 86400.0,
        mu_m3_s2=1.32712440018e20,
    )
    stats = prioritizer.score_candidate_legs([leg1, leg2])
    assert stats["n_legs"] == 2
    assert "tier0_max_dv_kms" in stats
    assert "tier0_sum_dv_kms" in stats
    assert "tier0_all_admitted" in stats
    # If both legs ran cleanly, max <= sum (and both should be finite for a
    # heliocentric short-period leg).
    if math.isfinite(stats["tier0_max_dv_kms"]) and math.isfinite(stats["tier0_sum_dv_kms"]):
        assert stats["tier0_max_dv_kms"] <= stats["tier0_sum_dv_kms"]


def test_legs_from_repeated_moon_candidate_pluto() -> None:
    """Reconstruct legs for a Pluto Hydra-Nix-Hydra 2-leg candidate.

    Matches the structure of repeated-moon-pluto-00000045 (review_queue line 1).
    The point of the test is to confirm the leg-reconstruction helper RUNS and
    produces N-1 PatchedConicLeg records; we do NOT assert on a specific phasing
    because the #264 close sweeps a 24-phase grid and picks the best by V_inf
    continuity, while our reconstruction picks the best by interior V_inf
    matching (a slightly different score). Both produce LEGS, which is the
    test's invariant.
    """
    legs = fp.legs_from_repeated_moon_candidate(
        primary="Pluto",
        sequence=["Hydra", "Nix", "Hydra"],
        n_rev=[0, 1],
        phase_samples=4,  # tiny grid for the smoke
        tof_resonance_grid=(1.0,),
    )
    assert legs is not None
    assert len(legs) == 2  # N-1 legs for an N=3-moon sequence
    for leg in legs:
        assert leg.label_from in {"Hydra", "Nix", "Charon"}
        assert leg.label_to in {"Hydra", "Nix", "Charon"}
        # Positions and velocities are SI.
        assert leg.r1_m.shape == (3,)
        assert leg.v1_m_s.shape == (3,)
        assert leg.dt_s > 0.0
        # Pluto mu in SI: 975.5 km^3/s^2 = 975.5 * (10^3)^3 m^3/s^2 = 9.755e11 m^3/s^2
        assert math.isclose(leg.mu_m3_s2, 9.755e11, rel_tol=1e-3)
        # Pluto-system Lambert: dt is in days for these small moons, dt < a few hundred days.
        assert leg.dt_s < 365.0 * 86400.0


def test_legs_from_repeated_moon_candidate_mismatched_nrev_raises() -> None:
    """n_rev length mismatch with n_legs raises ValueError."""
    with pytest.raises(ValueError, match="n_rev length"):
        fp.legs_from_repeated_moon_candidate(
            primary="Pluto",
            sequence=["Hydra", "Nix", "Hydra"],
            n_rev=[0],  # should be 2 for a 3-moon sequence
            phase_samples=2,
            tof_resonance_grid=(1.0,),
        )
