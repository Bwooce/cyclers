"""Tests for the two-tier accessibility prioritizer (#263).

The prioritizer composes the Braik-Ross energy-PRESERVING heading-fan reachable
scorer (#236/#247, ``reachable_network``) with the Zhou-Armellin energy-CHANGING
single-impulse reachable set (#239, ``reachable_impulsive``) into one API. The
two tiers are NOT redundant; each fills the other's blind spot. The headline
test (`test_tier2_lights_up_energy_changing_edge`) demonstrates that with the
tier-1 atlas built at a Hill-region-incompatible Jacobi constant the heading-fan
sees no overlap at all, while the tier-2 impulsive scorer (which is not bound to
a single C_J manifold) still confirms accessibility within a small ΔV budget --
the energy-changing complement value.

All EXPECTED values trace to: published source (Braik-Ross R21-S hard-access at
C_J=3.1294, Table 4 / Fig. 10), the existing energy-changing invariant of the
Zhou-Armellin impulse model (a 50 m/s impulse measurably changes C_J,
``test_apply_impulse_changes_jacobi_constant`` in ``test_reachable_impulsive``),
or pure composition mechanics (tier-1 threshold gating, identical-rep degeneracy).
The tests use the existing offline-recovered Braik-Ross subset fixture so no new
goldens are introduced.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

import cyclerfinder.search.reachable_network as rn
import cyclerfinder.search.reachable_representatives as rr
import cyclerfinder.search.two_tier_prioritizer as ttp

# ---------------------------------------------------------------------------
# Shared fixtures: a fast tier-1 grid and a fast prioritizer with the existing
# offline-recovered Braik-Ross subset. Recovery is ~1 s; tier-1 atlases at the
# coarse grid + short horizon are ~0.3 s each per rep; the suite stays under
# 90 s on a laptop. Settings are deliberately COARSE for speed -- this suite
# tests COMPOSITION, not solver quality (the gate suite already covers that).
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def system() -> object:
    return rr.braik_ross_system()


@pytest.fixture(scope="module")
def reps_subset() -> dict[str, rr.Representative]:
    """Subset of confirmed Braik-Ross offline-recovered representatives.

    Uses LL1, LL2, DPO, R21-S, C21 -- the smallest fast subset that includes
    one stable-resonant hard-access node (R21-S, Braik-Ross Table 4 bottom
    rank) and a cycler at a slightly different jacobi (C21, AAS-25-621
    ROSS_C21_JACOBI = 3.129389531).
    """
    sysm = rr.braik_ross_system()
    return {r.label: r for r in rr.recover_offline_set(sysm) if r.confirmed}


def _fast_prioritizer(system: object, *, tier1_min_overlap: float = 0.0) -> ttp.TwoTierPrioritizer:
    """Coarse-grid prioritizer tuned for fast tests (NOT for solver-quality)."""
    grid = rn.VoxelGrid(dx=0.05, dy=0.05, dtheta=math.radians(20.0))
    return ttp.TwoTierPrioritizer(
        system=system,  # type: ignore[arg-type]
        grid=grid,
        tier1_min_overlap=tier1_min_overlap,
        tier1_n_seeds=4,
        tier1_n_fan=3,
        tier1_delta_max=math.radians(30.0),
        tier1_horizon_frac=0.3,
        tier2_n_seeds=2,
        tier2_n_alpha=5,
        tier2_n_beta=9,
        tier2_n_dv_steps=2,
        tier2_horizon_frac=0.25,
        tier2_proximity_tol=0.08,  # ~30 800 km; coarse, matches the coarse grid
        tier2_n_target_samples=16,
    )


# ---------------------------------------------------------------------------
# (1) Reproduce-before-trust: tier-1 reproduces the Braik-Ross dominance
# direction already verified by the gate suite (the R21-S hard-access reading)
# on the recovered subset. We do not re-test the full Table-4 ranking (the gate
# already covers it); we only verify that wrapping tier 1 in the prioritizer
# does not alter the qualitative published structural result.
# ---------------------------------------------------------------------------


def test_tier1_reproduces_r21s_hard_access(
    system: object, reps_subset: dict[str, rr.Representative]
) -> None:
    """REPRODUCE-BEFORE-TRUST: tier-1 confirms R21-S as a hard-access node.

    Braik-Ross Table 4 / Fig. 10: the 2:1 stable resonant R21-S is the persistent
    hard-access family. On the recovered subset, the prioritizer's tier-1 mean
    proxy ΔV from R21-S to other nodes must be GREATER than the mean from LL1
    (an L1 Lyapunov, a Braik-Ross hub) -- i.e. R21-S sits at the bottom of
    accessibility from itself outward, consistent with the published reading.
    """
    prio = _fast_prioritizer(system)
    by = reps_subset
    others = ["LL1", "LL2", "DPO", "C21"]

    def _mean_tier1_kms(source: str) -> float:
        vals: list[float] = []
        for dest in others:
            if dest == source:
                continue
            out = prio.score_pair(
                by[source],
                by[dest],
                c_j=rr.C_J_BRAIK_ROSS,
                dv_budget_kms=1.5,
            )
            v = float(out["tier1_proxy_dv_kms"])  # type: ignore[arg-type]
            if math.isfinite(v):
                vals.append(v)
        if not vals:
            return math.inf
        return sum(vals) / len(vals)

    mean_from_ll1 = _mean_tier1_kms("LL1")
    mean_from_r21s = _mean_tier1_kms("R21-S")
    print(
        f"\nTier-1 mean proxy ΔV (km/s): LL1 -> others = {mean_from_ll1:.3f}, "
        f"R21-S -> others = {mean_from_r21s:.3f}"
    )
    # R21-S is the harder-to-leave node (Braik-Ross Table 4 R21-S bottom rank
    # under both strength and closeness). On the coarse subset the inequality is
    # qualitative: R21-S must NOT be the cheapest source.
    assert mean_from_r21s >= mean_from_ll1, (
        f"R21-S (hard-access) should not be cheaper than LL1 (hub): "
        f"R21-S mean={mean_from_r21s:.3f}, LL1 mean={mean_from_ll1:.3f}"
    )


# ---------------------------------------------------------------------------
# (2) HEADLINE: tier 2 lights up an energy-changing edge that tier 1 misses.
#
# Construction: evaluate the tier-1 atlas at a Hill-region-incompatible Jacobi
# constant (C_J = 3.7, tighter than any of our recovered reps' own energies).
# At this C_J the planar Hill region excludes both reps' state0 positions, so
# reachable_speed errors and the heading-fan atlas comes up empty -- tier 1
# returns proxy_dv = inf, zero overlap, zero score. Tier 2 is NOT bound to a
# C_J manifold (an impulse CHANGES C_J by construction; see
# ``test_apply_impulse_changes_jacobi_constant`` in
# ``tests/search/test_reachable_impulsive.py``), so it still reaches the
# target within a small ΔV budget. This IS the energy-changing complement.
# ---------------------------------------------------------------------------


def test_tier2_lights_up_energy_changing_edge(
    system: object, reps_subset: dict[str, rr.Representative]
) -> None:
    """Energy-changing edge: tier 2 reaches a pair tier 1 cannot bridge.

    With tier-1 evaluated at an out-of-Hill-region C_J, the heading-fan atlas
    returns no overlap (tier-1 ``accessible=False`` / score 0.0). Tier-2's
    single-impulse footprint CHANGES the Jacobi constant by construction and
    still reaches the target within the 1 km/s budget. The pair becomes
    ``accessible = True`` via ``dominant_tier = "tier2"`` -- the headline
    energy-changing complement.
    """
    prio = _fast_prioritizer(system)
    by = reps_subset
    # C_J = 3.7 is well above every Braik-Ross common-energy regime; the planar
    # Hill region collapses around the primaries and excludes our reps' state0.
    out = prio.score_pair(by["LL1"], by["DPO"], c_j=3.7, dv_budget_kms=1.0)
    # Tier 1: NO overlap at this incompatible C_J.
    tier1_score = out["tier1_heading_overlap"]
    assert tier1_score == 0.0, (
        f"Tier 1 should see no overlap at incompatible C_J=3.7; got score={tier1_score}"
    )
    assert math.isinf(float(out["tier1_proxy_dv_kms"])), (  # type: ignore[arg-type]
        f"Tier 1 proxy ΔV should be inf when atlas is empty; got {out['tier1_proxy_dv_kms']}"
    )
    assert out["tier1_admitted"] is False
    # Tier 2: lights up within budget -- the energy-changing complement.
    assert math.isfinite(float(out["tier2_impulsive_min_dv_kms"])), (  # type: ignore[arg-type]
        f"Tier 2 should find an impulsive bridge; got {out['tier2_impulsive_min_dv_kms']}"
    )
    assert out["tier2_within_budget"] is True
    # The pair is accessible ONLY because tier 2 caught what tier 1 missed.
    assert out["accessible"] is True
    assert out["dominant_tier"] == "tier2"


# ---------------------------------------------------------------------------
# (3) Composition mechanics: tier-1 threshold gating.
# ---------------------------------------------------------------------------


def test_tier1_threshold_gates_admission(
    system: object, reps_subset: dict[str, rr.Representative]
) -> None:
    """A high ``tier1_min_overlap`` rejects tier 1 admission; zero admits all.

    With ``tier1_min_overlap = 1e9`` no overlap score can ever exceed it (scores
    are in (0, 1]); ``tier1_admitted`` is uniformly False even when the pair
    DOES overlap. With ``tier1_min_overlap = 0.0`` any nonzero overlap admits.
    Tier 2 is unaffected by the threshold and continues to report its own
    verdict, so a STRICT-threshold pair can still be ``accessible`` via tier 2.
    """
    by = reps_subset
    prio_strict = _fast_prioritizer(system, tier1_min_overlap=1e9)
    prio_liberal = _fast_prioritizer(system, tier1_min_overlap=0.0)
    out_strict = prio_strict.score_pair(
        by["LL1"], by["LL2"], c_j=rr.C_J_BRAIK_ROSS, dv_budget_kms=1.0
    )
    out_liberal = prio_liberal.score_pair(
        by["LL1"], by["LL2"], c_j=rr.C_J_BRAIK_ROSS, dv_budget_kms=1.0
    )
    # Both produce the SAME underlying overlap score -- the threshold gates
    # admission, not the score itself.
    assert out_strict["tier1_heading_overlap"] == out_liberal["tier1_heading_overlap"]
    assert out_strict["tier1_admitted"] is False
    assert out_liberal["tier1_admitted"] is True
    # Tier-2 is unaffected by tier-1 threshold; results match.
    assert out_strict["tier2_impulsive_min_dv_kms"] == out_liberal["tier2_impulsive_min_dv_kms"]


# ---------------------------------------------------------------------------
# (4) API non-regression: identical reps + zero budget + tiny budget all return
# sensible (non-None) dicts.
# ---------------------------------------------------------------------------


def test_identical_reps_trivially_accessible(
    system: object, reps_subset: dict[str, rr.Representative]
) -> None:
    """Identical representatives short-circuit to a zero-cost accessible pair."""
    prio = _fast_prioritizer(system)
    by = reps_subset
    out = prio.score_pair(by["LL1"], by["LL1"], c_j=rr.C_J_BRAIK_ROSS, dv_budget_kms=1.0)
    assert isinstance(out, dict)
    assert out["accessible"] is True
    assert out["tier1_heading_overlap"] == 1.0
    assert out["tier1_proxy_dv_kms"] == 0.0
    assert out["tier2_impulsive_min_dv_kms"] == 0.0
    assert out["dominant_tier"] == "both"


def test_zero_budget_returns_sensible_dict(
    system: object, reps_subset: dict[str, rr.Representative]
) -> None:
    """A zero ΔV budget produces no tier-2 hit but does NOT crash.

    Tier 1 is independent of budget and still produces its overlap score; tier 2
    has no rungs to evaluate (``inf``), and the pair is accessible iff tier 1's
    proxy alone is within the (zero) budget -- which for a strictly-positive
    proxy is False. The dict has every documented field.
    """
    prio = _fast_prioritizer(system)
    by = reps_subset
    out = prio.score_pair(by["LL1"], by["LL2"], c_j=rr.C_J_BRAIK_ROSS, dv_budget_kms=0.0)
    assert out is not None
    for key in (
        "tier1_heading_overlap",
        "tier1_proxy_dv_kms",
        "tier1_admitted",
        "tier2_impulsive_min_dv_kms",
        "tier2_within_budget",
        "accessible",
        "dominant_tier",
    ):
        assert key in out, f"missing key {key}"
    assert math.isinf(float(out["tier2_impulsive_min_dv_kms"]))  # type: ignore[arg-type]
    # At zero budget nothing within budget: neither tier qualifies.
    assert out["tier2_within_budget"] is False
    assert out["dominant_tier"] in {"neither", "tier1"}


def test_rank_destinations_returns_sorted_list(
    system: object, reps_subset: dict[str, rr.Representative]
) -> None:
    """``rank_destinations`` returns a non-empty list sorted by cheapest bridge."""
    prio = _fast_prioritizer(system)
    by = reps_subset
    candidates = [by["LL2"], by["DPO"], by["R21-S"]]
    ranked = prio.rank_destinations(by["LL1"], candidates, c_j=rr.C_J_BRAIK_ROSS, dv_budget_kms=1.0)
    assert len(ranked) == 3

    def _key(d: dict[str, object]) -> float:
        return min(
            float(d["tier1_proxy_dv_kms"]),  # type: ignore[arg-type]
            float(d["tier2_impulsive_min_dv_kms"]),  # type: ignore[arg-type]
        )

    keys = [_key(d) for d in ranked]
    assert keys == sorted(keys), f"rank_destinations not sorted: {keys}"


# ---------------------------------------------------------------------------
# (5) Unit conversions are symmetric.
# ---------------------------------------------------------------------------


def test_kms_nd_round_trip() -> None:
    """``kms_to_nd`` and ``nd_to_kms`` are inverses (VU_MS ~ 1023.16 m/s)."""
    for dv in (0.0, 0.001, 0.1, 1.0, 10.0):
        assert ttp.nd_to_kms(ttp.kms_to_nd(dv)) == pytest.approx(dv, rel=1e-12)
    # The EM VU is ~1.023 km/s; 1 km/s should be ~0.977 nondimensional.
    assert ttp.kms_to_nd(1.0) == pytest.approx(1000.0 / rn.VU_MS, rel=1e-12)
    assert 0.95 < ttp.kms_to_nd(1.0) < 1.0


def test_repview_normalization_from_duck_typed_object() -> None:
    """A duck-typed object with ``state0`` + ``period`` is accepted by ``score_pair``.

    The prioritizer should not require a strict ``Representative`` instance --
    discovery campaigns may produce raw (state0, period) tuples wrapped in
    lightweight scan-output objects. Pure API mechanics; no propagation.
    """

    class Duck:
        label = "duck"
        state0 = np.array([0.85, 0.0, 0.0, 0.0, 0.55, 0.0])
        period = 6.0

    rv = ttp._as_repview(Duck(), "fallback")
    assert rv.label == "duck"
    assert rv.period == 6.0
    assert rv.state0.shape == (6,)
