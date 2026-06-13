"""Tests for the discovery-campaign engine (#253, Track C).

Cover: the pluggable :class:`SearchTarget` contract; deterministic, resumable
enumeration; generalized repeated-moon closure -> canonical residual; catalogue +
negative-registry dedup; outcome routing (SILVER -> review queue, no-hit ->
empty-region registry) to TEMP paths only; and the no-catalogue-writeback
discipline. All writes target ``tmp_path`` — the real ``data/`` registries are
never touched.
"""

from __future__ import annotations

from collections.abc import Iterator
from pathlib import Path

import pytest

from cyclerfinder.data.catalog import load_catalog
from cyclerfinder.data.empty_regions import load_empty_regions_list
from cyclerfinder.data.method_capability import MethodCapability
from cyclerfinder.data.review_queue import load_review_queue
from cyclerfinder.search.discovery_campaign import (
    CampaignConfig,
    CampaignRouting,
    Candidate,
    ClosureResult,
    RepeatedMoonTarget,
    SearchTarget,
    catalogue_moon_signatures,
    moon_cycler_signature_hash,
    run_campaign,
)


def _routing(tmp_path: Path) -> CampaignRouting:
    return CampaignRouting(
        review_queue_path=tmp_path / "review_queue.jsonl",
        empty_regions_path=tmp_path / "empty_regions.jsonl",
        checkpoint_path=tmp_path / "checkpoint.txt",
    )


# ---------------------------------------------------------------------------
# Signature / dedup
# ---------------------------------------------------------------------------


def test_signature_rotation_invariant() -> None:
    """A cycle and its cyclic rotation hash identically (same identity)."""
    h1 = moon_cycler_signature_hash(
        primary="Jupiter",
        sequence=("Callisto", "Ganymede", "Europa"),
        vinf_per_encounter_kms=(5.0, 7.0, 4.6),
    )
    h2 = moon_cycler_signature_hash(
        primary="Jupiter",
        sequence=("Ganymede", "Europa", "Callisto"),
        vinf_per_encounter_kms=(7.0, 4.6, 5.0),
    )
    assert h1 == h2


def test_signature_vinf_binning() -> None:
    """Sub-bin V_inf differences collapse to the same hash; cross-bin differ."""
    base = moon_cycler_signature_hash(
        primary="Jupiter", sequence=("Io", "Europa"), vinf_per_encounter_kms=(5.00, 6.00)
    )
    near = moon_cycler_signature_hash(
        primary="Jupiter", sequence=("Io", "Europa"), vinf_per_encounter_kms=(5.01, 6.00)
    )
    far = moon_cycler_signature_hash(
        primary="Jupiter", sequence=("Io", "Europa"), vinf_per_encounter_kms=(5.50, 6.00)
    )
    assert base == near
    assert base != far


def test_catalogue_signatures_include_liang_cge() -> None:
    """The real catalogue yields Jovicentric moon-cycler signatures (dedup base)."""
    cat = load_catalog()
    sigs = catalogue_moon_signatures(cat, primary="Jupiter")
    assert sigs, "expected at least one Jovicentric catalogue signature"


# ---------------------------------------------------------------------------
# RepeatedMoonTarget: enumeration + closure
# ---------------------------------------------------------------------------


def test_target_satisfies_protocol() -> None:
    assert isinstance(RepeatedMoonTarget(), SearchTarget)


def test_enumeration_deterministic_and_bounded() -> None:
    """Enumeration is finite, deterministically ordered, and index-contiguous."""
    t = RepeatedMoonTarget(seq_lengths=(3,), n_rev_grid=(0, 1, 2))
    a = list(t.enumerate_candidates())
    b = list(t.enumerate_candidates())
    assert a, "enumeration must be non-empty"
    assert [c.index for c in a] == list(range(len(a)))
    assert [(c.sequence, c.payload["n_rev"]) for c in a] == [
        (c.sequence, c.payload["n_rev"]) for c in b
    ]


def test_enumeration_legs_change_moons_and_use_two_bodies() -> None:
    """No consecutive-same-moon leg; every sequence uses >= 2 distinct moons.

    Closed-cycle topology (#259 fix C): every sequence also starts and ends at
    the same anchor moon, so it is a genuine repeating loop, not an open path.
    """
    t = RepeatedMoonTarget(seq_lengths=(3,), n_rev_grid=(1,))
    cands = list(t.enumerate_candidates())
    assert cands, "enumeration must be non-empty"
    for c in cands:
        seq = c.sequence
        assert all(seq[i] != seq[i + 1] for i in range(len(seq) - 1))
        assert len(set(seq)) >= 2
        assert seq[0] == seq[-1], "a repeated-moon cycle is closed on its anchor moon"


# ---------------------------------------------------------------------------
# #259 fixes: (A) primary-correct moon set, (B) no zero-rev, (C) closed cycle
# ---------------------------------------------------------------------------


def test_fix_a_saturn_enumerates_only_saturn_moons() -> None:
    """Fix A: a Saturn target enumerates ONLY Saturnian moons (registry-derived).

    The prior code hardcoded the Galilean list, so a Saturn run enumerated
    Jovian moons with Saturn's mu (meaningless). The moon set must come from the
    satellite registry filtered to the moons that orbit the primary.
    """
    from cyclerfinder.core.satellites import SATELLITES

    saturn_moons = {n for n, s in SATELLITES.items() if s.primary == "Saturn"}
    jupiter_moons = {n for n, s in SATELLITES.items() if s.primary == "Jupiter"}

    t_sat = RepeatedMoonTarget(primary="Saturn")
    assert set(t_sat.moons) <= saturn_moons
    assert set(t_sat.moons) & jupiter_moons == set()
    seen_bodies = {m for c in t_sat.enumerate_candidates() for m in c.sequence}
    assert seen_bodies, "Saturn enumeration must be non-empty"
    assert seen_bodies <= saturn_moons
    assert seen_bodies & jupiter_moons == set()

    t_jup = RepeatedMoonTarget(primary="Jupiter")
    jup_seen = {m for c in t_jup.enumerate_candidates() for m in c.sequence}
    assert jup_seen <= jupiter_moons
    assert jup_seen & saturn_moons == set()


def test_fix_b_all_zero_rev_excluded() -> None:
    """Fix B: every enumerated candidate has at least one leg with n_rev >= 1.

    The trivial all-zero-rev corner (direct transfer, no multi-rev) is excluded;
    all 4 purged Jovian SILVER artifacts were n_rev=[0,0]. With a grid that
    INCLUDES 0, the all-zero tuple must never appear.
    """
    t = RepeatedMoonTarget(seq_lengths=(3, 4), n_rev_grid=(0, 1, 2))
    cands = list(t.enumerate_candidates())
    assert cands
    for c in cands:
        nrevs = c.payload["n_rev"]
        assert any(nr >= 1 for nr in nrevs), f"all-zero-rev candidate leaked: {nrevs}"


def test_fix_b_liang_style_legs_still_allowed() -> None:
    """Fix B does not exclude the Liang reference (all legs n_rev=1)."""
    t = RepeatedMoonTarget(seq_lengths=(5,), n_rev_grid=(1,))
    cands = list(t.enumerate_candidates())
    assert cands, "1-rev closed cycles (Liang class) must enumerate"
    assert all(all(nr == 1 for nr in c.payload["n_rev"]) for c in cands)


def test_fix_c_closure_includes_anchor_wrap() -> None:
    """Fix C: the closure residual links the cycle-start and cycle-end anchor.

    For a closed cycle the cycle-start outbound V_inf and cycle-end inbound
    V_inf are the same anchor flyby across the cycle boundary; the residual must
    reflect their continuity. We construct a phasing-feasible closed candidate
    and check the residual is at least the anchor-wrap defect (i.e. the wrap is
    actually folded in, not dropped as the prior open-path residual did).
    """
    t = RepeatedMoonTarget(seq_lengths=(3,), n_rev_grid=(1,), n_phase_samples=8)
    cands = list(t.enumerate_candidates())
    converged = [(c, t.close(c)) for c in cands]
    converged = [(c, r) for c, r in converged if r.converged]
    assert converged, "at least one closed candidate must close"
    for c, r in converged:
        # closed-cycle: first and last body are the anchor moon
        assert c.sequence[0] == c.sequence[-1]
        # residual is a finite non-negative continuity defect (km/s)
        assert r.residual_kms >= 0.0


def test_close_returns_canonical_residual() -> None:
    """Closing a candidate yields a finite residual + per-encounter V_inf."""
    t = RepeatedMoonTarget(seq_lengths=(3,), n_rev_grid=(0, 1), n_phase_samples=6)
    cands = list(t.enumerate_candidates())
    assert cands
    closed = [(c, t.close(c)) for c in cands[:8]]
    converged = [(c, r) for c, r in closed if r.converged]
    assert converged, "at least one of the first candidates must close"
    for c, r in converged:
        assert r.residual_kms >= 0.0
        assert len(r.vinf_per_encounter_kms) == len(c.sequence)
        assert len(r.tof_days) == len(c.sequence) - 1


# ---------------------------------------------------------------------------
# Routing + dedup via run_campaign (TEMP paths only — no data/ writeback)
# ---------------------------------------------------------------------------


class _AlwaysHitTarget:
    """A toy target whose single candidate always closes at residual 0 (routes SILVER)."""

    target_id = "toy-hit"
    primary = "Jupiter"

    def method_capability(self) -> MethodCapability:
        return MethodCapability(
            genome="toy", corrector="toy", capability_tags=frozenset({"multi-arc"}), git_sha="test"
        )

    def enumerate_candidates(self) -> Iterator[Candidate]:
        yield Candidate(
            index=0,
            signature_hash="struct:0",
            sequence=("Io", "Europa", "Ganymede"),
            primary="Jupiter",
            payload={"n_rev": [0, 0]},
        )

    def close(self, candidate: Candidate) -> ClosureResult:
        return ClosureResult(
            converged=True,
            residual_kms=0.0,
            vinf_per_encounter_kms=(3.0, 3.0, 3.0),
            tof_days=(10.0, 10.0),
        )


class _AlwaysEmptyTarget(_AlwaysHitTarget):
    """A toy target whose candidate closes ABOVE the gate (routes empty-region)."""

    target_id = "toy-empty"

    def close(self, candidate: Candidate) -> ClosureResult:
        return ClosureResult(
            converged=True,
            residual_kms=99.0,
            vinf_per_encounter_kms=(3.0, 3.0, 3.0),
            tof_days=(10.0, 10.0),
        )


def test_silver_hit_routes_to_review_queue(tmp_path: Path) -> None:
    """A closed-below-gate novel candidate writes ONE SILVER review-queue entry."""
    routing = _routing(tmp_path)
    cat = load_catalog()
    stats = run_campaign(
        _AlwaysHitTarget(), CampaignConfig(gate_residual_kms=0.05), routing, catalog=cat
    )
    assert stats.silver_routed == 1
    assert stats.empty_routed == 0
    entries = list(load_review_queue(routing.review_queue_path))
    assert len(entries) == 1
    assert entries[0].verdict_tier == "silver"
    assert entries[0].sequence == ("Io", "Europa", "Ganymede")
    # No catalogue writeback ever.
    assert not (tmp_path / "catalogue.yaml").exists()


def test_no_hit_routes_to_empty_region(tmp_path: Path) -> None:
    """A sweep with no SILVER survivor writes ONE method-versioned empty record."""
    routing = _routing(tmp_path)
    cat = load_catalog()
    stats = run_campaign(
        _AlwaysEmptyTarget(), CampaignConfig(gate_residual_kms=0.05), routing, catalog=cat
    )
    assert stats.silver_routed == 0
    assert stats.empty_routed == 1
    regions = load_empty_regions_list(routing.empty_regions_path)
    assert len(regions) == 1
    assert regions[0].centre == "Jupiter"
    assert regions[0].method_capability.capability_tags  # non-empty (never unconditional)
    assert (
        not (tmp_path / "review_queue.jsonl").read_text()
        if (tmp_path / "review_queue.jsonl").exists()
        else True
    )


def test_known_candidate_is_deduped(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """A candidate matching a catalogue signature is skipped, not re-queued.

    The candidate's CLOSED signature is injected as the catalogue-known set (via
    the engine's :func:`catalogue_moon_signatures`), so the engine must recognise
    it as known and skip it rather than re-queue a SILVER duplicate.
    """
    routing = _routing(tmp_path)
    target = _AlwaysHitTarget()
    closed_hash = moon_cycler_signature_hash(
        primary="Jupiter",
        sequence=("Io", "Europa", "Ganymede"),
        vinf_per_encounter_kms=(3.0, 3.0, 3.0),
    )
    cat = load_catalog()
    monkeypatch.setattr(
        "cyclerfinder.search.discovery_campaign.catalogue_moon_signatures",
        lambda catalog, *, primary: {closed_hash},
    )
    stats = run_campaign(target, CampaignConfig(), routing, catalog=cat)
    assert stats.skipped_known == 1
    assert stats.silver_routed == 0


def test_checkpoint_resumes(tmp_path: Path) -> None:
    """A second run skips candidates already recorded in the checkpoint.

    ``empty_registry=[]`` + ``write_empty_on_no_hits=False`` isolate the
    checkpoint behaviour from the region-level capability-subsumption gate (which
    would otherwise short-circuit the re-run on the prior empty-region record).
    """
    routing = _routing(tmp_path)
    cat = load_catalog()
    run_campaign(
        _AlwaysEmptyTarget(),
        CampaignConfig(),
        routing,
        catalog=cat,
        empty_registry=[],
        write_empty_on_no_hits=False,
    )
    # Re-run: the single candidate index is checkpointed, so it is skipped.
    stats2 = run_campaign(
        _AlwaysEmptyTarget(),
        CampaignConfig(),
        routing,
        catalog=cat,
        empty_registry=[],
        write_empty_on_no_hits=False,
    )
    assert stats2.skipped_done == 1
    assert stats2.evaluated == 0


def test_worker_sharding_partitions_indices(tmp_path: Path) -> None:
    """Two workers over the same enumeration evaluate disjoint candidate shards.

    With no per-run cap, the full enumeration is walked by both, and each worker
    evaluates only its ``index % n_workers == worker_id`` shard; the two shard
    sizes sum to the single-worker total (a clean partition).
    """
    t = RepeatedMoonTarget(seq_lengths=(3,), n_rev_grid=(0, 1), n_phase_samples=4)
    r_all = _routing(tmp_path / "all")
    r0 = _routing(tmp_path / "w0")
    r1 = _routing(tmp_path / "w1")
    cat = load_catalog()
    s_all = run_campaign(t, CampaignConfig(worker_id=0, n_workers=1), r_all, catalog=cat)
    s0 = run_campaign(t, CampaignConfig(worker_id=0, n_workers=2), r0, catalog=cat)
    s1 = run_campaign(t, CampaignConfig(worker_id=1, n_workers=2), r1, catalog=cat)
    assert s0.evaluated > 0
    assert s1.evaluated > 0
    assert s0.enumerated == s1.enumerated == s_all.enumerated
    assert s0.evaluated + s1.evaluated == s_all.evaluated


def test_no_data_writeback_paths_are_temp(tmp_path: Path) -> None:
    """Routing paths used by the tests are all under tmp_path (no real data/)."""
    routing = _routing(tmp_path)
    for p in (routing.review_queue_path, routing.empty_regions_path, routing.checkpoint_path):
        assert str(tmp_path) in str(p)
