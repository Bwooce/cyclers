"""Tests for the CR3BP Jacobi-continuation family-search campaign (Phase 2).

Spec: ``docs/superpowers/specs/2026-06-12-cr3bp-family-continuation-search-design.md``
(campaign testing block), plus the #219 novelty-accounting correction:
  - a seeded family yields >= 1 member OR a clean EMPTY record;
  - a member matching a Ross-published one is labelled REPRODUCTION;
  - a non-reproduction member on a sourced-seed branch is KNOWN-FAMILY-CONTINUATION
    (the SAME published family at a neighbouring C), NEVER NOVEL-SILVER -- so the
    distinct-new-family count stays 0 for an all-sourced-seed campaign;
  - routed review-queue rows == known-family continuations + genuine novels
    (the counts must AGREE -- the 83-novel-vs-0-rows inconsistency regression);
  - the inertial cross-check runs on every kept member;
  - an off-family seed is abandoned (no fabricated members).

SOURCED-GOLDEN DISCIPLINE: the only published number used as an EXPECTED is Ross &
Roberts-Tsoukkas 2025 (AAS 25-621) -- the (3,3) family's C^stable (Table 3) as the
seed's own jacobi (the corrector enforces it). All other assertions are structural
(classification logic, gate plumbing, the off-family guard).

These tests run the REAL continuation + inertial gate on ONE seed (a few seconds);
the full 9-seed campaign is exercised by the script itself, not the suite.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import cast

import pytest

import cyclerfinder.core.cr3bp as cr3bp
import cyclerfinder.search.cr3bp_continuation as cc
import cyclerfinder.search.cr3bp_periodic as cp
import scripts.cr3bp_family_search as fs


def _ross_33_seed() -> fs.FamilySeed:
    return next(s for s in fs.ROSS_SEEDS if s.label == "ross-(3,3)")


@pytest.fixture(scope="module")
def ross_33_result() -> fs.SeedResult:
    # ONE real continuation + inertial-gate run, shared by the tests below (the
    # (3,3) family has a wide stable window, so it produces continuation members).
    return fs.run_seed(_ross_33_seed())


@pytest.mark.slow
def test_ross_seed_yields_members_and_runs_inertial_on_each(
    ross_33_result: fs.SeedResult,
) -> None:
    # A seeded family yields >= 1 member; the seed itself is a REPRODUCTION; every
    # kept member carries an inertial-crosscheck verdict (the gate ran on each).
    res = ross_33_result
    assert res.error == "", res.error
    assert res.n_members >= 1
    assert res.n_reproduction >= 1  # the published midpoint is a reproduction

    classifications = []
    for branch in res.branches:
        members = cast("list[dict[str, object]]", branch["members"])
        for m in members:
            classifications.append(m["classification"])
            # The inertial cross-check ran on EVERY kept member.
            assert m["inertial_verdict"] in {"PASS", "CHECK-FAILED", "INCONCLUSIVE"}
            assert "inertial_delta1_nd" in m
    assert any(c == "REPRODUCTION" for c in classifications)


@pytest.mark.slow
def test_continuation_of_sourced_seed_is_never_novel(ross_33_result: fs.SeedResult) -> None:
    # The #219 correction: every branch here is built by continuing a *published*
    # family representative, so a non-reproduction member is the SAME family at a
    # neighbouring C (KNOWN-FAMILY-CONTINUATION) -- NEVER a NOVEL-SILVER discovery.
    res = ross_33_result
    assert res.n_known_family_continuation >= 1  # (3,3) walks its wide stable window
    assert res.n_novel_silver == 0
    assert res.n_distinct_new_families == 0
    assert res.empty  # EMPTY-for-novelty: no genuinely-distinct new family
    for branch in res.branches:
        members = cast("list[dict[str, object]]", branch["members"])
        for m in members:
            assert m["classification"] != "NOVEL-SILVER"
            assert m["family_label"] == res.label  # provenance: the sourced family


@pytest.mark.slow
def test_seed_member_is_reproduction_of_itself() -> None:
    # The corrected seed (the published Ross point) must classify as REPRODUCTION
    # (within dedup tol of the published member), never as a discovery.
    seed = _ross_33_seed()
    system = cr3bp.CR3BPSystem(
        mu=seed.mu, primary=seed.primary, secondary=seed.secondary, l_km=seed.l_km, t_s=seed.t_s
    )
    corrected = cp.correct_symmetric_fixed_jacobi(
        system,
        seed.x0_seed,
        seed.jacobi,
        seed.period_guess,
        ydot0_sign=seed.ydot0_sign,
        half_crossings=seed.half_crossings,
        tol=1e-10,
    )
    member_state = corrected

    # Build a BranchMember-like for the dedup test via the gauntlet.
    ok, _reason, member = cc._run_gauntlet(
        system,
        member_state,
        period_floor=0.5 * corrected.period,
        period_ceiling=2.0 * corrected.period,
        max_speed_floor=cc.MAX_SPEED_FLOOR_ND,
        amplitude_floor=cc.AMPLITUDE_FLOOR_ND,
        jacobi_tol=cc.JACOBI_CONSERVATION_TOL,
        radau_closure_tol=1e-3,
        radau_jacobi_tol=1e-8,
        rtol=1e-12,
        atol=1e-12,
    )
    assert ok and member is not None
    assert fs._is_reproduction(member, seed)


@pytest.mark.slow
def test_offfamily_seed_is_abandoned_not_fabricated() -> None:
    # The Arenstorf seed uses half_crossings=1, which lands the symmetric corrector
    # on a DIFFERENT orbit (T~5 vs the figure-8's 17.06). The seed-period guard must
    # abandon it (SEED_OFF_FAMILY) and emit ZERO members -- never a fabricated one.
    res = fs.run_seed(fs.ARENSTORF_SEED)
    assert res.empty
    assert res.n_members == 0
    assert "off-family" in res.error


@pytest.mark.slow
def test_routed_members_reach_review_queue_and_counts_agree(
    ross_33_result: fs.SeedResult, tmp_path: Path
) -> None:
    # The 83-novel-vs-0-rows regression: every routed member (known-family
    # continuation + any genuine novel) MUST land in the review queue and the
    # written row count MUST equal the classification counts -- as SILVER, with
    # model provenance and the same-family novelty-scope flag, and nothing
    # touching the catalogue.
    res = ross_33_result
    expected = res.n_known_family_continuation + res.n_novel_silver
    assert expected >= 1  # (3,3) produces continuation members -- the test is not vacuous
    assert len(res.review_members) == expected
    queue = tmp_path / "review_queue.jsonl"
    orig = fs.REVIEW_QUEUE_PATH
    try:
        fs.REVIEW_QUEUE_PATH = queue
        n = fs.write_review_queue([res])
    finally:
        fs.REVIEW_QUEUE_PATH = orig
    assert n == expected
    rows = [json.loads(ln) for ln in queue.read_text().splitlines() if ln.strip()]
    assert len(rows) == n
    for row in rows:
        assert row["verdict_tier"] == "SILVER"
        assert row["model_assumption"] == "cr3bp"
        assert row["classification"] == "KNOWN-FAMILY-CONTINUATION"  # no novels occur here
        assert "novelty_scope" in row
        assert "NOT a discovery" in row["novelty_scope"]
        assert "literature-novelty claim" in row["novelty_scope"]
        # The inertial cross-check is among the recorded gates.
        assert any("inertial-crosscheck" in g for g in row["gates_passed"])


def test_empty_record_written_when_no_novel(tmp_path: Path) -> None:
    # A seed that produces no NOVEL-SILVER member yields a method-versioned EMPTY
    # record carrying the full search extent (no silent truncation).
    seed = _ross_33_seed()
    # Force an EMPTY-for-novelty result by hand: a SeedResult with reproductions only.
    res = fs.SeedResult(
        label=seed.label,
        kind=seed.kind,
        mu=seed.mu,
        primary=seed.primary,
        secondary=seed.secondary,
        seed_jacobi=seed.jacobi,
        seed_period=seed.period_guess,
        n_members=1,
        n_reproduction=1,
        n_novel_silver=0,
        n_stable_novel=0,
        branches=[
            {
                "direction": 1,
                "n_steps_taken": 0,
                "n_rejected": 0,
                "stop_reason": "jacobi_bound",
                "n_members": 1,
            }
        ],
    )
    empty = tmp_path / "empty_regions.jsonl"
    orig = fs.EMPTY_REGIONS_PATH
    try:
        fs.EMPTY_REGIONS_PATH = empty
        n = fs.write_empty_regions([res])
    finally:
        fs.EMPTY_REGIONS_PATH = orig
    assert n == 1
    rec = json.loads(empty.read_text().splitlines()[0])
    assert rec["result"]["novel_silver"] == 0
    assert rec["result"]["distinct_new_families"] == 0
    assert rec["method_capability"]["method"] == fs.METHOD_TAG
    assert rec["search_extent"]["jacobi_range"] == [fs.JACOBI_MIN, fs.JACOBI_MAX]
