"""Phase 6 Phase 4: method-capability partial order + re-sweep gate (Tasks 4.0a/4.0b)."""

from __future__ import annotations

from cyclerfinder.data.empty_regions import (
    DEFAULT_EMPTY_REGIONS_PATH,
    EmptyRegionReport,
    load_empty_regions_list,
)
from cyclerfinder.data.method_capability import MethodCapability, should_sweep, subsumes

SINGLE = MethodCapability(
    genome="single-ellipse free-return",
    corrector="ballistic_correct",
    capability_tags=frozenset({"ballistic", "patched-conic", "single-arc", "coplanar"}),
    git_sha="aaa",
)
MULTI = MethodCapability(
    genome="two-arc free-return chain",
    corrector="ballistic_correct",
    capability_tags=frozenset({"ballistic", "patched-conic", "multi-arc", "coplanar"}),
    git_sha="bbb",
)
BROKEN = MethodCapability(
    genome="single-ellipse free-return",
    corrector="ballistic_correct",
    capability_tags=frozenset({"ballistic", "patched-conic", "single-arc", "broken-plane"}),
    git_sha="ccc",
)


def test_subsumes_is_reflexive() -> None:
    assert subsumes(SINGLE, SINGLE) is True


def test_multi_arc_subsumes_single_ellipse() -> None:
    assert subsumes(MULTI, SINGLE) is True  # #163 ⊐ #137
    assert subsumes(SINGLE, MULTI) is False  # weaker does not subsume


def test_incomparable_methods_do_not_subsume() -> None:
    # coplanar-multi-arc vs broken-plane-single-arc: neither envelope ⊆ the other
    assert subsumes(MULTI, BROKEN) is False
    assert subsumes(BROKEN, MULTI) is False


# --- Task 4.0b: should_sweep gate -----------------------------------------


def _empty_record(*, region_id: str, method: MethodCapability) -> EmptyRegionReport:
    """A minimal EmptyRegionReport-shaped object carrying region_id + method."""
    return EmptyRegionReport(
        region_id=region_id,
        family="test",
        centre="Jupiter",
        topologies=(),
        method_capability=method,
        search_extent={"points_total": 1},
        prune_gates=("vilm",),
        result={},
        verdict="EMPTY",
        interpretation="",
        source_anchors="",
        run={},
    )


def test_weaker_method_skips_region_a_stronger_method_emptied() -> None:
    registry = [_empty_record(region_id="R", method=MULTI)]
    assert should_sweep(region_id="R", method=SINGLE, registry=registry) is False


def test_stronger_method_re_sweeps_region_a_weaker_method_emptied() -> None:
    registry = [_empty_record(region_id="R", method=SINGLE)]
    assert should_sweep(region_id="R", method=MULTI, registry=registry) is True


def test_incomparable_method_re_sweeps() -> None:
    registry = [_empty_record(region_id="R", method=MULTI)]
    assert should_sweep(region_id="R", method=BROKEN, registry=registry) is True


def test_no_prior_for_region_re_sweeps() -> None:
    assert should_sweep(region_id="UNSWEPT", method=SINGLE, registry=[]) is True


def test_equal_method_skips() -> None:
    registry = [_empty_record(region_id="R", method=SINGLE)]
    assert should_sweep(region_id="R", method=SINGLE, registry=registry) is False


# --- Task 8: leveraging ⊐ single-arc edge ----------------------------------------

_BALLISTIC = MethodCapability(
    genome="single-ellipse free-return (no-leveraging)",
    corrector="ballistic_correct",
    capability_tags=frozenset({"ballistic", "coplanar", "patched-conic", "single-arc"}),
    git_sha="061d42b",
)
_LEVERAGING = MethodCapability(
    genome="phase-full VILM endgame (leveraging)",
    corrector="solve_endgame",
    capability_tags=frozenset({"powered", "coplanar", "patched-conic", "leveraging"}),
    git_sha="deadbee",
)


def test_leveraging_subsumes_ballistic_no_leveraging() -> None:
    assert subsumes(_LEVERAGING, _BALLISTIC)


def test_ballistic_does_not_subsume_leveraging() -> None:
    assert not subsumes(_BALLISTIC, _LEVERAGING)


# --- #521-phase-1 migration proof: the real registry gates a real re-run ---------
#
# data/negative_results.yaml (prose, zero programmatic readers) was migrated into
# data/empty_regions.jsonl (this queryable schema) on 2026-07-02. This is the
# mechanical proof the migration actually works: an equal-or-weaker re-run of the
# #405/#411 SE<->EM patched-CR3BP closure search -- exactly the kind of re-run that
# produced #515-517 this week -- is now SKIPPED by should_sweep() reading the real
# on-disk registry, not a synthetic one.


def test_real_registry_skips_equal_capability_rerun_of_cross_system_se_em_search() -> None:
    registry = load_empty_regions_list(DEFAULT_EMPTY_REGIONS_PATH)
    proposed_rerun = MethodCapability(
        genome="patched-CR3BP SE<->EM connection matcher (planar)",
        corrector="correct_cross_cycle (theta-closure Newton)",
        capability_tags=frozenset(
            {"cr3bp", "patched-cr3bp", "coplanar", "sun-earth-moon", "planar"}
        ),
        git_sha="hypothetical-rerun",
    )
    assert (
        should_sweep(
            region_id="cross-system-se-em-l2-patched-cr3bp-2026-06-20",
            method=proposed_rerun,
            registry=registry,
        )
        is False
    )


def test_real_registry_still_sweeps_an_unrelated_region() -> None:
    registry = load_empty_regions_list(DEFAULT_EMPTY_REGIONS_PATH)
    proposed_rerun = MethodCapability(
        genome="patched-CR3BP SE<->EM connection matcher (planar)",
        corrector="correct_cross_cycle (theta-closure Newton)",
        capability_tags=frozenset(
            {"cr3bp", "patched-cr3bp", "coplanar", "sun-earth-moon", "planar"}
        ),
        git_sha="hypothetical-rerun",
    )
    assert (
        should_sweep(
            region_id="region-that-has-never-been-swept",
            method=proposed_rerun,
            registry=registry,
        )
        is True
    )
