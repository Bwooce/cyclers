"""Phase 4 — the novelty discovery loop (Forge plan §29-31).

These tests exercise :mod:`cyclerfinder.data.discover_novel`:

* :func:`cycler_from_closure` — the bridge from a real-ephemeris
  :class:`~cyclerfinder.search.correct.BallisticClosureResult` to a full
  :class:`~cyclerfinder.model.cycler.Cycler` (so the signature + Axis-A code-path
  machinery can consume a multi-arc closure).
* :func:`discover_novel` — construction-first discovery over a topology x scan
  grid, routed by :func:`~cyclerfinder.verify.gauntlet.run_gauntlet`.

GOLDEN DISCIPLINE: no EXPECTED side here is a value our own optimiser produced.
The shape/structure assertions (verdict tier routing, supersession surfacing,
match outcome) are logic predicates, not physics goldens. The one physics-touching
test (the Sanchez-regime closure → SILVER) asserts only *regime/feasibility/tier*,
never a sourced V_inf — mirroring tests/search/test_correct_sanchez_regime.py.
"""

from __future__ import annotations

from datetime import UTC, datetime, timedelta

import numpy as np
import pytest

from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.search.correct import ballistic_correct
from cyclerfinder.verify.gauntlet import VerdictTier

_J2000 = datetime(2000, 1, 1, 12, 0, 0, tzinfo=UTC)
_PERIOD_DAYS = (1.4612 + 2.8096) * 365.25  # 2-synodic Russell arcs (E-M-E-E)


# ---------------------------------------------------------------------------
# cycler_from_closure — the multi-arc bridge
# ---------------------------------------------------------------------------


@pytest.mark.slow
def test_cycler_from_closure_reconstructs_geometry() -> None:
    """A converged ballistic closure round-trips into a Cycler whose
    per-encounter V_inf magnitudes match the corrector's own vector magnitudes.

    EXPECTED side is the corrector's OWN output (an internal-consistency
    round-trip, NOT a sourced golden): the bridge must not distort the geometry.
    """
    from cyclerfinder.data.discover_novel import cycler_from_closure

    ephem = Ephemeris("astropy")
    t0 = ((datetime(2030, 3, 22, tzinfo=UTC) + timedelta(days=-20)) - _J2000).total_seconds()
    seq = ("E", "M", "E", "E")
    revs = (0, 0, 1)
    branch = ("single", "single", "low")
    closure = ballistic_correct(
        sequence=seq,
        per_leg_revs=revs,
        per_leg_branch=branch,
        t0_seed_sec=t0,
        tof_seed_days=(154.0, 379.0),
        period_sec=_PERIOD_DAYS * 86400.0,
        ephem=ephem,
        vinf_cap=14.0,
        slack_leg=2,
    )
    assert closure.converged

    cycler = cycler_from_closure(closure, seq, revs, branch, ephem)
    assert list(cycler.bodies) == list(seq)
    assert len(cycler.encounters) == len(seq)
    assert len(cycler.legs) == len(seq) - 1

    # Per-intermediate-encounter V_inf magnitude agrees with the corrector's
    # own per-encounter magnitude (the bridge re-solves the same Lambert legs).
    got = [
        max(float(np.linalg.norm(e.vinf_in)), float(np.linalg.norm(e.vinf_out)))
        for e in cycler.encounters
    ]
    # Corrector reports a per-encounter magnitude tuple; compare the maxima.
    assert max(got) == pytest.approx(max(closure.vinf_per_encounter_kms), abs=0.05)


# ---------------------------------------------------------------------------
# discover_novel — verdict routing (no live physics; synthetic closures)
# ---------------------------------------------------------------------------


def test_novel_unsourced_machine_confirmed_is_silver() -> None:
    """A closed, code-path-agreeing candidate whose signature matches NO
    catalogue row routes to SILVER (machine-confirmed, unsourced) — never GOLD.

    Uses the in-process classifier directly with synthesised inputs so the test
    is fast and asserts only the routing logic (golden-clean predicate test).
    """
    from cyclerfinder.data.discover_novel import classify_candidate_verdict

    verdict = classify_candidate_verdict(
        candidate_id="synthetic-novel-1",
        agreed=True,
        n_paths_available=2,
        n_paths_passed=2,
        match_outcome="novel",
        known_id=None,
        superseded_by=(),
        falsified=False,
    )
    assert verdict.tier is VerdictTier.SILVER
    assert verdict.provenance["match_status"] == "unmatched"


def test_novel_match_against_superseded_row_surfaces_chain() -> None:
    """R1 delta 3: a candidate matching a row that carries ``superseded_by`` must
    surface the supersession chain and never report a clean 'known'."""
    from cyclerfinder.data.discover_novel import classify_candidate_verdict

    verdict = classify_candidate_verdict(
        candidate_id="synthetic-superseded",
        agreed=True,
        n_paths_available=2,
        n_paths_passed=2,
        match_outcome="known",
        known_id="vem-emeeve-3syn",
        superseded_by=("vem-emeeve-3syn-realized",),
        falsified=False,
    )
    assert verdict.provenance["match_status"] == "superseded"
    assert verdict.provenance["superseded_by"] == ["vem-emeeve-3syn-realized"]


def test_bend_infeasible_closure_is_rejected_not_silver() -> None:
    """A converged-but-bend-INFEASIBLE closure is physically inadmissible and
    must route to REJECTED (a flyby that cannot deliver the required turn is an
    impossible bend), NEVER SILVER. evaluate_closure auto-falsifies it.
    """
    from cyclerfinder.data.discover_novel import evaluate_closure
    from cyclerfinder.search.correct import BallisticClosureResult

    ephem = Ephemeris("astropy")
    t0 = ((datetime(2030, 3, 22, tzinfo=UTC) + timedelta(days=-20)) - _J2000).total_seconds()
    seq = ("E", "M", "E", "E")
    revs = (0, 0, 1)
    branch = ("single", "single", "low")
    # A geometry that closes the residual but is NOT bend-feasible.
    real = ballistic_correct(
        sequence=seq,
        per_leg_revs=revs,
        per_leg_branch=branch,
        t0_seed_sec=t0,
        tof_seed_days=(154.0, 379.0),
        period_sec=_PERIOD_DAYS * 86400.0,
        ephem=ephem,
        vinf_cap=14.0,
        slack_leg=2,
    )
    infeasible = BallisticClosureResult(
        t0_sec=real.t0_sec,
        tof_days=real.tof_days,
        max_residual_kms=real.max_residual_kms,
        vinf_per_encounter_kms=real.vinf_per_encounter_kms,
        converged=True,
        bend_feasible=False,  # the inadmissibility
        vinf_cap_ok=True,
    )
    finding = evaluate_closure(
        infeasible,
        sequence=seq,
        per_leg_revs=revs,
        per_leg_branch=branch,
        period_k=2,
        ephem=ephem,
        slack_leg=2,
    )
    assert finding.verdict.tier is VerdictTier.REJECTED


def test_falsified_candidate_is_rejected() -> None:
    """A falsified candidate is REJECTED even if all other axes look clean."""
    from cyclerfinder.data.discover_novel import classify_candidate_verdict

    verdict = classify_candidate_verdict(
        candidate_id="synthetic-bogus",
        agreed=True,
        n_paths_available=2,
        n_paths_passed=2,
        match_outcome="novel",
        known_id=None,
        superseded_by=(),
        falsified=True,
    )
    assert verdict.tier is VerdictTier.REJECTED


def test_multi_arc_drops_inapplicable_construction_path() -> None:
    """For a multi-arc candidate the resonance-construction path (b) is treated
    as unavailable (cross-fidelity / structurally inapplicable), so a clean
    lamberthub + Kepler-reprop pair still yields 'agreed' on two paths.

    This is the predicate the e2e routing depends on: a multi-arc closure whose
    (a) and (c) paths pass must NOT be vetoed by a spurious path-(b) failure.
    """
    from cyclerfinder.data.discover_novel import _agreement_predicate
    from cyclerfinder.verify.agreement import (
        AgreementReport,
        ConstructionOptimiserPathResult,
        KeplerRepropPathResult,
        LamberthubPathResult,
    )

    lamberthub = LamberthubPathResult(available=True, per_leg=(), max_diff_mps=0.0, passed=True)
    # Path (b) RAN and FAILED — the spurious single-ellipse-vs-multi-arc veto.
    constr = ConstructionOptimiserPathResult(
        available=True,
        resonant_vinf_kms={"E": 5.0},
        cycler_vinf_kms={"E": 9.75},
        construction_max_diff_kms=1.76,
        optimiser_available=False,
        optimiser_vinf_kms={},
        optimiser_max_diff_kms=None,
        max_diff_kms=1.76,
        passed=False,
    )
    kepler = KeplerRepropPathResult(
        available=True, per_leg_residual_km=(0.001,), max_residual_km=0.001, passed=True
    )
    report = AgreementReport(
        lamberthub=lamberthub,
        construction_optimiser=constr,
        kepler_reprop=kepler,
        n_paths_available=3,
        n_paths_passed=2,
        agreed=False,  # the raw report vetoes on path (b)
    )
    # Single-ellipse: verbatim (still vetoed).
    assert _agreement_predicate(report, multi_arc=False) == (False, 3, 2)
    # Multi-arc: path (b) dropped -> (a)+(c) agree on two paths.
    assert _agreement_predicate(report, multi_arc=True) == (True, 2, 2)


# ---------------------------------------------------------------------------
# Sanchez-regime end-to-end: a real closure routes to SILVER
# ---------------------------------------------------------------------------


@pytest.mark.slow
def test_sanchez_regime_closure_routes_to_silver_e2e() -> None:
    """End-to-end: the demonstrated DE440 E-M-E-E near-ballistic closure, run
    through the novelty loop's per-candidate evaluator, yields a finding that is
    machine-confirmed (Axis A agrees) and unsourced -> SILVER (or BRONZE if a
    code path is unavailable). It must NOT be GOLD (no independent source) and
    must NOT be REJECTED (it closes feasibly).

    NON-GOLDEN: asserts only the tier band + that a real closure flows through
    the full bridge -> signature -> agreement -> gauntlet pipeline. No sourced
    V_inf is asserted (project memory feedback_golden_tests_sourced_only).
    """
    from cyclerfinder.data.discover_novel import evaluate_closure

    ephem = Ephemeris("astropy")
    t0 = ((datetime(2030, 3, 22, tzinfo=UTC) + timedelta(days=-20)) - _J2000).total_seconds()
    seq = ("E", "M", "E", "E")
    revs = (0, 0, 1)
    branch = ("single", "single", "low")
    closure = ballistic_correct(
        sequence=seq,
        per_leg_revs=revs,
        per_leg_branch=branch,
        t0_seed_sec=t0,
        tof_seed_days=(154.0, 379.0),
        period_sec=_PERIOD_DAYS * 86400.0,
        ephem=ephem,
        vinf_cap=14.0,
        slack_leg=2,
    )
    assert closure.converged and closure.bend_feasible

    finding = evaluate_closure(
        closure,
        sequence=seq,
        per_leg_revs=revs,
        per_leg_branch=branch,
        period_k=2,
        ephem=ephem,
    )
    # Machine-confirmed-or-thin, but never GOLD (no source) and never REJECTED
    # (it closes feasibly). Asserted via the tier set to keep both bounds explicit.
    assert finding.verdict.tier in {VerdictTier.SILVER, VerdictTier.BRONZE}
    # The signature is computed and the match outcome is recorded.
    assert finding.signature is not None
    assert finding.match_outcome in ("known", "probable-match-NEEDS-HUMAN", "novel")
