"""Tests for the precursor-MGA matcher (#302 / #289 Phase 4).

These tests exercise:

  1. Shape invariants on :class:`PrecursorMatch` and the closure record.
  2. Error handling on bad cycler_id / missing V_inf / unknown body codes.
  3. The matcher's filtering: only Earth-launched chains terminating at the
     cycler's first encounter body within the V_inf tolerance survive.
  4. The :func:`precursor_match_to_jsonl_record` shape (a downstream
     consumer reads this JSON, so the field set is part of the contract).
  5. A toy-catalogue end-to-end run on a synthetic 1-leg E->E "precursor"
     that closes trivially under the analytic-circular ephemeris.

The matcher's *physical* fidelity is validated by the Phase-1 / Phase-3
test surfaces (Tito free-return closure, Galileo VEEGA structural
reproduction).  These tests cover the wrapping layer's invariants and
contract, not the underlying dynamics.
"""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path

import pytest

from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.data.catalog import load_catalog
from cyclerfinder.genome.epoch_aware_genome import (
    EpochLockedClosure,
    EpochLockedTrajectory,
    close_epoch_locked,
)
from cyclerfinder.search.literature_check import SearchResult
from cyclerfinder.search.precursor_matcher import (
    PrecursorMatch,
    _epoch_alignment_score,
    _first_encounter_body_and_vinf,
    _signature_from_candidate,
    _terminal_vinf_kms,
    find_cycler_precursors,
    precursor_match_to_jsonl_record,
)

CATALOGUE_PATH = Path(__file__).resolve().parents[2] / "data" / "catalogue.yaml"


def _empty_search(_query: str) -> Sequence[SearchResult]:
    """Always-empty search — the literature_check returns 'inconclusive'."""
    return []


def _fresh_search(_query: str) -> Sequence[SearchResult]:
    """Always returns a single irrelevant hit — drives ``not-found``."""
    return [
        SearchResult(
            title="Unrelated topic about pizza toppings",
            url="https://example.com/pizza",
            snippet="Nothing to do with trajectories or gravity assists.",
        )
    ]


# --------------------------------------------------------------------------- #
# Catalogue-row helpers
# --------------------------------------------------------------------------- #


def test_first_encounter_body_and_vinf_extracts_seed() -> None:
    """The matcher reads the cycler's first encounter body + V_inf from YAML."""
    cat = load_catalog(CATALOGUE_PATH)
    aldrin = cat.by_id["aldrin-classic-em-k1-outbound"]
    body, vinf = _first_encounter_body_and_vinf(aldrin)
    assert body == "E"
    # Russell 2004 Table 3.4 cycler 1.0.1.-1: Earth V_inf = 6.5 km/s.
    assert vinf == pytest.approx(6.5, abs=1.0e-6)


def test_first_encounter_body_and_vinf_s1l1_seed() -> None:
    """S1L1 catalogue row's first encounter is Earth at 5.65 km/s."""
    cat = load_catalog(CATALOGUE_PATH)
    s1l1 = cat.by_id["s1l1-2syn-em-cpom"]
    body, vinf = _first_encounter_body_and_vinf(s1l1)
    assert body == "E"
    assert vinf == pytest.approx(5.65, abs=1.0e-6)


def test_first_encounter_body_and_vinf_raises_on_unknown_body() -> None:
    """Construct a fake CatalogueEntry whose first encounter body is unknown."""

    class _FakeEntry:
        id = "fake-row"
        vinf_kms_at_encounters = (("UnknownPlanet", 6.5),)

    with pytest.raises(ValueError, match="not in the matcher's body-code"):
        _first_encounter_body_and_vinf(_FakeEntry())  # type: ignore[arg-type]


def test_first_encounter_body_and_vinf_raises_on_missing_vinf() -> None:
    """Empty vinf_kms_at_encounters block is a hard error (no seed to match)."""

    class _FakeEntry:
        id = "fake-row"
        vinf_kms_at_encounters = ()

    with pytest.raises(ValueError, match="no vinf_kms_at_encounters"):
        _first_encounter_body_and_vinf(_FakeEntry())  # type: ignore[arg-type]


# --------------------------------------------------------------------------- #
# Epoch alignment scoring
# --------------------------------------------------------------------------- #


def test_epoch_alignment_score_inside_window_is_one() -> None:
    """Launch epoch inside the target window scores 1.0."""
    s = _epoch_alignment_score(
        "2030-06-15T00:00:00",
        ("2030-01-01T00:00:00", "2030-12-31T00:00:00"),
    )
    assert s == pytest.approx(1.0, abs=1.0e-6)


def test_epoch_alignment_score_far_outside_falls_to_zero() -> None:
    """A launch epoch one full synodic period before the window scores ~0."""
    s = _epoch_alignment_score(
        "2028-01-01T00:00:00",
        ("2030-06-01T00:00:00", "2030-07-01T00:00:00"),
    )
    # ~880 days before > 779.8 d => clamped to 0.
    assert s == pytest.approx(0.0, abs=1.0e-3)


def test_epoch_alignment_score_no_window_is_zero() -> None:
    """No target window => no signal, score is 0 (caller must not gate)."""
    s = _epoch_alignment_score("2030-06-15T00:00:00", None)
    assert s == 0.0


# --------------------------------------------------------------------------- #
# PrecursorMatch shape
# --------------------------------------------------------------------------- #


def _trivial_one_leg_trajectory_and_closure(
    eph: Ephemeris,
) -> tuple[EpochLockedTrajectory, EpochLockedClosure]:
    """Build a trivial 1-leg E->E synthetic precursor used by the shape tests.

    Closes by construction on the analytic-circular ephemeris (an E->E hop
    at a heliocentric resonance with TOF = 1 Earth year).  Returns
    (trajectory, closure) for use by the shape tests.
    """
    # 365.25-day E->E Hohmann at trivial V_inf (the analytic-circular
    # ephemeris's E orbit is exactly circular at 1.0 AU, so a 1-year arc
    # is an exact resonance with V_inf -> 0).
    traj = EpochLockedTrajectory(
        sequence=("E", "E"),
        leg_tofs_days=(365.25,),
        vinf_kms_at_encounters=(0.05, 0.05),
        launch_epoch_utc="2030-01-01T00:00:00",
        orbit_class="precursor_mga",
        n_returns=1,
        validity_window_start_utc="2030-01-01T00:00:00",
        validity_window_end_utc="2031-01-01T00:00:00",
        inserts_into="aldrin-classic-em-k1-outbound",
        notes="shape-test synthetic precursor (#302)",
    )
    closure = close_epoch_locked(
        traj,
        eph,
        closure_tol_kms=10.0,
        flyby_continuity_tol_kms=10.0,
        independent_cross_check=False,
        independent_tol_kms=10.0,
    )
    return traj, closure


def test_precursor_match_quality_score_weights() -> None:
    """``quality_score`` = vinf_match + closure + 2*flyby_continuity."""
    eph = Ephemeris(model="circular")
    traj, closure = _trivial_one_leg_trajectory_and_closure(eph)
    from cyclerfinder.search.literature_check import LiteratureCheckResult

    match = PrecursorMatch(
        candidate=traj,
        cycler_id="aldrin-classic-em-k1-outbound",
        cycler_seed_vinf_kms=6.5,
        vinf_match_residual_kms=0.1,
        epoch_alignment_score=0.5,
        closure=closure,
        literature_check=LiteratureCheckResult(
            status="not-found",
            citation=None,
            doi=None,
            confidence=0.0,
            query_trail=[],
        ),
    )
    expected = (
        match.vinf_match_residual_kms
        + closure.closure_residual_kms
        + 2.0 * closure.flyby_continuity_max_dv_kms
    )
    assert match.quality_score() == pytest.approx(expected, abs=1.0e-12)


def test_precursor_match_is_literature_fresh() -> None:
    """``is_literature_fresh`` only fires on the clean ``not-found`` verdict."""
    eph = Ephemeris(model="circular")
    traj, closure = _trivial_one_leg_trajectory_and_closure(eph)
    from cyclerfinder.search.literature_check import LiteratureCheckResult
    from cyclerfinder.search.literature_check import Status as LitStatus

    def _make_match(status: LitStatus) -> PrecursorMatch:
        return PrecursorMatch(
            candidate=traj,
            cycler_id="aldrin-classic-em-k1-outbound",
            cycler_seed_vinf_kms=6.5,
            vinf_match_residual_kms=0.1,
            epoch_alignment_score=0.0,
            closure=closure,
            literature_check=LiteratureCheckResult(
                status=status,
                citation=None,
                doi=None,
                confidence=0.0,
                query_trail=[],
            ),
        )

    assert _make_match("not-found").is_literature_fresh()
    assert not _make_match("published").is_literature_fresh()
    assert not _make_match("inconclusive").is_literature_fresh()


def test_precursor_match_jsonl_record_shape() -> None:
    """The JSONL record carries every field the doc + downstream gauntlet need."""
    eph = Ephemeris(model="circular")
    traj, closure = _trivial_one_leg_trajectory_and_closure(eph)
    from cyclerfinder.search.literature_check import LiteratureCheckResult

    match = PrecursorMatch(
        candidate=traj,
        cycler_id="aldrin-classic-em-k1-outbound",
        cycler_seed_vinf_kms=6.5,
        vinf_match_residual_kms=0.1,
        epoch_alignment_score=0.42,
        closure=closure,
        literature_check=LiteratureCheckResult(
            status="not-found",
            citation=None,
            doi=None,
            confidence=0.0,
            query_trail=[],
        ),
    )
    rec = precursor_match_to_jsonl_record(match)
    # Top-level fields the doc + caller depend on.
    for key in (
        "cycler_id",
        "cycler_seed_vinf_kms",
        "vinf_match_residual_kms",
        "epoch_alignment_score",
        "quality_score",
        "candidate",
        "closure",
        "literature_check",
        "is_literature_fresh",
    ):
        assert key in rec, f"JSONL record missing top-level key {key!r}"
    # Candidate sub-fields the gauntlet needs to reproduce the closure.
    for key in (
        "sequence",
        "leg_tofs_days",
        "vinf_kms_at_encounters",
        "launch_epoch_utc",
        "orbit_class",
        "inserts_into",
    ):
        assert key in rec["candidate"], f"candidate block missing {key!r}"
    # Closure sub-fields the doc cites.
    for key in (
        "closure_residual_kms",
        "flyby_continuity_max_dv_kms",
        "per_encounter_vinf_kms",
        "independent_check_residual_kms",
        "converged",
    ):
        assert key in rec["closure"], f"closure block missing {key!r}"
    # Literature check sub-fields the gauntlet flips on.
    for key in ("status", "citation", "doi", "confidence", "notes"):
        assert key in rec["literature_check"], f"literature_check block missing {key!r}"
    assert rec["is_literature_fresh"] is True


# --------------------------------------------------------------------------- #
# End-to-end smoke test
# --------------------------------------------------------------------------- #


def test_find_cycler_precursors_smoke_unknown_id_raises() -> None:
    cat = load_catalog(CATALOGUE_PATH)
    eph = Ephemeris(model="circular")
    with pytest.raises(ValueError, match="not in catalogue"):
        find_cycler_precursors(
            "this-cycler-id-does-not-exist",
            cat,
            eph,
            launch_window=("2030-01-01T00:00:00", "2030-02-01T00:00:00"),
        )


def test_find_cycler_precursors_smoke_returns_list() -> None:
    """A short, deliberately-narrow scan returns a (possibly empty) list.

    This smoke test runs the full pipeline against the Aldrin cycler row
    with a very narrow V_inf grid and short launch window — by design,
    we expect 0 or a handful of survivors (most chains won't close at
    this grain), and the test only verifies the return type and basic
    record invariants, not a specific count.

    Closure tolerance is wide (5 km/s) and flyby gate is wide (1 km/s) so
    a survivor doesn't require a publication-grade match — the matcher's
    SHAPE is being verified, not its dynamical fidelity.
    """
    cat = load_catalog(CATALOGUE_PATH)
    eph = Ephemeris(model="circular")

    matches = find_cycler_precursors(
        cycler_id="aldrin-classic-em-k1-outbound",
        catalogue=cat,
        ephemeris=eph,
        launch_window=("2030-01-01T00:00:00", "2030-02-01T00:00:00"),
        max_legs=2,
        intermediate_bodies=("V",),
        vinf_terminal_tol_kms=2.0,
        vinf_grid_kms=(5.0, 6.0, 7.0),
        tof_box_days_per_leg=(60.0, 400.0),
        epoch_step_days=30.0,
        tof_optimise=False,
        closure_tol_kms=5.0,
        flyby_continuity_tol_kms=1.0,
        independent_cross_check=False,
        independent_tol_kms=10.0,
        multi_shell=True,
        max_candidates_to_validate=20,
        literature_check_search=_empty_search,
    )
    assert isinstance(matches, list)
    for m in matches:
        assert isinstance(m, PrecursorMatch)
        assert m.cycler_id == "aldrin-classic-em-k1-outbound"
        # The terminal body must be Earth (the cycler's first encounter).
        assert m.candidate.sequence[-1] == "E"
        # The launch body is always Earth.
        assert m.candidate.sequence[0] == "E"
        # The closure must converge under the loose smoke gates.
        assert m.closure.closure_residual_kms <= 5.0
        # Literature check with empty search => inconclusive.
        assert m.literature_check.status == "inconclusive"


def test_signature_from_candidate_shape() -> None:
    """The literature-check signature carries the candidate's structure."""
    traj = EpochLockedTrajectory(
        sequence=("E", "V", "E", "M"),
        leg_tofs_days=(180.0, 180.0, 180.0),
        vinf_kms_at_encounters=(4.5, 4.5, 6.0, 6.5),
        launch_epoch_utc="2030-01-01T00:00:00",
        orbit_class="precursor_mga",
        n_returns=1,
        validity_window_start_utc="2030-01-01T00:00:00",
        validity_window_end_utc="2032-01-01T00:00:00",
        inserts_into="aldrin-classic-em-k1-outbound",
    )
    sig = _signature_from_candidate(traj)
    assert sig.primary == "Sun"
    assert sig.sequence == ("E", "V", "E", "M")
    assert sig.vinf_per_encounter_kms == (4.5, 4.5, 6.0, 6.5)


def test_terminal_vinf_kms_reads_last_encounter() -> None:
    eph = Ephemeris(model="circular")
    _, closure = _trivial_one_leg_trajectory_and_closure(eph)
    # The trivial 1-year E->E hop on circular ephemeris has |V_inf| ~ 0 at
    # both encounters; the terminal one is the second entry.
    assert _terminal_vinf_kms(closure) == pytest.approx(
        closure.per_encounter_vinf_kms[-1], abs=1.0e-12
    )
