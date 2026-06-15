"""Epoch-aware (window-bounded) genome data model + closure framework tests.

Phase 1 of task #289 (substrate for the four-class catalogue scope at schema
v4.7). These tests exercise the Tito 2018 Mars free-return as the canonical
end-to-end ``mga_tour`` closure, plus the class-invariant / window-search /
negative gates.

The Tito reproduction uses the catalogue row's ``vinf_kms_at_encounters`` as the
EXPECTED side, which traces to Tito 2013 Tables IV via the row's
``source_quotes`` (sourced-only golden-test discipline:
``feedback_golden_tests_sourced_only``). The ACTUAL side is computed by
:func:`close_epoch_locked` using DE440 — the documented residual budget is
<1.5% on any single ``|V_inf|``, which the row records is dominated by the
DE421->DE440 ephemeris difference accumulated over the 274-day return arc.

Slow path: the astropy ephemeris backend pulls DE440 Chebyshev kernels at
construction (a few hundred ms) and runs full JPL-grade body-state queries
per leg. A single closure costs ~0.5 s; the window-search test marches a
21-grid-point ±5-day box and runs ~10 s. Both are well below the
``--timeout=120`` repository CI ceiling.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pytest
import yaml  # type: ignore[import-untyped]

from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.genome.epoch_aware_genome import (
    EpochLockedClosure,
    EpochLockedTrajectory,
    close_epoch_locked,
    search_validity_window,
)

CATALOGUE_PATH = Path(__file__).resolve().parents[2] / "data" / "catalogue.yaml"
TITO_ID = "tito-2018-mars-free-return"

# --------------------------------------------------------------------------- #
# Catalogue fixture (Tito row).
# --------------------------------------------------------------------------- #


@pytest.fixture(scope="module")
def tito_row() -> dict[str, Any]:
    rows: list[dict[str, Any]] = yaml.safe_load(CATALOGUE_PATH.read_text())
    for row in rows:
        if row.get("id") == TITO_ID:
            return row
    raise AssertionError(
        f"Tito row {TITO_ID!r} missing from data/catalogue.yaml; "
        "this test depends on the schema v4.7 (#294) admission of the "
        "tito-2018-mars-free-return mga_tour row.",
    )


@pytest.fixture(scope="module")
def tito_trajectory(tito_row: dict[str, Any]) -> EpochLockedTrajectory:
    """Build the EpochLockedTrajectory from the (sourced) catalogue row.

    The Tito 2013 published epochs (Tables III, transcribed via the
    reproduction script) are:

      depart Earth : 2018-01-05 07:00:00.000 UTC
      flyby  Mars  : 2018-08-20 08:18:19.619 UTC
      return Earth : 2019-05-21 13:52:48.012 UTC

    The catalogue row stores the validity window endpoints but not the
    intermediate flyby UTC, so we supply the published flight times in days
    (Tables III: 227.05439374 d, 274.23227306 d). The published periapsis
    is 100 km (Tito Fig. 3 text), so the test sets a per-encounter override
    at the Mars node — the engineering default 300 km would tighten the
    bend cone enough to push the test past the 0.05 km/s continuity gate.
    """
    vinf_entries = tito_row["vinf_kms_at_encounters"]
    vinf_kms = tuple(float(entry["vinf_kms"]) for entry in vinf_entries)
    return EpochLockedTrajectory(
        sequence=("E", "M", "E"),
        leg_tofs_days=(227.05439374, 274.23227306),
        vinf_kms_at_encounters=vinf_kms,
        launch_epoch_utc=tito_row["validity_window"]["start"],
        orbit_class="mga_tour",
        n_returns=int(tito_row["n_returns"]),
        validity_window_start_utc=tito_row["validity_window"]["start"],
        validity_window_end_utc=tito_row["validity_window"]["end"],
        periapsis_altitudes_km=(None, 100.0, None),
        ephemeris="DE440",
        notes="Tito 2018 reproduction integration test",
    )


@pytest.fixture(scope="module")
def astropy_ephemeris() -> Ephemeris:
    """Single DE440-backed ephemeris instance reused across the module."""
    return Ephemeris("astropy")


@pytest.fixture(scope="module")
def tito_closure(
    tito_trajectory: EpochLockedTrajectory,
    astropy_ephemeris: Ephemeris,
) -> EpochLockedClosure:
    """Run the Tito closure once and share across the gate-level tests."""
    return close_epoch_locked(
        tito_trajectory,
        astropy_ephemeris,
        closure_tol_kms=0.5,
        flyby_continuity_tol_kms=0.05,
        independent_cross_check=True,
        independent_tol_kms=0.1,
    )


# --------------------------------------------------------------------------- #
# Test 1: Tito 2018 reproduction (the full Lambert + flyby + ephemeris loop).
# --------------------------------------------------------------------------- #


def test_tito_reproduction_closes_within_published_bounds(
    tito_closure: EpochLockedClosure,
    tito_row: dict[str, Any],
) -> None:
    """Tito 2018 closes on DE440 within the row-documented <1.5% on any |V_inf|.

    EXPECTED side: catalogue row vinf_kms_at_encounters (sourced from Tito
    2013 Tables IV via source_quotes). ACTUAL side: our closure on DE440.
    Per-encounter relative residual must stay below 1.5% — the row's
    documented worst is 1.42% on C3_arrival_Earth (i.e. ~0.71% on |V_inf|),
    dominated by the DE421->DE440 difference.
    """
    expected = [float(e["vinf_kms"]) for e in tito_row["vinf_kms_at_encounters"]]
    ours = list(tito_closure.per_encounter_vinf_kms)
    assert len(ours) == len(expected) == 3

    rel_residuals = [abs(o - e) / e for o, e in zip(ours, expected, strict=True)]
    worst_rel = max(rel_residuals)
    assert worst_rel < 0.015, (
        f"Tito |V_inf| relative residual {worst_rel:.4%} exceeds 1.5% gate; "
        f"per-encounter (ours, pub, rel%): "
        f"{list(zip(ours, expected, [f'{r:.4%}' for r in rel_residuals], strict=True))!r}"
    )

    # And the converged flag should be set under the row-grade tolerances.
    assert tito_closure.converged, (
        f"Tito closure did not converge: closure_residual_kms="
        f"{tito_closure.closure_residual_kms:.4f}, "
        f"flyby_continuity_max_dv_kms="
        f"{tito_closure.flyby_continuity_max_dv_kms:.4f}, "
        f"independent_check_residual_kms="
        f"{tito_closure.independent_check_residual_kms!r}",
    )


# --------------------------------------------------------------------------- #
# Test 2: Mars ballistic-continuity gate.
# --------------------------------------------------------------------------- #


def test_tito_mars_flyby_ballistic_continuity(tito_closure: EpochLockedClosure) -> None:
    """Tito's Mars flyby ballistic-continuity Delta V is below 0.05 km/s.

    The reproduction script reports 0.0105 km/s at the Mars flyby on DE440;
    the row docstring also asserts ballistic continuity is confirmed
    (required bend 33.4 deg < cone 34.2 deg at the published 100 km
    periapsis). The 0.05 km/s gate is Tito-grade -- mga_tour rows that fail
    this gate are powered, not ballistic.
    """
    assert tito_closure.flyby_continuity_max_dv_kms < 0.05, (
        f"Mars flyby continuity {tito_closure.flyby_continuity_max_dv_kms:.4f} "
        f"km/s exceeds 0.05 km/s ballistic gate"
    )


# --------------------------------------------------------------------------- #
# Test 3: Independent integrator cross-check.
# --------------------------------------------------------------------------- #


def test_tito_independent_cross_check_agrees(tito_closure: EpochLockedClosure) -> None:
    """The Lambert (BVP) and Kepler (IVP) propagators agree on the Tito conic.

    The orbit-closure discipline (``feedback_orbit_closure_discipline``)
    requires an independent solver cross-check on every closure. The
    universal-variable Lambert solver and the universal-variable Kepler
    propagator have different residual functions, different Newton loops,
    and different convergence criteria; their agreement on the same conic
    is the cross-check.
    """
    assert tito_closure.independent_check_residual_kms is not None, (
        "independent cross-check was requested but produced no residual"
    )
    assert tito_closure.independent_check_residual_kms < 0.1, (
        f"independent-check residual "
        f"{tito_closure.independent_check_residual_kms:.4f} km/s "
        f"exceeds 0.1 km/s gate"
    )


# --------------------------------------------------------------------------- #
# Test 4: Class invariants.
# --------------------------------------------------------------------------- #


def test_strict_cycler_class_rejected() -> None:
    """``orbit_class='cycler'`` is rejected by EpochLockedTrajectory.

    Strict cyclers are periodic and not epoch-locked; they belong in the
    Cycler class, not this data model.
    """
    with pytest.raises(ValueError, match=r"orbit_class.*'cycler'"):
        EpochLockedTrajectory(
            sequence=("E", "M", "E"),
            leg_tofs_days=(200.0, 200.0),
            vinf_kms_at_encounters=(6.0, 5.0, 9.0),
            launch_epoch_utc="2018-01-05T07:00:00Z",
            orbit_class="cycler",  # type: ignore[arg-type]
            n_returns=1,
            validity_window_start_utc="2018-01-05T07:00:00Z",
            validity_window_end_utc="2019-05-21T13:52:48Z",
        )


def test_precursor_mga_requires_inserts_into() -> None:
    """``precursor_mga`` must declare which cycler it inserts into."""
    with pytest.raises(ValueError, match=r"precursor_mga.*inserts_into"):
        EpochLockedTrajectory(
            sequence=("E", "M", "E"),
            leg_tofs_days=(200.0, 200.0),
            vinf_kms_at_encounters=(6.0, 5.0, 9.0),
            launch_epoch_utc="2018-01-05T07:00:00Z",
            orbit_class="precursor_mga",
            n_returns=1,
            validity_window_start_utc="2018-01-05T07:00:00Z",
            validity_window_end_utc="2019-05-21T13:52:48Z",
            # inserts_into deliberately omitted
        )


def test_precursor_mga_accepts_inserts_into() -> None:
    """A precursor_mga with a non-empty inserts_into id constructs cleanly."""
    traj = EpochLockedTrajectory(
        sequence=("E", "M", "E"),
        leg_tofs_days=(200.0, 200.0),
        vinf_kms_at_encounters=(6.0, 5.0, 9.0),
        launch_epoch_utc="2018-01-05T07:00:00Z",
        orbit_class="precursor_mga",
        n_returns=1,
        validity_window_start_utc="2018-01-05T07:00:00Z",
        validity_window_end_utc="2019-05-21T13:52:48Z",
        inserts_into="aldrin-1l1",
    )
    assert traj.inserts_into == "aldrin-1l1"


def test_n_returns_must_be_positive_integer() -> None:
    """``n_returns`` must be a positive integer."""
    common = dict(
        sequence=("E", "M", "E"),
        leg_tofs_days=(200.0, 200.0),
        vinf_kms_at_encounters=(6.0, 5.0, 9.0),
        launch_epoch_utc="2018-01-05T07:00:00Z",
        orbit_class="mga_tour",
        validity_window_start_utc="2018-01-05T07:00:00Z",
        validity_window_end_utc="2019-05-21T13:52:48Z",
    )
    with pytest.raises(ValueError, match="n_returns"):
        EpochLockedTrajectory(**common, n_returns=0)  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="n_returns"):
        EpochLockedTrajectory(**common, n_returns=-1)  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="n_returns"):
        EpochLockedTrajectory(**common, n_returns=1.5)  # type: ignore[arg-type]


def test_leg_tofs_shape_must_match_sequence() -> None:
    """``leg_tofs_days`` length must equal ``len(sequence) - 1``."""
    with pytest.raises(ValueError, match="leg_tofs_days"):
        EpochLockedTrajectory(
            sequence=("E", "M", "E"),
            leg_tofs_days=(200.0,),  # one short
            vinf_kms_at_encounters=(6.0, 5.0, 9.0),
            launch_epoch_utc="2018-01-05T07:00:00Z",
            orbit_class="mga_tour",
            n_returns=1,
            validity_window_start_utc="2018-01-05T07:00:00Z",
            validity_window_end_utc="2019-05-21T13:52:48Z",
        )


def test_vinf_shape_must_match_sequence() -> None:
    """``vinf_kms_at_encounters`` length must equal ``len(sequence)``."""
    with pytest.raises(ValueError, match="vinf_kms_at_encounters"):
        EpochLockedTrajectory(
            sequence=("E", "M", "E"),
            leg_tofs_days=(200.0, 200.0),
            vinf_kms_at_encounters=(6.0, 5.0),  # one short
            launch_epoch_utc="2018-01-05T07:00:00Z",
            orbit_class="mga_tour",
            n_returns=1,
            validity_window_start_utc="2018-01-05T07:00:00Z",
            validity_window_end_utc="2019-05-21T13:52:48Z",
        )


def test_unknown_body_code_rejected() -> None:
    """Unknown body codes in sequence raise ValueError."""
    with pytest.raises(ValueError, match="unknown body code"):
        EpochLockedTrajectory(
            sequence=("E", "X", "E"),
            leg_tofs_days=(200.0, 200.0),
            vinf_kms_at_encounters=(6.0, 5.0, 9.0),
            launch_epoch_utc="2018-01-05T07:00:00Z",
            orbit_class="mga_tour",
            n_returns=1,
            validity_window_start_utc="2018-01-05T07:00:00Z",
            validity_window_end_utc="2019-05-21T13:52:48Z",
        )


# --------------------------------------------------------------------------- #
# Test 5: Window search (Tito is a tight ±days window).
# --------------------------------------------------------------------------- #


def test_search_validity_window_tito_finds_nominal(
    tito_trajectory: EpochLockedTrajectory,
    astropy_ephemeris: Ephemeris,
) -> None:
    """search_validity_window finds the nominal Tito launch as a converged point.

    With the published flight times held fixed, slipping the launch by a
    day or two keeps the geometry close enough that closure can still pass
    tight tolerances (the Tito launch window is actually several days
    wide; the Patel-Longuski "two times every 15 years" phrasing describes
    the *recurrence between* windows, not the width of a single one).

    The discriminating assertion is therefore: the converged set is
    *non-empty*, *bounded* (much smaller than the 21-grid-point window),
    and *includes the nominal launch epoch*.
    """
    closures = search_validity_window(
        tito_trajectory,
        astropy_ephemeris,
        epoch_grid_step_days=1.0,
        epoch_grid_padding_days=10.0,
        closure_tol_kms=0.5,
        flyby_continuity_tol_kms=0.05,
        independent_cross_check=False,  # speed: skip the cross-check on grid sweep
        independent_tol_kms=0.1,
    )
    assert len(closures) >= 1, "expected at least one converged closure in the ±10-day Tito window"
    n_grid_points = 1 + 2 * 10  # ±10 days at 1 d step
    assert len(closures) < n_grid_points, (
        f"closure window should be bounded (a 'window' implies far-off epochs "
        f"don't close); got {len(closures)} / {n_grid_points} grid points "
        f"converged"
    )
    converged_epochs = [c.trajectory.launch_epoch_utc for c in closures]
    assert any(e.startswith("2018-01-05") for e in converged_epochs), (
        f"nominal Tito launch (2018-01-05) not among converged epochs: {converged_epochs!r}"
    )


# --------------------------------------------------------------------------- #
# Test 6: Negative — a bogus V_inf tuple should NOT close.
# --------------------------------------------------------------------------- #


def test_bogus_vinf_does_not_close(
    tito_trajectory: EpochLockedTrajectory,
    astropy_ephemeris: Ephemeris,
) -> None:
    """A trajectory carrying a deliberately wrong V_inf tuple fails the closure gate.

    The Lambert + ephemeris machinery doesn't know what the published
    V_inf is; it just computes one. The closure gate compares our computed
    V_inf against ``vinf_kms_at_encounters``. If we deliberately store a
    tuple far from the real answer, the closure residual blows past the
    0.5 km/s gate, so ``converged`` must be False.
    """
    # Each entry deliberately off by several km/s — far beyond the closure gate.
    bogus = EpochLockedTrajectory(
        sequence=tito_trajectory.sequence,
        leg_tofs_days=tito_trajectory.leg_tofs_days,
        vinf_kms_at_encounters=(1.0, 1.0, 1.0),
        launch_epoch_utc=tito_trajectory.launch_epoch_utc,
        orbit_class=tito_trajectory.orbit_class,
        n_returns=tito_trajectory.n_returns,
        validity_window_start_utc=tito_trajectory.validity_window_start_utc,
        validity_window_end_utc=tito_trajectory.validity_window_end_utc,
        periapsis_altitudes_km=tito_trajectory.periapsis_altitudes_km,
    )
    closure = close_epoch_locked(
        bogus,
        astropy_ephemeris,
        closure_tol_kms=0.5,
        flyby_continuity_tol_kms=0.05,
        independent_cross_check=False,
        independent_tol_kms=0.1,
    )
    assert not closure.converged, (
        "bogus V_inf tuple should not pass the closure gate; "
        f"got closure_residual_kms={closure.closure_residual_kms:.4f}"
    )
    # The closure residual is large by construction — Tito's |V_inf| are all
    # > 5 km/s and the bogus tuple is 1.0 km/s, so the worst gap is > 4 km/s.
    assert closure.closure_residual_kms > 4.0
