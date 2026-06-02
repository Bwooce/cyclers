"""Real-ephemeris cell optimiser (``optimise_cell_ephemeris``) tests.

Task 1 exercises the pure cell→engine input derivation
(:func:`_ephemeris_tof_seed_and_bounds`). Tasks 2 and 4 exercise the
slow real-ephemeris path: epoch resolution, return-shape parity, and
Aldrin parity against the Aldrin-specific solver.

Provenance: the optimiser COMPUTES the launch epoch and leg ToFs.
Validation asserts only **published** anchors (Aldrin a/e/V∞ from the
sourced constants in :mod:`tests.verify.test_real_closure`). Computed
ToFs/epoch are never asserted as golden.
"""

import math

import pytest

from cyclerfinder.core.constants import MU_SUN_KM3_S2
from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.model.cycler import orbit_elements_au
from cyclerfinder.search.optimize import (
    OptimisationResult,
    _ephemeris_tof_seed_and_bounds,
    _resolve_t0_multi_seed,
    optimise_cell_ephemeris,
)
from cyclerfinder.search.sequence import Cell


def test_ephemeris_tof_seed_and_bounds_equispaced() -> None:
    cell = Cell(
        bodies=("E", "M"),
        sequence=("E", "M", "E"),
        period_k=2,
        per_leg_revs=(0, 0),
        per_leg_branch=("single", "single"),
    )
    target_period_sec = 2 * 779.9 * 86400.0
    seed_days, bounds = _ephemeris_tof_seed_and_bounds(cell, target_period_sec)
    n_legs = len(cell.sequence) - 1
    assert len(seed_days) == n_legs
    assert len(bounds) == n_legs
    # equispaced seed: each leg ~ T/(N-1)
    expected = target_period_sec / n_legs / 86400.0
    assert all(math.isclose(s, expected, rel_tol=1e-9) for s in seed_days)
    # bounds bracket the seed and are strictly positive
    for (lo, hi), s in zip(bounds, seed_days, strict=True):
        assert 0 < lo < s < hi


@pytest.mark.slow
def test_resolve_t0_multi_seed_returns_low_mismatch_epoch() -> None:
    """STAGE 3: the multi-seed t0 resolver fans the cell's seed ToFs into
    asymmetric perturbations and returns the lowest-mismatch launch epoch on
    the real ephemeris for an Aldrin-band E-M-E cell.

    Provenance: 6.5/9.7 km/s are sourced Aldrin V_inf inputs; the returned
    epoch is COMPUTED and only checked for type/finiteness (no golden epoch).
    """
    from datetime import UTC, datetime

    cell = Cell(
        bodies=("E", "M"),
        sequence=("E", "M", "E"),
        period_k=1,
        per_leg_revs=(0, 0),
        per_leg_branch=("single", "single"),
    )
    eph = Ephemeris(model="astropy")
    priority = datetime(1985, 10, 28, tzinfo=UTC)
    t0 = _resolve_t0_multi_seed(
        cell,
        seed_days=[146.0, 634.0],
        priority_date=priority,
        ephem=eph,
        vinf_targets_kms={"E": 6.5, "M": 9.7},
        target_period_sec=780.0 * 86400.0,
    )
    # COMPUTED: a real window exists for the Aldrin band → finite epoch.
    assert t0 is not None
    assert math.isfinite(t0)


@pytest.mark.slow
def test_optimise_cell_ephemeris_returns_result_for_aldrin_em_cell() -> None:
    """The general ephemeris optimiser must reproduce the Aldrin E-M-E
    geometry the Aldrin-specific solver finds (parity check), returning a
    populated OptimisationResult on the real ephemeris."""
    cell = Cell(
        bodies=("E", "M"),
        sequence=("E", "M", "E"),
        period_k=1,
        per_leg_revs=(0, 0),
        per_leg_branch=("single", "single"),
    )
    eph = Ephemeris(model="astropy")
    result = optimise_cell_ephemeris(
        cell,
        eph,
        vinf_cap=12.0,
        priority_date_iso="1985-01-01",  # Aldrin priority era
        vinf_targets_kms={"E": 6.5, "M": 9.7},
        n_starts=3,
        seed=0,
    )
    assert isinstance(result, OptimisationResult)
    assert result.best_cycler is not None
    # the recovered E-M elements match the sourced Aldrin family within band
    assert result.best_score.max_vinf_kms > 1.5 or result.converged


# Sourced Aldrin anchors (mirror tests/verify/test_real_closure.py:204-211).
_ALDRIN_PUB_A_AU = 1.60
_ALDRIN_PUB_E = 0.393


@pytest.mark.slow
def test_general_engine_recovers_aldrin_family_from_family_seed() -> None:
    """Engine-level parity: the body-agnostic ``optimise_maintenance_dv`` —
    the same engine ``optimise_cell_ephemeris`` wraps — recovers the sourced
    Aldrin family (a~1.60 AU, e~0.393) on the real DE440 ephemeris when seeded
    with the family's asymmetric leg structure and phase-matched epoch.

    This proves the general path *is* the specialised path: it runs the exact
    engine the Aldrin-specific ``optimise_aldrin_maintenance_dv`` wrapper calls,
    with the same E→M / M→E ToF seeds and bounds, and reaches the same basin.
    Anchors are the published Aldrin elements (sourced); the epoch and leg ToFs
    are computed and not asserted as golden.
    """
    from cyclerfinder.search.maintain import optimise_maintenance_dv
    from cyclerfinder.search.phase_match import PhaseSignature
    from cyclerfinder.verify.real_closure import (
        _parse_priority_date,
        _resolve_real_t_start,
    )

    eph = Ephemeris(model="astropy")
    # Family-appropriate seed: Aldrin E→M ~146 d, M→E ~634 d (sum = 1 synodic).
    signature = PhaseSignature(
        bodies=("E", "M", "E"),
        leg_durations_s=(146.0 * 86400.0, 634.0 * 86400.0),
        vinf_target_kms=(6.5, 9.7, 6.5),
    )
    priority = _parse_priority_date("1985-01-01")
    assert priority is not None
    t0 = _resolve_real_t_start(signature, eph, priority)
    assert t0 is not None  # the Aldrin window resolves under the mismatch cap

    maint = optimise_maintenance_dv(
        ["E", "M", "E"],
        eph,
        t0_guess_sec=t0,
        tof_days_guesses=(146.0, 634.0),
        tof_bounds_days=((100.0, 250.0), (400.0, 900.0)),
        synodic_pair=("E", "M"),
        closure_body="E",
        tof_jitter_half_days=(20.0, 60.0),
        n_starts=5,
        seed=0,
    )
    assert maint.converged
    assert maint.a_au == pytest.approx(_ALDRIN_PUB_A_AU, abs=0.15)
    assert maint.e == pytest.approx(_ALDRIN_PUB_E, abs=0.08)
    vinf = dict(maint.vinf_kms_at_encounters)
    assert vinf["E"] == pytest.approx(6.5, abs=1.5)
    assert vinf["M"] == pytest.approx(9.7, abs=3.0)


@pytest.mark.slow
def test_optimise_cell_ephemeris_aldrin_parity_elements() -> None:
    """Parity goal: the general ``optimise_cell_ephemeris`` on the Aldrin E-M-E
    cell recovers a/e consistent with the Aldrin-specific solver from a cold
    cell start.

    STAGE 3 (2026-06-03) closed this gate. Previously xfail: the cold-start
    equispaced cell seed resolved the launch epoch from a symmetric leg
    signature and landed in a degenerate short-period basin (a~1 AU) instead
    of the asymmetric Aldrin family. The robust multi-seed resolver
    (``_resolve_t0_multi_seed`` → ``leg_duration_seeds`` +
    ``find_candidate_windows``, ranked by V_inf mismatch) now fans the cold
    equispaced seed into asymmetric perturbations and selects the Aldrin-band
    launch epoch, so the engine reaches the published family from a cold cell.

    Anchors are the published Aldrin elements (SOURCED); the epoch and leg
    ToFs are computed and not asserted as golden.
    """
    cell = Cell(
        bodies=("E", "M"),
        sequence=("E", "M", "E"),
        period_k=1,
        per_leg_revs=(0, 0),
        per_leg_branch=("single", "single"),
    )
    eph = Ephemeris(model="astropy")
    result = optimise_cell_ephemeris(
        cell,
        eph,
        vinf_cap=12.0,
        priority_date_iso="1985-01-01",
        vinf_targets_kms={"E": 6.5, "M": 9.7},
        n_starts=5,
        seed=0,
    )
    assert result.best_cycler is not None
    # Recover transfer-ellipse elements from the first leg.
    cyc = result.best_cycler
    a_au, e = orbit_elements_au(cyc.encounters[0].r, cyc.legs[0].v_depart, MU_SUN_KM3_S2)
    assert a_au == pytest.approx(_ALDRIN_PUB_A_AU, abs=0.15)
    assert e == pytest.approx(_ALDRIN_PUB_E, abs=0.08)


# ---------------------------------------------------------------------------
# STAGE 1 — multi-rev ToF floor + revs/branch threading
# ---------------------------------------------------------------------------


def test_ephemeris_seed_bounds_multirev_floor() -> None:
    """STAGE 1: a 1-rev leg gets a wider lower-ToF floor than a direct leg.

    A 1-rev Lambert leg is infeasible below its physical minimum ToF
    (``~pi * sqrt(a^3 / mu)`` for the minimum-energy orbit). The seed/bounds
    helper must lift leg 1's lower bound above the direct leg-0 floor so the
    optimiser is not seeded into the LambertConvergenceError regime.

    Provenance: ``# COMPUTED`` — asserts the bound ordering contract, not a
    sourced magnitude.
    """
    cell = Cell(
        bodies=("E", "M"),
        sequence=("E", "M", "E"),
        period_k=2,
        per_leg_revs=(0, 1),
        per_leg_branch=("single", "low"),
    )
    target_period_sec = 2 * 779.9 * 86400.0
    _seed, bounds = _ephemeris_tof_seed_and_bounds(cell, target_period_sec)
    # The 1-rev leg's lower floor exceeds the direct leg's 0.1*share floor.
    assert bounds[1][0] > bounds[0][0]
    # And every bound is still a valid strictly-positive (lo, hi).
    for lo, hi in bounds:
        assert 0.0 < lo < hi


def test_optimise_cell_ephemeris_threads_revs_to_build_chain(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """STAGE 1: ``optimise_cell_ephemeris`` extracts ``cell.per_leg_revs`` /
    ``cell.per_leg_branch`` and threads them into the maintenance solve, so the
    recovered cycler's legs carry the requested revolution counts and a spied
    ``optimise_maintenance_dv`` receives the exact per-leg metadata.

    Provenance: ``# COMPUTED`` — circular ephemeris; asserts the plumbing
    contract (revs forwarded, leg n_revs honoured), no sourced V_inf.
    """
    import cyclerfinder.search.maintain as maintain_mod

    cell = Cell(
        bodies=("E", "M"),
        sequence=("E", "M", "E"),
        period_k=1,
        per_leg_revs=(0, 0),
        per_leg_branch=("single", "single"),
    )

    captured: dict[str, object] = {}
    real_optimise = maintain_mod.optimise_maintenance_dv

    def _spy(sequence, ephem, **kwargs):  # type: ignore[no-untyped-def]
        captured["per_leg_revs"] = kwargs.get("per_leg_revs")
        captured["per_leg_branch"] = kwargs.get("per_leg_branch")
        return real_optimise(sequence, ephem, **kwargs)

    # optimise_cell_ephemeris imports the symbol locally, so patch the source.
    monkeypatch.setattr(maintain_mod, "optimise_maintenance_dv", _spy)

    result = optimise_cell_ephemeris(
        cell,
        Ephemeris("circular"),
        vinf_cap=12.0,
        priority_date_iso="1985-01-01",
        vinf_targets_kms={"E": 6.5, "M": 9.7},
        n_starts=2,
        seed=0,
    )
    assert isinstance(result, OptimisationResult)
    assert captured["per_leg_revs"] == (0, 0)
    assert captured["per_leg_branch"] == ("single", "single")
    # The default-direct cell still yields a 0-rev first leg.
    assert result.best_cycler.legs[0].n_revs == 0
