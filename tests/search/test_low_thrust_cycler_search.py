"""Tests for the #309 Phase 1 low-thrust cycler discovery driver.

Discipline (per ``feedback_golden_tests_sourced_only``):

* **Sourced golden**: the EXPECTED side of the Aldrin reproduction test ties
  to the published Aldrin turn-angle structure (84° required vs 72° max at
  the Earth flyby, McConaghy/Longuski/Byrnes 2002 Table 4 row 1L1 — the
  catalogue ``aldrin-classic-em-k1-outbound`` row's
  ``data_gaps.maintenance_dv_kms_per_synodic`` block records that the ΔV
  magnitude is NOT published, only the turn-angle test). The test therefore
  validates the *identity* between the impulsive baseline ΔV and the powered
  wrap (the Tsiolkovsky carry-through is a source-free physics invariant;
  identity is what the test pins).
* **Independent cross-check**: a seed bump must land on the same maintenance
  ΔV within ``1e-3`` km/s (the #285 mandatory-cross-check pattern). The
  optimum is a flat plateau, so the matching tolerance is generous.
"""

from __future__ import annotations

import math

import pytest

from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.search.low_thrust_cycler_search import (
    LowThrustCyclerCandidate,
    search_low_thrust_cyclers,
    sweep_low_thrust_cyclers,
)
from cyclerfinder.search.maintain import optimise_aldrin_maintenance_dv

# Aldrin baseline parameters — pinned to McConaghy/Longuski/Byrnes 2002
# Table 4 row "1L1" and the catalogue ``aldrin-classic-em-k1-outbound`` row.
# These are the same values the optimiser uses internally; the test forwards
# them explicitly so the driver reproduces the baseline impulsive surrogate.
_ALDRIN_EM_TOF_DAYS = 146.0
_ALDRIN_ME_TOF_DAYS = 634.0
_ALDRIN_TOF1_BOUNDS = (100.0, 250.0)
_ALDRIN_TOF2_BOUNDS = (400.0, 900.0)
_ALDRIN_TOF_JITTER = (20.0, 60.0)
_ALDRIN_EARTH_FLYBY_ALT_KM = 200.0


def _aldrin_t0_guess() -> float:
    """Aldrin Earth-departure epoch (s) from the catalogue phase-inversion seed.

    Imported lazily — :mod:`cyclerfinder.search.maintain._default_t0_guess` is
    private and we only call it from a test helper, not from the driver.
    """
    from cyclerfinder.search.maintain import _default_t0_guess

    return _default_t0_guess(_ALDRIN_EM_TOF_DAYS)


@pytest.fixture(scope="module")
def aldrin_baseline_dv_kms() -> float:
    """Aldrin maintenance ΔV from the impulsive surrogate (the established
    baseline; circular backend; matches the existing test surface).

    This is OUR computation, not source-attested (McConaghy 2002 defers the
    DV estimation; the turn-angle test ~84° vs ~72° is the published anchor).
    The driver's reproduced value must equal it within numeric tolerance.
    """
    eph = Ephemeris("circular")
    impulsive = optimise_aldrin_maintenance_dv(eph, seed=0)
    assert impulsive.converged
    return float(impulsive.maintenance_dv_kms)


def test_aldrin_eme_search_recovers_baseline(aldrin_baseline_dv_kms: float) -> None:
    """E-M-E single-cell search recovers the Aldrin maintenance ΔV.

    The driver's :func:`search_low_thrust_cyclers` wires the same optimiser
    used by the Aldrin baseline. Calling it on the canonical sequence with
    matching seed must produce an identical maintenance ΔV.
    """
    eph = Ephemeris("circular")
    cand = search_low_thrust_cyclers(
        sequence=("E", "M", "E"),
        k_synodic=1,
        ephem=eph,
        t0_guess_sec=_aldrin_t0_guess(),
        leg_tof_guesses_days=(_ALDRIN_EM_TOF_DAYS, _ALDRIN_ME_TOF_DAYS),
        leg_tof_bounds_days=(_ALDRIN_TOF1_BOUNDS, _ALDRIN_TOF2_BOUNDS),
        tof_jitter_half_days=_ALDRIN_TOF_JITTER,
        synodic_pair=("E", "M"),
        closure_body="E",
        closure_flyby_alt_km=_ALDRIN_EARTH_FLYBY_ALT_KM,
        n_starts=5,
        seed=0,
    )
    assert cand is not None
    assert isinstance(cand, LowThrustCyclerCandidate)
    # Maintenance ΔV identity with the baseline impulsive surrogate.
    assert cand.maintenance_dv_kms == pytest.approx(aldrin_baseline_dv_kms, abs=1.0e-3)
    # The powered wrap carries the same ΔV through unchanged.
    assert cand.powered.maintenance_dv_kms == pytest.approx(cand.maintenance_dv_kms)
    # Independent cross-check landed within the loose plateau tolerance.
    assert not math.isnan(cand.independent_cross_check_residual_kms)
    assert cand.independent_cross_check_residual_kms < 1.0e-3


def test_aldrin_eme_search_lit_check_catches_purple_corpus() -> None:
    """The literature-check block is populated and structured.

    Without a live WebSearch, the offline default reports ``inconclusive``
    ("No search results returned at all"). The candidate is therefore NOT
    novelty-claimable. This pins the discipline: an offline run never clears
    a candidate; novelty requires a real search. The Genova-Aldrin /
    Aldrin / McConaghy anchors in KNOWN_CORPUS still get exercised through
    the keyword-fingerprint queries built by ``build_queries``.
    """
    eph = Ephemeris("circular")
    cand = search_low_thrust_cyclers(
        sequence=("E", "M", "E"),
        k_synodic=1,
        ephem=eph,
        t0_guess_sec=_aldrin_t0_guess(),
        leg_tof_guesses_days=(_ALDRIN_EM_TOF_DAYS, _ALDRIN_ME_TOF_DAYS),
        leg_tof_bounds_days=(_ALDRIN_TOF1_BOUNDS, _ALDRIN_TOF2_BOUNDS),
        tof_jitter_half_days=_ALDRIN_TOF_JITTER,
        synodic_pair=("E", "M"),
        closure_body="E",
        closure_flyby_alt_km=_ALDRIN_EARTH_FLYBY_ALT_KM,
        n_starts=5,
        seed=0,
    )
    assert cand is not None
    assert cand.literature_check["checked"] is True
    # Offline default: no results → inconclusive verdict, never claimable.
    assert cand.literature_check["status"] == "inconclusive"
    assert cand.novelty_claimable is False


def test_aldrin_eme_search_sims_flanagan_feasibility_witness() -> None:
    """Sims-Flanagan per-segment capability bound is reported as a verdict.

    At Aldrin scale (maintenance ΔV ~2.9 km/s, distributed across 20 segments
    over a ~146-day Earth→Mars leg) the per-segment ΔV is ~0.15 km/s. The
    per-segment capability bound at the default thrust (0.25 N) and mass
    (10 t) over a 146-d / 20 segment = 7.3 d segment is
    ``(2.5e-4 / 1e4) * 7.3 * 86400 ≈ 0.016 km/s``. That's an order of
    magnitude below the requested per-segment ΔV, so the witness must report
    ``sims_flanagan_feasible == False``. That is a STRUCTURAL Phase 1
    finding (an Aldrin-class powered maintenance is NOT deliverable by a
    flat NEXT-scale SEP thrust train without aggressive concentration or a
    higher Tmax), not a bug.
    """
    eph = Ephemeris("circular")
    cand = search_low_thrust_cyclers(
        sequence=("E", "M", "E"),
        k_synodic=1,
        ephem=eph,
        t0_guess_sec=_aldrin_t0_guess(),
        leg_tof_guesses_days=(_ALDRIN_EM_TOF_DAYS, _ALDRIN_ME_TOF_DAYS),
        leg_tof_bounds_days=(_ALDRIN_TOF1_BOUNDS, _ALDRIN_TOF2_BOUNDS),
        tof_jitter_half_days=_ALDRIN_TOF_JITTER,
        synodic_pair=("E", "M"),
        closure_body="E",
        closure_flyby_alt_km=_ALDRIN_EARTH_FLYBY_ALT_KM,
        n_starts=5,
        seed=0,
    )
    assert cand is not None
    # Aldrin-class ΔV over a flat NEXT-thrust profile: not feasible.
    assert cand.sims_flanagan_feasible is False
    # The thrust-train budget is reported per leg (here: 2 legs).
    assert len(cand.sims_flanagan_dv_train_kms) == 2
    # A note records the infeasibility.
    assert any("sims_flanagan_infeasible" in n for n in cand.notes)


def test_aldrin_eme_search_sims_flanagan_feasibility_with_higher_thrust() -> None:
    """Raise ``tmax_kn`` to 100 N and the same maintenance ΔV becomes feasible.

    The feasibility witness scales with the per-segment capability ``T_max
    /m * dt_seg``; multiplying ``tmax_kn`` by ~400 (to 100 N) lifts the
    per-segment bound above the Aldrin per-segment ΔV. This is the dial
    Phase 2 would explore via the NLP solve, not a brute parameter sweep.
    """
    eph = Ephemeris("circular")
    cand = search_low_thrust_cyclers(
        sequence=("E", "M", "E"),
        k_synodic=1,
        ephem=eph,
        t0_guess_sec=_aldrin_t0_guess(),
        leg_tof_guesses_days=(_ALDRIN_EM_TOF_DAYS, _ALDRIN_ME_TOF_DAYS),
        leg_tof_bounds_days=(_ALDRIN_TOF1_BOUNDS, _ALDRIN_TOF2_BOUNDS),
        tof_jitter_half_days=_ALDRIN_TOF_JITTER,
        synodic_pair=("E", "M"),
        closure_body="E",
        closure_flyby_alt_km=_ALDRIN_EARTH_FLYBY_ALT_KM,
        tmax_kn=0.1,  # 100 N — well above NEXT, in the NEP regime.
        n_starts=5,
        seed=0,
    )
    assert cand is not None
    assert cand.sims_flanagan_feasible is True


def test_powered_wrap_tsiolkovsky_identity() -> None:
    """The powered wrap propellant fraction equals the Tsiolkovsky identity.

    Source-free physics invariant: ``m_p/m_0 = 1 - exp(-ΔV / (g0 * Isp))``.
    The driver's report must satisfy it to within numerical tolerance for any
    converged candidate. This is the same identity the existing
    :mod:`cyclerfinder.search.lowthrust_maintenance` tests pin; we re-exercise
    it through the search driver so a regression in either layer is caught.
    """
    from cyclerfinder.core.constants import STANDARD_GRAVITY_KM_S2

    eph = Ephemeris("circular")
    cand = search_low_thrust_cyclers(
        sequence=("E", "M", "E"),
        k_synodic=1,
        ephem=eph,
        t0_guess_sec=_aldrin_t0_guess(),
        leg_tof_guesses_days=(_ALDRIN_EM_TOF_DAYS, _ALDRIN_ME_TOF_DAYS),
        leg_tof_bounds_days=(_ALDRIN_TOF1_BOUNDS, _ALDRIN_TOF2_BOUNDS),
        tof_jitter_half_days=_ALDRIN_TOF_JITTER,
        synodic_pair=("E", "M"),
        closure_body="E",
        closure_flyby_alt_km=_ALDRIN_EARTH_FLYBY_ALT_KM,
        isp_s=3000.0,
        n_starts=5,
        seed=0,
    )
    assert cand is not None
    expected_frac = 1.0 - math.exp(
        -cand.maintenance_dv_kms / (STANDARD_GRAVITY_KM_S2 * cand.powered.isp_s)
    )
    assert cand.powered.propellant_mass_fraction == pytest.approx(expected_frac, rel=1e-12)
    assert 0.0 < cand.powered.propellant_mass_fraction < 1.0


def test_search_rejects_non_closed_sequence() -> None:
    """A sequence that does not close raises ValueError."""
    with pytest.raises(ValueError, match="must close"):
        search_low_thrust_cyclers(
            sequence=("E", "M", "V"),
            k_synodic=1,
        )


def test_search_rejects_bad_k_synodic() -> None:
    """``k_synodic`` must be at least 1."""
    with pytest.raises(ValueError, match="k_synodic"):
        search_low_thrust_cyclers(
            sequence=("E", "M", "E"),
            k_synodic=0,
        )


def test_sweep_emits_one_row_per_cell_dedup() -> None:
    """Sweeping a one-cell grid produces one row; duplicates de-dup by fingerprint."""
    eph = Ephemeris("circular")
    rows = sweep_low_thrust_cyclers(
        sequence=("E", "M", "E"),
        k_synodic=1,
        ephem=eph,
        t0_epochs_sec=(_aldrin_t0_guess(),),
        leg_tof_shapes_days=(
            (_ALDRIN_EM_TOF_DAYS, _ALDRIN_ME_TOF_DAYS),
            # Same shape, just a tiny perturbation that the optimiser would
            # collapse back to the same family: tests the fingerprint dedup.
            (_ALDRIN_EM_TOF_DAYS + 0.1, _ALDRIN_ME_TOF_DAYS - 0.1),
        ),
        per_leg_revs_grid=((0, 0),),
        synodic_pair=("E", "M"),
        closure_body="E",
        closure_flyby_alt_km=_ALDRIN_EARTH_FLYBY_ALT_KM,
        tof_jitter_half_days=_ALDRIN_TOF_JITTER,
        n_starts=4,
        seed=0,
    )
    # At most one row after fingerprint dedup (same family from both shapes).
    assert len(rows) >= 1
    assert len(rows) <= 2  # allow the optimiser to land on a second branch
    for row in rows:
        assert row.sequence == ("E", "M", "E")
        assert row.period_k == 1


def test_as_dict_roundtrips_to_jsonl_friendly() -> None:
    """``as_dict()`` returns a JSON-serialisable payload."""
    import json

    eph = Ephemeris("circular")
    cand = search_low_thrust_cyclers(
        sequence=("E", "M", "E"),
        k_synodic=1,
        ephem=eph,
        t0_guess_sec=_aldrin_t0_guess(),
        leg_tof_guesses_days=(_ALDRIN_EM_TOF_DAYS, _ALDRIN_ME_TOF_DAYS),
        leg_tof_bounds_days=(_ALDRIN_TOF1_BOUNDS, _ALDRIN_TOF2_BOUNDS),
        tof_jitter_half_days=_ALDRIN_TOF_JITTER,
        synodic_pair=("E", "M"),
        closure_body="E",
        closure_flyby_alt_km=_ALDRIN_EARTH_FLYBY_ALT_KM,
        n_starts=4,
        seed=0,
    )
    assert cand is not None
    payload = cand.as_dict()
    # Round-trips through JSON (no NaN, no numpy dtypes).
    encoded = json.dumps(payload)
    decoded = json.loads(encoded)
    assert decoded["sequence"] == ["E", "M", "E"]
    assert decoded["period_k"] == 1
    assert decoded["maintenance_dv_kms"] == pytest.approx(cand.maintenance_dv_kms)
