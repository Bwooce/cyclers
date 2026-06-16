"""Tests for #318 Phase 1 — multi-axis joint-search substrate.

Discipline:

* **Sourced golden anchor**: catalogue row ``aldrin-classic-em-k1-outbound`` is
  the regression target. At the joint-zero corner (powered=0, n_revs=0, z0=0,
  epoch=None) the joint driver must re-find the Aldrin cycler — same V_inf
  tuple, same maintenance ΔV as the established baseline.
* **Composition, not rewrite**: the joint driver wires the four existing axis
  modules and pins their already-tested behaviour. The tests are sanity gates
  on the COMPOSITION, not redefinitions of the axis verdicts (those have their
  own test files).
* **No catalogue writeback** — Phase 1 substrate only.
* **No novelty claims** — frame: "joint-axis Phase 1 substrate + small probe".
"""

from __future__ import annotations

import math
from typing import Any

import pytest

from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.search.low_thrust_cycler_search import (
    LowThrustCyclerCandidate,
)
from cyclerfinder.search.maintain import (
    _default_t0_guess,
    optimise_aldrin_maintenance_dv,
)
from cyclerfinder.search.multi_axis_search import (
    joint_axis_search,
)

# Aldrin params shared with tests/search/test_low_thrust_cycler_search.py.
_ALDRIN_EM_TOF_DAYS = 146.0
_ALDRIN_ME_TOF_DAYS = 634.0
_ALDRIN_TOF1_BOUNDS = (100.0, 250.0)
_ALDRIN_TOF2_BOUNDS = (400.0, 900.0)
_ALDRIN_TOF_JITTER = (20.0, 60.0)
_ALDRIN_EARTH_FLYBY_ALT_KM = 200.0


def _aldrin_joint_args() -> dict[str, Any]:
    """Canonical Aldrin-corner search args (epoch-blind, 0 powered, 0 revs, 0 z0)."""
    return {
        "primary_sequence": ("E", "M", "E"),
        "k_synodic": 1,
        "ephem": Ephemeris("circular"),
        "powered_budgets_kms": (0.0,),
        "n_revs_grid_per_leg": ((0,), (0,)),
        "z0_amplitudes_nondim": (0.0,),
        "launch_epoch_grid": (None,),
        "leg_tof_guesses_days": (_ALDRIN_EM_TOF_DAYS, _ALDRIN_ME_TOF_DAYS),
        "leg_tof_bounds_days": (_ALDRIN_TOF1_BOUNDS, _ALDRIN_TOF2_BOUNDS),
        "tof_jitter_half_days": _ALDRIN_TOF_JITTER,
        "synodic_pair": ("E", "M"),
        "closure_body": "E",
        "closure_flyby_alt_km": _ALDRIN_EARTH_FLYBY_ALT_KM,
        "n_starts": 5,
        "seed": 0,
        "t0_guess_sec_for_epoch_blind": _default_t0_guess(_ALDRIN_EM_TOF_DAYS),
    }


@pytest.fixture(scope="module")
def aldrin_baseline_dv_kms() -> float:
    """Established Aldrin baseline (impulsive turn-deficit surrogate)."""
    impulsive = optimise_aldrin_maintenance_dv(Ephemeris("circular"), seed=0)
    assert impulsive.converged
    return float(impulsive.maintenance_dv_kms)


# ---------------------------------------------------------------------------
# Test 1 — Joint-zero corner reduces to the existing #309 ballistic search
# ---------------------------------------------------------------------------


def test_joint_zero_corner_reduces_to_309_aldrin(aldrin_baseline_dv_kms: float) -> None:
    """At (powered=0, n_revs=0, z0=0, epoch=None) the joint driver must re-find
    the existing #309 result on the Aldrin tour bit-identically.

    This is the COMPOSITION regression: the joint driver MUST NOT change the
    behaviour of the ballistic + single-rev + planar + epoch-blind cell.
    """
    args = _aldrin_joint_args()
    candidates = list(joint_axis_search(**args))
    assert len(candidates) == 1
    cand = candidates[0]
    # Identity: same maintenance ΔV as the baseline (the #309 driver IS the
    # baseline in this corner, modulo the joint-axis tags).
    assert cand.powered_maintenance_dv_kms_per_synodic == pytest.approx(
        aldrin_baseline_dv_kms, abs=1.0e-3
    )
    # Axis coords echo the cell.
    assert cand.powered_budget_kms_requested == 0.0
    assert cand.n_revs_per_leg == (0, 0)
    assert cand.z0_amplitude_nondim == 0.0
    assert cand.launch_epoch_utc is None
    # The powered-axis record was attached.
    assert isinstance(cand.powered_axis_record, LowThrustCyclerCandidate)


# ---------------------------------------------------------------------------
# Test 2 — Aldrin reproduction at the joint-zero corner
# ---------------------------------------------------------------------------


def test_aldrin_reproduction_at_joint_zero_corner() -> None:
    """The joint driver re-finds the Aldrin V_inf tuple at the joint-zero corner.

    The published Aldrin V_inf values (Russell 2004 Table 3.4, cycler 1.0.1.-1)
    are 6.5 km/s at Earth and 9.7 km/s at Mars; the circular-coplanar
    optimum produced by the existing baseline lands near 6.52 and 9.73.
    """
    args = _aldrin_joint_args()
    candidates = list(joint_axis_search(**args))
    assert len(candidates) == 1
    cand = candidates[0]
    vinf = cand.vinf_tuple_kms
    assert len(vinf) == 3
    # Earth V_inf (first body and closing body) at the circular-coplanar Aldrin
    # optimum is ~6.52 km/s. Published is 6.5 km/s; tolerance 0.05 km/s.
    assert vinf[0] == pytest.approx(6.52, abs=0.05)
    assert vinf[-1] == pytest.approx(6.52, abs=0.05)
    # Mars V_inf at the same optimum is ~9.73 km/s; published 9.7 km/s.
    assert vinf[1] == pytest.approx(9.73, abs=0.05)


# ---------------------------------------------------------------------------
# Test 3 — Multi-rev grid expansion never reduces the surviving population
# ---------------------------------------------------------------------------


def test_multi_rev_grid_does_not_reduce_survivor_count() -> None:
    """At (powered=0, multi-rev grid, z0=0, epoch=None), the survivor count
    must be >= the direct-revs cell count.

    Multi-rev adds *additional* Lambert branches; it never removes the direct
    branch. So the joint sweep with revs=(0,1) per leg must produce >= the
    same candidates as revs=(0,) per leg.
    """
    args = _aldrin_joint_args()
    # Direct (revs=(0,) per leg).
    direct = list(joint_axis_search(**args))
    n_direct = len(direct)

    # Multi-rev (revs=(0, 1) per leg). The optimiser may or may not converge
    # on the rev=1 branches; the cell either produces a row or it doesn't.
    # What's pinned here: multi-rev never PRUNES the direct-cell row.
    args_mr = dict(args)
    args_mr["n_revs_grid_per_leg"] = ((0, 1), (0, 1))
    mr = list(joint_axis_search(**args_mr))
    n_mr = len(mr)
    assert n_mr >= n_direct
    # And the direct-cell row is among the multi-rev survivors (the optimiser
    # is deterministic at the same seed; the (0, 0) cell survives identically).
    direct_revs = [c for c in mr if c.n_revs_per_leg == (0, 0)]
    assert len(direct_revs) == 1
    assert direct_revs[0].powered_maintenance_dv_kms_per_synodic == pytest.approx(
        direct[0].powered_maintenance_dv_kms_per_synodic, abs=1.0e-3
    )


# ---------------------------------------------------------------------------
# Test 4 — 3D extension is recorded, not yet solved (Phase 1 contract)
# ---------------------------------------------------------------------------


def test_3d_extension_records_z0_amplitude_on_candidate() -> None:
    """At a non-zero z0 cell, the joint driver MUST record the amplitude on
    the candidate and tag the Phase-2 follow-up hook in notes.

    Phase 1 contract: the 2D Lambert engine is the production solve; the
    Axis-C amplitude is RECORDED for a Phase 2 follow-up that pipes the
    candidate through the #291 3D corrector.
    """
    args = _aldrin_joint_args()
    args["z0_amplitudes_nondim"] = (0.0, 1.0e-2)
    candidates = list(joint_axis_search(**args))
    # Both z0 cells should survive (the cell only differs in the Axis-C tag).
    assert len(candidates) == 2
    by_z0 = {c.z0_amplitude_nondim: c for c in candidates}
    assert 0.0 in by_z0
    assert pytest.approx(1.0e-2, abs=1e-12) == max(by_z0.keys())
    z3 = by_z0[max(by_z0.keys())]
    # Phase 1 follow-up hook tag in notes.
    assert any("axis_C_3d_request_recorded" in n for n in z3.notes)
    # The cycler verdict is the SAME on both cells (Phase 1 doesn't yet drive
    # the 3D corrector; the per-cell tag is the only difference).
    assert z3.powered_maintenance_dv_kms_per_synodic == pytest.approx(
        by_z0[0.0].powered_maintenance_dv_kms_per_synodic, abs=1e-6
    )


# ---------------------------------------------------------------------------
# Test 5 — Non-closed sequence is rejected
# ---------------------------------------------------------------------------


def test_non_closed_sequence_is_rejected() -> None:
    """A primary_sequence whose first != last body must raise."""
    with pytest.raises(ValueError, match="must close"):
        list(
            joint_axis_search(
                primary_sequence=("E", "M"),  # doesn't close
                k_synodic=1,
                ephem=Ephemeris("circular"),
                powered_budgets_kms=(0.0,),
                n_revs_grid_per_leg=((0,),),
                z0_amplitudes_nondim=(0.0,),
                launch_epoch_grid=(None,),
            )
        )


# ---------------------------------------------------------------------------
# Test 6 — k_synodic validation
# ---------------------------------------------------------------------------


def test_bad_k_synodic_is_rejected() -> None:
    """k_synodic < 1 must raise."""
    with pytest.raises(ValueError, match="k_synodic"):
        list(
            joint_axis_search(
                primary_sequence=("E", "M", "E"),
                k_synodic=0,
                ephem=Ephemeris("circular"),
            )
        )


# ---------------------------------------------------------------------------
# Test 7 — n_revs grid shape validation
# ---------------------------------------------------------------------------


def test_n_revs_grid_shape_mismatch_is_rejected() -> None:
    """A per-leg revs grid whose length doesn't match n_legs must raise."""
    with pytest.raises(ValueError, match="n_revs_grid_per_leg"):
        list(
            joint_axis_search(
                primary_sequence=("E", "M", "E"),  # n_legs = 2
                k_synodic=1,
                ephem=Ephemeris("circular"),
                n_revs_grid_per_leg=((0,),),  # length 1, not 2
            )
        )


# ---------------------------------------------------------------------------
# Test 8 — Candidate serialisation round-trips to JSONL-friendly dict
# ---------------------------------------------------------------------------


def test_candidate_as_dict_jsonl_friendly() -> None:
    """JointAxisCandidate.as_dict() must produce JSON-encodable primitives only."""
    import json

    args = _aldrin_joint_args()
    candidates = list(joint_axis_search(**args))
    assert len(candidates) == 1
    d = candidates[0].as_dict()
    # NaN -> None; tuples -> lists; floats / strs / bools / ints OK.
    json_text = json.dumps(d)  # raises if non-JSON-encodable
    assert isinstance(json_text, str)
    # The nested powered_axis_record dict is also serialised.
    assert isinstance(d["powered_axis_record"], dict)
    # No NaN sneaking through.
    assert all(not (isinstance(v, float) and math.isnan(v)) for v in d.values())


# ---------------------------------------------------------------------------
# Test 9 — Cartesian-product cell count (joint-axis enumeration)
# ---------------------------------------------------------------------------


def test_cartesian_product_cell_count() -> None:
    """The Cartesian product over (powered x revs x z0 x epoch) at the joint-
    zero corner enumerates all cells. With grids of size (1, 1*1, 2, 1) the
    driver must emit exactly 2 candidates (the only varying axis is z0)."""
    args = _aldrin_joint_args()
    args["z0_amplitudes_nondim"] = (0.0, 5e-3)
    candidates = list(joint_axis_search(**args))
    # 1 powered_budget x 1 revs tuple x 2 z0 x 1 epoch = 2 cells; each one
    # converges (same powered/n_revs/epoch, only the axis-C tag differs).
    assert len(candidates) == 2


# ---------------------------------------------------------------------------
# Test 10 — Powered budget annotation tracks the optimiser verdict
# ---------------------------------------------------------------------------


def test_powered_budget_above_optimum_no_advisory(aldrin_baseline_dv_kms: float) -> None:
    """If the cell's powered_budget_kms is ABOVE the optimum, no advisory note
    about budget exceedance is emitted. The Aldrin maintenance ΔV is ~1.29
    km/s; a budget of 5.0 km/s is comfortably above it."""
    args = _aldrin_joint_args()
    args["powered_budgets_kms"] = (5.0,)  # way above ~1.29 km/s optimum
    candidates = list(joint_axis_search(**args))
    assert len(candidates) == 1
    cand = candidates[0]
    assert cand.powered_budget_kms_requested == 5.0
    assert cand.powered_maintenance_dv_kms_per_synodic == pytest.approx(
        aldrin_baseline_dv_kms, abs=1.0e-3
    )
    # Optimum is below budget, so no exceedance note.
    assert not any("maintenance_dv_exceeds_cell_budget" in n for n in cand.notes)


def test_powered_budget_below_optimum_records_advisory() -> None:
    """If the cell's powered_budget_kms is BELOW the optimum, the cell is
    still emitted (Phase 1 records, does not gate) and notes call out the
    exceedance. The Aldrin maintenance ΔV is ~1.29 km/s; a budget of
    0.05 km/s is well below it."""
    args = _aldrin_joint_args()
    args["powered_budgets_kms"] = (0.05,)  # below ~1.29 km/s
    candidates = list(joint_axis_search(**args))
    assert len(candidates) == 1
    cand = candidates[0]
    assert cand.powered_budget_kms_requested == 0.05
    assert any("maintenance_dv_exceeds_cell_budget" in n for n in cand.notes)
