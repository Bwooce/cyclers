"""#318 Phase 2a — the real-ephemeris joint-cell evaluator + its positive control.

Kernel-free: cell construction + argument validation. Kernel-backed (skipped when JUP365
is absent, as CI has no 1.14 GB kernel): the positive control — the evaluator reproduces
the #223-validated Liang Member D CGCEC near-ballistic closure as a joint-cell and records
the four axis coordinates, including the real-ephemeris out-of-plane (broken-plane) extent
the Phase-1 circular-coplanar substrate could not represent.
"""

from __future__ import annotations

import os

import pytest

from cyclerfinder.nbody import jovian
from cyclerfinder.search.joint_cell import (
    JointCell,
    JointCellVerdict,
    evaluate_joint_cell,
    liang_member_d_cell,
)


def _resolve_kernel() -> str | None:
    env = jovian.jup365_kernel_path()
    if env is not None:
        return env
    fallback = os.path.expanduser("~/dev/references/kernels/jup365.bsp")
    return fallback if os.path.exists(fallback) else None


_KERNEL = _resolve_kernel()
_needs_kernel = pytest.mark.skipif(_KERNEL is None, reason="JUP365 kernel not furnished")


# --- kernel-free ---------------------------------------------------------------


def test_liang_cell_builder_shape() -> None:
    """The Liang positive-control cell is a well-formed closed 4-leg CGCEC cell."""
    cell = liang_member_d_cell()
    assert cell.sequence[0] == cell.sequence[-1]  # closed
    n_legs = len(cell.sequence) - 1
    assert n_legs == 4
    assert len(cell.n_revs) == len(cell.branches) == len(cell.tof_seed_days) == n_legs
    assert cell.primary == "Jupiter"
    assert cell.n_revs == (1, 1, 1, 1)


def test_evaluate_rejects_non_closed_sequence() -> None:
    """A non-closed sequence fails loud (never silently mis-evaluated)."""
    bad = JointCell(
        primary="Jupiter",
        sequence=("Callisto", "Ganymede", "Europa"),  # first != last
        epoch_iso="2033-09-25T00:00:00",
        n_revs=(1, 1),
        branches=("high", "low"),
        tof_seed_days=(30.0, 20.0),
    )
    with pytest.raises(ValueError, match="closed cycle"):
        evaluate_joint_cell(bad, ephem=None)  # type: ignore[arg-type]


def test_evaluate_rejects_axis_length_mismatch() -> None:
    """The Axis-B grids and ToF seed must match the leg count."""
    bad = JointCell(
        primary="Jupiter",
        sequence=("Callisto", "Ganymede", "Callisto"),  # 2 legs
        epoch_iso="2033-09-25T00:00:00",
        n_revs=(1,),  # only 1
        branches=("high", "low"),
        tof_seed_days=(30.0, 20.0),
    )
    with pytest.raises(ValueError, match="must each have 2 entries"):
        evaluate_joint_cell(bad, ephem=None)  # type: ignore[arg-type]


# --- kernel-backed positive control --------------------------------------------


@_needs_kernel
def test_joint_cell_reproduces_liang_member_d() -> None:
    """POSITIVE CONTROL: the joint-cell evaluator reproduces Liang Member D's closure.

    The evaluator must re-find the #223-validated near-ballistic CGCEC closure (sub-m/s
    cycle-1 defect, all flybys above the 100 km floor) AND record the four axis
    coordinates. Axis C (real out-of-plane extent) must be NON-ZERO — the broken-plane
    geometry the Phase-1 coplanar substrate discarded.
    """
    assert _KERNEL is not None
    ephem = jovian.JovianEphemeris(_KERNEL)
    v = evaluate_joint_cell(liang_member_d_cell(), ephem)
    assert isinstance(v, JointCellVerdict)

    # Near-ballistic closure (the #223 signature), feasible, ~100 d cycle.
    assert v.closure_defect_ms < 1.0, v.closure_defect_ms
    assert v.feasible, (v.min_alt_km, v.altitudes_km)
    assert 99.0 < v.cycle_tof_days < 101.0, v.cycle_tof_days
    assert v.min_alt_km > 100.0, v.min_alt_km

    # The four axis coordinates are recorded; Axis C is genuinely 3D (real ephemeris).
    assert v.axis_a_powered_dv_ms < 1.0
    assert v.axis_b_n_revs == (1, 1, 1, 1)
    assert v.axis_c_max_abs_z_km > 1.0e4, v.axis_c_max_abs_z_km  # real out-of-plane
    assert v.axis_d_epoch_iso == "2033-09-25T18:04:43"
    # CGCEC has 5 flyby nodes (Callisto x3, Ganymede, Europa) → 5 V∞ entries.
    assert len(v.vinf_kms) == 5
