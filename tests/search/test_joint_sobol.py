"""#318 Phase 2b — Sobol sampler tests + mandatory positive-control gate.

Kernel-free tests verify :func:`make_sobol_cells` structure and reproducibility.
The kernel-backed positive-control test proves the pipeline doesn't silently
discard the #223-validated Liang Member D closure
([[feedback_verify_gauntlet_with_positive_control]]).
"""

from __future__ import annotations

import os
from typing import Any

import pytest

from cyclerfinder.nbody import jovian
from cyclerfinder.search.joint_cell import evaluate_joint_cell, liang_member_d_cell
from cyclerfinder.search.joint_sobol import make_sobol_cells


def _resolve_kernel() -> str | None:
    env = jovian.jup365_kernel_path()
    if env is not None:
        return env
    fallback = os.path.expanduser("~/dev/references/kernels/jup365.bsp")
    return fallback if os.path.exists(fallback) else None


_KERNEL = _resolve_kernel()


# --- kernel-free -------------------------------------------------------


def test_make_sobol_cells_count() -> None:
    cells = make_sobol_cells(
        n_samples=8,
        sequence=list(jovian.CGCEC),
        epoch_window=("2033-01-01", "2035-01-01"),
        n_revs_range=(1, 2),
        tof_seed_range=(15.0, 40.0),
    )
    assert len(cells) == 8


def test_make_sobol_cells_structure() -> None:
    cells = make_sobol_cells(
        n_samples=16,
        sequence=list(jovian.CGCEC),
        epoch_window=("2033-01-01", "2035-01-01"),
        n_revs_range=(1, 2),
        tof_seed_range=(15.0, 40.0),
    )
    for cell in cells:
        assert cell.sequence == jovian.CGCEC
        assert cell.sequence[0] == cell.sequence[-1]
        assert len(cell.n_revs) == 4
        assert len(cell.branches) == 4
        assert len(cell.tof_seed_days) == 4
        assert all(1 <= r <= 2 for r in cell.n_revs), cell.n_revs
        assert all(b in ("high", "low") for b in cell.branches), cell.branches
        assert all(15.0 <= t <= 40.0 for t in cell.tof_seed_days), cell.tof_seed_days


def test_make_sobol_cells_epoch_range() -> None:
    from astropy.time import Time

    cells = make_sobol_cells(
        n_samples=8,
        sequence=list(jovian.CGCEC),
        epoch_window=("2033-01-01", "2035-01-01"),
        n_revs_range=(1, 1),
        tof_seed_range=(20.0, 35.0),
    )
    jd_lo = float(Time("2033-01-01", scale="tdb").jd)
    jd_hi = float(Time("2035-01-01", scale="tdb").jd)
    for cell in cells:
        jd = float(Time(cell.epoch_iso, scale="tdb").jd)
        assert jd_lo <= jd <= jd_hi + 1e-6, f"epoch {cell.epoch_iso!r} out of window"


def test_make_sobol_cells_rejects_open_sequence() -> None:
    with pytest.raises(ValueError, match="closed"):
        make_sobol_cells(
            n_samples=4,
            sequence=["Callisto", "Ganymede", "Europa"],  # not closed
            epoch_window=("2033-01-01", "2035-01-01"),
            n_revs_range=(1, 1),
            tof_seed_range=(20.0, 30.0),
        )


def test_make_sobol_cells_reproducible() -> None:
    """Same seed gives identical cells."""
    kw: dict[str, Any] = dict(
        sequence=list(jovian.CGCEC),
        epoch_window=("2033-01-01", "2035-01-01"),
        n_revs_range=(1, 1),
        tof_seed_range=(20.0, 35.0),
        seed=42,
    )
    a = make_sobol_cells(8, **kw)
    b = make_sobol_cells(8, **kw)
    assert [c.epoch_iso for c in a] == [c.epoch_iso for c in b]
    assert [c.n_revs for c in a] == [c.n_revs for c in b]
    assert [c.branches for c in a] == [c.branches for c in b]


def test_make_sobol_cells_different_seeds() -> None:
    """Different seeds give different epochs."""
    kw: dict[str, Any] = dict(
        sequence=list(jovian.CGCEC),
        epoch_window=("2033-01-01", "2035-01-01"),
        n_revs_range=(1, 1),
        tof_seed_range=(20.0, 35.0),
    )
    a = make_sobol_cells(8, seed=0, **kw)
    b = make_sobol_cells(8, seed=1, **kw)
    # Epochs should differ (different scramble)
    assert [c.epoch_iso for c in a] != [c.epoch_iso for c in b]


# --- kernel-backed positive control ------------------------------------


@pytest.mark.skipif(_KERNEL is None, reason="needs JUP365 kernel")
def test_sobol_pipeline_positive_control() -> None:
    """POSITIVE CONTROL: Liang Member D must survive the Sobol-pipeline prefilter.

    Proves the pipeline doesn't silently discard the #223-validated near-ballistic
    CGCEC closure ([[feedback_verify_gauntlet_with_positive_control]]). If this test
    fails, the Sobol sweep has no positive basis and must not be trusted.
    """
    assert _KERNEL is not None
    ephem = jovian.JovianEphemeris(_KERNEL)
    cell = liang_member_d_cell()
    v = evaluate_joint_cell(cell, ephem)

    assert v.feasible, (
        f"POSITIVE CONTROL FAILED: Liang Member D failed prefilter — "
        f"closure_defect={v.closure_defect_ms:.3f} m/s, min_alt={v.min_alt_km:.1f} km"
    )
    assert v.closure_defect_ms < 1.0, f"closure_defect_ms={v.closure_defect_ms}"
    assert v.min_alt_km > 100.0, f"min_alt_km={v.min_alt_km}"
    assert v.cycle_tof_days > 99.0, f"cycle_tof_days={v.cycle_tof_days}"
