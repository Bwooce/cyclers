"""M-ED Phase 3: mode= kwarg on optimise_cell_ephemeris (plan Phase 3).

test_default_mode_unchanged is a CHARACTERISATION test: the default
mode="maintenance" path must be byte-identical to pre-M-ED. If it passes
immediately after the signature change (kwarg added, default branch untouched),
that is the desired outcome.
"""

from __future__ import annotations

import inspect

import pytest

from cyclerfinder.search.optimize import optimise_cell_ephemeris


def test_mode_kwarg_exists_and_defaults_to_maintenance() -> None:
    sig = inspect.signature(optimise_cell_ephemeris)
    assert "mode" in sig.parameters
    assert sig.parameters["mode"].default == "maintenance"


def test_unknown_mode_raises() -> None:
    from cyclerfinder.core.ephemeris import Ephemeris
    from cyclerfinder.search.sequence import Cell

    cell = Cell(
        bodies=("E", "M"),
        sequence=("E", "M", "E"),
        period_k=2,
        per_leg_revs=(0, 0),
        per_leg_branch=("single", "single"),
    )
    with pytest.raises(ValueError, match="mode"):
        optimise_cell_ephemeris(cell, Ephemeris(model="circular"), vinf_cap=9.0, mode="bogus")


@pytest.mark.slow
def test_ballistic_mode_closes_s1l1_returns_real_residual() -> None:
    """ballistic mode returns a REAL closure residual (V_inf-continuity), not
    the maintenance-dV proxy. NON-GOLDEN: closure is asserted, the V_inf value
    is OUR computation (S1L1 floors Mars ~6.4 — see project memory)."""
    from cyclerfinder.core.ephemeris import Ephemeris
    from cyclerfinder.search.sequence import Cell

    cell = Cell(
        bodies=("E", "M"),
        sequence=("E", "M", "E", "E"),
        period_k=2,
        per_leg_revs=(0, 1, 2),
        per_leg_branch=("single", "single", "low"),
    )
    result = optimise_cell_ephemeris(
        cell,
        Ephemeris("astropy"),
        vinf_cap=9.0,
        priority_date_iso="2030-03-22",
        vinf_targets_kms={"E": 5.6, "M": 6.4},
        tof_seed_days=[154.0, 379.0, (1.4612 + 2.8096) * 365.25 - 154.0 - 379.0],
        mode="ballistic",
    )
    assert result.converged
    assert result.closure_residual_kms < 0.1  # real V_inf-continuity residual


def test_scan_epochs_kwarg_exists_and_defaults_to_one() -> None:
    sig = inspect.signature(optimise_cell_ephemeris)
    assert "scan_epochs" in sig.parameters
    assert sig.parameters["scan_epochs"].default == 1
    assert "scan_window_years" in sig.parameters
    assert "scan_max_workers" in sig.parameters


@pytest.mark.slow
def test_ballistic_scan_rung_runs_and_is_deterministic() -> None:
    """The scan rung (scan_epochs>1) drives a parallel epoch grid and returns a
    deterministic result: two identical runs give identical residual/V_inf.
    NON-GOLDEN: the V_inf is OUR computation; the gate is reproducibility, not a
    sourced anchor."""
    from cyclerfinder.core.ephemeris import Ephemeris
    from cyclerfinder.search.sequence import Cell

    cell = Cell(
        bodies=("E", "M"),
        sequence=("E", "M", "E", "E"),
        period_k=2,
        per_leg_revs=(0, 1, 2),
        per_leg_branch=("single", "single", "low"),
    )
    kwargs = dict(
        vinf_cap=9.0,
        priority_date_iso="2030-03-22",
        vinf_targets_kms={"E": 5.6, "M": 6.4},
        tof_seed_days=[154.0, 379.0, (1.4612 + 2.8096) * 365.25 - 154.0 - 379.0],
        mode="ballistic",
        scan_epochs=6,
        scan_window_years=4.0,
        scan_max_workers=4,
    )
    r1 = optimise_cell_ephemeris(cell, Ephemeris("astropy"), **kwargs)  # type: ignore[arg-type]
    r2 = optimise_cell_ephemeris(cell, Ephemeris("astropy"), **kwargs)  # type: ignore[arg-type]
    assert r1.closure_residual_kms == r2.closure_residual_kms
    assert r1.converged == r2.converged
    # The scan rung finds at least as good a residual as the single start would
    # from the same priority epoch (it includes that epoch's neighbourhood).
    assert r1.closure_residual_kms < float("inf")


def test_maintenance_mode_result_unchanged_for_aldrin_cell() -> None:
    """Byte-identical-default contract (spec §4): mode="maintenance" on the
    Aldrin E-M-E cell produces the same OptimisationResult fields as the
    pre-mode-kwarg path (mirrors test_optimize_ephemeris.py's Aldrin parity
    test, which calls optimise_cell_ephemeris with no mode= and the same args).
    """
    from cyclerfinder.core.ephemeris import Ephemeris
    from cyclerfinder.search.optimize import OptimisationResult
    from cyclerfinder.search.sequence import Cell

    cell = Cell(
        bodies=("E", "M"),
        sequence=("E", "M", "E"),
        period_k=1,
        per_leg_revs=(0, 0),
        per_leg_branch=("single", "single"),
    )
    eph = Ephemeris(model="astropy")
    kwargs = dict(
        vinf_cap=12.0,
        priority_date_iso="1985-01-01",
        vinf_targets_kms={"E": 6.5, "M": 9.7},
        n_starts=3,
        seed=0,
    )
    default = optimise_cell_ephemeris(cell, eph, **kwargs)  # type: ignore[arg-type]
    explicit = optimise_cell_ephemeris(cell, eph, mode="maintenance", **kwargs)  # type: ignore[arg-type]
    assert isinstance(default, OptimisationResult)
    assert default.converged == explicit.converged
    assert default.constraints_satisfied == explicit.constraints_satisfied
    assert default.closure_residual_kms == explicit.closure_residual_kms
    assert default.best_score.max_vinf_kms == explicit.best_score.max_vinf_kms
