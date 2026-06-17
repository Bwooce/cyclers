"""#347 Phase 2 P2.3 — smoke test for the sweep driver.

Tests the helper functions in :mod:`scripts.floquet_phase2_sweep` (the driver
module) that DON'T require a full sweep run. The full sweep is exercised by
the standalone driver via ``uv run python scripts/floquet_phase2_sweep.py``,
not by pytest — it takes ~10 minutes for the six-parent sweep.

Gates:

  1. The recovery helper ``_recover_parent_seed`` produces a SymmetricOrbit
     for the C32 anchor with period within 1% of the sourced 78.613 d.
  2. The recovery helper raises ValueError for an unknown label.
  3. The topology helper ``_topology_check`` returns (3, 2) for the C32
     anchor (independent reproduction of the Phase 1 P1.1 gate).
  4. The driver helper handles the `_down` suffix: C32_down recovers the
     same C32 IC but the driver picks direction=-1 (verified at the driver
     level — not by this smoke test).
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np

from cyclerfinder.genome.floquet_phase2_parents import SweepParent
from cyclerfinder.search.reachable_representatives import (
    C_J_BRAIK_ROSS,
    TU_DAYS,
    braik_ross_system,
)

# Load the driver module from scripts/ for direct helper-function access.
_REPO_ROOT = Path(__file__).resolve().parent.parent.parent
_DRIVER_PATH = _REPO_ROOT / "scripts" / "floquet_phase2_sweep.py"
_spec = importlib.util.spec_from_file_location("floquet_phase2_sweep_driver", _DRIVER_PATH)
assert _spec is not None and _spec.loader is not None
_driver = importlib.util.module_from_spec(_spec)
sys.modules["floquet_phase2_sweep_driver"] = _driver
_spec.loader.exec_module(_driver)


def _make_c32_parent() -> SweepParent:
    return SweepParent(
        label="C32",
        k1=3,
        k2=2,
        jacobi_anchor=C_J_BRAIK_ROSS,
        cj_window=(3.1294, 3.1544),
        dc=1e-4,
        n_steps=10,  # short for smoke test
        sourced_period_days=78.613,
        sourced_sigma_d=0.1583,
        notes="smoke test",
    )


def test_recover_parent_seed_c32_anchor() -> None:
    """The recovery helper produces a SymmetricOrbit at the C32 anchor."""
    parent = _make_c32_parent()
    seed = _driver._recover_parent_seed(parent)
    period_days = seed.period * TU_DAYS
    assert abs(period_days - 78.613) / 78.613 < 0.01, (
        f"recovered period {period_days:.3f} d differs from sourced 78.613 d by > 1%"
    )
    assert abs(seed.jacobi - C_J_BRAIK_ROSS) < 1e-6
    assert seed.converged


def test_recover_parent_seed_handles_down_suffix() -> None:
    """The `_down` suffix is stripped to find the base label C32."""
    parent = SweepParent(
        label="C32_down",
        k1=3,
        k2=2,
        jacobi_anchor=C_J_BRAIK_ROSS,
        cj_window=(3.10, 3.1294),
        dc=1e-4,
        n_steps=10,
        sourced_period_days=78.613,
        sourced_sigma_d=0.1583,
        notes="smoke test",
    )
    seed = _driver._recover_parent_seed(parent)
    period_days = seed.period * TU_DAYS
    # Same IC as the upward version.
    assert abs(period_days - 78.613) / 78.613 < 0.01


def test_recover_parent_seed_unknown_label_raises() -> None:
    """Unknown labels raise ValueError, not a silent default."""
    import pytest

    parent = SweepParent(
        label="UNKNOWN_FAMILY",
        k1=4,
        k2=3,
        jacobi_anchor=C_J_BRAIK_ROSS,
        cj_window=(3.0, 3.2),
        dc=1e-4,
        n_steps=10,
        sourced_period_days=80.0,
        sourced_sigma_d=0.5,
        notes="smoke test",
    )
    with pytest.raises(ValueError, match="unknown parent label"):
        _driver._recover_parent_seed(parent)


def test_topology_check_at_c32_anchor() -> None:
    """The C32 anchor returns winding topology (3, 2) (independent of Phase 1 P1.1)."""
    parent = _make_c32_parent()
    seed = _driver._recover_parent_seed(parent)
    state0 = np.array([seed.x0, 0.0, 0.0, 0.0, seed.ydot0, 0.0], dtype=np.float64)
    system = braik_ross_system()
    k1, k2, prograde = _driver._topology_check(system.mu, state0, seed.period)
    assert k1 == 3 and k2 == 2, f"C32 topology check returned ({k1}, {k2}) != (3, 2)"
    assert prograde is True


def test_topology_check_returns_sentinel_on_error() -> None:
    """A malformed input returns the (-1, -1, False) sentinel."""
    system = braik_ross_system()
    # Bad period (negative) triggers a ValueError inside winding_topology.
    k1, k2, prograde = _driver._topology_check(system.mu, np.zeros(6), -1.0)
    assert (k1, k2, prograde) == (-1, -1, False)
