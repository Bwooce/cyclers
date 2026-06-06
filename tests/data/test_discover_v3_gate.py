"""M-ED Phase 5: discover V3 branch is a real ballistic-closure gate, not a
stub-raise assumption (plan Phase 5 Task 5.0; spec §0 finding 3)."""

from __future__ import annotations

import inspect

from cyclerfinder.data import discover as discover_mod


def test_v3_branch_no_longer_assumes_notimplemented() -> None:
    src = inspect.getsource(discover_mod._auto_validate)
    # The stale comment / assumption must be gone.
    assert "raises until M6b lands" not in src
    assert 'mode="ballistic"' in src or "mode='ballistic'" in src
