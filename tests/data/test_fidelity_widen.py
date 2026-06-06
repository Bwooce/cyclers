"""M-3D Phase 4: circular-inclined Fidelity rung registered (plan §4)."""

from __future__ import annotations

from cyclerfinder.data.provenance import is_fidelity


def test_circular_inclined_is_a_known_fidelity() -> None:
    assert is_fidelity("circular-inclined")


def test_existing_fidelities_unchanged() -> None:
    for tier in ("circular-coplanar", "analytic-ephemeris", "real-de440"):
        assert is_fidelity(tier)
    assert not is_fidelity("bogus-rung")
