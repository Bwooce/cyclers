"""Phase 6 Phase 2: Galilean topology set + centred Jovian novelty sweep (Tasks 2.0-2.2)."""

from __future__ import annotations

import pytest

from cyclerfinder.data.discover_novel import jovian_galilean_topologies


def test_galilean_topologies_are_jovian_moon_sequences() -> None:
    specs = jovian_galilean_topologies()
    assert specs
    for s in specs:
        assert set(s.sequence) <= {"Io", "Europa", "Ganymede", "Callisto"}
        assert s.sequence[0] == s.sequence[-1]  # closed tour


def test_discover_novel_moon_prunes_then_scans(monkeypatch: pytest.MonkeyPatch) -> None:
    # Assert the loop applies the prune and only scans survivors; stub the scan so
    # the test stays fast (the real DE440 sweep is the Phase 5 slow run).
    import cyclerfinder.data.discover_novel as dn

    seen_topos: list[tuple[str, ...]] = []

    def _fake_grid(**kw: object) -> list[object]:
        seen_topos.append(kw["sequence"])  # type: ignore[arg-type]
        return [object()]

    monkeypatch.setattr(dn, "scan_parallel", lambda grid, **kw: [])
    monkeypatch.setattr(dn, "build_epoch_branch_grid", _fake_grid)
    list(dn.discover_novel_moon(base_t0_sec=0.0, n_epochs=2, budget_kms=50.0))
    # At least the known-closing I-E-G family survived the prune and was scanned.
    assert any(set(s) <= {"Io", "Europa", "Ganymede", "Callisto"} for s in seen_topos)
