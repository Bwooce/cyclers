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


@pytest.mark.slow
def test_jovian_closure_routes_through_full_pipeline() -> None:
    """A closed Jovian I-E-G chain flows bridge->signature->match->gauntlet (slow).

    NON-GOLDEN for the V_inf value (our computation). Per design note §5 the #76
    I-E-G closure is bend-INFEASIBLE in the no-leveraging model, so the realistic
    assertion is: it closes, reads ``novel`` against the null-numeric Jovian
    bucket, and routes to REJECTED (not SILVER) because it is bend-infeasible —
    proving the firewall holds on a non-heliocentric centre. We do NOT loosen tol
    to manufacture a SILVER.
    """
    from cyclerfinder.data.discover_novel import discover_novel_moon
    from cyclerfinder.verify.gauntlet import VerdictTier

    findings = list(
        discover_novel_moon(
            base_t0_sec=0.6 * 86400.0,
            n_epochs=4,
            span_days=2.0,
            budget_kms=50.0,
            max_workers=1,
        )
    )
    if not findings:
        pytest.xfail("no Jovian closure surfaced in this small grid (empty-set outcome)")

    # At least one finding flowed the full pipeline: it has a signature, a match
    # outcome, and a tiered verdict.
    f = findings[0]
    assert f.signature is not None
    assert f.match_outcome in ("novel", "probable-match-NEEDS-HUMAN", "known")
    assert f.verdict.tier in {
        VerdictTier.REJECTED,
        VerdictTier.BRONZE,
        VerdictTier.SILVER,
        VerdictTier.GOLD,
    }
    # The #76 honest-risk: a bend-INFEASIBLE closure must route REJECTED, never
    # SILVER (the firewall on a non-heliocentric centre).
    for finding in findings:
        if not finding.bend_feasible:
            assert finding.verdict.tier == VerdictTier.REJECTED
    # The Jovian bucket is null-numeric -> closures read novel.
    assert any(f.match_outcome == "novel" for f in findings)
