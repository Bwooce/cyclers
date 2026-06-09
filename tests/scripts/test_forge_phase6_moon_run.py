"""Phase 6 Phase 4: the Jovian moon-run orchestrator (plan Task 4.2, fast/stubbed).

Asserts (with a stubbed scan): (a) a barren sweep writes a valid method-versioned
empty-region report, and (b) re-running with the *same* method against that record
skips (should_sweep -> False). No real DE440 corrector runs here — the Phase 5
slow run exercises the actual sweep.
"""

from __future__ import annotations

from pathlib import Path

import scripts.forge_phase6_moon_run as run
from cyclerfinder.data.empty_regions import load_empty_regions_list


def test_barren_sweep_writes_valid_empty_region(tmp_path: Path, monkeypatch: object) -> None:
    # Stub the sweep to yield NOTHING (the EMPTY outcome — the expected base rate).
    monkeypatch.setattr(run, "discover_novel_moon", lambda **kw: iter(()))  # type: ignore[attr-defined]
    empty_path = tmp_path / "empty_regions.jsonl"
    queue_path = tmp_path / "review_queue.jsonl"

    rc = run.main(
        [
            "--epochs",
            "2",
            "--workers",
            "1",
            "--region-id",
            "test-region",
            "--empty-regions",
            str(empty_path),
            "--queue",
            str(queue_path),
        ]
    )
    assert rc == 0

    records = load_empty_regions_list(empty_path)
    assert len(records) == 1
    rec = records[0]
    assert rec.region_id == "test-region"
    assert rec.verdict.startswith("EMPTY")
    # Method-versioned: a non-empty capability envelope (validate_empty_region
    # would have raised on append otherwise).
    assert rec.method_capability.capability_tags
    assert rec.search_extent["points_total"] > 0
    assert rec.prune_gates


def test_re_sweep_gate_skips_same_method_on_recorded_empty(
    tmp_path: Path, monkeypatch: object
) -> None:
    monkeypatch.setattr(run, "discover_novel_moon", lambda **kw: iter(()))  # type: ignore[attr-defined]
    empty_path = tmp_path / "empty_regions.jsonl"
    queue_path = tmp_path / "review_queue.jsonl"
    argv = [
        "--epochs",
        "2",
        "--workers",
        "1",
        "--region-id",
        "test-region",
        "--empty-regions",
        str(empty_path),
        "--queue",
        str(queue_path),
    ]
    # First run records the empty region.
    assert run.main(argv) == 0
    n_after_first = len(load_empty_regions_list(empty_path))
    assert n_after_first == 1

    # Second run with the SAME (single-ellipse no-leveraging) method must SKIP
    # (should_sweep -> False): no new record appended.
    assert run.main(argv) == 0
    n_after_second = len(load_empty_regions_list(empty_path))
    assert n_after_second == 1  # skipped, no duplicate
