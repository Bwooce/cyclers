"""M7 writeback helper tests — spec §16.1 / §16.4 / §16.5, plan §3.3.

Covers :mod:`cyclerfinder.data.writeback`:

* **M7 GATE** — :func:`test_v2_writeback_populates_validation_block`
  (M6a ``StabilityReport`` → ``validation.gates.V2`` integration).
* V0/V1/V3 field copying + level promotion.
* ``record_rediscovery`` idempotency + attribution preservation.
* ``register_discovery`` status/source/attribution.
* ``serialise_entry_yaml`` round-trip.
"""

from __future__ import annotations

import yaml  # type: ignore[import-untyped]

from cyclerfinder.data.catalog import CatalogueEntry, _entry_from_yaml, load_catalog
from cyclerfinder.data.writeback import (
    apply_v0_v1_to_entry,
    apply_v2_to_entry,
    apply_v3_to_entry,
    record_rediscovery,
    register_discovery,
    serialise_entry_yaml,
)
from cyclerfinder.verify.crosscheck import LambertCrosscheckResult
from cyclerfinder.verify.propagate import StabilityReport
from cyclerfinder.verify.real_closure import RealClosureResult


def _entry() -> CatalogueEntry:
    """A real constructible catalogue entry for writeback exercises."""
    return load_catalog().by_id["s1l1-2syn-em-cpom"]


def _stability_report(*, stable: bool, max_drift_km: float) -> StabilityReport:
    return StabilityReport(
        cycler_id="s1l1-2syn-em-cpom",
        n_laps_propagated=3,
        max_drift_km=max_drift_km,
        max_drift_lap_index=2,
        per_lap_drift_km=(100.0, 200.0, max_drift_km),
        stable=stable,
        per_lap_dv=(0.0, 0.0, 0.0),
        total_tcm_dv=0.0,
        frame_used="dynamic",
    )


# ---------------------------------------------------------------------------
# M7 GATE — V2 writeback (M6a StabilityReport integration)
# ---------------------------------------------------------------------------


def test_v2_writeback_populates_validation_block() -> None:
    """M7 GATE — plan §4 item 10 / table row.

    ``apply_v2_to_entry`` copies the ``StabilityReport`` into the V2
    gate block; the input entry stays unmodified (frozen).
    """
    entry = _entry()
    report = _stability_report(stable=True, max_drift_km=12345.6)
    new_entry = apply_v2_to_entry(entry, report)

    v2 = new_entry.validation["gates"]["V2"]
    assert v2["pass"] is True
    assert v2["max_drift_km"] == 12345.6
    assert v2["n_laps"] == 3
    assert v2["per_lap_drift_km"] == [100.0, 200.0, 12345.6]
    assert v2["frame_used"] == report.frame_used

    # Input entry is frozen / untouched.
    assert "gates" not in entry.validation or "V2" not in entry.validation.get("gates", {})


def test_v2_writeback_promotes_validation_level() -> None:
    """A passing V2 promotes ``validation.level`` to ``"V2"``."""
    entry = _entry()
    new_entry = apply_v2_to_entry(entry, _stability_report(stable=True, max_drift_km=1.0))
    assert new_entry.validation["level"] == "V2"


def test_v2_writeback_failing_does_not_promote() -> None:
    """A failing V2 leaves the level at the implicit V0 floor."""
    entry = _entry()
    new_entry = apply_v2_to_entry(entry, _stability_report(stable=False, max_drift_km=9.9e9))
    assert new_entry.validation["gates"]["V2"]["pass"] is False
    assert new_entry.validation["level"] == "V0"


# ---------------------------------------------------------------------------
# V0 / V1 writeback
# ---------------------------------------------------------------------------


def _xcheck(passed: bool, max_diff_mps: float, leg_index: int = 0) -> LambertCrosscheckResult:
    return LambertCrosscheckResult(
        leg_index=leg_index,
        mine_v1_kms=(1.0, 2.0, 3.0),
        lamberthub_izzo_v1_kms=(1.0, 2.0, 3.0),
        lamberthub_gooding_v1_kms=(1.0, 2.0, 3.0),
        max_diff_mps=max_diff_mps,
        passed=passed,
    )


def test_v0_v1_writeback_all_pass_promotes_to_v1() -> None:
    entry = _entry()
    v1 = (_xcheck(True, 1e-5, 0), _xcheck(True, 2e-5, 1))
    new_entry = apply_v0_v1_to_entry(entry, True, v1)
    gates = new_entry.validation["gates"]
    assert gates["V0"]["pass"] is True
    assert gates["V1"]["pass"] is True
    assert gates["V1"]["max_diff_mps"] == 2e-5
    assert new_entry.validation["level"] == "V1"


def test_v0_v1_writeback_one_failing_leg_fails_v1() -> None:
    entry = _entry()
    v1 = (_xcheck(True, 1e-5, 0), _xcheck(False, 5.0, 1))
    new_entry = apply_v0_v1_to_entry(entry, True, v1)
    assert new_entry.validation["gates"]["V1"]["pass"] is False
    assert new_entry.validation["level"] == "V0"


def test_v0_v1_writeback_empty_v1_not_applicable() -> None:
    """An empty cross-check tuple is "V1 not applicable", not a pass."""
    entry = _entry()
    new_entry = apply_v0_v1_to_entry(entry, True, ())
    assert new_entry.validation["gates"]["V1"]["pass"] is False
    assert new_entry.validation["level"] == "V0"


# ---------------------------------------------------------------------------
# V3 writeback
# ---------------------------------------------------------------------------


def test_v3_writeback_populates_gate_and_metrics() -> None:
    entry = _entry()
    report = RealClosureResult(
        cycler_id="s1l1-2syn-em-cpom",
        n_cycles_propagated=2,
        max_drift_km=5000.0,
        per_cycle_drift_km=(1000.0, 5000.0),
        per_encounter_vinf_mismatch_kms=(0.0, 0.0),
        closes=True,
        v3_status="ok",
        horizon_tcm_mps=120.0,
        per_cycle_tcm_mps=(60.0, 60.0),
        frame_used="dynamic",
        t_start_sec=0.0,
    )
    new_entry = apply_v3_to_entry(entry, report)
    v3 = new_entry.validation["gates"]["V3"]
    assert v3["pass"] is True
    assert v3["horizon_tcm_mps"] == 120.0
    assert new_entry.validation["metrics"]["horizon_tcm_dv_mps"] == 120.0
    assert new_entry.validation["level"] == "V3"


# ---------------------------------------------------------------------------
# Rediscovery audit trail
# ---------------------------------------------------------------------------


def test_record_rediscovery_idempotent() -> None:
    entry = _entry()
    once = record_rediscovery(entry, "run-1", "cell-1", "2026-06-01")
    twice = record_rediscovery(once, "run-1", "cell-1", "2026-06-01")
    assert len(once.discovery["rediscoveries"]) == 1
    assert twice.discovery["rediscoveries"] == once.discovery["rediscoveries"]


def test_record_rediscovery_distinct_keys_append() -> None:
    entry = _entry()
    e1 = record_rediscovery(entry, "run-1", "cell-1", "2026-06-01")
    e2 = record_rediscovery(e1, "run-2", "cell-1", "2026-06-02")
    assert len(e2.discovery["rediscoveries"]) == 2


def test_record_rediscovery_preserves_attribution() -> None:
    entry = _entry()
    new_entry = record_rediscovery(entry, "run-1", "cell-1", "2026-06-01")
    assert new_entry.first_published == entry.first_published
    assert new_entry.priority_date == entry.priority_date


# ---------------------------------------------------------------------------
# Discovery registration
# ---------------------------------------------------------------------------


def test_register_discovery_sets_candidate_novel() -> None:
    skeleton = _entry()
    new_entry = register_discovery(skeleton, "run-9", "cell-9", "2026-06-01", "0.7.0")
    assert new_entry.source == "this-project"
    assert new_entry.our_status == "candidate-novel"
    assert new_entry.first_published is None
    assert new_entry.priority_date == "2026-06-01"
    assert new_entry.discovery_run == {
        "run_id": "run-9",
        "cell_id": "cell-9",
        "date": "2026-06-01",
        "finder_version": "0.7.0",
    }


# ---------------------------------------------------------------------------
# YAML round-trip
# ---------------------------------------------------------------------------


def test_serialise_entry_yaml_round_trip() -> None:
    entry = _entry()
    text = serialise_entry_yaml(entry)
    reloaded = _entry_from_yaml(yaml.safe_load(text))
    assert reloaded == entry


def test_serialise_entry_yaml_round_trip_after_v2() -> None:
    """A V2-updated entry serialises and reloads with its V2 block."""
    entry = apply_v2_to_entry(_entry(), _stability_report(stable=True, max_drift_km=42.0))
    text = serialise_entry_yaml(entry)
    reloaded = _entry_from_yaml(yaml.safe_load(text))
    assert reloaded.validation["gates"]["V2"]["max_drift_km"] == 42.0
    assert reloaded == entry
