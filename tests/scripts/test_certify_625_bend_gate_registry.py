"""#625 -- tests for the driver script that certifies every target entry.

Covers: (1) the Hohmann-floor reproduction for #571's Titan-pairs matches
that entry's own already-published numbers exactly (cross-check against
independently-sourced data, not just internal self-consistency); (2) the
data-grounded real-survivor reproduction for #609 Mars Phobos-Deimos matches
that entry's own recorded ``max_flyby_own_bend_deg_across_subgate_survivors``
(0.0159 deg) -- confirming this script's independent reproduction agrees
with the original production sweep's number; (3) end-to-end ``main()`` runs
clean and certifies the expected bodies while correctly reporting the
genuine Sylvia/ElektraBeta non-certifications (not silently dropping them);
(4) the write-back helper only touches targeted lines/keys on a throwaway
tempfile copy, never the real registry.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

mp = pytest.importorskip("mpmath", reason="mpmath is an optional 'interval' extra (task #610/#625)")

import scripts.certify_625_bend_gate_registry as driver  # noqa: E402


@pytest.fixture(autouse=True)
def _iv_precision() -> None:
    mp.mp.dps = 50
    mp.iv.dps = 50


def test_hohmann_floor_matches_571_own_recorded_numbers() -> None:
    """Cross-check against data/empty_regions.jsonl's own #571 entries (their
    Hohmann-floor numbers were independently computed by
    scripts/verify_571_gate_analytics.py at task #571 time) -- must match to
    displayed precision, confirming this script's reuse of the same helper
    reproduces the same numbers, not a divergent copy."""
    expected = {
        "Mimas": (4.5425, 0.0466),
        "Enceladus": (3.7086, 0.1704),
        "Tethys": (3.0566, 0.7930),
        "Dione": (2.3672, 2.2170),
    }
    iv = mp.iv
    for moon, (exp_vinf, exp_bend) in expected.items():
        result = driver.certify_hohmann_floor(iv, moon, "Titan", "Saturn", label=moon)
        assert result["box"]["vinf_kms"][0] == pytest.approx(exp_vinf, abs=1e-4)
        assert result["sup_bend_deg"] == pytest.approx(exp_bend, abs=1e-3)
        assert result["certified"] is True


def test_mars_phobos_deimos_reproduction_matches_609_entry() -> None:
    """The #609 entry itself records
    ``max_flyby_own_bend_deg_across_subgate_survivors: 0.0159`` (Phobos, the
    binding body, across ALL 52 sub-gate survivors from BOTH directions).
    This test reproduces ONE direction (52/2=26 survivors -- #625's own
    driver only needs one direction, see subgate_vinf_ranges_2moon's
    docstring) and must match that same figure."""
    data = driver.subgate_vinf_ranges_2moon("Phobos", "Deimos", primary="Mars")
    assert data["n_total"] == 256
    assert data["n_sub"] == 26
    iv = mp.iv
    result = driver.certify_data_grounded(
        iv, "Phobos", data["vinf"]["Phobos"], label="Phobos", vinf_hi_kms=50.0
    )
    assert result["sup_bend_deg"] == pytest.approx(0.0159, abs=1e-4)
    assert result["certified"] is True


def test_elektra_gamma_and_delta_certify_beta_does_not() -> None:
    """The genuine #625 finding for Elektra: Gamma/Delta's own bend never
    clears the gate across the real survivor range, but Beta's does at its
    own survivor minimum -- so system-level certification routes through
    Gamma/Delta, and Beta is honestly reported as NOT certified individually
    (matches the module docstring's own stated finding)."""
    data = driver.subgate_vinf_ranges_3moon(
        ("ElektraBeta", "ElektraGamma", "ElektraDelta"), primary="Elektra"
    )
    assert data["n_total"] == 96768
    assert data["n_sub"] == 16793
    iv = mp.iv
    beta = driver.certify_data_grounded(
        iv, "ElektraBeta", data["vinf"]["ElektraBeta"], label="Beta", vinf_hi_kms=1.0
    )
    gamma = driver.certify_data_grounded(
        iv, "ElektraGamma", data["vinf"]["ElektraGamma"], label="Gamma", vinf_hi_kms=1.0
    )
    delta = driver.certify_data_grounded(
        iv, "ElektraDelta", data["vinf"]["ElektraDelta"], label="Delta", vinf_hi_kms=1.0
    )
    assert beta["certified"] is False
    assert gamma["certified"] is True
    assert delta["certified"] is True


def test_main_runs_clean_and_reports_partial_607(capsys: pytest.CaptureFixture[str]) -> None:
    rc = driver.main([])  # dry run, no --write
    assert rc == 0
    out = capsys.readouterr().out
    assert "smallbody-multimoon-symmetric-closure-mass-limited-607-2026-07-16: PARTIAL" in out
    assert "Sylvia:Romulus: sup=5.1523 deg certified=False" in out
    assert "mars-phobos-deimos-symmetric-closure-609-2026-07-16: ALL bodies certified" in out


def test_write_results_only_touches_targeted_lines_and_preserves_extra_keys(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    fake_path = tmp_path / "empty_regions.jsonl"
    untouched_line = json.dumps({"region_id": "unrelated-entry", "extra_key": "preserve-me"})
    targeted_line = json.dumps(
        {"region_id": "target-entry", "extra_key": "also-preserve-me", "verdict": "EMPTY"}
    )
    fake_path.write_text(untouched_line + "\n" + targeted_line + "\n", encoding="utf-8")
    monkeypatch.setattr(driver, "DATA_PATH", fake_path)

    fake_result = {"label": "x", "sup_bend_deg": 0.1, "certified": True}
    driver._write_results({"target-entry": [fake_result]})

    lines = fake_path.read_text(encoding="utf-8").splitlines()
    assert len(lines) == 2
    untouched = json.loads(lines[0])
    targeted = json.loads(lines[1])
    assert untouched == {"region_id": "unrelated-entry", "extra_key": "preserve-me"}
    assert targeted["extra_key"] == "also-preserve-me"
    assert targeted["verdict"] == "EMPTY"
    assert targeted["bend_gate_certified_interval"]["task"] == 625
    assert targeted["bend_gate_certified_interval"]["per_body"] == [fake_result]


def test_write_results_raises_on_missing_region_id(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    fake_path = tmp_path / "empty_regions.jsonl"
    fake_path.write_text(json.dumps({"region_id": "present"}) + "\n", encoding="utf-8")
    monkeypatch.setattr(driver, "DATA_PATH", fake_path)
    with pytest.raises(SystemExit):
        driver._write_results({"absent-entry": [{"certified": True}]})
