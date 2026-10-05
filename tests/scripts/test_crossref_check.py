"""#964 -- scripts/crossref_check.py prints CONFIRMED only on a journal/volume match.

Offline: every test runs against Crossref responses recorded on 2026-10-05 in
``tests/scripts/fixtures/`` (trimmed to the fields the helper reads).
"""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path
from types import ModuleType
from typing import Any

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT = REPO_ROOT / "scripts" / "crossref_check.py"
FIX = Path(__file__).resolve().parent / "fixtures"


def _load() -> ModuleType:
    spec = importlib.util.spec_from_file_location("crossref_check", SCRIPT)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod  # dataclasses resolve their module via sys.modules
    spec.loader.exec_module(mod)
    return mod


cc = _load()


def _rec(name: str) -> dict[str, Any]:
    data: dict[str, Any] = json.loads((FIX / name).read_text())["message"]
    return data


HOLLISTER = "crossref-10.2514_3.29664.json"  # JSR 6(4):366-369, 1969
GILLESPIE_CONF = "crossref-10.2514_6.1966-37.json"  # AIAA 66-37 conference paper


def test_journal_match_confirms() -> None:
    v = cc.assess(
        _rec(HOLLISTER),
        cc.Expect(container="Journal of Spacecraft and Rockets", volume="6", year=1969),
    )
    assert v.status == "CONFIRMED", v.reasons
    assert v.doi == "10.2514/3.29664"


def test_volume_mismatch_is_not_confirmed() -> None:
    v = cc.assess(
        _rec(HOLLISTER), cc.Expect(container="Journal of Spacecraft and Rockets", volume="7")
    )
    assert v.status == "MISMATCH"
    assert any("volume" in r for r in v.reasons)


def test_conference_record_does_not_confirm_a_journal_citation() -> None:
    # The citation says JSR 4; the DOI Crossref ranks first is the AIAA conference paper.
    v = cc.assess(
        _rec(GILLESPIE_CONF), cc.Expect(container="Journal of Spacecraft and Rockets", volume="4")
    )
    assert v.status == "MISMATCH"
    assert any("conference" in r for r in v.reasons)


def test_no_container_or_volume_claimed_is_unchecked() -> None:
    assert cc.assess(
        _rec(HOLLISTER), cc.Expect(title="Periodic orbits for interplanetary flight")
    ).status == ("UNCHECKED")


def test_title_and_year_guards() -> None:
    bad_title = cc.assess(
        _rec(HOLLISTER), cc.Expect(volume="6", title="Venus swingby mission mode manned Mars")
    )
    assert bad_title.status == "MISMATCH"
    bad_year = cc.assess(_rec(HOLLISTER), cc.Expect(volume="6", year=1975))
    assert bad_year.status == "MISMATCH"


def test_missing_record_is_not_found() -> None:
    assert cc.assess(None, cc.Expect(volume="6")).status == "NOT_FOUND"


def test_search_picks_the_journal_version_not_the_first_hit() -> None:
    items = json.loads((FIX / "crossref-search-gillespie-ross.json").read_text())["message"][
        "items"
    ]
    assert items[0]["DOI"] == "10.2514/6.1966-37"  # Crossref ranks the conference paper first
    best, all_v = cc.pick_match(
        items, cc.Expect(container="Journal of Spacecraft", volume="4", year=1967)
    )
    assert best.status == "CONFIRMED"
    assert best.doi == "10.2514/3.28830"
    assert all_v[0].status == "MISMATCH"


def test_main_offline_exit_codes(
    monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    rec = _rec(HOLLISTER)
    monkeypatch.setattr(cc, "fetch_doi", lambda doi, fetch=None: rec)
    assert (
        cc.main(["doi", "10.2514/3.29664", "--container", "Journal of Spacecraft", "--volume", "6"])
        == 0
    )
    assert capsys.readouterr().out.startswith("CONFIRMED 10.2514/3.29664")
    assert cc.main(["doi", "10.2514/3.29664", "--volume", "9"]) == 1
