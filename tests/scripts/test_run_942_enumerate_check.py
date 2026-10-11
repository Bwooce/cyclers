"""run_942_enumerate.py check: the coverage validator reports a missing and a duplicated structure.

A synthetic two-shard layout for a small cell: every structure of the cell done once, then one
structure removed (missing), one repeated in the other shard (duplicate) and one zero line copied
(duplicated zero line).
"""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path
from typing import Any

import pytest

REPO = Path(__file__).resolve().parents[2]
SETTINGS = {
    "cell": "vm",
    "k": "1",
    "max_returns": "1,0",
    "transfer_revs": "0",
    "generic_revs": "1",
    "resonant_only": "",
    "visits": 1,
    "n_phase": 24,
    "n_split": 10,
    "n_refine": 60,
}


def _enum() -> Any:
    spec = importlib.util.spec_from_file_location(
        "run_942_enumerate_check", REPO / "scripts" / "run_942_enumerate.py"
    )
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    sys.modules["run_942_enumerate_check"] = mod
    spec.loader.exec_module(mod)
    return mod


def _write(d: Path, keys: list[str], zero_lines: list[dict[str, Any]]) -> None:
    d.mkdir(parents=True)
    (d / "settings.json").write_text(json.dumps(SETTINGS))
    (d / "structures.jsonl").write_text(
        "".join(json.dumps({"key": k, "k": 1, "n_zeros": 0, "statuses": []}) + "\n" for k in keys)
    )
    (d / "zeros.jsonl").write_text("".join(json.dumps(z) + "\n" for z in zero_lines))


def test_complete_layout_is_ok(tmp_path: Path) -> None:
    mod = _enum()
    keys = sorted(mod.expected_keys(SETTINGS))
    assert len(keys) >= 4
    _write(tmp_path / "s0", keys[0::2], [{"key": keys[0], "x_days": [1.0, 2.0]}])
    _write(tmp_path / "s1", keys[1::2], [{"key": keys[1], "x_days": [1.0, 2.0]}])
    summary, bad = mod.check_dirs([tmp_path / "s0", tmp_path / "s1"])
    assert bad == []
    assert f"expected {len(keys)} structures, done {len(keys)}" in summary


def test_missing_and_duplicated_structure_are_reported(tmp_path: Path) -> None:
    mod = _enum()
    keys = sorted(mod.expected_keys(SETTINGS))
    zero = {"key": keys[0], "x_days": [1.0, 2.0]}
    s0 = keys[0::2][1:]  # keys[0] is missing from shard 0 ...
    s1 = [*keys[1::2], keys[2]]  # ... and keys[2] is also in shard 1: done twice
    _write(tmp_path / "s0", s0, [zero, zero])  # duplicated zero line
    _write(tmp_path / "s1", s1, [])
    summary, bad = mod.check_dirs([tmp_path / "s0", tmp_path / "s1"])
    assert "1 structures not done" in bad
    assert "1 structures done more than once" in bad
    assert "1 duplicated zero lines" in bad
    assert f"done {len(keys) - 1}" in summary


def test_main_check_exits_nonzero_on_failure(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    mod = _enum()
    keys = sorted(mod.expected_keys(SETTINGS))
    _write(tmp_path / "s0", keys[1:], [])
    monkey_argv = ["run_942_enumerate.py", "check", str(tmp_path / "s0")]
    old = sys.argv
    sys.argv = monkey_argv
    try:
        with pytest.raises(SystemExit) as exc:
            mod.main()
    finally:
        sys.argv = old
    assert exc.value.code == 1
    assert "FAIL 1 structures not done" in capsys.readouterr().out
