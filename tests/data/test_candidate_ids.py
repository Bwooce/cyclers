"""#1071: stable candidate ids in the gauntlet JSONs. Pre-registration:
docs/notes/2026-10-12-1071-stable-candidate-ids-prereg.md."""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path
from typing import Any

import pytest

from cyclerfinder.data.candidate_ids import (
    candidate_cid,
    index_by_key,
    index_records,
    record_id,
    same_record,
)

REPO = Path(__file__).resolve().parents[2]


def test_same_key_distinct_cids() -> None:
    a = candidate_cid("ps", 3, [1.0, 2.0, 3.0])
    b = candidate_cid("ps", 4, [1.0, 2.0, 3.0])  # same dates, other position
    c = candidate_cid("ps", 3, [1.0, 2.0, 3.5])  # same position, other dates
    assert len({a, b, c}) == 3
    assert a.startswith("ps-0003-") and len(a.split("-")[2]) == 8
    # rounding to 6 decimals: a 1e-9 change does not move the cid
    assert candidate_cid("ps", 3, [1.0, 2.0, 3.0 + 1e-9]) == a


def test_two_records_with_one_key_are_distinguished() -> None:
    k = "k4|LCharon>Charon/1l|LCharon>Nix/0s|LNix>Charon/0s"
    r1 = {
        "key": k,
        "worst_ratio": 0.30,
        "x_days": [1.0, 2.0],
        "cid": candidate_cid("pn", 0, [1.0, 2.0]),
    }
    r2 = {
        "key": k,
        "worst_ratio": 0.31,
        "x_days": [1.0, 2.5],
        "cid": candidate_cid("pn", 1, [1.0, 2.5]),
    }
    assert r1["cid"] != r2["cid"]
    assert set(index_records([r1, r2])) == {r1["cid"], r2["cid"]}
    with pytest.raises(ValueError):
        index_by_key([r1, r2])  # the old dict-by-key would have kept only r2
    assert not same_record(r1, r2)
    assert same_record(r1, dict(r1))


def test_fallback_id_for_records_without_cid() -> None:
    k = "k3|LGanymede>Callisto/0s|LCallisto>Ganymede/0s"
    a = {"key": k, "flyby_table": [{"ratio": 0.5}, {"body": "x"}]}
    b = {"key": k, "flyby_table": [{"ratio": 0.6}]}
    assert record_id(a) == f"{k}@0.500000000" and record_id(b) != record_id(a)
    assert not same_record(a, b)


@pytest.mark.parametrize(
    "path",
    [
        "data/998_pluto_smallmoons/ps_gauntlet.json",
        "data/998_pluto_smallmoons/pk_gauntlet.json",
        "data/942_cell_ev_gauntlet.json",
        "data/973_gc/k4_gauntlet.json",
    ],
)
def test_committed_gauntlet_json_loads_one_to_one(path: str) -> None:
    """A committed file (no cid) matches its records one-to-one through the fallback id, and some
    of them do repeat a key (that is the defect)."""
    cands = json.loads((REPO / path).read_text())["candidates"]
    assert all("cid" not in c for c in cands)
    idx = index_records(cands)  # raises if two records share a fallback id
    assert len(idx) == len(cands)
    for i, c in enumerate(cands):
        assert idx[record_id(c)] is c
        assert cands[i] is c  # the committed files are also addressed by index
    if "998" in path or path.endswith("k4_gauntlet.json") or "ev_gauntlet" in path:
        assert len({c["key"] for c in cands}) < len(cands)  # keys repeat in these files


def _load(name: str, rel: str) -> Any:
    spec = importlib.util.spec_from_file_location(name, REPO / rel)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def test_writers_import_candidate_cid() -> None:
    """Both writers assign rec["cid"] from candidate_cid (new writes only)."""
    for rel in ("scripts/gauntlet_942.py", "scripts/gauntlet_998.py"):
        src = (REPO / rel).read_text()
        assert "from cyclerfinder.data.candidate_ids import candidate_cid" in src
        assert 'rec["cid"] = candidate_cid(args.cell, len(results), best["x_days"])' in src
