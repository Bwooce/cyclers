"""#1071: the #973 and #1065 gauntlet writers assign a distinct, stable ``cid`` to every candidate,
including two candidates that share a ``key`` (the defect).

Synthetic shard: two real committed candidates with the same key (different dates and V_inf) become
two passing zeros in a tiny ``zeros.jsonl`` and go through the writer's own ``gauntlet`` command.
"""

from __future__ import annotations

import collections
import importlib.util
import json
import sys
from pathlib import Path
from typing import Any

import pytest

from cyclerfinder.data.candidate_ids import candidate_cid

REPO = Path(__file__).resolve().parents[2]


def _load(name: str, rel: str) -> Any:
    spec = importlib.util.spec_from_file_location(name, REPO / rel)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def _duplicate_key_pair(rel: str) -> list[dict[str, Any]]:
    cands = json.loads((REPO / rel).read_text())["candidates"]
    by_key: dict[str, list[dict[str, Any]]] = collections.defaultdict(list)
    for c in cands:
        by_key[c["key"]].append(c)
    pair = next(v for v in by_key.values() if len(v) >= 2 and "|H" not in v[0]["key"])
    return pair[:2]


def _write_shard(tmp: Path, cands: list[dict[str, Any]], flybys_field: str) -> Path:
    shard = tmp / "shard"
    shard.mkdir()
    (shard / "structures.jsonl").write_text(
        json.dumps({"key": cands[0]["key"], "n_zeros": 2}) + "\n"
    )
    lines = []
    for c in cands:
        lines.append(
            {
                "key": c["key"],
                "k": c["k"],
                "x_days": c["x_days"],
                "status": "pass",
                "status_hm_floor": c.get("status_hm_floor", "pass"),
                "worst_ratio": c["worst_ratio"],
                "max_turn_deg": c["max_turn_deg"],
                "min_required_alt_km": c["min_required_alt_km"],
                "max_encounter_miss_km": 0.0,
                "vinf_kms": c.get("vinf_kms")
                or {"Charon": c["vinf_charon_kms"], c.get("moon", "Nix"): c["vinf_moon_kms"]},
                "r_min_km": c["r_min_km"],
                "r_max_km": c["r_max_km"],
                "flybys": c[flybys_field],
            }
        )
    (shard / "zeros.jsonl").write_text("".join(json.dumps(z) + "\n" for z in lines))
    return shard


def test_run_973_gauntlet_assigns_distinct_cids(tmp_path: Path) -> None:
    mod = _load("run_973_enumerate_1071", "scripts/run_973_enumerate.py")
    pair = _duplicate_key_pair("data/973_gc/k4_gauntlet.json")
    shard = _write_shard(tmp_path, pair, "flyby_table")
    out = tmp_path / "g.json"
    mod.cmd_gauntlet(["--cell", "gc", str(shard), "--out", str(out)])
    cands = json.loads(out.read_text())["candidates"]
    assert len(cands) == 2 and cands[0]["key"] == cands[1]["key"]
    assert cands[0]["cid"] != cands[1]["cid"]
    for i, c in enumerate(cands):
        assert c["cid"] == candidate_cid("gc", i, c["x_days"])


@pytest.mark.skipif(
    not (REPO / "scripts/run_1065_barycentric_enumerate.py").exists(), reason="no 1065"
)
def test_run_1065_gauntlet_assigns_distinct_cids(tmp_path: Path) -> None:
    mod = _load("run_1065_enumerate_1071", "scripts/run_1065_barycentric_enumerate.py")
    pair = _duplicate_key_pair("data/998_pluto_smallmoons/pn_gauntlet.json")
    for c in pair:
        c["moon"] = "Nix"
    shard = _write_shard(tmp_path, pair, "charon_flybys")
    out = tmp_path / "g.json"
    mod.cmd_gauntlet(["--cell", "pnb", str(shard), "--out", str(out)])
    cands = json.loads(out.read_text())["candidates"]
    assert len(cands) == 2 and cands[0]["key"] == cands[1]["key"]
    assert cands[0]["cid"] != cands[1]["cid"]
    for i, c in enumerate(cands):
        assert c["cid"] == candidate_cid("pnb", i, c["x_days"])
