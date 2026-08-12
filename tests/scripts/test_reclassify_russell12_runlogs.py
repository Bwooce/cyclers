"""#829 — the historical runlog re-classification rule.

``scripts/reclassify_russell12_runlogs.py`` applies the spec §14 V0 hard-constraint
gate to every persisted ``russell12-*`` runlog record. These tests pin the rule
itself (on synthetic records) and the committed re-classification of the real
runlogs, so a future edit cannot quietly un-retract a verdict.
"""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path
from types import ModuleType
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT = REPO_ROOT / "scripts" / "reclassify_russell12_runlogs.py"
ARTIFACT = REPO_ROOT / "data" / "runs" / "russell12-829-reclassification.json"


def _load() -> ModuleType:
    spec = importlib.util.spec_from_file_location("reclassify_russell12_runlogs", SCRIPT)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _rec(outcome: str, *, bend: bool | None = None, cap: bool | None = None) -> dict[str, Any]:
    audit: dict[str, Any] = {}
    if bend is not None:
        audit["bend_feasible"] = bend
    if cap is not None:
        audit["vinf_cap_ok"] = cap
    return {"row_id": "r", "outcome": outcome, "solver_audit": audit}


def test_bend_infeasible_closure_verdict_is_retracted() -> None:
    out = _load().reclassify_record(_rec("CLOSE-AND-MATCH", bend=False, cap=True))
    assert out["status"] == "RETRACTED"
    assert out["revised_outcome"] == "CLOSE-INADMISSIBLE"
    assert out["violations"] == ["bend"]


def test_vinf_cap_breach_is_also_a_retraction() -> None:
    out = _load().reclassify_record(_rec("CLOSE-OFF-ANCHOR", bend=True, cap=False))
    assert out["status"] == "RETRACTED"
    assert out["violations"] == ["vinf_cap"]


def test_admissible_closure_verdict_stands() -> None:
    out = _load().reclassify_record(_rec("CLOSE-AND-MATCH", bend=True, cap=True))
    assert out["status"] == "STANDS"
    assert out["revised_outcome"] == "CLOSE-AND-MATCH"


def test_no_close_is_not_reclassified() -> None:
    """NO-CLOSE asserts the ABSENCE of a closure — a bend-blind filter cannot
    have over-claimed it, so it is never touched."""
    out = _load().reclassify_record(_rec("NO-CLOSE"))
    assert out["status"] == "N/A"
    assert out["revised_outcome"] == "NO-CLOSE"


def test_unrecorded_flags_are_unknown_not_stands() -> None:
    """A record whose audit never stored the flags cannot be cleared — saying
    'STANDS' there would assert something the data does not support."""
    out = _load().reclassify_record(_rec("CLOSE-AND-MATCH"))
    assert out["status"] == "UNKNOWN"


def test_committed_reclassification_matches_the_runlogs() -> None:
    """The committed artifact is exactly what the script produces from the
    runlogs in the tree — it is a data claim, so it must be reproducible."""
    result = _load().reclassify_runlogs()
    committed = json.loads(ARTIFACT.read_text())
    assert result == committed


def test_820s_two_close_and_match_verdicts_survive_the_gate() -> None:
    """#820's headline grid verdicts were bend-FEASIBLE, so the #829 gate does
    not touch them; every other verdict in that run is retracted."""
    files = _load().reclassify_runlogs()["files"]
    run820 = files["russell12-circular-20260811T-820-reposed.jsonl"]
    stands = {e["row_id"] for e in run820 if e["status"] == "STANDS"}
    assert stands == {"russell-ch4-9.353Gg2", "russell-ch4-3.78Gg3"}
    assert all(e["status"] in {"STANDS", "RETRACTED"} for e in run820)
