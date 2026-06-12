"""Schema v4.5 field tests (spec §16.7.12).

Covers the additive top-level ``validation_level`` enum (V0-V5, the spec §14
gauntlet level) back-filled mechanically by
``scripts/backfill_validation_level.py``.

Both gates are exercised: the JSON Schema structural gate (enum membership) and
the Python semantic gate (``validate_validation_level``), which adds the rule
JSON Schema cannot express — a row may declare a level above V0 only when the
recorded mechanical evidence justifies it (the over-claim guard).
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import jsonschema  # type: ignore[import-untyped]
import pytest
import yaml  # type: ignore[import-untyped]

from cyclerfinder.data.validate import validate_validation_level
from scripts.backfill_validation_level import _LEVEL_BY_ID, preflight_registry_drift

CATALOGUE_PATH = Path(__file__).resolve().parent.parent.parent / "data" / "catalogue.yaml"
SCHEMA_PATH = Path(__file__).resolve().parent.parent.parent / "data" / "catalogue.schema.json"


def _load_schema() -> dict[str, Any]:
    return json.loads(SCHEMA_PATH.read_text())  # type: ignore[no-any-return]


def _load_rows() -> list[dict[str, Any]]:
    return yaml.safe_load(CATALOGUE_PATH.read_text())  # type: ignore[no-any-return]


def test_schema_version_is_at_least_4_5() -> None:
    # v4.5 fields remain present; the schema version has since advanced (>= 4.5).
    assert float(_load_schema()["version"]) >= 4.5


def test_valid_level_passes_schema() -> None:
    schema = _load_schema()
    jsonschema.validate([{"id": "ok", "validation_level": "V1"}], schema)  # must not raise


def test_unknown_level_fails_schema() -> None:
    schema = _load_schema()
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.validate([{"id": "bad", "validation_level": "V9"}], schema)


def test_live_catalogue_validates() -> None:
    jsonschema.validate(_load_rows(), _load_schema())


def test_semantic_gate_rejects_unjustified_level_above_v0() -> None:
    """A row claiming V1+ without recorded mechanical evidence is refused."""
    rows = [{"id": "some-random-row", "validation_level": "V1"}]
    errors = validate_validation_level(rows)
    assert errors, "expected an over-claim violation"
    assert any("validation_level" in e for e in errors)


def test_semantic_gate_rejects_unjustified_nested_validation_level() -> None:
    """C1 pin (#196): the NESTED ``validation.level`` field (the spelling the
    M7 writeback helpers maintain) is held to the same over-claim guard as the
    top-level ``validation_level`` tag — the historical bypass route."""
    rows = [{"id": "some-random-row", "validation": {"level": "V3"}}]
    errors = validate_validation_level(rows)
    assert errors, "expected an over-claim violation on nested validation.level"
    assert any("validation.level" in e for e in errors)


def test_semantic_gate_accepts_nested_v0_and_registered_nested_level() -> None:
    """Nested V0 needs no evidence; a registered nested (id, level) passes."""
    rows = [
        {"id": "anything", "validation": {"level": "V0"}},
        {"id": "aldrin-classic-em-k1-outbound", "validation": {"level": "V2"}},
    ]
    assert validate_validation_level(rows) == []


def test_semantic_gate_rejects_unknown_nested_level_enum() -> None:
    """A nested level outside V0..V5 is refused (enum check, nested spelling)."""
    rows = [{"id": "bad", "validation": {"level": "V9"}}]
    errors = validate_validation_level(rows)
    assert any("not one of V0..V5" in e for e in errors)


def test_semantic_gate_accepts_v0_anywhere() -> None:
    """V0 (internal-consistency floor) needs no evidence pointer."""
    rows = [{"id": "anything", "validation_level": "V0"}]
    assert validate_validation_level(rows) == []


def test_semantic_gate_accepts_known_justified_row() -> None:
    """The Aldrin outbound row's V1 is in the evidence registry."""
    rows = [{"id": "aldrin-classic-em-k1-outbound", "validation_level": "V1"}]
    assert validate_validation_level(rows) == []


def test_semantic_gate_accepts_known_justified_inbound_row() -> None:
    """The Aldrin inbound row's V1 (#125 Part 1) is in the evidence registry."""
    rows = [{"id": "aldrin-classic-em-k1-inbound", "validation_level": "V1"}]
    assert validate_validation_level(rows) == []


def test_semantic_gate_accepts_outbound_v2_powered() -> None:
    """The Aldrin outbound row's V2 (the 2026-06-07 §14 V2-powered class-split
    amendment) is in the evidence registry."""
    rows = [{"id": "aldrin-classic-em-k1-outbound", "validation_level": "V2"}]
    assert validate_validation_level(rows) == []


def test_live_v1_census_matches_recorded_evidence() -> None:
    """Mechanical §14 application. The powered Aldrin OUTBOUND is V2 (2026-06-07
    the §14 V2 class-split amendment: it clears the amended V2-powered gate — >=3
    consecutive in-family cycles each achieving their encounters with maintenance
    applied AND bounded intra-cycle drift). The Aldrin INBOUND is V1 (#125 Part 1's
    real-DE440 lamberthub + Kepler re-propagation evidence; it lands off-family at
    dV~0 so it is NOT V2-powered). Four Russell free-return rows whose single
    ellipse forms a closed, V_inf-continuous E->M->E cycler are V1 (#137 Part 1 +
    Part 3, circular like-for-like — the fourth, 9.353Gg2, promoted by the dense
    phase scan; they are NOT V2-ballistic — they are multi-arc single-ellipse
    slices, so no continuous >=3-lap trajectory exists). Six further Russell 2004
    Table 3.4 free-return rows are V1 (closer sweep 2026-06-08, #142 continuation —
    closed circular like-for-like single-ellipse arcs clearing §14 V1 mechanics;
    docs/notes/2026-06-08-closer-sweep-v1-candidates.md). Every other tagged row is
    V0; the only V2 is the Aldrin outbound. Two rows are V3: S1L1 (#167/#94 — the
    V3-ballistic type specimen) and russell-ch4-8.049gGf2 (#175/#170 — the first
    V3-powered, App-C #188 under its documented 420 m/s budget). Four further
    descriptor-bearing Russell Ch.4 rows are V1 (#181 ToF-fix 2026-06-10 — real-DE440
    closure via the joint (epoch, ToF) closer, emerged v∞ within 0.08 km/s of the
    sourced anchor, lamberthub + Kepler reprop; the 9.353Gg2/9.94Gg3 siblings were
    already V1 and carry the closure as added evidence;
    docs/notes/2026-06-10-tof-fix-closure-results.md). (V1=15, V2=1, V3=2.)"""
    rows = _load_rows()
    byid = {r["id"]: r.get("validation_level") for r in rows}
    assert byid.get("aldrin-classic-em-k1-outbound") == "V2"
    assert byid.get("aldrin-classic-em-k1-inbound") == "V1"
    above_v0 = {rid: lvl for rid, lvl in byid.items() if lvl not in (None, "V0")}
    assert above_v0 == {
        "aldrin-classic-em-k1-outbound": "V2",
        "aldrin-classic-em-k1-inbound": "V1",
        "russell-ch4-5.30gGf3": "V1",
        "russell-ch4-9.94Gg3": "V1",
        "russell-ch4-5.75ggF3": "V1",
        "russell-ch4-9.353Gg2": "V1",
        # closer sweep 2026-06-08 (#142 continuation): six Russell 2004 Table 3.4
        # rows whose single circular-coplanar ellipse closes to a V_inf-continuous
        # E->M->E cycler and clears §14 V1 mechanics like-for-like (emerged V∞
        # within 0.5 km/s of the sourced anchor, lamberthub + Kepler reprop, Mars
        # V_inf continuity intact). tests/search/test_closer_sweep_v1.py;
        # docs/notes/2026-06-08-closer-sweep-v1-candidates.md.
        "russell-ocampo-3.1.1+2": "V1",
        "russell-ocampo-3.1.3+0": "V1",
        "russell-ocampo-4.1.1-4": "V1",
        "russell-ocampo-4.1.2-2": "V1",
        "russell-ocampo-4.1.4-1": "V1",
        "russell-ocampo-4.6.3+0": "V1",
        # #181 ToF-fix (2026-06-10): four descriptor-bearing Russell Ch.4 rows
        # close on the real DE440 ephemeris via the joint (epoch, ToF) self-seed
        # closer (Stage-B coplanar-branch-ToF artifact corrected) and clear §14 V1
        # (emerged E/M v∞ within 0.08 km/s of the sourced anchor, lamberthub +
        # Kepler reprop). V3-CANDIDATES pending the multi-lap horizon-TCM.
        # docs/notes/2026-06-10-tof-fix-closure-results.md.
        "russell-ch4-3.78Gg3": "V1",
        "russell-ch4-6.44Gg3": "V1",
        "russell-ch4-3.64gGg3": "V1",
        "russell-ch4-5.30ggF3": "V1",
        # #167/#94: S1L1's FIRST V3 — corrected-topology closure CONFIRMED on the
        # real DE440 ephemeris (independent REBOUND/IAS15, all 7 Mars encounters in
        # 3x Mars-SOI at the published per-leg v∞, holding under Mars-perturbed
        # gravity; g-arcs sub-Mars). It is the V3-BALLISTIC type specimen (62 m/s
        # continuous TCM, 52% of the 120 m/s ballistic budget).
        # tests/nbody/test_s1l1_corrected_nbody.py;
        # docs/notes/2026-06-08-s1l1-corrected-closure-results.md.
        "russell-ch4-4.991gG2": "V3",
        # #175: the §14 V3 class-split. russell-ch4-8.049gGf2 (App-C #188) is the
        # first V3-POWERED — encounters confirmed on real DE440 by an independent
        # REBOUND/IAS15 integrator (#170), continuous TCM 163.6 m/s under its OWN
        # documented App-C ΔV (420 m/s; 0.39x budget), not the 120 m/s ballistic
        # bar. Encoded as the V3 enum value (powered class in _LEVEL_EVIDENCE text,
        # as V2-powered is under V2). The sibling #192 (russell-ch4-8.165Gfh-f2) is
        # NOT promoted (TCM 2040.6 m/s EXCEEDS its 1678 m/s budget, 1.22x).
        # docs/notes/2026-06-08-v3-powered-classsplit.md.
        "russell-ch4-8.049gGf2": "V3",
    }, above_v0
    # Exactly one row carries V2 today (the powered Aldrin outbound). Two rows carry
    # V3: S1L1 (#167/#94 — V3-ballistic) and russell-ch4-8.049gGf2 (#175/#170 — the
    # first V3-powered, App-C #188 under its documented 420 m/s budget).
    assert sum(1 for lvl in byid.values() if lvl == "V2") == 1
    assert sum(1 for lvl in byid.values() if lvl == "V3") == 2
    assert not any(lvl in ("V4", "V5") for lvl in byid.values())


def test_live_catalogue_validation_level_semantic_clean() -> None:
    """The whole live catalogue passes the validation_level semantic gate."""
    assert validate_validation_level(_load_rows()) == []


# ---------------------------------------------------------------------------
# M1 preflight (#196) — backfill registry must agree with _LEVEL_EVIDENCE
# ---------------------------------------------------------------------------


def test_backfill_preflight_rejects_fabricated_drift_entry() -> None:
    """M1 pin (#196): an (id, level>V0) present in the script's _LEVEL_BY_ID but
    absent from the validator's _LEVEL_EVIDENCE is drift — the preflight must
    name it so the script refuses to write."""
    drifted = {**_LEVEL_BY_ID, "fabricated-drift-row": "V2"}
    errors = preflight_registry_drift(drifted)
    assert errors, "expected the fabricated drift entry to be rejected"
    assert any("fabricated-drift-row" in e for e in errors)
    assert any("_LEVEL_EVIDENCE" in e for e in errors)


def test_backfill_preflight_allows_v0_floor_entries() -> None:
    """V0 (the internal-consistency floor) never needs an evidence pointer."""
    assert preflight_registry_drift({"any-row-at-the-floor": "V0"}) == []


def test_backfill_preflight_live_registries_agree() -> None:
    """The script's live _LEVEL_BY_ID and the validator's _LEVEL_EVIDENCE have
    NOT drifted — every level above V0 is evidence-backed."""
    assert preflight_registry_drift(_LEVEL_BY_ID) == []
