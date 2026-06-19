"""Frozen V0 gate for the #390 SPK-derived Voyager mga_tour catalogue rows.

The two Voyager rows admitted by #390 --
``voyager-1-jupiter-saturn-grand-tour`` (E-J-S) and ``voyager-2-grand-tour``
(E-J-S-U-N, the only four-giant-planet Grand Tour) -- carry per-encounter
hyperbolic-excess velocities (V_inf) DERIVED from each mission's archived NAIF
reconstructed spacecraft SPK at the flown flyby epochs
(:mod:`cyclerfinder.verify.mission_spk`; ``data/390_mission_vinf.jsonl``).

This module is the writeback gate: it freezes the contract that the V_inf
written into the catalogue rows are EXACTLY the values extracted by #390 (to the
display precision of three decimals), so a later edit to either the JSONL
extraction record or the catalogue row cannot silently drift apart. It is a
pure-data test (no network, no SPICE) -- the heavy SPK extraction itself is
gated by ``tests/verify/test_390_spk_vinf_extractor.py``.

Provenance: the V_inf are derived by our own code, so they are V0 (not a golden
sourced EXPECTED). What makes them admissible is the INDEPENDENT, sourced
closest-approach cross-check (each row's source_quotes cite the published CA
radius the SPK geometry reproduces to <1%); this test pins the V_inf-vs-JSONL
identity, not the V_inf-vs-published value (no published V_inf exists -- that is
exactly why #390 had to derive them).
"""

from __future__ import annotations

import json
from pathlib import Path

import yaml  # type: ignore[import-untyped]

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
CATALOGUE_PATH = REPO_ROOT / "data" / "catalogue.yaml"
VINF_JSONL_PATH = REPO_ROOT / "data" / "390_mission_vinf.jsonl"

# Map (catalogue row id, body code) -> (#390 JSONL mission, flyby_body) so the
# row V_inf can be checked against the extraction record it was sourced from.
_ROW_BODY_TO_JSONL: dict[tuple[str, str], tuple[str, str]] = {
    ("voyager-1-jupiter-saturn-grand-tour", "J"): ("Voyager 1", "Jupiter"),
    ("voyager-1-jupiter-saturn-grand-tour", "S"): ("Voyager 1", "Saturn"),
    ("voyager-2-grand-tour", "J"): ("Voyager 2", "Jupiter"),
    ("voyager-2-grand-tour", "S"): ("Voyager 2", "Saturn"),
    ("voyager-2-grand-tour", "U"): ("Voyager 2", "Uranus"),
    ("voyager-2-grand-tour", "N"): ("Voyager 2", "Neptune"),
}

# The exact bodies+sequence the two rows must carry (frozen).
_EXPECTED_SEQUENCE = {
    "voyager-1-jupiter-saturn-grand-tour": (["E", "J", "S"], "E-J-S"),
    "voyager-2-grand-tour": (["E", "J", "S", "U", "N"], "E-J-S-U-N"),
}

# Display-precision tolerance: the catalogue stores V_inf to 3 decimals.
_VINF_TOL_KMS = 5e-4


def _load_rows() -> dict[str, dict]:  # type: ignore[type-arg]
    rows = yaml.safe_load(CATALOGUE_PATH.read_text())
    return {r["id"]: r for r in rows}


def _load_jsonl() -> dict[tuple[str, str], dict]:  # type: ignore[type-arg]
    out: dict[tuple[str, str], dict] = {}  # type: ignore[type-arg]
    for line in VINF_JSONL_PATH.read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        rec = json.loads(line)
        out[(rec["mission"], rec["flyby_body"])] = rec
    return out


def test_voyager_rows_exist_and_are_mga_tour() -> None:
    """Both #390 Voyager rows are present and tagged mga_tour / multi-arc / V0."""
    rows = _load_rows()
    for rid in _EXPECTED_SEQUENCE:
        assert rid in rows, f"{rid} missing from catalogue"
        row = rows[rid]
        assert row["orbit_class"] == "mga_tour", f"{rid}: orbit_class != mga_tour"
        assert row["cycler_class"] == "multi-arc", f"{rid}: cycler_class != multi-arc"
        assert row["epoch_locked"] is True, f"{rid}: epoch_locked must be true"
        assert row["n_returns"] == 1, f"{rid}: n_returns must be 1"
        assert row["validation_level"] == "V0", f"{rid}: validation_level must be V0"
        assert row["orbit_source"] == "derived", f"{rid}: orbit_source must be derived"
        assert row["vinf_source"] == "derived", f"{rid}: vinf_source must be derived"


def test_voyager_rows_sequence_frozen() -> None:
    """bodies + sequence_canonical are exactly the flown itineraries."""
    rows = _load_rows()
    for rid, (bodies, seq) in _EXPECTED_SEQUENCE.items():
        row = rows[rid]
        assert row["bodies"] == bodies, f"{rid}: bodies {row['bodies']} != {bodies}"
        assert row["sequence_canonical"] == seq, (
            f"{rid}: sequence {row['sequence_canonical']!r} != {seq!r}"
        )


def test_voyager_row_vinf_matches_390_jsonl() -> None:
    """Every Voyager row V_inf equals the #390 SPK extraction record it cites."""
    rows = _load_rows()
    jsonl = _load_jsonl()
    checked = 0
    for (rid, body), (mission, flyby_body) in _ROW_BODY_TO_JSONL.items():
        row = rows[rid]
        encs = {e["body"]: e["vinf_kms"] for e in row["vinf_kms_at_encounters"]}
        assert body in encs, f"{rid}: no V_inf entry for body {body}"
        rec = jsonl.get((mission, flyby_body))
        assert rec is not None, f"#390 JSONL missing record for {mission} {flyby_body}"
        row_v = float(encs[body])
        jsonl_v = float(rec["vinf_kms"])
        assert abs(row_v - jsonl_v) <= _VINF_TOL_KMS, (
            f"{rid} {body}: catalogue V_inf {row_v} km/s != #390 JSONL "
            f"{jsonl_v} km/s ({mission} {flyby_body}) beyond {_VINF_TOL_KMS} km/s"
        )
        checked += 1
    assert checked == len(_ROW_BODY_TO_JSONL), "not every (row, body) pair was checked"


def test_voyager2_is_the_four_giant_planet_grand_tour() -> None:
    """Voyager 2 carries a V_inf at all four giant planets (J, S, U, N)."""
    rows = _load_rows()
    row = rows["voyager-2-grand-tour"]
    bodies_with_vinf = {e["body"] for e in row["vinf_kms_at_encounters"]}
    assert {"J", "S", "U", "N"} <= bodies_with_vinf, (
        f"voyager-2-grand-tour must visit all four giant planets; got {bodies_with_vinf}"
    )
