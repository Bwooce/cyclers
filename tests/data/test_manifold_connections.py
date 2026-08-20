"""Manifold-connection registry schema + validator tests (task #838 design,
#856 implementation).

The committed registry (``data/manifold_connections.yaml``) ships with the
two REAL, retroactively-populated ``#822``/``#828`` row-touching
connections (transcribed, not recomputed, from
``data/found/822_vaquero_em_free_transfer/results.json`` and
``vaquero_em_cycler_connections.OVERLAP_GRID_ICS``). The schema + validator
are additionally proven end-to-end here with a SYNTHETIC self-test entry
built in this file's own fixtures, referencing a REAL, already-catalogued
``data/catalogue.yaml`` row plus a fabricated uncatalogued counterpart --
mirrors ``test_cycler_networks.py``'s pattern exactly.
"""

from __future__ import annotations

import copy
from typing import Any

import pytest

from cyclerfinder.data.catalog import Catalog, load_catalog
from cyclerfinder.data.validate_connections import (
    CONNECTION_SCHEMA_PATH,
    CONNECTIONS_PATH,
    load_connection_schema,
    load_connections_raw,
    validate_connections,
    validate_connections_referential,
    validate_connections_schema,
    validate_connections_semantic,
)

# A real, already-catalogued V1 row (independently verified against the live
# catalogue, not hand-typed) -- reused as the self-test entry's row_ref.
_REAL_ROW_ID = "vaquero-31-c254-em-cycler-2013"


@pytest.fixture(scope="module")
def catalog() -> Catalog:
    return load_catalog()


def _build_self_test_connection() -> dict[str, Any]:
    """Build a synthetic self-test connection dict, well-formed end-to-end."""
    return {
        "id": "schema-self-test-connection-838",
        "kind": "heteroclinic",
        "model": {
            "type": "cr3bp",
            "system": "Earth-Moon planar CR3BP (schema self-test)",
            "mass_ratio": 0.01215058439469525,
        },
        "jacobi_constant": 2.54,
        "endpoints": [
            {
                "uncatalogued": {
                    "family": "Schema self-test uncatalogued family member",
                    "x0": 1.0905363960533268,
                    "ydot0": -0.8231863180949408,
                    "period_nd": 5.941227735609639,
                    "jacobi_constant": 2.54,
                    "lambda_max": 3.149761131186353,
                    "derivation": "Schema self-test only; not a real computation.",
                }
            },
            {
                "row_ref": _REAL_ROW_ID,
                "identity_evidence": (
                    "SCHEMA SELF-TEST ONLY -- not a real identity claim; exists solely "
                    "to exercise the schema/validator end-to-end."
                ),
            },
        ],
        "connection": {
            "branch_u": -1,
            "branch_s": 1,
            "k_u": 40,
            "k_s": 27,
            "epsilon": 1.0e-4,
            "crossing_x": 0.23322507758762506,
            "crossing_xdot": -1.785596430399187,
            "newton_residual": 3.336901283531185e-10,
        },
        "evidence": {
            "full_state_gap": 4.483722190016192e-06,
            "dv_kms": 0.0,
        },
        "evidence_class": (
            "SCHEMA SELF-TEST / INFRASTRUCTURE PROOF ONLY (task #838/#856). NOT a real "
            "discovered or adjudicated connection -- do not treat as a finding."
        ),
        "dv_kms": 0.0,
        "provenance": {
            "task_refs": ["#838", "#856"],
            "data": "schema self-test; no real data artifact",
            "module": "src/cyclerfinder/data/validate_connections.py",
            "commit": None,
            "notes": None,
        },
    }


# ---------------------------------------------------------------------------
# The committed registry passes the combined gate
# ---------------------------------------------------------------------------


def test_committed_registry_has_the_three_entries() -> None:
    connections = load_connections_raw(CONNECTIONS_PATH)
    ids = {c["id"] for c in connections}
    assert ids == {
        "em-vaquero-hetero-wu21c254-ws31c254-2026",
        "em-vaquero-hetero-wu21c266-ws31c266-2026",
        "em-kumar-hetero-wu31c254-ws21c254-2026",
    }


def test_committed_registry_passes_combined_gate(catalog: Catalog) -> None:
    connections = load_connections_raw(CONNECTIONS_PATH)
    errs = validate_connections(connections, catalog)
    assert errs == [], "validate_connections violations:\n" + "\n".join(errs)


def test_schema_file_loads_and_is_versioned() -> None:
    schema = load_connection_schema(CONNECTION_SCHEMA_PATH)
    assert schema["version"] == "1.0"
    assert schema["type"] == "array"


def test_committed_entries_each_have_exactly_one_catalogued_endpoint(catalog: Catalog) -> None:
    """Task #838's motivating case is explicitly HALF-catalogued -- pin that
    the two retroactively-populated entries actually are, not accidentally
    both-catalogued or zero-catalogued."""
    connections = load_connections_raw(CONNECTIONS_PATH)
    for conn in connections:
        row_refs = [ep["row_ref"] for ep in conn["endpoints"] if "row_ref" in ep]
        uncatalogued = [ep for ep in conn["endpoints"] if "uncatalogued" in ep]
        assert len(row_refs) == 1, conn["id"]
        assert len(uncatalogued) == 1, conn["id"]
        assert row_refs[0] in catalog.by_id, conn["id"]


def test_committed_entries_transcribed_values_match_source_artifact() -> None:
    """Guards against a transcription slip: re-load #822's own results.json
    and confirm the registry's numbers are not hand-typed drift."""
    import json
    from pathlib import Path

    results_path = Path("data/found/822_vaquero_em_free_transfer/results.json")
    raw = json.loads(results_path.read_text())
    sweep_by_c = {round(row["jacobi"], 2): row for row in raw["sweep"]}

    connections = {c["id"]: c for c in load_connections_raw(CONNECTIONS_PATH)}
    c254 = connections["em-vaquero-hetero-wu21c254-ws31c254-2026"]
    src254 = sweep_by_c[2.54]
    assert c254["connection"]["newton_residual"] == src254["connection"]["residual"]
    assert c254["connection"]["crossing_x"] == src254["connection"]["crossing_xv"][0]
    assert c254["connection"]["crossing_xdot"] == src254["connection"]["crossing_xv"][1]
    assert c254["evidence"]["full_state_gap"] == src254["evidence"]["full_state_gap"]

    c266 = connections["em-vaquero-hetero-wu21c266-ws31c266-2026"]
    src266 = sweep_by_c[2.66]
    assert c266["connection"]["newton_residual"] == src266["connection"]["residual"]
    assert c266["connection"]["crossing_x"] == src266["connection"]["crossing_xv"][0]
    assert c266["connection"]["crossing_xdot"] == src266["connection"]["crossing_xv"][1]


def test_kumar_entry_reproduces_fresh_not_just_transcribed() -> None:
    """The #854 entry's numbers are re-derived FRESH here (not read back from
    any cached results.json), seeded only from Kumar et al.'s own printed
    Table-5 state -- the same non-circularity check #854's own adjudication
    ran, repeated as a permanent ratchet."""
    import cyclerfinder.search.kumar_em_resonant_heteroclinics as keh

    system = keh.kumar_system()
    fresh = keh.reproduce_table5_intersection(system, 2.54)
    assert fresh.connection is not None

    connections = {c["id"]: c for c in load_connections_raw(CONNECTIONS_PATH)}
    entry = connections["em-kumar-hetero-wu31c254-ws21c254-2026"]

    assert entry["connection"]["newton_residual"] == pytest.approx(
        fresh.connection.residual, rel=1e-9
    )
    assert entry["connection"]["crossing_x"] == pytest.approx(
        fresh.connection.crossing_xv[0], rel=1e-9
    )
    assert entry["connection"]["crossing_xdot"] == pytest.approx(
        fresh.connection.crossing_xv[1], rel=1e-9
    )
    assert entry["connection"]["k_u"] == fresh.connection.k_u
    assert entry["connection"]["k_s"] == fresh.connection.k_s
    assert entry["connection"]["branch_u"] == fresh.connection.branch_u
    assert entry["connection"]["branch_s"] == fresh.connection.branch_s
    assert fresh.matched is True
    assert fresh.match_distance == pytest.approx(7.610531887910121e-07, rel=1e-6)


# ---------------------------------------------------------------------------
# Positive control: the synthetic self-test entry passes the FULL gate
# ---------------------------------------------------------------------------


def test_self_test_connection_passes_full_gate(catalog: Catalog) -> None:
    """Positive control: confirm the well-formed case clears every layer
    before trusting any negative-path assertion below."""
    conn = _build_self_test_connection()
    errs = validate_connections([conn], catalog)
    assert errs == [], "validate_connections violations on the self-test entry:\n" + "\n".join(errs)


# ---------------------------------------------------------------------------
# Negative controls: each layer actually catches its violation
# ---------------------------------------------------------------------------


def test_schema_rejects_zero_catalogued_endpoints() -> None:
    """The admission criterion: at least one endpoint must carry row_ref."""
    conn = _build_self_test_connection()
    conn["endpoints"][1] = {
        "uncatalogued": {
            "family": "Second uncatalogued node -- no row_ref at all",
            "derivation": "Schema self-test only.",
        }
    }
    errs = validate_connections_schema([conn])
    assert errs != [], "expected the schema's 'contains: row_ref' rule to fire"


def test_schema_rejects_row_ref_without_identity_evidence() -> None:
    conn = _build_self_test_connection()
    del conn["endpoints"][1]["identity_evidence"]
    errs = validate_connections_schema([conn])
    # oneOf failures report generically ("not valid under any of the given
    # schemas") rather than naming the missing field -- just confirm the
    # gate actually fires, which the positive control proves it wouldn't
    # do spuriously.
    assert errs != [], "expected the row_ref branch's required identity_evidence to fire"


def test_schema_rejects_uncatalogued_without_derivation() -> None:
    conn = _build_self_test_connection()
    del conn["endpoints"][0]["uncatalogued"]["derivation"]
    errs = validate_connections_schema([conn])
    assert errs != [], "expected the uncatalogued branch's required derivation to fire"


def test_schema_rejects_bad_kind_enum() -> None:
    conn = _build_self_test_connection()
    conn["kind"] = "not-a-valid-kind"
    errs = validate_connections_schema([conn])
    assert any("enum" in e or "not one of" in e for e in errs), errs


def test_schema_rejects_missing_required_top_level_fields() -> None:
    bad = {"id": "x", "kind": "heteroclinic"}
    errs = validate_connections_schema([bad])
    assert any("required" in e for e in errs), errs


def test_schema_rejects_wrong_endpoint_count() -> None:
    conn = _build_self_test_connection()
    conn["endpoints"] = [conn["endpoints"][0]]
    errs = validate_connections_schema([conn])
    assert errs != [], "expected minItems/maxItems=2 to fire"


def test_unresolvable_row_ref_flagged(catalog: Catalog) -> None:
    conn = _build_self_test_connection()
    conn["endpoints"][1]["row_ref"] = "not-a-real-catalogue-id"
    errs = validate_connections_referential([conn], catalog)
    assert any("not-a-real-catalogue-id" in e and "does not resolve" in e for e in errs), errs


def test_blank_identity_evidence_flagged_by_semantic_layer(catalog: Catalog) -> None:
    """A whitespace-only identity_evidence satisfies the schema's minLength:1
    but must still be caught -- guards the semantic layer's own value."""
    conn = _build_self_test_connection()
    conn["endpoints"][1]["identity_evidence"] = "   "
    errs = validate_connections_semantic([conn])
    assert any("blank/missing identity_evidence" in e for e in errs), errs


def test_duplicate_connection_id_flagged() -> None:
    conn = _build_self_test_connection()
    conn2 = copy.deepcopy(conn)
    errs = validate_connections_semantic([conn, conn2])
    assert any("duplicate manifold_connections id" in e for e in errs), errs


def test_unresolvable_reverse_of_flagged(catalog: Catalog) -> None:
    conn = _build_self_test_connection()
    conn["reverse_of"] = "not-a-real-connection-id"
    errs = validate_connections_referential([conn], catalog)
    assert any("not-a-real-connection-id" in e and "does not resolve" in e for e in errs), errs


def test_resolvable_reverse_of_clean(catalog: Catalog) -> None:
    conn = _build_self_test_connection()
    conn2 = _build_self_test_connection()
    conn2["id"] = "schema-self-test-connection-838-reverse"
    conn2["reverse_of"] = conn["id"]
    errs = validate_connections_referential([conn, conn2], catalog)
    assert errs == [], errs
