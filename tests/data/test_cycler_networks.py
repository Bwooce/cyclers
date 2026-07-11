"""Cycler-network registry schema + validator tests (task #570).

The committed registry (``data/cycler_networks.yaml``) ships GENUINELY
EMPTY in this pass -- no real, sourced network (Sanchez Net's or
otherwise) is populated (see that file's header + OUTSTANDING.md #570 /
#543). The schema + validator + taxi-cost cross-check are proven
end-to-end here with a SYNTHETIC self-test network built entirely in this
file's own fixtures, referencing REAL, already-catalogued
``data/catalogue.yaml`` rows. That network's ``source`` is ``"derived"``
and its ``notes`` explicitly states it is a schema self-test, never a
real discovered/sourced network -- it is intentionally NOT written into
the committed registry, so nobody can mistake it for a Sanchez-Net-style
finding later.
"""

from __future__ import annotations

import copy
from typing import Any

import pytest

from cyclerfinder.data.catalog import Catalog, load_catalog
from cyclerfinder.data.validate_networks import (
    NETWORK_SCHEMA_PATH,
    NETWORKS_PATH,
    _catalogue_entry_to_taxi_cycler,
    load_network_schema,
    load_networks_raw,
    validate_networks,
    validate_networks_referential,
    validate_networks_schema,
    validate_networks_semantic,
    validate_taxi_cost_cross_check,
)
from cyclerfinder.model.score import taxi_cost_kms

# Three REAL, already-catalogued, well-known V1+ rows (verified against the
# live catalogue, not hand-typed): the Aldrin classic outbound/inbound pair
# (V2/V1) and a Russell-Ocampo 2003/2004 short-transit variant (V1). All
# three carry a published Earth V_inf, so taxi_cost_kms('E') is non-zero
# for each -- a richer positive control than a body with no Earth encounter.
_MEMBER_IDS = (
    "aldrin-classic-em-k1-outbound",
    "aldrin-classic-em-k1-inbound",
    "russell-ocampo-2.5.1+0",
)


@pytest.fixture(scope="module")
def catalog() -> Catalog:
    return load_catalog()


def _build_self_test_network(catalog: Catalog) -> dict[str, Any]:
    """Build the synthetic self-test network dict.

    Every ``cost_kms`` value is computed by actually calling
    ``taxi_cost_kms`` on the referenced row's reconstructed
    :class:`Cycler` (via the same adapter the validator itself uses) --
    never hand-typed -- so this fixture cannot silently drift from the
    validator's own cross-check.
    """
    costs = []
    for cid in _MEMBER_IDS:
        entry = catalog.by_id[cid]
        cyc = _catalogue_entry_to_taxi_cycler(entry)
        cost = taxi_cost_kms(cyc, "E")
        costs.append({"cycler_id": cid, "cost_kms": cost, "taxi_body": "E"})

    return {
        "id": "schema-self-test-network-570",
        "name": "Schema self-test network (task #570 infrastructure proof)",
        "source": "derived",
        "member_cycler_ids": list(_MEMBER_IDS),
        "downlink_cadence": {
            "description": (
                "No real cadence claim -- this network exists only to prove the "
                "schema/validator/taxi-cost cross-check work end-to-end."
            ),
            "schedule": None,
        },
        "per_member_taxi_insertion_cost_kms": costs,
        "data_gaps": [
            {
                "path": "downlink_cadence.schedule",
                "kind": "not-applicable",
                "note": "synthetic self-test network; no real schedule exists to record.",
            }
        ],
        "first_published": None,
        "notes": (
            "SCHEMA SELF-TEST / INFRASTRUCTURE PROOF ONLY (task #570). NOT a real "
            "discovered or sourced cycler network -- do not treat as a Sanchez "
            "Net-style finding. Exists solely to exercise the cycler_network.schema.json "
            "+ validate_networks.py + taxi_cost_kms cross-check end-to-end in this test."
        ),
    }


# ---------------------------------------------------------------------------
# The committed registry ships empty and passes the combined gate
# ---------------------------------------------------------------------------


def test_committed_registry_is_empty() -> None:
    """Per the #570/#543 scoping decision, no real network is populated yet."""
    assert load_networks_raw(NETWORKS_PATH) == []


def test_committed_registry_passes_combined_gate(catalog: Catalog) -> None:
    networks = load_networks_raw(NETWORKS_PATH)
    errs = validate_networks(networks, catalog)
    assert errs == [], "validate_networks violations:\n" + "\n".join(errs)


def test_schema_file_loads_and_is_versioned() -> None:
    schema = load_network_schema(NETWORK_SCHEMA_PATH)
    assert schema["version"] == "1.0"
    assert schema["type"] == "array"


# ---------------------------------------------------------------------------
# Positive control: the synthetic self-test network passes the FULL gate
# ---------------------------------------------------------------------------


def test_self_test_network_passes_full_gate(catalog: Catalog) -> None:
    """Schema + referential + taxi-cost cross-check all pass end-to-end.

    This is the positive control: before trusting any negative-path
    assertion below, confirm the well-formed case clears every layer.
    """
    net = _build_self_test_network(catalog)
    errs = validate_networks([net], catalog)
    assert errs == [], "validate_networks violations on the self-test network:\n" + "\n".join(errs)


def test_self_test_network_is_labelled_derived_and_not_a_real_finding(catalog: Catalog) -> None:
    net = _build_self_test_network(catalog)
    assert net["source"] == "derived"
    assert "SCHEMA SELF-TEST" in net["notes"]
    assert "NOT a real discovered or sourced cycler network" in net["notes"]


def test_self_test_network_costs_are_nonzero_earth_insertion_costs(catalog: Catalog) -> None:
    """Sanity: all three chosen rows have a published Earth encounter, so
    their insertion cost is non-zero (exercises the real branch of
    taxi_cost_kms, not just the "no Earth encounter -> 0.0" fallback)."""
    net = _build_self_test_network(catalog)
    costs = {c["cycler_id"]: c["cost_kms"] for c in net["per_member_taxi_insertion_cost_kms"]}
    for cid in _MEMBER_IDS:
        assert costs[cid] > 0.0, f"{cid}: expected a non-zero Earth insertion cost"
    # Aldrin outbound/inbound are the same published cycler geometry
    # (Russell 2004 §3.8: "energy properties of inbound and outbound are
    # identical") -- their Earth V_inf, and therefore their insertion
    # cost, must match exactly.
    assert costs["aldrin-classic-em-k1-outbound"] == pytest.approx(
        costs["aldrin-classic-em-k1-inbound"]
    )


# ---------------------------------------------------------------------------
# Negative controls: each layer actually catches its violation
# ---------------------------------------------------------------------------


def test_unresolvable_member_id_flagged(catalog: Catalog) -> None:
    net = _build_self_test_network(catalog)
    net["member_cycler_ids"] = [*net["member_cycler_ids"], "not-a-real-catalogue-id"]
    errs = validate_networks_referential([net], catalog)
    assert any("not-a-real-catalogue-id" in e and "does not resolve" in e for e in errs), errs


def test_taxi_cost_drift_flagged(catalog: Catalog) -> None:
    """A stored cost_kms that disagrees with a fresh recomputation is a
    real validation failure, never silently tolerated."""
    net = _build_self_test_network(catalog)
    net["per_member_taxi_insertion_cost_kms"][0]["cost_kms"] += 5.0
    errs = validate_taxi_cost_cross_check([net], catalog)
    assert any("disagrees with a fresh taxi_cost_kms" in e for e in errs), errs


def test_taxi_cost_cross_check_is_a_real_recomputation_not_a_copy(catalog: Catalog) -> None:
    """Guards against a validator that merely echoes the stored value back
    (a no-op cross-check would never catch drift)."""
    net = _build_self_test_network(catalog)
    # Corrupt every cost_kms; the cross-check must flag ALL three, not zero.
    for item in net["per_member_taxi_insertion_cost_kms"]:
        item["cost_kms"] = -1.0
    errs = validate_taxi_cost_cross_check([net], catalog)
    assert len(errs) == len(_MEMBER_IDS), errs


def test_literature_source_without_first_published_flagged(catalog: Catalog) -> None:
    net = _build_self_test_network(catalog)
    net["source"] = "literature"
    net["first_published"] = None
    errs = validate_networks_semantic([net])
    assert any("requires a non-null first_published" in e for e in errs), errs


def test_literature_source_with_first_published_clean(catalog: Catalog) -> None:
    net = _build_self_test_network(catalog)
    net["source"] = "literature"
    net["first_published"] = {
        "authors": ["Sanchez Net, M."],
        "year": 2022,
        "title": "Earth-Mars Cycler Orbit fleet concept",
        "venue": "example venue",
    }
    errs = validate_networks_semantic([net])
    assert errs == [], errs


def test_duplicate_network_id_flagged(catalog: Catalog) -> None:
    net = _build_self_test_network(catalog)
    net2 = copy.deepcopy(net)
    errs = validate_networks_semantic([net, net2])
    assert any("duplicate cycler_networks id" in e for e in errs), errs


def test_cost_entry_cycler_id_not_a_member_flagged(catalog: Catalog) -> None:
    net = _build_self_test_network(catalog)
    net["per_member_taxi_insertion_cost_kms"].append(
        {"cycler_id": "not-a-member-id", "cost_kms": 1.0, "taxi_body": "E"}
    )
    errs = validate_networks_referential([net], catalog)
    assert any(
        "per_member_taxi_insertion_cost_kms references cycler_id 'not-a-member-id'" in e
        for e in errs
    ), errs


def test_schedule_cycler_id_not_a_member_flagged(catalog: Catalog) -> None:
    net = _build_self_test_network(catalog)
    net["downlink_cadence"] = {
        "description": "test",
        "schedule": [
            {
                "date": "2030-01-01T00:00:00",
                "cycler_id": "not-a-member-id",
                "me_transit_days": 150.0,
            }
        ],
    }
    errs = validate_networks_referential([net], catalog)
    assert any(
        "downlink_cadence.schedule references cycler_id 'not-a-member-id'" in e for e in errs
    ), errs


def test_schema_rejects_missing_required_fields() -> None:
    bad = {"name": "missing id/source/member_cycler_ids"}
    errs = validate_networks_schema([bad])
    assert any("required" in e for e in errs), errs


def test_schema_rejects_bad_source_enum() -> None:
    bad = {
        "id": "x",
        "name": "x",
        "source": "not-a-valid-source",
        "member_cycler_ids": ["aldrin-classic-em-k1-outbound"],
    }
    errs = validate_networks_schema([bad])
    assert any("enum" in e or "not one of" in e for e in errs), errs


def test_schema_rejects_empty_member_list() -> None:
    bad = {
        "id": "x",
        "name": "x",
        "source": "derived",
        "member_cycler_ids": [],
    }
    errs = validate_networks_schema([bad])
    assert any("non-empty" in e for e in errs), errs


def test_downlink_cadence_null_schedule_is_honest_not_an_error(catalog: Catalog) -> None:
    """A network with only a qualitative cadence description (no
    structured schedule) is valid -- matches this project's data_gap
    convention: absence records a known-unknown, not a defect."""
    net = _build_self_test_network(catalog)
    assert net["downlink_cadence"]["schedule"] is None
    errs = validate_networks([net], catalog)
    assert errs == [], errs
