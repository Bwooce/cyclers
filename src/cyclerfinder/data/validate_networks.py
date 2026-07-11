"""Cycler-network registry validator (task #570, #543-scoped first slice).

Mirrors :mod:`cyclerfinder.data.validate`'s pattern: every ``validate_*``
function returns a list of violation strings and never raises, so callers
(tests, CI, a future CLI) can collect every violation at once rather than
stopping at the first.

Layers
------
1. :func:`validate_networks_schema` -- JSON-Schema structural validation
   against ``data/cycler_network.schema.json`` (in-Python, via the
   ``jsonschema`` package -- unlike ``data/catalogue.yaml``, which is
   schema-checked out-of-process by the ``check-jsonschema`` pre-commit
   hook, this module also runs the schema check itself so
   :func:`validate_networks` is a single self-contained gate).
2. :func:`validate_networks_semantic` -- cross-field rules JSON Schema
   cannot express cheaply: entry-id uniqueness, and
   ``source == "literature"`` requiring a non-null ``first_published``.
3. :func:`validate_networks_referential` -- every ``member_cycler_ids``
   entry MUST resolve to a real, existing ``data/catalogue.yaml`` row id
   (never silently skipped), plus the same check for every
   ``cycler_id`` referenced inside ``downlink_cadence.schedule`` and
   ``per_member_taxi_insertion_cost_kms`` (those must also be declared
   members).
4. :func:`validate_taxi_cost_cross_check` -- for every
   ``per_member_taxi_insertion_cost_kms`` entry, reconstructs the
   referenced catalogue row as a minimal
   :class:`~cyclerfinder.model.cycler.Cycler` (see
   :func:`_catalogue_entry_to_taxi_cycler`) and calls
   :func:`cyclerfinder.model.score.taxi_cost_kms` fresh, then asserts
   the stored ``cost_kms`` matches within a tight tolerance. This is a
   cross-check against drift (the referenced row's encounters may
   change after the network entry was written), never a duplicate
   source of truth -- a mismatch is a real validation failure.

Design note carried over from task #570's scoping (verify independently,
don't just trust this comment): ``taxi_cost_kms`` computes a per-cycler
INSERTION cost (Earth -> that single cycler's own Earth encounter), NOT a
cost between two different cyclers. This registry therefore only ever
stores one insertion cost per member, keyed to that member's own row --
never an inter-cycler transfer cost (no such computation exists in this
codebase; see ``data/cycler_network.schema.json``'s top-level
description).
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any, Final

import numpy as np
import yaml  # type: ignore[import-untyped]

from cyclerfinder.data.catalog import Catalog, CatalogueEntry, load_catalog
from cyclerfinder.model.cycler import Cycler, Encounter
from cyclerfinder.model.score import taxi_cost_kms

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

NETWORKS_PATH: Final[Path] = (
    Path(__file__).resolve().parent.parent.parent.parent / "data" / "cycler_networks.yaml"
)
"""Resolved path to ``data/cycler_networks.yaml`` (mirrors
:data:`cyclerfinder.data.catalog.CATALOGUE_PATH`'s pattern)."""

NETWORK_SCHEMA_PATH: Final[Path] = (
    Path(__file__).resolve().parent.parent.parent.parent / "data" / "cycler_network.schema.json"
)
"""Resolved path to ``data/cycler_network.schema.json``."""

_TAXI_COST_TOL: Final[float] = 1e-9
"""Tolerance (both relative and absolute) for the stored-vs-fresh
``taxi_cost_kms`` cross-check. Tight because both sides compute the exact
same deterministic function of the exact same published V∞ magnitude --
any discrepancy above float round-off is drift, not noise."""

_DATA_GAP_KINDS: Final[frozenset[str]] = frozenset(
    {
        "unknown",
        "derive",
        "uncertain",
        "not-applicable",
        "unverified-provenance",
        "resolved",
        "conflict",
        "anomaly",
    }
)


# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------


def load_networks_raw(path: Path | str = NETWORKS_PATH) -> list[dict[str, Any]]:
    """Load ``data/cycler_networks.yaml`` as a list of raw dicts.

    Returns an empty list for a genuinely-empty (``[]`` or ``null``)
    registry -- the expected state for this task's pass (see that file's
    header comment).
    """
    raw = yaml.safe_load(Path(path).read_text())
    return list(raw) if raw else []


def load_network_schema(path: Path | str = NETWORK_SCHEMA_PATH) -> dict[str, Any]:
    """Load ``data/cycler_network.schema.json`` as a dict."""
    return dict(json.loads(Path(path).read_text()))


# ---------------------------------------------------------------------------
# Layer 1: JSON-Schema structural validation
# ---------------------------------------------------------------------------


def validate_networks_schema(
    networks: list[dict[str, Any]],
    *,
    schema: dict[str, Any] | None = None,
) -> list[str]:
    """Validate *networks* against ``data/cycler_network.schema.json``.

    Uses the ``jsonschema`` package directly (already a project dev
    dependency, see ``pyproject.toml``) rather than shelling out to
    ``check-jsonschema``, so this function is self-contained and callable
    from pytest without a subprocess.

    Returns every violation message (never raises); empty means clean.
    """
    from jsonschema import Draft202012Validator

    if schema is None:
        schema = load_network_schema()
    validator = Draft202012Validator(schema)
    errors: list[str] = []
    for err in sorted(validator.iter_errors(networks), key=lambda e: [str(p) for p in e.path]):
        loc = "/".join(str(p) for p in err.path) or "<root>"
        errors.append(f"schema: {loc}: {err.message}")
    return errors


# ---------------------------------------------------------------------------
# Layer 2: cross-field semantic rules JSON Schema cannot express
# ---------------------------------------------------------------------------


def validate_networks_semantic(networks: list[dict[str, Any]]) -> list[str]:
    """Cross-field rules: id uniqueness, and literature-source honesty.

    * Every ``id`` must be unique across the registry (JSON Schema has no
      "unique value of a property across array items" primitive).
    * ``source == "literature"`` requires a non-null ``first_published``
      block -- mirrors ``data/catalogue.schema.json``'s
      ``dv_band``-requires-``dv_band_source`` honesty-gate convention:
      every literature claim traces to a citation, never inferred.
    """
    errors: list[str] = []
    seen_ids: set[str] = set()
    for net in networks:
        net_id = str(net.get("id") or "<unknown>")
        if net_id in seen_ids:
            errors.append(f"duplicate cycler_networks id {net_id!r}")
        seen_ids.add(net_id)

        source = net.get("source")
        if source == "literature" and not net.get("first_published"):
            errors.append(
                f"{net_id}: source='literature' requires a non-null first_published block"
            )

        for gap in net.get("data_gaps") or []:
            if not isinstance(gap, dict):
                continue
            kind = gap.get("kind")
            if kind is not None and kind not in _DATA_GAP_KINDS:
                errors.append(f"{net_id}: data_gaps[].kind={kind!r} is not a recognised kind")

    return errors


# ---------------------------------------------------------------------------
# Layer 3: referential integrity against data/catalogue.yaml
# ---------------------------------------------------------------------------


def validate_networks_referential(
    networks: list[dict[str, Any]],
    catalog: Catalog,
) -> list[str]:
    """Every ``member_cycler_ids`` entry must resolve to a real catalogue row.

    Also checks that every ``cycler_id`` referenced inside
    ``downlink_cadence.schedule`` and ``per_member_taxi_insertion_cost_kms``
    is itself a declared member (a network cannot cite a downlink date or
    insertion cost for a cycler it doesn't claim as a member) -- this
    is stricter than the OUTSTANDING #570 minimum spec but a cheap,
    obviously-correct extension of the same referential-integrity rule.

    Never silently skips an unresolvable reference: every miss becomes a
    violation string.
    """
    errors: list[str] = []
    for net in networks:
        net_id = str(net.get("id") or "<unknown>")
        members_raw = net.get("member_cycler_ids") or []
        member_set: set[str] = set()
        for m in members_raw:
            if not isinstance(m, str) or not m:
                errors.append(
                    f"{net_id}: member_cycler_ids entries must be non-empty strings (got {m!r})"
                )
                continue
            member_set.add(m)
            if m not in catalog.by_id:
                errors.append(
                    f"{net_id}: member_cycler_ids entry {m!r} does not resolve to an "
                    f"existing data/catalogue.yaml row id"
                )

        schedule = ((net.get("downlink_cadence") or {}) or {}).get("schedule") or []
        for sched in schedule:
            if not isinstance(sched, dict):
                continue
            cid = sched.get("cycler_id")
            if cid is not None and cid not in member_set:
                errors.append(
                    f"{net_id}: downlink_cadence.schedule references cycler_id {cid!r} "
                    f"not listed in member_cycler_ids"
                )

        for item in net.get("per_member_taxi_insertion_cost_kms") or []:
            if not isinstance(item, dict):
                continue
            cid = item.get("cycler_id")
            if cid is not None and cid not in member_set:
                errors.append(
                    f"{net_id}: per_member_taxi_insertion_cost_kms references cycler_id "
                    f"{cid!r} not listed in member_cycler_ids"
                )

    return errors


# ---------------------------------------------------------------------------
# Layer 4: taxi-cost cross-check
# ---------------------------------------------------------------------------


def _catalogue_entry_to_taxi_cycler(entry: CatalogueEntry) -> Cycler:
    """Build the MINIMAL :class:`Cycler` sufficient for :func:`taxi_cost_kms`.

    ``taxi_cost_kms(cycler, taxi_body)`` reads exactly two things off each
    encounter -- ``.body`` and the NORM of ``.vinf_in`` (see
    ``model/score.py::taxi_cost_kms``: ``max(||enc.vinf_in||)`` over
    encounters whose body is ``taxi_body``). It never reads ``.t``,
    ``.r``, ``.v_planet``, ``.vinf_out``, or any :class:`Leg` field, and
    ``Cycler.period``/``Cycler.sense`` are likewise untouched by this
    function.

    ``data/catalog.py``'s :class:`CatalogueEntry` deliberately stores only
    the catalogue's PUBLISHED SCALAR ``vinf_kms_at_encounters`` magnitudes
    (a ``(body, float | None)`` tuple), not a full 3-vector spacecraft
    state -- no richer catalogue-row -> :class:`Cycler` adapter exists
    anywhere else in this codebase (searched: ``data/catalog.py``,
    ``verify/plausibility.py``, ``search/precursor_matcher.py`` all consume
    the scalar tuples directly, never a reconstructed :class:`Cycler`).

    This adapter places each published V∞ magnitude along an arbitrary
    fixed axis (direction is provably irrelevant here: ``taxi_cost_kms``
    only ever calls ``np.linalg.norm`` on the vector) to get the EXACT
    same norm ``taxi_cost_kms`` would compute from the real vector state.
    It is intentionally narrow -- scoped to the one function this
    validator cross-checks -- not a general-purpose catalogue-to-Cycler
    constructor (that would need real per-encounter geometry this
    dataclass does not carry).

    Encounters whose published V∞ is ``None`` (e.g. a CR3BP row's
    Jacobi-constant-only entries, see ``ross-rt-em-cycler-11-2025``) are
    skipped -- ``taxi_cost_kms`` cannot evaluate an unpublished magnitude,
    and skipping matches its own "no vinfs -> 0.0" empty-list behaviour.
    """
    encounters: list[Encounter] = []
    for body, vinf_kms in entry.vinf_kms_at_encounters:
        if vinf_kms is None:
            continue
        vinf_vec = np.array([float(vinf_kms), 0.0, 0.0])
        encounters.append(
            Encounter(
                body=body,
                t=0.0,
                r=np.zeros(3),
                v_planet=np.zeros(3),
                vinf_in=vinf_vec,
                vinf_out=vinf_vec,
            )
        )
    return Cycler(
        bodies=list(entry.bodies),
        period=0.0,  # unused by taxi_cost_kms
        encounters=encounters,
        legs=[],
        sense=entry.sense,
    )


def validate_taxi_cost_cross_check(
    networks: list[dict[str, Any]],
    catalog: Catalog,
) -> list[str]:
    """Re-derive every stored ``cost_kms`` fresh and assert it matches.

    Looks up each referenced catalogue row, reconstructs it as a minimal
    :class:`Cycler` (:func:`_catalogue_entry_to_taxi_cycler`), calls
    :func:`taxi_cost_kms` fresh, and compares against the stored value
    within :data:`_TAXI_COST_TOL`. A mismatch is a REAL validation
    failure (the row's encounters drifted since the network entry was
    written) -- never silently tolerated.
    """
    errors: list[str] = []
    for net in networks:
        net_id = str(net.get("id") or "<unknown>")
        for item in net.get("per_member_taxi_insertion_cost_kms") or []:
            if not isinstance(item, dict):
                errors.append(
                    f"{net_id}: per_member_taxi_insertion_cost_kms entry must be a dict "
                    f"(got {item!r})"
                )
                continue
            cycler_id = item.get("cycler_id")
            stored_cost = item.get("cost_kms")
            taxi_body = item.get("taxi_body")
            if cycler_id is None or stored_cost is None or taxi_body is None:
                errors.append(
                    f"{net_id}: per_member_taxi_insertion_cost_kms entry missing "
                    f"cycler_id/cost_kms/taxi_body (got {item!r})"
                )
                continue
            entry = catalog.by_id.get(str(cycler_id))
            if entry is None:
                # Already reported by validate_networks_referential; skip
                # re-flagging here rather than double-reporting the same
                # unresolvable id under a different message.
                continue
            cyc = _catalogue_entry_to_taxi_cycler(entry)
            fresh_cost = taxi_cost_kms(cyc, str(taxi_body))
            if not math.isclose(
                float(stored_cost), fresh_cost, rel_tol=_TAXI_COST_TOL, abs_tol=_TAXI_COST_TOL
            ):
                errors.append(
                    f"{net_id}: stored cost_kms={stored_cost!r} for cycler_id={cycler_id!r} "
                    f"taxi_body={taxi_body!r} disagrees with a fresh taxi_cost_kms "
                    f"recomputation ({fresh_cost!r}) -- drift, not tolerated (this is a "
                    f"cross-check against the live catalogue, not a duplicate source of truth)"
                )
    return errors


# ---------------------------------------------------------------------------
# Combined gate
# ---------------------------------------------------------------------------


def validate_networks(
    networks: list[dict[str, Any]],
    catalog: Catalog,
    *,
    schema: dict[str, Any] | None = None,
) -> list[str]:
    """Run all four layers over *networks*, returning every violation.

    Single combined entry point, mirroring
    :func:`cyclerfinder.data.validate.validate_catalogue`'s pattern.
    Never raises -- callers decide how to surface violations.
    """
    return (
        validate_networks_schema(networks, schema=schema)
        + validate_networks_semantic(networks)
        + validate_networks_referential(networks, catalog)
        + validate_taxi_cost_cross_check(networks, catalog)
    )


__all__ = [
    "NETWORKS_PATH",
    "NETWORK_SCHEMA_PATH",
    "load_network_schema",
    "load_networks_raw",
    "validate_networks",
    "validate_networks_referential",
    "validate_networks_schema",
    "validate_networks_semantic",
    "validate_taxi_cost_cross_check",
]


# ---------------------------------------------------------------------------
# Minimal CLI (mirrors data/catalog.py's ``_main`` pattern)
# ---------------------------------------------------------------------------


def _main(argv: list[str] | None = None) -> int:
    """``python -m cyclerfinder.data.validate_networks check`` — validate
    ``data/cycler_networks.yaml`` against ``data/catalogue.yaml`` and print
    every violation (empty output + exit 0 means clean)."""
    import argparse

    parser = argparse.ArgumentParser(prog="cyclerfinder.data.validate_networks")
    sub = parser.add_subparsers(dest="cmd", required=True)
    check = sub.add_parser("check", help="validate data/cycler_networks.yaml")
    check.add_argument("--networks", default=str(NETWORKS_PATH))
    check.add_argument("--catalogue", default=None)
    args = parser.parse_args(argv)

    if args.cmd == "check":
        networks = load_networks_raw(args.networks)
        catalog = load_catalog(args.catalogue) if args.catalogue else load_catalog()
        errors = validate_networks(networks, catalog)
        for e in errors:
            print(e)
        return 1 if errors else 0
    return 0


if __name__ == "__main__":
    raise SystemExit(_main())
