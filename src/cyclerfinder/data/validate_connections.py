"""Manifold-connection registry validator (task #838 design, #856 implementation).

Mirrors :mod:`cyclerfinder.data.validate_networks`'s pattern exactly (which
itself mirrors :mod:`cyclerfinder.data.validate`): every ``validate_*``
function returns a list of violation strings and never raises, so callers
(tests, CI, a future CLI) can collect every violation at once rather than
stopping at the first.

Layers
------
1. :func:`validate_connections_schema` -- JSON-Schema structural validation
   against ``data/manifold_connection.schema.json`` (in-Python via
   ``jsonschema``, so this module is a single self-contained gate; the
   ``check-jsonschema`` pre-commit hook additionally runs the same schema
   out-of-process, matching ``data/cycler_networks.yaml``'s convention).
   The schema's own ``endpoints`` ``contains`` constraint already enforces
   "at least one endpoint carries row_ref" structurally.
2. :func:`validate_connections_semantic` -- cross-field rules JSON Schema
   cannot express cheaply: entry-id uniqueness, and (defensively, since the
   schema's ``oneOf``/``minLength`` already cover most of this) that every
   ``row_ref`` endpoint's ``identity_evidence`` is genuinely non-blank
   (not just whitespace).
3. :func:`validate_connections_referential` -- every ``row_ref`` MUST
   resolve to a real, existing ``data/catalogue.yaml`` row id (never
   silently skipped), and every ``reverse_of`` MUST resolve to a real id
   *within this same registry*.

Scope note (see the schema's own top-level description for the full
rationale): this registry is for EXTRINSIC connections only -- a
transport statement BETWEEN two orbits, at least one of which is a
catalogued row. It never writes back to ``data/catalogue.yaml`` and never
changes a row's ``validation_level``/``our_status``/``orbit_class`` --
"provenance/audit only, not a promotion gate", the same rule already
stated on ``ccr4bp_provenance.connection`` and
``crnbp_provenance...seed_orbit_homoclinic.connection``.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Final

import yaml  # type: ignore[import-untyped]

from cyclerfinder.data.catalog import Catalog, load_catalog

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

CONNECTIONS_PATH: Final[Path] = (
    Path(__file__).resolve().parent.parent.parent.parent / "data" / "manifold_connections.yaml"
)
"""Resolved path to ``data/manifold_connections.yaml`` (mirrors
:data:`cyclerfinder.data.validate_networks.NETWORKS_PATH`'s pattern)."""

CONNECTION_SCHEMA_PATH: Final[Path] = (
    Path(__file__).resolve().parent.parent.parent.parent
    / "data"
    / "manifold_connection.schema.json"
)
"""Resolved path to ``data/manifold_connection.schema.json``."""


# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------


def load_connections_raw(path: Path | str = CONNECTIONS_PATH) -> list[dict[str, Any]]:
    """Load ``data/manifold_connections.yaml`` as a list of raw dicts.

    Returns an empty list for a genuinely-empty (``[]`` or ``null``)
    registry.
    """
    raw = yaml.safe_load(Path(path).read_text())
    return list(raw) if raw else []


def load_connection_schema(path: Path | str = CONNECTION_SCHEMA_PATH) -> dict[str, Any]:
    """Load ``data/manifold_connection.schema.json`` as a dict."""
    return dict(json.loads(Path(path).read_text()))


# ---------------------------------------------------------------------------
# Layer 1: JSON-Schema structural validation
# ---------------------------------------------------------------------------


def validate_connections_schema(
    connections: list[dict[str, Any]],
    *,
    schema: dict[str, Any] | None = None,
) -> list[str]:
    """Validate *connections* against ``data/manifold_connection.schema.json``.

    Returns every violation message (never raises); empty means clean.
    """
    from jsonschema import Draft202012Validator

    if schema is None:
        schema = load_connection_schema()
    validator = Draft202012Validator(schema)
    errors: list[str] = []
    for err in sorted(validator.iter_errors(connections), key=lambda e: [str(p) for p in e.path]):
        loc = "/".join(str(p) for p in err.path) or "<root>"
        errors.append(f"schema: {loc}: {err.message}")
    return errors


# ---------------------------------------------------------------------------
# Layer 2: cross-field semantic rules JSON Schema cannot express
# ---------------------------------------------------------------------------


def validate_connections_semantic(connections: list[dict[str, Any]]) -> list[str]:
    """Cross-field rules: id uniqueness, and non-blank identity evidence.

    * Every ``id`` must be unique across the registry (JSON Schema has no
      "unique value of a property across array items" primitive).
    * Every ``row_ref`` endpoint's ``identity_evidence`` must contain real
      text once whitespace is stripped -- a defensive check beyond the
      schema's ``minLength: 1`` (which a lone space would satisfy).
    """
    errors: list[str] = []
    seen_ids: set[str] = set()
    for conn in connections:
        conn_id = str(conn.get("id") or "<unknown>")
        if conn_id in seen_ids:
            errors.append(f"duplicate manifold_connections id {conn_id!r}")
        seen_ids.add(conn_id)

        for i, ep in enumerate(conn.get("endpoints") or []):
            if not isinstance(ep, dict):
                continue
            if "row_ref" in ep:
                evidence = ep.get("identity_evidence")
                if not isinstance(evidence, str) or not evidence.strip():
                    errors.append(
                        f"{conn_id}: endpoints[{i}] row_ref {ep.get('row_ref')!r} has "
                        f"blank/missing identity_evidence"
                    )

    return errors


# ---------------------------------------------------------------------------
# Layer 3: referential integrity
# ---------------------------------------------------------------------------


def validate_connections_referential(
    connections: list[dict[str, Any]],
    catalog: Catalog,
) -> list[str]:
    """Every ``row_ref`` must resolve to a real catalogue row; every
    ``reverse_of`` must resolve to a real id within this registry.

    Never silently skips an unresolvable reference: every miss becomes a
    violation string.
    """
    errors: list[str] = []
    all_ids = {str(c.get("id")) for c in connections if c.get("id")}

    for conn in connections:
        conn_id = str(conn.get("id") or "<unknown>")

        for i, ep in enumerate(conn.get("endpoints") or []):
            if not isinstance(ep, dict):
                continue
            row_ref = ep.get("row_ref")
            if row_ref is not None and row_ref not in catalog.by_id:
                errors.append(
                    f"{conn_id}: endpoints[{i}] row_ref {row_ref!r} does not resolve to an "
                    f"existing data/catalogue.yaml row id"
                )

        reverse_of = conn.get("reverse_of")
        if reverse_of is not None and reverse_of not in all_ids:
            errors.append(
                f"{conn_id}: reverse_of {reverse_of!r} does not resolve to an existing "
                f"data/manifold_connections.yaml entry id"
            )

    return errors


# ---------------------------------------------------------------------------
# Combined gate
# ---------------------------------------------------------------------------


def validate_connections(
    connections: list[dict[str, Any]],
    catalog: Catalog,
    *,
    schema: dict[str, Any] | None = None,
) -> list[str]:
    """Run all three layers over *connections*, returning every violation.

    Single combined entry point, mirroring
    :func:`cyclerfinder.data.validate_networks.validate_networks`'s
    pattern. Never raises -- callers decide how to surface violations.
    """
    return (
        validate_connections_schema(connections, schema=schema)
        + validate_connections_semantic(connections)
        + validate_connections_referential(connections, catalog)
    )


__all__ = [
    "CONNECTIONS_PATH",
    "CONNECTION_SCHEMA_PATH",
    "load_connection_schema",
    "load_connections_raw",
    "validate_connections",
    "validate_connections_referential",
    "validate_connections_schema",
    "validate_connections_semantic",
]


# ---------------------------------------------------------------------------
# Minimal CLI (mirrors validate_networks.py's ``_main`` pattern)
# ---------------------------------------------------------------------------


def _main(argv: list[str] | None = None) -> int:
    """``python -m cyclerfinder.data.validate_connections check`` --
    validate ``data/manifold_connections.yaml`` against
    ``data/catalogue.yaml`` and print every violation (empty output +
    exit 0 means clean)."""
    import argparse

    parser = argparse.ArgumentParser(prog="cyclerfinder.data.validate_connections")
    sub = parser.add_subparsers(dest="cmd", required=True)
    check = sub.add_parser("check", help="validate data/manifold_connections.yaml")
    check.add_argument("--connections", default=str(CONNECTIONS_PATH))
    check.add_argument("--catalogue", default=None)
    args = parser.parse_args(argv)

    if args.cmd == "check":
        connections = load_connections_raw(args.connections)
        catalog = load_catalog(args.catalogue) if args.catalogue else load_catalog()
        errors = validate_connections(connections, catalog)
        for e in errors:
            print(e)
        return 1 if errors else 0
    return 0


if __name__ == "__main__":
    raise SystemExit(_main())
