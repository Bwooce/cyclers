"""Schema-v4 semantic/cross-field validation gate (spec §16.7, Task 1.3).

Covers rules that JSON Schema (Task 1.4) cannot express: cross-field
semantics and census ratchets.  Returns a list of error strings (never
raises), so callers can collect all violations at once.

Rules enforced
--------------
1. ``cycler_class in {multi-arc, non-keplerian}`` ⇒ ``orbit_elements.a_au``
   and ``orbit_elements.e`` must be absent or null (no single-ellipse
   conic element at the top-level orbit block).
2. ``cycler_class == non-keplerian`` ⇒ ``primary != "Sun"`` (CR3BP /
   rotating-frame orbits are planet-centric, not heliocentric).
3. ``period.basis`` items, when present, must each carry both ``pair``
   and ``k`` keys.
"""

from __future__ import annotations

from typing import Any

_MULTI_CLASS = {"multi-arc", "non-keplerian"}


def validate_schema_invariants(rows: list[dict[str, Any]]) -> list[str]:
    """Validate schema-v4 semantic invariants across a list of raw YAML row dicts.

    Parameters
    ----------
    rows:
        List of raw row dicts (as loaded by ``yaml.safe_load`` on
        ``data/catalogue.yaml``).

    Returns
    -------
    list[str]
        All violation messages found.  Empty list means the data is clean.
        Does NOT raise — callers decide how to handle violations.
    """
    errors: list[str] = []
    for row in rows:
        rid = row.get("id") or "<unknown>"
        cls = str(row.get("cycler_class") or "single-ellipse")
        oe = row.get("orbit_elements") or {}
        period = row.get("period") or {}

        # Rule 1: multi-arc / non-keplerian must not carry top-level a/e
        if cls in _MULTI_CLASS:
            a_au = oe.get("a_au")
            e = oe.get("e")
            if a_au is not None:
                errors.append(
                    f"{rid}: cycler_class={cls!r} must not have orbit_elements.a_au "
                    f"(got {a_au!r}); top-level a/e is only valid for single-ellipse"
                )
            if e is not None:
                errors.append(
                    f"{rid}: cycler_class={cls!r} must not have orbit_elements.e "
                    f"(got {e!r}); top-level a/e is only valid for single-ellipse"
                )

        # Rule 2: non-keplerian implies primary != "Sun"
        if cls == "non-keplerian":
            primary = str(row.get("primary") or "Sun")
            if primary == "Sun":
                errors.append(
                    f"{rid}: cycler_class=non-keplerian requires a non-Sun primary "
                    f"(CR3BP orbit is planet-centric); got primary={primary!r}"
                )

        # Rule 3: period.basis items must have pair + k
        basis = period.get("basis")
        if basis is not None:
            for i, item in enumerate(basis):
                if not isinstance(item, dict):
                    errors.append(
                        f"{rid}: period.basis[{i}] must be a dict with 'pair' and 'k', "
                        f"got {type(item).__name__!r}"
                    )
                    continue
                if "pair" not in item:
                    errors.append(f"{rid}: period.basis[{i}] is missing required key 'pair'")
                if "k" not in item:
                    errors.append(f"{rid}: period.basis[{i}] is missing required key 'k'")

    return errors


__all__ = ["validate_schema_invariants"]
