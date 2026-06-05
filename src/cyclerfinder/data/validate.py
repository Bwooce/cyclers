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
4. (v4.2, spec §16.7.9) ``trajectory.segments[].tof_days_bounds``, when
   present, must be exactly two positive numbers with ``min <= max``.  It
   is deliberately NOT required to contain ``tof_days``.
5. (v4.2, spec §16.7.9) ``source_ephemeris``, when present, must be a
   non-empty string.

Dispatch helper
---------------
:func:`anchors_for` returns which expected-output anchor categories apply
to a given catalogue entry dict, dispatching by ``cycler_class`` per
spec §16.7.5.  Callers should use this rather than hard-coding class
checks inline.
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

        # Rule 4 (v4.2): trajectory.segments[].tof_days_bounds shape/order.
        # Note: tof_days (when present) is deliberately NOT required to lie
        # inside the bounds — e.g. the Aldrin outbound segment carries
        # tof_days=146 (circular-coplanar idealization) while Rogers et al.
        # 2012 Table 4 publishes a STOUR range of 161-172 d; both are sourced
        # framings of the same leg, so forcing containment would reject valid data.
        trajectory = row.get("trajectory") or {}
        segments = trajectory.get("segments")
        if segments is not None:
            for i, seg in enumerate(segments):
                if not isinstance(seg, dict):
                    continue
                bounds = seg.get("tof_days_bounds")
                if bounds is None:
                    continue
                seg_id = seg.get("id", i)
                if not isinstance(bounds, (list, tuple)) or len(bounds) != 2:
                    errors.append(
                        f"{rid}: trajectory.segments[{seg_id}].tof_days_bounds must be "
                        f"exactly 2 numbers [min, max] (got {bounds!r})"
                    )
                    continue
                lo, hi = bounds
                if (
                    not isinstance(lo, (int, float))
                    or not isinstance(hi, (int, float))
                    or isinstance(lo, bool)
                    or isinstance(hi, bool)
                ):
                    errors.append(
                        f"{rid}: trajectory.segments[{seg_id}].tof_days_bounds items must "
                        f"both be numbers (got {bounds!r})"
                    )
                    continue
                if lo <= 0 or hi <= 0:
                    errors.append(
                        f"{rid}: trajectory.segments[{seg_id}].tof_days_bounds values must "
                        f"both be > 0 (got {bounds!r})"
                    )
                if lo > hi:
                    errors.append(
                        f"{rid}: trajectory.segments[{seg_id}].tof_days_bounds must have "
                        f"min <= max (got min={lo!r}, max={hi!r})"
                    )

        # Rule 5 (v4.2): source_ephemeris, when present, must be a non-empty string.
        if "source_ephemeris" in row:
            se = row.get("source_ephemeris")
            if se is not None and (not isinstance(se, str) or se.strip() == ""):
                errors.append(
                    f"{rid}: source_ephemeris must be a non-empty string when present (got {se!r})"
                )

    return errors


def anchors_for(entry: dict[str, Any]) -> dict[str, bool]:
    """Return which expected-output anchor categories apply to *entry*.

    Dispatches by ``cycler_class`` per spec §16.7.5.  The returned dict
    maps anchor-category names to ``True`` when that category applies and
    ``False`` when it must not be applied for the given class.

    Anchor categories
    -----------------
    ``"vinf"``
        V∞ multiset check — applicable to all classes (where
        ``vinf_kms_at_encounters`` is populated).
    ``"a_e"``
        Semi-major axis + eccentricity identity — applies only to
        ``single-ellipse``.  Must NOT be applied to multi-arc or
        non-keplerian entries (they have no single conic; the
        ``construct_resonant_cycler`` constructor must not be invoked for
        them).
    ``"period"``
        Period (years / k) check — applies to ``single-ellipse`` and
        ``multi-arc``; not meaningful for ``non-keplerian`` entries whose
        period is expressed in non-dimensional CR3BP time.
    ``"invariants"``
        Cycle-level identity (aphelion_ratio / turn_ratio /
        transit_times_days) — applies only to ``multi-arc``.
    ``"cr3bp"``
        CR3BP identity triple (jacobi_constant / period_nd /
        stability_index) — applies only to ``non-keplerian``.

    Parameters
    ----------
    entry:
        A raw YAML row dict (as loaded by ``yaml.safe_load`` on
        ``data/catalogue.yaml``).

    Returns
    -------
    dict[str, bool]
        Mapping of anchor-category names to applicability flags.
    """
    cls = str(entry.get("cycler_class") or "single-ellipse")
    return {
        "vinf": True,
        "a_e": cls == "single-ellipse",
        "period": cls in ("single-ellipse", "multi-arc"),
        "invariants": cls == "multi-arc",
        "cr3bp": cls == "non-keplerian",
    }


__all__ = ["anchors_for", "validate_schema_invariants"]
