"""Publication-layer plausibility predicates — refuse-with-reason (task #127).

The 55.32 km/s near-miss lesson: a degenerate maintenance-ΔV solve slid onto an
off-family basin (55.3 km/s) and *reached a publication surface* — it was caught
only by a manual cross-check, not by any automated gate. This module is the
automated gate that should have been there: a small, pure set of per-quantity
predicates with **refuse-with-reason** semantics, called mechanically by every
publication surface (the site exporter, the validation-level writeback, the
review-queue append) so an implausible value is recorded-and-refused instead of
silently written.

Two kinds of bar, kept explicitly distinct (the discipline the brief demands):

* **Engineering conventions** — operational bars, not physics. A maintenance ΔV
  above ~3 km/s is "implausibly large for a TCM budget" by the project's own
  ``dv < 3.0`` convention (``tests/search/test_maintain.py:153``); it is an
  *engineering* judgement about what a station-keeping budget should look like,
  not a law of nature. A high-energy *cycler* could legitimately want more — the
  bar guards the maintenance-ΔV *quantity kind*, framed as a budget.

* **Physics ceilings** — hard limits. The elliptic-periodicity V_inf ceiling
  (:data:`~cyclerfinder.core.constants.VINF_CEILING_KMS`) is physics: a periodic
  heliocentric orbit *cannot* have V_inf above ``v_esc_sun(r_B) + v_B``. A value
  above it is impossible, full stop.

The V_inf predicate layers three bars (most-permissive sourced bar first, then
the data-layer invariant, then the physics ceiling), so an over-ceiling value is
attributed to physics and an over-distribution-but-sub-ceiling value to the
sourced catalogue distribution.

Pure: depends only on the constants registry and a small live-catalogue read for
the sourced V_inf distribution maximum (computed once, cached).
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum
from functools import lru_cache
from pathlib import Path
from typing import Any

import yaml  # type: ignore[import-untyped]

from cyclerfinder.core.constants import VINF_CEILING_KMS

# --- Engineering conventions (bars, NOT physics) ---------------------------

MAINTENANCE_DV_CONVENTION_KMS: float = 3.0
"""Maintenance-ΔV engineering ceiling (km/s). CONVENTION, not physics.

A per-synodic TCM / station-keeping budget above this is implausibly large for a
maintenance quantity — the project's own sanity bar
(``tests/search/test_maintain.py:153``: ``dv < 3.0``). It exists to catch the
degenerate-basin artifact (the 55.3 km/s off-family solve), NOT to bound a
high-energy mission. Not a match against any published value.
"""

# --- Data-layer invariant (mirrors data/validate.py) -----------------------

VINF_DATA_INVARIANT_KMS: float = 50.0
"""Catalogue V_inf unit-error invariant (km/s); mirrors ``data/validate.py``.

The data layer refuses any catalogue V_inf >= 50 km/s as a likely unit error
(m/s entered in a km/s field lands at ~10^3). Sits above all real catalogue data
(max ~20.3) yet far below the 1000x unit-error class. We re-assert it here so the
publication layer enforces the SAME bar the data loader does.
"""

# --- Sourced-distribution headroom -----------------------------------------

VINF_DISTRIBUTION_HEADROOM_KMS: float = 5.0
"""Documented headroom (km/s) added to the live-catalogue V_inf maximum.

CONVENTION. The sourced distribution currently tops out at ~20.3 km/s
(Russell-Ocampo). A *published* V_inf meaningfully above the sourced maximum is
suspect — but the catalogue grows, so we add fixed headroom rather than pinning
to the exact current max. A value within (max, max+headroom] passes with a note;
above max+headroom it is refused as out-of-distribution.
"""

DEFAULT_CATALOGUE_PATH: Path = (
    Path(__file__).resolve().parent.parent.parent.parent / "data" / "catalogue.yaml"
)


class QuantityKind(StrEnum):
    """The publication quantities this module knows how to vet."""

    MAINTENANCE_DV_KMS = "maintenance_dv_kms"
    VINF_KMS = "vinf_kms"


@dataclass(frozen=True)
class PlausibilityVerdict:
    """Outcome of a plausibility check: publishable-or-not, with a reason.

    ``ok`` True means the value may be published; ``reason`` then carries an
    informational note (which bars it cleared, or a sourced-headroom caveat).
    ``ok`` False means REFUSE: ``reason`` is the human-readable refusal cause to
    record in place of the value.
    """

    ok: bool
    reason: str


@lru_cache(maxsize=4)
def _sourced_vinf_max_kms(catalogue_path: str) -> float:
    """Largest sourced V_inf in the live catalogue (km/s), cached per path.

    Reads ``vinf_kms_at_encounters[].vinf_kms`` across every row. Returns 0.0 for
    an unreadable / empty catalogue (the headroom bar then degrades to the
    data-invariant bar, never to a free pass — the ceiling still applies).
    """
    p = Path(catalogue_path)
    try:
        rows = yaml.safe_load(p.read_text(encoding="utf-8"))
    except (OSError, yaml.YAMLError):
        return 0.0
    if not isinstance(rows, list):
        return 0.0
    best = 0.0
    for row in rows:
        for enc in row.get("vinf_kms_at_encounters") or []:
            v = enc.get("vinf_kms")
            if isinstance(v, (int, float)) and v > best:
                best = float(v)
    return best


def sourced_vinf_max_kms(catalogue_path: Path | str | None = None) -> float:
    """Public accessor for the live-catalogue sourced V_inf maximum (km/s)."""
    path = str(catalogue_path) if catalogue_path is not None else str(DEFAULT_CATALOGUE_PATH)
    return _sourced_vinf_max_kms(path)


def _check_maintenance_dv(value: float, context: dict[str, Any]) -> PlausibilityVerdict:
    """Maintenance-ΔV predicate: finite, non-negative, below the engineering bar."""
    if value != value or value in (float("inf"), float("-inf")):  # NaN / inf
        return PlausibilityVerdict(False, f"maintenance ΔV is not finite: {value}")
    if value < 0.0:
        return PlausibilityVerdict(False, f"maintenance ΔV is negative: {value} km/s")
    if value > MAINTENANCE_DV_CONVENTION_KMS:
        return PlausibilityVerdict(
            False,
            f"maintenance ΔV {value:.4g} km/s exceeds the engineering bar "
            f"{MAINTENANCE_DV_CONVENTION_KMS} km/s (convention, not physics: a TCM "
            f"budget this large signals a degenerate / off-family solve, the 55.32 "
            f"km/s near-miss class) — refusing to publish.",
        )
    return PlausibilityVerdict(
        True,
        f"maintenance ΔV {value:.4g} km/s within the {MAINTENANCE_DV_CONVENTION_KMS} "
        f"km/s engineering bar.",
    )


def _check_vinf(value: float, context: dict[str, Any]) -> PlausibilityVerdict:
    """V_inf predicate: physics ceiling, data invariant, sourced distribution.

    ``context`` may carry ``body`` (str) to apply that body's elliptic-periodicity
    ceiling; absent a body, the most-permissive ceiling (Mercury) is used so the
    physics bar is never stricter than the unknown body could justify.
    ``catalogue_path`` overrides the sourced-distribution source.
    """
    if value != value or value in (float("inf"), float("-inf")):
        return PlausibilityVerdict(False, f"V_inf is not finite: {value}")
    if value < 0.0:
        return PlausibilityVerdict(False, f"V_inf is negative: {value} km/s")

    body = context.get("body")
    if body is not None and body in VINF_CEILING_KMS:
        ceiling = VINF_CEILING_KMS[body]
        ceiling_label = f"the {body} elliptic-periodicity ceiling"
    else:
        # No body: use the most permissive physics ceiling so we never refuse a
        # value that SOME body could physically support.
        ceiling = max(VINF_CEILING_KMS.values())
        ceiling_label = "the most-permissive (Mercury) elliptic-periodicity ceiling"

    # Physics ceiling first — a breach here is impossible, not merely unusual.
    if value > ceiling:
        return PlausibilityVerdict(
            False,
            f"V_inf {value:.4g} km/s exceeds {ceiling_label} ({ceiling:.4g} km/s): "
            f"PHYSICALLY IMPOSSIBLE for a periodic heliocentric cycler "
            f"(|V_inf| <= v_esc_sun(r_B) + v_B) — refusing to publish.",
        )
    # Data-layer unit-error invariant.
    if value >= VINF_DATA_INVARIANT_KMS:
        return PlausibilityVerdict(
            False,
            f"V_inf {value:.4g} km/s breaches the data-layer invariant "
            f"{VINF_DATA_INVARIANT_KMS} km/s (likely m/s-in-a-km/s-field unit error) "
            f"— refusing to publish.",
        )
    # Sourced-distribution bar (most permissive; only refuses far-outliers).
    sourced_max = sourced_vinf_max_kms(context.get("catalogue_path"))
    distribution_bar = sourced_max + VINF_DISTRIBUTION_HEADROOM_KMS
    if value > distribution_bar:
        return PlausibilityVerdict(
            False,
            f"V_inf {value:.4g} km/s exceeds the sourced-distribution bar "
            f"{distribution_bar:.4g} km/s (live-catalogue max {sourced_max:.4g} + "
            f"{VINF_DISTRIBUTION_HEADROOM_KMS} km/s headroom) — refusing to publish "
            f"a value far outside the sourced distribution.",
        )
    note = (
        f"V_inf {value:.4g} km/s within sourced distribution (<= {sourced_max:.4g} km/s)."
        if value <= sourced_max
        else (
            f"V_inf {value:.4g} km/s above the sourced max {sourced_max:.4g} km/s but "
            f"within the {VINF_DISTRIBUTION_HEADROOM_KMS} km/s documented headroom."
        )
    )
    return PlausibilityVerdict(True, note)


_CHECKERS = {
    QuantityKind.MAINTENANCE_DV_KMS: _check_maintenance_dv,
    QuantityKind.VINF_KMS: _check_vinf,
}


def check_publishable(
    quantity_kind: QuantityKind | str,
    value: float,
    context: dict[str, Any] | None = None,
) -> PlausibilityVerdict:
    """Vet a single quantity for publication; refuse-with-reason.

    Parameters
    ----------
    quantity_kind:
        Which predicate to apply (a :class:`QuantityKind` or its string value).
    value:
        The numeric value about to be published.
    context:
        Optional per-quantity context. For ``vinf_kms``: ``body`` (str) selects
        the body's physics ceiling, ``catalogue_path`` overrides the sourced
        distribution source.

    Returns
    -------
    PlausibilityVerdict
        ``ok=True`` -> publishable (``reason`` is an informational note);
        ``ok=False`` -> REFUSE, record ``reason`` in place of the value.
    """
    try:
        kind = QuantityKind(quantity_kind)
    except ValueError as exc:
        raise ValueError(
            f"unknown quantity_kind {quantity_kind!r}; expected one of "
            f"{[k.value for k in QuantityKind]}"
        ) from exc
    return _CHECKERS[kind](float(value), context or {})


__all__ = [
    "MAINTENANCE_DV_CONVENTION_KMS",
    "VINF_DATA_INVARIANT_KMS",
    "VINF_DISTRIBUTION_HEADROOM_KMS",
    "PlausibilityVerdict",
    "QuantityKind",
    "check_publishable",
    "sourced_vinf_max_kms",
]
