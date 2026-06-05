"""Source + fidelity registry and validation-tier classifier.

Phase 0 (data-validation-hardening, Tasks 1 + 4). This module turns the
free-text provenance carried in catalogue ``note:`` prose into a small,
machine-checkable vocabulary so a *test* — not just a human reader — can
distinguish a real cross-source check from a circular same-source one, and
so the S1L1-class **model-fidelity mismatch** (comparing a coplanar V∞
against a real-ephemeris V∞) is representable and refusable.

Two pieces:

* :data:`SOURCE_REGISTRY` — short stable keys → full citations, plus the
  two pseudo-sources ``derived`` and ``computed``.
* :data:`Fidelity` — the three model-fidelity tiers a quantity can be
  stated at, ordered coplanar → analytic-ephemeris → real-DE440.
* :func:`classify_validation` — pure function mapping
  ``(orbit_source, vinf_source, same_fidelity)`` to a :class:`Tier`. It is
  source-independent logic (no catalogue read), so it is golden-clean and
  unit-testable over synthetic inputs.

The catalogue back-fill that *attaches* these keys to YAML fields is a
separate, additive task (plan Task 3); this module only defines and
classifies.
"""

from __future__ import annotations

import enum
from typing import Literal

# ---------------------------------------------------------------------------
# Source registry (Task 1)
# ---------------------------------------------------------------------------

SOURCE_REGISTRY: dict[str, str] = {
    "rogers-2012-t1": (
        "Rogers, B.A., Hughes, K.M., Longuski, J.M. (2012). 'Establishing "
        "Cycler Trajectories Between Earth and Mars.' Acta Astronautica — "
        "Table 1 (ballistic Earth-Mars cycler orbital elements)."
    ),
    "russell-2004-t34": (
        "Russell, R.P., Ocampo, C.A. (2004). 'Systematic method for "
        "constructing Earth-Mars cyclers...' JGCD — Table 3/4 region "
        "(cycler V_inf and geometry tabulation)."
    ),
    "russell-2004-t39_311": (
        "Russell, R.P., Ocampo, C.A. (2004), Table 3.9 / family 3.11 "
        "(Earth-Mars cycler V_inf at encounters)."
    ),
    "russell-2004-t49_413": (
        "Russell, R.P., Ocampo, C.A. (2004), Table 4.9 / family 4.13 "
        "(Earth-Mars cycler V_inf at encounters)."
    ),
    "mcconaghy-2002": (
        "McConaghy, T.T., Longuski, J.M., et al. (2002). Earth-Mars cycler "
        "trajectory analysis (S_nL_m taxonomy and V_inf tabulation)."
    ),
    "mcconaghy-2006": (
        "McConaghy, T.T., et al. (2006). Follow-on Earth-Mars cycler "
        "analysis (orbital elements / encounter geometry)."
    ),
    "spec-9": (
        "Project specification §9 — anchor synodic periods and reference "
        "V_inf values used as in-repo sourced anchors."
    ),
    "hollister-1970-t3": (
        "Hollister, W.M. (1970). 'Castles in space.' Astronautica Acta — "
        "Table 3 (Earth-Venus / Earth-Mars cycler families)."
    ),
    "friedlander-1986": (
        "Friedlander, A.L., et al. (1986). Circulating (cycler) transportation "
        "systems for Mars (orbital elements and V_inf)."
    ),
    # Pseudo-sources: not a citation but an in-row provenance class.
    "derived": (
        "Derived in-repo from other sourced fields on the same row "
        "(e.g. a = (perihelion + aphelion) / 2). NOT an independent source — "
        "must not be paired with its own inputs to claim cross-validation."
    ),
    "computed": (
        "Computed by this repository's own code (optimiser / constructor). "
        "A pseudo-source: never an EXPECTED value in a golden check, and "
        "never counts as an independent corroborating source."
    ),
}
"""Short stable provenance keys → full citation strings.

``derived`` and ``computed`` are pseudo-sources, not citations: ``derived``
marks a value computed from other fields on the same row; ``computed`` marks
a value produced by this repo's own code. Neither may stand in for an
independent source when classifying a cross-validated quantity (golden
discipline)."""

Fidelity = Literal["circular-coplanar", "analytic-ephemeris", "real-de440"]
"""Model-fidelity tier a quantity is stated at, ascending in realism.

* ``circular-coplanar`` — planets on circular, coplanar orbits at mean sma
  (the :func:`~cyclerfinder.search.resonant_construct.construct_resonant_cycler`
  model).
* ``analytic-ephemeris`` — analytic/mean-element ephemeris (eccentric,
  inclined orbits; e.g. the Jones AAS 17-577 multisets).
* ``real-de440`` — full numerical ephemeris (JPL DE440) closure.

Cross-fidelity comparison is the S1L1 5.65-vs-4.99 bug class: a value at one
tier must never be compared against, or corroborated by, a value at another
tier without an explicit fidelity flag."""

_FIDELITIES: frozenset[str] = frozenset({"circular-coplanar", "analytic-ephemeris", "real-de440"})
"""Runtime set of the :data:`Fidelity` literal values (``Literal`` is not
introspectable at runtime, so the membership set is kept explicitly)."""

_PSEUDO_SOURCES: frozenset[str] = frozenset({"derived", "computed"})
"""Provenance keys that are NOT independent citations."""


def is_fidelity(value: str) -> bool:
    """Return True iff *value* is one of the three :data:`Fidelity` tiers."""
    return value in _FIDELITIES


def is_registry_key(value: str) -> bool:
    """Return True iff *value* is a known :data:`SOURCE_REGISTRY` key."""
    return value in SOURCE_REGISTRY


def is_independent_source(value: str) -> bool:
    """Return True iff *value* is a real citation (a registry key that is
    not a pseudo-source).

    ``derived`` / ``computed`` return ``False`` — they can never count as an
    independent corroborating source (golden discipline).
    """
    return value in SOURCE_REGISTRY and value not in _PSEUDO_SOURCES


# ---------------------------------------------------------------------------
# Validation-tier classifier (Task 4)
# ---------------------------------------------------------------------------


class Tier(enum.Enum):
    """How strongly a reproduction is validated, per plan Task 4.

    A quantity (the row's orbit + V∞) is:

    * :attr:`CROSS_VALIDATED` — its orbit and its V∞ trace to two
      *different independent* sources at the *same* fidelity. The strongest
      tier: an independent second source corroborates, and no cross-fidelity
      mismatch can hide in the comparison.
    * :attr:`CONSISTENCY_CHECKED` — both are set but trace to the *same*
      source. Catches transcription/round-trip error only, not a source
      mistake; cannot be promoted to cross-validated without a 2nd source.
    * :attr:`UNVALIDATED` — anything else (a source missing, a pseudo-source
      standing in for a real one, or a cross-fidelity comparison).
    """

    CROSS_VALIDATED = "cross_validated"
    CONSISTENCY_CHECKED = "consistency_checked"
    UNVALIDATED = "unvalidated"


def classify_validation(
    orbit_source: str | None,
    vinf_source: str | None,
    *,
    same_fidelity: bool,
) -> Tier:
    """Classify a reproduction's validation strength (pure function).

    Parameters
    ----------
    orbit_source:
        Provenance key for the row's orbit elements, or ``None`` if untagged.
    vinf_source:
        Provenance key for the row's V∞ multiset, or ``None`` if untagged.
    same_fidelity:
        Whether the orbit and V∞ are stated at the *same* model fidelity.
        A cross-fidelity comparison can never be cross-validated (the S1L1
        bug class), so ``False`` here caps the result at
        :attr:`Tier.UNVALIDATED`.

    Returns
    -------
    Tier
        * :attr:`Tier.CROSS_VALIDATED` iff both sources are set, are
          *independent* (real citations, not pseudo-sources), are *different*,
          and ``same_fidelity`` is ``True``.
        * :attr:`Tier.CONSISTENCY_CHECKED` iff both are set, are real
          citations, are the *same* source, and ``same_fidelity`` is ``True``.
        * :attr:`Tier.UNVALIDATED` otherwise.

    Notes
    -----
    Source-independent by construction (no catalogue read), so this is
    golden-clean: it *enforces* the independence rule rather than asserting a
    computed value. A pseudo-source (``derived``/``computed``) on either side
    can never reach a validated tier — it is not an independent citation.
    """
    if orbit_source is None or vinf_source is None:
        return Tier.UNVALIDATED
    # A pseudo-source can never validate: it is not an independent citation.
    if not is_independent_source(orbit_source) or not is_independent_source(vinf_source):
        return Tier.UNVALIDATED
    # Cross-fidelity comparison is never trustworthy (S1L1 mismatch class).
    if not same_fidelity:
        return Tier.UNVALIDATED
    if orbit_source != vinf_source:
        return Tier.CROSS_VALIDATED
    return Tier.CONSISTENCY_CHECKED


__all__ = [
    "SOURCE_REGISTRY",
    "Fidelity",
    "Tier",
    "classify_validation",
    "is_fidelity",
    "is_independent_source",
    "is_registry_key",
]
