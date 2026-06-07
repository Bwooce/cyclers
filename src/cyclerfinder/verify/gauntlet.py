"""Phase 3 — the gauntlet combiner: Axes A-D → a tiered :class:`ValidationVerdict`.

The four cross-validation axes each already exist as code and emit their own
frozen report:

* **Axis A** — code-path agreement
  (:func:`cyclerfinder.verify.agreement.crosscheck_code_paths` → ``AgreementReport``):
  ``>= 2`` independent in-house code paths must agree on V∞/geometry.
* **Axis B** — fidelity persistence
  (:func:`cyclerfinder.verify.fidelity.fidelity_persistence` → ``PersistenceReport``):
  a tracked quantity must persist across the fidelity ladder, or shift in a
  *documented* direction.
* **Axis C** — provenance / corroboration
  (:func:`cyclerfinder.data.provenance.classify_validation` /
  :func:`~cyclerfinder.data.provenance.score_corroboration`): is the row backed
  by an *independent* second source?
* **Axis D** — falsification / independent-tool guards: a deliberately bogus
  candidate must be refuted; an adversarial panel that refutes kills the verdict.

:func:`run_gauntlet` is a **pure combiner**. It does not recompute any physics
and never invents a value — it folds the axis reports (which the caller computes
with the existing machinery) into one of four tiers per the Forge spec's
governing principle: *independence enforced at every tier*.

Combination rules (Forge spec §16.7 governing principle)
--------------------------------------------------------
Let A = Axis A ``agreed``; B = no ``SHIFTS_UNDOCUMENTED`` among the supplied
persistence reports; C_source = an *independent* corroborating source exists
(Axis C is ``CROSS_VALIDATED`` or ``STRONGLY_SOURCED``); D_ok = not falsified.

* **REJECTED** — ``not D_ok`` (falsified / adversarially refuted), OR an axis
  that *ran and failed* in a disqualifying way: Axis A available-but-not-agreed
  (a failing available path vetoes — see ``AgreementReport.agreed``), or a
  ``SHIFTS_UNDOCUMENTED`` persistence report, or a ``DISPUTED`` corroboration
  whose dispute is a cross-fidelity mismatch (the S1L1 bug class). Falsification
  dominates: a refuted candidate is REJECTED even if other axes look clean.
* **GOLD** — ``A and B and C_source and D_ok``: machine-confirmed by ``>= 2``
  agreeing paths, fidelity-robust, AND independently sourced. The only tier that
  may be called a credible reproduction without further human work.
* **SILVER** — ``A and B and D_ok`` but ``not C_source``: machine-confirmed but
  *unsourced* (or single-sourced). Explicitly **capped pending human review** —
  this is where a novel discovery lands; it is never auto-promoted to GOLD.
* **BRONZE** — none of the failing conditions, but the machine-confirmation
  predicate (A) is not met because Axis A was *unavailable* (fewer than two
  paths ran) rather than failing. A weak, non-refuted signal: some evidence,
  not enough to assert agreement. (If A were available-but-failing the verdict
  would be REJECTED, above.)

Supersession-aware matching (Forge R1 delta 3)
----------------------------------------------
A verdict on a candidate that matched a catalogue row carrying ``superseded_by``
must carry the supersession chain in :attr:`ValidationVerdict.provenance`
(``superseded_by``) — it must never report a clean "known" match against an
invalidated premise. :func:`run_gauntlet` records the chain when the caller
supplies it via ``superseded_by``; it does not downgrade the tier on that basis
(the superseding rows may themselves be GOLD), but it makes the invalidated
premise visible in every verdict's provenance.

Golden discipline
-----------------
This module asserts no computed physics value. Each axis input is a
*classification* the caller computed elsewhere (agreement predicate, persistence
class, provenance tier, falsification flag). The combiner only applies the
independence rules above. ``confidence`` is an ordinal label derived from the
tier, not a fabricated probability.
"""

from __future__ import annotations

import enum
from dataclasses import dataclass, field
from typing import Any

from cyclerfinder.data.provenance import Corroboration, Tier
from cyclerfinder.verify.agreement import AgreementReport
from cyclerfinder.verify.fidelity import PersistenceClass, PersistenceReport


class VerdictTier(enum.Enum):
    """The four gauntlet tiers, strongest → weakest, with REJECTED separate.

    * :attr:`GOLD` — machine-confirmed AND independently sourced; credible.
    * :attr:`SILVER` — machine-confirmed but unsourced; capped pending human
      review (the novel-discovery holding tier).
    * :attr:`BRONZE` — partial, non-refuted evidence (Axis A unavailable, not
      failing); a weak signal.
    * :attr:`REJECTED` — falsified, adversarially refuted, or an axis ran and
      failed disqualifyingly.
    """

    GOLD = "gold"
    SILVER = "silver"
    BRONZE = "bronze"
    REJECTED = "rejected"


# Ordinal confidence label per tier (NOT a fabricated probability; an ordinal).
_TIER_CONFIDENCE: dict[VerdictTier, str] = {
    VerdictTier.GOLD: "high",
    VerdictTier.SILVER: "medium",
    VerdictTier.BRONZE: "low",
    VerdictTier.REJECTED: "none",
}


@dataclass(frozen=True)
class ValidationVerdict:
    """Frozen gauntlet verdict — the Phase 3 combiner output.

    Attributes
    ----------
    tier:
        The :class:`VerdictTier`.
    confidence:
        Ordinal label (``"high"`` / ``"medium"`` / ``"low"`` / ``"none"``)
        derived from the tier. Not a probability.
    axis_results:
        Per-axis machine-readable summary of what each axis contributed
        (``"A"`` / ``"B"`` / ``"C"`` / ``"D"`` keys). Carries the booleans the
        combiner branched on plus the raw class labels, so a verdict is fully
        auditable from this field alone.
    provenance:
        Audit trail: ``candidate_id``, the matched ``known_id`` (if any), the
        ``superseded_by`` chain (Forge R1 delta 3), and any free-form
        ``notes``. Never a clean "known" when ``superseded_by`` is non-empty.
    """

    tier: VerdictTier
    confidence: str
    axis_results: dict[str, Any] = field(default_factory=dict)
    provenance: dict[str, Any] = field(default_factory=dict)


def _axis_a_summary(report: AgreementReport | None) -> tuple[dict[str, Any], bool, bool]:
    """Return (summary, available, failing) for Axis A.

    * ``available`` — at least two paths ran (the machine-confirmation
      predicate can be evaluated).
    * ``failing`` — an available Axis-A path *ran and failed its tolerance*
      (``n_paths_passed < n_paths_available``). This is the disqualifying case
      (REJECTED). A single passing path that simply does not reach the ``>= 2``
      machine-confirmation predicate is NOT failing — it is thin coverage
      (``available=False`` with ``failing=False`` ⇒ BRONZE). The veto is a path
      that contradicts the others, never merely a shortage of paths.
    """
    if report is None:
        return ({"available": False, "agreed": None}, False, False)
    available = report.n_paths_available >= 2
    # "Failing" means an available path ran and did NOT pass its tolerance — a
    # contradicting witness, which vetoes. (Distinct from too-few-paths.)
    failing = report.n_paths_passed < report.n_paths_available
    summary = {
        "available": available,
        "agreed": report.agreed,
        "n_paths_available": report.n_paths_available,
        "n_paths_passed": report.n_paths_passed,
    }
    return (summary, available, failing)


def _axis_b_summary(
    reports: tuple[PersistenceReport, ...],
) -> tuple[dict[str, Any], bool, bool]:
    """Return (summary, ran, undocumented_shift) for Axis B.

    * ``ran`` — at least one persistence report was supplied.
    * ``undocumented_shift`` — any supplied report is ``SHIFTS_UNDOCUMENTED``
      (the red-flag cross-fidelity confusion class → REJECTED).
    """
    ran = len(reports) > 0
    undocumented = any(r.classification is PersistenceClass.SHIFTS_UNDOCUMENTED for r in reports)
    summary = {
        "ran": ran,
        "classifications": [r.classification.value for r in reports],
        "undocumented_shift": undocumented,
    }
    return (summary, ran, undocumented)


def _axis_c_summary(
    tier: Tier | None,
    corroboration: Corroboration | None,
    *,
    corroboration_cross_fidelity: bool,
) -> tuple[dict[str, Any], bool, bool]:
    """Return (summary, has_independent_source, disputed_cross_fidelity).

    * ``has_independent_source`` — Axis C says an independent second source
      corroborates: ``Tier.CROSS_VALIDATED`` or
      ``Corroboration.STRONGLY_SOURCED``.
    * ``disputed_cross_fidelity`` — ``Corroboration.DISPUTED`` whose spread is
      (at least partly) a cross-fidelity mismatch (the S1L1 bug class →
      REJECTED). A pure single-fidelity numeric dispute does not reject (it is
      surfaced, not auto-killed — a human resolves it), but it does not count
      as an independent source either.
    """
    has_source = tier is Tier.CROSS_VALIDATED or corroboration is Corroboration.STRONGLY_SOURCED
    disputed_cross_fid = corroboration is Corroboration.DISPUTED and corroboration_cross_fidelity
    summary = {
        "tier": tier.value if tier is not None else None,
        "corroboration": corroboration.value if corroboration is not None else None,
        "has_independent_source": has_source,
        "disputed_cross_fidelity": disputed_cross_fid,
    }
    return (summary, has_source, disputed_cross_fid)


def run_gauntlet(
    candidate_id: str,
    *,
    agreement: AgreementReport | None = None,
    persistence_reports: tuple[PersistenceReport, ...] = (),
    provenance_tier: Tier | None = None,
    corroboration: Corroboration | None = None,
    corroboration_cross_fidelity: bool = False,
    falsified: bool = False,
    known_id: str | None = None,
    superseded_by: tuple[str, ...] = (),
    notes: str | None = None,
) -> ValidationVerdict:
    """Combine Axes A-D into a tiered :class:`ValidationVerdict` (pure combiner).

    Parameters
    ----------
    candidate_id:
        Identifier of the candidate / cell under verdict (carried into
        provenance for the audit trail).
    agreement:
        Axis A :class:`~cyclerfinder.verify.agreement.AgreementReport`, or
        ``None`` if Axis A was not run.
    persistence_reports:
        Axis B :class:`~cyclerfinder.verify.fidelity.PersistenceReport` tuple
        (one per tracked quantity); empty if Axis B was not run.
    provenance_tier:
        Axis C :class:`~cyclerfinder.data.provenance.Tier`
        (``classify_validation`` output), or ``None``.
    corroboration:
        Axis C :class:`~cyclerfinder.data.provenance.Corroboration`
        (``score_corroboration`` output), or ``None``.
    corroboration_cross_fidelity:
        Whether the Axis C corroboration dispute spans more than one fidelity
        (``CorroborationScore.cross_fidelity``). Only meaningful when
        ``corroboration is DISPUTED``.
    falsified:
        Axis D: ``True`` iff a falsification / independent-tool guard refuted
        the candidate (an adversarial-panel majority-refute). Dominates: forces
        REJECTED.
    known_id:
        The catalogue row this candidate matched, if any (provenance).
    superseded_by:
        The ``superseded_by`` chain of the matched row (Forge R1 delta 3). When
        non-empty, the verdict's provenance carries the chain and the match is
        never reported as a clean "known".
    notes:
        Free-form audit note (provenance).

    Returns
    -------
    ValidationVerdict
        The combined frozen verdict.
    """
    a_summary, a_available, a_failing = _axis_a_summary(agreement)
    b_summary, _b_ran, b_undocumented = _axis_b_summary(persistence_reports)
    c_summary, c_has_source, c_disputed_cross_fid = _axis_c_summary(
        provenance_tier,
        corroboration,
        corroboration_cross_fidelity=corroboration_cross_fidelity,
    )

    d_summary = {"falsified": falsified}

    axis_results: dict[str, Any] = {
        "A": a_summary,
        "B": b_summary,
        "C": c_summary,
        "D": d_summary,
    }

    machine_confirmed = bool(agreement is not None and agreement.agreed and a_available)

    # --- Tier decision (REJECTED dominates) ---------------------------------
    if falsified or a_failing or b_undocumented or c_disputed_cross_fid:
        tier = VerdictTier.REJECTED
    elif machine_confirmed and c_has_source:
        tier = VerdictTier.GOLD
    elif machine_confirmed:
        tier = VerdictTier.SILVER
    else:
        tier = VerdictTier.BRONZE

    provenance: dict[str, Any] = {
        "candidate_id": candidate_id,
        "known_id": known_id,
        "superseded_by": list(superseded_by),
    }
    if superseded_by:
        # Never report a clean "known" against an invalidated premise.
        provenance["match_status"] = "superseded"
    elif known_id is not None:
        provenance["match_status"] = "known"
    else:
        provenance["match_status"] = "unmatched"
    if notes is not None:
        provenance["notes"] = notes

    return ValidationVerdict(
        tier=tier,
        confidence=_TIER_CONFIDENCE[tier],
        axis_results=axis_results,
        provenance=provenance,
    )


def validate_verdict(verdict: ValidationVerdict) -> None:
    """Re-derive the tier from a verdict's own ``axis_results`` and refuse it if
    they do not justify the claimed tier (the fabrication / drift guard).

    :func:`run_gauntlet` produces self-consistent verdicts; this validator is
    the *teeth*: it independently checks that the recorded per-axis booleans
    actually entail the stated tier and confidence. A hand-built or tampered
    verdict that over-claims (e.g. ``GOLD`` without an independent source, or a
    ``REJECTED`` with no refuting cause) is rejected. It enforces the same
    governing rules :func:`run_gauntlet` applies, never invents physics, and is
    safe to run over a persisted ledger record reconstructed into a verdict.

    Raises
    ------
    ValueError
        If the verdict is internally inconsistent (the tier is not the one the
        ``axis_results`` imply, the confidence label does not match the tier,
        or the supersession provenance is malformed).
    """
    a = verdict.axis_results.get("A", {})
    b = verdict.axis_results.get("B", {})
    c = verdict.axis_results.get("C", {})
    d = verdict.axis_results.get("D", {})

    # Re-derive the predicates the combiner branched on.
    a_available = bool(a.get("available"))
    a_agreed = a.get("agreed")
    # available-but-failing: a path ran but did NOT pass its tolerance (a
    # contradicting witness). Mirrors _axis_a_summary's `failing` rule.
    # Review finding I6 (2026-06-07): guessing missing path counts (the old
    # `2 if a_available else 0` default) could make this validator disagree
    # with the combiner's own verdict for hand-built axis dicts. A validator
    # whose purpose is catching fabricated verdicts must be strict on
    # ambiguity: an AVAILABLE Axis A must carry explicit counts.
    if a_available and ("n_paths_available" not in a or "n_paths_passed" not in a):
        raise ValueError(
            "invalid verdict: axis_results['A'] is available but missing "
            "explicit n_paths_available/n_paths_passed (refusing to guess; "
            "see _axis_a_summary for the producing contract)"
        )
    n_paths_available = int(a.get("n_paths_available", 0))
    n_paths_passed = int(a.get("n_paths_passed", 0))
    a_failing = n_paths_passed < n_paths_available
    b_undocumented = bool(b.get("undocumented_shift"))
    c_has_source = bool(c.get("has_independent_source"))
    c_disputed_cross_fid = bool(c.get("disputed_cross_fidelity"))
    falsified = bool(d.get("falsified"))

    machine_confirmed = a_available and a_agreed is True

    if falsified or a_failing or b_undocumented or c_disputed_cross_fid:
        expected = VerdictTier.REJECTED
    elif machine_confirmed and c_has_source:
        expected = VerdictTier.GOLD
    elif machine_confirmed:
        expected = VerdictTier.SILVER
    else:
        expected = VerdictTier.BRONZE

    # Tier-specific over-claim guards, phrased for diagnostics.
    if verdict.tier is VerdictTier.GOLD and not c_has_source:
        raise ValueError(
            "invalid verdict: GOLD requires an independent source "
            "(axis C has_independent_source must be True)"
        )
    if verdict.tier is VerdictTier.REJECTED and expected is not VerdictTier.REJECTED:
        raise ValueError(
            "invalid verdict: REJECTED requires a refuting cause "
            "(falsified, a failing Axis A, an undocumented shift, or a "
            "cross-fidelity dispute) in axis_results"
        )
    if verdict.tier is not expected:
        raise ValueError(
            f"invalid verdict: tier {verdict.tier.value!r} does not match the "
            f"tier {expected.value!r} implied by axis_results"
        )
    if verdict.confidence != _TIER_CONFIDENCE[verdict.tier]:
        raise ValueError(
            f"invalid verdict: confidence {verdict.confidence!r} does not match "
            f"tier {verdict.tier.value!r} (expected {_TIER_CONFIDENCE[verdict.tier]!r})"
        )

    # Supersession provenance must be a list; a non-empty chain forbids a clean
    # "known" match status (Forge R1 delta 3).
    chain = verdict.provenance.get("superseded_by", [])
    if not isinstance(chain, list):
        raise ValueError("invalid verdict: provenance.superseded_by must be a list")
    if chain and verdict.provenance.get("match_status") == "known":
        raise ValueError("invalid verdict: a superseded match must not report match_status 'known'")


__all__ = [
    "ValidationVerdict",
    "VerdictTier",
    "run_gauntlet",
    "validate_verdict",
]
