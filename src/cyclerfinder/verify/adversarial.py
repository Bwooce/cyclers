"""Phase 5 — the adversarial panel (the Forge novelty gate's teeth).

For each candidate finding, run **N independent verifiers**; a *majority refute*
kills the candidate (forces a REJECTED outcome downstream). Each verifier is a
genuinely independent re-check of a different failure mode:

1. **Perturbed-seed re-solve** — re-run the N-arc corrector from a *perturbed*
   ``(t0, ToF)`` seed and require it to re-close to a chain with a comparable
   peak V_inf (family robustness). A spurious single-seed root that vanishes
   under perturbation is refuted.
2. **Falsification probe** — assert the candidate's own claimed self-consistency
   actually holds: converged (residual below the corrector tolerance),
   bend-feasible, and within the V_inf cap. A fabricated candidate (perturbed
   V_inf, impossible bend, non-converged) fails this immediately.
3. **Re-closure feasibility** — an independent fresh corrector run from the
   *reported* seed must itself converge and be bend-feasible. A candidate whose
   reported geometry does not actually close under an independent solve is
   refuted.

Majority rule: ``n_refuted >= ceil(n_verifiers / 2)`` -> ``majority_refute``.

GOLDEN DISCIPLINE: every probe is a physical self-consistency / robustness check
(closure residual, bend feasibility, V_inf cap, seed-perturbation persistence).
No probe asserts a sourced or self-computed *target* value — they test that the
candidate is what it claims to be, not that it equals a particular number.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any

import numpy as np

from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.search.correct import BallisticClosureResult, ballistic_correct

DAY_S = 86400.0


@dataclass(frozen=True)
class PanelResult:
    """Outcome of an adversarial panel over one candidate (frozen, audit-ready).

    Attributes
    ----------
    n_verifiers:
        How many independent verifiers ran.
    n_refuted:
        How many refuted the candidate.
    majority_refute:
        ``n_refuted >= ceil(n_verifiers / 2)`` — the kill condition.
    verifier_verdicts:
        Per-verifier ``"confirmed"`` / ``"refuted"`` labels (in run order).
    notes:
        Free-form diagnostic.
    """

    n_verifiers: int
    n_refuted: int
    majority_refute: bool
    verifier_verdicts: tuple[str, ...]
    notes: str = ""

    def as_dict(self) -> dict[str, Any]:
        """Serialise for the review-queue / ledger audit trail."""
        return {
            "n_verifiers": self.n_verifiers,
            "n_refuted": self.n_refuted,
            "majority_refute": self.majority_refute,
            "verifier_verdicts": list(self.verifier_verdicts),
            "notes": self.notes,
        }


def _falsification_probe(
    closure: BallisticClosureResult,
    *,
    vinf_cap: float,
) -> bool:
    """Probe 2: the candidate's own claimed self-consistency must hold.

    Returns ``True`` (confirmed) iff converged, bend-feasible, and within cap.
    A fabricated / non-converged candidate fails here.
    """
    if not closure.converged:
        return False
    if not closure.bend_feasible:
        return False
    vinf = closure.vinf_per_encounter_kms
    if not vinf:
        return False
    return max(vinf) <= vinf_cap + 1e-9


def _reclose_probe(
    closure: BallisticClosureResult,
    *,
    sequence: tuple[str, ...],
    per_leg_revs: tuple[int, ...],
    per_leg_branch: tuple[str, ...],
    period_sec: float,
    slack_leg: int,
    vinf_cap: float,
    ephem: Ephemeris,
) -> bool:
    """Probe 3: a fresh corrector run from the reported seed must re-close
    feasibly. Independent of whatever produced ``closure``."""
    free_tofs = [t for i, t in enumerate(closure.tof_days) if i != slack_leg]
    redo = ballistic_correct(
        sequence=sequence,
        per_leg_revs=per_leg_revs,
        per_leg_branch=per_leg_branch,
        t0_seed_sec=closure.t0_sec,
        tof_seed_days=free_tofs,
        period_sec=period_sec,
        ephem=ephem,
        vinf_cap=vinf_cap,
        slack_leg=slack_leg,
    )
    return redo.converged and redo.bend_feasible and redo.vinf_cap_ok


def _perturbed_seed_probe(
    closure: BallisticClosureResult,
    *,
    sequence: tuple[str, ...],
    per_leg_revs: tuple[int, ...],
    per_leg_branch: tuple[str, ...],
    period_sec: float,
    slack_leg: int,
    vinf_cap: float,
    ephem: Ephemeris,
    rng: np.random.Generator,
    epoch_jitter_days: float = 25.0,
    tof_jitter_days: float = 8.0,
    vinf_tol_kms: float = 1.5,
) -> bool:
    """Probe 1: a perturbed-seed re-solve must re-find the family.

    Jitter the launch epoch + free ToFs, re-run the corrector, and require a
    converged, bend-feasible re-closure whose peak V_inf is within
    ``vinf_tol_kms`` of the candidate's (same family). A brittle single-seed root
    that vanishes under perturbation is refuted.
    """
    if not closure.vinf_per_encounter_kms:
        return False
    target_peak = max(closure.vinf_per_encounter_kms)
    free_tofs = [t for i, t in enumerate(closure.tof_days) if i != slack_leg]
    jittered_t0 = closure.t0_sec + float(rng.uniform(-epoch_jitter_days, epoch_jitter_days)) * DAY_S
    jittered_tofs = [
        max(1.0, t + float(rng.uniform(-tof_jitter_days, tof_jitter_days))) for t in free_tofs
    ]
    redo = ballistic_correct(
        sequence=sequence,
        per_leg_revs=per_leg_revs,
        per_leg_branch=per_leg_branch,
        t0_seed_sec=jittered_t0,
        tof_seed_days=jittered_tofs,
        period_sec=period_sec,
        ephem=ephem,
        vinf_cap=vinf_cap,
        slack_leg=slack_leg,
    )
    if not (redo.converged and redo.bend_feasible):
        return False
    redo_peak = max(redo.vinf_per_encounter_kms) if redo.vinf_per_encounter_kms else float("inf")
    return abs(redo_peak - target_peak) <= vinf_tol_kms


def adversarial_panel(
    closure: BallisticClosureResult,
    *,
    sequence: tuple[str, ...],
    per_leg_revs: tuple[int, ...],
    per_leg_branch: tuple[str, ...],
    period_sec: float,
    slack_leg: int,
    vinf_cap: float,
    ephem: Ephemeris,
    n_verifiers: int = 3,
    seed: int = 0,
) -> PanelResult:
    """Run an N-verifier adversarial panel over one candidate closure.

    The three probe *kinds* (falsification, re-closure, perturbed-seed) are each
    independent re-checks. With ``n_verifiers == 3`` each kind runs once; with
    more verifiers the perturbed-seed probe is repeated with fresh jitter (the
    only stochastic kind), broadening the robustness sweep. A candidate is killed
    iff a *majority* of verifiers refute it.

    Returns
    -------
    PanelResult
        With ``majority_refute`` set per the majority rule.
    """
    if n_verifiers < 1:
        raise ValueError("n_verifiers must be >= 1")
    rng = np.random.default_rng(seed)
    verdicts: list[str] = []

    def _record(confirmed: bool) -> None:
        verdicts.append("confirmed" if confirmed else "refuted")

    # Verifier 1 — falsification probe (claimed self-consistency).
    _record(_falsification_probe(closure, vinf_cap=vinf_cap))

    # Verifier 2 — independent re-closure from the reported seed (only meaningful
    # for a candidate that claims to converge; a fabricated non-converged
    # candidate is refuted here too without burning a full solve when the
    # falsification probe already failed badly).
    if n_verifiers >= 2:
        if closure.converged:
            _record(
                _reclose_probe(
                    closure,
                    sequence=sequence,
                    per_leg_revs=per_leg_revs,
                    per_leg_branch=per_leg_branch,
                    period_sec=period_sec,
                    slack_leg=slack_leg,
                    vinf_cap=vinf_cap,
                    ephem=ephem,
                )
            )
        else:
            _record(False)

    # Verifiers 3..N — perturbed-seed robustness (repeated with fresh jitter).
    for _ in range(max(0, n_verifiers - 2)):
        if closure.converged:
            _record(
                _perturbed_seed_probe(
                    closure,
                    sequence=sequence,
                    per_leg_revs=per_leg_revs,
                    per_leg_branch=per_leg_branch,
                    period_sec=period_sec,
                    slack_leg=slack_leg,
                    vinf_cap=vinf_cap,
                    ephem=ephem,
                    rng=rng,
                )
            )
        else:
            _record(False)

    n_refuted = sum(1 for v in verdicts if v == "refuted")
    majority = n_refuted >= math.ceil(n_verifiers / 2)
    return PanelResult(
        n_verifiers=n_verifiers,
        n_refuted=n_refuted,
        majority_refute=majority,
        verifier_verdicts=tuple(verdicts),
    )


__all__ = [
    "PanelResult",
    "adversarial_panel",
]
