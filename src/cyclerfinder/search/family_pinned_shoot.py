"""Family-pinned penalty homotopy closer (#388).

Drives the full STM multiple-shooting corrector (``nbody.shooter.shoot``) toward
the PUBLISHED V∞ family by ramping a V∞-anchor penalty weight from a calibrated
``W`` down to zero. Each rung warm-starts from the previous rung's corrected
states (the ``continuation.continuation_correct`` ladder pattern). The final
``weight == 0`` rung is the verdict solve: its emerged V∞ is the recorded value
(the penalty is a basin-selector ramped to zero, never the recorded number — the
golden-discipline requirement of ``feedback_golden_tests_sourced_only``).

See ``docs/superpowers/specs/2026-06-22-family-pinned-homotopy-closer-design.md``.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass
from typing import TYPE_CHECKING

from cyclerfinder.nbody.shooter import ShootingSeed, ShootResult, _seed_with_states, shoot

if TYPE_CHECKING:
    from cyclerfinder.core.ephemeris import Ephemeris


@dataclass(frozen=True)
class FamilyPinnedResult:
    """Result of a family-pinned penalty homotopy run.

    ``final`` is the unpenalized (``weight == 0``) ShootResult — the verdict solve
    whose emerged V∞ is the recorded value. ``trace`` is one
    ``(weight, defect_norm, vinf_per_encounter_kms)`` tuple per ladder rung.
    ``anchor_retention_kms`` is the change in the best per-anchor V∞ residual from
    the FIRST penalized rung to the final ``weight == 0`` rung — small means V∞
    held near the published family as the penalty lifted; large means it snapped
    off-anchor (the stronger characterized negative).
    """

    final: ShootResult
    final_weight: float
    trace: list[tuple[float, float, list[float]]]
    anchor_retention_kms: float
    vinf_anchors: dict[str, float]


def _best_anchor_residual(
    vinf_per_encounter_kms: Sequence[float], vinf_anchors: Mapping[str, float]
) -> float:
    """Worst-over-anchors of the best-matching encounter V∞ (km/s)."""
    worst = 0.0
    for anchor in vinf_anchors.values():
        best = min((abs(v - float(anchor)) for v in vinf_per_encounter_kms), default=float("inf"))
        worst = max(worst, best)
    return worst


def family_pinned_shoot(
    seed: ShootingSeed,
    *,
    ephem: Ephemeris,
    bodies: Sequence[str],
    vinf_anchors: Mapping[str, float],
    weight_ladder: Sequence[float] = (40.0, 10.0, 2.5, 0.5, 0.0),
    accuracy: float = 1e-9,
    max_nfev: int = 100,
    max_wall_sec: float = 30.0,
    progress: Callable[[str, int, float, float], None] | None = None,
) -> FamilyPinnedResult:
    """Ramp the V∞-anchor penalty from ``weight_ladder[0]`` down to (a final) 0.

    Each rung runs ``shoot(jacobian="stm", vinf_anchors=..., vinf_weight=w)``,
    warm-started from the previous rung's corrected states. ``weight_ladder`` MUST
    end at 0.0 (asserted) so the verdict solve is unpenalized. Returns a
    :class:`FamilyPinnedResult`.
    """
    ladder = [float(w) for w in weight_ladder]
    if not ladder or ladder[-1] != 0.0:
        raise ValueError("weight_ladder must be non-empty and end at 0.0 (the verdict rung)")
    anchors = {str(k): float(v) for k, v in vinf_anchors.items()}

    cur = seed
    trace: list[tuple[float, float, list[float]]] = []
    last: ShootResult | None = None
    first_penalized_resid: float | None = None

    for w in ladder:
        res = shoot(
            cur,
            ephem=ephem,
            bodies=bodies,
            accuracy=accuracy,
            max_nfev=max_nfev,
            max_wall_sec=max_wall_sec,
            jacobian="stm",
            vinf_anchors=anchors,
            vinf_weight=w,
            progress=progress,
        )
        vinf = list(res.vinf_per_encounter_kms)
        trace.append((w, res.defect_norm, vinf))
        if w > 0.0 and first_penalized_resid is None:
            first_penalized_resid = _best_anchor_residual(vinf, anchors)
        # warm-start the next rung from this rung's corrected states
        cur = _seed_with_states(cur, res.corrected_states)
        last = res

    assert last is not None
    final_resid = _best_anchor_residual(list(last.vinf_per_encounter_kms), anchors)
    retention = (
        abs(final_resid - first_penalized_resid)
        if first_penalized_resid is not None
        else float("nan")
    )
    return FamilyPinnedResult(
        final=last,
        final_weight=ladder[-1],
        trace=trace,
        anchor_retention_kms=retention,
        vinf_anchors=anchors,
    )


__all__ = ["FamilyPinnedResult", "family_pinned_shoot"]
