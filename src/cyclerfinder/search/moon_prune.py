"""Incremental VILM/Tisserand/bend prune gate for moon-pair legs (Forge Phase 6).

The Ceriotti 2010 incremental-pruning import (design note
``docs/notes/2026-06-08-forge-phase6-discovery-design.md`` §2; mining note
``docs/notes/2026-06-07-ceriotti-2010-mga-global-opt-mining.md`` §7): promote the
post-hoc per-leg feasibility test into a *level gate* that prunes a moon-pair leg
BEFORE the next leg's grid is built.

The per-leg criterion is NOT the ΔV — it is the geometric/physics feasibility
already computed elsewhere:

* :func:`cyclerfinder.search.vilm.vilm_dv_floor` — the admissible VILM ΔV-floor
  (escape+capture insertion cost; #76 deviation, *not* the no-GA quadrature),
* :func:`cyclerfinder.search.tisserand.linkable` — the Jovicentric constant-V∞
  contour intersection (centre-aware via ``mu=PRIMARIES[primary]``),
* :func:`cyclerfinder.search.correct._max_bend_deg` — bend feasibility (can the
  flyby turn at all).

The objective is a sum of non-negative per-leg terms (Bellman), so a per-leg gate
never discards a globally-feasible candidate: one dead leg kills the prefix.
Every verdict carries a ``reason`` string so the prune is *recorded* (the
Phase-4 empty-region report consumes the reason), never silent.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from cyclerfinder.core.constants import AU_KM
from cyclerfinder.core.satellites import PRIMARIES, SATELLITES
from cyclerfinder.search.correct import _max_bend_deg
from cyclerfinder.search.tisserand import linkable
from cyclerfinder.search.vilm import vilm_dv_floor

if TYPE_CHECKING:
    from cyclerfinder.data.discover_novel import TopologySpec

# Padding factor on the moon-pair SMA window handed to ``linkable``'s
# ``a_range_au`` (the contour-intersection search box in AU). The two moons'
# about-primary SMAs bracket the transfer; pad by 25% on each side so a transfer
# ellipse whose apo/peri straddle the pair is not boxed out.
_A_RANGE_PAD = 1.25


def _a_range_au(moon_a: str, moon_b: str) -> tuple[float, float]:
    """``(a_min, a_max)`` in AU spanning the two moons' about-primary SMAs, padded.

    ``linkable`` filters contour branches to this box; the moons' SMAs are in km
    (about the primary), so convert to AU and pad so the transfer ellipse is not
    boxed out at the edges.
    """
    a_a = SATELLITES[moon_a].sma_km
    a_b = SATELLITES[moon_b].sma_km
    lo = min(a_a, a_b) / _A_RANGE_PAD
    hi = max(a_a, a_b) * _A_RANGE_PAD
    return lo / AU_KM, hi / AU_KM


def moon_leg_admissible(
    moon_a: str,
    moon_b: str,
    *,
    vinf_kms: float,
    budget_kms: float,
    primary: str,
) -> tuple[bool, str]:
    """Is the single moon-pair leg ``moon_a`` -> ``moon_b`` admissible?

    Combines the three sourced per-leg gates and returns ``(admissible, reason)``.
    The leg is admissible iff ALL three pass:

    1. ``vilm_dv_floor(moon_a, moon_b) <= budget_kms`` — the VILM leveraging
       ΔV-floor is within budget (#76 escape+capture admissible floor).
    2. ``linkable(...)`` truthy — the Jovicentric constant-V∞ contours intersect
       (a flyby sequence between the two moons at this energy is possible).
    3. ``_max_bend_deg(vinf_kms, moon_b) > 0.0`` — the flyby can turn at all.

    The ``reason`` records WHY: the first failing gate, or ``"admissible"``.
    """
    floor = vilm_dv_floor(moon_a, moon_b)
    if floor > budget_kms:
        return False, f"vilm floor {floor:.3f} km/s > budget {budget_kms:.3f} km/s"

    mu = PRIMARIES[primary]
    if not linkable(
        moon_a,
        moon_b,
        vinf_kms,
        a_range_au=_a_range_au(moon_a, moon_b),
        mu=mu,
    ):
        return False, f"not linkable (Jovicentric contours disjoint at vinf {vinf_kms:.3f})"

    bend = _max_bend_deg(vinf_kms, moon_b)
    if bend <= 0.0:
        return False, f"max bend {bend:.4f} deg <= 0 at {moon_b} (vinf {vinf_kms:.3f})"

    return True, f"admissible (vilm floor {floor:.3f}<=budget, linkable, bend {bend:.3f} deg)"


def prune_topology_legs(
    spec: TopologySpec,
    *,
    vinf_seed_kms: float,
    budget_kms: float,
    primary: str,
) -> tuple[bool, tuple[str, ...]]:
    """Incremental per-leg box gate over a topology (the Ceriotti level-gate).

    Walk ``spec.sequence`` leg by leg; for each consecutive moon pair call
    :func:`moon_leg_admissible`. The topology survives only if EVERY leg is
    admissible (back-pruning: one dead leg kills the prefix). Returns
    ``(survives, per_leg_reasons)`` with one reason per leg — the prune is
    recorded, not silent, so an empty-region report can quote why a topology was
    dropped.
    """
    reasons: list[str] = []
    survives = True
    for moon_a, moon_b in zip(spec.sequence[:-1], spec.sequence[1:], strict=False):
        ok, reason = moon_leg_admissible(
            moon_a,
            moon_b,
            vinf_kms=vinf_seed_kms,
            budget_kms=budget_kms,
            primary=primary,
        )
        reasons.append(f"{moon_a}->{moon_b}: {reason}")
        if not ok:
            survives = False
    return survives, tuple(reasons)


__all__ = ["moon_leg_admissible", "prune_topology_legs"]
