# src/cyclerfinder/search/leveraging_chain.py
"""Resonant-hop V∞ descent chain — the VILM endgame orchestrator (#465).

Where a single retarget leg (``DsmReleg`` / ``LowThrustReleg``) sheds a leg's
whole V∞ defect in ONE impulse (paying the single-VILM *maximum*, Eq.14), a
multi-rev leveraging endgame walks V∞ DOWN across many revolutions with a CHAIN
of resonant gravity-assist legs, each shedding a small slice (Campagnola-Russell
"The Endgame Problem" Part-1, digest ``2026-06-05-endgame-tisserand-mining.md``).
The total approaches the Eq.(13) quadrature *minimum* — roughly an order of
magnitude cheaper than the single-impulse shed.

This module is the **chain orchestrator** only: it chooses the per-hop resonance
``(n:m)`` by a greedy "zigzag low" descent (mining note p.11) and sums the per-hop
apse-burn ΔV. Each hop IS one existing :func:`cyclerfinder.search.leveraging_leg.
evaluate_leveraging_leg` call (#179, the phase-full apse VILM with its Γ-floor
cross-check). No new cost model, no new optimiser — the primitive and the floor
both already exist and are golden.

FIDELITY CAVEAT (inherited from ``leveraging_leg``): a converged chain is a
geometrically valid, Γ-floor-respecting V∞-lowering sequence — a conservative
lower bound on the realised endgame ΔV, not yet a phasing-closed trajectory (the
real return-to-moon phasing is an n-body / downstream concern). The chain reports
``converged=False`` rather than fabricating a ΔV when no feasible hop advances.
"""

from __future__ import annotations

from dataclasses import dataclass

from cyclerfinder.search import vilm
from cyclerfinder.search.leveraging_leg import LeveragingLegResult, evaluate_leveraging_leg

# Resonance grid the greedy descent searches at each hop. (n moon revs : m SC
# revs); the small integers span the leveraging graph the paper uses (Part-1
# Fig. 5 resonances are all low-order). Kept modest so the per-hop search is cheap
# and the realised chain stays near the continuous quadrature.
_RESONANCES: tuple[tuple[int, int], ...] = tuple((n, m) for n in range(1, 13) for m in range(1, 13))

# Default V∞ step per hop (km/s). Smaller → closer to the continuous floor but
# more hops; 0.05 reproduces the published Europa endgame (~137 m/s vs the 128 m/s
# continuous floor and 154 m/s published 3-VILM discrete) inside the finite-chain
# band.
_DEFAULT_STEP_KMS: float = 0.05


@dataclass(frozen=True)
class ChainResult:
    """A resonant-hop V∞ descent at one moon. CONSTRAINED vs EMERGED separated.

    ``converged=True`` asserts only that every hop is a geometrically valid,
    Γ-floor-respecting V∞-lowering leg and the chain reached the requested target
    within tol — the same conservative-lower-bound caveat as
    :class:`cyclerfinder.search.leveraging_leg.LeveragingLegResult`. It does NOT
    assert phasing closure of the whole multi-rev sequence (an n-body concern).
    """

    total_dv_kms: float  # CONSTRAINED — Σ per-hop apse-burn ΔV
    hops: tuple[LeveragingLegResult, ...]  # the realised hop sequence (EMERGED)
    vinf_end_kms: float  # EMERGED — V∞ the chain actually reached
    total_revs: int  # EMERGED — Σ hop n_moon_revs (the ToF proxy)
    converged: bool


def _infeasible(vinf_end_kms: float, hops: tuple[LeveragingLegResult, ...]) -> ChainResult:
    """A stalled descent — never a fabricated ΔV (orbit-closure discipline)."""
    return ChainResult(
        total_dv_kms=float("nan"),
        hops=hops,
        vinf_end_kms=vinf_end_kms,
        total_revs=sum(h.resonance[0] for h in hops),
        converged=False,
    )


def walk_vinf_down(
    moon: str,
    vinf_hi_kms: float,
    vinf_lo_kms: float,
    *,
    exterior: bool = True,
    max_hops: int = 200,
    max_revs: int = 2000,
    step_kms: float = _DEFAULT_STEP_KMS,
    vinf_tol_kms: float = 1.0e-3,
) -> ChainResult:
    """Walk V∞ from ``vinf_hi_kms`` down to ``vinf_lo_kms`` via resonant hops.

    Greedy "zigzag low" descent: at each step target ``v - step`` (clamped at the
    low target) and pick the CHEAPEST feasible
    :func:`cyclerfinder.search.leveraging_leg.evaluate_leveraging_leg` hop whose
    achieved V∞ steps strictly toward ``vinf_lo_kms`` and which respects its Γ
    floor. Sum the per-hop ΔV. Stop at ``vinf_lo`` (within ``vinf_tol_kms``), or —
    reporting INFEASIBLE (``converged=False``) — at ``max_hops`` / ``max_revs``, or
    when no feasible hop advances.

    Leveraging is only physical/efficient above the moon's efficiency threshold
    V̄∞ (Eq.9, Part-1 Table 3): a requested ``vinf_lo`` below V̄∞ is structurally
    unreachable by leveraging and returns infeasible immediately (chaining walks
    V∞ within the efficient regime, it cannot manufacture sub-threshold V∞).
    """
    if vinf_lo_kms >= vinf_hi_kms - vinf_tol_kms:
        # No walk needed (target at or above the natural arrival V∞).
        return ChainResult(
            total_dv_kms=0.0,
            hops=(),
            vinf_end_kms=vinf_hi_kms,
            total_revs=0,
            converged=True,
        )

    vbar = vilm.min_vinf_for_vilm(moon, exterior=exterior)
    if vinf_lo_kms < vbar - vinf_tol_kms:
        # Below the efficiency threshold — leveraging cannot reach it.
        return _infeasible(vinf_hi_kms, ())

    hops: list[LeveragingLegResult] = []
    total_dv = 0.0
    total_revs = 0
    v = vinf_hi_kms
    while v - vinf_lo_kms > vinf_tol_kms:
        v_target = max(vinf_lo_kms, v - step_kms)
        best: LeveragingLegResult | None = None
        for n, m in _RESONANCES:
            leg = evaluate_leveraging_leg(
                moon=moon,
                n_moon_revs=n,
                m_sc_revs=m,
                vinf_in_kms=v,
                vinf_out_target_kms=v_target,
                exterior=exterior,
            )
            if not (leg.converged and leg.gamma_floor_ok):
                continue
            # The hop must STRICTLY advance toward the low target.
            if leg.vinf_out_kms >= v - 1.0e-9:
                continue
            if best is None or leg.dv_dsm_kms < best.dv_dsm_kms:
                best = leg
        if best is None:
            return _infeasible(v, tuple(hops))
        hops.append(best)
        total_dv += best.dv_dsm_kms
        total_revs += best.resonance[0]
        v = best.vinf_out_kms
        if len(hops) >= max_hops or total_revs >= max_revs:
            # Pathological deep walk — a feasibility kill, not a fabricated close.
            if v - vinf_lo_kms > vinf_tol_kms:
                return _infeasible(v, tuple(hops))
            break

    return ChainResult(
        total_dv_kms=total_dv,
        hops=tuple(hops),
        vinf_end_kms=v,
        total_revs=total_revs,
        converged=True,
    )
