"""Powered moon-tour driver — close a cycle with a swappable releg backend (#449).

The driver loops the legs of a tour skeleton with a chosen
:class:`cyclerfinder.search.releg_solver.Releg` backend and returns a
:class:`PoweredCycleVerdict`. The closure is the leg-swap of the ballistic
``_close_one_phasing`` (``discovery_campaign``): instead of a near-zero ballistic
continuity residual, a *powered* leg DELIVERS a budgeted ΔV so each interior
flyby is V_inf-continuous AFTER the maneuver, and the cycle is scored on total
delivered ΔV vs the powered dv-band.

Cheapest-first prefilter (the C21 lesson, design §4.2)
------------------------------------------------------
Before any (expensive) DSM solve the driver runs the existing sourced
VILM/Tisserand/bend prefilter
(:func:`cyclerfinder.search.moon_prune.moon_leg_admissible`). If any leg is
unbridgeable at *every* probed V_inf — the Uranian disjoint-Tisserand-contour
case — the skeleton is reported EMPTY (``prefilter_skipped=True``) WITHOUT running
the powered solve. This reproduces the structural negative
``uranus-neptune-regular-moon-endgame-vilm-2026-06-23`` honestly: a powered leg
cannot bridge contours that are disjoint at all V_inf, so the genome does not
fabricate a bridge.

Continuity model
----------------
A flyby rotates V_inf but preserves its magnitude, so a closed cycle is
V_inf-continuous iff every leg departs AND arrives at a common V_inf magnitude
``T`` at the shared flyby body. The driver picks ``T`` (minimising total ΔV over a
probe band) and asks each powered leg to PIN its departure V_inf magnitude to
``T`` and RETARGET its arrival V_inf to ``T`` — so continuity holds by
construction (the closed-cycle wrap included). The total delivered ΔV is the sum
of the per-leg DSM impulses; it is classified against the powered dv-band
(:mod:`cyclerfinder.verify.dv_band_acceptance`).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from cyclerfinder.core.satellites import PRIMARIES, SATELLITES
from cyclerfinder.data.empty_regions import EmptyRegionReport
from cyclerfinder.data.method_capability import MethodCapability
from cyclerfinder.search.discovery_campaign import DAY_S, _mean_motion_rad_day, _moon_state
from cyclerfinder.search.moon_prune import moon_leg_admissible
from cyclerfinder.search.releg_solver import Releg, RelegResult

# V_inf magnitudes (km/s) the driver tries as the common continuity target ``T``
# and the prefilter linkability probe. Spans the registry's positive-control link
# (V_inf=4) up to the upper screen the Uranus probe used (15 km/s).
_VINF_PROBE_KMS: tuple[float, ...] = (4.0, 5.0, 6.0, 8.0, 10.0, 12.0, 15.0)

# Closed-cycle continuity acceptance (the ballistic gate the powered retarget must
# beat AFTER the maneuver, ``CampaignConfig.gate_residual_kms``).
_CONTINUITY_GATE_KMS: float = 0.05

# Prefilter budget (km/s): the VILM ΔV-floor admissibility ceiling. The powered
# band's per-cycle sanity ceiling (3.5 km/s/cycle) is the natural budget — a leg
# whose irreducible escape+capture floor already exceeds it is unbridgeable even
# by an optimal powered leg.
_PREFILTER_BUDGET_KMS: float = 3.5


@dataclass(frozen=True)
class PoweredCycleVerdict:
    """Outcome of closing one moon-tour cycle with a powered releg backend.

    Attributes
    ----------
    feasible:
        ``True`` iff every leg closed and the post-retarget continuity residual
        is below the ballistic gate. ``False`` on a prefilter skip or an
        infeasible leg.
    total_dv_kms:
        Sum of the per-leg delivered ΔV over the cycle (km/s). ``inf`` when
        prefiltered or infeasible.
    per_leg_dv_kms:
        Per-leg delivered ΔV (km/s); empty on a prefilter skip.
    continuity_residual_kms:
        Worst V_inf-magnitude continuity defect at any interior flyby + the
        closed-cycle wrap, AFTER the powered retarget (km/s). ``inf`` when
        prefiltered/infeasible.
    target_vinf_kms:
        The common flyby V_inf magnitude the legs were pinned/retargeted to
        (the continuity target minimising total ΔV). ``nan`` when prefiltered.
    dv_band:
        The measured dv-band classification of the total ΔV (from
        :func:`cyclerfinder.verify.dv_band_acceptance.classify_dv_band`), or
        ``None`` when prefiltered/infeasible.
    prefilter_skipped:
        ``True`` iff the VILM/linkability prefilter marked the skeleton
        unbridgeable and NO powered solve was run (the structural negative).
    prefilter_reasons:
        Per-leg prefilter reason strings (recorded, never silent).
    """

    feasible: bool
    total_dv_kms: float
    per_leg_dv_kms: tuple[float, ...]
    continuity_residual_kms: float
    target_vinf_kms: float
    dv_band: str | None
    prefilter_skipped: bool
    prefilter_reasons: tuple[str, ...]


def _leg_admissible_any_vinf(moon_a: str, moon_b: str, *, primary: str) -> tuple[bool, str]:
    """Is the leg admissible at ANY probed V_inf? Returns ``(ok, reason)``.

    A leg is bridgeable if it is admissible at at least one V_inf in the probe
    band; it is unbridgeable (the structural case) only if it fails at EVERY
    probed V_inf. The reason is the last failure (or the first success).
    """
    last_reason = ""
    for vinf in _VINF_PROBE_KMS:
        ok, reason = moon_leg_admissible(
            moon_a, moon_b, vinf_kms=vinf, budget_kms=_PREFILTER_BUDGET_KMS, primary=primary
        )
        if ok:
            return True, f"{moon_a}->{moon_b}: {reason}"
        last_reason = reason
    return False, f"{moon_a}->{moon_b}: {last_reason}"


def _moon_states(
    sequence: tuple[str, ...],
    leg_tofs_days: tuple[float, ...],
    phasing: dict[str, float],
    mu: float,
) -> list[tuple[np.ndarray, np.ndarray]]:
    """Planet-frame (pos, vel) at each encounter (mirrors ``_close_one_phasing``)."""
    epochs = [0.0]
    for tof in leg_tofs_days:
        epochs.append(epochs[-1] + tof)
    states: list[tuple[np.ndarray, np.ndarray]] = []
    for moon, t in zip(sequence, epochs, strict=True):
        sma = SATELLITES[moon].sma_km
        n_rad_day = _mean_motion_rad_day(mu, sma)
        states.append(_moon_state(phasing[moon], n_rad_day, t, sma, mu))
    return states


def _close_at_target(
    sequence: tuple[str, ...],
    states: list[tuple[np.ndarray, np.ndarray]],
    leg_tofs_days: tuple[float, ...],
    n_revs: tuple[int, ...],
    releg: Releg,
    mu: float,
    target: float | None,
) -> tuple[bool, float, float, tuple[float, ...]] | None:
    """Close every leg at a common flyby V_inf ``target``; return the cycle stats.

    Returns ``(feasible, total_dv, continuity_residual, per_leg_dv)`` or ``None``
    if any leg is infeasible. ``target`` pins each leg's departure V_inf magnitude
    and retargets its arrival V_inf, so continuity holds by construction; when
    ``target`` is ``None`` (the ballistic backend) the legs run unconstrained.
    """
    n_legs = len(sequence) - 1
    per_leg_dv: list[float] = []
    results: list[RelegResult] = []
    for k in range(n_legs):
        r_a, v_a = states[k]
        r_b, v_b = states[k + 1]
        res = releg.solve(
            r_a,
            v_a,
            r_b,
            v_b,
            leg_tofs_days[k] * DAY_S,
            mu,
            n_rev=n_revs[k],
            vinf_target_in=target,
            vinf_depart_mag=target,
        )
        if not res.feasible:
            return None
        results.append(res)
        per_leg_dv.append(res.dv_kms)

    # Continuity defect at every interior flyby + the closed-cycle anchor wrap
    # (the same definition as ``_close_one_phasing`` fix C), AFTER the retarget.
    worst = 0.0
    for k in range(1, n_legs):
        worst = max(worst, abs(results[k - 1].vinf_in - results[k].vinf_out))
    worst = max(worst, abs(results[0].vinf_out - results[-1].vinf_in))
    total_dv = float(sum(per_leg_dv))
    return True, total_dv, worst, tuple(per_leg_dv)


def close_powered_cycle(
    *,
    primary: str,
    sequence: tuple[str, ...],
    leg_tofs_days: tuple[float, ...],
    n_revs: tuple[int, ...],
    releg: Releg,
    phasing: dict[str, float],
    dv_band: str | None = "powered_dsm",
) -> PoweredCycleVerdict:
    """Close one moon-tour cycle with a powered (or ballistic) releg backend.

    Pipeline (design §4):

    1. **VILM/linkability prefilter.** For each consecutive moon pair, is the leg
       admissible at ANY probed V_inf? If any leg is unbridgeable at EVERY probed
       V_inf (disjoint Tisserand contours — the Uranus case), return EMPTY with
       ``prefilter_skipped=True`` and NO powered solve.
    2. **Powered close.** For each candidate common flyby V_inf ``T`` in the probe
       band, pin every leg's departure V_inf to ``T`` and retarget its arrival to
       ``T`` (continuity by construction); keep the ``T`` minimising total ΔV.
    3. **Classify.** Bin the total ΔV against the powered dv-band.

    ``dv_band`` is the catalogue band the caller expects (default
    ``"powered_dsm"``); it is carried for the audit trail. The verdict's
    ``dv_band`` field is the MEASURED classification of the achieved ΔV.

    Returns a :class:`PoweredCycleVerdict`. A closed cycle is necessary-not-
    sufficient: a hit must still clear the V2 moontour gauntlet + lit-novelty.
    """
    if len(sequence) < 2 or sequence[0] != sequence[-1]:
        raise ValueError("sequence must be a closed cycle (first == last, >= 1 leg)")
    n_legs = len(sequence) - 1
    if len(leg_tofs_days) != n_legs or len(n_revs) != n_legs:
        raise ValueError(
            f"leg_tofs_days/n_revs must each have {n_legs} entries, got "
            f"{len(leg_tofs_days)}/{len(n_revs)}"
        )

    # --- Stage 1: cheap prefilter (no DSM solve). ---
    reasons: list[str] = []
    bridgeable = True
    for k in range(n_legs):
        ok, reason = _leg_admissible_any_vinf(sequence[k], sequence[k + 1], primary=primary)
        reasons.append(reason)
        if not ok:
            bridgeable = False
    if not bridgeable:
        return PoweredCycleVerdict(
            feasible=False,
            total_dv_kms=math.inf,
            per_leg_dv_kms=(),
            continuity_residual_kms=math.inf,
            target_vinf_kms=math.nan,
            dv_band=None,
            prefilter_skipped=True,
            prefilter_reasons=tuple(reasons),
        )

    mu = PRIMARIES[primary]
    states = _moon_states(sequence, leg_tofs_days, phasing, mu)

    # --- Stage 2: powered close at the cheapest common flyby V_inf. ---
    # A powered backend (DsmReleg / LowThrustReleg) retargets each leg to a common
    # flyby V_inf so continuity holds by construction; the ballistic backend has no
    # retarget knob, so it runs once unconstrained. A backend is powered iff it is
    # NOT the regression-locked BallisticReleg — so LowThrustReleg swaps in for
    # DsmReleg with no driver change (the #449 Task 7 swap contract).
    from cyclerfinder.search.releg_solver import BallisticReleg

    targets: tuple[float | None, ...]
    targets = (None,) if isinstance(releg, BallisticReleg) else _VINF_PROBE_KMS

    # Select the common flyby V_inf ``T`` that CLOSES continuity (residual below
    # the gate) at minimum ΔV; if no probed ``T`` closes, keep the one with the
    # smallest continuity residual (so the verdict reports the honest best-effort
    # defect, which the caller's gate then fails). Sort key: (does-not-close,
    # continuity-if-not-closed, total-ΔV) — a closing target always beats a
    # non-closing one, and among closing targets the cheapest ΔV wins.
    best: tuple[bool, float, float, tuple[float, ...]] | None = None
    best_key: tuple[int, float, float] | None = None
    best_target = math.nan
    for target in targets:
        out = _close_at_target(sequence, states, leg_tofs_days, n_revs, releg, mu, target)
        if out is None:
            continue
        _, total_dv_t, continuity_t, _ = out
        closes = continuity_t < _CONTINUITY_GATE_KMS
        key = (0 if closes else 1, 0.0 if closes else continuity_t, total_dv_t)
        if best_key is None or key < best_key:
            best_key = key
            best = out
            best_target = target if target is not None else math.nan
    if best is None:
        return PoweredCycleVerdict(
            feasible=False,
            total_dv_kms=math.inf,
            per_leg_dv_kms=(),
            continuity_residual_kms=math.inf,
            target_vinf_kms=math.nan,
            dv_band=None,
            prefilter_skipped=False,
            prefilter_reasons=tuple(reasons),
        )

    _, total_dv, continuity, per_leg_dv = best
    feasible = continuity < _CONTINUITY_GATE_KMS
    band: str | None = None
    if feasible:
        from cyclerfinder.verify.dv_band_acceptance import classify_dv_band

        # One cycle's deterministic ΔV (m/s), classified on its own basis.
        band = classify_dv_band(total_dv * 1000.0, n_cycles=1)
    return PoweredCycleVerdict(
        feasible=feasible,
        total_dv_kms=total_dv,
        per_leg_dv_kms=per_leg_dv,
        continuity_residual_kms=continuity,
        target_vinf_kms=best_target,
        dv_band=band,
        prefilter_skipped=False,
        prefilter_reasons=tuple(reasons),
    )


# ---------------------------------------------------------------------------
# Capability-subsumption re-stamp (design §4.7 / honesty contract). A powered
# re-test that ALSO finds a region empty (the Uranus disjoint-contour case)
# produces a SUBSUMING method+version record suitable for appending to
# empty_regions.jsonl — but the writeback itself is a campaign-issue action, NOT
# a build action. This module only BUILDS the record (the test asserts its
# shape); it never writes it.
# ---------------------------------------------------------------------------

# The powered-releg capability envelope: it strictly subsumes the ballistic
# repeated-moon lane AND the phase-full VILM endgame lane (it adds the in-leg
# one-DSM impulse DOF on top of the leveraging genome). Tags are drawn from the
# method_capability partial order (``data.method_capability._CAPABILITY_EDGES``):
# ``one-dsm-per-leg`` ⊐ ``single-arc``, ``powered`` ⊐ ``ballistic``,
# ``leveraging`` ⊐ ``single-arc``.
_RELEG_CAPABILITY_TAGS = frozenset(
    {"powered", "one-dsm-per-leg", "leveraging", "multi-arc", "patched-conic", "coplanar"}
)


def powered_releg_method_capability(*, git_sha: str = "uncommitted") -> MethodCapability:
    """The DSM powered-releg capability envelope (capability-subsumption record).

    The genome adds a one-DSM-per-leg in-leg impulse to the leveraging moon-tour
    lane, so it strictly subsumes both the ballistic repeated-moon lane and the
    phase-full VILM endgame lane (it reaches every tour they reach, plus the
    powered ones). ``git_sha`` is the producing commit (reproducibility, NOT the
    capability — that is the tag set).
    """
    return MethodCapability(
        genome="DSM powered releg (one-DSM-per-leg V_inf-retarget moon tour, #449)",
        corrector="releg_moontour.close_powered_cycle (DsmReleg, VILM-floor prefilter)",
        capability_tags=_RELEG_CAPABILITY_TAGS,
        git_sha=git_sha,
    )


def build_powered_empty_restamp(
    *,
    region_id: str,
    family: str,
    centre: str,
    sequence: tuple[str, ...],
    verdict: PoweredCycleVerdict,
    git_sha: str = "uncommitted",
    run_date: str = "",
) -> EmptyRegionReport:
    """Build (NOT write) a capability-subsumption re-stamp for a powered-empty region.

    When a powered re-test of a moon-tour dead region ALSO finds it empty — the
    Uranian disjoint-Tisserand-contour case (``verdict.prefilter_skipped`` true),
    or a feasible close whose ΔV exceeds the powered band — the negative is
    re-stamped with the SUBSUMING powered method+version (a powered emptiness is a
    stronger, more final negative than a ballistic one; design §1.2 honesty
    caveat). This builds the :class:`EmptyRegionReport` for that re-stamp; the
    ACTUAL append to ``empty_regions.jsonl`` is a campaign-issue action, not done
    here (the #449 plan ships the capability + golden; the campaign spends it).

    The record is :func:`validate_empty_region`-valid (bounded extent, recorded
    prune gates, non-empty capability tags).
    """
    if verdict.prefilter_skipped:
        reason = "all legs unbridgeable: disjoint Tisserand/resonance contours at all probed V_inf"
        n_probed = len(_VINF_PROBE_KMS)
    else:
        reason = (
            f"powered close exceeds the powered dv-band ceiling "
            f"(total {verdict.total_dv_kms:.3f} km/s/cycle)"
        )
        n_probed = len(_VINF_PROBE_KMS)
    return EmptyRegionReport(
        region_id=region_id,
        family=family,
        centre=centre,
        topologies=({"sequence": list(sequence), "period_k": 1},),
        method_capability=powered_releg_method_capability(git_sha=git_sha),
        search_extent={
            "vinf_seeds_probed_kms": list(_VINF_PROBE_KMS),
            "points_total": n_probed,
            "ephem_model": "circular",
            "center": centre,
        },
        prune_gates=(
            "linkable (Tisserand/resonance contour overlap)",
            "vilm_dv_floor<=budget",
            "powered dv-band ceiling",
        ),
        result={
            "prefilter_skipped": verdict.prefilter_skipped,
            "total_dv_kms": verdict.total_dv_kms,
            "continuity_residual_kms": verdict.continuity_residual_kms,
            "prefilter_reasons": list(verdict.prefilter_reasons),
        },
        verdict=(
            "EMPTY (STRUCTURAL, powered re-test) — disjoint contours; powered releg "
            "cannot bridge what is unlinkable at all V_inf"
            if verdict.prefilter_skipped
            else "EMPTY (powered) — closes only above the powered dv-band ceiling"
        ),
        interpretation=reason,
        source_anchors=(
            "capability-subsumption re-stamp of the prior ballistic/VILM negative (#449)"
        ),
        run={"date": run_date, "git_sha": git_sha, "method": "releg_moontour (DSM)"},
    )
