"""#473 resonance-lock moon-tour search (run-and-report) — CORRECTED.

Background
---------
#468 found 10 tours that CLOSE IN-BAND (cheap single-window powered close
< 3.5 km/s/cycle). #470 tried to admit them and got 0/10 against the strict
50,000 km V2 drift floor. STEP 0 of #473 re-adjudicated those 10 by drift-series
SHAPE and found 9/10 are bounded-oscillating (not monotonic-divergence) — but
with envelopes of 800k-3.7M km, far larger than the #339 quasi-cycler reference
envelope (~530k km). The #470 closes are GEOMEAN-ToF (no resonance lock).

The #339 existence proof (Umbriel-Oberon-Umbriel, SILVER_TOF=14.94 d/leg) is a
quasi_cycler because its cycle is RESONANCE-LOCKED to the encounter synodic
period: after one full cycle the encounter moons have advanced ~5 whole synodic
revolutions and return near their cycle-0 geometry -> bounded ~86k-530k km drift.

This search imposes the RESONANCE-LOCK CONDITION DURING the search: it sweeps
ToFs where the per-leg span is a (half-)integer multiple of the leg's encounter
synodic period, and scores each candidate by the BOUNDED-vs-DIVERGENT drift
SHAPE over n_cycles=10 (NOT by the strict floor only — that was the #470 bug),
AND records the drift ENVELOPE (max drift). A quasi_cycler survivor = bounded
oscillating shape with an envelope at or below a #339-referenced ceiling.
Candidates that ALSO pass the strict 50k floor are flagged as STRICT cyclers
(the bigger prize). The #465 MultiRevLeveragingReleg provides leg dV.

Positive control: run by scripts/_posctl_473_check.py (HARD GATE, already PASSED
— #339 reproduces its bounded-oscillating-and-returns 86k-530k km signature).

Outputs (per-unit JSONL runlog, one flushed line per candidate):
  out/resonance_lock_473_runlog.jsonl
  out/resonance_lock_473_summary.json
"""

from __future__ import annotations

import functools
import json
import math
import subprocess
import sys
import time
from collections.abc import Iterator
from dataclasses import dataclass, field
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

from cyclerfinder.core.satellites import PRIMARIES, SATELLITES
from cyclerfinder.data.validation.v2_moontour import (
    V2_MOONTOUR_DRIFT_FLOOR_KMS,
    run_v2_moontour,
)
from cyclerfinder.parallel import ParallelSweepConfig, parallel_sweep
from cyclerfinder.search.cache_warm import warm_moon_leg_caches
from cyclerfinder.search.discovery_campaign import _mean_motion_rad_day
from cyclerfinder.search.releg_solver import MultiRevLeveragingReleg
from cyclerfinder.search.vilm import vilm_dv_floor

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
from _drift_shape_473 import BOUNDED_SHAPES, classify_drift_shape  # noqa: E402
from _relative_drift_473 import (  # noqa: E402
    R_MARGIN,
    R_REF,
    normalising_moon,
    relative_drift_ratio,
)

RUNLOG = REPO / "out" / "resonance_lock_473_runlog.jsonl"
SUMMARY = REPO / "out" / "resonance_lock_473_summary.json"

DAY_S = 86400.0
N_CYCLES_SEARCH = 10  # score on the multi-cycle drift SHAPE over 10 cycles

# #477 admissible ΔV-floor pre-filter (branch-and-bound). Before the expensive
# per-candidate leveraging-chain solve (run_v2_moontour with the
# MultiRevLeveragingReleg, whose walk_vinf_down dominates the cost), screen each
# candidate's legs with the CHEAP, #472-memoized admissible lower bound
# cyclerfinder.search.vilm.vilm_dv_floor — the #76 escape+capture insertion cost,
# which is admissible (<= EVERY routing, no-GA and with-GA alike: leverage and
# finite phasing only ADD to it). A candidate is pruned iff ANY leg's irreducible
# VILM floor exceeds the powered-band per-cycle budget ceiling PRUNE_BUDGET_KMS:
# the per-cycle ΔV objective is a Bellman sum of non-negative per-leg costs, so
# one over-budget leg proves the whole tour out-of-band (one dead leg kills the
# prefix). Only the ΔV-floor is used — NOT the constant-V∞ linkable/bend gates of
# moon_prune.moon_leg_admissible: those are V∞-specific and (at the 6 km/s common
# target) wrongly report the wide-spaced Galilean/Saturnian/Uranian pairs as
# "contours disjoint", which would drop real candidates (the EGC + #339 controls).
# The ΔV-floor is V∞-band-independent and provably admissible, so it is the only
# sound branch-and-bound bound here.
#
# POSITIVE CONTROL (#477): the 2 known #470 flips (Galilean IEG, EGC) have per-leg
# floors <= 1.57 km/s < PRUNE_BUDGET_KMS, so the gate provably does NOT drop them
# — guarded in tests/search/test_resonance_lock_473_prefilter.py.
PRUNE_BUDGET_KMS = 3.5  # powered_dsm band per-cycle ceiling (V2-powered sanity max)


def _prefilter_cell(sk: Skeleton) -> tuple[bool, str]:
    """Cheap #477 admissible ΔV-floor pre-filter — survives iff EVERY leg's
    irreducible VILM floor is within the per-cycle powered budget.

    Returns ``(survives, reason)``. Self-legs (A->A phasing loops, no inter-moon
    transfer) carry no VILM floor and are skipped (always admissible). The reason
    records the FIRST failing leg (so the prune is auditable in the runlog) or
    ``"admissible"``.
    """
    for k in range(len(sk.sequence) - 1):
        a, b = sk.sequence[k], sk.sequence[k + 1]
        if a == b:
            continue
        floor = vilm_dv_floor(a, b)
        if floor > PRUNE_BUDGET_KMS:
            return False, f"{a}->{b}: vilm floor {floor:.3f} km/s > budget {PRUNE_BUDGET_KMS} km/s"
    return True, "admissible"


# RELATIVE (scale-invariant) envelope — see scripts/_relative_drift_473.py.
# A survivor must be bounded-in-shape AND have a RELATIVE drift ratio
# (max_drift / SMA_of_outermost_moon) no worse than the #339-calibrated bar
# R_REF~0.91 (margin R_MARGIN=1.0). An ABSOLUTE-km ceiling was WRONG: 530k km is
# Uranian-scale (Oberon SMA 583,511 km); applied to Jovian/Saturnian systems it
# mis-scales (Ganymede 1.07M, Callisto 1.88M, Titan 1.22M legitimately drift more
# in absolute km in a bigger system). Strict cyclers (max drift <= 50k) are
# flagged separately as the bigger prize.
SILVER_339_ENVELOPE_KMS = 530_000.0


def _git_sha() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"], cwd=REPO, text=True
        ).strip()
    except Exception:
        return "uncommitted"


def _now() -> str:
    return datetime.now(UTC).astimezone().isoformat(timespec="seconds")


def _period_days(moon: str, mu: float) -> float:
    s = SATELLITES[moon]
    return 2.0 * math.pi * math.sqrt(s.sma_km**3 / mu) / DAY_S


def _synodic_days(moon_a: str, moon_b: str, mu: float) -> float:
    na = _mean_motion_rad_day(mu, SATELLITES[moon_a].sma_km)
    nb = _mean_motion_rad_day(mu, SATELLITES[moon_b].sma_km)
    if abs(na - nb) < 1e-12:
        return math.inf
    return 2.0 * math.pi / abs(na - nb)


@dataclass(frozen=True)
class Skeleton:
    system: str
    sequence: tuple[str, ...]
    label: str


# The #468 sequences (which closed in-band but failed #470's strict floor) PLUS
# the #339 positive control as an in-search sentinel.
SKELETONS: tuple[Skeleton, ...] = (
    Skeleton("Uranus", ("Umbriel", "Oberon", "Umbriel"), "POSCTL Umbriel-Oberon (#339)"),
    Skeleton("Jupiter", ("Io", "Europa", "Ganymede", "Io"), "Galilean IEG"),
    Skeleton("Jupiter", ("Europa", "Ganymede", "Callisto", "Europa"), "Galilean EGC"),
    Skeleton("Jupiter", ("Io", "Europa", "Io"), "Jovian IE pair"),
    Skeleton("Jupiter", ("Ganymede", "Callisto", "Ganymede"), "Jovian GC pair"),
    Skeleton("Jupiter", ("Callisto", "Ganymede", "Europa", "Callisto"), "Galilean CGE (Liang)"),
    Skeleton("Saturn", ("Titan", "Rhea", "Dione", "Titan"), "Saturnian TRD"),
    Skeleton("Saturn", ("Rhea", "Dione", "Tethys", "Rhea"), "Saturnian RDT"),
    Skeleton("Saturn", ("Dione", "Tethys", "Enceladus", "Dione"), "Saturnian DTE"),
    Skeleton("Saturn", ("Rhea", "Dione", "Rhea"), "Saturnian RD pair"),
    Skeleton("Saturn", ("Tethys", "Enceladus", "Tethys"), "Saturnian TE pair"),
)

# Resonance-lock grid: per-leg ToF = mult * synodic period of the leg's two
# encounter moons (so the encounter geometry re-phases each cycle — the #339
# mechanism). #339's lock is T_leg ~= 2.5 * T_syn. Sweep half-integer AND integer
# multiples to catch both 2-moon (A-B-A) and 3-moon (A-B-C-A) locks.
SYNODIC_MULTIPLES: tuple[float, ...] = (
    0.5,
    1.0,
    1.5,
    2.0,
    2.5,
    3.0,
    3.5,
    4.0,
    4.5,
    5.0,
    6.0,
    7.0,
    8.0,
)
# Phase offset of the second distinct moon relative to the first (deg). #339's
# basin floor is rel_offset = 180 deg; sweep the circle.
PHASE_OFFSETS_DEG: tuple[float, ...] = tuple(float(d) for d in range(0, 360, 30))
# Anchor longitude of the first moon (absolute rotation of the whole tour).
# The rendezvous drift is ||r_final_k - r_final_0||; rotating the WHOLE tour by
# phase0 rotates r_final_k and r_final_0 equally, leaving the drift NORM exactly
# invariant (verified empirically: identical max_drift to the decimal across
# phase0 in {0,90,180,270} for a Jovian probe). phase0 is therefore a redundant
# axis for this metric — a single value suffices and cuts the grid 4x.
PHASE0_DEG: tuple[float, ...] = (0.0,)


@dataclass
class CandidateResult:
    label: str
    system: str
    sequence: tuple[str, ...]
    leg_tofs_days: tuple[float, ...]
    syn_multiple: float
    phase0_deg: float
    rel_offset_deg: float
    n_cycles_completed: int
    max_drift_kms: float
    min_drift_kms: float
    rel_drift_ratio: float
    norm_moon: str
    drift_shape: str
    per_cycle_drift_kms: tuple[float, ...]
    max_resid_kms: float
    powered_dv_kms: float | None
    is_bounded_quasi: bool
    is_strict_cycler: bool
    extra: dict[str, Any] = field(default_factory=dict)


def _leg_tofs_resonance_locked(
    sequence: tuple[str, ...], mu: float, mult: float
) -> tuple[float, ...]:
    out: list[float] = []
    for k in range(len(sequence) - 1):
        a, b = sequence[k], sequence[k + 1]
        t = mult * _period_days(a, mu) if a == b else mult * _synodic_days(a, b, mu)
        out.append(t)
    return tuple(out)


@dataclass(frozen=True)
class Cell:
    sk_index: int
    tofs: tuple[float, ...]
    syn_multiple: float
    phase0_deg: float
    rel_offset_deg: float


def _candidates(sk: Skeleton, sk_index: int) -> Iterator[Cell]:
    mu = PRIMARIES[sk.system]
    for mult in SYNODIC_MULTIPLES:
        tofs = _leg_tofs_resonance_locked(sk.sequence, mu, mult)
        if any(t <= 0.0 or not math.isfinite(t) for t in tofs):
            continue
        # Cap total cycle period so the search stays in a sane window.
        if sum(tofs) > 400.0:
            continue
        for ph0 in PHASE0_DEG:
            for rel in PHASE_OFFSETS_DEG:
                yield Cell(sk_index, tofs, mult, ph0, rel)


def _eval_cell(cell: Cell) -> CandidateResult:
    """Top-level closure (pickle-safe) — run the OFFICIAL gauntlet at 10 cycles."""
    sk = SKELETONS[cell.sk_index]
    n_legs = len(sk.sequence) - 1
    verdict = run_v2_moontour(
        candidate_id=f"{sk.label}-m{cell.syn_multiple}-ph{cell.phase0_deg}-rel{cell.rel_offset_deg}",
        sequence=sk.sequence,
        vinf_tuple_kms=(6.0,) * len(sk.sequence),  # carried only for audit
        leg_tofs_days=cell.tofs,
        rel_offset_deg=cell.rel_offset_deg,
        system=None,
        n_cycles=N_CYCLES_SEARCH,
        n_revs=tuple(0 for _ in range(n_legs)),
        phase0_deg=cell.phase0_deg,
        releg=MultiRevLeveragingReleg(),
        dv_band="powered_dsm",
    )
    drifts = tuple(float(c.rendezvous_drift_kms) for c in verdict.per_cycle)
    cls = classify_drift_shape(list(drifts[1:]))
    completed_all = verdict.n_cycles_completed >= N_CYCLES_SEARCH
    is_strict = bool(verdict.passes_v2)
    ratio = relative_drift_ratio(verdict.max_drift_kms, sk.sequence)
    is_bounded_quasi = bool(
        completed_all
        and cls["shape"] in BOUNDED_SHAPES
        and not is_strict
        and ratio <= R_MARGIN  # RELATIVE drift bar (replaces the absolute-km ceiling)
    )
    return CandidateResult(
        label=sk.label,
        system=sk.system,
        sequence=sk.sequence,
        leg_tofs_days=cell.tofs,
        syn_multiple=cell.syn_multiple,
        phase0_deg=cell.phase0_deg,
        rel_offset_deg=cell.rel_offset_deg,
        n_cycles_completed=verdict.n_cycles_completed,
        max_drift_kms=verdict.max_drift_kms,
        min_drift_kms=float(cls["min"]) if cls.get("min") is not None else math.inf,
        rel_drift_ratio=ratio,
        norm_moon=normalising_moon(sk.sequence),
        drift_shape=str(cls["shape"]),
        per_cycle_drift_kms=drifts,
        max_resid_kms=verdict.max_closure_residual_kms,
        powered_dv_kms=verdict.powered_total_dv_kms,
        is_bounded_quasi=is_bounded_quasi,
        is_strict_cycler=is_strict,
    )


def _ser(r: CandidateResult) -> dict[str, Any]:
    return {
        "label": r.label,
        "system": r.system,
        "sequence": list(r.sequence),
        "leg_tofs_days": [round(t, 3) for t in r.leg_tofs_days],
        "syn_multiple": r.syn_multiple,
        "phase0_deg": r.phase0_deg,
        "rel_offset_deg": r.rel_offset_deg,
        "n_cycles_completed": r.n_cycles_completed,
        "max_drift_kms": round(r.max_drift_kms, 1),
        "min_drift_kms": round(r.min_drift_kms, 1),
        "rel_drift_ratio": round(r.rel_drift_ratio, 4),
        "norm_moon": r.norm_moon,
        "drift_shape": r.drift_shape,
        "per_cycle_drift_kms": [round(d, 1) for d in r.per_cycle_drift_kms],
        "max_resid_kms": round(r.max_resid_kms, 4),
        "powered_dv_kms": (round(r.powered_dv_kms, 4) if r.powered_dv_kms is not None else None),
        "is_bounded_quasi": r.is_bounded_quasi,
        "is_strict_cycler": r.is_strict_cycler,
    }


def main() -> None:
    git_sha = _git_sha()
    RUNLOG.parent.mkdir(parents=True, exist_ok=True)
    runlog = RUNLOG.open("w", encoding="utf-8")

    # Build per-skeleton cell lists. We sweep ONE SKELETON AT A TIME and flush
    # the runlog after each, so (a) a process kill loses at most one skeleton's
    # work, not the whole 6864-cell run (the prior single-batch run was killed
    # exit-144 at 1120/6864 with NOTHING written), and (b) live results stay
    # bounded to one skeleton (~hundreds of cells) instead of accumulating all.
    cells_by_sk: list[list[Cell]] = [list(_candidates(sk, i)) for i, sk in enumerate(SKELETONS)]
    n_total = sum(len(c) for c in cells_by_sk)
    print(f"[{_now()}] {n_total} resonance-locked candidates across {len(SKELETONS)} skeletons")

    # #477 admissible ΔV-floor pre-filter (branch-and-bound). The per-leg VILM
    # floor / linkable / bend gate depends ONLY on the skeleton's moon-pairs (not
    # the phasing), so we screen ONCE per skeleton and skip the whole skeleton's
    # cells when a leg is inadmissible — pruning the expensive leveraging-chain
    # solve on candidates the cheap lower bound already proves out-of-band.
    prefilter: dict[int, tuple[bool, str]] = {
        i: _prefilter_cell(sk) for i, sk in enumerate(SKELETONS)
    }
    n_pruned_cells = sum(len(cells_by_sk[i]) for i, (ok, _) in prefilter.items() if not ok)
    for i, sk in enumerate(SKELETONS):
        ok, reason = prefilter[i]
        if not ok:
            print(f"[{_now()}] PRUNED skeleton {sk.label}: {reason} ({len(cells_by_sk[i])} cells)")
    print(
        f"[{_now()}] pre-filter: {n_pruned_cells}/{n_total} cells pruned cheaply, "
        f"{n_total - n_pruned_cells} survive to the leveraging-chain solve"
    )

    # #474 cross-worker cache: fork under multiprocessing with the parent's #472
    # caches pre-warmed over the campaign's leg/V∞/budget domain (COW-shared).
    vinf_grid = (4.0, 5.0, 6.0, 8.0)
    budget_grid = (3.5,)
    cfg = ParallelSweepConfig(
        backend="multiprocessing",
        prewarm=functools.partial(
            warm_moon_leg_caches, vinf_grid=vinf_grid, budget_grid=budget_grid
        ),
        verbose=5,
    )

    best_per_sk: dict[str, CandidateResult] = {}
    quasi_survivors: list[CandidateResult] = []
    strict_survivors: list[CandidateResult] = []
    item_id = 0
    n_ok = 0
    n_failed = 0
    t0 = time.time()
    for sk_index, sk in enumerate(SKELETONS):
        cells = cells_by_sk[sk_index]
        if not cells:
            continue
        ok, reason = prefilter[sk_index]
        if not ok:
            # Record the cheap prune for every cell (auditable runlog), skip solve.
            for _cell in cells:
                item_id += 1
                runlog.write(
                    json.dumps(
                        {
                            "item_id": item_id,
                            "sub_step": sk.label,
                            "result": "pruned_prefilter",
                            "prune_reason": reason,
                            "ts": _now(),
                            "k_of_N": f"{item_id}/{n_total}",
                            "elapsed_s": round(time.time() - t0, 1),
                        },
                        ensure_ascii=True,
                    )
                    + "\n"
                )
            runlog.flush()
            continue
        ts0 = time.time()
        sweep = parallel_sweep(cells, _eval_cell, config=cfg)
        n_ok += sweep.n_succeeded
        n_failed += sweep.n_failed
        print(
            f"[{_now()}] {sk.label}: {sweep.n_succeeded}/{sweep.n_cells} ok, "
            f"{sweep.n_failed} failed, {time.time() - ts0:.1f}s "
            f"(cum {time.time() - t0:.1f}s)"
        )
        for res in sweep.results:
            item_id += 1
            if res is None:
                continue
            prev = best_per_sk.get(sk.label)
            # "best" = smallest RELATIVE drift ratio (scale-invariant) for this skeleton.
            if prev is None or res.rel_drift_ratio < prev.rel_drift_ratio:
                best_per_sk[sk.label] = res
            if res.is_strict_cycler:
                strict_survivors.append(res)
            elif res.is_bounded_quasi:
                quasi_survivors.append(res)
            runlog.write(
                json.dumps(
                    {
                        "item_id": item_id,
                        "sub_step": sk.label,
                        "result": (
                            "strict_cycler"
                            if res.is_strict_cycler
                            else ("bounded_quasi" if res.is_bounded_quasi else res.drift_shape)
                        ),
                        "residual": round(res.max_drift_kms, 1),
                        "ts": _now(),
                        "k_of_N": f"{item_id}/{n_total}",
                        "elapsed_s": round(time.time() - t0, 1),
                        "detail": _ser(res),
                    },
                    ensure_ascii=True,
                )
                + "\n"
            )
        runlog.flush()  # checkpoint after each skeleton
    elapsed = time.time() - t0
    print(f"[{_now()}] sweep done: {n_ok}/{n_total} ok, {n_failed} failed, {elapsed:.1f}s")
    runlog.flush()
    runlog.close()

    # Rank survivors by RELATIVE drift ratio (closest to / below #339 first).
    quasi_survivors.sort(key=lambda r: r.rel_drift_ratio)
    strict_survivors.sort(key=lambda r: r.max_drift_kms)

    for label, r in best_per_sk.items():
        tag = (
            "STRICT-CYCLER"
            if r.is_strict_cycler
            else ("BOUNDED-QUASI" if r.is_bounded_quasi else f"{r.drift_shape}")
        )
        print(
            f"  {label:26s} best ratio={r.rel_drift_ratio:6.3f} "
            f"(env={r.max_drift_kms:>12,.0f} km norm={r.norm_moon} "
            f"mult={r.syn_multiple}, ph0={r.phase0_deg}, rel={r.rel_offset_deg}) {tag}"
        )

    summary = {
        "git_sha": git_sha,
        "run_date": "2026-06-26",
        "n_cycles_search": N_CYCLES_SEARCH,
        "strict_floor_kms": V2_MOONTOUR_DRIFT_FLOOR_KMS,
        "criterion": "shape=bounded AND rel_drift_ratio<=R_MARGIN",
        "r_ref": R_REF,
        "r_margin": R_MARGIN,
        "r_ref_provenance": "#339 peak 530000 km / Oberon SMA 583511 km = 0.9083",
        "silver_339_envelope_kms": SILVER_339_ENVELOPE_KMS,
        "n_candidates": n_total,
        "prefilter_budget_kms": PRUNE_BUDGET_KMS,
        "n_pruned_prefilter_cells": n_pruned_cells,
        "pruned_skeletons": {
            SKELETONS[i].label: reason for i, (ok, reason) in prefilter.items() if not ok
        },
        "best_per_skeleton": {k: _ser(v) for k, v in best_per_sk.items()},
        "n_strict_cyclers": len(strict_survivors),
        "n_bounded_quasi_survivors": len(quasi_survivors),
        "strict_cyclers": [_ser(s) for s in strict_survivors],
        "bounded_quasi_survivors": [_ser(s) for s in quasi_survivors[:50]],
        "elapsed_s": round(elapsed, 1),
    }
    SUMMARY.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(
        f"\n[{_now()}] strict cyclers: {len(strict_survivors)}; "
        f"bounded quasi (ratio<={R_MARGIN}): {len(quasi_survivors)}"
    )
    print(f"summary -> {SUMMARY}")


if __name__ == "__main__":
    main()
