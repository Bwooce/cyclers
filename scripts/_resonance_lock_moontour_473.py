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

from cyclerfinder.core.satellites import PRIMARIES, SATELLITES
from cyclerfinder.data.validation.v2_moontour import (
    V2_MOONTOUR_DRIFT_FLOOR_KMS,
    run_v2_moontour,
)
from cyclerfinder.parallel import ParallelSweepConfig, parallel_sweep
from cyclerfinder.search.cache_warm import warm_moon_leg_caches
from cyclerfinder.search.discovery_campaign import _mean_motion_rad_day
from cyclerfinder.search.releg_solver import MultiRevLeveragingReleg

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
from _drift_shape_473 import BOUNDED_SHAPES, classify_drift_shape  # noqa: E402

RUNLOG = REPO / "out" / "resonance_lock_473_runlog.jsonl"
SUMMARY = REPO / "out" / "resonance_lock_473_summary.json"

DAY_S = 86400.0
N_CYCLES_SEARCH = 10  # score on the multi-cycle drift SHAPE over 10 cycles

# #339 quasi-cycler reference envelope is max drift ~530k km. A resonance-lock
# survivor must be bounded-in-shape AND have an envelope no worse than ~1.5x the
# #339 reference (we use 800k km as the quasi_cycler ceiling so we don't admit a
# "bounded-but-enormous" excursion as if it were #339-quality). Strict cyclers
# (max drift <= 50k) are flagged separately as the bigger prize.
QUASI_ENVELOPE_CEILING_KMS = 800_000.0
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
PHASE0_DEG: tuple[float, ...] = (0.0, 90.0, 180.0, 270.0)


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
    drift_shape: str
    per_cycle_drift_kms: tuple[float, ...]
    max_resid_kms: float
    powered_dv_kms: float | None
    is_bounded_quasi: bool
    is_strict_cycler: bool
    extra: dict = field(default_factory=dict)


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
    is_bounded_quasi = bool(
        completed_all
        and cls["shape"] in BOUNDED_SHAPES
        and not is_strict
        and verdict.max_drift_kms <= QUASI_ENVELOPE_CEILING_KMS
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
        drift_shape=str(cls["shape"]),
        per_cycle_drift_kms=drifts,
        max_resid_kms=verdict.max_closure_residual_kms,
        powered_dv_kms=verdict.powered_total_dv_kms,
        is_bounded_quasi=is_bounded_quasi,
        is_strict_cycler=is_strict,
    )


def _ser(r: CandidateResult) -> dict:
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

    # Build the full cell list (all skeletons) for one parallel sweep.
    cells: list[Cell] = []
    for i, sk in enumerate(SKELETONS):
        cells.extend(_candidates(sk, i))
    n_total = len(cells)
    print(f"[{_now()}] {n_total} resonance-locked candidates across {len(SKELETONS)} skeletons")

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

    t0 = time.time()
    sweep = parallel_sweep(cells, _eval_cell, config=cfg)
    elapsed = time.time() - t0
    print(
        f"[{_now()}] sweep done: {sweep.n_succeeded}/{sweep.n_cells} ok, "
        f"{sweep.n_failed} failed, {elapsed:.1f}s"
    )
    if sweep.notes:
        print(f"  notes: {sweep.notes}")

    best_per_sk: dict[str, CandidateResult] = {}
    quasi_survivors: list[CandidateResult] = []
    strict_survivors: list[CandidateResult] = []
    for item_id, (cell, res) in enumerate(zip(cells, sweep.results, strict=True), start=1):
        if res is None:
            continue
        sk = SKELETONS[cell.sk_index]
        prev = best_per_sk.get(sk.label)
        # "best" = smallest envelope among bounded-shape candidates, else smallest drift.
        if prev is None or res.max_drift_kms < prev.max_drift_kms:
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
                    "elapsed_s": round(elapsed, 1),
                    "detail": _ser(res),
                },
                ensure_ascii=True,
            )
            + "\n"
        )
    runlog.flush()
    runlog.close()

    # Rank quasi survivors by envelope (closest to / below #339 reference first).
    quasi_survivors.sort(key=lambda r: r.max_drift_kms)
    strict_survivors.sort(key=lambda r: r.max_drift_kms)

    for label, r in best_per_sk.items():
        tag = (
            "STRICT-CYCLER"
            if r.is_strict_cycler
            else ("BOUNDED-QUASI" if r.is_bounded_quasi else f"{r.drift_shape}")
        )
        print(
            f"  {label:26s} best env={r.max_drift_kms:>12,.0f} km "
            f"(mult={r.syn_multiple}, ph0={r.phase0_deg}, rel={r.rel_offset_deg}) {tag}"
        )

    summary = {
        "git_sha": git_sha,
        "run_date": "2026-06-26",
        "n_cycles_search": N_CYCLES_SEARCH,
        "strict_floor_kms": V2_MOONTOUR_DRIFT_FLOOR_KMS,
        "quasi_envelope_ceiling_kms": QUASI_ENVELOPE_CEILING_KMS,
        "silver_339_envelope_kms": SILVER_339_ENVELOPE_KMS,
        "n_candidates": n_total,
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
        f"bounded quasi (env<=800k): {len(quasi_survivors)}"
    )
    print(f"summary -> {SUMMARY}")


if __name__ == "__main__":
    main()
