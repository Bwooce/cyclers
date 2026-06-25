"""#436 — direct-e>0 ER3BP discovery runner (CR3BP-INDEPENDENT seeds).

#432 (continuation-from-CR3BP) and #435 (generated high-e Sun-planet seeds) both
START from a circular-problem family and continue in ``e``. By construction such
seeds CANNOT discover an "e>0-only" family — one with NO circular 3-body limit —
because they trace back to a CR3BP ancestor.

This runner closes that mechanism gap. It builds a blind symmetric x-axis-crossing
initial-condition grid placed DIRECTLY at the target eccentricity (no CR3BP
ancestry), forward-converges each seed, and then reverse-continues every converged
orbit toward e=0 to test for a circular limit:

  * ``cr3bp_continuous`` — the family continues back to e~0 (a CR3BP ancestor
    exists -> NOT an e>0-only discovery).
  * ``e_only_candidate`` — the family dies at ``death_e > floor`` (no CR3BP limit
    -> candidate e>0-only family; adjudicated against the literature corpus and
    probed with the saddle-center branch-switcher).
  * ``inconclusive`` — dies at ``death_e <= floor`` (numerical noise near e=0).

Systems probed at their real orbital eccentricity:
  * Earth-Moon   e=0.0549  (mu=0.012155)
  * Sun-Mercury  e=0.206
  * Sun-Mars     e=0.093
  * Sun-Pluto    e=0.249

Report-only — NO catalogue writeback. Per the literature_check discipline a
"not-found" literature status is NECESSARY-NOT-SUFFICIENT for novelty; the
V0-V5 gauntlet still governs.

Usage::

    uv run python scripts/run_436_direct_er3bp.py            # full 4-system 12x12
    CYCLERS_436_SMOKE=1 uv run python scripts/run_436_direct_er3bp.py   # 3x3 Earth-Moon only
"""

from __future__ import annotations

import contextlib
import json
import math
import os
import signal
import time
from collections.abc import Iterator, Sequence
from pathlib import Path

from joblib import Parallel, delayed

from cyclerfinder.core.cr3bp import cr3bp_system
from cyclerfinder.core.er3bp import ER3BPSystem
from cyclerfinder.genome.er3bp_branching import branch_at_saddle_center_er3bp
from cyclerfinder.genome.er3bp_periodic import ER3BPPeriodicOrbit
from cyclerfinder.search.er3bp_direct_seeding import (
    DirectEr3bpSeed,
    classify_no_cr3bp_limit,
    converge_direct_seed,
    direct_e_seed_grid,
)
from cyclerfinder.search.literature_check import (
    KNOWN_CORPUS,
    CandidateSignature,
    SearchResult,
    check_literature,
)

_DATA_DIR = Path(__file__).resolve().parents[1] / "data"

# Deliverable grid density. PARAMETERISED so the smoke run can shrink it to 3x3.
_N_X = 12
_N_YDOT = 12

# Blind symmetric-IC band: x0 spans the rotating frame between the primaries'
# vicinity and just inside the secondary; ydot0 spans prograde/retrograde.
_X0_RANGE = (0.1, 0.95)
_YDOT0_RANGE = (-4.0, 4.0)

# Both a full-period (2*pi) and a half-period (pi) symmetric guess are tried so a
# family whose true period is pi is not missed by a 2*pi-only span.
_PERIOD_FS: tuple[float, ...] = (2.0 * math.pi, math.pi)

# Earth-Moon mu per the plan (the #432/#435 Earth-Moon value); Sun-planet mu is
# read from cr3bp_system. Real orbital eccentricities.
_EARTH_MOON_MU = 0.012155
_SYSTEMS: tuple[tuple[str, str, float], ...] = (
    ("Earth", "Moon", 0.0549),
    ("Sun", "Mercury", 0.206),
    ("Sun", "Mars", 0.093),
    ("Sun", "Pluto", 0.249),
)


# Per-seed wall-clock budget (seconds). A blind grid hits pathological ICs (deep
# in a primary's well) where the single-shooting corrector grinds its full
# max_iter * backtrack ladder over a stiff propagation, stalling for minutes. The
# corrector exposes no max_iter through converge_direct_seed, so the runner caps
# each seed's converge+classify with a SIGALRM budget: a seed that blows the
# budget is treated as a non-convergence (skipped), not a campaign-wide hang.
_SEED_BUDGET_S = 20

# Parallel workers for the independent per-seed converge+classify tasks. Leave a
# few cores free (concurrent #434 sweep + BLAS/system overhead).
_N_JOBS = int(os.environ.get("CYCLERS_436_NJOBS", "12"))

# Soft per-seed wall-clock threshold (seconds) for pathological-state visibility.
# A seed slower than this is logged as productive-but-slow and NEVER terminated
# by the runner (the SIGALRM _SEED_BUDGET_S is a separate, intentional hard cap
# inside the worker; this constant only governs the parent's warning line).
_SLOW_UNIT_WARN_S = 1800.0


class _SeedBudgetError(Exception):
    """Raised when a single seed's converge+classify exceeds its wall-clock budget."""


@contextlib.contextmanager
def _seed_budget(seconds: int) -> Iterator[None]:
    """SIGALRM-based per-seed wall-clock guard (Linux main-thread only)."""

    def _handler(signum: int, frame: object) -> None:
        raise _SeedBudgetError

    old = signal.signal(signal.SIGALRM, _handler)
    signal.alarm(seconds)
    try:
        yield
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, old)


# --- Offline literature search (reused verbatim from run_435) -----------------
def _offline_corpus_hits() -> list[SearchResult]:
    hits: list[SearchResult] = []
    for anchor in KNOWN_CORPUS:
        bodies_str = " ".join(sorted(anchor.body_set))
        snippet = (
            f"{anchor.citation}. Authors: {', '.join(anchor.authors)}. "
            f"Keywords: {', '.join(anchor.keywords)}. Bodies: {bodies_str}. "
            f"Primary: {anchor.primary}. Subject: cycler trajectory mission design."
        )
        hits.append(
            SearchResult(
                title=anchor.name + " (cycler trajectory)",
                url=anchor.doi or f"https://example.org/{anchor.name.replace(' ', '_')}",
                snippet=snippet,
            )
        )
    return hits


_OFFLINE_HITS = _offline_corpus_hits()


def _tokenise(s: str) -> set[str]:
    return {t for t in "".join(c.lower() if c.isalnum() else " " for c in s).split() if len(t) > 2}


def offline_search(query: str) -> Sequence[SearchResult]:
    q_terms = _tokenise(query)
    out: list[tuple[int, SearchResult]] = []
    for r in _OFFLINE_HITS:
        text_terms = _tokenise(r.title + " " + r.snippet)
        overlap = len(q_terms & text_terms)
        if overlap >= 2:
            out.append((overlap, r))
    out.sort(key=lambda t: t[0], reverse=True)
    return [r for _, r in out]


def _print_progress(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def _build_system(primary: str, secondary: str, e: float) -> ER3BPSystem:
    """Build the ER3BP system: literal Earth-Moon mu, else cr3bp_system mu."""
    if (primary, secondary) == ("Earth", "Moon"):
        mu = _EARTH_MOON_MU
    else:
        mu = cr3bp_system(primary, secondary).mu
    return ER3BPSystem(mu=mu, e=e, primary_name=primary, secondary_name=secondary)


def _adjudicate_candidate(
    seed: DirectEr3bpSeed,
    orbit: ER3BPPeriodicOrbit,
) -> dict[str, object]:
    """Adjudicate an e_only_candidate: literature check + branch-switch probe."""
    # Structural fingerprint: the secondary body is the ER3BP perturber; the
    # primary is the central body. No resonance/vinf structure is asserted for a
    # raw periodic-orbit candidate.
    primary = seed.system.primary_name
    secondary = seed.system.secondary_name
    sig = CandidateSignature(
        primary=primary,
        sequence=(primary[0], secondary[0]),
        topology_label=frozenset({"resonant"}),
    )
    lit = check_literature(sig, search=offline_search)

    branch_orbit, branch_info = branch_at_saddle_center_er3bp(
        seed.system,
        orbit.state0,
        orbit.period_f / 2.0,
    )
    branched = branch_orbit is not None
    return {
        "literature_status": lit.status,
        "literature_citation": lit.citation,
        "branch_switched": branched,
        "branch_info": {k: str(v) for k, v in branch_info.items()},
    }


def _eval_one_seed(system: ER3BPSystem, seed: DirectEr3bpSeed) -> dict[str, object]:
    """Converge + classify ONE seed (the parallel work unit).

    Runs in a joblib worker process; the SIGALRM per-seed budget fires in the
    worker's own main thread. Returns a small picklable result dict; the parent
    aggregates the tally and writes records (so no shared mutable state).
    """
    t_unit = time.time()
    try:
        with _seed_budget(_SEED_BUDGET_S):
            orbit = converge_direct_seed(seed)
            if orbit is None:
                return {
                    "outcome": "nonconverged",
                    "label": seed.label,
                    "wall_s": time.time() - t_unit,
                }
            cls = classify_no_cr3bp_limit(orbit, system)
    except _SeedBudgetError:
        return {"outcome": "timeout", "label": seed.label, "wall_s": time.time() - t_unit}
    status = str(cls["status"])
    rec: dict[str, object] = {
        "label": seed.label,
        "primary": system.primary_name,
        "secondary": system.secondary_name,
        "mu": system.mu,
        "target_e": seed.target_e,
        "period_f_guess": seed.period_f,
        "x0": float(seed.state0[0]),
        "ydot0": float(seed.state0[4]),
        "status": status,
        "min_e": cls.get("min_e"),
        "death_e": cls.get("death_e"),
        "corrector_residual": orbit.corrector_residual,
        "source": seed.source,
    }
    if status == "e_only_candidate":
        rec.update(_adjudicate_candidate(seed, orbit))
    return {
        "outcome": "converged",
        "status": status,
        "record": rec,
        "label": seed.label,
        "wall_s": time.time() - t_unit,
    }


def _process_grid(
    system: ER3BPSystem,
    seeds: Sequence[DirectEr3bpSeed],
    progress: dict[str, object],
) -> tuple[list[dict[str, object]], dict[str, int]]:
    """Converge + classify a grid IN PARALLEL across _N_JOBS workers.

    Streams results in completion order (return_as="generator_unordered") and
    writes every converged record to the shared open file handle in
    ``progress["f"]`` IMMEDIATELY (with flush), so a kill mid-run preserves all
    completed seeds. ``progress`` carries the cross-grid unit counter ``k``, the
    campaign total ``n_total``, the campaign start ``t0`` and the file handle.
    """
    f = progress["f"]
    n_total = int(progress["n_total"])  # type: ignore[call-overload]
    t0 = float(progress["t0"])  # type: ignore[arg-type]

    gen = Parallel(n_jobs=_N_JOBS, backend="loky", return_as="generator_unordered")(
        delayed(_eval_one_seed)(system, seed) for seed in seeds
    )
    records: list[dict[str, object]] = []
    tally = {
        "grid_size": len(seeds),
        "converged": 0,
        "cr3bp_continuous": 0,
        "e_only_candidate": 0,
        "inconclusive": 0,
        "timed_out": 0,
    }
    while True:
        try:
            res = next(gen)
        except StopIteration:
            break
        except Exception as exc:  # a worker propagated; log + keep draining
            progress["k"] = int(progress["k"]) + 1  # type: ignore[call-overload]
            _print_progress(
                f"ERROR: unit {progress['k']}/{n_total} raised in worker "
                f"({exc!r}); logged, continuing with remaining seeds"
            )
            continue
        progress["k"] = int(progress["k"]) + 1  # type: ignore[call-overload]
        k = int(progress["k"])  # type: ignore[call-overload]
        try:
            wall_s = float(res.get("wall_s", float("nan")))  # type: ignore[union-attr]
            label = str(res.get("label", "?"))  # type: ignore[union-attr]
            outcome = res["outcome"]
            n_recs = 0
            if outcome == "timeout":
                tally["timed_out"] += 1
            elif outcome == "nonconverged":
                pass
            else:
                tally["converged"] += 1
                status = str(res["status"])
                tally[status] = tally.get(status, 0) + 1
                rec = res["record"]
                assert isinstance(rec, dict)
                if status == "e_only_candidate":
                    _print_progress(
                        f"  *** e_only_candidate: {rec['label']} "
                        f"({system.primary_name}-{system.secondary_name}) "
                        f"x0={float(rec['x0']):.4f} ydot0={float(rec['ydot0']):.4f} "  # type: ignore[arg-type]
                        f"death_e={rec.get('death_e')}"
                    )
                records.append(rec)
                # Incremental write: flush each converged record immediately so a
                # kill mid-run preserves everything already completed.
                f.write(json.dumps(rec) + "\n")  # type: ignore[union-attr]
                f.flush()  # type: ignore[union-attr]
                n_recs = 1

            elapsed = time.time() - t0
            avg = elapsed / k
            eta = (n_total - k) * avg
            _print_progress(
                f"unit {k}/{n_total} done | {label} | {outcome} | {n_recs} recs | "
                f"wall {wall_s:.1f}s | elapsed {elapsed:.1f}s | ETA {eta:.1f}s"
            )
            if wall_s > _SLOW_UNIT_WARN_S:
                _print_progress(
                    f"WARNING: unit {label} took {wall_s:.1f}s "
                    f"(>{_SLOW_UNIT_WARN_S}s) — productive-but-slow, NOT terminated"
                )
        except Exception as exc:  # one bad unit must not lose the good ones
            _print_progress(
                f"ERROR: unit {k}/{n_total} post-processing raised ({exc!r}); "
                "logged, continuing with remaining seeds"
            )
    return records, tally


def main() -> None:
    t0 = time.time()
    smoke = os.environ.get("CYCLERS_436_SMOKE") == "1"
    n_x = 3 if smoke else _N_X
    n_ydot = 3 if smoke else _N_YDOT
    systems = _SYSTEMS[:1] if smoke else _SYSTEMS

    _print_progress(
        f"#436 direct-e>0 ER3BP discovery "
        f"(grid={n_x}x{n_ydot}, periods={[round(p, 4) for p in _PERIOD_FS]}, "
        f"systems={len(systems)}, smoke={smoke})"
    )

    all_records: list[dict[str, object]] = []
    overall = {
        "grid_size": 0,
        "converged": 0,
        "cr3bp_continuous": 0,
        "e_only_candidate": 0,
        "inconclusive": 0,
        "timed_out": 0,
    }
    candidate_records: list[dict[str, object]] = []

    # Pre-build every (system, period_f) grid up front so the campaign-wide unit
    # total (n_total) is known before streaming — needed for the per-unit ETA.
    grids: list[tuple[ER3BPSystem, str, str, float, float, Sequence[DirectEr3bpSeed]]] = []
    for primary, secondary, e in systems:
        system = _build_system(primary, secondary, e)
        for period_f in _PERIOD_FS:
            seeds = direct_e_seed_grid(
                system,
                _X0_RANGE,
                _YDOT0_RANGE,
                n_x,
                n_ydot,
                period_f,
                is_half_period_residual=True,
            )
            grids.append((system, primary, secondary, e, period_f, seeds))
    n_total = sum(len(seeds) for *_, seeds in grids)

    # Open the output JSONL ONCE; records stream in and flush per-unit (a kill
    # mid-run preserves all completed seeds).
    # Smoke runs write to a *_smoke.jsonl path so a smoke-test never clobbers the
    # real campaign output under the canonical path.
    out_path = _DATA_DIR / ("er3bp_direct_436_smoke.jsonl" if smoke else "er3bp_direct_436.jsonl")
    out_path.parent.mkdir(parents=True, exist_ok=True)

    sys_tally_by_pair: dict[tuple[str, str], dict[str, int]] = {}
    with out_path.open("w") as f:
        progress: dict[str, object] = {"f": f, "k": 0, "n_total": n_total, "t0": t0}
        last_pair: tuple[str, str] | None = None
        for system, primary, secondary, e, period_f, seeds in grids:
            pair = (primary, secondary)
            if pair != last_pair:
                _print_progress(f"=== {primary}-{secondary} (mu={system.mu:.6g}, e={e}) ===")
                sys_tally_by_pair[pair] = {
                    "grid_size": 0,
                    "converged": 0,
                    "cr3bp_continuous": 0,
                    "e_only_candidate": 0,
                    "inconclusive": 0,
                    "timed_out": 0,
                }
                last_pair = pair
            sys_tally = sys_tally_by_pair[pair]
            _print_progress(
                f"{primary}-{secondary}: period_f={period_f:.4f} grid -> "
                f"{len(seeds)} seeds; converging..."
            )
            records, tally = _process_grid(system, seeds, progress)
            for key in sys_tally:
                sys_tally[key] += tally.get(key, 0)
            all_records.extend(records)
            candidate_records.extend(r for r in records if r["status"] == "e_only_candidate")

    for pair, sys_tally in sys_tally_by_pair.items():
        primary, secondary = pair
        _print_progress(
            f"{primary}-{secondary} tally: grid={sys_tally['grid_size']} "
            f"converged={sys_tally['converged']} "
            f"cr3bp_continuous={sys_tally['cr3bp_continuous']} "
            f"e_only_candidate={sys_tally['e_only_candidate']} "
            f"inconclusive={sys_tally['inconclusive']} "
            f"timed_out={sys_tally['timed_out']}"
        )
        for key in overall:
            overall[key] += sys_tally[key]

    _print_progress(f"Wrote {out_path.relative_to(_DATA_DIR.parent)} ({len(all_records)} records)")

    _print_progress(
        f"OVERALL tally: grid={overall['grid_size']} converged={overall['converged']} "
        f"cr3bp_continuous={overall['cr3bp_continuous']} "
        f"e_only_candidate={overall['e_only_candidate']} "
        f"inconclusive={overall['inconclusive']} "
        f"timed_out={overall['timed_out']}"
    )
    if candidate_records:
        _print_progress(f"!!! {len(candidate_records)} e_only_candidate(s) FLAGGED:")
        for rec in candidate_records:
            _print_progress(
                f"  -> {rec['label']} ({rec['primary']}-{rec['secondary']}) "
                f"death_e={rec.get('death_e')} lit={rec.get('literature_status')} "
                f"branch_switched={rec.get('branch_switched')}"
            )
    else:
        _print_progress(
            "No e_only_candidate found at this grid: every converged family has a "
            "CR3BP limit or is numerical noise (registry-grade negative)."
        )

    _print_progress(f"Campaign complete in {time.time() - t0:.1f}s")


if __name__ == "__main__":
    main()
