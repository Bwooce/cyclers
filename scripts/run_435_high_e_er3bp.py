"""#435 — high-e Sun-planet ER3BP discovery runner.

#432's ER3BP discovery campaign was seed-limited to Earth-Moon ICs (the Broucke
IC is Earth-Moon-μ-specific), so its high-e probes could only push a synthetic
Earth-Moon family — it could NOT reach the high-departure regime of genuinely
eccentric Sun-planet systems. This runner closes that gap: it GENERATES
rotating-frame CR3BP seeds (small-amplitude L1 Lyapunov + DRO) directly at the
high-e Sun-planet systems and feeds them to the already-built #432 ER3BP
pipeline (``continue_and_monitor`` / ``adjudicate_trace``), continuing each to
the body's real orbital eccentricity.

Systems probed (real eccentricities):
  * Sun-Mercury  e=0.206
  * Sun-Mars     e=0.093
  * Sun-Pluto    e=0.249

Report-only — NO catalogue writeback. Per the literature_check discipline a
"not-found" literature status is NECESSARY-NOT-SUFFICIENT for novelty; the
V0-V5 gauntlet still governs.

Usage::

    uv run python scripts/run_435_high_e_er3bp.py
"""

from __future__ import annotations

import json
import time
from collections.abc import Sequence
from pathlib import Path

from cyclerfinder.core.cr3bp import cr3bp_system
from cyclerfinder.core.er3bp import ER3BPSystem
from cyclerfinder.search.cr3bp_seed_generator import dro_seed, lyapunov_seed
from cyclerfinder.search.er3bp_discovery import (
    Er3bpContinuationTrace,
    Er3bpSeed,
    adjudicate_trace,
    continue_and_monitor,
)
from cyclerfinder.search.literature_check import KNOWN_CORPUS, SearchResult

_DATA_DIR = Path(__file__).resolve().parents[1] / "data"

# Fine stepping: the secant predictor dies on step size below ~n_steps=30.
# 60 keeps the predictor well inside its convergence basin.
_N_STEPS = 60
assert _N_STEPS >= 30, "continuation secant predictor dies on coarse steps"

# High-e Sun-planet systems probed at their real orbital eccentricity.
_SYSTEMS: tuple[tuple[str, str, float], ...] = (
    ("Sun", "Mercury", 0.206),
    ("Sun", "Mars", 0.093),
    ("Sun", "Pluto", 0.249),
)

_SEED_SOURCE = (
    "#435 generated CR3BP Lyapunov/DRO seed (lagrange linear seed + fixed-Jacobi corrector)"
)


# --- Offline literature search (reused verbatim from run_432) -----------------
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


def _record_for_seed(
    seed: Er3bpSeed,
    trace: Er3bpContinuationTrace,
    lit_status: str,
    lit_citation: str | None,
) -> dict[str, object]:
    return {
        "label": seed.label,
        "primary": seed.system.primary_name,
        "secondary": seed.system.secondary_name,
        "mu": seed.system.mu,
        "target_e": trace.target_e,
        "outcome": trace.outcome,
        "e_max_reached": trace.e_max_reached,
        "e_star": trace.e_star,
        "source": seed.source,
        "literature_status": lit_status,
        "literature_citation": lit_citation,
        "steps": [
            {
                "e": st.e,
                "corrector_residual": st.corrector_residual,
                "stability_tag": st.stability_tag,
                "on_unit_circle": st.on_unit_circle,
            }
            for st in trace.steps
        ],
    }


def _process_seed(seed: Er3bpSeed) -> tuple[dict[str, object], Er3bpContinuationTrace]:
    trace = continue_and_monitor(seed, n_steps=_N_STEPS)
    verdict = adjudicate_trace(trace, seed, literature_check_search=offline_search)
    record = _record_for_seed(seed, trace, verdict.status, verdict.citation)
    return record, trace


def _outcome_breakdown(traces: Sequence[Er3bpContinuationTrace]) -> dict[str, int]:
    counts = {"survives": 0, "dies": 0, "bifurcates": 0}
    for t in traces:
        counts[t.outcome] = counts.get(t.outcome, 0) + 1
    return counts


def _generate_seeds(primary: str, secondary: str, target_e: float) -> list[Er3bpSeed]:
    """Generate L1-Lyapunov (mandatory) + DRO (best-effort) seeds for a system."""
    cr3bp = cr3bp_system(primary, secondary)
    system = ER3BPSystem.from_cr3bp(cr3bp, target_e)
    seeds: list[Er3bpSeed] = []

    st, period_f = lyapunov_seed(cr3bp, point="L1")
    seeds.append(
        Er3bpSeed(
            label=f"{secondary}-L1-lyapunov",
            system=system,
            state0=st,
            period_f=period_f,
            is_half_period_residual=True,
            target_e=target_e,
            source=_SEED_SOURCE,
        )
    )

    try:
        st_dro, period_dro = dro_seed(cr3bp)
    except Exception as exc:  # DRO may not converge for some μ; degrade gracefully
        _print_progress(
            f"  {secondary}: DRO seed did not converge ({type(exc).__name__}: {exc}); "
            "skipping DRO, continuing with Lyapunov only."
        )
    else:
        seeds.append(
            Er3bpSeed(
                label=f"{secondary}-dro",
                system=system,
                state0=st_dro,
                period_f=period_dro,
                is_half_period_residual=True,
                target_e=target_e,
                source=_SEED_SOURCE,
            )
        )

    return seeds


def main() -> None:
    t0 = time.time()
    _print_progress(f"#435 high-e Sun-planet ER3BP discovery (n_steps={_N_STEPS})")

    all_records: list[dict[str, object]] = []
    all_traces: list[Er3bpContinuationTrace] = []

    for primary, secondary, target_e in _SYSTEMS:
        _print_progress(f"=== {primary}-{secondary} (target e={target_e}) ===")
        seeds = _generate_seeds(primary, secondary, target_e)
        _print_progress(f"{primary}-{secondary}: {len(seeds)} seed(s) generated")

        sys_records: list[dict[str, object]] = []
        sys_traces: list[Er3bpContinuationTrace] = []
        for i, seed in enumerate(seeds, 1):
            _print_progress(
                f"{primary}-{secondary} [{i}/{len(seeds)}] continuing {seed.label} -> e={target_e}"
            )
            rec, trace = _process_seed(seed)
            sys_records.append(rec)
            sys_traces.append(trace)
            _print_progress(
                f"  {seed.label}: outcome={trace.outcome} e_max={trace.e_max_reached:.4f} "
                f"e_star={trace.e_star} lit={rec['literature_status']}"
            )

        counts = _outcome_breakdown(sys_traces)
        _print_progress(
            f"{primary}-{secondary} outcome breakdown: survives={counts['survives']} "
            f"dies={counts['dies']} bifurcates={counts['bifurcates']} (n={len(sys_traces)})"
        )
        all_records.extend(sys_records)
        all_traces.extend(sys_traces)

    # Write the per-seed JSONL.
    out_path = _DATA_DIR / "er3bp_discovery_435_highE.jsonl"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w") as f:
        for rec in all_records:
            f.write(json.dumps(rec) + "\n")
    _print_progress(f"Wrote {out_path.relative_to(_DATA_DIR.parent)} ({len(all_records)} records)")

    # Overall report.
    counts = _outcome_breakdown(all_traces)
    _print_progress(
        f"OVERALL outcome breakdown: survives={counts['survives']} "
        f"dies={counts['dies']} bifurcates={counts['bifurcates']} (n={len(all_traces)})"
    )
    for rec, t in zip(all_records, all_traces, strict=True):
        if t.outcome == "bifurcates":
            _print_progress(
                f"  BIFURCATION: {rec['label']} ({rec['primary']}-{rec['secondary']}) "
                f"e_star={t.e_star} target_e={t.target_e}"
            )
    lit = {"published": 0, "not-found": 0, "inconclusive": 0}
    for rec in all_records:
        status = str(rec["literature_status"])
        lit[status] = lit.get(status, 0) + 1
    _print_progress(
        f"OVERALL literature breakdown: published={lit['published']} "
        f"not-found={lit['not-found']} inconclusive={lit['inconclusive']}"
    )

    _print_progress(f"Campaign complete in {time.time() - t0:.1f}s")


if __name__ == "__main__":
    main()
