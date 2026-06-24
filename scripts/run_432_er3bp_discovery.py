"""#432 — ER3BP discovery campaign runner (Phase A / Phase B).

Continues rotating-frame CR3BP cycler families into the elliptic restricted
three-body problem (e>0), Floquet-monitoring each continuation step for survival,
death, or an elliptic<->hyperbolic bifurcation, and adjudicates each survivor's
structural signature against the published literature (offline corpus).

Phase A — the guaranteed seed floor:
  * the sourced Broucke-1969 Earth-Moon family (``standard_family_seeds``), and
  * any CR3BP-class catalogue seed for which a rotating-frame IC exists
    (``catalogue_cr3bp_seeds``; currently the Koblick NRHO table, Earth-Moon),
  each continued to the REAL Earth-Moon eccentricity e=0.0549.

Phase B — the standard floor family re-run at higher e to probe where the
Earth-Moon Broucke family bifurcates or dies. High-e Sun-planet systems
(Sun-Mercury e=0.206, Sun-Mars e=0.093) are SEED-LIMITED: the Broucke IC is
Earth-Moon-mu-specific, so a Sun-Mercury continuation would need a Sun-Mercury
CR3BP seed IC we do NOT have. Phase B therefore attempts only systems with an
encoded seed (Earth-Moon at synthetic high-e targets) and LOGS the seed-limited
Sun-planet systems explicitly rather than fabricating ICs.

Report-only — NO catalogue writeback. Per the
:mod:`cyclerfinder.search.literature_check` discipline preamble a "not-found"
literature status is NECESSARY-NOT-SUFFICIENT for novelty; the V0-V5 gauntlet
still governs.

Usage::

    uv run python scripts/run_432_er3bp_discovery.py
"""

from __future__ import annotations

import json
import time
from collections.abc import Sequence
from pathlib import Path

from cyclerfinder.search.er3bp_discovery import (
    Er3bpContinuationTrace,
    Er3bpSeed,
    adjudicate_trace,
    catalogue_cr3bp_seeds,
    continue_and_monitor,
    standard_family_seeds,
)
from cyclerfinder.search.literature_check import KNOWN_CORPUS, SearchResult

_DATA_DIR = Path(__file__).resolve().parents[1] / "data"

# Real Earth-Moon orbital eccentricity — the Phase-A continuation target.
_EM_REAL_E = 0.0549

# Fine stepping: the secant predictor dies on step size below ~n_steps=30
# (Δe<=0.0018). 60 keeps the predictor well inside its convergence basin.
_N_STEPS = 60

# Phase B synthetic high-e probes for the Earth-Moon Broucke floor family — push
# past the real e=0.0549 to find where the EM family bifurcates / dies.
_PHASE_B_EM_TARGETS = (0.10, 0.15)

# Phase B Sun-planet systems we WOULD probe at their real high e but CANNOT: the
# Broucke IC is Earth-Moon-mu-specific, so each needs a Sun-planet CR3BP seed IC
# we do not have. Logged as seed-limited; NO ICs fabricated.
_PHASE_B_SEED_LIMITED = (
    ("Sun", "Mercury", 0.206),
    ("Sun", "Mars", 0.093),
)


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


def _report_phase(
    phase: str,
    records: Sequence[dict[str, object]],
    traces: Sequence[Er3bpContinuationTrace],
) -> None:
    counts = _outcome_breakdown(traces)
    _print_progress(
        f"{phase} outcome breakdown: survives={counts['survives']} "
        f"dies={counts['dies']} bifurcates={counts['bifurcates']} "
        f"(n={len(traces)})"
    )
    for rec, t in zip(records, traces, strict=True):
        if t.outcome == "bifurcates":
            _print_progress(
                f"  {phase} BIFURCATION: {rec['label']} e_star={t.e_star} target_e={t.target_e}"
            )
    lit = {"published": 0, "not-found": 0, "inconclusive": 0}
    for rec in records:
        status = str(rec["literature_status"])
        lit[status] = lit.get(status, 0) + 1
    _print_progress(
        f"{phase} literature breakdown: published={lit['published']} "
        f"not-found={lit['not-found']} inconclusive={lit['inconclusive']}"
    )


def _write_jsonl(path: Path, records: Sequence[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as f:
        for rec in records:
            f.write(json.dumps(rec) + "\n")
    _print_progress(f"Wrote {path.relative_to(_DATA_DIR.parent)} ({len(records)} records)")


def _run_phase_a() -> None:
    _print_progress("=== Phase A: seed floor continued to real Earth-Moon e ===")
    seeds: list[Er3bpSeed] = list(standard_family_seeds(target_e=_EM_REAL_E))
    cat_seeds, n_skipped = catalogue_cr3bp_seeds(target_e=_EM_REAL_E)
    seeds.extend(cat_seeds)
    _print_progress(
        f"Phase A seeds: {len(seeds)} "
        f"(standard floor + {len(cat_seeds)} catalogue, {n_skipped} skipped no-IC)"
    )

    records: list[dict[str, object]] = []
    traces: list[Er3bpContinuationTrace] = []
    for i, seed in enumerate(seeds, 1):
        _print_progress(f"Phase A [{i}/{len(seeds)}] continuing {seed.label} -> e={_EM_REAL_E}")
        rec, trace = _process_seed(seed)
        records.append(rec)
        traces.append(trace)
        _print_progress(
            f"  {seed.label}: outcome={trace.outcome} e_max={trace.e_max_reached:.4f} "
            f"e_star={trace.e_star} lit={rec['literature_status']}"
        )

    _write_jsonl(_DATA_DIR / "er3bp_discovery_phaseA.jsonl", records)
    _report_phase("Phase A", records, traces)


def _run_phase_b() -> None:
    _print_progress("=== Phase B: floor family at high-e (seed-limited) ===")
    for primary, secondary, real_e in _PHASE_B_SEED_LIMITED:
        _print_progress(
            f"Phase B SEED-LIMITED: {primary}-{secondary} (real e={real_e}) — no "
            f"{primary}-{secondary} CR3BP seed IC available; the Broucke IC is "
            "Earth-Moon-mu-specific. Skipping; ICs NOT fabricated."
        )

    records: list[dict[str, object]] = []
    traces: list[Er3bpContinuationTrace] = []
    for target_e in _PHASE_B_EM_TARGETS:
        floor = standard_family_seeds(target_e=target_e)
        for seed in floor:
            label = f"{seed.label}-highE{target_e:g}"
            probe = Er3bpSeed(
                label=label,
                system=seed.system,
                state0=seed.state0,
                period_f=seed.period_f,
                is_half_period_residual=seed.is_half_period_residual,
                target_e=target_e,
                source=seed.source + f" (Phase B high-e probe target_e={target_e})",
            )
            _print_progress(f"Phase B continuing {label} -> e={target_e}")
            rec, trace = _process_seed(probe)
            records.append(rec)
            traces.append(trace)
            _print_progress(
                f"  {label}: outcome={trace.outcome} e_max={trace.e_max_reached:.4f} "
                f"e_star={trace.e_star} lit={rec['literature_status']}"
            )

    _write_jsonl(_DATA_DIR / "er3bp_discovery_phaseB.jsonl", records)
    _report_phase("Phase B", records, traces)


def main() -> None:
    t0 = time.time()
    _print_progress(f"#432 ER3BP discovery campaign (n_steps={_N_STEPS})")
    _run_phase_a()
    _run_phase_b()
    _print_progress(f"Campaign complete in {time.time() - t0:.1f}s")


if __name__ == "__main__":
    main()
