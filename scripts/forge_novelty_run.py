"""Forge Phase 5 orchestrator — adversarial novelty run + human gate.

WORKFLOW-TOOL CHOICE (plan §33-34). The plan asks for "a Workflow-tool script
(.claude/workflows/forge-novelty.md or .js) per the Workflow tool conventions IF
determinable from the repo, OTHERWISE a Python orchestrator with the same
fan-out + per-finding adversarial-panel structure." No ``.claude/workflows/``
convention exists in this repo (verified 2026-06-06: the directory is absent and
no workflow artifact is referenced anywhere in src/ scripts/ data/), so this is
the Python orchestrator alternative the plan sanctions.

Structure:

* **Fan-out** — :func:`cyclerfinder.data.discover_novel.discover_novel` drives
  the E-M multi-arc space (topology x epoch scan grid) in parallel
  (16-core ``scan_parallel``). This is the empirically-demonstrated bend-feasible
  frontier; VEM single-ellipse novelty is empirically nil (#110) and is NOT
  scanned.
* **Per-finding adversarial panel** —
  :func:`cyclerfinder.verify.adversarial.adversarial_panel` runs N independent
  verifiers (falsification probe, independent re-closure, perturbed-seed
  robustness); a majority-refute kills the candidate.
* **Human gate** — surviving SILVER candidates are appended to
  ``data/review_queue.jsonl`` (NON-catalogue; no auto-promotion).

Run::

    uv run python scripts/forge_novelty_run.py --epochs 16 --workers 16 \
        --queue data/review_queue.jsonl --report /tmp/forge_novelty.txt

This is a *slow* run (real DE440 corrector x epoch grid x adversarial re-solves).
Capture stdout to a file; the per-finding lines are the verbatim result.
"""

from __future__ import annotations

import argparse
import sys
from datetime import UTC, datetime

from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.data.discover_novel import (
    NoveltyFinding,
    discover_novel,
    em_multiarc_topologies,
)
from cyclerfinder.data.review_queue import (
    ReviewQueueEntry,
    append_review_entry,
    validate_review_entry,
)
from cyclerfinder.verify.adversarial import PanelResult, adversarial_panel
from cyclerfinder.verify.gauntlet import VerdictTier

_J2000 = datetime(2000, 1, 1, 12, 0, 0, tzinfo=UTC)


def _base_t0_sec(date_iso: str) -> float:
    dt = datetime.fromisoformat(date_iso)
    if dt.tzinfo is None:
        dt = dt.replace(tzinfo=UTC)
    return (dt - _J2000).total_seconds()


def _panel_for_finding(
    finding: NoveltyFinding,
    *,
    ephem: Ephemeris,
    n_verifiers: int,
    vinf_cap: float,
) -> PanelResult:
    """Run the adversarial panel for one finding by reconstructing its closure.

    The finding carries its exact topology + seed (sequence / per_leg_revs /
    per_leg_branch / t0 / slack_leg), so the panel re-runs the corrector on the
    identical structure — no lossy id-parsing.
    """
    from cyclerfinder.search.correct import ballistic_correct

    free_tofs = [t for i, t in enumerate(finding.tof_days) if i != finding.slack_leg]
    period_sec = sum(finding.tof_days) * 86400.0
    closure = ballistic_correct(
        sequence=finding.sequence,
        per_leg_revs=finding.per_leg_revs,
        per_leg_branch=finding.per_leg_branch,
        t0_seed_sec=finding.t0_sec,
        tof_seed_days=free_tofs,
        period_sec=period_sec,
        ephem=ephem,
        vinf_cap=vinf_cap,
        slack_leg=finding.slack_leg,
    )
    return adversarial_panel(
        closure,
        sequence=finding.sequence,
        per_leg_revs=finding.per_leg_revs,
        per_leg_branch=finding.per_leg_branch,
        period_sec=period_sec,
        slack_leg=finding.slack_leg,
        vinf_cap=vinf_cap,
        ephem=ephem,
        n_verifiers=n_verifiers,
    )


def _finding_to_queue_entry(
    finding: NoveltyFinding,
    panel: PanelResult,
) -> ReviewQueueEntry:
    assert finding.signature is not None
    return ReviewQueueEntry(
        candidate_id=finding.candidate_id,
        signature_hash=finding.signature.hash,
        verdict_tier=finding.verdict.tier.value,
        match_outcome=finding.match_outcome,
        known_id=finding.known_id,
        superseded_by=finding.superseded_by,
        vinf_per_encounter_kms=finding.vinf_per_encounter_kms,
        tof_days=finding.tof_days,
        bend_feasible=finding.bend_feasible,
        max_vinf_kms=finding.max_vinf_kms,
        sequence=finding.sequence,
        period_k=finding.period_k,
        model_assumption="analytic-ephemeris",
        verdict_audit=finding.verdict.axis_results,
        panel=panel.as_dict(),
        t_added=datetime.now(UTC).isoformat(),
    )


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description="Forge Phase 5 adversarial novelty run")
    p.add_argument("--date", default="2030-03-02", help="centre launch date (ISO)")
    p.add_argument("--epochs", type=int, default=16, help="epoch grid points per topology")
    p.add_argument("--span-days", type=float, default=1280.0, help="launch-epoch span")
    p.add_argument("--vinf-cap", type=float, default=14.0)
    p.add_argument("--workers", type=int, default=16)
    p.add_argument("--n-verifiers", type=int, default=3)
    p.add_argument("--queue", default="data/review_queue.jsonl")
    p.add_argument("--report", default=None, help="optional path to mirror stdout summary")
    p.add_argument(
        "--no-queue-write", action="store_true", help="dry run: do not write the review queue"
    )
    args = p.parse_args(argv)

    ephem = Ephemeris("astropy")
    base_t0 = _base_t0_sec(args.date)
    topologies = em_multiarc_topologies()

    lines: list[str] = []

    def emit(s: str) -> None:
        print(s, flush=True)
        lines.append(s)

    emit("=== Forge Phase 5 — adversarial novelty run ===")
    emit(f"date={args.date} epochs={args.epochs} span_days={args.span_days} ")
    emit(
        f"topologies={len(topologies)} (E-M multi-arc) "
        f"vinf_cap={args.vinf_cap} workers={args.workers}"
    )
    emit("")

    n_silver = n_bronze = n_rejected = n_known = n_novel = n_supersede = 0
    n_queued = n_panel_killed = 0

    for finding in discover_novel(
        ephem=ephem,
        topologies=topologies,
        base_t0_sec=base_t0,
        n_epochs=args.epochs,
        span_days=args.span_days,
        vinf_cap=args.vinf_cap,
        max_workers=args.workers,
        distinct_only=True,
    ):
        tier = finding.verdict.tier
        if tier is VerdictTier.SILVER:
            n_silver += 1
        elif tier is VerdictTier.BRONZE:
            n_bronze += 1
        elif tier is VerdictTier.REJECTED:
            n_rejected += 1
        if finding.match_outcome == "known":
            n_known += 1
        elif finding.match_outcome == "novel":
            n_novel += 1
        if finding.superseded_by:
            n_supersede += 1

        sig = finding.signature
        emit(
            f"[{tier.value.upper()}] {finding.candidate_id} "
            f"match={finding.match_outcome} known_id={finding.known_id} "
            f"superseded_by={list(finding.superseded_by)}"
        )
        emit(
            f"    vinf_per_enc={[round(v, 3) for v in finding.vinf_per_encounter_kms]} "
            f"max_vinf={finding.max_vinf_kms:.3f} bend_feasible={finding.bend_feasible} "
            f"tof_days={[round(t, 1) for t in finding.tof_days]}"
        )
        if sig is not None:
            emit(
                f"    sig={sig.hash} seq={sig.sequence_canonical} k={sig.period_k} "
                f"vinf_multiset={[(b, round(v, 2)) for b, v in sig.vinf_multiset_binned]}"
            )

        # Per-finding adversarial panel (only for non-rejected machine-confirmed
        # candidates worth promoting; a BRONZE/SILVER both get panelled).
        if tier in (VerdictTier.SILVER, VerdictTier.GOLD):
            panel = _panel_for_finding(
                finding, ephem=ephem, n_verifiers=args.n_verifiers, vinf_cap=args.vinf_cap
            )
            emit(
                f"    PANEL n_verifiers={panel.n_verifiers} n_refuted={panel.n_refuted} "
                f"majority_refute={panel.majority_refute} verdicts={list(panel.verifier_verdicts)}"
            )
            if panel.majority_refute:
                n_panel_killed += 1
                emit("    -> KILLED by adversarial panel (majority refute)")
                continue
            entry = _finding_to_queue_entry(finding, panel)
            try:
                validate_review_entry(entry)
            except ValueError as exc:
                emit(f"    -> NOT queued: {exc}")
                continue
            if not args.no_queue_write:
                append_review_entry(args.queue, entry, validate=False)
            n_queued += 1
            emit(f"    -> QUEUED for human review ({args.queue})")
        emit("")

    emit("")
    emit("=== SUMMARY ===")
    emit(f"SILVER={n_silver} BRONZE={n_bronze} REJECTED={n_rejected}")
    emit(f"match: known={n_known} novel={n_novel} superseded={n_supersede}")
    emit(f"panel: queued={n_queued} killed={n_panel_killed}")

    if args.report:
        from pathlib import Path

        Path(args.report).write_text("\n".join(lines) + "\n", encoding="utf-8")
        print(f"(report mirrored to {args.report})", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
