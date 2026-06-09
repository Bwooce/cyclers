"""Forge Phase 6 orchestrator — Jovian VILM moon-system novelty sweep.

The first-campaign run: the planet-centric Galilean Jovian VILM-gated space
(design note ``docs/notes/2026-06-08-forge-phase6-discovery-design.md`` §5; plan
Task 4.2). Mirrors :mod:`scripts.forge_novelty_run` (fan-out -> per-finding
adversarial panel -> human queue), swapping the topology set for
:func:`cyclerfinder.data.discover_novel.jovian_galilean_topologies` +
:func:`cyclerfinder.data.discover_novel.discover_novel_moon`, and adding two
Phase-6 first-class pieces:

* **The re-sweep gate (design §6b).** BEFORE sweeping, query
  ``data/empty_regions.jsonl`` via
  :func:`cyclerfinder.data.method_capability.should_sweep`: SKIP the region iff a
  prior >=-capable method already emptied it. The proposed method here is the
  single-ellipse no-leveraging :class:`MethodCapability`.
* **The empty-region emit (design §6).** If the sweep yields zero promotable
  candidates (the EXPECTED EMPTY outcome — the #76 honest-risk generalises), emit
  a bounded + method-versioned :class:`EmptyRegionReport` carrying the actual
  search extent + best-achieved V_inf + the V_inf-floor gap, so a later
  multi-arc/n-body/low-thrust method re-sweeps per §6b.

HONEST-RISK (binding): a bend-INFEASIBLE closure is auto-REJECTED inside
:func:`evaluate_closure` (never SILVER); a barren sweep is a SUCCESS, recorded as
a rigorous bounded negative. We do NOT loosen tol/budget/bend-cap to manufacture a
survivor.

Run::

    uv run python scripts/forge_phase6_moon_run.py --epochs 64 --workers 16 \
        --empty-regions data/empty_regions.jsonl --queue data/review_queue.jsonl \
        --report /tmp/forge_phase6_jovian.txt
"""

from __future__ import annotations

import argparse
import subprocess
import sys
import time
from datetime import UTC, datetime

from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.core.satellites import PRIMARIES
from cyclerfinder.data.discover_novel import (
    NoveltyFinding,
    discover_novel_moon,
    jovian_galilean_permutation_topologies,
    jovian_galilean_topologies,
    saturnian_titan_tour_topologies,
)
from cyclerfinder.data.empty_regions import (
    EmptyRegionReport,
    append_empty_region,
    load_empty_regions_list,
)
from cyclerfinder.data.method_capability import MethodCapability, should_sweep
from cyclerfinder.data.review_queue import (
    ReviewQueueEntry,
    append_review_entry,
    validate_review_entry,
)
from cyclerfinder.verify.adversarial import PanelResult, adversarial_panel
from cyclerfinder.verify.gauntlet import VerdictTier

# The single-ellipse no-leveraging method-capability descriptor for THIS campaign
# (design §6a). The Jovian centred scan is a single ballistic arc per leg,
# patched-conic, coplanar-circular about the primary — no DSM, no multi-arc, no
# n-body. A later more-capable method (multi-arc / n-body / low-thrust) will read
# this and re-sweep per the §6b gate.
_CAMPAIGN_TAGS = frozenset({"ballistic", "patched-conic", "single-arc", "coplanar"})

# Topology-set registry (Forge Phase 6 #178). Each entry maps a --topology-set
# choice to (enumerator, region-id prefix). The dated region id is built per-run.
# The default keeps the original first-campaign Jovian-Galilean behaviour identical.
_TOPOLOGY_SETS: dict[str, tuple[object, str]] = {
    "jovian-galilean": (jovian_galilean_topologies, "jovian-IEG-vilm"),
    "jovian-permutations": (jovian_galilean_permutation_topologies, "jovian-perm-vilm"),
    "saturnian-titan": (saturnian_titan_tour_topologies, "saturnian-titan-vilm"),
}


def _git_sha() -> str:
    try:
        return subprocess.check_output(["git", "rev-parse", "--short", "HEAD"], text=True).strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        return "unknown"


def _campaign_method(center: str) -> MethodCapability:
    return MethodCapability(
        genome=f"single-ellipse free-return ({center}-centric, no-leveraging)",
        corrector=f"ballistic_correct (mu_central=PRIMARIES[{center}])",
        capability_tags=_CAMPAIGN_TAGS,
        git_sha=_git_sha(),
    )


def _panel_for_finding(
    finding: NoveltyFinding,
    *,
    ephem: Ephemeris,
    n_verifiers: int,
    vinf_cap: float,
    mu_central: float,
) -> PanelResult:
    """Adversarial panel for one Jovian finding (centred ephemeris + mu_central)."""
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
        mu_central=mu_central,
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
        mu_central=mu_central,
    )


def _finding_to_queue_entry(finding: NoveltyFinding, panel: PanelResult) -> ReviewQueueEntry:
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
        model_assumption="circular-coplanar",
        verdict_audit=finding.verdict.axis_results,
        panel=panel.as_dict(),
        t_added=datetime.now(UTC).isoformat(),
        literature_check=None,  # machine default: not yet checked (design §3c)
    )


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description="Forge Phase 6 VILM moon-system sweep")
    p.add_argument(
        "--topology-set",
        default="jovian-galilean",
        choices=sorted(_TOPOLOGY_SETS),
        help="which Phase 6 topology set to sweep (default: jovian-galilean, first campaign)",
    )
    p.add_argument("--center", default="Jupiter")
    p.add_argument("--region-id", default=None, help="empty-region id (defaults to dated)")
    p.add_argument("--epochs", type=int, default=64, help="epoch grid points per topology")
    p.add_argument("--span-days", type=float, default=8.0)
    p.add_argument("--budget-kms", type=float, default=50.0, help="VILM prune budget")
    p.add_argument("--vinf-cap", type=float, default=14.0)
    p.add_argument(
        "--vinf-floor-kms",
        type=float,
        default=6.0,
        help="target V_inf floor (the bend-feasible band; for the gap report)",
    )
    p.add_argument("--workers", type=int, default=16)
    p.add_argument("--n-verifiers", type=int, default=3)
    p.add_argument("--empty-regions", default="data/empty_regions.jsonl")
    p.add_argument("--queue", default="data/review_queue.jsonl")
    p.add_argument("--report", default=None, help="optional path to mirror stdout summary")
    p.add_argument("--no-write", action="store_true", help="dry run: write neither artefact")
    args = p.parse_args(argv)

    center = args.center
    enumerator, region_prefix = _TOPOLOGY_SETS[args.topology_set]
    region_id = args.region_id or f"{region_prefix}-{datetime.now(UTC).date().isoformat()}"
    method = _campaign_method(center)
    mu_central = PRIMARIES[center]

    lines: list[str] = []

    def emit(s: str) -> None:
        print(s, flush=True)
        lines.append(s)

    emit("=== Forge Phase 6 — VILM moon-system novelty sweep ===")
    emit(f"topology_set={args.topology_set} center={center} region_id={region_id}")
    emit(f"epochs={args.epochs} budget_kms={args.budget_kms}")
    emit(f"method={method.genome} tags={sorted(method.capability_tags)} git_sha={method.git_sha}")

    # --- The re-sweep gate (design §6b) -----------------------------------
    registry = load_empty_regions_list(args.empty_regions)
    if not should_sweep(region_id=region_id, method=method, registry=registry):
        emit(
            f"-> SKIP: a prior >=-capable method already emptied region {region_id!r} "
            f"(should_sweep=False). Nothing new to learn; not re-sweeping."
        )
        return 0

    ephem = Ephemeris(model="circular", center=center)
    topologies = enumerator()  # type: ignore[operator]
    base_t0 = 0.6 * 86400.0  # the #76 converging seed window (NON-GOLDEN)
    t_start = time.monotonic()

    n_silver = n_bronze = n_rejected = n_novel = n_closed = 0
    n_queued = n_panel_killed = 0
    best_max_vinf = 0.0

    for finding in discover_novel_moon(
        base_t0_sec=base_t0,
        topologies=topologies,
        center=center,
        budget_kms=args.budget_kms,
        n_epochs=args.epochs,
        span_days=args.span_days,
        vinf_cap=args.vinf_cap,
        max_workers=args.workers,
        distinct_only=True,
    ):
        n_closed += 1
        tier = finding.verdict.tier
        if tier is VerdictTier.SILVER:
            n_silver += 1
        elif tier is VerdictTier.BRONZE:
            n_bronze += 1
        elif tier is VerdictTier.REJECTED:
            n_rejected += 1
        if finding.match_outcome == "novel":
            n_novel += 1
        best_max_vinf = max(best_max_vinf, finding.max_vinf_kms)

        emit(
            f"[{tier.value.upper()}] {finding.candidate_id} match={finding.match_outcome} "
            f"max_vinf={finding.max_vinf_kms:.3f} bend_feasible={finding.bend_feasible}"
        )

        if tier in (VerdictTier.SILVER, VerdictTier.GOLD):
            panel = _panel_for_finding(
                finding,
                ephem=ephem,
                n_verifiers=args.n_verifiers,
                vinf_cap=args.vinf_cap,
                mu_central=mu_central,
            )
            emit(
                f"    PANEL n_refuted={panel.n_refuted}/{panel.n_verifiers} "
                f"majority_refute={panel.majority_refute}"
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
            if not args.no_write:
                append_review_entry(args.queue, entry, validate=False)
            n_queued += 1
            emit(f"    -> QUEUED for human review (literature_check=None) ({args.queue})")

    wall_s = time.monotonic() - t_start
    # search_extent: the ACTUAL grid covered (epochs x topologies x branch points).
    points_total = args.epochs * len(topologies)
    gap_kms = max(0.0, best_max_vinf - args.vinf_floor_kms) if n_closed else float("nan")

    emit("")
    emit("=== SUMMARY ===")
    emit(
        f"closed={n_closed} novel={n_novel} SILVER={n_silver} "
        f"BRONZE={n_bronze} REJECTED={n_rejected}"
    )
    emit(f"panel: queued={n_queued} killed={n_panel_killed}")
    emit(
        f"best_max_vinf={best_max_vinf:.3f} km/s "
        f"vinf_floor_target={args.vinf_floor_kms} gap={gap_kms}"
    )
    emit(f"wall_s={wall_s:.1f}")

    # --- The empty-region emit (design §6) --------------------------------
    if n_queued == 0:
        report = EmptyRegionReport(
            region_id=region_id,
            family=f"planet-centric moon system ({center})",
            centre=center,
            topologies=tuple(
                {
                    "sequence": list(s.sequence),
                    "per_leg_revs": list(s.per_leg_revs),
                    "period_k": s.period_k,
                }
                for s in topologies
            ),
            method_capability=method,
            search_extent={
                "n_epochs": args.epochs,
                "span_days": args.span_days,
                "n_topologies": len(topologies),
                "points_total": points_total,
                "ephem_model": "circular",
                "center": center,
            },
            prune_gates=(
                "vilm_dv_floor<=budget",
                "linkable(Jovicentric)",
                "max_bend_deg feasibility",
            ),
            result={
                "closed": n_closed,
                "novel": n_novel,
                "silver_survivors": n_silver,
                "bend_feasible_queued": n_queued,
                "best_max_vinf_kms": best_max_vinf,
                "vinf_floor_target_kms": args.vinf_floor_kms,
                "gap_kms": gap_kms,
                "rejected": n_rejected,
                "panel_killed": n_panel_killed,
            },
            verdict="EMPTY — no bend-feasible ballistic closure promoted below the V_inf floor",
            interpretation=(
                "no-leveraging corrector closes a higher-V_inf family; VILM gating did not "
                "surface a bend-feasible Laplace-resonant tour promotable past SILVER"
            ),
            source_anchors=(
                "none populated in the (circular-coplanar, Jupiter) bucket; the Jovian "
                "catalogue rows are family-seed null-numeric"
            ),
            run={
                "date": datetime.now(UTC).date().isoformat(),
                "cores": args.workers,
                "git_sha": method.git_sha,
                "wall_s": round(wall_s, 1),
            },
        )
        if not args.no_write:
            append_empty_region(args.empty_regions, report)
        emit(f"-> EMPTY-REGION report emitted ({args.empty_regions}) region_id={region_id}")
    else:
        emit(f"-> {n_queued} SILVER survivor(s) queued; NO empty-region report (region not barren)")

    if args.report:
        from pathlib import Path

        Path(args.report).write_text("\n".join(lines) + "\n", encoding="utf-8")
        print(f"(report mirrored to {args.report})", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
