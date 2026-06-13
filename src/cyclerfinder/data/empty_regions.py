"""The first-class negative: the empty-region report artefact (Forge Phase 6).

The empty-set lesson (Hughes/Jones: a thorough sweep that finds nothing is a
publishable scientific result, not a failed run) captured structurally. Today the
empty-set findings (#110/#120/#122) live as prose in ``data/OUTSTANDING.md``;
this gives them a machine-readable, *bounded*, *reproducible*, *method-versioned*
report (design note ``docs/notes/2026-06-08-forge-phase6-discovery-design.md``
§6).

The mirror of the human-review queue: SILVER survivors go to
``data/review_queue.jsonl``; empty regions go to ``data/empty_regions.jsonl``.
Both are non-catalogue artefacts; neither auto-promotes anything
(:func:`is_catalogue_source` returns ``False`` by contract).

The bar for a negative to count (enforced by :func:`validate_empty_region`):

* ``search_extent["points_total"]`` present and > 0 — an unbounded negative is a
  silently-dropped negative,
* ``prune_gates`` non-empty — without them an empty set may be an over-pruning
  artefact,
* ``method_capability`` carries a non-empty tag set — an *unconditional* "empty"
  claim is forbidden; "empty" is never unconditional (design §6a).
"""

from __future__ import annotations

import json
from collections.abc import Iterator, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from cyclerfinder.data.method_capability import MethodCapability

DEFAULT_EMPTY_REGIONS_PATH: Path = (
    Path(__file__).resolve().parent.parent.parent.parent / "data" / "empty_regions.jsonl"
)


@dataclass(frozen=True)
class EmptyRegionReport:
    """One swept region that yielded zero promotable candidates (design §6).

    Carries enough to make the negative reproducible (the grid + git SHA),
    bounded (what was NOT covered), and method-versioned (the capability
    envelope, so a later more-capable method knows to re-sweep).
    """

    region_id: str
    family: str
    centre: str
    topologies: tuple[dict[str, Any], ...]
    method_capability: MethodCapability
    search_extent: dict[str, Any]
    prune_gates: tuple[str, ...]
    result: dict[str, Any]
    verdict: str
    interpretation: str
    source_anchors: str
    run: dict[str, Any]


def is_catalogue_source() -> bool:
    """The empty-region log is NON-catalogue (golden discipline).

    Returns ``False`` always: an empty-region record is a negative-result audit
    trail and never feeds the catalogue.
    """
    return False


def validate_empty_region(report: EmptyRegionReport) -> None:
    """Refuse a negative that is not first-class (design §6).

    Raises :class:`ValueError` if the report is unbounded (no positive
    ``points_total``), un-prunable-to-audit (empty ``prune_gates``), or makes an
    unconditional "empty" claim (empty ``method_capability.capability_tags``).
    """
    points_total = report.search_extent.get("points_total")
    if not points_total or points_total <= 0:
        raise ValueError(
            f"empty-region report {report.region_id!r} has no bounded search extent "
            f"(search_extent['points_total']={points_total!r}); an unbounded negative "
            f"is a silently-dropped negative."
        )
    if not report.prune_gates:
        raise ValueError(
            f"empty-region report {report.region_id!r} records no prune gates; cannot "
            f"distinguish a real empty set from an over-pruning artefact."
        )
    if not report.method_capability.capability_tags:
        raise ValueError(
            f"empty-region report {report.region_id!r} carries no method-capability "
            f"tags; 'empty' is never unconditional (design §6a)."
        )


def _to_payload(report: EmptyRegionReport) -> dict[str, Any]:
    """Serialise a report to a JSON-able dict (method_capability flattened)."""
    mc = report.method_capability
    return {
        "region_id": report.region_id,
        "family": report.family,
        "centre": report.centre,
        "topologies": list(report.topologies),
        "method_capability": {
            "genome": mc.genome,
            "corrector": mc.corrector,
            "capability_tags": sorted(mc.capability_tags),
            "git_sha": mc.git_sha,
        },
        "search_extent": report.search_extent,
        "prune_gates": list(report.prune_gates),
        "result": report.result,
        "verdict": report.verdict,
        "interpretation": report.interpretation,
        "source_anchors": report.source_anchors,
        "run": report.run,
    }


def _from_payload(payload: dict[str, Any]) -> EmptyRegionReport:
    """Inverse of :func:`_to_payload` — rebuild a frozen report from JSON."""
    mc_raw = payload["method_capability"]
    method = MethodCapability(
        genome=mc_raw["genome"],
        corrector=mc_raw["corrector"],
        capability_tags=frozenset(mc_raw["capability_tags"]),
        git_sha=mc_raw["git_sha"],
    )
    return EmptyRegionReport(
        region_id=payload["region_id"],
        family=payload["family"],
        centre=payload["centre"],
        # topologies / prune_gates / source_anchors are optional metadata absent
        # from the earliest (#219) entries; default to empty so the whole file
        # round-trips. Pinned by test_load_real_empty_regions_file.
        topologies=tuple(payload.get("topologies", ())),
        method_capability=method,
        search_extent=payload["search_extent"],
        prune_gates=tuple(payload.get("prune_gates", ())),
        result=payload["result"],
        verdict=payload["verdict"],
        interpretation=payload["interpretation"],
        source_anchors=payload.get("source_anchors", ""),
        run=payload["run"],
    )


def append_empty_region(
    path: Path | str,
    report: EmptyRegionReport,
    *,
    validate: bool = True,
) -> None:
    """Append one validated report to the empty-region log (creating parent dirs)."""
    if validate:
        validate_empty_region(report)
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    line = json.dumps(_to_payload(report), ensure_ascii=True)
    with p.open("a", encoding="utf-8") as fh:
        fh.write(line + "\n")


def load_empty_regions(path: Path | str) -> Iterator[EmptyRegionReport]:
    """Yield every :class:`EmptyRegionReport` from the log (read-only)."""
    p = Path(path)
    if not p.exists():
        return
    for raw in p.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line:
            continue
        yield _from_payload(json.loads(line))


def load_empty_regions_list(path: Path | str) -> Sequence[EmptyRegionReport]:
    """Eager :func:`load_empty_regions` (a registry for :func:`should_sweep`)."""
    return list(load_empty_regions(path))


__all__ = [
    "DEFAULT_EMPTY_REGIONS_PATH",
    "EmptyRegionReport",
    "append_empty_region",
    "is_catalogue_source",
    "load_empty_regions",
    "load_empty_regions_list",
    "validate_empty_region",
]
