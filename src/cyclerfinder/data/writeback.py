"""Catalogue-record writeback helpers (M7) — spec §16.1 / §16.4 / §16.5.

Pure ``entry -> entry'`` functions that merge a validation-gate result
or a discovery/rediscovery record into a :class:`CatalogueEntry`,
returning a fresh frozen instance (immutability preserved per the M3 /
M6a pattern). The corresponding ``raw`` YAML projection is updated in
lock-step so :func:`serialise_entry_yaml` can round-trip the entry
without re-deriving the YAML shape.

Spec / plan references
----------------------
* Spec §16.1 — ``validation.gates`` schema (V0-V5 per-gate blocks).
* Spec §16.4 — rediscovery audit trail (attribution never overwritten).
* Spec §16.5 — ``source: this-project`` discovery registration.
* Plan: ``docs/phases/m7-catalogue-novelty-matching/plan.md`` §3.3.

CI safety
---------
None of these helpers touch ``data/catalogue.yaml``. They operate on
in-memory :class:`CatalogueEntry` instances; on-disk writeback is an
operator-driven step (plan §5 risk #6).
"""

from __future__ import annotations

from typing import Any

import yaml  # type: ignore[import-untyped]

from cyclerfinder.data.catalog import CatalogueEntry, _replace_entry
from cyclerfinder.data.validate import has_level_evidence
from cyclerfinder.verify.crosscheck import LambertCrosscheckResult
from cyclerfinder.verify.propagate import StabilityReport
from cyclerfinder.verify.real_closure import RealClosureResult

_LEVELS: tuple[str, ...] = ("V0", "V1", "V2", "V3", "V4", "V5")


def _assert_level_evidence(cycler_id: str, level: str) -> None:
    """Refuse to persist a level above V0 without registered evidence.

    The C1 guard (task #196): the apply_* helpers write the NESTED
    ``validation.level`` field, which historically routed around the
    over-claim guard (``validate_validation_level`` read only the
    top-level ``validation_level`` tag). Both ends are now closed — the
    validator also reads the nested field, and this assertion stops an
    unregistered promotion at the writeback source. ``V0`` (the
    internal-consistency floor) never needs evidence.
    """
    if not has_level_evidence(cycler_id, level):
        raise ValueError(
            f"{cycler_id}: refusing to write validation level {level!r} above V0 — "
            f"({cycler_id!r}, {level!r}) has no recorded mechanical evidence in "
            f"_LEVEL_EVIDENCE (cyclerfinder.data.validate). Register the evidence "
            f"pointer first; when in doubt, V0."
        )


def _gate_passed(name: str, gate: object) -> bool:
    """Whether a gate block counts as passed.

    A gate passes when its block carries an explicit ``pass: True``.
    ``V0`` is treated as implicitly passing when its block is present
    but lacks an explicit flag (any validated entry has satisfied the
    internal-consistency floor).
    """
    if not isinstance(gate, dict):
        return False
    if "pass" in gate:
        return bool(gate["pass"])
    return name == "V0"


def _highest_passing_level(gates: dict[str, Any]) -> str:
    """Return the largest ``V0..V5`` gate whose block passes.

    Scans cheapest-first and keeps the highest passing level. Defaults
    to ``"V0"`` (the implicit internal-consistency floor) when no gate
    passes — any entry reaching writeback has at least been
    constructed.
    """
    best = "V0"
    for name in _LEVELS:
        gate = gates.get(name)
        if gate is not None and _gate_passed(name, gate):
            best = name
    return best


def _with_validation(entry: CatalogueEntry, new_validation: dict[str, Any]) -> CatalogueEntry:
    """Replace ``entry.validation`` and mirror it into ``entry.raw``."""
    new_raw = {**entry.raw, "validation": new_validation}
    return _replace_entry(entry, validation=new_validation, raw=new_raw)


def apply_v0_v1_to_entry(
    entry: CatalogueEntry,
    v0_result: bool,
    v1_result: tuple[LambertCrosscheckResult, ...],
) -> CatalogueEntry:
    """Merge V0 (internal consistency) + V1 (Lambert cross-check) gates.

    ``v0_result`` is the internal-consistency pass flag (e.g.
    ``OptimisationResult.constraints_satisfied``). ``v1_result`` is the
    tuple returned by
    :func:`cyclerfinder.verify.crosscheck.crosscheck_cycler`; V1 passes
    only when it is non-empty AND every leg passes (an all-multi-rev
    cycler returns an empty tuple, which is "V1 not applicable", not a
    pass).
    """
    new_validation = dict(entry.validation)
    gates = dict(new_validation.get("gates", {}))
    gates["V0"] = {"pass": bool(v0_result)}
    if v1_result:
        max_diff_mps = max(r.max_diff_mps for r in v1_result)
        v1_pass = all(r.passed for r in v1_result)
    else:
        max_diff_mps = 0.0
        v1_pass = False
    gates["V1"] = {"pass": v1_pass, "max_diff_mps": max_diff_mps}
    new_validation["gates"] = gates
    level = _highest_passing_level(gates)
    _assert_level_evidence(entry.id, level)
    new_validation["level"] = level
    return _with_validation(entry, new_validation)


def apply_v2_to_entry(entry: CatalogueEntry, report: StabilityReport) -> CatalogueEntry:
    """Merge the M6a :class:`StabilityReport` into ``validation.gates.V2``.

    Spec §16.1 carries ``V2: {max_drift_km}``; M7 adds ``pass``,
    ``n_laps``, ``per_lap_drift_km``, and ``frame_used`` as additive
    diagnostic fields (consumers ignore unknown fields per the v2
    additive-fields convention).
    """
    new_validation = dict(entry.validation)
    gates = dict(new_validation.get("gates", {}))
    gates["V2"] = {
        "pass": report.stable,
        "max_drift_km": report.max_drift_km,
        "n_laps": report.n_laps_propagated,
        "per_lap_drift_km": list(report.per_lap_drift_km),
        "frame_used": report.frame_used,
    }
    new_validation["gates"] = gates
    level = _highest_passing_level(gates)
    _assert_level_evidence(entry.id, level)
    new_validation["level"] = level
    return _with_validation(entry, new_validation)


def apply_v3_to_entry(entry: CatalogueEntry, report: RealClosureResult) -> CatalogueEntry:
    """Merge the M6b TCM-budget result into ``validation.gates.V3`` + metrics.

    M7 only persists the M6b output; the optimisation that produces the
    TCM budget is M6b's responsibility. Until M6b populates non-zero
    ``horizon_tcm_mps`` the V3 block records the (zero) budget and the
    real-closure pass flag.
    """
    new_validation = dict(entry.validation)
    gates = dict(new_validation.get("gates", {}))
    gates["V3"] = {
        "pass": report.closes,
        "max_drift_km": report.max_drift_km,
        "n_cycles": report.n_cycles_propagated,
        "horizon_tcm_mps": report.horizon_tcm_mps,
        "per_cycle_tcm_mps": list(report.per_cycle_tcm_mps),
        "frame_used": report.frame_used,
    }
    new_validation["gates"] = gates
    level = _highest_passing_level(gates)
    _assert_level_evidence(entry.id, level)
    new_validation["level"] = level
    metrics = dict(new_validation.get("metrics", {}))
    metrics["horizon_tcm_dv_mps"] = report.horizon_tcm_mps
    metrics["horizon_laps"] = report.n_cycles_propagated
    new_validation["metrics"] = metrics
    return _with_validation(entry, new_validation)


def record_rediscovery(
    entry: CatalogueEntry,
    run_id: str,
    cell_id: str,
    date: str,
) -> CatalogueEntry:
    """Append a spec §16.4 rediscovery record (idempotent).

    A ``(run_id, cell_id)`` pair already present is not re-appended.
    Attribution (``first_published`` / ``priority_date``) is never
    touched — the rediscovering run only extends the audit trail.
    """
    rediscoveries = list(entry.discovery.get("rediscoveries", []))
    key = (run_id, cell_id)
    if any((r.get("run_id"), r.get("cell_id")) == key for r in rediscoveries):
        return entry
    rediscoveries.append({"run_id": run_id, "cell_id": cell_id, "date": date})
    new_discovery = {**entry.discovery, "rediscoveries": rediscoveries}
    new_raw = {**entry.raw, "discovery": new_discovery}
    return _replace_entry(entry, discovery=new_discovery, raw=new_raw)


def register_discovery(
    entry_skeleton: CatalogueEntry,
    run_id: str,
    cell_id: str,
    date: str,
    finder_version: str,
) -> CatalogueEntry:
    """Promote a candidate to a ``source: this-project`` entry (spec §16.5).

    Sets ``source="this-project"``, ``our_status="candidate-novel"``,
    ``priority_date=date``, ``first_published=None``, and a
    ``discovery_run`` provenance block. The on-disk catalogue writer is
    the consumer; this helper only shapes the entry.
    """
    discovery_run = {
        "run_id": run_id,
        "cell_id": cell_id,
        "date": date,
        "finder_version": finder_version,
    }
    new_raw = {
        **entry_skeleton.raw,
        "source": "this-project",
        "our_status": "candidate-novel",
        "priority_date": date,
        "discovery_run": discovery_run,
    }
    new_raw.pop("first_published", None)
    return _replace_entry(
        entry_skeleton,
        source="this-project",
        our_status="candidate-novel",
        priority_date=date,
        first_published=None,
        discovery_run=discovery_run,
        raw=new_raw,
    )


def serialise_entry_yaml(entry: CatalogueEntry) -> str:
    """Produce the YAML block for ``entry`` (spec §16.1 record shape).

    Dumps the entry's ``raw`` projection — kept in lock-step with the
    typed fields by the writeback helpers above — so round-tripping an
    unmodified entry reproduces it field-for-field (modulo comments,
    which PyYAML does not preserve).
    """
    dumped: str = yaml.safe_dump(
        entry.raw,
        sort_keys=False,
        allow_unicode=True,
        default_flow_style=False,
    )
    return dumped


__all__ = [
    "apply_v0_v1_to_entry",
    "apply_v2_to_entry",
    "apply_v3_to_entry",
    "record_rediscovery",
    "register_discovery",
    "serialise_entry_yaml",
]
