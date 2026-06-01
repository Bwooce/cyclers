"""M6b catalogue loader (test infrastructure).

Reads ``data/catalogue.yaml`` and filters to entries in scope for the
M6b real-ephemeris closure verification gate (spec §14 V2-real). Entries
that fall outside M6b's scope (cr3bp, analytic-ephemeris, non-Sun
primary, non-ballistic) are excluded at load time so callers do not
have to re-implement the filter at every test site.

This loader is **test-only**; M7's full ``data/catalog.py`` module (per
spec §16.1) supersedes it. The duplication is intentional: M6b ships
ahead of the catalogue subpackage and the test surface should not depend
on a production-grade loader that does not yet exist.

References
----------
* Plan: ``docs/phases/m6b-real-ephemeris-closure/plan.md`` §1.3, §3.2,
  §4.5.
* Spec §16.1 schema v2 (``model_assumption``, ``trajectory_regime``,
  ``primary`` fields).
* Spec §12.2 (idealised → real-ephemeris representation framework — the
  loader returns only ``circular-coplanar`` entries as the V1 seed
  population).
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Final

import yaml  # type: ignore[import-untyped]

_REPO_ROOT = Path(__file__).resolve().parent.parent.parent
_CATALOGUE_PATH = _REPO_ROOT / "data" / "catalogue.yaml"


M6B_REGRESSION_IDS: Final[tuple[str, ...]] = (
    "aldrin-classic-em-k1-outbound",
    "aldrin-classic-em-k1-inbound",
    "mcconaghy-2006-em-k2",
    "russell-ocampo-2.1.1+2-case2",
    "russell-ocampo-2.5.1+0",
)
"""The 5-entry regression set per plan §3.2.

Hand-picked literature-anchored E-M cyclers most likely to close on real
ephemeris with single-rev Lambert. Entries that turn out to require
multi-rev Lambert or have incomplete leg data are surfaced via the
``EXPECTED_SKIPS`` registry in :mod:`cyclerfinder.verify.real_closure`.
"""


def _entry_in_scope(entry: dict[str, Any]) -> bool:
    """Return True iff the catalogue entry is in M6b's scope.

    Per plan §3.2:

    * ``model_assumption in (None, "circular-coplanar")`` — accepts
      entries that omit the field per schema-v2 defaults and explicit
      circular-coplanar entries; rejects ``cr3bp`` and
      ``analytic-ephemeris``.
    * ``trajectory_regime in (None, "ballistic")`` — rejects future
      low-thrust entries.
    * ``primary in (None, "Sun")`` — heliocentric only; rejects
      Earth-primary lunar / Jovian / Saturnian entries since
      :func:`cyclerfinder.search.phase_match.find_real_windows` only
      supports the Sun.
    """
    model_assumption = entry.get("model_assumption")
    if model_assumption not in (None, "circular-coplanar"):
        return False
    trajectory_regime = entry.get("trajectory_regime")
    if trajectory_regime not in (None, "ballistic"):
        return False
    primary = entry.get("primary")
    return primary in (None, "Sun")


def _passes_v1_filter(entry: dict[str, Any]) -> bool:
    """Return True iff the entry passes V1 OR has no V1 gate record.

    M6a's V1 gate writes ``validation.gates.V1.pass`` per spec §16.1.
    At M6b authorship time M6a may not yet have written the field for
    every entry; treat absent V1 as pass-through so the loader is
    forward-compatible.
    """
    validation = entry.get("validation")
    if validation is None:
        return True
    gates = validation.get("gates") if isinstance(validation, dict) else None
    if not isinstance(gates, dict):
        return True
    v1 = gates.get("V1")
    if not isinstance(v1, dict):
        return True
    passed = v1.get("pass")
    if passed is None:
        return True
    return bool(passed)


def load_m6b_entries(catalogue_path: Path | None = None) -> list[dict[str, Any]]:
    """Load and filter the M6b regression-candidate catalogue entries.

    Reads ``data/catalogue.yaml`` (or the override path), applies the
    plan §3.2 scope filter, and returns the raw entry dicts. Callers
    turn each entry into a :class:`~cyclerfinder.model.cycler.Cycler`
    via
    :func:`cyclerfinder.verify.real_closure.construct_real_ephemeris_cycler`.

    Parameters
    ----------
    catalogue_path:
        Override path to the catalogue YAML. ``None`` (default) reads
        the repository's ``data/catalogue.yaml``.

    Returns
    -------
    list[dict[str, Any]]
        Filtered catalogue entries, in the order they appear in the
        YAML.
    """
    path = catalogue_path or _CATALOGUE_PATH
    raw = yaml.safe_load(path.read_text())
    if not isinstance(raw, list):
        raise ValueError(
            f"expected catalogue YAML to be a list of entries, got {type(raw).__name__}"
        )
    out: list[dict[str, Any]] = []
    for entry in raw:
        if not isinstance(entry, dict):
            continue
        if "id" not in entry:
            continue
        if not _entry_in_scope(entry):
            continue
        if not _passes_v1_filter(entry):
            continue
        out.append(entry)
    return out


__all__ = ["M6B_REGRESSION_IDS", "load_m6b_entries"]
