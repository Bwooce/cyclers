"""Russell free_return_arcs[] descriptor parser (spec §16.7.7).

Maps the catalogue's ``free_return_arcs[]`` Earth-Earth leg descriptors onto the
per-leg ``(n_revs, branch)`` topology + asymmetric ToF seeds that the N-arc
corrector (search/correct.py) consumes. Isolated here so the corrector stays
catalogue-agnostic (spec §4).

Field semantics (Russell 2004 pp.126-127, spec §16.7.7):
  * ``arc_type`` -- generic (g/G), half-rev (h/H), full-rev (f/F).
  * ``tof_years`` -- Earth-Earth leg ToF in years; g/h arcs only (null for f/F).
  * ``resonance`` -- ``M:N`` resonant orbit; full-rev arcs only (M = spacecraft
    revs, N = Earth years).
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any

from cyclerfinder.core.constants import DAYS_PER_JULIAN_YEAR


def arc_to_leg_topology(arc_type: str, *, resonance: str | None) -> tuple[int, str]:
    """Per-leg ``(n_revs, branch)`` for a single Earth-Earth return arc.

    generic / half-rev returns are direct E-E legs (``n_revs=0, branch="single"``).
    A full-rev return is a resonant multi-rev loop: the spacecraft does ``M`` revs
    (``M:N``), so ``n_revs=M, branch="low"`` (spec §16.7.7).
    """
    if arc_type in ("generic", "half-rev"):
        return (0, "single")
    if arc_type == "full-rev":
        if resonance is None:
            raise ValueError("full-rev arc requires an M:N resonance")
        m_revs = int(resonance.split(":")[0])
        return (m_revs, "low")
    raise ValueError(f"unknown arc_type {arc_type!r}")


def arc_tof_seed_days(arc_type: str, *, tof_years: float | None, resonance: str | None) -> float:
    """ToF seed (days) for a single Earth-Earth return arc.

    For g/h arcs the seed is the sourced ``tof_years * DAYS_PER_JULIAN_YEAR``
    (spec §16.7.7). For f/F arcs ``tof_years`` is null; the seed is derived from
    the ``M:N`` resonance: ``M`` spacecraft revs over ``N`` Earth years means the
    Earth-Earth resonant interval is approximately ``N`` Earth years
    (``N * DAYS_PER_JULIAN_YEAR``). This is a *seed only*, refined by the
    corrector.
    """
    if arc_type in ("generic", "half-rev"):
        if tof_years is None:
            raise ValueError(f"{arc_type} arc requires tof_years")
        return tof_years * DAYS_PER_JULIAN_YEAR
    if arc_type == "full-rev":
        if resonance is None:
            raise ValueError("full-rev arc requires an M:N resonance")
        n_years = int(resonance.split(":")[1])
        return n_years * DAYS_PER_JULIAN_YEAR
    raise ValueError(f"unknown arc_type {arc_type!r}")


def parse_free_return_arcs(
    arcs: Sequence[Mapping[str, Any]],
) -> tuple[tuple[int, ...], tuple[str, ...], tuple[float, ...]]:
    """Map a catalogue ``free_return_arcs[]`` list (one arc per Earth-Earth leg)
    onto the three per-leg tuples ``(per_leg_revs, per_leg_branch, tof_seed_days)``
    that the corrector consumes.

    The S1L1 descriptor ``g(1.4612,...) G(2.8096,...)`` yields two generic arcs:
    revs ``(0, 0)``, branches ``("single", "single")``, seeds
    ``[1.4612 yr, 2.8096 yr]`` in days (matching the prototype's pinned arcs,
    ``correct_s1l1_twoarc.py:40``).
    """
    revs: list[int] = []
    branches: list[str] = []
    seeds: list[float] = []
    for arc in arcs:
        arc_type = arc["arc_type"]
        resonance = arc.get("resonance")
        n_revs, branch = arc_to_leg_topology(arc_type, resonance=resonance)
        seed = arc_tof_seed_days(arc_type, tof_years=arc.get("tof_years"), resonance=resonance)
        revs.append(n_revs)
        branches.append(branch)
        seeds.append(seed)
    return tuple(revs), tuple(branches), tuple(seeds)
