"""Russell free_return_arcs[] descriptor parser (spec §16.7.7).

Maps the catalogue's ``free_return_arcs[]`` Earth-Earth leg descriptors onto the
per-leg ``(n_revs, branch)`` topology + asymmetric ToF seeds that the N-arc
corrector (search/correct.py) consumes. Isolated here so the corrector stays
catalogue-agnostic (spec §4).

Field semantics (Russell 2004 pp.126-127, spec §16.7.7; McConaghy, Russell &
Longuski 2005, JSR 42(4) DOI 10.2514/1.8123, "Full-Revolution Transfers"):
  * ``arc_type`` -- generic (g/G), half-rev (h/H), full-rev (f/F).
  * ``tof_years`` -- Earth-Earth leg ToF in years; g/h arcs only (null for f/F).
  * ``resonance`` -- ``M:N`` resonant orbit; full-rev arcs only. Per the
    primary source: "The first parameter, M, is the number of Earth
    revolutions, and so for Earth-Earth transfers, M also equals the transfer
    time of flight in years. The second parameter, N, is the number of
    spacecraft revolutions" (a = a_E (M/N)^(2/3), their Eq. (14)).
    [#794, 2026-08-09: this module previously read M as spacecraft revs and N
    as years -- reversed vs. BOTH primary sources; harmless for the 1:1 arcs
    but wrong for 3:2 / 2:1. Fixed.]
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any

from cyclerfinder.core.constants import DAYS_PER_JULIAN_YEAR


def arc_to_leg_topology(arc_type: str, *, resonance: str | None) -> tuple[int, str]:
    """Per-leg ``(n_revs, branch)`` for a single Earth-Earth return arc.

    generic / half-rev returns are direct E-E legs (``n_revs=0, branch="single"``).
    A full-rev return is a resonant multi-rev loop: in ``M:N``, ``N`` (after the
    colon) is the spacecraft revolution count (McConaghy/Russell/Longuski 2005:
    "The second parameter, N, is the number of spacecraft revolutions"; Russell
    2004 p.126: "the number following the colon ... represent[s] the number of
    revolutions by the spacecraft"), so ``n_revs=N, branch="low"``.
    """
    if arc_type in ("generic", "half-rev"):
        return (0, "single")
    if arc_type == "full-rev":
        if resonance is None:
            raise ValueError("full-rev arc requires an M:N resonance")
        n_revs = int(resonance.split(":")[1])
        return (n_revs, "low")
    raise ValueError(f"unknown arc_type {arc_type!r}")


def arc_tof_seed_days(arc_type: str, *, tof_years: float | None, resonance: str | None) -> float:
    """ToF seed (days) for a single Earth-Earth return arc.

    For g/h arcs the seed is the sourced ``tof_years * DAYS_PER_JULIAN_YEAR``
    (spec §16.7.7). For f/F arcs ``tof_years`` is null; the ToF is exactly the
    ``M`` (before the colon) of the ``M:N`` resonance in Earth years
    (McConaghy/Russell/Longuski 2005: "M also equals the transfer time of
    flight in years"; a full-rev arc departs and re-meets Earth at the same
    point after exactly M Earth years / N spacecraft revs). This is a *seed
    only*, refined by the corrector.
    """
    if arc_type in ("generic", "half-rev"):
        if tof_years is None:
            raise ValueError(f"{arc_type} arc requires tof_years")
        return tof_years * DAYS_PER_JULIAN_YEAR
    if arc_type == "full-rev":
        if resonance is None:
            raise ValueError("full-rev arc requires an M:N resonance")
        m_years = int(resonance.split(":")[0])
        return m_years * DAYS_PER_JULIAN_YEAR
    raise ValueError(f"unknown arc_type {arc_type!r}")


def designated_arc_index(arcs: Sequence[Mapping[str, Any]]) -> int:
    """Index of the DESIGNATED arc: the UPPERCASE letter in Russell's own
    leg-descriptor notation (Russell 2004 SS4.8 pp.125-127: "The transit times
    and Mars v-inf are calculated using the designated transit leg, as
    indicated by an uppercase descriptor letter"; established for this
    catalogue by #794, first coded by #820). NOT always ``arcs[0]``: e.g. it is
    ``arcs[1]`` (G) for russell-ch4-5.30gGf3 and ``arcs[2]`` (F) for the
    ggF/gfF-pattern rows -- callers must never read ``free_return_arcs``
    positionally to find the designated leg (#849's own defect class)."""
    ups = [i for i, a in enumerate(arcs) if str(a.get("raw_descriptor") or "")[:1].isupper()]
    if len(ups) != 1:
        raise ValueError(
            f"expected exactly one designated (uppercase) arc, found {len(ups)} in "
            f"{[a.get('raw_descriptor') for a in arcs]}"
        )
    return ups[0]


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
