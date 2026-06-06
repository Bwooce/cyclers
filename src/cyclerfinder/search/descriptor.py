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
