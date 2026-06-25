"""#450 Task 4: global fixed-point enumerator over the (x, xdot) section domain.

Sweeps ``residual(s, n)`` over a box, isolates sub-tolerance cells, dedups into
basins, and emits coarse candidate ICs ``(x0, xdot0, c_target, n)``. Backend-
agnostic (consumes the SectionMap interface).
"""

from __future__ import annotations

import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.genome.da_hotm_backend import SamplingSectionMap
from cyclerfinder.genome.da_hotm_enumerator import DomainBox, enumerate_fixed_points


def _em() -> cr3bp.CR3BPSystem:
    return cr3bp.cr3bp_system("Earth", "Moon")


def test_enumerator_recovers_n1_dro_section_point() -> None:
    """On a coarse box containing the n=1 EM DRO at C=3.00022, emit a candidate
    near the published DRO section point (x0~0.885, xdot0~0)."""
    system = _em()
    backend = SamplingSectionMap(system, c_target=3.00022)
    box = DomainBox(x_lo=0.84, x_hi=0.92, xdot_lo=-0.06, xdot_hi=0.06)
    cands = enumerate_fixed_points(backend, box, n=1, residual_tol=2e-2, grid=(33, 25))
    assert cands, "no candidate emitted for the n=1 DRO box"
    # Some emitted candidate is within tol of the paper's n=1 DRO (x0=0.88500968).
    near = [c for c in cands if abs(c.x0 - 0.88500968) < 0.02 and abs(c.xdot0) < 0.02]
    assert near, [(round(c.x0, 4), round(c.xdot0, 4), c.residual) for c in cands]


def test_candidates_are_deduplicated_into_basins() -> None:
    """A continuum of near-residual cells collapses to one representative each."""
    system = _em()
    backend = SamplingSectionMap(system, c_target=3.00022)
    box = DomainBox(x_lo=0.86, x_hi=0.91, xdot_lo=-0.04, xdot_hi=0.04)
    cands = enumerate_fixed_points(
        backend, box, n=1, residual_tol=3e-2, grid=(21, 17), dedup_radius=0.02
    )
    # The single DRO basin in this tight box yields a small number of reps, not
    # every sub-tol cell.
    assert 1 <= len(cands) <= 3, [(round(c.x0, 4), round(c.xdot0, 4)) for c in cands]


def test_candidate_carries_ic_tuple() -> None:
    system = _em()
    backend = SamplingSectionMap(system, c_target=3.00022)
    box = DomainBox(x_lo=0.86, x_hi=0.91, xdot_lo=-0.03, xdot_hi=0.03)
    cands = enumerate_fixed_points(backend, box, n=1, residual_tol=3e-2, grid=(17, 13))
    c = cands[0]
    assert isinstance(c.x0, float)
    assert isinstance(c.xdot0, float)
    assert c.c_target == 3.00022
    assert c.n == 1
    assert c.residual >= 0.0


def test_empty_box_yields_no_candidates() -> None:
    """A box with no sub-tolerance cell emits nothing (a legitimate negative)."""
    system = _em()
    backend = SamplingSectionMap(system, c_target=3.00022)
    # An off-family box at a tight tolerance: no fixed point here.
    box = DomainBox(x_lo=0.60, x_hi=0.62, xdot_lo=0.30, xdot_hi=0.32)
    cands = enumerate_fixed_points(backend, box, n=1, residual_tol=1e-3, grid=(9, 9))
    assert cands == []
