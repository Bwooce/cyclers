"""M8 multi-body period dispatch tests (spec §3, §9).

The only sourced VEM number we assert against is the beat period
3 x E-M ~ 4 x E-V ~ 6.40 yr (spec §9 line 160). Cross-checked against the
already-passing resonance gate tests/search/test_resonance.py:57-64
(top tuple (4,3) -> 6.406 yr) — same physics, exercised here through the
optimiser's period resolver.
"""

from __future__ import annotations

import pytest

from cyclerfinder.core.constants import DAYS_PER_JULIAN_YEAR, SECONDS_PER_DAY
from cyclerfinder.search.optimize import _target_period_sec
from cyclerfinder.search.sequence import Cell


def _vem_cell(period_k: int = 3, basis: tuple[str, str] | None = ("E", "M")) -> Cell:
    # Anchor-pair convention (plan §2): the catalogue's sourced EMEEVE sequence
    # E-M-E-E-V-E (catalogue.yaml:1777) with its sourced anchor pair (E,M) and
    # sourced k=3. period_k is NOT rewritten — the basis tells the resolver how
    # to interpret it. The spec §13.8 E-V-M-E-M-E string is illustrative only.
    seq = ("E", "M", "E", "E", "V", "E")
    n_legs = len(seq) - 1
    return Cell(
        bodies=("V", "E", "M"),
        sequence=seq,
        period_k=period_k,
        per_leg_revs=(0,) * n_legs,
        per_leg_branch=("single",) * n_legs,
        period_basis=basis,
    )


def test_vem_anchor_pair_period_is_6_41yr() -> None:
    """EMEEVE cell with basis (E,M) and k=3 resolves to 3*T_syn(E,M) ~ 6.41 yr,
    the sourced beat (spec §9, catalogue.yaml:1782)."""
    t_sec = _target_period_sec(_vem_cell(period_k=3, basis=("E", "M")))
    t_yr = t_sec / SECONDS_PER_DAY / DAYS_PER_JULIAN_YEAR
    assert t_yr == pytest.approx(6.41, abs=0.02), f"got {t_yr:.4f} yr"


def test_vem_no_basis_falls_back_to_natural_beat() -> None:
    """With period_basis=None and >=3 bodies, the resolver falls back to the
    natural beat multi_body_beat_days(...)[0] ~ 6.406 yr (spec §9)."""
    t_sec = _target_period_sec(_vem_cell(period_k=1, basis=None))
    t_yr = t_sec / SECONDS_PER_DAY / DAYS_PER_JULIAN_YEAR
    assert t_yr == pytest.approx(6.406, abs=0.01), f"got {t_yr:.4f} yr"


def test_vem_anchor_period_scales_with_period_k() -> None:
    """period_k scales the anchor-pair period linearly (the higher branches,
    spec §3 line 41); doubling k doubles T."""
    t1 = _target_period_sec(_vem_cell(period_k=3, basis=("E", "M")))
    t2 = _target_period_sec(_vem_cell(period_k=6, basis=("E", "M")))
    assert t2 == pytest.approx(2.0 * t1, rel=1e-12)


def test_two_body_period_path_unchanged() -> None:
    """The 2-body fast path is byte-identical to the pre-M8 single-pair formula."""
    from cyclerfinder.search.resonance import synodic_period_days

    cell = Cell(
        bodies=("E", "M"),
        sequence=("E", "M", "E"),
        period_k=2,
        per_leg_revs=(0, 0),
        per_leg_branch=("single", "single"),
    )
    expected = synodic_period_days("E", "M") * 2 * SECONDS_PER_DAY
    assert _target_period_sec(cell) == expected  # exact, not approx
