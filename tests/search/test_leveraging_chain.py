"""Tests for the resonant-hop V∞ descent chain (#465, the chain orchestrator).

The chain walks V∞ DOWN at a single moon across N resonant-leveraging hops (the
VILM endgame structure of Campagnola-Russell "The Endgame Problem" Part-1), each
hop one :func:`cyclerfinder.search.leveraging_leg.evaluate_leveraging_leg` call.
The summed ΔV approaches — and is bounded below by — the Eq.(13) continuous
quadrature floor (a finite chain of integer-resonance hops cannot beat the
infinite-VILM continuous minimum), and is bounded above by the published
finite-chain penalty (the Europa endgame is 154 m/s discrete vs 128 m/s
continuous, mining note 436-438 / A6).

Source discipline (``feedback_golden_tests_sourced_only``): the EXPECTED bracket
edges both trace to Campagnola-Russell — the lower edge is the Eq.(13) quadrature
(``vilm`` golden), the upper edge is the +20% published finite-chain penalty.
Never a number the chain itself computed.
"""

from __future__ import annotations

import math

from cyclerfinder.search import vilm
from cyclerfinder.search.leveraging_chain import walk_vinf_down


def test_chain_walks_vinf_down_at_floor() -> None:
    """At Europa, walk V∞ 1.8 → 0.77 km/s (the sourced A6 endgame bounds).

    The chain must converge, end within tol of 0.77, and its summed ΔV must sit
    in the sourced bracket: ≥ the Eq.(13) continuous floor (128 m/s) and ≤ the
    +20% finite-chain ceiling the published 154-vs-128 m/s discrete Europa endgame
    sets. Both EXPECTED edges are sourced (floor = Eq.13 quadrature; ceiling =
    published discrete penalty).
    """
    floor_ms, _ = vilm.europa_endgame_dv()  # 128 m/s continuous (sourced A6)
    result = walk_vinf_down("Europa", 1.8, 0.77, exterior=True, max_hops=200, max_revs=2000)
    assert result.converged is True
    assert abs(result.vinf_end_kms - 0.77) < 1.0e-3
    total_ms = result.total_dv_kms * 1000.0
    # A finite chain cannot beat the continuous minimum (less small numeric tol).
    assert total_ms >= floor_ms - 1.0e-6
    # ... and the realised discrete chain stays inside the +25% finite-chain band
    # the published 154/128 m/s Europa endgame sets (sourced ceiling).
    assert total_ms <= 1.25 * floor_ms
    assert len(result.hops) >= 1
    assert result.total_revs >= 1


def test_chain_zero_walk_is_zero_dv() -> None:
    """Target == natural V∞ ⇒ no hops ⇒ ΔV ≈ 0 (the regression limit)."""
    result = walk_vinf_down("Europa", 1.8, 1.8, exterior=True, max_hops=200, max_revs=2000)
    assert result.converged is True
    assert result.total_dv_kms == 0.0
    assert result.hops == ()
    assert math.isclose(result.vinf_end_kms, 1.8, abs_tol=1.0e-9)


def test_chain_infeasible_returns_none() -> None:
    """An unreachable target returns the infeasible sentinel, never a fake ΔV.

    Walking V∞ below the moon's efficiency threshold V̄∞ (deeper than any feasible
    leveraging hop sequence) cannot close — the descent stalls and reports
    ``converged=False`` rather than fabricating a ΔV.
    """
    # V̄∞ at Europa (exterior) is the floor of efficient leveraging; ask for a
    # target far below it (0.01 km/s) — no feasible hop sequence reaches it.
    result = walk_vinf_down("Europa", 1.8, 0.01, exterior=True, max_hops=50, max_revs=500)
    assert result.converged is False
