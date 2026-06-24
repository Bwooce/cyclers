"""Golden tests for the e=0 five-MMR resonant-PO seed base (#440 Phase 1).

Each interior MMR's circular (e=0) symmetric periodic orbit must converge at
mu=0.001 and land at its published semi-major axis a1 (Antoniadou & Libert 2018,
DOI 10.1007/s10569-018-9834-8) within a finite-amplitude tolerance, ON the
correct mean-motion ratio n = p/q. The SOURCED assertions are the published a1
gap locations and the physical resonance ratio n = a_helio**-1.5 = p/q — NOT an
unsourced initial condition. Closure (crossing_residual) must be < 1e-9.

The 3/2 (gate-validated) is the most eccentric member and the tightest test of
the recipe: its osculating-a tracks a1 to ~0.5% while its instantaneous crossing
radius is ~4% off — the test asserts on the osculating a (a_helio), the correct
resonance element.
"""

from __future__ import annotations

import pytest

import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.search.er3bp_isolated_seeds import (
    MMR_SEMI_MAJOR_AXES,
    all_mmr_seeds,
    resonant_po_seed,
)

# Finite-amplitude tolerances (the 3/2 lands ~0.5% off a1; give headroom).
A1_REL_TOL = 0.015  # 1.5% on the heliocentric (osculating) semi-major axis
N_REL_TOL = 0.02  # 2% on the mean-motion ratio n = a_helio**-1.5 vs p/q
CLOSURE_TOL = 1e-9


@pytest.mark.parametrize(("p", "q", "a1"), MMR_SEMI_MAJOR_AXES, ids=lambda v: str(v))
def test_each_mmr_converges_on_resonance(p: int, q: int, a1: float) -> None:
    system = cr3bp.CR3BPSystem(mu=0.001, primary="Sun", secondary="planet", l_km=1.0, t_s=1.0)
    seed = resonant_po_seed(system, p, q, a1)

    assert seed.converged, f"{p}:{q} did not converge (res={seed.crossing_residual:.2e})"
    assert seed.crossing_residual < CLOSURE_TOL, (
        f"{p}:{q} closure {seed.crossing_residual:.2e} >= {CLOSURE_TOL:.0e}"
    )
    # Sourced: lands at the published a1 gap location.
    assert seed.a_helio == pytest.approx(a1, rel=A1_REL_TOL), (
        f"{p}:{q} a_helio {seed.a_helio:.6f} not within {A1_REL_TOL:.1%} of a1={a1}"
    )
    # Sourced physical assertion: the mean-motion ratio IS the resonance p/q.
    assert seed.n_ratio == pytest.approx(p / q, rel=N_REL_TOL), (
        f"{p}:{q} n_ratio {seed.n_ratio:.4f} not within {N_REL_TOL:.1%} of {p / q:.4f}"
    )


def test_three_two_is_gate_validated_member() -> None:
    """The 3/2 must converge to the gate's published member (a~0.7596, n~1.51)."""
    system = cr3bp.CR3BPSystem(mu=0.001, primary="Sun", secondary="planet", l_km=1.0, t_s=1.0)
    seed = resonant_po_seed(system, 3, 2, 0.763143)
    assert seed.converged
    assert seed.crossing_residual < 1e-12  # gate: 1.6e-15
    assert seed.a_helio == pytest.approx(0.7596, abs=2e-3)
    assert seed.n_ratio == pytest.approx(1.51, abs=2e-2)
    assert seed.period == pytest.approx(11.93, abs=0.1)


def test_all_mmr_seeds_returns_all_five_converged() -> None:
    seeds = all_mmr_seeds(mu=0.001)
    assert len(seeds) == len(MMR_SEMI_MAJOR_AXES)
    labels = {s.label for s in seeds}
    assert labels == {"3:2", "5:2", "3:1", "4:1", "5:1"}
    for s in seeds:
        assert s.converged, f"{s.label} failed to converge in all_mmr_seeds"
        assert s.crossing_residual < CLOSURE_TOL
        assert s.n_ratio == pytest.approx(s.p / s.q, rel=N_REL_TOL)
        # IC is a perpendicular x-axis crossing: y=z=xdot=zdot=0.
        assert s.state0[1] == 0.0 and s.state0[2] == 0.0
        assert s.state0[3] == 0.0 and s.state0[5] == 0.0
