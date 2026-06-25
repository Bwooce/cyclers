"""#450 Task 3: SectionMap interface + SamplingSectionMap backend.

The SectionMap reduces the planar CR3BP return map to ``(x, xdot)`` on the
Poincare section ``Sigma = {y=0, ydot>=0, 0<x<1-mu}`` at a fixed Jacobi constant
(ydot recovered from C, sign = +1 for the ydot>=0 branch). The SamplingSectionMap
is the brute-force float-propagator realization (the reachable_impulsive.py
precedent) used as the validation oracle.
"""

from __future__ import annotations

import numpy as np
import pytest

import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.genome.da_hotm_backend import SamplingSectionMap, SectionPoint


def _em() -> cr3bp.CR3BPSystem:
    return cr3bp.cr3bp_system("Earth", "Moon")


def test_single_rev_matches_direct_propagation() -> None:
    """single_rev(s) equals a direct propagate to the first ydot>=0 y=0 return."""
    system = _em()
    backend = SamplingSectionMap(system, c_target=3.00022)
    s = SectionPoint(x=0.807357887647950, xdot=-0.0956081545978604)
    out = backend.single_rev(s)
    # Reconstruct the lifted IC and propagate directly to the same crossing.
    state0 = backend.lift(s)
    arc = cr3bp.propagate(system, state0, out.t)
    assert abs(arc.state_f[0] - out.point.x) < 1e-9
    assert abs(arc.state_f[3] - out.point.xdot) < 1e-9
    assert abs(arc.state_f[1]) < 1e-8  # back on the y=0 section


def test_compose_equals_sequential_single_revs() -> None:
    """compose(s, n) equals n chained single_rev returns (point + total time)."""
    system = _em()
    backend = SamplingSectionMap(system, c_target=3.00022)
    s = SectionPoint(x=0.85, xdot=-0.05)
    n = 3
    # Manual chain.
    cur = s
    total_t = 0.0
    for _ in range(n):
        step = backend.single_rev(cur)
        cur = step.point
        total_t += step.t
    composed = backend.compose(s, n)
    assert abs(composed.point.x - cur.x) < 1e-9
    assert abs(composed.point.xdot - cur.xdot) < 1e-9
    assert abs(composed.t - total_t) < 1e-9


def test_residual_small_at_known_dro_section_point() -> None:
    """A DRO section point near C~3.0002 returns near itself under n=1.

    JPL DRO at C~3.000114 is x0=0.88494, vy0=0.47062 (mining note JPL
    triangulation). At C=3.00022 the n=1 DRO is a near-fixed point of P^1; the
    single-rev section residual is small (a true fixed point would be ~0).
    """
    system = _em()
    backend = SamplingSectionMap(system, c_target=3.00022)
    s = SectionPoint(x=0.88500968, xdot=0.0)  # paper's n=1 DRO (Table 3)
    r = backend.residual(s, n=1)
    assert r < 5e-2, r  # coarse: a near-fixed point, refined later by the corrector


def test_lift_recovers_jacobi_constant() -> None:
    """The lifted full state has the requested Jacobi constant (ydot from C)."""
    system = _em()
    backend = SamplingSectionMap(system, c_target=3.00022)
    s = SectionPoint(x=0.807357887647950, xdot=-0.0956081545978604)
    state0 = backend.lift(s)
    assert abs(cr3bp.jacobi_constant(state0, system.mu) - 3.00022) < 1e-10
    assert state0[4] > 0.0  # ydot >= 0 branch


def test_lift_raises_when_jacobi_infeasible() -> None:
    """An (x, xdot) where the Jacobi radicand is negative raises ValueError."""
    system = _em()
    backend = SamplingSectionMap(system, c_target=3.00022)
    # Large xdot at this x drives the radicand negative.
    s = SectionPoint(x=0.5, xdot=5.0)
    with pytest.raises(ValueError):
        backend.lift(s)


def test_residual_returns_inf_for_infeasible_point() -> None:
    """An infeasible section point has an infinite (non-finite) residual, not a crash."""
    system = _em()
    backend = SamplingSectionMap(system, c_target=3.00022)
    s = SectionPoint(x=0.5, xdot=5.0)
    r = backend.residual(s, n=1)
    assert not np.isfinite(r)
