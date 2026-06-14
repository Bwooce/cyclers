"""Tests for the asymmetric (general) CR3BP periodic-orbit corrector (#249).

The corrector drops the symmetric perpendicular-crossing assumption (``xdot0 =
0``) and finds a y=0 return-map fixed point in the free variables ``(x0,
xdot0)`` at fixed Jacobi. Symmetric orbits fall out as the ``xdot0 ~= 0``
special case -- the reproduce-before-trust gate.

Convention notes / conditioning
--------------------------------
For a symmetric member the perpendicular half-period crossing is at y=0 crossing
index ``H`` (the symmetric corrector's ``half_crossings``); the FULL state
returns to the IC at crossing ``2H``. So the general corrector closes at
``half_crossings = 2H`` and its period equals the symmetric period ``T``.

The multi-crossing return arc is ill-conditioned (condition number ~1e7 for the
C11a member), so the corrector's quadratic basin is narrow. The reproduce test
therefore seeds the general corrector from the symmetric corrector's CONVERGED
``x0`` (the "known symmetric member"), and the perturbation test uses a small
perturbation inside the basin. Both are honest: the deliverable is recovering a
known orbit and confirming Jacobi/closure, not a wide capture radius.
"""

from __future__ import annotations

import numpy as np
import pytest

import cyclerfinder.core.cr3bp as cr3bp
import cyclerfinder.search.cr3bp_periodic as cp
from cyclerfinder.search.cr3bp_general_periodic import (
    GeneralPeriodicOrbit,
    _ydot0_general,
    correct_general_periodic,
)

# Earth-Moon Ross system (mu, scales) and the confirmed C11a member.
_SYS = cr3bp.CR3BPSystem(
    mu=1.2150584270572e-2, primary="Earth", secondary="Moon", l_km=384400.0, t_s=375699.8
)
_C = 3.1294
_C11A_X0_SEED = -0.8116407  # rounded literal; refined by the symmetric corrector
_C11A_HALF_CROSSINGS_SYM = 3  # perpendicular half-period crossing index
_C11A_HALF_CROSSINGS_GEN = 2 * _C11A_HALF_CROSSINGS_SYM  # full-return crossing index


def _symmetric_c11a() -> cp.SymmetricOrbit:
    sym = cp.correct_symmetric_fixed_jacobi(
        _SYS, _C11A_X0_SEED, _C, 9.69, ydot0_sign=-1.0, half_crossings=_C11A_HALF_CROSSINGS_SYM
    )
    assert sym.converged
    return sym


def test_reproduce_symmetric_c11a_as_xdot0_zero() -> None:
    """Reproduce-before-trust: seeded at the known symmetric C11a member with
    xdot0_guess=0, the general corrector must converge to the SAME orbit and
    report xdot0 ~= 0 (the symmetric special case)."""
    sym = _symmetric_c11a()
    g = correct_general_periodic(
        _SYS,
        sym.x0,
        0.0,
        _C,
        sym.period,
        half_crossings=_C11A_HALF_CROSSINGS_GEN,
        ydot0_sign=-1.0,
    )
    assert isinstance(g, GeneralPeriodicOrbit)
    assert g.converged
    # Recovers the symmetric special case: xdot0 ~= 0 (asymmetry vanishes).
    assert abs(g.xdot0) < 1e-8
    assert g.asymmetry < 1e-8
    assert g.y0 == 0.0
    # Same orbit as the symmetric corrector (tight).
    assert abs(g.x0 - sym.x0) < 1e-9
    assert abs(g.period - sym.period) < 1e-7
    # Held Jacobi (independent recompute) and an independent full-period closure.
    assert abs(g.jacobi - _C) < 1e-9
    state0 = np.array([g.x0, 0.0, 0.0, g.xdot0, g.ydot0, 0.0])
    assert abs(cr3bp.jacobi_constant(state0, _SYS.mu) - _C) < 1e-9
    assert g.closure_residual < 1e-7


def test_small_asymmetric_perturbation_holds_jacobi_and_closes() -> None:
    """A small xdot0 perturbation stays converged; Jacobi is held (independent
    recompute) and an independent full-period re-propagation closes."""
    sym = _symmetric_c11a()
    g = correct_general_periodic(
        _SYS,
        sym.x0,
        1e-8,  # inside the (narrow) quadratic basin
        _C,
        sym.period,
        half_crossings=_C11A_HALF_CROSSINGS_GEN,
        ydot0_sign=-1.0,
    )
    assert g.converged
    # Jacobi held: re-derive C from the corrected IC independently.
    state0 = np.array([g.x0, 0.0, 0.0, g.xdot0, g.ydot0, 0.0])
    assert abs(cr3bp.jacobi_constant(state0, _SYS.mu) - _C) < 1e-9
    # Independent full-period re-propagation closes.
    assert np.isfinite(g.closure_residual)
    assert g.closure_residual < 1e-7
    # The auto-closing ydot component matches at the crossing.
    assert g.ydot_residual < 1e-8


def test_infeasible_jacobi_returns_cleanly() -> None:
    """A Jacobi constant unattainable at the seed (negative radicand) must be
    handled cleanly: the algebraic helper raises, and the corrector returns
    converged=False without raising."""
    # The algebraic ydot0 solve raises on a negative radicand.
    with pytest.raises(ValueError, match="negative radicand"):
        _ydot0_general(-0.81, 0.0, 100.0, _SYS.mu, -1.0)
    # The corrector swallows it and returns a clean non-converged result.
    g = correct_general_periodic(
        _SYS, -0.81, 0.0, 100.0, 9.69, half_crossings=_C11A_HALF_CROSSINGS_GEN, ydot0_sign=-1.0
    )
    assert not g.converged
    assert np.isnan(g.ydot0)
    assert np.isnan(g.jacobi)
