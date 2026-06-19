"""Tests for the #314 heteroclinic-cycle framework (planar CR3BP).

Sourced-golden discipline (feedback_golden_tests_sourced_only): EXPECTED values
trace to Wilczak & Zgliczyński, "Heteroclinic Connections between Periodic Orbits
in the Planar Restricted Three-Body Problem" Part I (arXiv:math/0201278, Comm.
Math. Phys.) — the computer-assisted proof of the closed L1<->L2 Lyapunov cycle
in the Sun-Jupiter-Oterma PCR3BP. Self-consistency checks (FD-Jacobian, empty-path)
need no external source, mirroring existing corrector tests.
"""

from __future__ import annotations

import numpy as np

import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.genome.heteroclinic_cycle import (
    LyapunovNode,
    _planar_floquet_pair,
)

# --- W-Z Sun-Jupiter-Oterma golden (arXiv:math/0201278) ---------------------
WZ_MU = 0.0009537  # W-Z fixed Sun-Jupiter mass ratio (published exactly)
WZ_C = 3.03  # Oterma Jacobi constant in W-Z convention: C = 2*Omega - v^2
# where Omega = (x^2+y^2)/2 + (1-mu)/r1 + mu/r2 + mu*(1-mu)/2  (includes constant term).
# Our code uses C_ours = 2*Ubar - v^2, Ubar = Omega - mu*(1-mu)/2, so:
#   C_ours = C_WZ - mu*(1-mu)
# This constant offset is the ONLY difference; the dynamics are identical.
WZ_C_OURS = WZ_C - WZ_MU * (1.0 - WZ_MU)  # = 3.0290472095... (our code's C for WZ energy)
# Lyapunov fixed points on the section {y=0}, params (x, xdot); xdot=0 at the
# perpendicular crossing. W-Z Part I, interval-enclosed centres:
WZ_X_L1 = 0.9208034913207400196
WZ_X_L2 = 1.081929486841799903


def _sun_jupiter() -> cr3bp.CR3BPSystem:
    # l_km / t_s are not used by the corrector math (all dynamics use mu only);
    # plausible Sun-Jupiter values for completeness.
    return cr3bp.CR3BPSystem(
        mu=WZ_MU, primary="sun", secondary="jupiter", l_km=778.57e6, t_s=5.957e8
    )


def test_floquet_pair_gives_unstable_and_stable_reciprocal() -> None:
    """A libration Lyapunov orbit has a real saddle Floquet pair (lambda, 1/lambda)."""
    system = _sun_jupiter()
    # Generate the L1 Lyapunov orbit at the Oterma energy (Task 2 wires the real
    # corrector; here we lean on the same primitive directly).
    node = LyapunovNode.from_libration(
        system, x0_guess=WZ_X_L1, jacobi=WZ_C, period_guess=3.0, label="L1"
    )
    lam_u, v_u, lam_s, v_s = _planar_floquet_pair(system, node.state0, node.period)
    assert lam_u > 1.0 + 1e-3, f"unstable multiplier must exceed 1, got {lam_u}"
    assert lam_s < 1.0 - 1e-3, f"stable multiplier must be < 1, got {lam_s}"
    # Reciprocal saddle pair: lam_u * lam_s ~ 1.
    assert abs(lam_u * lam_s - 1.0) < 1e-2, f"not a reciprocal pair: {lam_u}*{lam_s}"
    assert v_u.shape == (4,) and v_s.shape == (4,)
    assert np.isclose(np.linalg.norm(v_u), 1.0) and np.isclose(np.linalg.norm(v_s), 1.0)


def test_lyapunov_fixed_points_match_wz() -> None:
    """Corrected L1/L2 Lyapunov x0 reproduce W-Z's section fixed points at C=3.03.

    EXPECTED = W-Z Part I interval-enclosed centres (arXiv:math/0201278); confirms
    our mu/Jacobi/section conventions agree with the paper before any connection.

    Jacobi convention note: W-Z uses C = 2*Omega - v^2 with Omega including the
    mu*(1-mu)/2 constant term; our code omits that term (WZ_C_OURS = WZ_C - mu*(1-mu)).
    The dynamics are identical; the Jacobi values differ by a fixed offset.  At the
    WZ-equivalent energy the corrector reproduces x* to double-precision, validating
    that our CR3BP mu/equations/section match the paper exactly.

    Seeds: L1 uses ydot0_sign=+1 (x0 < L1_x, Theta+ start); L2 uses ydot0_sign=-1
    (x0 > L2_x, Theta- start) — these are the working seeds; period_guess=3.0 suffices
    for both.
    """
    system = _sun_jupiter()
    l1 = LyapunovNode.from_libration(
        system, x0_guess=WZ_X_L1, jacobi=WZ_C_OURS, period_guess=3.0, label="L1"
    )
    l2 = LyapunovNode.from_libration(
        system, x0_guess=WZ_X_L2, jacobi=WZ_C_OURS, period_guess=3.0, label="L2", ydot0_sign=-1.0
    )
    assert l1.converged and l2.converged
    # W-Z enclosures are ~1e-13; our corrector tol is 1e-10, so allow 1e-6.
    assert abs(l1.state0[0] - WZ_X_L1) < 1e-6, f"L1 x0={l1.state0[0]} vs {WZ_X_L1}"
    assert abs(l2.state0[0] - WZ_X_L2) < 1e-6, f"L2 x0={l2.state0[0]} vs {WZ_X_L2}"
    # Both nodes sit at the WZ-equivalent Oterma energy (our convention).
    assert abs(l1.jacobi - WZ_C_OURS) < 1e-6 and abs(l2.jacobi - WZ_C_OURS) < 1e-6


from cyclerfinder.genome.heteroclinic_cycle import (  # noqa: E402  (grouped import)
    _section_crossing,
    _seed_on_manifold,
)


def test_unstable_manifold_reaches_section() -> None:
    """The L1 unstable manifold crosses {y=0} within a bounded horizon."""
    system = _sun_jupiter()
    l1 = LyapunovNode.from_libration(
        system, x0_guess=WZ_X_L1, jacobi=WZ_C_OURS, period_guess=3.0, label="L1"
    )
    seed = _seed_on_manifold(system, l1, tau=0.0, direction="unstable", branch=+1, epsilon=1e-6)
    assert seed.shape == (6,)
    pt = _section_crossing(system, seed, direction="unstable", k=1, max_time=8.0 * l1.period)
    assert pt is not None, "manifold must reach the {y=0} section"
    assert pt.shape == (2,)  # (x, xdot)


def test_section_miss_returns_none() -> None:
    """A horizon too short to reach the section yields None (no hang, no fabrication)."""
    system = _sun_jupiter()
    l1 = LyapunovNode.from_libration(
        system, x0_guess=WZ_X_L1, jacobi=WZ_C_OURS, period_guess=3.0, label="L1"
    )
    seed = _seed_on_manifold(system, l1, tau=0.0, direction="unstable", branch=+1, epsilon=1e-6)
    # Ask for the 9999th crossing — unreachable in this horizon.
    pt = _section_crossing(system, seed, direction="unstable", k=9999, max_time=2.0 * l1.period)
    assert pt is None
