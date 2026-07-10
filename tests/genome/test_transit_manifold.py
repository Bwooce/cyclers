"""Transit-vs-non-transit branch positive control (#547).

Pins the first from-first-principles, textbook-validated positive control for the
Conley-McGehee / Koon-Lo-Marsden-Ross transit-branch classification in the planar
CR3BP -- the sub-problem task #534 flagged but never closed for the
``qp_tori``/``qp_torus_heteroclinic`` linking-number method family.

The Earth-Moon L1 Lyapunov orbit at ``C ~ 3.1869`` sits between ``C_L2 = 3.1722``
and ``C_L1 = 3.1883``: the L1 neck is open, the L2 neck is closed, so the Moon
realm is bounded (the classic KLMR L1-gateway setting). Its ``+`` unstable branch
is a genuine TRANSIT trajectory (threads the neck into the Moon realm, crosses
``x = 1 - mu``, closely approaches the Moon); its ``-`` branch is NON-TRANSIT
(stays interior toward Earth, never crosses). Values are not taken from any
value our own code computed against itself -- the transit/non-transit dichotomy
and the energy ordering are the published KLMR picture; the test asserts the
QUALITATIVE dichotomy plus loose quantitative bounds.
"""

from __future__ import annotations

import numpy as np

import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.genome.transit_manifold import (
    classify_unstable_branch,
    libration_lyapunov,
)
from cyclerfinder.search.reachable_representatives import lagrange_collinear_x

# Earth-Moon collinear Jacobi thresholds (independent of any orbit our code
# corrects): C_L1 = 3.188341..., C_L2 = 3.172160... (verified in-session).
EM_C_L1 = 3.188341106545848
EM_C_L2 = 3.172160451379475


def _em_system() -> cr3bp.CR3BPSystem:
    return cr3bp.cr3bp_system("Earth", "Moon")


def test_l1_lyapunov_is_a_genuine_saddle_at_gateway_energy() -> None:
    """The corrected small L1 Lyapunov is a strong saddle at C in (C_L2, C_L1)."""
    system = _em_system()
    node = libration_lyapunov(system, "L1", amplitude=0.005)
    assert node.converged
    # Energy is above the L1 threshold (lower C) but below L2's -- neck open at
    # L1, closed at L2: the classic bounded-Moon-realm gateway setting.
    assert EM_C_L2 < node.jacobi < EM_C_L1, f"C={node.jacobi} outside gateway band"
    # x0 stayed near L1 (did not slide onto an unrelated family).
    x_l1 = lagrange_collinear_x(system.mu, "L1")
    assert abs(node.state0[0] - x_l1) < 0.02
    # Strongly unstable in-plane saddle (Floquet unstable eigenvalue >> 1).
    from cyclerfinder.genome.heteroclinic_cycle import _planar_floquet_pair

    lam_u, _v_u, lam_s, _v_s = _planar_floquet_pair(system, node.state0, node.period)
    assert lam_u > 100.0, f"expected a strong saddle, got lam_u={lam_u}"
    assert abs(lam_u * lam_s - 1.0) < 1e-2, "monodromy saddle pair must be reciprocal"


def test_transit_and_nontransit_branches_are_distinct() -> None:
    """One unstable branch transits into the Moon realm; the other does not."""
    system = _em_system()
    mu = system.mu
    surface_x = 1.0 - mu
    node = libration_lyapunov(system, "L1", amplitude=0.005)

    plus = classify_unstable_branch(system, node, +1, surface_x=surface_x, t_max=8.0)
    minus = classify_unstable_branch(system, node, -1, surface_x=surface_x, t_max=8.0)

    # Exactly one branch is a transit branch (the KLMR dichotomy).
    assert plus.transits != minus.transits, (
        f"expected one transit + one non-transit branch; got "
        f"+:{plus.transits} (nx={plus.n_crossings}) -:{minus.transits} (nx={minus.n_crossings})"
    )
    transit = plus if plus.transits else minus
    nontransit = minus if plus.transits else plus

    # The transit branch reaches the Moon realm and closely approaches the Moon.
    assert transit.n_crossings >= 1
    assert transit.x_max > surface_x, "transit branch must cross into the Moon realm"
    assert transit.min_secondary_distance < 0.05, (
        f"transit branch should pass near the Moon, min dist={transit.min_secondary_distance}"
    )
    assert np.isfinite(transit.first_crossing_time)

    # The non-transit branch stays interior (never reaches the section) and swings
    # back toward Earth (x well below the secondary).
    assert nontransit.n_crossings == 0
    assert nontransit.x_max < surface_x, "non-transit branch must not cross the section"
    assert nontransit.x_min < node.state0[0], "non-transit branch swings toward Earth"
    assert nontransit.min_secondary_distance > transit.min_secondary_distance


def test_classify_rejects_bad_branch() -> None:
    system = _em_system()
    node = libration_lyapunov(system, "L1", amplitude=0.005)
    try:
        classify_unstable_branch(system, node, 0, surface_x=1.0 - system.mu, t_max=8.0)
    except ValueError:
        pass
    else:  # pragma: no cover
        raise AssertionError("branch=0 must raise ValueError")
