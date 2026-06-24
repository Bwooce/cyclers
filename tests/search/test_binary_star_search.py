"""Tests for the direct fixed-mu binary-star cycler search (#255).

The load-bearing test is the *reproduce-before-search* gate: the winding-number
topology classifier must label the published Earth-Moon (3,1) and (1,1) Ross
cyclers correctly (and as prograde) before it can be trusted to classify
binary-star candidates. The EXPECTED side here is the published (k1,k2) label
of a sourced member, not a value our own code produced.
"""

from __future__ import annotations

import numpy as np
import pytest

import cyclerfinder.core.cr3bp as cr3bp
import cyclerfinder.search.cr3bp_periodic as cp
from cyclerfinder.search.binary_star_search import (
    collinear_lpoints,
    topology_3d,
    winding_topology,
    z_oscillation_count,
)

ROSS_MU = 1.2150584270572e-2


def _system(mu: float) -> cr3bp.CR3BPSystem:
    return cr3bp.CR3BPSystem(mu=mu, primary="Earth", secondary="Moon", l_km=384400.0, t_s=375699.8)


def test_lpoints_earth_moon_ordering() -> None:
    l1, l2, l3 = collinear_lpoints(ROSS_MU)
    # P1 at -mu ~ -0.012, P2 at 1-mu ~ 0.988
    assert -ROSS_MU < l1 < 1.0 - ROSS_MU  # L1 between the primaries
    assert l2 > 1.0 - ROSS_MU  # L2 beyond the secondary
    assert l3 < -ROSS_MU  # L3 beyond the primary
    # sourced Earth-Moon collinear points (Szebehely): L1~0.8369, L2~1.1557
    assert l1 == pytest.approx(0.8369, abs=2e-3)
    assert l2 == pytest.approx(1.1557, abs=2e-3)


@pytest.mark.parametrize(
    ("label", "x0", "c", "hc", "k1", "k2"),
    [
        ("(3,1)", -0.3209891696, 3.161784147013429, 3, 3, 1),
        ("(1,1)", -0.7682140805, 3.151175879508174, 3, 1, 1),
    ],
)
def test_winding_classifier_reproduces_known_em_members(
    label: str, x0: float, c: float, hc: int, k1: int, k2: int
) -> None:
    """The classifier must reproduce the published (k1,k2) labels of the held
    Earth-Moon Ross cyclers, prograde -- the reproduce-before-search gate."""
    system = _system(ROSS_MU)
    o = cp.correct_symmetric_fixed_jacobi(
        system, x0, c, 12.0, ydot0_sign=-1.0, half_crossings=hc, tol=1e-11
    )
    assert o.converged
    state0 = np.array([o.x0, 0.0, 0.0, 0.0, o.ydot0, 0.0])
    topo = winding_topology(ROSS_MU, state0, o.period)
    assert (topo.k1, topo.k2) == (k1, k2), f"{label}: got ({topo.k1},{topo.k2})"
    assert topo.prograde, f"{label}: published cyclers are prograde"
    assert topo.reaches_secondary, f"{label}: a cycler must reach the secondary realm"


def test_z_oscillation_count_conventions() -> None:
    # planar: z identically 0 -> no crossings
    assert z_oscillation_count(np.zeros(100)) == 0
    # one full sine period sin(2*pi*t) over [0,1) has two sign changes
    t = np.linspace(0, 1, 1001, endpoint=False)
    assert z_oscillation_count(np.sin(2 * np.pi * t)) == 2


def test_z_oscillation_count_ignores_exact_zero_samples() -> None:
    # exact-zero samples must not be double-counted; periodic convention closes
    # the loop, so + 0 - over a period crosses going down then back up: 2.
    assert z_oscillation_count(np.array([1.0, 0.0, -1.0])) == 2
    # touching zero and returning to the same sign is zero crossings
    assert z_oscillation_count(np.array([1.0, 0.0, 1.0])) == 0
    # z_center offset shifts the reference plane; signs (-, 0, +) over a period
    # cross up then wrap back down: 2.
    assert z_oscillation_count(np.array([1.0, 2.0, 3.0]), z_center=2.0) == 2


def test_topology_3d_planar_orbit_has_kz_zero() -> None:
    """A genuinely planar (z=0) CR3BP orbit must have k_z == 0 while still
    reproducing the published planar (k1, k2) label."""
    system = _system(ROSS_MU)
    o = cp.correct_symmetric_fixed_jacobi(
        system, -0.3209891696, 3.161784147013429, 12.0, ydot0_sign=-1.0, half_crossings=3, tol=1e-11
    )
    assert o.converged
    state0 = np.array([o.x0, 0.0, 0.0, 0.0, o.ydot0, 0.0])
    topo = topology_3d(ROSS_MU, state0, o.period)
    assert (topo.k1, topo.k2) == (3, 1)
    assert topo.k_z == 0
