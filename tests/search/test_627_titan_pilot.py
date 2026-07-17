"""Task #627: Ross-Roberts-Tsoukkas mu-continuation Titan-pilot test suite.

Covers the two genuinely new, reusable pieces of code this pilot adds:

1. :func:`real_binary_kk_sweep.mu_step_to_system_tracking_c_l1` -- the
   downward-continuation driver that interleaves C-walks (tracking C below
   C_L1(mu)) with each mu step. A FAST regression exercises it over a short,
   already-validated mu range (self-consistency: converged, periodic,
   independent-Radau closes, Jacobi conserved). A SLOW regression (marked
   ``@pytest.mark.slow``, ~7 min) reproduces the actual #627 pilot finding:
   continuing the Ross-RT 2026 Table-I (1,1) mu=0.001 anchor all the way down
   to Saturn-Titan's real mu=2.37e-4 DOES converge to a genuine, stable
   periodic orbit, but that orbit does NOT preserve the (1,1)
   secondary-reaching topology -- a clean, reproducible engineering negative,
   not a numerics failure (see the #627 OUTSTANDING.md bullet + the pilot
   script's own run log for the full writeup).

2. :func:`perimoon_passage.find_perimoon_passage` -- the perimoon-passage
   geometry helper used to judge encounter-relevance (gate b). Verified
   against the ALREADY-ADMITTED, sourced-family Pluto-Charon (3,2) stable
   cycler (#494/#504): its closest approach to Charon must be a small
   positive altitude (a genuine close flyby, not a collision and not a
   distant non-encounter), and refining past the coarse grid must not move
   the answer by more than the coarse grid spacing (internal self-consistency,
   not a sourced golden -- this project's own computed geometry, checked for
   correctness of the SEARCH, not asserted against a published number).
"""

from __future__ import annotations

import numpy as np
import pytest

import cyclerfinder.core.cr3bp as cr3bp
import cyclerfinder.search.cr3bp_periodic as cp
from cyclerfinder.core.cr3bp import cr3bp_system
from cyclerfinder.search.binary_star_search import winding_topology
from cyclerfinder.search.perimoon_passage import find_perimoon_passage
from cyclerfinder.search.pluto_charon_kk_sweep import make_pluto_charon_system
from cyclerfinder.search.real_binary_kk_sweep import mu_step_to_system_tracking_c_l1

# ---------------------------------------------------------------------------
# 1a. mu_step_to_system_tracking_c_l1 -- fast self-consistency regression
# ---------------------------------------------------------------------------


@pytest.mark.slow
def test_627_c_tracking_short_hop_is_self_consistent() -> None:
    """A short, well-behaved mu decrease yields a genuine, continuous periodic orbit.

    Starting from the Ross-RT 2026 Table-I (3,3) anchor (mu=0.012150584270572)
    and stepping down to mu=0.0115 (a small, already-well-behaved decrease --
    NOT the fold region #627 found immediately adjacent to the (1,1) mu=0.01215
    anchor), the landed member must be a genuine perpendicular-crossing
    periodic orbit: small crossing residual, independent-Radau closure, exact
    Jacobi self-consistency, and continuity from the anchor (no discontinuous
    jump in x0/period).
    """
    target = cr3bp.CR3BPSystem(mu=0.0115, primary="P1", secondary="P2", l_km=1.0, t_s=1.0)
    anchor_x0 = -0.322477620583087
    anchor_jacobi = 3.183379082910527
    anchor_period = 19.503763587070285

    landed = mu_step_to_system_tracking_c_l1(
        0.012150584270572,
        target,
        anchor_x0,
        anchor_jacobi,
        anchor_period,
        hc=7,
        sign=-1.0,
        n_steps=10,
        c_margin=0.005,
        tol=1e-10,
    )
    assert landed is not None, "short, well-behaved mu-decrease unexpectedly failed to converge"
    assert landed.converged
    assert landed.crossing_residual < 1e-8

    # Independent-Radau cross-check + exact Jacobi self-consistency.
    state0 = np.array([landed.x0, 0.0, 0.0, 0.0, landed.ydot0, 0.0])
    po = cp.PeriodicOrbit(
        state0=state0,
        period=landed.period,
        jacobi=landed.jacobi,
        converged=True,
        closure_residual=landed.crossing_residual,
    )
    ok_cc, dj = cp.crosscheck_periodic(target, po, closure_tol=1e-6, jacobi_tol=1e-8)
    assert ok_cc, f"independent-Radau crosscheck failed (dj={dj:.2e})"

    c_check = cr3bp.jacobi_constant(state0, target.mu)
    assert abs(c_check - landed.jacobi) < 1e-9

    # Continuity: x0/period must not jump discontinuously over this small hop.
    assert abs(landed.x0 - anchor_x0) < 0.05
    assert abs(landed.period - anchor_period) < 0.5 * anchor_period


def test_627_c_tracking_reduces_to_corrector_at_zero_hop() -> None:
    """n_steps continuation with anchor_mu == target.mu reproduces the anchor.

    A degenerate zero-length mu-hop must just re-run the ordinary fixed-Jacobi
    corrector at the SAME mu and land back on (within corrector tolerance) the
    Table-I anchor itself -- the trivial case every continuation driver must
    get right before trusting it over a real range.
    """
    mu = 0.012150584270572
    target = cr3bp.CR3BPSystem(mu=mu, primary="P1", secondary="P2", l_km=1.0, t_s=1.0)
    x0_pub = -0.322477620583087
    c_pub = 3.183379082910527
    t_pub = 19.503763587070285

    # c_margin=0.001 is smaller than this anchor's own natural gap to C_L1
    # (0.0049), so no C-walk is triggered even by this function's proactive
    # margin-tracking -- the pure zero-hop case.
    landed = mu_step_to_system_tracking_c_l1(
        mu, target, x0_pub, c_pub, t_pub, hc=7, sign=-1.0, n_steps=5, c_margin=0.001, tol=1e-10
    )
    assert landed is not None
    assert abs(landed.x0 - x0_pub) < 1e-6
    assert abs(landed.period - t_pub) < 1e-4
    assert abs(landed.jacobi - c_pub) < 1e-12


# ---------------------------------------------------------------------------
# 1b. The #627 pilot finding itself (slow, ~7 min -- full mu=0.001 -> Titan)
# ---------------------------------------------------------------------------


@pytest.mark.slow
def test_627_titan_pilot_11_continues_but_loses_secondary_reach() -> None:
    """Reproduces the #627 pilot's actual (1,1) finding at Saturn-Titan mu.

    Continuing the Ross-RT 2026 Table-I Rep-1 (1,1) anchor (mu=0.001, the
    paper's own floor) down to the REAL Saturn-Titan mass ratio via
    mu_step_to_system_tracking_c_l1 (200 steps, c_margin=0.02) DOES converge
    to a genuine, linearly-STABLE periodic orbit at the target mu -- but this
    is a clean engineering negative, not a discovery: the landed orbit's
    winding topology is (1,0), not the sought (1,1), i.e. it does NOT reach
    Titan's realm at all (``reaches_secondary`` is False). This is the
    reproducible basis for #627's OUTSTANDING.md verdict; not a golden
    (no published expected value exists at this mu) -- a regression on this
    project's OWN computed, self-consistent result.
    """
    titan = cr3bp_system("Saturn", "Titan")
    landed = mu_step_to_system_tracking_c_l1(
        0.001,
        titan,
        -0.647047499999966,
        3.031605708907296,
        14.774502790974823,
        hc=1,
        sign=-1.0,
        n_steps=200,
        c_margin=0.02,
        tol=1e-10,
    )
    assert landed is not None, "continuation to Titan mu unexpectedly failed to converge"
    assert landed.converged
    assert landed.crossing_residual < 1e-8

    state0 = np.array([landed.x0, 0.0, 0.0, 0.0, landed.ydot0, 0.0])
    topo = winding_topology(titan.mu, state0, landed.period)
    nu, _lam = cp.barden_stability(titan, landed, rtol=1e-13, atol=1e-13)

    # The genuine, reproducible #627 finding: numerically stable...
    assert abs(nu) < 1.0
    # ...but NOT the sought (1,1) secondary-reaching topology.
    assert not topo.reaches_secondary
    assert (topo.k1, topo.k2) != (1, 1)


# ---------------------------------------------------------------------------
# 2. find_perimoon_passage -- verified against the admitted PC (3,2) cycler
# ---------------------------------------------------------------------------


def test_627_perimoon_passage_on_admitted_pluto_charon_cycler() -> None:
    """Perimoon-passage geometry on the already-admitted PC (3,2) stable cycler.

    Uses the #494-derived stable (3,2) member at Pluto-Charon mu=0.10851
    (docs/notes/2026-06-30-494-phase2-3-mu-family-pluto-charon-verdict.md):
    since this orbit's topology is sourced as ``reaches_secondary=True``, its
    closest approach to Charon must be a genuine close passage: a small
    POSITIVE altitude (not below Charon's surface -- this is a real periodic
    orbit, not a collision course) and well inside the CR3BP length scale (not
    a "closest approach" at the far side of the domain, which would signal a
    bug in the search, not a real periapsis).
    """
    system = make_pluto_charon_system()
    x0 = -0.693189765944
    c = 3.579222016200
    t_guess = 11.8366755503

    orbit = cp.correct_symmetric_fixed_jacobi(
        system, x0, c, t_guess, ydot0_sign=-1.0, half_crossings=6, tol=1e-11, rtol=1e-13, atol=1e-13
    )
    assert orbit.converged
    state0 = np.array([orbit.x0, 0.0, 0.0, 0.0, orbit.ydot0, 0.0])
    topo = winding_topology(system.mu, state0, orbit.period)
    assert topo.reaches_secondary, "test fixture orbit must reach the secondary"

    # Charon isn't in SATELLITES (it's a co-orbiting binary, not a "moon" in that
    # registry's sense); Charon mean radius (Nimmo et al. 2017) used directly,
    # only as a sanity floor for the below-surface check.
    charon_radius_km = 606.0

    passage = find_perimoon_passage(
        system, state0, orbit.period, charon_radius_km, rtol=1e-12, atol=1e-12
    )

    # Genuine close passage: positive altitude, well inside one CR3BP length unit.
    assert not passage.below_surface
    assert 0.0 < passage.r2_km < system.l_km
    assert passage.speed_rel_kms > 0.0

    # Internal self-consistency: a finer coarse grid must not move the answer
    # by more than a couple of the ORIGINAL grid's own sample spacings.
    passage_fine = find_perimoon_passage(
        system, state0, orbit.period, charon_radius_km, n_coarse=12000, rtol=1e-12, atol=1e-12
    )
    coarse_dt = orbit.period / 4000.0
    assert abs(passage_fine.t_periapsis - passage.t_periapsis) < 3.0 * coarse_dt
    assert abs(passage_fine.r2_km - passage.r2_km) < 50.0
