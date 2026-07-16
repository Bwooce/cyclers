"""Tests for the generalized L1/L2 halo-family-at-Jacobi-constant builder (#613).

`#613` needed a version of `#555`'s/`#612`'s C=3.15-hardcoded halo builders
that accepts an arbitrary target Jacobi constant, to map the L1/L2 quasi-halo
rotation number vs. ``C`` near the L1 planar-Lyapunov -> halo bifurcation
(`#555` independently located this at C~3.1745). This file pins numbers
reproduced LIVE in this session (2026-07-16) by running the module directly --
not copied from a docstring or a prior task's report.

Key headline finding this task's mapping produced (see `data/OUTSTANDING.md`'s
`#613` bullet for the full writeup): the L1 quasi-halo rotation number is NOT
flat vs. ``C`` (only flat vs. amplitude at FIXED ``C``, per `#555`) -- it rises
monotonically from 0 at the bifurcation to ~0.49 near C~3.02, crossing Owen &
Baresi's target 0.2739 at C~3.076. At that SAME ``C``, the L2 family's natural
rotation number is ~0.23 -- nowhere near O&B's 0.02163 (which only holds in a
narrow window near the L2 family's OWN bifurcation, C~3.152). These tests pin
that mismatch directly, confirming the module reproduces it.
"""

from __future__ import annotations

import math

import pytest

import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.genome.qp_tori import correct_qp_torus
from cyclerfinder.search.halo_family_at_jacobi import (
    best_resonance_k,
    halo_center_pair,
    l1_halo_at_jacobi,
    l2_halo_at_jacobi,
    linear_rotation_number,
)

MU = 0.012153643
SYS = cr3bp.CR3BPSystem(mu=MU, primary="Earth", secondary="Moon", l_km=384400.0, t_s=382981.0)

# The C~3.076 crossing this task located: L1's natural rotation number equals
# Owen & Baresi's L1 target 0.2739 here (NOT at C=3.15, per #555).
L1_MATCH_C = 3.07626
OWEN_BARESI_L1_TARGET = 0.2739
OWEN_BARESI_L2_TARGET = 0.02163


def test_l1_halo_at_jacobi_reproduces_612_known_c315_member() -> None:
    """Regression floor: at C=3.15 this module must reproduce `#612`'s
    independently-built halo (x0=0.82401, z0=-0.05453, T=2.76109, live-verified
    in `#612`'s own diagnosis) to within the corrector's own tolerance."""
    r = l1_halo_at_jacobi(SYS, 3.15)
    assert r is not None
    assert r.converged
    assert r.jacobi == pytest.approx(3.15, abs=2e-3)
    assert r.z0 == pytest.approx(-0.0545, abs=2e-3)
    assert pytest.approx(2.761, abs=2e-2) == r.T_TU


def test_l1_halo_at_jacobi_toward_lower_c_direction_is_probed_correctly() -> None:
    """The L1 family's continuation direction (which sign of dx0 increases C)
    must be probed, not assumed: targeting a C BELOW the ~3.145 seed must
    still converge (this is the direction #555/#612 never exercised -- they
    only ever continued TOWARD the bifurcation, i.e. C increasing)."""
    r = l1_halo_at_jacobi(SYS, L1_MATCH_C)
    assert r is not None
    assert r.converged
    assert r.jacobi == pytest.approx(L1_MATCH_C, abs=1e-3)
    # A genuinely larger-amplitude halo than the C=3.15 member (further from
    # the bifurcation -> bigger z0).
    assert abs(r.z0) > 0.10


def test_l1_natural_rotation_number_crosses_owen_baresi_target_near_c_3076() -> None:
    """The headline mapping result: the L1 quasi-halo's natural (linear
    Floquet-phase) rotation number at C~3.076 is close to O&B's L1 target
    0.2739 -- confirmed against an ACTUAL small-amplitude GMOS torus build
    (not just the linear estimate), matching this task's live exploration."""
    r = l1_halo_at_jacobi(SYS, L1_MATCH_C)
    assert r is not None and r.converged
    state0, cands = halo_center_pair(SYS, r)
    assert len(cands) == 2
    linear_rot = linear_rotation_number(cands)
    assert linear_rot is not None
    assert linear_rot == pytest.approx(OWEN_BARESI_L1_TARGET, rel=5e-3)

    # Confirm with an actual small-amplitude GMOS torus (not just the linear
    # estimate) -- the task's "cheap confirmation" step.
    phi = abs(math.atan2(cands[0].imag, cands[0].real))
    k = best_resonance_k(phi)
    torus = correct_qp_torus(
        SYS,
        state0,
        r.T_TU,
        (cands[0], cands[1]),
        k=k,
        n_trans=4,
        initial_torus_amplitude=5e-4,
        tol=1e-7,
        max_iter=40,
    )
    assert torus.invariance_residual < 1e-4
    ratio = abs(torus.omega_trans / torus.omega_long)
    assert ratio == pytest.approx(OWEN_BARESI_L1_TARGET, rel=5e-3)


def test_l2_halo_at_jacobi_reproduces_owen_baresi_l2_regime_at_c315() -> None:
    """Regression floor: at C=3.15 the L2 family's natural rotation number
    sits in O&B's ~0.0216 regime (per #555's independent finding), confirmed
    live via this module."""
    r = l2_halo_at_jacobi(SYS, 3.15)
    assert r is not None and r.converged
    _, cands = halo_center_pair(SYS, r)
    assert len(cands) == 2
    rot = linear_rotation_number(cands)
    assert rot is not None
    assert rot == pytest.approx(0.02, abs=5e-3)


def test_l2_natural_rotation_number_does_not_match_at_the_l1_matched_c() -> None:
    """The decisive negative half of this task's finding: at the SAME C where
    L1 crosses O&B's 0.2739 (~3.076), the L2 family's natural rotation number
    is FAR from O&B's L2 target 0.02163 (off by more than an order of
    magnitude) -- the convenient L2 match `#555` found at C=3.15 does NOT
    generalize to this new C. There is therefore no single Jacobi constant
    where this codebase's per-family natural rotation numbers match BOTH O&B
    targets simultaneously."""
    r = l2_halo_at_jacobi(SYS, L1_MATCH_C)
    assert r is not None and r.converged
    assert r.jacobi == pytest.approx(L1_MATCH_C, abs=1e-3)
    _, cands = halo_center_pair(SYS, r)
    assert len(cands) == 2
    rot = linear_rotation_number(cands)
    assert rot is not None
    # Nowhere near 0.02163 -- off by more than an order of magnitude.
    assert abs(rot - OWEN_BARESI_L2_TARGET) > 0.15


def test_best_resonance_k_picks_nearest_primitive_root() -> None:
    """A phase very close to 2*pi/5 should resolve to k=5."""
    phi = 2 * math.pi / 5 + 1e-4
    assert best_resonance_k(phi) == 5
