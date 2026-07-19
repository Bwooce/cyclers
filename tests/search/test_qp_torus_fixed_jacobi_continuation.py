"""Tests for the free-rho fixed-Jacobi torus continuation + halo rotation
profiler (#615).

`#615` tested two untested Owen & Baresi (2024) reproduction hypotheses,
motivated by a direct re-read of the paper:

* **H1** -- the paper's Section 2.3 says quasi-periodic tori at a FIXED Jacobi
  constant form a 1-parameter family in rotation number, so `#612`'s hard rho
  PIN may have hidden a branch that drifts to O&B's L1 target 0.2739 at large
  amplitude. :func:`continue_qp_torus_fixed_jacobi` walks that genuine free-rho
  branch (pseudo-arclength, Jacobi held by a constraint row, rho free). Result:
  at C=3.15 the rotation number is confined to |rho| ~[0.057, 0.075] and moves
  AWAY from 0.2739 as amplitude grows -- a decisive negative that also refines
  `#612`'s "flat" claim (rho does vary; it is just nowhere near 0.2739).

* **H2** -- the paper's Fig. 5 shows an L2 NRHO, so the §4.1.1 orbits might be
  large-amplitude NRHOs, not the small-amplitude near-bifurcation halos every
  prior task assumed. :func:`halo_family_rotation_profile` shows the halo
  family's Jacobi constant is CAPPED at its planar-Lyapunov bifurcation
  (~3.15), reached at ~zero amplitude; NRHOs (deep |z0|) sit at C well below
  3.15. So no NRHO exists at C=3.15 -- the orbit-type hypothesis is excluded.

All pinned numbers were reproduced LIVE in this session (2026-07-16) by running
the module directly. The rotation numbers are energy-pinned/stable; the
residual floors are truncation/integrator dependent, pinned with generous
bounds (DOP853 stepping is not bit-reproducible across libm versions).
"""

from __future__ import annotations

import itertools
import math

import numpy as np
import pytest
from numpy.typing import NDArray

import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.genome.qp_tori import QPTorus, correct_qp_torus
from cyclerfinder.search.bifurcation_detector import floquet_multipliers, monodromy
from cyclerfinder.search.halo_family_at_jacobi import l1_halo_at_jacobi, l2_halo_at_jacobi
from cyclerfinder.search.nrho_continuation import SymmetricNRHO
from cyclerfinder.search.qp_torus_fixed_jacobi_continuation import (
    continue_qp_torus_fixed_jacobi,
    halo_family_rotation_profile,
    torus_jacobi_and_gradient,
)

MU = 0.012153643
SYS = cr3bp.CR3BPSystem(mu=MU, primary="Earth", secondary="Moon", l_km=384400.0, t_s=382981.0)

OWEN_BARESI_L1_TARGET = 0.2739
OWEN_BARESI_L2_TARGET = 0.02163


def _best_k(phi: float) -> int:
    bk, bd = 5, math.inf
    for kk in range(3, 81):
        for j in range(1, kk):
            if math.gcd(j, kk) == 1 and abs(phi - 2 * math.pi * j / kk) < bd:
                bd, bk = abs(phi - 2 * math.pi * j / kk), kk
    return bk


def _center_pair(r: SymmetricNRHO) -> tuple[NDArray[np.float64], list[complex]]:
    s0 = np.array([r.x0, 0.0, r.z0, 0.0, r.ydot0, 0.0])
    eigs = floquet_multipliers(monodromy(SYS, s0, r.T_TU))
    cands = [complex(e) for e in eigs if abs(e - 1.0) > 1e-3 and abs(e.imag) > 1e-4]
    return s0, cands


def _l1_gmos_seed_at_315() -> QPTorus:
    """Small-amplitude GMOS L1 quasi-halo torus at C=3.15 (the H1 bootstrap)."""
    l1 = l1_halo_at_jacobi(SYS, 3.15)
    assert l1 is not None and l1.converged
    s0, cands = _center_pair(l1)
    assert len(cands) >= 2
    k = _best_k(abs(math.atan2(cands[0].imag, cands[0].real)))
    return correct_qp_torus(
        SYS,
        s0,
        l1.T_TU,
        (cands[0], cands[1]),
        k=k,
        n_trans=4,
        initial_torus_amplitude=5e-4,
        tol=1e-8,
        max_iter=40,
    )


# ---------------------------------------------------------------------------
# Unit: torus Jacobi + analytic gradient.
# ---------------------------------------------------------------------------


def test_torus_jacobi_matches_cr3bp_and_gradient_matches_fd() -> None:
    """``torus_jacobi_and_gradient`` returns the codebase-convention Jacobi
    constant of ``u(0, 0)`` (matching ``cr3bp.jacobi_constant``) and an analytic
    gradient that agrees with a central finite difference."""
    n1, n2 = 3, 2
    a1, a2 = 2 * n1 + 1, 2 * n2 + 1
    # Coeffs whose (theta1,theta2)=(0,0) evaluation is a chosen state: put the
    # state on the constant mode, and sprinkle small cos/sin content elsewhere.
    rng = np.random.default_rng(7)
    coeffs = 1e-3 * rng.normal(size=(6, a1, a2))
    state = np.array([0.86, 0.02, -0.03, 0.01, 0.18, -0.02])
    # Reference angle sums the constant+cos block; zero those extra cos entries
    # and set the constant so u(0,0) equals `state` exactly.
    coeffs[:, 1 : n1 + 1, :] = 0.0
    coeffs[:, :, 1 : n2 + 1] = 0.0
    for c in range(6):
        coeffs[c, 0, 0] = state[c]

    c_val, grad = torus_jacobi_and_gradient(coeffs, n1, n2, MU)
    assert c_val == pytest.approx(cr3bp.jacobi_constant(state, MU), abs=1e-12)

    from cyclerfinder.search.variational_qp_torus import _pack, _unpack

    z = _pack(coeffs, 1.7, 0.3)
    eps = 1e-6
    fd = np.zeros_like(z)
    for i in range(z.size):
        zp, zm = z.copy(), z.copy()
        zp[i] += eps
        zm[i] -= eps
        cp, _ = torus_jacobi_and_gradient(_unpack(zp, n1, n2)[0], n1, n2, MU)
        cm, _ = torus_jacobi_and_gradient(_unpack(zm, n1, n2)[0], n1, n2, MU)
        fd[i] = (cp - cm) / (2 * eps)
    rel = np.max(np.abs(grad - fd)) / max(np.max(np.abs(fd)), 1e-30)
    assert rel < 1e-6, f"jacobi gradient rel error {rel:.2e}"
    # The frequency columns do not affect the Jacobi constant.
    assert grad[-1] == 0.0 and grad[-2] == 0.0


# ---------------------------------------------------------------------------
# H1: free-rho, fixed-C continuation stays far from 0.2739.
# ---------------------------------------------------------------------------


def test_h1_free_rho_continuation_stays_far_from_owen_baresi_l1_target() -> None:
    """The decisive H1 negative: walking the genuine free-rho branch at FIXED
    C=3.15 (rho NOT pinned), the rotation number is confined near the parent
    halo's ~0.074 and DECREASES in magnitude as amplitude grows -- moving AWAY
    from O&B's 0.2739, never toward it. Confirms the paper's 1-parameter family
    is real (rho varies) but its range at C=3.15 does not contain 0.2739."""
    gmos = _l1_gmos_seed_at_315()
    # GMOS seed sits at the energy-pinned L1 rotation number ~0.074 (#555/#612).
    assert abs(gmos.omega_trans / gmos.omega_long) == pytest.approx(0.07402, abs=1e-3)

    steps = continue_qp_torus_fixed_jacobi(
        SYS, gmos, 3.15, n1=10, n2=4, ds=0.01, n_steps=8, max_nfev=120
    )
    assert len(steps) == 8

    # Energy held fixed and the invariance PDE stays converged the whole way.
    for s in steps:
        assert abs(s.jacobi - 3.15) < 1e-4
        assert s.residual_rms < 1e-5

    rots = [s.rotation_number for s in steps]
    amps = [s.transverse_amplitude for s in steps]
    mags = [abs(r) for r in rots]
    # Assert the physics on the FORWARD-MARCHING branch, up to where the
    # arclength continuation folds back near the fragile boundary. #635 pinned
    # the Neimark-Sacker eigenvector PHASE at source (+45-degree convention),
    # which makes the GMOS seed -- and hence this whole continuation -- platform-
    # independent, but as a side effect the seed's arclength tangent now folds
    # the 8-step ds=0.01 march at step ~7 (amp reaches ~0.066 before turning
    # back) instead of marching monotonically past 0.07 as it did under the old
    # BLAS-phase-dependent raw eigenvector. The terminal fold is a seed-
    # parametrization artifact, NOT physics: over the forward-marching prefix the
    # rotation number is identical to before (0.0748 -> ~0.068 as amplitude
    # grows). So assert the decisive negative on the forward branch, robust to
    # where the fold lands. (#635.)
    i_max = int(np.argmax(amps))
    assert i_max >= 4, f"continuation should march several steps before folding (got {i_max})"
    fwd_amps = amps[: i_max + 1]
    fwd_mags = mags[: i_max + 1]
    # Amplitude genuinely grows into the large-amplitude branch (not the pin).
    assert amps[0] < 0.02 and max(amps) > 0.06
    assert all(b > a for a, b in itertools.pairwise(fwd_amps)), "amplitude must increase"
    # rho is reported with the canonical POSITIVE sign (#632: the Neimark-Sacker
    # eigenvalue is pinned to its positive-imaginary representative so the
    # rotation-number sign is reproducible cross-platform); its MAGNITUDE
    # decreases monotonically along the forward branch as amplitude grows.
    assert all(r > 0 for r in rots)
    assert all(b < a for a, b in itertools.pairwise(fwd_mags)), "|rho| must decrease with amp"
    # Energy-pinned endpoints of the forward branch: 0.0748 -> ~0.068.
    assert mags[0] == pytest.approx(0.07484, abs=1e-3)
    assert min(fwd_mags) < 0.070
    # rho genuinely varies (a real 1-parameter family, not a degenerate pin).
    assert max(mags) - min(mags) > 5e-3
    # The whole family stays confined near 0.074, a factor of >3 below O&B's L1
    # target -- the decisive negative: its range at C=3.15 does not contain 0.2739.
    assert all(0.06 < m < 0.08 for m in mags)
    assert max(mags) < 0.10
    assert abs(max(mags) - OWEN_BARESI_L1_TARGET) > 0.15


# ---------------------------------------------------------------------------
# H2: NRHOs do not exist at C=3.15 (family Jacobi ceiling = bifurcation).
# ---------------------------------------------------------------------------


def test_h2_l1_family_jacobi_capped_at_bifurcation_nrhos_below_315() -> None:
    """H2 (L1): marching the L1 halo family into large amplitude, the Jacobi
    constant is CAPPED at ~3.15 (its planar-Lyapunov bifurcation ceiling,
    reached at ~zero amplitude); deep NRHO members (|z0| > 0.2) sit at C well
    below 3.05. So no large-amplitude NRHO exists at C=3.15 -- O&B's C=3.15
    quasi-halos cannot be NRHOs. (En route the rotation number rises to a max
    ~0.498 near C~3.02, crossing 0.2739 there, NOT at C=3.15 -- consistent with
    `#613`.)"""
    seed = l1_halo_at_jacobi(SYS, 3.15)
    assert seed is not None and seed.converged
    prof = halo_family_rotation_profile(SYS, seed, dx0=1e-3, n_steps=90)
    assert len(prof) > 40

    jacobis = [p.jacobi for p in prof]
    # No member exceeds ~3.15 -- the family's energy ceiling is the bifurcation.
    assert max(jacobis) < 3.151
    assert max(jacobis) == pytest.approx(3.15, abs=2e-3)

    # The family genuinely reaches deep NRHO amplitude, and there at low energy.
    deepest = max(prof, key=lambda p: abs(p.z0))
    assert abs(deepest.z0) > 0.20
    assert deepest.jacobi < 3.05

    # Rotation number climbs to a large maximum far from the C=3.15 value 0.074.
    rots = [p.rotation_number for p in prof if p.rotation_number is not None]
    assert max(rots) > 0.40


def test_h2_l2_rotation_number_minimised_at_the_c315_seed() -> None:
    """H2 (L2): the L2 family's rotation number is MINIMISED at the small-
    amplitude C=3.15 seed (~0.0198, in O&B's 0.02163 regime); marching to larger
    amplitude only INCREASES it (to >0.3). So O&B's L2 target lives at the
    near-bifurcation small-amplitude member, not on an NRHO."""
    seed = l2_halo_at_jacobi(SYS, 3.15)
    assert seed is not None and seed.converged
    prof = halo_family_rotation_profile(SYS, seed, dx0=-1e-3, n_steps=70)
    assert len(prof) > 30

    rots = [p.rotation_number for p in prof if p.rotation_number is not None]
    # Seed rotation number near O&B's L2 target, and it is the family minimum.
    assert prof[0].rotation_number == pytest.approx(OWEN_BARESI_L2_TARGET, abs=3e-3)
    assert prof[0].rotation_number == pytest.approx(min(rots), abs=1e-4)
    # Larger amplitude drives the rotation number far ABOVE the target (an
    # order of magnitude up; how far the march gets before the near-Moon region
    # slows the standard propagator varies, so pin a robust lower bound).
    assert max(rots) > 0.20
    # Energy ceiling again at the bifurcation.
    assert max(p.jacobi for p in prof) < 3.151
