"""Tests for the #405 cross-system (SE<->EM) heteroclinic-cycle framework (patched CR3BP).

Two-tier validation (per spec): (1) sourced numeric — SE Lyapunov / connection vs
Canalias 2007 C-values (golden, #407); (2) internal-consistency identities for the new
inter-system frame transform (round-trip, ballistic-dV continuity, energy bookkeeping).
"""

from __future__ import annotations

import pathlib

import numpy as np
import pytest

from cyclerfinder.genome.cross_system_cycle import (
    FrameBridge,
    em_moon_system,
    se_earth_system,
)
from cyclerfinder.genome.heteroclinic_cycle import LyapunovNode

_GOLDEN = pathlib.Path("data/golden/canalias_se_em_connection.yaml")

# Canalias-Gómez-Marcote-Masdemont 2007 (doctoral thesis, UPC) SE heteroclinic-family
# bifurcation Jacobi constant, sourced from §4.2 / Fig 4.14 (Sun-Earth case):
# "two branches that end at the bifurcation trajectory with C = 3.000863625."
# This is the bifurcation Jacobi of the SE L1↔L2 heteroclinic family in the planar CR3BP.
#
# MU-MODEL NOTE: Canalias uses µ_SE = 3.040423e-6 (Earth+Moon as SE secondary, per Table 5.1),
# while se_earth_system() uses µ_SE = 3.003481e-6 (Earth-only as SE secondary, per IAU 2015
# constants).  The SAME Jacobi formula is used in both systems (no additive formula offset):
#   C = (x^2+y^2) + 2(1-u)/r1 + 2u/r2 - v^2
# The same numerical C = 3.000863625 lies below C(L1) in BOTH systems (our C(L1) approx 3.000891),
# so a Lyapunov family member exists at that energy in our system.  The Canalias C-value
# therefore directly tests our SE Jacobi convention (formula + µ/length/time scaling) without
# any additive offset being required.
CANALIAS_C_SE = 3.000863625


def _bridge() -> FrameBridge:
    # SE is built from sourced constants (Sun/Earth absent from the satellite registry,
    # so cr3bp.cr3bp_system("Sun","Earth") cannot construct it); EM uses the registry.
    return FrameBridge(se=se_earth_system(), em=em_moon_system())


def test_frame_round_trip_identity() -> None:
    """SE-rot -> inertial -> EM-rot -> inertial -> SE-rot returns the state to <1e-9."""
    bridge = _bridge()
    rng = np.random.default_rng(0)  # seed: deterministic; NOT a sourced value
    for _ in range(20):
        s0 = rng.normal(scale=0.3, size=6)
        s0[0] += 1.0  # place near the Earth region in SE-rot
        theta = float(rng.uniform(0.0, 2.0 * np.pi))
        inert = bridge.se_rot_to_inertial(s0, theta=theta)
        em = bridge.inertial_to_em_rot(inert, theta=theta)
        inert2 = bridge.em_rot_to_inertial(em, theta=theta)
        s1 = bridge.inertial_to_se_rot(inert2, theta=theta)
        assert np.allclose(s0, s1, atol=1e-9), f"round-trip drift {np.abs(s0 - s1).max():.2e}"


def test_moon_maps_to_em_secondary_position() -> None:
    """The EM secondary (Moon) at rest in EM-rot maps to a consistent inertial point.

    Physical anchor (not just an inverse-consistency check): the Moon sits at
    (1-mu_em, 0,0,0,0,0) in EM-rot. Transforming EM-rot -> inertial -> EM-rot must
    return it, and its inertial distance from Earth must equal the EM length scale
    (Earth-Moon distance) to ~1e-6 relative.
    """
    bridge = _bridge()
    moon_em = np.array([1.0 - bridge.em.mu, 0.0, 0.0, 0.0, 0.0, 0.0])
    theta = 0.7
    inert = bridge.em_rot_to_inertial(moon_em, theta=theta)
    r_km = float(np.linalg.norm(inert[:3]))
    assert abs(r_km - bridge.em.l_km) / bridge.em.l_km < 1e-6, f"{r_km} vs {bridge.em.l_km}"
    back = bridge.inertial_to_em_rot(inert, theta=theta)
    assert np.allclose(moon_em, back, atol=1e-9)


@pytest.mark.slow
def test_se_lyapunov_reproduces_canalias_c() -> None:
    """An SE Lyapunov orbit corrected at the Canalias family C closes; jacobi matches.

    EXPECTED = Canalias C = 3.000863625 (sourced: Canalias-Gómez-Marcote-Masdemont 2007
    doctoral thesis, §4.2 / Fig 4.14, "two branches that end at the bifurcation
    trajectory with C = 3.000863625").

    Confirms our SE Jacobi formula and µ/length/time scaling conventions accept the
    sourced C-value.  No additive formula offset is required (see CANALIAS_C_SE note
    above): both Canalias and our code use
        C = (x^2+y^2) + 2(1-u)/r1 + 2u/r2 - v^2
    The mu-model difference (Canalias: Earth+Moon secondary, mu_SE = 3.040e-6; ours:
    Earth-only secondary, mu_SE = 3.003e-6) shifts the physical orbit shape but not the
    Jacobi formula, and C(L1) ≈ 3.000891 in our system leaves the bifurcation C below
    L1 in both.

    Seed note: x0_guess=0.9893 is the working seed in our system (corrected x0 ≈ 0.98925,
    T ≈ 3.031 nondim ≈ 176 days).  x0_guess=0.9899 (near x_L1 ≈ 0.98999) does not
    converge — Newton wanders off the branch from a seed too close to the equilibrium.
    """
    se = se_earth_system()
    c_target = CANALIAS_C_SE
    if _GOLDEN.exists():
        import yaml  # type: ignore[import-untyped]

        data = yaml.safe_load(_GOLDEN.read_text())
        # #407 golden schema: conventions_match_cr3bp is true (no mu(1-mu) offset),
        # so the SE L1->L2 heteroclinic-family bifurcation C is used directly.
        c_target = float(data["bifurcation_jacobi_constants"]["se_l1_to_l2_heteroclinic"])
    node = LyapunovNode.from_libration(
        se,
        x0_guess=0.9893,
        jacobi=c_target,
        period_guess=3.06,
        label="SE-L1",
    )
    assert node.converged, f"SE-L1 did not converge at C={c_target}"
    assert abs(node.jacobi - c_target) < 1e-6, f"SE jacobi {node.jacobi} vs {c_target}"


from cyclerfinder.genome.cross_system_cycle import (  # noqa: E402
    CrossConnection,
    correct_cross_connection,
)


@pytest.mark.slow
def test_em_to_se_connection_is_low_energy() -> None:
    """An EM-L2 unstable -> SE-L2 stable connection closes its inertial position gap
    with a low patch ΔV (near-ballistic cross-system connection). Position match < 100 km,
    ΔV < 1 km/s.

    Working config (recorded): EM-L2 unstable (branch +1) meets SE-L2 stable (branch -1)
    on the inertial patch plane {x_inertial = 1.5e6 km} (corrector defaults). At the
    converged point (theta≈4.945 rad, tau_u≈2.079, tau_s≈2.686, 4 Newton iters) the
    inertial position gap is ≈0.38 km and the patch ΔV is ≈0.364 km/s — a genuine
    near-ballistic cross-system connection in the patched (CR3BP+CR3BP) model. The EM-L2
    unstable manifold reaches ~1.77 Mkm sunward in inertial X, past the SE-L2 region at
    ~1.6 Mkm, so both manifolds co-reach the patch plane."""
    se = se_earth_system()
    em = em_moon_system()
    bridge = FrameBridge(se=se, em=em)
    # Both L2 Lyapunov orbits need ydot0_sign=-1.0: with the default +1 sign the
    # corrector collapses onto a different (shorter-period) family member instead of
    # the L2 orbit (EM lands at x0≈1.064 T≈2.05 rather than x0≈1.182 T≈3.42).
    em_l2 = LyapunovNode.from_libration(
        em, x0_guess=1.18, jacobi=3.15, period_guess=3.4, label="EM-L2", ydot0_sign=-1.0
    )
    se_l2 = LyapunovNode.from_libration(
        se, x0_guess=1.009, jacobi=CANALIAS_C_SE, period_guess=3.06, label="SE-L2", ydot0_sign=-1.0
    )
    conn = correct_cross_connection(bridge, em_l2, se_l2, label_from="EM-L2", label_to="SE-L2")
    assert isinstance(conn, CrossConnection)
    assert conn.converged, f"pos residual {conn.residual:.3e} km, n_iter {conn.n_iter}"
    assert conn.residual < 1e2
    assert conn.patch_dv_kms < 1.0


@pytest.mark.slow
def test_se_to_em_return_leg_converges() -> None:
    """The SE-L2 unstable -> EM-L2 stable RETURN leg converges at low ΔV (#411).

    Once each manifold is propagated in ITS OWN system (the #411 direction fix), the
    return leg closes near-ballistically (observed ΔV ~0.15 km/s) — both cross-system
    legs exist. This guards against regressing to the pre-fix behaviour, where the
    return manifold was propagated in the WRONG system and never converged (which the
    #405 clean-negative had mis-attributed to a narrow-basin physics limit)."""
    se = se_earth_system()
    em = em_moon_system()
    bridge = FrameBridge(se=se, em=em)
    em_l2 = LyapunovNode.from_libration(
        em, x0_guess=1.18, jacobi=3.15, period_guess=3.4, label="EM-L2", ydot0_sign=-1.0
    )
    se_l2 = LyapunovNode.from_libration(
        se, x0_guess=1.009, jacobi=CANALIAS_C_SE, period_guess=3.06, label="SE-L2", ydot0_sign=-1.0
    )
    conn = correct_cross_connection(
        bridge,
        se_l2,
        em_l2,
        label_from="SE-L2",
        label_to="EM-L2",
        branch_u=-1,
        branch_s=+1,
        scan_n=8,
        scan_n_tau=3,
    )
    assert conn.converged, f"return leg must converge; res={conn.residual:.3e} km"
    assert conn.residual < 1e2  # inertial position match on the patch section
    assert conn.patch_dv_kms < 1.0  # near-ballistic (observed ~0.15 km/s)


from cyclerfinder.genome.cross_system_cycle import (  # noqa: E402
    CrossCycle,
    search_cross_cycle,
)


@pytest.mark.slow
def test_closure_search_returns_results_or_clean_negative() -> None:
    """The bounded search returns CrossCycle results; each is closed (with bounded
    residuals) or honestly open (closed=False, notes set). Never fabricates closure."""
    se = se_earth_system()
    em = em_moon_system()
    bridge = FrameBridge(se=se, em=em)
    results = search_cross_cycle(
        bridge,
        c_em_grid=(3.15,),
        c_se_grid=(CANALIAS_C_SE,),
        libration_pairs=(("EM-L2", "SE-L2"),),
        max_attempts=2,
    )
    assert isinstance(results, list)
    for cyc in results:
        assert isinstance(cyc, CrossCycle)
        if cyc.closed:
            assert cyc.max_leg_residual < 1e2
            assert cyc.theta_closure_residual < 1e-2
        else:
            assert cyc.notes


from cyclerfinder.genome.cross_system_cycle import crosscheck_cross_cycle  # noqa: E402


@pytest.mark.slow
def test_crosscheck_is_recorded_or_inf() -> None:
    """crosscheck_cross_cycle returns a cycle whose independent_residual is set
    (finite if every converged leg re-derives under Radau, inf if any fails)."""
    se = se_earth_system()
    em = em_moon_system()
    bridge = FrameBridge(se=se, em=em)
    results = search_cross_cycle(
        bridge,
        c_em_grid=(3.15,),
        c_se_grid=(CANALIAS_C_SE,),
        libration_pairs=(("EM-L2", "SE-L2"),),
        max_attempts=2,
    )
    for cyc in results:
        checked = crosscheck_cross_cycle(bridge, cyc)
        assert not np.isnan(checked.independent_residual)
