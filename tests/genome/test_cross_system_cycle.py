"""Tests for the #405 cross-system (SE<->EM) heteroclinic-cycle framework (patched CR3BP).

Two-tier validation (per spec): (1) sourced numeric — SE Lyapunov / connection vs
Canalias 2007 C-values (golden, #407); (2) internal-consistency identities for the new
inter-system frame transform (round-trip, ballistic-dV continuity, energy bookkeeping).
"""

from __future__ import annotations

import numpy as np

from cyclerfinder.genome.cross_system_cycle import (
    FrameBridge,
    em_moon_system,
    se_earth_system,
)


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
