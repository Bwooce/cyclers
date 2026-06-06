"""M-3D Phase 1: vector-omega to_rotating is a bit-exact superset (plan §1).

COPLANAR-LIMIT GOLDEN GATE: with omega_vec = (0,0,omega) the vector form must be
numpy.array_equal to the scalar to_rotating across the M3 golden inputs. EXPECTED
side = the scalar function's own output (a refactor-equivalence regression pin,
design §6(b)), never a sourced value.
"""

from __future__ import annotations

import numpy as np

from cyclerfinder.core.frames import (
    synodic_omega,
    to_rotating,
    to_rotating_omega_vec,
)

_RS = [
    np.array([1.496e8, 0.0, 0.0]),
    np.array([0.0, 2.279e8, 0.0]),
    np.array([1.0e8, -1.0e8, 3.0e7]),  # genuine z-component
]
_VS = [
    np.array([0.0, 29.78, 0.0]),
    np.array([-24.07, 0.0, 0.0]),
    np.array([10.0, 10.0, 1.0]),
]


def test_vector_omega_z_only_matches_scalar_bit_for_bit() -> None:
    omega = synodic_omega("E")
    for omega_val in (0.0, omega, -omega, 2.0 * omega):
        for r, v in zip(_RS, _VS, strict=True):
            for t_sec in (0.0, 1.0e6, -3.7e6):
                r_s, v_s = to_rotating(r, v, t_sec, omega_val)
                r_x, v_x = to_rotating_omega_vec(r, v, t_sec, np.array([0.0, 0.0, omega_val]))
                assert np.array_equal(r_s, r_x)
                assert np.array_equal(v_s, v_x)


def test_vector_omega_roundtrip_identity_general_axis() -> None:
    """from(to(x)) == x for a TILTED omega vector (general 3-D frame)."""
    from cyclerfinder.core.frames import from_rotating_omega_vec

    omega_vec = np.array([0.02, -0.015, 0.11]) * synodic_omega("E")
    for r, v in zip(_RS, _VS, strict=True):
        for t_sec in (0.0, 5.0e5, -2.0e6):
            r_rot, v_rot = to_rotating_omega_vec(r, v, t_sec, omega_vec)
            r_back, v_back = from_rotating_omega_vec(r_rot, v_rot, t_sec, omega_vec)
            r_mag = float(np.linalg.norm(r))
            v_mag = float(np.linalg.norm(v))
            omega_mag = float(np.linalg.norm(omega_vec))
            assert float(np.linalg.norm(r_back - r)) / r_mag < 1e-10
            # Velocity reconstruction picks up an omega*|r| cross-term error.
            v_scale = max(v_mag, omega_mag * r_mag)
            assert float(np.linalg.norm(v_back - v)) / v_scale < 1e-10


def test_vector_omega_z_only_inverse_matches_scalar_bit_for_bit() -> None:
    from cyclerfinder.core.frames import from_rotating, from_rotating_omega_vec

    omega = synodic_omega("E")
    for r, v in zip(_RS, _VS, strict=True):
        r_s, v_s = from_rotating(r, v, 1.0e6, omega)
        r_x, v_x = from_rotating_omega_vec(r, v, 1.0e6, np.array([0.0, 0.0, omega]))
        assert np.array_equal(r_s, r_x)
        assert np.array_equal(v_s, v_x)
