"""B-plane / powered-flyby kernel — Jones AAS 17-577 Eqs. 1-5 (#142, Phase C).

SELF-CONSISTENCY checks only. Jones gives NO worked numeric example of Eqs. 1-5
(deep-dive §7 "not extractable"), so there is no source-traced EXPECTED value to
assert against. These tests pin the kernel's internal consistency (orthonormal
frame, monotonic r_p, ballistic limits, published-tolerance gate) — never a golden.
"""

from __future__ import annotations

import numpy as np

from cyclerfinder.core.constants import PLANETS
from cyclerfinder.nbody.bplane import (
    MAX_FLYBY_ALT_KM,
    MIN_FLYBY_ALT_KM,
    VINF_MISMATCH_TOL_KMS,
    bplane_angle_rad,
    bplane_frame,
    flyby_altitude_km,
    interior_flyby_feasible,
    periapsis_radius_km,
    tangential_dv_kms,
    turn_angle_rad,
)


def test_turn_angle_zero_for_aligned_asymptotes() -> None:
    v = np.array([3.0, 1.0, -0.5])
    assert turn_angle_rad(v, v) == 0.0


def test_turn_angle_ninety_degrees() -> None:
    a = np.array([4.0, 0.0, 0.0])
    b = np.array([0.0, 4.0, 0.0])
    assert turn_angle_rad(a, b) == np.pi / 2


def test_bplane_frame_is_orthonormal_right_handed() -> None:
    frame = bplane_frame(np.array([2.0, -1.0, 0.7]))
    for u in (frame.s_hat, frame.t_hat, frame.r_hat):
        assert np.isclose(np.linalg.norm(u), 1.0)  # unit vectors
    # Mutual orthogonality.
    assert abs(float(np.dot(frame.s_hat, frame.t_hat))) < 1e-12
    assert abs(float(np.dot(frame.s_hat, frame.r_hat))) < 1e-12
    assert abs(float(np.dot(frame.t_hat, frame.r_hat))) < 1e-12
    # Right-handed: Ŝ x T̂ = R̂.
    assert np.allclose(np.cross(frame.s_hat, frame.t_hat), frame.r_hat)


def test_bplane_frame_handles_pole_aligned_asymptote() -> None:
    # Ŝ ∥ k̂: must still return an orthonormal triad (fallback branch).
    frame = bplane_frame(np.array([0.0, 0.0, 5.0]))
    assert abs(float(np.dot(frame.s_hat, frame.t_hat))) < 1e-12
    assert np.isclose(np.linalg.norm(frame.t_hat), 1.0)


def test_bplane_angle_in_range() -> None:
    theta = bplane_angle_rad(np.array([3.0, 0.5, 0.0]), np.array([2.5, 1.5, 0.4]))
    assert -2.0 * np.pi <= theta <= np.pi


def test_periapsis_radius_larger_for_smaller_bend() -> None:
    # Eq. 2 left side is monotonically decreasing in r_p: a gentler bend needs a
    # LARGER periapsis. Equal-magnitude asymptotes so only the angle changes.
    v_in = np.array([4.0, 0.0, 0.0])
    gentle = np.array([4.0 * np.cos(0.2), 4.0 * np.sin(0.2), 0.0])
    sharp = np.array([4.0 * np.cos(0.8), 4.0 * np.sin(0.8), 0.0])
    rp_gentle = periapsis_radius_km(v_in, gentle, "E")
    rp_sharp = periapsis_radius_km(v_in, sharp, "E")
    assert rp_gentle > rp_sharp > 0.0


def test_periapsis_radius_infinite_for_no_bend() -> None:
    v = np.array([3.0, 0.0, 0.0])
    assert periapsis_radius_km(v, v, "E") == float("inf")


def test_tangential_dv_reduces_to_magnitude_change_when_ballistic() -> None:
    # Pure magnitude change, no bend -> r_p infinite -> Δv = |v_out| - |v_in|.
    v_in = np.array([3.0, 0.0, 0.0])
    v_out = np.array([3.5, 0.0, 0.0])
    assert np.isclose(tangential_dv_kms(v_in, v_out, "E"), 0.5)


def test_tangential_dv_zero_for_identical_asymptotes() -> None:
    v = np.array([4.0, 1.0, 0.0])
    assert np.isclose(tangential_dv_kms(v, v, "M"), 0.0)


def test_flyby_altitude_consistent_with_periapsis_radius() -> None:
    v_in = np.array([4.0, 0.0, 0.0])
    v_out = np.array([4.0 * np.cos(0.5), 4.0 * np.sin(0.5), 0.0])
    rp = periapsis_radius_km(v_in, v_out, "E")
    alt = flyby_altitude_km(v_in, v_out, "E")
    assert np.isclose(alt, rp - PLANETS["E"].radius_eq_km)


def test_interior_feasibility_rejects_large_vinf_mismatch() -> None:
    # Mismatch above the published 200 m/s tolerance -> infeasible.
    v_in = np.array([4.0, 0.0, 0.0])
    v_out = np.array([4.0, 0.5, 0.0])  # |Δ| = 0.5 km/s >> 0.2
    assert np.linalg.norm(v_out - v_in) > VINF_MISMATCH_TOL_KMS
    assert not interior_flyby_feasible(v_in, v_out, "E")


def test_published_tolerances_have_expected_values() -> None:
    # Pin the sourced constants (Jones Eq. 7, deep-dive §4) so a silent drift trips.
    assert VINF_MISMATCH_TOL_KMS == 0.200
    assert MIN_FLYBY_ALT_KM == 100.0
    assert MAX_FLYBY_ALT_KM == 100_000.0
