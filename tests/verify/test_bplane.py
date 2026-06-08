r"""Tests for the Jones AAS 17-577 B-plane flyby-targeting kernel (GMAT V4 #171).

Golden discipline (binding). Jones tabulates NO worked
:math:`(v_\infty^-, v_\infty^+) \to (\theta_B, r_p, BdotR, BdotT)` example
(deep-dive §7), so every assertion here is either a **mathematical identity**
(frame orthonormality, the B-vector definition) or a **self-consistency
round-trip** (the computed goal reproduces the intended turn) — NEVER a value
sourced from an evaluator's own output. The only sourced inputs are the
constructed :math:`v_\infty` vectors (the caller's published nodes in production).
"""

from __future__ import annotations

import numpy as np
import pytest

from cyclerfinder.core.constants import PLANETS, SAFE_PERIHELION_KM
from cyclerfinder.core.flyby import max_bend
from cyclerfinder.verify.bplane import (
    BPlaneTarget,
    bplane_frame,
    bplane_target,
    bplane_target_for,
)

_MU_MARS = PLANETS["M"].mu_km3_s2
_RP_MIN_MARS = SAFE_PERIHELION_KM["M"]


def _rotate(v: np.ndarray, axis: np.ndarray, angle: float) -> np.ndarray:
    """Rotate ``v`` by ``angle`` about the component of ``axis`` perpendicular to ``v``.

    Projecting the axis perpendicular to ``v`` guarantees the angle between ``v`` and
    the result is exactly ``angle`` (a rotation about an axis with a component along
    ``v`` would turn ``v`` by less). This lets the tests construct a (v_-, v_+) pair
    with a known, exact turn angle.
    """
    v_hat = v / np.linalg.norm(v)
    perp = axis - np.dot(axis, v_hat) * v_hat
    k = perp / np.linalg.norm(perp)
    rotated = (
        v * np.cos(angle)
        + np.cross(k, v) * np.sin(angle)
        + k * np.dot(k, v) * (1.0 - np.cos(angle))
    )
    return np.asarray(rotated, dtype=np.float64)


# --- Task 0.1: frame orthonormality (Jones Eq.4) -----------------------------


def test_bplane_frame_orthonormal() -> None:
    """(S, T, R) are unit-length, mutually orthogonal, right-handed (Jones Eq.4)."""
    vinf_minus = np.array([4.0, -2.0, 1.5], dtype=np.float64)
    s_hat, t_hat, r_hat = bplane_frame(vinf_minus)

    for vec in (s_hat, t_hat, r_hat):
        assert np.linalg.norm(vec) == pytest.approx(1.0, abs=1e-12)

    assert float(np.dot(s_hat, t_hat)) == pytest.approx(0.0, abs=1e-12)
    assert float(np.dot(s_hat, r_hat)) == pytest.approx(0.0, abs=1e-12)
    assert float(np.dot(t_hat, r_hat)) == pytest.approx(0.0, abs=1e-12)

    # S = v_inf_minus_hat; R = S x T (right-handed).
    np.testing.assert_allclose(s_hat, vinf_minus / np.linalg.norm(vinf_minus), atol=1e-12)
    np.testing.assert_allclose(np.cross(s_hat, t_hat), r_hat, atol=1e-12)
    # T is perpendicular to the pole (T . k == 0 by construction, T = (S x k)/|.|).
    assert float(np.dot(t_hat, np.array([0.0, 0.0, 1.0]))) == pytest.approx(0.0, abs=1e-12)


def test_bplane_frame_rejects_pole_parallel() -> None:
    """A v_inf along the pole has no equatorial B-plane frame."""
    with pytest.raises(ValueError, match="parallel to the pole"):
        bplane_frame(np.array([0.0, 0.0, 3.0]))


def test_bplane_frame_rejects_zero() -> None:
    with pytest.raises(ValueError, match="non-zero"):
        bplane_frame(np.array([0.0, 0.0, 0.0]))


# --- Task 0.2: target round-trip (self-consistency, NOT a golden) ------------


def test_bplane_target_roundtrip() -> None:
    r"""The computed B-vector reproduces the intended turn.

    CONSTRUCT a feasible equal-magnitude (v_-, v_+) pair (turn strictly inside the
    Mars bend cone). Per Jones, "the flyby bends v_inf such that the projection of
    v_inf^+ onto the B-plane is along -B". Assert the in-plane (T,R) projection of
    v_hat^+ is antiparallel to B_hat to < 1e-9. Self-consistency round-trip.
    """
    vmag = 5.0
    vinf_minus = np.array([vmag, 0.0, 0.0], dtype=np.float64)
    # A feasible turn: comfortably inside the cone for this v_inf at Mars.
    cone = max_bend(_MU_MARS, _RP_MIN_MARS, vmag)
    turn = 0.5 * cone
    # Rotate within a tilted plane so theta_B is non-trivial (not axis-aligned).
    axis = np.array([0.3, 0.4, 1.0], dtype=np.float64)
    vinf_plus = _rotate(vinf_minus, axis, turn)
    assert np.linalg.norm(vinf_plus) == pytest.approx(vmag, abs=1e-12)

    tgt = bplane_target(vinf_minus, vinf_plus, _MU_MARS, _RP_MIN_MARS)

    # Recovered turn matches the intended turn (Eq.1 reuse).
    assert tgt.turn_rad == pytest.approx(turn, abs=1e-12)
    assert tgt.feasible is True

    # B-vector in the B-plane: B = (B.T) T + (B.R) R == |B| (cos theta_B T + sin theta_B R).
    b_vec = tgt.bdot_t_km * tgt.t_hat + tgt.bdot_r_km * tgt.r_hat
    b_vec_polar = tgt.b_mag_km * (
        np.cos(tgt.theta_b_rad) * tgt.t_hat + np.sin(tgt.theta_b_rad) * tgt.r_hat
    )
    np.testing.assert_allclose(b_vec, b_vec_polar, atol=1e-9)
    b_hat = b_vec / np.linalg.norm(b_vec)

    # Projection of v_hat^+ onto the B-plane (drop the S component).
    vplus_hat = vinf_plus / np.linalg.norm(vinf_plus)
    in_plane = vplus_hat - np.dot(vplus_hat, tgt.s_hat) * tgt.s_hat
    in_plane_hat = in_plane / np.linalg.norm(in_plane)

    # Jones: the in-plane projection of v^+ is along -B.
    assert float(np.dot(in_plane_hat, b_hat)) == pytest.approx(-1.0, abs=1e-9)


def test_bplane_target_solves_eq2() -> None:
    r"""The returned r_p satisfies Jones Eq.2 for the requested turn."""
    vinf_minus = np.array([4.5, 1.0, 0.7], dtype=np.float64)
    vinf_plus = np.array([3.9, 2.0, 1.1], dtype=np.float64)  # unequal |v_inf| allowed
    tgt = bplane_target(vinf_minus, vinf_plus, _MU_MARS, _RP_MIN_MARS)

    vm = float(np.linalg.norm(vinf_minus))
    vp = float(np.linalg.norm(vinf_plus))
    lhs = np.arcsin(_MU_MARS / (_MU_MARS + tgt.rp_km * vm * vm)) + np.arcsin(
        _MU_MARS / (_MU_MARS + tgt.rp_km * vp * vp)
    )
    assert lhs == pytest.approx(tgt.turn_rad, abs=1e-9)


def test_bplane_target_b_mag_definition() -> None:
    """|B| = r_p sqrt(1 + 2 mu / (r_p v_-^2)) and BdotR/BdotT decompose it."""
    vinf_minus = np.array([6.0, -1.0, 0.5], dtype=np.float64)
    vinf_plus = _rotate(vinf_minus, np.array([0.2, 1.0, 0.3]), 0.3)
    tgt = bplane_target(vinf_minus, vinf_plus, _MU_MARS, _RP_MIN_MARS)

    vm = float(np.linalg.norm(vinf_minus))
    expected_b = tgt.rp_km * np.sqrt(1.0 + 2.0 * _MU_MARS / (tgt.rp_km * vm * vm))
    assert tgt.b_mag_km == pytest.approx(expected_b, rel=1e-12)
    assert np.hypot(tgt.bdot_r_km, tgt.bdot_t_km) == pytest.approx(tgt.b_mag_km, rel=1e-12)


# --- Task 0.3: infeasible turn is flagged ------------------------------------


def test_bplane_charges_infeasible_turn() -> None:
    """A turn beyond the bend cone is flagged infeasible (powered TCM needed)."""
    vmag = 9.0  # fast flyby: narrow cone
    cone = max_bend(_MU_MARS, _RP_MIN_MARS, vmag)
    vinf_minus = np.array([vmag, 0.0, 0.0], dtype=np.float64)
    # Demand MORE turn than the cone allows.
    turn = min(np.pi - 1e-3, 2.0 * cone + 0.2)
    vinf_plus = _rotate(vinf_minus, np.array([0.1, 0.2, 1.0]), turn)
    tgt = bplane_target(vinf_minus, vinf_plus, _MU_MARS, _RP_MIN_MARS)

    assert tgt.turn_rad > cone
    assert tgt.feasible is False


def test_bplane_target_for_resolves_registry() -> None:
    """The planet-aware wrapper matches the explicit-mu call for Mars."""
    vinf_minus = np.array([5.2, 0.7, 1.44], dtype=np.float64)
    vinf_plus = _rotate(vinf_minus, np.array([0.4, 0.1, 1.0]), 0.25)
    a = bplane_target_for("M", vinf_minus, vinf_plus)
    b = bplane_target(vinf_minus, vinf_plus, _MU_MARS, _RP_MIN_MARS)
    assert isinstance(a, BPlaneTarget)
    assert a.bdot_r_km == pytest.approx(b.bdot_r_km, rel=1e-12)
    assert a.bdot_t_km == pytest.approx(b.bdot_t_km, rel=1e-12)
    assert a.rp_km == pytest.approx(b.rp_km, rel=1e-12)
