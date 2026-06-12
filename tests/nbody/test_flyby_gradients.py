"""Tests for the analytic flyby-continuity gradients (Ellison 2018, #189).

CONSISTENCY TESTS, NOT GOLDENS: Ellison, Conway, Englander & Ozimek (JGCD 41(7),
2018, doi:10.2514/1.G003077) publish no unit-level numeric gradient values (no
worked Appendix example with numbers — mining note §6), so the validation here
is FD-vs-analytic agreement with central differences across randomized states —
the paper's own recommended verification pattern for derivative code (Sec. VI).
Both sides are our own computations; nothing here is a sourced EXPECTED value.

The one cross-implementation anchor is internal: for equal-magnitude v∞ pairs
the Ellison Eq. 4 closed-form periapsis radius must agree with the *independent*
Jones Eq. 2 bisection solver (``nbody/bplane.py:periapsis_radius_km``), since
``asin(μ/(μ + r_p v∞²)) = δ/2`` inverts exactly to ``r_p = (μ/v∞²)(1/sin(δ/2)-1)``
when both asymptote magnitudes are equal. Also a consistency check (two in-repo
implementations), flagged as such.
"""

from __future__ import annotations

from collections.abc import Callable

import numpy as np
import pytest
from numpy.typing import NDArray

from cyclerfinder.nbody.bplane import periapsis_radius_km
from cyclerfinder.nbody.flyby_gradients import (
    flyby_altitude_gradient,
    flyby_continuity_gradients,
    vinf_continuity_gradient,
)

Vec3 = NDArray[np.float64]

# Central-difference relative agreement floor for the analytic gradients. CD
# truncation is O(h²); with a scale-aware step the agreement is comfortably
# inside 1e-6 relative (the #189 acceptance threshold).
_REL_TOL = 1e-6


def _central_diff(f: Callable[[Vec3], float], v: Vec3) -> Vec3:
    """Central-difference gradient of scalar ``f`` at ``v`` (scale-aware step)."""
    h = 1e-6 * max(1.0, float(np.linalg.norm(v)))
    g = np.zeros(3, dtype=np.float64)
    for j in range(3):
        e = np.zeros(3, dtype=np.float64)
        e[j] = h
        g[j] = (f(v + e) - f(v - e)) / (2.0 * h)
    return g


def _rel_err(analytic: Vec3, fd: Vec3) -> float:
    denom = max(float(np.linalg.norm(analytic)), 1e-30)
    return float(np.linalg.norm(analytic - fd)) / denom


def _random_vinf_pairs(
    rng: np.random.Generator, n: int, *, min_turn_deg: float = 5.0, max_turn_deg: float = 175.0
) -> list[tuple[Vec3, Vec3]]:
    """Random v∞ pairs (1-10 km/s) with a non-degenerate turn angle."""
    pairs: list[tuple[Vec3, Vec3]] = []
    while len(pairs) < n:
        v_in = rng.uniform(1.0, 10.0) * _random_unit(rng)
        v_out = rng.uniform(1.0, 10.0) * _random_unit(rng)
        cos_d = float(np.dot(v_in, v_out) / (np.linalg.norm(v_in) * np.linalg.norm(v_out)))
        delta_deg = float(np.degrees(np.arccos(max(-1.0, min(1.0, cos_d)))))
        if min_turn_deg < delta_deg < max_turn_deg:
            pairs.append((v_in, v_out))
    return pairs


def _random_unit(rng: np.random.Generator) -> Vec3:
    v = rng.normal(size=3)
    return np.asarray(v / np.linalg.norm(v), dtype=np.float64)


class TestVinfContinuityGradient:
    def test_value_is_magnitude_defect(self) -> None:
        v_in = np.array([3.0, 0.0, 0.0])
        v_out = np.array([0.0, 4.0, 0.0])
        c, _, _ = vinf_continuity_gradient(v_in, v_out)
        assert c == pytest.approx(1.0)

    def test_gradients_are_signed_unit_vectors(self) -> None:
        v_in = np.array([1.0, 2.0, -2.0])
        v_out = np.array([-3.0, 0.0, 4.0])
        _, d_in, d_out = vinf_continuity_gradient(v_in, v_out)
        np.testing.assert_allclose(d_in, -v_in / 3.0, rtol=1e-15)
        np.testing.assert_allclose(d_out, v_out / 5.0, rtol=1e-15)

    @pytest.mark.parametrize("seed", [1, 2, 3])
    def test_fd_agreement(self, seed: int) -> None:
        # CONSISTENCY test (FD vs analytic; no sourced golden exists — see module
        # docstring).
        rng = np.random.default_rng(seed)
        for v_in, v_out in _random_vinf_pairs(rng, 5):
            _, d_in, d_out = vinf_continuity_gradient(v_in, v_out)

            def c_of_vin(v: Vec3, vo: Vec3 = v_out) -> float:
                return vinf_continuity_gradient(v, vo)[0]

            def c_of_vout(v: Vec3, vi: Vec3 = v_in) -> float:
                return vinf_continuity_gradient(vi, v)[0]

            assert _rel_err(d_in, _central_diff(c_of_vin, v_in)) < _REL_TOL
            assert _rel_err(d_out, _central_diff(c_of_vout, v_out)) < _REL_TOL

    def test_zero_vinf_raises(self) -> None:
        with pytest.raises(ValueError, match="zero-magnitude"):
            vinf_continuity_gradient(np.zeros(3), np.array([1.0, 0.0, 0.0]))


class TestFlybyAltitudeGradient:
    @pytest.mark.parametrize("body", ["E", "V", "M"])
    @pytest.mark.parametrize("seed", [11, 12])
    def test_fd_agreement_random_states(self, body: str, seed: int) -> None:
        # CONSISTENCY test (FD-vs-analytic, central differences, rel tol 1e-6):
        # the Ellison Appendix publishes formulas but no numeric example, so the
        # validation is the paper's own Sec. VI cross-check pattern, never a
        # sourced golden.
        rng = np.random.default_rng(seed)
        worst = 0.0
        for v_in, v_out in _random_vinf_pairs(rng, 8):
            _, d_in, d_out, _, _, active = flyby_altitude_gradient(v_in, v_out, body)
            assert active

            def c_of_vin(v: Vec3, vo: Vec3 = v_out) -> float:
                return flyby_altitude_gradient(v, vo, body)[0]

            def c_of_vout(v: Vec3, vi: Vec3 = v_in) -> float:
                return flyby_altitude_gradient(vi, v, body)[0]

            fd_in = _central_diff(c_of_vin, v_in)
            fd_out = _central_diff(c_of_vout, v_out)
            worst = max(worst, _rel_err(d_in, fd_in), _rel_err(d_out, fd_out))
        assert worst < _REL_TOL

    def test_constraint_value_matches_eq4(self) -> None:
        # Direct re-evaluation of Eq. 4 with the same inputs (formula audit).
        from cyclerfinder.core.constants import PLANETS

        v_in = np.array([4.0, 1.0, 0.5])
        v_out = np.array([1.0, 4.2, -0.3])
        body = "E"
        c, _, _, r_p, delta, active = flyby_altitude_gradient(v_in, v_out, body)
        assert active
        p = PLANETS[body]
        xi = float(np.dot(v_out, v_out))
        expected_rp = (p.mu_km3_s2 / xi) * (1.0 / np.sin(0.5 * delta) - 1.0)
        assert r_p == pytest.approx(expected_rp, rel=1e-14)
        assert c == pytest.approx(expected_rp - (p.radius_eq_km + p.safe_alt_km), rel=1e-12)

    def test_equal_magnitude_periapsis_matches_jones_bisection(self) -> None:
        # CONSISTENCY cross-check between two independent in-repo
        # implementations: Ellison Eq. 4 closed form (this module) vs the Jones
        # Eq. 2 bisection solver (nbody/bplane.py). For ‖v∞⁻‖ = ‖v∞⁺‖ the Jones
        # equation inverts exactly to the Ellison form. Not a sourced golden.
        rng = np.random.default_rng(21)
        for body in ("E", "M"):
            for _ in range(4):
                mag = rng.uniform(2.0, 8.0)
                v_in = mag * _random_unit(rng)
                # Random direction with a moderate turn relative to v_in.
                v_out = mag * _random_unit(rng)
                delta_deg = np.degrees(
                    np.arccos(np.clip(np.dot(v_in, v_out) / (mag * mag), -1.0, 1.0))
                )
                if not 10.0 < delta_deg < 170.0:
                    continue
                _, _, _, r_p, _, _ = flyby_altitude_gradient(v_in, v_out, body)
                r_p_jones = periapsis_radius_km(v_in, v_out, body)
                assert r_p == pytest.approx(r_p_jones, abs=1e-4, rel=1e-6)

    def test_h_safe_override_shifts_constraint_only(self) -> None:
        v_in = np.array([3.0, 1.0, 0.0])
        v_out = np.array([1.0, 3.0, 0.2])
        c0, d_in0, d_out0, rp0, _, _ = flyby_altitude_gradient(v_in, v_out, "E")
        c1, d_in1, d_out1, rp1, _, _ = flyby_altitude_gradient(v_in, v_out, "E", h_safe_km=1000.0)
        assert rp1 == rp0  # periapsis geometry unchanged
        from cyclerfinder.core.constants import PLANETS

        assert c0 - c1 == pytest.approx(1000.0 - PLANETS["E"].safe_alt_km, rel=1e-12)
        np.testing.assert_array_equal(d_in0, d_in1)  # gradient independent of h_safe
        np.testing.assert_array_equal(d_out0, d_out1)

    def test_zero_bend_is_inactive(self) -> None:
        v = np.array([5.0, 0.0, 0.0])
        c, d_in, d_out, r_p, _, active = flyby_altitude_gradient(v, v.copy(), "E")
        assert not active
        assert c == float("inf")
        assert r_p == float("inf")
        np.testing.assert_array_equal(d_in, np.zeros(3))
        np.testing.assert_array_equal(d_out, np.zeros(3))

    def test_antiparallel_cusp_returns_finite_value_zero_gradient(self) -> None:
        v = np.array([5.0, 0.0, 0.0])
        c, d_in, d_out, r_p, delta, active = flyby_altitude_gradient(v, -v, "E")
        assert active
        assert delta == pytest.approx(np.pi)
        assert np.isfinite(c)
        assert r_p == pytest.approx(0.0)  # full reversal needs r_p -> 0
        np.testing.assert_array_equal(d_in, np.zeros(3))
        np.testing.assert_array_equal(d_out, np.zeros(3))

    def test_zero_vinf_raises(self) -> None:
        with pytest.raises(ValueError, match="zero-magnitude"):
            flyby_altitude_gradient(np.array([1.0, 0.0, 0.0]), np.zeros(3), "E")


class TestCombined:
    def test_combined_matches_parts(self) -> None:
        v_in = np.array([2.0, 3.0, -1.0])
        v_out = np.array([-1.0, 3.5, 0.5])
        g = flyby_continuity_gradients(v_in, v_out, "V")
        c_v, dv_in, dv_out = vinf_continuity_gradient(v_in, v_out)
        c_a, da_in, da_out, r_p, delta, active = flyby_altitude_gradient(v_in, v_out, "V")
        assert g.c_vinf_kms == c_v
        assert g.c_altitude_km == c_a
        assert g.r_periapse_km == r_p
        assert g.turn_angle_rad == delta
        assert g.altitude_active == active
        np.testing.assert_array_equal(g.d_cvinf_d_vinf_in, dv_in)
        np.testing.assert_array_equal(g.d_cvinf_d_vinf_out, dv_out)
        np.testing.assert_array_equal(g.d_calt_d_vinf_in, da_in)
        np.testing.assert_array_equal(g.d_calt_d_vinf_out, da_out)
