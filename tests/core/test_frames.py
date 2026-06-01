"""Rotating-frame transform tests.

These tests address spec §10's "closure-frame correctness" risk. The two
load-bearing checks:

1. ``test_round_trip_identity`` — ``from_rotating(to_rotating(...))`` is the
   identity to ~1e-13 relative tolerance at 1-AU position scale (this is
   the float64 machine-precision floor for arithmetic on ~1.5e8 km values).
   This catches any algebraic mistake in either function.
2. ``test_circular_orbit_is_stationary_in_its_own_frame`` — a body on a
   prograde circular orbit at radius ``a`` with speed ``sqrt(mu/a)``, viewed
   in the frame rotating at its own mean motion, has zero rotating-frame
   velocity. This catches sign errors in the ``omega x r`` correction
   that would silently survive the round-trip identity.

Tolerance note
--------------
Plan §4.1 names ``TOL_ROUND_TRIP_KM = 1e-10``. At 1 AU position scale
(~1.5e8 km) this would require ~6.7e-19 relative precision — below the
float64 machine epsilon (~2.2e-16). Plan §3.2 prose says "identity to ~1e-13
rel" which IS the achievable floor and the constraint actually being
defended; the §4.1 absolute value is an authoring inconsistency. We use the
§3.2 relative form. To get a 1e-10 km absolute bound we'd need a separate
test at sub-km position scale, which would not exercise the trig terms in
their real regime; the relative form is the more meaningful test.

Plan: ``docs/phases/m3-model-construct/plan.md`` §3.2, §4.1.
"""

from __future__ import annotations

from math import pi, sqrt

import numpy as np
import pytest

from cyclerfinder.core.constants import AU_KM, MU_SUN_KM3_S2, PLANETS, SECONDS_PER_DAY
from cyclerfinder.core.frames import from_rotating, synodic_omega, to_rotating

# Tolerances — module-level so loosening is a one-line change.
# Round-trip identity is expressed as a *relative* error against the input
# magnitude (cf. plan §3.2: "identity to ~1e-13 rel"). See module docstring
# for why an absolute km tolerance is not the right shape of test here.
TOL_ROUND_TRIP_REL: float = 1.0e-13
TOL_CIRC_STATIONARY_KMS: float = 1.0e-9
TOL_OMEGA_REL: float = 1.0e-12


def test_round_trip_identity() -> None:
    """``from_rotating ∘ to_rotating`` is the identity to ~1e-13 relative.

    Random states at 1-AU position scale, ~30 km/s velocity scale; several
    values of ``t`` and ``omega`` (including negative omega for retrograde
    frames). This is the primary spec §10 risk-mitigation test.

    The velocity reconstruction picks up an extra ``omega * |r|`` error term
    (from the ``omega x r`` add-back step propagating the rotated-velocity
    error); we bound that explicitly rather than collapsing it into the
    relative form.
    """
    rng = np.random.default_rng(seed=20260531)
    omega_earth = synodic_omega("E")
    # Mix omega scales: zero, Earth, retrograde, Mars-like.
    omega_set = [0.0, omega_earth, -omega_earth, 0.5 * omega_earth, 2.0 * omega_earth]
    # Mix t scales: 0, fraction of year, multi-year.
    t_set_days = [0.0, 50.0, 146.0, 365.0, 800.0, 2500.0]
    for omega in omega_set:
        for t_days in t_set_days:
            t_sec = t_days * SECONDS_PER_DAY
            for _trial in range(8):
                # 1-AU-scale position, ~30 km/s velocity (Earth-orbit scale).
                r = AU_KM * rng.uniform(-1.5, 1.5, size=3).astype(np.float64)
                v = 30.0 * rng.uniform(-1.0, 1.0, size=3).astype(np.float64)
                r_rot, v_rot = to_rotating(r, v, t_sec, omega)
                r_back, v_back = from_rotating(r_rot, v_rot, t_sec, omega)
                r_mag = float(np.linalg.norm(r))
                v_mag = float(np.linalg.norm(v))
                dr = float(np.linalg.norm(r_back - r))
                dv = float(np.linalg.norm(v_back - v))
                assert dr / r_mag < TOL_ROUND_TRIP_REL, (
                    f"omega={omega}, t={t_sec}, |dr|/|r|={dr / r_mag}"
                )
                # Velocity floor: own scale + cross-term from omega*|r| error.
                vel_scale = max(v_mag, abs(omega) * r_mag)
                assert dv / vel_scale < TOL_ROUND_TRIP_REL, (
                    f"omega={omega}, t={t_sec}, |dv|/scale={dv / vel_scale}"
                )


def test_omega_zero_is_identity() -> None:
    """At ``omega = 0`` the transform reduces to identity for any ``t``."""
    rng = np.random.default_rng(seed=42)
    for t_days in (0.0, 100.0, 1000.0):
        t_sec = t_days * SECONDS_PER_DAY
        r = AU_KM * rng.uniform(-1.5, 1.5, size=3).astype(np.float64)
        v = 30.0 * rng.uniform(-1.0, 1.0, size=3).astype(np.float64)
        r_rot, v_rot = to_rotating(r, v, t_sec, 0.0)
        assert np.allclose(r_rot, r, atol=0.0, rtol=0.0)
        assert np.allclose(v_rot, v, atol=0.0, rtol=0.0)


def test_circular_orbit_is_stationary_in_its_own_frame() -> None:
    """A body on a circular orbit, viewed in the frame at its own mean motion,
    has zero rotating-frame velocity.

    This is the physics-anchored sanity test that catches ``omega x r`` sign
    errors invisible to the round-trip identity. Setup: prograde circular
    orbit at radius ``a`` (Earth's heliocentric distance), with speed
    ``sqrt(mu/a)`` perpendicular to the radius vector; mean motion
    ``n = sqrt(mu/a^3)``. In the frame rotating at +n, the body sits at a
    fixed point in the rotating x-y plane and ``v_rot`` is the zero vector.
    """
    a = PLANETS["E"].sma_au * AU_KM  # km
    mu = MU_SUN_KM3_S2
    n = sqrt(mu / a**3)  # rad/s, mean motion
    v_circ = sqrt(mu / a)  # km/s

    # Place the body at several phases of its orbit and propagate.
    for theta0 in (0.0, 0.5, 1.0, 2.5, 3.0):
        # Inertial state at t = (theta0 / n).
        t_sec = theta0 / n
        # Body at angle (n * t_sec) = theta0:
        r_inertial = np.array(
            [a * np.cos(theta0), a * np.sin(theta0), 0.0],
            dtype=np.float64,
        )
        v_inertial = np.array(
            [-v_circ * np.sin(theta0), v_circ * np.cos(theta0), 0.0],
            dtype=np.float64,
        )
        # Transform to the rotating frame at the body's own mean motion.
        _r_rot, v_rot = to_rotating(r_inertial, v_inertial, t_sec, n)
        assert np.linalg.norm(v_rot) < TOL_CIRC_STATIONARY_KMS, (
            f"theta0={theta0}, |v_rot|={np.linalg.norm(v_rot)}"
        )


def test_synodic_omega_earth_matches_constants() -> None:
    """``synodic_omega('E')`` agrees with ``PLANETS['E'].mean_motion_deg_day``."""
    expected = PLANETS["E"].mean_motion_deg_day * (pi / 180.0) / SECONDS_PER_DAY
    got = synodic_omega("E")
    assert abs(got - expected) / expected < TOL_OMEGA_REL


def test_synodic_omega_anchors_on_named_body() -> None:
    """``synodic_omega(body)`` returns that body's own mean motion (the frame
    anchor), not Earth's — Mars and Venus differ from Earth."""
    for code in ("M", "V"):
        expected = PLANETS[code].mean_motion_deg_day * (pi / 180.0) / SECONDS_PER_DAY
        assert abs(synodic_omega(code) - expected) / expected < TOL_OMEGA_REL
    assert synodic_omega("M") != synodic_omega("E")
    assert synodic_omega("V") != synodic_omega("E")


def test_synodic_omega_unknown_body_raises_keyerror() -> None:
    """A body code not in :data:`PLANETS` surfaces as ``KeyError``."""
    with pytest.raises(KeyError):
        synodic_omega("X")
