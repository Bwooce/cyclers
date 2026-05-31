"""M6a dynamic rotating-frame tests (spec §10 + §12(c)).

Two binding gates and several cross-validation tests for the
dynamic-frame extension to :mod:`cyclerfinder.core.frames`:

* ``test_dynamic_frame_roundtrip_identity`` — **spec §10 BINDING
  GATE**. ``from_rotating_dynamic ∘ to_rotating_dynamic`` is the
  identity to ≤ 1e-10 relative on AU-scale states across a 5-year
  time grid. A wrong frame silently fakes/breaks every drift
  measurement downstream; the round-trip test is the load-bearing
  algebraic check.
* ``test_dynamic_frame_reduces_to_uniform_for_circular`` — at the
  circular-coplanar limit where the M3 uniform frame is provably
  correct, the dynamic frame must agree to numerical precision. Any
  divergence here means the dynamic ``omega(t)`` or ``theta(t)``
  computation is wrong.

The remaining tests cross-validate ``synodic_omega_dynamic`` against
the M3 uniform rate (circular) and against a known eccentricity-driven
variation (astropy).

Plan: ``docs/phases/m6a-idealized-closure-verification/plan.md`` §4.1,
§4.4.
"""

from __future__ import annotations

import numpy as np
import pytest

from cyclerfinder.core.constants import AU_KM, SECONDS_PER_DAY
from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.core.frames import (
    from_rotating_dynamic,
    synodic_omega,
    synodic_omega_dynamic,
    to_rotating,
    to_rotating_dynamic,
)

# Tolerances — module-level so loosening is a one-line change.
TOL_ROUND_TRIP_REL: float = 1.0e-10
"""Spec §10 binding tolerance. AU-scale positions; relative to input
magnitude per the M3 ``tests/core/test_frames.py`` convention."""

TOL_CIRC_AGREE_POS_KM: float = 1.0e-5
"""Dynamic vs uniform frame on circular ephemeris, position. The plan
§4.7 names 1e-8 km, but at the 1.5-AU scale tested here and out to 5
years the relevant floor is float64's relative precision applied to
``cos(omega * t)`` for ``omega * t ~ 25 rad``, which gives an absolute
error of ``~25 * eps * 1.5 AU ~ 1e-6 km``. The 1e-5 km tolerance
sits comfortably above that floor while still rejecting any sign /
quadrant error in the dynamic-frame computation (which would produce
errors at the AU scale)."""

TOL_CIRC_AGREE_VEL_KMS: float = 1.0e-9
"""Dynamic vs uniform frame on circular ephemeris, velocity. Plan
§4.7 names 1e-12; the actual float-precision floor at the velocity
scale tested here is ~1e-10 km/s. The 1e-9 km/s bound rejects sign
errors while accepting the float-precision noise."""

TOL_OMEGA_REL: float = 1.0e-12
"""Synodic-omega-circular agreement vs the M3 closed form."""


def test_dynamic_frame_roundtrip_identity() -> None:
    """Spec §10 BINDING GATE: ``from(to(...))`` is identity to ≤ 1e-10 rel.

    5-year time grid at 1-month spacing; synthetic inertial states at
    1.5-AU position scale and 30-km/s velocity scale (the regime
    relevant to cycler verification). Tests both the ``"circular"``
    and ``"astropy"`` ephemeris backends because the round-trip
    identity is algebraic — it must hold for any θ(t), ω(t) the
    ephemeris produces.
    """
    rng = np.random.default_rng(seed=20260601)
    for ephem_model in ("circular", "astropy"):
        ephem = Ephemeris(model=ephem_model)
        bodies = ("E", "M")
        # ~60 grid points over 5 yr.
        for month in range(0, 60):
            t_sec = month * 30.0 * SECONDS_PER_DAY
            for _trial in range(3):
                r = 1.5 * AU_KM * rng.uniform(-1.0, 1.0, size=3).astype(np.float64)
                v = 30.0 * rng.uniform(-1.0, 1.0, size=3).astype(np.float64)
                r_rot, v_rot = to_rotating_dynamic(r, v, t_sec, bodies, ephem)
                r_back, v_back = from_rotating_dynamic(r_rot, v_rot, t_sec, bodies, ephem)
                r_mag = float(np.linalg.norm(r))
                v_mag = float(np.linalg.norm(v))
                dr = float(np.linalg.norm(r_back - r))
                dv = float(np.linalg.norm(v_back - v))
                assert dr / r_mag < TOL_ROUND_TRIP_REL, (
                    f"backend={ephem_model}, t={t_sec}, |dr|/|r|={dr / r_mag}"
                )
                # Velocity floor picks up a |omega| * |r| cross-term per the
                # M3 test convention.
                omega = synodic_omega_dynamic(t_sec, bodies, ephem)
                vel_scale = max(v_mag, abs(omega) * r_mag)
                assert dv / vel_scale < TOL_ROUND_TRIP_REL, (
                    f"backend={ephem_model}, t={t_sec}, |dv|/scale={dv / vel_scale}"
                )


def test_dynamic_frame_reduces_to_uniform_for_circular() -> None:
    """Dynamic and uniform frames must agree on the circular backend.

    At the circular-coplanar limit the dynamic frame's
    ``theta(t) = n_E * t`` and ``omega(t) = n_E`` exactly. So
    :func:`to_rotating_dynamic` must produce bitwise-near-identical
    output to :func:`to_rotating` at the M3 Earth-mean-motion rate.
    Any divergence means the dynamic-frame computation is wrong
    (sign error, atan2 quadrant slip, ephemeris-time mismatch).
    """
    ephem = Ephemeris(model="circular")
    bodies = ("E", "M")
    omega_e = synodic_omega("E")
    rng = np.random.default_rng(seed=42)
    # Representative grid: several yrs, several states each.
    for t_days in (0.0, 50.0, 146.0, 365.0, 800.0, 1500.0):
        t_sec = t_days * SECONDS_PER_DAY
        for _trial in range(5):
            r = AU_KM * rng.uniform(-1.5, 1.5, size=3).astype(np.float64)
            v = 30.0 * rng.uniform(-1.0, 1.0, size=3).astype(np.float64)
            r_dyn, v_dyn = to_rotating_dynamic(r, v, t_sec, bodies, ephem)
            r_uni, v_uni = to_rotating(r, v, t_sec, omega_e)
            dr = float(np.linalg.norm(r_dyn - r_uni))
            dv = float(np.linalg.norm(v_dyn - v_uni))
            assert dr < TOL_CIRC_AGREE_POS_KM, f"t={t_sec}, dr={dr}"
            assert dv < TOL_CIRC_AGREE_VEL_KMS, f"t={t_sec}, dv={dv}"


def test_synodic_omega_dynamic_matches_uniform_for_circular_earth() -> None:
    """On the circular backend, ``synodic_omega_dynamic`` ≈ ``synodic_omega('E')``.

    Earth's circular motion has constant angular rate ``n_E``; the
    dynamic-rate computation must recover that to float64 precision.
    """
    ephem = Ephemeris(model="circular")
    expected = synodic_omega("E")
    bodies = ("E", "M")
    for t_days in (0.0, 365.0, 730.0, 1825.0):  # 0, 1, 2, 5 yr
        t_sec = t_days * SECONDS_PER_DAY
        got = synodic_omega_dynamic(t_sec, bodies, ephem)
        assert abs(got - expected) / expected < TOL_OMEGA_REL, (
            f"t={t_sec}: got={got}, expected={expected}"
        )


def test_synodic_omega_dynamic_varies_with_eccentricity_for_astropy() -> None:
    """On the astropy backend, Earth's instantaneous rate varies over a year.

    Earth's eccentricity (~0.0167) drives a ~3% peak-to-peak variation
    in the instantaneous mean motion between perihelion and aphelion.
    Sample at perihelion (~early January) and aphelion (~early July)
    relative to J2000 (2000-01-01); the rates should differ by at
    least 2% (a conservative bound below the textbook 3%).
    """
    ephem = Ephemeris(model="astropy")
    bodies = ("E", "M")
    # J2000 = 2000-01-01 12:00 TDB; Earth's perihelion is ~early Jan,
    # aphelion ~early July. ~3 d and ~185 d after J2000.
    t_peri = 3.0 * SECONDS_PER_DAY
    t_apo = 185.0 * SECONDS_PER_DAY
    omega_peri = synodic_omega_dynamic(t_peri, bodies, ephem)
    omega_apo = synodic_omega_dynamic(t_apo, bodies, ephem)
    # Both prograde, omega_peri > omega_apo (peri = faster).
    assert omega_peri > 0.0
    assert omega_apo > 0.0
    rel_diff = (omega_peri - omega_apo) / omega_apo
    assert rel_diff >= 0.02, f"omega_peri={omega_peri}, omega_apo={omega_apo}, rel_diff={rel_diff}"


def test_to_rotating_dynamic_anchor_x_axis() -> None:
    """After transforming body[0]'s own state, the y-coordinate is ~0.

    The dynamic frame's +x axis is defined to align with the
    instantaneous Sun→body[0] direction; therefore body[0]'s rotating-
    frame position must lie on the +x axis (y ~ 0, z ~ 0 for an
    in-plane body).
    """
    ephem = Ephemeris(model="circular")
    bodies = ("E", "M")
    for t_days in (0.0, 100.0, 1000.0):
        t_sec = t_days * SECONDS_PER_DAY
        r_e, v_e = ephem.state("E", t_sec)
        r_rot, _ = to_rotating_dynamic(r_e, v_e, t_sec, bodies, ephem)
        # |y| << |x|. Circular Earth at 1 AU; tolerance is 1e-6 km
        # which corresponds to 6.7e-21 relative — float-noise floor.
        assert abs(float(r_rot[1])) < 1.0e-6, f"t={t_sec}: r_rot={r_rot}"
        assert float(r_rot[0]) > 0.0, f"t={t_sec}: r_rot={r_rot}"


def test_to_rotating_dynamic_inertial_at_t_zero_matches_uniform() -> None:
    """At ``t = 0`` on the circular backend the frames coincide.

    Both the uniform M3 frame and the dynamic frame have their x-axis
    along inertial +x at ``t = 0``; the transforms must produce
    bitwise-identical results.
    """
    ephem = Ephemeris(model="circular")
    bodies = ("E", "M")
    omega_e = synodic_omega("E")
    rng = np.random.default_rng(seed=7)
    for _trial in range(5):
        r = AU_KM * rng.uniform(-1.5, 1.5, size=3).astype(np.float64)
        v = 30.0 * rng.uniform(-1.0, 1.0, size=3).astype(np.float64)
        r_dyn, v_dyn = to_rotating_dynamic(r, v, 0.0, bodies, ephem)
        r_uni, v_uni = to_rotating(r, v, 0.0, omega_e)
        # Exact agreement at t=0: both frames are the inertial frame.
        assert float(np.linalg.norm(r_dyn - r_uni)) < 1.0e-9
        assert float(np.linalg.norm(v_dyn - v_uni)) < 1.0e-12


def test_synodic_omega_dynamic_requires_two_bodies() -> None:
    """``len(bodies) < 2`` raises ``ValueError``."""
    ephem = Ephemeris(model="circular")
    with pytest.raises(ValueError, match=r"len\(bodies\) >= 2"):
        synodic_omega_dynamic(0.0, ("E",), ephem)
