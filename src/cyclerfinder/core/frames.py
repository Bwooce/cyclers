"""Uniform synodic rotating-frame transforms.

A right-handed heliocentric inertial state :math:`(r, v)` is mapped into a
**uniform** rotating frame about the heliocentric z-axis at constant angular
rate :math:`\\omega` (rad/s). The frame's x-axis is aligned with the inertial
x-axis at :math:`t = 0`. The transform is:

.. math::

    r_{\\text{rot}} &= R(-\\omega t) \\, r_{\\text{inertial}} \\\\
    v_{\\text{rot}} &= R(-\\omega t) \\, (v_{\\text{inertial}}
                       - \\boldsymbol{\\omega} \\times r_{\\text{inertial}})

where :math:`R(\\theta)` is the right-handed rotation matrix about
:math:`+z` and :math:`\\boldsymbol{\\omega} = \\omega \\hat z`. The inverse
``from_rotating`` is the algebraic inverse, not a matrix inverse: a left
rotation by :math:`+\\omega t`, then add back :math:`\\boldsymbol{\\omega}
\\times r`.

Scope (M3 — uniform frame only)
-------------------------------

This module implements ONLY the uniform rotating frame. It is **exact** for
the M1 circular-coplanar :class:`~cyclerfinder.core.ephemeris.Ephemeris`
(every planet rides a circle in the ecliptic plane at its mean motion). It
is **NOT** correct for a real ephemeris (eccentric, inclined orbits): in
that regime the spec §12(c) dynamic frame is required — a non-uniform
rotating frame anchored to instantaneous Sun-Earth (or Sun-pair)
geometry, with tolerant ("geometric breathing") closure checks. The
dynamic frame is an M6 deliverable; :func:`synodic_omega` raises
``NotImplementedError`` for Venus as a placeholder reminder that the
Venus-anchored frames belong to the VEM campaign (M8) rather than M3.

References
----------

* Vallado, D. A., *Fundamentals of Astrodynamics and Applications*, 4th
  ed., Microcosm Press, 2013, §3.4 (rotating frames) and §3.5 (Coriolis /
  centripetal cross-term in velocity transformation).
* spec.md §3 (synodic background), §10 (closure-frame correctness risk),
  §12(c) (dynamic ephemeris frame, deferred to M6).
* Plan: ``docs/phases/m3-model-construct/plan.md`` §3.2.
"""

from __future__ import annotations

from math import cos, pi, sin

import numpy as np
from numpy.typing import NDArray

from cyclerfinder.core.constants import PLANETS, SECONDS_PER_DAY

Vec3 = NDArray[np.float64]  # shape (3,), dtype float64


def to_rotating(
    r_inertial: Vec3,
    v_inertial: Vec3,
    t: float,
    omega_rad_per_s: float,
) -> tuple[Vec3, Vec3]:
    """Inertial heliocentric -> uniform rotating frame at angular rate ``omega`` about z.

    Parameters
    ----------
    r_inertial, v_inertial:
        Position (km) and velocity (km/s) in the heliocentric inertial frame.
        Length-3 float64 arrays. The z-axis is the ecliptic normal.
    t:
        Time since the rotating frame's reference epoch (s). The rotating
        frame coincides with the inertial frame at ``t = 0``.
    omega_rad_per_s:
        Rotation rate of the frame about :math:`+\\hat z` (rad/s). A
        positive value rotates the frame's x-axis counter-clockwise.

    Returns
    -------
    (r_rot, v_rot):
        New length-3 float64 arrays (no in-place mutation of the inputs).

    Notes
    -----
    The transform composes a passive rotation of the position by
    :math:`-\\omega t` with the Coriolis-style velocity correction
    :math:`v - \\boldsymbol{\\omega} \\times r`, then rotates the corrected
    velocity by the same angle. See module docstring for derivation.
    """
    theta = omega_rad_per_s * t
    c, s = cos(theta), sin(theta)

    # Passive rotation by -theta (i.e. R(-theta)) about +z.
    # R(-theta) = [[ cos t,  sin t, 0],
    #              [-sin t,  cos t, 0],
    #              [   0,      0,  1]]
    rx, ry, rz = float(r_inertial[0]), float(r_inertial[1]), float(r_inertial[2])
    vx, vy, vz = float(v_inertial[0]), float(v_inertial[1]), float(v_inertial[2])

    # Velocity correction: v - omega x r, with omega = omega_rad_per_s * z_hat.
    # omega x r = (-omega * r_y, +omega * r_x, 0).
    vx_corr = vx - (-omega_rad_per_s * ry)
    vy_corr = vy - (omega_rad_per_s * rx)
    vz_corr = vz  # z-component of omega x r is zero by construction.

    r_rot = np.array(
        [c * rx + s * ry, -s * rx + c * ry, rz],
        dtype=np.float64,
    )
    v_rot = np.array(
        [c * vx_corr + s * vy_corr, -s * vx_corr + c * vy_corr, vz_corr],
        dtype=np.float64,
    )
    return r_rot, v_rot


def from_rotating(
    r_rot: Vec3,
    v_rot: Vec3,
    t: float,
    omega_rad_per_s: float,
) -> tuple[Vec3, Vec3]:
    """Inverse of :func:`to_rotating`.

    Algebraically: rotate the position by ``+omega * t`` to recover the
    inertial position, rotate the velocity likewise to recover the
    ``v - omega x r`` quantity, then add back ``omega x r_inertial`` to
    obtain the inertial velocity.

    Parameters
    ----------
    r_rot, v_rot:
        State in the rotating frame (km, km/s). Length-3 float64.
    t:
        Time since the frame epoch (s).
    omega_rad_per_s:
        Rotation rate used in the forward transform (rad/s).

    Returns
    -------
    (r_inertial, v_inertial):
        Reconstructed inertial state. ``from_rotating(to_rotating(...))``
        is the identity to ~1e-13 relative for inputs at 1 AU scale.
    """
    theta = omega_rad_per_s * t
    c, s = cos(theta), sin(theta)

    rx_r, ry_r, rz_r = float(r_rot[0]), float(r_rot[1]), float(r_rot[2])
    vx_r, vy_r, vz_r = float(v_rot[0]), float(v_rot[1]), float(v_rot[2])

    # R(+theta) about +z: [[cos t, -sin t, 0], [sin t, cos t, 0], [0, 0, 1]].
    rx_i = c * rx_r - s * ry_r
    ry_i = s * rx_r + c * ry_r
    rz_i = rz_r

    # Rotate velocity by +theta, then add back omega x r_inertial.
    vx_corr_i = c * vx_r - s * vy_r
    vy_corr_i = s * vx_r + c * vy_r
    vz_corr_i = vz_r

    # omega x r_inertial = (-omega * r_y, omega * r_x, 0).
    vx_i = vx_corr_i + (-omega_rad_per_s * ry_i)
    vy_i = vy_corr_i + (omega_rad_per_s * rx_i)
    vz_i = vz_corr_i

    r_inertial = np.array([rx_i, ry_i, rz_i], dtype=np.float64)
    v_inertial = np.array([vx_i, vy_i, vz_i], dtype=np.float64)
    return r_inertial, v_inertial


def synodic_omega(body: str) -> float:
    """Angular rate (rad/s) of the synodic rotating frame for an Earth-``body`` pair.

    For Earth-Mars cyclers (the M3 case) this is **Earth's** mean motion: the
    spacecraft's repeating geometry is naturally described in the frame that
    follows Earth's longitude. For Earth-Venus or VEM cyclers the caller chooses
    which body anchors the frame (typically the slowest member of the pair set);
    M3 only needs the E-M case.

    Parameters
    ----------
    body:
        One-letter planet code (``"E"``, ``"M"``, or ``"V"``).

    Returns
    -------
    float
        Angular rate in rad/s.

    Raises
    ------
    NotImplementedError
        For ``"V"``: Venus-anchored frames are part of the VEM campaign (M8),
        and the M3 single-frame primitive should not pretend to support them
        silently.
    KeyError
        For body codes not in :data:`PLANETS`.
    """
    if body == "V":
        raise NotImplementedError(
            "Venus-anchored synodic frame is part of the M8 VEM campaign; "
            "M3 only supports Earth-anchored frames for E-M cyclers."
        )
    if body not in ("E", "M"):
        # Forward to PLANETS lookup so unknown codes surface as KeyError.
        _ = PLANETS[body]
    # Earth's mean motion (deg/day) -> rad/s. Both "E" and "M" use Earth's
    # mean motion because the spacecraft repeats its geometry in Earth's
    # rotating frame for an E-M cycler.
    n_deg_day = PLANETS["E"].mean_motion_deg_day
    return n_deg_day * (pi / 180.0) / SECONDS_PER_DAY
