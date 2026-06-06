"""Synodic rotating-frame transforms (uniform + dynamic).

This module hosts two coexisting rotating-frame transforms:

* **Uniform (M3).** Constant angular rate :math:`\\omega` about the
  heliocentric :math:`+\\hat z` axis. Frame x-axis coincides with the
  inertial x-axis at :math:`t = 0`. Exact for the M1 circular-coplanar
  :class:`~cyclerfinder.core.ephemeris.Ephemeris` (every planet rides a
  circle in the ecliptic plane at its mean motion). M3/M5 consumers
  (closure residual, idealized-mode optimiser) use this transform.

* **Dynamic (M6a).** Non-uniform rotating frame anchored to instantaneous
  Sun-body[0] geometry per spec §12(c). The frame x-axis at time
  :math:`t` aligns with the instantaneous Sun→body[0] direction; the
  angular rate is the instantaneous synodic rate of the body[0]→body[1]
  line. Real-ephemeris verification (M6a binding gate, M6b TCM
  optimiser, M7 V2 batch validation) uses this transform. For a
  circular-coplanar ephemeris it agrees with the uniform frame to
  numerical precision (see :func:`to_rotating_dynamic`'s docstring and
  :func:`tests.core.test_frames_dynamic.test_dynamic_frame_reduces_to_uniform_for_circular`).

The uniform-frame transform is:

.. math::

    r_{\\text{rot}} &= R(-\\omega t) \\, r_{\\text{inertial}} \\\\
    v_{\\text{rot}} &= R(-\\omega t) \\, (v_{\\text{inertial}}
                       - \\boldsymbol{\\omega} \\times r_{\\text{inertial}})

where :math:`R(\\theta)` is the right-handed rotation matrix about
:math:`+z` and :math:`\\boldsymbol{\\omega} = \\omega \\hat z`. The
inverse ``from_rotating`` is the algebraic inverse, not a matrix inverse:
a left rotation by :math:`+\\omega t`, then add back
:math:`\\boldsymbol{\\omega} \\times r`.

The dynamic-frame transform replaces the constant
:math:`\\theta = \\omega t` with the instantaneous frame angle
:math:`\\theta(t) = \\mathrm{atan2}(r_{b_0}(t)_y, r_{b_0}(t)_x)` read
directly from the ephemeris (not integrated from :math:`\\omega(t)`).
Reading :math:`\\theta(t)` directly makes the round-trip identity
algebraically exact — per spec §10's binding correctness risk, a wrong
frame would silently fake/break every drift measurement. The Coriolis
velocity correction uses the instantaneous :math:`\\omega(t)` from
:func:`synodic_omega_dynamic`.

Frame anchor (dynamic frame)
----------------------------
``bodies[0]`` defines the +x direction; ``bodies[1]`` enters only via
:func:`synodic_omega_dynamic` for the Coriolis-correction angular rate.
For ``len(bodies) >= 3`` the additional bodies are ignored at the
frame-definition level; multi-body VEM frames belong to M8.

Scope
-----

* :func:`to_rotating` / :func:`from_rotating` / :func:`synodic_omega` —
  M3 uniform-frame surface. Bit-stable; M3/M5 consumers must keep getting
  identical results.
* :func:`to_rotating_dynamic` / :func:`from_rotating_dynamic` /
  :func:`synodic_omega_dynamic` — M6a dynamic-frame extension. Spec
  §12(c) compliant. Round-trip identity to ≤ 1e-10 rel (the spec §10
  binding unit test in ``tests/core/test_frames_dynamic.py``).

References
----------

* Vallado, D. A., *Fundamentals of Astrodynamics and Applications*, 4th
  ed., Microcosm Press, 2013, §3.4 (rotating frames) and §3.5 (Coriolis /
  centripetal cross-term in velocity transformation).
* spec.md §3 (synodic background), §10 (closure-frame correctness risk),
  §12(c) (dynamic ephemeris frame).
* Plans: ``docs/phases/m3-model-construct/plan.md`` §3.2 (uniform),
  ``docs/phases/m6a-idealized-closure-verification/plan.md`` §3.1
  (dynamic extension).
"""

from __future__ import annotations

from math import atan2, cos, pi, sin

import numpy as np
from numpy.typing import NDArray

from cyclerfinder.core.constants import PLANETS, SECONDS_PER_DAY
from cyclerfinder.core.ephemeris import Ephemeris

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


def _rodrigues_rotate(vec: Vec3, axis_hat: Vec3, angle: float) -> Vec3:
    """Rotate ``vec`` about the unit ``axis_hat`` by ``angle`` (rad), Rodrigues.

    Right-handed rotation: ``v cos a + (k x v) sin a + k (k . v)(1 - cos a)``
    for unit axis ``k``. Module-internal helper for the vector-omega frame.
    """
    c, s = cos(angle), sin(angle)
    kx, ky, kz = float(axis_hat[0]), float(axis_hat[1]), float(axis_hat[2])
    vx, vy, vz = float(vec[0]), float(vec[1]), float(vec[2])
    # k x v
    cross_x = ky * vz - kz * vy
    cross_y = kz * vx - kx * vz
    cross_z = kx * vy - ky * vx
    # k . v
    dot = kx * vx + ky * vy + kz * vz
    return np.array(
        [
            vx * c + cross_x * s + kx * dot * (1.0 - c),
            vy * c + cross_y * s + ky * dot * (1.0 - c),
            vz * c + cross_z * s + kz * dot * (1.0 - c),
        ],
        dtype=np.float64,
    )


def to_rotating_omega_vec(
    r_inertial: Vec3,
    v_inertial: Vec3,
    t: float,
    omega_vec: Vec3,
) -> tuple[Vec3, Vec3]:
    """Inertial heliocentric -> uniform rotating frame at vector rate ``omega_vec``.

    A strict superset of :func:`to_rotating`: the scalar form's ``omega * z_hat``
    is replaced by an explicit angular-velocity vector ``omega_vec`` so the frame
    spins about the unit axis ``omega_hat`` at rate ``|omega_vec|`` and the
    Coriolis correction ``v - omega_vec x r`` is a full 3-D cross product.

    Coplanar-limit gate (binding, plan §1): when
    ``omega_vec = (0, 0, omega)`` this returns arrays that are
    ``numpy.array_equal`` to ``to_rotating(r, v, t, omega)`` bit-for-bit. That
    is guaranteed by delegating the pure-z case directly to the scalar
    :func:`to_rotating`; only a genuinely tilted ``omega_vec`` takes the general
    Rodrigues path.

    Parameters
    ----------
    r_inertial, v_inertial:
        Inertial heliocentric state, km and km/s. Length-3 float64.
    t:
        Time since the rotating frame's reference epoch (s).
    omega_vec:
        Angular-velocity vector of the frame (rad/s), length-3. The frame
        rotates about ``omega_vec / |omega_vec|`` at rate ``|omega_vec|``.

    Returns
    -------
    (r_rot, v_rot):
        Two length-3 float64 arrays (no in-place mutation of the inputs).
    """
    wx, wy, wz = float(omega_vec[0]), float(omega_vec[1]), float(omega_vec[2])
    # Pure-z omega: delegate to the scalar form so the coplanar-limit gate is
    # bit-for-bit (numpy.array_equal), not merely close.
    if wx == 0.0 and wy == 0.0:
        return to_rotating(r_inertial, v_inertial, t, wz)

    omega_mag = (wx * wx + wy * wy + wz * wz) ** 0.5
    omega_hat = np.array([wx / omega_mag, wy / omega_mag, wz / omega_mag], dtype=np.float64)
    theta = omega_mag * t

    # Coriolis correction: v - omega_vec x r.
    rx, ry, rz = float(r_inertial[0]), float(r_inertial[1]), float(r_inertial[2])
    vx, vy, vz = float(v_inertial[0]), float(v_inertial[1]), float(v_inertial[2])
    wxr_x = wy * rz - wz * ry
    wxr_y = wz * rx - wx * rz
    wxr_z = wx * ry - wy * rx
    v_corr = np.array([vx - wxr_x, vy - wxr_y, vz - wxr_z], dtype=np.float64)

    # Passive rotation by -theta about omega_hat.
    r_rot = _rodrigues_rotate(np.asarray(r_inertial, dtype=np.float64), omega_hat, -theta)
    v_rot = _rodrigues_rotate(v_corr, omega_hat, -theta)
    return r_rot, v_rot


def synodic_omega(anchor_body: str) -> float:
    """Angular rate (rad/s) of the synodic rotating frame anchored on ``anchor_body``.

    The synodic frame follows the longitude of a chosen *anchor* body, so its
    angular rate is that body's mean motion. For an Earth-Mars cycler the
    natural anchor is **Earth** (``synodic_omega("E")``): the spacecraft's
    geometry repeats in the frame that co-rotates with Earth's longitude. For
    Earth-Venus or VEM cyclers the caller passes whichever body anchors the
    frame (commonly the cycler's home/first-encounter body).

    Parameters
    ----------
    anchor_body:
        One-letter planet code in :data:`PLANETS` (e.g. ``"E"``, ``"M"``,
        ``"V"``). The frame rate is this body's mean motion.

    Returns
    -------
    float
        Angular rate in rad/s.

    Raises
    ------
    KeyError
        For body codes not in :data:`PLANETS`.
    """
    n_deg_day = PLANETS[anchor_body].mean_motion_deg_day
    return n_deg_day * (pi / 180.0) / SECONDS_PER_DAY


# ---------------------------------------------------------------------------
# Dynamic frame (M6a, spec §12(c))
# ---------------------------------------------------------------------------


def synodic_omega_dynamic(
    t_sec: float,
    bodies: tuple[str, ...],
    ephem: Ephemeris,
) -> float:
    """Instantaneous synodic angular rate (rad/s) for the body pair anchor.

    For ``len(bodies) >= 2`` returns the instantaneous rate at which the
    Sun→body[0] line is rotating in the inertial frame, computed from
    ``bodies[0]``'s heliocentric :math:`(r, v)`. This is the
    **non-uniform** synodic rate per spec §12(c) — the value depends on
    ``t_sec`` because real orbits are eccentric.

    The angular rate is obtained from the planar formula
    :math:`\\omega = (r \\times v)_z / |r|^2`, evaluated on the
    body's heliocentric state at ``t_sec``. The :math:`z` component
    picks the prograde-positive sign convention.

    For ``bodies[0]`` on a circular orbit (the M1
    ``Ephemeris("circular")`` backend) this equals
    :func:`synodic_omega(bodies[0]) <synodic_omega>` to floating-point
    precision at every ``t_sec``. M6a unit-tests this degeneracy in
    ``test_synodic_omega_dynamic_matches_uniform_for_circular_earth``.

    Parameters
    ----------
    t_sec:
        Seconds since the ephemeris reference epoch.
    bodies:
        Body codes; ``bodies[0]`` defines the frame anchor. Length ≥ 2.
        Lengths > 2 use only ``bodies[0]`` for the frame angle.
    ephem:
        Heliocentric state provider.

    Returns
    -------
    Instantaneous frame rotation rate (rad/s), signed positive for
    prograde.

    Raises
    ------
    ValueError
        If ``len(bodies) < 2``.
    """
    if len(bodies) < 2:
        raise ValueError(
            f"synodic_omega_dynamic requires len(bodies) >= 2; got {bodies!r}",
        )
    r, v = ephem.state(bodies[0], t_sec)
    # In-plane angular rate: (r cross v)_z / |r|^2. Positive for
    # prograde motion in the J2000 ecliptic (the convention used
    # throughout cyclerfinder).
    rxv_z = float(r[0] * v[1] - r[1] * v[0])
    r_sq = float(r[0] * r[0] + r[1] * r[1] + r[2] * r[2])
    if r_sq == 0.0:
        raise ValueError(
            f"synodic_omega_dynamic: bodies[0]={bodies[0]!r} has zero "
            f"heliocentric radius at t_sec={t_sec}",
        )
    return rxv_z / r_sq


def _frame_angle_dynamic(
    t_sec: float,
    bodies: tuple[str, ...],
    ephem: Ephemeris,
) -> float:
    """Dynamic-frame angle ``theta(t) = atan2(r_b0_y, r_b0_x)``, rad.

    Read directly from the ephemeris rather than integrated from
    :func:`synodic_omega_dynamic`. Per plan §3.1.2 this makes the
    round-trip identity algebraically exact (the forward and inverse
    transforms evaluate the same ``atan2`` and recover the inputs to
    float64 precision).

    Module-internal helper.
    """
    if len(bodies) < 2:
        raise ValueError(
            f"_frame_angle_dynamic requires len(bodies) >= 2; got {bodies!r}",
        )
    r, _ = ephem.state(bodies[0], t_sec)
    return atan2(float(r[1]), float(r[0]))


def to_rotating_dynamic(
    r_inertial: Vec3,
    v_inertial: Vec3,
    t_sec: float,
    bodies: tuple[str, ...],
    ephem: Ephemeris,
) -> tuple[Vec3, Vec3]:
    """Inertial heliocentric → dynamic rotating frame anchored to Sun→body[0].

    Per spec §12(c) the rotating frame at time ``t_sec`` has its x-axis
    aligned with the instantaneous Sun→body[0] direction; the angular
    rate is the instantaneous synodic rate from
    :func:`synodic_omega_dynamic`.

    The frame angle ``θ(t) = atan2(r_b0(t)_y, r_b0(t)_x)`` is read from
    the ephemeris (not integrated), making the round-trip identity
    exact rather than just numerically close (plan §3.1.2).

    The velocity transform applies the instantaneous Coriolis
    correction ``v - omega(t) x r`` before rotating; this matches the
    M3 uniform frame's structure but uses ``omega(t_sec)`` rather than
    a constant.

    Parameters
    ----------
    r_inertial, v_inertial:
        Inertial heliocentric state, km and km/s. Length-3 float64.
    t_sec:
        Seconds since the ephemeris reference epoch.
    bodies:
        Body codes; ``bodies[0]`` defines the frame anchor. Length ≥ 2.
    ephem:
        Heliocentric state provider.

    Returns
    -------
    (r_rot, v_rot):
        Two length-3 float64 arrays (km, km/s). See
        :func:`from_rotating_dynamic` for the inverse.

    Raises
    ------
    ValueError
        If ``len(bodies) < 2`` or body[0]'s heliocentric radius vanishes.
    """
    theta = _frame_angle_dynamic(t_sec, bodies, ephem)
    omega = synodic_omega_dynamic(t_sec, bodies, ephem)
    c, s = cos(theta), sin(theta)

    rx, ry, rz = float(r_inertial[0]), float(r_inertial[1]), float(r_inertial[2])
    vx, vy, vz = float(v_inertial[0]), float(v_inertial[1]), float(v_inertial[2])

    # Coriolis correction: v - omega(t) x r, with omega along +z.
    vx_corr = vx - (-omega * ry)
    vy_corr = vy - (omega * rx)
    vz_corr = vz

    # Passive rotation by -theta about +z.
    r_rot = np.array(
        [c * rx + s * ry, -s * rx + c * ry, rz],
        dtype=np.float64,
    )
    v_rot = np.array(
        [c * vx_corr + s * vy_corr, -s * vx_corr + c * vy_corr, vz_corr],
        dtype=np.float64,
    )
    return r_rot, v_rot


def from_rotating_dynamic(
    r_rot: Vec3,
    v_rot: Vec3,
    t_sec: float,
    bodies: tuple[str, ...],
    ephem: Ephemeris,
) -> tuple[Vec3, Vec3]:
    """Inverse of :func:`to_rotating_dynamic`.

    Exact inverse because ``θ(t)`` is read from the ephemeris (no
    numerical integration). Round-trip identity:
    ``from_rotating_dynamic(*to_rotating_dynamic(r, v, t, bodies, ephem),
    t, bodies, ephem) == (r, v)`` to ≤ 1e-10 relative on AU-scale
    states. The M6a gate test
    ``test_dynamic_frame_roundtrip_identity`` asserts this.

    Parameters
    ----------
    r_rot, v_rot:
        Rotating-frame state, km and km/s. Length-3 float64.
    t_sec:
        Seconds since the ephemeris reference epoch (the same value
        used in the forward transform).
    bodies:
        Body codes; ``bodies[0]`` defines the frame anchor. Length ≥ 2.
    ephem:
        Heliocentric state provider.

    Returns
    -------
    (r_inertial, v_inertial):
        Reconstructed inertial state, km and km/s.
    """
    theta = _frame_angle_dynamic(t_sec, bodies, ephem)
    omega = synodic_omega_dynamic(t_sec, bodies, ephem)
    c, s = cos(theta), sin(theta)

    rx_r, ry_r, rz_r = float(r_rot[0]), float(r_rot[1]), float(r_rot[2])
    vx_r, vy_r, vz_r = float(v_rot[0]), float(v_rot[1]), float(v_rot[2])

    # R(+theta) about +z.
    rx_i = c * rx_r - s * ry_r
    ry_i = s * rx_r + c * ry_r
    rz_i = rz_r

    vx_corr_i = c * vx_r - s * vy_r
    vy_corr_i = s * vx_r + c * vy_r
    vz_corr_i = vz_r

    # Add back omega(t) x r_inertial.
    vx_i = vx_corr_i + (-omega * ry_i)
    vy_i = vy_corr_i + (omega * rx_i)
    vz_i = vz_corr_i

    r_inertial = np.array([rx_i, ry_i, rz_i], dtype=np.float64)
    v_inertial = np.array([vx_i, vy_i, vz_i], dtype=np.float64)
    return r_inertial, v_inertial
