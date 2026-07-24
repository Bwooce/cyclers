"""CCR4BP-to-real-ephemeris consistency check (`#704`).

`#701` (commits ``2aa13e6``/``727eccb``, `#702` fix ``2d8d563``) found a
ghost-guard-verified homoclinic connection of a resonant quasi-periodic torus
in the Uranus-Umbriel-Titania CCR4BP -- an IDEALIZED model: Umbriel and
Titania on CONCENTRIC CIRCULAR COPLANAR orbits, real physical masses. Real
Umbriel (e=0.0041) and Titania (e=0.0021) are both eccentric and mutually
inclined by a fraction of a degree to Uranus's equatorial (Laplace) plane --
neither is modelled by the idealized CCR4BP at all.

This module asks the load-bearing vetting-chain question `#704` was
dispatched to answer: does `#701`'s idealized connection approximately
survive when Umbriel's and Titania's REAL (SPICE-sourced, URA111) ephemeris
replaces the circular-coplanar approximation?  No prior real-ephemeris
validation capability for a torus/manifold object exists anywhere in this
codebase (Uranian or otherwise) -- `#312`'s own V0-V5 gauntlet
(``data/validation/v0_uranus.py``...``v4_uranus_strict.py``) is built
entirely around discrete Lambert-arc point-to-point tour closures, a
structurally different object from a continuous-flow-time manifold
connection. This module is new, additive capability, built specifically for
that different object shape.

Method
------
1. Take `#701`'s own converged unstable-manifold DEPARTURE state (before any
   flow -- the torus point perturbed by ``eps`` along the extracted CLV
   eigenvector, in the idealized CCR4BP's own Uranus-Umbriel synodic
   rotating frame, nondimensional units) and its own converged
   stable-manifold TARGET state (``state_s``, the point the idealized
   connection's own residual drives to near-coincide with the unstable
   arrival -- see :mod:`cyclerfinder.search.ccr4bp_heteroclinic_search`).
2. Convert the departure state to a REAL physical state (position, velocity)
   at a chosen real epoch, using an INSTANTANEOUS "osculating rotating
   frame" built from Umbriel's REAL SPICE state at that epoch (real
   x-axis = toward Umbriel, real z-axis = Umbriel's own instantaneous orbit
   normal, real angular rate = Umbriel's own instantaneous ``h/r^2`` --
   :func:`osculating_frame`). This is the natural real-ephemeris analogue of
   "Umbriel fixed on the rotating frame's x-axis," evaluated pointwise
   rather than assumed constant (Umbriel's real angular rate is NOT
   constant -- e=0.0041 -- so this differs at every epoch and phase).
3. Propagate that real physical state FORWARD for the SAME elapsed flow
   time `#701`'s connection used (``t_u``, converted from nondimensional
   time units to seconds), under a real force model: Uranus point-mass
   central term + Umbriel's and Titania's REAL SPICE-sourced third-body
   perturbations (:func:`real_nbody_rhs`) in the Uranus-centred J2000
   inertial frame (NOT the idealized rotating frame, whose own definition
   assumes circular coplanar moons and therefore cannot represent the real
   trajectory).
4. Convert the idealized connection's OWN target state (``state_s``) to the
   SAME real/physical frame at the CORRESPONDING real epoch (``epoch0 +
   t_u`` -- the epoch the real propagation actually reaches), using the
   SAME osculating-frame construction (now built from Umbriel's real state
   at that later epoch).
5. Compare the propagated real endpoint to the converted target: the
   position/velocity gap is the headline number
   (:class:`ConsistencyCheckResult`).

Sourcing / reuse discipline
----------------------------
* SPICE kernel paths and the FURNSH/spkezr mechanics are imported UNMODIFIED
  from :mod:`cyclerfinder.data.validation.v4_uranus_strict` (`#312`'s own V4
  gauntlet) -- the exact same kernel-loading mechanism, not re-derived
  (``_spice_furnsh_all``, ``_moon_state_spice``, ``_ephemeris_time_seconds``,
  ``DEFAULT_LSK_PATH``/``DEFAULT_PCK_PATH``/``DEFAULT_URA_PATH``). Reusing a
  module's private helpers across files is an already-precedented pattern in
  this arc (see :mod:`cyclerfinder.search.ccr4bp_heteroclinic_search`'s own
  docstring on why, and its own re-use of ``hs._L_KM`` from test code).
* The third-body Battin acceleration formula is imported UNMODIFIED from
  :mod:`cyclerfinder.data.validation.v4_uranus` (``_third_body_acceleration_kms2``),
  the SAME formula `#312`'s own V4/V4-strict gauntlets use.
* :mod:`cyclerfinder.core.ccr4bp`, :mod:`cyclerfinder.search.variational_ccr4bp_torus`,
  :mod:`cyclerfinder.search.ccr4bp_whisker`,
  :mod:`cyclerfinder.search.ccr4bp_manifold_globalize`, and
  :mod:`cyclerfinder.search.ccr4bp_heteroclinic_search` are all consumed
  read-only through their existing public API -- this module does NOT modify
  any of `#689`-`#694`'s or `#701`'s code.

Force-model scope and honest limitations
-----------------------------------------
* **No Uranus J2.** The idealized CCR4BP itself has no J2 term either -- to
  isolate exactly the delta this task asks about (circular-coplanar ->
  real-ephemeris moon positions), the real force model here changes ONLY
  the moon-position source, not the primary's gravity model. J2 is a real,
  separate effect (`#332`'s V4 fallback quantifies it at ~5e-5 of the
  central acceleration at Umbriel's SMA) that is NOT included and would need
  a separate sensitivity pass to characterize on top of this result.
* **No Sun.** Checked, not assumed: at Uranus's ~19 AU heliocentric
  distance, the Sun's tidal acceleration on a spacecraft ~3e5 km from
  Uranus is ~3.4e-12 km/s^2 (``2*GM_sun*r_sc/d_sun^3``), roughly 5e-8 of the
  central Uranus acceleration at that distance (~6.4e-5 km/s^2) -- about
  three orders of magnitude below the ALREADY-negligible J2 term `#332`'s V4
  fallback documented and omitted. Negligible on this connection's ~13-day
  timescale.
* **No other Uranian moons** (Miranda, Ariel, Oberon) as perturbers -- the
  idealized CCR4BP itself has only Umbriel (base) and Titania (perturber);
  adding the others would test a DIFFERENT question (whether the whole
  5-body real system is self-consistent) than this task's own scope (does
  THIS idealized 3-body-plus-forcing model's own connection survive real
  Umbriel/Titania ephemeris).
* **Uranus-centred propagation frame, but the barycentre offset IS
  corrected, not neglected.** The idealized CCR4BP's own origin is the
  Uranus-Umbriel barycentre (Uranus fixed at nondim ``(-mu, 0, 0)``); the
  physical offset from Uranus's own centre is tiny (``mu * L_KM ~ 3.9 km``
  for ``mu ~ 1.469e-5``), but this system's own strong hyperbolicity
  amplifies even that few-km offset by roughly an order of magnitude over a
  several-TU flow -- confirmed empirically (this module's own reduction-test
  positive control, :func:`nondim_state_to_inertial`'s docstring). Every
  nondim<->inertial conversion therefore takes ``mu`` explicitly and
  re-origins to Uranus's own centre before scaling/rotating -- NOT treated
  as an interchangeable-with-barycentre approximation.
* **Point-mass moons, no patched-conic softening inside a Hill sphere.**
  Unlike `#312`'s V4 gauntlet (a between-encounter propagator that never
  flies THROUGH a moon's Hill sphere by Lambert-geometry construction), this
  connection's own trajectory passes near Titania's orbital radius by
  design (`#701`'s own base-orbit note). :func:`propagate_real` reports the
  closest approach to each perturber during the propagated arc
  (``closest_approach_km``) so a caller can judge whether the point-mass
  model broke down at any point, rather than silently trusting a possibly
  unphysical close pass.
* **Single chosen epoch + a discrete epoch scan**, not a continuous duty
  cycle. This is a NEW methodology (unlike `#312`'s own already-validated
  Lambert-arc gauntlet), so the scan resolution is a pragmatic, stated
  choice, not a proof of the full epoch-dependence structure.

No catalogue writeback; this module produces a CONSISTENCY-CHECK result, not
a vetted discovery, and does not implement the (separate, later) schema/
writeback tasks.
"""

from __future__ import annotations

import math
from collections.abc import Callable
from dataclasses import dataclass, field

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import solve_ivp

from cyclerfinder.core.satellites import PRIMARIES, SATELLITES
from cyclerfinder.data.validation.v4_uranus import _third_body_acceleration_kms2
from cyclerfinder.data.validation.v4_uranus_strict import (
    DEFAULT_LSK_PATH,
    DEFAULT_PCK_PATH,
    DEFAULT_URA_PATH,
    _ephemeris_time_seconds,
    _moon_state_spice,
    _spice_furnsh_all,
)

__all__ = [
    "DEFAULT_LSK_PATH",
    "DEFAULT_PCK_PATH",
    "DEFAULT_URA_PATH",
    "PRIMARIES",
    "ConsistencyCheckResult",
    "MoonStateFn",
    "check_connection_survives_real_ephemeris",
    "idealized_moon_state_fn",
    "inertial_state_to_nondim",
    "nondim_state_to_inertial",
    "osculating_frame",
    "propagate_real",
    "real_nbody_rhs",
    "spice_moon_state_fn",
    "tu_to_seconds",
]

#: Signature every moon-state source (real SPICE or idealized circular
#: substitute) must satisfy: ``(moon_name, et_seconds) -> (r_km, v_km_s)``,
#: Uranus-centred, in whatever frame the caller is working in (J2000 for the
#: real path, the same idealized-inertial frame for the reduction-test
#: substitute -- see :func:`idealized_moon_state_fn`).
MoonStateFn = Callable[[str, float], tuple[NDArray[np.float64], NDArray[np.float64]]]


def spice_moon_state_fn(
    moon_name: str, et_seconds: float
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Real SPICE-sourced moon state (thin wrapper on the reused `#312` helper).

    Requires the caller to have already FURNSH'd the kernels (see
    :func:`check_connection_survives_real_ephemeris`, which does this in a
    try/finally around the whole check -- matching `#312`'s own
    ``run_v4_uranus_strict`` discipline).
    """
    return _moon_state_spice(moon_name, et_seconds)


# --------------------------------------------------------------------------- #
# Rotating <-> inertial frame conversion
# --------------------------------------------------------------------------- #


def osculating_frame(
    r_km: NDArray[np.float64], v_km_s: NDArray[np.float64]
) -> tuple[NDArray[np.float64], float]:
    """Instantaneous "osculating rotating frame" built from a body's state.

    Returns ``(R, omega_inst)`` where ``R`` is the 3x3 rotation matrix whose
    COLUMNS are the rotating frame's ``(x_hat, y_hat, z_hat)`` basis vectors
    expressed in the inertial frame (``x_hat`` toward the body, ``z_hat``
    along its instantaneous orbit-normal, ``y_hat = z_hat x x_hat``), and
    ``omega_inst`` is the body's own instantaneous angular rate about that
    ``z_hat`` axis, rad/s (``|h| / |r|^2`` -- exact for planar motion at any
    eccentricity, from the standard areal-velocity relation ``h = r^2
    dtheta/dt``; NOT assumed constant, unlike a circular-orbit mean motion).

    This is the "as-if-Umbriel-defines-the-rotating-frame" construction,
    evaluated pointwise at a single instant from the body's REAL state
    rather than assuming a fixed circular rate -- the natural real-ephemeris
    analogue of the idealized CCR4BP's own rotating frame.
    """
    r = np.asarray(r_km, dtype=np.float64)
    v = np.asarray(v_km_s, dtype=np.float64)
    r_norm = float(np.linalg.norm(r))
    if r_norm <= 0.0:
        raise ValueError("osculating_frame: zero-norm position vector")
    x_hat = r / r_norm
    h = np.cross(r, v)
    h_norm = float(np.linalg.norm(h))
    if h_norm <= 0.0:
        raise ValueError("osculating_frame: zero angular momentum (degenerate orbit)")
    z_hat = h / h_norm
    y_hat = np.cross(z_hat, x_hat)
    rot = np.column_stack([x_hat, y_hat, z_hat])
    omega_inst = h_norm / (r_norm * r_norm)
    return rot, float(omega_inst)


def tu_to_seconds(l_km: float, v_unit_km_s: float) -> float:
    """Nondimensional CCR4BP time unit (1 TU = 1/n1) in seconds.

    ``n1 = v_unit_km_s / l_km`` (the base Uranus-Umbriel mean motion the
    CCR4BP nondimensionalises time by -- see
    :mod:`cyclerfinder.core.ccr4bp_umbriel_titania`), so ``1 TU = l_km /
    v_unit_km_s`` seconds.
    """
    return l_km / v_unit_km_s


def nondim_state_to_inertial(
    state6_nondim: NDArray[np.float64],
    r_moon_km: NDArray[np.float64],
    v_moon_km_s: NDArray[np.float64],
    l_km: float,
    v_unit_km_s: float,
    *,
    mu: float = 0.0,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Convert an idealized CCR4BP rotating-frame nondim 6-state to a real
    inertial, URANUS-CENTRED (not barycentre-centred) physical state, using
    the osculating frame built from ``(r_moon_km, v_moon_km_s)`` (the base
    moon's -- Umbriel's -- real state at the target epoch; see
    :func:`osculating_frame`).

    ``state6_nondim = (x, y, z, vx, vy, vz)`` in Uranus-Umbriel synodic
    nondim units (length = ``l_km``, velocity = ``v_unit_km_s``), where the
    ORIGIN is the Uranus-Umbriel barycentre by the CR3BP/CCR4BP's own
    convention (Uranus fixed at nondim ``(-mu, 0, 0)``, Umbriel at
    ``(1-mu, 0, 0)``, both frame-fixed). ``mu`` re-origins this to
    URANUS-CENTRED before scaling/rotating (``mu=0.0`` recovers the naive
    barycentre~=Uranus approximation; pass the system's own ``mu`` -- e.g.
    ``1.469e-5`` for Umbriel-Titania -- for an exact reduction, which
    matters here because this system's own strong hyperbolicity can amplify
    even the tiny (``~mu*l_km``, a few km) barycentre offset by orders of
    magnitude over a several-TU flow -- confirmed empirically by this
    module's own reduction-test positive control, which is why ``mu`` is a
    real, non-default parameter rather than a documented-away approximation).
    Velocity needs NO offset correction: both primaries are FIXED in the
    rotating frame by construction, so a constant position shift has zero
    time derivative there.

    Returns ``(r_km, v_km_s)``.
    """
    rot, omega_inst = osculating_frame(r_moon_km, v_moon_km_s)
    r_rot_km = l_km * (np.asarray(state6_nondim[:3], dtype=np.float64) + np.array([mu, 0.0, 0.0]))
    v_rot_km_s = v_unit_km_s * np.asarray(state6_nondim[3:], dtype=np.float64)
    omega_cross_r = np.array(
        [-omega_inst * r_rot_km[1], omega_inst * r_rot_km[0], 0.0], dtype=np.float64
    )
    r_inertial = rot @ r_rot_km
    v_inertial = rot @ (v_rot_km_s + omega_cross_r)
    return r_inertial, v_inertial


def inertial_state_to_nondim(
    r_km: NDArray[np.float64],
    v_km_s: NDArray[np.float64],
    r_moon_km: NDArray[np.float64],
    v_moon_km_s: NDArray[np.float64],
    l_km: float,
    v_unit_km_s: float,
    *,
    mu: float = 0.0,
) -> NDArray[np.float64]:
    """Inverse of :func:`nondim_state_to_inertial` (used only for the
    round-trip regression test)."""
    rot, omega_inst = osculating_frame(r_moon_km, v_moon_km_s)
    r_rot_km = rot.T @ np.asarray(r_km, dtype=np.float64)
    v_rot_km_s_total = rot.T @ np.asarray(v_km_s, dtype=np.float64)
    omega_cross_r = np.array(
        [-omega_inst * r_rot_km[1], omega_inst * r_rot_km[0], 0.0], dtype=np.float64
    )
    v_rot_km_s = v_rot_km_s_total - omega_cross_r
    state6 = np.concatenate([r_rot_km / l_km - np.array([mu, 0.0, 0.0]), v_rot_km_s / v_unit_km_s])
    return np.asarray(state6, dtype=np.float64)


# --------------------------------------------------------------------------- #
# Real (or idealized-substitute) N-body force model
# --------------------------------------------------------------------------- #


def real_nbody_rhs(
    t_s: float,
    y: NDArray[np.float64],
    *,
    mu_uranus: float,
    et0: float,
    perturber_moons: tuple[str, ...],
    moon_state_fn: MoonStateFn,
) -> NDArray[np.float64]:
    """Uranus-centred inertial RHS: central point-mass + third-body
    perturbers (Umbriel, Titania by default), no J2, no Sun (see module
    docstring for the quantified justification). ``moon_state_fn`` supplies
    each perturber's Uranus-centred state at absolute ephemeris time ``et0 +
    t_s`` -- either :func:`spice_moon_state_fn` (real) or
    :func:`idealized_moon_state_fn` (the reduction-test substitute).
    """
    r_sc = y[:3]
    v_sc = y[3:]
    r_norm = float(np.linalg.norm(r_sc))
    a_central = -mu_uranus * r_sc / (r_norm**3) if r_norm > 0.0 else np.zeros(3)
    a_total = np.array(a_central, dtype=np.float64)
    et = et0 + t_s
    for moon in perturber_moons:
        r_moon, _v_moon = moon_state_fn(moon, et)
        a_total += _third_body_acceleration_kms2(r_sc, r_moon, mu_body=SATELLITES[moon].mu_km3_s2)
    return np.concatenate([v_sc, a_total])


@dataclass(frozen=True)
class RealPropagationResult:
    """Result of :func:`propagate_real`."""

    r_f_km: NDArray[np.float64]
    v_f_km_s: NDArray[np.float64]
    success: bool
    closest_approach_km: dict[str, float]
    """Minimum spacecraft-perturber distance (km) achieved anywhere along
    the propagated arc, per perturber name -- a diagnostic for whether the
    point-mass model's own Hill-sphere assumption plausibly held (see module
    docstring)."""


def propagate_real(
    r0_km: NDArray[np.float64],
    v0_km_s: NDArray[np.float64],
    tof_s: float,
    *,
    et0: float,
    mu_uranus: float,
    perturber_moons: tuple[str, ...],
    moon_state_fn: MoonStateFn,
    rtol: float = 1e-12,
    atol: float = 1e-6,
    n_diag_samples: int = 400,
) -> RealPropagationResult:
    """Propagate a Uranus-centred inertial state forward under
    :func:`real_nbody_rhs` for ``tof_s`` seconds starting at absolute
    ephemeris time ``et0``."""
    y0 = np.concatenate(
        [np.asarray(r0_km, dtype=np.float64), np.asarray(v0_km_s, dtype=np.float64)]
    )

    def _rhs(t_s: float, y: NDArray[np.float64]) -> NDArray[np.float64]:
        return real_nbody_rhs(
            t_s,
            y,
            mu_uranus=mu_uranus,
            et0=et0,
            perturber_moons=perturber_moons,
            moon_state_fn=moon_state_fn,
        )

    sol = solve_ivp(
        _rhs,
        (0.0, float(tof_s)),
        y0,
        method="DOP853",
        rtol=rtol,
        atol=atol,
        dense_output=True,
    )
    closest: dict[str, float] = dict.fromkeys(perturber_moons, float("inf"))
    if sol.success and sol.sol is not None:
        ts = np.linspace(0.0, float(tof_s), n_diag_samples)
        for t in ts:
            y = sol.sol(t)
            r_sc = y[:3]
            for moon in perturber_moons:
                r_moon, _ = moon_state_fn(moon, et0 + t)
                d = float(np.linalg.norm(r_moon - r_sc))
                if d < closest[moon]:
                    closest[moon] = d
    if not sol.success:
        return RealPropagationResult(
            r_f_km=np.full(3, np.nan),
            v_f_km_s=np.full(3, np.nan),
            success=False,
            closest_approach_km=closest,
        )
    yf = sol.y[:, -1]
    return RealPropagationResult(
        r_f_km=yf[:3], v_f_km_s=yf[3:], success=True, closest_approach_km=closest
    )


def idealized_moon_state_fn(
    system: object,
    l_km: float,
    v_unit_km_s: float,
    mu: float,
    theta_gan0: float,
    a_gan: float,
    omega_gan: float,
) -> MoonStateFn:
    """Build a :data:`MoonStateFn` from the IDEALIZED circular-coplanar
    CCR4BP model -- used ONLY as the reduction-test substitute (see module
    docstring / the test module's positive control): feeding this in place
    of :func:`spice_moon_state_fn` must make :func:`propagate_real` +
    :func:`nondim_state_to_inertial` reduce EXACTLY (integrator tolerance)
    to a direct idealized-rotating-frame :func:`cyclerfinder.core.ccr4bp.propagate_ccr4bp`
    call transformed to inertial via the SAME frame construction.

    Returns positions URANUS-CENTRED (matching :func:`spice_moon_state_fn`'s
    own ``observer="URANUS"`` convention, NOT the CCR4BP's own
    barycentre-centred nondim convention where Umbriel sits at ``(1-mu, 0,
    0)``): Umbriel is exactly ``l_km`` from Uranus (the primary-separation
    IS the length unit, by definition -- the ``(1-mu)`` factor is a
    barycentre-relative radius, not a Uranus-relative one). Titania's own
    barycentre-relative circular position (radius ``a_gan`` at rotating-frame
    angle ``theta_gan0 + omega_gan*t_nondim``) is shifted by the SAME
    ``+mu`` barycentre-to-Uranus offset :func:`nondim_state_to_inertial`
    applies to the spacecraft state, so both bodies are consistently
    Uranus-centred. (An earlier version of this function returned
    barycentric positions unshifted -- caught by this module's own
    reduction-test positive control, which is EXACTLY the TDD check this
    bug-class exists to catch: a ``~mu``-relative, ``~1e-10 km/s^2``
    systematic acceleration mismatch that integrates to a multi-km position
    error after a few TU, quadratic in elapsed time as expected for a
    constant spurious acceleration offset.)

    Both Umbriel and Titania move on circles about Uranus in an inertial
    frame that itself rotates at the frame's own constant rate ``n1 =
    v_unit_km_s/l_km`` -- i.e. Umbriel's inertial angle is ``n1*t_s`` and
    Titania's is ``n1*t_s + (theta_gan0+omega_gan*n1*t_s)``.
    """
    n1 = v_unit_km_s / l_km  # rad/s, base Uranus-Umbriel mean motion

    def _fn(moon_name: str, et_seconds: float) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        # et_seconds is treated as "seconds since the idealized t=0 epoch"
        # by the caller (see check_connection_survives_real_ephemeris's
        # reduction-test wiring) -- this substitute does not consult SPICE
        # at all, so its absolute epoch has no real-world meaning.
        t_s = et_seconds
        theta_frame = n1 * t_s  # inertial angle of the rotating frame's x-axis
        # Barycentre-to-Uranus offset (+mu, 0, 0) in rotating-frame nondim
        # units, expressed in inertial km/km-s at this instant -- the SAME
        # shift nondim_state_to_inertial applies to the spacecraft.
        mu_shift_r = l_km * mu * np.array([math.cos(theta_frame), math.sin(theta_frame), 0.0])
        mu_shift_v = l_km * mu * n1 * np.array([-math.sin(theta_frame), math.cos(theta_frame), 0.0])
        if moon_name == "Umbriel":
            # Uranus-Umbriel separation IS l_km by definition (the length
            # unit) -- Umbriel's Uranus-centred radius is exactly l_km, no
            # (1-mu) factor (that factor is barycentre-relative).
            r_umbriel_km = l_km
            theta = theta_frame
            r = np.array([r_umbriel_km * math.cos(theta), r_umbriel_km * math.sin(theta), 0.0])
            v = np.array(
                [-r_umbriel_km * n1 * math.sin(theta), r_umbriel_km * n1 * math.cos(theta), 0.0]
            )
            return r, v
        if moon_name == "Titania":
            r_titania_km = l_km * a_gan
            theta_rot = theta_gan0 + omega_gan * (n1 * t_s)  # rotating-frame angle
            theta = theta_frame + theta_rot  # inertial angle
            omega_titania_inertial = n1 * (1.0 + omega_gan)
            r = (
                np.array([r_titania_km * math.cos(theta), r_titania_km * math.sin(theta), 0.0])
                + mu_shift_r
            )
            v = (
                np.array(
                    [
                        -r_titania_km * omega_titania_inertial * math.sin(theta),
                        r_titania_km * omega_titania_inertial * math.cos(theta),
                        0.0,
                    ]
                )
                + mu_shift_v
            )
            return r, v
        raise ValueError(f"idealized_moon_state_fn: unsupported moon {moon_name!r}")

    return _fn


# --------------------------------------------------------------------------- #
# Top-level driver
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class ConsistencyCheckResult:
    """Headline result of one real-ephemeris consistency check at one epoch."""

    epoch0_utc: str
    et0: float
    l_km: float
    v_unit_km_s: float
    t_u_tu: float
    t_u_seconds: float
    r0_km: NDArray[np.float64]
    v0_km_s: NDArray[np.float64]
    r_f_km: NDArray[np.float64]
    v_f_km_s: NDArray[np.float64]
    r_target_km: NDArray[np.float64]
    v_target_km_s: NDArray[np.float64]
    pos_gap_km: float
    vel_gap_km_s: float
    propagation_success: bool
    closest_approach_km: dict[str, float]
    idealized_off_torus_km_for_scale: float
    """`#701`'s OWN ``corrected_off_torus_km`` for this candidate (the
    idealized model's own "genuinely off the torus" gate, 1000 km) --
    reported alongside for direct scale comparison, not used in this
    module's own computation."""
    notes: str = field(default="")


def check_connection_survives_real_ephemeris(
    epoch0_utc: str,
    departure_state6_nondim: NDArray[np.float64],
    target_state6_nondim: NDArray[np.float64],
    t_u_tu: float,
    l_km: float,
    v_unit_km_s: float,
    mu: float,
    *,
    perturber_moons: tuple[str, ...] = ("Umbriel", "Titania"),
    idealized_off_torus_km_for_scale: float = float("nan"),
    kernel_paths: tuple[str, str, str] | None = None,
    rtol: float = 1e-12,
    atol: float = 1e-6,
) -> ConsistencyCheckResult:
    """Run one real-ephemeris consistency check at a chosen real epoch.

    ``departure_state6_nondim`` is `#701`'s own unstable-manifold departure
    state (nondim, idealized rotating frame, BEFORE flowing --
    ``t_flow=0`` point). ``target_state6_nondim`` is `#701`'s own converged
    stable-manifold state (``RefinedConnection.state_s``) -- the point the
    idealized connection's own residual drives to near-coincide with the
    unstable arrival. ``t_u_tu`` is the SAME elapsed unstable-branch flow
    time (nondimensional TU) `#701`'s own connection used. ``mu`` is the
    system's own Umbriel mass ratio (``system.mu``, ``1.469e-5`` for
    Umbriel-Titania) -- the Uranus-vs-barycentre re-origining
    :func:`nondim_state_to_inertial` needs (see its own docstring on why
    this is NOT negligible for this system).

    FURNSHes the SPICE kernels for the duration of this call only (kclear'd
    in a finally, matching `#312`'s own ``run_v4_uranus_strict`` discipline
    -- does not pollute the SPICE kernel pool across calls).
    """
    if kernel_paths is None:
        kernel_paths = (str(DEFAULT_LSK_PATH), str(DEFAULT_PCK_PATH), str(DEFAULT_URA_PATH))

    mu_uranus = PRIMARIES["Uranus"]
    t_u_seconds = t_u_tu * tu_to_seconds(l_km, v_unit_km_s)

    import spiceypy as spice

    spice.kclear()
    try:
        _spice_furnsh_all(kernel_paths)
        et0 = _ephemeris_time_seconds(epoch0_utc)

        r_umbriel0, v_umbriel0 = spice_moon_state_fn("Umbriel", et0)
        r0_km, v0_km_s = nondim_state_to_inertial(
            departure_state6_nondim, r_umbriel0, v_umbriel0, l_km, v_unit_km_s, mu=mu
        )

        prop = propagate_real(
            r0_km,
            v0_km_s,
            t_u_seconds,
            et0=et0,
            mu_uranus=mu_uranus,
            perturber_moons=perturber_moons,
            moon_state_fn=spice_moon_state_fn,
            rtol=rtol,
            atol=atol,
        )

        et_target = et0 + t_u_seconds
        r_umbriel_target, v_umbriel_target = spice_moon_state_fn("Umbriel", et_target)
        r_target_km, v_target_km_s = nondim_state_to_inertial(
            target_state6_nondim, r_umbriel_target, v_umbriel_target, l_km, v_unit_km_s, mu=mu
        )
    finally:
        spice.kclear()

    if prop.success:
        pos_gap_km = float(np.linalg.norm(prop.r_f_km - r_target_km))
        vel_gap_km_s = float(np.linalg.norm(prop.v_f_km_s - v_target_km_s))
        notes = "converged"
    else:
        pos_gap_km = float("inf")
        vel_gap_km_s = float("inf")
        notes = "real-ephemeris propagation failed"

    return ConsistencyCheckResult(
        epoch0_utc=epoch0_utc,
        et0=et0,
        l_km=l_km,
        v_unit_km_s=v_unit_km_s,
        t_u_tu=t_u_tu,
        t_u_seconds=t_u_seconds,
        r0_km=r0_km,
        v0_km_s=v0_km_s,
        r_f_km=prop.r_f_km,
        v_f_km_s=prop.v_f_km_s,
        r_target_km=r_target_km,
        v_target_km_s=v_target_km_s,
        pos_gap_km=pos_gap_km,
        vel_gap_km_s=vel_gap_km_s,
        propagation_success=prop.success,
        closest_approach_km=prop.closest_approach_km,
        idealized_off_torus_km_for_scale=idealized_off_torus_km_for_scale,
        notes=notes,
    )
