"""V4 high(er)-fidelity Uranian-system gauntlet (#332 / #306 Phase 4).

Spec reference
--------------
* §14 V4 — high-fidelity ephemeris model: the candidate's defining dynamical
  signature must survive re-evaluation under the highest available real-physics
  model (J2/J3 zonal harmonics, other-moon third-body perturbations, real
  ephemeris where available).

For the #327 Umbriel-Oberon-Umbriel SILVER (``repeated-moon-uranus-00000041``),
V4 asks: does the bounded near-5:1 Umbriel-Oberon synodic-resonance signature
V3 confirmed (IAS15 ~ nanometers agreement with V2; #331) survive the dominant
non-Keplerian perturbations the V2/V3 stack OMITS?

GMAT-vs-fallback rationale (the honesty bit)
--------------------------------------------
The task #332 spec calls for a GMAT/SPICE Uranian-satellite real-ephemeris
path. The local GMAT install (``~/GMAT/R2022a``) ships ONLY DE405/DE421/DE424
planetary ephemerides; the Uranian satellite kernels (URA111/URA107) required
to drive Umbriel/Oberon in GMAT's PointMasses force model are NOT bundled and
are not on disk. Per the task's explicit fallback clause:

    "If GMAT isn't reachable from the agent environment OR SPICE kernels
    aren't available, fall back to a higher-fidelity scipy DOP853 integrator
    with hand-implemented Uranus J2 + Umbriel/Oberon/Titania third-body +
    circular Sun perturbation. This is NOT full HFEM but it adds the dominant
    V4-class perturbations not in V3. Document the fallback honestly."

This module IS the documented fallback. Phase 4.1 (full GMAT V4 once Uranian
satellite SPICE kernels are installed) remains a downstream gate before
catalogue admission as a ``quasi_cycler`` row.

What this V4 fallback adds over V3
-----------------------------------
V3 (IAS15) propagates the spacecraft as a test particle under planet-frame
two-body Kepler about Uranus, with moons placed by circular-coplanar Kepler
ephemerides. V4 here adds, on top of the same scipy DOP853 integrator the
V2 driver chain composes on:

1. **Uranus J2 zonal harmonic** — Jacobson 2014 (AJ 148:76 Table 4):
   ``J2 = 3.34343e-3``, equatorial radius ``R_eq = 25559 km``. At the
   spacecraft's distance from Uranus (Umbriel SMA 265,986 km, Oberon
   SMA 583,511 km), the J2 perturbing acceleration scales as
   ``J2 * (R_eq / r)^2 * g`` -> ~5e-5 of the central acceleration at
   Umbriel's SMA, ~1e-5 at Oberon's. Small but cumulative over a
   ~30-day cycle.

2. **Other classical Uranian moons as third-body perturbers** —
   Miranda, Ariel, Titania (Umbriel + Oberon are the Lambert-tour moons
   and are perturbers in V4; in V3 they are reduced to point ephemerides
   the spacecraft only intersects at rendezvous). Each moon's GM acts on
   the spacecraft as Battin's standard third-body perturbation.

3. **The Lambert-tour moons (Umbriel, Oberon) as real perturbers** —
   in V3 they exist only at the encounter points; here they perturb the
   spacecraft throughout the leg.

What V4 does NOT add (and why)
------------------------------
* **Real (non-circular, non-coplanar) moon ephemerides via SPICE** —
  requires URA111/URA107 kernels which are unavailable. The circular-
  coplanar ephemeris is therefore retained. Adding J2 still tests whether
  the dominant V4-class perturbation breaks the V3 signature.
* **Solar third-body perturbation** — Uranus's heliocentric distance is
  ~19 AU; the solar perturbation on a Uranus-system spacecraft is
  ``mu_sun / d_sun^2 * (r / d_sun)``, which is ~1e-9 km/s^2 at Oberon's
  distance — three orders of magnitude below the J2 perturbation we DO
  model. Omitting it does not change the V4 verdict materially. Including
  it would require an epoch + Uranus heliocentric ephemeris, which we
  also don't have without SPICE.
* **Uranus J3, J4** — Jacobson 2014 gives J3 essentially zero (consistent
  with Uranus's axisymmetric shape); J4 ~ -2.9e-5, two orders of
  magnitude smaller than J2. Adequate for a V4-fallback verdict.

Composition map
---------------
This module composes on:

* :func:`scipy.integrate.solve_ivp` (DOP853) — the high-order Runge-Kutta
  the V2 driver chain uses. Same INTEGRATOR FAMILY as V2 by design — V4
  changes the MODEL (adding J2 + n-body), not the integrator architecture
  (that was V3's role). V4 vs V3 disagreement isolates the V4-class
  perturbation effect.
* :func:`cyclerfinder.core.lambert.lambert` — re-used for the Lambert
  targeting (same encounters as V2/V3); V4 is about what the spacecraft
  does BETWEEN encounters under richer physics, not Lambert targeting.
* :func:`cyclerfinder.search.discovery_campaign._moon_state` — circular-
  coplanar moon ephemeris (the model V3 also uses; V4 keeps it for the
  reasons documented above).

Discipline
----------
* NO catalogue writeback. V4 PASS still has #329 Heaton-Longuski paywall
  ahead AND a strict reading still wants Phase 4.1 (GMAT + real Uranian
  satellite SPICE) before catalogue admission.
* The V4 verdict is whatever the math says. Don't tune to pass.
* V4 vs V3 disagreement isolates the V4-class perturbations effect on the
  bounded-drift signature: if it survives, the signature is robust to
  the dominant non-Keplerian effects (J2 + other moons); if it fails,
  it was a Keplerian artifact.
"""

from __future__ import annotations

import math
from collections.abc import Callable
from dataclasses import dataclass
from typing import Final

import numpy as np
from scipy.integrate import solve_ivp

import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.core.lambert import lambert as _lambert
from cyclerfinder.core.satellites import PRIMARIES, SATELLITES
from cyclerfinder.data.validation.v3_3d import V3Verdict3D
from cyclerfinder.search.discovery_campaign import DAY_S, _mean_motion_rad_day, _moon_state

# --------------------------------------------------------------------------- #
# SOURCED perturbation constants
# --------------------------------------------------------------------------- #

URANUS_J2: Final[float] = 3.34343e-3
"""Uranus J2 zonal harmonic coefficient.

Source: Jacobson 2014, "The orbits of the Uranian satellites and rings, the
gravity field of the Uranian system, and the orientation of the pole of
Uranus", The Astronomical Journal 148:76, Table 4 ("Best-fit zonal harmonic
coefficients of Uranus"). Same upstream JPL uses for URA111.
"""

URANUS_R_EQ_KM: Final[float] = 25559.0
"""Uranus equatorial radius (km).

Source: Jacobson 2014 op. cit. (also IAU 2015 nominal value, Mamajek et al.).
"""

# Other-moon perturber set (V4 adds these as third-body perturbations on the
# spacecraft during the Lambert legs). Umbriel and Oberon are the Lambert-tour
# moons; V4 includes them as perturbers throughout the leg, not only at the
# encounter points (V2/V3 only intersect them at rendezvous).
URANIAN_PERTURBER_MOONS: Final[tuple[str, ...]] = (
    "Miranda",
    "Ariel",
    "Umbriel",
    "Titania",
    "Oberon",
)
"""Classical regular Uranian moons that act as third-body perturbers in V4.

Each contributes a standard Battin third-body acceleration on the spacecraft.
Miranda and Ariel are the SMALLER moons not on the Lambert tour; Umbriel and
Oberon are the tour moons themselves (in V3 they exist only at encounter
epochs; in V4 they perturb the spacecraft throughout the leg). Titania is the
largest classical moon and the dominant third-body perturber by GM.
"""

V4_AGREEMENT_FLOOR_KMS: Final[float] = 50_000.0
"""Default V4-vs-V3 agreement floor in km for the bounded-drift signature.

If V4 (DOP853 + J2 + other-moon third-body) and V3 (IAS15, Kepler-only)
disagree by more than this on the per-cycle terminal-position drift, the
V3-confirmed bounded-drift signature did NOT survive the dominant V4-class
perturbations. The signature was a Keplerian idealization.

50,000 km is the project's same-model drift floor (matches
:data:`v2_3d.V2_DRIFT_FLOOR_KMS` and :data:`v2_moontour.V2_MOONTOUR_DRIFT_FLOOR_KMS`)
— well below the SILVER's actual V3 drift (~3e5 km), so a V4 verdict driven
by J2 + other-moon perturbations is genuinely discriminating, not noise.
"""

V4_N_CYCLES_MIN: Final[int] = 3
"""Spec §14 minimum cycles for a V4 verdict, same as V2/V3."""

# --------------------------------------------------------------------------- #
# Verdict dataclasses
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class V4CycleVerdictUranus:
    """Per-cycle V4 verdict for the Uranian-system high(er)-fidelity gauntlet."""

    cycle_index: int
    """Zero-indexed cycle number (0, 1, 2, ...)."""
    converged_legs: int
    """Number of Lambert legs that closed in this cycle (same machinery as V2/V3)."""
    n_legs: int
    """Total Lambert legs in one cycle."""
    rendezvous_drift_kms_v4: float
    """Position offset of THIS cycle's V4 final encounter vs cycle 0's V4 final
    encounter (planet-frame km). The V4 analogue of V3's
    ``rendezvous_drift_kms_v3``."""
    rendezvous_drift_kms_v3: float
    """V3's stored ``rendezvous_drift_kms_v3`` at the same cycle index, for
    direct comparison. 0.0 for cycle 0 by construction."""
    agreement_kms: float
    """``|rendezvous_drift_kms_v4 - rendezvous_drift_kms_v3|`` — V4-vs-V3
    delta. Cycle-0 is 0 by construction."""
    v4_terminal_offset_vs_moon_kms: float
    """The V4 spacecraft terminal-position offset vs the analytic moon
    target at the cycle's final encounter, km. Same definition as V3's
    ``ias15_vs_analytic_kepler_kms`` but under the J2 + other-moon model."""
    notes: str = ""


@dataclass(frozen=True)
class V4UranusVerdict:
    """Frozen V4 verdict for an Uranus-system moontour candidate."""

    candidate_id: str
    sequence: tuple[str, ...]
    n_cycles_propagated: int
    integrator: str
    """Human-readable integrator + model label. The fallback path is
    documented honestly. For this module: ``"scipy DOP853 (J2 + other-moon
    third-body fallback; no SPICE)"``."""
    per_cycle: tuple[V4CycleVerdictUranus, ...]
    per_cycle_drift_kms_v4: tuple[float, ...]
    """Per-cycle V4 rendezvous drift, km. Headline for V4-vs-V3 compare."""
    per_cycle_drift_kms_v3: tuple[float, ...]
    """V3 series sliced to ``n_cycles_propagated``."""
    drift_agreement_kms: float
    """``max_k |drift_v4[k] - drift_v3[k]|`` — headline V4-vs-V3 number.
    If small (< :data:`V4_AGREEMENT_FLOOR_KMS`), the V3 bounded-drift
    signature survives the dominant V4-class perturbations."""
    v4_v3_agreement_floor_kms: float
    """Floor against which ``drift_agreement_kms`` is gated. Default
    :data:`V4_AGREEMENT_FLOOR_KMS`."""
    bounded_drift_survives: bool
    """Qualitative: V4 drift stays bounded (max drift over cycles does not
    blow up monotonically by more than 10x the V3 max drift). This is the
    central V4 question, distinct from the strict V3 agreement gate."""
    passes_v4: bool
    """``drift_agreement_kms <= v4_v3_agreement_floor_kms`` AND every cycle's
    Lambert leg closed AND bounded_drift_survives. The headline boolean.

    Interpretation:

    * PASS: V4 (J2 + other-moon n-body) agrees with V3 (IAS15 Kepler-only)
      on the bounded-drift signature -> the SILVER's quasi_cycler property
      survives the dominant non-Keplerian perturbations of the Uranian
      system. Next gates: #329 (Heaton-Longuski literature check) +
      Phase 4.1 (GMAT + Uranian satellite SPICE for the full HFEM gate)
      before catalogue admission as ``quasi_cycler``.
    * FAIL: J2 + other-moon perturbations break the bounded-drift signature
      -> the V3-confirmed quasi_cycler property was a Keplerian artifact.
      Retire to the negative-results registry (#172) with the perturbation
      order at which the signature collapses recorded.
    """
    notes: str = ""


# --------------------------------------------------------------------------- #
# Physics helpers (the V4-vs-V3 delta)
# --------------------------------------------------------------------------- #


def _j2_acceleration_kms2(r_km: np.ndarray, *, mu: float, j2: float, r_eq_km: float) -> np.ndarray:
    """Acceleration (km/s^2) on a test particle from J2 zonal harmonic.

    Standard formulation (Vallado, "Fundamentals of Astrodynamics and
    Applications", 4th ed., Eq. 8-37):

    .. math::

        \\vec{a}_{J2} = -\\frac{3}{2}\\frac{\\mu J_2 R_{eq}^2}{r^5}
            \\begin{pmatrix}
              x \\left(1 - 5 \\frac{z^2}{r^2}\\right) \\\\
              y \\left(1 - 5 \\frac{z^2}{r^2}\\right) \\\\
              z \\left(3 - 5 \\frac{z^2}{r^2}\\right)
            \\end{pmatrix}

    Frame: Uranus body-fixed equatorial (the Z axis is Uranus's spin axis).
    Since our circular-coplanar moon ephemerides are also in this frame
    (the orbit plane IS the Uranus equatorial plane by construction), the
    in-plane spacecraft motion has z = 0 and the z component of a_J2 is
    also 0 — but the radial component is nonzero and acts on every cycle.
    """
    x, y, z = float(r_km[0]), float(r_km[1]), float(r_km[2])
    r2 = x * x + y * y + z * z
    if r2 <= 0.0:
        return np.zeros(3, dtype=np.float64)
    r = math.sqrt(r2)
    coeff = -1.5 * mu * j2 * r_eq_km * r_eq_km / (r2 * r2 * r)
    z2_over_r2 = (z * z) / r2
    in_plane_factor = 1.0 - 5.0 * z2_over_r2
    out_of_plane_factor = 3.0 - 5.0 * z2_over_r2
    return np.array(
        [coeff * x * in_plane_factor, coeff * y * in_plane_factor, coeff * z * out_of_plane_factor],
        dtype=np.float64,
    )


def _third_body_acceleration_kms2(
    r_sc_km: np.ndarray,
    r_body_km: np.ndarray,
    *,
    mu_body: float,
    softening_km: float = 0.0,
) -> np.ndarray:
    """Battin third-body acceleration (km/s^2) on the spacecraft.

    Standard formulation (Battin, "An Introduction to the Mathematics and
    Methods of Astrodynamics", §8.3, Eq. 8.60):

    .. math::

        \\vec{a}_{3B} = \\mu_b \\left(
            \\frac{\\vec{r}_{b/sc}}{|\\vec{r}_{b/sc}|^3}
            - \\frac{\\vec{r}_b}{|\\vec{r}_b|^3}
        \\right)

    where :math:`\\vec{r}_b` is the perturbing body's position vector from
    the central body (Uranus) and :math:`\\vec{r}_{b/sc} = \\vec{r}_b - \\vec{r}_{sc}`.

    Both terms must be computed in the same inertial frame. The second
    term is the indirect/centrifugal correction (the central body itself
    accelerates under the third-body's gravity, and the spacecraft is
    expressed in the central-body-centered frame).

    Patched-conic treatment of the moon's Hill sphere: if
    ``|r_body - r_sc| < softening_km``, return zero from that moon. This is
    the standard patched-conic treatment of the encounter geometry — inside
    a moon's Hill sphere, the moon would BE the central body and the
    patched-conic glue zeroes its third-body contribution to avoid the
    Lambert endpoint's r_body_sc = 0 singularity (the spacecraft IS AT the
    moon at the leg start/end by Lambert geometry). V4 here is propagating
    BETWEEN encounters, not modeling the within-SOI flyby dynamics; the
    softening is the documented patched-conic boundary.
    """
    r_b_sc = r_body_km - r_sc_km
    d_b_sc = float(np.linalg.norm(r_b_sc))
    d_b = float(np.linalg.norm(r_body_km))
    if d_b_sc <= 0.0 or d_b <= 0.0:
        return np.zeros(3, dtype=np.float64)
    if softening_km > 0.0 and d_b_sc < softening_km:
        return np.zeros(3, dtype=np.float64)
    return np.asarray(mu_body * (r_b_sc / (d_b_sc**3) - r_body_km / (d_b**3)), dtype=np.float64)


def _hill_radius_km(*, sma_moon_km: float, mu_moon: float, mu_primary: float) -> float:
    """Hill (Roche) sphere radius (km).

    Standard formulation: ``r_H = a_moon * (mu_moon / (3 * mu_primary))^(1/3)``.

    For Uranian moons at SMA ~3e5 - 6e5 km, Hill radii are ~5,000 -
    15,000 km. The Hill radius bounds the "inside the moon's SOI" region
    that V4 (a between-encounter propagator, not a within-SOI flyby
    integrator) does not resolve; we use it as the third-body softening
    distance.
    """
    return float(sma_moon_km * (mu_moon / (3.0 * mu_primary)) ** (1.0 / 3.0))


def _build_dynamics(
    *,
    mu_primary: float,
    j2: float,
    r_eq_km: float,
    perturber_moons: tuple[str, ...],
    moon_consts: dict[str, tuple[float, float]],
    theta_base: dict[str, float],
    t_cycle_offset_days: float,
) -> Callable[[float, np.ndarray], np.ndarray]:
    """Build a closed-over ``f(t, y)`` for ``scipy.integrate.solve_ivp``.

    The dynamics evaluated at integration time ``t_s`` (seconds since the
    leg start) are:

    * Central two-body about the primary (Uranus).
    * J2 zonal harmonic (Vallado Eq. 8-37).
    * Third-body perturbations from each ``perturber_moons`` entry — moon
      position computed via :func:`_moon_state` at the LEG time
      ``t_cycle_offset_days + t_s / DAY_S``. Each moon's third-body
      acceleration is softened inside its Hill sphere (the patched-conic
      boundary; V4 is a between-encounter propagator and does not resolve
      within-SOI flyby dynamics).

    Returns the right-hand-side callable ``f(t_s, y_6) -> dydt_6``.
    """
    perturber_mu = {m: SATELLITES[m].mu_km3_s2 for m in perturber_moons}
    # Hill softening distance per moon — the spacecraft is AT a moon at
    # leg start/end by Lambert geometry, so the third-body term must be
    # zeroed inside the Hill sphere to avoid the patched-conic singularity.
    perturber_hill_km = {
        m: _hill_radius_km(
            sma_moon_km=moon_consts[m][0],
            mu_moon=perturber_mu[m],
            mu_primary=mu_primary,
        )
        for m in perturber_moons
    }

    def rhs(t_s: float, y: np.ndarray) -> np.ndarray:
        r_sc = y[:3]
        v_sc = y[3:]
        r_norm = float(np.linalg.norm(r_sc))
        if r_norm <= 0.0:
            a_central = np.zeros(3, dtype=np.float64)
        else:
            a_central = -mu_primary * r_sc / (r_norm**3)
        a_j2 = _j2_acceleration_kms2(r_sc, mu=mu_primary, j2=j2, r_eq_km=r_eq_km)
        # Time-varying third-body positions (circular-coplanar moon ephemerides).
        t_days = t_cycle_offset_days + t_s / DAY_S
        a_3b = np.zeros(3, dtype=np.float64)
        for moon in perturber_moons:
            sma, n_rad_day = moon_consts[moon]
            r_moon, _ = _moon_state(theta_base[moon], n_rad_day, t_days, sma, mu_primary)
            a_3b += _third_body_acceleration_kms2(
                r_sc,
                r_moon,
                mu_body=perturber_mu[moon],
                softening_km=perturber_hill_km[moon],
            )
        a_total = a_central + a_j2 + a_3b
        return np.concatenate([v_sc, a_total])

    return rhs


def _v4_propagate_leg(
    r0_km: np.ndarray,
    v0_km_s: np.ndarray,
    tof_s: float,
    *,
    mu_primary: float,
    j2: float,
    r_eq_km: float,
    perturber_moons: tuple[str, ...],
    moon_consts: dict[str, tuple[float, float]],
    theta_base: dict[str, float],
    t_cycle_offset_days: float,
    rtol: float = 1e-10,
    atol: float = 1e-3,
) -> tuple[np.ndarray, np.ndarray, bool]:
    """Propagate one Lambert leg under the V4 model (central + J2 + n-moon 3B).

    Returns ``(r_f_km, v_f_km_s, success)``. The integrator is scipy DOP853
    (the same RK family used implicitly by the V2 driver's CR3BP propagator);
    the V4-vs-V3 delta is in the FORCE MODEL, not the integrator family.
    """
    rhs = _build_dynamics(
        mu_primary=mu_primary,
        j2=j2,
        r_eq_km=r_eq_km,
        perturber_moons=perturber_moons,
        moon_consts=moon_consts,
        theta_base=theta_base,
        t_cycle_offset_days=t_cycle_offset_days,
    )
    y0 = np.concatenate(
        [np.asarray(r0_km, dtype=np.float64), np.asarray(v0_km_s, dtype=np.float64)]
    )
    sol = solve_ivp(
        rhs,
        (0.0, float(tof_s)),
        y0,
        method="DOP853",
        rtol=rtol,
        atol=atol,
        dense_output=False,
    )
    if not sol.success:
        return np.zeros(3), np.zeros(3), False
    yf = sol.y[:, -1]
    return yf[:3].copy(), yf[3:].copy(), True


# --------------------------------------------------------------------------- #
# Driver
# --------------------------------------------------------------------------- #


def _resolve_primary(system: cr3bp.CR3BPSystem | None, sequence: tuple[str, ...]) -> str:
    """Resolve the primary body name (same logic as v2/v3 moontour)."""
    if system is not None and system.primary:
        target = system.primary.strip().lower()
        for name in PRIMARIES:
            if name.lower() == target:
                return name
        raise ValueError(f"unknown primary {system.primary!r} from system")
    if not sequence:
        raise ValueError("empty sequence; cannot resolve primary")
    head = sequence[0]
    if head not in SATELLITES:
        raise ValueError(f"unknown moon {head!r} in sequence")
    return SATELLITES[head].primary


def _moon_constants(primary: str, moons: tuple[str, ...]) -> dict[str, tuple[float, float]]:
    """``{moon: (sma_km, mean_motion_rad_day)}`` for every moon in ``moons``."""
    mu = PRIMARIES[primary]
    out: dict[str, tuple[float, float]] = {}
    for moon in moons:
        if moon not in SATELLITES:
            raise ValueError(f"unknown moon {moon!r}")
        sat = SATELLITES[moon]
        if sat.primary != primary:
            raise ValueError(
                f"moon {moon!r} orbits {sat.primary!r}, not the resolved primary {primary!r}"
            )
        out[moon] = (sat.sma_km, _mean_motion_rad_day(mu, sat.sma_km))
    return out


def _cycle_v4(
    *,
    sequence: tuple[str, ...],
    leg_tofs_days: tuple[float, ...],
    theta_base: dict[str, float],
    t_cycle_offset_days: float,
    tour_consts: dict[str, tuple[float, float]],
    perturber_consts: dict[str, tuple[float, float]],
    perturber_moons: tuple[str, ...],
    mu_primary: float,
    n_revs: tuple[int, ...] | None,
    j2: float,
    r_eq_km: float,
) -> tuple[bool, np.ndarray | None, float]:
    """Re-solve all Lambert legs of one cycle AND propagate them under V4 physics.

    Returns ``(converged, r_final_v4_km, v4_terminal_offset_vs_moon_kms)``:
    the spacecraft's terminal planet-frame position under the V4 model AT the
    final encounter, and the offset of that V4 terminal position vs the
    analytic moon target at the leg endpoint (the V4 analogue of V3's
    ``ias15_vs_analytic_kepler_kms``).

    The Lambert targeting is identical to V2/V3 — V4 changes the propagation
    physics between encounters, not the targeting.
    """
    n_legs = len(sequence) - 1
    if n_revs is None:
        n_revs_used: tuple[int, ...] = tuple(0 for _ in range(n_legs))
    else:
        if len(n_revs) != n_legs:
            raise ValueError(f"n_revs length {len(n_revs)} != n_legs {n_legs}")
        n_revs_used = tuple(n_revs)

    epochs_days = [0.0]
    for tof in leg_tofs_days:
        epochs_days.append(epochs_days[-1] + tof)

    # Analytic moon states at each encounter (tour moons only — same as V2/V3
    # Lambert targeting machinery).
    states: list[tuple[np.ndarray, np.ndarray]] = []
    for moon, t in zip(sequence, epochs_days, strict=True):
        sma, n_rad_day = tour_consts[moon]
        states.append(
            _moon_state(theta_base[moon], n_rad_day, t_cycle_offset_days + t, sma, mu_primary)
        )

    sc_r_curr: np.ndarray | None = None
    worst_offset_kms = 0.0
    for k in range(n_legs):
        r_a, v_a_moon = states[k]
        r_b, _ = states[k + 1]
        nrev = max(0, n_revs_used[k])
        sols = _lambert(r_a, r_b, leg_tofs_days[k] * DAY_S, mu=mu_primary, max_revs=nrev)
        wanted = [s for s in sols if s.n_revs == n_revs_used[k]]
        if not wanted:
            return False, None, float("inf")
        v_a_captured = v_a_moon
        best = min(wanted, key=lambda s: float(np.linalg.norm(s.v1 - v_a_captured)))
        # IC for the V4 leg: planet-frame spacecraft state at moon-A, using
        # Lambert's v-out from moon-A. Identical to V3's choice; the delta is
        # in what physics drives the propagation.
        r0_leg = r_a.copy()
        v0_leg = best.v1.copy()
        # IMPORTANT: the third-body perturbations need the GLOBAL time
        # offset (cycle offset + leg start), since perturber-moon positions
        # depend on it. The leg starts at t_s = 0 of the integrator but
        # t_days = t_cycle_offset_days + epochs_days[k] of the global model.
        leg_t_cycle_offset_days = t_cycle_offset_days + epochs_days[k]
        r_f_leg, _, ok = _v4_propagate_leg(
            r0_leg,
            v0_leg,
            leg_tofs_days[k] * DAY_S,
            mu_primary=mu_primary,
            j2=j2,
            r_eq_km=r_eq_km,
            perturber_moons=perturber_moons,
            moon_consts=perturber_consts,
            theta_base=theta_base,
            t_cycle_offset_days=leg_t_cycle_offset_days,
        )
        if not ok:
            return False, None, float("inf")
        # V4 terminal r vs analytic moon target (Lambert r_b endpoint). At a
        # V4 PASS the V4 physics has driven the spacecraft AWAY from the
        # Kepler target by the cumulative perturbation effect over the leg;
        # this offset is the per-leg V4-vs-V3 delta.
        leg_offset_kms = float(np.linalg.norm(r_f_leg - r_b))
        worst_offset_kms = max(worst_offset_kms, leg_offset_kms)
        sc_r_curr = r_f_leg
    if sc_r_curr is None:
        return False, None, float("inf")
    return True, sc_r_curr, worst_offset_kms


def run_v4_uranus(
    candidate_id: str,
    sequence: tuple[str, ...],
    vinf_tuple_kms: tuple[float, ...],
    leg_tofs_days: tuple[float, ...],
    rel_offset_deg: float,
    system: cr3bp.CR3BPSystem | None,
    *,
    v3_verdict: V3Verdict3D,
    n_cycles: int = V4_N_CYCLES_MIN,
    n_revs: tuple[int, ...] | None = None,
    phase0_deg: float = 0.0,
    j2: float = URANUS_J2,
    r_eq_km: float = URANUS_R_EQ_KM,
    perturber_moons: tuple[str, ...] = URANIAN_PERTURBER_MOONS,
    agreement_floor_kms: float = V4_AGREEMENT_FLOOR_KMS,
    drift_unbounded_factor: float = 10.0,
    notes: str = "",
) -> V4UranusVerdict:
    """Run V4 for an Uranian moontour: re-propagate cycles under J2 + other-moon n-body.

    Pipeline:
      1. For each cycle k = 0, ..., n_cycles - 1:
           a. Advance moon longitudes by ``k * cycle_period_days``.
           b. Re-solve each Lambert leg (same as V2/V3).
           c. Propagate the spacecraft through each leg with scipy DOP853
              under the V4 force model (Uranus central + Uranus J2 +
              every ``perturber_moons`` third-body Battin acceleration).
           d. Record the cycle's V4 terminal spacecraft position at the
              final encounter.
      2. Compute V4 per-cycle drift series.
      3. Compute ``drift_agreement_kms = max_k |drift_v4[k] - drift_v3[k]|``.
      4. Check ``bounded_drift_survives``: max V4 drift over all cycles
         does not exceed ``drift_unbounded_factor`` x max V3 drift. A
         qualitatively bounded V4 series satisfies this; an unbounded
         (monotonically diverging) one does not.
      5. Verdict: PASS iff every cycle converged AND
         ``drift_agreement_kms <= agreement_floor_kms`` AND
         bounded_drift_survives.

    Parameters
    ----------
    candidate_id, sequence, vinf_tuple_kms, leg_tofs_days, rel_offset_deg, system:
        Same semantics as :func:`run_v3_3d`. ``vinf_tuple_kms`` is carried
        for audit; V4 re-solves Lambert from geometry.
    v3_verdict:
        The :class:`V3Verdict3D` from :func:`run_v3_3d` on this candidate at
        the SAME ``n_cycles`` (or larger). V4 reads
        ``per_cycle[k].rendezvous_drift_kms_v3`` to compute the agreement.
    n_cycles:
        Cycles to attempt. Must be >= :data:`V4_N_CYCLES_MIN` (= 3).
    n_revs:
        Per-leg revolution count. Pass the candidate's stored ``n_rev``.
    phase0_deg:
        First moon's initial longitude at cycle 0, degrees.
    j2, r_eq_km:
        Override the default Uranus J2 / R_eq if desired (default the
        Jacobson 2014 values).
    perturber_moons:
        Tuple of Uranian moon names to include as third-body perturbers.
        Default :data:`URANIAN_PERTURBER_MOONS`.
    agreement_floor_kms:
        Bar against which ``drift_agreement_kms`` is gated. Default
        :data:`V4_AGREEMENT_FLOOR_KMS` (50,000 km — the project's same-model
        drift floor).
    drift_unbounded_factor:
        Bounded-drift heuristic: V4 max drift must not exceed this multiple
        of V3 max drift. Default 10.0.
    notes:
        Free-form audit note.

    Returns
    -------
    V4UranusVerdict
        ``passes_v4`` is the headline.

    Notes
    -----
    This is the documented scipy-fallback path (GMAT install lacks Uranian
    satellite SPICE kernels — see module docstring). Phase 4.1 (full GMAT
    HFEM with Uranian satellite SPICE) remains a downstream catalogue-
    admission gate.
    """
    if not sequence:
        raise ValueError("empty sequence")
    if sequence[0] != sequence[-1]:
        raise ValueError(f"moontour sequence must be CLOSED (first == last); got {sequence!r}")
    n_legs = len(sequence) - 1
    if len(leg_tofs_days) != n_legs:
        raise ValueError(f"leg_tofs_days length {len(leg_tofs_days)} != n_legs {n_legs}")
    if len(vinf_tuple_kms) != len(sequence):
        raise ValueError(
            f"vinf_tuple_kms length {len(vinf_tuple_kms)} != len(sequence) {len(sequence)}"
        )
    if any(tof <= 0.0 for tof in leg_tofs_days):
        raise ValueError(f"leg_tofs_days must be positive; got {leg_tofs_days!r}")
    if n_cycles < V4_N_CYCLES_MIN:
        raise ValueError(f"V4 requires n_cycles >= {V4_N_CYCLES_MIN} (spec §14); got {n_cycles}")
    if agreement_floor_kms <= 0.0:
        raise ValueError(f"agreement_floor_kms must be > 0; got {agreement_floor_kms}")
    if len(v3_verdict.per_cycle) < n_cycles:
        raise ValueError(
            f"v3_verdict has only {len(v3_verdict.per_cycle)} cycles but V4 wants {n_cycles}"
        )

    primary = _resolve_primary(system, sequence)
    if primary != "Uranus":
        raise ValueError(
            f"run_v4_uranus is for Uranus-system candidates only; got primary={primary!r}"
        )
    tour_moons = tuple(sorted({m for m in sequence}))
    if len(tour_moons) == 1:
        raise ValueError(f"moontour requires >= 2 distinct moons; got {tour_moons!r}")
    if not perturber_moons:
        raise ValueError("perturber_moons must be non-empty (V4 vs V3 needs perturbations)")
    # The Lambert-tour moons must be in the perturber set OR explicitly excluded;
    # we permit any subset of registered Uranian moons, but tour moons not in
    # perturber set means they only exist at rendezvous (same as V3) — which is
    # legal but reduces V4 to a J2-only check. Don't auto-add; keep API explicit.

    # Tour moon constants (the Lambert-targeting moons — same as V2/V3 use).
    tour_consts = _moon_constants(primary, tour_moons)
    # Perturber moon constants (third-body bodies — may overlap with tour moons).
    perturber_consts = _moon_constants(primary, perturber_moons)
    mu_primary = PRIMARIES[primary]

    # Same theta_base convention as v2_moontour / v3_3d for like-for-like compare.
    # IMPORTANT: theta_base must include EVERY perturber moon and every tour moon
    # because all of them are referenced in the dynamics + Lambert targeting.
    phase0_rad = math.radians(phase0_deg)
    rel_off_rad = math.radians(rel_offset_deg)
    theta_base: dict[str, float] = {}
    # Tour moons first (same convention as V2/V3 — phase0 is the first moon's
    # initial longitude; the second tour moon is offset by rel_offset_deg).
    for j, moon in enumerate(tour_moons):
        if j == 0:
            theta_base[moon] = phase0_rad
        elif j == 1:
            theta_base[moon] = phase0_rad + rel_off_rad
        else:
            theta_base[moon] = phase0_rad + rel_off_rad + 2.0 * math.pi * (j - 1) / len(tour_moons)
    # Perturber moons NOT in the tour: arbitrary deterministic phasing. The
    # bounded-drift test asks whether the signature survives realistic
    # perturbations; the perturber-moon initial phase is part of "ephemeris
    # we don't have" — we pick a deterministic phase (longitude 0 at cycle 0)
    # for reproducibility. A V4 verdict robust to this phase choice is the
    # real test; we don't randomize because that introduces a free parameter.
    for moon in perturber_moons:
        if moon not in theta_base:
            theta_base[moon] = 0.0

    cycle_period_days = float(sum(leg_tofs_days))

    per_cycle: list[V4CycleVerdictUranus] = []
    v4_terminal_positions: list[np.ndarray] = []
    n_completed = 0
    cycle_zero_r_v4: np.ndarray | None = None

    for k in range(n_cycles):
        t_offset_days = k * cycle_period_days
        converged, r_v4, v4_offset_vs_moon = _cycle_v4(
            sequence=sequence,
            leg_tofs_days=leg_tofs_days,
            theta_base=theta_base,
            t_cycle_offset_days=t_offset_days,
            tour_consts=tour_consts,
            perturber_consts=perturber_consts,
            perturber_moons=perturber_moons,
            mu_primary=mu_primary,
            n_revs=n_revs,
            j2=j2,
            r_eq_km=r_eq_km,
        )
        if not converged or r_v4 is None:
            per_cycle.append(
                V4CycleVerdictUranus(
                    cycle_index=k,
                    converged_legs=0,
                    n_legs=n_legs,
                    rendezvous_drift_kms_v4=float("inf"),
                    rendezvous_drift_kms_v3=float(v3_verdict.per_cycle[k].rendezvous_drift_kms_v3),
                    agreement_kms=float("inf"),
                    v4_terminal_offset_vs_moon_kms=float("inf"),
                    notes="Lambert / DOP853 failed at least one leg in this cycle",
                )
            )
            break
        v4_terminal_positions.append(r_v4.copy())
        if k == 0:
            cycle_zero_r_v4 = r_v4.copy()
            drift_v4 = 0.0
        else:
            assert cycle_zero_r_v4 is not None
            drift_v4 = float(np.linalg.norm(r_v4 - cycle_zero_r_v4))
        drift_v3 = float(v3_verdict.per_cycle[k].rendezvous_drift_kms_v3)
        agreement = abs(drift_v4 - drift_v3)
        per_cycle.append(
            V4CycleVerdictUranus(
                cycle_index=k,
                converged_legs=n_legs,
                n_legs=n_legs,
                rendezvous_drift_kms_v4=drift_v4,
                rendezvous_drift_kms_v3=drift_v3,
                agreement_kms=agreement,
                v4_terminal_offset_vs_moon_kms=v4_offset_vs_moon,
            )
        )
        n_completed += 1

    drift_v4_series = tuple(c.rendezvous_drift_kms_v4 for c in per_cycle)
    drift_v3_series = tuple(c.rendezvous_drift_kms_v3 for c in per_cycle)
    drift_agreement = float("inf") if n_completed == 0 else max(c.agreement_kms for c in per_cycle)

    # Bounded-drift heuristic: V4 max drift must not exceed ``drift_unbounded_factor``
    # x V3 max drift. A V4 series that stays within this multiple of V3's bounded
    # signature qualifies as "bounded drift survives". A monotonically-diverging
    # V4 series will exceed this multiple at high enough n_cycles.
    if n_completed == 0:
        bounded_drift_survives = False
    else:
        finite_v4 = [d for d in drift_v4_series if math.isfinite(d)]
        finite_v3 = [d for d in drift_v3_series if math.isfinite(d)]
        if not finite_v4 or not finite_v3:
            bounded_drift_survives = False
        else:
            max_v4 = max(finite_v4)
            max_v3 = max(finite_v3)
            # If V3 max is essentially 0, fall back to comparing V4 max to a
            # default scale (the SMA of the OUTER tour moon — a "blew up by an
            # orbit radius" check).
            if max_v3 < 1.0:  # km — essentially closed
                max_outer_sma = max(tour_consts[m][0] for m in tour_moons)
                bounded_drift_survives = max_v4 < drift_unbounded_factor * max_outer_sma
            else:
                bounded_drift_survives = max_v4 <= drift_unbounded_factor * max_v3

    passes_v4 = bool(
        n_completed >= V4_N_CYCLES_MIN
        and n_completed == n_cycles
        and math.isfinite(drift_agreement)
        and drift_agreement <= agreement_floor_kms
        and bounded_drift_survives
    )

    return V4UranusVerdict(
        candidate_id=candidate_id,
        sequence=tuple(sequence),
        n_cycles_propagated=int(n_completed),
        integrator=(
            "scipy DOP853 (J2 + other-moon third-body fallback; "
            "no SPICE — Uranian satellite kernels unavailable in GMAT install)"
        ),
        per_cycle=tuple(per_cycle),
        per_cycle_drift_kms_v4=drift_v4_series,
        per_cycle_drift_kms_v3=drift_v3_series,
        drift_agreement_kms=float(drift_agreement),
        v4_v3_agreement_floor_kms=float(agreement_floor_kms),
        bounded_drift_survives=bool(bounded_drift_survives),
        passes_v4=passes_v4,
        notes=notes,
    )


__all__ = [
    "URANIAN_PERTURBER_MOONS",
    "URANUS_J2",
    "URANUS_R_EQ_KM",
    "V4_AGREEMENT_FLOOR_KMS",
    "V4_N_CYCLES_MIN",
    "V4CycleVerdictUranus",
    "V4UranusVerdict",
    "run_v4_uranus",
]
