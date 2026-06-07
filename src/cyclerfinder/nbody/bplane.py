"""B-plane targeting + powered-flyby kernel (Jones AAS 17-577 Eqs. 1-5, #142).

The targeting structure the Phase C shooter was missing. The method deep-dive
(``docs/notes/2026-06-07-jones-aas17-577-method-deepdive.md`` §2.2 / §5) flagged the
B-plane frame Ŝ,T̂,R̂ + θ_B (their Eqs. 4-5) as **DOES NOT MAP — Phase C shooter
needs this**, and named it the single most implementable finding. This module
codes those five equations verbatim, with the published tolerances as the gate
thresholds, so the shooter's flyby evaluation matches Jones's broad-search filter.

Equations (Jones, Hernandez & Jesick, AAS 17-577, pp. 5-6, transcribed in the
deep-dive note; corroborated by Russell 2004 Eq. 5.5 for the powered-SOI Δv):

  Eq. 1  turn angle           δ = ∠(v∞⁻, v∞⁺)
  Eq. 2  periapsis radius     asin(μ/(μ + r_p v∞⁻²)) + asin(μ/(μ + r_p v∞⁺²)) = δ
  Eq. 3  periapsis speeds     v_p± = sqrt(v∞±² + 2μ/r_p);  tangential Δv = v_p⁺ - v_p⁻
  Eq. 4  B-plane frame        Ŝ = v̂∞⁻;  T̂ = (Ŝ x k̂)/‖Ŝ x k̂‖;  R̂ = Ŝ x T̂  (k̂ = pole)
  Eq. 5  B-plane angle        θ_B = atan2(v̂∞⁺·R̂, v̂∞⁺·T̂) - π

**Published tolerances** (deep-dive §4 — parameter anchors, NOT goldens; Jones
gives no worked numeric example of Eqs. 1-5, so any unit check here is a
self-consistency check, never a sourced golden):

  - interior-flyby v∞-mismatch tolerance  Δv∞^max ∈ [100, 200] m/s  (Eq. 7)
  - flyby altitude window                 100 km - 100,000 km        (Eq. 7)

GOLDEN DISCIPLINE: nothing in this module is an EXPECTED value. Everything is a
computed quantity the shooter feeds into its residual / feasibility gate.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from cyclerfinder.core.constants import PLANETS

Vec3 = NDArray[np.float64]

# Published interior-flyby v∞-mismatch tolerance (Jones Eq. 7, deep-dive §2.3):
# "Δv∞^max between 100 and 200 m/s". We adopt the *looser* published bound (200
# m/s) as the feasibility gate (Jones: "velocity increments below 200 m/sec are
# permitted since experience has shown these can be differentially corrected in
# high-fidelity dynamics to be entirely ballistic"). km/s.
VINF_MISMATCH_TOL_KMS = 0.200

# Published flyby altitude window (Jones Eq. 7 / deep-dive §2.1): 100 km .. 1e5 km.
MIN_FLYBY_ALT_KM = 100.0
MAX_FLYBY_ALT_KM = 100_000.0

# Body spin-pole unit vector. Jones uses k̂ = (0,0,1) (body-centered equatorial,
# Eq. 4). In our heliocentric J2000-ecliptic working frame this is the ecliptic
# +Z; the B-plane angle is frame-relative and the shooter only uses θ_B as a
# targeting parameter, so the ecliptic pole is a consistent, documented choice.
_K_HAT = np.array([0.0, 0.0, 1.0])


def turn_angle_rad(vinf_in: Vec3, vinf_out: Vec3) -> float:
    """Eq. 1 — the bend δ between incoming and outgoing v∞ asymptotes (rad)."""
    a = np.asarray(vinf_in, dtype=np.float64)
    b = np.asarray(vinf_out, dtype=np.float64)
    na = float(np.linalg.norm(a))
    nb = float(np.linalg.norm(b))
    if na <= 0.0 or nb <= 0.0:
        return 0.0
    cos_d = float(np.dot(a, b) / (na * nb))
    cos_d = max(-1.0, min(1.0, cos_d))
    return float(np.arccos(cos_d))


def periapsis_radius_km(vinf_in: Vec3, vinf_out: Vec3, body: str) -> float:
    """Eq. 2 — periapsis radius r_p delivering the turn δ for a powered flyby.

    Solves ``asin(μ/(μ + r_p v∞⁻²)) + asin(μ/(μ + r_p v∞⁺²)) = δ`` for ``r_p`` by
    bisection on ``[0, r_outer]``. The left side is monotonically *decreasing* in
    ``r_p`` (larger periapsis → weaker bend), so a unique root exists for any
    achievable ``δ`` in ``(0, δ_max]``. Subsurface solutions are allowed here
    (Jones removes them later via the altitude filter); ``r_p`` may be < R_body.

    Returns ``+inf`` when ``δ`` is ~0 (no bend → any radius works; the limiting
    r_p → ∞). Returns ``0.0`` when even ``r_p = 0`` cannot deliver ``δ`` (the bend
    exceeds the body's maximum possible turn).
    """
    mu = PLANETS[body].mu_km3_s2
    delta = turn_angle_rad(vinf_in, vinf_out)
    if delta <= 1e-12:
        return float("inf")
    v_in2 = float(np.dot(vinf_in, vinf_in))
    v_out2 = float(np.dot(vinf_out, vinf_out))

    def lhs(r_p: float) -> float:
        term_in = np.arcsin(mu / (mu + r_p * v_in2)) if (mu + r_p * v_in2) > 0 else np.pi / 2
        term_out = np.arcsin(mu / (mu + r_p * v_out2)) if (mu + r_p * v_out2) > 0 else np.pi / 2
        return float(term_in + term_out)

    # At r_p = 0 the bend is maximal (π); if even that is below δ, infeasible.
    if lhs(0.0) < delta:
        return 0.0
    # Bracket: grow the upper bound until lhs falls below δ.
    hi = 1.0
    for _ in range(200):
        if lhs(hi) < delta:
            break
        hi *= 2.0
    else:
        return float("inf")
    lo = 0.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if lhs(mid) > delta:
            lo = mid
        else:
            hi = mid
        if hi - lo < 1e-6:
            break
    return 0.5 * (lo + hi)


def tangential_dv_kms(vinf_in: Vec3, vinf_out: Vec3, body: str) -> float:
    """Eq. 3 — tangential powered-flyby Δv = v_p⁺ - v_p⁻ at the Eq. 2 periapsis.

    ``v_p± = sqrt(v∞±² + 2μ/r_p)`` (energy at periapsis). The tangential maneuver
    magnitude is ``|v_p⁺ - v_p⁻|`` (a fast sub-optimal filter; Jones §2.2). When
    the turn is ballistic (r_p → ∞), this reduces to ``|v∞⁺ - v∞⁻|`` (the pure
    magnitude change), as the 2μ/r_p term vanishes.
    """
    mu = PLANETS[body].mu_km3_s2
    r_p = periapsis_radius_km(vinf_in, vinf_out, body)
    v_in = float(np.linalg.norm(vinf_in))
    v_out = float(np.linalg.norm(vinf_out))
    if not np.isfinite(r_p) or r_p <= 0.0:
        # No finite bending radius: degenerate to the magnitude change.
        return abs(v_out - v_in)
    vp_in = float(np.sqrt(v_in * v_in + 2.0 * mu / r_p))
    vp_out = float(np.sqrt(v_out * v_out + 2.0 * mu / r_p))
    return abs(vp_out - vp_in)


@dataclass(frozen=True)
class BPlaneFrame:
    """Body-centered B-plane orthonormal triad (Jones Eq. 4)."""

    s_hat: Vec3
    t_hat: Vec3
    r_hat: Vec3


def bplane_frame(vinf_in: Vec3, *, k_hat: Vec3 = _K_HAT) -> BPlaneFrame:
    """Eq. 4 — the B-plane frame Ŝ,T̂,R̂ from the incoming asymptote.

    ``Ŝ = v̂∞⁻``;  ``T̂ = (Ŝ x k̂)/‖Ŝ x k̂‖``;  ``R̂ = Ŝ x T̂``. ``k̂`` is the body
    spin pole (Jones: (0,0,1) body-centered equatorial). Degenerate when Ŝ ∥ k̂;
    we fall back to the x-axis for T̂'s seed cross-product so the triad stays
    orthonormal.
    """
    s = np.asarray(vinf_in, dtype=np.float64)
    ns = float(np.linalg.norm(s))
    if ns <= 0.0:
        raise ValueError("incoming v∞ has zero magnitude; B-plane frame undefined")
    s_hat = s / ns
    cross = np.cross(s_hat, k_hat)
    nc = float(np.linalg.norm(cross))
    if nc <= 1e-12:
        # Ŝ ∥ pole: pick an arbitrary in-plane seed orthogonal to Ŝ.
        seed = np.array([1.0, 0.0, 0.0])
        if abs(float(np.dot(s_hat, seed))) > 0.9:
            seed = np.array([0.0, 1.0, 0.0])
        cross = np.cross(s_hat, seed)
        nc = float(np.linalg.norm(cross))
    t_hat = cross / nc
    r_hat = np.cross(s_hat, t_hat)
    return BPlaneFrame(s_hat=s_hat, t_hat=t_hat, r_hat=r_hat)


def bplane_angle_rad(vinf_in: Vec3, vinf_out: Vec3, *, k_hat: Vec3 = _K_HAT) -> float:
    """Eq. 5 — the B-plane angle θ_B targeting the outgoing asymptote (rad).

    ``θ_B = atan2(v̂∞⁺·R̂, v̂∞⁺·T̂) - π`` (atan2 range (-π, π]). The flyby bends v∞
    so the projection of v∞⁺ onto the B-plane lies along the -B vector.
    """
    frame = bplane_frame(vinf_in, k_hat=k_hat)
    b = np.asarray(vinf_out, dtype=np.float64)
    nb = float(np.linalg.norm(b))
    if nb <= 0.0:
        return 0.0
    b_hat = b / nb
    return float(
        np.arctan2(float(np.dot(b_hat, frame.r_hat)), float(np.dot(b_hat, frame.t_hat))) - np.pi
    )


def flyby_altitude_km(vinf_in: Vec3, vinf_out: Vec3, body: str) -> float:
    """Periapsis altitude (r_p - R_body) of the powered flyby, km. May be < 0."""
    r_p = periapsis_radius_km(vinf_in, vinf_out, body)
    if not np.isfinite(r_p):
        return float("inf")
    return r_p - PLANETS[body].radius_eq_km


def interior_flyby_feasible(
    vinf_in: Vec3,
    vinf_out: Vec3,
    body: str,
    *,
    vinf_mismatch_tol_kms: float = VINF_MISMATCH_TOL_KMS,
) -> bool:
    """Eq. 7 — the published interior-flyby feasibility gate.

    ``‖v∞⁺ - v∞⁻‖ < Δv∞^max`` (magnitude continuity, default 200 m/s) AND the
    powered-flyby periapsis altitude lies in the published 100 km - 100,000 km
    window. Both thresholds are sourced (Jones Eq. 7); the inputs are our computed
    v∞ vectors. Pure feasibility predicate (no golden EXPECTED).
    """
    mismatch = float(np.linalg.norm(np.asarray(vinf_out) - np.asarray(vinf_in)))
    if mismatch >= vinf_mismatch_tol_kms:
        return False
    alt = flyby_altitude_km(vinf_in, vinf_out, body)
    return MIN_FLYBY_ALT_KM < alt < MAX_FLYBY_ALT_KM


__all__ = [
    "MAX_FLYBY_ALT_KM",
    "MIN_FLYBY_ALT_KM",
    "VINF_MISMATCH_TOL_KMS",
    "BPlaneFrame",
    "bplane_angle_rad",
    "bplane_frame",
    "flyby_altitude_km",
    "interior_flyby_feasible",
    "periapsis_radius_km",
    "tangential_dv_kms",
    "turn_angle_rad",
]
