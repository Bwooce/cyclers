"""Analytic flyby-continuity constraint gradients (Ellison et al. 2018, #189).

Closed-form gradients of the two patched-conic flyby-continuity constraints used
by two-sided-shooting transcriptions, w.r.t. the Cartesian components of the
incoming/outgoing v∞ vectors:

  Eq. 3   v∞-magnitude continuity   c_v∞ = ‖v∞⁺‖ - ‖v∞⁻‖ = 0
  Eq. 4   minimum-altitude turn     c_alt = r_periapse - (r_body + h_safe) ≥ 0
          with  r_periapse = (μ/‖v∞⁺‖²)·[1/sin(δ/2) - 1]
  Eq. 5   turn angle                δ = acos( v∞⁻·v∞⁺ / (‖v∞⁻‖‖v∞⁺‖) )

Source: D. H. Ellison, B. A. Conway, J. A. Englander, M. T. Ozimek, "Analytic
Gradient Computation for Bounded-Impulse Trajectory Models Using Two-Sided
Shooting," *Journal of Guidance, Control, and Dynamics*, Vol. 41, No. 7, 2018,
pp. 1449-1462, doi:10.2514/1.G003077 — constraints Eqs. 3-5, gradients Appendix
Eqs. A1-A6 (mining note ``docs/notes/2026-06-10-ellison-2018-analytic-gradients-
mining.md`` §4; this is the "immediately portable regardless of path" item).

TRANSCRIPTION NOTE (scaling): this module implements the gradient of Eq. 4
exactly as printed (c_alt in km). The published Appendix expressions match this
gradient term by term — same numerator vectors (A1-A3: ``gamma v∞⁺ - φ v∞⁻``;
A4-A6: ``ξ v∞⁻ - φ v∞⁺``), same ``cos(acos alpha / 2)``, ``(alpha - 1)``,
``[1 - φ²/(gammaβ)]^{1/2}`` factors and ``gamma^{3/2}β^{3/2}`` / ``psi^{1/2}ξ^{5/2}``
powers — but each printed term carries an additional ``1/r_periapse`` (A1-A3,
first terms of A4-A6) or ``1/r_flyby`` (final terms of A4-A6) factor, i.e. the
appendix publishes the gradient of an r-scaled (nondimensionalised) form of the
constraint, with the scaling radius printed inconsistently between terms. The
unscaled gradient implemented here is the dimensionally consistent ∂c_alt/∂v∞
(km per km/s = s) and is validated against central differences — the paper's own
recommended verification pattern for derivative code (Sec. VI).

Both constraints are purely **algebraic in the v∞ vectors** — no STM chain, no
propagation enters these gradients (the structural point of the two-sided
transcription: both phase endpoints are pinned to bodies through v∞ decision
variables). They plug into any corrector that carries v∞ vectors as decision
variables (the Path-B MGAnDSMs re-transcription of the DSM lane; mining note
§7). They are provided as separate functions — an analytic OPTION alongside the
existing finite-difference machinery, not a silent behaviour change to any
solver.

Frame: any inertial Cartesian frame, applied consistently to both vectors (the
constraints are frame-invariant scalars). Units: km, km/s, km³/s².

GOLDEN DISCIPLINE: Ellison publishes no unit-level numeric gradient values (no
worked example with numbers; mining note §6), so the validation tests are
FD-vs-analytic CONSISTENCY checks, never sourced goldens.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from cyclerfinder.core.constants import PLANETS

Vec3 = NDArray[np.float64]

# Below this value of sin(δ) the turn geometry is degenerate for differentiation:
# δ ≈ 0 (no bend: r_periapse → ∞, the altitude constraint is inactive and its
# gradient is zero in the limit) or δ ≈ π (head-on reversal: δ(v∞) has a
# non-smooth cusp; the directional derivative is unbounded). Both cases return
# zero altitude gradients, with ``altitude_active`` distinguishing them.
_DEGENERATE_SIN_DELTA = 1e-12


@dataclass(frozen=True)
class FlybyContinuityGradients:
    """Flyby-continuity constraint values + analytic v∞ gradients (Eqs. 3-5, A1-A6).

    ``c_vinf_kms`` is the Eq. 3 magnitude-continuity defect (km/s; zero when the
    unpowered flyby conserves v∞ magnitude). ``c_altitude_km`` is the Eq. 4
    minimum-altitude constraint value (km; feasible when ≥ 0). The ``d_*`` fields
    are the gradients w.r.t. the Cartesian components of v∞⁻ (``vinf_in``) and
    v∞⁺ (``vinf_out``). ``r_periapse_km`` is the Eq. 4 periapsis radius
    delivering the turn (``inf`` when the bend is negligible);
    ``turn_angle_rad`` is δ (Eq. 5). ``altitude_active`` is False when δ ≈ 0
    (no periapsis bound; c_altitude is ``+inf`` and its gradients are zero).
    """

    c_vinf_kms: float
    d_cvinf_d_vinf_in: Vec3
    d_cvinf_d_vinf_out: Vec3
    c_altitude_km: float
    d_calt_d_vinf_in: Vec3
    d_calt_d_vinf_out: Vec3
    r_periapse_km: float
    turn_angle_rad: float
    altitude_active: bool


def vinf_continuity_gradient(vinf_in: Vec3, vinf_out: Vec3) -> tuple[float, Vec3, Vec3]:
    """Eq. 3 constraint + analytic gradient: ``c = ‖v∞⁺‖ - ‖v∞⁻‖``.

    Returns ``(c, ∂c/∂v∞⁻, ∂c/∂v∞⁺)``. The gradients are the (signed) unit
    vectors ``∂c/∂v∞⁺ = v̂∞⁺``, ``∂c/∂v∞⁻ = -v̂∞⁻`` (Ellison Appendix intro:
    "Gradient of Eq. (3) is trivial"). Raises ``ValueError`` on a zero-magnitude
    input (the norm is not differentiable at the origin).
    """
    v_in = np.asarray(vinf_in, dtype=np.float64)
    v_out = np.asarray(vinf_out, dtype=np.float64)
    n_in = float(np.linalg.norm(v_in))
    n_out = float(np.linalg.norm(v_out))
    if n_in <= 0.0 or n_out <= 0.0:
        raise ValueError("v-infinity continuity gradient undefined for zero-magnitude v-infinity")
    return n_out - n_in, -v_in / n_in, v_out / n_out


def flyby_altitude_gradient(
    vinf_in: Vec3,
    vinf_out: Vec3,
    body: str,
    *,
    h_safe_km: float | None = None,
) -> tuple[float, Vec3, Vec3, float, float, bool]:
    """Eq. 4 constraint + analytic gradient (Appendix Eqs. A1-A6, unscaled).

    ``c_alt = (μ/‖v∞⁺‖²)·[1/sin(δ/2) - 1] - (r_body + h_safe)`` with δ from
    Eq. 5; feasible (the turn is ballistically realizable above the safe
    altitude) when ``c_alt ≥ 0``. ``h_safe_km`` defaults to the body's
    ``safe_alt_km`` from :data:`~cyclerfinder.core.constants.PLANETS`.

    Returns ``(c_alt, ∂c/∂v∞⁻, ∂c/∂v∞⁺, r_periapse, δ, active)``.

    Closed form (scalars per the Appendix: ``φ = v∞⁻·v∞⁺``, ``gamma = psi = ‖v∞⁻‖²``,
    ``β = ξ = ‖v∞⁺‖²``, ``alpha = eps = φ/√(gammaξ) = cos δ``; ``s = sin(δ/2)`` so that
    ``1 - alpha = 2s²`` and ``√(1-alpha²) = sin δ``):

      ∂c/∂v∞⁻ = k · (gamma v∞⁺ - φ v∞⁻) / (gamma^{3/2} ξ^{3/2})                 (A1-A3)
      ∂c/∂v∞⁺ = -(2μ/ξ²)(1/s - 1) v∞⁺
                + k · (ξ v∞⁻ - φ v∞⁺) / (gamma^{1/2} ξ^{5/2})               (A4-A6)
      k = μ cos(δ/2) / (2 s² sin δ)

    (the printed A1-A6 are these expressions divided by r_periapse / r_flyby —
    their NLP constraint scaling; see the module docstring).

    Degenerate geometry: δ ≈ 0 returns ``c_alt = +inf``, zero gradients,
    ``active=False`` (no periapsis bound). δ ≈ π returns the finite constraint
    value with zero gradients (δ(v∞) has a cusp there; documented, not raised).
    Raises ``ValueError`` on a zero-magnitude input.
    """
    p = PLANETS[body]
    mu = p.mu_km3_s2
    r_flyby = p.radius_eq_km + (p.safe_alt_km if h_safe_km is None else float(h_safe_km))

    v_in = np.asarray(vinf_in, dtype=np.float64)
    v_out = np.asarray(vinf_out, dtype=np.float64)
    gamma = float(np.dot(v_in, v_in))  # ‖v∞⁻‖² (gamma in A1-A3, psi in A4-A6)
    xi = float(np.dot(v_out, v_out))  # ‖v∞⁺‖² (β in A1-A3, ξ in A4-A6)
    if gamma <= 0.0 or xi <= 0.0:
        raise ValueError("flyby altitude gradient undefined for zero-magnitude v-infinity")
    phi = float(np.dot(v_in, v_out))
    alpha = phi / float(np.sqrt(gamma * xi))  # cos δ
    alpha = max(-1.0, min(1.0, alpha))
    delta = float(np.arccos(alpha))
    sin_delta = float(np.sqrt(max(0.0, 1.0 - alpha * alpha)))
    zero3 = np.zeros(3, dtype=np.float64)

    if sin_delta <= _DEGENERATE_SIN_DELTA and alpha > 0.0:
        # δ ≈ 0: no bend — any periapsis works (r_p → ∞), constraint inactive.
        return float("inf"), zero3, zero3.copy(), float("inf"), delta, False

    s = float(np.sin(0.5 * delta))
    r_p = (mu / xi) * (1.0 / s - 1.0)
    c_alt = r_p - r_flyby

    if sin_delta <= _DEGENERATE_SIN_DELTA:
        # δ ≈ π: head-on reversal — δ(v∞) is non-smooth (cusp); the constraint
        # value is finite but its gradient is unbounded. Zero by convention.
        return c_alt, zero3, zero3.copy(), r_p, delta, True

    # k = ∂r_p/∂δ · (-1/sin δ) applied through ∂δ/∂alpha = -1/√(1-alpha²):
    # ∂r_p/∂δ = -(μ/ξ)·cos(δ/2)/(2s²) < 0 (sharper bend → lower periapsis).
    k = mu * float(np.cos(0.5 * delta)) / (2.0 * s * s * sin_delta)
    d_in = k * (gamma * v_out - phi * v_in) / (gamma**1.5 * xi**1.5)
    d_out = (-2.0 * mu / (xi * xi)) * (1.0 / s - 1.0) * v_out + k * (xi * v_in - phi * v_out) / (
        float(np.sqrt(gamma)) * xi**2.5
    )
    return c_alt, d_in, d_out, r_p, delta, True


def flyby_continuity_gradients(
    vinf_in: Vec3,
    vinf_out: Vec3,
    body: str,
    *,
    h_safe_km: float | None = None,
) -> FlybyContinuityGradients:
    """Both flyby-continuity constraints + analytic gradients (Eqs. 3-5, A1-A6).

    The combined evaluation a v∞-decision-variable corrector consumes: one call
    per interior flyby yields the two constraint rows and their exact Jacobian
    entries w.r.t. the six v∞ components — no finite-difference re-propagation.
    """
    c_v, dv_in, dv_out = vinf_continuity_gradient(vinf_in, vinf_out)
    c_a, da_in, da_out, r_p, delta, active = flyby_altitude_gradient(
        vinf_in, vinf_out, body, h_safe_km=h_safe_km
    )
    return FlybyContinuityGradients(
        c_vinf_kms=c_v,
        d_cvinf_d_vinf_in=dv_in,
        d_cvinf_d_vinf_out=dv_out,
        c_altitude_km=c_a,
        d_calt_d_vinf_in=da_in,
        d_calt_d_vinf_out=da_out,
        r_periapse_km=r_p,
        turn_angle_rad=delta,
        altitude_active=active,
    )


__all__ = [
    "FlybyContinuityGradients",
    "flyby_altitude_gradient",
    "flyby_continuity_gradients",
    "vinf_continuity_gradient",
]
