# src/cyclerfinder/search/leveraging_leg.py
"""Phase-full VILM single-leg evaluator (spec 2026-06-09, Component 1).

The phase-FULL counterpart of the phase-free :mod:`cyclerfinder.search.vilm`
(Campagnola & Russell, "The Endgame Problem"). A VILM leg departs a moon M with
V∞_in on an orbit resonant with M, applies a deep-space impulse at the apse, and
returns to M with a changed V∞_out — the apse impulse IS the leveraging maneuver.

Canonical units about the primary (see plan "Canonical units"). Coplanar, circular
moon. Pure: math/scipy + core.satellites + search.vilm only.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from cyclerfinder.core.satellites import PRIMARIES, SATELLITES
from cyclerfinder.search import vilm


def v_m_kms(moon: str) -> float:
    """Moon circular velocity about its primary, km/s (the canonical V_M)."""
    sat = SATELLITES[moon]
    return math.sqrt(PRIMARIES[sat.primary] / sat.sma_km)


def resonant_sma(n_moon_revs: int, m_sc_revs: int) -> float:
    """Canonical SC semimajor axis for an n:m resonance (a = (n/m)**(2/3), a_M=1)."""
    return float((n_moon_revs / m_sc_revs) ** (2.0 / 3.0))


def tisserand_vinf(*, a: float, e: float) -> float:
    """Adimensional V∞ at the moon for an SC orbit (a, e). nan if no real V∞."""
    val = 3.0 - 1.0 / a - 2.0 * math.sqrt(a * (1.0 - e * e))
    if abs(val) < 1e-15:
        return 0.0
    return math.sqrt(val) if val > 0.0 else float("nan")


def eccentricity_from_vinf(*, a: float, vinf: float) -> float:
    """Eccentricity of the orbit with semimajor axis ``a`` and V∞ ``vinf``.

    Inverts :func:`tisserand_vinf`. nan if no real moon-crossing orbit.
    """
    h = (3.0 - 1.0 / a - vinf * vinf) / 2.0
    if h < 0.0:
        return float("nan")
    ratio = (h * h) / a
    if ratio > 1.0:
        return float("nan")
    e = math.sqrt(1.0 - ratio)
    # The leg is a moon flyby at r=a_M=1, so the orbit MUST cross r=1
    # (periapsis <= 1 <= apoapsis). An orbit that never reaches the moon hosts no
    # encounter and no leg — do NOT remove this check (orbit-closure-discipline).
    if not (a * (1.0 - e) <= 1.0 + 1e-12 and a * (1.0 + e) + 1e-12 >= 1.0):
        return float("nan")
    return e


@dataclass(frozen=True)
class LeveragingLegResult:
    """One phase-full VILM leg. CONSTRAINED vs EMERGED separated (golden rule).

    FIDELITY CAVEAT (this rung): ``converged=True`` asserts only that the leg is a
    geometrically valid V∞-SHAPING burn (real near-root, bound post-burn orbit,
    both orbits cross the moon) and that ΔV ≥ the Γ floor. It does NOT assert
    phasing closure — ``resonance_residual`` is REPORTED as a diagnostic but is
    deliberately not enforced here (the post-burn orbit need not land on an
    integer resonance). True return-to-moon phasing is confirmed downstream by the
    n-body step (:func:`cyclerfinder.data.discover_novel.endgame_route_to_nbody_request`).
    So a converged leg / route is a conservative V∞-lowering lower bound, not yet a
    phasing-closed endgame design.
    """

    dv_dsm_kms: float  # CONSTRAINED — the apse leveraging burn
    vinf_out_kms: float  # EMERGED — achieved excess speed at return
    vinf_in_kms: float
    resonance: tuple[int, int]  # (n_moon_revs, m_sc_revs)
    apse_radius_km: float  # EMERGED — burn-point radius about the primary
    exterior: bool  # apoapsis (True) vs periapsis (False) burn
    moon: str
    resonance_residual: float  # phasing DIAGNOSTIC (rad), NOT enforced — see caveat
    converged: bool
    gamma_floor_ok: bool = False


def _nan_leg(moon: str, n: int, m: int, vinf_in: float, exterior: bool) -> LeveragingLegResult:
    """A non-converged result (infeasible geometry) — returned, never raised."""
    return LeveragingLegResult(
        dv_dsm_kms=float("nan"),
        vinf_out_kms=float("nan"),
        vinf_in_kms=vinf_in,
        resonance=(n, m),
        apse_radius_km=float("nan"),
        exterior=exterior,
        moon=moon,
        resonance_residual=float("nan"),
        converged=False,
        gamma_floor_ok=False,
    )


def gamma_floor_kms(*, moon: str, vinf_lo_kms: float, vinf_hi_kms: float, exterior: bool) -> float:
    """Γ-quadrature analytic-minimum ΔV (km/s) for one leg's V∞ step.

    The independent cross-check: ∫ V∞/Γ dV∞ over [lo, hi], adimensional, scaled by
    V_M. Reuses the phase-free :mod:`cyclerfinder.search.vilm` (a different code
    path — closed-form quadrature vs the apse-DSM near-root solve). A realised leg
    ΔV below this is non-physical.
    """
    v_m = v_m_kms(moon)
    lo, hi = sorted((vinf_lo_kms / v_m, vinf_hi_kms / v_m))
    return vilm._quadrature_dv_adim(lo, hi, exterior=exterior) * v_m


def evaluate_leveraging_leg(
    *,
    moon: str,
    n_moon_revs: int,
    m_sc_revs: int,
    vinf_in_kms: float,
    vinf_out_target_kms: float,
    exterior: bool,
    epoch_sec: float = 0.0,
) -> LeveragingLegResult:
    """Evaluate one phase-full VILM leg at ``moon`` (canonical core, km/s out).

    Pre-burn orbit: resonant a = (n/m)**(2/3); e from V∞_in (Tisserand). A
    tangential burn at the fixed apse radius R collapses the Tisserand relation to
    a quadratic in the apse speed ``v``::

        V∞**2 = v**2 - 2*R*v + (3 - 2/R)

    so the post-burn apse speed solving for the target V∞_out is
    ``v' = R ± sqrt(R**2 - (3 - 2/R - V∞_out**2))``. We take the **near root**
    (continuous from the pre-burn ``v``) — that is the small, *leveraged* burn; the
    far root is the unphysical flip-to-the-other-orbit solution. ΔV = |v' - v|.
    Returns ``converged=False`` on any infeasible geometry (no real root, unbound
    post-burn orbit, or an orbit that no longer crosses the moon) rather than
    raising. ``epoch_sec`` is carried for later trajectory realisation; the ΔV/V∞
    are phase-free (the paper notes the leverage ΔV is ~phase-independent).
    """
    v_m = v_m_kms(moon)
    vin = vinf_in_kms / v_m
    vout_t = vinf_out_target_kms / v_m

    a = resonant_sma(n_moon_revs, m_sc_revs)
    e = eccentricity_from_vinf(a=a, vinf=vin)
    if math.isnan(e):
        return _nan_leg(moon, n_moon_revs, m_sc_revs, vinf_in_kms, exterior)

    # Burn at the apoapsis (exterior VILM) or periapsis (interior); R held fixed.
    r_apse = a * (1.0 + e) if exterior else a * (1.0 - e)
    v_apse = math.sqrt(2.0 / r_apse - 1.0 / a)

    # V∞_out**2 = v'**2 - 2 R v' + (3 - 2/R); solve, take the near (leveraged) root.
    disc = r_apse * r_apse - (3.0 - 2.0 / r_apse - vout_t * vout_t)
    if disc < 0.0:
        # Target V∞ unreachable from this apse in a single leg.
        return _nan_leg(moon, n_moon_revs, m_sc_revs, vinf_in_kms, exterior)
    s = math.sqrt(disc)
    root_lo, root_hi = r_apse - s, r_apse + s
    v_apse_p = root_lo if abs(root_lo - v_apse) <= abs(root_hi - v_apse) else root_hi
    if v_apse_p <= 0.0:
        return _nan_leg(moon, n_moon_revs, m_sc_revs, vinf_in_kms, exterior)

    inv_aprime = 2.0 / r_apse - v_apse_p * v_apse_p
    if inv_aprime <= 0.0:
        # Post-burn orbit unbound about the primary — not a closing leg.
        return _nan_leg(moon, n_moon_revs, m_sc_revs, vinf_in_kms, exterior)
    aprime = 1.0 / inv_aprime
    eprime = abs(r_apse / aprime - 1.0)
    if eprime >= 1.0:
        return _nan_leg(moon, n_moon_revs, m_sc_revs, vinf_in_kms, exterior)
    if not (aprime * (1.0 - eprime) <= 1.0 + 1e-12 and aprime * (1.0 + eprime) + 1e-12 >= 1.0):
        # Post-burn orbit no longer crosses the moon — no return encounter.
        return _nan_leg(moon, n_moon_revs, m_sc_revs, vinf_in_kms, exterior)
    vout = tisserand_vinf(a=aprime, e=eprime)
    if math.isnan(vout):
        return _nan_leg(moon, n_moon_revs, m_sc_revs, vinf_in_kms, exterior)

    dv_kms = abs(v_apse_p - v_apse) * v_m
    # Phasing diagnostic: post-burn period vs nearest integer re-encounter (the
    # real closure check is the n-body step; this is a cheap proxy).
    revs = aprime**1.5  # canonical T_SC / T_M
    residual = abs(revs - round(revs)) * 2.0 * math.pi
    floor = gamma_floor_kms(
        moon=moon, vinf_lo_kms=vout * v_m, vinf_hi_kms=vinf_in_kms, exterior=exterior
    )
    sma_km = SATELLITES[moon].sma_km
    return LeveragingLegResult(
        dv_dsm_kms=dv_kms,
        vinf_out_kms=vout * v_m,
        vinf_in_kms=vinf_in_kms,
        resonance=(n_moon_revs, m_sc_revs),
        apse_radius_km=r_apse * sma_km,
        exterior=exterior,
        moon=moon,
        resonance_residual=residual,
        converged=True,
        gamma_floor_ok=dv_kms >= floor - 1e-9,
    )
