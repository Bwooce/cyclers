"""VILM (V-infinity-leveraging maneuver) feasibility / ΔV-floor layer — Tier-1.

Transcribes the **phase-free** half of Campagnola & Russell, "The Endgame
Problem" Part 1 (Leveraging Graph): the n:m_K± leg taxonomy, the Γ function
(Eq. 25, Appendix A), the V̄∞-efficiency root (Eq. 9), and the theoretical-minimum
ΔV quadrature (Eq. 13). All velocities are *adimensional* (normalised by the
minor body's circular velocity V_M = sqrt(mu_primary / a_M)); a dimensional ΔV /
V∞ is recovered by multiplying by V_M.

GATE SCOPE (plan Phase 5, binding): the numeric outputs are gated by the
transcribed Endgame Part-1 Tables 1-3 (mining note A1-A3) + the worked Europa
scalar (A6). The two flagged suspect Part-2 Table 1 cells are NEVER goldens. The
Γ shape + the T = 3 - v∞² identity are physics-invariant (validated by
round-trip/algebra in Phase 4), not table-gated.

Equation references are to the mining note transcription,
``docs/notes/2026-06-05-endgame-tisserand-mining.md``.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from scipy.optimize import brentq

from cyclerfinder.core.satellites import PRIMARIES, SATELLITES

# ---------------------------------------------------------------------------
# n:m_K± leg taxonomy (mining note lines 98-130) — physics-invariant
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class VilmLeg:
    """One V∞-leveraging leg classified by the paper's n:m_K± scheme (pp.3-4).

    resonance:
        ``(n, m)`` — n (m) is the approximate number of minor-body (spacecraft)
        revolutions during the VILM.
    body:
        The minor body (moon) the leg leverages at.
    exterior:
        Exterior (ΔV_AB at apocenter, r_A > a_M) vs Interior (at pericenter).
    k_revs:
        K — number of full spacecraft revolutions on the H-B arc (default 0).
    long_transfer:
        H- (long-transfer) vs H+ (short-transfer) encounter point.
    """

    resonance: tuple[int, int]
    body: str
    exterior: bool
    k_revs: int = 0
    long_transfer: bool = True


def classify_vilm_leg(
    *,
    n: int,
    m: int,
    body: str,
    exterior: bool,
    k_revs: int = 0,
    long_transfer: bool = True,
) -> VilmLeg:
    """Classify a VILM leg into the n:m_K± taxonomy (mining note 98-130)."""
    return VilmLeg(
        resonance=(n, m),
        body=body,
        exterior=exterior,
        k_revs=k_revs,
        long_transfer=long_transfer,
    )


# ---------------------------------------------------------------------------
# Adimensional VILM physics (V_c = velocities normalised by V_M)
# ---------------------------------------------------------------------------


def gamma(vinf: float, *, exterior: bool) -> float:
    """Γ^(E,I)(V∞) — phase-free leverage gain, Eq. (25) (mining note 88-91).

    ``Γ^(E,I)(V∞) = V∞ (V∞³ ± 3V∞² - V∞ ∓ 7) / (V∞³ ± 3V∞² + V∞ ∓ 1)`` with the
    upper sign for the Exterior VILM. Adimensional (V∞ normalised by V_M).
    """
    s = 1.0 if exterior else -1.0
    num = vinf * (vinf**3 + s * 3.0 * vinf**2 - vinf - s * 7.0)
    den = vinf**3 + s * 3.0 * vinf**2 + vinf - s * 1.0
    return num / den


def _vc_adim(moon: str) -> float:
    """Adimensional flyby circular velocity V_c = sqrt(mu_M / r_π) / V_M.

    r_π = the moon's radius + its (paper-matching) safe altitude; V_M the moon's
    circular velocity about the primary. This is the V_c parameter of the Eq. (9)
    proposition (mining note 147-149).
    """
    sat = SATELLITES[moon]
    v_m = math.sqrt(PRIMARIES[sat.primary] / sat.sma_km)
    r_pi = sat.radius_eq_km + sat.safe_alt_km
    return math.sqrt(sat.mu_km3_s2 / r_pi) / v_m


def _vinf_of_vpi(vpi: float, vc: float) -> float:
    """V∞ from pericenter speed V_π, Eq. (1)/(10): V∞ = sqrt(V_π² - 2 V_c²)."""
    val = vpi * vpi - 2.0 * vc * vc
    return math.sqrt(val) if val > 0.0 else float("nan")


def _vbar_vinf_adim(moon: str, *, exterior: bool) -> float:
    """Adimensional V̄∞ — the Eq. (9) efficiency-threshold root for ``moon``.

    The VILM strategy is efficient iff V∞L > V̄∞, where V̄∞ = sqrt(V̄_π² - 2V_c²)
    and V̄_π is the root of f(V_π) = Γ∘V∞(V_π; V_c) - V_π (mining note 147-149,
    Eq. 9). Solved numerically with a bracketed brentq on V_π.
    """
    vc = _vc_adim(moon)
    lo = math.sqrt(2.0) * vc + 1e-9

    def f(vpi: float) -> float:
        vinf = _vinf_of_vpi(vpi, vc)
        if math.isnan(vinf):
            return -vpi
        return gamma(vinf, exterior=exterior) - vpi

    # f(lo) ~ Γ(0) - lo = -lo < 0; f grows past the root. Expand the upper
    # bracket until a sign change (the root sits within ~2 V_c of lo).
    hi = lo + 0.01
    f_lo = f(lo)
    while f(hi) * f_lo > 0.0 and hi < lo + 5.0:
        hi += 0.01
    vpi_bar = brentq(f, lo, hi, xtol=1e-12)
    return _vinf_of_vpi(vpi_bar, vc)


def min_vinf_for_vilm(moon: str, *, exterior: bool = True) -> float:
    """Minimum *dimensional* V∞ (km/s) at which a VILM at ``moon`` is efficient.

    The Eq. (9) root V-inf-bar (Part-1 Table 3 V-inf-bar E/I column). Adimensional
    root times V_M. ``exterior=True`` reproduces the E column, ``False`` the I.
    """
    sat = SATELLITES[moon]
    v_m = math.sqrt(PRIMARIES[sat.primary] / sat.sma_km)
    return _vbar_vinf_adim(moon, exterior=exterior) * v_m
