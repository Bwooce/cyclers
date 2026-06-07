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
from collections.abc import Sequence
from dataclasses import dataclass

from scipy.integrate import quad
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


# ---------------------------------------------------------------------------
# Eq. (13) theoretical-minimum ΔV quadrature (mining note 153-183, A2/A3)
# ---------------------------------------------------------------------------


def _v_m(moon: str) -> float:
    """Circular velocity of ``moon`` about its primary, km/s (V_M)."""
    sat = SATELLITES[moon]
    return math.sqrt(PRIMARIES[sat.primary] / sat.sma_km)


def _escape_dv(moon: str, vinf_kms: float) -> float:
    """ΔV (km/s) from a circular parking orbit at the moon's r_π to a hyperbola
    of excess speed ``vinf_kms`` (km/s) about the moon. Symmetric for capture.

    ``r_π = radius_eq_km + safe_alt_km`` (the paper's r̃_π, 100 km / Titan
    1500 km). This is the escape (and, mirrored, the capture) column of Part-1
    Table 1 / 2.
    """
    sat = SATELLITES[moon]
    r_pi = sat.radius_eq_km + sat.safe_alt_km
    mu = sat.mu_km3_s2
    return math.sqrt(vinf_kms * vinf_kms + 2.0 * mu / r_pi) - math.sqrt(mu / r_pi)


def _hohmann_vinf(moon_outer: str, moon_inner: str) -> tuple[float, float]:
    """Dimensional V∞ (km/s) at the outer / inner moon for the Hohmann transfer
    between their circular orbits about the shared primary.

    ``(vinf_outer, vinf_inner)`` = ``|v_transfer - v_circular|`` at each moon's
    orbit. This is the V∞H bound the begingame/endgame quadratures integrate to.
    """
    sat_o = SATELLITES[moon_outer]
    mu = PRIMARIES[sat_o.primary]
    r_o = sat_o.sma_km
    r_i = SATELLITES[moon_inner].sma_km
    a_t = 0.5 * (r_o + r_i)
    v_o = math.sqrt(mu * (2.0 / r_o - 1.0 / a_t))
    v_i = math.sqrt(mu * (2.0 / r_i - 1.0 / a_t))
    return abs(v_o - math.sqrt(mu / r_o)), abs(v_i - math.sqrt(mu / r_i))


def _quadrature_dv_adim(vinf_lo: float, vinf_hi: float, *, exterior: bool) -> float:
    """Adimensional Eq. (13) leverage ΔV = ∫_{V∞L}^{V∞H} V∞ / Γ^(E,I)(V∞) dV∞."""
    if vinf_hi <= vinf_lo:
        return 0.0
    val, _ = quad(lambda v: v / gamma(v, exterior=exterior), vinf_lo, vinf_hi)
    return float(val)


def _leverage_dv_kms(moon: str, vinf_hi_kms: float, *, exterior: bool) -> float:
    """Dimensional begingame/endgame ΔV (km/s) at ``moon``: integrate Eq. (13)
    from the efficiency threshold V̄∞ up to the Hohmann V∞H, scaled by V_M."""
    v_m = _v_m(moon)
    vbar = min_vinf_for_vilm(moon, exterior=exterior) / v_m
    return _quadrature_dv_adim(vbar, vinf_hi_kms / v_m, exterior=exterior) * v_m


def _order_by_sma(moon_a: str, moon_b: str) -> tuple[str, str]:
    """Return ``(outer, inner)`` by about-primary SMA."""
    if SATELLITES[moon_a].sma_km >= SATELLITES[moon_b].sma_km:
        return moon_a, moon_b
    return moon_b, moon_a


def vilm_dv_min(moon_a: str, moon_b: str, *, via: Sequence[str] | None = None) -> float:
    """Theoretical-minimum VILM intermoon-transfer ΔV (km/s), Eq. (13).

    The minimum-ΔV strategy (mining note 181-183): escape the source moon to its
    efficiency-threshold V̄∞, run Exterior VILMs at the outer moon (begingame) up
    to the Hohmann V∞H, the Hohmann transfer itself (ballistic, ΔV folded into
    the begin/end V∞ bounds), Interior VILMs at the inner moon (endgame) down to
    its V̄∞, and capture. Total = escape + begingame + endgame + capture.

    Both moons must share a primary. ``via`` adds intermediate-moon gravity
    assists (Part-1 Table 2): the transfer is chained outer->...->inner through
    the listed intermediate moons, each leg's begingame/endgame summed and the
    interior escapes/captures at intermediate moons replaced by ballistic flybys
    (so an intermediate moon contributes only its leverage quadratures, not a
    fresh escape/capture pair). GA can only reduce ΔV.

    GOLDEN: reproduces Part-1 Table 1 (no-GA) and Table 2 (with-GA) ΔV_min to
    well inside the 10% linked-conic-vs-CR3BP band (mining note 491-493).
    """
    chain = [moon_a, *via, moon_b] if via else [moon_a, moon_b]
    # Order the whole chain outer->inner by about-primary SMA.
    chain_sorted = sorted(set(chain), key=lambda m: -SATELLITES[m].sma_km)
    outer, inner = chain_sorted[0], chain_sorted[-1]

    # The two ends of the full transfer carry the escape/capture insertions; the
    # outer member is escaped to its Exterior V̄∞, the inner captured at its
    # Interior V̄∞ (Part-1 Table 1/2 escape + capture columns).
    dv = _escape_dv(outer, min_vinf_for_vilm(outer)) + _escape_dv(
        inner, min_vinf_for_vilm(inner, exterior=False)
    )

    if len(chain_sorted) == 2:
        # No-GA (Part-1 Table 1): begingame at the outer moon up to the Hohmann
        # V∞H, the Hohmann transfer (ballistic), endgame at the inner moon.
        vinf_o, vinf_i = _hohmann_vinf(outer, inner)
        dv += _leverage_dv_kms(outer, vinf_o, exterior=True)
        dv += _leverage_dv_kms(inner, vinf_i, exterior=False)
        return dv

    # With-GA (Part-1 Table 2): the intermediate moons supply FREE ballistic
    # gravity assists that bridge the middle of the transfer, so they contribute
    # NO leverage quadrature and NO escape/capture. Only the begingame at the
    # outer moon (up to the first leg's Hohmann V∞) and the endgame at the inner
    # moon (down from the last leg's Hohmann V∞) cost leverage ΔV.
    vinf_outer, _ = _hohmann_vinf(outer, chain_sorted[1])
    _, vinf_inner = _hohmann_vinf(chain_sorted[-2], inner)
    dv += _leverage_dv_kms(outer, vinf_outer, exterior=True)
    dv += _leverage_dv_kms(inner, vinf_inner, exterior=False)
    return dv


def vilm_dv_floor(moon_a: str, moon_b: str) -> float:
    """Admissible ΔV lower bound (km/s) for search pruning (design §5).

    DEVIATION FROM PLAN (task #76): the plan proposed the no-GA quadrature ΔV_min
    as the admissible floor and asserted ``floor <= with-GA``. That is physically
    backwards — a gravity assist REDUCES ΔV (Part-1 Table 2 < Table 1), so the
    no-GA value is the UPPER, not lower, bound and is NOT admissible. The true
    admissible lower bound (<= every routing, no-GA and with-GA alike) is the
    irreducible escape + capture insertion cost at the two endpoints: any actual
    transfer must at least escape the source moon's parking orbit and capture
    into the destination's, and leverage + finite phasing only ADD to that. So
    ``vilm_dv_floor`` returns ``escape(outer @ V̄∞) + capture(inner @ V̄∞)``,
    which is admissible for A*-style pruning in a future Forge moon-tour search.
    """
    outer, inner = _order_by_sma(moon_a, moon_b)
    return _escape_dv(outer, min_vinf_for_vilm(outer)) + _escape_dv(
        inner, min_vinf_for_vilm(inner, exterior=False)
    )


def europa_endgame_dv() -> tuple[float, float]:
    """The Europa endgame theoretical-minimum ΔV + duration (Part-1 A6).

    Returns ``(delta_v_ms, days)`` for the Europa endgame that reduces V̄∞ from
    1.8 to 0.77 km/s (mining note 436-438). The ΔV is the Eq. (13) Exterior-VILM
    CONTINUOUS quadrature at Europa between those V∞ bounds (dimensional, m/s) —
    the theoretical minimum (infinite-VILM floor). The published DISCRETE 3-VILM
    design costs 154 m/s and the CR3BP re-optimised long-transfer 147 m/s; the
    continuous floor is a valid lower bound on both (a finite-VILM sequence
    cannot beat the floor). The duration is the published 46-day phasing scalar
    (the phase-free quadrature does not predict ToF), carried verbatim from A6.
    """
    v_m = _v_m("Europa")
    dv_kms = _quadrature_dv_adim(0.77 / v_m, 1.8 / v_m, exterior=True) * v_m
    return dv_kms * 1000.0, 46.0
