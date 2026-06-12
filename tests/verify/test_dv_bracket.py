"""Tests for the per-leg ΔV bracket diagnostic (``verify/dv_bracket.py``).

Golden discipline (binding): the EXPECTED side of every golden traces to a
published value, never to a number our own code produced. The sourced golden is
Şaloğlu & Taheri (2025) §5.1 (the self-contained geocentric Table 1 problem);
its μ_Earth is not printed in the paper, so μ is PINNED FIRST by reproducing the
published two-impulse base ΔV from the Table-1 elements with a Curtis-standard GM
BEFORE any tolerance on the bracket is set (see ``test_pin_mu_earth_*`` below).

Source
------
K. Şaloğlu, E. Taheri, "Classification and Feasibility Assessment of Infinitely
Many Iso-Impulse Three-Dimensional Trajectories," *J. Astronautical Sciences*
(2025), DOI 10.1007/s40295-025-00528-0; arXiv:2501.01583. Table 1 (p. 18) and
the two-impulse base-solution ΔV breakdown (p. 18). Mining note:
``docs/notes/2026-06-10-saloglu-2025-iso-impulse-3d-mining.md`` §5.1.

Published geocentric Table 1 (co-axial, ω = Ω = 0):
* Initial orbit: a = 7000 km, e = 0.02, i = 60°
* Target orbit:  a = 105000 km, e = 0.3, i = 12°   (48° plane change)
* Two-impulse base solution: ΔV_total = 3.9618011 km/s
  (Δv₁ = 2.8246140, Δv₂ = 1.1371871 km/s); departs initial-orbit perigee,
  arrives target apogee, plane change concentrated at the apogee impulse.
"""

from __future__ import annotations

import numpy as np
import pytest

from cyclerfinder.verify.dv_bracket import (
    BracketVerdict,
    bracket_leg_dv,
    iso_dv_split_feasible,
    lambert_ceiling_dv_kms,
)

# --- Şaloğlu & Taheri 2025, Table 1 (geocentric) ---------------------------
A1_KM, E1, I1_DEG = 7000.0, 0.02, 60.0
A2_KM, E2, I2_DEG = 105000.0, 0.3, 12.0
PUBLISHED_BASE_DV_KMS = 3.9618011  # two-impulse base ΔV_total (p. 18)

# Curtis-standard geocentric GM (Curtis, *Orbital Mechanics for Engineering
# Students*, Table A.2): 398600 km³/s²; the 4-figure value the paper's regime is
# consistent with (verified below). The trailing digits are solver-dependent.
MU_EARTH_CURTIS = 398600.0


def _bitangential_base_dv(mu: float) -> tuple[float, float, float]:
    """Bi-tangential coaxial transfer with the plane change at apogee.

    Departs the initial orbit at its perigee onto a transfer ellipse whose
    apogee is the target apogee, then circularises... no - matches the target
    orbit at its apogee with the full 48° plane change folded into the apogee
    burn. This is the *specific construction the paper names* for the two-impulse
    base (perigee departure / apogee arrival, plane change at the second burn);
    with FREE true anomalies the Eq.-1 optimiser can only do better, so this is a
    tight UPPER estimate of the published base value - exactly what is needed to
    pin μ (the published 3.9618 must sit just below it).
    """
    di = np.radians(I1_DEG - I2_DEG)
    rp1 = A1_KM * (1.0 - E1)  # initial perigee radius
    ra2 = A2_KM * (1.0 + E2)  # target apogee radius
    at = 0.5 * (rp1 + ra2)  # transfer ellipse semi-major axis

    v_init_p = np.sqrt(mu * (2.0 / rp1 - 1.0 / A1_KM))
    v_t_p = np.sqrt(mu * (2.0 / rp1 - 1.0 / at))
    dv1 = abs(v_t_p - v_init_p)  # in-plane only

    v_t_a = np.sqrt(mu * (2.0 / ra2 - 1.0 / at))
    v_tgt_a = np.sqrt(mu * (2.0 / ra2 - 1.0 / A2_KM))
    dv2 = float(np.sqrt(v_t_a**2 + v_tgt_a**2 - 2.0 * v_t_a * v_tgt_a * np.cos(di)))
    return float(dv1), float(dv2), float(dv1 + dv2)


# ---------------------------------------------------------------------------
# Pin μ FIRST: reproduce the published base ΔV regime before any tolerance.
# ---------------------------------------------------------------------------


def test_pin_mu_earth_bitangential_brackets_published_base() -> None:
    """μ is Curtis-standard: the bi-tangential base sits just ABOVE the published.

    The paper's named two-impulse base construction (perigee→apogee, plane change
    at apogee) with FIXED true anomalies is an upper estimate of its phase-FREE
    optimum. Computed with Curtis GM it gives ≈3.96594 km/s - within ~4.2 m/s
    ABOVE the published 3.9618011 km/s (the optimiser's true-anomaly freedom buys
    the small remainder). This both PINS μ ≈ 398600 km³/s² and validates the
    published number's physical regime; a wrong μ (e.g. a solar GM) misses by
    orders of magnitude.
    """
    _dv1, _dv2, total = _bitangential_base_dv(MU_EARTH_CURTIS)
    # The phase-free optimum cannot exceed the fixed-anomaly bi-tangential value.
    assert total >= PUBLISHED_BASE_DV_KMS
    # ...and the gap is small (true-anomaly freedom only), confirming Curtis μ.
    assert total - PUBLISHED_BASE_DV_KMS < 0.010  # < 10 m/s


def test_pin_mu_earth_wrong_mu_is_rejected() -> None:
    """A non-Curtis μ (solar GM) does NOT reproduce the published base regime."""
    from cyclerfinder.core.constants import MU_SUN_KM3_S2

    _dv1, _dv2, total = _bitangential_base_dv(MU_SUN_KM3_S2)
    assert abs(total - PUBLISHED_BASE_DV_KMS) > 1.0  # wildly off


# ---------------------------------------------------------------------------
# Sourced golden: the certificate ordering base ≤ dsm ≤ lambert.
# ---------------------------------------------------------------------------


def test_bracket_ordering_holds_on_published_base() -> None:
    """A DSM solution inside [published base, a Lambert ceiling] is WITHIN_BRACKET.

    Uses the SOURCED published base floor (3.9618011 km/s) as the lower bound.
    The Lambert ceiling here is a representative feasible upper bound above the
    base (a real per-leg ceiling comes from ``lambert_ceiling_dv_kms``); a DSM
    value between them must order-check clean.
    """
    base = PUBLISHED_BASE_DV_KMS
    ceiling = base + 0.50  # a feasible two-impulse ceiling above the base floor
    dsm = base + 0.20  # a solution comfortably inside the band
    res = bracket_leg_dv(base, dsm, ceiling)
    assert res.verdict is BracketVerdict.WITHIN_BRACKET
    assert 0.0 < res.headroom_fraction < 1.0
    # Headroom is measured toward the floor: (ceiling - dsm)/(ceiling - base).
    assert res.headroom_fraction == pytest.approx((ceiling - dsm) / (ceiling - base))


def test_bracket_flags_near_lambert_ceiling() -> None:
    """A DSM solution hugging the Lambert ceiling is flagged for unexploited headroom."""
    base = PUBLISHED_BASE_DV_KMS
    ceiling = base + 0.50
    dsm = ceiling - 0.01  # only 2% of the span below the ceiling
    res = bracket_leg_dv(base, dsm, ceiling, near_ceiling_frac=0.05)
    assert res.verdict is BracketVerdict.NEAR_LAMBERT_CEILING
    assert res.headroom_fraction < 0.05


def test_bracket_flags_below_base_floor() -> None:
    """A DSM ΔV below the phase-free floor is impossible → BELOW_BASE_FLOOR (bug)."""
    base = PUBLISHED_BASE_DV_KMS
    res = bracket_leg_dv(base, base - 0.10, base + 0.50)
    assert res.verdict is BracketVerdict.BELOW_BASE_FLOOR


def test_bracket_flags_above_lambert_ceiling() -> None:
    """A DSM ΔV above the ballistic Lambert ceiling → ABOVE_LAMBERT_CEILING (bug)."""
    base = PUBLISHED_BASE_DV_KMS
    ceiling = base + 0.50
    res = bracket_leg_dv(base, ceiling + 0.10, ceiling)
    assert res.verdict is BracketVerdict.ABOVE_LAMBERT_CEILING


def test_bracket_degenerate_band_is_within() -> None:
    """When Lambert ≈ base (no slack) the leg is WITHIN_BRACKET with zero headroom."""
    base = PUBLISHED_BASE_DV_KMS
    res = bracket_leg_dv(base, base, base)
    assert res.verdict is BracketVerdict.WITHIN_BRACKET
    assert res.headroom_fraction == 0.0


# ---------------------------------------------------------------------------
# Lambert ceiling utility - construction invariants.
# ---------------------------------------------------------------------------


def test_lambert_ceiling_ballistic_arc_is_nonnegative_and_finite() -> None:
    """The Lambert ceiling is a finite, non-negative two-impulse ΔV.

    A heliocentric Earth→Mars-like geometry: positions ~1 and ~1.5 AU apart over
    a ~200-day arc; the ceiling sums the departure and arrival burns against the
    imposed boundary velocities. Construction invariant only (no sourced
    magnitude): finite, ≥ 0.
    """
    from cyclerfinder.core.constants import AU_KM, MU_SUN_KM3_S2, SECONDS_PER_DAY

    r1 = np.array([AU_KM, 0.0, 0.0])
    r2 = np.array([0.0, 1.5 * AU_KM, 0.0])
    # Boundary velocities: rough circular speeds at each radius.
    v1_before = np.array([0.0, np.sqrt(MU_SUN_KM3_S2 / AU_KM), 0.0])
    v2_after = np.array([-np.sqrt(MU_SUN_KM3_S2 / (1.5 * AU_KM)), 0.0, 0.0])
    tof = 200.0 * SECONDS_PER_DAY
    dv = lambert_ceiling_dv_kms(r1, v1_before, r2, v2_after, tof, mu=MU_SUN_KM3_S2)
    assert np.isfinite(dv)
    assert dv >= 0.0


def test_lambert_ceiling_rejects_nonpositive_tof() -> None:
    r1 = np.array([1.0e8, 0.0, 0.0])
    r2 = np.array([0.0, 1.0e8, 0.0])
    v = np.zeros(3)
    with pytest.raises(ValueError, match="TOF must be positive"):
        lambert_ceiling_dv_kms(r1, v, r2, v, 0.0)


# ---------------------------------------------------------------------------
# Eq.-16 iso-ΔV split-feasibility probe - sourced golden (Şaloğlu 2023 §5.1)
# ---------------------------------------------------------------------------
#
# Earth-Dionysus benchmark (Şaloğlu, Taheri & Landau 2023, Sec. III.A; mining
# note ``docs/notes/2026-06-10-saloglu-2023-iso-impulse-mining.md`` §5.1):
#   T_E (initial-orbit period, = T(alpha=0)) = 365.25 d
#   T_pf (phase-free-arc period, = T(alpha=1)) = 1161.47 d
#   surplus TOF = 2490.48 d
# Published feasible families include M=0, N₁=3 with phasing period 830.16 d
# (Fig. 6a) - which lies inside (365.25, 1161.47). The §III.E revolution-count
# bound is n_p < TOF/T(alpha=0) = 2490.48/365.25 = 6.82 ⇒ n_p ≤ 6.

_DAY = 86400.0
_ED_T_E_S = 365.25 * _DAY
_ED_T_PF_S = 1161.47 * _DAY
_ED_SURPLUS_S = 2490.48 * _DAY


def test_iso_split_earth_dionysus_3rev_feasible_matches_published_period() -> None:
    """SOURCED: Earth-Dionysus N₁=3 split is feasible; period matches Fig. 6a.

    The published 3-impulse family (M=0, N₁=3) has phasing period 830.16 d, which
    the Eq.-16 probe must (a) report feasible and (b) reproduce as the required
    period - a number the paper prints, so it is a sourced golden.
    """
    res = iso_dv_split_feasible(_ED_SURPLUS_S, _ED_T_E_S, _ED_T_PF_S, n_phasing_revs=3)
    assert res.feasible is True
    # Published phasing period for the (M=0, N₁=3) family: 830.16 d.
    assert res.required_phasing_period_s / _DAY == pytest.approx(830.16, abs=0.01)


def test_iso_split_earth_dionysus_revcount_bound_matches_published() -> None:
    """SOURCED: the §III.E phasing-rev upper bound is n_p ≤ 6 (2490.48/365.25)."""
    res = iso_dv_split_feasible(_ED_SURPLUS_S, _ED_T_E_S, _ED_T_PF_S, n_phasing_revs=1)
    assert res.max_phasing_revs == 6  # floor(6.82)


def test_iso_split_high_revcount_exceeds_window_infeasible() -> None:
    """Too many phasing revs shrinks the required period below T(alpha=0) → infeasible.

    With Σ N_k = 7 the required period is 2490.48/7 = 355.78 d < T_E = 365.25 d,
    i.e. below the window floor - the §III.E reason the count caps at 6.
    """
    res = iso_dv_split_feasible(_ED_SURPLUS_S, _ED_T_E_S, _ED_T_PF_S, n_phasing_revs=7)
    assert res.feasible is False
    assert res.required_phasing_period_s < _ED_T_E_S


def test_iso_split_short_leg_infeasible() -> None:
    """A short (< 1 rev) leg admits no phasing orbit (note §3): not improvable by split.

    Surplus far below the initial-orbit period: no phasing orbit fits, so the
    iso-ΔV split mechanism is unavailable (and a single interior impulse is not
    leaving anchor-type savings on the table).
    """
    res = iso_dv_split_feasible(0.3 * _ED_T_E_S, _ED_T_E_S, _ED_T_PF_S, n_phasing_revs=1)
    assert res.feasible is False
    assert res.max_phasing_revs == 0


def test_iso_split_rejects_zero_revs() -> None:
    with pytest.raises(ValueError, match="n_phasing_revs must be >= 1"):
        iso_dv_split_feasible(_ED_SURPLUS_S, _ED_T_E_S, _ED_T_PF_S, n_phasing_revs=0)
