"""Independent integration confirmation of the Oberth periapsis flyby cost.

Task #154. Tests for :mod:`cyclerfinder.verify.flyby_integrate`, the
direct-numerical-integration check (option a) of
:func:`cyclerfinder.core.flyby.dv_powered_flyby_periapsis`.

Test layers
-----------
1. **Harness self-validation** — the integrated turn of an *unpowered*
   hyperbola reproduces the analytic ``max_bend`` cone, and the
   asymptote-baseline reconstruction reproduces ``dv_from_turn_deficit``.
   These validate the integrator + asymptote extraction before it is used to
   judge the Oberth formula (task §3).
2. **Aldrin-point agreement gate** — at V∞=6.86 km/s, required turn 93.01°
   (Earth ``rp_min`` from constants), the integrated Δv agrees with the
   formula to ≤1% relative (it agrees to ~1e-9 in practice).
3. **Sweep gate** (slow) — V∞ in {2,4,6.86,9,12} km/s by deficit in
   {5°,15°,24.52°,40°}; every point agrees to ≤1% relative, and the high-V∞
   crossover (Oberth exceeds asymptote) is confirmed to exist.

Golden discipline: the integration *is* the independent check (a different
numerical method on the same two-body physics, like the SPICE cross-check
pattern). No invented expected Δv values — the only expectations are the
model-agreement tolerance and physically-mandated facts (zero within the
cone, turn achieved == required, crossover existence).
"""

from __future__ import annotations

from math import degrees, radians, sqrt

import pytest

from cyclerfinder.core.constants import PLANETS, SAFE_PERIHELION_KM
from cyclerfinder.core.flyby import dv_from_turn_deficit, max_bend
from cyclerfinder.verify.flyby_integrate import (
    asymptote_baseline_dv,
    confirm_oberth_point,
    integrate_turn_angle,
)

MU_E = PLANETS["E"].mu_km3_s2
RP_E = SAFE_PERIHELION_KM["E"]

# Agreement gate: ≤1% relative on Δv. Observed worst case across the full
# sweep is ~3.4e-5 (a single near-escape steep-turn point); everywhere else
# ≤2e-7. 1e-2 is a deliberately loose publishable gate; the realised numbers
# are reported in the confirmation note.
REL_DV_GATE = 1.0e-2

# A reduced far field keeps the fast tests quick while still resolving the
# asymptote to <1e-6 rad even at the slowest-converging low-V∞ corner (the
# direction error scales like (2 mu / rp) / (vinf**2 * N), worst at vinf=2);
# the slow sweep uses the module default (1e5). Integration cost is ~constant
# in N (adaptive DOP853), so this only tightens accuracy, not runtime.
FAST_FAR_FIELD = 2.0e4


# ---------------------------------------------------------------------------
# Layer 1: harness self-validation (asymptote case).
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("vinf", [2.0, 4.0, 6.86, 9.0, 12.0])
def test_integrated_cone_matches_max_bend(vinf: float) -> None:
    """Unpowered hyperbola's integrated turn == analytic ``max_bend`` cone.

    This validates the integrator and the asymptote-direction extraction: if
    the harness can reproduce the closed-form cone to <1e-5 rad, its verdict
    on the powered case is trustworthy.
    """
    vp_in = sqrt(vinf * vinf + 2.0 * MU_E / RP_E)
    turn = integrate_turn_angle(vp_in, RP_E, MU_E, FAST_FAR_FIELD)
    expected = max_bend(MU_E, RP_E, vinf)
    assert turn == pytest.approx(expected, abs=1.0e-5)


@pytest.mark.parametrize("vinf", [2.0, 6.86, 12.0])
@pytest.mark.parametrize("deficit_deg", [5.0, 24.52, 40.0])
def test_asymptote_baseline_matches_formula(vinf: float, deficit_deg: float) -> None:
    """Harness asymptote reconstruction == ``dv_from_turn_deficit`` (task §3).

    The baseline burn is a magnitude-preserving rotation of the V∞ vector at
    infinity; building it from explicit vectors and a subtraction must equal
    the closed-form ``2 V∞ sin(deficit/2)``. Trivial agreement here proves the
    harness's notion of "rotate the asymptote" is the same physics the formula
    encodes — validating the harness itself.
    """
    delta_max = max_bend(MU_E, RP_E, vinf)
    delta_required = delta_max + radians(deficit_deg)
    harness = asymptote_baseline_dv(vinf, delta_required, delta_max)
    formula = dv_from_turn_deficit(vinf, delta_required, delta_max)
    assert harness == pytest.approx(formula, rel=1.0e-12, abs=1.0e-12)


def test_within_cone_is_zero() -> None:
    """A required turn inside the cone costs exactly 0 in both models."""
    vinf = 6.86
    delta_max = max_bend(MU_E, RP_E, vinf)
    delta_required = delta_max - radians(5.0)
    res = confirm_oberth_point(vinf, delta_required, MU_E, RP_E, FAST_FAR_FIELD)
    assert res.dv_integrated == 0.0
    assert res.dv_formula == 0.0


# ---------------------------------------------------------------------------
# Layer 2: Aldrin-point agreement gate.
# ---------------------------------------------------------------------------


def test_aldrin_point_agreement() -> None:
    """Aldrin Earth flyby: integrated Δv agrees with the formula to ≤1%.

    V∞=6.86 km/s, required turn 93.01° (the #151/#148 Aldrin maintenance
    geometry). The integration root-solves for the periapsis speed that turns
    by 93.01° and charges ``2|vp_in - vp_target|``; this must match
    ``dv_powered_flyby_periapsis`` for the same geometry, and the achieved turn
    must equal the required turn.
    """
    vinf = 6.86
    delta_required = radians(93.01)
    res = confirm_oberth_point(vinf, delta_required, MU_E, RP_E)

    # Turn actually delivered == required (root-solve converged).
    assert degrees(res.turn_achieved_rad) == pytest.approx(93.01, abs=1.0e-4)
    # Asymptote extraction is essentially exact at the default far field.
    assert res.asymptote_residual_rad < 1.0e-6
    # The headline confirmation: integrated vs formula Δv.
    assert res.rel_dv_error <= REL_DV_GATE
    # Sanity: the Oberth cost is well below the asymptote baseline here
    # (favorable regime — the deficit is far below the crossover).
    delta_max = max_bend(MU_E, RP_E, vinf)
    asym = asymptote_baseline_dv(vinf, delta_required, delta_max)
    assert res.dv_integrated < asym


# ---------------------------------------------------------------------------
# Layer 3: sweep gate + crossover (slow).
# ---------------------------------------------------------------------------


@pytest.mark.slow
@pytest.mark.parametrize("vinf", [2.0, 4.0, 6.86, 9.0, 12.0])
@pytest.mark.parametrize("deficit_deg", [5.0, 15.0, 24.52, 40.0])
def test_sweep_agreement(vinf: float, deficit_deg: float) -> None:
    """Full grid: integrated Δv agrees with the formula to ≤1% everywhere."""
    delta_max = max_bend(MU_E, RP_E, vinf)
    delta_required = delta_max + radians(deficit_deg)
    assert delta_required < radians(180.0)
    res = confirm_oberth_point(vinf, delta_required, MU_E, RP_E)
    assert degrees(res.turn_achieved_rad) == pytest.approx(degrees(delta_required), abs=1.0e-3)
    assert res.rel_dv_error <= REL_DV_GATE


@pytest.mark.slow
def test_high_vinf_crossover_exists() -> None:
    """Confirm the #151 high-V∞ crossover: Oberth exceeds asymptote.

    The note asserts the periapsis maneuver loses its advantage at high V∞
    (the ballistic cone is so narrow that opening it demands a large doubled
    magnitude excursion). We confirm *by integration* that at V∞=12 km/s with
    a small deficit the integrated Oberth cost exceeds the asymptote baseline,
    while at the Aldrin V∞=6.86 it does not — i.e. a crossover exists between
    the two regimes.
    """
    delta_max_hi = max_bend(MU_E, RP_E, 12.0)
    hi = confirm_oberth_point(12.0, delta_max_hi + radians(5.0), MU_E, RP_E)
    asym_hi = asymptote_baseline_dv(12.0, delta_max_hi + radians(5.0), delta_max_hi)
    assert hi.dv_integrated > asym_hi  # Oberth LOSES at high V∞ (crossover side)

    delta_max_lo = max_bend(MU_E, RP_E, 6.86)
    lo = confirm_oberth_point(6.86, delta_max_lo + radians(24.52), MU_E, RP_E)
    asym_lo = asymptote_baseline_dv(6.86, delta_max_lo + radians(24.52), delta_max_lo)
    assert lo.dv_integrated < asym_lo  # Oberth WINS at Aldrin (favorable side)
