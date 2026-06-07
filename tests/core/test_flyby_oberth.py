"""Tests for :func:`cyclerfinder.core.flyby.dv_powered_flyby_periapsis`.

The Oberth-credited periapsis powered-flyby cost (Takao 2025 Eq. 11 / Russell
2004 Eq. 5.5 family). These are *mathematics* (label-mechanics) checks of the
function's stated properties — not goldens:

* (a) zero deficit -> exactly ``0.0``;
* (b) Oberth cost ``<=`` asymptote-rotation cost for the same deficit, strict
  for nonzero deficit, **in the deep-well regime** (asserted over a bounded
  sweep restricted to that regime — the guarantee is not universal, see the
  function docstring and ``docs/notes/2026-06-07-oberth-flyby-recost.md``);
* (c) monotone non-decreasing in the required turn;
* (d) consistent units/frames with :func:`dv_from_turn_deficit`.
"""

from __future__ import annotations

from math import pi, radians

import numpy as np
import pytest

from cyclerfinder.core.constants import PLANETS, SAFE_PERIHELION_KM
from cyclerfinder.core.flyby import (
    dv_from_turn_deficit,
    dv_powered_flyby_periapsis,
    max_bend,
)

MU_EARTH: float = PLANETS["E"].mu_km3_s2
RP_EARTH: float = PLANETS["E"].radius_eq_km + 200.0  # sourced Aldrin 200 km flyby
MU_MARS: float = PLANETS["M"].mu_km3_s2
RP_MARS: float = SAFE_PERIHELION_KM["M"]


# --- (a) zero deficit -> zero cost ----------------------------------------


def test_zero_deficit_is_exactly_zero() -> None:
    vinf = 6.0
    dmax = max_bend(MU_EARTH, RP_EARTH, vinf)
    # required turn within the cone -> no maneuver.
    assert dv_powered_flyby_periapsis(vinf, 0.5 * dmax, dmax, MU_EARTH, RP_EARTH) == 0.0
    # exactly at the cone edge -> still zero.
    assert dv_powered_flyby_periapsis(vinf, dmax, dmax, MU_EARTH, RP_EARTH) == 0.0


def test_zero_vinf_zero_cost() -> None:
    # vinf == 0 -> cone is pi, no finite deficit possible; cost is 0.
    dmax = max_bend(MU_EARTH, RP_EARTH, 0.0)
    assert dv_powered_flyby_periapsis(0.0, dmax, dmax, MU_EARTH, RP_EARTH) == 0.0


# --- (c) monotone non-decreasing in required turn -------------------------


@pytest.mark.parametrize("vinf", [2.0, 4.0, 6.0, 6.8615])
def test_monotone_in_required_turn(vinf: float) -> None:
    dmax = max_bend(MU_EARTH, RP_EARTH, vinf)
    turns = np.linspace(dmax, min(pi - 1e-3, dmax + radians(80.0)), 60)
    costs = [dv_powered_flyby_periapsis(vinf, t, dmax, MU_EARTH, RP_EARTH) for t in turns]
    diffs = np.diff(costs)
    assert np.all(diffs >= -1e-9), f"non-monotone: min step {diffs.min():.3e}"
    # strictly positive once past the cone.
    assert costs[-1] > 0.0


# --- (b) Oberth <= asymptote in the deep-well regime ----------------------
#
# The two models are only ordered Oberth <= asymptote when the well dominates
# (2 mu / rp >> vinf^2). For an Earth 200 km flyby that holds up to
# vinf ~ 6.9 km/s across the whole deficit range; for Mars (shallower well) the
# threshold is lower. We sweep the deep-well regime and assert the strict
# ordering, then separately document (test below) that it is NOT universal.


@pytest.mark.parametrize(
    ("mu", "rp", "vinf_max"),
    [
        (MU_EARTH, RP_EARTH, 6.5),  # Earth deep-well band (covers Aldrin 6.86 nearby)
        (MU_MARS, RP_MARS, 2.0),  # Mars shallow well -> tighter band
    ],
)
def test_oberth_le_asymptote_deep_well_sweep(mu: float, rp: float, vinf_max: float) -> None:
    n_checked = 0
    for vinf in np.linspace(0.5, vinf_max, 12):
        dmax = max_bend(mu, rp, float(vinf))
        for extra_deg in (1.0, 5.0, 10.0, 24.5, 40.0):
            dreq = dmax + radians(extra_deg)
            if dreq >= pi - 1e-3:
                continue
            oberth = dv_powered_flyby_periapsis(float(vinf), dreq, dmax, mu, rp)
            asym = dv_from_turn_deficit(float(vinf), dreq, dmax)
            assert oberth < asym + 1e-9, (
                f"Oberth not <= asymptote at vinf={vinf:.2f}, extra={extra_deg} deg: "
                f"oberth={oberth:.4f} asym={asym:.4f}"
            )
            # strict: a genuine Oberth credit for nonzero deficit.
            assert oberth < asym, f"no strict credit at vinf={vinf:.2f}, extra={extra_deg}"
            n_checked += 1
    assert n_checked > 0


def test_oberth_advantage_not_universal_high_vinf() -> None:
    # Documented caveat: at high vinf (shallow well relative to vinf) the
    # periapsis maneuver can EXCEED the asymptote rotation. This is the honest
    # boundary of property (b), not a regression.
    vinf = 12.0
    dmax = max_bend(MU_EARTH, RP_EARTH, vinf)
    dreq = dmax + radians(10.0)
    oberth = dv_powered_flyby_periapsis(vinf, dreq, dmax, MU_EARTH, RP_EARTH)
    asym = dv_from_turn_deficit(vinf, dreq, dmax)
    assert oberth > asym, "expected high-vinf Oberth disadvantage as documented"


# --- (d) units / frame consistency + Aldrin point -------------------------


def test_aldrin_point_value() -> None:
    # The recovered Aldrin closure: vinf 6.8615 km/s, required turn 93.01 deg,
    # 200 km Earth flyby. Asymptote = 2.9138 km/s (the regression anchor);
    # Oberth = 1.9336 km/s. Mathematics of the two formulas at one point.
    vinf = 6.861549580708436
    dreq = radians(93.01070843411748)
    dmax = max_bend(MU_EARTH, RP_EARTH, vinf)
    asym = dv_from_turn_deficit(vinf, dreq, dmax)
    oberth = dv_powered_flyby_periapsis(vinf, dreq, dmax, MU_EARTH, RP_EARTH)
    assert asym == pytest.approx(2.9138, abs=1e-3)
    assert oberth == pytest.approx(1.9336, abs=1e-3)
    assert oberth < asym


def test_invalid_inputs_raise() -> None:
    dmax = max_bend(MU_EARTH, RP_EARTH, 6.0)
    with pytest.raises(ValueError):
        dv_powered_flyby_periapsis(-1.0, radians(90.0), dmax, MU_EARTH, RP_EARTH)
    with pytest.raises(ValueError):
        dv_powered_flyby_periapsis(6.0, radians(90.0), dmax, -MU_EARTH, RP_EARTH)
    with pytest.raises(ValueError):
        dv_powered_flyby_periapsis(6.0, radians(90.0), dmax, MU_EARTH, -RP_EARTH)


def test_required_turn_at_or_above_pi_falls_back() -> None:
    # No finite periapsis target exists for dreq >= pi; fall back to the
    # asymptote cost rather than returning a spuriously small number.
    vinf = 3.0
    dmax = max_bend(MU_EARTH, RP_EARTH, vinf)
    val = dv_powered_flyby_periapsis(vinf, pi, dmax, MU_EARTH, RP_EARTH)
    assert val == pytest.approx(dv_from_turn_deficit(vinf, pi, dmax))
