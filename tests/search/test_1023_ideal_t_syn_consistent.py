"""#1023: the ideal Galilean model's synodic period must be the one its smas are built for.

``resonant_conic.ideal_moon_smas`` sets the Europa and Ganymede smas so that every moon advances
2*pi*k + Delta (Delta = 5.2 deg, Hernandez et al. 2017 p.3) in one synodic period, i.e.
T_syn = (2*pi + Delta) / n_Ganymede. Over that period all three moons advance by exactly Delta (a
rigid rotation); over the coded ``ideal_t_syn`` (the ideal Ganymede period) they do not.
"""

from __future__ import annotations

import math

import pytest

from cyclerfinder.search.resonant_conic import (
    IDEAL_DELTA_RAD,
    MU_JUPITER_KM3_S2,
    ideal_moon_smas,
    ideal_t_syn,
)


def _advance_mod_2pi(a_km: float, t_s: float) -> float:
    n = math.sqrt(MU_JUPITER_KM3_S2 / a_km**3)
    return (n * t_s) % (2.0 * math.pi)


def test_coded_ideal_t_syn_is_not_a_rigid_rotation() -> None:
    adv = {m: _advance_mod_2pi(a, ideal_t_syn()) for m, a in ideal_moon_smas().items()}
    assert max(adv.values()) - min(adv.values()) > math.radians(10.0)


@pytest.mark.xfail(strict=True, reason="#1023: ideal_t_syn_consistent added in the next commit")
def test_consistent_t_syn_advances_every_moon_by_delta() -> None:
    import cyclerfinder.search.resonant_conic as rc

    t = rc.ideal_t_syn_consistent()  # type: ignore[attr-defined]
    for a in ideal_moon_smas().values():
        assert _advance_mod_2pi(a, t) == pytest.approx(IDEAL_DELTA_RAD, abs=1e-12)
    assert t / 86400.0 == pytest.approx(7.10536, abs=1e-4)
