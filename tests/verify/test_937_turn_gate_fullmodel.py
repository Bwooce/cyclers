"""#937: the demanded-turn gate on full-model orbits, with published positive controls.

The #899 review reported that Casoliva's published 7-3b/c orbit "demands 32.1 degrees against
19.1 available" and so fails the #888 gate. That number is the angle between two Moon-relative
velocities expressed in the ROTATING axes of two instants a day apart (the axes turn 13.2 degrees
per day), compared with the turn of the osculating hyperbola at the pass radius rather than with
the bend available at the altitude floor. :mod:`cyclerfinder.verify.turn_gate_fullmodel` measures
the turn in inertial axes; these tests hold it to published orbits and to closed forms.

Sources of every state (printed values only):

* Casoliva, Mondelo, Villac, Mease, Gomez, Lizy-Destrez 2010 Table 3 (JGCD 33(5) p.1630), vendored
  in :mod:`cyclerfinder.search.earth_moon_resonant_families`, corrected at the paper's mu by
  :func:`second_species_continuation.table3_reference` (moves the printed crossing by at most
  3e-10, ``tests/search/test_second_species_continuation.py``); Casoliva et al. 2008 Table 2 seed
  32a (mu = 1e-6) from the same module.
* Leiva & Briozzo 2005 (CMDA 91:357, Sect. 4.3 p368) three-body seed, and Leiva & Briozzo 2008
  (CMDA 101:225) Table 1 (p234), Table 5 (p241); transcriptions as in
  ``tests/core/test_leiva_briozzo_2005_more.py`` and
  ``tests/core/test_leiva_briozzo_2008_tables.py``.
* Oshima 2022 (ASR 70:1325) Table 2 row 1 (p1332) and Table 1 constants (p1326), as in
  ``tests/core/test_bcr4bp_oshima_2022.py``.
* Lantoine & Russell 2011 (J. Astronaut. Sci. 58(3)) Table 4 (p353), as in
  ``tests/core/test_cr3bp_lantoine_russell_2011.py``.

The expected outcome of a published ballistic orbit is "not rejected": ``pass`` where the pass is
a hyperbola, ``indeterminate`` where theory says the patched conic does not apply (a slow,
temporarily captured pass), ``no_encounter`` where the orbit never enters the sphere of influence,
and ``fail`` only where the periapsis is inside the body by the paper's own numbers.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from cyclerfinder.core import cr3bp, qbcp
from cyclerfinder.core.bcr4bp import BCR4BPSystem
from cyclerfinder.core.satellites import SATELLITES
from cyclerfinder.search import earth_moon_resonant_families as emrf
from cyclerfinder.search import second_species_continuation as ssc
from cyclerfinder.verify import turn_gate_fullmodel as tgf
from cyclerfinder.verify.turn_gate import (
    TIDAL_BAND_FACTOR,
    Encounter,
    demanded_turn_gate,
    evaluate_encounter,
)

MOON = SATELLITES["Moon"]
EM = cr3bp.cr3bp_system("Earth", "Moon")


def _em_model(
    mu: float, radius_km: float = MOON.radius_eq_km, floor_km: float = MOON.safe_alt_km
) -> tgf.RotatingModel:
    return tgf.cr3bp_model(
        mu, l_km=384400.0, t_s=EM.t_s, body="Moon", radius_km=radius_km, alt_floor_km=floor_km
    )


def _casoliva(designation: str) -> tgf.FullModelTurnReport:
    ref = ssc.table3_reference(designation)
    return tgf.full_model_pass_turns(_em_model(ssc.CASOLIVA_MU_2010), ref.nodes[0], ref.period)


# ---------------------------------------------------------------------------
# Instrument control with a known answer: a pure Moon-centred hyperbola
# ---------------------------------------------------------------------------


def _kepler_only_model(mu: float, radius_km: float) -> tgf.RotatingModel:
    """Rotating frame (rate 1) in which the spacecraft feels only the body (fixed at
    ``c = (1 - mu, 0, 0)``, itself on a circle): its body-relative motion is an exact two-body
    hyperbola, so the inertial turn is known in closed form and the non-body part is zero."""
    c = np.array([1.0 - mu, 0.0, 0.0])
    z = np.array([0.0, 0.0, 1.0])

    def rhs(t: float, s: np.ndarray) -> np.ndarray:
        r, v = s[:3], s[3:]
        rho = r - c
        a = -mu * rho / np.linalg.norm(rho) ** 3 - c
        a = a - 2.0 * np.cross(z, v) - np.cross(z, np.cross(z, r))
        return np.concatenate([v, a])

    return tgf.RotatingModel(
        "kepler", mu, rhs, 384400.0, EM.t_s, "Moon", radius_km, MOON.safe_alt_km
    )


def _hyperbola_start(
    mu: float, rp: float, vinf: float, r0: float
) -> tuple[np.ndarray, float, float]:
    """Rotating-frame state at radius r0 on the incoming branch (t = 0), the time to periapsis,
    and the eccentricity."""
    e = 1.0 + rp * vinf**2 / mu
    p = rp * (1.0 + e)
    f = -math.acos((p / r0 - 1.0) / e)
    rho = r0 * np.array([math.cos(f), math.sin(f), 0.0])
    w = math.sqrt(mu / p) * np.array([-math.sin(f), e + math.cos(f), 0.0])
    c = np.array([1.0 - mu, 0.0, 0.0])
    v_rot = w - np.cross([0.0, 0.0, 1.0], rho)
    big_f = math.acosh((e + math.cos(f)) / (1.0 + e * math.cos(f)))
    a = mu / vinf**2
    t_peri = (e * math.sinh(big_f) - big_f) / math.sqrt(mu / a**3)
    return np.concatenate([c + rho, v_rot]), t_peri, e


def _window_turn_closed_form(e: float, rp: float, r_w: float) -> float:
    """Angle between the velocities at r = r_w before and after periapsis: 2 (f - gamma), with
    cos f = (p / r - 1) / e and tan gamma = e sin f / (1 + e cos f) (flight-path angle)."""
    p = rp * (1.0 + e)
    f = math.acos((p / r_w - 1.0) / e)
    gamma = math.atan2(e * math.sin(f), 1.0 + e * math.cos(f))
    return math.degrees(2.0 * (f - gamma))


def test_pure_hyperbola_turn_matches_the_closed_form() -> None:
    """Known answer: the window turn of an exact two-body hyperbola, viewed from the rotating
    frame for the 1.6 days it spends inside the sphere, equals the closed form, the non-body part
    is zero and the osculating turn is 2 asin(1/e) (the asymptotic limit). The rotating-axes angle
    differs from it by the frame's rotation over the window."""
    mu = 0.0121505
    rp, vinf = 0.01, 0.5
    model = _kepler_only_model(mu, MOON.radius_eq_km)
    s0, t_peri, e = _hyperbola_start(mu, rp, vinf, 1.2 * model.soi)
    rep = tgf.full_model_pass_turns(model, s0, 2.5 * t_peri, periodic=False)
    (p,) = rep.passes
    assert p.complete and p.n_periapses == 1
    assert p.rp_km == pytest.approx(rp * 384400.0, rel=1e-9)
    assert p.e_osc == pytest.approx(e, rel=1e-9)
    assert p.turn_osc_deg == pytest.approx(math.degrees(2.0 * math.asin(1.0 / e)), abs=1e-8)
    assert p.demanded_turn_deg == pytest.approx(
        _window_turn_closed_form(e, rp, model.soi), abs=1e-7
    )
    assert abs(p.rest_part_deg) < 1e-7
    assert p.body_part_deg == pytest.approx(p.demanded_turn_deg, abs=1e-7)
    frame_rotation = math.degrees(p.t_out - p.t_in)
    assert abs(p.rotating_axes_turn_deg - p.demanded_turn_deg) == pytest.approx(
        frame_rotation, abs=1e-6
    )
    assert p.status == "pass"


def test_pure_hyperbola_below_the_surface_fails() -> None:
    """Negative control in the full-model path: the same hyperbola with its periapsis at
    1153 km from the centre (inside the Moon) fails, and so does the gate on its window vectors
    (the required periapsis is inside the body)."""
    mu = 0.0121505
    rp, vinf = 0.003, 0.5
    model = _kepler_only_model(mu, MOON.radius_eq_km)
    s0, t_peri, _ = _hyperbola_start(mu, rp, vinf, 1.2 * model.soi)
    rep = tgf.full_model_pass_turns(model, s0, 2.5 * t_peri, periodic=False)
    (p,) = rep.passes
    assert p.status == "fail"
    assert "inside the body" in p.reason
    assert p.gate is not None and not p.gate.turn_feasible
    assert rep.status == "fail"


# ---------------------------------------------------------------------------
# Synthetic negative controls on the patched-conic gate itself
# ---------------------------------------------------------------------------


def test_turn_beyond_the_floor_bend_still_fails() -> None:
    """At the Moon, 1 km/s, 100 km floor: e = 1 + 1837.4 / 4902.8 = 1.375, so the largest bend is
    2 asin(1/1.375) = 93.3 degrees; a 150 degree demand fails with status "fail" (the tidal band,
    3 x 1.8 degrees, is far from the 57 degree deficit)."""
    v = 1.0
    a = math.radians(150.0)
    enc = Encounter.for_body("Moon", [v, 0.0, 0.0], [v * math.cos(a), v * math.sin(a), 0.0])
    e = evaluate_encounter(enc)
    assert e.available_bend_deg == pytest.approx(
        math.degrees(2.0 * math.asin(1.0 / (1.0 + 1837.4 * v * v / MOON.mu_km3_s2))), abs=1e-9
    )
    assert not e.turn_feasible
    assert e.status == "fail"
    assert demanded_turn_gate([enc]).status == "fail"


def test_reversal_is_a_rejection_not_a_non_encounter() -> None:
    """A 180 degree demand at 1.8 km/s at Oberon (the scale of the withdrawn Uranian rows) fails:
    a near-reversal is a demanded turn like any other, not "no encounter"."""
    v = 1.8
    e = evaluate_encounter(Encounter.for_body("Oberon", [v, 0.0, 0.0], [-v, 0.0, 0.0]))
    assert e.demanded_turn_deg == pytest.approx(180.0)
    assert not e.turn_feasible
    assert e.status == "fail"


def test_slow_pass_is_indeterminate_and_keeps_its_feasibility_flag() -> None:
    """At 0.1 km/s the Moon's tidal turn scale (n r_soi / v)^2 is about 3 radians: the patched
    conic has no meaning there (the slow-pass regime), so status is "indeterminate" while
    ``turn_feasible`` keeps the plain comparison."""
    v = 0.1
    a = math.radians(170.0)
    e = evaluate_encounter(
        Encounter.for_body("Moon", [v, 0.0, 0.0], [v * math.cos(a), v * math.sin(a), 0.0])
    )
    assert math.radians(e.tidal_turn_scale_deg) > 1.0
    assert e.status == "indeterminate"
    assert e.turn_feasible == (e.demanded_turn_deg <= e.available_bend_deg)


def test_marginal_verdict_inside_the_tidal_band_is_indeterminate() -> None:
    """A demand 0.5 degree above the bend at 0.6 km/s (tidal scale about 5 degrees) is not a
    rejection: the part of the turn the patched conic leaves out could flip it."""
    v = 0.6
    bend = math.degrees(2.0 * math.asin(1.0 / (1.0 + 1837.4 * v * v / MOON.mu_km3_s2)))
    a = math.radians(bend + 0.5)
    e = evaluate_encounter(
        Encounter.for_body("Moon", [v, 0.0, 0.0], [v * math.cos(a), v * math.sin(a), 0.0])
    )
    assert not e.turn_feasible
    assert e.status == "indeterminate"


# ---------------------------------------------------------------------------
# Published positive controls
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("designation", ["7-3a", "7-3b", "7-3c"])
def test_casoliva_7_3_lunar_passes_pass(designation: str) -> None:
    """Every lunar pass of the published 7-3a, 7-3b and 7-3c orbits (ballistic by publication)
    passes in inertial axes (ratio below 1, gate status "pass"), and the non-Moon part of each
    turn stays inside the gate's tidal band (TIDAL_BAND_FACTOR tidal turn scales)."""
    rep = _casoliva(designation)
    assert rep.status == "pass"
    assert rep.passes
    for p in rep.passes:
        assert p.status == "pass"
        assert p.gate is not None and p.gate.turn_feasible and p.gate.status == "pass"
        assert p.gate.ratio < 1.0
        assert abs(p.rest_part_deg) < TIDAL_BAND_FACTOR * p.tidal_turn_scale_deg
        assert p.body_part_deg + p.rest_part_deg == pytest.approx(p.demanded_turn_deg, abs=1e-6)


def test_casoliva_7_3b_review_error_reproduced() -> None:
    """Regression for the #899 measurement (which reported 32.1 demanded against 19.1, a
    cross-reference only): at the closest 7-3b pass the rotating-axes angle is the inertial turn
    plus the frame's rotation over the window (the pass turns in the frame's sense), it exceeds
    the osculating hyperbola's turn at the pass radius, so that comparison 'fails' a published
    orbit; the inertial turn fits the floor bend."""
    rep = _casoliva("7-3b")
    p = min(rep.passes, key=lambda q: q.rp_km)
    frame_rotation = math.degrees(p.t_out - p.t_in)
    assert p.rotating_axes_turn_deg == pytest.approx(p.demanded_turn_deg + frame_rotation, abs=1e-6)
    assert p.rotating_axes_turn_deg > p.turn_osc_deg  # the review's false failure
    assert p.gate is not None
    assert p.demanded_turn_deg < p.gate.available_bend_deg
    assert p.status == "pass"


@pytest.mark.parametrize("designation", ["1-2b", "1-2c", "1-2d", "1-2e", "2-1a", "2-1b", "3-2c"])
def test_casoliva_rows_without_a_lunar_encounter(designation: str) -> None:
    """These published orbits never enter the Moon's sphere of influence (66,100 km): the gate
    has nothing to judge, which is not evidence for or against them."""
    rep = _casoliva(designation)
    assert rep.status == "no_encounter"
    assert rep.min_distance_km > 66_100.0


def test_leiva_briozzo_2005_seed_passes() -> None:
    """The printed three-body seed of Leiva & Briozzo 2005 (Sect. 4.3): its one lunar pass is a
    hyperbola above the floor and passes."""
    mu = qbcp.qbcp_default().mu
    period = qbcp.qbcp_default().sun_period_tu
    seed = np.array([1.02379270, 0.0, 0.0, 0.0, -1.91110553, 0.0])
    rep = tgf.full_model_pass_turns(_em_model(mu), seed, period)
    assert rep.status == "pass"
    (p,) = rep.passes
    assert p.e_osc > 1.0


LB2008_MU = 0.0121505482
LB2008_SECTION_X = 0.836915310
LB2008_SUN_PERIOD = 2.0 * math.pi / 0.925195985520347
# Table 1 (p234): number, h, y, ydot, p, q (paper frame; sign-reversed below).
LB2008_TABLE_1 = {
    "053_2": (-1.58753537, 0.0348725952, -0.0149583173, 5, 2),
    "013": (-1.58740571, -0.0399746624, -0.0441622383, 4, 1),
    "032B_1": (-1.58703219, -0.0280775876, 0.00719954961, 5, 1),
}


def _lb2008(name: str) -> tgf.FullModelTurnReport:
    h, y, ydot, p, q = LB2008_TABLE_1[name]
    mu = LB2008_MU
    x, y, ydot = LB2008_SECTION_X, -y, -ydot
    r1 = math.hypot(x + mu, y)
    r2 = math.hypot(x - 1.0 + mu, y)
    xdot = math.sqrt(2.0 * h + x * x + y * y + 2.0 * (1.0 - mu) / r1 + 2.0 * mu / r2 - ydot**2)
    state = np.array([x, y, 0.0, xdot, ydot, 0.0])
    return tgf.full_model_pass_turns(_em_model(mu), state, p / q * LB2008_SUN_PERIOD)


@pytest.mark.parametrize("name", ["053_2", "013"])
def test_leiva_briozzo_2008_slow_passes_are_indeterminate(name: str) -> None:
    """The Table 1 orbits (Jacobi constant 3.17) pass the Moon slowly (under 0.4 km/s at the
    sphere): at periapsis the osculating orbit is an ellipse, so no hyperbola and no patched-conic
    turn exists. The gate must return "indeterminate", never "fail", for these published orbits."""
    rep = _lb2008(name)
    assert rep.status == "indeterminate"
    assert rep.passes
    for p in rep.passes:
        assert p.status == "indeterminate"
        assert p.e_osc < 1.0 or p.n_periapses > 1 or not p.complete


def test_leiva_briozzo_2008_032b_passes_inside_the_moon() -> None:
    """Orbit 032B_1 of Table 1 comes within 433 km of the Moon's centre; the paper's own Table 5
    gives 725 km for its quasi-bicircular continuation 032B_1_t4, also inside the Moon
    (1737.4 km). A point-mass orbit through the body is not flyable: "fail" is correct here."""
    rep = _lb2008("032B_1")
    assert rep.status == "fail"
    assert rep.min_distance_km < MOON.radius_eq_km
    assert any("inside the body" in p.reason for p in rep.passes)


def test_oshima_2022_never_enters_the_sphere() -> None:
    """Oshima 2022 Table 2 row 1 (bicircular model): closest lunar distance about 84,000 km."""
    omega = 0.925195985
    system = BCR4BPSystem(
        mu=0.0121506683, mu_sun=328900.541, a_sun_nondim=388.811143, omega_sun_nondim=omega
    )
    model = tgf.bcr4bp_model(
        system,
        l_km=384400.0,
        t_s=EM.t_s,
        body="Moon",
        radius_km=MOON.radius_eq_km,
        alt_floor_km=MOON.safe_alt_km,
    )
    state = np.array([-1.107328855, 0.0, 0.0, 0.0, 2.056256309, -0.156230339])
    rep = tgf.full_model_pass_turns(model, state, 2.0 * math.pi / omega)
    assert rep.status == "no_encounter"
    assert rep.min_distance_km > 66_100.0


def test_lantoine_russell_resonant_orbit_never_enters_the_sphere() -> None:
    """Lantoine & Russell 2011 Table 4, the Ganymede 3:4 resonant orbit: closest approach about
    38,500 km, outside Ganymede's sphere of influence (24,300 km)."""
    mu = 7.8037e-5
    sysj = cr3bp.cr3bp_system("Jupiter", "Ganymede")
    gan = SATELLITES["Ganymede"]
    model = tgf.cr3bp_model(
        mu,
        l_km=sysj.l_km,
        t_s=sysj.t_s,
        body="Ganymede",
        radius_km=gan.radius_eq_km,
        alt_floor_km=gan.safe_alt_km,
    )
    state = np.array([0.9639250025, 0.0, 0.0, 0.0, -0.037537693295765, 0.0])
    rep = tgf.full_model_pass_turns(model, state, 19.1527202833)
    assert rep.status == "no_encounter"
    assert rep.min_distance_km > model.soi * model.l_km


def test_deep_pass_matches_ks_and_the_hyperbola() -> None:
    """Casoliva 2008 seed 32a at mu = 1e-6 (a point-mass Moon: radius 1 m, no floor): its lunar pass
    at 9.5e-5 lunar distances is propagated again with the KS-regularised integrator, which agrees
    at the window exit to 1e-10; the window turn is the osculating hyperbola's turn to 0.01 degree
    and the non-Moon part is under 1e-3 degree."""
    seed = ssc.table2_seed("32a")
    rep = tgf.full_model_pass_turns(
        _em_model(ssc.SEED_MU, radius_km=1e-3, floor_km=0.0), seed.state_project(), seed.period
    )
    (p,) = rep.passes
    assert p.rp_km / 384400.0 < tgf.KS_CHECK_RP
    assert p.ks_exit_mismatch < 1e-10
    assert p.demanded_turn_deg == pytest.approx(p.turn_osc_deg, abs=0.01)
    assert abs(p.rest_part_deg) < 1e-3
    assert p.status == "pass"


def test_table3_designations_are_the_vendored_ones() -> None:
    """The controls above cover every catalogued Table 3 row plus the printed 1-2b."""
    covered = {"7-3a", "7-3b", "7-3c", "1-2b", "1-2c", "1-2d", "1-2e", "2-1a", "2-1b", "3-2c"}
    assert set(emrf.TABLE3_VALID_DESIGNATIONS) <= covered
