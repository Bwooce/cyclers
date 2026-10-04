"""Bruno 1981 Tables I, III, IV; Hitzl & Henon 1977b Table 1; Perko 1981 example; the
collision velocity, demanded turn and indeterminacy flags (#906).

Expected values are printed numbers (fixtures transcribed from the digests), never
values computed by this module.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path

import pytest
from scipy.optimize import minimize_scalar

from cyclerfinder.search import second_species_arcs as m

PI = math.pi
FX = Path(__file__).parent / "fixtures" / "second_species_arcs"


def _read(name: str) -> list[dict[str, str]]:
    with open(FX / name) as f:
        return list(csv.DictReader(f))


# --- Bruno 1981 Table I: |W| and V at e = 1 ---------------------------------------------------
TABLE1 = _read("bruno1981_table1.csv")


def test_table1_has_all_printed_rows() -> None:
    # the digest says 23 rows; the table as transcribed has 25 (a = 0.5 and a = infinity included)
    assert len(TABLE1) == 25


@pytest.mark.parametrize("r", TABLE1, ids=lambda r: f"a={r['a']}")
def test_bruno_table1_w_and_v_at_e_equal_one(r: dict[str, str]) -> None:
    a = float("inf") if r["a"] == "infinity" else float(r["a"])
    if math.isinf(a):
        w, v = (math.sqrt(3.0) / math.sqrt(2.0) - 1.0) / 3.0, math.sqrt(3.0)
    else:
        w, v = m.w_e1(a)
    if r["abs_W"] == "infinity":
        assert math.isinf(w)
    else:
        assert w == pytest.approx(float(r["abs_W"]), abs=1e-6)
    if r["V"]:
        assert v == pytest.approx(float(r["V"]), abs=1e-6)


def test_table1_limit_gamma_and_threshold() -> None:
    """gamma = (sqrt 3/sqrt 2 - 1)/3 = 0.0749149 and |W| = 1 near a = 0.57045 (eq. 16)."""
    assert (math.sqrt(3.0) / math.sqrt(2.0) - 1.0) / 3.0 == pytest.approx(0.0749149, abs=1e-7)
    assert m.w_e1(0.57045)[0] == pytest.approx(1.0, abs=2e-3)


# --- Bruno 1981 Table III: the acceptable e = 1 arcs ------------------------------------------
TABLE3 = _read("bruno1981_table3.csv")


def _e1_below(a_max: float) -> list[m.SArc]:
    out: list[m.SArc] = []
    for k in range(1, 5):
        for eps2 in (1, -1):
            out.extend(a for a in m.e1_arcs(k, eps2) if a.a < a_max)
    return out


def test_table3_has_eleven_arcs_below_a_0_725() -> None:
    """Eleven arcs with a < 0.725 in tau/pi <= 4, eta/pi <= 6.6: the table's eleven rows
    (the e = 1 solve also finds a = 0.74280 which Bruno excludes by the 0.725 bound, C34)."""
    arcs = [a for a in _e1_below(0.75) if a.a < 0.725]
    assert len(arcs) == len(TABLE3) == 11


@pytest.mark.parametrize("r", TABLE3, ids=lambda r: f"{r['family']}-{r['a']}")
def test_bruno_table3_a_and_sign_of_w(r: dict[str, str]) -> None:
    printed = float(r["a"])
    cands = [a for a in _e1_below(0.725) if abs(a.a - printed) < 2e-5]
    assert len(cands) == 1
    arc = cands[0]
    if (r["family"], r["a"]) == ("C25", "0.57888"):
        pytest.skip(
            "printed a differs from the documented 0.57889; see the strict xfail test below"
        )
    assert arc.a == pytest.approx(printed, abs=5.2e-6)
    cd = m.collision_data(arc)
    assert (cd.w_eq12 > 0) == (r["sgn_w"] == "+")


@pytest.mark.xfail(strict=True, reason="Table III C25: printed 0.57888, documented as 0.57889")
def test_bruno_table3_c25_printed_0_57888() -> None:
    (arc,) = [a for a in _e1_below(0.725) if abs(a.a - 0.57888) < 2e-5]
    assert arc.a == pytest.approx(0.57888, abs=5.2e-6)


def test_bruno_table3_c25_documented_value() -> None:
    (arc,) = [a for a in _e1_below(0.725) if abs(a.a - 0.57889) < 2e-5]
    assert arc.a == pytest.approx(0.57889, abs=5.2e-6)


def test_table3_exterior_orbits_are_the_three_positive_w() -> None:
    """Text p.266: the orbits of three families (C24, C25, C36) are exterior (W > 0)."""
    pos = sorted(r["family"] for r in TABLE3 if r["sgn_w"] == "+")
    assert pos == ["C24", "C25", "C36"]


# --- Bruno 1981 Table IV: A0 segment, sidereal-form W -----------------------------------------
TABLE4 = _read("bruno1981_table4.csv")


def _table4_arc(r: dict[str, str]) -> m.SArc:
    return m.arc_elements(float(r["tau_pi"]) * PI, float(r["eta_pi"]) * PI, -1)


TABLE4_FINITE = [r for r in TABLE4 if r["a"] != "infinity" and r["eta_pi"] != "0.43"]


@pytest.mark.parametrize("r", TABLE4_FINITE, ids=lambda r: f"tau{r['tau_pi']}")
def test_bruno_table4_rows(r: dict[str, str]) -> None:
    """a, e, rho = a (1 - e), C and the Table IV W column (the sidereal-speed form)."""
    arc = _table4_arc(r)
    cd = m.collision_data(arc)
    assert arc.a == pytest.approx(float(r["a"]), abs=5e-4 + 5e-5 * float(r["a"]))
    assert arc.e == pytest.approx(float(r["e"]), abs=3e-5)
    assert arc.a * (1 - arc.e) == pytest.approx(float(r["rho"]), abs=2e-5)
    assert arc.jacobi == pytest.approx(float(r["C"]), abs=6e-5)
    assert cd.w_sidereal == pytest.approx(float(r["W"]), abs=1.5e-5)


def test_bruno_table4_first_row_is_the_parabolic_orbit() -> None:
    r = TABLE4[0]
    par = m.parabolic_arc()
    assert par.x0 == pytest.approx(-float(r["rho"]), abs=6e-6)
    assert par.jacobi == pytest.approx(float(r["C"]), abs=6e-6)
    w = (par.speed / math.sqrt(2.0) - 1.0) / par.speed**2  # a = infinity: sidereal speed sqrt 2
    assert w == pytest.approx(float(r["W"]), abs=1.5e-5)


@pytest.mark.xfail(
    strict=True, reason="Table IV row tau/pi 0.2358 prints eta/pi 0.43; it closes at 0.43381"
)
def test_bruno_table4_printed_row_eta_0_43_does_not_close() -> None:
    r = next(r for r in TABLE4 if r["eta_pi"] == "0.43")
    arc = _table4_arc(r)
    assert arc.a == pytest.approx(float(r["a"]), abs=5e-4)


def test_bruno_table4_row_closes_at_eta_0_43381() -> None:
    r = next(r for r in TABLE4 if r["eta_pi"] == "0.43")
    arc = m.arc_elements(float(r["tau_pi"]) * PI, 0.43381 * PI, -1)
    assert arc.a == pytest.approx(float(r["a"]), abs=1e-4)
    assert arc.e == pytest.approx(float(r["e"]), abs=5e-5)
    assert arc.jacobi == pytest.approx(float(r["C"]), abs=3e-5)
    assert m.collision_data(arc).w_sidereal == pytest.approx(float(r["W"]), abs=3e-5)


def test_w_forms_agree_at_e_equal_one_and_differ_below() -> None:
    """Eq. 12 and the sidereal form coincide only at e = 1 (digest, answer 2)."""
    (arc,) = [a for a in m.e1_arcs(1, 1) if abs(a.eta / PI - 1.36836) < 1e-3]
    cd = m.collision_data(arc)
    assert abs(cd.w_eq12) == pytest.approx(cd.w_sidereal, rel=1e-9)
    r = next(r for r in TABLE4 if r["tau_pi"] == "0.3")
    cd2 = m.collision_data(_table4_arc(r))
    assert abs(cd2.w_eq12 - cd2.w_sidereal) > 0.1  # 0.555 against 0.2407


@pytest.mark.xfail(strict=True, reason="Table IV W column is not Bruno eq. 12 (digest 2.4)")
def test_bruno_table4_w_column_is_not_eq12() -> None:
    r = next(r for r in TABLE4 if r["tau_pi"] == "0.3")
    assert m.collision_data(_table4_arc(r)).w_eq12 == pytest.approx(float(r["W"]), abs=1.5e-5)


def test_bruno_max_c_on_a0() -> None:
    """Text p.267: the maximum Jacobi constant on the A0 segment is C = -0.39913 (orbit
    a = 1.41019, e = 0.88445) and the critical function S vanishes there."""

    def c_of_tau(tp: float) -> float:
        (eta,) = m.find_etas(tp * PI, -1, eta_max=0.6 * PI)
        return -m.arc_elements(tp * PI, eta, -1).jacobi

    res = minimize_scalar(c_of_tau, bounds=(0.2, 0.23), method="bounded", options={"xatol": 1e-10})
    tp = float(res.x)
    (eta,) = m.find_etas(tp * PI, -1, eta_max=0.6 * PI)
    arc = m.arc_elements(tp * PI, eta, -1)
    assert arc.jacobi == pytest.approx(-0.39913, abs=1e-5)
    assert arc.a == pytest.approx(1.41019, abs=2e-4)
    assert arc.e == pytest.approx(0.88445, abs=2e-5)
    assert abs(m.critical_function_s(tp * PI, eta, -1)) < 1e-4


# --- Hitzl & Henon 1977b Table 1: seven critical orbits ---------------------------------------
HH = _read("hitzl_henon1977b_table1.csv")


def _sg(s: str) -> int:
    return 1 if s == "+" else -1


@pytest.mark.parametrize("r", HH, ids=lambda r: r["name"])
def test_hitzl_henon_table1_critical_orbits(r: dict[str, str]) -> None:
    """a, e, x0, x1, C, T from the 5-decimal (tau/pi, eta/pi) and signs, tolerance about 1e-3;
    S = 0 at each critical orbit (S has size O(1) away from them)."""
    tau, eta = float(r["tau_pi"]) * PI, float(r["eta_pi"]) * PI
    s0, s1, s2 = _sg(r["s0"]), _sg(r["s1"]), _sg(r["s2"])
    arc = m.arc_elements(tau, eta, s0 * s2, eps1=s1)
    assert (arc.eps, arc.eps1, arc.eps2) == (s0, s1, s2)
    for col, val in (
        ("a", arc.a),
        ("e", arc.e),
        ("x0", arc.x0),
        ("x1", arc.x1),
        ("C", arc.jacobi),
    ):
        assert val == pytest.approx(float(r[col]), abs=1e-3), col
    assert 2.0 * PI * float(r["tau_pi"]) == pytest.approx(float(r["T"]), abs=2e-4)
    assert abs(m.critical_function_s(tau, eta, s0 * s2)) < 1e-3


def test_s_is_not_small_at_a_non_critical_orbit() -> None:
    """Negative control: S at the A0 row tau/pi = 0.3 (not critical) is O(1)."""
    (eta,) = m.find_etas(0.3 * PI, -1, eta_max=0.6 * PI)
    assert abs(m.critical_function_s(0.3 * PI, eta, -1)) > 0.5


def test_hitzl_a0_minus_1_demanded_turn() -> None:
    """Digest 1977b 9.4 item 3: demanded turn 65.094 degrees and V = 1.8437 at A0(-1)."""
    r = HH[0]
    arc = m.arc_elements(float(r["tau_pi"]) * PI, float(r["eta_pi"]) * PI, -1, eps1=-1)
    cd = m.collision_data(arc)
    assert cd.turn_deg == pytest.approx(65.094, abs=0.02)
    assert cd.speed == pytest.approx(1.8437, abs=2e-4)
    assert not cd.turn_indeterminate


# --- Perko 1981 (2, 1) resonance bifurcation, C = -0.406767 -----------------------------------
@pytest.mark.parametrize(
    ("tau_pi", "eta_pi", "sigma"),
    [(0.203581, 0.366926, -1), (1.796418, 0.633073, 1), (2.203581, 1.366926, 1)],
    ids=["A0-arc", "B1-arc", "E1-ellipse"],
)
def test_perko_bifurcation_orbit(tau_pi: float, eta_pi: float, sigma: int) -> None:
    arc = m.arc_elements(tau_pi * PI, eta_pi * PI, sigma)
    assert arc.jacobi == pytest.approx(-0.406767, abs=3e-6)
    assert arc.speed == pytest.approx(1.845743, abs=3e-6)
    # the retrograde ellipse of mean motion 1/2: a^(-3/2) = 1/2
    assert arc.a ** (-1.5) == pytest.approx(0.5, abs=1e-5)
    assert arc.eps1 == -1
    # the two durations add to 2 pi (1.999999 pi printed)
    assert pytest.approx(2.0, abs=2e-6) == 0.203581 + 1.796418


# --- collision velocity, demanded turn, indeterminate flags (#906) ----------------------------
def test_v1_squared_plus_v2_squared_is_three_minus_c_on_every_table4_row() -> None:
    """Brjuno 1978 4.C: V1^2 + V2^2 = 3 - C (1e-4 on the printed rows)."""
    for r in TABLE4_FINITE:
        arc = _table4_arc(r)
        assert arc.v1**2 + arc.v2**2 == pytest.approx(3.0 - float(r["C"]), abs=1e-4)


def test_turn_formula_sin_half_delta_is_v1_over_v() -> None:
    r = TABLE4_FINITE[3]
    arc = _table4_arc(r)
    cd = m.collision_data(arc)
    delta = math.radians(cd.turn_deg)
    assert math.sin(delta / 2.0) == pytest.approx(abs(arc.v1) / cd.speed, rel=1e-12)
    assert math.cos(delta) == pytest.approx((arc.v2**2 - arc.v1**2) / cd.speed**2, abs=1e-12)


def test_tangent_ellipse_has_zero_turn_flag() -> None:
    """Type II (eqs. 40-41): V1 = 0, so the demanded turn is exactly 0 and flagged."""
    arc = m.tangent_arc(2, 1, eps1=1)
    cd = m.collision_data(arc)
    assert "tangent_resonance" in cd.indeterminate
    assert cd.turn_deg == pytest.approx(0.0, abs=1e-6)
    assert cd.jacobi == pytest.approx(2.97093, abs=6e-5)  # Table 2 row 2.0, eps1 = +
    assert m.collision_data(m.tangent_arc(2, 1, eps1=-1)).jacobi == pytest.approx(
        -1.71101, abs=6e-5
    )


def test_circular_orbit_flag_row() -> None:
    """Table 2 row tau/pi 0.5: a = 1, e = 0, C = -1, V = 2; zero turn, flagged."""
    arc = m.arc_elements(0.5 * PI, 0.5 * PI, -1)
    cd = m.collision_data(arc)
    assert "circular" in cd.indeterminate
    assert cd.turn_deg == pytest.approx(0.0, abs=1e-6)
    assert cd.jacobi == pytest.approx(-1.0, abs=1e-12)
    assert cd.speed == pytest.approx(2.0, abs=1e-12)
    assert cd.w_sidereal == pytest.approx(0.25, abs=1e-12)  # printed Table IV last row
    assert math.isinf(cd.w_eq12)  # eq. 12 has V1 = 0 here, "an infinite W" (digest 2.4)


def test_radial_arc_flag() -> None:
    (arc,) = [a for a in m.e1_arcs(1, 1) if abs(a.eta / PI - 1.36836) < 1e-3]
    assert "radial" in m.collision_data(arc).indeterminate


def test_half_turn_and_zero_speed_flags_on_closed_forms() -> None:
    # V2 = 0: eps1 sqrt(a (1 - e^2)) = 1, e.g. e = 0.6, a = 1/(1 - 0.36)
    e = 0.6
    a = 1.0 / (1.0 - e * e)
    arc = m.SArc(1.0, 2.0, 1, 1, 1, 1, a, e, 0.0, 0.0)
    cd = m.collision_data(arc)
    assert "half_turn" in cd.indeterminate
    # type IV*: a = 1, e = 0, direct: C = 3, V = 0
    iv = m.SArc(1.0, 2.0, 1, 1, 1, 1, 1.0, 0.0, 0.0, 0.0)
    cd4 = m.collision_data(iv)
    assert cd4.jacobi == pytest.approx(m.JUNCTION_CONSTANTS["type_IV_star_C"])
    assert "zero_relative_speed" in cd4.indeterminate
    # type III: a = 1, e = 0, retrograde: C = -1, V = 2, turn 0
    iii = m.SArc(1.0, 2.0, 1, 1, -1, 1, 1.0, 0.0, 0.0, 0.0)
    cd3 = m.collision_data(iii)
    assert cd3.jacobi == pytest.approx(m.JUNCTION_CONSTANTS["type_III_C"])
    assert cd3.speed == pytest.approx(2.0)


def test_type_ii_junction_closed_form_matches_henon_rows() -> None:
    """C = 1/a + 2 eps1 sqrt(2 - 1/a) at the tangent points of Henon Tables 2 and 3 (tau/pi=2)."""
    # (p, q) = (1, 1): N = 2, a = 2^(-2/3)?  a = (p/(p+q))^(2/3) = 0.62996 is Table 7 row 1
    j = m.type_ii_junction(1, 1, 1)
    assert j["a"] == pytest.approx(0.62996, abs=6e-6)
    assert j["e"] == pytest.approx(0.58740, abs=6e-6)
    assert j["C"] == pytest.approx(2.87208, abs=6e-5)  # Table 7 first row
    assert j["V"] == pytest.approx(0.35766, abs=6e-5)
    assert m.type_ii_junction(1, 1, -1)["C"] == pytest.approx(0.30272, abs=6e-5)  # Table 7 row 7
    assert j["turn_deg"] == 0.0


def test_elliptic_collision_data_not_implemented() -> None:
    arc = m.arc_elements(1.0, 2.0, -1, eps1=1, e_p=0.5)
    with pytest.raises(NotImplementedError):
        _ = arc.jacobi
    with pytest.raises(NotImplementedError):
        m.collision_data(arc)
