"""Brjuno 1978 parts II and III: Tables I and II, the e* membership test (Theorem 2.2),
the asymmetric arcs T_N and the junction closed forms.

Expected values are printed numbers (fixtures transcribed from the digests); the e* test
is checked against Henon 1968's printed tables, which it did not use.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path

import pytest

from cyclerfinder.search import second_species_arcs as m

PI = math.pi
FX = Path(__file__).parent / "fixtures" / "second_species_arcs"


def _read(name: str) -> list[dict[str, str]]:
    with open(FX / name) as f:
        return list(csv.DictReader(f))


# --- Table I: 97 synodic orbits of the families A0 .. C35 -------------------------------------
TABLE1 = _read("brjuno1978b_table1.csv")
TANGENT_ORBITS = {"6", "11", "14", "18", "22", "26", "27", "35", "41", "49", "104", "109"}
MISPRINTS = {"53", "61"}


def _family_indices(fam: str) -> tuple[str, int, int]:
    if fam[0] == "A":
        return "A", int(fam[1:]), 0
    if fam[0] == "B":
        return "B", 0, int(fam[1:])
    return "C", int(fam[1]), int(fam[2])


def _signs(r: dict[str, str]) -> tuple[float, float, int, int, int]:
    at, et = float(r["a_tilde"]), float(r["e_tilde"])
    eps = 1 if at > 0 else -1
    eps1 = 1 if et > 0 else -1
    eps2 = 1 if abs(et) > 1 else -1
    return abs(at), abs(abs(et) - 1.0), eps, eps1, eps2


def test_table1_has_97_rows() -> None:
    assert len(TABLE1) == 97


def test_table1_sign_of_a_tilde_is_the_family_sigma() -> None:
    """Brjuno's eps = sgn a-tilde equals sigma: A -> -, B -> +, C_(p,p+q) -> (-1)^q
    (the parity table p.40 puts q even in omega_1, omega_2 where a-tilde > 0)."""
    for r in TABLE1:
        k, j, kk = _family_indices(r["family"])
        want = m.family_sigma(k, j, kk)
        assert (1 if float(r["a_tilde"]) > 0 else -1) == want, r["orbit"]


GENERIC = [r for r in TABLE1 if r["orbit"] not in TANGENT_ORBITS | MISPRINTS | {"2"}]
PRINTED_WRONG = [r for r in TABLE1 if r["orbit"] in MISPRINTS]
TANGENT = [r for r in TABLE1 if r["orbit"] in TANGENT_ORBITS]


def _nearest_tau(r: dict[str, str]) -> float:
    a, e, eps, eps1, eps2 = _signs(r)
    br = m.tau_branches_from_a_e(a, e, eps, eps1, eps2)
    tp = float(r["tau_pi"])
    return min((abs(t / PI - tp) for _, t in br), default=9.0)


@pytest.mark.parametrize("r", GENERIC, ids=lambda r: f"orbit{r['orbit']}-{r['family']}")
def test_brjuno_table1_tau_from_a_and_e(r: dict[str, str]) -> None:
    """tau/pi from (a-tilde, e-tilde) through eqs. 2.3 and 3.9 (3e-4 on 82 rows)."""
    assert _nearest_tau(r) < 3e-4


@pytest.mark.parametrize("r", TANGENT, ids=lambda r: f"orbit{r['orbit']}-{r['family']}")
def test_brjuno_table1_tangent_rows_to_three_e_minus_3(r: dict[str, str]) -> None:
    """Type II rows (integer tau/pi, |cos eta| = 1): the root is ill conditioned, 3e-3."""
    assert _nearest_tau(r) < 3e-3


def test_table1_circular_row_is_the_a0_circle() -> None:
    """Orbit 2: a-tilde = e-tilde = -1, tau/pi = 0.5: Henon's circular row (a = 1, e = 0)."""
    r = next(r for r in TABLE1 if r["orbit"] == "2")
    arc = m.arc_elements(float(r["tau_pi"]) * PI, 0.5 * PI, -1)
    assert arc.a == pytest.approx(abs(float(r["a_tilde"])), abs=1e-12)
    assert arc.e == pytest.approx(0.0, abs=1e-12)
    assert arc.eps1 == (1 if float(r["e_tilde"]) > 0 else -1)


@pytest.mark.parametrize(
    ("orbit", "check"), [("53", "a"), ("61", "tau")], ids=["orbit53-a", "orbit61-tau"]
)
@pytest.mark.xfail(strict=True, reason="printed value is wrong (digest 2.1)")
def test_brjuno_table1_each_misprint_strict(orbit: str, check: str) -> None:
    r = next(r for r in TABLE1 if r["orbit"] == orbit)
    assert _nearest_tau(r) < 3e-4


def test_brjuno_table1_corrected_orbits_53_and_61() -> None:
    """The recomputed values: orbit 53 a = 1.66394 (printed e-tilde 1.69049 and tau/pi 4.56 agree
    with it); orbit 61 tau/pi = 3.56002 (the leading 3 was dropped in the print)."""
    r53 = next(r for r in TABLE1 if r["orbit"] == "53")
    a, e, eps, eps1, eps2 = _signs(r53)
    br = m.tau_branches_from_a_e(1.66394, e, eps, eps1, eps2)
    assert min(abs(t / PI - float(r53["tau_pi"])) for _, t in br) < 3e-4
    r61 = next(r for r in TABLE1 if r["orbit"] == "61")
    a, e, eps, eps1, eps2 = _signs(r61)
    br = m.tau_branches_from_a_e(a, e, eps, eps1, eps2)
    assert min(abs(t / PI - 3.56) for _, t in br) < 1e-4


# --- e* membership test (Theorem 2.2) --------------------------------------------------------
HENON = _read("henon1968_tables2_9.csv")
HENON_SIGMA = {"A0": -1, "A1": -1, "A2": -1, "B1": 1, "B2": 1, "C12": -1, "C23": -1, "C24": 1}
HENON_FAMILY = {
    "A0": ("A", 0, 0),
    "A1": ("A", 1, 0),
    "A2": ("A", 2, 0),
    "B1": ("B", 0, 1),
    "B2": ("B", 0, 2),
    "C12": ("C", 1, 2),
    "C23": ("C", 2, 3),
    "C24": ("C", 2, 4),
}
HENON_DEFECT_ROW = ("3", "1.80000", "1.24723")


def _henon_determinate() -> list[tuple[dict[str, str], m.SArc]]:
    out = []
    for r in HENON:
        if r["eta_pi"].startswith("0.0000"):
            continue
        eta = float(r["eta_pi"]) * PI
        if abs(math.sin(eta)) < 1e-12:
            continue
        arc = m.arc_elements(float(r["tau_pi"]) * PI, eta, HENON_SIGMA[r["family"]])
        if arc.eps1 == 0 or arc.e < 1e-6 or arc.e > 1 - 1e-6 or abs(arc.a - 1) >= arc.a * arc.e:
            continue
        out.append((r, arc))
    return out


def test_e_star_reproduces_every_determinate_henon_row() -> None:
    """Digest: 276 determinate rows (202 with a > 1, 74 with a < 1) lie on the Theorem 2.2
    characteristic of their family; the only exception is Table 3's misprinted row."""
    rows = _henon_determinate()
    fails = []
    n_gt = n_lt = 0
    for r, arc in rows:
        k, j, kk = HENON_FAMILY[r["family"]]
        if not m.e_star_matches(arc.a, arc.e, arc.eps1, k, j, kk):
            fails.append((r["table"], r["tau_pi"], r["eta_pi"]))
        elif arc.a > 1:
            n_gt += 1
        else:
            n_lt += 1
    assert fails == [HENON_DEFECT_ROW]
    assert (n_gt, n_lt) == (202, 74)


def test_e_star_rejects_a_wrong_family_label() -> None:
    """Negative control: relabelled families (A1 as A3, A0 as A3, B1 as B2, C23 as C24)
    match under a fifth of the rows.  (A1 and A2 share the omega_4 form of Theorem 2.2,
    so that pair is not a discriminating control.)"""
    rows = _henon_determinate()
    for fam, label in (
        ("A1", ("A", 3, 0)),
        ("A0", ("A", 3, 0)),
        ("B1", ("B", 0, 2)),
        ("C23", ("C", 2, 4)),
    ):
        sel = [arc for r, arc in rows if r["family"] == fam]
        wrong = sum(m.e_star_matches(arc.a, arc.e, arc.eps1, *label) for arc in sel)
        assert wrong < 0.2 * len(sel), (fam, wrong, len(sel))


def test_e_star_on_table1_rows_lies_on_the_family_characteristic() -> None:
    """Of the 91 Table I rows with a collision off the tangent points, 90 sit on their
    family's Theorem 2.2 curve; orbit 53 (printed a-tilde wrong) is the one that does not."""
    ok: list[str] = []
    off: list[str] = []
    for r in TABLE1:
        a, e, _eps, eps1, _eps2 = _signs(r)
        if e < 1e-9 or abs(a - 1.0) >= a * e * (1 + 1e-9):
            continue  # tangent rows: the radicand of e* vanishes
        k, j, kk = _family_indices(r["family"])
        (ok if m.e_star_matches(a, e, eps1, k, j, kk) else off).append(r["orbit"])
    assert off == ["53"]
    assert len(ok) == 90


# --- Table II: the curve f (P = 0), five checks plus the e* values ----------------------------
TABLE2 = _read("brjuno1978b_table2.csv")


def test_table2_has_28_rows() -> None:
    assert len(TABLE2) == 28


@pytest.mark.parametrize("r", TABLE2, ids=lambda r: f"a={r['a']}")
def test_brjuno_table2_row(r: dict[str, str]) -> None:
    a, e = float(r["a"]), float(r["e"])
    ni, es = float(r["ninv"]), float(r["e_star"])
    x, phi = float(r["x"]), float(r["phi"])
    # (i) P(a, e) = 0 (eq. 2.16') to the printed digits
    assert abs(m.curve_f_p(a, e)) < 3e-5
    # (ii) N^-1 = a^(3/2)
    assert a**1.5 == pytest.approx(ni, abs=1.5e-5)
    # (iii) x = 1/(N - 1), N = a^(-3/2), from the printed N^-1 (5 digits limits the last rows)
    xx = 1.0 / (1.0 / ni - 1.0)
    assert xx == pytest.approx(x, rel=3 * 5e-6 / (ni * (1.0 - ni)))
    # (iv) phi = y - (x - z0) with y = 2 arccos(e*)/(pi (1 - N^-1))
    y = 2.0 * math.acos(es) / (PI * (1.0 - ni))
    assert y - (xx - m.Z0) == pytest.approx(phi, abs=4e-5 + 6 * 5e-6 / (1.0 - ni))
    # (v) parametric form a = (c^2 + 1)/(c^3 - c + 2), 1 - e^2 = c^2/a, c^2 = a (1 - e^2)
    c = math.sqrt(a * (1.0 - e * e))
    a_par, e_par = m.curve_f_from_c(c)
    assert a_par == pytest.approx(a, abs=2e-5)
    assert e_par == pytest.approx(e, abs=2e-5)
    # (vi) the e* reconstruction (eps1 = -1, a < 1) reproduces the printed e* column
    assert m.e_star_abs(a, e, -1) == pytest.approx(es, abs=3e-5)


def test_table2_asymptotic_constants_and_eq_2_32() -> None:
    """z0 = 0.5469181607, L1 = -0.4428283507, L2 = -0.3525286384 and eq. 2.32 at the last row."""
    z0 = 1.0 / (2.0 * math.sqrt(2.0) - 1.0)
    l1 = -(16.0 / (9.0 * PI)) * (3.0 / 8.0) ** 0.25
    l2 = -(39.0 / (45.0 * PI)) * (8.0 / 3.0) ** 0.25
    assert z0 == pytest.approx(0.5469181607, abs=1e-10)
    assert l1 == pytest.approx(-0.4428283507, abs=1e-10)
    assert l2 == pytest.approx(-0.3525286384, abs=1e-10)
    x = 1650.0017
    phi = 0.5 + z0 + l1 * x**0.25 + l2 * x**-0.25
    assert phi == pytest.approx(-1.83099, abs=5e-4)  # O(x^(-3/4)) remainder, 3e-4 in the digest
    assert max(float(r["phi"]) for r in TABLE2) == pytest.approx(0.18381, abs=1e-12)


def test_curve_f_near_a_equal_one_expansion() -> None:
    """Eq. 2.33: a - 1 = -(1/4) e^4 - (5/16) e^6 at the printed row a = 0.99960, e = 0.19802."""
    e = 0.19802
    assert -(e**4) / 4.0 - 5.0 * e**6 / 16.0 == pytest.approx(0.99960 - 1.0, abs=5e-6)


# --- Type III expansion: printed 1/2 is wrong, 1/12 is right (digest) -------------------------
def test_type_iii_expansion_coefficient_is_one_twelfth_not_one_half() -> None:
    """Near tau = eta = pi/2 on A0, (a - 1)/e^4 -> 1/12; Henon Table 2 rows tau/pi 0.4 and 0.6
    (e = 0.31255 and 0.31069) bracket it, while the printed 1/2 would give a - 1 near 0.005."""
    rows = {r["tau_pi"]: r for r in HENON if r["table"] == "2"}
    for key in ("0.40000", "0.60000"):
        r = rows[key]
        ratio = (float(r["a"]) - 1.0) / float(r["e"]) ** 4
        assert abs(ratio - 1.0 / 12.0) < 0.06  # 0.128 and 0.061, the fifth-order term is +- 0.1 e
        assert ratio < 0.2


@pytest.mark.xfail(strict=True, reason="Brjuno II p.44 prints (a - 1)/e^4 -> 1/2; it is 1/12")
def test_type_iii_expansion_printed_half() -> None:
    rows = {r["tau_pi"]: r for r in HENON if r["table"] == "2"}
    for key in ("0.40000", "0.60000"):
        r = rows[key]
        assert (float(r["a"]) - 1.0) / float(r["e"]) ** 4 == pytest.approx(0.5, rel=0.3)


# --- T_N asymmetric arcs -----------------------------------------------------------------------
@pytest.mark.parametrize(
    ("p", "q", "n_num", "n_den"),
    [
        (1, 2, 3, 1),
        (1, 1, 2, 1),
        (3, 2, 5, 3),
        (2, 1, 3, 2),
        (3, 1, 4, 3),
        (1, 0, 1, 1),
        (2, -1, 1, 2),
    ],
)
def test_brjuno_p_q_to_n_map(p: int, q: int, n_num: int, n_den: int) -> None:
    """Printed (p, q) <-> N = 3, 2, 5/3, 3/2, 4/3, 1, 1/2 (eq. 2.7-2.8, p.14)."""
    assert (p + q) / p == pytest.approx(n_num / n_den)


def test_t_n_exists_below_two_root_two() -> None:
    assert m.t_n_exists(1, 1)  # N = 2
    assert m.t_n_exists(3, 2)  # N = 5/3
    assert not m.t_n_exists(1, 2)  # N = 3: the intersection with Delta_3 is empty
    with pytest.raises(ValueError):
        m.t_n_arc(1, 2, 1, 0.5)


@pytest.mark.parametrize(("p", "q"), [(1, 1), (3, 2), (2, 1), (3, 1), (2, -1)])
@pytest.mark.parametrize("eps1", [1, -1])
def test_t_n_arcs_are_collision_arcs_of_the_resonant_ellipse(p: int, q: int, eps1: int) -> None:
    """Every admissible eta gives an ellipse of mean motion N through P2's circle at r = 1 with
    equal departure and arrival speed (turn 0) and V1^2 + V2^2 = 3 - C (eq. 3.8, 4.C)."""
    n_checked = 0
    for _k, lo, hi in m.t_n_eta_intervals(p, q, 3):
        for frac in (0.1, 0.5, 0.9):
            eta = lo + frac * (hi - lo)
            if abs(math.cos(eta)) < 1e-6:
                continue
            arc = m.t_n_arc(p, q, eps1, eta)
            n_checked += 1
            assert arc.a ** (-1.5) == pytest.approx((p + q) / p, rel=1e-12)
            # eq. 2.2: the collision point is at distance 1 from P1
            assert arc.a * (1.0 - arc.eps2 * arc.e * math.cos(eta)) == pytest.approx(1.0, abs=1e-12)
            assert arc.v1**2 + arc.v2**2 == pytest.approx(3.0 - arc.jacobi, abs=1e-12)
            assert arc.turn_deg == 0.0
            # time between successive collisions at the same point: (p + q) revolutions of P3,
            # p revolutions of P2
            assert arc.duration / (TWO_PI * arc.a**1.5) == pytest.approx(p + q, rel=1e-12)
            assert arc.duration / TWO_PI == pytest.approx(p)
    assert n_checked > 0


TWO_PI = 2.0 * PI


def test_t_n_eta_outside_the_intervals_is_rejected() -> None:
    p, q = 1, 1  # a = 0.62996, |1 - 1/a| = 0.5874, J_0 = [-0.9437, 0.9437]
    with pytest.raises(ValueError):
        m.t_n_arc(p, q, 1, 1.2)


def test_t_n_at_n_equal_one_is_the_circle_with_the_moon() -> None:
    """N = 1 (p, q) = (1, 0): a = 1 and eps2 e = 0 for cos eta != 0; the arc is the circle
    (type III for eps1 = -1: C = -1; type IV* for eps1 = +1: C = 3)."""
    arc_m = m.t_n_arc(1, 0, -1, 0.3)
    arc_p = m.t_n_arc(1, 0, 1, 0.3)
    assert arc_m.jacobi == pytest.approx(m.JUNCTION_CONSTANTS["type_III_C"])
    assert arc_p.jacobi == pytest.approx(m.JUNCTION_CONSTANTS["type_IV_star_C"])
    assert arc_p.speed == pytest.approx(0.0, abs=1e-12)


def test_t_n_joins_the_s_family_at_the_type_ii_point() -> None:
    """At the type II point the S arcs are the same orbit described k times: the T_N arc with
    eta = j pi has a = (p/(p+q))^(2/3), the tangent ellipse of Henon eq. 41 (here (1, 1)
    against Table 7's first row)."""
    arc = m.t_n_arc(1, 1, 1, 0.0)
    assert arc.a == pytest.approx(0.62996, abs=6e-6)
    assert arc.e == pytest.approx(0.58740, abs=6e-6)
    assert arc.jacobi == pytest.approx(2.87208, abs=6e-5)  # Table 7 row 1


# --- junction constants (Theorem 3.1 text) ----------------------------------------------------
def test_junction_values_are_isolated_in_c() -> None:
    """Type II junctions have distinct C (algebraic) apart from the conjugate pair, and
    differ from the type III (-1) and type IV* (3) values (Theorem 3.1 at small p, q)."""
    seen: dict[float, tuple[int, int, int]] = {}
    for p in range(1, 7):
        for q in range(-p + 1, 4 * p):
            if math.gcd(p, abs(q)) != 1 or (p + q) / p >= 2.0 * math.sqrt(2.0) or p + q <= 0:
                continue
            if (p + q) / p == 1.0:
                continue
            for eps1 in (1, -1):
                c = round(m.type_ii_junction(p, q, eps1)["C"], 9)
                assert c not in seen, ((p, q, eps1), seen[c])
                seen[c] = (p, q, eps1)
                assert abs(c - (-1.0)) > 1e-6
                assert abs(c - 3.0) > 1e-6
