"""Henon 2001 chapter 18 (R-arcs, R-orbits) with Devaney's counts, and the Gomez & Olle 1986
elliptic extension (timing equation, parabolic arc, eq. 21 existence rule).

Expected values are printed numbers (Henon 2001 Tables 18.2 and 18.3, the p.191 orbit and
eq. 18.49; Gomez & Olle eqs. 19 and 21 and the values quoted in the text).
"""

from __future__ import annotations

import csv
import itertools
import math
from pathlib import Path

import numpy as np
import pytest
from scipy.optimize import fsolve

from cyclerfinder.search import second_species_arcs as m

PI = math.pi
FX = Path(__file__).parent / "fixtures" / "second_species_arcs"


def _read(name: str) -> list[dict[str, str]]:
    with open(FX / name) as f:
        return list(csv.DictReader(f))


# --- R-arcs: Table 18.2 ------------------------------------------------------------------------
T182 = _read("henon2001_table18_2.csv")


def test_table_18_2_has_31_values() -> None:
    assert len(T182) == 31
    # 2^(n-2) per n among sequences beginning R+: 1 + 2 + 4 + 8 + 16 = 31 for n = 2 .. 6
    assert sum(2 ** (int(r["n"]) - 2) for r in T182[:0]) == 0
    per_n = {n: sum(1 for r in T182 if int(r["n"]) == n) for n in range(2, 7)}
    assert per_n == {2: 1, 3: 2, 4: 4, 5: 8, 6: 16}


@pytest.mark.parametrize("r", T182, ids=lambda r: f"n{r['n']}-{r['signs_y1_to_y_nm1']}")
def test_henon_table_18_2_r_arc_y1(r: dict[str, str]) -> None:
    """y_1 of each R-arc to the printed 8 digits (sign code = sign of y_1 .. y_(n-1))."""
    n = int(r["n"])
    signs = [1 if c == "+" else -1 for c in r["signs_y1_to_y_nm1"]]
    assert len(signs) == n - 1
    y = m.r_arc(signs)
    assert y[0] == pytest.approx(float(r["y1"]), abs=6e-9)
    assert m.r_arc_residual(y) < 1e-10


def test_r_arc_closed_forms_n2_n3() -> None:
    """18.11-18.13: n = 2: y_1 = 1/sqrt 2; n = 3: y_1 = y_2 = 1 and y_1 = -y_2 = 1/sqrt 3."""
    assert m.r_arc([1])[0] == pytest.approx(1.0 / math.sqrt(2.0), abs=1e-12)
    y = m.r_arc([1, 1])
    assert y == pytest.approx([1.0, 1.0], abs=1e-12)
    y = m.r_arc([1, -1])
    assert y == pytest.approx([1.0 / math.sqrt(3.0), -1.0 / math.sqrt(3.0)], abs=1e-12)


def test_r_arc_sequences_beginning_r_minus_are_sign_reversals() -> None:
    """Property E' (18.3): reversing every sign of the code reverses the arc."""
    for code in itertools.product((1, -1), repeat=4):
        assert m.r_arc([-c for c in code]) == pytest.approx(-m.r_arc(code), abs=1e-11)


@pytest.mark.parametrize("n", range(2, 11))
def test_r_arc_count_is_two_to_the_n_minus_one(n: int) -> None:
    """Devaney / Henon: exactly 2^(n-1) R-arcs, one per sign code (hard assertion inside)."""
    arcs = m.r_arcs(n)
    assert len(arcs) == 2 ** (n - 1)
    assert max(m.r_arc_residual(y) for y in arcs.values()) < 1e-9
    assert all(np.all(np.sign(y) == np.array(code)) for code, y in arcs.items())


# --- R-orbits: Table 18.3 ----------------------------------------------------------------------
T183 = _read("henon2001_table18_3.csv")


def test_table_18_3_has_21_rows() -> None:
    assert len(T183) == 21


@pytest.mark.parametrize("r", T183, ids=lambda r: f"n{r['n']}-{r['y0']}-{r['y1']}")
def test_henon_table_18_3_r_orbit(r: dict[str, str]) -> None:
    n = int(r["n"])
    printed = [float(r[f"y{i}"]) for i in range(n)]
    y = m.r_orbit([1 if v > 0 else -1 for v in printed])
    assert y == pytest.approx(printed, abs=1e-6)
    assert m.r_orbit_residual(y) < 1e-9


def test_r_orbit_n7_unsymmetric_orbit() -> None:
    """p.191: the first orbit with no symmetry, +++-+--: 9 printed digits for y_0, y_1 and 6 for
    the others."""
    y = m.r_orbit([1, 1, 1, -1, 1, -1, -1])
    # the book's 9th digit is off by 2.4e-9 for y_0 (computed 0.8801420916)
    assert y[0] == pytest.approx(0.880142094, abs=5e-9)
    assert y[1] == pytest.approx(1.302150709, abs=5e-9)
    assert y[2:] == pytest.approx([0.956199, -0.435560, 0.468576, -0.761411, -0.678047], abs=1e-6)


def test_r_orbit_n5_symmetric_polynomial_roots() -> None:
    """Eq. 18.49: 5 y^6 - 20 y^4 + 17 y^2 - 4 = 0; its positive roots are the y_0 = y_1 of the
    three symmetric n = 5 orbits of Table 18.3 (0.652966, 0.799673, 1.712938)."""
    roots = np.roots([5.0, 0.0, -20.0, 0.0, 17.0, 0.0, -4.0])
    pos = sorted(float(r.real) for r in roots if abs(r.imag) < 1e-9 and r.real > 0)
    assert pos == pytest.approx([0.652966, 0.799673, 1.712938], abs=1e-6)
    rows = [r for r in T183 if r["n"] == "5" and float(r["y0"]) > 0]
    assert sorted(float(r["y0"]) for r in rows) == pytest.approx(pos, abs=1e-6)
    solved = sorted(
        float(m.r_orbit(code)[0])
        for code in ([1, 1, 1, -1, 1], [1, 1, -1, 1, -1], [1, 1, -1, -1, -1])
    )
    assert solved == pytest.approx(pos, abs=1e-9)


@pytest.mark.parametrize("n", range(2, 11))
def test_r_orbit_count_is_two_to_the_n_minus_two(n: int) -> None:
    """Devaney Corollary B / Henon 18.2: exactly 2^n - 2 R-orbits, one per sign code except
    all + and all -; each satisfies sum 1/y_i = 0 (18.84) and the recurrence."""
    orbits = m.r_orbits(n)
    assert len(orbits) == 2**n - 2
    for code, y in orbits.items():
        assert np.all(np.sign(y) == np.array(code))
        assert abs(float(np.sum(1.0 / y))) < 1e-8
    assert max(m.r_orbit_residual(y) for y in orbits.values()) < 1e-9


def test_r_orbit_n6_count_decomposition_62() -> None:
    """The book's count for n = 6: 9 listed families x 6 shifts = 54, plus 2 sub-period-2 and 6
    sub-period-3 solutions, total 62 = 2^6 - 2."""
    orbits = m.r_orbits(6)
    ys = list(orbits.values())
    sub2 = [y for y in ys if np.allclose(y, np.roll(y, 2), atol=1e-9)]
    sub3 = [y for y in ys if np.allclose(y, np.roll(y, 3), atol=1e-9)]
    assert len(sub2) == 2
    assert len(sub3) == 6
    full = [
        y
        for y in ys
        if not (
            np.allclose(y, np.roll(y, 2), atol=1e-9) or np.allclose(y, np.roll(y, 3), atol=1e-9)
        )
    ]
    assert len(full) == 54
    classes = {tuple(sorted(tuple(np.sign(np.roll(y, k))) for k in range(6))) for y in full}
    assert len(classes) == 9  # the nine printed families (Table 18.3 lists nine for n = 6)


def test_r_orbit_n5_has_thirty_solutions_six_families() -> None:
    orbits = m.r_orbits(5)
    assert len(orbits) == 30
    classes = {
        tuple(sorted(tuple(np.sign(np.roll(y, k))) for k in range(5))) for y in orbits.values()
    }
    assert len(classes) == 6


def test_r_orbit_rejects_all_equal_sign_codes() -> None:
    with pytest.raises(ValueError):
        m.r_orbit([1, 1, 1])
    with pytest.raises(ValueError):
        m.r_orbit([-1, -1])


def test_r_orbit_extended_precision_polish_agrees_with_double() -> None:
    """n = 13: Newton in 50 digits changes the double-precision orbit by under 1e-10 (so double
    precision is adequate here), and the 50-digit residual vanishes."""
    code = [1, 1, 1, 1, 1, 1, -1, 1, 1, -1, -1, 1, 1]
    y = m.r_orbit(code)
    polished = m.refine_r_orbit_mp(list(y), dps=50)
    diff = max(abs(float(p) - float(v)) for p, v in zip(polished, y, strict=True))
    assert diff < 1e-10
    import mpmath as mp

    mp.mp.dps = 50
    yy = [mp.mpf(p) for p in polished]
    n = len(yy)
    res = max(abs(yy[i - 1] - 2 * yy[i] + yy[(i + 1) % n] + 1 / yy[i]) for i in range(n))
    assert res < mp.mpf(10) ** -40


@pytest.mark.slow
def test_r_region_counts_to_n_equal_14() -> None:
    """Full sweep beyond the default range (several minutes at most)."""
    for n in (11, 12, 13, 14):
        arcs = m.r_arcs(n)
        assert len(arcs) == 2 ** (n - 1)
        orbits = m.r_orbits(n)
        assert len(orbits) == 2**n - 2
        assert max(m.r_orbit_residual(y) for y in orbits.values()) < 1e-8


# --- Gomez & Olle 1986: elliptic extension -----------------------------------------------------
def _parabola_first_principles(e_p: float, eps_p: int, eps: int) -> float:
    """tau/pi of the parabolic arc from the three collision equations directly: P2 on its
    Kepler ellipse (position and time from the eccentric anomaly), P3 on x = eps (p/2)(1 - s^2),
    y = eps eps1 p s, t = p^(3/2)(s/2 + s^3/6), eps1 = -1."""

    def p2(f: float) -> tuple[float, float, float]:
        nu = f if eps_p == 1 else f + PI  # true anomaly from pericentre
        r = (1.0 - e_p**2) / (1.0 + e_p * math.cos(nu))
        ecc = math.atan2(math.sqrt(1.0 - e_p**2) * math.sin(nu), e_p + math.cos(nu))
        if eps_p == -1:
            ecc %= 2.0 * PI  # unwrap: nu = f + pi lies in (pi, 2 pi)
        t = ecc - e_p * math.sin(ecc)
        if eps_p == -1:
            t -= PI  # time measured from apocentre
        return eps_p * r * math.cos(f), eps_p * r * math.sin(f), t

    def eqs(v: list[float]) -> list[float]:
        tau, pp, sig = v
        x, y, t = p2(tau)
        return [
            eps * (pp / 2.0) * (1.0 - sig**2) - x,
            eps * -1.0 * pp * sig - y,
            pp**1.5 * (sig / 2.0 + sig**3 / 6.0) - t,
        ]

    best = None
    for tau0 in np.linspace(0.05, 3.0, 25):
        for sg0 in (1.0, 3.0):
            sol, _info, ier, _msg = fsolve(eqs, [tau0, 0.2, sg0], full_output=True, xtol=1e-13)
            good = ier == 1 and 0.0 < sol[0] < PI and sol[1] > 0 and sol[2] > 0
            if good and max(abs(v) for v in eqs(list(sol))) < 1e-10:
                best = float(sol[0]) / PI
            break
        if best is not None:
            break
    assert best is not None
    return best


def test_parabolic_arc_circular_limit_value() -> None:
    assert m.parabolic_tau(0.0, 1) / PI == pytest.approx(0.16393, abs=6e-6)


def test_parabolic_arc_eps_p_minus_one_printed_value() -> None:
    """e_p = 0.5, eps_p = -1, eps = +1: printed tau/pi = 0.1055."""
    assert m.parabolic_tau(0.5, -1) / PI == pytest.approx(0.1055, abs=1e-4)
    assert m.parabolic_tau(0.5, -1) / PI == pytest.approx(
        _parabola_first_principles(0.5, -1, 1), abs=1e-9
    )


def test_parabolic_arc_eps_p_plus_one_matches_first_principles() -> None:
    """e_p = 0.5, eps_p = +1: eq. 19 against a direct solve of the collision equations
    (independent of the closed form); the printed 0.2318 is the strict xfail below."""
    ours = m.parabolic_tau(0.5, 1) / PI
    assert ours == pytest.approx(_parabola_first_principles(0.5, 1, -1), abs=1e-9)
    assert ours == pytest.approx(0.2311, abs=1e-4)


@pytest.mark.xfail(
    strict=True,
    reason="Gomez & Olle p.43 print tau/pi = 0.2318 (Fig. 5: 0.2317); the equations give 0.231149",
)
def test_parabolic_arc_printed_0_2318() -> None:
    assert m.parabolic_tau(0.5, 1) / PI == pytest.approx(0.2318, abs=1e-4)


@pytest.mark.parametrize(
    ("i", "bound"),
    [(1, 1.539), (2, 16.0), (3, 4.617), (4, 32.0), (5, 7.695)],
)
def test_eq21_bounds_at_e_p_half(i: int, bound: float) -> None:
    """Eq. 21 numbers printed at e_p = 0.5 (eps_p = +1): i = 1..5."""
    assert m.c_family_max_j(i, 0.5, 1) == pytest.approx(bound, abs=6e-4 * bound)


def test_eq21_prefactors_at_e_p_half_and_098() -> None:
    assert pytest.approx(1.539, abs=1e-3) == (2.0 / 1.5) ** 1.5  # 1.5396 printed as 1.539
    assert pytest.approx(8.0, abs=1e-12) == (2.0 / 0.5) ** 1.5
    assert pytest.approx(1.015, abs=5e-4) == (2.0 / 1.98) ** 1.5
    assert (
        pytest.approx(1000.0, rel=1e-12) == (2.0 / 0.02) ** 1.5
    )  # printed "1.000", a lost separator


def test_eq21_families_that_appear_and_disappear_at_e_p_half() -> None:
    """p.44: at e_p = 0.5, eps_p = +1 the families C12, C35, C36 disappear (i odd) and C26 ..
    C2,16 appear (i even); at e_p = 0 Henon's 1 < j/i <= 2 sqrt 2 holds."""
    assert not m.c_family_exists(1, 2, 0.5, 1)
    assert not m.c_family_exists(3, 5, 0.5, 1)
    assert not m.c_family_exists(3, 6, 0.5, 1)
    assert m.c_family_exists(3, 4, 0.5, 1)
    assert m.c_family_exists(2, 6, 0.5, 1)
    assert m.c_family_exists(2, 16, 0.5, 1)
    assert not m.c_family_exists(2, 17, 0.5, 1)
    for i, j in itertools.product(range(1, 9), range(1, 30)):
        assert m.c_family_exists(i, j) == (1.0 < j / i <= 2.0 * math.sqrt(2.0) + 1e-12)


def test_eq21_fig7_eps_p_minus_one_at_e_p_half() -> None:
    """Fig. 7 (eps_p = -1, e_p = 0.5): i = 1 allows j <= 8 (C12 .. C17 drawn), i = 2 allows
    j <= 3.08 (only C23 drawn)."""
    assert m.c_family_max_j(1, 0.5, -1) == pytest.approx(8.0)
    assert m.c_family_max_j(2, 0.5, -1) == pytest.approx(3.08, abs=1e-2)
    assert [j for i, j in m.c_family_indices(2, 0.5, -1) if i == 2] == [3]
    assert max(j for i, j in m.c_family_indices(1, 0.5, -1) if i == 1) == 8


def test_first_c_families_at_e_p_098() -> None:
    """p.44: the first odd-i family is C67,68 (eps_p = +1) and the first even-i family is
    C68,69 (eps_p = -1), the latter with the rounded 1.015 of the text."""
    assert m.first_c_family(0.98, 1, 1) == (67, 68)
    assert m.first_c_family(0.98, -1, 0, factor=1.015) == (68, 69)


def test_first_even_c_family_with_exact_eq21_is_c66_67() -> None:
    """Exact eq. 21 puts the first even-i family at C66,67 (66 * 1.01523 = 67.005), not the
    printed C68,69: the printed value follows from the rounded factor 1.015."""
    assert m.first_c_family(0.98, -1, 0) == (66, 67)


@pytest.mark.xfail(
    strict=True, reason="printed C68,69 needs the rounded 1.015; exact eq. 21 gives C66,67"
)
def test_first_even_c_family_printed_c68_69_exact() -> None:
    assert m.first_c_family(0.98, -1, 0) == (68, 69)


@pytest.mark.parametrize(("i", "j"), [(2, 5), (3, 4), (2, 6), (1, 1)])
def test_tangent_ellipse_formula_agrees_with_the_timing_equation_at_e_p_half(
    i: int, j: int
) -> None:
    """The double point (i pi, j pi) of eq. 17 at e_p = 0.5 carries the ellipse a = (i/j)^(2/3),
    e = |1 - (j/i)^(2/3)(1 + (-1)^(i+1) eps_p e_p)| (text after eq. 20): the arcs of eq. 17 near
    the point (offset 2e-3 in tau) have those elements to a few 1e-3."""
    sigma = (-1) ** (i + j)
    a_pred, e_pred = m.tangent_ellipse(i, j, 0.5, 1)
    tau = i * PI + 2e-3
    roots = m.find_etas(
        tau, sigma, eta_min=(j - 0.3) * PI, eta_max=(j + 0.3) * PI, e_p=0.5, eps_p=1
    )
    near = [x for x in roots if abs(x - j * PI) < 0.2]
    assert near
    arc = m.arc_elements(tau, near[0], sigma, e_p=0.5, eps_p=1)
    assert arc.a == pytest.approx(a_pred, abs=5e-3)
    assert arc.e == pytest.approx(e_pred, abs=5e-3)


def test_circular_limit_of_the_elliptic_timing_equation() -> None:
    """Eq. 17 at e_p -> 0 is Henon's eq. 30: Table 2 row tau/pi 0.17 closes to the print."""
    assert abs(m.timing_residual(0.17 * PI, 0.16734 * PI, -1, e_p=1e-9)) < 1e-4
    assert abs(m.timing_residual(0.17 * PI, 0.16734 * PI, -1, e_p=0.0)) < 1e-4
    arc = m.arc_elements(0.17 * PI, 0.16734 * PI, -1, e_p=1e-9)
    assert arc.a == pytest.approx(6.92689, abs=5e-4)
    assert arc.e == pytest.approx(0.98922, abs=1e-5)


def test_elliptic_extension_has_no_jacobi_integral() -> None:
    arc = m.arc_elements(0.17 * PI, 0.16734 * PI, -1, eps1=-1, e_p=0.5)
    with pytest.raises(NotImplementedError):
        _ = arc.v1
    with pytest.raises(NotImplementedError):
        _ = arc.jacobi
    with pytest.raises(NotImplementedError):
        m.collision_data(arc)
