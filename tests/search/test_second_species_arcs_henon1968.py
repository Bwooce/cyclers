"""Henon 1968 Tables 1-9, the parabolic arc and the closed-form special solutions.

Every EXPECTED value here is a number printed in the paper (transcribed in
``fixtures/second_species_arcs/henon1968_*.csv`` from
``docs/notes/2026-10-04-digest-henon-1968-consecutive-collision-orbits.md``);
nothing is a value our own code computed.  The three printed defects the digest
found are strict expected failures; the corrected quantities are checked against
the printed values of the other columns of the same row.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path

import pytest

from cyclerfinder.search import second_species_arcs as m

PI = math.pi
FX = Path(__file__).parent / "fixtures" / "second_species_arcs"

SIGMA = {"A0": -1, "A1": -1, "A2": -1, "B1": 1, "B2": 1, "C12": -1, "C23": -1, "C24": 1}

# Rows whose printed values contradict the paper's own equations (digest section 2.3).
DEFECTS = {
    ("3", "1.80000", "1.24723"),  # eta/pi printed 1.24723, equations give 1.14723
    ("5", "3.70000", "1.28122"),  # C printed 1.31608, V = 1.08808 gives 1.81607
    ("6", "1.90000", "2.10810"),  # x1 printed 0.85698, eq. 32 gives -0.85698
}


def _rows() -> list[dict[str, str]]:
    with open(FX / "henon1968_tables2_9.csv") as f:
        return list(csv.DictReader(f))


ROWS = _rows()


def _key(r: dict[str, str]) -> tuple[str, str, str]:
    return (r["table"], r["tau_pi"], r["eta_pi"])


def _arc_for(r: dict[str, str], eps1: int | None = None) -> m.SArc:
    tau_pi, eta_pi = float(r["tau_pi"]), float(r["eta_pi"])
    sigma = SIGMA[r["family"]]
    if abs(math.sin(eta_pi * PI)) < 1e-12:
        return m.tangent_arc(round(tau_pi), round(eta_pi), eps1=eps1 or 1)
    return m.arc_elements(tau_pi * PI, eta_pi * PI, sigma, eps1=eps1)


def _senses(arc_eps1: int, r: dict[str, str]) -> list[int]:
    """Senses to try: the inferred one, or both on a tangent-ellipse double point."""
    if abs(math.sin(float(r["eta_pi"]) * PI)) < 1e-12:
        return [1, -1]
    return [arc_eps1]


def _tol(printed: float) -> float:
    # printed to 5 decimals; a and x1 are large near tau/pi = 0.17, hence the relative part
    return 1.0e-3 + 1.0e-4 * abs(printed)


def _compare_row(r: dict[str, str]) -> list[str]:
    """Names of columns that disagree with the printed row (best sense chosen)."""
    base = _arc_for(r)
    worst: list[str] | None = None
    for e1 in _senses(base.eps1, r):
        arc = _arc_for(r, e1)
        bad = []
        for col, val in (
            ("a", arc.a),
            ("e", arc.e),
            ("x0", arc.x0),
            ("x1", arc.x1),
            ("C", arc.jacobi),
            ("V", arc.speed),
        ):
            p = float(r[col])
            if abs(val - p) > _tol(p):
                bad.append(col)
        if worst is None or len(bad) < len(worst):
            worst = bad
    assert worst is not None
    return worst


PARABOLIC_ROW = [r for r in ROWS if r["eta_pi"].startswith("0.00000")]
ELLIPTIC_ROWS = [r for r in ROWS if not r["eta_pi"].startswith("0.00000")]
GOOD_ROWS = [r for r in ELLIPTIC_ROWS if _key(r) not in DEFECTS]


def _by_key(key: tuple[str, str, str]) -> dict[str, str]:
    return next(r for r in ELLIPTIC_ROWS if _key(r) == key)


def test_row_counts_match_the_paper() -> None:
    per_table: dict[str, int] = {}
    for r in ROWS:
        per_table[r["table"]] = per_table.get(r["table"], 0) + 1
    assert per_table == {"2": 42, "3": 58, "4": 39, "5": 65, "6": 47, "7": 26, "8": 31, "9": 16}
    assert sum(per_table.values()) == 324


@pytest.mark.parametrize("r", GOOD_ROWS, ids=lambda r: f"T{r['table']}-{r['tau_pi']}-{r['eta_pi']}")
def test_henon_tables_2_to_9_row_closes_from_tau_eta(r: dict[str, str]) -> None:
    """a, e, x0, x1, C and V of every printed row follow from (tau, eta) and the family sign
    through eqs. 29, 32, 33; the timing equation residual is small at the printed point."""
    assert _compare_row(r) == []
    tau, eta = float(r["tau_pi"]) * PI, float(r["eta_pi"]) * PI
    if abs(math.sin(eta)) > 1e-12:  # at a double point F has a vanishing gradient, not a root test
        # printed (tau, eta) carry 5 decimals; F changes by O(1) per radian away from the curve
        assert abs(m.timing_residual(tau, eta, SIGMA[r["family"]])) < 3.0e-4 * (1.0 + abs(tau))


@pytest.mark.parametrize("r", GOOD_ROWS, ids=lambda r: f"T{r['table']}-{r['tau_pi']}-{r['eta_pi']}")
def test_henon_row_lies_on_the_timing_curve(r: dict[str, str]) -> None:
    """Moving the printed (tau, eta) onto F = 0 (along the better-conditioned axis)
    changes it by less than the printed rounding, except at two ill-conditioned
    near-vertical-tangent rows."""
    tau, eta = float(r["tau_pi"]) * PI, float(r["eta_pi"]) * PI
    if abs(math.sin(eta)) < 1e-12 and abs(math.sin(tau)) < 1e-12:
        pytest.skip("double point: F has zero gradient (eq. 40-43 cover it)")
    t1, h1 = m.refine_point(tau, eta, SIGMA[r["family"]])
    shift = max(abs(t1 - tau), abs(h1 - eta)) / PI
    # The curve is nearly tangent to both axes at the sharp bends near (i + 1/2, i + 1/2)
    # (eq. 45), where a 5-decimal print cannot fix the point better than about 0.01.
    limit = 1.0e-2 if (r["table"], r["tau_pi"]) == ("4", "2.45500") else 2.0e-4
    assert shift < limit


def test_corrected_defect_rows_are_consistent_with_the_other_columns() -> None:
    """Each corrected quantity reproduces the printed values of the rest of the row."""
    # Table 3, tau/pi = 1.8: the root of eq. 30 is eta/pi 1.14723; the other columns are on the row
    row = _by_key(("3", "1.80000", "1.24723"))
    arc = m.arc_elements(1.8 * PI, 1.14723 * PI, -1)
    for col, val in (
        ("a", arc.a),
        ("e", arc.e),
        ("x0", arc.x0),
        ("x1", arc.x1),
        ("C", arc.jacobi),
        ("V", arc.speed),
    ):
        assert abs(val - float(row[col])) < _tol(float(row[col])), col
    # Table 5, tau/pi = 3.7: printed V = 1.08808 means C = 3 - V^2 = 1.81607 (the smooth sequence)
    row5 = _by_key(("5", "3.70000", "1.28122"))
    assert abs((3.0 - float(row5["V"]) ** 2) - 1.81607) < 1e-4
    arc5 = _arc_for(row5)
    assert abs(arc5.jacobi - 1.81607) < 1e-3
    # Table 6, tau/pi = 1.9: x1 = -0.85698 by eq. 32 (a dropped minus sign in the print)
    row6 = _by_key(("6", "1.90000", "2.10810"))
    arc6 = _arc_for(row6)
    assert abs(arc6.x1 - (-float(row6["x1"]))) < _tol(0.85698)


def test_the_three_printed_defects_are_exactly_these() -> None:
    """Over all 324 rows only the three digest-listed rows disagree with the equations."""
    failing = [_key(r) for r in ELLIPTIC_ROWS if _compare_row(r)]
    assert set(failing) == DEFECTS


# --- strict expected failures: printed values compared as printed --------------------------
@pytest.mark.xfail(
    strict=True, reason="Table 3 tau/pi 1.8: printed eta/pi 1.24723 (equations give 1.14723)"
)
def test_table3_printed_eta_1_24723() -> None:
    row = _by_key(("3", "1.80000", "1.24723"))
    assert _compare_row(row) == []


@pytest.mark.xfail(
    strict=True, reason="Table 5 tau/pi 3.7: printed C 1.31608 (V = 1.08808 gives 1.81607)"
)
def test_table5_printed_c_1_31608() -> None:
    row = _by_key(("5", "3.70000", "1.28122"))
    assert _compare_row(row) == []


@pytest.mark.xfail(
    strict=True, reason="Table 6 tau/pi 1.9: printed x1 +0.85698 (eq. 32 gives -0.85698)"
)
def test_table6_printed_x1_plus_0_85698() -> None:
    row = _by_key(("6", "1.90000", "2.10810"))
    assert _compare_row(row) == []


# --- parabolic arc (eqs. 16-22) -------------------------------------------------------------
def test_parabolic_arc_matches_printed_eqs_21_22() -> None:
    par = m.parabolic_arc()
    assert par.tau / PI == pytest.approx(0.16393, abs=6e-6)
    assert par.x0 == pytest.approx(-0.06485, abs=6e-6)
    assert par.jacobi == pytest.approx(-0.72028, abs=6e-6)
    assert par.speed == pytest.approx(1.92880, abs=6e-6)
    assert par.sigma_param == pytest.approx(3.7973, abs=1e-4)
    # the first row of Table 2 is the same orbit
    row = PARABOLIC_ROW[0]
    assert par.tau / PI == pytest.approx(float(row["tau_pi"]), abs=6e-6)
    assert par.x0 == pytest.approx(float(row["x0"]), abs=6e-6)
    assert par.jacobi == pytest.approx(float(row["C"]), abs=6e-6)
    assert par.speed == pytest.approx(float(row["V"]), abs=6e-6)


def test_parabolic_tau_with_eps_p_minus_one_is_the_rotated_figure() -> None:
    """Gomez & Olle: at e_p = 0 and eps_p = -1 the situation is the same turned by pi."""
    assert m.parabolic_tau(0.0, -1) == pytest.approx(m.parabolic_tau(0.0, 1), abs=1e-12)


# --- Table 1, hyperbolic orbits (eqs. 5-14) ------------------------------------------------
def _hyp_rows() -> list[dict[str, str]]:
    with open(FX / "henon1968_table1_hyperbolic.csv") as f:
        return list(csv.DictReader(f))


@pytest.mark.parametrize("r", _hyp_rows(), ids=lambda r: f"tau{r['tau_pi']}")
def test_table1_hyperbolic_rows(r: dict[str, str]) -> None:
    tau = float(r["tau_pi"]) * PI
    eta = float(r["eta"])
    # rows 0.01 and 0.16 are ill conditioned (digest 2.1): a ~ 1e-3 and a ~ 9.35
    loose = r["tau_pi"] in {"0.01000", "0.16000"}
    h = m.hyperbolic_arc(tau, eta)
    assert h.a == pytest.approx(float(r["a"]), abs=2e-4 if loose else 2e-5)
    assert h.e == pytest.approx(float(r["e"]), abs=2e-5)
    assert h.x0 == pytest.approx(float(r["x0"]), abs=2e-5)
    assert h.speed == pytest.approx(float(r["V"]), abs=2e-4 if loose else 2e-5)
    assert h.jacobi == pytest.approx(
        float(r["C"]), abs=6e-3 if loose else max(2e-4, 2e-6 * abs(float(r["C"])))
    )
    # eq. 6, third: tau = a^(3/2) (e sinh eta - eta)
    assert h.a**1.5 * (h.e * math.sinh(eta) - eta) == pytest.approx(tau, abs=1e-6)
    # solving eq. 10 for eta from tau alone reproduces the printed eta
    assert m.hyperbolic_arc(tau).eta == pytest.approx(eta, abs=2e-5)


# --- eq. 41 tangent ellipses and e = 1 rows ------------------------------------------------
@pytest.mark.parametrize(
    ("i", "j", "a", "e"),
    [
        (2, 1, 1.58740, 0.37004),  # Table 2 row 2.0
        (4, 1, 2.51984, 0.60315),  # Table 2 last row / Table 3 first row
        (2, 3, 0.76314, 0.31037),  # Table 8 first row
        (1, 2, 0.62996, 0.58740),  # Table 7 first row
    ],
)
def test_tangent_ellipse_eq41(i: int, j: int, a: float, e: float) -> None:
    ta, te = m.tangent_ellipse(i, j)
    assert ta == pytest.approx(a, abs=6e-6)
    assert te == pytest.approx(e, abs=6e-6)


def test_tangent_ellipse_exists_only_below_two_root_two() -> None:
    """Eq. 42: e <= 1 requires j/i <= 2 sqrt 2."""
    assert m.tangent_ellipse(5, 14)[1] <= 1.0  # 14/5 = 2.8 < 2.828
    with pytest.raises(ValueError):
        m.tangent_ellipse(1, 3)


def test_b1_e_equal_one_interior_orbit() -> None:
    """Table 5 e = 1 row at tau/pi = 1: eta/pi 1.36836, a = 0.71333 (Bruno's interior orbit)."""
    (arc,) = [a for a in m.e1_arcs(1, 1) if abs(a.eta / PI - 1.36836) < 1e-3]
    assert arc.eta / PI == pytest.approx(1.36836, abs=6e-6)
    assert arc.a == pytest.approx(0.71333, abs=6e-6)
    assert arc.e == 1.0


@pytest.mark.parametrize(("tau_pi", "eta_pi"), [(0.17, 0.16734), (0.30, 0.48269)])
def test_a0_eta_root_from_tau(tau_pi: float, eta_pi: float) -> None:
    """Digest 2.4: the root of eq. 30 (sigma = -1) at tau/pi 0.17 and 0.3 against the print."""
    roots = [x / PI for x in m.find_etas(tau_pi * PI, -1, eta_max=0.6 * PI)]
    assert any(abs(x - eta_pi) < 6e-6 for x in roots)


def test_a0_family_is_the_only_hyperbolic_continuation() -> None:
    """Eq. 11: solutions of the hyperbolic equation exist only for eps = eps' = -1; the
    A0 hyperbolic arc at the table's tau/pi = 0.16 joins the parabolic orbit at tau/pi 0.16393."""
    h = m.hyperbolic_arc(0.16 * PI)
    assert h.eta < 0.5  # eta -> 0 at the parabolic orbit
    assert m.parabolic_arc().tau / PI > 0.16


# --- enumerator completeness: every printed row is found by find_etas / enumerate_arcs -----------
def _groups_by_tau_and_sigma() -> dict[tuple[float, int], list[tuple[tuple[str, str, str], float]]]:
    groups: dict[tuple[float, int], list[tuple[tuple[str, str, str], float]]] = {}
    for r in ELLIPTIC_ROWS:
        if _key(r) in DEFECTS:
            continue
        tp, ep = float(r["tau_pi"]), float(r["eta_pi"])
        if abs(math.sin(ep * PI)) < 1e-9 and abs(math.sin(tp * PI)) < 1e-9:
            continue  # double point (tangent ellipse): F has a vanishing gradient there
        groups.setdefault((tp, SIGMA[r["family"]]), []).append((_key(r), ep))
    return groups


def test_enumerator_finds_every_printed_row() -> None:
    """For each distinct printed (tau/pi, sigma) of Tables 2-9, every printed eta/pi appears
    among the roots of the timing equation within 2e-4 (298 rows in 125 groups).  This
    includes Table 8 tau/pi 2.43883, the printed turning point, where two roots (2.50301 and
    2.50352) sit inside one 0.0017 pi grid cell."""
    groups = _groups_by_tau_and_sigma()
    assert len(groups) == 125
    n = 0
    missing = []
    for (tp, sg), lst in groups.items():
        roots = [x / PI for x in m.find_etas(tp * PI, sg, eta_max=7.0 * PI)]
        for key, ep in lst:
            n += 1
            if min((abs(x - ep) for x in roots), default=9.0) > 2e-4:
                missing.append(key)
    assert n == 298
    assert missing == []


def test_enumerator_resolves_two_roots_inside_one_grid_cell() -> None:
    """Near the C23 turning point (tau/pi about 2.43884) the two roots at eta/pi 2.5030 and
    2.5035 are 5e-4 pi apart, closer than the grid spacing 1.7e-3 pi; the printed Table 8 row
    tau/pi 2.43883, eta/pi 2.50300 is the lower one."""
    roots = [x / PI for x in m.find_etas(2.43883 * PI, -1, eta_max=2.6 * PI, eta_min=2.4 * PI)]
    assert len(roots) == 2
    assert roots[1] - roots[0] < 1.7e-3
    assert roots[0] == pytest.approx(2.50300, abs=2e-5)
    assert m.find_etas(2.4390 * PI, -1, eta_max=2.6 * PI, eta_min=2.4 * PI) == []


def test_enumerate_arcs_returns_arcs_with_both_signs() -> None:
    """enumerate_arcs at tau/pi 0.17 returns the single A0 arc (Table 2 row 2) and, for sigma = +1,
    the B-type arcs of the same duration, none of which is the coincident circle."""
    arcs = m.enumerate_arcs(0.17 * PI, eta_max=7.0 * PI)
    a0 = [a for a in arcs if a.sigma == -1 and a.eta / PI < 0.5]
    assert len(a0) == 1
    assert a0[0].a == pytest.approx(6.92689, abs=5e-4)
    assert all(not (a.sigma == 1 and abs(a.a - 1.0) < 1e-9 and a.e < 1e-9) for a in arcs)
