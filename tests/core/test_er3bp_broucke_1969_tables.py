"""#933: Broucke 1969 (JPL TR 32-1360) printed orbit tables closed in ``core.er3bp``.

Source: R. A. Broucke, "Periodic Orbits in the Elliptic Restricted Three-Body Problem", JPL
Technical Report 32-1360 (15 July 1969), NTRS 19700005781. The rows are transcribed in the digests
``docs/notes/2026-10-04-digest-broucke-1969-elliptic-periodic-orbits-part-a.md`` and ``-part-b.md``
and stored in ``tests/core/data/broucke_1969_tables.csv`` (printed columns only; the digests' own
residual columns are not stored, so no expected value here was computed by project code).

Tables covered (846 transcribed rows; the two e = 1 rows are skipped):

* Tables 14 to 18 (part B, 431 rows): families 8P, 8A, 11P, 11A (mu = 0.5, start (x0, 0, 0, ydot0))
  and 10P (start (0, y0, xdot0, 0), quarter orbit). Rows were read from page images directly.
* Table 12 (7P, 131 rows) and Table 13 (7A, 120 rows), mu = 0.012155 exactly (Earth-Moon).
* Table 9 (9P, mu = 0.5, 11 rows with e < 1; inertial) and Table 5 (12A, mu = 0.5, 151 rows with
  e < 1; inertial, apoapsis start, period 4 pi, no printed end state: only the perpendicular
  crossing at the half period f = pi -> 3 pi is tested). Both use the conversion printed in
  Broucke Eqs. 62 and 64 and derived in digest part A section 3.1.

Not tested: the rectilinear tables (4, 6, 7, 8, 10, 11, and the e = 1 rows of Tables 5 and 9)
because core.er3bp cannot start at the collision (division by r2^3 and 1 + e cos f = 0 at e = 1) and
a rectilinear inertial integrator is a separate piece of work (#928). Their known defects (Table 6
truncates Table 7 to six decimals and numbers orbits 1 and 2 the other way round; Table 11 rows 9
and 10 are out of order; Tables 8, 10, 11 use the mirror of the frame of eq. 36c,d) are therefore
recorded here, not exercised. Table 19 (69 collision rows) needs the Birkhoff-regularised
propagator (#928) and is likewise not tested.

Provenance caveat for part A (Tables 5, 9, 12, 13): the digits were obtained by per-cell optical
recognition with majority voting, then checked by integration and column smoothness, and flagged
rows were corrected from the page image; a seventh-digit error on a weakly sensitive row could
survive. Part B was read from the images directly. This is why part A passes with a rounding bound
and not with a tighter independent value.

Tolerance. A printed row has seven decimals, so each printed start component carries up to 5e-8 of
rounding. The closure residual (y and x' at the end, plus x1 and ydot1 against the printed end
state) is therefore expected to be bounded by the state transition matrix of the half revolution
applied to that rounding, plus 5e-8 for the printed end state. Each row's tolerance is TOL_FACTOR
(= 2) times that bound, computed with the project's own variational equations (a scale, never a
value being tested). With this, the 8P rows with e >= 0.815, whose residual grows with e to 1e-2,
pass because the bound grows with them (rounding amplification, as the digest argues), while a
real digit error of 1e-6 or more on an insensitive row fails. Measured ratios residual / bound
were at most 1.5 for every row of Tables 9, 12 to 18 except Table 17 row 93 (e = 0.895, the
printed close approach to a primary where the linearised bound is not valid) and Table 5.
Table 5 rows near e = 1 close only to about 1e-4 (the report states that "five to six places"
are required for the inertial continuation, so its printed digits are not all significant): its
tolerance is the fixed 1e-4 and Table 17 row 93 gets 2e-4, both listed in OVERRIDES.

Known printed defects encoded here: Table 18 rows 55 and 56 have YDOT1 and X1 swapped (the printed
rows are strict expected failures, the corrected rows pass); Table 15 rows 75 to 77 are out of
eccentricity order (asserted, and each row still closes at its printed e); the 8P e = 0 row equals
the 8A e = 0 row, likewise 11P and 11A.

Marker choice. The repo marks tests longer than about 10 s ``slow`` (skipped by default,
``-m slow`` runs them). The full sweep of the 842 testable rows (the two faulty Table 18 rows excluded) takes about 25 s, so it is
``@pytest.mark.slow``; the default suite runs a sample (both ends of each table, every 10th row,
every flagged row, every 8P row with e >= 0.815 within the sample) in about 5 s.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path

import numpy as np
import pytest

from cyclerfinder.core.er3bp import ER3BPSystem, propagate_er3bp

DATA = Path(__file__).parent / "data" / "broucke_1969_tables.csv"
ROUND = 0.5e-7  # half a unit in the seventh printed decimal
TOL_FACTOR = 2.0
# (table, nr) -> absolute tolerance replacing TOL_FACTOR * bound, with the reason in the docstring.
OVERRIDES: dict[tuple[str, int], float] = {("17", 93): 2e-4}
TABLE5_TOL = 1e-4


def _load() -> list[dict[str, str]]:
    with DATA.open() as fh:
        return list(csv.DictReader(fh))


ROWS = _load()
ROWS_E_LT_1 = [r for r in ROWS if float(r["e"]) < 1.0]


def _row(table: str, nr: int) -> dict[str, str]:
    return next(r for r in ROWS if r["table"] == table and int(r["nr"]) == nr)


def closure(
    row: dict[str, str], override_end: tuple[float, float] | None = None
) -> tuple[float, float]:
    """Return (largest residual, rounding bound) for one printed row."""
    mu = float(row["mu"])
    e = float(row["e"])
    a = float(row["start_a"])
    b = float(row["start_b"])
    fam = row["family"]
    frame = row["frame"]
    sysm = ER3BPSystem(mu=mu, e=e, primary_name="m1", secondary_name="m2")
    p = 1.0 - e * e
    sp = math.sqrt(p)
    rp, ra = 1.0 - e, 1.0 + e
    if frame == "pulsating":
        if row["kind"] == "x0ydot0":
            s0 = [a, 0.0, 0.0, 0.0, b, 0.0]
            cols = (0, 4)
        else:  # 10P: start on the y axis
            s0 = [0.0, a, 0.0, b, 0.0, 0.0]
            cols = (1, 3)
        dx = (ROUND, ROUND)
        f0 = 0.0 if fam.endswith("P") else math.pi
        span = math.pi
    elif frame == "inertial_peri":
        s0 = [a / rp, 0.0, 0.0, 0.0, -a / rp + rp * b / sp, 0.0]
        dx = (ROUND / rp, ROUND / rp + ROUND * rp / sp)
        cols = (0, 4)
        f0, span = 0.0, math.pi
    elif frame == "inertial_apo":
        s0 = [-a / ra, 0.0, 0.0, 0.0, a / ra - ra * b / sp, 0.0]
        dx = (ROUND / ra, ROUND / ra + ROUND * ra / sp)
        cols = (0, 4)
        f0, span = math.pi, 2.0 * math.pi
    else:
        raise ValueError(frame)
    _, y, phi = propagate_er3bp(
        np.array(s0), (f0, f0 + span), sysm, rtol=1e-13, atol=1e-13, with_stm=True
    )
    end = y[:, -1]
    res = [abs(end[1]), abs(end[3])]
    if row["x1"] != "":
        x1p, yd1p = (
            (float(row["x1"]), float(row["ydot1"])) if override_end is None else override_end
        )
        if frame == "inertial_peri":
            res += [abs(-ra * end[0] - x1p), abs(-(sp / ra) * (end[4] + end[0]) - yd1p)]
        else:
            res += [abs(end[0] - x1p), abs(end[4] - yd1p)]
    bound = 0.0
    for k in (1, 3):
        bound = max(bound, abs(phi[k, cols[0]]) * dx[0] + abs(phi[k, cols[1]]) * dx[1])
    out_round = ROUND * (1.0 if frame == "pulsating" else 3.0)
    for k in (0, 4):
        bound = max(bound, abs(phi[k, cols[0]]) * dx[0] + abs(phi[k, cols[1]]) * dx[1] + out_round)
    return max(res), bound


def _tolerance(row: dict[str, str], bound: float) -> float:
    if row["table"] == "5":
        return TABLE5_TOL
    return OVERRIDES.get((row["table"], int(row["nr"])), TOL_FACTOR * bound)


def _check(row: dict[str, str]) -> None:
    res, bound = closure(row)
    tol = _tolerance(row, bound)
    assert res <= tol, (
        f"table {row['table']} row {row['nr']} e={row['e']}: residual {res:.3e} > tol {tol:.3e}"
    )


def _sample() -> list[tuple[str, int]]:
    chosen: set[tuple[str, int]] = set()
    for table in ("5", "9", "12", "13", "14", "15", "16", "17", "18"):
        rs = [r for r in ROWS_E_LT_1 if r["table"] == table]
        nrs = [int(r["nr"]) for r in rs]
        chosen.update((table, n) for n in (nrs[:2] + nrs[-2:]))
        chosen.update((table, n) for n in nrs if n % 10 == 0 and (table != "5" or n % 20 == 0))
    # flagged rows
    chosen.update(("15", n) for n in (74, 75, 76, 77, 78))  # eccentricity order slip
    chosen.add(("17", 93))  # close approach, e = 0.895
    chosen.update(("14", n) for n in (92, 93, 95, 100, 104, 110, 114, 118))  # 8P, e >= 0.815 tail
    chosen.update(("18", n) for n in (53, 54, 57, 58))  # neighbours of the printer fault
    chosen.discard(("18", 55))
    chosen.discard(("18", 56))
    return sorted(chosen, key=lambda t: (int(t[0]), t[1]))


SAMPLE = _sample()


def test_row_counts_match_table_20_and_printed_tables() -> None:
    counts = {t: sum(1 for r in ROWS if r["table"] == t) for t in {r["table"] for r in ROWS}}
    assert counts == {
        "5": 152,
        "9": 12,
        "12": 131,
        "13": 120,
        "14": 118,
        "15": 92,
        "16": 52,
        "17": 93,
        "18": 76,
    }


def test_first_rows_of_p_and_a_families_are_identical() -> None:
    # Tables 14 and 15, and 16 and 17, start from the same e = 0 orbit.
    for p_table, a_table in (("14", "15"), ("16", "17")):
        p, a = _row(p_table, 1), _row(a_table, 1)
        for key in ("e", "start_a", "start_b", "x1", "ydot1"):
            assert p[key] == a[key]


def test_table_15_rows_75_to_77_are_out_of_eccentricity_order_as_printed() -> None:
    e = [float(_row("15", n)["e"]) for n in range(74, 79)]
    assert e == pytest.approx([0.67, 0.70, 0.68, 0.69, 0.71])
    # Every other row of the table is in increasing order; each of 75 to 77 still closes at its
    # printed e (sampled below), so this is a print-order slip and not a digit error.
    others = [float(r["e"]) for r in ROWS if r["table"] == "15" and not 75 <= int(r["nr"]) <= 77]
    assert others == sorted(others)


@pytest.mark.parametrize(("table", "nr"), SAMPLE, ids=[f"T{t}-row{n}" for t, n in SAMPLE])
def test_broucke_sample_row_closes_in_er3bp(table: str, nr: int) -> None:
    _check(_row(table, nr))


def test_wrong_apse_start_fails_by_order_one() -> None:
    """Negative control: an apoapsis family started at the periapsis epoch misses by order 1."""
    row = dict(_row("15", 30))
    row["family"] = "8P"  # same state, wrong start epoch f = 0
    res, bound = closure(row)
    assert res > 0.1 > 1000 * bound


def test_8p_tail_residual_is_rounding_amplification_not_a_defect() -> None:
    low, _ = closure(_row("14", 60))
    high, bound_high = closure(_row("14", 118))
    assert low < 1e-6
    assert high > 1e-3  # e = 0.975: residual grows to about 1e-2
    assert high <= TOL_FACTOR * bound_high  # and is inside the rounding-amplification bound
    assert bound_high > 1e3 * ROUND


def _corrected_table_18_rows() -> tuple[dict[str, str], dict[str, str]]:
    r55, r56 = dict(_row("18", 55)), dict(_row("18", 56))
    # Printed: row 55 YDOT1 = -1.0908355 (the size of X1) and row 56 X1 = -0.4151629 (the size
    # of YDOT1); the two values are swapped between the rows.
    r55["ydot1"], r56["x1"] = r56["x1"], r55["ydot1"]
    return r55, r56


def test_table_18_swap_is_the_printed_defect_and_corrected_rows_close() -> None:
    r55, r56 = _corrected_table_18_rows()
    assert float(_row("18", 55)["ydot1"]) == pytest.approx(-1.0908355)
    assert float(_row("18", 56)["x1"]) == pytest.approx(-0.4151629)
    for row in (r55, r56):
        res, bound = closure(row)
        assert res <= TOL_FACTOR * bound


@pytest.mark.xfail(
    strict=True,
    reason="Broucke Table 18 rows 55 and 56 (p.80) print YDOT1 and X1 swapped between the rows "
    "(row 55 YDOT1 = -1.0908355 is the size of X1; row 56 X1 = -0.4151629 is the size of YDOT1, "
    "and row 56 also prints mass ratio 0). Integrating the printed starts gives the swapped "
    "values to better than 1e-6 (digest part B section 6 item 1). The printed rows must fail.",
)
@pytest.mark.parametrize("nr", [55, 56])
def test_table_18_printed_rows_55_56_as_printed(nr: int) -> None:
    _check(_row("18", nr))


@pytest.mark.slow
@pytest.mark.parametrize(
    ("table", "nr"),
    [
        (r["table"], int(r["nr"]))
        for r in ROWS_E_LT_1
        if not (r["table"] == "18" and r["nr"] in ("55", "56"))
    ],
)
def test_broucke_full_table_sweep_closes_in_er3bp(table: str, nr: int) -> None:
    _check(_row(table, nr))
