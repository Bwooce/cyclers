"""#935: printed tables of the Dvorak & Henrard volume (CMDA 56, 1993) against ``core.er3bp``.

Source and digest:
``docs/notes/2026-10-05-digest-dvorak-henrard-1993-cmda-56-selected-chapters.md``, section 13
Table A only (printed values); Table B (the digest agent's own runs) is not used as an expected
value anywhere here.

(a) Hadjidemetriou, "Resonant motion in the restricted three body problem", CMDA 56:201-219,
pp. 211, 216, 217 (Sun-Jupiter, mu = 0.00095387535). All 43 printed rows are stored in
``data/hadjidemetriou_cmda56_table_a.csv`` (39 distinct states: the e' = 0 rows of II_c, II_e and
IV_e repeat those of I_c, I_e and III_e, and the digest counts them in its tiers). Printed ydot0
is relative to a frame rotating at the
INSTANTANEOUS angular rate of the Sun-Jupiter line at t = 0 (with the mean rate 1 the elliptic rows
fail by 0.3 to 1.2), x0 is in Jupiter-orbit-radius units, the primaries are at perihelion (I) or
aphelion (II) at t = 0. That is Broucke's Eq. 62 rotating frame, so the pulsating state of
``core.er3bp`` is x = x0 / r, ydot = r * ydot0 / sqrt(1 - e'^2) with r = 1 -/+ e' (the apse has
r' = 0). Symmetric orbits start perpendicular to the x axis, so the half-period conditions are
y(pi) = 0 and x'(pi) = 0.

Tiers. A printed row is judged by the distance to the nearest exact symmetric orbit (a two-parameter
fsolve of the half-period conditions in (x0, ydot0) at fixed e', started at the printed row):
the 19 ``golden`` rows reproduce to the printed digits, distance <= 2e-5 (the class boundary in the
digest, section 5.3); the other rows are only a documented band: the digest reports an exact
orbit within 6e-3 of each, which is the digest agent's computation and not a printed value, so the
test uses the cheap full-period closure residual with a loose band of 0.1 (digest range 1e-3 to
6e-2 for these rows) as a sanity check, not as a golden. The e' = "0.010" row of 2:1 I_e (p.211)
is a misprint for 0.100 (closure scan, section 5.3): at the printed e' it fails the band by a
factor of 27 (strict expected failure), at 0.100 it closes to the neighbouring rows' level.
The printed 4.700496 (p.215) against 4.700478 (p.217) for the III_e e' = 0 orbit is an unresolved
slip in print; the table value 4.700478 is stored and used, and listed in the coarse tier.

(b) Hagel & Trenkler, "A computer aided analysis of the Sitnikov problem", CMDA 56:81-98, Table I
(p.86) and Table V (p.98): the linearised monodromy of the z motion at the barycentre, mu = 0.5,
R = [[r1, r2], [r3, r4]] = STM[2,2], STM[2,5], STM[5,2], STM[5,5] of ``propagate_er3bp`` over one
period, starting at periapsis f = 0 for e > 0 and at apoapsis f = pi for the table's negative e
(its convention: negative e means the primaries start at greatest separation). Tolerances as the
digest section 8.2 (it rounds the largest miss to 3e-4; it is 3.3e-4, e = -0.60 r3, printed
-1.2706 against -1.27027, so the test uses 4e-4): 4e-4 absolute for |e| <= 0.6, 3e-3 at |e| = 0.8.
Printed slips are strict expected failures: Tr R at e = -0.40; r4 at e = +-0.80 and Tr R at
e = -0.80; every column at e = +-0.99 (the source's RK4 with step 2 pi / 200 is not converged there). Table V's first row
(exact linear frequency Q(e) at T(0) = 0: 2.84802, 2.91273, 3.04723, 3.34258, 4.96398 at
e = 0.2 .. 0.99) gives cos(2 pi Q) = r1, an independent check that also holds at e = 0.99,
where Table I fails. Eq. (26), Q = sqrt(8) [1 + (24/121) e^2], does not reproduce Table I.

Marker choice. Repo convention: ``@pytest.mark.slow`` (skipped by default) for runs over about
10 s. The full 19 golden rows take about 6 s and the band rows about 5 s, so the default suite
runs 8 golden rows, the misprint pair and 6 band rows (about 6 s) and the remaining golden and band
rows are ``slow``.
"""

from __future__ import annotations

import csv
import math
import warnings
from pathlib import Path

import numpy as np
import pytest
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve

from cyclerfinder.core.er3bp import ER3BPSystem, er3bp_eom, propagate_er3bp

DATA = Path(__file__).parent / "data" / "hadjidemetriou_cmda56_table_a.csv"
MU_SUN_JUPITER = 0.00095387535
GOLDEN_TOL = 2e-5
BAND = 0.1


def _rows() -> list[dict[str, str]]:
    with DATA.open() as fh:
        return list(csv.DictReader(fh))


ROWS = _rows()


def _rid(r: dict[str, str]) -> str:
    return f"{r['resonance']}-{r['family']}-ep{r['eprime_printed']}"


def _to_pulsating(r: dict[str, str], e: float) -> tuple[float, float, float]:
    rr = 1.0 - e if r["phase"] == "peri" else 1.0 + e
    sp = math.sqrt(1.0 - e * e)
    f0 = 0.0 if r["phase"] == "peri" else math.pi
    return rr, sp, f0


def _end(xp: float, yp: float, e: float, f0: float, span: float) -> np.ndarray:
    sol = solve_ivp(
        er3bp_eom,
        (f0, f0 + span),
        [xp, 0.0, 0.0, 0.0, yp, 0.0],
        args=(MU_SUN_JUPITER, e),
        method="DOP853",
        rtol=1e-12,
        atol=1e-12,
    )
    return np.asarray(sol.y[:, -1])


def distance_to_exact_orbit(r: dict[str, str]) -> float:
    e = float(r["eprime_run"])
    rr, sp, f0 = _to_pulsating(r, e)
    x0, yd0 = float(r["x0"]), float(r["ydot0"])

    def half_period_residual(v: np.ndarray) -> list[float]:
        end = _end(v[0] / rr, rr * v[1] / sp, e, f0, math.pi)
        return [float(end[1]), float(end[3])]  # y and x' vanish at the perpendicular crossing

    with warnings.catch_warnings():
        # fsolve may report "not making good progress" once it sits at the xtol floor; only the
        # distance to the converged point is used, and it is checked against GOLDEN_TOL.
        warnings.simplefilter("ignore", RuntimeWarning)
        sol = fsolve(half_period_residual, [x0, yd0], xtol=1e-13)
    return float(max(abs(sol[0] - x0), abs(sol[1] - yd0)))


def full_period_closure(r: dict[str, str], e: float | None = None) -> float:
    e = float(r["eprime_run"]) if e is None else e
    rr, sp, f0 = _to_pulsating(r, e)
    xp, yp = float(r["x0"]) / rr, rr * float(r["ydot0"]) / sp
    end = _end(xp, yp, e, f0, 2.0 * math.pi)
    start = np.array([xp, 0.0, 0.0, 0.0, yp, 0.0])
    return float(np.max(np.abs(end - start)))


GOLDEN = [r for r in ROWS if r["tier"] == "golden"]
COARSE = [r for r in ROWS if r["tier"] == "coarse"]
GOLDEN_DEFAULT = [
    r
    for r in GOLDEN
    if _rid(r)
    in {
        "2:1-II_e-ep0.075",
        "4:1-I_c-ep0",
        "4:1-I_c-ep0.35",
        "4:1-II_c-ep0.03",
        "4:1-I_e-ep0.080",
        "4:1-I_e-ep0.090",
        "4:1-II_e-ep0.03",
        "4:1-III_e-ep0.15",
    }
]
GOLDEN_SLOW = [r for r in GOLDEN if r not in GOLDEN_DEFAULT]
COARSE_DEFAULT = [
    r
    for r in COARSE
    if _rid(r)
    in {
        "2:1-I_e-ep0",
        "2:1-II_e-ep0.048",
        "4:1-II_c-ep0.16",
        "4:1-III_e-ep0.08",
        "4:1-IV_e-ep0.020",
        "4:1-IV_e-ep0.041",
    }
]
COARSE_SLOW = [r for r in COARSE if r not in COARSE_DEFAULT]


def test_hadjidemetriou_table_a_row_counts() -> None:
    tiers = [r["tier"] for r in ROWS]
    assert len(ROWS) == 43
    assert tiers.count("golden") == 19  # the 19 rows the digest says reproduce to printed digits
    assert tiers.count("misprint") == 1
    assert len(GOLDEN_DEFAULT) == 8 and len(COARSE_DEFAULT) == 6


@pytest.mark.parametrize("row", GOLDEN_DEFAULT, ids=_rid)
def test_hadjidemetriou_golden_row_is_an_exact_orbit_to_printed_digits(row: dict[str, str]) -> None:
    assert distance_to_exact_orbit(row) <= GOLDEN_TOL


@pytest.mark.slow
@pytest.mark.parametrize("row", GOLDEN_SLOW, ids=_rid)
def test_hadjidemetriou_golden_rows_full(row: dict[str, str]) -> None:
    assert distance_to_exact_orbit(row) <= GOLDEN_TOL


@pytest.mark.parametrize("row", COARSE_DEFAULT, ids=_rid)
def test_hadjidemetriou_other_row_closes_within_band(row: dict[str, str]) -> None:
    assert full_period_closure(row) <= BAND


@pytest.mark.slow
@pytest.mark.parametrize("row", COARSE_SLOW, ids=_rid)
def test_hadjidemetriou_other_rows_full(row: dict[str, str]) -> None:
    assert full_period_closure(row) <= BAND


def _misprint_row() -> dict[str, str]:
    return next(r for r in ROWS if r["tier"] == "misprint")


@pytest.mark.xfail(
    strict=True,
    reason="Hadjidemetriou p.211, 2:1 I_e row printed e' = 0.010 is a misprint for 0.100: at e' = "
    "0.010 the printed (x0, ydot0) = (0.186246, 2.785537) has no orbit nearby (closure residual "
    "about 2.7, digest section 5.3); only e' = 0.100 closes.",
)
def test_hadjidemetriou_misprinted_eprime_row_as_printed() -> None:
    assert full_period_closure(_misprint_row(), e=0.010) <= BAND


def test_hadjidemetriou_misprinted_row_closes_at_corrected_eprime() -> None:
    row = _misprint_row()
    assert float(row["eprime_printed"]) == 0.010 and float(row["eprime_run"]) == 0.100
    assert full_period_closure(row) <= 3e-3  # digest scan: 8.6e-4 at e' = 0.100


def test_hadjidemetriou_frame_convention_mean_rate_fails() -> None:
    """With the mean rate (omega = 1, i.e. ydot0 used as the inertial offset), e' > 0 rows fail."""
    row = next(r for r in GOLDEN if _rid(r) == "4:1-I_c-ep0.35")
    e = 0.35
    # Wrong: the printed ydot0 as a pulsating-frame velocity (omega = 1, no r or sqrt(p) scaling).
    wrong = _end(float(row["x0"]), float(row["ydot0"]), e, 0.0, 2.0 * math.pi)
    start = np.array([float(row["x0"]), 0.0, 0.0, 0.0, float(row["ydot0"]), 0.0])
    assert float(np.max(np.abs(wrong - start))) > 0.05  # digest: 0.3 to 1.2 with omega = 1
    assert full_period_closure(row) < 1e-2


# ---- Hagel & Trenkler Table I (p.86) -------------------------------------------------------

SQRT8 = math.sqrt(8.0)
# e: (r1, r2, r3, r4, Tr R) as printed.
TABLE_I: dict[float, tuple[float, float, float, float, float]] = {
    -0.99: (0.9325, -0.0135, 9.6365, 0.9325, 1.8650),
    -0.80: (-0.5498, 0.1373, -5.0817, -0.5438, 1.0995),
    -0.60: (0.9563, 0.0675, -1.2706, 0.9563, 1.9125),
    -0.40: (0.8534, -0.1455, 1.8667, 0.8534, 1.6117),
    -0.20: (0.5777, -0.2606, 2.5563, 0.5777, 1.1554),
    0.20: (0.5777, -0.3131, 2.1280, 0.5777, 1.1554),
    0.40: (0.8534, -0.2139, 1.2703, 0.8534, 1.7068),
    0.60: (0.9563, 0.1271, -0.6733, 0.9563, 1.9125),
    0.80: (-0.5498, 0.3812, -1.8303, -0.5450, -1.0995),
    0.99: (0.9325, -0.1007, 0.4458, 0.9325, 1.8650),
}
COLUMNS = ("r1", "r2", "r3", "r4", "TrR")
# (e, column) -> reason, for printed slips (strict expected failures).
SLIPS: dict[tuple[float, str], str] = {
    (-0.40, "TrR"): "printed Tr R = 1.6117 at e = -0.40 is a slip (the table's own text says Q "
    "is symmetric in e; +0.40 prints 1.7068)",
    (-0.80, "r4"): "printed r4 = -0.5438 at e = -0.80 (det R = 1 and equal diagonals need r4 = r1)",
    (0.80, "r4"): "printed r4 = -0.5450 at e = +0.80 (equal diagonals need r4 = r1)",
    (-0.80, "TrR"): "printed Tr R = +1.0995 at e = -0.80 has the wrong sign (Tr R is symmetric "
    "in e and is -1.0989 at +0.80)",
}
for _e in (-0.99, 0.99):
    for _c in COLUMNS:
        SLIPS[(_e, _c)] = (
            "Hagel-Trenkler Table I rows e = +-0.99 come from RK4 with step 2 pi / 200, not "
            "converged near |e| = 1; Table V's Q(0.99) = 4.96398 agrees with the project instead"
        )


def monodromy(e: float) -> tuple[float, float, float, float]:
    f0 = 0.0 if e > 0 else math.pi
    sys_ = ER3BPSystem(mu=0.5, e=abs(e), primary_name="m1", secondary_name="m2")
    _, _, phi = propagate_er3bp(
        np.zeros(6), (f0, f0 + 2.0 * math.pi), sys_, rtol=1e-13, atol=1e-13, with_stm=True
    )
    return float(phi[2, 2]), float(phi[2, 5]), float(phi[5, 2]), float(phi[5, 5])


def _params() -> list[object]:
    out: list[object] = []
    for e in sorted(TABLE_I):
        for k, col in enumerate(COLUMNS):
            marks = []
            if (e, col) in SLIPS:
                marks.append(pytest.mark.xfail(strict=True, reason=SLIPS[(e, col)]))
            out.append(pytest.param(e, k, col, marks=marks, id=f"e{e:+.2f}-{col}"))
    return out


@pytest.mark.parametrize(("e", "k", "col"), _params())
def test_sitnikov_table_i_against_er3bp_monodromy(e: float, k: int, col: str) -> None:
    r1, r2, r3, r4 = monodromy(e)
    computed = (r1, r2, r3, r4, r1 + r4)[k]
    tol = 4e-4 if abs(e) <= 0.6 else 3e-3
    assert computed == pytest.approx(TABLE_I[e][k], abs=tol)


def test_sitnikov_e0_closed_form_row() -> None:
    r1, r2, r3, r4 = monodromy(1e-12)
    w = 2.0 * math.pi * SQRT8
    assert r1 == pytest.approx(math.cos(w), abs=1e-6)
    assert r2 == pytest.approx(math.sin(w) / SQRT8, abs=1e-6)
    assert r3 == pytest.approx(-8.0 * r2, abs=1e-6)
    assert r1 + r4 == pytest.approx(0.9461, abs=5e-5)  # printed Tr R at e = 0


@pytest.mark.parametrize(
    ("e", "q_printed"),
    [(0.2, 2.84802), (0.4, 2.91273), (0.6, 3.04723), (0.8, 3.34258), (0.99, 4.96398)],
)
def test_sitnikov_table_v_frequency_matches_monodromy(e: float, q_printed: float) -> None:
    """cos(2 pi Q) = r1 with Q the printed Table V value (first row, T(0) = 0)."""
    r1, _, _, r4 = monodromy(e)
    assert math.cos(2.0 * math.pi * q_printed) == pytest.approx(r1, abs=1e-4)
    assert r1 == pytest.approx(r4, abs=1e-9)


@pytest.mark.xfail(
    strict=True,
    reason="Hagel-Trenkler eq. (26), Q = sqrt(8) [1 + (24/121) e^2], does not reproduce Table I "
    "(2.8509 against 2.8480 at e = 0.20; 2.9182 against 2.9127 at e = 0.40); the second-order "
    "coefficient is 21/124 (digest section 8.3, transposed digits).",
)
@pytest.mark.parametrize(("e", "q_table_i"), [(0.20, 2.8480), (0.40, 2.9127)])
def test_sitnikov_eq26_as_printed_reproduces_table_i(e: float, q_table_i: float) -> None:
    assert SQRT8 * (1.0 + (24.0 / 121.0) * e * e) == pytest.approx(q_table_i, abs=5e-4)
