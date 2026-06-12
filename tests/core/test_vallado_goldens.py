"""Published golden values from Vallado, D. A., TR-91-6 (1991), Appendix E.

Source: Vallado, D. A., *Methods of Astrodynamics, A Computer Approach*,
Version 3.0, USAFA Technical Report TR-91-6, 1991 (DTIC AD-A239 662).
All Appendix E cases are in Earth canonical units (DU, TU, mu = 1), so they
exercise ``lambert`` / ``propagate`` / ``_stumpff`` unit-independently.

Golden discipline: every EXPECTED value below is transcribed from the printed
report (table/example and page cited per case); none was produced by our code.
Each wired cell was reproduced against the repo implementation before wiring
(load-bearing-numbers verification pass #203, PDF-adjudicated 2026-06-12).

Tolerance: the report's author states cross-machine floating-point differences
appear "usually in the 5th or 6th decimal place" (TR-91-6 p. iv), so the
source's own precision is ~1e-5 — used here as the assert tolerance (achieved
agreement is ~1e-6..1e-9 per case; see per-case comments).

DO-NOT-USE — three TR-91-6 cells are PDF-confirmed SOURCE-PRINT DEFECTS and
must never be wired as goldens (adjudication record:
docs/notes/2026-06-12-load-bearing-numbers-verification.md):

1. Stumpff row z = -39.47842 (p. E-12): printed C = 5.83559577,
   S = 0.97444596; power series AND closed cosh/sinh forms both give
   C = 6.75677528, S = 1.05406777, and no nearby z reproduces both printed
   values. The published cell is wrong.
2. Kepler case BMW App. D.3-3 hyperbolic propagated state (p. E-14): printed
   R_y = -0.1172043; repo ``propagate``, stepped propagation, and an
   independent RK4 all give R_y = -0.11822049 (delta ~1e-3, far above the
   report's own ~1e-5 floor). The printed a/e reproduce, so only the
   propagated state is defective.
3. Interplanetary Hohmann table Mars row (p. E-36): printed Mars distance
   227.8e6 km is inconsistent with its own printed outputs, which match
   277.8e6 km exactly (227.8 <-> 277.8 digit transposition). Not wired here
   (out of scope for lambert/kepler/stumpff) — recorded so it is never wired.
"""

from __future__ import annotations

import numpy as np
import pytest

from cyclerfinder.core._stumpff import stumpff_c, stumpff_s
from cyclerfinder.core.kepler import propagate
from cyclerfinder.core.lambert import lambert

# The report's stated cross-machine precision floor (TR-91-6 p. iv).
_TOL = 1.0e-5

# ---------------------------------------------------------------------------
# Stumpff C(z) and S(z) — FindCandS table, TR-91-6 p. E-12 (ref BMW p. 209)
# ---------------------------------------------------------------------------

# (z, C(z), S(z)) exactly as printed (8 decimal places). The fifth printed row
# z = -39.47842 is a SOURCE-PRINT DEFECT (see module docstring) and is
# deliberately absent. Reproduction (#203): all four rows agree to <= 4e-9.
_STUMPFF_ROWS = [
    (0.0, 0.50000000, 0.16666667),
    (0.57483, 0.47650300, 0.16194146),
    (39.47842, 0.00000000, 0.02533029),
    (50.0, 0.00589304, 0.01799504),
]


@pytest.mark.parametrize(("z", "c_expected", "s_expected"), _STUMPFF_ROWS)
def test_stumpff_tr916_table(z: float, c_expected: float, s_expected: float) -> None:
    """C(z) and S(z) match the printed table to its 8-decimal precision."""
    # 1e-8 = one unit in the last printed decimal place.
    assert abs(stumpff_c(z) - c_expected) < 1.0e-8
    assert abs(stumpff_s(z) - s_expected) < 1.0e-8


# ---------------------------------------------------------------------------
# Lambert (Gauss problem) — GAUSS test cases, TR-91-6 pp. E-20 ... E-27
# ---------------------------------------------------------------------------

# Each entry: (id, r1 DU, r2 DU, dt TU, prograde flag, printed v1, printed v2).
# The report labels cases "short way"/"long way"; the repo solver's
# short/long selection is geometry-driven via the ``prograde`` flag — each
# case below matched exactly one flag to < 1e-6 (#203). Velocities DU/TU,
# printed to 7 decimals. Reproduction (#203 + re-run before wiring): worst
# per-component |dv| <= 9.4e-7 across all six cases.
_LAMBERT_CASES = [
    # BMW App. D.4-1, long way (p. E-20): worst |dv| = 9.3e-8.
    (
        "D.4-1_long",
        (0.5, 0.6, 0.7),
        (0.0, -1.0, 0.0),
        20.0,
        True,
        (-0.1229814, 1.1921622, -0.1721740),
        (0.6698699, 0.4804848, 0.9378179),
    ),
    # BMW App. D.4-2, short way (p. E-20): worst |dv| = 9.4e-7.
    (
        "D.4-2_short",
        (0.3, 0.7, 0.4),
        (0.6, -1.4, 0.8),
        5.0,
        False,
        (0.7326124, -0.1048188, 0.9768165),
        (-0.3438450, -0.1048188, -0.4584600),
    ),
    # BMW App. D.4-6, short way, vectors almost 180 deg apart (p. E-22):
    # worst |dv| = 1.0e-7.
    (
        "D.4-6_short",
        (-0.4, 0.6, -1.201),
        (0.2, -0.3, 0.6),
        5.0,
        True,
        (0.2551050, -0.3826576, -0.5738817),
        (-0.7292157, 1.0938236, 0.4920219),
    ),
    # BMW p. 275 prob 5.11b, short way, planar (p. E-24): worst |dv| = 1.2e-7.
    (
        "5.11b_short",
        (1.2, 0.0, 0.0),
        (0.0, 2.0, 0.0),
        10.0,
        True,
        (0.7497686, 0.7090867, 0.0),
        (-0.4254520, -0.4661339, 0.0),
    ),
    # BMW p. 275 prob 5.11e, long way, planar (p. E-25): worst |dv| = 1.0e-7.
    (
        "5.11e_long",
        (2.0, 0.0, 0.0),
        (-2.0, -0.2, 0.0),
        20.0,
        True,
        (0.3083363, 0.7157383, 0.0),
        (0.3778475, -0.6779535, 0.0),
    ),
    # A423 test case, short way, 35 TU (p. E-27): worst |dv| = 4.7e-8.
    (
        "A423_35TU_short",
        (1.05, 0.0, 0.0),
        (0.0, 0.9, 0.0),
        35.0,
        True,
        (1.1418714, 0.5381541, 0.0),
        (-0.6278465, -1.2315637, 0.0),
    ),
]


@pytest.mark.parametrize(
    ("r1", "r2", "dt", "prograde", "v1_expected", "v2_expected"),
    [case[1:] for case in _LAMBERT_CASES],
    ids=[case[0] for case in _LAMBERT_CASES],
)
def test_lambert_tr916_golden(
    r1: tuple[float, float, float],
    r2: tuple[float, float, float],
    dt: float,
    prograde: bool,
    v1_expected: tuple[float, float, float],
    v2_expected: tuple[float, float, float],
) -> None:
    """Single-rev Lambert (mu = 1) reproduces the printed terminal velocities."""
    sols = lambert(
        np.array(r1, dtype=np.float64),
        np.array(r2, dtype=np.float64),
        dt,
        mu=1.0,
        prograde=prograde,
    )
    assert len(sols) == 1
    sol = sols[0]
    assert float(np.max(np.abs(sol.v1 - np.array(v1_expected)))) < _TOL
    assert float(np.max(np.abs(sol.v2 - np.array(v2_expected)))) < _TOL


# ---------------------------------------------------------------------------
# Kepler propagation — KEPLER test cases, TR-91-6 pp. E-13 ... E-17
# ---------------------------------------------------------------------------

# Each entry: (id, Ro DU, Vo DU/TU, Dt TU, printed R, printed V, rel tol).
# BMW App. D.3-3 (p. E-14) is deliberately absent — SOURCE-PRINT DEFECT in the
# propagated R_y (see module docstring). Reproduction (#203 + re-run before
# wiring): relative errors per case noted below.
_KEPLER_CASES = [
    # BMW example p. 210, ellipse (TR p. E-13): dr 7.5e-7, dv 1.9e-6 rel.
    (
        "BMW_p210_ellipse",
        (1.0, 0.0, 0.0),
        (0.0, 0.0, 1.1),
        2.0,
        (-0.3206670, 0.0, 1.2364349),
        (-0.8799766, -0.0, -0.0373113),
        _TOL,
    ),
    # BMW App. D.3-4, backward multi-rev, Dt = -20 TU (p. E-15): the report's
    # "elliptical orbit with multi-revs and backwards propagation" case.
    # dr 2.0e-8, dv 1.6e-7 rel.
    (
        "D.3-4_backward_multirev",
        (0.5, 0.7, 0.8),
        (0.0, 0.1, 0.9),
        -20.0,
        (0.0401556, 0.2664818, 1.9566242),
        (-0.2291452, -0.2755040, 0.0410620),
        _TOL,
    ),
    # Kaplan example p. 307, hyperbolic first guess (TR p. E-16):
    # dr 2.4e-7, dv 4.3e-7 rel.
    (
        "Kaplan_p307_hyperbolic",
        (1.5679, 0.0, 0.0),
        (0.0, 1.1638, 0.0),
        13.386,
        (-4.8259941, 7.3013686, -0.0),
        (-0.4571857, 0.3135849, -0.0),
        _TOL,
    ),
    # BMW problem p. 225 #4.18, radial escape, 219.6 TU (TR p. E-17): the
    # near-rectilinear escape geometry reproduces only loosely (dr 2.4e-5,
    # dv 3.7e-5 rel — #203 verdict "PASS (loose; radial escape)"), so this
    # case carries a documented 1e-4 relative tolerance instead of the 1e-5
    # report floor.
    (
        "p225_4.18_radial_escape",
        (0.2, 0.0, 0.0),
        (3.1622770, 0.0, 0.0),
        219.6,
        (60.1009021, 0.0, 0.0),
        (0.1824184, 0.0, 0.0),
        1.0e-4,
    ),
]


@pytest.mark.parametrize(
    ("r0", "v0", "dt", "r_expected", "v_expected", "rel_tol"),
    [case[1:] for case in _KEPLER_CASES],
    ids=[case[0] for case in _KEPLER_CASES],
)
def test_kepler_tr916_golden(
    r0: tuple[float, float, float],
    v0: tuple[float, float, float],
    dt: float,
    r_expected: tuple[float, float, float],
    v_expected: tuple[float, float, float],
    rel_tol: float,
) -> None:
    """Universal-variable propagation (mu = 1) reproduces the printed states."""
    r, v = propagate(
        np.array(r0, dtype=np.float64),
        np.array(v0, dtype=np.float64),
        dt,
        mu=1.0,
    )
    r_exp = np.array(r_expected)
    v_exp = np.array(v_expected)
    assert float(np.linalg.norm(r - r_exp) / np.linalg.norm(r_exp)) < rel_tol
    assert float(np.linalg.norm(v - v_exp) / np.linalg.norm(v_exp)) < rel_tol
