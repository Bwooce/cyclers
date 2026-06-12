"""INDEPENDENT-ORACLE CROSS-CHECKS against Vallado's published code release.

These are NOT print-goldens (see test_vallado_goldens.py for those). Every
EXPECTED value below is the output of a SECOND PUBLISHED IMPLEMENTATION:

    Vallado, Fundamentals of Astrodynamics and Applications, 5th ed. code
    release (github.com/CelesTrak/fundamentals-of-astrodynamics),
    software/python/src/valladopy — astro/twobody/utils.py (findc2c3),
    astro/twobody/kepler.py (kepler), astro/iod/lambert.py (universal),
    retrieved 2026-06-13.

The oracle values were computed once from that release (AGPL-3.0; no code
from it is copied here — only its numeric outputs, cited as facts) and are
hard-coded as literals so this suite does NOT depend on the reference clone.
Generating inputs are recorded per case. valladopy works in km/km/s/s with
Earth MU = 398600.4415 km^3/s^2 hard-wired; canonical-unit cases (mu = 1)
were scaled through DU = RE = 6378.1363 km, TU = sqrt(DU^3/MU), VU = DU/TU
before/after calling the oracle. Stumpff C/S correspond to valladopy c2/c3.

TR-91-6 DEFECT-CELL EXTERNAL ADJUDICATION (#220, run 2026-06-13) — the three
source-print defects recorded in test_vallado_goldens.py were re-run through
Vallado's own current code:

1. Stumpff z = -39.47842: findc2c3 gives C = 6.756775284482481,
   S = 1.0540677718375635 — agrees with OUR adjudicated values to 0 ulp and
   REFUTES the 1991 print (5.83559577 / 0.97444596). CONFIRMS #203.
2. Kepler BMW App. D.3-3 hyperbolic (Ro=(0.3,1,0), Vo=(3,0,0), Dt=5, mu=1):
   their kepler gives R = (13.962281215332403, -0.11822048981640433, 0) —
   R_y agrees with OUR -0.11822049 to 1.8e-15 and REFUTES the 1991 print
   (-0.1172043). CONFIRMS #203. (Wired below as an oracle-backed case, which
   rescues the defect cell's inputs for regression coverage.)
3. Mars-Hohmann 227.8 <-> 277.8 transposition (TR-91-6 p. E-36): the 5th-ed
   code release contains NO interplanetary-Hohmann example (matlab examples
   stop at ch11; no ch12 scripts; no Mars heliocentric-distance constant in
   valladopy) and none of the six errata PDFs mentions 227.8/277.8 —
   EXTERNALLY UNADJUDICABLE from this release; the #203 internal-consistency
   adjudication (printed outputs match 277.8e6 km exactly) stands.

Implementation note observed while generating oracle values: valladopy's
findc2c3 switches to the z ~ 0 limit only for |z| <= 1e-10, so at z = 1e-9
its trig branch suffers catastrophic cancellation (returns 0.50000004...
vs the series value 0.5 - z/24 ~ 0.49999999996). No near-zero oracle row is
wired for that reason; our _stumpff series branch is the better behaved one
there.

Observed agreement at the wired cases is machine precision (<= 5.2e-15 per
component); asserts use 1e-12 absolute-per-component (Lambert/Kepler) and
1e-12 relative (Stumpff) for headroom.
"""

from __future__ import annotations

import numpy as np
import pytest

from cyclerfinder.core._stumpff import stumpff_c, stumpff_s
from cyclerfinder.core.kepler import propagate
from cyclerfinder.core.lambert import lambert

# Oracle-agreement tolerance (observed <= 5.2e-15; see module docstring).
_TOL = 1.0e-12

# ---------------------------------------------------------------------------
# Stumpff C(z), S(z) vs valladopy findc2c3 (c2 = C, c3 = S)
# ---------------------------------------------------------------------------

# (z, oracle C, oracle S) — oracle: findc2c3(z), full repr precision.
# Our values matched the oracle to 0 ulp on all four rows.
_STUMPFF_ORACLE_ROWS = [
    # Deep hyperbolic, well past the TR-91-6 table's most-negative good row.
    (-120.0, 238.35883931017943, 21.751512890965863),
    # The TR-91-6 p. E-12 defect-cell z: oracle REFUTES the 1991 print and
    # confirms our adjudicated C/S (see module docstring, item 1).
    (-39.47842, 6.756775284482481, 1.0540677718375635),
    # Mildly hyperbolic, between table rows.
    (-2.0, 0.5890917783042855, 0.1841494360042954),
    # Large positive z, beyond the table's z = 50 (sqrt(z) ~ 15.8 rad).
    (250.0, 0.007978625588375858, 0.004026118072480247),
]


@pytest.mark.parametrize(("z", "c_oracle", "s_oracle"), _STUMPFF_ORACLE_ROWS)
def test_stumpff_vs_vallado_code(z: float, c_oracle: float, s_oracle: float) -> None:
    """C(z)/S(z) agree with Vallado's findc2c3 to ~machine precision."""
    assert abs(stumpff_c(z) - c_oracle) < _TOL * max(1.0, abs(c_oracle))
    assert abs(stumpff_s(z) - s_oracle) < _TOL * max(1.0, abs(s_oracle))


# ---------------------------------------------------------------------------
# Kepler propagation vs valladopy kepler (canonical mu = 1 via DU/TU scaling)
# ---------------------------------------------------------------------------

# (id, Ro DU, Vo DU/TU, Dt TU, oracle R DU, oracle V DU/TU)
# Oracle call: kepler(ro*DU, vo*VU, dt*TU); outputs divided by DU / VU.
_KEPLER_ORACLE_CASES = [
    # BMW App. D.3-3 hyperbolic — the TR-91-6 p. E-14 defect-cell inputs,
    # now oracle-backed (module docstring item 2). Ours agreed to 1.8e-15.
    (
        "D.3-3_hyperbolic_oracle",
        (0.3, 1.0, 0.0),
        (3.0, 0.0, 0.0),
        5.0,
        (13.962281215332403, -0.11822048981640433, 0.0),
        (2.6779022951452025, -0.237538756730562, 0.0),
    ),
    # New regime vs the #187 golden set: inclined ellipse, multi-rev forward
    # propagation (~8 revs). Ours agreed to 5.7e-15.
    (
        "inclined_ellipse_multirev",
        (1.0, 0.05, 0.1),
        (0.1, 1.05, 0.2),
        50.0,
        (0.9012785519028897, 0.6427329705938947, 0.1987949539173345),
        (-0.3924298766887529, 0.8796086159292749, 0.12425339592452701),
    ),
    # Golden-input reuse (BMW p. 210 ellipse) at full oracle precision —
    # ties the print-golden case to the code oracle. Ours agreed to 4.9e-16.
    (
        "BMW_p210_ellipse_oracle",
        (1.0, 0.0, 0.0),
        (0.0, 0.0, 1.1),
        2.0,
        (-0.3206678684492108, 0.0, 1.236434486125314),
        (-0.8799780238144452, 0.0, -0.03731220212764113),
    ),
]


@pytest.mark.parametrize(
    ("r0", "v0", "dt", "r_oracle", "v_oracle"),
    [case[1:] for case in _KEPLER_ORACLE_CASES],
    ids=[case[0] for case in _KEPLER_ORACLE_CASES],
)
def test_kepler_vs_vallado_code(
    r0: tuple[float, float, float],
    v0: tuple[float, float, float],
    dt: float,
    r_oracle: tuple[float, float, float],
    v_oracle: tuple[float, float, float],
) -> None:
    """Universal-variable propagation agrees with Vallado's kepler."""
    r, v = propagate(
        np.array(r0, dtype=np.float64),
        np.array(v0, dtype=np.float64),
        dt,
        mu=1.0,
    )
    assert float(np.max(np.abs(r - np.array(r_oracle)))) < _TOL
    assert float(np.max(np.abs(v - np.array(v_oracle)))) < _TOL


# ---------------------------------------------------------------------------
# Lambert vs valladopy universal (canonical mu = 1 via DU/TU scaling)
# ---------------------------------------------------------------------------

# (id, r1 DU, r2 DU, dt TU, prograde, dm, oracle v1, oracle v2)
# Oracle call: universal(r1*DU, [0,0,1], r2*DU, dt*TU, dm, LOW, nrev=0,
# kbi=0, tol=1e-8, n_iter=60); outputs divided by VU. `dm` records which
# DirectionOfMotion (SHORT/LONG) matched our geometry-driven prograde branch;
# the other dm gives the complementary arc (checked, O(1) different).
_LAMBERT_ORACLE_CASES = [
    # Golden-input reuse (BMW p. 275 prob 5.11b) at full oracle precision.
    # dm=SHORT. Ours agreed to 4.4e-16.
    (
        "5.11b_planar_short",
        (1.2, 0.0, 0.0),
        (0.0, 2.0, 0.0),
        10.0,
        True,
        (0.7497684879725941, 0.7090867634500305, 0.0),
        (-0.4254520580700183, -0.4661337825925817, 0.0),
    ),
    # Golden-input reuse (BMW App. D.4-1, long way, 3D). dm=LONG.
    # Ours agreed to 2.2e-15.
    (
        "D.4-1_3d_long",
        (0.5, 0.6, 0.7),
        (0.0, -1.0, 0.0),
        20.0,
        True,
        (-0.12298143871958415, 1.192162120874134, -0.17217401420741776),
        (0.6698699236688174, 0.48048470742678623, 0.9378178931363442),
    ),
    # New regime vs the #187 golden set: HYPERBOLIC fast transfer
    # (dt = 0.5 TU; v1^2/2 - 1/r1 > 0). dm=SHORT. Ours agreed to 2.2e-15.
    (
        "hyperbolic_fast_short",
        (1.0, 0.0, 0.0),
        (0.0, 1.2, 0.3),
        0.5,
        True,
        (-1.7475966249092854, 2.5530993758370037, 0.6382748439592509),
        (-2.1275828131975034, 2.184458625110377, 0.5461146562775943),
    ),
]


@pytest.mark.parametrize(
    ("r1", "r2", "dt", "prograde", "v1_oracle", "v2_oracle"),
    [case[1:] for case in _LAMBERT_ORACLE_CASES],
    ids=[case[0] for case in _LAMBERT_ORACLE_CASES],
)
def test_lambert_vs_vallado_code(
    r1: tuple[float, float, float],
    r2: tuple[float, float, float],
    dt: float,
    prograde: bool,
    v1_oracle: tuple[float, float, float],
    v2_oracle: tuple[float, float, float],
) -> None:
    """Single-rev Lambert agrees with Vallado's universal-variable solver."""
    sols = lambert(
        np.array(r1, dtype=np.float64),
        np.array(r2, dtype=np.float64),
        dt,
        mu=1.0,
        prograde=prograde,
    )
    assert len(sols) == 1
    sol = sols[0]
    assert float(np.max(np.abs(sol.v1 - np.array(v1_oracle)))) < _TOL
    assert float(np.max(np.abs(sol.v2 - np.array(v2_oracle)))) < _TOL
