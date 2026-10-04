"""#896: printed Earth-Moon L1 and L2 numbers of the circular restricted three-body problem.

Part 1, L1 linear frequencies. Jorba, A., Jorba-Cusco, M. and Rosales, J. J. (2020), "The vicinity
of the Earth-Moon L1 point in the bicircular problem", Celestial Mechanics and Dynamical Astronomy
132:11, DOI 10.1007/s10569-019-9940-2, filed in the private paper corpus as
jorba-jorba-cusco-rosales-2020-vicinity-earth-moon-l1-point-bicircular-problem-cmda-132-11-doi-10.1007-s10569-019-9940-2.pdf.
Section 3.2, page 13, compares the Floquet frequencies of the bicircular L1 replacement orbit
(tested in ``test_bcr4bp_sun_sense.py``) with "the frequencies related to the equilibrium point L1
in the RTBP (2.33438585628816 and 2.2688310655411, respectively)": the planar (in-plane elliptic)
and the vertical frequency of the linearisation at L1, in units where the Earth-Moon mean motion
is 1. Table 1 (page 3) gives the mass ratio as mu = 0.012150581623433623.

Measured with the linearisation of ``core.cr3bp`` (``cr3bp_stm_eom`` with the identity as STM) at
an L1 found to machine precision:

- at the Table 1 mu the planar and vertical frequencies are 2.3044e-9 and 2.3545e-9 BELOW the
  printed values, about 1e-9 relative. This is the "about 1e-8" agreement of the digest.
- The two misses are the same mass-ratio error. Both frequencies grow with mu (d omega/d mu = 7.802
  planar, 7.974 vertical), and both misses are removed by the same change of mu, +2.9527e-10 (the
  two estimates agree to 4e-16 in mu). At mu = 1/(1 + 81.300585) = 0.012150581918706896 the model
  gives 2.33438585628800 and 2.26883106554112: residuals -1.6e-13 (planar, in double precision;
  1.6e-14 in 40-digit arithmetic) and +2e-14 (vertical, within the last printed digit).
- Table 1's double is itself exactly 1/(1 + 81.300587) (equal to the last bit). INFERRED: the
  printed RTBP frequencies were computed with an Earth/Moon mass ratio of 81.300585, two units in
  the sixth decimal from the 81.300587 behind Table 1. The paper does not say this; it is the
  single value that reproduces both printed frequencies to their printed precision. It is not a
  rounding of mu, a loose L1 root, or a different quantity: the quantities are the linearisation
  frequencies, and they match to 1e-13 once the mass ratio is changed.

Asserted: both frequencies within 1e-12 at mu = 1/(1 + 81.300585); at the Table 1 mu, both
misses within 2e-11 of the miss predicted from d omega/d mu and the mu difference; control: the
project's default Earth-Moon mu (0.01215058439469525, from ``cr3bp_system("Earth", "Moon")``)
misses both by about +1.9e-8, more than eight times the Table 1 miss.

Part 2, L2 bifurcation periods. Singh, J., Park, B. and Howell, K. C. (2026), "Evolution of L2
orbit families and bifurcations within intermediary-fidelity models", AAS/AIAA paper AAS 26-654,
filed in the private paper corpus as
singh-park-howell-2026-evolution-l2-orbit-families-bifurcations-intermediary-fidelity-models-qbcp-er3bp-AAS-26-654.pdf.
The CR3BP rows of Table 2 (p12, L2 Lyapunov family at the halo bifurcation, out-of-plane pair
"Center -> Saddle", 14.8319 days), Table 3 (p12, L2 Lyapunov family at the axial bifurcation,
"Saddle -> Center", 18.7183 days) and Table 4 (p14, L2 vertical family at the axial bifurcation,
"Center -> Saddle", 19.2033 days). The paper prints neither mu nor its length and time units. The
project's Earth-Moon system (``cr3bp_system("Earth", "Moon")``: mu = 0.01215058439469525,
384,400 km, time unit 4.342480 days) is used; the digest infers the same nominal units from the
tables' p:q brackets, and the test below pins it from those brackets. Measured on the vertical
family, the bifurcation period moves by 2.35e-6 day per 1e-7 change of mu, so the unknown mu does
not matter at four decimals (closing the vertical miss below would need mu changed by about 4e-6).

Method, written here on ``core.cr3bp`` (``cr3bp_eom``, ``cr3bp_stm_eom``, ``propagate``):
single-shooting perpendicular-crossing families with natural continuation, then a root of the
stability function of the bifurcating pair.
- Lyapunov: (x0, 0, 0, 0, vy0, 0) on the far side of L2, half period to the next y = 0 crossing,
  vx = 0 there; continued in x0 from the linear seed. For a planar orbit the (z, vz) block of the
  monodromy matrix decouples, and with B the (z, vz) block of the half-period STM its trace minus 2
  is 4 B[z, vz] B[vz, z] (negative: centre; positive: saddle). Checked against the full-period
  monodromy below.
- Vertical: (x0, 0, 0, 0, vy0, vz0) at the node of the figure eight, quarter period to the next
  y = 0 crossing, vx = vz = 0 there; continued in vz0 from a seed at vz0 = 0.2 (the seed is only a
  start for Newton). The pair is tracked by trace invariants of the full-period monodromy M:
  s1 + s2 = tr M - 2, s1^2 + s2^2 = tr M^2 + 2 (s = lambda + 1/lambda per pair), so that the
  small-|s| pair is found without picking eigenvalues near +1, where the bifurcating pair meets
  the trivial pair.

Measured: halo bifurcation 14.831874 d (printed 14.8319, -2.6e-5 d); Lyapunov axial bifurcation
18.718299 d (printed 18.7183, -7e-7 d); vertical axial bifurcation 19.203197 d (printed 19.2033,
-1.03e-4 d). Each stability change is in the printed direction. The vertical value does not move
with the integrator: 19.2031974 d at rtol = atol = 1e-12, at rtol 1e-13 / atol 1e-14, and with
``stm_mode="fixed_path"`` for the monodromy. The stability function is checked once against the
eigenvalues of the monodromy at a member far from the root (vz0 = 0.4, s - 2 = -0.3707).

The time unit is pinned by the paper itself: each of the six ER3BP bracket ends of Tables 2 to 4
is a p:q orbit of period q P_sys / p with P_sys = 2 pi time units, so P_sys = p P / q, with
(p / q) 5e-5 day of rounding. The six intervals intersect in [27.284585, 27.284615] day; the
project's 2 pi t_s = 27.284606 day lies inside, which pins the paper's time unit to within -7.8e-7
and +3.2e-7 relative of the project's. Matching 19.2033 would need a unit at least 2.7e-6 longer,
which those brackets exclude. So the vertical miss survives both the rounding and the unknown units.
It is asserted at the printed precision (half a unit in the fourth decimal, 5e-5 day) and kept as a
strict xfail: a finding, about one unit in the last printed digit. Cause not determined (a coarse
bifurcation estimate in the paper's continuation, which samples discrete family members, would be
one explanation; INFERRED, not stated in the paper).
Control: a time unit from the sidereal month (27.321661 d / 2 pi = 4.348371 d) moves the three
periods by 0.020 to 0.026 day, far outside the printed digits.
"""

from __future__ import annotations

import dataclasses

import numpy as np
import pytest
from numpy.typing import NDArray
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from cyclerfinder.core.cr3bp import (
    CR3BPSystem,
    cr3bp_eom,
    cr3bp_stm_eom,
    cr3bp_system,
    propagate,
)

FloatArray = NDArray[np.float64]

# Jorba, Jorba-Cusco & Rosales 2020, section 3.2, p13: RTBP L1 planar and vertical frequencies.
_JORBA_L1_PLANAR = 2.33438585628816
_JORBA_L1_VERTICAL = 2.2688310655411
# Table 1, p3.
_JORBA_TABLE1_MU = 0.012150581623433623
# INFERRED (module docstring): the mass ratio behind the printed frequencies.
_MU_FROM_FREQUENCIES = 1.0 / (1.0 + 81.300585)


def _l1_x(mu: float) -> float:
    def ax(x: float) -> float:
        return float(cr3bp_eom(0.0, np.array([x, 0.0, 0.0, 0.0, 0.0, 0.0]), mu)[3])

    return float(brentq(ax, 0.5, 1.0 - mu - 1e-3, xtol=1e-15, rtol=8.9e-16))


def _l1_frequencies(mu: float) -> tuple[float, float]:
    """Planar and vertical frequency of the linearisation of ``core.cr3bp`` at L1."""
    y42 = np.concatenate([[_l1_x(mu), 0.0, 0.0, 0.0, 0.0, 0.0], np.eye(6).ravel()])
    jac = cr3bp_stm_eom(0.0, y42, mu)[6:].reshape(6, 6)
    eigs = np.linalg.eigvals(jac)
    vertical = float(np.sqrt(-jac[5, 2]))
    imaginary = [abs(e.imag) for e in eigs if abs(e.real) < 1e-9 and abs(e.imag) > 0.0]
    planar = [w for w in imaginary if abs(w - vertical) > 1e-3]
    assert len(imaginary) == 4 and len(planar) == 2  # saddle x centre x centre
    return float(planar[0]), vertical


def test_jorba_2020_l1_frequencies_at_the_mass_ratio_81_300585() -> None:
    planar, vertical = _l1_frequencies(_MU_FROM_FREQUENCIES)
    assert abs(planar - _JORBA_L1_PLANAR) < 1e-12
    assert abs(vertical - _JORBA_L1_VERTICAL) < 1e-12


def test_jorba_2020_table1_mu_misses_both_frequencies_by_one_mass_ratio_error() -> None:
    assert _JORBA_TABLE1_MU == 1.0 / (1.0 + 81.300587)
    planar, vertical = _l1_frequencies(_JORBA_TABLE1_MU)
    miss_p, miss_v = planar - _JORBA_L1_PLANAR, vertical - _JORBA_L1_VERTICAL
    assert -2.4e-9 < miss_p < -2.2e-9 and -2.45e-9 < miss_v < -2.25e-9
    h = 1e-7
    hi = _l1_frequencies(_JORBA_TABLE1_MU + h)
    lo = _l1_frequencies(_JORBA_TABLE1_MU - h)
    slope_p, slope_v = (hi[0] - lo[0]) / (2 * h), (hi[1] - lo[1]) / (2 * h)
    d_mu = _JORBA_TABLE1_MU - _MU_FROM_FREQUENCIES
    assert abs(miss_p - slope_p * d_mu) < 2e-11
    assert abs(miss_v - slope_v * d_mu) < 2e-11


def test_project_default_earth_moon_mu_is_discriminated() -> None:
    mu = cr3bp_system("Earth", "Moon").mu
    planar, vertical = _l1_frequencies(mu)
    assert planar - _JORBA_L1_PLANAR > 1.5e-8
    assert vertical - _JORBA_L1_VERTICAL > 1.5e-8


# ---------------------------------------------------------------------------------------------
# Part 2: Singh, Park & Howell 2026, CR3BP rows of Tables 2 to 4.

_SINGH_HALO_BIFURCATION_D = 14.8319  # Table 2, p12, L2 Lyapunov, centre to saddle
_SINGH_LYAPUNOV_AXIAL_BIFURCATION_D = 18.7183  # Table 3, p12, L2 Lyapunov, saddle to centre
_SINGH_VERTICAL_AXIAL_BIFURCATION_D = 19.2033  # Table 4, p14, L2 vertical, centre to saddle
_HALF_UNIT_4TH_DECIMAL_D = 5e-5
_SIDEREAL_MONTH_D = 27.321661
# ER3BP p:q bracket ends of Tables 2 to 4 (period in days, p, q), Singh et al. 2026 p12, p14.
_SINGH_ER3BP_BRACKET_ENDS = (
    (14.7792, 24, 13),
    (14.8825, 11, 6),
    (18.7007, 89, 61),
    (18.7094, 35, 24),
    (19.1730, 37, 26),
    (19.3266, 24, 17),
)

_EM = cr3bp_system("Earth", "Moon")
_DAY_PER_UNIT = _EM.t_s / 86400.0


def _to_crossing(state6: FloatArray, sense: float, mu: float) -> tuple[float, FloatArray]:
    """Time and state+STM (42) at the first y = 0 crossing after t = 0 (y moving with ``sense``)."""

    def rhs(t: float, y: FloatArray) -> FloatArray:
        return cr3bp_stm_eom(t, y, mu)

    def event(t: float, y: FloatArray) -> float:
        return float(y[1])

    event.terminal = True  # type: ignore[attr-defined]
    event.direction = -sense  # type: ignore[attr-defined]
    y0 = np.concatenate([state6, np.eye(6).ravel()])
    sol = solve_ivp(
        rhs,
        (0.0, 10.0),
        y0,
        method="DOP853",
        rtol=1e-12,
        atol=1e-12,
        events=event,
    )
    assert sol.t_events is not None and sol.y_events is not None
    assert sol.t_events[0].size == 1
    return float(sol.t_events[0][0]), np.asarray(sol.y_events[0][0], dtype=np.float64)


def _crossing_jacobian(
    yf: FloatArray, cols: tuple[int, ...], rows: tuple[int, ...], mu: float
) -> FloatArray:
    """d(final state rows)/d(initial state cols) at the y = 0 crossing, time free."""
    phi, s = yf[6:].reshape(6, 6), yf[:6]
    f = cr3bp_eom(0.0, s, mu)
    jac = np.empty((len(rows), len(cols)))
    for j, c in enumerate(cols):
        col = phi[:, c] - f * phi[1, c] / s[4]
        jac[:, j] = col[list(rows)]
    return jac


def _lyapunov(x0: float, vy0: float, mu: float) -> tuple[float, float, FloatArray]:
    """L2 Lyapunov orbit through (x0, 0) with vx = 0: (vy0, period, half-period state+STM)."""
    for _ in range(25):
        t, yf = _to_crossing(np.array([x0, 0.0, 0.0, 0.0, vy0, 0.0]), np.sign(vy0), mu)
        if abs(yf[3]) < 1e-12:
            return vy0, 2.0 * t, yf
        vy0 -= yf[3] / _crossing_jacobian(yf, (4,), (3,), mu)[0, 0]
    raise RuntimeError(f"Lyapunov corrector did not converge at x0 = {x0}")


def _lyapunov_out_of_plane(yf: FloatArray) -> float:
    """Full-period (z, vz) monodromy trace minus 2, from the half-period STM."""
    phi = yf[6:].reshape(6, 6)
    return 4.0 * float(phi[2, 5] * phi[5, 2])


def _vertical(vz0: float, x0: float, vy0: float, mu: float) -> tuple[float, float, float]:
    """L2 vertical orbit through (x0, 0, 0) with velocity (0, vy0, vz0): (x0, vy0, period)."""
    for _ in range(30):
        t, yf = _to_crossing(np.array([x0, 0.0, 0.0, 0.0, vy0, vz0]), np.sign(vy0), mu)
        res = np.array([yf[3], yf[5]])
        if float(np.max(np.abs(res))) < 1e-11:
            return x0, vy0, 4.0 * t
        dx = np.linalg.solve(_crossing_jacobian(yf, (0, 4), (3, 5), mu), -res)
        x0, vy0 = x0 + float(dx[0]), vy0 + float(dx[1])
    raise RuntimeError(f"vertical corrector did not converge at vz0 = {vz0}")


def _small_pair_s_minus_2(system: CR3BPSystem, state6: FloatArray, period: float) -> float:
    """s - 2 for the non-trivial pair with the smaller |s| (s = lambda + 1/lambda)."""
    stm = propagate(system, state6, period, with_stm=True).stm
    assert stm is not None
    p = float(np.trace(stm)) - 2.0
    sum_sq = float(np.trace(stm @ stm)) + 2.0
    q = (p * p - sum_sq) / 2.0
    s_large = (p + np.sign(p) * np.sqrt(max(p * p - 4.0 * q, 0.0))) / 2.0
    return float(q / s_large) - 2.0


@dataclasses.dataclass(frozen=True)
class _Bifurcation:
    period_units: float
    before: float  # stability function a little below the root (in the family's period order)
    after: float


@pytest.fixture(scope="module")
def l2_bifurcations() -> dict[str, _Bifurcation]:
    mu = _EM.mu

    def ax(x: float) -> float:
        return float(cr3bp_eom(0.0, np.array([x, 0.0, 0.0, 0.0, 0.0, 0.0]), mu)[3])

    x_l2 = float(brentq(ax, 1.0, 1.3, xtol=1e-15))
    r1, r2 = x_l2 + mu, x_l2 - 1.0 + mu
    c2 = (1.0 - mu) / r1**3 + mu / r2**3
    wp = np.sqrt(
        ((2.0 - c2) + np.sqrt((2.0 - c2) ** 2 + 4.0 * (c2 - 1.0) * (1.0 + 2.0 * c2))) / 2.0
    )
    kappa = (wp * wp + 1.0 + 2.0 * c2) / (2.0 * wp)
    out: dict[str, _Bifurcation] = {}

    # Lyapunov family, continued outwards in x0 from the linear seed.
    xs: list[float] = []
    vys: list[float] = []
    gs: list[float] = []
    x0, vy0 = x_l2 + 0.003, float(-kappa * 0.003 * wp)
    while True:
        vy0, period, yf = _lyapunov(x0, vy0, mu)
        xs.append(x0)
        vys.append(vy0)
        gs.append(_lyapunov_out_of_plane(yf))
        if period * _DAY_PER_UNIT > 19.0:
            break
        vy0 = 2 * vys[-1] - vys[-2] if len(vys) > 1 else vy0
        x0 += 0.002

    def g_lyap(x: float) -> float:
        return _lyapunov_out_of_plane(_lyapunov(x, float(np.interp(x, xs, vys)), mu)[2])

    roots = [
        float(brentq(g_lyap, xs[i], xs[i + 1], xtol=1e-13))
        for i in range(len(gs) - 1)
        if gs[i] * gs[i + 1] < 0.0
    ]
    assert len(roots) == 2  # halo, then axial
    for name, xr in zip(("halo", "lyapunov_axial"), roots, strict=True):
        vy, period, yf = _lyapunov(xr, float(np.interp(xr, xs, vys)), mu)
        out[name] = _Bifurcation(period, g_lyap(xr - 1e-4), g_lyap(xr + 1e-4))
        if name == "halo":
            # Instrument check: the half-period formula against the full-period monodromy.
            stm = propagate(_EM, np.array([xr, 0.0, 0.0, 0.0, vy, 0.0]), period, with_stm=True).stm
            assert stm is not None
            assert abs(stm[2, 2] + stm[5, 5] - 2.0) < 1e-9

    # Vertical family, continued in vz0 from a seed at vz0 = 0.2 (period about 15.58 d).
    vzs: list[float] = []
    xvs: list[float] = []
    vyvs: list[float] = []
    hs: list[float] = []
    vz0, x0, vy0 = 0.2, 1.1453, -0.0253
    while True:
        x0, vy0, period = _vertical(vz0, x0, vy0, mu)
        vzs.append(vz0)
        xvs.append(x0)
        vyvs.append(vy0)
        state = np.array([x0, 0.0, 0.0, 0.0, vy0, vz0])
        hs.append(_small_pair_s_minus_2(_EM, state, period))
        if len(vzs) == 41:  # vz0 = 0.4, far from the root: check against the eigenvalues
            stm = propagate(_EM, state, period, with_stm=True).stm
            assert stm is not None
            s_vals = [float((e + 1.0 / e).real) for e in np.linalg.eigvals(stm)]
            s_small = min(s_vals, key=lambda v: abs(v - 2.0) if abs(v - 2.0) > 1e-3 else 1e9)
            assert abs(hs[-1] - (s_small - 2.0)) < 1e-7
        if hs[-1] > 0.0:
            break
        if len(vzs) > 1:
            x0, vy0 = 2 * xvs[-1] - xvs[-2], 2 * vyvs[-1] - vyvs[-2]
        vz0 += 0.005

    def vert(vz: float) -> tuple[FloatArray, float]:
        x, vy, period = _vertical(
            vz, float(np.interp(vz, vzs, xvs)), float(np.interp(vz, vzs, vyvs)), mu
        )
        return np.array([x, 0.0, 0.0, 0.0, vy, vz]), period

    def h_vert(vz: float) -> float:
        return _small_pair_s_minus_2(_EM, *vert(vz))

    vr = float(brentq(h_vert, vzs[-2], vzs[-1], xtol=1e-13))
    out["vertical_axial"] = _Bifurcation(vert(vr)[1], h_vert(vr - 1e-4), h_vert(vr + 1e-4))
    return out


_PERIOD_CASES = [
    pytest.param("halo", _SINGH_HALO_BIFURCATION_D, id="table2_halo"),
    pytest.param("lyapunov_axial", _SINGH_LYAPUNOV_AXIAL_BIFURCATION_D, id="table3_lyapunov_axial"),
    pytest.param(
        "vertical_axial",
        _SINGH_VERTICAL_AXIAL_BIFURCATION_D,
        id="table4_vertical_axial",
        marks=pytest.mark.xfail(
            strict=True,
            reason="#896: core.cr3bp puts the L2 vertical-family axial bifurcation at 19.203197 d, "
            "1.03e-4 d below the printed 19.2033 (two half-units of the fourth decimal); not "
            "integrator tolerance or mu, and the paper's own ER3BP brackets pin its time unit to "
            "1e-6, excluding the 2.7e-6 longer unit a match would need",
            raises=AssertionError,
        ),
    ),
]


@pytest.mark.parametrize(("name", "printed_days"), _PERIOD_CASES)
def test_singh_2026_cr3bp_l2_bifurcation_period(
    l2_bifurcations: dict[str, _Bifurcation], name: str, printed_days: float
) -> None:
    period_d = l2_bifurcations[name].period_units * _DAY_PER_UNIT
    assert abs(period_d - printed_days) < _HALF_UNIT_4TH_DECIMAL_D


@pytest.mark.parametrize(
    ("name", "change"),
    [
        ("halo", "centre_to_saddle"),  # Table 2
        ("lyapunov_axial", "saddle_to_centre"),  # Table 3
        ("vertical_axial", "centre_to_saddle"),  # Table 4
    ],
)
def test_singh_2026_stability_change_is_in_the_printed_direction(
    l2_bifurcations: dict[str, _Bifurcation], name: str, change: str
) -> None:
    """Family ordered by increasing period; negative: centre, positive: saddle."""
    b = l2_bifurcations[name]
    if change == "centre_to_saddle":
        assert b.before < 0.0 < b.after
    else:
        assert b.before > 0.0 > b.after


def test_singh_2026_er3bp_brackets_pin_the_time_unit() -> None:
    """The paper's ER3BP p:q bracket ends give P_sys = 2 pi time units in days, within rounding."""
    lo = max(p * (period - 5e-5) / q for period, p, q in _SINGH_ER3BP_BRACKET_ENDS)
    hi = min(p * (period + 5e-5) / q for period, p, q in _SINGH_ER3BP_BRACKET_ENDS)
    assert lo < 2.0 * np.pi * _DAY_PER_UNIT < hi
    # A unit long enough to put the vertical bifurcation at 19.2033 is outside the window.
    assert 2.0 * np.pi * _DAY_PER_UNIT * (1.0 + 2.7e-6) > hi


def test_singh_2026_sidereal_month_time_unit_is_discriminated(
    l2_bifurcations: dict[str, _Bifurcation],
) -> None:
    sidereal_unit_d = _SIDEREAL_MONTH_D / (2.0 * np.pi)
    for name, printed in (
        ("halo", _SINGH_HALO_BIFURCATION_D),
        ("lyapunov_axial", _SINGH_LYAPUNOV_AXIAL_BIFURCATION_D),
        ("vertical_axial", _SINGH_VERTICAL_AXIAL_BIFURCATION_D),
    ):
        assert abs(l2_bifurcations[name].period_units * sidereal_unit_d - printed) > 0.015
