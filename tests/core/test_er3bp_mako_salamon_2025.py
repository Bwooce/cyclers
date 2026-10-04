"""Published checks for ``core/er3bp.py``: the weak-stability sweeps over the Earth's true anomaly
in Mako & Salamon (2025).

Source: Z. Mako and J. Salamon, "A note on the weak stability transition region in the planar
elliptic restricted three-body problem", Celestial Mechanics and Dynamical Astronomy 137:34 (2025),
DOI 10.1007/s10569-025-10264-0. Filed in the private paper corpus as
mako-salamon-2025-weak-stability-transition-region-planar-elliptic-restricted-three-body-problem-
cmda-137-34-doi-10.1007-s10569-025-10264-0.pdf. The paper has no tables; the numbers used here are
printed in the text:

- Sun-Earth system, mu = 3.003158242e-6, e = 0.0167 (Section 5, p8); collision with the Earth when
  the distance reaches R_E = 6378 km (p8); speed grid step step_v2 = 0.0063 km/s (p8).
- Figure 7 sweep (p13): alpha = 45 deg, r2 = 381 829 km, v2 = 1.20446 km/s, direct. "If f0 in
  [0, 177 deg] ... weakly stable, but if f0 in [178, 215 deg] then the orbit will be weak unstable,
  and for f0 in [216, 360 deg) we get weak stable orbits again." With e = 0 the orbit is weakly
  stable for every f0.
- Figure 8 sweep (p15): alpha = 135 deg, r2 = 12 754 km, v2 = 7.63689 km/s, direct; at this distance
  v_c = 5.496891 km/s and v_e = 7.773778 km/s. Weakly stable for f0 in [0, 118 deg], weakly unstable
  for [119, 259 deg], weakly stable for [260, 360 deg); with e = 0 all orbits are weakly stable.
  Both points are said to be "a stable point of WSTR ... near to the first boundary curve".

The classifier is the paper's Definition 1 (p3-4) and Definition 3 (p5), implemented here on
``er3bp_eom``. P3 starts at the periapsis of an osculating ellipse about P2, at physical distance r2
on the half-line l(f0, alpha), the half-line from P2 at angle alpha to the P1P2 axis (measured from
the direction away from P1, counter-clockwise, Fig. 9). Its velocity relative to P2 in the fixed
frame P2xy is perpendicular to that half-line, of magnitude v2, counter-clockwise (direct,
p_theta > 0, A.5). The orbit is weakly stable if it makes a full cycle around P2 and returns to the
half-line (theta in P2xy, theta = theta-tilde + f by A4, has advanced by 2 pi, with theta-dot > 0)
with Kepler energy about P2 (A5) negative. It is weakly unstable if the energy at that return is
positive, if it collides (distance below R_E), or if it does not return. The paper does not
quantify "going near to P1" or give a time limit; here the integration stops at 0.2 AU from the
Earth or after one revolution of the primaries (delta f = 2 pi). Moving the cut to 0.05 AU, or the
limit to 4 pi, changed no classification in either integer-degree sweep.

Speed scale. The paper's km/s are the two-body speeds times (1 - e): the printed v_c and v_e at
12 754 km are sqrt(GM/r) and sqrt(2 GM/r) times 0.98327 (the digest's section 3.4; the paper does
not say why). The speeds are therefore used as ratios to the paper's own circular speed, v_c(r) =
5.496891 km/s * sqrt(12 754 km / r): v2/v_c = 1.198910 for Figure 7 and 1.389311 for Figure 8. Only
the ratio and r2/a (a = 1 au = 149 597 870.7 km) enter the dimensionless problem.

Results (measured with this file's classifier, DOP853, rtol 1e-11; rtol 1e-13 gives the same
classification over the whole Figure 7 sweep):

- Figure 7: weakly stable at every integer f0 except 327 to 331 deg, where the Kepler energy at the
  return is positive. The first unstable speed of the model, as a ratio to v_c, varies smoothly with
  f0 between 1.1989 (near f0 = 330 deg) and 1.2055 (near 150 deg): the printed point 1.198910 sits
  2e-6 above the model's lowest boundary speed and 5e-4 below its boundary at f0 = 0, which is the
  paper's own speed resolution (eps_v2 = step_v2/10, p10). So the model agrees that the point lies
  on the first boundary, but its least stable f0 is near 330 deg (the start just before perihelion,
  where the Hill radius is smallest), not 178 to 215 deg (aphelion). The printed unstable band is
  not reproduced and five points of the second stable band are unstable.
- Figure 8: weakly stable at every integer f0. At v2/v_c = 1.389 (osculating eccentricity 0.930) the
  two-body apoapsis is 3.5e5 km, a quarter of the Hill radius (1.47e6 km), and the Sun cannot unbind
  the orbit in one 9-day revolution. The model's first non-stable speed is about 0.998 v_e
  (collision after a perturbed passage), for every f0. The orbits drawn in Figure 8 reach several
  1e6 km, which is what the model gives near 0.998 v_e, not at the printed speed: the printed speed
  is probably not the speed the figure was computed with. The printed unstable band [119, 259 deg]
  is not reproduced.
- The variant "1.22426 km/s" for the Figure 7 speed (the digest's section 3.4, from the p15 Moon
  remark) gives v2/v_c = 1.2186, which escapes at every f0 in the model, so it does not help.
- Two other readings were tried and do not change the conclusion: the return section taken as the
  half-line rotating with the primaries (theta-tilde advancing by 2 pi, as the proof of
  Proposition 1, p6, can be read), and the periapsis condition A7 used with its printed sign
  (p_r = +e r sin f / (1 + e cos f), which gives a radial component; the A4 relations need the
  opposite sign). The model's boundary speed versus f0 is the same to 5e-4 in all four readings.

The failing checks are strict xfails. The definitional controls (a circular start is weakly stable,
a start at 1.2 v_e, the top of the paper's speed grid, is weakly unstable) and a check of the time
variable (the first return of a circular low orbit after its Kepler period) show the classifier is
not stuck on one answer.
"""

from __future__ import annotations

import math

import numpy as np
import pytest
from numpy.typing import NDArray
from scipy.integrate import solve_ivp

import cyclerfinder.core.er3bp as er3bp

FloatArray = NDArray[np.float64]

# Section 5, p8.
_MU = 3.003158242e-6
_ECC = 0.0167
_R_EARTH_KM = 6378.0
_STEP_V2_KMS = 0.0063
# The Earth's semimajor axis (the paper says only "a"); IAU 2012 au.
_AU_KM = 149_597_870.7

# p15: circular and escape speed at 12 754 km, in the paper's km/s scale.
_VC_REF_KMS = 5.496891
_VE_REF_KMS = 7.773778
_R_REF_KM = 12_754.0

# Figure 7 (p13) and Figure 8 (p15) launch points: (alpha deg, r2 km, v2 km/s).
_FIG7 = (45.0, 381_829.0, 1.20446)
_FIG8 = (135.0, 12_754.0, 7.63689)

_ESCAPE_KM = 0.2 * _AU_KM


def _paper_vc_kms(r_km: float) -> float:
    """The paper's circular speed at ``r_km``, scaled from its printed value at 12 754 km."""
    return _VC_REF_KMS * math.sqrt(_R_REF_KM / r_km)


def _sep(f: float, e: float) -> float:
    """Primary separation R(f) in units of a (A1)."""
    return (1.0 - e * e) / (1.0 + e * math.cos(f))


def _dfdt(f: float, e: float) -> float:
    """df/dt with G(m1 + m2) = 1 and a = 1 (Section 2, p4)."""
    return float((1.0 + e * math.cos(f)) ** 2 / (1.0 - e * e) ** 1.5)


def _initial_state(
    f0: float, alpha_deg: float, r_km: float, v_over_vc: float, e: float
) -> FloatArray:
    """Pulsating-rotating state [x, y, z, x', y', z'] (primes d/df) for a direct start at the
    periapsis of an osculating ellipse about P2, speed ``v_over_vc`` times the circular speed."""
    alpha = math.radians(alpha_deg)
    u = np.array([math.cos(alpha), math.sin(alpha)])
    w = np.array([-u[1], u[0]])
    d = r_km / _AU_KM
    sep = _sep(f0, e)
    dsep = sep * e * math.sin(f0) / (1.0 + e * math.cos(f0))
    rho = (d / sep) * u
    v_inertial = v_over_vc * math.sqrt(_MU / d) * w  # rotating axes, time unit 1/n
    rho_prime = (v_inertial / _dfdt(f0, e) - dsep * rho) / sep - np.array([-rho[1], rho[0]])
    return np.array([1.0 - _MU + rho[0], rho[1], 0.0, rho_prime[0], rho_prime[1], 0.0])


def _kepler_energy(f: float, state: FloatArray, e: float) -> tuple[float, float]:
    """Kepler energy about P2 (A5, written in Cartesian form) and the physical distance."""
    rho = np.array([state[0] - 1.0 + _MU, state[1]])
    rho_prime = np.array([state[3], state[4]])
    sep = _sep(f, e)
    dsep = sep * e * math.sin(f) / (1.0 + e * math.cos(f))
    vel = _dfdt(f, e) * (dsep * rho + sep * (rho_prime + np.array([-rho[1], rho[0]])))
    dist = sep * float(np.linalg.norm(rho))
    return 0.5 * float(vel @ vel) - _MU / dist, dist


def _classify(
    f0_deg: float, alpha_deg: float, r_km: float, v_over_vc: float, e: float = _ECC
) -> tuple[str, float]:
    """Return (label, true anomaly at return) with label "stable", "unstable", "collision",
    "escape" or "no return"."""
    f0 = math.radians(f0_deg)
    state0 = np.concatenate((_initial_state(f0, alpha_deg, r_km, v_over_vc, e), [0.0]))

    def rhs(f: float, s: FloatArray) -> FloatArray:
        d6 = er3bp.er3bp_eom(f, s[:6], _MU, e)
        px, py = s[0] - 1.0 + _MU, s[1]
        theta_rate = (px * s[4] - py * s[3]) / (px * px + py * py) + 1.0  # theta = theta~ + f
        return np.concatenate((d6, [theta_rate]))

    def ev_return(f: float, s: FloatArray) -> float:
        return float(s[6] - 2.0 * math.pi)

    def ev_collision(f: float, s: FloatArray) -> float:
        return _sep(f, e) * math.hypot(float(s[0]) - 1.0 + _MU, float(s[1])) - _R_EARTH_KM / _AU_KM

    def ev_escape(f: float, s: FloatArray) -> float:
        return _sep(f, e) * math.hypot(float(s[0]) - 1.0 + _MU, float(s[1])) - _ESCAPE_KM / _AU_KM

    ev_return.terminal = True  # type: ignore[attr-defined]
    ev_return.direction = 1  # type: ignore[attr-defined]
    ev_collision.terminal = True  # type: ignore[attr-defined]
    ev_collision.direction = -1  # type: ignore[attr-defined]
    ev_escape.terminal = True  # type: ignore[attr-defined]
    ev_escape.direction = 1  # type: ignore[attr-defined]

    sol = solve_ivp(
        rhs,
        (f0, f0 + 2.0 * math.pi),
        state0,
        method="DOP853",
        rtol=1e-11,
        atol=1e-14,
        events=(ev_return, ev_collision, ev_escape),
    )
    t_events, y_events = sol.t_events, sol.y_events
    assert t_events is not None and y_events is not None
    if t_events[1].size:
        return "collision", math.nan
    if t_events[2].size:
        return "escape", math.nan
    if not t_events[0].size:
        return "no return", math.nan
    f_ret = float(t_events[0][0])
    energy, _ = _kepler_energy(f_ret, y_events[0][0][:6], e)
    return ("stable" if energy < 0.0 else "unstable"), f_ret


_K7 = _FIG7[2] / _paper_vc_kms(_FIG7[1])
_K8 = _FIG8[2] / _VC_REF_KMS


@pytest.fixture(scope="module")
def fig7_sweep() -> dict[int, str]:
    return {f0: _classify(f0, _FIG7[0], _FIG7[1], _K7)[0] for f0 in range(360)}


@pytest.fixture(scope="module")
def fig8_sweep() -> dict[int, str]:
    return {f0: _classify(f0, _FIG8[0], _FIG8[1], _K8)[0] for f0 in range(360)}


def test_speed_scale_is_internally_consistent() -> None:
    """The printed v_e / v_c at 12 754 km is sqrt(2) (to 4e-8), so the two printed speeds share one
    scale and the ratio v2/v_c is meaningful. Ratios: Figure 7 1.198910, Figure 8 1.389311
    (v2/v_e = 0.982391)."""
    assert math.isclose(_VE_REF_KMS / _VC_REF_KMS, math.sqrt(2.0), rel_tol=1e-6)
    assert math.isclose(_K7, 1.198910, abs_tol=1e-6)
    assert math.isclose(_K8, 1.389311, abs_tol=1e-6)


def test_first_return_of_a_low_circular_orbit_comes_after_its_kepler_period() -> None:
    """Instrument check of the frame and time variable: a circular start at 2 R_E returns to its
    half-line after one Kepler period, converted to true anomaly at perihelion, to 1e-3. Measured
    relative difference 1.7e-7 (the return is weakly stable)."""
    label, f_ret = _classify(0.0, 135.0, _R_REF_KM, 1.0)
    assert label == "stable"
    d = _R_REF_KM / _AU_KM
    period = 2.0 * math.pi * math.sqrt(d**3 / _MU)
    expected = period * _dfdt(0.0, _ECC)
    assert f_ret == pytest.approx(expected, rel=1e-3)


@pytest.mark.parametrize("fig", [_FIG7, _FIG8], ids=["fig7", "fig8"])
@pytest.mark.parametrize("f0", [0.0, 180.0])
def test_definitional_controls(fig: tuple[float, float, float], f0: float) -> None:
    """At the paper's two launch points: a circular start (e3 = 0) is weakly stable, and a start
    at 1.2 v_e, the top of the paper's speed grid (p8; Kepler energy positive, unstable by the
    two-body rule of p3), is weakly unstable."""
    alpha, r_km, _ = fig
    assert _classify(f0, alpha, r_km, 1.0)[0] == "stable"
    assert _classify(f0, alpha, r_km, 1.2 * math.sqrt(2.0))[0] != "stable"


@pytest.mark.parametrize(("fig", "k"), [(_FIG7, _K7), (_FIG8, _K8)], ids=["fig7", "fig8"])
@pytest.mark.parametrize("f0", [0.0, 90.0, 180.0, 270.0])
def test_circular_problem_is_stable_for_every_f0(
    fig: tuple[float, float, float], k: float, f0: float
) -> None:
    """p13 and p15: with e = 0 both launch points are weakly stable for all f0. This does not
    discriminate here: the model is also weakly stable at these f0 with e = 0.0167."""
    assert _classify(f0, fig[0], fig[1], k, e=0.0)[0] == "stable"


def _band(sweep: dict[int, str], lo: int, hi: int) -> dict[int, str]:
    return {f0: sweep[f0] for f0 in range(lo, hi + 1)}


def test_fig7_first_stable_band(fig7_sweep: dict[int, str]) -> None:
    """p13: weakly stable for f0 in [0, 177 deg]. Measured: stable at all 178 integer degrees."""
    band = _band(fig7_sweep, 0, 177)
    assert all(v == "stable" for v in band.values()), band


@pytest.mark.xfail(
    strict=True,
    reason="#896: the model is weakly stable at every integer f0 in the printed unstable band "
    "[178, 215] deg (Figure 7, p13); its least stable f0 is near 330 deg instead",
)
def test_fig7_unstable_band(fig7_sweep: dict[int, str]) -> None:
    """p13: weakly unstable for f0 in [178, 215 deg]. Measured: stable at all 38 points."""
    band = _band(fig7_sweep, 178, 215)
    assert all(v != "stable" for v in band.values()), band


def test_fig7_second_stable_band_outside_327_331(fig7_sweep: dict[int, str]) -> None:
    """p13: weakly stable for f0 in [216, 360 deg). Measured: stable at 216-326 and 332-359."""
    band = {f0: v for f0, v in _band(fig7_sweep, 216, 359).items() if not 327 <= f0 <= 331}
    assert all(v == "stable" for v in band.values()), band


@pytest.mark.xfail(
    strict=True,
    reason="#896: at f0 = 327-331 deg the model's Kepler energy at the return is positive "
    "(the printed speed is 2e-6 above the model's lowest boundary speed there), where the paper "
    "prints weakly stable (Figure 7, p13)",
)
def test_fig7_second_stable_band_327_331(fig7_sweep: dict[int, str]) -> None:
    band = _band(fig7_sweep, 327, 331)
    assert all(v == "stable" for v in band.values()), band


def test_fig8_first_stable_band(fig8_sweep: dict[int, str]) -> None:
    """p15: weakly stable for f0 in [0, 118 deg]. Measured: stable at all 119 integer degrees."""
    band = _band(fig8_sweep, 0, 118)
    assert all(v == "stable" for v in band.values()), band


@pytest.mark.xfail(
    strict=True,
    reason="#896: the model is weakly stable at every integer f0 in the printed unstable band "
    "[119, 259] deg (Figure 8, p15); at v2/v_c = 1.389 the Kepler apoapsis is 3.5e5 km, a quarter "
    "of the Hill radius",
)
def test_fig8_unstable_band(fig8_sweep: dict[int, str]) -> None:
    """p15: weakly unstable for f0 in [119, 259 deg]. Measured: stable at all 141 points."""
    band = _band(fig8_sweep, 119, 259)
    assert all(v != "stable" for v in band.values()), band


def test_fig8_second_stable_band(fig8_sweep: dict[int, str]) -> None:
    """p15: weakly stable for f0 in [260, 360 deg). Measured: stable at all 100 integer degrees."""
    band = _band(fig8_sweep, 260, 359)
    assert all(v == "stable" for v in band.values()), band


# Each printed transition at the paper's 1 deg resolution: the last degree of one band and the
# first of the next. Rows whose unstable side the model does not reproduce are strict xfails.
_XF7 = pytest.mark.xfail(strict=True, reason="#896: model weakly stable here (Figure 7, p13)")
_XF8 = pytest.mark.xfail(strict=True, reason="#896: model weakly stable here (Figure 8, p15)")


@pytest.mark.parametrize(
    ("fig", "f0", "printed"),
    [
        pytest.param("7", 177, "stable", id="fig7-177"),
        pytest.param("7", 178, "unstable", marks=_XF7, id="fig7-178"),
        pytest.param("7", 215, "unstable", marks=_XF7, id="fig7-215"),
        pytest.param("7", 216, "stable", id="fig7-216"),
        pytest.param("8", 118, "stable", id="fig8-118"),
        pytest.param("8", 119, "unstable", marks=_XF8, id="fig8-119"),
        pytest.param("8", 259, "unstable", marks=_XF8, id="fig8-259"),
        pytest.param("8", 260, "stable", id="fig8-260"),
    ],
)
def test_printed_transitions(
    fig: str,
    f0: int,
    printed: str,
    fig7_sweep: dict[int, str],
    fig8_sweep: dict[int, str],
) -> None:
    sweep = fig7_sweep if fig == "7" else fig8_sweep
    assert (sweep[f0] == "stable") == (printed == "stable"), sweep[f0]


@pytest.mark.parametrize(
    ("fig", "vc_kms"),
    [
        pytest.param(_FIG7, _paper_vc_kms(_FIG7[1]), id="fig7"),
        pytest.param(
            _FIG8,
            _VC_REF_KMS,
            marks=pytest.mark.xfail(
                strict=True,
                reason="#896: one grid step above the printed Figure 8 speed the model is still "
                "weakly stable; its first non-stable speed is about 0.998 v_e (p15 point is "
                "0.982 v_e)",
            ),
            id="fig8",
        ),
    ],
)
def test_launch_point_is_within_one_speed_step_of_the_boundary(
    fig: tuple[float, float, float], vc_kms: float
) -> None:
    """Both points are "a stable point of WSTR ... near to the first boundary curve" (p13, p15).
    At f0 = 0 (stated stable) the printed speed is weakly stable and the next speed of the
    paper's grid, v2 + step_v2 = v2 + 0.0063 km/s (p8), is not. Measured for Figure 7: stable,
    then unstable (the model's boundary at f0 = 0 is v2/v_c = 1.1995, printed 1.19891)."""
    alpha, r_km, v2 = fig
    assert _classify(0.0, alpha, r_km, v2 / vc_kms)[0] == "stable"
    assert _classify(0.0, alpha, r_km, (v2 + _STEP_V2_KMS) / vc_kms)[0] != "stable"
