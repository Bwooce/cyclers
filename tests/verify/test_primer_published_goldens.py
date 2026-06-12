"""Published golden values for the verify-layer primer / multi-impulse machinery.

Two independent published sources are wired here (golden-test discipline: the
EXPECTED side of every assertion is a value printed in the publication, never a
number our own code produced; reproduction was confirmed at printed precision
before each cell was wired):

* Iorfida, E., Palmer, P. L., & Roberts, M. (2016). "Geometric Approach to the
  Perpendicular Thrust Case for Trajectory Optimization." *Journal of Guidance,
  Control, and Dynamics* 39(5), 1059-1068. DOI 10.2514/1.G001525.
  The out-of-plane primer component ``p_h`` on a conic coast has the closed
  form ``p_h(nu) = A x(nu) + B y(nu)`` (their Eq. 19) with parallel straight
  ``p_h``-isolines ``y = m x + D p_h`` (Eq. 24). Their Sec. V.C worked example
  (Table 5, p. 1067) prints ``m, D`` and the line-ellipse intersection roots
  ``x_1, x_2`` for both boundary-condition cases, plus optimality verdicts —
  an independent ANALYTIC cross-check of our STM-integrated primer
  (:mod:`cyclerfinder.verify.primer`), exercised here end-to-end.

* Shakouri, A., Kiani, M., & Pourtakdoust, S. H. (2019). "A New Shape-Based
  Multiple-Impulse Strategy for Coplanar Orbital Maneuvers."
  arXiv:1905.04543v1 [math.OC] (preprint submitted to Acta Astronautica).
  Tables 2 (p. 14) and 3 (p. 18) print optimal two-impulse transfer costs
  (``J_c`` = sum of impulse norms, their Eq. 26; ``J_m`` = max impulse norm,
  their Eq. 27) with their optimal transfer times — independent goldens for
  our Lambert solver + Kepler propagation on multi-impulse coplanar cases.

mu resolution (Shakouri): the paper never states the Earth gravitational
parameter. Reproduction was run with both 398600.0 and 398600.4418 km^3/s^2;
every wired cell agrees with the published 4-decimal value under BOTH (the
difference enters below 1e-5 km/s). 398600.0 is used here.

DO-NOT-USE record (cells that FAILED reproduction; per golden-test discipline
they are documented, not wired — a clean negative is a valid outcome):

* Shakouri Table 2, "1-Impulse" row (J_c = 2.6305 km/s, t_f = 2631 / 3463 s):
  the initial (a=13756 km, e=0.5, omega=10 deg) and final (circular 13756 km)
  orbits intersect at true anomaly 120/240 deg on the initial orbit, where our
  independently computed single-impulse cost is 2.7864 km/s at BOTH points
  (symmetric geometry) — not 2.6305. The coast times from theta_12 = 270 deg
  to the intersections are 4039 / 14623 s — neither matches 2631 / 3463 s.
  No interpretation tried (either omega sign convention, either coast
  direction, time-from-perigee) reproduces the printed row.
* Shakouri Table 2, "2-Impulse (Remark 2.7)" row, t_f = 2315 s cell only: the
  closed-form transfer is a half revolution of the intermediate ellipse
  (a_2 = 10317 km), whose half-period is 5214 s. The J_c / J_m cells of the
  same row DO reproduce exactly and are wired below.
* Shakouri Tables 2-3, "2-Impulse Lambert (b)" rows (t_f and theta_f free):
  the feasible set is ambiguous ("allows the spacecraft to remain in its
  initial orbit before applying the first impulse", p. 12). With the
  departure point fixed at theta_12 our optimum is 1.93 km/s (case 1, vs
  printed 1.4677); with the departure anomaly also free our optima fall
  BELOW the printed values (1.337 vs 1.4677; 0.782 vs 0.7831; case 2: 2.508
  vs 2.5604, 1.283 vs 1.3344). Irreproducible as printed.
* Shakouri Table 3, Lambert (a) "J_c* = 5.1176 / J_m* = 7.9455" labels: as
  printed the row is internally impossible (a max impulse norm cannot exceed
  the sum of impulse norms; 7.9455 > 5.1176). The two VALUES both reproduce
  with the labels exchanged: 5.1176 is the min-over-t_f of the MAX impulse
  (at its printed t_f = 2894 s), and 7.9455 is the min-over-t_f of the SUM —
  see the wired tests below. The printed t_f = 3570 s for the 7.9455 cell
  does NOT reproduce (our minimiser of the sum sits at t_f ~ 5126 s):
  DO-NOT-USE for that t_f cell.
* Shakouri Tables 2-3, "Shape-Based" (SBGM) 2-/3-impulse rows: not attempted.
  There is no SBGM implementation in this project to validate; reproducing
  those rows would test a re-implementation of the paper's own method, not
  project code.
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy.optimize import minimize_scalar

from cyclerfinder.core.kepler import coe_to_rv, propagate
from cyclerfinder.core.lambert import lambert
from cyclerfinder.verify.primer import (
    PrimerVerdict,
    _coast_stm,
    _solve_primer_rate,
    primer_on_coast,
)

# ---------------------------------------------------------------------------
# Iorfida et al. 2016 — Sec. V.C worked example (Table 5, p. 1067)
# ---------------------------------------------------------------------------
# Canonical units mu = 1 DU^3/TU^2 (their Fig. 2 convention). Transfer ellipse
# e = 0.3, a = 1.0 DU, nu_0 = 30 deg, nu_f = 100 deg (body text, p. 1067; the
# Table 5 caption prints a = 0.1 DU but the printed r_0, r_f only reproduce
# with a = 1.0 — and m, D, x are a-independent anyway, their Eq. 23).

MU_CANONICAL = 1.0
IORFIDA_E = 0.3
IORFIDA_A = 1.0
IORFIDA_NU0 = np.radians(30.0)
IORFIDA_NUF = np.radians(100.0)
IORFIDA_L = 1.0 - IORFIDA_E**2  # their l = 1 - e^2 (semi-latus rectum, a = 1)

# Published values are printed to 2 decimals; half-ulp tolerance.
TOL_2DP = 0.005

ZHAT = np.array([0.0, 0.0, 1.0])


def _tof_between(nu0: float, nuf: float, a: float, e: float, mu: float) -> float:
    """Elliptic time of flight between true anomalies (Kepler's equation)."""

    def mean_anom(nu: float) -> float:
        ecc = 2.0 * np.arctan(np.sqrt((1.0 - e) / (1.0 + e)) * np.tan(nu / 2.0))
        return float(ecc - e * np.sin(ecc))

    dm = mean_anom(nuf) - mean_anom(nu0)
    if dm < 0.0:
        dm += 2.0 * np.pi
    return float(dm / np.sqrt(mu / a**3))


def _iorfida_primer(
    p_hf_sign: float, n_samples: int = 400
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """STM-propagated primer on the Sec. V.C transfer arc, pure out-of-plane BCs.

    Returns ``(times, ref_states, p)`` where ``p`` is the (n, 3) primer
    history recovered by OUR machinery (:func:`_coast_stm` +
    :func:`_solve_primer_rate`), pinned to ``p(0) = +z_hat`` and
    ``p(T) = p_hf_sign * z_hat`` (their Eq. 25 boundary conditions).
    """
    r0, v0 = coe_to_rv(IORFIDA_A, IORFIDA_E, IORFIDA_NU0, MU_CANONICAL)
    tof = _tof_between(IORFIDA_NU0, IORFIDA_NUF, IORFIDA_A, IORFIDA_E, MU_CANONICAL)
    times, ref, stms = _coast_stm(r0, v0, tof, MU_CANONICAL, n_samples=n_samples)
    phi = stms[-1]
    pdot0, ill, _rcond = _solve_primer_rate(phi[:3, :3], phi[:3, 3:], ZHAT, p_hf_sign * ZHAT)
    assert ill is False  # 70-degree arc: the two-point BVP is well-posed
    state0 = np.concatenate([ZHAT, pdot0])
    p = np.array([stms[k][:3, :] @ state0 for k in range(times.shape[0])])
    return times, ref, p


def _recover_m_d(ref: np.ndarray, p: np.ndarray) -> tuple[float, float, float]:
    """Fit our numeric p_h to the published form p_h = A x + B y (Eq. 19).

    The p_h-isolines are then ``y = m x + D p_h`` (Eq. 24) with ``m = -A/B``
    and ``D = 1/B``. Returns ``(m, D, max_fit_residual)``.
    """
    x, y, ph = ref[:, 0], ref[:, 1], p[:, 2]
    coeffs, _res, _rank, _sv = np.linalg.lstsq(np.column_stack([x, y]), ph, rcond=None)
    a_c, b_c = float(coeffs[0]), float(coeffs[1])
    resid = float(np.max(np.abs(a_c * x + b_c * y - ph)))
    return -a_c / b_c, 1.0 / b_c, resid


def test_iorfida_example_endpoints_match_published() -> None:
    """Construction gate: r_0 = [0.63, 0.36, 0], r_f = [-0.17, 0.95, 0] DU (p. 1067)."""
    r0, _v0 = coe_to_rv(IORFIDA_A, IORFIDA_E, IORFIDA_NU0, MU_CANONICAL)
    rf, _vf = coe_to_rv(IORFIDA_A, IORFIDA_E, IORFIDA_NUF, MU_CANONICAL)
    assert r0 == pytest.approx([0.63, 0.36, 0.0], abs=TOL_2DP)
    assert rf == pytest.approx([-0.17, 0.95, 0.0], abs=TOL_2DP)


@pytest.mark.parametrize("p_hf_sign", [1.0, -1.0])
def test_iorfida_ph_decouples_and_is_axby(p_hf_sign: float) -> None:
    """Structural (their Eqs. 15, 18-19): p_h decouples and p_h = A x + B y.

    With pure out-of-plane boundary conditions the in-plane primer pair must
    stay identically zero (the gravity-gradient eigendecomposition decouples
    the h-channel), and the out-of-plane component must match the perifocal
    position components in the analytic two-coefficient form — on OUR
    STM-integrated primer, not a re-derivation.
    """
    _times, ref, p = _iorfida_primer(p_hf_sign)
    assert float(np.max(np.abs(p[:, :2]))) < 1.0e-12
    _m, _d, resid = _recover_m_d(ref, p)
    assert resid < 1.0e-9


def test_iorfida_table5_same_sign_m_d_roots() -> None:
    """Table 5 (p. 1067), p_h0 = p_hf case: m = -0.74, D = 0.82, x = 0.08, -1.29.

    ``m, D`` are recovered from OUR STM-propagated primer profile; the roots
    are the intersections of the limiting line ``p_h = -1`` (``y = m x - D``,
    their Sec. V.A.1) with the transfer ellipse via their Eqs. 33-34 evaluated
    on our recovered ``(m, D)``.
    """
    _times, ref, p = _iorfida_primer(+1.0)
    m, d, _resid = _recover_m_d(ref, p)
    assert m == pytest.approx(-0.74, abs=TOL_2DP)  # Table 5, m (p_h0 = p_hf)
    assert d == pytest.approx(0.82, abs=TOL_2DP)  # Table 5, D

    # Eq. 34 on the p_h = -1 line (intercept -D): (l+m^2)x^2 + 2(mD+le)x + (D^2-l^2) = 0.
    el, ec = IORFIDA_L, IORFIDA_E
    dd = -d
    quad_a = el + m * m
    quad_b = m * dd + el * ec
    disc = quad_b * quad_b - (dd * dd - el * el) * (m * m + el)
    x1 = (-quad_b + np.sqrt(disc)) / quad_a
    x2 = (-quad_b - np.sqrt(disc)) / quad_a
    assert x1 == pytest.approx(0.08, abs=TOL_2DP)  # Table 5, x_1
    assert x2 == pytest.approx(-1.29, abs=TOL_2DP)  # Table 5, x_2


def test_iorfida_table5_opposite_sign_m_d_roots() -> None:
    """Table 5 (p. 1067), p_h0 != p_hf case: m = 2.85, D = -1.42, x = 0.21, -0.79.

    Roots via their Eqs. 40-41 (existence discriminants Delta_+/-) and Eq. 43,
    evaluated on the ``(m, D)`` recovered from OUR primer profile. Each Eq. 43
    pair's other root must coincide with the known endpoint abscissa (x_0 or
    x_f) — asserted as an internal consistency gate.
    """
    _times, ref, p = _iorfida_primer(-1.0)
    m, d, _resid = _recover_m_d(ref, p)
    assert m == pytest.approx(2.85, abs=TOL_2DP)  # Table 5, m (p_h0 != p_hf)
    assert d == pytest.approx(-1.42, abs=TOL_2DP)  # Table 5, D

    el, ec = IORFIDA_L, IORFIDA_E
    delta_plus = -d * d + el * m * m + 2.0 * ec * d * m + el  # Eq. 40
    delta_minus = -d * d + el * m * m - 2.0 * ec * d * m + el  # Eq. 41
    assert delta_plus > 0.0 and delta_minus > 0.0  # existence condition, Eq. 42

    quad_a = el + m * m
    x_0 = (-(m * d + el * ec) + np.sqrt(el * delta_plus)) / quad_a  # Eq. 43
    x_1 = (-(m * d + el * ec) - np.sqrt(el * delta_plus)) / quad_a
    x_f = (-(-m * d + el * ec) + np.sqrt(el * delta_minus)) / quad_a
    x_2 = (-(-m * d + el * ec) - np.sqrt(el * delta_minus)) / quad_a
    assert x_1 == pytest.approx(0.21, abs=TOL_2DP)  # Table 5, x_1
    assert x_2 == pytest.approx(-0.79, abs=TOL_2DP)  # Table 5, x_2
    # The duplicate roots are the endpoints themselves (P_0 and P_f on the ellipse).
    assert x_0 == pytest.approx(ref[0, 0], abs=1.0e-6)
    assert x_f == pytest.approx(ref[-1, 0], abs=1.0e-6)


def test_iorfida_published_verdicts() -> None:
    """Published optimality verdicts (p. 1067) via primer_on_coast.

    Same-sign case: "the transfer is either non-optimal or sub-optimal in both
    the anti-clockwise and clockwise directions" — our (anticlockwise) coast
    must violate the |p| <= 1 necessary condition (IMPROVABLE). Opposite-sign
    case: "the clockwise transfer is sub-optimal, whereas the anti-clockwise
    transfer is optimal" — our coast must satisfy it (OPTIMAL). Qualitative
    verdicts only; no unpublished magnitude is asserted.
    """
    r0, v0 = coe_to_rv(IORFIDA_A, IORFIDA_E, IORFIDA_NU0, MU_CANONICAL)
    tof = _tof_between(IORFIDA_NU0, IORFIDA_NUF, IORFIDA_A, IORFIDA_E, MU_CANONICAL)

    same = primer_on_coast(r0, v0, ZHAT, ZHAT, tof, mu=MU_CANONICAL, n_samples=400)
    assert same.verdict is PrimerVerdict.IMPROVABLE_ADD_IMPULSE
    assert same.max_primer_magnitude > 1.0

    opp = primer_on_coast(r0, v0, ZHAT, -ZHAT, tof, mu=MU_CANONICAL, n_samples=400)
    assert opp.verdict is PrimerVerdict.OPTIMAL_NECESSARY_CONDITIONS_MET
    assert opp.max_primer_magnitude <= 1.0 + 1.0e-6


# ---------------------------------------------------------------------------
# Shakouri et al. 2019 — Tables 2 (p. 14) and 3 (p. 18)
# ---------------------------------------------------------------------------
# Their planar-orbit convention (Eq. 1, p. 4): r = a(1-e^2)/(1+e cos(theta+omega))
# with theta the inertial polar angle, so the true anomaly is nu = theta+omega
# and the perifocal frame is rotated by -omega. mu: see module docstring.

MU_SHAKOURI = 398600.0

# Case study 1 inputs (Sec. 4.1, p. 13): eccentric LEO -> circular LEO.
CASE1_INITIAL = (13756.0, 0.5, 10.0)  # a (km), e, omega (deg)
CASE1_FINAL = (13756.0, 0.0, 60.0)
CASE1_THETA12 = 270.0  # departure polar angle (deg)
CASE1_THETA34 = 30.0  # arrival polar angle (deg)

# Case study 2 inputs (Sec. 4.2, p. 16): near-circular LEO -> Molniya.
CASE2_INITIAL = (6644.4, 0.01, 60.0)
CASE2_FINAL = (26562.0, 0.74105, 30.0)
CASE2_THETA12 = 45.0
CASE2_THETA34 = 15.0

# Published costs carry 4 decimals; reproduction agreed to <= 1e-4 km/s on
# every wired cell, so 5e-4 absorbs the print rounding plus the unstated mu.
TOL_COST = 5.0e-4


def _shakouri_state(
    elements: tuple[float, float, float], theta_deg: float
) -> tuple[np.ndarray, np.ndarray]:
    """Inertial state on an orbit at polar angle theta (their Eq. 1 convention)."""
    a_km, e, omega_deg = elements
    nu = np.radians(theta_deg + omega_deg)
    return coe_to_rv(a_km, e, nu, MU_SHAKOURI, arg_peri_rad=np.radians(-omega_deg))


def _two_impulse_costs(
    r1: np.ndarray,
    v1_orbit: np.ndarray,
    r2: np.ndarray,
    v2_orbit: np.ndarray,
    tof_s: float,
) -> tuple[float, float]:
    """(J_c, J_m) of the fixed-endpoint two-impulse Lambert transfer (Eqs. 26-27)."""
    sol = lambert(r1, r2, tof_s, mu=MU_SHAKOURI)[0]
    dv1 = float(np.linalg.norm(np.asarray(sol.v1) - v1_orbit))
    dv2 = float(np.linalg.norm(v2_orbit - np.asarray(sol.v2)))
    return dv1 + dv2, max(dv1, dv2)


def _case_states(
    case: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    if case == 1:
        r1, v1 = _shakouri_state(CASE1_INITIAL, CASE1_THETA12)
        r2, v2 = _shakouri_state(CASE1_FINAL, CASE1_THETA34)
    else:
        r1, v1 = _shakouri_state(CASE2_INITIAL, CASE2_THETA12)
        r2, v2 = _shakouri_state(CASE2_FINAL, CASE2_THETA34)
    return r1, v1, r2, v2


def test_shakouri_case1_lambert_a_min_sum_golden() -> None:
    """Table 2 (p. 14), 2-Impulse Lambert (a), CE row: J_c* = 4.4539 km/s at t_f = 3750 s.

    Fixed endpoints (theta_12 = 270 deg on the initial orbit, theta_34 = 30 deg
    on the final orbit, a 120-deg prograde transfer); t_f is the paper's only
    free variable, so OUR Lambert cost evaluated at their printed optimal t_f
    must reproduce the printed optimal cost — and be locally stationary there.
    """
    r1, v1, r2, v2 = _case_states(1)
    j_c, _j_m = _two_impulse_costs(r1, v1, r2, v2, 3750.0)
    assert j_c == pytest.approx(4.4539, abs=TOL_COST)  # Table 2, J_c*
    # Local-minimum sanity (not a golden): the printed t_f is their optimiser's.
    assert j_c < _two_impulse_costs(r1, v1, r2, v2, 3600.0)[0]
    assert j_c < _two_impulse_costs(r1, v1, r2, v2, 3900.0)[0]


def test_shakouri_case1_lambert_a_min_max_golden() -> None:
    """Table 2 (p. 14), 2-Impulse Lambert (a), MI row: J_m* = 2.2989 km/s at t_f = 3184 s."""
    r1, v1, r2, v2 = _case_states(1)
    _j_c, j_m = _two_impulse_costs(r1, v1, r2, v2, 3184.0)
    assert j_m == pytest.approx(2.2989, abs=TOL_COST)  # Table 2, J_m*
    assert j_m < _two_impulse_costs(r1, v1, r2, v2, 3050.0)[1]
    assert j_m < _two_impulse_costs(r1, v1, r2, v2, 3320.0)[1]


def test_shakouri_case1_remark27_closed_form_golden() -> None:
    """Table 2 (p. 14), 2-Impulse (Remark 2.7) row: J_c = 1.5210, J_m = 0.9878 km/s.

    The closed-form generalized Hohmann (their Remark 2.7 / Eq. 16, p. 7):
    first tangent impulse at the perigee of the initial orbit onto the
    intermediate ellipse with e_2 = (a_4 - a_1(1-e_1)) / (a_4 + a_1(1-e_1)),
    second tangent impulse half a revolution later onto the circular target.
    The transfer coast runs through OUR universal-variable propagator; arrival
    on the target radius is asserted structurally.

    DO-NOT-USE: the row's t_f = 2315 s cell. The transfer is a half revolution
    of the a_2 = 10317 km intermediate ellipse, whose half-period is 5214 s;
    no interpretation reproduces 2315 s (see module docstring).
    """
    a_1, e_1, omega_1 = CASE1_INITIAL
    a_4 = CASE1_FINAL[0]
    r_p1 = a_1 * (1.0 - e_1)
    e_2 = (a_4 - r_p1) / (a_4 + r_p1)  # Eq. 16
    a_2 = r_p1 / (1.0 - e_2)

    # Both perigees sit at polar angle -omega_1 (nu = 0, shared apse line).
    _r_pi, v_peri_initial = _shakouri_state(CASE1_INITIAL, -omega_1)
    r_peri, v_peri_transfer = _shakouri_state((a_2, e_2, omega_1), -omega_1)
    dv1 = float(np.linalg.norm(v_peri_transfer - v_peri_initial))

    half_period_s = float(np.pi * np.sqrt(a_2**3 / MU_SHAKOURI))
    r_apo, v_apo = propagate(r_peri, v_peri_transfer, half_period_s, MU_SHAKOURI)
    r_apo_norm = float(np.linalg.norm(r_apo))
    # Structural: the coast arrives ON the circular target orbit.
    assert r_apo_norm == pytest.approx(a_4, abs=1.0)
    ang = np.arctan2(r_apo[1], r_apo[0])
    v_circ = np.sqrt(MU_SHAKOURI / r_apo_norm) * np.array([-np.sin(ang), np.cos(ang), 0.0])
    dv2 = float(np.linalg.norm(v_circ - v_apo))

    assert dv1 + dv2 == pytest.approx(1.5210, abs=TOL_COST)  # Table 2, J_c
    assert max(dv1, dv2) == pytest.approx(0.9878, abs=TOL_COST)  # Table 2, J_m


def test_shakouri_case2_lambert_a_min_max_golden() -> None:
    """Table 3 (p. 18), 2-Impulse Lambert (a): min-max impulse 5.1176 km/s at t_f = 2894 s.

    LABEL NOTE (see module docstring): Table 3 prints this value in the CE row
    (J_c*), but as printed the row pair is impossible (J_m* = 7.9455 > J_c* =
    5.1176, and a max impulse norm can never exceed the sum). The value
    reproduces as the MAX-impulse optimum: our J_m at the printed t_f = 2894 s
    equals 5.1176, and our t_f-minimiser of J_m sits at the printed t_f. The
    geometry is the 330-deg prograde transfer (theta: 45 deg -> 15 deg).
    """
    r1, v1, r2, v2 = _case_states(2)
    _j_c, j_m = _two_impulse_costs(r1, v1, r2, v2, 2894.0)
    assert j_m == pytest.approx(5.1176, abs=2.0e-4)  # Table 3, printed 5.1176
    res = minimize_scalar(
        lambda t: _two_impulse_costs(r1, v1, r2, v2, float(t))[1],
        bounds=(2000.0, 4000.0),
        method="bounded",
    )
    assert res.fun == pytest.approx(5.1176, abs=2.0e-4)
    assert res.x == pytest.approx(2894.0, abs=10.0)  # Table 3, t_f


def test_shakouri_case2_lambert_a_min_sum_golden() -> None:
    """Table 3 (p. 18), 2-Impulse Lambert (a): min-sum cost 7.9455 km/s (label swapped).

    The printed "J_m* = 7.9455" reproduces as the t_f-free MINIMUM of the SUM
    of impulses on the same fixed-endpoint geometry (see module docstring and
    the label note above; their own discussion, p. 19, has Lambert (a) as the
    "worst (highest) costs", matching this as the CE optimum).

    DO-NOT-USE: the printed t_f = 3570 s for this cell. Our minimiser sits at
    t_f ~ 5126 s; our J_c(3570 s) = 8.2459, nowhere near 7.9455.
    """
    r1, v1, r2, v2 = _case_states(2)
    res = minimize_scalar(
        lambda t: _two_impulse_costs(r1, v1, r2, v2, float(t))[0],
        bounds=(4000.0, 7000.0),
        method="bounded",
    )
    assert res.fun == pytest.approx(7.9455, abs=TOL_COST)  # Table 3, printed 7.9455
