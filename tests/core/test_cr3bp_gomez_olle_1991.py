"""#896 (j): printed second-species periodic orbits of Gomez & Olle 1991 in ``core.cr3bp`` and
``core.er3bp``.

Sources:

* G. Gomez and M. Olle, "Second-species solutions in the circular and elliptic restricted
  three-body problem. I. Existence and asymptotic approximation", Celest. Mech. Dyn. Astron. 52,
  107-146 (1991), DOI 10.1007/BF00049446. Filed in the private paper corpus as
  gomez-olle-1991-second-species-solutions-circular-elliptic-restricted-three-body-problem-I-existence-asymptotic-cmda-52-107-doi-10.1007-BF00049446.pdf.
* G. Gomez and M. Olle, "Second-species solutions in the circular and elliptic restricted
  three-body problem. II. Numerical explorations", Celest. Mech. Dyn. Astron. 52, 147-166 (1991),
  DOI 10.1007/BF00049447. Filed in the private paper corpus as
  gomez-olle-1991-second-species-solutions-circular-elliptic-restricted-three-body-problem-II-numerical-explorations-cmda-52-147-doi-10.1007-BF00049447.pdf.

Printed values used (Part II, all at mu = 1e-6, "in what follows mu = 10^-6", p. 160; every
digit re-read from 300 dpi renderings of the pages on 2026-10-04):

* Fig. 8 caption, p. 159, circular problem: (a) type 2, T = 6.273357979998948,
  x = -2.176283412498242, ydot = 1.63874359406601, C = 2.969728009802066; (b) type 3,
  T = 6.294070698585960, x = -2.173709228554383, ydot = 1.634970253275049, C = 2.971971477141059.
  T is the half period (p. 154, "T: half-period of the SPS").
* p. 163, circular problem, the starting orbits of Section 4 "with half period T = k pi":
  A0 x = -2.22979002878, ydot = 2.22978812061 (T = pi); A1 x = -4.30865819030,
  ydot = 4.30865753635 (T = 3 pi); B1 x = 3.348331642398, ydot = -3.3483307149 (T = 2 pi).
* Captions of Figs. 14, 16 and 18, pp. 162, 164, 165, elliptic problem, e_m = 0.5:
  A_{0,1}: x = -4.898301112912180, ydot = 4.898298973960655 (half period pi);
  A_{1,3}: x = -8.884269706188229, ydot = 8.884269157978011 (3 pi);
  B_{1,2}: x = 6.467038249660822, ydot = -6.467037966687732 (2 pi). The half period j pi is
  given by the family subscript (p. 163, "the second subindex j gives the value of the
  half-period (T = j pi)"); Part I p. 109: symmetric periodic orbits of the elliptic problem
  cross the axis perpendicularly with the small primary at pericentre or apocentre and have
  period 2 k pi. The primaries are at pericentre at t = 0 (Part II p. 148).

Conventions (read or confirmed here):

* Frame. Big primary at -mu, small at 1 - mu, as in ``core.cr3bp``: the printed Fig. 8 C is
  reproduced to 1e-14 only in this frame. Part II p. 154 writes the Levi-Civita distance as
  |x + 1 - mu| (small primary at -(1 - mu)), which conflicts with p. 152 and with the numbers;
  for the elliptic orbits a frame with the big primary at the origin (x - mu) fails (control).
* Jacobi constant. The formula is not printed. The printed C equals ``jacobi_constant`` plus
  mu(1 - mu) (Szebehely's form, Part I ref. [16]) to 1e-14; without the constant it misses by
  1e-6 (control).
* Elliptic independent variable. Not stated in words. Part I p. 141, Remark 2, prints the
  synodic initial condition of a rectilinear solution as x_s = r0/(1 - e_m), xdot_s =
  -r0/(1 - e_m): a point at rest in the inertial frame at pericentre has synodic velocity -x
  only if the derivative is with respect to the true anomaly f (with respect to time it would
  be -x sqrt(1 + e_m)/(1 - e_m)^(3/2) = -3.46 x). ``core.er3bp`` uses f, the pulsating frame and
  primaries at -mu, 1 - mu; the time reading fails (control).
* Broucke types (p. 152 sketches): type 2 has the second crossing P'' between the primaries,
  type 3 beyond the small primary. The p. 155 text ("type 3: the particle passes on the left of
  the small primary") says the same in the mirrored frame of p. 154. The p. 152 rule
  sign(x(P'') - 1 + mu) = sign(epsilon'' sin eta) with epsilon'' = -1, 0 < eta < pi (Part I
  p. 109) puts the second crossing of the nearly rectilinear orbits between the primaries.

First-order closest approach (Part I): Delta = mu K0 with K0 = C0/(v1^2 sqrt(|x_m(t1)|^2 v1^2
- C0^2)), C0 = sqrt(1 - e_m^2) for a rectilinear generating ellipse (Lemma 7, p. 123); e1 =
[1 + (v_inf^2 Delta/mu)^2]^(1/2) and r_p = Delta ((e1 - 1)/(e1 + 1))^(1/2) (eq. (9), p. 122),
v_inf = |v1| to leading order. The generating ellipse is read off the printed x through Remark 2,
p. 141 (x = r0/(1 - e_m), r0 = -2 epsilon a0): a0 = |x|(1 - e_m)/2, falling radially onto the
small primary at t1 = k pi, where |x_m| = 1 + e_m (k odd, apocentre) or 1 - e_m (k even), the
small primary moves at vis-viva speed perpendicular to x_m and |v1|^2 = (2/|x_m| - 1/a0) + v_m^2.
For e_m = 0 this is Theorem 10's |v1|^2 = 3 - 1/a0 (p. 136). The prediction is in inertial units;
the pulsating-frame distance is r_p/|x_m|.

Integrator. Measured on 2026-10-04, with every quantity located by an event (y = 0, or the
periapsis condition (r - r_m) . v = 0) rather than at a fixed time:

* Circular, the three computed starting members (below), plain ``cr3bp_eom`` DOP853: periapsis
  angle |y_p|/r_p = 2.5e-6, 4.3e-6, 7.9e-6 at rtol = atol = 1e-12 and 8.3e-7, 5.1e-7, 5.1e-6 at
  1e-13, r_p within 1.6e-5 (1e-12) and 2.4e-6 (1e-13) of the regularised value; the r2-Sundman
  integration at 1e-13 gives 5.8e-7, 1.4e-7, 4.8e-6.
* Elliptic, Figs. 14, 16, 18, plain ``er3bp_eom`` DOP853: xdot/ydot at the crossing 8e-7,
  1.2e-5, 3.9e-5 at 1e-12 and 8e-8, 1.2e-6, 3.7e-6 at 1e-13; r2-Sundman 5.8e-7, 4.1e-6, 2.7e-5
  at 1e-12 and 5e-8, 6e-7, 2.4e-6 at 1e-13.

So the plain integrators resolve these 1e-7 passages to about the level of the r2-Sundman one; the
tests use dt/ds = r2 at 1e-13 (``core.cr3bp_regularized.sundman_rhs`` for the circular problem; the
same factor applied to ``core.er3bp.er3bp_eom`` here, a reparametrisation of the project's
equations, not a different model) because the plain integration at 1e-13 sometimes takes 1e5
function evaluations through the passage. What does NOT resolve the passage is a residual taken at a
FIXED time: near periapsis d(xdot)/dt ~ v_p^2/r_p ~ 1e8, so a 1e-12 error in the passage timing
moves xdot at f = k pi by 1e-4. Measured xdot at exactly f = k pi from the printed elliptic states
(plain, 1e-12): 4.7e-5, 2.4e-3, 5.3e-4. The periapsis angle and the crossing time are the
well-conditioned forms of the same condition. The Fig. 8 orbits stay 1.6e-3 from the small primary
and need no regularisation. The project has no Levi-Civita integrator, which is what the paper used.

Measured on 2026-10-04:

==========  ==============================================  ============================
orbit       half-period crossing                            closest approach / first order
==========  ==============================================  ============================
Fig. 8(a)   t - T = -4.3e-11, xdot = 1.2e-10, x'' = 0.998102 1.897e-3 (O(mu^nu), Table I)
Fig. 8(b)   t - T = +1e-11, xdot = 4.0e-11, x'' = 1.001634   1.635e-3
Fig. 14     f_p - pi = 1.0e-13, |y_p|/r_p = 9e-8             2.21646e-7 / 2.21641e-7
Fig. 16     f_p - 3 pi = 5.6e-13, |y_p|/r_p = 1.1e-6         9.5171e-8 / 9.5170e-8
Fig. 18     f_p - 2 pi = 4.6e-13, |y_p|/r_p = 4.2e-6         1.17134e-7 / 1.17135e-7
==========  ==============================================  ============================

Closure after the full period: Fig. 8(a) 1.3e-9, Fig. 8(b) 4.2e-10 (plain DOP853); Fig. 14
9.3e-7 (regularised). Figs. 16 and 18 return to 4.8e-5 and 5.6e-4: the two 1e-7 passages
amplify the half-period residual, so only Fig. 14 is asserted.

Finding (strict xfail). The three starting orbits of p. 163 do not cross perpendicularly at T = k pi
in this model. From the printed state the crossing near k pi has xdot = -1.10, -1.26, -1.20 (at t -
k pi = 1.35e-6, 1.33e-6, 1.46e-6) and the orbit passes 9.9e-6, 2.4e-5, 3.5e-5 from the small primary
(the paper's regular SPSSS pass at O(mu) = 1e-7). The printed values have 11 to 12 decimals; the
miss changes by about 1 to 3 per unit change of ydot, so rounding does not explain it, and neither
does mu = 1e-5 or 1e-4, the big primary at the origin, or the mirrored frame (each tried). Solving
for the member with periapsis on the axis at t = k pi (Newton on (x, ydot), regularised) gives A0
(-2.229784142403, 2.229782234233), A1 (-4.308649107687316, 4.308648453731817), B1 (3.348316277269,
-3.348315336354): offsets (dx, dydot) = (+5.89e-6, -5.89e-6), (+9.08e-6, -9.08e-6), (-1.537e-5,
+1.538e-5) from the printed values, with x + ydot (the inertial velocity at t = 0) equal to the
printed x + ydot to within the printed digits for A0, 5.5e-12 for A1 and 1.3e-8 for B1 (about 250
times B1's 10-decimal rounding, so weaker there). x + ydot is Theorem 10's family parameter dv
(inertial ydot at t = 0, pp. 136-137). INFERRED: the printed family parameter matches the T = k pi
member while the printed x is off by 6e-6 to 1.5e-5. These computed members (not printed values)
pass the small primary at the first-order distance the paper predicts (1.81067e-7, 1.12374e-7,
1.28520e-7 against 1.81064e-7, 1.12373e-7, 1.28516e-7), on the side the p. 152 rule gives.

ER3BP corrector check (``genome.er3bp_periodic.correct_er3bp_periodic``, used by
``genome.er3bp_continuation``; ``search.er3bp_periodic`` holds only a coordinate converter and a
monodromy helper). Its symmetric mode (free x, ydot; residual y, xdot at f = period_f, the half
period) is Part I's condition (5) for the elliptic problem. At its default tolerance 1e-10 it raises
ConvergenceError on all three printed elliptic orbits (line search stalls at residual 5.3e-8,
9.8e-6, 4.4e-6): its residual is xdot at the fixed f = k pi, ill-conditioned at a 1e-7 passage (see
Integrator; the printed state itself has residual 4.7e-5, 2.4e-3, 5.3e-4 there). With tol = 1e-5 it
accepts Figs. 14 and 18 after moving the printed state by 8.8e-12 and 1.8e-11 (Fig. 16 also, by
4.4e-11, but takes 7 s and is left out). Its own Radau full-period check returns 4.2e-5 and 2.2e-2,
above its 1e-5 independent tolerance, which only logs a warning: a continuation through such orbits
would accept them without notice. Since #930 that check raises ClosureError (the tests below
pass ``independent_tol=None`` to study the crossing condition alone, and assert the rejection).
"""

from __future__ import annotations

import functools
import math
from dataclasses import dataclass

import numpy as np
import pytest
from numpy.typing import NDArray
from scipy.integrate import solve_ivp

from cyclerfinder.core.cr3bp import cr3bp_eom, jacobi_constant
from cyclerfinder.core.cr3bp_regularized import sundman_rhs
from cyclerfinder.core.er3bp import ER3BPSystem, er3bp_eom
from cyclerfinder.genome.er3bp_periodic import ClosureError, correct_er3bp_periodic

FloatArray = NDArray[np.float64]

MU = 1.0e-6
E_M = 0.5
X_M = 1.0 - MU  # small primary, project frame
REG_TOL = 1e-13
PI = math.pi


@dataclass(frozen=True)
class Fig8:
    t_half: str
    x: str
    ydot: str
    c: str
    broucke_type: int


FIG8 = {
    "fig8a": Fig8(
        "6.273357979998948", "-2.176283412498242", "1.63874359406601", "2.969728009802066", 2
    ),
    "fig8b": Fig8(
        "6.294070698585960", "-2.173709228554383", "1.634970253275049", "2.971971477141059", 3
    ),
}

# Part II p. 163: (x, ydot, k) with half period k pi, as printed (do not pad).
STARTING = {
    "A0": ("-2.22979002878", "2.22978812061", 1),
    "A1": ("-4.30865819030", "4.30865753635", 3),
    "B1": ("3.348331642398", "-3.3483307149", 2),
}

# Members with periapsis on the axis at t = k pi, found on 2026-10-04 by a regularised Newton
# solve from the printed states (computed here, NOT printed; see the module docstring).
STARTING_MEMBER = {
    "A0": (-2.229784142403, 2.229782234233),
    "A1": (-4.308649107687316, 4.308648453731817),
    "B1": (3.348316277269, -3.348315336354),
}

# Part II captions of Figs. 14, 16, 18 (e_m = 0.5): (x, ydot, k), half period k pi.
ELLIPTIC = {
    "fig14_A01": ("-4.898301112912180", "4.898298973960655", 1),
    "fig16_A13": ("-8.884269706188229", "8.884269157978011", 3),
    "fig18_B12": ("6.467038249660822", "-6.467037966687732", 2),
}


def _first_order_rp(x: float, e_m: float, k: int, mu: float = MU) -> float:
    """Part I first-order closest approach, in the synodic (pulsating) unit; see the docstring."""
    a0 = abs(x) * (1.0 - e_m) / 2.0
    r_m = 1.0 + e_m if k % 2 else 1.0 - e_m
    v_m2 = 2.0 / r_m - 1.0  # vis-viva of the small primary, a = 1
    v_r2 = 2.0 / r_m - 1.0 / a0  # radial fall from apocentre 2 a0 to r_m
    v1_2 = v_r2 + v_m2
    c0 = math.sqrt(1.0 - e_m * e_m)
    k0 = c0 / (v1_2 * math.sqrt(r_m * r_m * v1_2 - c0 * c0))
    delta = mu * k0
    e1 = math.sqrt(1.0 + (v1_2 * delta / mu) ** 2)
    rp = delta * math.sqrt((e1 - 1.0) / (e1 + 1.0))
    return rp / r_m


def _er3bp_reg_rhs(s: float, y: FloatArray, mu: float, e: float) -> FloatArray:
    """dX/ds with ds = df/r2 (r2-Sundman reparametrisation of ``core.er3bp.er3bp_eom``)."""
    r2 = math.sqrt((y[0] - 1.0 + mu) ** 2 + y[1] ** 2 + y[2] ** 2)
    out = np.empty(7)
    out[:6] = er3bp_eom(float(y[6]), y[:6], mu, e) * r2
    out[6] = r2
    return out


def _cr3bp_reg_rhs(s: float, y: FloatArray, mu: float, e: float) -> FloatArray:
    return sundman_rhs(s, y, mu, "r2")


@dataclass(frozen=True)
class Passage:
    cross: FloatArray  # (x, y, z, xdot, ydot, zdot, time) at the y = 0 crossing nearest k pi
    peri: FloatArray  # 7-state at the periapsis about the small primary nearest k pi
    end: FloatArray  # 7-state at time 2 k pi (only when full=True)

    @property
    def rp(self) -> float:
        return math.hypot(self.peri[0] - X_M, self.peri[1])


@functools.cache
def _passage(
    x: float, ydot: float, k: int, elliptic: bool, full: bool = False, mu: float = MU
) -> Passage:
    """Regularised integration from (x, 0, 0, 0, ydot, 0) at time 0 to k pi (2 k pi if full)."""
    rhs = _er3bp_reg_rhs if elliptic else _cr3bp_reg_rhs
    t_half = k * PI
    t_end = 2.0 * t_half if full else t_half + 0.05

    def y_cross(s: float, y: FloatArray, *_a: object) -> float:
        return float(y[1])

    def periapsis(s: float, y: FloatArray, *_a: object) -> float:
        return float((y[0] - 1.0 + mu) * y[3] + y[1] * y[4])

    periapsis.direction = 1.0  # type: ignore[attr-defined]

    def stop(s: float, y: FloatArray, *_a: object) -> float:
        return float(y[6] - t_end)

    stop.terminal = True  # type: ignore[attr-defined]
    sol = solve_ivp(
        rhs,
        (0.0, 1e4),
        np.array([x, 0.0, 0.0, 0.0, ydot, 0.0, 0.0]),
        args=(mu, E_M if elliptic else 0.0),
        method="DOP853",
        rtol=REG_TOL,
        atol=REG_TOL,
        events=[y_cross, periapsis, stop],
    )
    assert sol.y_events is not None
    near = [y for y in sol.y_events[0] if abs(y[6] - t_half) < 0.05]
    peris = [y for y in sol.y_events[1] if abs(y[6] - t_half) < 0.05]
    assert near and peris, "no axis crossing or periapsis within 0.05 of k pi"
    cross = min(near, key=lambda y: abs(y[6] - t_half))
    peri = min(peris, key=lambda y: math.hypot(y[0] - 1.0 + mu, y[1]))
    end = sol.y_events[2][0] if full else np.full(7, np.nan)
    return Passage(cross, peri, end)


# ---------------------------------------------------------------------------------------------
# Circular problem, Fig. 8 (p. 159)
# ---------------------------------------------------------------------------------------------


@pytest.mark.parametrize("key", sorted(FIG8))
def test_fig8_printed_jacobi_constant_includes_mu_one_minus_mu(key: str) -> None:
    """Measured: |C_project + mu(1 - mu) - C_printed| = 9e-15 (a), 5e-16 (b); without the constant
    the difference is 1.0e-6 (control)."""
    o = FIG8[key]
    s = np.array([float(o.x), 0.0, 0.0, 0.0, float(o.ydot), 0.0])
    c_project = jacobi_constant(s, MU)
    assert abs(c_project + MU * (1.0 - MU) - float(o.c)) < 1e-13
    assert abs(c_project - float(o.c)) > 5e-7  # control: the convention without the constant


@functools.cache
def _fig8_run(key: str) -> tuple[FloatArray, FloatArray, float, FloatArray]:
    """Plain ``cr3bp_eom`` DOP853 at 1e-13 over the full period: (first y = 0 crossing after t = 0,
    the crossing nearest T, the closest approach to the small primary, the state at 2T)."""
    o = FIG8[key]
    s0 = np.array([float(o.x), 0.0, 0.0, 0.0, float(o.ydot), 0.0])
    t_half = float(o.t_half)

    def y_cross(t: float, y: FloatArray, *_a: object) -> float:
        return float(y[1])

    sol = solve_ivp(
        cr3bp_eom,
        (0.0, 2.0 * t_half),
        s0,
        args=(MU,),
        method="DOP853",
        rtol=1e-13,
        atol=1e-13,
        events=y_cross,
        dense_output=True,
    )
    assert sol.t_events is not None and sol.y_events is not None and sol.sol is not None
    ts, ys = sol.t_events[0], sol.y_events[0]
    i = int(np.argmin(np.abs(ts - t_half)))
    tt = np.linspace(t_half - 0.05, t_half + 0.05, 20001)
    yy = sol.sol(tt)
    r_min = float(np.hypot(yy[0] - X_M, yy[1]).min())
    return np.append(ys[i], ts[i]), s0, r_min, sol.y[:, -1]


@pytest.mark.parametrize("key", sorted(FIG8))
def test_fig8_perpendicular_crossing_at_printed_half_period(key: str) -> None:
    """Measured: crossing time minus printed T = -4.3e-11 (a), +1.0e-11 (b); xdot there 1.2e-10,
    4.0e-11; second crossing at x = 0.998102 (a, between the primaries: type 2) and 1.001634
    (b, beyond the small primary: type 3); closest approach 1.897e-3 and 1.635e-3, i.e. O(mu^nu)
    near the A0 / E21 bifurcation (Table I), not O(mu); closure after 2T 1.3e-9, 4.3e-10."""
    o = FIG8[key]
    cross, s0, r_min, end = _fig8_run(key)
    assert abs(cross[6] - float(o.t_half)) < 1e-9
    assert abs(cross[3]) < 1e-8
    assert abs(cross[1]) < 1e-12
    if o.broucke_type == 2:
        assert -MU < cross[0] < X_M
    else:
        assert cross[0] > X_M
    assert float(np.abs(end - s0).max()) < 1e-8
    assert 1e-4 < r_min < 1e-2


# ---------------------------------------------------------------------------------------------
# Circular problem, the starting orbits of Section 4 (p. 163)
# ---------------------------------------------------------------------------------------------


@pytest.mark.xfail(
    strict=True,
    raises=AssertionError,
    reason=(
        "#896: the printed A0, A1, B1 starting orbits (p. 163) do not cross perpendicularly at "
        "T = k pi: xdot = -1.10, -1.26, -1.20 there and the closest approach is 9.9e-6, 2.4e-5, "
        "3.5e-5 instead of O(mu); the members that do are 6e-6 to 1.5e-5 away in x"
    ),
)
@pytest.mark.parametrize("key", sorted(STARTING))
def test_starting_orbit_perpendicular_at_k_pi(key: str) -> None:
    x, ydot, k = STARTING[key]
    p = _passage(float(x), float(ydot), k, elliptic=False)
    assert abs(p.cross[3]) < 1e-3, f"xdot at the crossing near k pi = {p.cross[3]:.3e}"
    assert p.rp < 1e-6, f"closest approach {p.rp:.3e}"


@pytest.mark.parametrize("key", sorted(STARTING_MEMBER))
def test_starting_member_first_order_closest_approach(key: str) -> None:
    """The computed T = k pi member nearest each printed starting orbit (module docstring): its
    periapsis lies on the axis at t = k pi, on the side of the p. 152 rule, at the Part I first-
    order distance (from the PRINTED x) within 1 percent.

    Measured: t_p - k pi = 9e-13, -9e-14, 2.8e-13; |y_p|/r_p = 5.8e-7, 1.4e-7, 4.8e-6; r_p =
    1.81067e-7, 1.12374e-7, 1.28520e-7 against the prediction 1.81064e-7, 1.12373e-7, 1.28516e-7
    (agreement 1.5e-5, 8.9e-6, 2.9e-5). The
    1 percent bound is the first-order theory's, set after seeing the numbers."""
    x_m, ydot_m = STARTING_MEMBER[key]
    x_p, _, k = STARTING[key]
    p = _passage(x_m, ydot_m, k, elliptic=False)
    assert abs(p.peri[6] - k * PI) < 1e-9
    assert abs(p.peri[1]) / p.rp < 1e-4
    assert p.peri[0] < X_M  # sign(x'' - 1 + mu) = sign(epsilon'' sin eta) = -1
    predicted = _first_order_rp(float(x_p), 0.0, k)
    assert abs(p.rp / predicted - 1.0) < 1e-2


# ---------------------------------------------------------------------------------------------
# Elliptic problem, e_m = 0.5, Figs. 14, 16, 18
# ---------------------------------------------------------------------------------------------


@pytest.mark.parametrize("key", sorted(ELLIPTIC))
def test_elliptic_orbit_periapsis_on_axis_at_k_pi(key: str) -> None:
    """Measured (regularised, 1e-13): f_p - k pi = 1.0e-13, 5.6e-13, 4.6e-13; |y_p|/r_p = 9e-8,
    1.1e-6, 4.2e-6; xdot/ydot at the crossing 5e-8, 6e-7, 2.4e-6; r_p (pulsating) = 2.21646e-7,
    9.5171e-8, 1.17134e-7 against the first-order 2.21641e-7, 9.5170e-8, 1.17135e-7 (agreement
    2.4e-5, 1.1e-5, -9.1e-6; the 1 percent bound was set after seeing these)."""
    x, ydot, k = ELLIPTIC[key]
    p = _passage(float(x), float(ydot), k, elliptic=True)
    assert abs(p.peri[6] - k * PI) < 1e-10
    assert abs(p.cross[6] - k * PI) < 1e-10
    assert abs(p.peri[1]) / p.rp < 1e-4
    assert abs(p.cross[3] / p.cross[4]) < 1e-4
    assert p.peri[0] < X_M  # p. 152 rule, epsilon'' = -1
    predicted = _first_order_rp(float(x), E_M, k)
    assert abs(p.rp / predicted - 1.0) < 1e-2


def test_elliptic_fig14_closes_after_full_period() -> None:
    """Measured: max |X(2 pi) - X(0)| = 9.3e-7 (regularised, 1e-13). Figs. 16 and 18 give 4.8e-5
    and 5.6e-4 (two 1e-7 passages amplify the half-period residual) and are not asserted."""
    x, ydot, k = ELLIPTIC["fig14_A01"]
    p = _passage(float(x), float(ydot), k, elliptic=True, full=True)
    s0 = np.array([float(x), 0.0, 0.0, 0.0, float(ydot), 0.0])
    assert abs(p.end[6] - 2.0 * PI) < 1e-12
    assert float(np.abs(p.end[:6] - s0).max()) < 1e-5


@pytest.mark.parametrize("key", sorted(ELLIPTIC))
def test_elliptic_controls_time_derivative_and_origin_at_big_primary_fail(key: str) -> None:
    """Controls. (i) ydot read as d/dt (divided by df/dt = sqrt(1 + e)/(1 - e)^(3/2) = 3.46 at
    pericentre): measured, none of the three orbits crosses the axis or has a periapsis about the
    small primary within 0.05 of k pi. (ii) The printed x measured from the
    big primary (x - mu in the project frame): measured xdot/ydot at the crossing 0.51, 1.16,
    1.11 and closest approach 1.9e-6, 5.4e-6, 3.1e-6."""
    x, ydot, k = ELLIPTIC[key]
    fac = math.sqrt(1.0 + E_M) / (1.0 - E_M) ** 1.5
    with pytest.raises(AssertionError, match="no axis crossing or periapsis"):
        _passage(float(x), float(ydot) / fac, k, elliptic=True)
    p_o = _passage(float(x) - MU, float(ydot), k, elliptic=True)
    assert abs(p_o.cross[3] / p_o.cross[4]) > 0.05
    assert p_o.rp > 1e-6


# ---------------------------------------------------------------------------------------------
# genome.er3bp_periodic corrector (used by genome.er3bp_continuation)
# ---------------------------------------------------------------------------------------------


@pytest.mark.parametrize("key", ["fig14_A01", "fig18_B12"])
def test_er3bp_corrector_accepts_printed_orbit(key: str) -> None:
    """Symmetric half-period mode (free x, ydot; residual y, xdot at f = period_f = k pi, the
    printed half period) at tol = 1e-5, the level its fixed-time residual reaches at a 1e-7
    passage. Measured: converges after moving the printed state by 8.8e-12 (Fig. 14) and 1.8e-11
    (Fig. 18). At the default tol = 1e-10 it raises ConvergenceError (module docstring)."""
    x, ydot, k = ELLIPTIC[key]
    s0 = np.array([float(x), 0.0, 0.0, 0.0, float(ydot), 0.0])
    system = ER3BPSystem(MU, E_M, "P0", "P1")
    orbit = correct_er3bp_periodic(system, s0, k * PI, tol=1e-5, independent_tol=None)
    assert float(np.abs(orbit.state0 - s0).max()) < 1e-9
    # The Radau full-period check fails for both (4.2e-5 and 2.2e-2 measured before #930, against
    # the 1e-5 default bound): the corrector now refuses to return them.
    assert orbit.independent_residual > 1e-5
    with pytest.raises(ClosureError):
        correct_er3bp_periodic(system, s0, k * PI, tol=1e-5)
