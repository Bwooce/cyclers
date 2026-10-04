"""#896 (c): printed periodic orbits of the second species close in ``core.cr3bp``.

Sources (both at mu = 1e-4, C_J = 2.8, encounter circle of radius mu^alpha, alpha = 2/5):

* J. Font, A. Nunes and C. Simo, "Consecutive quasi-collisions in the planar circular RTBP",
  Nonlinearity 15, 115-142 (2002), DOI 10.1088/0951-7715/15/1/306. Captions of Figures 8, 9 and
  10, pp. 138-139: three periodic orbits, (phi, psi) and "stability parameter", no period.
  Filed in the private paper corpus as
  font-nunes-simo-2002-consecutive-quasi-collisions-planar-circular-RTBP-nonlinearity-15-115-doi-10.1088-0951-7715-15-1-306.pdf.
* J. Font, A. Nunes and C. Simo, "A numerical study of the orbits of second species of the planar
  circular RTBP", Celest. Mech. Dyn. Astron. 103, 143-162 (2009), DOI 10.1007/s10569-008-9176-z.
  Tables 3 and 4, p. 155: five periodic orbits, (phi, psi), period and stability parameter.
  Filed in the private paper corpus as
  font-nunes-simo-2009-numerical-study-orbits-second-species-planar-circular-RTBP-cmda-103-143-doi-10.1007-s10569-008-9176-z.pdf.

Convention, read from the papers (not inferred):

* Frame (2009 eq. (1.1), p. 145; 2002 eq. (1)): synodic, big primary E at (mu, 0), small primary
  M at (mu - 1, 0), equations ``xddot - 2 ydot = Omega_x``. ``core.cr3bp`` has E at (-mu, 0), M at
  (1 - mu, 0) and the same Coriolis sign, so the paper frame is the project frame turned by pi
  (X = -x, Y = -y): same sense of rotation, every angle shifted by pi.
* Jacobi constant (2009 eq. (1.2)): C_J = 2 Omega - v^2 with Omega including the constant
  mu(1 - mu)/2. ``core.cr3bp.jacobi_constant`` omits it, so C_J = 2.8 is a project value of
  2.8 - mu(1 - mu) = 2.79990001.
* psi is "the angle between the x-axis and the synodic velocity" (2009 p. 151; also the mu = 0
  initial condition (xdot, ydot) = (v cos psi, v sin psi) on p. 146).
* phi is the position angle on the circle of radius mu^alpha centred on M, from the +x axis,
  counter-clockwise: the marked initial points in 2002 Figures 8(b), 9(b), 10(b) sit at
  (x - x_M, y) of about (0.023, 0.012), (0.013, 0.022) and (0.015, 0.020), polar angles 0.48,
  1.06, 0.93, the printed phi. The points are exits from the disk (2009 p. 145: "the set of
  initial conditions exiting the circle").

What is asserted (``_shoot``): the return map on the circle is the composition of the outer arcs
(exit to next entry) and the passages inside the disk (entry to next exit). The orbits are very
unstable (printed stability parameters 2.6e4 to 2.9e13), so one double-precision integration of
a full period cannot show closure for most of them. Instead a multiple-shooting Newton solve with
a node at every crossing of the circle (2 per encounter, the number of encounters being the
number of printed symbols) finds the fixed point of the return map nearest the printed point,
with section-projected variational (``cr3bp_stm_eom``) Jacobians. Measured on 2026-10-04
(DOP853, rtol = atol = 1e-13):

=========  ======================  ===============  =====================  ===================
orbit      |fixed pt - printed|    |T - printed T|  trace / printed - 1    closest approach
=========  ======================  ===============  =====================  ===================
2002 F8    1.5e-14                 (not printed)    -0.9 (exponent, xfail)  1.2e-3
2002 F9    1.2e-13                 (not printed)    3.8e-7                 5.5e-4
2002 F10   1.7e-13                 (not printed)    3.4e-8                 2.5e-4
2009 12a   2.7e-12                 2.3e-12          1.6e-9                 4.2e-5
2009 12c   3.3e-14                 2.5e-14          7.0e-10                1.4e-3
2009 12e   4.2e-14                 6.0e-13          -2.4e-9                3.0e-5
2009 12g   8.9e-14                 1.4e-14          -3.5e-9                3.8e-5
2009 12i   5.7e-11                 2.8e-11          -2.7e-8                8.7e-6
=========  ======================  ===============  =====================  ===================

(closest approach to M in units of the primaries' separation; the circle radius is 2.5e-2.)
One integration of a period from the printed state (the literal check) misses by 3.2e-9 (12a),
4.5e-7 (Fig. 8), 7.7e-7 (Fig. 9), 6.6e-7 (12c) in the state norm.

The 12g phi: the digest of 2026-10-04 transcribed it as 0.98857090907035220406767336927899 (32
decimals); the table as printed (p. 155, re-read at 400 dpi) reads 0.988570907035220406767336927899
(30 decimals, like its psi). The digested string closes only to 2.0e-9 (an inserted "09"); the
printed one closes to 1e-12. The test uses the printed string.

The stability parameter is not defined in either paper. It agrees with the trace of the
linearised return map (lambda + 1/lambda, coordinate free at a fixed point), INFERRED.

Integrator error is bounded, not eliminated, by (i) Jacobi conservation at every crossing and
(ii) a Sundman (``core.cr3bp_regularized.sundman_rhs``, dt/ds = r2) recomputation of every
passage inside the disk: about 1e-11 for passages farther than 1e-4 from M, about 1e-9 for the
closer ones (12i comes within 8.7e-6; there both integrators drift by 6e-10 in C and differ by
2.3e-10 in the exit angles). These two bounds (``_passage_tolerance``) were set after measuring;
the physical bounds (fixed point 1e-10, period 1e-10, trace 1e-6) were set before. The project
has no Levi-Civita integrator, which is what the papers used.

Controls: the most plausible wrong readings (the digest's recipe applied literally in the
project frame; C_J without mu(1 - mu); angles clockwise; psi from the radial direction) miss
the start by more than 1e-3 after one period where
the right reading misses by less than |stability parameter| x 1e-11.
"""

from __future__ import annotations

import functools
import math
from collections.abc import Callable
from dataclasses import dataclass

import numpy as np
import pytest
from numpy.typing import NDArray
from scipy.integrate import solve_ivp

from cyclerfinder.core.cr3bp import cr3bp_eom, cr3bp_stm_eom, jacobi_constant
from cyclerfinder.core.cr3bp_regularized import sundman_rhs

FloatArray = NDArray[np.float64]

MU = 1.0e-4
C_J = 2.8
ALPHA = 0.4
RADIUS = MU**ALPHA  # encounter circle, 0.02512 (2009 p. 146; 2002 p. 133)
X_M = 1.0 - MU  # small primary in the project frame
C_PROJECT = C_J - MU * (1.0 - MU)  # the paper's C_J without its constant mu(1 - mu)
TOL = 1e-13
TWO_PI = 2.0 * math.pi


@dataclass(frozen=True)
class Printed:
    """One printed periodic orbit, digits as printed."""

    source: str
    phi: str
    psi: str
    symbols: tuple[tuple[str, int, int, int], ...]  # (X, p, q, s); 2002 orbits use X = "R"
    stability: float
    period: str | None


ORBITS: dict[str, Printed] = {
    # 2002, Figure 8 caption, p. 138: (4,3,1) n odd, (2,1,1) n even.
    "2002_fig8": Printed(
        "2002 Fig. 8",
        "0.4840458051695259",
        "0.3657582929540782",
        (("R", 4, 3, 1), ("R", 2, 1, 1)),
        2.338645e7,
        None,
    ),
    # 2002, Figure 9 caption, p. 138: (2,1,1), (2,1,-1); symmetric.
    "2002_fig9": Printed(
        "2002 Fig. 9",
        "1.060031115468069",
        "0.9672594502381595",
        (("R", 2, 1, 1), ("R", 2, 1, -1)),
        1.132360e6,
        None,
    ),
    # 2002, Figure 10 caption, p. 139: (2,1,1), (2,1,-1), (3,2,1).
    "2002_fig10": Printed(
        "2002 Fig. 10",
        "0.91622091785178612688",
        "0.94492534193201301129",
        (("R", 2, 1, 1), ("R", 2, 1, -1), ("R", 3, 2, 1)),
        -4.228471e9,
        None,
    ),
    # 2009, Table 3, p. 155.
    "2009_12a": Printed(
        "2009 Table 3, Fig. 12a",
        "2.685003594282268",
        "2.639352410113041",
        (("S", 1, 1, -1),),
        0.2639321981e05,
        "7.933918152222289",
    ),
    "2009_12c": Printed(
        "2009 Table 3, Fig. 12c",
        "3.0141823898790555595",
        "3.0316423709511530133",
        (("S", 1, 2, -1), ("S", 3, 1, 1)),
        0.1946971773492e08,
        "27.003867331650326275",
    ),
    "2009_12e": Printed(
        "2009 Table 3, Fig. 12e",
        "2.5457836425596942403",
        "2.4204845762785093917",
        (("S", 2, 1, 1), ("R", 3, 2, -1)),
        -0.1201986472185e09,
        "22.310316405350279103",
    ),
    # 2009, Table 4, p. 155. phi as printed (the 2026-10-04 digest has an extra "09").
    "2009_12g": Printed(
        "2009 Table 4, Fig. 12g",
        "0.988570907035220406767336927899",
        "0.954760770266642116216411782118",
        (("R", 2, 1, 1), ("S", 2, 3, 1), ("R", 1, 3, 1)),
        0.2867317722679597242e14,
        "45.8085897638589254831031186215",
    ),
    "2009_12i": Printed(
        "2009 Table 4, Fig. 12i",
        "3.27824276059703769776947746511",
        "3.21924695055804868387324783072",
        (("R", 1, 2, -1), ("S", 1, 2, -1), ("R", 1, 2, 1)),
        0.9608334138976017932e11,
        "40.8236298146304851361874699996",
    ),
}

# The phi of 12g as transcribed in docs/notes/2026-10-04-digest-font-nunes-simo-2009-...md.
DIGEST_12G_PHI = "0.98857090907035220406767336927899"


def _wrap(a: float) -> float:
    return (a + math.pi) % TWO_PI - math.pi


def _state(phi: float, psi: float, reading: str = "paper") -> FloatArray:
    """Project-frame state on the encounter circle for printed (phi, psi) under ``reading``.

    ``paper`` is the convention read from the papers (module docstring). The others are the
    wrong readings used as controls.
    """
    c_target = C_PROJECT
    if reading == "clockwise":
        phi, psi = -phi, -psi
    elif reading == "psi_from_radial":
        psi = phi + psi
    elif reading == "no_mu_term":
        c_target = C_J
    if reading == "digest_literal":
        # Angles about M at (1 - mu, 0) in the project frame, no rotation by pi.
        pos_angle, vel_angle = phi, psi
    else:
        pos_angle, vel_angle = phi + math.pi, psi + math.pi
    s = np.array([X_M + RADIUS * math.cos(pos_angle), RADIUS * math.sin(pos_angle), 0, 0, 0, 0])
    speed = math.sqrt(jacobi_constant(s, MU) - c_target)
    s[3] = speed * math.cos(vel_angle)
    s[4] = speed * math.sin(vel_angle)
    return s


def _angles(s: FloatArray) -> FloatArray:
    """Inverse of ``_state`` (paper reading): (phi, psi) of a state on the circle."""
    return np.array([math.atan2(-s[1], -(s[0] - X_M)), math.atan2(-s[4], -s[3])], dtype=np.float64)


def _angles_jacobian(s: FloatArray) -> FloatArray:
    dx, dy = s[0] - X_M, s[1]
    rho2 = dx * dx + dy * dy
    v2 = s[3] ** 2 + s[4] ** 2
    a = np.zeros((2, 6))
    a[0, 0], a[0, 1] = -dy / rho2, dx / rho2
    a[1, 3], a[1, 4] = -s[4] / v2, s[3] / v2
    return a


def _state_jacobian(z: FloatArray) -> FloatArray:
    """d state / d (phi, psi) on the Jacobi level (central differences of a smooth map)."""
    h = 1e-6
    b = np.zeros((6, 2))
    for j in range(2):
        dz = np.zeros(2)
        dz[j] = h
        b[:, j] = (_state(*(z + dz)) - _state(*(z - dz))) / (2 * h)
    return b


def _circle_event(direction: float) -> Callable[[float, FloatArray], float]:
    def ev(t: float, y: FloatArray, *_args: object) -> float:
        return float(math.hypot(y[0] - X_M, y[1]) - RADIUS)

    ev.direction = direction  # type: ignore[attr-defined]
    ev.terminal = True  # type: ignore[attr-defined]
    return ev


@dataclass(frozen=True)
class Leg:
    end: FloatArray  # (phi, psi) at the next crossing
    time: float
    jac: FloatArray  # 2x2 section-to-section Jacobian
    r2_min: float
    jacobi_end: float
    end_state: FloatArray


def _leg(z: FloatArray, kind: str) -> Leg:
    """From a node on the circle: an exit runs to the next entry, an entry to the next exit.

    The start is exactly on the circle, so the event direction is chosen to exclude t = 0: an
    exit (r increasing) looks for a decreasing crossing and vice versa.
    """
    s0 = _state(float(z[0]), float(z[1]))
    y0 = np.concatenate([s0, np.eye(6).ravel()])
    sol = solve_ivp(
        cr3bp_stm_eom,
        (0.0, 40.0),
        y0,
        args=(MU,),
        method="DOP853",
        rtol=TOL,
        atol=TOL,
        events=_circle_event(-1.0 if kind == "exit" else 1.0),
    )
    assert sol.t_events is not None and sol.y_events is not None
    assert sol.t_events[0].size == 1, f"no crossing of the circle on a {kind} leg"
    ye = sol.y_events[0][0]
    se = ye[:6]
    phi_t = ye[6:].reshape(6, 6)
    f = cr3bp_eom(0.0, se, MU)
    grad = np.zeros(6)
    rho = math.hypot(se[0] - X_M, se[1])
    grad[0], grad[1] = (se[0] - X_M) / rho, se[1] / rho
    proj = np.eye(6) - np.outer(f, grad) / float(grad @ f)
    jac = _angles_jacobian(se) @ proj @ phi_t @ _state_jacobian(z)
    r2 = np.hypot(sol.y[0] - X_M, sol.y[1])
    return Leg(
        _angles(se), float(sol.t_events[0][0]), jac, float(r2.min()), jacobi_constant(se, MU), se
    )


def _crossings(s0: FloatArray, t_end: float) -> list[tuple[float, FloatArray]]:
    sol = solve_ivp(
        cr3bp_eom,
        (0.0, t_end),
        s0,
        args=(MU,),
        method="DOP853",
        rtol=TOL,
        atol=TOL,
        events=lambda t, y, *_a: math.hypot(y[0] - X_M, y[1]) - RADIUS,
    )
    assert sol.t_events is not None and sol.y_events is not None
    return [(float(t), y) for t, y in zip(sol.t_events[0], sol.y_events[0], strict=True)]


@dataclass(frozen=True)
class Shot:
    nodes: FloatArray  # (2n, 2): exit, entry, exit, entry, ...
    legs: tuple[Leg, ...]
    residual: float

    @property
    def period(self) -> float:
        return sum(leg.time for leg in self.legs)

    @property
    def monodromy(self) -> FloatArray:
        m = np.eye(2)
        for leg in self.legs:
            m = leg.jac @ m
        return m


@functools.cache
def _shoot(key: str) -> Shot:
    orb = ORBITS[key]
    phi0, psi0 = float(orb.phi), float(orb.psi)
    n = len(orb.symbols)
    t_guess = float(orb.period) if orb.period else TWO_PI * sum(q for _, _, q, _ in orb.symbols)
    # Initial nodes: forward crossings for the first n // 2 encounters, backward for the rest,
    # so no guess is taken more than half a period from the printed point.
    s0 = _state(phi0, psi0)
    nf = n // 2
    fw = [c for c in _crossings(s0, 0.75 * t_guess) if c[0] > 1e-9][: 2 * nf]
    bw = [c for c in _crossings(s0, -0.75 * t_guess) if c[0] < -1e-9][: 2 * (n - nf) - 1]
    nodes = [np.array([phi0, psi0])] + [_angles(y) for _, y in fw + bw[::-1]]
    assert len(nodes) == 2 * n
    z = np.array(nodes)
    m = 2 * n
    kinds = ["exit" if k % 2 == 0 else "entry" for k in range(m)]
    best: Shot | None = None
    for _ in range(8):
        legs = tuple(_leg(z[k], kinds[k]) for k in range(m))
        f = np.zeros(2 * m)
        jac = np.zeros((2 * m, 2 * m))
        for k, leg in enumerate(legs):
            kn = (k + 1) % m
            f[2 * k] = _wrap(leg.end[0] - z[kn, 0])
            f[2 * k + 1] = _wrap(leg.end[1] - z[kn, 1])
            jac[2 * k : 2 * k + 2, 2 * k : 2 * k + 2] = leg.jac
            jac[2 * k : 2 * k + 2, 2 * kn : 2 * kn + 2] -= np.eye(2)
        res = float(np.max(np.abs(f)))
        shot = Shot(z.copy(), legs, res)
        if best is None or res < best.residual:
            best = shot
        if res < 1e-11:
            break
        z = z + np.linalg.solve(jac, -f).reshape(m, 2)
    assert best is not None
    return best


def _passage_tolerance(leg: Leg) -> float:
    """Error budget of one leg of the 1e-13 DOP853 integration.

    Measured: Jacobi drift and plain-versus-Sundman disagreement are below 1e-12 for every leg
    that stays farther than 1e-4 from M, and reach 1.9e-11 (12a, 4.2e-5), 2.6e-11 (12g, 3.8e-5),
    1.8e-11 (12e, 3.0e-5) and 6e-10 (12i, 8.7e-6) on the closest passages, for both integrators
    and both 1e-13 and 1e-14. That is the double-precision floor of an integration that is not
    Levi-Civita regularised (the papers used Levi-Civita and extended precision), not a model
    difference; it is why the fixed-point bound is 1e-10 and not 1e-12.
    """
    return 1e-11 if leg.r2_min > 1e-4 else 1e-9


def _direct_miss(key: str, reading: str) -> float:
    """State-space miss after one period of a single integration from the printed state.

    The period is the time of the n-th outward crossing of the circle (n = number of printed
    symbols); infinity if the orbit does not make n outward crossings within 1.5 periods.
    """
    orb = ORBITS[key]
    n = len(orb.symbols)
    t_guess = float(orb.period) if orb.period else TWO_PI * sum(q for _, _, q, _ in orb.symbols)
    s0 = _state(float(orb.phi), float(orb.psi), reading)
    exits = []
    for t, y in _crossings(s0, 1.5 * t_guess):
        if t > 1e-9 and (y[0] - X_M) * y[3] + y[1] * y[4] > 0:
            exits.append(y)
    if len(exits) < n:
        return math.inf
    return float(np.linalg.norm(exits[n - 1][:6] - s0))


# --------------------------------------------------------------------------------------------


@pytest.mark.parametrize("key", list(ORBITS))
def test_constructed_state_has_the_printed_jacobi_constant(key: str) -> None:
    """Paper eq. (1.2) in the paper frame gives 2.8; core.cr3bp gives 2.8 - mu(1 - mu)."""
    orb = ORBITS[key]
    s = _state(float(orb.phi), float(orb.psi))
    x, y, vx, vy = -s[0], -s[1], -s[3], -s[4]  # back to the paper frame
    r1 = math.hypot(x - MU, y)
    r2 = math.hypot(x - MU + 1.0, y)
    assert abs(r2 - RADIUS) < 1e-15
    omega = (1 - MU) / r1 + MU / r2 + 0.5 * (x * x + y * y) + 0.5 * MU * (1 - MU)
    assert 2 * omega - (vx * vx + vy * vy) == pytest.approx(C_J, abs=1e-14)
    assert jacobi_constant(s, MU) == pytest.approx(2.8 - 1e-4 * (1 - 1e-4), abs=1e-14)
    # Leaving the disk (2009 p. 145; the strips lie in |phi - psi| < pi/2).
    assert abs(_wrap(float(orb.phi) - float(orb.psi))) < math.pi / 2


@pytest.mark.parametrize("key", list(ORBITS))
def test_printed_point_is_the_fixed_point_of_the_return_map(key: str) -> None:
    """The fixed point of the n-encounter return map is within 1e-10 rad of the printed point.

    1e-10 is about 100 times the shooting residual floor (1e-12, set by the 1e-13 integrator
    tolerance at the passage inside the disk); the wrong readings are off by order 1.
    """
    orb = ORBITS[key]
    shot = _shoot(key)
    assert shot.residual < 1e-10
    assert abs(_wrap(shot.nodes[0, 0] - float(orb.phi))) < 1e-10
    assert abs(_wrap(shot.nodes[0, 1] - float(orb.psi))) < 1e-10
    for leg in shot.legs:  # Jacobi conserved at every crossing
        assert abs(leg.jacobi_end - C_PROJECT) < _passage_tolerance(leg)
    # The section map is conjugate to an area-preserving one: loop determinant 1.
    assert float(np.prod([np.linalg.det(leg.jac) for leg in shot.legs])) == pytest.approx(
        1.0, abs=1e-5
    )
    # Every resonant (T) arc spends q turns of the frame outside the disk, to O(mu^alpha)
    # (2002 Definition 1, p. 116). Outer legs are the even-numbered ones.
    outer = [shot.legs[k].time / TWO_PI for k in range(0, len(shot.legs), 2)]
    for x_type, _p, q, _s in orb.symbols:
        if x_type != "R":
            continue
        match = [i for i, turns in enumerate(outer) if abs(turns - q) < 2 * RADIUS]
        assert match, f"no outer arc of {q} turns: {outer}"
        outer.pop(match[0])


@pytest.mark.parametrize("key", [k for k, o in ORBITS.items() if o.period])
def test_period_matches_the_printed_period(key: str) -> None:
    """Sum of the leg times of the fixed point equals the printed period (2009 Tables 3, 4)."""
    shot = _shoot(key)
    assert shot.period == pytest.approx(float(ORBITS[key].period or "nan"), abs=1e-10)


def test_digested_12g_phi_is_not_the_printed_one() -> None:
    """The digest's 12g phi (an inserted '09') sits 2.0e-9 from the fixed point."""
    shot = _shoot("2009_12g")
    off = abs(_wrap(shot.nodes[0, 0] - float(DIGEST_12G_PHI)))
    assert 1.9e-9 < off < 2.2e-9


_FIG8_EXPONENT = pytest.mark.xfail(
    strict=True,
    reason="#896: 2002 Fig. 8 caption prints 2.338645E+7; the trace is 2.338645e6 (all seven "
    "printed digits agree, the exponent differs by one; Figs. 9 and 10 agree with exponent)",
)


@pytest.mark.parametrize(
    "key",
    [pytest.param(k, marks=_FIG8_EXPONENT) if k == "2002_fig8" else k for k in ORBITS],
)
def test_stability_parameter_is_the_trace_of_the_return_map(key: str) -> None:
    """Printed stability parameter = trace of the linearised n-encounter return map (INFERRED).

    Measured trace / printed - 1: 3.8e-7 (Fig. 9, printed to 7 digits), 3.4e-8 (Fig. 10),
    1.6e-9 (12a), 7.0e-10 (12c), -2.4e-9 (12e), -3.5e-9 (12g), -2.7e-8 (12i). Fig. 8: trace
    2.3386449e6 against the printed 2.338645E+7.
    """
    shot = _shoot(key)
    trace = float(np.trace(shot.monodromy))
    assert trace == pytest.approx(ORBITS[key].stability, rel=1e-6)
    # Half the trace is off by a factor 2: the definition is the trace, not the half-trace.
    assert abs(0.5 * trace / ORBITS[key].stability - 1.0) > 0.4


def test_fig8_trace_matches_the_caption_with_exponent_6() -> None:
    """2002 Fig. 8: the trace agrees with the caption's seven digits 2.338645 at exponent 6.

    The caption prints 2.338645E+7 (p. 138, re-read at 300 dpi). Read as E+6, the exponent
    corrected reading, it agrees to the same 1e-6 as Figs. 9 and 10 (measured -4.7e-8).
    """
    trace = float(np.trace(_shoot("2002_fig8").monodromy))
    assert trace == pytest.approx(2.338645e6, rel=1e-6)


@pytest.mark.parametrize("key", list(ORBITS))
def test_passages_agree_with_regularized_integration(key: str) -> None:
    """Each passage inside the disk, redone with the Sundman (dt/ds = r2) form, ends at the same
    (phi, psi) and time, to 1e-11 (1e-9 within 1e-4 of M, see ``_passage_tolerance``)."""
    shot = _shoot(key)
    for k in range(1, len(shot.legs), 2):  # entry legs
        z = shot.nodes[k]
        y0 = np.append(_state(float(z[0]), float(z[1])), 0.0)
        sol = solve_ivp(
            sundman_rhs,
            (0.0, 1e5),
            y0,
            args=(MU, "r2"),
            method="DOP853",
            rtol=TOL,
            atol=TOL,
            events=_circle_event(1.0),
        )
        assert sol.y_events is not None
        ye = sol.y_events[0][0]
        plain = shot.legs[k]
        tol = _passage_tolerance(plain)
        assert np.max(np.abs([_wrap(a) for a in _angles(ye[:6]) - plain.end])) < tol
        assert abs(ye[6] - plain.time) < tol
        assert abs(jacobi_constant(ye[:6], MU) - C_PROJECT) < tol
        assert plain.r2_min > 5e-6  # closest approach, about 8.7e-6 for 12i


_DIRECT = ["2009_12a", "2002_fig9", "2002_fig8", "2009_12c"]


@pytest.mark.parametrize("key", _DIRECT)
def test_single_integration_returns_after_one_period(key: str) -> None:
    """Literal closure: one integration from the printed state returns to it after n encounters.

    Only for the four least unstable orbits; the miss is the stability parameter times the
    initial error (1e-16 rounding plus about 1e-13 integration error per encounter), so the bound
    is |stability parameter| x 1e-11. The others close through the shooting fixed point above.
    """
    assert _direct_miss(key, "paper") < abs(ORBITS[key].stability) * 1e-11


@pytest.mark.parametrize(
    "reading", ["digest_literal", "no_mu_term", "clockwise", "psi_from_radial"]
)
@pytest.mark.parametrize("key", ["2009_12a", "2002_fig9"])
def test_wrong_readings_do_not_close(key: str, reading: str) -> None:
    """Control: each wrong convention misses the start by more than 1e-3 after one period."""
    assert _direct_miss(key, reading) > 1e-3
