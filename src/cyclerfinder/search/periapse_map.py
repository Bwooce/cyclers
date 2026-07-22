"""Periapsis Poincaré ("periapse") maps for the planet-moon planar CR3BP (task #683).

The one classical seedless discovery map this codebase genuinely lacked
(survey-verified in `#679`'s strategy pass: no periapsis/apsis Poincaré-map
machinery existed anywhere in ``search/`` before this module). Instead of
shooting for exact periodic orbits (this project's correctors) or measuring
statistical transport (``set_oriented_transfer_operator.py`` / ``#664``/``#685``'s
GAIO pipeline), a periapse map fixes the Jacobi constant, takes the surface of
section AT periapsis with respect to the *secondary* (the moon), and reads the
map's own geometric lobe structure -- capture / escape / impact / transit --
directly, with NO Monte-Carlo noise and NO shooting correction.

Method and conventions -- Davis & Howell 2011
---------------------------------------------
This module follows, and is validated against, the periapse-map lineage of

    Davis, D. C. & Howell, K. C., "Trajectory evolution in the multi-body
    problem with applications in the Saturnian system," *Acta Astronautica*
    69 (2011) 1038-1049, doi:10.1016/j.actaastro.2011.07.007

(itself building on Villac & Scheeres 2003, JGCD 26(2) 224-232, who first
introduced the periapsis Poincaré map, and Davis & Howell 2012, JGCD 35(1)).
The Davis-Howell 2011 paper is the one that applies the map DIRECTLY to the
Saturn-Titan system with explicit published Jacobi values, so it is the
positive control this module reproduces. See the corpus digest
``docs/notes/2026-07-22-683-digest-davis-howell-2011-periapse-maps.md``.

Key definitions taken verbatim from that paper (its Sec. 2.4, "Periapsis
Poincaré maps"):

* **Jacobi integral** (paper Eq. 2, identical to this project's
  :func:`cyclerfinder.core.cr3bp.jacobi_constant` convention -- verified: the
  paper's ``J = x^2 + y^2 + 2 mu / r + 2(1-mu)/d - v^2`` with ``r`` = distance
  to the secondary ``P2`` and ``d`` = distance to the primary ``P1`` is exactly
  ``core.cr3bp``'s ``C``).

* **Periapsis (the section)**: a state satisfying ``rdot = 0`` (radial rate
  with respect to ``P2``) with ``rddot > 0`` (a local *minimum* of the
  distance to the secondary). Outside the ``rdot = 0`` contour it is an
  apoapsis (``rddot < 0``).

* **Map coordinates**: the periapse *position* in configuration space, plotted
  relative to ``P2`` in units of the Hill radius ``r_H = (mu/3)^(1/3)``:
  ``x_p = (x - (1-mu)) / r_H``, ``y_p = y / r_H``. The periapse angle
  ``omega_r = atan2(y, x - (1-mu))`` and periapse radius ``r_p`` are the polar
  form of the same point (paper Fig. 1).

* **Parametrisation** (paper Sec. 2.2): for a chosen periapse position and a
  fixed Jacobi constant, the velocity MAGNITUDE is set by the Jacobi integral
  and the velocity DIRECTION is perpendicular to the ``P2``-relative radius
  (the apse condition). The prograde choice fixes the remaining sign. Hence
  each in-contour periapse position defines a single planar prograde
  trajectory.

* **Fate after propagation** (paper Sec. 3): propagate a periapse initial
  condition forward and classify the outcome -- ``IMPACT`` (passes on or within
  the radius of ``P2``), ``ESCAPE_L1`` / ``ESCAPE_L2`` (an ``x``-coordinate
  lying more than ``escape_margin = 0.01`` nondimensional units beyond ``L1`` or
  ``L2`` respectively), or ``CAPTURED`` (remains near ``P2`` for the whole
  propagation). This is the paper's "initial condition map" (its Figs. 3, 4).

Propagation is this project's own validated ``core.cr3bp`` DOP853 propagator
(``cr3bp_eom``) end to end -- no new integrator, matching every other search
module here.
"""

from __future__ import annotations

import math
from collections.abc import Callable
from dataclasses import dataclass
from enum import Enum

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from cyclerfinder.core.constants import PLANETS
from cyclerfinder.core.cr3bp import (
    CR3BPSystem,
    cr3bp_eom,
    cr3bp_system,
    jacobi_constant,
)
from cyclerfinder.core.satellites import SATELLITES

_EventFn = Callable[[float, NDArray[np.float64], float], float]


class PeriapseFate(Enum):
    """Outcome of propagating a periapse initial condition (Davis-Howell 2011,
    Sec. 3): the four regions of the paper's initial-condition map."""

    IMPACT = "impact"
    ESCAPE_L1 = "escape_l1"
    ESCAPE_L2 = "escape_l2"
    CAPTURED = "captured"


@dataclass(frozen=True)
class PeriapseMapSystem:
    """A planet-moon CR3BP prepared for periapse-map work.

    Bundles the validated :class:`CR3BPSystem` with the Hill radius, the two
    collinear libration points bracketing the secondary (``L1`` interior,
    ``L2`` exterior) and their Jacobi constants, and the secondary's physical
    radius in nondimensional length units (for the impact test).
    """

    system: CR3BPSystem
    r_hill: float
    x_l1: float  # rotating-frame x of L1 (interior collinear point)
    x_l2: float  # rotating-frame x of L2 (exterior collinear point)
    c_l1: float  # Jacobi constant at L1
    c_l2: float  # Jacobi constant at L2
    secondary_radius_nd: float  # secondary body radius / l_km

    @property
    def mu(self) -> float:
        return self.system.mu

    @property
    def x_secondary(self) -> float:
        return 1.0 - self.system.mu


def _dudx_on_axis(x: float, mu: float) -> float:
    """d(effective potential)/dx along the x-axis (y=z=0). Its zeros are the
    collinear libration points."""
    r1 = abs(x + mu)
    r2 = abs(x - 1.0 + mu)
    return x - (1.0 - mu) * (x + mu) / r1**3 - mu * (x - 1.0 + mu) / r2**3


def collinear_l1_l2(mu: float) -> tuple[float, float]:
    """Rotating-frame x-locations of the interior (``L1``) and exterior
    (``L2``) collinear libration points bracketing the secondary at
    ``x = 1-mu``. Located as the two zeros of :func:`_dudx_on_axis` on either
    side of the secondary."""
    x_sec = 1.0 - mu
    x_l1 = brentq(_dudx_on_axis, x_sec - 0.5, x_sec - 1e-9, args=(mu,))
    x_l2 = brentq(_dudx_on_axis, x_sec + 1e-9, x_sec + 0.5, args=(mu,))
    return float(x_l1), float(x_l2)


def _jacobi_on_axis(x: float, mu: float) -> float:
    """Jacobi constant of the (zero-velocity) state at ``(x, 0)``."""
    return jacobi_constant(np.array([x, 0.0, 0.0, 0.0, 0.0, 0.0]), mu)


def _secondary_radius_km(primary: str, secondary: str) -> float:
    """Physical equatorial radius (km) of the secondary body, resolving from
    the moon registry (moon secondaries) or the planet registry (Sun-primary
    planet secondaries, e.g. Sun-Saturn)."""
    if primary == "Sun":
        for planet in PLANETS.values():
            if planet.name == secondary:
                return planet.radius_eq_km
        raise KeyError(f"no Sun-primary planet named {secondary!r}")
    return SATELLITES[secondary].radius_eq_km


def build_periapse_map_system(primary: str, secondary: str) -> PeriapseMapSystem:
    """Assemble a :class:`PeriapseMapSystem` from this project's registries.

    Reuses :func:`cyclerfinder.core.cr3bp.cr3bp_system` for ``mu`` and the
    length/time scales, and the moon registry for the secondary's physical
    radius (for the impact test)."""
    system = cr3bp_system(primary, secondary)
    mu = system.mu
    r_hill = (mu / 3.0) ** (1.0 / 3.0)
    x_l1, x_l2 = collinear_l1_l2(mu)
    c_l1 = _jacobi_on_axis(x_l1, mu)
    c_l2 = _jacobi_on_axis(x_l2, mu)
    secondary_radius_nd = _secondary_radius_km(primary, secondary) / system.l_km
    return PeriapseMapSystem(
        system=system,
        r_hill=r_hill,
        x_l1=x_l1,
        x_l2=x_l2,
        c_l1=c_l1,
        c_l2=c_l2,
        secondary_radius_nd=secondary_radius_nd,
    )


def secondary_radial_rate(state6: NDArray[np.float64], mu: float) -> float:
    """Radial rate ``rdot`` of the spacecraft with respect to the secondary
    ``P2`` (fixed at ``(1-mu, 0)`` in the rotating frame). Zero at an apse."""
    dx = float(state6[0]) - (1.0 - mu)
    dy = float(state6[1])
    r2 = math.hypot(dx, dy)
    if r2 == 0.0:
        return 0.0
    return (dx * float(state6[3]) + dy * float(state6[4])) / r2


def secondary_radial_accel(state6: NDArray[np.float64], mu: float) -> float:
    """Radial acceleration ``rddot`` of the spacecraft w.r.t. the secondary.

    ``rddot = (|v_rel|^2 + rel . a - rdot^2) / r2`` where ``a`` is the full
    rotating-frame acceleration (``cr3bp_eom``) -- the secondary is stationary
    in the rotating frame, so the relative acceleration IS ``a``. At an apse
    (``rdot = 0``) this reduces to ``(|v_rel|^2 + rel . a) / r2``; ``rddot > 0``
    marks a periapsis (local minimum of ``r2``), ``rddot < 0`` an apoapsis."""
    dx = float(state6[0]) - (1.0 - mu)
    dy = float(state6[1])
    r2 = math.hypot(dx, dy)
    if r2 == 0.0:
        return math.inf
    vx, vy = float(state6[3]), float(state6[4])
    acc = cr3bp_eom(0.0, state6, mu)
    ax, ay = float(acc[3]), float(acc[4])
    rdot = (dx * vx + dy * vy) / r2
    v_rel_sq = vx * vx + vy * vy
    rel_dot_a = dx * ax + dy * ay
    return (v_rel_sq + rel_dot_a - rdot * rdot) / r2


def is_periapsis(state6: NDArray[np.float64], mu: float) -> bool:
    """True iff ``state6`` is at (or effectively at) a periapsis w.r.t. the
    secondary: ``rddot > 0`` (a local minimum of distance to ``P2``)."""
    return secondary_radial_accel(state6, mu) > 0.0


def construct_periapse_state(
    x: float,
    y: float,
    jacobi: float,
    mu: float,
    *,
    prograde: bool = True,
) -> NDArray[np.float64] | None:
    """Full planar CR3BP 6-state at a periapse position ``(x, y)`` and fixed
    Jacobi constant (Davis-Howell 2011, Sec. 2.2 parametrisation).

    The velocity magnitude ``v`` is set by the Jacobi integral
    (``v^2 = x^2 + y^2 + 2(1-mu)/r1 + 2 mu/r2 - J``); the direction is
    perpendicular to the ``P2``-relative radius (the apse condition). The
    ``prograde`` sense (counter-clockwise angular motion about ``P2`` in the
    rotating frame) fixes the remaining sign: ``v_hat = (-dy, dx)/r2``.

    Returns ``None`` if the point is outside the zero-velocity contour at this
    energy (``v^2 < 0`` -- no real velocity solves the Jacobi constraint there).
    A returned state is guaranteed on the section (``rdot = 0``) but is NOT
    guaranteed to be a periapsis rather than an apoapsis; call
    :func:`is_periapsis` to keep only the in-contour (``rddot > 0``) points, as
    the paper does.
    """
    r1 = math.hypot(x + mu, y)
    dx = x - (1.0 - mu)
    dy = y
    r2 = math.hypot(dx, dy)
    if r1 == 0.0 or r2 == 0.0:
        return None
    v_sq = x * x + y * y + 2.0 * (1.0 - mu) / r1 + 2.0 * mu / r2 - jacobi
    if v_sq < 0.0:
        return None
    v = math.sqrt(v_sq)
    # Prograde (ccw about P2): velocity perpendicular to the P2-relative radius.
    vx = v * (-dy / r2)
    vy = v * (dx / r2)
    if not prograde:
        vx, vy = -vx, -vy
    return np.array([x, y, 0.0, vx, vy, 0.0], dtype=np.float64)


def periapse_map_coords(
    state6: NDArray[np.float64], sys_or_mu: PeriapseMapSystem | float
) -> tuple[float, float, float, float]:
    """Map a periapse state to ``(x_p, y_p, r_p, omega_r)``.

    ``x_p, y_p`` are the periapse position relative to ``P2`` in Hill-radius
    units (the paper's Fig.-1 map axes); ``r_p`` is the periapse radius in
    nondimensional length units; ``omega_r = atan2(y, x-(1-mu))`` is the
    periapse angle (radians). Accepts either a :class:`PeriapseMapSystem`
    (so ``x_p, y_p`` are scaled by its Hill radius) or a bare ``mu`` (then
    ``x_p, y_p`` are returned in raw nondimensional units, i.e. ``r_H = 1``)."""
    if isinstance(sys_or_mu, PeriapseMapSystem):
        mu = sys_or_mu.mu
        r_h = sys_or_mu.r_hill
    else:
        mu = sys_or_mu
        r_h = 1.0
    dx = float(state6[0]) - (1.0 - mu)
    dy = float(state6[1])
    r_p = math.hypot(dx, dy)
    omega_r = math.atan2(dy, dx)
    return dx / r_h, dy / r_h, r_p, omega_r


def _periapsis_event(t: float, y: NDArray[np.float64], mu: float) -> float:
    """``rdot * r2`` (the numerator of the radial rate) w.r.t. the secondary.
    Zero at every apse; ``direction=+1`` selects the increasing crossing =
    periapsis (local minimum of distance to ``P2``)."""
    dx = float(y[0]) - (1.0 - mu)
    dy = float(y[1])
    return dx * float(y[3]) + dy * float(y[4])


_periapsis_event.terminal = True  # type: ignore[attr-defined]
_periapsis_event.direction = 1.0  # type: ignore[attr-defined]


def _make_impact_event(radius_nd: float) -> _EventFn:
    def event(t: float, y: NDArray[np.float64], mu: float) -> float:
        dx = float(y[0]) - (1.0 - mu)
        dy = float(y[1])
        return math.hypot(dx, dy) - radius_nd

    event.terminal = True  # type: ignore[attr-defined]
    event.direction = -1.0  # type: ignore[attr-defined]
    return event


def _make_escape_l1_event(x_l1: float, margin: float) -> _EventFn:
    def event(t: float, y: NDArray[np.float64], mu: float) -> float:
        return float(y[0]) - (x_l1 - margin)

    event.terminal = True  # type: ignore[attr-defined]
    event.direction = -1.0  # type: ignore[attr-defined]  # x decreasing past L1
    return event


def _make_escape_l2_event(x_l2: float, margin: float) -> _EventFn:
    def event(t: float, y: NDArray[np.float64], mu: float) -> float:
        return float(y[0]) - (x_l2 + margin)

    event.terminal = True  # type: ignore[attr-defined]
    event.direction = 1.0  # type: ignore[attr-defined]  # x increasing past L2
    return event


@dataclass(frozen=True)
class SegmentResult:
    """Outcome of one propagated segment between apse/impact/escape events."""

    kind: str  # "periapsis" | "impact" | "escape_l1" | "escape_l2" | "none"
    state: NDArray[np.float64] | None  # state at the firing event (None if kind=="none")
    t: float  # nondim time reached (event time, or t_max if kind=="none")


def _advance_to_event(
    state6: NDArray[np.float64],
    pmsys: PeriapseMapSystem,
    *,
    escape_margin: float,
    t_max: float,
    kickoff_dt: float,
    rtol: float,
    atol: float,
) -> SegmentResult:
    """Propagate from a periapse state to the FIRST of: the next periapsis,
    an impact, an L1 escape, or an L2 escape (whichever event fires earliest in
    time), or exhaust ``t_max``.

    A tiny fixed ``kickoff_dt`` step is taken WITHOUT event monitoring first so
    the starting periapsis (which sits exactly on the ``rdot=0`` section) is not
    spuriously self-detected as the return -- the same guard the sibling
    ``quasi_hilda_positive_control.poincare_first_return`` uses.
    """
    mu = pmsys.mu
    kick = solve_ivp(
        cr3bp_eom,
        (0.0, kickoff_dt),
        np.asarray(state6, dtype=np.float64),
        args=(mu,),
        method="DOP853",
        rtol=1e-13,
        atol=1e-13,
    )
    if not kick.success:
        return SegmentResult(kind="none", state=None, t=0.0)
    y0 = kick.y[:, -1]

    events: list[_EventFn] = [
        _periapsis_event,
        _make_impact_event(pmsys.secondary_radius_nd),
        _make_escape_l1_event(pmsys.x_l1, escape_margin),
        _make_escape_l2_event(pmsys.x_l2, escape_margin),
    ]
    kinds = ["periapsis", "impact", "escape_l1", "escape_l2"]
    sol = solve_ivp(
        cr3bp_eom,
        (0.0, t_max),
        y0,
        args=(mu,),  # type: ignore[call-overload]
        method="DOP853",
        rtol=rtol,
        atol=atol,
        events=events,
        dense_output=False,
    )
    if not sol.success:
        return SegmentResult(kind="none", state=None, t=float(sol.t[-1]) + kickoff_dt)

    # Find the earliest-firing event across all four channels.
    best_t = math.inf
    best_kind = "none"
    best_state: NDArray[np.float64] | None = None
    for k, kind in enumerate(kinds):
        te = sol.t_events[k]
        if te.size == 0:
            continue
        if float(te[0]) < best_t:
            best_t = float(te[0])
            best_kind = kind
            best_state = np.asarray(sol.y_events[k][0], dtype=np.float64)
    if best_kind == "none":
        return SegmentResult(kind="none", state=None, t=t_max + kickoff_dt)
    return SegmentResult(kind=best_kind, state=best_state, t=best_t + kickoff_dt)


def next_periapse(
    state6: NDArray[np.float64],
    pmsys: PeriapseMapSystem,
    *,
    t_max: float = 4.0 * math.pi,
    kickoff_dt: float = 1e-6,
    rtol: float = 1e-11,
    atol: float = 1e-11,
) -> NDArray[np.float64] | None:
    """Propagate to the next periapsis w.r.t. the secondary and return that
    state, IGNORING impact/escape (pure apsidal return map). Returns ``None``
    if no periapsis is reached within ``t_max``. Used by the long-term
    itinerary tracker, where impact/escape are handled by the caller."""
    mu = pmsys.mu
    kick = solve_ivp(
        cr3bp_eom,
        (0.0, kickoff_dt),
        np.asarray(state6, dtype=np.float64),
        args=(mu,),
        method="DOP853",
        rtol=1e-13,
        atol=1e-13,
    )
    if not kick.success:
        return None
    sol = solve_ivp(
        cr3bp_eom,
        (0.0, t_max),
        kick.y[:, -1],
        args=(mu,),  # type: ignore[call-overload]
        method="DOP853",
        rtol=rtol,
        atol=atol,
        events=[_periapsis_event],
        dense_output=False,
    )
    if not sol.success or sol.t_events[0].size == 0:
        return None
    return np.asarray(sol.y_events[0][0], dtype=np.float64)


def classify_fate(
    state6: NDArray[np.float64],
    pmsys: PeriapseMapSystem,
    *,
    max_revs: int = 1,
    escape_margin: float = 0.01,
    t_max_per_rev: float = 4.0 * math.pi,
    rtol: float = 1e-11,
    atol: float = 1e-11,
) -> tuple[PeriapseFate, list[NDArray[np.float64]]]:
    """Classify a periapse initial condition's fate over up to ``max_revs``
    revolutions (Davis-Howell 2011 initial-condition map, Figs. 3/4).

    Propagates periapse-to-periapse; returns as soon as an impact or escape
    event fires, otherwise ``CAPTURED`` after ``max_revs`` complete periapsis
    returns. Also returns the list of periapse states visited (each is exactly
    on the section). ``escape_margin`` (default 0.01, the paper's own value) is
    the nondimensional distance beyond ``L1``/``L2`` that defines escape."""
    periapses: list[NDArray[np.float64]] = []
    state = np.asarray(state6, dtype=np.float64)
    for _ in range(max_revs):
        seg = _advance_to_event(
            state,
            pmsys,
            escape_margin=escape_margin,
            t_max=t_max_per_rev,
            kickoff_dt=1e-6,
            rtol=rtol,
            atol=atol,
        )
        if seg.kind == "impact":
            return PeriapseFate.IMPACT, periapses
        if seg.kind == "escape_l1":
            return PeriapseFate.ESCAPE_L1, periapses
        if seg.kind == "escape_l2":
            return PeriapseFate.ESCAPE_L2, periapses
        if seg.kind == "none" or seg.state is None:
            # No further periapsis found within the time budget and no
            # impact/escape -- treat as captured (bounded near P2).
            return PeriapseFate.CAPTURED, periapses
        periapses.append(seg.state)
        state = seg.state
    return PeriapseFate.CAPTURED, periapses


def collect_titan_periapses(
    state6: NDArray[np.float64],
    pmsys: PeriapseMapSystem,
    *,
    t_total: float,
    rtol: float = 1e-10,
    atol: float = 1e-10,
) -> tuple[list[tuple[float, NDArray[np.float64]]], bool]:
    """Collect EVERY periapsis w.r.t. the secondary over a single continuous
    propagation of total nondimensional time ``t_total`` (one ``solve_ivp``
    call, so no per-revolution time-budget truncation).

    Returns ``(events, impacted)`` where ``events`` is the time-ordered list of
    ``(t, state)`` at each periapsis (radial minimum of distance to ``P2``),
    and ``impacted`` is ``True`` if the propagation terminated on an impact
    (closest approach within the secondary radius). Because the energy is fixed
    below the escape threshold the trajectory stays bound to the Saturn-Titan
    system, so no far-field cut-off is needed -- a trajectory that leaves
    Titan's vicinity through a gateway is tracked continuously (its closest
    approaches to Titan while it orbits the primary are recorded as periapses
    with a large ``r_p``), and a genuine re-capture reappears as later periapses
    with small ``r_p`` -- exactly the capture -> excursion -> re-capture signal
    the #683 search needs, with no truncation blind spot.
    """
    mu = pmsys.mu

    def peri_event(t: float, y: NDArray[np.float64], mu: float) -> float:
        dx = float(y[0]) - (1.0 - mu)
        dy = float(y[1])
        return dx * float(y[3]) + dy * float(y[4])

    peri_event.terminal = False  # type: ignore[attr-defined]
    peri_event.direction = 1.0  # type: ignore[attr-defined]

    impact_event = _make_impact_event(pmsys.secondary_radius_nd)
    sol = solve_ivp(
        cr3bp_eom,
        (0.0, t_total),
        np.asarray(state6, dtype=np.float64),
        args=(mu,),  # type: ignore[call-overload]
        method="DOP853",
        rtol=rtol,
        atol=atol,
        events=[peri_event, impact_event],
        dense_output=False,
    )
    if not sol.success:
        # Return whatever periapses were found before the failure.
        pass
    events: list[tuple[float, NDArray[np.float64]]] = []
    te = sol.t_events[0]
    ye = sol.y_events[0]
    for k in range(te.size):
        events.append((float(te[k]), np.asarray(ye[k], dtype=np.float64)))
    impacted = sol.t_events[1].size > 0
    return events, impacted
