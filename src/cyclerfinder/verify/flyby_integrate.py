"""Independent numerical confirmation of the Oberth periapsis powered-flyby cost.

Task #154. This module is the *independent check* (option a) for the
Oberth-credited periapsis cost model
:func:`cyclerfinder.core.flyby.dv_powered_flyby_periapsis` introduced under
task #151 and re-cost in ``docs/notes/2026-06-07-oberth-flyby-recost.md``.

Method (deliberately different from the formula under test, same physics)
------------------------------------------------------------------------
The formula computes the periapsis maneuver analytically from the hyperbolic
bend relation ``sin(delta/2) = 1/(1 + rp * vinf**2 / mu)``. Here we instead
**numerically integrate** the two-body equations of motion (the flyby body's
``mu`` only) with :func:`scipy.integrate.solve_ivp`, extract the realised
turn angle by propagating far enough up- and down-stream that the velocity
direction has converged to its asymptote, and **root-solve** (1-D, over the
periapsis speed) for the maneuver that delivers the required turn. The Δv we
report is then a function only of integrated trajectories and a Brent
root-find — it never calls :func:`dv_powered_flyby_periapsis`. Agreement
between the two is therefore an independent confirmation, in the same spirit
as the SPICE ephemeris cross-check: a second method, identical physics, no
shared formula.

The physical model being confirmed (note §"The two models"):

* The residual turn is supplied at periapsis. To make the *ballistic* cone
  open to the full required turn, the spacecraft must run the hyperbola at a
  lower periapsis speed ``vp_target`` (equivalently a lower excess speed). A
  slower hyperbola bends more.
* The closure must preserve ``|V_inf|``: slow into the widened cone on the
  inbound side, restore the original speed on the outbound side. Two
  tangential periapsis impulses, each of magnitude ``|vp_in - vp_target|``,
  hence ``ΔV = 2 |vp_in - vp_target|``.

So the independent task is: find, by integration, the periapsis speed
``vp_target_integrated`` whose hyperbola turns by exactly ``delta_required``;
then the integrated cost is ``2 |vp_in - vp_target_integrated|`` with
``vp_in = sqrt(vinf**2 + 2 mu / rp)`` (also confirmed by integrating the
unpowered incoming hyperbola and reading its turn == ``delta_max``).

Harness self-validation (note §3 of the task)
--------------------------------------------
:func:`integrate_turn_angle` integrating an *unpowered* hyperbola at the
incoming speed must reproduce :func:`cyclerfinder.core.flyby.max_bend` (the
analytic cone), and the asymptote-rotation baseline
:func:`cyclerfinder.core.flyby.dv_from_turn_deficit` is recovered trivially
as ``2 V_inf sin(deficit/2)`` (a rotation of the asymptote vector at
infinity, no integration needed). Both validate that the integration frame /
asymptote-extraction is correct before it is used to judge the Oberth model.

All quantities in km, km/s, km**3/s**2, radians — matching ``core/flyby.py``.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import cos, sin, sqrt
from typing import Final

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from cyclerfinder.core.flyby import max_bend

Vec3 = NDArray[np.float64]

# Far-field radius (in units of rp) at which the velocity direction is taken
# as the asymptote. At r = N * rp the speed exceeds vinf by a factor
# sqrt(1 + (2 mu / rp) / (vinf**2 * N)); the *direction* error vs the true
# asymptote falls off like 1/r. N = 1e5 keeps the asymptote-direction error
# below ~1e-6 rad for the Earth flybys swept here (documented + asserted by
# the harness-self-validation test against max_bend).
_FAR_FIELD_RP_MULTIPLE: Final[float] = 1.0e5


@dataclass(frozen=True)
class IntegratedFlyby:
    """Result of the independent integration for one (vinf, delta_req) point.

    Attributes
    ----------
    dv_integrated:
        Δv (km/s) the integration + root-solve find for the equal-|V_inf|
        periapsis closure: ``2 |vp_in - vp_target_integrated|``.
    dv_formula:
        :func:`cyclerfinder.core.flyby.dv_powered_flyby_periapsis` prediction
        for the identical geometry (the value under test).
    turn_achieved_rad:
        Turn angle the integrated powered hyperbola actually realises, read at
        the far field. Should equal ``delta_required`` to the root-solve
        tolerance.
    turn_required_rad:
        The ``delta_required`` target.
    vp_in:
        Incoming periapsis speed ``sqrt(vinf**2 + 2 mu / rp)`` (km/s).
    vp_target_integrated:
        Periapsis speed whose integrated hyperbola turns by ``delta_required``
        (km/s), found by the root-solve.
    rel_dv_error:
        ``|dv_integrated - dv_formula| / dv_formula`` (0.0 if formula is 0).
    asymptote_residual_rad:
        Residual of the unpowered incoming hyperbola's integrated turn vs the
        analytic ``max_bend`` cone — the documented asymptote-extraction error.
    """

    dv_integrated: float
    dv_formula: float
    turn_achieved_rad: float
    turn_required_rad: float
    vp_in: float
    vp_target_integrated: float
    rel_dv_error: float
    asymptote_residual_rad: float


def integrate_turn_angle(
    vp: float,
    rp: float,
    mu: float,
    far_field_rp_multiple: float = _FAR_FIELD_RP_MULTIPLE,
) -> float:
    """Numerically integrate a periapsis state and return the realised turn.

    The hyperbola is started at periapsis on the +x axis with the velocity
    purely tangential (+y), speed ``vp``. We propagate *outbound* until the
    radius reaches ``far_field_rp_multiple * rp`` and read the outgoing
    velocity direction; the incoming asymptote is the mirror of the outgoing
    one across the periapsis (symmetry of the unperturbed two-body
    hyperbola), so the total asymptote-to-asymptote turn is
    ``pi - 2 * angle(periapsis_velocity, outgoing_asymptote_velocity)`` ...
    but we avoid relying on that identity and instead integrate the *full*
    arc: backward to far-field (incoming asymptote) and forward to far-field
    (outgoing asymptote), then take the angle between the two velocity
    vectors. This makes no symmetry assumption and so also self-checks the
    integrator.

    Parameters
    ----------
    vp:
        Speed at periapsis, km/s. Must give a hyperbola (``vp**2 > 2 mu / rp``).
    rp:
        Periapsis radius, km.
    mu:
        Gravitational parameter, km**3/s**2.
    far_field_rp_multiple:
        Radius (in units of ``rp``) at which the velocity direction is read as
        the asymptote.

    Returns
    -------
    float
        Asymptote-to-asymptote turn angle, radians, in ``[0, pi]``.
    """
    if vp * vp <= 2.0 * mu / rp:
        raise ValueError("vp does not give a hyperbolic (escape) trajectory")

    r_far = far_field_rp_multiple * rp
    # Periapsis state: position +x, velocity +y (counter-clockwise hyperbola).
    y0 = np.array([rp, 0.0, 0.0, vp])

    # Conservative time span: well beyond the time to reach r_far at ~vinf.
    vinf = sqrt(vp * vp - 2.0 * mu / rp)
    t_max = 10.0 * r_far / max(vinf, 1.0e-3)

    def _rhs(_t: float, y: NDArray[np.float64]) -> NDArray[np.float64]:
        x, yy, vx, vy = y
        r = sqrt(x * x + yy * yy)
        a = -mu / (r * r * r)
        return np.array([vx, vy, a * x, a * yy])

    def _hit_far(_t: float, y: NDArray[np.float64]) -> float:
        return sqrt(y[0] * y[0] + y[1] * y[1]) - r_far

    _hit_far.terminal = True  # type: ignore[attr-defined]
    _hit_far.direction = 1.0  # type: ignore[attr-defined]

    out = solve_ivp(
        _rhs,
        (0.0, t_max),
        y0,
        method="DOP853",
        events=_hit_far,
        rtol=1.0e-12,
        atol=1.0e-9,
    )
    back = solve_ivp(
        _rhs,
        (0.0, -t_max),
        y0,
        method="DOP853",
        events=_hit_far,
        rtol=1.0e-12,
        atol=1.0e-9,
    )
    if out.t_events is None or len(out.t_events[0]) == 0:
        raise RuntimeError("outbound integration never reached far field")
    if back.t_events is None or len(back.t_events[0]) == 0:
        raise RuntimeError("inbound integration never reached far field")

    v_out = out.y[2:4, -1]
    v_in = back.y[2:4, -1]
    # ``solve_ivp`` stores the *instantaneous* velocity at each sample, so the
    # vector at the far point reached by backward integration is already the
    # inbound direction of travel (motion approaching the body) — no sign flip.
    # The deflection (turn) angle is the angle between this inbound asymptote
    # velocity and the outbound asymptote velocity. Validated against
    # :func:`max_bend` to 1e-8 rad by the harness-self-validation test.
    cos_turn = float(np.dot(v_in, v_out) / (np.linalg.norm(v_in) * np.linalg.norm(v_out)))
    return float(np.arccos(min(1.0, max(-1.0, cos_turn))))


def asymptote_baseline_dv(vinf: float, delta_required: float, delta_max: float) -> float:
    """Independent reconstruction of the asymptote-rotation baseline.

    The baseline maneuver rotates the ``V_inf`` *asymptote vector at infinity*
    by the deficit ``delta_required - delta_max`` while preserving its
    magnitude ``vinf``. The Δv of rotating a vector of length ``vinf`` by an
    angle ``phi`` (chord length) is exactly ``2 vinf sin(phi/2)``. We build it
    here from the vector geometry (two explicit ``vinf`` vectors and a
    subtraction) rather than calling
    :func:`cyclerfinder.core.flyby.dv_from_turn_deficit`, so the test that the
    two agree validates the harness's notion of "rotate the asymptote".

    Returns ``0.0`` when within the cone.
    """
    deficit = max(0.0, delta_required - delta_max)
    if deficit == 0.0:
        return 0.0
    v_before = np.array([vinf, 0.0])
    v_after = np.array([vinf * cos(deficit), vinf * sin(deficit)])
    return float(np.linalg.norm(v_after - v_before))


def confirm_oberth_point(
    vinf: float,
    delta_required: float,
    mu: float,
    rp: float,
    far_field_rp_multiple: float = _FAR_FIELD_RP_MULTIPLE,
) -> IntegratedFlyby:
    """Independently confirm :func:`dv_powered_flyby_periapsis` at one point.

    Steps (all independent of the formula under test):

    1. ``vp_in = sqrt(vinf**2 + 2 mu / rp)`` and confirm, by integrating the
       *unpowered* incoming hyperbola, that its turn equals the analytic cone
       ``max_bend`` (records the residual — the asymptote-extraction error).
    2. Root-solve (Brent) over periapsis speed ``vp`` in
       ``(escape, vp_in]`` for the ``vp`` whose integrated hyperbola turns by
       ``delta_required``. A *slower* periapsis speed turns more, so a
       deficit (required > cone) is met by ``vp < vp_in``.
    3. Integrated cost = ``2 |vp_in - vp_target_integrated|`` (slow in,
       restore out), the equal-|V_inf| closure.
    4. Compare to :func:`dv_powered_flyby_periapsis` for the identical
       geometry.

    Returns
    -------
    IntegratedFlyby
    """
    delta_max = max_bend(mu, rp, vinf)

    # Lazy import: the formula under test, only used for the comparison field,
    # never inside the integrated path.
    from cyclerfinder.core.flyby import dv_powered_flyby_periapsis

    dv_formula = dv_powered_flyby_periapsis(vinf, delta_required, delta_max, mu, rp)

    vp_in = sqrt(vinf * vinf + 2.0 * mu / rp)

    # (1) self-validate: unpowered incoming hyperbola turn vs analytic cone.
    turn_cone = integrate_turn_angle(vp_in, rp, mu, far_field_rp_multiple)
    asymptote_residual = abs(turn_cone - delta_max)

    if delta_required <= delta_max:
        # Within the cone: no maneuver. Integrated cost is exactly 0.
        return IntegratedFlyby(
            dv_integrated=0.0,
            dv_formula=dv_formula,
            turn_achieved_rad=turn_cone,
            turn_required_rad=delta_required,
            vp_in=vp_in,
            vp_target_integrated=vp_in,
            rel_dv_error=0.0,
            asymptote_residual_rad=asymptote_residual,
        )

    # (2) root-solve for the periapsis speed whose turn == delta_required.
    v_escape = sqrt(2.0 * mu / rp)

    def _turn_deficit(vp: float) -> float:
        return integrate_turn_angle(vp, rp, mu, far_field_rp_multiple) - delta_required

    # Bracket: at vp_in the turn is delta_max < delta_required so deficit < 0;
    # just above escape the turn approaches pi so deficit > 0. Nudge off the
    # escape singularity for a well-posed bracket.
    lo = v_escape * (1.0 + 1.0e-6)
    hi = vp_in
    f_lo = _turn_deficit(lo)
    f_hi = _turn_deficit(hi)
    if f_lo * f_hi > 0.0:
        raise RuntimeError(
            f"no bracket for turn root-solve: f(lo)={f_lo}, f(hi)={f_hi} "
            f"(delta_required={delta_required}, delta_max={delta_max})"
        )

    vp_target = float(brentq(_turn_deficit, lo, hi, xtol=1.0e-9, rtol=1.0e-12))
    turn_achieved = integrate_turn_angle(vp_target, rp, mu, far_field_rp_multiple)

    dv_integrated = 2.0 * abs(vp_in - vp_target)
    rel_err = abs(dv_integrated - dv_formula) / dv_formula if dv_formula > 0.0 else 0.0

    return IntegratedFlyby(
        dv_integrated=dv_integrated,
        dv_formula=dv_formula,
        turn_achieved_rad=turn_achieved,
        turn_required_rad=delta_required,
        vp_in=vp_in,
        vp_target_integrated=vp_target,
        rel_dv_error=rel_err,
        asymptote_residual_rad=asymptote_residual,
    )
