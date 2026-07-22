"""Sun-Earth PCR3BP glue for the #685 own-system application of the #664
set-oriented transfer-operator (GAIO) pipeline (`#661` shortlist item 2).

This is the system-specific glue analogous to
``search/quasi_hilda_positive_control.py`` (the Sun-Jupiter Dellnitz-2005
positive control), but pointed for the FIRST time at one of THIS project's own
flagship systems -- Earth-Mars transport -- searching for a genuinely novel
almost-invariant transport SET rather than reproducing a published one. It is
built entirely on top of the generic, model-agnostic GAIO primitives in
``search/set_oriented_transfer_operator.py`` (``BoxGrid``,
``build_transition_matrix``, ``almost_invariant_sets_spectral``,
``almost_invariance_ratio``, ``transport_probability``) -- those are reused
verbatim, never rebuilt -- and reuses this project's own validated
``core.cr3bp`` DOP853 propagator end to end (no new integrator).

Why Sun-Earth PCR3BP (and NOT a literal 3-body "Sun-Earth-Mars")
----------------------------------------------------------------
The CR3BP admits exactly ONE secondary. The Dellnitz-2005 control's own
approach -- reused here -- is to take a two-body CR3BP (there Sun-Jupiter) and
define the "crossing" region for a THIRD reference body (there Mars) purely
GEOMETRICALLY, via osculating heliocentric two-body elements computed from the
state; Mars is never a body in the model, only a semimajor-axis reference
line. This module does the same, one system in: Sun + Earth as the CR3BP
primaries, with Mars entering ONLY as the osculating-aphelion reference
``a_Mars`` that defines the "Mars-reaching" region ``Q``.

Sun-Earth (rather than Sun-Mars) is chosen deliberately: Earth is the dominant
perturber for Earth-Mars cyclers and the operationally important encounters
happen at Earth; the length unit is Earth's own SMA, so nondimensional ``x``
is (to a part in 3e-6) the heliocentric distance in AU, and Mars sits at a
clean ``a_Mars ~= 1.5237`` length units. This is a genuinely DIFFERENT-strength
question from `#681`'s recent Sun-Mars WSB clean negative: `#681` searched for a
STRICT, repeating BALLISTIC-CAPTURE chain (a periodic capture<->escape claim);
this asks the fundamentally weaker, PROBABILISTIC question of whether an
almost-invariant transport SET of Earth-neighborhood transfer orbits leaks,
with measurable residence/transport statistics, into Mars-reaching orbits --
not redundant with that negative.

Model / section / energy
------------------------
Planar CR3BP, Sun + Earth primaries (``core.cr3bp.cr3bp_system("Sun",
"Earth")``, ``mu ~= 3.0e-6``). Poincare section: ``y = 0``, ``ydot < 0``,
``x > 0`` -- the EXTERIOR branch (spacecraft on Earth's side of the Sun,
heliocentric distance ``~= x``), the natural analog of Dellnitz's own
interior-realm ``x < 0`` branch, flipped outward because Earth-Mars transport
lives OUTSIDE Earth's orbit. Energy is fixed at the Jacobi constant of the
minimum-energy Earth-Mars (Hohmann) transfer ellipse (perihelion at Earth's
SMA, aphelion at Mars's SMA) -- :data:`SUN_EARTH_MARS_TRANSFER_C`, computed
below. This energy sits just BELOW ``C_L2`` for Sun-Earth (the L2 neck is
OPEN), which is deliberate and the OPPOSITE of Dellnitz's closed-neck choice:
Earth-Mars transport REQUIRES the exterior realm to communicate with Earth's
vicinity, so the physically relevant regime is exactly the open-neck one. Any
almost-invariant set found here is therefore a genuine metastable trap in a
LEAKY domain, not an artifact of an artificially sealed energy well.

The regions R and Q
-------------------
Both are defined geometrically on the section from osculating heliocentric
two-body elements (:func:`osculating_elements_at_section`), exactly as the
Dellnitz control defined its Mars-crossing line, and are constructed
DISJOINT (verified directly, see the tests) so any measured R->Q transport is
a genuine dynamical effect, not geometric overlap:

* ``Q`` (Mars-reaching): a bound (``e < 1``) osculating orbit whose aphelion
  ``r_a = a(1+e)`` sits in the Mars-encounter annulus
  (``a_Mars <= r_a <= 1.70``; see :data:`MARS_REACHING_R_A_MAX`). The outward
  analog of Dellnitz's "osculating periapsis <= Mars SMA" Mars-crossing line,
  banded to Mars's actual orbit so it stays a COMPACT, cycler-meaningful
  "aphelion at Mars" region rather than an unbounded beyond-Mars half-plane
  dominated by near-parabolic escape orbits.
* ``R`` (Earth-neighborhood sub-Mars transfer region): a bound orbit confined
  to the Earth-to-sub-Mars annulus -- perihelion at or above Earth's orbit
  (``r_p >= 0.95``) AND aphelion strictly short of Mars (``r_a <= 1.45``, a
  deliberate gap below ``a_Mars ~= 1.5237`` for clean disjointness). These are
  the Earth-coupled transfer ellipses that a resonant-pumping cycler mechanism
  would have to lift across the Mars line to make a genuine Earth-Mars cycler.

R->Q transport thus measures precisely the Earth-Mars "pumping" question:
under the real Sun-Earth CR3BP dynamics, does an Earth-neighborhood transfer
orbit get pumped across the Mars-reaching line, and with what
residence/transport statistics? A large almost-invariant set coincident with
R plus negligible R->Q transport is a clean negative (Earth-neighborhood
orbits are metastably trapped short of Mars); a substantial almost-invariant
set with non-negligible transport would be a genuinely novel finding to report
for adjudication (NOT written to the catalogue here -- the open schema
question of how a metastable SET earns a catalogue row is inherited unresolved
from `#664`).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import solve_ivp

from cyclerfinder.core.cr3bp import CR3BPSystem, cr3bp_eom, cr3bp_system

# Earth / Mars heliocentric semimajor axes (au), Standish & Williams J2000 --
# duplicated as bare floats from ``core/constants.py``'s own ``_EARTH_SMA_AU`` /
# ``_MARS_SMA_AU`` (Sun-Earth is served from the planet registry for the CR3BP
# system itself, but these two floats keep this module's Mars-reference /
# length-unit API self-contained, exactly as ``quasi_hilda_positive_control``
# duplicates the Mars/Jupiter SMAs for the same reason).
_EARTH_SMA_AU = 1.00000261
_MARS_SMA_AU = 1.52371034


def sun_earth_system() -> CR3BPSystem:
    """Sun-Earth CR3BP system, built from this project's own registry."""
    return cr3bp_system("Sun", "Earth")


def a_mars_nondim() -> float:
    """Mars's heliocentric SMA in Sun-Earth CR3BP length units (Earth's own
    SMA about the Sun); ``~= 1.5237``."""
    return _MARS_SMA_AU / _EARTH_SMA_AU


def hohmann_transfer_elements() -> tuple[float, float]:
    """``(a, e)`` of the minimum-energy Earth->Mars Hohmann transfer ellipse in
    nondimensional (Earth-SMA) length units: perihelion at Earth's orbit
    (``r_p = 1``), aphelion at Mars's orbit (``r_a = a_Mars``)."""
    r_p = 1.0
    r_a = a_mars_nondim()
    a = 0.5 * (r_p + r_a)
    e = (r_a - r_p) / (r_a + r_p)
    return a, e


def _tisserand_jacobi(a: float, e: float) -> float:
    """Planar Tisserand-parameter approximation to the Jacobi constant of a
    heliocentric osculating orbit ``(a, e)`` (cos i = 1):
    ``C ~= 1/a + 2*sqrt(a(1-e^2))``. Exact in the ``mu -> 0`` limit and, at
    Sun-Earth's ``mu ~= 3e-6``, within the section's own energy-shell
    thickness -- the section state is subsequently built to reproduce THIS
    ``C`` exactly under ``core.cr3bp.jacobi_constant`` (see the tests), so the
    approximation only sets which energy shell we work on, not the shell's
    internal consistency."""
    return 1.0 / a + 2.0 * math.sqrt(a * (1.0 - e * e))


#: Jacobi constant of the minimum-energy Earth-Mars Hohmann transfer ellipse --
#: the physically-motivated energy shell this task works on (see module
#: docstring). Computed, not hand-typed, from the Earth/Mars SMAs. ``~= 2.990``,
#: which sits just below ``C_L2`` for Sun-Earth (open L2 neck -- the transport
#: regime).
SUN_EARTH_MARS_TRANSFER_C = _tisserand_jacobi(*hohmann_transfer_elements())


def zero_velocity_v(x: float, mu: float, c_target: float) -> float:
    """``V(x) = x^2 + 2(1-mu)/r1 + 2mu/r2 - C`` at the section (``y = 0``).

    The accessible ``xdot`` range at this ``x`` is ``|xdot| <= sqrt(V(x))``
    when ``V(x) >= 0``; the point is inaccessible (imaginary ``ydot``) when
    ``V(x) < 0``. Identical in form to the Sun-Jupiter control's own
    ``zero_velocity_v`` -- the CR3BP zero-velocity relation is system-agnostic.
    """
    r1 = abs(x + mu)
    r2 = abs(x - 1.0 + mu)
    return x * x + 2.0 * (1.0 - mu) / r1 + 2.0 * mu / r2 - c_target


def accessible_xdot_max(x: float, mu: float, c_target: float) -> float | None:
    """``sqrt(V(x))``, or ``None`` if ``x`` is outside the accessible
    zero-velocity band at this energy."""
    v = zero_velocity_v(x, mu, c_target)
    if v < 0.0:
        return None
    return math.sqrt(v)


def section_state6(x: float, xdot: float, mu: float, c_target: float) -> NDArray[np.float64] | None:
    """Full CR3BP rotating-frame 6-state at section point ``(x, xdot)``, or
    ``None`` if ``(x, xdot)`` is off the energy manifold (``ydot^2 < 0``).
    ``ydot`` is taken NEGATIVE (the section branch, ``ydot < 0``)."""
    v = zero_velocity_v(x, mu, c_target) - xdot * xdot
    if v < 0.0:
        return None
    ydot = -math.sqrt(v)
    return np.array([x, 0.0, 0.0, xdot, ydot, 0.0], dtype=np.float64)


def _y_zero_event(t: float, y: NDArray[np.float64], mu: float) -> float:
    return float(y[1])


_y_zero_event.terminal = True  # type: ignore[attr-defined]
_y_zero_event.direction = -1.0  # type: ignore[attr-defined]  # only ydot<0 crossings


def poincare_first_return_exterior(
    state6: NDArray[np.float64],
    mu: float,
    *,
    t_segment: float = 2.0 * math.pi,
    max_segments: int = 40,
    rtol: float = 1e-10,
    atol: float = 1e-10,
    kickoff_dt: float = 1e-6,
) -> NDArray[np.float64] | None:
    """Propagate ``state6`` forward to the first ``y=0, ydot<0, x>0`` return.

    This is the EXTERIOR-branch (``x > 0``) analog of the Sun-Jupiter control's
    own ``poincare_first_return`` (which accepts ``x < 0``); the only
    difference is the accepted-crossing side. It carries the SAME #664
    section-crossing ``t=0`` self-detection fix VERBATIM in intent: every state
    passed in sits EXACTLY on the section (``y = 0``) by construction, so a
    naively event-monitored integration would spuriously self-detect the START
    point as its own return (the very first bracket ``[g(0)=0, g(h)<0]`` reads
    as a valid ``ydot<0`` crossing at ``t=0``). The fix is identical: advance a
    tiny FIXED, event-free ``kickoff_dt`` step first (four+ orders of magnitude
    below the shortest osculating return time reachable in this domain, which
    is ``~2*pi`` for a near-Earth orbit), landing strictly off the section,
    then hand off to the event-monitored segment.

    A crossing with ``x <= 0`` (the interior branch, on the far side of the
    Sun) is not accepted -- the search continues from that crossing into a
    fresh kick-off-then-monitor segment. Returns ``None`` if no qualifying
    ``x > 0`` crossing is found within ``max_segments`` segments, or if the
    integrator fails.
    """
    y = np.asarray(state6, dtype=np.float64).copy()
    for _ in range(max_segments):
        kick = solve_ivp(
            cr3bp_eom,
            (0.0, kickoff_dt),
            y,
            args=(mu,),
            method="DOP853",
            rtol=1e-13,
            atol=1e-13,
        )
        if not kick.success:
            return None
        y_kicked = kick.y[:, -1]
        sol = solve_ivp(
            cr3bp_eom,
            (0.0, t_segment),
            y_kicked,
            args=(mu,),  # type: ignore[call-overload]
            method="DOP853",
            rtol=rtol,
            atol=atol,
            events=[_y_zero_event],
            dense_output=False,
        )
        if not sol.success:
            return None
        if sol.t_events[0].size == 0:
            y = sol.y[:, -1]
            continue
        candidate = sol.y_events[0][0]
        if candidate[0] > 0.0:
            return np.asarray(candidate, dtype=np.float64)
        y = np.asarray(candidate, dtype=np.float64)
    return None


def section_map_xxdot(
    x: float,
    xdot: float,
    mu: float,
    c_target: float,
    **kwargs: float,
) -> NDArray[np.float64] | None:
    """The exterior-branch Poincare return map on section coordinates
    ``(x, xdot)``.

    Returns the image ``(x_new, xdot_new)`` (shape ``(2,)``), or ``None`` if
    ``(x, xdot)`` is off the energy manifold or the propagation fails to find a
    qualifying ``x>0`` return. Suitable as the ``map_fn`` for
    :func:`cyclerfinder.search.set_oriented_transfer_operator.build_transition_matrix`
    when partially applied over ``(mu, c_target)``.
    """
    state0 = section_state6(x, xdot, mu, c_target)
    if state0 is None:
        return None
    t_segment = float(kwargs.get("t_segment", 2.0 * math.pi))
    max_segments = int(kwargs.get("max_segments", 40))
    result = poincare_first_return_exterior(
        state0, mu, t_segment=t_segment, max_segments=max_segments
    )
    if result is None:
        return None
    return np.array([result[0], result[3]], dtype=np.float64)


@dataclass(frozen=True)
class OsculatingElements:
    a: float
    e: float
    r_p: float
    r_a: float


def osculating_elements_at_section(
    x: float, xdot: float, mu: float, c_target: float
) -> OsculatingElements | None:
    """Osculating heliocentric two-body elements ``(a, e, r_p, r_a)`` at
    section point ``(x, xdot)`` -- the general-conic vis-viva/angular-momentum
    formulas applied to the instantaneous heliocentric position/velocity,
    exactly as the Sun-Jupiter control does it (this is the paper's own
    Mars-crossing-line definition, evaluated directly from the state).

    Sun sits at ``(-mu, 0)`` in the rotating frame, so ``r_vec = (x + mu, y)``;
    heliocentric velocity uses the rotating-to-inertial relation (frame rate 1,
    aligned at ``t=0``): ``vx_inertial = xdot - y``, ``vy_inertial = ydot + x``.
    For a bound orbit (``e < 1``) ``r_a = a(1+e)``; for an unbound orbit
    (``e >= 1``, no aphelion) ``r_a = +inf``. Returns ``None`` if ``(x, xdot)``
    is off the energy manifold.
    """
    state = section_state6(x, xdot, mu, c_target)
    if state is None:
        return None
    _, y, _, xdot_, ydot, _ = (float(v) for v in state)
    gm = 1.0 - mu
    r_vec = np.array([x + mu, y])
    v_vec = np.array([xdot_ - y, ydot + x])
    h = r_vec[0] * v_vec[1] - r_vec[1] * v_vec[0]
    r = float(np.linalg.norm(r_vec))
    v2 = float(v_vec @ v_vec)
    eps = v2 / 2.0 - gm / r
    p_slr = h * h / gm
    ecc_sq = 1.0 + 2.0 * eps * h * h / (gm * gm)
    ecc = math.sqrt(max(0.0, ecc_sq))
    r_p = p_slr / (1.0 + ecc)
    a = -gm / (2.0 * eps) if abs(eps) > 1e-13 else math.inf
    r_a = a * (1.0 + ecc) if ecc < 1.0 and math.isfinite(a) else math.inf
    return OsculatingElements(a=a, e=ecc, r_p=r_p, r_a=r_a)


#: Upper aphelion bound (nondim, Earth-SMA units) for the Mars-reaching region
#: Q. Mars's own aphelion is ``a_Mars(1+e_Mars) ~= 1.524*1.093 ~= 1.666``; the
#: cap ``1.70`` admits orbits whose aphelion sits anywhere on Mars's actual
#: orbit (perihelion to aphelion) while EXCLUDING near-parabolic escape orbits
#: whose aphelion runs to many AU (a direct scan of this energy shell shows the
#: unbanded ``r_a >= a_Mars`` set is dominated by such near-escape orbits, with
#: median ``r_a ~= 1.8`` and a long tail past 100 -- physically these are not
#: Earth-Mars transfer/cycler orbits, whose aphelion is AT Mars). Banding Q
#: this way keeps it a COMPACT, cycler-meaningful target region, much closer to
#: Dellnitz's own compact Mars-crossing region than an "everything beyond Mars"
#: half-plane would be.
MARS_REACHING_R_A_MAX = 1.70


def mars_reaching_indicator(
    x: float,
    xdot: float,
    mu: float,
    c_target: float,
    *,
    r_a_max: float = MARS_REACHING_R_A_MAX,
) -> bool | None:
    """``True`` iff the osculating orbit at this section point is BOUND
    (``e < 1``) with its aphelion in the Mars-encounter annulus
    (``a_Mars <= r_a <= r_a_max``, default ``r_a_max = 1.70``) -- the outward
    analog of the Sun-Jupiter control's Mars-crossing line, banded to Mars's
    actual orbit (see :data:`MARS_REACHING_R_A_MAX`) so it is a compact,
    cycler-meaningful "aphelion at Mars" region rather than an unbounded
    beyond-Mars half-plane. ``None`` if ``(x, xdot)`` is off the energy
    manifold.
    """
    elems = osculating_elements_at_section(x, xdot, mu, c_target)
    if elems is None:
        return None
    return bool(elems.e < 1.0 and a_mars_nondim() <= elems.r_a <= r_a_max)


def earth_neighborhood_region_indicator(
    x: float,
    xdot: float,
    mu: float,
    c_target: float,
    *,
    r_p_min: float = 0.95,
    r_a_max: float = 1.45,
) -> bool | None:
    """``True`` iff the osculating orbit at this section point is a BOUND
    Earth-to-sub-Mars transfer ellipse: perihelion at/above Earth's orbit
    (``r_p >= r_p_min``, default 0.95) AND aphelion strictly short of Mars
    (``r_a <= r_a_max``, default 1.45 -- a deliberate gap below
    ``a_Mars ~= 1.5237`` so this region is geometrically DISJOINT from
    :func:`mars_reaching_indicator`, verified directly in the tests).

    These are the Earth-coupled transfer orbits confined to the
    Earth-to-sub-Mars annulus -- the candidate reservoir a resonant-pumping
    Earth-Mars cycler mechanism would have to lift across the Mars line.
    ``None`` if ``(x, xdot)`` is off the energy manifold.
    """
    elems = osculating_elements_at_section(x, xdot, mu, c_target)
    if elems is None:
        return None
    return bool(elems.e < 1.0 and elems.r_p >= r_p_min and elems.r_a <= r_a_max)
