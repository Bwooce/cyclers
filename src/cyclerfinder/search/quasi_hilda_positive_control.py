"""Sun-Jupiter PCR3BP glue for the Dellnitz et al. 2005 quasi-Hilda positive
control (task #664, `#661` shortlist item 2).

Reproduces the model/method of: Dellnitz, M., Junge, O., Lo, M. W., Marsden,
J. E., Padberg, K., Preis, R., Ross, S. D. and Thiere, B., "Transport of
Mars-Crossing Asteroids from the Quasi-Hilda Region," Phys. Rev. Lett. 94
(2005) 231102. See ``docs/notes/2026-06-30-digest-dellnitz2005-mars-transport-gaio.md``
for the full corpus digest this module is built against.

Model (paper Sec. 2, digest Sec. 2)
------------------------------------
Planar CR3BP, Sun + Jupiter primaries (this project's own
``core.cr3bp.cr3bp_system("Sun", "Jupiter")`` -- reused verbatim, no new
integrator). Poincare section: ``y = 0``, ``ydot < 0``, ``x < 0`` (the
"interior realm" branch -- conjunctions on the opposite side of the Sun from
Jupiter). Energy ``E = -1.52`` (just below ``E_L1 = -1.5199``), which the
paper chose specifically because at this level the L1 neck is CLOSED (no
direct interior<->Jupiter-realm connection) while Jupiter's perturbation
remains dynamically significant and the 3:2 Hilda resonance island is
present. The section coordinate ``xdot`` is the ROTATING-frame ``dx/dt``
(paper's own coordinate choice); ``ydot`` at each section point is SOLVED
from the Jacobi-constant energy constraint, not sampled independently.

Reconstructing the paper's exact box-covering domain and its "quasi-Hilda
island" shape honestly
------------------------------------------------------------------------------
The corpus digest (see module docstring reference above) captures the
paper's MODEL, METHOD, and sourced numerical RESULTS in full, but Fig. 1's
box-covering domain and the exact pixel-space shape of the 3:2 resonance
island are only described qualitatively ("sideways U-shaped region on the
left of the SOS") -- the digest has no extracted coordinate data for them,
and this module does not invent any. Instead, every geometric object used
below is reconstructed independently from first-principles orbital mechanics
at the SAME energy level, self-consistently verified against this project's
own :func:`cyclerfinder.core.cr3bp.jacobi_constant`:

* The box-covering domain is the accessible zero-velocity-curve lens on the
  section at ``x < 0`` closest to the origin (the true "interior realm,"
  ``a < a_Jupiter``) -- see :func:`accessible_xdot_max`. A direct numerical
  scan (see this task's own working notes) shows the ``x < 0`` half of the
  energy manifold at ``E = -1.52`` splits into two disjoint accessible
  bands separated by a forbidden gap straddling the L3 vicinity (``x`` in
  roughly ``(-1.117, -0.892)``); only the INNER band (``x`` in roughly
  ``(-0.892, 0)``, containing ``a < a_Jupiter`` orbits including the 3:2
  Hilda resonance at ``a/a_Jupiter = (2/3)**(2/3) ~= 0.7631``) is the
  physically relevant one, matching the paper's own "interior realm"
  description exactly.

* The "quasi-Hilda region" ``R`` and "Mars-crossing region" ``Q`` are both
  defined GEOMETRICALLY on the section via osculating two-body elements
  computed directly from the full (rotating-frame, real-mu) state at each
  section point -- exactly the paper's own definition of the Mars-crossing
  line ("periapsis ``r_p = a(1-e)`` equal to Mars's semimajor axis") applied
  identically to both regions, not merely to ``Q``. ``R`` is centered on the
  3:2 resonance semimajor axis with a moderate-eccentricity ceiling (chosen
  so ``R`` and ``Q`` are geometrically DISJOINT on the section, verified
  directly -- any measured R-to-Q transport is then necessarily a genuine
  dynamical effect, not geometric overlap).

Everything above is independently re-derivable from the paper's own stated
setup (model, section, energy, Mars-crossing-line definition) plus standard
two-body osculating-element formulas; NONE of it depends on reading pixel
coordinates off Fig. 1. This is flagged explicitly, per this task's own
instructions, as the qualitative-reproduction path taken because the digest
is a text-only extraction of a 4-page PRL letter with no machine-readable
figure data.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import solve_ivp

from cyclerfinder.core.cr3bp import CR3BPSystem, cr3bp_eom, cr3bp_system

# Mars / Jupiter heliocentric semimajor axes (au), Standish & Williams J2000 --
# same source (`core/constants.py`'s own `_MARS_SMA_AU` / `_JUPITER_SMA_AU`)
# duplicated here as bare floats to keep this module's public API
# (`SUN_JUPITER_SYSTEM`, `a_mars_nondim`) independent of the planet-registry
# lookup path used for full ephemeris planets (Sun-Jupiter is not one of this
# project's registered CR3BP pairs elsewhere).
_MARS_SMA_AU = 1.52371034
_JUPITER_SMA_AU = 5.20288700

#: Paper's own energy choice (Sec. "The energy we consider..."): just below
#: the L1 equilibrium-point energy E_L1 = -1.5199.
DELLNITZ_2005_ENERGY = -1.52

#: 3:2 mean-motion resonance semimajor axis, in units of Jupiter's own SMA
#: (Kepler's third law ratio, independent of AU normalization).
HILDA_RESONANCE_A_NONDIM = (2.0 / 3.0) ** (2.0 / 3.0)


def sun_jupiter_system() -> CR3BPSystem:
    """Sun-Jupiter CR3BP system, built from this project's own registry."""
    return cr3bp_system("Sun", "Jupiter")


def energy_to_jacobi_constant(energy: float, mu: float) -> float:
    """Convert Dellnitz et al.'s Hamiltonian energy ``E`` to this project's
    Jacobi-constant convention ``C`` (``core.cr3bp.jacobi_constant``).

    Derivation: with canonical momenta ``px = xdot - y``, ``py = ydot + x``,
    Dellnitz's Eq. 1 Hamiltonian reduces algebraically to
    ``H = (v^2 - (x^2+y^2))/2 - (1-mu)/r1 - mu/r2 - mu(1-mu)/2``, while this
    project's ``jacobi_constant`` computes
    ``C = (x^2+y^2) + 2(1-mu)/r1 + 2mu/r2 - v^2``. Comparing term-by-term:
    ``C = -2H - mu(1-mu)``.
    """
    return -2.0 * energy - mu * (1.0 - mu)


def a_mars_nondim() -> float:
    """Mars's heliocentric SMA in Sun-Jupiter CR3BP length units (Jupiter's
    own SMA about the Sun)."""
    return _MARS_SMA_AU / _JUPITER_SMA_AU


def zero_velocity_v(x: float, mu: float, c_target: float) -> float:
    """``V(x) = x^2 + 2(1-mu)/r1 + 2mu/r2 - C`` at the section (``y=0``).

    The accessible ``xdot`` range at this ``x`` is ``|xdot| <=
    sqrt(V(x))`` when ``V(x) >= 0``; the point is inaccessible (would require
    imaginary ``ydot``) when ``V(x) < 0``.
    """
    r1 = abs(x + mu)
    r2 = abs(x - 1.0 + mu)
    return x * x + 2.0 * (1.0 - mu) / r1 + 2.0 * mu / r2 - c_target


def accessible_xdot_max(x: float, mu: float, c_target: float) -> float | None:
    """``sqrt(V(x))``, or ``None`` if ``x`` is outside the accessible zero-
    velocity-curve band at this energy (no real ``ydot`` solves the section
    constraint for any ``xdot`` at this ``x``)."""
    v = zero_velocity_v(x, mu, c_target)
    if v < 0.0:
        return None
    return math.sqrt(v)


def section_state6(x: float, xdot: float, mu: float, c_target: float) -> NDArray[np.float64] | None:
    """Full CR3BP rotating-frame 6-state at section point ``(x, xdot)``, or
    ``None`` if this ``(x, xdot)`` pair is outside the energy manifold
    (``ydot^2 < 0``). ``ydot`` is taken NEGATIVE (the paper's own section
    branch, ``ydot < 0``)."""
    v = zero_velocity_v(x, mu, c_target) - xdot * xdot
    if v < 0.0:
        return None
    ydot = -math.sqrt(v)
    return np.array([x, 0.0, 0.0, xdot, ydot, 0.0], dtype=np.float64)


def _y_zero_event(t: float, y: NDArray[np.float64], mu: float) -> float:
    return float(y[1])


_y_zero_event.terminal = True  # type: ignore[attr-defined]
_y_zero_event.direction = -1.0  # type: ignore[attr-defined]  # only ydot<0 crossings


def poincare_first_return(
    state6: NDArray[np.float64],
    mu: float,
    *,
    t_segment: float = 2.0 * math.pi,
    max_segments: int = 40,
    rtol: float = 1e-10,
    atol: float = 1e-10,
    kickoff_dt: float = 1e-6,
) -> NDArray[np.float64] | None:
    """Propagate ``state6`` forward to the first ``y=0, ydot<0, x<0`` return.

    Every state passed in (the very first call's ``state6``, and every
    ``y = candidate`` re-entry point below) sits EXACTLY on the section
    (``y=0``) by construction. ``solve_ivp``'s event machinery brackets a
    crossing between consecutive accepted steps by testing the event
    function's sign at each step boundary; starting a monitored integration
    exactly AT a zero of the event function makes the very first bracket
    ``[0.0, g(h_1)]`` -- with ``g(h_1) < 0`` since ``ydot < 0`` here -- read as
    a valid "decreasing through zero" crossing AT ``t=0`` itself, spuriously
    self-detecting the START point as the return (confirmed directly: without
    the kick-off below, this function always returned its own input
    unchanged). The standard fix is to advance a tiny FIXED step
    (``kickoff_dt``, far smaller than any physically plausible return time in
    this project's box domain -- the shortest osculating period reachable
    there is ``~0.03`` nondimensional time units, four orders of magnitude
    above ``kickoff_dt``) WITHOUT event monitoring first, landing strictly off
    the section, and only then hand off to the event-monitored segment.

    Integrates each post-kick-off segment for nondimensional duration
    ``t_segment`` (default: one Jupiter orbital period) with a TERMINAL
    ``ydot<0`` section-crossing event; a crossing with ``x >= 0`` (the OTHER
    branch of the section, on Jupiter's side of the Sun) is not accepted --
    the search continues from that crossing point into a fresh
    kick-off-then-monitor segment. Returns ``None`` if no qualifying crossing
    is found within ``max_segments`` segments, or if the underlying
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
        if candidate[0] < 0.0:
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
    """The Poincare return map on section coordinates ``(x, xdot)``.

    Returns the image ``(x_new, xdot_new)`` (shape ``(2,)``), or ``None`` if
    ``(x, xdot)`` is off the energy manifold or the propagation fails to find
    a qualifying return (see :func:`poincare_first_return`). Suitable as the
    ``map_fn`` for
    :func:`cyclerfinder.search.set_oriented_transfer_operator.build_transition_matrix`
    when partially applied over ``(mu, c_target)``.
    """
    state0 = section_state6(x, xdot, mu, c_target)
    if state0 is None:
        return None
    t_segment = float(kwargs.get("t_segment", 2.0 * math.pi))
    max_segments = int(kwargs.get("max_segments", 40))
    result = poincare_first_return(state0, mu, t_segment=t_segment, max_segments=max_segments)
    if result is None:
        return None
    return np.array([result[0], result[3]], dtype=np.float64)


@dataclass(frozen=True)
class OsculatingElements:
    a: float
    e: float
    r_p: float


def osculating_elements_at_section(
    x: float, xdot: float, mu: float, c_target: float
) -> OsculatingElements | None:
    """Osculating heliocentric two-body elements ``(a, e, r_p)`` at section
    point ``(x, xdot)`` -- the general-conic vis-viva/angular-momentum
    formulas (valid for ellipse, parabola or hyperbola alike), applied to the
    instantaneous (rotating-frame-aligned-with-inertial-at-this-instant)
    heliocentric position/velocity. This is EXACTLY the paper's own
    definition of the Mars-crossing line (``r_p = a(1-e)``), just evaluated
    directly from the state rather than backed out from an assumed orbital-
    element parametrization -- so it applies unchanged to any section point,
    not only ones constructed from a chosen ``(a, e, omega)``.

    Heliocentric position: Sun sits at ``(-mu, 0)`` in the rotating frame, so
    ``r_vec = (x + mu, y)``; heliocentric velocity uses the standard
    rotating-to-inertial relation (frame angular rate 1, aligned at ``t=0``):
    ``vx_inertial = xdot - y``, ``vy_inertial = ydot + x``. Returns ``None``
    if ``(x, xdot)`` is off the energy manifold.
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
    return OsculatingElements(a=a, e=ecc, r_p=r_p)


def mars_crossing_indicator(x: float, xdot: float, mu: float, c_target: float) -> bool | None:
    """``True`` iff the osculating periapsis at this section point is inside
    Mars's own heliocentric SMA (the paper's Mars-crossing-line criterion).
    ``None`` if ``(x, xdot)`` is off the energy manifold."""
    elems = osculating_elements_at_section(x, xdot, mu, c_target)
    if elems is None:
        return None
    return elems.r_p <= a_mars_nondim()


def quasi_hilda_region_indicator(
    x: float,
    xdot: float,
    mu: float,
    c_target: float,
    *,
    a_lo: float = 0.70,
    a_hi: float = 0.85,
    e_max: float = 0.35,
) -> bool | None:
    """``True`` iff the osculating ``(a, e)`` at this section point sits in a
    window around the 3:2 Hilda resonance (default ``a`` in ``[0.70, 0.85]``
    x Jupiter's SMA, bracketing ``HILDA_RESONANCE_A_NONDIM ~= 0.7631``) at
    moderate eccentricity (default ``e <= 0.35``, comfortably below what any
    orbit at ``a`` in this range needs to become Mars-crossing -- verified
    directly: this default window is geometrically DISJOINT from
    :func:`mars_crossing_indicator`'s region at ``E = -1.52`` for the Sun-
    Jupiter system, over a direct grid check).

    ``None`` if ``(x, xdot)`` is off the energy manifold.
    """
    elems = osculating_elements_at_section(x, xdot, mu, c_target)
    if elems is None:
        return None
    return bool(a_lo <= elems.a <= a_hi and elems.e <= e_max)
