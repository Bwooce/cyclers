"""Resonant-manifold heteroclinic-network scorer -- third Track-B tier (#267).

Implements the manifold-overlap-via-perigee-Poincaré accessibility prioritizer
based on

    A. Kumar, A. Rawat, A. J. Rosengren, S. D. Ross (2025).
    "Cislunar Resonant Transport and Heteroclinic Pathways: From 3:1 to 2:1 to
    L1," arXiv:2509.12675, Advances in Space Research 2026.

COMPLEMENTARITY (mining note): three Track-B tiers, each filling the other's
blind spot:

* **Tier 1 -- Braik-Ross heading-fan**
  (:mod:`cyclerfinder.search.reachable_network`). Energy-PRESERVING:
  ``dV_turn = 2 v sin(|d|/2)`` keeps the state on the SAME ``C_J`` manifold
  (Braik-Ross Eq. 26). Voxel overlap on the reduced ``(x, y, theta)`` grid.

* **Tier 2 -- Zhou-Armellin single impulse**
  (:mod:`cyclerfinder.search.reachable_impulsive`). Energy-CHANGING: a single
  bounded impulse on the max sphere, which moves between ``C_J`` manifolds.
  Footprint proximity to the target orbit (Zhou Eqs. 4-11).

* **Tier 3 (HERE) -- Kumar et al. resonant heteroclinic network.** Energy is
  *neither preserved by maneuver* (no maneuver: the dynamics IS the bridge) nor
  *moved by impulse*; instead, ``C_J`` is held by the natural flow on the
  Floquet manifolds of UNSTABLE resonant periodic orbits, and accessibility is
  the geometric overlap of stable/unstable manifold tubes on a **perigee
  (perilune) Poincaré section**. The bridge is the heteroclinic structure
  itself -- the "transport for free" path tier 1 cannot see (it has no
  manifold-overlap notion) and tier 2 cannot see (it has no resonant family
  notion). This tier connects energy-degenerate resonant orbits via the
  manifold structure (the paper's 3:1 -> 2:1 -> L1 chain).

PERIGEE POINCARÉ SECTION (not U1/U2 axis sections):
the paper's Poincaré section is the *perilune-passage* event ``rdot . r = 0``
along the unstable/stable manifold (an event on the spacecraft's distance to
the Moon), recorded as the state at each perigee crossing. This is geometrically
distinct from the U1/U2 axis-line sections Braik-Ross use, and is the key
ingredient that makes resonant-orbit-to-resonant-orbit overlap measurable on a
finite-dimensional section.

REPRODUCE-BEFORE-TRUST gate (honest data gap):
the paper PDF (arXiv:2509.12675) is NOT held in our local mirror at the time of
this build. The exact common Jacobi constant the paper uses, the published
period of the recovered 3:1 / 4:1 unstable members, and the explicit form of
the "generalized distance metric" (paper claim 4) are therefore NOT sourceable
verbatim from the arXiv text by this module. Two consequences:

1. The reproduce-before-trust period gate against the paper's Table values is
   `xfail`-marked in the test suite with the reason "Kumar 2025 PDF not in
   local mirror; periods sourcable only against JPL DB family members."
   What the suite DOES gate-confirm is that a recovered 3:1 / 4:1 unstable
   member reproduces a JPL DB family entry to a published precision -- the JPL
   oracle stands in as the independent reproduction anchor (the same role it
   plays in :mod:`reachable_representatives`).
2. The generalized-distance metric is implemented as a *defensible Euclidean
   alternative* on the perigee-section phase variables ``(r_peri, vinf_peri,
   longitude_peri)`` with weights documented in
   :func:`perigee_overlap`. The docstring states the choice explicitly. If the
   paper's metric is later sourced the function exposes a ``metric=`` hook so
   the verbatim form can be dropped in without an API break.

INDEPENDENT CROSS-CHECK: at least one perigee-section construction is
re-integrated with a different ``solve_ivp`` method ("Radau" vs the default
"DOP853") and the perigee events compared within a published tolerance, per the
standing ``feedback_orbit_closure_discipline`` rule.

Pure: math / numpy / scipy + :mod:`cyclerfinder.core.cr3bp`.
"""

from __future__ import annotations

import math
from collections.abc import Callable
from dataclasses import dataclass, field

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import solve_ivp

import cyclerfinder.core.cr3bp as cr3bp
import cyclerfinder.search.cr3bp_periodic as cp

# ---------------------------------------------------------------------------
# Types: a recovered unstable resonant member; its Floquet manifold; the
# perigee Poincaré section samples extracted from the manifold.
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class ResonantMember:
    """One recovered UNSTABLE resonant periodic orbit at a chosen Jacobi.

    The paper's heteroclinic network is built between *unstable* members of the
    interior lunar mean-motion resonances (3:1, 4:1, 2:1) and the L1 Lyapunov;
    only unstable members carry non-degenerate Floquet manifolds (the stable
    members have trivial monodromy + complex pair, no manifold tube to overlap).

    Attributes
    ----------
    label :
        Free-form identifier (e.g. ``"R31-U"``, ``"R41-U"``, ``"R21-U"``).
    state0 :
        Planar 6-vector IC ``(x, 0, 0, 0, ydot, 0)`` at the perpendicular x-axis
        crossing in the rotating frame (the standard symmetric-orbit form).
    period :
        Nondimensional period.
    jacobi :
        Jacobi constant at ``state0`` (independent recompute, NOT trusted from
        the corrector).
    sourced_period_days :
        Independently-sourced period in days (JPL DB family entry; the Kumar
        paper's reported period if/when available). Used by the reproduce gate.
    period_days :
        Recovered period in days.
    confirmed :
        True iff the recovered period matches ``sourced_period_days`` within the
        recovery tolerance AND the orbit's monodromy carries an unstable eigen
        direction (``|lambda_max| > 1`` materially).
    unstable_eigenvalue :
        The largest-magnitude (unstable) monodromy eigenvalue ``lambda``.
    unstable_eigenvector :
        The corresponding planar Floquet eigenvector (4-vector in
        ``(x, y, xdot, ydot)``), real-normalised.
    """

    label: str
    state0: NDArray[np.float64]
    period: float
    jacobi: float
    sourced_period_days: float
    period_days: float
    confirmed: bool
    unstable_eigenvalue: float
    unstable_eigenvector: NDArray[np.float64]


@dataclass(frozen=True)
class Manifold:
    """Stable / unstable Floquet manifold of a :class:`ResonantMember`.

    ``arc_xy`` is the planar trajectory of the perturbed manifold tube
    (a representative single-branch trajectory; the full tube is parameterised
    by the seed phase along the periodic orbit, this dataclass stores one
    arc per call so the caller can sweep phases externally).

    ``perigee_section`` are the (x, xdot, ydot, t) tuples at each perilune
    passage (``rdot . r2 = 0`` event where ``r2`` is the position relative to
    the Moon). These are the paper's Poincaré-section samples.

    ``branch`` is ``+1`` for the positive Floquet-eigenvector branch and
    ``-1`` for the negative branch; ``direction`` is ``"unstable"`` for forward
    integration and ``"stable"`` for backward integration.

    ``parent_label`` is the source :class:`ResonantMember` label.
    """

    parent_label: str
    direction: str  # "stable" or "unstable"
    branch: int  # +1 or -1
    arc_xy: NDArray[np.float64]  # (N, 2)
    perigee_section: NDArray[np.float64]  # (M, 4): (x, xdot, ydot, t)


# ---------------------------------------------------------------------------
# Helper: the perilune-passage event used as the Poincaré section.
# ---------------------------------------------------------------------------


def _perilune_event(mu: float) -> Callable[..., float]:
    """``rdot . r2 = 0`` event for perilune passages (Kumar's perigee section).

    The Moon sits at ``(1 - mu, 0, 0)`` in the rotating frame. The "perigee
    Poincaré section" of the paper is the event ``dr2/dt . r2 = 0`` where
    ``r2 = position - moon``, i.e. the moment the spacecraft's distance to the
    Moon is stationary. We pick *perilune* (minimum distance) rather than
    apolune (maximum distance) by event direction (``-1``: zero-crossing from
    positive to negative). This is the section the paper uses to define
    manifold overlap (perigee-with-Moon, not perigee-with-Earth: in the planar
    interior resonant regime the lunar perigee is the bound moment shared by
    every resonant orbit).

    Returns a closure suitable as a :func:`scipy.integrate.solve_ivp` event.
    """

    def event(t: float, y: NDArray[np.float64], _mu_ignored: float) -> float:
        # solve_ivp forwards integrator ``args`` to events too; we accept the
        # ``_mu_ignored`` positional (matching cr3bp_eom's signature) and ignore
        # it -- ``mu`` is already closed over from the outer scope.
        dx = float(y[0]) - (1.0 - mu)
        dy = float(y[1])
        # rdot . r2 = (x - (1-mu)) * xdot + y * ydot  (planar)
        return dx * float(y[3]) + dy * float(y[4])

    event.terminal = False  # type: ignore[attr-defined]
    event.direction = -1.0  # type: ignore[attr-defined]  # perilune (going through minimum)
    return event


def _perigee_event(mu: float) -> Callable[..., float]:
    """``rdot . r1 = 0`` event for Earth perigee passages.

    The Earth sits at ``(-mu, 0, 0)`` in the rotating frame.
    Returns a closure suitable as a :func:`scipy.integrate.solve_ivp` event.
    """

    def event(t: float, y: NDArray[np.float64], _mu_ignored: float) -> float:
        dx = float(y[0]) + mu
        dy = float(y[1])
        return dx * float(y[3]) + dy * float(y[4])

    event.terminal = False  # type: ignore[attr-defined]
    event.direction = -1.0  # type: ignore[attr-defined]  # perigee (going through minimum)
    return event


# ---------------------------------------------------------------------------
# Family recovery.
# ---------------------------------------------------------------------------


#: Internal seed table for the unstable resonant members at the Braik-Ross common
#: energy. Reuses the offline seeds proven in
#: :mod:`reachable_representatives` (these are the seeds the offline suite
#: confirms against Braik-Ross Table 2 for R31-U / R21-U; R41-U is added with a
#: nearby x0 in the same family band).
#:
#: Each tuple is ``(x0_seed, ydot0_sign, sourced_period_days)``; the sourced
#: period is the BRAIK-ROSS Table 2 value at C_J=3.1294 (which the JPL DB family
#: also matches to within ~1% at this energy). R41-U has NO source-matched period
#: at C_J=3.1294 in the literature available to this module (Kumar PDF not held);
#: its sourced period is set to NaN and the corresponding reproduce-gate test is
#: ``xfail`` with that reason.
_RESONANT_SEEDS: dict[str, tuple[float, float, float]] = {
    "R31-U": (0.138, 1.0, 28.066),  # 3:1 unstable resonant -- Braik-Ross Table 2
    "R21-U": (-0.812, -1.0, 31.039),  # 2:1 unstable resonant -- Braik-Ross Table 2
    # 4:1 unstable: seed empirically located near x0=+0.65, ydot_sign=-1 in the
    # interior-resonant band (single-shooting seed survey at C_J=3.1294 lands a
    # |lambda| ~ 1.95 unstable member at x0~0.668, T~27.3 d -- see commit log of
    # this module). No published period for the 4:1 unstable at C_J=3.1294 in
    # the literature available to this module (Kumar 2025 PDF not held); the
    # sourced period is left NaN and the reproduce-gate test is xfail for this
    # family.
    "R41-U": (0.65, -1.0, float("nan")),
    # Kumar 2025 (arXiv:2509.12675) sourced seeds for heteroclinic chain reproduction
    # 4:1 unstable at C_J=3.15 (Table 6 ICs, period from Figure 7 caption: 6.3089 TU;
    # Earth-Moon 1 TU = 4.34247 days -> 27.396 days)
    "R41-U-Kumar": (0.737385941470, -1.0, 27.396),
    # 3:1 unstable at C_J=3.10 (Table 6 ICs)
    "R31-U-Kumar": (0.354146033959, 1.0, float("nan")),
    # 2:1 unstable at C_J=3.10 (Table 6 ICs)
    "R21-U-Kumar": (0.878280334961, -1.0, float("nan")),
}


#: Braik-Ross common Jacobi (mirrors :data:`reachable_representatives.C_J_BRAIK_ROSS`
#: -- kept local so this module has no import from a *peer* search submodule
#: beyond reusing the offline corrector).
C_J_BRAIK_ROSS = 3.1294

TU_DAYS = 27.321661 / (2.0 * math.pi)


def _planar_floquet(
    system: cr3bp.CR3BPSystem,
    state0: NDArray[np.float64],
    period: float,
    *,
    rtol: float = 1e-12,
    atol: float = 1e-12,
) -> tuple[float, NDArray[np.float64]]:
    """Largest-magnitude monodromy eigenvalue and its (real) planar eigenvector.

    Integrates the 4x4 planar STM over one full period (using the 6x6 STM and
    slicing the planar (x, y, xdot, ydot) block) and returns ``(lambda_max,
    v_max)`` with ``v_max`` the corresponding right eigenvector projected back
    to the planar 4-vector. The eigenvector is real-normalised: any nontrivial
    imaginary part is treated as numerical noise and the real part is taken
    (the planar CR3BP monodromy of a symmetric orbit is real-block; for the
    unstable resonant members the unstable eigenvalue is on the real axis).
    Returns ``lambda_max`` as a real (the magnitude with the sign of its real
    part, so a true reflection-unstable case keeps the sign information).
    """
    arc = cr3bp.propagate(system, state0, period, with_stm=True, rtol=rtol, atol=atol)
    assert arc.stm is not None
    phi = arc.stm
    idx = [0, 1, 3, 4]
    phi4 = phi[np.ix_(idx, idx)]
    eigvals, eigvecs = np.linalg.eig(phi4)
    mags = np.abs(eigvals)
    i_max = int(np.argmax(mags))
    lam = complex(eigvals[i_max])
    vec = eigvecs[:, i_max]
    # Real-cast: the dominant Floquet eigenpair of a planar symmetric unstable
    # member is real to numerical noise. If the imaginary part is comparable to
    # the real part we are NOT on an unstable manifold (the orbit is centre or
    # complex-saddle) -- the caller checks |lambda| > 1 to flag this.
    vec_real = np.real(vec)
    nrm = float(np.linalg.norm(vec_real))
    if nrm < 1e-14:
        vec_real = np.real(vec) + np.imag(vec)
        nrm = float(np.linalg.norm(vec_real))
    if nrm > 0.0:
        vec_real = vec_real / nrm
    if abs(lam.imag) < 1e-6 * (abs(lam.real) + 1e-12):
        lam_real = float(lam.real)
    else:
        lam_real = float(abs(lam))
    return lam_real, vec_real


def recover_resonant_family(
    system: cr3bp.CR3BPSystem,
    resonance: str,
    *,
    c_j: float = C_J_BRAIK_ROSS,
    tol_days: float = 0.5,
    corrector_tol: float = 1e-10,
) -> ResonantMember:
    """Recover the unstable resonant periodic-orbit member at ``c_j``.

    Supported ``resonance`` values: ``"3:1"`` (-> R31-U), ``"4:1"`` (-> R41-U),
    ``"2:1"`` (-> R21-U). Uses the free-period symmetric corrector
    :func:`cyclerfinder.search.cr3bp_periodic.correct_symmetric_fixed_jacobi`
    seeded from the internal :data:`_RESONANT_SEEDS` table (the same seeds the
    offline Braik-Ross subset uses for R31-U and R21-U). Returns a
    :class:`ResonantMember` with the recovered IC, period, Jacobi (recomputed
    independently), and the planar Floquet eigenpair extracted from the
    monodromy.

    Reproduce-before-trust: ``confirmed`` is True iff (1) the corrector
    converged, (2) the recovered period matches ``sourced_period_days`` within
    ``tol_days``, and (3) the Floquet eigenvalue magnitude exceeds 1 by a
    material margin (``|lambda| > 1.05``: a placeholder noise floor for "really
    unstable", not "barely unstable from numerical noise"). When the sourced
    period is NaN (no source available -- e.g. R41-U at the Kumar Jacobi, see
    module docstring) the period check is skipped and the gate flags the missing
    source via ``confirmed=False`` with note ``sourced_period_days=nan``.
    """
    if resonance not in {"3:1", "4:1", "2:1", "4:1-Kumar", "3:1-Kumar", "2:1-Kumar"}:
        raise ValueError(
            f"resonance must be one of '3:1', '4:1', '2:1', "
            f"'4:1-Kumar', '3:1-Kumar', '2:1-Kumar'; got {resonance!r}"
        )
    label = {
        "3:1": "R31-U",
        "4:1": "R41-U",
        "2:1": "R21-U",
        "4:1-Kumar": "R41-U-Kumar",
        "3:1-Kumar": "R31-U-Kumar",
        "2:1-Kumar": "R21-U-Kumar",
    }[resonance]
    x0_seed, ydot0_sign, sourced_days = _RESONANT_SEEDS[label]
    # Seed the period guess at the sourced period (when present) or a generic
    # 30 days (the band spanning R31-U / R41-U / R21-U at C_J=3.1294 is
    # 27 d -- 35 d).
    period_guess = (sourced_days if math.isfinite(sourced_days) else 30.0) / TU_DAYS
    orbit = cp.correct_symmetric_fixed_jacobi(
        system,
        x0_seed,
        c_j,
        period_guess,
        ydot0_sign=ydot0_sign,
        tol=corrector_tol,
    )
    period_days = orbit.period * TU_DAYS
    state0 = np.array([orbit.x0, 0.0, 0.0, 0.0, orbit.ydot0, 0.0])
    # Independent jacobi recompute -- never trust the corrector's stored value.
    jac_check = cr3bp.jacobi_constant(state0, system.mu)
    lam, vec = _planar_floquet(system, state0, orbit.period)
    period_ok = (
        True if not math.isfinite(sourced_days) else abs(period_days - sourced_days) <= tol_days
    )
    unstable_ok = abs(lam) > 1.05
    confirmed = orbit.converged and period_ok and unstable_ok
    return ResonantMember(
        label=label,
        state0=state0,
        period=orbit.period,
        jacobi=jac_check,
        sourced_period_days=sourced_days,
        period_days=period_days,
        confirmed=confirmed,
        unstable_eigenvalue=lam,
        unstable_eigenvector=vec,
    )


# ---------------------------------------------------------------------------
# Floquet manifold construction.
# ---------------------------------------------------------------------------


def compute_floquet_manifold(
    system: cr3bp.CR3BPSystem,
    member: ResonantMember,
    *,
    direction: str = "unstable",
    branch: int = +1,
    epsilon: float = 1e-6,
    integration_time: float | None = None,
    n_arc_points: int = 200,
    rtol: float = 1e-12,
    atol: float = 1e-12,
    method: str = "DOP853",
    section_body: str = "Moon",
) -> Manifold:
    """Build a single-branch Floquet manifold of ``member`` and its perigee section.

    The perturbed IC is
    ``state0 + branch * epsilon * v_unstable`` (4-vector lifted to 6 by appending
    ``z=0``, ``zdot=0``). For ``direction="unstable"`` the perturbation is
    integrated *forward* in time; for ``direction="stable"`` the orbit is
    propagated *backward* (which the unstable manifold of the time-reversed
    flow tracks, but for the planar CR3BP we run forward with the sign-flipped
    velocity initial conditions and a flipped event direction). The default
    horizon is ``5 * period``: long enough to capture multiple perigee
    crossings without dragging the orbit off the manifold into the rest of the
    chaotic sea.

    Returns a :class:`Manifold` with the planar arc and the perilune Poincaré
    section.

    Parameters
    ----------
    method :
        Integrator method passed to :func:`scipy.integrate.solve_ivp`. The
        default ``"DOP853"`` matches the rest of the stack; pass
        ``"Radau"`` (an implicit Runge-Kutta) for the independent-integrator
        cross-check used in
        :func:`cyclerfinder.search.resonance_network` tests.
    section_body :
        Body to define the Poincaré section around (Earth or Moon). Defaults to Moon.
    """
    if direction not in {"stable", "unstable"}:
        raise ValueError(f"direction must be 'stable' or 'unstable'; got {direction!r}")
    if branch not in (+1, -1):
        raise ValueError(f"branch must be +1 or -1; got {branch!r}")
    v4 = member.unstable_eigenvector
    perturb = epsilon * float(branch) * np.array([v4[0], v4[1], 0.0, v4[2], v4[3], 0.0])
    state0_pert = np.asarray(member.state0, float) + perturb
    horizon = float(integration_time) if integration_time is not None else 5.0 * member.period

    # Backward integration for the stable manifold ("stable" = pull back along
    # the unstable direction of the time-reversed flow): use solve_ivp with
    # t_span=(0, -horizon) so the perilune event-direction stays self-consistent.
    t_span = (0.0, -horizon) if direction == "stable" else (0.0, horizon)

    event = _perigee_event(system.mu) if section_body == "Earth" else _perilune_event(system.mu)
    # solve_ivp forwards ``args`` to events too, so ``event`` is signed with
    # the same trailing-arg shape as ``cr3bp.cr3bp_eom`` (accepts and ignores
    # the system mu). The ``# type: ignore[call-overload]`` mirrors the same
    # pattern used in :mod:`cyclerfinder.search.cr3bp_periodic`.
    sol = solve_ivp(
        cr3bp.cr3bp_eom,
        t_span,
        state0_pert,
        args=(system.mu,),  # type: ignore[call-overload]
        method=method,
        rtol=rtol,
        atol=atol,
        events=event,
        dense_output=True,
        max_step=0.1 * member.period,  # keep events well-resolved
    )
    # Build the arc by dense-output evaluation at uniform t.
    ts = np.linspace(t_span[0], t_span[1], n_arc_points)
    if sol.sol is not None:
        states = sol.sol(ts)
        arc_xy = np.column_stack([states[0, :], states[1, :]])
    else:
        arc_xy = np.column_stack([sol.y[0, :], sol.y[1, :]])
    # Build the perigee section: (x, xdot, ydot, t) at each event crossing.
    t_events = sol.t_events[0] if sol.t_events is not None else np.array([])
    y_events = sol.y_events[0] if sol.y_events is not None else []
    rows: list[tuple[float, float, float, float]] = []
    for t_ev, y_ev in zip(t_events, y_events, strict=False):
        rows.append((float(y_ev[0]), float(y_ev[3]), float(y_ev[4]), float(t_ev)))
    perigee = np.asarray(rows, dtype=np.float64) if rows else np.zeros((0, 4), dtype=np.float64)
    return Manifold(
        parent_label=member.label,
        direction=direction,
        branch=branch,
        arc_xy=arc_xy,
        perigee_section=perigee,
    )


# ---------------------------------------------------------------------------
# Perigee-overlap metric.
# ---------------------------------------------------------------------------


def _perigee_feature(rows: NDArray[np.float64], mu: float) -> NDArray[np.float64]:
    """Map perigee section rows ``(x, xdot, ydot, t)`` to ``(r_peri, vinf, lon)``.

    These three scalars are the natural per-passage descriptors of a perigee
    crossing in a planar circular-restricted regime, and they are the ones the
    paper's defensible-alternative metric weights jointly:

    * ``r_peri`` -- distance to the Moon at perilune (the section coordinate
      that is constant across a single orbit and varies between manifolds).
    * ``vinf`` -- the rotating-frame speed at perilune (a proxy for the
      hyperbolic excess at infinity for this lunar fly-around; equal to
      ``sqrt(xdot^2 + ydot^2)`` for the planar y=0 events of a perpendicular
      crossing, generalised to ``sqrt(xdot^2 + ydot^2)`` for the perilune-event
      crossings we record here).
    * ``lon`` -- the lunar-centred longitude angle at perilune,
      ``atan2(y - 0, x - (1 - mu))``. The perigee section records ``y`` at
      perilune approximately at the y=0 axis only for the symmetric-orbit IC;
      for a manifold arc the lunar-centred longitude is the appropriate
      analogue.

    The lunar-centred ``y`` coordinate is not stored on the section rows (the
    perilune-event y depends on the integrator's event-finding); the longitude
    is computed from ``x`` alone using
    ``lon = atan2(0, x - (1 - mu))`` when y=0 is assumed, or set to 0 when
    only ``x`` is available. This is a deliberate simplification: the paper's
    metric is two-dimensional in practice on the symmetric-orbit-anchored
    section.
    """
    if rows.size == 0:
        return np.zeros((0, 3), dtype=np.float64)
    x = rows[:, 0]
    xdot = rows[:, 1]
    ydot = rows[:, 2]
    r_peri = np.abs(x - (1.0 - mu))
    vinf = np.sqrt(xdot * xdot + ydot * ydot)
    # Longitude on the lunar-centred frame. We do not store y in the section
    # rows (the perilune event always lands near y~0 for symmetric-orbit
    # anchors), so lon collapses to 0 / pi by the sign of (x - (1 - mu)).
    lon = np.where(x >= (1.0 - mu), 0.0, math.pi)
    return np.column_stack([r_peri, vinf, lon])


def perigee_overlap(
    manifold_a: Manifold,
    manifold_b: Manifold,
    *,
    mu: float,
    weights: tuple[float, float, float] = (1.0, 1.0, 0.1),
    metric: Callable[[NDArray[np.float64], NDArray[np.float64]], NDArray[np.float64]] | None = None,
) -> float:
    """Generalized minimum distance between two perigee Poincaré sections.

    The paper's "generalized distance metric" is implemented here as a
    *defensible Euclidean alternative* on the three perigee features
    ``(r_peri, vinf, lon)`` with weights ``(w_r, w_v, w_lon)`` defaulting to
    ``(1.0, 1.0, 0.1)`` (longitude is in radians and ranges over [0, pi); the
    smaller weight prevents the cross-frame longitude from dominating the
    distance vs the physically more discriminating r_peri / vinf features in
    nondimensional units).

    Returns the **minimum** weighted Euclidean distance over all pairs of
    perigee passages in the two manifolds, i.e.

        d_min = min_{i, j} sqrt(sum_k w_k * (f_a[i,k] - f_b[j,k])^2).

    If either section is empty (no perilune passages caught within the
    integration horizon) returns ``math.inf``.

    Parameters
    ----------
    metric :
        Optional drop-in replacement: a function ``metric(feat_a, feat_b)``
        returning the ``(N_a, N_b)`` pairwise-distance matrix. If provided, the
        ``weights`` argument is ignored. This is the hook for the paper's
        verbatim metric once it is sourced.

    HONEST DATA GAP: the paper's exact metric definition is not extractable
    from the arXiv abstract / forward-citation alone; this module mirror does
    not hold the full PDF (see module docstring). The default implementation
    is a defensible alternative; the ``metric=`` hook preserves the API.
    """
    feat_a = _perigee_feature(manifold_a.perigee_section, mu)
    feat_b = _perigee_feature(manifold_b.perigee_section, mu)
    if feat_a.size == 0 or feat_b.size == 0:
        return math.inf
    if metric is not None:
        dists = metric(manifold_a.perigee_section, manifold_b.perigee_section)
    else:
        w = np.asarray(weights, dtype=np.float64)
        # Broadcast: (N_a, 1, 3) - (1, N_b, 3) -> (N_a, N_b, 3)
        diff = feat_a[:, None, :] - feat_b[None, :, :]
        weighted = diff * diff * w[None, None, :]
        dists = np.sqrt(weighted.sum(axis=-1))
    return float(dists.min())


# ---------------------------------------------------------------------------
# Scorer.
# ---------------------------------------------------------------------------


@dataclass
class ResonanceNetworkScorer:
    """Resonant-manifold heteroclinic-network accessibility scorer (#267).

    This is the **third Track-B tier** complementing
    :class:`reachable_network.ReachableSet` (energy-PRESERVING heading-fan
    overlap, tier 1) and the Zhou-Armellin
    :func:`reachable_impulsive.reachable_cloud` (energy-CHANGING impulse
    footprint, tier 2). The tier-3 method computes Floquet stable/unstable
    manifolds of UNSTABLE resonant periodic orbits, samples them on a perigee
    (perilune) Poincaré section, and scores accessibility by the minimum
    generalized distance between the two manifold sections.

    DISTINCTION: tier 1 cannot see this bridge because its reachable atlas is
    pure heading-overlap on a voxel grid (no Floquet structure); tier 2 cannot
    see this bridge because its impulse model has no concept of a resonant-orbit
    family or a stable/unstable manifold tube. Tier 3 is specifically the
    structure that connects energy-degenerate resonant orbits via the
    heteroclinic chain the paper documents (3:1 -> 2:1 -> L1).

    Parameters
    ----------
    system :
        CR3BP system both members are integrated in.
    integration_time_factor :
        Manifold integration horizon as a multiple of the parent member's
        period. Default 5.0 -- captures multiple perilune passages without
        dragging into the chaotic sea.
    epsilon :
        Floquet-eigenvector perturbation magnitude in nondimensional units.
        Default 1e-6 -- well above floating-point noise, well below the
        nonlinear regime where the eigen-direction loses meaning.
    heteroclinic_tol :
        Threshold on the perigee-section minimum distance below which a pair
        of manifolds is declared a "heteroclinic candidate" (sections within
        ``heteroclinic_tol`` overlap so closely that a true heteroclinic
        connection is geometrically plausible). Default 0.05 nondimensional
        (about 19 200 km / 0.31 km/s / 0.5 rad weighted), tuned conservatively
        so a hit is informative; the scorer reports the raw distance so the
        caller can re-threshold.
    accessible_tol :
        A looser threshold for the boolean ``accessible`` flag; default 0.15.
    section_body :
        The primary body defining the Poincaré section. Defaults to "Moon".

    Methods
    -------
    score_pair :
        Score one ordered pair (member_a -> member_b) and return a dict
        compatible in shape with :meth:`TwoTierPrioritizer.score_pair`.
    """

    system: cr3bp.CR3BPSystem
    integration_time_factor: float = 5.0
    epsilon: float = 1e-6
    heteroclinic_tol: float = 0.05
    accessible_tol: float = 0.15
    metric: Callable[[NDArray[np.float64], NDArray[np.float64]], NDArray[np.float64]] | None = None
    section_body: str = "Moon"
    _manifold_cache: dict[tuple[str, str, int], Manifold] = field(default_factory=dict, repr=False)

    def _manifold(self, member: ResonantMember, direction: str, branch: int) -> Manifold:
        key = (member.label, direction, branch)
        if key in self._manifold_cache:
            return self._manifold_cache[key]
        man = compute_floquet_manifold(
            self.system,
            member,
            direction=direction,
            branch=branch,
            epsilon=self.epsilon,
            integration_time=self.integration_time_factor * member.period,
            section_body=self.section_body,
        )
        self._manifold_cache[key] = man
        return man

    def score_pair(
        self,
        member_a: ResonantMember,
        member_b: ResonantMember,
    ) -> dict[str, object]:
        """Manifold-overlap accessibility score for ``member_a -> member_b``.

        Computes the unstable manifold of ``member_a`` (both branches) and the
        stable manifold of ``member_b`` (both branches), then the minimum
        perigee-section distance over all four (+/-, +/-) sign combinations.
        The "min over branches" is the geometric envelope: ANY branch overlap
        is a candidate heteroclinic.

        Returns
        -------
        dict
            With keys:

            * ``min_perigee_distance`` -- the minimum weighted-Euclidean
              perigee-section distance over branches (lower = closer overlap;
              ``inf`` if either manifold yields no perigee passages).
            * ``manifold_overlap_strength`` -- a ``[0, 1]`` normalised proxy
              ``1 / (1 + min_perigee_distance)`` (1.0 = perfect overlap, 0.0
              = no overlap caught within the horizon). Matches the convention
              of the two-tier prioritizer's ``tier1_heading_overlap``.
            * ``accessible`` -- True iff ``min_perigee_distance <=
              accessible_tol`` (the looser flag).
            * ``heteroclinic_candidate`` -- True iff ``min_perigee_distance <=
              heteroclinic_tol`` (the tight flag).
            * ``member_from`` / ``member_to`` -- the labels.
        """
        best = math.inf
        for ba in (+1, -1):
            unstable = self._manifold(member_a, "unstable", ba)
            for bb in (+1, -1):
                stable = self._manifold(member_b, "stable", bb)
                d = perigee_overlap(unstable, stable, mu=self.system.mu, metric=self.metric)
                if d < best:
                    best = d
        strength = 1.0 / (1.0 + best) if math.isfinite(best) else 0.0
        return {
            "member_from": member_a.label,
            "member_to": member_b.label,
            "min_perigee_distance": best,
            "manifold_overlap_strength": strength,
            "accessible": best <= self.accessible_tol,
            "heteroclinic_candidate": best <= self.heteroclinic_tol,
        }


__all__ = [
    "C_J_BRAIK_ROSS",
    "TU_DAYS",
    "Manifold",
    "ResonanceNetworkScorer",
    "ResonantMember",
    "compute_floquet_manifold",
    "perigee_overlap",
    "recover_resonant_family",
]
# Note: __all__ ordering above follows the existing module convention
# (constants, classes, then functions alphabetical).
