"""Circular-restricted three-body problem (CR3BP) dynamics core (spec 2026-06-10).

Nondimensional rotating frame; mu = m2/(m1+m2); primary at (-mu,0,0), secondary at
(1-mu,0,0). See the plan's "CR3BP equations" block. Pure: math/numpy/scipy +
core.satellites only.
"""

from __future__ import annotations

import math
import time
from collections.abc import Callable
from dataclasses import dataclass
from typing import TYPE_CHECKING, Literal

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import solve_ivp

from cyclerfinder.core.constants import AU_KM, MU_SUN_KM3_S2, PLANETS, PlanetData
from cyclerfinder.core.satellites import PRIMARIES, SATELLITES

if TYPE_CHECKING:
    # Only exposed under scipy's private `_ivp` submodule at runtime (not
    # `scipy.integrate.OdeResult`) -- imported solely for the type annotation
    # below, never executed, so this stays robust to scipy internals moving.
    from scipy.integrate._ivp.ivp import OdeResult

# STM integration modes. See ``propagate(with_stm=True, stm_mode=...)`` for the
# semantics; ``"variable"`` is the legacy variable-step variational path,
# ``"fixed_path"`` is the Pellegrini-Russell 2016 mitigation (record the state
# only path's accepted step grid, replay state+STM along that pre-scheduled
# grid so the step size has no IC dependence).
StmMode = Literal["variable", "fixed_path"]


class CR3BPCloseEncounterError(RuntimeError):
    """A #652 ``solve_ivp`` termination event fired: the trajectory dynamically
    evolved into a genuine close encounter with the SECONDARY body partway
    through STM-augmented propagation.

    Subclasses ``RuntimeError`` deliberately (not a bare new exception type) so
    every existing ``except RuntimeError`` caller of ``propagate(with_stm=True)``
    (``correct_periodic`` and its many downstream callers across ``search/`` and
    ``genome/``) keeps working unchanged -- this is a MORE SPECIFIC failure than
    the pre-existing generic ``RuntimeError`` raised when ``sol.success`` is
    False (e.g. a collision trajectory near the PRIMARY), distinguishable by
    type and by its own honestly-worded message, never silently conflated with
    it.

    Background (#651, #652): the augmented state+STM variational equations
    (:func:`cr3bp_stm_eom`) contain second-derivative terms that scale as
    ``~1/r2^5`` near the secondary body -- far more singular than the state
    equations' own ``~1/r2^2`` -- so a close approach to the secondary that only
    develops DURING integration (not visible from the initial condition alone)
    can collapse the adaptive DOP853 step size and grind for minutes+ without
    ever tripping ``sol.success=False``. This event-based guard fires FAST,
    before that grind, by monitoring distance-to-secondary directly (the exact
    quantity driving the singularity) rather than relying on a coarse
    wall-clock/step-count budget.
    """


class CR3BPPropagationTimeoutError(RuntimeError):
    """#652 defense-in-depth wall-clock backstop fired during STM-augmented
    CR3BP propagation, WITHOUT the close-secondary-encounter event
    (:class:`CR3BPCloseEncounterError`) having fired first.

    Subclasses ``RuntimeError`` for the same existing-caller-compatibility
    reason as :class:`CR3BPCloseEncounterError`. This is the coarser,
    secondary guard -- it exists in case some OTHER pathological case (not the
    near-secondary-encounter #651 diagnosed) also stalls the integrator; it is
    NOT the primary mechanism (that is the close-encounter event, which is
    diagnostic) and should be rare in practice. Its budget is generous (see
    :data:`_PROPAGATION_WALLCLOCK_BUDGET_S`) so it never fires for real work.
    """


# #652: distance from the SECONDARY body (nondimensional CR3BP length units --
# 1.0 is the secondary's own SMA about the primary) below which #651's bounded
# SIGALRM diagnostics confirmed the augmented state+STM DOP853 integration's
# adaptive step size collapses into a minutes+ grind. Empirically calibrated
# against #651's own confirmed pathological draw (Sun-Earth mu, model rng seed
# 651, n=100 draw index 77 -- see tests/core/test_cr3bp_close_encounter.py):
# that trajectory's true closest approach bottoms out at r2 ~= 2.18e-7, and a
# direct sweep of candidate thresholds on that exact case showed the grind
# ramps up steeply between r2 ~= 9e-6 (terminates in ~17ms) and r2 ~= 4e-6
# (already ~7s and climbing) -- i.e. the pathological regime lives BELOW
# ~1e-5. This threshold sits about two orders of magnitude above the observed
# grind onset (comfortable margin to fire fast and reliably) while staying
# roughly 4x below the smallest physically-realistic close approach this
# project's own systems registry supports (Sun-Earth LEO/surface grazing,
# r2 ~= 4.3e-5) and two-plus orders of magnitude below the tightest
# legitimate lunar flyby this project models (Earth-Moon LLO, r2 ~= 0.0048) --
# so it does not clip any known intentional close-approach trajectory.
_CLOSE_ENCOUNTER_R2_THRESHOLD = 1e-5

# #652 defense-in-depth backstop only (see CR3BPPropagationTimeoutError) --
# generous relative to a normal propagate() call (milliseconds), so it should
# never fire for real work; it exists purely to bound worst-case wall time for
# a pathological case the close-encounter event above does not anticipate.
_PROPAGATION_WALLCLOCK_BUDGET_S = 30.0


def _make_close_encounter_event(
    threshold: float = _CLOSE_ENCOUNTER_R2_THRESHOLD,
) -> Callable[[float, NDArray[np.float64], float], float]:
    """Build a ``solve_ivp`` terminal event: distance from the current state's
    position (``y[:3]`` -- valid whether ``y`` is a bare 6-state or a 6+36
    augmented state+STM vector) to the secondary body crossing below
    ``threshold``. See #652 / :class:`CR3BPCloseEncounterError`. Follows this
    project's existing ``events=``+``.terminal``/``.direction`` convention
    (see ``search/cr3bp_periodic.py``'s own ``_y_event``).
    """

    def event(t: float, y: NDArray[np.float64], mu: float) -> float:
        r2 = math.sqrt((float(y[0]) - 1.0 + mu) ** 2 + float(y[1]) ** 2 + float(y[2]) ** 2)
        return r2 - threshold

    event.terminal = True  # type: ignore[attr-defined]
    event.direction = -1.0  # type: ignore[attr-defined] # only fire approaching, not departing
    return event


def _make_wallclock_budget_event(
    budget_s: float = _PROPAGATION_WALLCLOCK_BUDGET_S,
) -> Callable[[float, NDArray[np.float64], float], float]:
    """Build a ``solve_ivp`` terminal event: #652 defense-in-depth wall-clock
    backstop. See :class:`CR3BPPropagationTimeoutError`. Stateful (captures
    its own start time at construction) -- build ONE instance per
    :func:`propagate` call and reuse it across every ``solve_ivp`` invocation
    within that call (the fixed-path mode makes several) so the budget covers
    the WHOLE call, not just one sub-interval.
    """
    start = time.monotonic()

    def event(t: float, y: NDArray[np.float64], mu: float) -> float:
        return budget_s - (time.monotonic() - start)

    event.terminal = True  # type: ignore[attr-defined]
    event.direction = -1.0  # type: ignore[attr-defined]
    return event


def _stm_termination_events() -> list[Callable[[float, NDArray[np.float64], float], float]]:
    """Build a fresh #652 event pair for one :func:`propagate` call (the
    events read ``mu`` from ``solve_ivp``'s own ``args`` at call time, not from
    here)."""
    return [_make_close_encounter_event(), _make_wallclock_budget_event()]


def _raise_for_stm_event_termination(
    sol: OdeResult[np.float64], *, mu: float, t_requested: float, context: str
) -> None:
    """If a #652 termination event fired, raise the specific, honestly-worded
    exception for it. A no-op if ``sol`` ran to completion or failed for the
    existing, unrelated ``sol.success=False`` reason (that case is still
    handled by each caller's own existing generic ``RuntimeError``, unchanged).
    Uses ``getattr`` (not direct attribute access) so a caller's own
    monkeypatched/stubbed ``solve_ivp`` replacement (e.g.
    ``test_propagate_raises_on_integrator_failure``'s ``_FailedSol``, which has
    no ``t_events`` attribute at all) degrades to a no-op here instead of an
    unrelated ``AttributeError`` -- a real ``scipy.integrate.OdeResult``
    always carries ``t_events`` (``None`` when no ``events=`` were passed).
    """
    t_events = getattr(sol, "t_events", None)
    if t_events is None:
        return
    if len(t_events) > 0 and len(t_events[0]) > 0:
        t_evt = float(t_events[0][0])
        y_events = sol.y_events
        assert y_events is not None
        y_evt = y_events[0][0]
        r2 = math.sqrt(
            (float(y_evt[0]) - 1.0 + mu) ** 2 + float(y_evt[1]) ** 2 + float(y_evt[2]) ** 2
        )
        raise CR3BPCloseEncounterError(
            f"{context}: terminated at t={t_evt:.6g} (of requested t={t_requested:.6g}) -- "
            f"distance to secondary crossed below {_CLOSE_ENCOUNTER_R2_THRESHOLD:.3g} "
            f"(nondimensional; reached {r2:.3g}). This is a genuine close encounter with "
            "the secondary body that develops DURING propagation, not visible from the "
            "initial condition alone -- the augmented state+STM variational equations "
            "become increasingly singular near the secondary and would otherwise collapse "
            "the adaptive step size into a minutes+ grind rather than a clean "
            "sol.success=False (#651, #652). Treat this as a precision-floor "
            "non-convergence for this initial condition, not a generic solver failure."
        )
    if len(t_events) > 1 and len(t_events[1]) > 0:
        t_evt = float(t_events[1][0])
        raise CR3BPPropagationTimeoutError(
            f"{context}: exceeded the {_PROPAGATION_WALLCLOCK_BUDGET_S:.0f}s #652 "
            f"defense-in-depth wall-clock budget at t={t_evt:.6g} (of requested "
            f"t={t_requested:.6g}) WITHOUT the close-secondary-encounter event firing "
            "first -- an unanticipated pathological case, not the near-secondary-encounter "
            "#651 diagnosed (that case is caught fast by CR3BPCloseEncounterError instead)."
        )


def _r1_r2(x: float, y: float, z: float, mu: float) -> tuple[float, float]:
    r1 = math.sqrt((x + mu) ** 2 + y * y + z * z)
    r2 = math.sqrt((x - 1.0 + mu) ** 2 + y * y + z * z)
    return r1, r2


def jacobi_constant(state6: NDArray[np.float64], mu: float) -> float:
    """C = (x^2+y^2) + 2(1-mu)/r1 + 2mu/r2 - v^2 (conserved)."""
    x, y, z, vx, vy, vz = (float(v) for v in state6)
    r1, r2 = _r1_r2(x, y, z, mu)
    return float(
        (x * x + y * y) + 2.0 * (1.0 - mu) / r1 + 2.0 * mu / r2 - (vx * vx + vy * vy + vz * vz)
    )


def jacobi_gap_dv_min(speed: float, delta_c: float) -> float:
    """Minimum single-impulse |dv| to change the Jacobi constant by ``delta_c``.

    An impulse at fixed position changes C only through the ``-v^2`` term, so
    ``delta_c = v0^2 - vf^2`` where ``vf`` is the post-impulse speed. Any
    impulse achieving the gap satisfies ``|dv| >= |vf - v0|`` (triangle
    inequality), with equality for a tangential burn — hence

        dv_min = |sqrt(speed^2 - delta_c) - speed|.

    This is the Jacobi-gap minimum-DV technique of Cuevas del Valle, Urrutxua &
    Solano-Lopez 2026 ("Fuel-optimal Rendezvous in the CR3BP via MPC and
    Proximal Operators", CEAS EuroGNC 2026, paper CEAS-GNC-2026-012, Sec. 7.2):
    the Jacobi energy integral floors the l2 DV-requirement to bridge the
    energy gap between two CR3BP states. The bound is exact (tight) for a
    single impulse applied at speed ``speed``; it is a rigorous lower bound for
    any single-impulse transfer-cost estimate at that state. Units: any
    consistent set — ``delta_c`` carries speed^2 units; nondimensional in our
    usage (multiply by ``l_km / t_s`` to dimensionalise the result).

    Raises
    ------
    ValueError
        If ``speed`` is negative, or ``delta_c > speed**2`` (the required
        post-impulse speed squared would be negative: no single impulse at this
        state can raise C past the v = 0 ceiling).
    """
    if speed < 0.0:
        raise ValueError(f"jacobi_gap_dv_min: negative speed {speed}")
    vf_sq = speed * speed - delta_c
    if vf_sq < 0.0:
        raise ValueError(
            f"jacobi_gap_dv_min: delta_c={delta_c} exceeds speed^2={speed * speed}; "
            "C cannot be raised past the zero-velocity ceiling by one impulse"
        )
    return abs(math.sqrt(vf_sq) - speed)


def cr3bp_eom(t: float, state6: NDArray[np.float64], mu: float) -> NDArray[np.float64]:
    """Rotating-frame equations of motion (state [x,y,z,vx,vy,vz])."""
    x, y, z, vx, vy, vz = (float(v) for v in state6)
    r1, r2 = _r1_r2(x, y, z, mu)
    r1c, r2c = r1**3, r2**3
    ax = x + 2.0 * vy - (1.0 - mu) * (x + mu) / r1c - mu * (x - 1.0 + mu) / r2c
    ay = y - 2.0 * vx - (1.0 - mu) * y / r1c - mu * y / r2c
    az = -(1.0 - mu) * z / r1c - mu * z / r2c
    return np.array([vx, vy, vz, ax, ay, az], dtype=np.float64)


@dataclass(frozen=True)
class CR3BPSystem:
    mu: float
    primary: str
    secondary: str
    l_km: float  # characteristic length (secondary SMA about primary)
    t_s: float  # characteristic time = sqrt(l^3 / (G(m1+m2)))


@dataclass(frozen=True)
class CR3BPArc:
    state_f: NDArray[np.float64]
    stm: NDArray[np.float64] | None
    t: float


def cr3bp_stm_eom(t: float, y42: NDArray[np.float64], mu: float) -> NDArray[np.float64]:
    """State (6) + flattened 6x6 STM (36) variational EOM.

    Implements the classical variational equations: Phi'(t) = A(X(t)) Phi(t)
    where A is the Jacobian of the CR3BP RHS at the current state X(t).

    Pellegrini-Russell 2016 caveat
    ------------------------------
    The variational equations capture the smooth-flow sensitivity but, when
    integrated alongside the state by a *variable-step* integrator, miss a
    per-step term

        f(X_i) (partial delta_t / partial X|_X_i) Phi_i

    (eq. 17, Pellegrini & Russell, JGCD 2016, DOI 10.2514/1.G001920, p. 3) --
    the step size delta_t selected by the adaptive controller depends on the
    initial condition through the local error estimate, and this dependence
    contributes to the STM but is not propagated by the variational equations.
    The bias is small (~ epsilon^(4/5) for an RKF(7)8 controller per eq. 19,
    same page) but is most pronounced for *highly sensitive* orbits where the
    STM has large magnitude. See :func:`propagate` for the ``stm_mode``
    parameter that switches to a pre-scheduled fixed-path replay (Pellegrini-
    Russell Conclusion 2, p. 14) -- it pins the step grid to the state-only
    pass so the IC-dependence vanishes.

    Reference
    ---------
    Pellegrini, E. & Russell, R.P. (2016), "On the Computation and Accuracy of
    Trajectory State Transition Matrices", *Journal of Guidance, Control, and
    Dynamics*, DOI 10.2514/1.G001920.
    """
    s = y42[:6]
    phi = y42[6:].reshape(6, 6)
    x, y, z = float(s[0]), float(s[1]), float(s[2])
    r1, r2 = _r1_r2(x, y, z, mu)
    r1c, r2c = r1**3, r2**3
    r1f, r2f = r1**5, r2**5
    om1 = 1.0 - mu
    # Pseudo-potential second derivatives (Uxx etc.) for the A matrix.
    uxx = (
        1 - om1 / r1c - mu / r2c + 3 * om1 * (x + mu) ** 2 / r1f + 3 * mu * (x - 1 + mu) ** 2 / r2f
    )
    uyy = 1 - om1 / r1c - mu / r2c + 3 * om1 * y * y / r1f + 3 * mu * y * y / r2f
    uzz = -om1 / r1c - mu / r2c + 3 * om1 * z * z / r1f + 3 * mu * z * z / r2f
    uxy = 3 * om1 * (x + mu) * y / r1f + 3 * mu * (x - 1 + mu) * y / r2f
    uxz = 3 * om1 * (x + mu) * z / r1f + 3 * mu * (x - 1 + mu) * z / r2f
    uyz = 3 * om1 * y * z / r1f + 3 * mu * y * z / r2f
    mat_a = np.zeros((6, 6))
    mat_a[0:3, 3:6] = np.eye(3)
    mat_a[3, 0], mat_a[3, 1], mat_a[3, 2] = uxx, uxy, uxz
    mat_a[4, 0], mat_a[4, 1], mat_a[4, 2] = uxy, uyy, uyz
    mat_a[5, 0], mat_a[5, 1], mat_a[5, 2] = uxz, uyz, uzz
    mat_a[3, 4], mat_a[4, 3] = 2.0, -2.0  # Coriolis
    dphi = mat_a @ phi
    return np.concatenate([cr3bp_eom(t, s, mu), dphi.reshape(36)])


def propagate(
    system: CR3BPSystem,
    state6: NDArray[np.float64],
    t: float,
    *,
    with_stm: bool = False,
    rtol: float = 1e-12,
    atol: float = 1e-12,
    stm_mode: StmMode = "variable",
) -> CR3BPArc:
    """Propagate a state (and optionally the STM) for nondimensional time ``t``.

    Parameters
    ----------
    system, state6, t, with_stm, rtol, atol :
        Standard CR3BP propagation inputs. ``rtol`` / ``atol`` are the DOP853
        local-error tolerances (default 1e-12 / 1e-12 -- the project standard).
    stm_mode :
        Selects the STM-integration method when ``with_stm=True``. Ignored
        otherwise.

        * ``"variable"`` (default; backward-compatible): integrate the augmented
          state + STM together with scipy's variable-step DOP853. This is the
          legacy path. Subject to the Pellegrini-Russell 2016 small-bias
          caveat: the adaptive step size delta_t depends on the IC through the
          local error estimate, contributing a term

              f(X_i) (partial delta_t / partial X|_X_i) Phi_i

          (eq. 17, p. 3) that the variational equations do NOT capture. The
          bias scales as ~ epsilon^(4/5) for an RKF(7)8 controller (eq. 19,
          p. 3) and is most pronounced for highly sensitive orbits with large
          ||Phi||.

        * ``"fixed_path"``: implements Pellegrini-Russell's "fixed-path"
          mitigation (§III.A.2, p. 6; Conclusion 2, p. 14). Step 1: run a
          state-only DOP853 propagation at the same tolerances, *recording*
          its accepted step grid ``[t_0, t_1, ..., t_N=t]``. Step 2: replay
          the augmented state + STM along that pre-scheduled grid, taking
          exactly one DOP853 step per sub-interval (h_i = t_{i+1} - t_i fed
          as ``first_step=max_step=h_i``). Because the grid is fixed in
          advance from the unperturbed trajectory, ``partial delta_t /
          partial X = 0``, the eq. 17 contamination term vanishes, and the
          STM is unbiased relative to the achieved state path. Cost: a
          state-only pass plus a per-sub-interval solver call (~3x the
          variable-step path's wall time for our typical orbits).

    Returns
    -------
    CR3BPArc :
        Always carries ``state_f`` at time ``t``; ``stm`` is the 6x6 monodromy
        block when ``with_stm=True``, else ``None``.

    Raises
    ------
    RuntimeError
        If the integrator fails (``sol.success`` is False) -- e.g. a collision
        trajectory driving the step size below floating-point resolution near a
        primary. Returning ``sol.y[:, -1]`` in that case would silently hand back
        the state where the integrator gave up, not the state at time ``t``.
    CR3BPCloseEncounterError (with_stm=True only; subclasses RuntimeError)
        #652: a close-secondary-encounter ``solve_ivp`` termination event
        fired -- the trajectory dynamically evolved into a genuine close
        approach to the SECONDARY body partway through STM-augmented
        propagation (undetectable from the initial condition alone), which
        would otherwise collapse the adaptive step size into a minutes+ grind
        without ever reaching ``sol.success=False`` (the failure mode above is
        near the PRIMARY and unrelated). Fires fast, before the grind, by
        design; see the class docstring for detail.
    CR3BPPropagationTimeoutError (with_stm=True only; subclasses RuntimeError)
        #652 defense-in-depth wall-clock backstop -- fired WITHOUT the
        close-encounter event above also firing, i.e. an unanticipated
        pathological case. Should be rare; see the class docstring.
    ValueError
        If ``stm_mode`` is not one of the supported literals.

    Reference
    ---------
    Pellegrini, E. & Russell, R.P. (2016), "On the Computation and Accuracy of
    Trajectory State Transition Matrices", *Journal of Guidance, Control, and
    Dynamics*, DOI 10.2514/1.G001920. The eq. 17 derivation (p. 3) plus the
    PO2-test-case Fig. 13 (p. 13) demonstrate that variable-step variational
    can underperform CSD by several orders of magnitude on highly sensitive
    orbits at loose tolerance; fixed-path closes most of that gap without a
    complex-number integrator refactor.
    """
    if stm_mode not in ("variable", "fixed_path"):
        raise ValueError(
            f"propagate: stm_mode must be 'variable' or 'fixed_path', got {stm_mode!r}"
        )
    if with_stm:
        if stm_mode == "variable":
            return _propagate_with_stm_variable(system, state6, t, rtol=rtol, atol=atol)
        return _propagate_with_stm_fixed_path(system, state6, t, rtol=rtol, atol=atol)
    sol = solve_ivp(
        cr3bp_eom,
        (0.0, t),
        np.asarray(state6, float),
        args=(system.mu,),
        rtol=rtol,
        atol=atol,
        method="DOP853",
    )
    if not sol.success:
        raise RuntimeError(f"CR3BP propagation failed at t={sol.t[-1]}: {sol.message}")
    return CR3BPArc(state_f=sol.y[:, -1], stm=None, t=t)


def _propagate_with_stm_variable(
    system: CR3BPSystem,
    state6: NDArray[np.float64],
    t: float,
    *,
    rtol: float,
    atol: float,
) -> CR3BPArc:
    """Legacy variable-step augmented (state+STM) DOP853 propagation.

    #652: guarded by the close-secondary-encounter + wall-clock-budget
    termination events (see :func:`_stm_termination_events`) -- see
    :class:`CR3BPCloseEncounterError` / :class:`CR3BPPropagationTimeoutError`
    for what each raises and why.
    """
    y0 = np.concatenate([np.asarray(state6, float), np.eye(6).reshape(36)])
    sol = solve_ivp(
        cr3bp_stm_eom,
        (0.0, t),
        y0,
        args=(system.mu,),  # type: ignore[call-overload]
        rtol=rtol,
        atol=atol,
        method="DOP853",
        dense_output=False,
        events=_stm_termination_events(),
    )
    _raise_for_stm_event_termination(
        sol, mu=system.mu, t_requested=t, context="CR3BP STM propagation (variable-step)"
    )
    if not sol.success:
        raise RuntimeError(f"CR3BP STM propagation failed at t={sol.t[-1]}: {sol.message}")
    yf = sol.y[:, -1]
    return CR3BPArc(state_f=yf[:6], stm=yf[6:].reshape(6, 6), t=t)


def _propagate_with_stm_fixed_path(
    system: CR3BPSystem,
    state6: NDArray[np.float64],
    t: float,
    *,
    rtol: float,
    atol: float,
) -> CR3BPArc:
    """Pellegrini-Russell fixed-path STM propagation.

    Step 1: state-only DOP853 to record the accepted-step grid. Step 2: replay
    the augmented (state, STM) system along that recorded grid, one DOP853
    step per sub-interval, with ``first_step = max_step = h_i`` (and the
    same ``rtol`` / ``atol`` as the state-only pass).

    Pinning the step grid in advance to the unperturbed state trajectory kills
    the ``partial delta_t / partial X`` term in eq. 17 (Pellegrini-Russell
    2016, p. 3): the per-step size is now an INPUT, not an IC-dependent
    output. The cost is a separate state-only pass and per-sub-interval
    solve_ivp instantiation; the benefit is an unbiased STM aligned exactly
    with the state path.

    The single-step DOP853 (no internal subdivision) gives the most faithful
    reading of Pellegrini-Russell §III.A.2 ("fixed-path feature"): the same
    8th-order scheme that produced the recorded grid replays each step
    exactly once at the same h_i, so the LOCAL truncation error per step is
    consistent with the original state-only pass.

    #652: this path is guarded by the SAME close-secondary-encounter +
    wall-clock-budget termination events as
    :func:`_propagate_with_stm_variable` (see :func:`_stm_termination_events`),
    applied to BOTH the state-only pre-pass and every per-sub-interval
    augmented replay call, sharing one event pair (and therefore one combined
    wall-clock budget) across the whole call. This mode is EXPOSED TO THE SAME
    #651-diagnosed risk despite its state-only pre-pass: empirically, the
    state-only pass alone does not hang on #651's confirmed pathological draw
    (it reaches the same close approach in well under a second, since the
    state equations' own singularity is only ``~1/r2^2``), but the
    per-sub-interval AUGMENTED replay in step 2 does, on that exact draw --
    confirmed directly (not inferred) before adding this guard. Sharing the
    close-encounter event with the cheap state-only pass actually catches
    most pathological cases EARLIER (and more cheaply) than variable-step
    mode would, before step 2's expensive replay is ever reached.
    """
    state_arr = np.asarray(state6, dtype=np.float64)
    if state_arr.shape != (6,):
        raise ValueError(
            "propagate(stm_mode='fixed_path'): state6 must be a 6-vector, "
            f"got shape {state_arr.shape}"
        )
    events = _stm_termination_events()
    # Step 1: record state-only step grid.
    sol_state = solve_ivp(
        cr3bp_eom,
        (0.0, t),
        state_arr,
        args=(system.mu,),  # type: ignore[call-overload]
        rtol=rtol,
        atol=atol,
        method="DOP853",
        events=events,
    )
    _raise_for_stm_event_termination(
        sol_state, mu=system.mu, t_requested=t, context="CR3BP fixed-path state-only pre-pass"
    )
    if not sol_state.success:
        raise RuntimeError(
            "CR3BP fixed-path state-only pre-pass failed at "
            f"t={sol_state.t[-1]}: {sol_state.message}"
        )
    grid = np.asarray(sol_state.t, dtype=np.float64)
    if grid.size < 2:
        raise RuntimeError(
            "CR3BP fixed-path: state-only DOP853 returned a degenerate grid "
            f"(size {grid.size}); cannot replay STM"
        )
    # Step 2: replay augmented system along the recorded grid.
    y = np.concatenate([state_arr, np.eye(6).reshape(36)])
    for i in range(grid.size - 1):
        h = float(grid[i + 1] - grid[i])
        if h <= 0.0:
            raise RuntimeError(f"CR3BP fixed-path: non-positive recorded step h={h} at index {i}")
        sol_step = solve_ivp(
            cr3bp_stm_eom,
            (grid[i], grid[i + 1]),
            y,
            args=(system.mu,),  # type: ignore[call-overload]
            rtol=rtol,
            atol=atol,
            method="DOP853",
            first_step=h,
            max_step=h,
            events=events,
        )
        _raise_for_stm_event_termination(
            sol_step,
            mu=system.mu,
            t_requested=t,
            context=f"CR3BP fixed-path STM replay at sub-interval [{grid[i]}, {grid[i + 1]}]",
        )
        if not sol_step.success:
            raise RuntimeError(
                "CR3BP fixed-path STM replay failed at sub-interval "
                f"[{grid[i]}, {grid[i + 1]}]: {sol_step.message}"
            )
        y = sol_step.y[:, -1]
    return CR3BPArc(state_f=y[:6], stm=y[6:].reshape(6, 6), t=t)


def cr3bp_system(primary: str, secondary: str) -> CR3BPSystem:
    """Build a CR3BPSystem from the registry: mu = GM2/G(m1+m2), scales from SMA.

    ``PRIMARIES[primary]`` is the JPL *system* GM (gm_de440 planetary constants;
    see ``core/satellites.py``) — it ALREADY includes the satellites' GM. The
    pair total ``G(m1+m2)`` is therefore the system GM itself; the historical
    ``PRIMARIES[primary] + gm2`` double-counted the secondary (#212 Part A:
    Earth-Moon mu came out 0.0120047, -1.2% vs the canonical 0.0121505843,
    and t_s was -0.6%). For Earth-Moon the system GM is exactly G(Earth+Moon);
    for the Jupiter/Saturn pairs the OTHER moons' mass stays folded into the
    primary term (<= 2.4e-4 of the system GM, dominated by Titan/the Galileans)
    — the best decomposition available from system GMs, and far below the
    SILVER members' reported precision (quantified in tests/core/test_cr3bp.py).

    Sun-primary heliocentric pairs (e.g. ``cr3bp_system("Sun", "Earth")``) are
    served from the planet registry (``core.constants.PLANETS``) instead of the
    moon registry: the secondary's GM and heliocentric SMA come from ``PLANETS``
    and the primary GM is ``MU_SUN_KM3_S2`` (IAU 2015 nominal solar GM). The
    convention is identical — pair GM = G(m1+m2), ``l_km`` = secondary SMA about
    primary, ``t_s`` = sqrt(l_km^3 / G(m1+m2)). This is the single source of
    truth for the #405/#411 cross-system (Sun-Earth <-> Earth-Moon) work.
    """
    if primary == "Sun":
        planet = _planet_by_name(secondary)
        gm2 = planet.mu_km3_s2
        gm_pair = MU_SUN_KM3_S2 + gm2  # G(Sun + planet)
        mu = gm2 / gm_pair
        l_km = planet.sma_au * AU_KM
        t_s = float(math.sqrt(l_km**3 / gm_pair))
        return CR3BPSystem(mu=mu, primary=primary, secondary=secondary, l_km=l_km, t_s=t_s)
    gm2 = SATELLITES[secondary].mu_km3_s2
    gm_pair = PRIMARIES[primary]  # system GM == G(m1 + m2 [+ other moons])
    mu = gm2 / gm_pair
    l_km = SATELLITES[secondary].sma_km
    t_s = float(math.sqrt(l_km**3 / gm_pair))
    return CR3BPSystem(mu=mu, primary=primary, secondary=secondary, l_km=l_km, t_s=t_s)


def _planet_by_name(name: str) -> PlanetData:
    """Resolve a full English planet name (e.g. ``"Earth"``) to its PLANETS entry.

    ``PLANETS`` is keyed by short code (``"E"``, ``"M"``, ...); cross-system callers
    name the secondary by its English name, matching the ``(primary, secondary)``
    string convention of the moon-registry path. Raises KeyError with the valid
    names if no match (parallels ``SATELLITES[secondary]`` failing for a bad moon).
    """
    for planet in PLANETS.values():
        if planet.name == name:
            return planet
    valid = sorted(p.name for p in PLANETS.values())
    raise KeyError(f"no Sun-primary planet named {name!r}; valid secondaries: {valid}")
