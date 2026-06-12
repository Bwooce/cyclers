"""Natural-parameter (Jacobi) continuation of CR3BP symmetric periodic orbits.

Spec: ``docs/superpowers/specs/2026-06-12-cr3bp-family-continuation-search-design.md``
(Phase 1). Walks a family of perpendicular-x-axis-crossing symmetric orbits by
stepping the Jacobi constant ``C`` and re-converging the fixed-Jacobi symmetric
single-shooting corrector (Ross & Roberts-Tsoukkas 2025,
``correct_symmetric_fixed_jacobi``) at each step. The ``x0`` of the next member is
*predicted* by a secant extrapolation from the last two converged members
(natural-parameter / pseudo-tangent in ``(C, x0)``), then corrected.

Every kept member must pass the full convergence gauntlet (spec constraint 3):
  - corrector converged (perpendicular re-crossing residual < ``tol``);
  - period in bounds (no collapse, no run-past): ``period_floor <= T`` and
    ``T <= period_ceiling`` relative to the seed period;
  - not an equilibrium: ``max|v|`` over the propagated period >= ``max_speed_floor``
    AND position amplitude >= ``amplitude_floor`` (a libration point trivially
    closes for any period -- the "it closed!" danger signal);
  - distinct from every prior member on this branch (dedup on ``(C, x0)``);
  - Jacobi conserved over the period: ``|C(T) - C(0)| <= jacobi_tol``;
  - independent-Radau cross-check (``crosscheck_periodic``) agrees (re-closure +
    Jacobi drift) -- a different integrator than the DOP853 used by the corrector.

Stability ``nu = 1/2 (lambda + 1/lambda)`` is computed for every kept member via
the Barden half-period monodromy (``barden_stability``); ``|nu| < 1`` is linearly
stable.

The walk STOPS at the first of: the Jacobi bound (``min_jacobi`` / ``max_jacobi``);
a turning point / fold (the corrector fails to converge, the Jacobi radicand goes
negative -- ``ydot0_from_jacobi`` raises -- or the predicted/corrected step
reverses sign in ``x0``, i.e. the branch folded back over itself); or a gauntlet
rejection that is not recoverable by re-prediction. Each stop reason is recorded.

Model: pure planar CR3BP (PCR3BP). Pure (math/numpy/scipy + the CR3BP core).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import StrEnum

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import solve_ivp

import cyclerfinder.core.cr3bp as cr3bp
import cyclerfinder.search.cr3bp_periodic as cp


class StopReason(StrEnum):
    """Why a continuation branch terminated (recorded per direction)."""

    JACOBI_BOUND = "jacobi_bound"  # reached min_jacobi / max_jacobi
    MAX_STEPS = "max_steps"  # n_steps exhausted (not a physical edge)
    FOLD_RADICAND = "fold_radicand"  # ydot0_from_jacobi radicand < 0 (turning point)
    FOLD_REVERSAL = "fold_reversal"  # x0 step reversed sign (branch folded back)
    NO_CONVERGE = "no_converge"  # corrector failed at the stepped C
    GAUNTLET_REJECT = "gauntlet_reject"  # converged but failed the gauntlet
    TOPOLOGY_JUMP = "topology_jump"  # period jumped vs the previous member (off-family)
    SEED_OFF_FAMILY = "seed_off_family"  # the corrected seed does not match its known period


@dataclass(frozen=True)
class BranchMember:
    """One converged, gauntlet-passing member of a continuation branch.

    All quantities are in-model (pure CR3BP). ``state0`` is the perpendicular
    x-axis-crossing IC ``(x0, 0, 0, 0, ydot0, 0)`` (nondimensional).
    """

    state0: NDArray[np.float64]
    x0: float
    ydot0: float
    period: float
    jacobi: float
    nu: float
    abs_lambda: float
    crossing_residual: float  # |xdot(t_half)| from the corrector
    radau_djacobi: float  # independent-Radau Jacobi drift over the period
    max_speed_nd: float  # equilibrium gate diagnostic
    amplitude_nd: float  # equilibrium gate diagnostic
    n_iter: int
    stable: bool  # |nu| < 1


@dataclass
class FamilyBranch:
    """The ordered branch produced by continuing one seed in one C-direction."""

    seed_label: str
    direction: int  # +1 (increasing C) or -1 (decreasing C)
    d_jacobi: float
    half_crossings: int
    ydot0_sign: float
    members: list[BranchMember] = field(default_factory=list)
    stop_reason: StopReason = StopReason.MAX_STEPS
    n_steps_taken: int = 0
    n_rejected: int = 0


# Default gauntlet thresholds (nondimensional). Mirror the #182 moontour gates
# (max|v| and amplitude floors, dedup tol) and the spec's Jacobi-conservation
# bound (<= 1e-10).
MAX_SPEED_FLOOR_ND = 1e-6
AMPLITUDE_FLOOR_ND = 1e-6
JACOBI_CONSERVATION_TOL = 1e-10
DEDUP_X0_TOL = 1e-9
DEDUP_C_TOL = 1e-12


def _gate_metrics(
    system: cr3bp.CR3BPSystem,
    state0: NDArray[np.float64],
    period: float,
    *,
    rtol: float,
    atol: float,
) -> tuple[float, float]:
    """max|v| and position-amplitude over one propagated period (equilibrium gate)."""
    # Sample the period to capture the extrema (endpoints alone miss the swing).
    n = 128
    sol = solve_ivp(
        cr3bp.cr3bp_eom,
        (0.0, period),
        np.asarray(state0, float),
        args=(system.mu,),
        method="DOP853",
        rtol=rtol,
        atol=atol,
        t_eval=np.linspace(0.0, period, n),
    )
    if not sol.success:
        raise RuntimeError(f"equilibrium-gate propagation failed: {sol.message}")
    states = sol.y
    speeds = np.linalg.norm(states[3:6, :], axis=0)
    disp = np.linalg.norm(states[0:3, :] - states[0:3, :1], axis=0)
    return float(np.max(speeds)), float(np.max(disp))


def _run_gauntlet(
    system: cr3bp.CR3BPSystem,
    orbit: cp.SymmetricOrbit,
    *,
    period_floor: float,
    period_ceiling: float,
    max_speed_floor: float,
    amplitude_floor: float,
    jacobi_tol: float,
    radau_closure_tol: float,
    radau_jacobi_tol: float,
    rtol: float,
    atol: float,
) -> tuple[bool, str, BranchMember | None]:
    """Run the full convergence gauntlet on a corrected symmetric orbit.

    Returns ``(ok, reason, member)``: on pass ``ok=True`` and ``member`` is the
    populated :class:`BranchMember`; on rejection ``ok=False`` and ``reason``
    names the failing gate (``member`` is ``None``).
    """
    if not orbit.converged:
        return False, "not_converged", None
    if not (period_floor <= orbit.period <= period_ceiling):
        return False, "period_bounds", None

    state0 = np.array([orbit.x0, 0.0, 0.0, 0.0, orbit.ydot0, 0.0])

    # Equilibrium gate (max|v| + position amplitude over the period).
    max_speed, amplitude = _gate_metrics(system, state0, orbit.period, rtol=rtol, atol=atol)
    if max_speed < max_speed_floor or amplitude < amplitude_floor:
        return False, "equilibrium", None

    # Jacobi conservation over the period (in-model invariant).
    arc = cr3bp.propagate(system, state0, orbit.period, with_stm=False, rtol=rtol, atol=atol)
    dj = abs(cr3bp.jacobi_constant(arc.state_f, system.mu) - orbit.jacobi)
    if dj > jacobi_tol:
        return False, "jacobi_conservation", None

    # Independent-Radau cross-check (different integrator than DOP853).
    po = cp.PeriodicOrbit(
        state0=state0,
        period=orbit.period,
        jacobi=orbit.jacobi,
        converged=orbit.converged,
        closure_residual=orbit.crossing_residual,
    )
    radau_ok, radau_dj = cp.crosscheck_periodic(
        system, po, closure_tol=radau_closure_tol, jacobi_tol=radau_jacobi_tol
    )
    if not radau_ok:
        return False, "radau_crosscheck", None

    nu, lam = cp.barden_stability(system, orbit, rtol=rtol, atol=atol)
    member = BranchMember(
        state0=state0,
        x0=orbit.x0,
        ydot0=orbit.ydot0,
        period=orbit.period,
        jacobi=orbit.jacobi,
        nu=float(nu),
        abs_lambda=float(abs(lam)),
        crossing_residual=orbit.crossing_residual,
        radau_djacobi=float(radau_dj),
        max_speed_nd=max_speed,
        amplitude_nd=amplitude,
        n_iter=orbit.n_iter,
        stable=abs(float(nu)) < 1.0,
    )
    return True, "", member


def continue_family(
    system: cr3bp.CR3BPSystem,
    seed: cp.SymmetricOrbit,
    *,
    direction: int,
    d_jacobi: float,
    n_steps: int,
    min_jacobi: float,
    max_jacobi: float,
    half_crossings: int,
    ydot0_sign: float,
    seed_label: str = "",
    corrector_tol: float = 1e-10,
    period_floor_frac: float = 0.5,
    period_ceiling_frac: float = 2.0,
    period_step_frac: float = 0.10,
    max_speed_floor: float = MAX_SPEED_FLOOR_ND,
    amplitude_floor: float = AMPLITUDE_FLOOR_ND,
    jacobi_tol: float = JACOBI_CONSERVATION_TOL,
    radau_closure_tol: float = 1e-3,
    radau_jacobi_tol: float = 1e-8,
    dedup_x0_tol: float = DEDUP_X0_TOL,
    dedup_c_tol: float = DEDUP_C_TOL,
    rtol: float = 1e-12,
    atol: float = 1e-12,
) -> FamilyBranch:
    """Natural-parameter (Jacobi) continuation of one symmetric-orbit family.

    Starting from ``seed`` (a converged
    :class:`~cyclerfinder.search.cr3bp_periodic.SymmetricOrbit`),
    step the Jacobi constant by ``direction * abs(d_jacobi)`` and re-converge the
    fixed-Jacobi symmetric corrector at each step, predicting the next ``x0`` by a
    secant extrapolation from the last two converged members. Each member is run
    through the full gauntlet (closure / period-bounds / equilibrium / dedup /
    Jacobi-conservation / independent-Radau) and assigned a Barden ``nu``.

    Parameters
    ----------
    direction:
        ``+1`` walks toward larger ``C``, ``-1`` toward smaller ``C``.
    d_jacobi:
        Jacobi step magnitude (its sign is taken from ``direction``).
    n_steps:
        Maximum number of continuation steps (no silent cap -- the count and the
        stop reason are returned).
    min_jacobi, max_jacobi:
        Hard Jacobi bounds; the walk stops when the next ``C`` would leave
        ``[min_jacobi, max_jacobi]``.
    half_crossings, ydot0_sign:
        The crossing-index and velocity sign that place the corrector on this
        family's branch (held fixed across the walk -- they identify the family).
    period_floor_frac, period_ceiling_frac:
        Period bounds as fractions of the seed period (collapse / run-past gate).
    period_step_frac:
        Period-continuity gate: a continuation member whose period differs from the
        PREVIOUS member by more than this fraction is a topology jump (the corrector
        landed on a different perpendicular crossing / a different family), not a
        continuation step -- the branch stops (``TOPOLOGY_JUMP``). A genuine family
        varies its period smoothly with C; a near-doubling is the "it closed!"
        danger signal in disguise.

    Returns
    -------
    FamilyBranch
        The ordered list of gauntlet-passing members (seed first), the stop
        reason, the number of steps taken, and the rejection count.
    """
    if direction not in (1, -1):
        raise ValueError(f"direction must be +1 or -1, got {direction}")
    step = float(direction) * abs(float(d_jacobi))
    seed_period = float(seed.period)
    period_floor = period_floor_frac * seed_period
    period_ceiling = period_ceiling_frac * seed_period

    branch = FamilyBranch(
        seed_label=seed_label,
        direction=direction,
        d_jacobi=step,
        half_crossings=half_crossings,
        ydot0_sign=ydot0_sign,
    )

    # Seed member: it must itself pass the gauntlet (the seed is a sourced member,
    # but we validate it through the same gates so the branch is uniform).
    ok, _reason, seed_member = _run_gauntlet(
        system,
        seed,
        period_floor=period_floor,
        period_ceiling=period_ceiling,
        max_speed_floor=max_speed_floor,
        amplitude_floor=amplitude_floor,
        jacobi_tol=jacobi_tol,
        radau_closure_tol=radau_closure_tol,
        radau_jacobi_tol=radau_jacobi_tol,
        rtol=rtol,
        atol=atol,
    )
    if ok and seed_member is not None:
        branch.members.append(seed_member)

    # (C, x0, T) history for the secant predictor + period-continuity gate.
    c_hist: list[float] = [float(seed.jacobi)]
    x0_hist: list[float] = [float(seed.x0)]
    prev_period = float(seed.period)
    last_x0_delta = 0.0  # for the fold-reversal test

    c_curr = float(seed.jacobi)
    for k in range(1, n_steps + 1):
        c_next = c_curr + step
        if c_next < min_jacobi or c_next > max_jacobi:
            branch.stop_reason = StopReason.JACOBI_BOUND
            break

        # Predict x0 at c_next by secant extrapolation in (C, x0).
        if len(c_hist) >= 2 and abs(c_hist[-1] - c_hist[-2]) > 0.0:
            slope = (x0_hist[-1] - x0_hist[-2]) / (c_hist[-1] - c_hist[-2])
            x0_pred = x0_hist[-1] + slope * (c_next - c_hist[-1])
        else:
            x0_pred = x0_hist[-1]

        # Correct at the stepped C. A negative Jacobi radicand at the predicted
        # x0 is a turning point (the family does not extend past here).
        try:
            orbit = cp.correct_symmetric_fixed_jacobi(
                system,
                x0_pred,
                c_next,
                seed_period,
                ydot0_sign=ydot0_sign,
                half_crossings=half_crossings,
                tol=corrector_tol,
                rtol=rtol,
                atol=atol,
            )
        except ValueError:
            branch.stop_reason = StopReason.FOLD_RADICAND
            break

        branch.n_steps_taken = k

        if not orbit.converged:
            branch.stop_reason = StopReason.NO_CONVERGE
            break

        # Period-continuity gate: a genuine continuation member's period varies
        # smoothly with C. A jump (the corrector landed on a different
        # perpendicular crossing => a different family/topology) is rejected; the
        # branch stops rather than emitting the off-family orbit.
        if abs(orbit.period - prev_period) > period_step_frac * prev_period:
            branch.stop_reason = StopReason.TOPOLOGY_JUMP
            break

        # Dedup against every prior member on this branch (same (C, x0) = one orbit).
        is_dup = any(
            abs(orbit.x0 - m.x0) < dedup_x0_tol and abs(orbit.jacobi - m.jacobi) < dedup_c_tol
            for m in branch.members
        )
        if is_dup:
            # The corrector slid back onto an existing member: the branch folded.
            branch.stop_reason = StopReason.FOLD_REVERSAL
            break

        # Fold-reversal: the corrected x0 stepped opposite to the running trend.
        x0_delta = orbit.x0 - x0_hist[-1]
        if last_x0_delta != 0.0 and x0_delta * last_x0_delta < 0.0:
            branch.stop_reason = StopReason.FOLD_REVERSAL
            break

        ok, _reason, member = _run_gauntlet(
            system,
            orbit,
            period_floor=period_floor,
            period_ceiling=period_ceiling,
            max_speed_floor=max_speed_floor,
            amplitude_floor=amplitude_floor,
            jacobi_tol=jacobi_tol,
            radau_closure_tol=radau_closure_tol,
            radau_jacobi_tol=radau_jacobi_tol,
            rtol=rtol,
            atol=atol,
        )
        if not ok or member is None:
            branch.n_rejected += 1
            branch.stop_reason = StopReason.GAUNTLET_REJECT
            break

        branch.members.append(member)
        c_hist.append(orbit.jacobi)
        x0_hist.append(orbit.x0)
        prev_period = orbit.period
        last_x0_delta = x0_delta
        c_curr = c_next
    else:
        branch.stop_reason = StopReason.MAX_STEPS

    return branch
