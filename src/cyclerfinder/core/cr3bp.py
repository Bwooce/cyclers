"""Circular-restricted three-body problem (CR3BP) dynamics core (spec 2026-06-10).

Nondimensional rotating frame; mu = m2/(m1+m2); primary at (-mu,0,0), secondary at
(1-mu,0,0). See the plan's "CR3BP equations" block. Pure: math/numpy/scipy +
core.satellites only.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Literal

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import solve_ivp

from cyclerfinder.core.satellites import PRIMARIES, SATELLITES

# STM integration modes. See ``propagate(with_stm=True, stm_mode=...)`` for the
# semantics; ``"variable"`` is the legacy variable-step variational path,
# ``"fixed_path"`` is the Pellegrini-Russell 2016 mitigation (record the state
# only path's accepted step grid, replay state+STM along that pre-scheduled
# grid so the step size has no IC dependence).
StmMode = Literal["variable", "fixed_path"]


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
    """Legacy variable-step augmented (state+STM) DOP853 propagation."""
    y0 = np.concatenate([np.asarray(state6, float), np.eye(6).reshape(36)])
    sol = solve_ivp(
        cr3bp_stm_eom,
        (0.0, t),
        y0,
        args=(system.mu,),
        rtol=rtol,
        atol=atol,
        method="DOP853",
        dense_output=False,
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
    """
    state_arr = np.asarray(state6, dtype=np.float64)
    if state_arr.shape != (6,):
        raise ValueError(
            "propagate(stm_mode='fixed_path'): state6 must be a 6-vector, "
            f"got shape {state_arr.shape}"
        )
    # Step 1: record state-only step grid.
    sol_state = solve_ivp(
        cr3bp_eom,
        (0.0, t),
        state_arr,
        args=(system.mu,),
        rtol=rtol,
        atol=atol,
        method="DOP853",
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
            args=(system.mu,),
            rtol=rtol,
            atol=atol,
            method="DOP853",
            first_step=h,
            max_step=h,
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
    """
    gm2 = SATELLITES[secondary].mu_km3_s2
    gm_pair = PRIMARIES[primary]  # system GM == G(m1 + m2 [+ other moons])
    mu = gm2 / gm_pair
    l_km = SATELLITES[secondary].sma_km
    t_s = float(math.sqrt(l_km**3 / gm_pair))
    return CR3BPSystem(mu=mu, primary=primary, secondary=secondary, l_km=l_km, t_s=t_s)
