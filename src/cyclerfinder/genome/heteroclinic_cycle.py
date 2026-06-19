"""#314 Heteroclinic-cycle framework (planar CR3BP).

Discovers and certifies CLOSED HETEROCLINIC CYCLES: chains O_1 -> O_2 -> ... -> O_1
of transversal invariant-manifold connections among equal-energy unstable orbits.
This is a new closure definition -- periodic-up-to-rotation (recurrence to the same
orbit, phase along it free) -- distinct from the strict state(T)=state(0) periodicity
every other genome assumes.

Validated against Wilczak & Zgliczynski's computer-assisted proof of the closed
L1<->L2 Lyapunov cycle in the Sun-Jupiter-Oterma PCR3BP (arXiv:math/0201278). See
docs/superpowers/specs/2026-06-19-314-heteroclinic-cycle-framework-design.md.

Reuses core/cr3bp (propagate + STM + Jacobi) and search/cr3bp_periodic (Lyapunov
orbit corrector); replicates the focused Floquet manifold-seeding pattern from
search/resonance_network.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import solve_ivp

import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.search.cr3bp_periodic import correct_symmetric_fixed_jacobi

_PLANAR_IDX = [0, 1, 3, 4]  # (x, y, xdot, ydot) block of the 6x6 STM


@dataclass(frozen=True)
class LyapunovNode:
    """A libration-point Lyapunov orbit serving as a cycle node.

    ``state0`` is the full 6-vector IC ``(x0, 0, 0, 0, ydot0, 0)`` on the section
    {y=0}; ``unstable_eigvec`` / ``stable_eigvec`` are the planar 4-vectors
    (x, y, xdot, ydot) of the Floquet saddle pair, used to seed the manifolds.
    """

    label: str
    state0: NDArray[np.float64]
    period: float
    jacobi: float
    unstable_eigvec: NDArray[np.float64]
    stable_eigvec: NDArray[np.float64]
    converged: bool

    @classmethod
    def from_libration(
        cls,
        system: cr3bp.CR3BPSystem,
        *,
        x0_guess: float,
        jacobi: float,
        period_guess: float,
        label: str,
        ydot0_sign: float = 1.0,
        tol: float = 1e-10,
    ) -> LyapunovNode:
        """Correct a Lyapunov orbit at fixed Jacobi and extract its Floquet pair."""
        orbit = correct_symmetric_fixed_jacobi(
            system, x0_guess, jacobi, period_guess, ydot0_sign=ydot0_sign, tol=tol
        )
        state0 = np.array([orbit.x0, 0.0, 0.0, 0.0, orbit.ydot0, 0.0], dtype=np.float64)
        jac = cr3bp.jacobi_constant(state0, system.mu)
        _lu, v_u, _ls, v_s = _planar_floquet_pair(system, state0, orbit.period)
        return cls(
            label=label,
            state0=state0,
            period=orbit.period,
            jacobi=jac,
            unstable_eigvec=v_u,
            stable_eigvec=v_s,
            converged=orbit.converged,
        )


def _real_unit(v: NDArray[np.complex128]) -> NDArray[np.float64]:
    """Real-cast + unit-normalise a (possibly complex) eigenvector (4-vector).

    Raises ``ValueError`` if the result is numerically zero even after the
    real+imag fallback: that means the eigenvector has no usable real direction
    (the orbit is not a genuine saddle), and returning a zero vector would
    silently produce a zero manifold perturbation downstream.
    """
    vr = np.real(v)
    n = float(np.linalg.norm(vr))
    if n < 1e-14:
        vr = np.real(v) + np.imag(v)
        n = float(np.linalg.norm(vr))
    if n < 1e-14:
        raise ValueError("eigenvector is numerically zero — orbit is not a saddle")
    return (vr / n).astype(np.float64)


def _planar_floquet_pair(
    system: cr3bp.CR3BPSystem,
    state0: NDArray[np.float64],
    period: float,
    *,
    rtol: float = 1e-12,
    atol: float = 1e-12,
) -> tuple[float, NDArray[np.float64], float, NDArray[np.float64]]:
    """Return ``(|lam_u|, v_u, |lam_s|, v_s)`` for the planar monodromy saddle pair.

    Integrates the 6x6 STM over one period, slices the planar (x, y, xdot, ydot)
    block, and returns the largest- and smallest-magnitude eigenvalues with their
    real-normalised eigenvectors (the unstable and stable Floquet directions).
    Mirrors ``search/resonance_network._planar_floquet`` but returns BOTH ends of
    the reciprocal pair (the connection needs unstable-of-A and stable-of-B).

    NB: eigenvalues are returned as MAGNITUDES (sign stripped). This is correct
    for a Lyapunov saddle (positive real reciprocal pair, the W-Z target). A
    caller working with reflection/period-doubling orbits whose saddle eigenvalue
    is negative real needs the signed value and should use
    ``search/resonance_network._planar_floquet`` directly instead.
    """
    arc = cr3bp.propagate(system, state0, period, with_stm=True, rtol=rtol, atol=atol)
    assert arc.stm is not None
    phi4 = arc.stm[np.ix_(_PLANAR_IDX, _PLANAR_IDX)]
    eigvals, eigvecs = np.linalg.eig(phi4)
    mags = np.abs(eigvals)
    i_u = int(np.argmax(mags))
    i_s = int(np.argmin(mags))
    return (
        float(mags[i_u]),
        _real_unit(eigvecs[:, i_u]),
        float(mags[i_s]),
        _real_unit(eigvecs[:, i_s]),
    )


def _seed_on_manifold(
    system: cr3bp.CR3BPSystem,
    node: LyapunovNode,
    *,
    tau: float,
    direction: str,
    branch: int,
    epsilon: float,
    rtol: float = 1e-12,
    atol: float = 1e-12,
) -> NDArray[np.float64]:
    """Manifold seed at phase ``tau`` along ``node``.

    Propagates the orbit to phase ``tau`` (state + STM), transports the chosen
    Floquet eigenvector with ``Phi(tau)``, renormalises, and ε-perturbs:
    ``state(tau) + branch * epsilon * v_hat(tau)``. ``direction`` selects the
    unstable (forward) or stable (backward) eigenvector.
    """
    if direction not in {"stable", "unstable"}:
        raise ValueError(f"direction must be 'stable' or 'unstable'; got {direction!r}")
    if branch not in (+1, -1):
        raise ValueError(f"branch must be +1 or -1; got {branch!r}")
    tau = float(tau) % float(node.period)
    if tau <= 0.0:
        state_tau = np.asarray(node.state0, float)
        phi = np.eye(6)
    else:
        arc = cr3bp.propagate(system, node.state0, tau, with_stm=True, rtol=rtol, atol=atol)
        assert arc.stm is not None
        state_tau = arc.state_f
        phi = arc.stm
    v4 = node.unstable_eigvec if direction == "unstable" else node.stable_eigvec
    v6 = np.array([v4[0], v4[1], 0.0, v4[2], v4[3], 0.0], dtype=np.float64)
    v_tau = phi @ v6
    n = float(np.linalg.norm(v_tau))
    if n > 0.0:
        v_tau = v_tau / n
    return (state_tau + float(branch) * float(epsilon) * v_tau).astype(np.float64)


def _section_crossing(
    system: cr3bp.CR3BPSystem,
    seed: NDArray[np.float64],
    *,
    direction: str,
    k: int,
    max_time: float,
    section_y: float = 0.0,
    ydot_sign: int | None = None,
    rtol: float = 1e-12,
    atol: float = 1e-12,
    method: str = "DOP853",
) -> NDArray[np.float64] | None:
    """Integrate ``seed`` to the k-th qualifying crossing of {y=section_y}.

    Forward in time for ``direction="unstable"``, backward for ``"stable"``.
    ``ydot_sign`` (if given) restricts to the Theta+/Theta- half of the section
    (sign of ydot at the crossing). Returns the section point ``(x, xdot)`` at the
    k-th crossing (1-based), or ``None`` if fewer than ``k`` qualifying crossings
    occur within ``max_time`` (bounded — never hangs, never fabricates a crossing).
    """
    if direction not in {"stable", "unstable"}:
        raise ValueError(f"direction must be 'stable' or 'unstable'; got {direction!r}")
    horizon = abs(float(max_time))
    t_span = (0.0, horizon) if direction == "unstable" else (0.0, -horizon)

    def _y_event(t: float, y: NDArray[np.float64], _mu: float) -> float:
        return float(y[1] - section_y)

    _y_event.terminal = False  # type: ignore[attr-defined]
    _y_event.direction = 0.0  # type: ignore[attr-defined]

    sol = solve_ivp(
        cr3bp.cr3bp_eom,
        t_span,
        np.asarray(seed, float),
        args=(system.mu,),  # type: ignore[call-overload]
        method=method,
        rtol=rtol,
        atol=atol,
        events=_y_event,
        max_step=horizon / 500.0,  # fine enough that no {y=0} crossing is stepped over
    )
    t_events = sol.t_events[0] if sol.t_events is not None else np.array([])
    y_events = sol.y_events[0] if sol.y_events is not None else []
    t_floor = 1e-6 * horizon
    count = 0
    for t_ev, y_ev in zip(t_events, y_events, strict=False):
        if abs(float(t_ev)) <= t_floor:
            continue  # skip the t~0 root at the seed
        if ydot_sign is not None and int(np.sign(float(y_ev[4]))) != ydot_sign:
            continue
        count += 1
        if count == k:
            return np.array([float(y_ev[0]), float(y_ev[3])], dtype=np.float64)
    return None


@dataclass(frozen=True)
class HeteroclinicConnection:
    """A certified (or attempted) Wu(A) -> Ws(B) connection on {y=0}."""

    orbit_from: str
    orbit_to: str
    jacobi: float
    tau_u: float
    tau_s: float
    k_u: int
    k_s: int
    crossing_xv: NDArray[np.float64]  # (x, xdot) at the matched crossing
    residual: float
    converged: bool
    n_iter: int
    notes: str = ""


def _connection_residual(
    system: cr3bp.CR3BPSystem,
    a: LyapunovNode,
    b: LyapunovNode,
    *,
    tau_u: float,
    tau_s: float,
    k_u: int,
    k_s: int,
    epsilon: float,
    branch_u: int,
    branch_s: int,
    max_time: float,
    ydot_sign_u: int | None,
    ydot_sign_s: int | None,
    method: str = "DOP853",
) -> tuple[NDArray[np.float64] | None, NDArray[np.float64] | None]:
    """Return ``(residual2, crossing_xv)`` or ``(None, None)`` if a leg misses the section."""
    seed_u = _seed_on_manifold(
        system, a, tau=tau_u, direction="unstable", branch=branch_u, epsilon=epsilon
    )
    p_u = _section_crossing(
        system,
        seed_u,
        direction="unstable",
        k=k_u,
        max_time=max_time,
        ydot_sign=ydot_sign_u,
        method=method,
    )
    seed_s = _seed_on_manifold(
        system, b, tau=tau_s, direction="stable", branch=branch_s, epsilon=epsilon
    )
    p_s = _section_crossing(
        system,
        seed_s,
        direction="stable",
        k=k_s,
        max_time=max_time,
        ydot_sign=ydot_sign_s,
        method=method,
    )
    if p_u is None or p_s is None:
        return None, None
    return (p_u - p_s).astype(np.float64), p_u


def _scan_starts(
    resid: Callable[[float, float], tuple[NDArray[np.float64] | None, NDArray[np.float64] | None]],
    period_u: float,
    period_s: float,
    *,
    n: int = 20,
    fixed_u: float | None = None,
    fixed_s: float | None = None,
) -> tuple[float, float]:
    """Coarse grid over the FREE phase(s) -> ``(tau_u, tau_s)`` of least residual norm.

    The Newton target is codimension-1; a blind centre start may sit outside the
    basin. This grids the free phase(s) and returns the cell whose section gap is
    smallest so Newton starts near a genuine intersection. Falls back to the centre
    (or the pinned value) if every cell misses.

    If ``fixed_u``/``fixed_s`` is supplied, that phase is held fixed and only the
    free phase is gridded -- so a caller who pins exactly one phase still gets a
    sensible start for the other (rather than a joint scan whose pinned component is
    then overwritten). Cost is O(n^2) residual evaluations with both phases free,
    O(n) with one pinned; each evaluation is two manifold propagations, so raising
    ``n`` increases cost quadratically.
    """
    us = [float(fixed_u)] if fixed_u is not None else [period_u * (i + 0.5) / n for i in range(n)]
    ss = [float(fixed_s)] if fixed_s is not None else [period_s * (j + 0.5) / n for j in range(n)]
    best = float("inf")
    best_tu = float(fixed_u) if fixed_u is not None else 0.5 * period_u
    best_ts = float(fixed_s) if fixed_s is not None else 0.5 * period_s
    for tu in us:
        for ts in ss:
            res, _ = resid(tu, ts)
            if res is None:
                continue
            rn = float(np.linalg.norm(res))
            if rn < best:
                best, best_tu, best_ts = rn, tu, ts
    return best_tu, best_ts


def correct_connection(
    system: cr3bp.CR3BPSystem,
    orbit_from: LyapunovNode,
    orbit_to: LyapunovNode,
    *,
    k_u: int = 3,
    k_s: int = 4,
    epsilon: float = 1e-6,
    branch_u: int = -1,
    branch_s: int = +1,
    tau_u0: float | None = None,
    tau_s0: float | None = None,
    ydot_sign_u: int | None = None,
    ydot_sign_s: int | None = None,
    tol: float = 1e-7,
    max_iter: int = 40,
    fd_step: float = 1e-6,
    max_time_factor: float = 8.0,
    jacobi_tol: float = 1e-6,
    scan_n: int = 20,
) -> HeteroclinicConnection:
    """Certify Wu(orbit_from) ∩ Ws(orbit_to) on {y=0} by 2-D Newton on (tau_u, tau_s).

    Residual = section-plane gap (Δx, Δxdot). Free vars = manifold phases.
    Jacobian finite-differenced (2x2); Newton step damped by backtracking. Raises
    ``ValueError`` on an energy mismatch. A leg that never reaches the section ->
    ``converged=False`` with a diagnostic note (never a fabricated closure).

    When ``tau_u0``/``tau_s0`` are not supplied, a coarse ``scan_n``-by-``scan_n`` grid
    over the two phases seeds Newton at the cell of least section gap (the codim-1
    intersection rarely lies at a blind centre start).

    Default branch/crossing selection (W-Z Oterma L1->L2)
    -----------------------------------------------------
    The defaults ``branch_u=-1, branch_s=+1, k_u=3, k_s=4`` pick the *neck-facing*
    branch of each manifold: from L1 (interior) and L2 (exterior) only those
    branches reach the L1-L2 neck where the connection lives. With the opposite
    (orbit-hugging) branches the first few crossings cluster on the source orbit
    and the two section curves stay in disjoint x-ranges (no intersection). The
    crossing indices ``(3, 4)`` are the lowest-index pair whose Wu(L1) and Ws(L2)
    section curves cross transversally in (x, xdot); the meeting point lands at
    x~=0.9588 in the neck (near the W-Z Part I crossing 0.95792). Other valid
    transversal connections exist at higher (k_u, k_s) -- e.g. (4, 3) at x~=1.040,
    (5, 4), (5, 5) -- reflecting the rich symbolic dynamics of the W-Z system;
    pass them explicitly to certify a specific one. The exact W-Z crossing-by-
    crossing match is Task 8.
    """
    if max_iter < 1:
        raise ValueError(f"max_iter must be >= 1; got {max_iter}")
    if abs(orbit_from.jacobi - orbit_to.jacobi) > jacobi_tol:
        raise ValueError(
            f"connection requires equal Jacobi (energy): "
            f"{orbit_from.label} C={orbit_from.jacobi:.6f} vs "
            f"{orbit_to.label} C={orbit_to.jacobi:.6f}"
        )
    max_time = max_time_factor * max(orbit_from.period, orbit_to.period)

    def _resid(
        tu: float, ts: float
    ) -> tuple[NDArray[np.float64] | None, NDArray[np.float64] | None]:
        return _connection_residual(
            system,
            orbit_from,
            orbit_to,
            tau_u=tu,
            tau_s=ts,
            k_u=k_u,
            k_s=k_s,
            epsilon=epsilon,
            branch_u=branch_u,
            branch_s=branch_s,
            max_time=max_time,
            ydot_sign_u=ydot_sign_u,
            ydot_sign_s=ydot_sign_s,
        )

    if tau_u0 is not None and tau_s0 is not None:
        tau_u = float(tau_u0)
        tau_s = float(tau_s0)
    else:
        # Grid only the free phase(s); a pinned phase is held fixed during the scan.
        tau_u, tau_s = _scan_starts(
            _resid,
            orbit_from.period,
            orbit_to.period,
            n=scan_n,
            fixed_u=tau_u0,
            fixed_s=tau_s0,
        )

    res, xv = _resid(tau_u, tau_s)
    n_iter = 0
    for n_iter in range(1, max_iter + 1):
        if res is None:
            return HeteroclinicConnection(
                orbit_from.label,
                orbit_to.label,
                orbit_from.jacobi,
                tau_u,
                tau_s,
                k_u,
                k_s,
                np.zeros(2),
                float("inf"),
                False,
                n_iter,
                notes="manifold leg did not reach the section",
            )
        rn = float(np.linalg.norm(res))
        if rn < tol:
            break
        jac = np.zeros((2, 2), dtype=np.float64)
        ok = True
        for j, (du, ds) in enumerate([(fd_step, 0.0), (0.0, fd_step)]):
            rp, _ = _resid(tau_u + du, tau_s + ds)
            if rp is None:
                ok = False
                break
            jac[:, j] = (rp - res) / fd_step
        if not ok:
            return HeteroclinicConnection(
                orbit_from.label,
                orbit_to.label,
                orbit_from.jacobi,
                tau_u,
                tau_s,
                k_u,
                k_s,
                xv if xv is not None else np.zeros(2),
                rn,
                False,
                n_iter,
                notes="FD-Jacobian probe left the manifold's section branch",
            )
        try:
            step = np.linalg.solve(jac, -res)
        except np.linalg.LinAlgError:
            step, *_ = np.linalg.lstsq(jac, -res, rcond=None)
        alpha = 1.0
        improved = False
        for _ in range(20):
            tu_t = (tau_u + alpha * float(step[0])) % orbit_from.period
            ts_t = (tau_s + alpha * float(step[1])) % orbit_to.period
            res_t, xv_t = _resid(tu_t, ts_t)
            if res_t is not None and float(np.linalg.norm(res_t)) < rn:
                tau_u, tau_s, res, xv = tu_t, ts_t, res_t, xv_t
                improved = True
                break
            alpha *= 0.5
        if not improved:
            break
    final_rn = float(np.linalg.norm(res)) if res is not None else float("inf")
    converged = res is not None and final_rn < tol
    return HeteroclinicConnection(
        orbit_from=orbit_from.label,
        orbit_to=orbit_to.label,
        jacobi=orbit_from.jacobi,
        tau_u=tau_u,
        tau_s=tau_s,
        k_u=k_u,
        k_s=k_s,
        crossing_xv=xv if xv is not None else np.zeros(2),
        residual=final_rn,
        converged=converged,
        n_iter=n_iter,
        notes="" if converged else "did not reach tol",
    )
