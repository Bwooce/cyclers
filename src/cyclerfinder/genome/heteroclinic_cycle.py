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
