"""CR3BP multiple-shooting periodic-orbit corrector (task #687).

Extends this project's single-shooting convention
(:func:`cyclerfinder.search.cr3bp_periodic.correct_periodic` -- STM-based Newton
with a min-norm least-squares step) to ``N`` patch points, reusing
:func:`cyclerfinder.core.cr3bp.propagate` (``with_stm=True``) verbatim as the
segment propagator. NOT a new correction scheme: same variational-STM Newton, same
full-state periodicity constraint, just segmented.

Why multiple shooting is needed here
-------------------------------------
Single shooting integrates the STM across the WHOLE orbit in one arc. For a long,
strongly-hyperbolic trajectory -- e.g. a Saturn-Titan capture -> escape ->
re-capture transit whose closest approach to Saturn is ~40 Titan-Hill-radii out
(task #683's seed #217) -- the full-arc monodromy norm is astronomical
(empirically ``~1e18`` for that itinerary: a product of four ``~1e4-1e5`` per-arc
STMs). A Newton step built from a ``1e18``-norm Jacobian is numerically useless,
and single shooting on that seed diverges to a spurious, off-energy far-field
orbit (Jacobi wanders from ``3.015`` to ``-11.6``). Breaking the trajectory at the
periapse nodes keeps each segment's STM ``~1e4-1e5`` and the Newton system
conditioned -- the textbook reason multiple shooting exists.

Formulation
-----------
Free variables ``z`` (``7N`` of them): the ``N`` patch-point 6-states
``X_0 .. X_{N-1}`` (``6N``) plus the ``N`` segment durations ``tau_0 .. tau_{N-1}``.
Constraints ``F`` (``6N`` of them): full-state continuity around a periodic loop,

    F_i = phi(X_i, tau_i) - X_{(i+1) mod N} = 0,   i = 0 .. N-1

where ``phi`` is the CR3BP flow. The ``i``-th 6-row block of the Jacobian is

    dF_i/dX_i        = Phi_i           (segment STM, 6x6)
    dF_i/dX_{i+1}    = -I              (6x6, wraps mod N)
    dF_i/dtau_i      = f(phi(X_i,tau_i))   (CR3BP RHS at the segment end, 6x1)

The system is underdetermined (``7N`` unknowns, ``6N`` constraints -- the extra
``N`` freedoms are the along-flow time-shift/redistribution of the nodes), so each
step is the min-norm Levenberg-Marquardt solution
``dz = J^T (J J^T + lm I)^{-1} (-F)`` (the ``correct_periodic`` min-norm idea with
a small ``lm`` ridge for the near-rank-deficient hyperbolic case), globalised by a
backtracking line search so a bad full step cannot blow the iterate up. Converged
iff ``||F|| < tol`` with all segment times positive. Jacobi is conserved to
integrator precision by the dynamics; it is reported, not constrained (matching
``correct_periodic``).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

import cyclerfinder.core.cr3bp as cr3bp


@dataclass(frozen=True)
class MultipleShootingOrbit:
    """Result of a multiple-shooting periodicity correction.

    ``nodes``/``seg_times`` are the corrected patch points and the segment
    durations between them (``nodes[i]`` flows to ``nodes[(i+1) % N]`` over
    ``seg_times[i]``). ``period`` is ``sum(seg_times)``; ``closure_residual`` is
    the L2 norm of the stacked full-state continuity defect ``F`` (position AND
    velocity at every node) -- the honest analogue of
    :attr:`cyclerfinder.search.cr3bp_periodic.PeriodicOrbit.closure_residual`.
    """

    nodes: list[NDArray[np.float64]]
    seg_times: list[float]
    period: float
    jacobi: float
    converged: bool
    closure_residual: float
    n_iter: int


def _residual_and_jacobian(
    system: cr3bp.CR3BPSystem,
    nodes: list[NDArray[np.float64]],
    seg_times: list[float],
    *,
    rtol: float,
    atol: float,
) -> tuple[NDArray[np.float64], NDArray[np.float64], list[NDArray[np.float64]]]:
    """Stacked continuity residual ``F`` (6N), its Jacobian ``J`` (6N x 7N), and
    the per-segment STMs. Reuses :func:`cr3bp.propagate` (``with_stm=True``).

    Propagation failures (integrator ``sol.success=False``, or the #652
    close-secondary-encounter / wall-clock events, all ``RuntimeError``
    subclasses) propagate out to the caller, which treats them as a
    non-convergence for this iterate rather than a crash.
    """
    n = len(nodes)
    mu = system.mu
    f_vec = np.zeros(6 * n, dtype=np.float64)
    stms: list[NDArray[np.float64]] = []
    f_ends: list[NDArray[np.float64]] = []
    for i in range(n):
        arc = cr3bp.propagate(system, nodes[i], seg_times[i], with_stm=True, rtol=rtol, atol=atol)
        assert arc.stm is not None
        stms.append(arc.stm)
        f_ends.append(cr3bp.cr3bp_eom(0.0, arc.state_f, mu))
        f_vec[6 * i : 6 * i + 6] = arc.state_f - nodes[(i + 1) % n]
    jac = np.zeros((6 * n, 7 * n), dtype=np.float64)
    eye6 = np.eye(6)
    for i in range(n):
        j_next = (i + 1) % n
        jac[6 * i : 6 * i + 6, 6 * i : 6 * i + 6] = stms[i]
        jac[6 * i : 6 * i + 6, 6 * j_next : 6 * j_next + 6] -= eye6
        jac[6 * i : 6 * i + 6, 6 * n + i] = f_ends[i]
    return f_vec, jac, stms


def correct_multiple_shooting(
    system: cr3bp.CR3BPSystem,
    nodes_guess: list[NDArray[np.float64]],
    seg_times_guess: list[float],
    *,
    tol: float = 1e-10,
    max_iter: int = 80,
    lm: float = 1e-8,
    rtol: float = 1e-12,
    atol: float = 1e-12,
    min_seg_time: float = 1e-4,
    max_period_factor: float = 10.0,
) -> MultipleShootingOrbit:
    """Multiple-shooting periodicity correction (see module docstring).

    Parameters
    ----------
    nodes_guess, seg_times_guess :
        Initial patch points (each a planar or spatial CR3BP 6-state) and the
        segment durations between successive nodes (``len(seg_times) ==
        len(nodes) == N``; the last segment closes ``nodes[N-1] -> nodes[0]``).
    tol :
        Convergence threshold on ``||F||`` (full-state continuity, L2).
    lm :
        Levenberg-Marquardt ridge for the min-norm step (stabilises the
        near-rank-deficient normal equations of a strongly-hyperbolic orbit).
    max_period_factor :
        Abort if the total period exceeds this multiple of the initial guess
        (divergence guard, mirroring ``correct_periodic``'s ``max_period_factor``).

    Returns
    -------
    MultipleShootingOrbit
        ``converged`` is ``True`` only if ``||F|| < tol`` with every segment
        time positive. A propagation failure at any iterate (e.g. a #652
        close-secondary encounter developing mid-arc) ends the iteration and
        returns the best iterate so far with ``converged=False``.
    """
    if len(nodes_guess) != len(seg_times_guess):
        raise ValueError(
            f"nodes ({len(nodes_guess)}) and seg_times ({len(seg_times_guess)}) length mismatch"
        )
    if len(nodes_guess) == 0:
        raise ValueError("need at least one node")
    n = len(nodes_guess)
    nodes = [np.asarray(x, dtype=np.float64).copy() for x in nodes_guess]
    seg = [float(t) for t in seg_times_guess]
    period_guess = sum(seg)
    period_hi = max_period_factor * abs(period_guess)

    try:
        f_vec, jac, _ = _residual_and_jacobian(system, nodes, seg, rtol=rtol, atol=atol)
    except RuntimeError:
        return _package(system, nodes, seg, float("inf"), 0)
    res = float(np.linalg.norm(f_vec))

    n_iter = 0
    for n_iter in range(1, max_iter + 1):  # noqa: B007 -- returned as the iteration count
        if res < tol:
            break
        # Min-norm Levenberg-Marquardt step: dz = J^T (J J^T + lm I)^{-1} (-F).
        gram = jac @ jac.T + lm * np.eye(6 * n)
        try:
            dz = jac.T @ np.linalg.solve(gram, -f_vec)
        except np.linalg.LinAlgError:
            dz, *_ = np.linalg.lstsq(jac, -f_vec, rcond=None)
        # Backtracking line search: shrink the step until ||F|| strictly
        # decreases and every segment time stays positive.
        alpha = 1.0
        accepted = False
        for _ in range(40):
            trial_nodes = [nodes[i] + alpha * dz[6 * i : 6 * i + 6] for i in range(n)]
            trial_seg = [seg[i] + alpha * float(dz[6 * n + i]) for i in range(n)]
            if any(t <= min_seg_time for t in trial_seg) or sum(trial_seg) > period_hi:
                alpha *= 0.5
                continue
            try:
                f_try, jac_try, _ = _residual_and_jacobian(
                    system, trial_nodes, trial_seg, rtol=rtol, atol=atol
                )
            except RuntimeError:
                alpha *= 0.5
                continue
            res_try = float(np.linalg.norm(f_try))
            if res_try < res:
                nodes, seg, f_vec, jac, res = trial_nodes, trial_seg, f_try, jac_try, res_try
                accepted = True
                break
            alpha *= 0.5
        if not accepted:
            # No downhill step exists from here: the residual surface has no
            # descent toward zero near this iterate (a genuine stall, not a
            # crash) -- the corrector's honest "no nearby periodic orbit" signal.
            break

    converged = res < tol and all(t > 0.0 for t in seg)
    return _package(system, nodes, seg, res, n_iter, converged=converged)


def _package(
    system: cr3bp.CR3BPSystem,
    nodes: list[NDArray[np.float64]],
    seg: list[float],
    res: float,
    n_iter: int,
    *,
    converged: bool = False,
) -> MultipleShootingOrbit:
    return MultipleShootingOrbit(
        nodes=nodes,
        seg_times=seg,
        period=sum(seg),
        jacobi=cr3bp.jacobi_constant(nodes[0], system.mu),
        converged=converged,
        closure_residual=res,
        n_iter=n_iter,
    )


def monodromy(
    system: cr3bp.CR3BPSystem,
    orbit: MultipleShootingOrbit,
    *,
    rtol: float = 1e-12,
    atol: float = 1e-12,
) -> NDArray[np.float64]:
    """Full-period monodromy matrix as the ordered product of the segment STMs
    ``Phi = Phi_{N-1} @ ... @ Phi_1 @ Phi_0`` (each ``Phi_i`` from
    :func:`cr3bp.propagate` over ``seg_times[i]``). Its eigenvalues are the
    orbit's Floquet multipliers."""
    mat = np.eye(6)
    for x, t in zip(orbit.nodes, orbit.seg_times, strict=True):
        arc = cr3bp.propagate(system, x, t, with_stm=True, rtol=rtol, atol=atol)
        assert arc.stm is not None
        mat = arc.stm @ mat
    return mat


def floquet_multipliers(
    system: cr3bp.CR3BPSystem,
    orbit: MultipleShootingOrbit,
    *,
    rtol: float = 1e-12,
    atol: float = 1e-12,
) -> NDArray[np.complex128]:
    """Floquet multipliers (eigenvalues of the :func:`monodromy` matrix). A
    genuine periodic orbit has a trivial multiplier pair near ``1``; the largest
    ``|lambda|`` is the linear instability rate over one period."""
    eigs = np.linalg.eigvals(monodromy(system, orbit, rtol=rtol, atol=atol))
    return np.asarray(eigs, dtype=np.complex128)
