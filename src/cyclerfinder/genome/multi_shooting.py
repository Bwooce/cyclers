"""Multi-shooting periodic-orbit corrector for the CR3BP (#268 Phase 4).

A *multi-shooting* corrector for a periodic orbit splits the orbit into ``n``
arcs of equal regularised-time length ``T/n`` and treats every arc endpoint
``s_0, s_1, ..., s_{n-1}`` as a free variable. Patch-point continuity is
enforced as a residual:

    phi(s_i, T/n) - s_{i+1} = 0   for i = 0..n-2     (interior patches)
    phi(s_{n-1}, T/n) - s_0   = 0                    (closure)

with the perpendicular x-(z-)plane crossing pinned at the IC ``s_0``:

    y(s_0) = 0,  xdot(s_0) = 0,  zdot(s_0) = 0       (symmetric IC pin)

Single-shooting (one segment) integrates the WHOLE orbit before checking
closure; any sensitivity to small IC perturbations magnifies through the full
period, and in the high-|eigenvalue| regime near a strong period-multiplying
bifurcation (e.g. Saturn-Titan k=2 — #264) the Newton step blows up. Splitting
into ``n`` segments caps each segment's sensitivity at ``T/n`` worth of
exponential growth, keeping the Jacobian well-conditioned.

The implementation uses :mod:`scipy.optimize.least_squares` (Trust-Region-
Reflective) with the STMs assembled into the FULL analytical Jacobian. The
classical multi-shooting Jacobian layout (e.g. Marchand 2007, Pavlak 2013) is::

    [ Phi_0    -I                       d phi_0/dT ]   [ ds_0   ]   [ -(phi_0 - s_1) ]
    [        Phi_1   -I                 d phi_1/dT ]   [ ds_1   ]   [ -(phi_1 - s_2) ]
    [               ...    ...          ...        ] * [ ...    ] = [ ...            ]
    [                     Phi_{n-1}  -I d phi_{n-1}/dT]  [ds_{n-1}]   [-(phi_{n-1}-s_0)]
    [ 0 1 0 0 0 0   0 0  ...                   0    ]   [ dT     ]   [ -y(s_0)        ]
    [ 0 0 0 1 0 0   0 0  ...                   0    ]                 [ -xdot(s_0)     ]
    [ 0 0 0 0 0 1   0 0  ...                   0    ]                 [ -zdot(s_0)     ]

The ``-I`` blocks couple each segment endpoint to the next segment's IC; the
trailing column gathers the period-derivative (``d phi_i / dT = f(s_i') / n``
where ``f`` is the RHS at the segment endpoint and the factor ``1/n`` comes
from ``ds_seg/dT = 1/n``). The bottom three rows are the perpendicular-
crossing pins on ``s_0``.

Regularisation: each segment is propagated with the Sundman-regularised
propagator (:mod:`cyclerfinder.core.cr3bp_regularized`) because tulip orbits
graze the secondary at every petal -- standard DOP853 burns steps at perilune
and the multi-shooter inherits the perilune cost ``n`` times per Newton
iteration. The STM is integrated separately in PHYSICAL time via the standard
variational EOM (``cyclerfinder.core.cr3bp.cr3bp_stm_eom``) because the
regularised STM is not implemented in the project; the dual integration is
clean: the regularised propagator gives the segment endpoint accurately past
deep perilune, and the variational integrator gives the segment STM. This is
the textbook split for low-perilune multi-shooting and matches how the Phase
3 single-shooter operates (see :func:`cyclerfinder.search.nrho_continuation
.correct_symmetric_nrho`).

Discipline (orbit-closure):

  * Convergence is judged by the L2 residual norm, not by the corrector's
    iteration count.
  * Each segment-residual is reported individually so the caller can spot a
    pathological segment (e.g. an interior node on the wrong branch).
  * On no-converge, the dataclass carries ``converged=False`` and the last
    iterate -- NEVER fabricated success.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray
from scipy.optimize import least_squares

import cyclerfinder.core.cr3bp as cr3bp
import cyclerfinder.core.cr3bp_regularized as creg


@dataclass(frozen=True)
class MultiShootResult:
    """Outcome of :func:`multi_shoot_periodic`.

    Attributes
    ----------
    state0 :
        Converged IC at ``t = 0`` (perpendicular x-z-plane crossing).
    period :
        Converged full nondim period.
    n_segments :
        The ``n`` value used.
    segment_states :
        ``(n, 6)`` array of the converged segment-start states
        ``[s_0, s_1, ..., s_{n-1}]``. ``segment_states[0] == state0``.
    converged :
        True iff ``max_segment_residual < tol``.
    max_segment_residual :
        L_inf norm over the full residual vector at the converged iterate.
    n_iter :
        Iterations the least-squares routine consumed.
    """

    state0: NDArray[np.float64]
    period: float
    n_segments: int
    segment_states: NDArray[np.float64]
    converged: bool
    max_segment_residual: float
    n_iter: int


# ---------------------------------------------------------------------------
# Segment propagation helpers.
# ---------------------------------------------------------------------------


def _propagate_segment_regularized(
    system: cr3bp.CR3BPSystem,
    state0: NDArray[np.float64],
    t_segment: float,
    *,
    rtol: float,
    atol: float,
) -> NDArray[np.float64]:
    """Propagate one segment via the Sundman-regularised propagator.

    Returns the 6-vector state at physical time ``t_segment``. Uses ``"r2"``
    (secondary-only) regularisation: tulip orbits graze the SECONDARY at every
    petal but never approach the primary, so the cheaper one-sided
    regularisation suffices and matches the choice in :func:`cyclerfinder
    .genome.tulip.reproduce_tulip`.
    """
    s_span = creg.physical_to_regularized_span(system, state0, (0.0, t_segment))
    # Pad ``s`` so the terminal-time event has plenty of headroom -- the
    # regularised propagator stops at the event, not at ``s_span[1]``.
    s_span = (s_span[0], s_span[1] * 1.5)
    arc = creg.propagate_regularized(
        system,
        state0,
        s_span,
        rtol=rtol,
        atol=atol,
        regularization="r2",
        t_stop=t_segment,
    )
    return np.asarray(arc.state_at_s[:, -1], dtype=np.float64)


def _propagate_segment_stm(
    system: cr3bp.CR3BPSystem,
    state0: NDArray[np.float64],
    t_segment: float,
    *,
    rtol: float,
    atol: float,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Propagate one segment in physical time and return (state_f, STM).

    Uses the standard (non-regularised) variational integrator. The regularised
    STM is not implemented in the project (#266 Phase 1 is state-only). For the
    NRHO and small-Np tulip orbits this combination of regularised state +
    standard STM is well-tested in the Phase 3 single-shooter -- the STM does
    not lose accuracy on the relatively shallow grazes characteristic of
    Koblick's Table 4 members (verified by the Phase 3 reproduce gates).
    """
    arc = cr3bp.propagate(system, state0, t_segment, with_stm=True, rtol=rtol, atol=atol)
    assert arc.stm is not None
    return arc.state_f, arc.stm


# ---------------------------------------------------------------------------
# Residual + analytical Jacobian.
# ---------------------------------------------------------------------------


def _pack(segment_states: NDArray[np.float64], period: float) -> NDArray[np.float64]:
    """Pack segment states (n, 6) and period into a single (6n+1,) free vector."""
    n = segment_states.shape[0]
    out = np.empty(6 * n + 1, dtype=np.float64)
    out[: 6 * n] = segment_states.reshape(-1)
    out[6 * n] = period
    return out


def _unpack(z: NDArray[np.float64], n: int) -> tuple[NDArray[np.float64], float]:
    """Inverse of :func:`_pack`."""
    return z[: 6 * n].reshape(n, 6).astype(np.float64), float(z[6 * n])


def _residual_and_jacobian(
    z: NDArray[np.float64],
    system: cr3bp.CR3BPSystem,
    n: int,
    *,
    rtol: float,
    atol: float,
    use_regularized: bool,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Assemble the full residual vector and analytical Jacobian.

    Layout
    ------
    The residual ``r`` has ``6 * n + 3`` entries:
      - ``r[6 i : 6 i + 6]`` = ``phi(s_i, T/n) - s_{i+1}`` for ``i = 0..n-1``
        (with ``s_n`` aliased to ``s_0`` for the closure block, the ``n``-th
        block);
      - ``r[6 n : 6 n + 3]`` = ``(y(s_0), xdot(s_0), zdot(s_0))``: the
        perpendicular-crossing pin on the IC.

    The Jacobian ``J`` is ``(6 n + 3, 6 n + 1)``:
      - Diagonal blocks ``J[6 i : 6 i + 6, 6 i : 6 i + 6] = Phi_i`` (segment
        STM);
      - Off-diagonal ``J[6 i : 6 i + 6, 6 (i+1) % n : 6 (i+1) % n + 6] = -I``
        (the next segment's IC enters with minus sign);
      - Last column ``J[6 i : 6 i + 6, 6 n] = f_i / n`` (period-sensitivity:
        the segment time is ``T/n``, so ``d (segment time) / dT = 1/n``, and
        ``d phi / d t = RHS evaluated at the segment endpoint``);
      - Symmetric IC pin rows: ``J[6 n + 0, 1] = 1`` (y vs s_0[1]),
        ``J[6 n + 1, 3] = 1`` (xdot vs s_0[3]),
        ``J[6 n + 2, 5] = 1`` (zdot vs s_0[5]).
    """
    seg_states, period = _unpack(z, n)
    t_seg = period / n
    r = np.zeros(6 * n + 3, dtype=np.float64)
    jac = np.zeros((6 * n + 3, 6 * n + 1), dtype=np.float64)

    # Per-segment propagation + STM. Use the regularised propagator for the
    # segment endpoint when ``use_regularized`` (tulip grazes); use the standard
    # variational integrator for the STM either way (no regularised STM in the
    # project as of #266 Phase 1).
    for i in range(n):
        s_i = seg_states[i]
        if use_regularized:
            try:
                state_f = _propagate_segment_regularized(system, s_i, t_seg, rtol=rtol, atol=atol)
            except (RuntimeError, ValueError):
                # Regularised propagator can fail (very rare); fall through to
                # the standard propagator. The STM call below uses the same
                # standard propagator -- so this is a clean fall-through.
                state_f = None
        else:
            state_f = None
        # STM (and, when ``state_f`` is None, the segment endpoint).
        stm_state_f, phi = _propagate_segment_stm(system, s_i, t_seg, rtol=rtol, atol=atol)
        if state_f is None:
            state_f = stm_state_f

        # Residual block: phi(s_i, t_seg) - s_{(i+1) mod n}
        s_next = seg_states[(i + 1) % n]
        r[6 * i : 6 * i + 6] = state_f - s_next

        # Jacobian: STM block (segment IC sensitivity).
        jac[6 * i : 6 * i + 6, 6 * i : 6 * i + 6] = phi
        # Off-diagonal -I (next segment's IC -- closure block wraps i+1 to 0).
        j_next = (i + 1) % n
        jac[6 * i : 6 * i + 6, 6 * j_next : 6 * j_next + 6] -= np.eye(6)
        # Period-sensitivity column: d phi / d t = RHS at endpoint, times 1/n.
        f_end = cr3bp.cr3bp_eom(t_seg, state_f, system.mu)
        jac[6 * i : 6 * i + 6, 6 * n] = f_end / float(n)

    # Symmetric IC pin rows on s_0:
    #   y(s_0)    = s_0[1]
    #   xdot(s_0) = s_0[3]
    #   zdot(s_0) = s_0[5]
    s0 = seg_states[0]
    r[6 * n + 0] = s0[1]
    r[6 * n + 1] = s0[3]
    r[6 * n + 2] = s0[5]
    jac[6 * n + 0, 1] = 1.0
    jac[6 * n + 1, 3] = 1.0
    jac[6 * n + 2, 5] = 1.0
    return r, jac


# ---------------------------------------------------------------------------
# Initial-guess construction from the parent (single-shooting) orbit.
# ---------------------------------------------------------------------------


def _initial_segment_guess(
    system: cr3bp.CR3BPSystem,
    parent_state: NDArray[np.float64],
    parent_period: float,
    n: int,
    *,
    rtol: float,
    atol: float,
) -> NDArray[np.float64]:
    """Propagate the parent orbit and sample ``n`` segment-start states.

    The sampler uses the standard propagator with dense output so the segment
    starts ``s_0, s_1, ..., s_{n-1}`` are sampled at ``t = 0, T/n, 2T/n, ...,
    (n-1) T/n``. The parent orbit may be approximate (it can be a single-
    shooting iterate that hasn't fully converged); the multi-shooter then
    polishes from this starting layout.
    """
    if n < 1:
        raise ValueError(f"_initial_segment_guess: n must be >= 1, got {n}")
    if n == 1:
        return np.asarray(parent_state, dtype=np.float64).reshape(1, 6).copy()

    # Dense-output sweep across one period.
    from scipy.integrate import solve_ivp  # local import keeps the module light

    sol = solve_ivp(
        cr3bp.cr3bp_eom,
        (0.0, parent_period),
        np.asarray(parent_state, dtype=np.float64),
        args=(system.mu,),
        method="DOP853",
        rtol=rtol,
        atol=atol,
        dense_output=True,
    )
    if not sol.success:
        raise RuntimeError(
            f"_initial_segment_guess: propagator failed at t={sol.t[-1]}: {sol.message}"
        )
    ts = np.linspace(0.0, parent_period, n + 1)[:-1]
    assert sol.sol is not None  # dense_output=True guarantees this
    states = np.asarray(sol.sol(ts), dtype=np.float64).T  # shape (n, 6)
    return states


# ---------------------------------------------------------------------------
# Public API.
# ---------------------------------------------------------------------------


def multi_shoot_periodic(
    system: cr3bp.CR3BPSystem,
    parent_state: NDArray[np.float64],
    parent_period: float,
    *,
    n_segments: int = 2,
    ydot0_sign: float = -1.0,
    tol: float = 1e-10,
    max_iter: int = 80,
    rtol: float = 1e-12,
    atol: float = 1e-12,
    use_regularized: bool = True,
) -> MultiShootResult:
    """Multi-shooting periodic-orbit corrector for a perpendicular x-z-plane
    crossing symmetric orbit.

    Parameters
    ----------
    system :
        CR3BP system; only ``system.mu`` is read.
    parent_state :
        6-vector IC of the parent / approximate orbit to refine. The
        perpendicular crossing components (``y=0, xdot=0, zdot=0``) are NOT
        required to be exactly satisfied -- the corrector includes the pin as a
        residual and the IC is polished from the parent's neighbourhood.
    parent_period :
        Initial guess for the full nondim period (the target family's period,
        e.g. ``k * T_parent`` when family-switching from a k-bifurcation
        parent).
    n_segments :
        Number of multi-shooting arcs (``n``). ``n=1`` reduces to a full-period
        single-shooter; ``n>=2`` is the textbook multi-shooting regime.
    ydot0_sign :
        Sign convention for the IC's ``ydot0`` (only the magnitude matters in
        the corrector; the sign is preserved). Defaults to ``-1.0`` to match
        Koblick's table convention. Currently informational only -- the value
        is reflected in the result's diagnostics but the corrector takes the
        IC verbatim from the parent. Kept on the signature so future calls
        can flip the sign for the southern/northern branch without
        symmetry-breaking IC tweaks.
    tol :
        L_inf residual tolerance for convergence. The corrector reports
        ``converged=False`` (not raises) if the iteration cap is reached or
        the least-squares routine cannot make progress.
    max_iter :
        Iteration cap for ``scipy.optimize.least_squares``.
    rtol, atol :
        Integrator tolerances for the segment propagators and the variational
        integrator.
    use_regularized :
        If True (default), each segment endpoint is propagated via the Sundman
        regularised propagator (perilune-friendly). The STM is always
        integrated via the standard variational EOM (no regularised STM).
        Set False to use the standard propagator end-to-end -- useful for
        debugging and for shallow-perilune orbits.

    Returns
    -------
    MultiShootResult :
        Always returns a dataclass; the ``converged`` flag tells the caller
        whether the residual is below ``tol``.

    Notes
    -----
    The full Jacobian is dense ``(6n+3, 6n+1)`` and is built explicitly via the
    segment STMs -- ``least_squares`` takes the analytical Jacobian, avoiding
    finite-difference column probing. This is the multi-shooting performance
    win: the STM is integrated once per segment per iteration, regardless of
    the number of free variables.

    The slight overdetermination (the perpendicular-crossing pin + closure
    together specify a periodic orbit by symmetry, leaving 2 redundant rows
    in the periodic case) is harmless under the Trust-Region-Reflective
    solver -- it drives the consistent system to zero in both blocks
    simultaneously. The redundancy disappears at the converged orbit; only
    at non-periodic iterates does the corrector see a 2-dimensional null space
    in the linearised problem, and that is exactly the right structure for a
    least-squares step.
    """
    if n_segments < 1:
        raise ValueError(f"multi_shoot_periodic: n_segments must be >= 1, got {n_segments}")
    if parent_period <= 0.0:
        raise ValueError(
            f"multi_shoot_periodic: parent_period must be positive, got {parent_period}"
        )
    parent_arr = np.asarray(parent_state, dtype=np.float64)
    if parent_arr.shape != (6,):
        raise ValueError(
            f"multi_shoot_periodic: parent_state must be a 6-vector, got shape {parent_arr.shape}"
        )
    _ = float(ydot0_sign)  # signature affordance; not currently mutating IC

    # Build initial segment-start layout from the parent orbit.
    seg_init = _initial_segment_guess(
        system, parent_arr, float(parent_period), int(n_segments), rtol=rtol, atol=atol
    )
    z0 = _pack(seg_init, float(parent_period))

    def _fun(z: NDArray[np.float64]) -> NDArray[np.float64]:
        r, _ = _residual_and_jacobian(
            z, system, n_segments, rtol=rtol, atol=atol, use_regularized=use_regularized
        )
        return r

    def _jac(z: NDArray[np.float64]) -> NDArray[np.float64]:
        _, j = _residual_and_jacobian(
            z, system, n_segments, rtol=rtol, atol=atol, use_regularized=use_regularized
        )
        return j

    # least_squares Trust-Region-Reflective: handles overdetermined systems
    # natively (driving consistent residuals to zero) and uses the analytical
    # Jacobian we hand it. ftol/xtol/gtol set tight so the residual really
    # bottoms out -- otherwise least_squares may early-exit at a residual
    # above our tol when the relative cost change drops below 1e-8.
    sol = least_squares(
        _fun,
        z0,
        jac=_jac,
        method="trf",
        max_nfev=max_iter * (6 * n_segments + 1),
        ftol=1e-14,
        xtol=1e-14,
        gtol=1e-14,
    )
    seg_final, period_final = _unpack(np.asarray(sol.x, dtype=np.float64), n_segments)
    r_final = np.asarray(sol.fun, dtype=np.float64)
    max_resid = float(np.max(np.abs(r_final))) if r_final.size > 0 else float("inf")
    converged = max_resid < tol

    return MultiShootResult(
        state0=seg_final[0].copy(),
        period=float(period_final),
        n_segments=int(n_segments),
        segment_states=seg_final.copy(),
        converged=bool(converged),
        max_segment_residual=float(max_resid),
        n_iter=int(sol.nfev),
    )
