"""ER3BP fixed-period saddle-center branch-switcher (#436 Task 3).

CONTINGENT INFRASTRUCTURE — HONEST CAVEAT
-----------------------------------------
This module builds the ER3BP analogue of the CR3BP saddle-center branch-switcher
:func:`cyclerfinder.genome.asymmetric_branch.branch_at_saddle_center`. It is meant
to be invoked when an e-continuation flags a genuine elliptic<->hyperbolic Floquet
transition (a saddle-center bifurcation): perturb the parent IC along the marginal
eigenvector and re-converge to land on the new branch.

HOWEVER: #432 and #435 found ZERO bifurcations in the ER3BP e-continuations
investigated so far (the Broucke 7P family stays hyperbolic — no stability regime
flip along e). So this switcher has NO current target. It is built correctly and
tested for SAFE behaviour (it must not fabricate a fake branch on a non-bifurcating
parent), but it is contingent infrastructure: it only does useful work once a real
saddle-center is discovered in the e-continuation.

KEY DIFFERENCE FROM THE CR3BP SWITCHER
--------------------------------------
The ER3BP is non-autonomous: the true anomaly ``f`` appears explicitly in the EOMs,
so strict periodicity requires the period in ``f`` to be a multiple of ``2*pi``.
The period is therefore FIXED (locked). Unlike the CR3BP switcher — which frees the
period in a 7-variable corrector (free vars including T) — the ER3BP re-convergence
uses :func:`cyclerfinder.genome.er3bp_periodic.correct_er3bp_periodic` with the SAME
``period_f`` and a 2-variable symmetric corrector (free (x, y') -> residual (y, x')
at the half-period). The period is never a free variable.

PERIOD CONVENTION (caller contract)
-----------------------------------
``parent_period_f`` is the integration span the corrector expects — i.e. the
HALF-period (``pi`` for a full ``2*pi`` orbit) when the half-period symmetric
residual is used. This matches the span passed by
:func:`cyclerfinder.search.er3bp_discovery.continue_and_monitor` and
:func:`cyclerfinder.search.er3bp_direct_seeding.converge_direct_seed`. The monodromy
is computed over the FULL period (``2 * parent_period_f`` here), since the
saddle-center Floquet structure lives on the full-period monodromy.
"""

from __future__ import annotations

import math

import numpy as np
from numpy.typing import NDArray

import cyclerfinder.core.er3bp as er3bp
from cyclerfinder.genome.asymmetric_branch import _select_saddle_center_eigenvector
from cyclerfinder.genome.er3bp_periodic import (
    IDX_X,
    IDX_XDOT,
    IDX_Y,
    IDX_YDOT,
    ER3BPPeriodicOrbit,
    correct_er3bp_periodic,
)
from cyclerfinder.search.er3bp_floquet import er3bp_monodromy


def branch_at_saddle_center_er3bp(
    system: er3bp.ER3BPSystem,
    parent_state0: NDArray[np.float64],
    parent_period_f: float,
    *,
    epsilon: float = 1e-3,
    tol: float = 1e-10,
    max_iter: int = 80,
    rtol: float = 1e-12,
    atol: float = 1e-12,
    independent_tol: float = 1e-5,
) -> tuple[ER3BPPeriodicOrbit | None, dict[str, object]]:
    """Branch off an ER3BP symmetric family at a saddle-center bifurcation.

    Pipeline (mirrors :func:`asymmetric_branch.branch_at_saddle_center`, but with
    the period FIXED — see module docstring):

      1. Compute the parent's FULL-period monodromy (over ``2 * parent_period_f``).
      2. Pick the marginal eigenvector ``v`` (the secondary non-trivial pair's
         closest-to-+1 right-eigenvector) via the shared
         :func:`asymmetric_branch._select_saddle_center_eigenvector`,
         Gram-Schmidt-orthogonalised against the parent's flow tangent
         ``f(0, state0)`` so the perturbation does not slide along the parent.
      3. Perturb the parent IC by ``+/- eps * v`` over a small eps/sign ladder.
      4. Re-converge each perturbed IC via
         :func:`er3bp_periodic.correct_er3bp_periodic` with the SAME ``period_f``
         (the half-period symmetric corrector — the period is never freed).
      5. Return the first converged orbit (corrector residual < ``tol`` AND
         independent closure < ``independent_tol``), or ``None`` with a reason.

    Parameters
    ----------
    system :
        ER3BP system (``mu`` and ``e`` are read).
    parent_state0 :
        6-vector parent IC at ``f = 0``. Should be a symmetric family member at a
        flagged saddle-center; the marginal eigenvector is computed at THIS state.
    parent_period_f :
        The corrector's integration span — the HALF-period (``pi`` for a full
        ``2*pi`` orbit). See module docstring "PERIOD CONVENTION".
    epsilon :
        Perturbation amplitude along the eigenvector (nondim state units).
    tol, max_iter, rtol, atol, independent_tol :
        Forwarded to :func:`correct_er3bp_periodic`. A converged orbit must
        satisfy corrector residual < ``tol`` AND independent closure
        < ``independent_tol``.

    Returns
    -------
    tuple[ER3BPPeriodicOrbit | None, dict] :
        ``(orbit, {"epsilon", "sign", "eigenvalue"})`` on success, else
        ``(None, {"reason": ...})`` with reason ``"no marginal eigenvector"`` or
        ``"no perturbation converged"``.

    Notes
    -----
    Discipline (orbit-closure): ``None`` is returned on no-converge — never a
    "near-converged" orbit relabelled as passing. On a non-bifurcating parent the
    switcher must not invent a false discovery: the perturbed IC either fails to
    converge (returns ``None``) or re-converges back onto the parent (an orbit
    within tolerance of ``parent_state0``), not onto a fabricated distinct family.
    """
    parent_state0 = np.asarray(parent_state0, dtype=np.float64).copy()
    if parent_state0.shape != (6,):
        raise ValueError(f"parent_state0 must have shape (6,); got {parent_state0.shape}")
    if parent_period_f <= 0.0 or not math.isfinite(parent_period_f):
        raise ValueError(f"parent_period_f must be > 0 finite; got {parent_period_f}")
    if epsilon <= 0.0 or not math.isfinite(epsilon):
        raise ValueError(f"epsilon must be > 0 finite; got {epsilon}")

    # The saddle-center Floquet structure lives on the FULL-period monodromy. The
    # caller passes the half-period integration span, so the full period is 2x.
    full_period_f = 2.0 * parent_period_f
    mono = er3bp_monodromy(parent_state0, full_period_f, system)

    # Parent flow tangent at the IC: f(0, state0). This is the time-translation
    # (trivial) monodromy eigenvector; projecting it out keeps the perturbation
    # off the parent's own orbit. EOM signature: er3bp_eom(f, state, mu, e).
    parent_tangent = er3bp.er3bp_eom(0.0, parent_state0, system.mu, system.e)

    pick = _select_saddle_center_eigenvector(mono, parent_tangent=parent_tangent)
    if pick is None:
        return None, {"reason": "no marginal eigenvector"}
    v, lam = pick

    # Epsilon ladder mirrors branch_at_saddle_center: caller's epsilon first (both
    # signs), then a small set of nearby amplitudes (dedup). The saddle-center
    # basin is marginal; the ladder only ADDS fallback attempts after the
    # requested amplitude is exhausted.
    epsilon_ladder = [epsilon]
    for fallback in (7e-4, 3e-4, 1e-3, 1e-4, 1.5e-3, 5e-4):
        if all(abs(fallback - e) > 1e-15 for e in epsilon_ladder):
            epsilon_ladder.append(fallback)

    for eps in epsilon_ladder:
        for sign in (+1, -1):
            perturbed = parent_state0 + sign * eps * v
            try:
                orbit = correct_er3bp_periodic(
                    system,
                    perturbed,
                    parent_period_f,
                    free_vars=(IDX_X, IDX_YDOT),
                    residual_indices=(IDX_Y, IDX_XDOT),
                    is_half_period_residual=True,
                    tol=tol,
                    max_iter=max_iter,
                    rtol=rtol,
                    atol=atol,
                    independent_tol=independent_tol,
                )
            except Exception:
                # ConvergenceError / SingularJacobian / propagation failure — try
                # the next sign/amplitude rather than fabricating a result.
                continue
            if orbit.corrector_residual < tol and orbit.independent_residual < independent_tol:
                return orbit, {"epsilon": eps, "sign": sign, "eigenvalue": complex(lam)}

    return None, {"reason": "no perturbation converged"}
