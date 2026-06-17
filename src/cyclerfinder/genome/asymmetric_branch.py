"""Asymmetric / 3D branch corrector at a saddle-center bifurcation (#347 Phase 1).

Mirrors :mod:`cyclerfinder.genome.family_switch` (which handles period-multiplying
k>=2 bifurcations and ALWAYS re-corrects via the SYMMETRIC perpendicular-x-axis-
crossing corrector :func:`cyclerfinder.search.nrho_continuation.correct_symmetric_nrho`)
but routes the post-perturbation IC through
:func:`cyclerfinder.search.cr3bp_general_periodic_3d.correct_general_periodic_3d`
in **full-asymmetric mode**: free vars (x0, y0, z0, xdot0, ydot0, zdot0, T),
residual (x - x0, y - y0, z - z0, xdot - xdot0, ydot - ydot0, zdot - zdot0) at T.

This is the corrector for the k=1 saddle-center / pitchfork bifurcation (the
codimension-1 lambda=+1 bifurcation in Hamiltonian-with-symmetry systems —
RTR2026 p.5, classical theory in Golubitsky-Stewart-Marsden). The pre-existing
:func:`cyclerfinder.genome.family_switch.switch_family` always re-corrects via
the symmetric corrector, which is the right move for period-doubling
(k=2 → 2T) but the WRONG move for a saddle-center / pitchfork: the new family
may break the time-reversal symmetry (the residual ``(y, xdot, zdot)|T/2 = 0``
no longer holds at the branched orbit). Single-shooting on the full 6D state
closure is the textbook fix.

The 3D variant (z != 0) is consistent with the (3,2) Earth-Moon C32 saddle-
center: the marginal eigenvectors at i=124 (C=3.14180, post-bifurcation) of
the P1.2 walk artifact point in the (z, zdot) direction — components 2 and 5
of the right-eigenvector — with the in-plane components (x, y, xdot, ydot)
near zero. The branched family is a 3D family with z amplitude.

Discipline:

  * NO catalogue writeback in Phase 1 (the orbit emitted here is "discovery
    candidate", not a sourced row).
  * The marginal eigenvector is freshly computed at the bifurcation point's
    PARENT family member (the last member where the secondary pair is still
    complex on the unit circle). Branching off the POST-bifurcation member
    would seed the corrector on the already-bifurcated geometry; the
    pre-bifurcation member is the "parent" in the branch-switching sense.
  * Independent topology cross-check via
    :func:`cyclerfinder.search.binary_star_search.winding_topology`.
  * Per ``feedback_orbit_closure_discipline``: "it closed!" is the danger
    signal. The corrector residual is necessary; the topology cross-check
    is the orthogonal gate that catches a residual-passing orbit on the
    wrong basin.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.search.bifurcation_detector import monodromy
from cyclerfinder.search.binary_star_search import Topology, winding_topology
from cyclerfinder.search.cr3bp_general_periodic_3d import (
    FREE_VARS_FULL_ASYMMETRIC,
    RESIDUAL_FULL_STATE_AT_T,
    Periodic3DOrbit,
    correct_general_periodic_3d,
)


@dataclass(frozen=True)
class BranchedOrbit:
    """A branched orbit from an asymmetric / 3D corrector landing.

    Attributes
    ----------
    parent_state0, parent_period :
        The pre-bifurcation parent IC + period (the marginal eigenvector was
        computed at this state).
    parent_topology :
        (k1, k2) winding-topology of the parent at perturbation time.
    eigenvalue_used :
        The eigenvalue closest to +1 among the SECONDARY non-trivial pair at
        the parent (the bifurcation direction).
    eigenvector_used :
        The corresponding length-6 real right-eigenvector (unit-normalised).
    epsilon :
        The perturbation amplitude actually applied.
    sign :
        +1 or -1 — the eigenvector is sign-ambiguous; the corrector may have
        landed on either side of the parent. Both signs are tried; this
        records which sign produced the converged orbit.
    branched_orbit :
        The :class:`Periodic3DOrbit` returned by
        :func:`cyclerfinder.search.cr3bp_general_periodic_3d.correct_general_periodic_3d`.
    branched_topology :
        Winding topology of the converged branched orbit.
    topology_changed :
        ``branched_topology.k1 != parent_topology.k1 or
        branched_topology.k2 != parent_topology.k2``. The Phase 1 exit
        criterion (Section 4 of the design doc).
    """

    parent_state0: NDArray[np.float64]
    parent_period: float
    parent_topology: Topology
    eigenvalue_used: complex
    eigenvector_used: NDArray[np.float64]
    epsilon: float
    sign: int
    branched_orbit: Periodic3DOrbit
    branched_topology: Topology
    topology_changed: bool


def _select_saddle_center_eigenvector(
    monodromy_matrix: NDArray[np.float64],
) -> tuple[NDArray[np.float64], complex] | None:
    """Pick the marginal eigenvector at a saddle-center bifurcation.

    The saddle-center / pitchfork at lambda=+1 has its bifurcation direction
    along the eigenvector(s) of the monodromy whose eigenvalue is closest to
    +1 AMONG THE NON-TRIVIAL PAIR. We exclude the two trivial-pair eigenvalues
    (closest to +1; the energy + time-translation pair) and the two strongly
    unstable saddle-pair eigenvalues (largest |log|lambda||); the remaining 2
    are the secondary pair, and the secondary eigenvalue closest to +1 is the
    bifurcation direction.

    Returns
    -------
    (v, lam) :
        ``v`` = unit-normalised real right-eigenvector (length 6).
        ``lam`` = the matching eigenvalue.
        ``None`` if no secondary eigenvalue can be cleanly identified.

    Notes
    -----
    Mirrors the structure of
    :func:`cyclerfinder.genome.family_switch._select_period_multiplying_eigenvector`
    but tuned for k=1 (saddle-center). For a complex eigenvalue ``a + bi`` the
    eigenvector is generically complex; the perturbation direction is its
    real part. For a real eigenvalue the eigenvector is real.
    """
    eigvals, eigvecs = np.linalg.eig(monodromy_matrix)
    if eigvals.shape[0] < 4:
        return None
    # Step 1: exclude 2 trivial-pair eigenvalues by argsort(|lam - 1|).
    dists_to_one = np.abs(eigvals - 1.0)
    trivial_idx = set(int(i) for i in np.argsort(dists_to_one)[:2])
    # Step 2: exclude the 2 primary-saddle eigenvalues by largest |log|lam||.
    log_mags = np.array(
        [abs(math.log(max(abs(complex(e)), 1e-300))) for e in eigvals],
        dtype=np.float64,
    )
    primary_candidates = [i for i in range(eigvals.shape[0]) if i not in trivial_idx]
    primary_candidates.sort(key=lambda i: -log_mags[i])
    primary_idx = set(primary_candidates[:2])
    # Remaining indices are the secondary pair.
    secondary_idx = [
        i for i in range(eigvals.shape[0]) if i not in trivial_idx and i not in primary_idx
    ]
    if not secondary_idx:
        return None
    # Pick the secondary eigenvalue closest to +1.
    secondary_idx.sort(key=lambda i: dists_to_one[i])
    idx = secondary_idx[0]
    lam = complex(eigvals[idx])
    v_complex = eigvecs[:, idx]
    v_real = np.real(v_complex).astype(np.float64)
    norm = float(np.linalg.norm(v_real))
    if norm < 1e-14:
        v_real = np.imag(v_complex).astype(np.float64)
        norm = float(np.linalg.norm(v_real))
    if norm < 1e-14:
        return None
    v_real /= norm
    return v_real, lam


def branch_at_saddle_center(
    system: cr3bp.CR3BPSystem,
    parent_state0: NDArray[np.float64],
    parent_period: float,
    *,
    epsilon: float = 1e-3,
    tol: float = 1e-10,
    max_iter: int = 80,
    rtol: float = 1e-12,
    atol: float = 1e-12,
    independent_tol: float = 1e-6,
    require_monotone_decrease: bool = False,
) -> BranchedOrbit | None:
    """Branch off the parent symmetric family at a saddle-center bifurcation.

    Pipeline:

      1. Compute the parent's monodromy + Floquet eigenvectors.
      2. Pick the marginal eigenvector ``v`` (the secondary non-trivial pair's
         closest-to-+1 right-eigenvector — see
         :func:`_select_saddle_center_eigenvector`).
      3. Compute the parent's winding topology for the cross-check.
      4. Perturb the parent IC by ``+epsilon * v`` (and by ``-epsilon * v`` if
         the +sign converges to the parent's own topology / does not converge).
      5. Hand to
         :func:`cyclerfinder.search.cr3bp_general_periodic_3d.correct_general_periodic_3d`
         in **full-asymmetric mode**: free vars = (x0, y0, z0, xdot0, ydot0,
         zdot0, T); residual = full 6D state closure at T. Period guess = the
         parent's period (the saddle-center spawns a family with NEAR the
         parent's period, not k * parent's period — the k=1 special case is
         period-preserving).
      6. Compute the branched orbit's winding topology; compare to the parent's.

    Parameters
    ----------
    system :
        CR3BP system; only ``system.mu`` is read.
    parent_state0 :
        6-vector parent IC. The parent should be a SYMMETRIC family member
        (typically the bifurcation point's pre-bifurcation neighbour; the
        marginal eigenvector is freshly computed at THIS state).
    parent_period :
        Full nondim period of the parent.
    epsilon :
        Perturbation amplitude along the eigenvector (nondim state units).
        Too small: corrector lands back on the parent (the planar manifold is
        invariant if the eigenvector is in the (z, zdot) direction). Too
        large: the corrector diverges. Default 1e-3 matches family_switch's
        default; the design doc Section 5 budgets a sweep in P1.4 (this
        function is called per-epsilon).
    tol, max_iter, rtol, atol, independent_tol :
        Forwarded to :func:`correct_general_periodic_3d`.
    require_monotone_decrease :
        Forwarded to :func:`correct_general_periodic_3d`. Default ``False``
        (the damped-Newton pattern that matches ``correct_symmetric_nrho``).
        The strict monotone-decrease line search rejects the transient steps
        the damped Newton needs to settle into the quadratic basin at the
        post-bifurcation parent's basin; with the default the corrector
        converges from the perturbed C32 i=124 seed at eps=5e-4 to residual
        9.4e-12 in ~30 iterations. Set ``True`` for very ill-conditioned long-
        arc closures where transient residual increase would escape the basin.

    Returns
    -------
    BranchedOrbit | None :
        The converged + topology-checked branched orbit, or ``None`` if neither
        sign of the eigenvector produced a converged orbit.

    Notes
    -----
    Discipline (orbit-closure):

      * "Converged" means corrector residual < ``tol`` AND independent
        closure check (Radau re-propagation) < ``independent_tol``. This is
        :class:`Periodic3DOrbit`'s built-in compound gate; see its docstring.
      * Topology check is INDEPENDENT of the corrector (it runs `solve_ivp`
        with DOP853 over the closed orbit and counts windings). A residual-
        passing orbit whose topology matches the parent EXACTLY is not
        necessarily a failure — the saddle-center may bifurcate into two
        copies of the parent's (k1, k2) with different stability character.
        The :attr:`BranchedOrbit.topology_changed` flag is informational, not
        a hard gate.
      * ``None`` is returned on no-converge — never a "near-converged" orbit
        relabelled as passing.
    """
    parent_state0 = np.asarray(parent_state0, dtype=np.float64).copy()
    if parent_state0.shape != (6,):
        raise ValueError(f"parent_state0 must have shape (6,); got {parent_state0.shape}")
    if parent_period <= 0.0 or not math.isfinite(parent_period):
        raise ValueError(f"parent_period must be > 0 finite; got {parent_period}")
    if epsilon <= 0.0 or not math.isfinite(epsilon):
        raise ValueError(f"epsilon must be > 0 finite; got {epsilon}")

    mono = monodromy(system, parent_state0, parent_period, rtol=rtol, atol=atol)
    pick = _select_saddle_center_eigenvector(mono)
    if pick is None:
        return None
    v, lam = pick

    parent_topo = winding_topology(system.mu, parent_state0, parent_period)

    for sign in (+1, -1):
        perturbed = parent_state0 + sign * epsilon * v
        result = correct_general_periodic_3d(
            system,
            perturbed,
            parent_period,
            free_vars=FREE_VARS_FULL_ASYMMETRIC,
            residual_indices=RESIDUAL_FULL_STATE_AT_T,
            is_half_period_residual=False,
            tol=tol,
            max_iter=max_iter,
            rtol=rtol,
            atol=atol,
            independent_tol=independent_tol,
            require_monotone_decrease=require_monotone_decrease,
        )
        if not result.converged:
            continue
        # Independent topology cross-check (winding numbers).
        try:
            branched_topo = winding_topology(system.mu, result.state0, result.T_TU)
        except (RuntimeError, ValueError):
            continue
        topology_changed = bool(
            branched_topo.k1 != parent_topo.k1 or branched_topo.k2 != parent_topo.k2
        )
        return BranchedOrbit(
            parent_state0=parent_state0,
            parent_period=parent_period,
            parent_topology=parent_topo,
            eigenvalue_used=lam,
            eigenvector_used=v,
            epsilon=epsilon,
            sign=sign,
            branched_orbit=result,
            branched_topology=branched_topo,
            topology_changed=topology_changed,
        )
    return None
