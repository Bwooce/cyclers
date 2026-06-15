"""Family-switching corrector at a period-multiplying bifurcation (#266 Phase 3).

When a one-parameter NRHO family passes through a period-multiplying
bifurcation (a Floquet multiplier crosses a primitive k-th root of unity),
a NEW family with k times the parent period branches off at that point. The
classical genome-finding move is:

  1. Identify the bifurcation point (the parent family member where a
     multiplier sits on / near a primitive k-th root). The
     :mod:`cyclerfinder.search.bifurcation_detector` module provides the
     bracket; the parent member is one of ``bifurcation.members``.

  2. Compute the right-eigenvector ``v`` of the monodromy at the bifurcation
     point corresponding to the period-multiplying multiplier. This vector
     points along the branching family in state space.

  3. Perturb the parent IC ``s0`` along ``v`` by a small amplitude;
     simultaneously multiply the period by ``k``.

  4. Re-converge to a true periodic orbit at the new (kT, perturbed s0) guess
     -- this lands on the new family.

This module is the minimum-viable family-switcher for k=2 period-doubling. It
ALWAYS goes through the single-shooting symmetric corrector
:func:`cyclerfinder.search.nrho_continuation.correct_symmetric_nrho` first
because:

  - Single-shooting is the simpler, fewer-moving-parts approach;
  - For the L2 Southern NRHO period-doubling (k=2) the literature
    (Howell-Breakwell 1984; Koblick 2023 AMOSTECH) reports that single-shooting
    succeeds when seeded from the bifurcation point along the eigenvector;
  - Multi-shooting is escalated only when single-shooting fails (Phase 4).

Discipline (orbit-closure):

  * The switched family member is INDEPENDENTLY cross-checked by counting
    petals via :func:`cyclerfinder.genome.tulip.petal_count`. ``petal_count
    == k`` is the topological gate; the corrector's residual is the closure
    gate.
  * On no-converge, ``None`` is returned (NEVER a near-converged result
    relabeled as "passing" -- see the project memory entry
    ``feedback_orbit_closure_discipline``).
"""

from __future__ import annotations

import math

import numpy as np

import cyclerfinder.core.cr3bp as cr3bp
import cyclerfinder.genome.tulip as tulip
from cyclerfinder.search.bifurcation_detector import BifurcationPoint
from cyclerfinder.search.nrho_continuation import (
    SymmetricNRHO,
    correct_symmetric_nrho,
)

# ---------------------------------------------------------------------------
# Eigenvector selection.
# ---------------------------------------------------------------------------


def _select_period_multiplying_eigenvector(
    monodromy_matrix: np.ndarray,
    k: int,
) -> tuple[np.ndarray, complex, float] | None:
    """Pick the eigenvector of the monodromy whose eigenvalue is closest to a
    primitive k-th root of unity.

    Returns
    -------
    (v, lam, dist) :
        The eigenvector ``v`` (length-6 real vector), the matching eigenvalue
        ``lam`` (complex), and the distance from ``lam`` to its nearest
        primitive k-th root of unity. ``None`` if no eigenvalue is within a
        loose tolerance (``0.5``) of any primitive k-th root -- the bifurcation
        bracket may have been mislabeled.

    Notes
    -----
    For k=2 (period doubling) the primitive 2nd root of unity is ``-1``. The
    standard period-doubling multiplier sits NEAR ``-1`` at the bifurcation
    point; for the actual EIGENVECTOR we expect a real one (since both the
    monodromy and the multiplier are real). For higher k the multiplier is
    complex but the eigenvector is generally complex; we take the REAL part of
    the (complex-conjugate) eigenvector pair as the perturbation direction
    (the real subspace spanned by the complex pair).
    """
    eigvals, eigvecs = np.linalg.eig(monodromy_matrix)
    # Primitive k-th roots: e^(2 pi i j / k) for j coprime with k, 0 < j < k.
    roots: list[complex] = []
    for j in range(1, k):
        if math.gcd(j, k) == 1:
            roots.append(complex(math.cos(2 * math.pi * j / k), math.sin(2 * math.pi * j / k)))
    if not roots:
        return None
    best_dist = float("inf")
    best_idx = -1
    best_root = complex(1.0)
    for i, lam in enumerate(eigvals):
        for r in roots:
            d = abs(complex(lam) - r)
            if d < best_dist:
                best_dist = d
                best_idx = i
                best_root = r
    if best_idx < 0 or best_dist > 0.5:
        return None
    v_complex = eigvecs[:, best_idx]
    # Real direction in state space: for a real eigenvalue (k=2) the
    # eigenvector is real (up to a phase); for complex eigenvalues, take
    # the real part of the (normalised) eigenvector -- this spans the
    # 2-D real invariant subspace shared by the conjugate pair.
    v_real = np.real(v_complex).astype(np.float64)
    norm = float(np.linalg.norm(v_real))
    if norm < 1e-14:
        # Real part vanished; fall back to imaginary part.
        v_real = np.imag(v_complex).astype(np.float64)
        norm = float(np.linalg.norm(v_real))
    if norm < 1e-14:
        return None
    v_real /= norm
    del best_root  # not currently emitted
    return v_real, complex(eigvals[best_idx]), float(best_dist)


# ---------------------------------------------------------------------------
# Family switch.
# ---------------------------------------------------------------------------


def switch_family(
    parent: SymmetricNRHO,
    bifurcation: BifurcationPoint,
    system: cr3bp.CR3BPSystem,
    *,
    k: int,
    eigenvector_step: float = 1e-3,
    tol: float = 1e-9,
    max_iter: int = 80,
    rtol: float = 1e-12,
    atol: float = 1e-12,
    verify_petal_count: bool = True,
) -> SymmetricNRHO | None:
    """Branch off the parent family at the bifurcation, landing on a kT family.

    1. Read the parent's monodromy from ``parent.monodromy``; require it.
    2. Pick the eigenvector ``v`` of the monodromy whose eigenvalue is closest
       to a primitive k-th root of unity (see
       :func:`_select_period_multiplying_eigenvector`).
    3. Perturb the parent IC's ``(z0, ydot0)`` along the corresponding
       components of ``v`` by ``eigenvector_step``.
    4. Re-correct via the symmetric corrector
       (:func:`cyclerfinder.search.nrho_continuation.correct_symmetric_nrho`)
       with the period guess ``k * parent.T_TU`` and the fixed ``x0`` of the
       parent.
    5. (Optional) Independently verify ``petal_count(result) == k`` via the
       topological classifier.

    Parameters
    ----------
    parent :
        A converged parent family member; must carry a non-None ``monodromy``.
    bifurcation :
        The bracketed bifurcation. Currently only ``k`` is used (the eigenvector
        is freshly computed at the parent for numerical clarity).
    system :
        CR3BP system; only ``system.mu`` is read.
    k :
        Period-multiplying integer (k=2 doubles the period). Should match
        ``bifurcation.k`` -- callers may override for testing.
    eigenvector_step :
        Amplitude of the IC perturbation along the eigenvector (in nondim
        state-vector units). Small enough that the perturbed IC sits near the
        parent in state space; large enough that the Newton corrector is pulled
        off the parent family and onto the new k-branch.
    tol :
        Closure tolerance forwarded to the corrector.
    max_iter, rtol, atol :
        Forwarded to the corrector.
    verify_petal_count :
        If True, the result is independently checked via
        :func:`cyclerfinder.genome.tulip.petal_count` -- and ``None`` is
        returned if ``petal_count != k`` (topological mismatch). The corrector's
        residual gate is necessary but not sufficient: the corrector may
        converge to a *different* periodic orbit (e.g. a kT orbit on a sibling
        family) that happens to satisfy the perpendicular crossing residual.
        The petal_count classifier is the independent cross-check.

    Returns
    -------
    SymmetricNRHO | None :
        The converged k-period family member or ``None`` on failure.

    Notes
    -----
    Single-shooting first. Multi-shooting (one segment per petal; patch-point
    continuity) is the textbook upgrade when single-shooting fails on high-k
    bifurcations; this module deliberately ships only the single-shooting
    path, with multi-shooting deferred to Phase 4 per the task spec.
    """
    if k < 2:
        raise ValueError(f"switch_family: k must be >= 2, got {k}")
    if parent.monodromy is None:
        raise ValueError(
            "switch_family: parent.monodromy is None; correct the parent with "
            "with_monodromy=True (default) before family-switching."
        )

    pick = _select_period_multiplying_eigenvector(parent.monodromy, k)
    if pick is None:
        return None
    v, _lam, _dist = pick
    # Defensive: parent state vector layout.
    # state = (x, y, z, xdot, ydot, zdot). For the symmetric IC we hold y=0,
    # xdot=0, zdot=0; the eigenvector's z (index 2) and ydot (index 4)
    # components drive the (z0, ydot0) IC perturbation.
    dz0 = float(v[2]) * eigenvector_step
    dydot0 = float(v[4]) * eigenvector_step

    # Perturb and multiply the period.
    z0_new = parent.z0 + dz0
    ydot0_new = parent.ydot0 + dydot0
    period_guess = float(k) * parent.T_TU

    member = correct_symmetric_nrho(
        system,
        parent.x0,
        z0_new,
        ydot0_new,
        period_guess,
        tol=tol,
        max_iter=max_iter,
        rtol=rtol,
        atol=atol,
        with_monodromy=True,
    )
    if not member.converged:
        # Try the opposite sign of the eigenvector (eigenvectors are sign-
        # ambiguous; the branch may lie on the OTHER side).
        z0_new = parent.z0 - dz0
        ydot0_new = parent.ydot0 - dydot0
        member = correct_symmetric_nrho(
            system,
            parent.x0,
            z0_new,
            ydot0_new,
            period_guess,
            tol=tol,
            max_iter=max_iter,
            rtol=rtol,
            atol=atol,
            with_monodromy=True,
        )
    if not member.converged:
        return None

    # Independent topological gate: the petal count must equal k. Otherwise
    # the corrector landed on a kT orbit of the PARENT family (a k-fold cover
    # of the parent), not the new k-branch.
    if verify_petal_count:
        s0 = np.array([member.x0, 0.0, member.z0, 0.0, member.ydot0, 0.0], dtype=np.float64)
        try:
            n = tulip.petal_count(s0, member.T_TU, system, rtol=rtol, atol=atol)
        except RuntimeError:
            return None
        if n != k:
            return None
    return member
