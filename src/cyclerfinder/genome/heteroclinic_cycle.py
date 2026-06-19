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
    """Real-cast + unit-normalise a (possibly complex) eigenvector (4-vector)."""
    vr = np.real(v)
    n = float(np.linalg.norm(vr))
    if n < 1e-14:
        vr = np.real(v) + np.imag(v)
        n = float(np.linalg.norm(vr))
    return (vr / n).astype(np.float64) if n > 0.0 else vr.astype(np.float64)


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
