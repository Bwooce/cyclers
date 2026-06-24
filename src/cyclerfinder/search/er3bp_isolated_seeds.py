"""e=0 (circular RTBP) symmetric resonant periodic-orbit seed base (#440 Phase 1).

Generalises the GO-gate-validated 3/2 recipe (``correct_symmetric_fixed_jacobi``
at mu=0.001) to all five interior mean-motion resonances (MMRs) tabulated by
Antoniadou & Libert 2018 (CMDA, DOI 10.1007/s10569-018-9834-8). These circular
(e=0) members are the seed base for the later isolated-family hunt (Phase 2,
which continues each into the pulsating eccentricity e via the #437 fold-aware
continuator). This module is PURE CR3BP — it uses the fixed-Jacobi symmetric
CR3BP corrector, NOT a free-period corrector (the latter is a family-selection
trap that drifts to the exterior 1:2 / L4-L5 family).

Validated recipe (verbatim from the GO gate, plan
``2026-06-25-440-isolated-er3bp-family-hunt-plan.md``):

- Build ``CR3BPSystem(mu=0.001, ...)``.
- For interior MMR p/q with semi-major axis a1: seed ``x0 = -mu + a1``,
  ``ydot0 = sqrt((1-mu)/a1) - a1`` (inertial circular speed minus synodic
  rotation), ``jacobi = jacobi_constant([x0,0,0,0,ydot0,0], mu)``.
- Converge with ``correct_symmetric_fixed_jacobi(system, x0, jacobi,
  period_guess=T0, ydot0_sign=+1, half_crossings=1, tol=1e-11)`` where
  ``T0 = 2*pi / (a1**-1.5 - 1)`` (the synodic/resonant period; 4*pi for 3/2).
  T0 is the linearized bifurcation period — the converged member sits at a
  nearby, finite-amplitude-shifted T, which is expected; do NOT pin T.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.search.cr3bp_periodic import (
    SymmetricOrbit,
    correct_symmetric_fixed_jacobi,
)

# Interior-MMR semi-major axes (Antoniadou & Libert 2018; see the plan + digest
# ``2026-06-25-digest-antoniadou-libert-2018-circ-ellip-continuation.md``).
# Cross-check: for p/q interior resonance a1 = (q/p)**(2/3).
MMR_SEMI_MAJOR_AXES: tuple[tuple[int, int, float], ...] = (
    (3, 2, 0.763143),
    (5, 2, 0.5428),
    (3, 1, 0.4807),
    (4, 1, 0.3968),
    (5, 1, 0.3419),
)


@dataclass(frozen=True)
class ResonantSeed:
    """A recovered e=0 symmetric resonant periodic-orbit seed at mu.

    ``state0`` is the perpendicular-x-axis-crossing IC ``[x0,0,0,0,ydot0,0]``;
    ``a_helio`` is the osculating heliocentric semi-major axis at the IC (the
    physically-correct resonance measure — for the eccentric 3/2 the crossing
    radius differs from a by ~4%, but the energy-based osculating a tracks a1 to
    ~0.5%); ``n_ratio = a_helio**-1.5`` is the mean-motion ratio that must equal
    the resonance p/q.
    """

    label: str
    p: int
    q: int
    a1: float
    state0: np.ndarray
    period: float
    jacobi: float
    a_helio: float
    n_ratio: float
    converged: bool
    crossing_residual: float


def _osculating_a_helio(x0: float, ydot0: float, mu: float) -> float:
    """Osculating heliocentric semi-major axis of the IC ``[x0,0,0,0,ydot0,0]``.

    The primary (star) sits at ``-mu``, so the heliocentric radius is
    ``r = x0 + mu``. In the synodic (rotating, omega=1) frame the velocity is
    purely ``ydot0`` in y; the inertial heliocentric tangential speed adds the
    frame rotation ``omega * r``, giving ``v_inert = ydot0 + r``. With
    ``mu_helio = 1 - mu`` the vis-viva energy ``E = v^2/2 - mu_helio/|r|`` fixes
    ``a = -mu_helio / (2 E)``. This is the resonance-defining mean orbit element,
    not the instantaneous crossing radius.
    """
    r = x0 + mu
    v_inert = ydot0 + r  # omega = 1 in nondimensional CR3BP
    mu_helio = 1.0 - mu
    energy = 0.5 * v_inert * v_inert - mu_helio / abs(r)
    return float(-mu_helio / (2.0 * energy))


def resonant_po_seed(
    system: cr3bp.CR3BPSystem,
    p: int,
    q: int,
    a1: float,
    *,
    tol: float = 1e-11,
) -> ResonantSeed:
    """Converge the e=0 symmetric interior-MMR periodic orbit for one resonance.

    Applies the GO-gate recipe (module docstring) for resonance ``p/q`` at
    semi-major axis ``a1`` and returns the recovered :class:`ResonantSeed`.
    ``.converged`` mirrors the corrector's crossing-residual convergence flag;
    callers should inspect it rather than assume success.
    """
    mu = system.mu
    x0 = -mu + a1
    ydot0_seed = math.sqrt((1.0 - mu) / a1) - a1
    seed_state = np.array([x0, 0.0, 0.0, 0.0, ydot0_seed, 0.0])
    jacobi = cr3bp.jacobi_constant(seed_state, mu)
    t0 = 2.0 * math.pi / (a1**-1.5 - 1.0)

    orbit: SymmetricOrbit = correct_symmetric_fixed_jacobi(
        system,
        x0,
        jacobi,
        period_guess=t0,
        ydot0_sign=1.0,
        half_crossings=1,
        tol=tol,
    )

    state0 = np.array([orbit.x0, 0.0, 0.0, 0.0, orbit.ydot0, 0.0])
    a_helio = _osculating_a_helio(orbit.x0, orbit.ydot0, mu)
    n_ratio = a_helio**-1.5
    return ResonantSeed(
        label=f"{p}:{q}",
        p=p,
        q=q,
        a1=a1,
        state0=state0,
        period=orbit.period,
        jacobi=orbit.jacobi,
        a_helio=a_helio,
        n_ratio=n_ratio,
        converged=orbit.converged,
        crossing_residual=orbit.crossing_residual,
    )


def all_mmr_seeds(*, mu: float = 0.001, tol: float = 1e-11) -> list[ResonantSeed]:
    """Converge all five interior-MMR e=0 resonant POs at ``mu``.

    Does NOT raise on a single non-convergence: a failing MMR is returned with
    ``.converged == False`` so the caller can see exactly which resonances
    landed. Builds one ``CR3BPSystem`` at ``mu`` (unit scales — pure CR3BP
    dynamics; the physical l_km/t_s are irrelevant to the nondimensional
    resonant geometry).
    """
    system = cr3bp.CR3BPSystem(mu=mu, primary="Sun", secondary="planet", l_km=1.0, t_s=1.0)
    seeds: list[ResonantSeed] = []
    for p, q, a1 in MMR_SEMI_MAJOR_AXES:
        seeds.append(resonant_po_seed(system, p, q, a1, tol=tol))
    return seeds
