"""Generalized L1/L2 halo-family-member builder, parameterized by Jacobi constant (#613).

Motivation
----------
`#548`->`#553`->`#555`->`#556`->`#612` all tested the Owen & Baresi
latitudinal-frequency reproduction question (L1 target rotation number
0.2739, L2 target 0.02163) at a SINGLE, hardcoded Jacobi constant, C=3.15.
`#555` proved the L1 quasi-halo rotation number at C=3.15 is energy-pinned
near 0.074 (flat in amplitude) -- NOT 0.2739 -- but never tested whether a
DIFFERENT C (still inside the physically relevant L1 halo regime, below its
Neimark-Sacker-adjacent planar-to-halo bifurcation at C~3.1745) gives a
different rotation number. Every halo-family builder in this codebase before
this module (`tests/search/test_variational_qp_torus.py`'s private
``_l1_halo_at_315``/``_l2_halo_at_315``, ``scripts/run_555_owen_baresi_c315_final.py``'s
``l1_halo_at_c``/``l2_halo_at_c``) was hardcoded to continue FROM a C=3.15-ish
seed TOWARD C=3.15/the bifurcation (a single continuation direction). This
module generalizes those to an ARBITRARY target Jacobi constant, in EITHER
direction (toward or away from the bifurcation), by picking the natural-
parameter continuation sign automatically.

Method (unchanged from the prior scripts/tests -- productized, not redesigned)
-------------------------------------------------------------------------
Both builders use :func:`cyclerfinder.search.nrho_continuation.correct_symmetric_nrho`
(fixed-``x0`` perpendicular-crossing symmetric corrector) in natural-parameter
continuation on ``x0``, tracking the Jacobi constant at each converged member,
and stop when the family's Jacobi constant brackets the target (then refine
with a few finer sub-steps). The L1 seed is the same near-C=3.145 halo IC used
throughout `#548`-`#612`; the L2 seed is built via the `#553`-scoped
planar-Lyapunov -> halo bifurcation SEED GENERATOR
(:func:`cyclerfinder.search.cr3bp_periodic.correct_symmetric_fixed_jacobi`)
at a reference Jacobi constant near the L2 family's own small-amplitude
regime, then continued in ``x0`` to the target (the small-``z0`` step-off
that works AT its own bifurcation does NOT work far from it -- confirmed
empirically in this task: step-off from a Lyapunov built directly at a
far-from-bifurcation target C collapses back to the planar orbit).

``halo_center_pair`` and ``best_resonance_k`` extract the Neimark-Sacker-type
complex-conjugate Floquet pair and its nearest-primitive-root resonance
order -- the same pattern as the private test/script helpers, generalized to
take the built halo directly.

Discipline
----------
* No catalogue writeback from this module; it is a family-continuation
  utility, not a discovery result.
* The returned Floquet-phase rotation number is a LINEAR (infinitesimal-
  amplitude) estimate -- cheap, used for broad mapping. Confirming a specific
  candidate Jacobi constant against a genuine target frequency requires
  building an actual small-amplitude torus with
  :func:`cyclerfinder.genome.qp_tori.correct_qp_torus` (not duplicated here).
"""

from __future__ import annotations

import math

import numpy as np
from numpy.typing import NDArray

import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.search.bifurcation_detector import floquet_multipliers, monodromy
from cyclerfinder.search.cr3bp_periodic import correct_symmetric_fixed_jacobi
from cyclerfinder.search.nrho_continuation import SymmetricNRHO, correct_symmetric_nrho

# The C~3.145 L1 halo IC used throughout #548-#612 (near, but not at, the
# genuine C=3.15 target -- the continuation below walks it to whatever
# target_jacobi the caller asks for).
_L1_SEED_X0 = 0.82431
_L1_SEED_Z0 = -0.06058
_L1_SEED_YDOT0 = 0.17154
_L1_SEED_T = 2.7647

# The #553-scoped L2 planar-Lyapunov seed guesses (near-libration family;
# T~3.42 branch, NOT the large-amplitude NRHO branch).
_L2_LYAP_GUESSES = ((1.1225, 3.43), (1.1204, 3.42), (1.1180, 3.40), (1.1250, 3.45))
_L2_REFERENCE_JACOBI = 3.15
_L2_SEED_Z0_LADDER = (0.006, 0.008, 0.012, 0.016, 0.02)


def _continue_to_jacobi(
    system: cr3bp.CR3BPSystem,
    seed: SymmetricNRHO,
    target_jacobi: float,
    *,
    dx0: float,
    max_steps: int,
    refine_steps: int,
) -> SymmetricNRHO | None:
    """Natural-parameter (``x0``) continuation of a symmetric-halo family until
    the Jacobi constant brackets ``target_jacobi``, then refine.

    Direction is chosen by PROBING the local sign of ``dC/dx0`` (one trial
    step of ``+abs(dx0)``) rather than assumed: this task found the L1 and L2
    families have OPPOSITE conventions (decreasing ``x0`` increases Jacobi for
    the L1 family; decreasing ``x0`` DECREASES Jacobi for the L2 family, the
    reverse), so a single hardcoded sign is wrong for one of the two families.
    """
    if abs(seed.jacobi - target_jacobi) < 1e-9:
        return seed
    probe = correct_symmetric_nrho(
        system, seed.x0 + abs(dx0), seed.z0, seed.ydot0, seed.T_TU, with_monodromy=False
    )
    if not probe.converged:
        return None
    dc_dx0_positive = probe.jacobi > seed.jacobi
    want_increase = target_jacobi > seed.jacobi
    direction = 1.0 if dc_dx0_positive == want_increase else -1.0
    step = direction * abs(dx0)
    x0, z0, ydot0, t_period = seed.x0, seed.z0, seed.ydot0, seed.T_TU
    prev_jacobi = seed.jacobi
    for _ in range(max_steps):
        x0 += step
        member = correct_symmetric_nrho(
            system, float(x0), z0, ydot0, t_period, with_monodromy=False
        )
        if not member.converged:
            return None
        if (prev_jacobi - target_jacobi) * (member.jacobi - target_jacobi) <= 0:
            return _refine_to_jacobi(
                system, x0 - step, x0, z0, ydot0, t_period, target_jacobi, refine_steps
            )
        prev_jacobi = member.jacobi
        z0, ydot0, t_period = member.z0, member.ydot0, member.T_TU
    return None


def _refine_to_jacobi(
    system: cr3bp.CR3BPSystem,
    x0_lo: float,
    x0_hi: float,
    z0: float,
    ydot0: float,
    t_period: float,
    target_jacobi: float,
    refine_steps: int,
) -> SymmetricNRHO | None:
    """Sub-divide the bracket ``[x0_lo, x0_hi]`` into ``refine_steps`` finer
    natural-parameter steps, returning as soon as the finer stepping brackets
    ``target_jacobi`` (or the last converged member if it never re-brackets --
    this happens right at a pitchfork bifurcation, where ``dC/dx0`` diverges
    and a single ORIGINAL-scale step can jump past the target by more than the
    whole fine-stepping range spans; the returned member is then the closest
    achievable approximation at this ``dx0`` resolution, not an exact hit)."""
    dx0_fine = (x0_hi - x0_lo) / max(refine_steps, 1)
    x0 = x0_lo
    prev_jacobi = correct_symmetric_nrho(
        system, float(x0), z0, ydot0, t_period, with_monodromy=False
    ).jacobi
    member: SymmetricNRHO | None = None
    for _ in range(refine_steps):
        x0 += dx0_fine
        member = correct_symmetric_nrho(
            system, float(x0), z0, ydot0, t_period, with_monodromy=False
        )
        if not member.converged:
            return None
        if (prev_jacobi - target_jacobi) * (member.jacobi - target_jacobi) <= 0:
            return member
        prev_jacobi = member.jacobi
        z0, ydot0, t_period = member.z0, member.ydot0, member.T_TU
    return member


def l1_halo_at_jacobi(
    system: cr3bp.CR3BPSystem,
    target_jacobi: float,
    *,
    dx0: float = 3e-5,
    max_steps: int = 3000,
    refine_steps: int = 10,
) -> SymmetricNRHO | None:
    """Build the EM L1 halo family member at ``target_jacobi`` (generalizes
    the C=3.15-hardcoded ``_l1_halo_at_315``/``l1_halo_at_c`` helpers in
    `#555`'s script and `#612`'s test file to an arbitrary target, in either
    continuation direction).

    Returns ``None`` if the seed itself fails to converge, or the family's
    continuation loses convergence before bracketing ``target_jacobi`` (e.g.
    ``target_jacobi`` sits past the family's physical extent).
    """
    seed = correct_symmetric_nrho(
        system, _L1_SEED_X0, _L1_SEED_Z0, _L1_SEED_YDOT0, _L1_SEED_T, with_monodromy=False
    )
    if not seed.converged:
        return None
    return _continue_to_jacobi(
        system, seed, target_jacobi, dx0=dx0, max_steps=max_steps, refine_steps=refine_steps
    )


def _l2_halo_at_reference(
    system: cr3bp.CR3BPSystem, reference_jacobi: float
) -> SymmetricNRHO | None:
    """The `#553`-scoped planar-Lyapunov -> halo bifurcation seed generator,
    at ``reference_jacobi`` (must be close enough to the L2 family's own
    bifurcation, C~3.1521, for the small-``z0`` step-off ladder to land on the
    genuine small-amplitude 3D halo -- confirmed to work at 3.15, confirmed in
    this task NOT to work at e.g. 3.076, far from the bifurcation)."""
    lyap = None
    for xg, tg in _L2_LYAP_GUESSES:
        cand = correct_symmetric_fixed_jacobi(
            system,
            x0_guess=xg,
            jacobi=reference_jacobi,
            period_guess=tg,
            ydot0_sign=+1.0,
            half_crossings=1,
            tol=1e-10,
            x0_bounds=(1.10, 1.20),
        )
        if cand.converged and abs(cand.period - 3.42) < 0.15:
            lyap = cand
            break
    if lyap is None:
        return None
    for z0_seed in _L2_SEED_Z0_LADDER:
        member = correct_symmetric_nrho(
            system,
            float(lyap.x0),
            z0_seed,
            float(lyap.ydot0),
            float(lyap.period),
            with_monodromy=False,
        )
        if (
            member.converged
            and abs(member.z0) > 3e-3
            and abs(member.jacobi - reference_jacobi) < 3e-3
        ):
            return member
    return None


def l2_halo_at_jacobi(
    system: cr3bp.CR3BPSystem,
    target_jacobi: float,
    *,
    reference_jacobi: float = _L2_REFERENCE_JACOBI,
    dx0: float = 3e-5,
    max_steps: int = 3000,
    refine_steps: int = 10,
) -> SymmetricNRHO | None:
    """Build the EM L2 halo family member at ``target_jacobi``.

    First builds the genuine small-amplitude 3D halo at ``reference_jacobi``
    (default 3.15, where the bifurcation step-off is known to work -- see
    `#555`), then continues in ``x0`` to ``target_jacobi`` if it differs.
    Returns ``None`` if either stage fails.
    """
    seed = _l2_halo_at_reference(system, reference_jacobi)
    if seed is None:
        return None
    return _continue_to_jacobi(
        system, seed, target_jacobi, dx0=dx0, max_steps=max_steps, refine_steps=refine_steps
    )


def halo_center_pair(
    system: cr3bp.CR3BPSystem, halo: SymmetricNRHO
) -> tuple[NDArray[np.float64], list[complex]]:
    """Return ``(state0, candidates)``: the halo's IC and the Floquet
    multipliers that are neither the trivial ``(1,1)`` pair nor real (i.e. the
    Neimark-Sacker-type complex-conjugate "center" pair used to seed a GMOS
    quasi-halo torus). ``candidates`` has length 0 or 2 (a conjugate pair)."""
    state0 = np.array([halo.x0, 0.0, halo.z0, 0.0, halo.ydot0, 0.0])
    mono = monodromy(system, state0, halo.T_TU)
    eigs = floquet_multipliers(mono)
    cands = [complex(e) for e in eigs if abs(e - 1.0) > 1e-3 and abs(e.imag) > 1e-4]
    return state0, cands


def best_resonance_k(phi: float, *, k_max: int = 80) -> int:
    """Nearest primitive ``k``-th-root-of-unity order to angle ``phi``
    (radians) -- the ``k`` :func:`cyclerfinder.genome.qp_tori.correct_qp_torus`
    needs to seed its Neimark-Sacker eigenvector search."""
    best_k, best_dist = 5, math.inf
    for k in range(3, k_max + 1):
        for j in range(1, k):
            if math.gcd(j, k) == 1:
                dist = abs(phi - 2 * math.pi * j / k)
                if dist < best_dist:
                    best_dist, best_k = dist, k
    return best_k


def linear_rotation_number(candidates: list[complex]) -> float | None:
    """Infinitesimal-amplitude rotation-number ESTIMATE: ``|arg(lambda)| / 2pi``
    of the first center-pair candidate. Cheap (no torus build) -- the
    quantity this task's broad ``C``-mapping sweeps, confirmed at specific
    candidate ``C`` values against an actual small-amplitude
    ``correct_qp_torus`` build (not duplicated here). Returns ``None`` if no
    center pair is available (halo has no genuine unit-circle complex pair)."""
    if not candidates:
        return None
    e = candidates[0]
    return abs(math.atan2(e.imag, e.real)) / (2 * math.pi)


__all__ = [
    "best_resonance_k",
    "halo_center_pair",
    "l1_halo_at_jacobi",
    "l2_halo_at_jacobi",
    "linear_rotation_number",
]
