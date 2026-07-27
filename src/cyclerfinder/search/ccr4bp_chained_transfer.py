"""CCR4BP chained-transfer glue: Callisto-L1 <-> Ganymede <-> Europa-L2 (`#715`).

Second independent positive control for the CCR4BP heteroclinic-search
pipeline (`#689`-`#694`), grounded in Aryan & Fitzgerald (2024), "Four Body
Invariant Structures And Chaos Analysis For Jovian Multi-Moon Ballistic
Transfers," AAS 24-103 (digested at `#710`,
``docs/notes/2026-07-26-710-digest-aryan-fitzgerald-2024-jovian-pccfbp.md``).

That paper's Section 4 (see the digest) reports a genuine, PUBLISHED
transit trajectory: a Jupiter-Callisto-Ganymede L1 UNSTABLE manifold reaches
a Ganymede rendezvous after ~60 days, then a Jupiter-Europa-Ganymede L2
STABLE manifold (backward-propagated) reaches Europa after ~74 days -- a
genuine two-hop, Ganymede-mediated manifold CHAIN across TWO different
``CCR4BPSystem`` instances (base=Callisto/perturber=Ganymede, then
base=Europa/perturber=Ganymede), not a single-system tube-to-tube near-
coincidence the way `#694`'s own ``ccr4bp_heteroclinic_search.py`` is scoped
(``ManifoldTube`` pairs within ONE system).

What this module does
----------------------
1. :func:`lyapunov_seed_at_jacobi` -- a GENERIC (any CR3BP ``mu``) planar
   Lyapunov-orbit seed at a caller-specified Jacobi constant, built entirely
   by REUSING two already-existing, unmodified project functions:
   :func:`cyclerfinder.search.periapse_map.collinear_l1_l2` (L1/L2 x-location)
   and :func:`cyclerfinder.search.cr3bp_periodic.correct_symmetric_fixed_jacobi`
   (Ross & Roberts-Tsoukkas fixed-Jacobi symmetric corrector). Every prior
   CCR4BP torus this project has built (`#690`'s JEG 3:4, `#696`'s Io-
   Ganymede 4:1, ...) seeds from a DISTANT spacecraft:moon resonant orbit,
   not a genuine L1/L2 Lyapunov orbit -- this is the first CCR4BP torus in
   this project built from an actual libration-point seed, needed here
   because the paper's own claim is specifically about L1/L2 manifolds
   (a distant resonant-orbit substitute was tried first during this task's
   own development and found to sit at a fixed "safe" radius the manifold
   provably cannot leave within any reasonable propagation window -- see the
   `#715` commit history / test module docstrings for that ruled-out
   approach).
2. :func:`ganymede_synodic_position` -- Ganymede's (x, y) position in a
   CCR4BP system's OWN base-moon synodic frame (re-implements the same
   circle parametrisation as ``core.ccr4bp._ganymede_position``, which is
   private to that module; duplicated locally rather than imported, mirroring
   ``ccr4bp_manifold_globalize.py``'s own precedent of not treating another
   module's private helpers as a supported cross-module contract -- cross-
   checked for exact equality in this module's own test file).
3. :func:`scan_ganymede_approach` -- scans an already-globalized
   :class:`~cyclerfinder.search.ccr4bp_manifold_globalize.ManifoldTube`'s full
   ``(theta2, t)`` grid for its closest planar approach to Ganymede's
   synodic position (evaluated at the SAME absolute time each grid state was
   sampled at), returning a :class:`GanymedeApproach`. Works for EITHER leg
   (Callisto-side unstable tube or Europa-side stable tube) since both
   systems carry Ganymede as their own perturber with its own tracked
   synodic position.
4. :func:`to_physical_km_days` -- converts a :class:`GanymedeApproach` to
   physical units, given the caller-supplied length unit (km, the relevant
   system's OWN base-moon SMA) and base-moon orbital period (days) -- kept as
   explicit caller-supplied conversion factors (not re-deriving `#694`'s own
   Europa-hardcoded ``_L_KM``/``_v_unit_km_s``) so this module works
   correctly for BOTH the Callisto-based and Europa-based systems without
   inheriting that documented generality gap.

Honest result (read the full account in `#715`'s own commit message / final
report): steps 1-3 all succeed cleanly on BOTH systems, at the paper's own
exact quoted Jacobi constants (Callisto C=3.0044, Europa C=3.0034) -- a
genuine, independently-verified SECOND positive control for the CCR4BP torus
corrector specifically. The CHAIN itself (a genuine close Ganymede
encounter on either leg) was NOT found within this task's own search scope
(a modest theta2-phase grid, propagated up to ~90 days) -- see
``tests/search/test_ccr4bp_chained_transfer.py`` for the pinned, honestly-
reported negative bound on both legs' closest approach. This is consistent
with the `#710` digest's own repeated caution that the paper's specific
transit result is a carefully-located trajectory in a large search space,
"not an independently-sourced target this project could grade against
without re-deriving their exact seed" -- reproducing the SPECIFIC transport
route would need either a much finer/adaptive phase search (the paper's own
Poincare-map methodology) or substantially more compute than this task
spent, not a methodological wall in the pipeline itself (the pipeline's own
mechanics -- torus convergence, manifold globalization, Ganymede-distance
evaluation -- all run cleanly on this new, previously-unexercised
INTERIOR-perturber system).

Scope
-----
Pure post-hoc glue module. Does NOT modify ``core/ccr4bp.py``,
``search/variational_ccr4bp_torus.py``, ``search/ccr4bp_manifold_globalize.py``,
``search/ccr4bp_whisker.py``, ``search/cr3bp_periodic.py``, or
``search/periapse_map.py`` -- only calls their public functions. No catalogue
writeback (a capability-proof + discovery-attempt module, not a vetted
result).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.core.ccr4bp import CCR4BPSystem
from cyclerfinder.search.ccr4bp_manifold_globalize import ManifoldTube
from cyclerfinder.search.cr3bp_periodic import SymmetricOrbit, correct_symmetric_fixed_jacobi
from cyclerfinder.search.periapse_map import collinear_l1_l2


@dataclass(frozen=True)
class LyapunovSeed:
    """A converged planar Lyapunov-orbit seed at a fixed Jacobi constant."""

    state0: NDArray[np.float64]  # (x0, 0, 0, 0, ydot0, 0)
    period: float
    jacobi: float
    l_point: str  # "L1" or "L2"
    orbit: SymmetricOrbit  # the underlying corrector result, for diagnostics


def lyapunov_seed_at_jacobi(
    mu: float,
    target_jacobi: float,
    l_point: str,
    *,
    x0_offset: float = 0.005,
    period_guess: float = 3.0,
    ydot0_sign: float = 1.0,
    tol: float = 1e-10,
    x0_bracket_halfwidth: float = 0.3,
    max_iter: int = 60,
) -> LyapunovSeed:
    """Planar Lyapunov-orbit seed at ``l_point`` ("L1" or "L2") with Jacobi
    constant ``target_jacobi``, generic over ``mu`` (works for ANY CR3BP
    mass ratio, not just the pairs this project has already built).

    Thin wrapper around :func:`cyclerfinder.search.periapse_map.collinear_l1_l2`
    (locates L1/L2) and
    :func:`cyclerfinder.search.cr3bp_periodic.correct_symmetric_fixed_jacobi`
    (the actual corrector) -- reused UNMODIFIED. ``x0_offset`` is the initial
    displacement from the bare libration point (a small-to-moderate-amplitude
    Lyapunov seed; the corrector walks it to the exact ``target_jacobi``
    regardless of the seed amplitude being only approximate).
    ``ydot0_sign``/``period_guess`` select which family branch the corrector
    lands on when more than one converges near the same seed (the L1/L2
    Lyapunov family generally has more than one topologically-distinct
    member reachable from a given seed region -- callers building a NEW
    system should verify the returned orbit is the intended simple,
    single-loop small-to-moderate-amplitude branch, e.g. by checking
    ``period`` against the expected order-few nondimensional-time-unit scale,
    not a resonant multiple of it).

    Raises
    ------
    RuntimeError
        If the corrector fails to converge to ``target_jacobi`` from the
        given seed.
    """
    if l_point not in ("L1", "L2"):
        raise ValueError(f"l_point must be 'L1' or 'L2', got {l_point!r}")
    l1_x, l2_x = collinear_l1_l2(mu)
    l_x = l1_x if l_point == "L1" else l2_x
    sign = -1.0 if l_point == "L1" else 1.0
    x0_guess = l_x + sign * x0_offset
    # primary/secondary strings are cosmetic (log_outcome metadata only);
    # l_km/t_s are unused by the corrector's own math (mu is the only system
    # field the corrector actually reads).
    system = cr3bp.CR3BPSystem(mu=mu, primary="_generic", secondary="_generic", l_km=1.0, t_s=1.0)
    orbit = correct_symmetric_fixed_jacobi(
        system,
        x0_guess,
        target_jacobi,
        period_guess,
        ydot0_sign=ydot0_sign,
        tol=tol,
        max_iter=max_iter,
        x0_bounds=(l_x - x0_bracket_halfwidth, l_x + x0_bracket_halfwidth),
    )
    if not orbit.converged:
        raise RuntimeError(
            f"lyapunov_seed_at_jacobi: corrector did not converge "
            f"(l_point={l_point}, mu={mu}, target_jacobi={target_jacobi}, "
            f"x0_guess={x0_guess}, crossing_residual={orbit.crossing_residual:.3e})"
        )
    state0 = np.array([orbit.x0, 0.0, 0.0, 0.0, orbit.ydot0, 0.0], dtype=np.float64)
    return LyapunovSeed(
        state0=state0, period=orbit.period, jacobi=orbit.jacobi, l_point=l_point, orbit=orbit
    )


def ganymede_synodic_position(t: float, system: CCR4BPSystem) -> tuple[float, float]:
    """Ganymede's ``(x, y)`` position in ``system``'s own base-moon synodic
    frame at absolute nondimensional time ``t``.

    Re-implements the SAME two-line circle parametrisation as
    ``core.ccr4bp._ganymede_position`` (private to that module -- see module
    docstring for why this is duplicated rather than imported). Cross-checked
    for exact equality against that private function in
    ``tests/search/test_ccr4bp_chained_transfer.py``.
    """
    theta = system.theta_gan0 + system.omega_gan * t
    return system.a_gan * math.cos(theta), system.a_gan * math.sin(theta)


@dataclass(frozen=True)
class GanymedeApproach:
    """Closest planar approach of a :class:`ManifoldTube` to Ganymede."""

    dmin: float  # nondimensional (system's own base-moon SMA units)
    theta2: float  # torus phase the closest-approach departure came from
    t_flow: float  # elapsed flow time since departure (>= 0), tube's own convention
    r_min: float  # minimum spacecraft radius from the system origin over the SAME scan
    n_valid_departures: int  # how many theta2 departures actually propagated (diagnostic)


def scan_ganymede_approach(tube: ManifoldTube, system: CCR4BPSystem) -> GanymedeApproach:
    """Scan ``tube``'s full ``(theta2, t)`` grid for its closest planar
    approach to Ganymede's synodic position in ``system`` (the SAME system
    the tube's own torus was built from -- Ganymede is that system's
    perturber, so its position is a pure function of absolute time via
    :func:`ganymede_synodic_position`).

    ``tube``'s states are sampled at ``theta1_section=0`` departure (the
    :func:`~cyclerfinder.search.ccr4bp_manifold_globalize.globalize_manifold_tube`
    convention), so absolute time is ``sign_dir * t_samples[j]`` (forward for
    ``branch="unstable"``, backward for ``branch="stable"``) -- matches
    :class:`~cyclerfinder.search.ccr4bp_manifold_globalize.ManifoldTube`'s own
    documented convention exactly.
    """
    sign_dir = 1.0 if tube.branch == "unstable" else -1.0
    dmin = math.inf
    theta2_min = 0.0
    t_min = 0.0
    r_min = math.inf
    n_valid = 0
    for i, theta2 in enumerate(tube.theta2_grid):
        if not tube.valid[i]:
            continue
        n_valid += 1
        for j, t in enumerate(tube.t_samples):
            st = tube.states[i, j]
            if np.any(np.isnan(st)):
                continue
            t_abs = sign_dir * float(t)
            gx, gy = ganymede_synodic_position(t_abs, system)
            d = math.hypot(float(st[0]) - gx, float(st[1]) - gy)
            r = math.hypot(float(st[0]), float(st[1]))
            r_min = min(r_min, r)
            if d < dmin:
                dmin = d
                theta2_min = float(theta2)
                t_min = float(t)
    return GanymedeApproach(
        dmin=dmin, theta2=theta2_min, t_flow=t_min, r_min=r_min, n_valid_departures=n_valid
    )


def to_physical_km_days(
    approach: GanymedeApproach, l_km: float, base_moon_period_days: float
) -> tuple[float, float]:
    """Convert a :class:`GanymedeApproach` to physical ``(dmin_km, t_days)``.

    ``l_km`` is the system's OWN base-moon semi-major axis (its length unit);
    ``base_moon_period_days`` is that same base moon's orbital period --
    caller-supplied explicitly rather than hardcoded, so this works for
    BOTH the Callisto-based and Europa-based systems without inheriting
    `#694`'s own documented Europa-only ``_L_KM`` hardcode (see module
    docstring).
    """
    tu_days = base_moon_period_days / (2.0 * math.pi)
    return approach.dmin * l_km, approach.t_flow * tu_days
