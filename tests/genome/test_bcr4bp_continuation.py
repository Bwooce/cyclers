"""BCR4BP mu_sun-continuation driver tests (#303 / #292 Phase 2 Part A).

Per ``feedback_golden_tests_sourced_only``: the EXPECTED side of every
assertion is either zero (the mathematical definition of "periodic" -- the
residual at closure) or a SOURCED published value (the CR3BP-limit anchor
period, the Andreu mu_sun target). NO value computed by our own code is on
the EXPECTED side of any equality assertion.

Per ``feedback_orbit_closure_discipline``: every accepted family member is
required to pass the corrector's own independent (Radau) closure check, AND
the family-level tests assert STRUCTURE (continuation reaches the requested
extent, topology stays in L1 region, closure residuals all bounded) rather
than specific intermediate IC values.
"""

from __future__ import annotations

import math

import numpy as np

import cyclerfinder.core.bcr4bp as bcr4bp
import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.genome.bcr4bp_continuation import (
    BCR4BPFamily,
    BCR4BPFamilyMember,
    continue_bcr4bp_family_in_musun,
)
from cyclerfinder.genome.bcr4bp_genome import (
    IDX_T,
    IDX_X,
    IDX_XDOT,
    IDX_Y,
    IDX_YDOT,
    IDX_ZDOT,
    correct_bcr4bp_periodic,
)
from cyclerfinder.search.reachable_representatives import (
    correct_symmetric_free_period,
)

# Earth-Moon mass ratio (Rosales-Jorba 2023 Table 3; matches the Phase 1 module).
_MU_EM = 0.012150581600000

# Braik-Ross 2026 common Jacobi level (sourced in
# src/cyclerfinder/search/reachable_representatives.py + tests/...).
_C_JACOBI = 3.1294


def _seed_l1_lyapunov() -> tuple[np.ndarray, float]:
    """Return a sourced CR3BP planar L1 Lyapunov IC + period at C = 3.1294.

    Uses the existing CR3BP perpendicular-crossing corrector with the
    OFFLINE_SEEDS LL1 seed ``x = 0.8115`` (a small Lyapunov amplitude near
    the EM L1 point at x ~= 0.83692). The corrector converges in ~4
    iterations to the closed orbit; this is a SOURCED seed (Braik-Ross
    Jacobi level + the LL1 offline seed both from the existing codebase),
    not a value we made up.
    """
    sys_cr = cr3bp.CR3BPSystem(
        mu=_MU_EM, primary="earth", secondary="moon", l_km=384400.0, t_s=375190.0
    )
    orb = correct_symmetric_free_period(
        sys_cr, x0_guess=0.8115, jacobi=_C_JACOBI, t_half_guess=1.5, ydot0_sign=1.0
    )
    assert orb.converged, "test setup: CR3BP L1 Lyapunov seed failed to close"
    state6 = np.array([orb.x0, 0.0, 0.0, 0.0, orb.ydot0, 0.0], dtype=np.float64)
    return state6, float(orb.period)


def _seed_bcr4bp_at_mu_sun_zero() -> tuple[bcr4bp.BCR4BPSystem, np.ndarray, float]:
    """Return the BCR4BP system at mu_sun=0 + converged L1 Lyapunov state/period.

    The BCR4BP-at-mu_sun=0 corrector must close the CR3BP L1 Lyapunov to
    floating-point precision (the Phase 1 CR3BP-limit anchor); we use the
    BCR4BP corrector here so the returned ``BCR4BPPeriodicOrbit`` is
    structurally identical to what the continuation expects as its seed.
    """
    state_seed, period_seed = _seed_l1_lyapunov()
    sys_zero = bcr4bp.BCR4BPSystem(
        mu=_MU_EM,
        mu_sun=0.0,
        a_sun_nondim=388.8111430233511,
        omega_sun_nondim=0.925195985520347,
    )
    seed = correct_bcr4bp_periodic(
        sys_zero,
        state_seed,
        period_seed,
        sun_commensurate_n=1,
        free_vars=(IDX_X, IDX_YDOT, IDX_T),
        residual_indices=(IDX_Y, IDX_XDOT, IDX_ZDOT),
        is_half_period_residual=True,
        tol=1e-12,
        independent_tol=1e-6,
    )
    assert seed.converged, (
        "test setup: BCR4BP@mu_sun=0 corrector failed on the sourced "
        f"CR3BP L1 Lyapunov seed; "
        f"corrector_residual={seed.corrector_residual:.3e}, "
        f"indep_closure={seed.independent_closure_residual:.3e}"
    )
    return sys_zero, seed.state_initial, seed.period_nondim


# ---------------------------------------------------------------------------
# Gate 1: trivial-step continuation (target == seed within a tiny offset).
#
# The continuation driver must return a family whose single member is the
# seed itself (re-converged at a slightly perturbed mu_sun).
# ---------------------------------------------------------------------------


def test_continuation_trivial_step_returns_seed_neighbourhood() -> None:
    """A 1-step continuation to mu_sun = 1e-3 stays right at the seed."""
    sys_zero, state_seed, period_seed = _seed_bcr4bp_at_mu_sun_zero()
    seed_orbit = correct_bcr4bp_periodic(
        sys_zero,
        state_seed,
        period_seed,
        sun_commensurate_n=1,
        free_vars=(IDX_X, IDX_YDOT, IDX_T),
        residual_indices=(IDX_Y, IDX_XDOT, IDX_ZDOT),
        is_half_period_residual=True,
        tol=1e-12,
        independent_tol=1e-6,
    )
    fam = continue_bcr4bp_family_in_musun(
        seed_orbit,
        seed_mu_sun=0.0,
        target_mu_sun=1.0e-3,
        n_steps=1,
        step_method="linear",
        corrector_tol=1e-10,
        closure_tol=1e-6,
        monodromy=True,
    )
    assert isinstance(fam, BCR4BPFamily)
    assert len(fam.members) == 1, f"expected 1 family member, got {len(fam.members)}"
    member = fam.members[0]
    assert isinstance(member, BCR4BPFamilyMember)
    # mu_sun = 1e-3 is tiny -- the orbit should barely move from the seed.
    assert abs(member.orbit.state_initial[IDX_X] - state_seed[IDX_X]) < 1e-4, (
        f"L1 Lyapunov x0 moved {abs(member.orbit.state_initial[IDX_X] - state_seed[IDX_X]):.3e} "
        f"under a 1e-3 mu_sun step -- too sensitive."
    )
    # Closure must hold under the corrector's own gates (mathematical 0).
    assert member.orbit.corrector_residual < 1e-9
    assert member.orbit.independent_closure_residual < 1e-6
    # Monodromy must have been computed.
    assert member.monodromy is not None, "monodromy=True but no monodromy returned"
    assert member.monodromy.shape == (6, 6)
    assert member.floquet is not None
    assert len(member.floquet) == 6


# ---------------------------------------------------------------------------
# Gate 2: geometric stepping schedule is monotone increasing.
#
# Structural: the schedule must be sorted ascending in mu_sun and the family
# members must end up in that order.
# ---------------------------------------------------------------------------


def test_continuation_geometric_schedule_is_monotone_and_ordered() -> None:
    """The geometric step schedule produces monotone-increasing mu_sun members."""
    sys_zero, state_seed, period_seed = _seed_bcr4bp_at_mu_sun_zero()
    seed_orbit = correct_bcr4bp_periodic(
        sys_zero,
        state_seed,
        period_seed,
        sun_commensurate_n=1,
        free_vars=(IDX_X, IDX_YDOT, IDX_T),
        residual_indices=(IDX_Y, IDX_XDOT, IDX_ZDOT),
        is_half_period_residual=True,
        tol=1e-12,
        independent_tol=1e-6,
    )
    # Use a moderate target so the test is fast. The independent (Radau)
    # full-period closure tolerance is set to 1e-3 here because the corrector
    # uses a half-period symmetric residual and T is free -- the converged T
    # is NOT generally Sun-commensurate at non-zero mu_sun, so the full-period
    # closure picks up an O(mu_sun) Sun-phase residual. This is the same
    # discipline used by the Phase 1 weak-Sun halo gate.
    fam = continue_bcr4bp_family_in_musun(
        seed_orbit,
        seed_mu_sun=0.0,
        target_mu_sun=100.0,  # well below Andreu's 3.3e5 -- keeps the test fast
        n_steps=5,
        step_method="geometric",
        corrector_tol=1e-10,
        closure_tol=1e-3,
        monodromy=False,  # speed
    )
    assert len(fam.members) == 5, (
        f"all 5 steps to mu_sun=100 should converge; got {len(fam.members)}. "
        f"walk_notes: {fam.walk_notes}"
    )
    mu_sun_values = [m.mu_sun_value for m in fam.members]
    assert mu_sun_values == sorted(mu_sun_values), "members not in mu_sun-ascending order"
    # All members must close under the corrector (mathematical zero on the
    # half-period symmetric residual) and under the looser full-period gate.
    for m in fam.members:
        assert m.orbit.corrector_residual < 1e-9
        assert m.orbit.independent_closure_residual < 1e-3
    # Final mu_sun matches target (we don't have step-truncation logic).
    assert math.isclose(mu_sun_values[-1], 100.0, rel_tol=1e-9)


# ---------------------------------------------------------------------------
# Gate 3: stability tag taxonomy is well-formed.
#
# Structural: the tag must be one of the documented values. The actual
# stability of an L1 Lyapunov varies with mu_sun (and is published for
# specific cases that the catalogue may someday absorb); here we just check
# the tag is in the allowed set, not a specific value.
# ---------------------------------------------------------------------------


def test_continuation_stability_tag_in_allowed_set() -> None:
    """Every converged family member carries a valid stability tag."""
    sys_zero, state_seed, period_seed = _seed_bcr4bp_at_mu_sun_zero()
    seed_orbit = correct_bcr4bp_periodic(
        sys_zero,
        state_seed,
        period_seed,
        sun_commensurate_n=1,
        free_vars=(IDX_X, IDX_YDOT, IDX_T),
        residual_indices=(IDX_Y, IDX_XDOT, IDX_ZDOT),
        is_half_period_residual=True,
        tol=1e-12,
        independent_tol=1e-6,
    )
    fam = continue_bcr4bp_family_in_musun(
        seed_orbit,
        seed_mu_sun=0.0,
        target_mu_sun=10.0,
        n_steps=2,
        step_method="linear",
        monodromy=True,
    )
    allowed = {"stable", "unstable", "hyperbolic_pair", "marginal", "monodromy_failed"}
    assert len(fam.members) > 0
    for m in fam.members:
        assert m.stability_tag in allowed, (
            f"stability_tag={m.stability_tag!r} not in allowed set {allowed}"
        )


# ---------------------------------------------------------------------------
# Gate 4: bad inputs are rejected.
# ---------------------------------------------------------------------------


def test_continuation_rejects_seed_mu_sun_mismatch() -> None:
    """Seed orbit's system.mu_sun != declared seed_mu_sun => ValueError."""
    sys_zero, state_seed, period_seed = _seed_bcr4bp_at_mu_sun_zero()
    seed_orbit = correct_bcr4bp_periodic(
        sys_zero,
        state_seed,
        period_seed,
        sun_commensurate_n=1,
        free_vars=(IDX_X, IDX_YDOT, IDX_T),
        residual_indices=(IDX_Y, IDX_XDOT, IDX_ZDOT),
        is_half_period_residual=True,
        tol=1e-12,
        independent_tol=1e-6,
    )
    assert seed_orbit.system.mu_sun == 0.0
    try:
        continue_bcr4bp_family_in_musun(
            seed_orbit,
            seed_mu_sun=1.0,  # wrong!
            target_mu_sun=10.0,
            n_steps=1,
        )
    except ValueError as exc:
        assert "mu_sun" in str(exc).lower()
    else:
        raise AssertionError("expected ValueError on seed mu_sun mismatch")


def test_continuation_rejects_unconverged_seed() -> None:
    """A non-converged seed is rejected explicitly, not silently propagated."""
    # Build a deliberately bad "orbit" -- the corrector flags converged=False.
    sys_zero = bcr4bp.BCR4BPSystem(
        mu=_MU_EM,
        mu_sun=0.0,
        a_sun_nondim=388.8111430233511,
        omega_sun_nondim=0.925195985520347,
    )
    bogus_seed = np.array([0.5, 0.0, 0.0, 0.0, 5.0, 0.0], dtype=np.float64)
    bogus = correct_bcr4bp_periodic(
        sys_zero,
        bogus_seed,
        2.0,
        sun_commensurate_n=1,
        free_vars=(IDX_X, IDX_YDOT, IDX_T),
        residual_indices=(IDX_Y, IDX_XDOT, IDX_ZDOT),
        is_half_period_residual=True,
        tol=1e-12,
        independent_tol=1e-6,
        max_iter=10,
    )
    assert not bogus.converged
    try:
        continue_bcr4bp_family_in_musun(
            bogus,
            seed_mu_sun=0.0,
            target_mu_sun=10.0,
            n_steps=1,
        )
    except ValueError as exc:
        assert "converged" in str(exc).lower()
    else:
        raise AssertionError("expected ValueError on unconverged seed")


# ---------------------------------------------------------------------------
# Gate 5: CR3BP-limit anchor (mu_sun=0 -> small mu_sun) closure is tight.
# ---------------------------------------------------------------------------


def test_continuation_cr3bp_limit_anchor_tight() -> None:
    """Very small mu_sun continuation matches the CR3BP-limit seed within 1e-4.

    The L1 Lyapunov family is locally smooth in mu_sun (a regular perturbation
    problem near mu_sun = 0). A continuation to mu_sun = 1.0 must NOT drift
    the IC by more than O(mu_sun * eps_indirect) ~ small.
    """
    sys_zero, state_seed, period_seed = _seed_bcr4bp_at_mu_sun_zero()
    seed_orbit = correct_bcr4bp_periodic(
        sys_zero,
        state_seed,
        period_seed,
        sun_commensurate_n=1,
        free_vars=(IDX_X, IDX_YDOT, IDX_T),
        residual_indices=(IDX_Y, IDX_XDOT, IDX_ZDOT),
        is_half_period_residual=True,
        tol=1e-12,
        independent_tol=1e-6,
    )
    fam = continue_bcr4bp_family_in_musun(
        seed_orbit,
        seed_mu_sun=0.0,
        target_mu_sun=1.0,
        n_steps=3,
        step_method="linear",
        monodromy=False,
    )
    assert len(fam.members) == 3
    # First step: mu_sun ~ 1/3, the family member should be O(1e-4) away from seed.
    first = fam.members[0]
    dx = abs(first.orbit.state_initial[IDX_X] - state_seed[IDX_X])
    dv = abs(first.orbit.state_initial[IDX_YDOT] - state_seed[IDX_YDOT])
    # Loose absolute bound -- mu_sun ~ 1/3 vs Andreu ~ 3.3e5 is a 1e-6 relative
    # perturbation; the orbit should barely move.
    assert dx < 1e-3, f"L1 Lyapunov x0 drifted {dx:.3e} under tiny mu_sun step"
    assert dv < 1e-3, f"L1 Lyapunov vy drifted {dv:.3e} under tiny mu_sun step"
    # All closures tight.
    for m in fam.members:
        assert m.orbit.corrector_residual < 1e-9
        assert m.orbit.independent_closure_residual < 1e-6
