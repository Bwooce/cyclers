"""#649: cheap Hill-radius/scaled-energy coordinate-transform pilot for #628's
generative seed model's cross-mu problem.

`#642` proved `#628`'s Earth-Moon-trained generative model's raw
``(state0, period)`` output, fed directly to ``correct_periodic`` at a
foreign mass ratio mu, collapses onto degenerate L4/L5 equilibria at BOTH
tested cross-mu targets (mu=0.001 and Sun-Earth mu~3.0e-6) -- a STRUCTURAL
failure (the model has no mu-conditioning), not a sampling artifact. `#645`'s
shortlist item 4 proposed testing a CHEAP, non-ML coordinate fix before
accepting that cross-mu transfer is permanently out of reach: instead of
feeding the model's raw Earth-Moon-shaped output straight to the target-mu
corrector, first re-express it in a mu-INDEPENDENT coordinate, then invert
that coordinate at the target mu to build a target-mu-appropriate seed.

**The transform, precisely (the specific choices this task's own bullet asks
to be explicit about):**

1. **Energy**: `#629` established rho = (C-3)/(C_L1(mu)-3) as a genuinely
   useful mu-invariant "how deep in the ballistic-capture corridor" energy
   coordinate (`#629`'s own Titan mu-continuation work found sourced
   Ross-Roberts-Tsoukkas anchors sit at a near-invariant rho ~0.79-0.80
   across a 12x mu range). Compute the RAW decoded guess's own Jacobi
   constant C_guess = jacobi_constant(state0_guess, mu_train) (NOT the
   model's separately-generated ``jacobi`` feature column, which is not
   guaranteed to be numerically consistent with the decoded state0 -- the
   guess's OWN dynamically-computed C is the honest input), then
   rho = (C_guess - 3) / (C_L1(mu_train) - 3). At the target mu, invert:
   C_target = 3 + rho * (C_L1(mu_target) - 3) -- literally #629's own
   formula, solved for C instead of rho.

2. **Length scale**: this module's own :func:`hill_radius` is the EXACT
   numeric distance from the secondary (x=1-mu) to L1 (a root of
   ``collinear_lpoints``), not the classical mu^(1/3) leading-order
   asymptotic approximation -- reusing the project's own already-validated
   L1 solver rather than a cruder closed form. Position is decomposed as an
   offset from the secondary (``state0[:3] - (1-mu_train, 0, 0)``), scaled by
   ``hill_radius(mu_target) / hill_radius(mu_train)``, then re-anchored to
   the target system's own secondary location ``(1-mu_target, 0, 0)``. This
   is the direct generalization of `#629`'s energy-scaling idea to length:
   "how many Hill radii from the secondary" is treated as the mu-invariant
   quantity, exactly parallel to "how far through the C_L1 energy corridor".

3. **Velocity magnitude + direction**: raw velocity is scaled by the SAME
   Hill-radius ratio (Hill's equations have O(1) natural frequencies in
   units of the mean motion regardless of mu to leading order, so a velocity
   already expressed in Hill-normalized units transfers via the same length
   ratio -- no separate velocity scale). This Hill-scaled velocity vector's
   DIRECTION is kept, but its MAGNITUDE is then corrected so the constructed
   state's actual Jacobi constant (computed at mu_target, at the scaled
   position) exactly equals C_target from step 1 -- solving
   ``v_needed^2 = 2*Omega(pos_target; mu_target) - C_target`` for the speed
   along the preserved direction. This is necessary because position-scaling
   and energy-scaling are two independent operations on a nonlinear
   potential; simply scaling velocity by the Hill ratio like position does
   NOT in general reproduce the target rho exactly, and hitting the target
   rho *exactly* is the entire point of the transform (this task's own
   step 3: "construct a seed with the SAME rho ... by solving for the
   target-mu Jacobi constant that gives that rho"). When
   ``v_needed^2 <= 0`` (the Hill-scaled position sits outside the target
   system's own zero-velocity surface at C_target -- an under-determined /
   physically-unrealizable combination for that particular draw), the
   transform honestly returns ``None`` rather than fabricating a complex or
   negative-energy velocity; the caller counts this as a construction
   failure, exactly like any other seed the corrector rejects.

4. **Period**: left UNCHANGED. Hill's equations (the mu->0 limit CR3BP
   reduces to near a secondary, after Hill-normalizing length) are already
   nondimensionalized in units of 1/(mean motion) -- the SAME time unit CR3BP
   itself uses at any mu -- so a period already expressed in that common time
   unit does not, to leading (Hill) order, need an additional rescale the
   way length does. This is a documented, honestly-labeled APPROXIMATION
   (leading-order Hill scaling, not exact for mu gaps as large as
   0.01215 -> 3e-6), not a claim of exactness; the refinement corrector
   downstream is what absorbs any residual mismatch (or fails to, which is
   itself part of what this pilot measures).

None of steps 1-4 retrain or otherwise touch `#628`'s model -- this module
only reinterprets its raw output through a fixed, mu-parametrized coordinate
change before handing the result to the EXISTING
``cr3bp_periodic.correct_periodic`` corrector, exactly like `#624`'s original
(untransformed) protocol did.

See ``scripts/run_649_coordinate_fix_pilot.py`` for the evaluation that
reuses `#624`'s exact protocol (same N, same mu targets, same
``is_physically_sane``) to test whether this transform rescues any of the
cross-mu value `#642` found the raw approach lacks.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from cyclerfinder.core.cr3bp import jacobi_constant
from cyclerfinder.search.binary_star_search import collinear_lpoints


def jacobi_at_l1(mu: float) -> float:
    """Jacobi constant AT the L1 libration point for mass-ratio ``mu``.

    Same quantity `#629`'s own ``pluto_charon_kk_sweep._c_l1`` computes,
    reimplemented locally (not imported) so this module's dependency surface
    stays to the two generic, non-private helpers below
    (``collinear_lpoints``, ``jacobi_constant``) rather than reaching into
    another script's underscore-prefixed private helper.
    """
    l1, _l2, _l3 = collinear_lpoints(mu)
    return float(jacobi_constant(np.array([l1, 0.0, 0.0, 0.0, 0.0, 0.0]), mu))


def hill_radius(mu: float) -> float:
    """Exact numeric distance from the secondary (x=1-mu) to L1.

    This module's own mu-dependent length scale -- the direct, exact
    (not mu^(1/3)-asymptotic) analog of the classical Hill radius, reusing
    the project's own already-validated ``collinear_lpoints`` root-finder.
    """
    l1, _l2, _l3 = collinear_lpoints(mu)
    return float((1.0 - mu) - l1)


def rho_scaled_energy(jacobi: float, mu: float) -> float:
    """rho = (C - 3) / (C_L1(mu) - 3) -- `#629`'s scaled-energy coordinate."""
    return (jacobi - 3.0) / (jacobi_at_l1(mu) - 3.0)


@dataclass(frozen=True)
class TransformedSeed:
    """One coordinate-transformed seed, plus the intermediate quantities that
    produced it (kept for diagnostics/reporting, not just the final state)."""

    state0: NDArray[np.float64]
    period: float
    rho: float
    scale: float
    c_guess: float
    c_target: float


def transform_seed_to_target_mu(
    state0_guess: NDArray[np.float64] | list[float],
    period_guess: float,
    mu_train: float,
    mu_target: float,
) -> TransformedSeed | None:
    """Re-express a raw generated ``(state0_guess, period_guess)`` (decoded at
    ``mu_train``) as a seed appropriate for ``mu_target``, via the rho-matched
    Hill-radius-scaled transform documented in this module's own docstring.

    Returns ``None`` if the transform cannot construct a real-valued velocity
    at the target mu for this particular draw (see step 3 above) -- an
    honest construction failure, not a crash. The caller should treat this
    exactly like any other rejected seed (counted in the attempt total, not
    converged/sane).

    When ``mu_target == mu_train`` (bit-for-bit), returns the IDENTITY
    transform (unchanged state0/period) -- the self-consistency case this
    module's own test suite checks explicitly.
    """
    state0_guess = np.asarray(state0_guess, dtype=np.float64)
    if state0_guess.shape != (6,):
        raise ValueError(f"state0_guess must be a 6-vector, got shape {state0_guess.shape}")

    c_guess = jacobi_constant(state0_guess, mu_train)
    c_l1_train = jacobi_at_l1(mu_train)
    rho = (c_guess - 3.0) / (c_l1_train - 3.0)

    if mu_target == mu_train:
        return TransformedSeed(
            state0=state0_guess.copy(),
            period=period_guess,
            rho=rho,
            scale=1.0,
            c_guess=c_guess,
            c_target=c_guess,
        )

    c_l1_target = jacobi_at_l1(mu_target)
    c_target = 3.0 + rho * (c_l1_target - 3.0)

    r_h_train = hill_radius(mu_train)
    r_h_target = hill_radius(mu_target)
    if r_h_train <= 0.0 or r_h_target <= 0.0:
        raise ValueError(
            f"non-positive Hill radius (train={r_h_train}, target={r_h_target}) "
            f"-- mu_train={mu_train}, mu_target={mu_target}"
        )
    scale = r_h_target / r_h_train

    sec_train_x = 1.0 - mu_train
    sec_target_x = 1.0 - mu_target
    offset = state0_guess[:3] - np.array([sec_train_x, 0.0, 0.0])
    pos_target = np.array([sec_target_x, 0.0, 0.0]) + offset * scale
    vel_scaled = state0_guess[3:6] * scale

    x, y, z = pos_target
    r1 = math.sqrt((x + mu_target) ** 2 + y * y + z * z)
    r2 = math.sqrt((x - 1.0 + mu_target) ** 2 + y * y + z * z)
    omega = 0.5 * (x * x + y * y) + (1.0 - mu_target) / r1 + mu_target / r2
    v2_needed = 2.0 * omega - c_target
    if v2_needed <= 0.0:
        return None

    v_scaled_norm = float(np.linalg.norm(vel_scaled))
    if v_scaled_norm < 1e-12:
        return None

    v_dir = vel_scaled / v_scaled_norm
    v_mag = math.sqrt(v2_needed)
    vel_target = v_dir * v_mag

    state0_target = np.concatenate([pos_target, vel_target])
    return TransformedSeed(
        state0=state0_target,
        period=period_guess,
        rho=rho,
        scale=scale,
        c_guess=c_guess,
        c_target=c_target,
    )
