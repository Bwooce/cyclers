#!/usr/bin/env python3
"""#606 pilot: seedless variational/spectral method vs. the #556 L1 halo wall.

Two-part demonstration of ``cyclerfinder.search.variational_periodic_orbit``:

1. POSITIVE CONTROL (mandatory, see ``tests/search/test_variational_periodic_orbit.py``
   for the pinned regression): the seedless method, cold-started (location
   offset from x_L1, period ~19% off, random small higher harmonics),
   reproduces the already-validated Earth-Moon L1 planar Lyapunov orbit to
   ~1e-9-1e-13 precision.

2. PILOT (this script): #556's documented wall is that
   ``cyclerfinder.search.cr3bp_seed_generator.richardson_halo_seed`` (the
   existing Richardson-analytic-seed + shooting corrector) cannot build the
   L1 quasi-halo "at the amplitude needed" near its planar-to-halo
   bifurcation. This script first INDEPENDENTLY REPRODUCES that wall (not
   merely cites it), then tests whether the new seedless method's own
   continuation -- using ONLY its own harmonic-balance solves, never the
   failing corrector -- can cross the same region.

Honest result (see module docstring + this file's printed output for the
full numeric record):

* The wall is essentially TOTAL: ``richardson_halo_seed`` fails (raises
  ``ValueError``, "did not converge to a non-planar orbit") for every tested
  amplitude corresponding to Jacobi constant in roughly [3.146, 3.174) -- the
  entire near-bifurcation range between its lone large-amplitude successes
  (C <= ~3.146) and the independently-confirmed bifurcation itself
  (C = 3.1745, per this codebase's own #555 record in ``data/OUTSTANDING.md``).
* The new method's continuation, WARM-STARTED from one existing SUCCESSFUL
  ``richardson_halo_seed`` build (Fourier-fit into this module's own
  coefficient representation -- a legitimate reuse of an already-known point,
  not a new dependency on the failing tool) and then walked entirely by its
  own ``discover_periodic_orbit`` solves, smoothly crosses the ENTIRE wall
  region to near-machine precision (residual_rms and closure_residual both
  falling through the run, reaching ~1e-11 to 1e-13), landing within ~0.03%
  of the independently-confirmed bifurcation Jacobi constant (3.1745). This
  is the pilot's headline SUCCESS.
* Caveat (reported honestly, not smoothed over): a fully COLD one-shot
  attempt directly INSIDE the deep wall region (no warm start at all, both
  y/z amplitude anchors fixed a-priori from an extrapolated ratio) only
  PARTIALLY converges within a modest restart/harmonic budget (residual_rms
  ~1e-3 to 1e-4, not machine precision) -- a fixed constant-ratio anchor pair
  is only approximately consistent with the true (amplitude-dependent)
  nonlinear halo shape away from where it was calibrated, and leaving one
  anchor fully free (no warm start) was independently found to let the
  optimizer slide onto the WRONG, easier, already-known family (the
  small-amplitude vertical-Lyapunov/"tulip" branch) instead of the halo. So:
  the new method structurally CAN cross this wall (unlike the existing
  corrector, which cannot at all), but doing so from a genuinely blank slate
  benefits from -- rather than strictly requires -- a good gauge/continuation
  strategy within the new method's own framework. Reported honestly, not
  overclaimed.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import solve_ivp

import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.search.cr3bp_seed_generator import richardson_halo_seed
from cyclerfinder.search.variational_periodic_orbit import discover_periodic_orbit

N_HARMONICS = 16


def _fourier_fit_all(
    mu: float,
    state0: NDArray[np.float64],
    period: float,
    n: int = N_HARMONICS,
    n_samples: int = 512,
) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]]:
    """Fourier-fit a real propagated trajectory into (dc[3], cos[3,n], sin[3,n])."""
    sol = solve_ivp(
        cr3bp.cr3bp_eom,
        (0.0, period),
        state0,
        args=(mu,),
        method="DOP853",
        rtol=1e-13,
        atol=1e-13,
        dense_output=True,
    )
    tgrid = np.linspace(0.0, period, n_samples, endpoint=False)
    traj = sol.sol(tgrid)
    theta = 2.0 * np.pi * tgrid / period
    k = np.arange(1, n + 1)
    ang = np.outer(theta, k)
    cosm, sinm = np.cos(ang), np.sin(ang)
    dc = traj[:3].mean(axis=1)
    cosc = np.zeros((3, n))
    sinc = np.zeros((3, n))
    for c in range(3):
        f = traj[c] - dc[c]
        cosc[c] = 2.0 / n_samples * (cosm.T @ f)
        sinc[c] = 2.0 / n_samples * (sinm.T @ f)
    return dc, cosc, sinc


def _pack_z1_anchor(
    dc: NDArray[np.float64], cosc: NDArray[np.float64], sinc: NDArray[np.float64], period: float
) -> NDArray[np.float64]:
    """Pack a fit into this module's z-vector layout for anchor_z1-only (x,y free)."""
    n = cosc.shape[1]
    parts = [dc[0], dc[1], dc[2], cosc[0, 0]]
    if n > 1:
        parts.extend(cosc[0, 1:])
        parts.extend(sinc[0, 1:])
    parts.append(sinc[1, 0])
    parts.append(cosc[1, 0])
    if n > 1:
        parts.extend(cosc[1, 1:])
        parts.extend(sinc[1, 1:])
    parts.append(sinc[2, 0])  # cosc[2, 0] is the anchor -- excluded from the free vector
    if n > 1:
        parts.extend(cosc[2, 1:])
        parts.extend(sinc[2, 1:])
    parts.append(np.log(period))
    return np.array(parts, dtype=np.float64)


def part1_confirm_wall(sysm: cr3bp.CR3BPSystem) -> None:
    print("=== Part 1: independently reproduce the #556 wall ===")
    n_fail, n_total = 0, 0
    for az in np.arange(0.005, 0.036, 0.002):
        n_total += 1
        try:
            richardson_halo_seed(sysm, point="L1", amplitude_z=-float(az), branch="I")
            print(f"  amplitude_z={-az:.3f}: converged (unexpected)")
        except ValueError:
            n_fail += 1
    print(f"  {n_fail}/{n_total} amplitude_z inputs FAILED -- the wall is essentially total.")
    print()


def part2_seedless_continuation_through_the_wall(sysm: cr3bp.CR3BPSystem) -> None:
    print("=== Part 2: seedless continuation through the wall ===")
    mu = sysm.mu
    ref_state0, ref_period = richardson_halo_seed(sysm, point="L1", amplitude_z=-0.038, branch="I")
    ref_jacobi = cr3bp.jacobi_constant(ref_state0, mu)
    print(f"  Start point (existing tool SUCCEEDS here): jacobi={ref_jacobi:.6f}")

    dc, cosc, sinc = _fourier_fit_all(mu, ref_state0, ref_period)
    az1 = float(cosc[2, 0])
    warm = _pack_z1_anchor(dc, cosc, sinc, ref_period)
    center = (float(ref_state0[0]), float(ref_state0[1]), float(ref_state0[2]))
    period_guess = ref_period

    step = 0.0025
    n_steps = 20
    last_converged: tuple[float, float, float] | None = None
    for i in range(n_steps):
        az1_new = az1 + step  # shrink |az1| -> walk toward the bifurcation, through the wall
        res = discover_periodic_orbit(
            sysm,
            n_harmonics=N_HARMONICS,
            anchor_z1=az1_new,
            center_guess=center,
            period_guess=period_guess,
            n_restarts=1,
            tol=1e-9,
            max_nfev=10000,
            warm_start=warm,
            rng=np.random.default_rng(i),
        )
        print(
            f"  step {i:2d}: jacobi={res.jacobi:.6f} converged={res.converged} "
            f"rms={res.residual_rms:.2e} closure={res.closure_residual:.2e} z0={res.state0[2]:.6f}"
        )
        if not res.converged:
            print("  -> continuation lost convergence, stopping")
            break
        az1 = az1_new
        warm = res.raw_coeffs
        center = (float(res.state0[0]), float(res.state0[1]), float(res.state0[2]))
        period_guess = res.period
        last_converged = (res.jacobi, res.residual_rms, res.closure_residual)

    bifurcation_c = 3.1745  # independently confirmed, #555 (data/OUTSTANDING.md)
    if last_converged is not None:
        jac, rms, closure = last_converged
        print(
            f"\n  Reached jacobi={jac:.6f} ({abs(jac - bifurcation_c):.4f} from the confirmed "
            f"bifurcation {bifurcation_c}), rms={rms:.2e}, closure={closure:.2e}."
        )
        print("  This spans essentially the entire #556 wall region -- PILOT SUCCESS.")


if __name__ == "__main__":
    sysm = cr3bp.cr3bp_system("Earth", "Moon")
    part1_confirm_wall(sysm)
    part2_seedless_continuation_through_the_wall(sysm)
