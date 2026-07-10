#!/usr/bin/env python3
"""Search runner for Task #546: Uranian-system CR3BP-level torus-connection band screen.

Per data/OUTSTANDING.md #546 ("decoupled Uranian-system CR3BP-level band screen,
un-gated from #538/#544"): point the existing #522/#534/#536 linking-number
heteroclinic-connection screen at the Uranian moon system, the one system that
has ever produced a confirmed novel catalogue row (#312, Umbriel-Oberon
quasi_cycler). Genuine Jacobi-constant band (not a single point).

CRITICAL SCOPE NOTE (established this session, not assumed): unlike Earth-Moon
(#534, Howell 1984 / Koblick 2023 published NRHO seeds) and Jupiter-Europa
(#536, in-repo literature seeds), there is NO published halo/NRHO seed table
for any Uranian moon anywhere in this codebase or (as far as this session's
search found) the literature -- these would be the first CR3BP halo/quasi-halo
families ever computed for this system. A naive linear vertical-mode seed at
fixed x0=x_L collapses to the trivial degenerate (zero-amplitude) solution
(empirically confirmed in this session's scratch testing). Instead, this
script bootstraps genuine finite-amplitude PLANAR LYAPUNOV orbits directly
from the collinear libration point using the already-validated, already-used-
elsewhere FIXED-JACOBI corrector (`cr3bp_periodic.correct_symmetric_fixed_jacobi`,
Ross & Roberts-Tsoukkas 2025 AAS 25-621 Eqs 9-12) with x0 free (not fixed),
which converges cleanly to nonzero-amplitude family members. Their intrinsic
out-of-plane (vertical) Floquet pair supplies the Neimark-Sacker seed for
`genome.qp_tori.correct_qp_torus` -- the same corrector run_534/536 use, just
fed a differently-sourced (but equally legitimate, physically standard)
parent orbit. This sidesteps needing any literature NRHO table or hand-typed
Richardson-approximation coefficients (neither of which exist for these
systems), at the cost of building QUASI-LYAPUNOV (near-planar) tori rather
than 3D quasi-halo tori -- a genuine capability extension, not a shortcut.

POSITIVE-CONTROL FINDING (this session): re-ran scripts/run_534_torus_connection.py
(Earth-Moon L1<->L2, the qp_tori/qp_torus_heteroclinic method's own best
same-method precedent) to completion. ALL THREE crossing searches (first-crossing,
and both copies of the advanced ydot<0 search) returned zero sign changes /
all-zero linking numbers -- i.e. #534 does NOT itself produce a positive
connection when actually run to completion, matching its own undercommitted
"Phase 2" notes in OUTSTANDING.md (the transit-vs-non-transit branch problem
was flagged as an open sub-problem, never resolved). #536 (Jupiter-Europa) is
ALSO a documented zero-connection run. So this method family is 0-for-2 on
every prior application, with NO validated positive anywhere in its history --
this session could not manufacture one either (task #312 is a structurally
different object, a Lambert/patched-conic multi-arc quasi-cycler, not a
libration-point torus-heteroclinic connection, so it cannot serve this role).
Per this project's own "verify a gauntlet with a positive control before
trusting 0/N" discipline, a Uranian 0/N from this method is NOT a validated
clean negative -- it is uninterpretable until the method itself is shown to
find at least one genuine connection somewhere. This script's results are
reported on that basis: exploratory data, not a certified empty-region claim.
"""

from __future__ import annotations

import math
import pathlib
from dataclasses import dataclass
from typing import Any

import numpy as np
from numpy.typing import NDArray
from scipy.optimize import brentq

import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.data.method_capability import MethodCapability
from cyclerfinder.data.preflight import preflight_search
from cyclerfinder.genome.qp_tori import correct_qp_torus
from cyclerfinder.genome.qp_torus_heteroclinic import build_manifold_grids, scan_linking_number
from cyclerfinder.search.bifurcation_detector import floquet_multipliers, monodromy
from cyclerfinder.search.cr3bp_periodic import correct_symmetric_fixed_jacobi

_REGION_ID = "uranian-moon-l1-l2-torus-linking-number-2026-07-10"
_METHOD = MethodCapability(
    genome=(
        "Uranian-moon (Miranda/Ariel/Umbriel/Oberon) L1/L2 quasi-Lyapunov torus "
        "manifold crossing grids (8x8 long/lat), Jacobi-constant band sweep, "
        "planar-Lyapunov-seeded (no literature NRHO table exists for this system)"
    ),
    corrector=(
        "correct_symmetric_fixed_jacobi (Ross & Roberts-Tsoukkas 2025) bootstrap "
        "-> correct_qp_torus -> build_manifold_grids + scan_linking_number"
    ),
    capability_tags=frozenset(
        {"cr3bp", "qp-torus", "heteroclinic", "linking-number", "uranian", "planar-lyapunov-seed"}
    ),
    git_sha="working-tree",
)


def l1_l2_x(mu: float) -> tuple[float, float]:
    """Root-find the collinear L1/L2 x-positions (elementary; no literature
    coefficients, just the zero of the rotating-frame x-acceleration at
    y=z=vx=vy=vz=0, which is the definition of a collinear libration point).
    """

    def f(x: float) -> float:
        r1 = x + mu
        r2 = x - 1.0 + mu
        return x - (1.0 - mu) * r1 / abs(r1) ** 3 - mu * r2 / abs(r2) ** 3

    x_sec = 1.0 - mu
    x_l1 = brentq(f, x_sec - 0.5, x_sec - 1e-10)
    x_l2 = brentq(f, x_sec + 1e-10, x_sec + 0.5)
    return x_l1, x_l2


def in_plane_frequency(x_l: float, mu: float, h: float = 1e-6) -> float | None:
    """Numeric (finite-difference) in-plane oscillatory frequency at a
    collinear point -- avoids trusting a memorized closed-form coefficient.
    """

    def accel_xy(x: float, y: float) -> NDArray[np.float64]:
        state = np.array([x, y, 0.0, 0.0, 0.0, 0.0])
        return cr3bp.cr3bp_eom(0.0, state, mu)[[3, 4]]

    jac = np.zeros((2, 2))
    for j, (dx, dy) in enumerate(((h, 0.0), (0.0, h))):
        plus = accel_xy(x_l + dx, dy)
        minus = accel_xy(x_l - dx, -dy)
        jac[:, j] = (plus - minus) / (2 * h)
    a4 = np.zeros((4, 4))
    a4[0, 2] = 1.0
    a4[1, 3] = 1.0
    a4[2, 0:2] = jac[0, :]
    a4[2, 3] = 2.0
    a4[3, 0:2] = jac[1, :]
    a4[3, 2] = -2.0
    eigvals = np.linalg.eigvals(a4)
    omega_candidates = [abs(e.imag) for e in eigvals if abs(e.real) < 1e-6 and abs(e.imag) > 1e-6]
    return max(omega_candidates) if omega_candidates else None


@dataclass
class LyapunovSeed:
    label: str
    x0: float
    ydot0: float
    period: float
    jacobi: float
    converged: bool
    floquet_pair: tuple[complex, complex] | None
    k_guess: int


def bootstrap_lyapunov_seed(
    system: cr3bp.CR3BPSystem,
    x_l: float,
    x0_guess_offset: float,
    c_target: float,
    t_guess: float,
    *,
    branch: str,
) -> LyapunovSeed | None:
    """Bootstrap a finite-amplitude planar Lyapunov orbit at a target Jacobi
    constant near a collinear point, then check for a genuine off-real-axis
    (Neimark-Sacker candidate) Floquet pair near the unit circle.

    ``branch`` is ``"L1"`` or ``"L2"`` -- used as a topology guard. The
    free-x0 fixed-Jacobi corrector (``correct_symmetric_fixed_jacobi``) has no
    knowledge of "which side of the secondary" it should stay on; at larger
    amplitude it can wander across and re-converge onto the OTHER collinear
    point's own family (empirically observed in this session's exploratory
    testing: an L2 seed collapsed onto the L1 family, identical x0/T, at
    C=2.999524 for Umbriel). Reject any result whose x0 landed on the wrong
    side of the secondary (x = 1 - mu) as a non-convergence to the intended
    branch, not a valid seed.
    """
    res = correct_symmetric_fixed_jacobi(
        system, x0_guess=x_l + x0_guess_offset, jacobi=c_target, period_guess=t_guess, tol=1e-10
    )
    if not res.converged:
        return None
    x_sec = 1.0 - system.mu
    if branch == "L1" and res.x0 >= x_sec:
        return None
    if branch == "L2" and res.x0 <= x_sec:
        return None
    state0 = np.array([res.x0, 0.0, 0.0, 0.0, res.ydot0, 0.0])
    mono = monodromy(system, state0, res.period, rtol=1e-13, atol=1e-13)
    eigs = floquet_multipliers(mono)
    # NS candidate threshold: imag part must be well clear of floating-point
    # noise around the always-present (1,1) trivial pair (empirically, noise
    # sits at ~1e-14..1e-5; genuine NS pairs found in this session's testing
    # had |imag| ~ 0.7-0.99), so 1e-3 cleanly separates the two.
    complex_pairs = [
        (e, ec)
        for e in eigs
        for ec in eigs
        if abs(e.imag) > 1e-3 and abs(ec - np.conj(e)) < 1e-6 and abs(abs(e) - 1.0) < 0.05
    ]
    if not complex_pairs:
        return LyapunovSeed("", res.x0, res.ydot0, res.period, res.jacobi, True, None, 0)
    lam_a, lam_b = complex_pairs[0]
    arg = abs(math.atan2(lam_a.imag, lam_a.real))
    k_guess = max(3, min(20, round(2 * math.pi / arg))) if arg > 1e-6 else 4
    return LyapunovSeed(
        "", res.x0, res.ydot0, res.period, res.jacobi, True, (lam_a, lam_b), k_guess
    )


def screen_moon(moon: str, c_targets: list[float], n_grid: int = 8) -> dict[str, Any]:
    """Screen one Uranus-<moon> CR3BP system for L1<->L2 torus heteroclinic
    connections across the given Jacobi-constant band.
    """
    sysm = cr3bp.cr3bp_system("Uranus", moon)
    mu = sysm.mu
    x_l1, x_l2 = l1_l2_x(mu)
    c_l1 = cr3bp.jacobi_constant(np.array([x_l1, 0, 0, 0, 0, 0]), mu)
    c_l2 = cr3bp.jacobi_constant(np.array([x_l2, 0, 0, 0, 0, 0]), mu)
    omega_l1 = in_plane_frequency(x_l1, mu)
    omega_l2 = in_plane_frequency(x_l2, mu)
    print(
        f"\n=== Uranus-{moon}: mu={mu:.6e}  xL1={x_l1:.7f} (C={c_l1:.6f})  "
        f"xL2={x_l2:.7f} (C={c_l2:.6f}) ==="
    )
    results = []
    surface_x = 1.0 - mu
    for c_target in c_targets:
        if c_target >= min(c_l1, c_l2):
            print(f"  C={c_target:.6f}: skipped (>= min(C_L1,C_L2), no finite-amplitude family)")
            continue
        # Offset AWAY from the secondary on each L-point's own side: L1 sits
        # below x=1-mu so a negative offset moves further from the secondary
        # (deeper into L1's own family); L2 sits above x=1-mu so the offset
        # must be positive for the same reason. Using the same-signed offset
        # for both (this session's first attempt used -2e-3 for both) let the
        # L2 corrector wander across x=1-mu into L1's own family at larger
        # amplitude -- exactly the failure the branch guard above now also
        # catches defensively.
        seed_l1 = bootstrap_lyapunov_seed(
            sysm,
            x_l1,
            -2e-3,
            c_target,
            t_guess=2 * math.pi / omega_l1 if omega_l1 else 3.0,
            branch="L1",
        )
        seed_l2 = bootstrap_lyapunov_seed(
            sysm,
            x_l2,
            2e-3,
            c_target,
            t_guess=2 * math.pi / omega_l2 if omega_l2 else 3.0,
            branch="L2",
        )
        if seed_l1 is None or seed_l2 is None:
            print(f"  C={c_target:.6f}: Lyapunov bootstrap FAILED to converge (L1 or L2)")
            continue
        if seed_l1.floquet_pair is None or seed_l2.floquet_pair is None:
            print(
                f"  C={c_target:.6f}: no near-unit-circle complex Floquet pair (no NS torus seed)"
            )
            continue
        print(
            f"  C={c_target:.6f}: L1 x0={seed_l1.x0:.6f} T={seed_l1.period:.5f} "
            f"lam={seed_l1.floquet_pair[0]:.4f} k~{seed_l1.k_guess} | "
            f"L2 x0={seed_l2.x0:.6f} T={seed_l2.period:.5f} "
            f"lam={seed_l2.floquet_pair[0]:.4f} k~{seed_l2.k_guess}"
        )
        try:
            torus_l1 = correct_qp_torus(
                sysm,
                np.array([seed_l1.x0, 0.0, 0.0, 0.0, seed_l1.ydot0, 0.0]),
                seed_l1.period,
                seed_l1.floquet_pair,
                k=seed_l1.k_guess,
                n_trans=8,
                initial_torus_amplitude=5e-4,
                tol=1e-5,
            )
            torus_l2 = correct_qp_torus(
                sysm,
                np.array([seed_l2.x0, 0.0, 0.0, 0.0, seed_l2.ydot0, 0.0]),
                seed_l2.period,
                seed_l2.floquet_pair,
                k=seed_l2.k_guess,
                n_trans=8,
                initial_torus_amplitude=5e-4,
                tol=1e-5,
            )
        except Exception as exc:
            print(f"  C={c_target:.6f}: torus correction raised {exc!r}")
            continue
        print(
            f"    torus L1 residual={torus_l1.invariance_residual:.2e} "
            f"torus L2 residual={torus_l2.invariance_residual:.2e}"
        )
        stable_grid, unstable_grid = build_manifold_grids(
            torus_l2,
            torus_l1,
            n_long=n_grid,
            n_lat=n_grid,
            eps=1e-5,
            surface_x=surface_x,
            t_max=15.0,
            stable_sign=-1.0,
            unstable_sign=1.0,
        )
        n_cross_u = int(np.sum(np.isfinite(unstable_grid.endpoints[:, :, 0])))
        n_cross_s = int(np.sum(np.isfinite(stable_grid.endpoints[:, :, 0])))
        n_pts = n_grid * n_grid
        print(f"    crossings: unstable={n_cross_u}/{n_pts} stable={n_cross_s}/{n_pts}")
        n_sign_changes = 0
        z_s = stable_grid.endpoints[:, :, 2]
        z_u = unstable_grid.endpoints[:, :, 2]
        z_s_finite = z_s[np.isfinite(z_s)]
        z_u_finite = z_u[np.isfinite(z_u)]
        if len(z_s_finite) > 0 and len(z_u_finite) > 0:
            overlap_min = max(z_s_finite.min(), z_u_finite.min())
            overlap_max = min(z_s_finite.max(), z_u_finite.max())
            if overlap_min < overlap_max:
                d_values = np.linspace(overlap_min, overlap_max, 30)
                link_result = scan_linking_number(
                    stable_grid,
                    unstable_grid,
                    scanning_component="z",
                    curve_components=("y", "ydot", "zdot"),
                    d_values=d_values,
                )
                sign_changes = link_result.sign_change_locations()
                n_sign_changes = len(sign_changes)
                print(
                    f"    z overlap [{overlap_min:.5f},{overlap_max:.5f}] "
                    f"sign_changes={sign_changes}"
                )
            else:
                print("    no z-overlap between stable/unstable crossing grids")
        else:
            print("    no finite crossings on one or both grids")
        results.append(
            {
                "moon": moon,
                "C_target": c_target,
                "n_cross_u": n_cross_u,
                "n_cross_s": n_cross_s,
                "n_grid_pts": n_pts,
                "n_sign_changes": n_sign_changes,
            }
        )
    return {"moon": moon, "mu": mu, "C_L1": c_l1, "C_L2": c_l2, "points": results}


def main() -> None:
    preflight_search(
        task_no=546,
        region_id=_REGION_ID,
        method=_METHOD,
        script_path=pathlib.Path(__file__),
        n_points=8 * 8 * 2 * 3 * 3,  # grid_pts * (stable+unstable) * ~3 C-values * ~3 moons
        override_reason=(
            "exploratory first-pass band screen with a brand-new (this-session) "
            "planar-Lyapunov torus bootstrap for a system with no prior timing "
            "history; grid deliberately kept small (8x8, not 16x16) and the "
            "Jacobi band deliberately kept short (a handful of points, not "
            "dozens) specifically to bound wall-clock risk while the bootstrap's "
            "own convergence rate at each (moon, C) point is still unknown; a "
            "timing pilot would cost the same wall-clock as just running the "
            "bounded sweep itself. See #546 OUTSTANDING entry for the full "
            "positive-control caveat this run is reported under."
        ),
    )

    all_results = []
    # Umbriel first: directly adjacent to #312 (Umbriel-Oberon quasi_cycler).
    for moon in ["Umbriel", "Oberon", "Ariel"]:
        sysm = cr3bp.cr3bp_system("Uranus", moon)
        mu = sysm.mu
        x_l1, x_l2 = l1_l2_x(mu)
        c_l1 = cr3bp.jacobi_constant(np.array([x_l1, 0, 0, 0, 0, 0]), mu)
        c_l2 = cr3bp.jacobi_constant(np.array([x_l2, 0, 0, 0, 0, 0]), mu)
        c_ref = min(c_l1, c_l2)
        c_targets = [c_ref - frac * abs(c_ref) for frac in (2e-4, 5e-4, 8e-4, 1.1e-3, 1.5e-3, 2e-3)]
        result = screen_moon(moon, c_targets, n_grid=8)
        all_results.append(result)

    print("\n\n===== SUMMARY =====")
    for r in all_results:
        print(f"{r['moon']}: mu={r['mu']:.4e} C_L1={r['C_L1']:.6f} C_L2={r['C_L2']:.6f}")
        for p in r["points"]:
            print(
                f"  C={p['C_target']:.6f}: cross u={p['n_cross_u']}/{p['n_grid_pts']} "
                f"s={p['n_cross_s']}/{p['n_grid_pts']} sign_changes={p['n_sign_changes']}"
            )


if __name__ == "__main__":
    main()
