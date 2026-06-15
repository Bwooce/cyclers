"""Tulip-orbit genome reproduce-before-trust gates (#266 Phase 2).

Four mandatory gates per the orbit-closure discipline:

  1. **Sourced reproduce (Np=2)** -- the pumpkyn ``getTulip.m`` Np=2 butterfly
     IC closes within published-precision tolerance under our Sundman-
     regularised propagator, with the correct petal count and a Jacobi value
     in Koblick's documented band (~3.05 ± 0.05).

  2. **Independent integrator cross-check** -- the SAME IC closes under
     ``solve_ivp(method="Radau")`` (an implicit Runge-Kutta integrator,
     algorithmically independent from the DOP853 used by our propagators).
     Two independent integrators agreeing on the same IC verifies the IC, not
     just one integrator.

  3. **Petal classifier on Np=2** -- the topological classifier labels the
     Np=2 IC as 2 petals.

  4. **Petal classifier on Np=1 (parent NRHO)** -- after the NRHO seed is
     re-corrected to a true periodic orbit, the classifier labels it as 1
     petal.

The IC table itself only carries two rows (Np=1, Np=2) until the Koblick 2023
AMOSTECH paper is digitised. The discovery discipline forbids manufacturing
ICs; see the project memory entries `feedback_orbit_closure_discipline` and
`feedback_published_rounded_values_are_display`.
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp

import cyclerfinder.core.cr3bp as cr3bp
import cyclerfinder.search.cr3bp_periodic as cp
from cyclerfinder.genome.tulip import (
    KOBLICK_2023_TABLE4,
    koblick_system,
    petal_count,
    reproduce_tulip,
)

# ---------------------------------------------------------------------------
# Gate 1: sourced reproduce (Np = 2).
# ---------------------------------------------------------------------------


def test_koblick_np2_reproduces() -> None:
    """Pumpkyn Np=2 butterfly IC closes, has 2 petals, and lies in Koblick's
    Jacobi band.

    Closure tolerance is 1e-4 nondim per the task spec -- published precision
    is the limiter, not solver accuracy. Empirically the pumpkyn IC closes to
    < 1e-9 under DOP853 rtol/atol=1e-12, so the actual residual is far below
    the gate; the loose gate is the principled choice (relaxed to published
    precision; do NOT tighten without re-correcting the IC).
    """
    result = reproduce_tulip(np_target=2, n_periods=1)
    assert result.closed, (
        f"Np=2 closure residual {result.closure_residual:.3e} > 1e-4 -- "
        "pumpkyn IC failed the reproduce gate. This is a faithful negative; "
        "do NOT tune the tolerance to pass."
    )
    # Empirical floor: should be far below the loose 1e-4 gate.
    assert result.closure_residual < 1e-6, (
        f"observed residual {result.closure_residual:.3e} larger than the "
        "regression floor (1e-6) -- the propagator may have regressed."
    )
    assert result.n_petals_observed == 2, (
        f"observed {result.n_petals_observed} petals on Np=2 butterfly seed -- "
        "petal classifier broken or IC family-mismatched."
    )
    # Koblick (per the task description) reports the butterfly tulip Jacobi
    # band as approximately 3.05 +/- 0.05. The pumpkyn IC sits at ~3.058 -- in
    # the band.
    assert 3.00 < result.jacobi_constant < 3.10, (
        f"Jacobi {result.jacobi_constant:.4f} outside Koblick's [3.00, 3.10] "
        "band -- IC family-mismatched."
    )


# ---------------------------------------------------------------------------
# Gate 2: independent integrator cross-check.
# ---------------------------------------------------------------------------


def test_independent_cross_check_pumpkyn() -> None:
    """The pumpkyn Np=2 IC closes under solve_ivp(method='Radau').

    Radau is an implicit Runge-Kutta integrator independent of the DOP853
    (explicit Runge-Kutta with Dormand-Prince coefficients) used by the
    project's propagators -- agreement between the two on the same IC
    verifies the IC, not just one integrator's numerical fingerprint.

    See the project memory `feedback_orbit_closure_discipline`: independent
    cross-check is MANDATORY before trusting a sourced IC. The pumpkyn IC IS
    the cross-check anchor (sourced to Coorbital's public-domain MATLAB code,
    independent of any reading of the Koblick 2023 AMOSTECH paper).
    """
    row = KOBLICK_2023_TABLE4[2]
    state0 = np.array(
        [
            float(row["x0"]),  # type: ignore[arg-type]
            float(row["y0"]),  # type: ignore[arg-type]
            float(row["z0"]),  # type: ignore[arg-type]
            float(row["xdot0"]),  # type: ignore[arg-type]
            float(row["ydot0"]),  # type: ignore[arg-type]
            float(row["zdot0"]),  # type: ignore[arg-type]
        ]
    )
    period = float(row["T_TU"])  # type: ignore[arg-type]
    sysm = koblick_system()

    sol = solve_ivp(
        cr3bp.cr3bp_eom,
        (0.0, period),
        state0,
        args=(sysm.mu,),
        method="Radau",
        rtol=1e-11,
        atol=1e-11,
    )
    assert sol.success, f"Radau integrator failed: {sol.message}"

    state_f = sol.y[:, -1]
    closure = float(np.linalg.norm(state_f - state0))
    assert closure < 1e-4, (
        f"Radau closure residual {closure:.3e} > 1e-4 -- pumpkyn IC failed "
        "the independent cross-check. This is the orbit-closure discipline's "
        "hard stop; do NOT relax to pass."
    )

    # Sanity floor: Radau at rtol=1e-11 should hit ~1e-9 or better on this IC.
    assert closure < 1e-6, f"Radau closure {closure:.3e} above the regression floor (1e-6)."

    # Jacobi conservation across one period under the independent integrator.
    c0 = cr3bp.jacobi_constant(state0, sysm.mu)
    cf = cr3bp.jacobi_constant(state_f, sysm.mu)
    assert abs(cf - c0) < 1e-8, (
        f"Jacobi drift {abs(cf - c0):.3e} under Radau -- integrator failing "
        "to conserve the integral, IC may be off-family."
    )


# ---------------------------------------------------------------------------
# Gate 3: petal classifier labels Np=2 as 2.
# ---------------------------------------------------------------------------


def test_petal_classifier_labels_np2_as_2() -> None:
    """Classifier returns N_petals = 2 for the pumpkyn Np=2 butterfly IC."""
    row = KOBLICK_2023_TABLE4[2]
    state0 = np.array(
        [
            float(row["x0"]),  # type: ignore[arg-type]
            float(row["y0"]),  # type: ignore[arg-type]
            float(row["z0"]),  # type: ignore[arg-type]
            float(row["xdot0"]),  # type: ignore[arg-type]
            float(row["ydot0"]),  # type: ignore[arg-type]
            float(row["zdot0"]),  # type: ignore[arg-type]
        ]
    )
    period = float(row["T_TU"])  # type: ignore[arg-type]
    sysm = koblick_system()
    n = petal_count(state0, period, sysm)
    assert n == 2, (
        f"petal_count returned {n} on the Np=2 butterfly IC -- expected 2. "
        "Topology classifier or IC mismatched."
    )


# ---------------------------------------------------------------------------
# Gate 4: petal classifier labels (corrected) NRHO as 1.
# ---------------------------------------------------------------------------


def test_petal_classifier_labels_np1_as_1() -> None:
    """Classifier returns N_petals = 1 for the parent NRHO after correction.

    The Np=1 seed in the table is the widely-published NASA Gateway 9:2 NRHO
    starting guess (x, z, ydot) = (1.0213, -0.1824, -0.1031), T ~ 1.5111. It
    is **not** a converged periodic orbit at full precision -- the same seed
    is used by tests/search/test_reachable_impulsive.py with the corrector to
    refine it. We do the same here: re-correct, then run the classifier on
    the converged orbit.
    """
    row = KOBLICK_2023_TABLE4[1]
    state0 = np.array(
        [
            float(row["x0"]),  # type: ignore[arg-type]
            float(row["y0"]),  # type: ignore[arg-type]
            float(row["z0"]),  # type: ignore[arg-type]
            float(row["xdot0"]),  # type: ignore[arg-type]
            float(row["ydot0"]),  # type: ignore[arg-type]
            float(row["zdot0"]),  # type: ignore[arg-type]
        ]
    )
    period_guess = float(row["T_TU"])  # type: ignore[arg-type]
    sysm = koblick_system()
    orbit = cp.correct_periodic(sysm, state0, period_guess)
    assert orbit.converged, (
        f"NRHO corrector failed (residual {orbit.closure_residual:.3e}) -- "
        "test gate cannot proceed."
    )
    n = petal_count(orbit.state0, orbit.period, sysm)
    assert n == 1, f"petal_count returned {n} on the corrected Np=1 NRHO -- expected 1."


# ---------------------------------------------------------------------------
# Smoke gate: reproduce_tulip respects np_periods >= 1 and produces real
# Floquet eigenvalues with non-trivial structure.
# ---------------------------------------------------------------------------


def test_reproduce_tulip_floquet_structure_on_np2() -> None:
    """The Np=2 butterfly's Floquet multipliers exhibit the reciprocal-pair
    Hamiltonian structure (every multiplier has a 1/lambda partner) and a near-
    unit-modulus profile.

    This is a structural sanity gate -- not a numeric reproduction of Koblick's
    nu=1.144 value (which would require the original Koblick IC, not the
    pumpkyn cross-check). The HONEST report: the pumpkyn Np=2 IC's planar
    multiplier pair is (~1.004, ~0.996) -- products to ~1 (reciprocal) and the
    Barden stability index nu = 0.5(lambda + 1/lambda) sits at ~1.000009, i.e.
    on the stability boundary. The Koblick paper reports nu=1.144 for the
    same orbit family; if our IC sits at a slightly different family member
    (very plausible at 15 digits of IC precision), nu will differ. We assert
    only the structural property (reciprocal pairs, near-unit modulus) here.
    """
    result = reproduce_tulip(np_target=2)
    eigs = result.monodromy_eigs
    assert eigs.shape == (6,)

    # Reciprocal-pair structure: for every eigenvalue, 1/lambda is also a
    # multiplier within numerical tolerance. Verified by sorting eigs by
    # magnitude and pairing top with bottom, etc.
    mags = np.abs(eigs)
    # Sort by magnitude ascending.
    order = np.argsort(mags)
    sorted_eigs = eigs[order]
    # Three reciprocal pairs (in the planar-orbit case the trivial pair is
    # included). |sorted[0]| * |sorted[5]| ~ 1.
    products = [
        float(abs(sorted_eigs[0]) * abs(sorted_eigs[5])),
        float(abs(sorted_eigs[1]) * abs(sorted_eigs[4])),
        float(abs(sorted_eigs[2]) * abs(sorted_eigs[3])),
    ]
    for p in products:
        assert abs(p - 1.0) < 1e-3, (
            f"Floquet pairs failed reciprocal product = 1 check: products={products}"
        )

    # Near-unit-modulus profile: every multiplier sits within 1% of |.|=1 for
    # this near-bifurcation orbit.
    for m in mags:
        assert 0.99 < m < 1.01, (
            f"Floquet multiplier modulus {m:.4f} outside [0.99, 1.01] -- "
            "Np=2 butterfly should be near-marginal, not strongly hyperbolic."
        )
