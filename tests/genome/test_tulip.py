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


# ---------------------------------------------------------------------------
# Phase 3 end-to-end gate: continuation + family-switching lands Np=2.
# ---------------------------------------------------------------------------


def test_find_tulip_via_continuation_lands_np2() -> None:
    """End-to-end Phase 3 reproduce gate: from the Koblick Np=1 seed,
    continuation + family-switching lands a Np=2 family member.

    Discipline:
      - The switched orbit must converge (residual below the corrector tol).
      - The petal count must equal 2 (independent topological cross-check).
      - The period must be within 10% of the Koblick paper Np=2 period
        (T = tau0 = 2.756 TU). Note: our family-switching corrector lands on
        the SAME-family Np=2 orbit but at the BIFURCATION POINT x0 (which is
        NOT the paper's fixed x0=1.0237 -- it's whatever x0 the multiplier
        crossed -1 at). So we compare PERIOD and JACOBI within precision,
        NOT raw IC.
      - The Jacobi must be in Koblick's documented Np=2 band ~3.05 +- 0.05.

    On failure: file as HONEST DIAGNOSTIC via xfail; do NOT tune to pass.
    """
    from cyclerfinder.genome.tulip import (
        KOBLICK_2023_TABLE4_PAPER,
        find_tulip_via_continuation,
    )

    result = find_tulip_via_continuation(np_target=2, d_x0=5e-4, n_steps_max=40)
    assert result.success, (
        f"find_tulip_via_continuation failed: reason={result.reason}; "
        f"branch_members={len(result.branch_members)}, "
        f"bifurcation={'yes' if result.bifurcation is not None else 'no'}"
    )
    switched = result.switched
    assert switched is not None
    # Topological check.
    s0 = np.array([switched.x0, 0, switched.z0, 0, switched.ydot0, 0])
    n = petal_count(s0, switched.T_TU, koblick_system())
    assert n == 2, f"switched orbit petal_count={n}, expected 2"
    # Period within +-10% of Koblick Np=2 tau0 (the FULL period, per the
    # Phase 3 finding that tau0 is the full period in Table 4).
    t_target = float(KOBLICK_2023_TABLE4_PAPER[2]["tau0"])
    period_ratio = switched.T_TU / t_target
    assert 0.90 < period_ratio < 1.10, (
        f"switched T={switched.T_TU:.4f} not within 10% of Koblick Np=2 "
        f"T={t_target:.4f} (ratio={period_ratio:.4f})"
    )
    # Jacobi in Koblick's Np=2 band (~3.05 +- 0.05).
    assert 3.00 < switched.jacobi < 3.10, (
        f"switched Jacobi {switched.jacobi:.4f} outside Koblick Np=2 band [3.00, 3.10]"
    )


# ---------------------------------------------------------------------------
# Phase 4 backward-compatibility gate: multi_shooting=True at Earth-Moon Np=2
# matches the Phase 3 single-shooting outcome to within the corrector tol.
# ---------------------------------------------------------------------------


def test_find_tulip_via_continuation_with_multi_shooting_matches_single_shooting_at_np2() -> None:
    """multi_shooting=True at the Earth-Moon k=2 bifurcation lands the same
    Np=2 family as single-shooting (period within 1%, Jacobi within 0.01).

    Phase 4's multi-shooter is an ESCALATION, not a replacement: at cases where
    single-shooting succeeds, both paths must land on the same family member.
    The gate is loose (1% on period, 0.01 on Jacobi) because the multi-shooter's
    free-variable layout differs from single-shooting's, so the two correctors
    can land on slightly different family-curve points -- but both must be in
    the same Np=2 family.
    """
    from cyclerfinder.genome.tulip import find_tulip_via_continuation

    # Reuse the standard Phase 3 settings to anchor the single-shooting path.
    ss = find_tulip_via_continuation(np_target=2, d_x0=5e-4, n_steps_max=40)
    ms = find_tulip_via_continuation(
        np_target=2,
        d_x0=5e-4,
        n_steps_max=40,
        multi_shooting=True,
    )
    assert ss.success, f"Phase 3 single-shooting did not succeed: {ss.reason}"
    assert ms.success, f"Phase 4 multi-shooting did not succeed: {ms.reason}"
    assert ss.switched is not None
    assert ms.switched is not None
    # Both should land in the same Np=2 family: period within 1%, Jacobi within
    # 0.01 of each other.
    period_diff = abs(ms.switched.T_TU - ss.switched.T_TU)
    assert period_diff / ss.switched.T_TU < 0.01, (
        f"multi-shooting period {ms.switched.T_TU:.4f} differs from single-shooting "
        f"{ss.switched.T_TU:.4f} by {period_diff:.3e} (> 1%)."
    )
    jacobi_diff = abs(ms.switched.jacobi - ss.switched.jacobi)
    assert jacobi_diff < 0.01, (
        f"multi-shooting Jacobi {ms.switched.jacobi:.4f} differs from single-shooting "
        f"{ss.switched.jacobi:.4f} by {jacobi_diff:.3e} (> 0.01)."
    )
    # Petal counts must both be 2.
    s0_ss = np.array([ss.switched.x0, 0, ss.switched.z0, 0, ss.switched.ydot0, 0])
    s0_ms = np.array([ms.switched.x0, 0, ms.switched.z0, 0, ms.switched.ydot0, 0])
    n_ss = petal_count(s0_ss, ss.switched.T_TU, koblick_system())
    n_ms = petal_count(s0_ms, ms.switched.T_TU, koblick_system())
    assert n_ss == 2 and n_ms == 2, (
        f"Phase 4 vs Phase 3 petal_count mismatch: ss={n_ss}, ms={n_ms} (both expected 2)"
    )


# ---------------------------------------------------------------------------
# #322 regression gates: 3D-topology gate rejects planar Np-petal collapse.
# ---------------------------------------------------------------------------


def test_322_marsphobos_paperseed_np4_does_not_claim_tulip() -> None:
    """Regression for #322: at Mars-Phobos (mu ~ 1.65e-8) the Koblick Np=4
    paper IC collapses to planar (z0 -> ~1e-14) when corrected, and
    petal_count returns 2 -- matching ``np_target=2`` and TRIGGERING the
    Tier A success gate UNDER THE OLD CODE (see ``data/scan_313_*.jsonl``).

    Post-fix expected behavior: ``find_tulip_at_system`` MUST NOT return
    ``success=True`` for this case. The result must carry
    ``topology_verdict == "planar Np-petal collapse"`` AND ``success=False``,
    since the orbit is not a genuine 3D tulip.
    """
    from cyclerfinder.genome.tulip import (
        KOBLICK_2023_TABLE4_PAPER,
        find_tulip_at_system,
    )

    system = cr3bp.cr3bp_system("Mars", "Phobos")
    result = find_tulip_at_system(
        system,
        np_target=2,
        seed_row=KOBLICK_2023_TABLE4_PAPER[4],
        try_direct_seed=True,
        multi_shooting=False,  # avoid the multishoot RuntimeError at very small mu
        n_steps_max=4,  # keep Tier B fallback brief; we only care about Tier A verdict
    )
    assert result is not None, "result should not be None when seed corrects"
    # The fix MUST refuse to claim success on the planar-collapse case.
    assert not result.success, (
        f"#322 regression: find_tulip_at_system claimed success on the "
        f"Mars-Phobos planar-collapse case (reason={result.reason}). "
        "The petal_count gate alone IS NOT SUFFICIENT to certify 3D tulip "
        "topology."
    )
    # And it must specifically diagnose the planar collapse.
    assert result.topology_verdict == "planar Np-petal collapse", (
        f"#322 regression: expected topology_verdict='planar Np-petal collapse', "
        f"got '{result.topology_verdict}' (reason={result.reason})."
    )
    # Sanity: max_abs_z should be small (below the threshold).
    assert result.max_abs_z is not None
    from cyclerfinder.genome.tulip import TULIP_Z_AMPLITUDE_FLOOR_NONDIM

    assert result.max_abs_z < TULIP_Z_AMPLITUDE_FLOOR_NONDIM, (
        f"planar-collapse case had max|z|={result.max_abs_z:.3e}, "
        f"expected < {TULIP_Z_AMPLITUDE_FLOOR_NONDIM:.3e}."
    )


def test_322_earth_moon_genuine_3d_tulip_still_admitted() -> None:
    """Regression for #322: the genuine 3D tulip at Earth-Moon (the existing
    Tier B end-to-end gate target) MUST continue to be admitted under the
    fixed code. Belt-and-braces: this test asserts ``success=True`` AND
    ``topology_verdict='3D tulip'`` AND ``max_abs_z >= floor`` for the
    Earth-Moon Np=2 case.
    """
    from cyclerfinder.genome.tulip import (
        TULIP_Z_AMPLITUDE_FLOOR_NONDIM,
        find_tulip_at_system,
    )

    system = koblick_system()
    # Force the Tier B path -- Tier A from the Koblick Np=1 seed at
    # Earth-Moon mu produces Np=1 (not Np=2), so the canonical EM workflow
    # exercises the Tier B continuation + family-switch.
    result = find_tulip_at_system(
        system,
        np_target=2,
        try_direct_seed=False,
        multi_shooting=False,
        d_x0=5e-4,
        n_steps_max=40,
    )
    assert result is not None
    assert result.success, (
        f"#322 regression: the canonical Earth-Moon Np=2 tulip MUST still "
        f"be admitted. reason={result.reason}, verdict={result.topology_verdict}."
    )
    assert result.topology_verdict == "3D tulip", (
        f"#322 regression: Earth-Moon Np=2 mis-classified as '{result.topology_verdict}'."
    )
    assert result.max_abs_z is not None
    assert result.max_abs_z >= TULIP_Z_AMPLITUDE_FLOOR_NONDIM, (
        f"Earth-Moon Np=2 max|z|={result.max_abs_z:.3e} below floor "
        f"{TULIP_Z_AMPLITUDE_FLOOR_NONDIM:.3e} -- the gate would reject a "
        "real tulip."
    )


def test_322_threshold_sanity_at_half_and_double() -> None:
    """Numeric sanity check for the chosen z-amplitude floor.

    Constructs a SYNTHETIC orbit that would pass the gate at z0 = 2*floor and
    fail at z0 = 0.5*floor. We use a CR3BP propagation from a perpendicular-
    crossing IC at chosen z0; with the symmetric corrector NOT involved (no
    z0 collapse), the orbit's z(t) excursion is bounded by an order-of-mag
    factor of z0 over short time scales. Both directions of the test are
    asserted:

      * z0 = 0.5 * floor -> ``is_three_dimensional`` returns False.
      * z0 = 2.0 * floor -> ``is_three_dimensional`` returns True.
    """
    from cyclerfinder.genome.tulip import (
        TULIP_Z_AMPLITUDE_FLOOR_NONDIM,
        is_three_dimensional,
    )

    sysm = koblick_system()
    # Use the Koblick Np=2 paper IC family at Earth-Moon, scaled in z. The IC
    # is at the perpendicular crossing (y0 = xdot0 = zdot0 = 0); we vary z0
    # and ydot0 stays the published value. (The orbit will NOT be closed at
    # these adjusted z0, but we only need ``out_of_plane_amplitude`` over one
    # period, not a closed orbit, for the gate-threshold test.)
    period = 2.756426
    floor = TULIP_Z_AMPLITUDE_FLOOR_NONDIM

    # Below-floor case: should return is_3d=False (short-circuits on z0 check).
    state_below = np.array([1.023731, 0.0, 0.5 * floor, 0.0, -0.082095, 0.0])
    is_3d_below, _ = is_three_dimensional(state_below, period, sysm)
    assert not is_3d_below, (
        f"is_three_dimensional returned True for z0={0.5 * floor:.3e} "
        f"(half the floor {floor:.3e}) -- threshold check broken."
    )

    # Above-floor case: must return is_3d=True. The Koblick Np=2 family has
    # max|z| comparable to z0, so 2*floor in z0 produces max|z| comfortably
    # above the floor.
    state_above = np.array([1.023731, 0.0, 2.0 * floor, 0.0, -0.082095, 0.0])
    is_3d_above, max_abs_z_above = is_three_dimensional(state_above, period, sysm)
    assert is_3d_above, (
        f"is_three_dimensional returned False for z0={2.0 * floor:.3e} "
        f"(twice the floor {floor:.3e}), max|z|={max_abs_z_above:.3e} -- "
        "threshold check overly strict."
    )
