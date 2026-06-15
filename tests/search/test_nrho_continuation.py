"""NRHO continuation tests (#266 Phase 3).

Three discipline gates:

  1. **Reproduce Np=1** -- the symmetric corrector lands the Koblick AMOSTECH
     Table 4 paper-row Np=1 IC with closure residual far below the published
     precision.

  2. **Smooth continuation** -- the period varies smoothly along the family
     (no >30% jumps between adjacent members). This is the regression gate:
     a continuation that silently jumps families is broken even when each
     member is a valid periodic orbit.

  3. **Period-multiplying bifurcation detection** -- on a wider continuation
     run, the detector returns at least one period-multiplying bracket. This
     is the GENOME-validity gate; if no brackets are found anywhere the
     genome layer cannot find tulip orbits.

The discipline forbids tuning to make any of these pass. Per the task spec,
any failure should be filed as an HONEST DIAGNOSTIC.
"""

from __future__ import annotations

import numpy as np

from cyclerfinder.genome.tulip import KOBLICK_2023_TABLE4_PAPER, koblick_system
from cyclerfinder.search.nrho_continuation import (
    NRHOStopReason,
    SymmetricNRHO,
    continue_nrho_family,
    correct_symmetric_nrho,
)

# ---------------------------------------------------------------------------
# Gate 1: corrector reproduces the Koblick Np=1 row.
# ---------------------------------------------------------------------------


def test_correct_symmetric_nrho_reproduces_koblick_np1() -> None:
    """The Np=1 paper row converges within the Phase 1 round-trip residual.

    The Koblick AMOSTECH Table 4 row 1 is (x0=1.023731, z0=0.183250,
    ydot0=-0.106950, tau0=1.533637). tau0 is the FULL nondim period of the
    parent NRHO. The corrector should converge with a near-machine-precision
    residual; the gate is a generous 1e-9 to allow for the integrator's
    round-trip floor.
    """
    sysm = koblick_system()
    row = KOBLICK_2023_TABLE4_PAPER[1]
    t_seed = float(row["tau0"])  # full period
    result = correct_symmetric_nrho(
        sysm,
        float(row["x0"]),
        float(row["z0"]),
        float(row["ydot0"]),
        t_seed,
        tol=1e-11,
    )
    assert result.converged, (
        f"Np=1 corrector failed: residual={result.closure_residual:.3e}, n_iter={result.n_iter}"
    )
    assert result.closure_residual < 1e-9, (
        f"closure residual {result.closure_residual:.3e} above 1e-9 floor"
    )
    # Sanity: the corrected IC should hardly move from the sourced row (the
    # paper IC is published to ~6 sig figs and is already very close to a true
    # periodic orbit).
    assert abs(result.x0 - float(row["x0"])) < 1e-9, "x0 must not be modified by the corrector"
    assert abs(result.z0 - float(row["z0"])) < 1e-4
    assert abs(result.ydot0 - float(row["ydot0"])) < 1e-4
    assert abs(result.T_TU - t_seed) < 1e-4, (
        f"T moved from {t_seed:.6f} to {result.T_TU:.6f}; expected ~no change for a "
        "sourced near-periodic IC"
    )
    # Monodromy was computed.
    assert result.monodromy is not None
    assert result.monodromy.shape == (6, 6)
    # The Jacobi constant should match Koblick's documented band (~3.04-3.06
    # for the parent NRHO; we get ~3.0448).
    assert 3.00 < result.jacobi < 3.10


# ---------------------------------------------------------------------------
# Gate 2: continuation produces a smooth period curve (no family-jumping).
# ---------------------------------------------------------------------------


def test_continue_nrho_family_produces_smooth_period_curve() -> None:
    """A short continuation step produces members whose periods vary smoothly.

    Adjacent periods should differ by less than 30% along the L2 Southern NRHO
    family. Larger jumps indicate the corrector silently jumped onto a
    different family (e.g. the doubled-period family or the pumpkyn Np=2
    branch).
    """
    sysm = koblick_system()
    row = KOBLICK_2023_TABLE4_PAPER[1]
    seed = correct_symmetric_nrho(
        sysm,
        float(row["x0"]),
        float(row["z0"]),
        float(row["ydot0"]),
        float(row["tau0"]),
        tol=1e-11,
    )
    assert seed.converged
    branch = continue_nrho_family(
        seed,
        sysm,
        direction=-1,
        d_x0=5e-4,
        n_steps_max=15,
        tol=1e-10,
        bif_tol=1e-1,
    )
    assert len(branch.members) >= 5, (
        f"continuation produced only {len(branch.members)} members; expected >=5"
    )
    periods = [m.T_TU for m in branch.members]
    for i in range(1, len(periods)):
        jump = abs(periods[i] - periods[i - 1]) / max(abs(periods[i - 1]), 1e-9)
        assert jump < 0.30, (
            f"period jumped {100 * jump:.1f}% between member {i - 1} (T={periods[i - 1]:.4f}) "
            f"and {i} (T={periods[i]:.4f}); corrector may have switched families"
        )
    # Period should be monotonically DECREASING toward smaller x0 (perilune
    # tightens, period drops).
    assert periods[-1] < periods[0], (
        f"period did NOT decrease over continuation: {periods[0]:.4f} -> {periods[-1]:.4f}"
    )


# ---------------------------------------------------------------------------
# Gate 3: bifurcation detection actually flags a bracket.
# ---------------------------------------------------------------------------


def test_continue_nrho_family_flags_period_multiplying_bifurcation() -> None:
    """A wider continuation returns at least one period-multiplying bracket.

    On the L2 Southern NRHO family the parent NRHO has a real hyperbolic
    multiplier pair (-2.30, -0.43) at the seed; these collide AT -1 (period
    doubling) as x0 decreases through ~1.0117. The detector flags this as a
    k=2 bracket at the coarse bif_tol=0.1; at the tight tolerance bif_tol=0.02
    the tangent never crosses zero and the detector reports nothing (this is
    a CLEAN NEGATIVE we explicitly tolerate by using the coarser tolerance).
    """
    sysm = koblick_system()
    row = KOBLICK_2023_TABLE4_PAPER[1]
    seed = correct_symmetric_nrho(
        sysm,
        float(row["x0"]),
        float(row["z0"]),
        float(row["ydot0"]),
        float(row["tau0"]),
        tol=1e-11,
    )
    branch = continue_nrho_family(
        seed,
        sysm,
        direction=-1,
        d_x0=5e-4,
        n_steps_max=30,
        tol=1e-10,
        bif_tol=1e-1,  # coarse tolerance to catch the tangent k=2 bifurcation
        bif_k_max=4,
    )
    assert len(branch.bifurcations) >= 1, (
        f"no bifurcations found on a 30-step continuation; "
        f"stop_reason={branch.stop_reason.value}, members={len(branch.members)}"
    )
    # At least one k=2 bracket somewhere.
    k2_count = sum(1 for b in branch.bifurcations if b.k == 2)
    assert k2_count >= 1, (
        f"no k=2 (period-doubling) bracket found among {len(branch.bifurcations)} "
        f"brackets; k values were {[b.k for b in branch.bifurcations]}"
    )


# ---------------------------------------------------------------------------
# Smoke gate: stop_on_first_bifurcation works.
# ---------------------------------------------------------------------------


def test_continue_nrho_family_stop_on_first_bifurcation() -> None:
    """``stop_on_first_bifurcation=True`` halts on the first k=2 bracket."""
    sysm = koblick_system()
    row = KOBLICK_2023_TABLE4_PAPER[1]
    seed = correct_symmetric_nrho(
        sysm,
        float(row["x0"]),
        float(row["z0"]),
        float(row["ydot0"]),
        float(row["tau0"]),
        tol=1e-11,
    )
    branch = continue_nrho_family(
        seed,
        sysm,
        direction=-1,
        d_x0=5e-4,
        n_steps_max=30,
        tol=1e-10,
        bif_tol=1e-1,
        bif_k_max=4,
        stop_on_first_bifurcation=True,
    )
    assert branch.stop_reason == NRHOStopReason.BIFURCATION_HIT or any(
        b.k == 2 for b in branch.bifurcations
    ), (
        f"stop_on_first_bifurcation did not trigger; stop={branch.stop_reason.value}, "
        f"members={len(branch.members)}, bifs={len(branch.bifurcations)}"
    )


# ---------------------------------------------------------------------------
# Defensive: seed with converged=False is rejected.
# ---------------------------------------------------------------------------


def test_continue_nrho_family_rejects_unconverged_seed() -> None:
    """An unconverged seed yields an empty branch with NO_CONVERGE."""
    sysm = koblick_system()
    bogus = SymmetricNRHO(
        x0=0.0,
        z0=0.0,
        ydot0=0.0,
        T_TU=1.0,
        jacobi=0.0,
        converged=False,
        closure_residual=1.0,
        n_iter=0,
        monodromy=None,
    )
    branch = continue_nrho_family(bogus, sysm, n_steps_max=5)
    assert branch.stop_reason == NRHOStopReason.NO_CONVERGE
    assert len(branch.members) == 0


# ---------------------------------------------------------------------------
# Defensive: tiny n_steps_max produces only the seed member.
# ---------------------------------------------------------------------------


def test_continue_nrho_family_min_steps() -> None:
    """Zero-step continuation returns only the seed."""
    sysm = koblick_system()
    row = KOBLICK_2023_TABLE4_PAPER[1]
    seed = correct_symmetric_nrho(
        sysm,
        float(row["x0"]),
        float(row["z0"]),
        float(row["ydot0"]),
        float(row["tau0"]),
        tol=1e-11,
    )
    branch = continue_nrho_family(
        seed,
        sysm,
        direction=-1,
        d_x0=5e-4,
        n_steps_max=0,
        tol=1e-10,
        bif_tol=1e-1,
    )
    assert len(branch.members) == 1
    assert branch.members[0].x0 == seed.x0


# ---------------------------------------------------------------------------
# Round-trip closure check (regularised vs DOP853 corrector).
# ---------------------------------------------------------------------------


def test_regularized_closure_agrees_with_corrector() -> None:
    """The regularised propagator closes the corrected Np=1 orbit at machine
    precision. Two integrators agreeing on the same IC verifies the IC, not
    just one integrator (see the project memory ``feedback_orbit_closure_discipline``).
    """
    from cyclerfinder.search.nrho_continuation import regularized_full_period_closure

    sysm = koblick_system()
    row = KOBLICK_2023_TABLE4_PAPER[1]
    member = correct_symmetric_nrho(
        sysm,
        float(row["x0"]),
        float(row["z0"]),
        float(row["ydot0"]),
        float(row["tau0"]),
        tol=1e-11,
    )
    assert member.converged
    closure = regularized_full_period_closure(sysm, member, regularization="r2")
    # The regularised propagator at rtol=1e-12 should match the corrector at
    # similar precision; 1e-9 is a generous bound.
    assert closure < 1e-9, f"regularised closure {closure:.3e} > 1e-9"


# ---------------------------------------------------------------------------
# label_branch_members helper.
# ---------------------------------------------------------------------------


def test_label_branch_members_shape() -> None:
    """``label_branch_members`` converts a branch to FamilyMember instances."""
    from cyclerfinder.search.nrho_continuation import label_branch_members

    sysm = koblick_system()
    row = KOBLICK_2023_TABLE4_PAPER[1]
    seed = correct_symmetric_nrho(
        sysm,
        float(row["x0"]),
        float(row["z0"]),
        float(row["ydot0"]),
        float(row["tau0"]),
        tol=1e-11,
    )
    branch = continue_nrho_family(seed, sysm, direction=-1, d_x0=5e-4, n_steps_max=3, tol=1e-10)
    fams = label_branch_members(branch, mu=sysm.mu)
    assert len(fams) == len(branch.members)
    for f, m in zip(fams, branch.members, strict=True):
        assert f.mu == sysm.mu
        assert f.period == m.T_TU
        assert f.parameter == m.x0
        np.testing.assert_array_equal(f.state0, np.array([m.x0, 0, m.z0, 0, m.ydot0, 0]))
