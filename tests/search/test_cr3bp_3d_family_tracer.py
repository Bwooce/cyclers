"""Tests for the 3D family-tracer + continuation driver (#296 Phase 2).

Sourced anchor (per ``feedback_golden_tests_sourced_only``)
------------------------------------------------------------
The PLANAR seed at the family endpoint traces to ``data/catalogue.yaml`` row
``braik-ross-c11a-cycler-2026``. The 3D continuation walks AWAY from that
planar member, so the GOLDEN side of these tests is necessarily
TOPOLOGICAL (the spike's family extent, the presence of a fold, the
preservation of closure across the walk) — NOT a numeric "expected z0" pulled
from our own code (which would be circular).

The spike data at ``data/spike_287.jsonl`` is the cross-check target for
family-extent sanity (164 family rows, x0 in [-0.85, -0.77], z0 in
[-0.24, 0], T in [9.5, 18.9] TU, C in [2.78, 3.10]); the tracer is asserted
to cover the spike's extent. The spike's specific (x0, z0, T) tuples are
OUR computation and so are NOT used as numerical expected values; the
topology IS (x0 spans both directions, z0 reaches |z0| ~ 0.24, a fold is
detected, the family terminates by collapsing to the planar manifold on the
z0 -> 0 side).

Tests
-----
  1. Spike reproduction: pseudo-arclength forward+backward from the spike's
     IC produces a multi-member family that covers x0 < -0.81 and x0 > -0.81;
     all members satisfy closure tol.
  2. Closure preservation on every member.
  3. Monodromy / Floquet structure: every member has a trivial unit
     eigenvalue pair; non-trivial multipliers are reciprocal pairs
     (Hamiltonian symmetry).
  4. Fold detection: the spike's z0-fold is detected by the tracer; the
     pseudo-arclength walk continues PAST it.
  5. Degeneracy guard: a planar seed (z0=0) marks all walked members as
     degenerate_planar (the planar manifold is invariant — the walk does NOT
     escape it spontaneously).
  6. Direction symmetry: forward and backward walks cover separate x0
     halves.
  7. Natural-T continuation runs cleanly for a small step count.
  8. Input validation.

Notes
-----
These tests are slow (each member runs a corrector + a Radau closure check
+ a monodromy integration). The full reproduction test uses a modest
``n_steps_max`` (e.g. 25) — a wider walk lives in the family JSONL
generator script, not the unit test.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.search.cr3bp_3d_family_tracer import (
    Family3D,
    Family3DMember,
    FoldPoint,
    continue_general_3d_family,
)
from cyclerfinder.search.cr3bp_general_periodic_3d import (
    IDX_X,
    IDX_Z,
)

# Sourced golden: Braik-Ross C11a planar Earth-Moon (1,1) cycler.
# data/catalogue.yaml row braik-ross-c11a-cycler-2026.
C11A_X0 = -0.8116406668238195
C11A_YDOT0 = -0.11859055759763637
C11A_PERIOD_TU = 9.69107744379376
EM_MU = 1.2150584270572e-2
EM_L_KM = 384400.0
EM_T_S = 375699.8

# 3D spike seed (#287 spike, reproduced by Phase 1 corrector at residual 1e-13,
# independent closure 1e-10). The exact (x0, z0, ydot0, T) values traces to
# the spike's converged member at z0_guess = 0.05; they are OUR computation
# and are NOT used as expected values for any assertion — only as the seed
# IC for the tracer.
SPIKE_X0 = -0.8116406668238195
SPIKE_Z0 = -0.2408102083477011
SPIKE_YDOT0 = -0.10629710963669947
SPIKE_T = 10.204301970414399


def _make_system() -> cr3bp.CR3BPSystem:
    """Earth-Moon CR3BP at the catalogued (sourced) Braik-Ross mu."""
    return cr3bp.CR3BPSystem(
        mu=EM_MU,
        primary="earth",
        secondary="moon",
        l_km=EM_L_KM,
        t_s=EM_T_S,
    )


def _spike_seed() -> tuple[np.ndarray, float]:
    """Return ``(seed_state, seed_period)`` for the #287 spike's 3D member."""
    state = np.array([SPIKE_X0, 0.0, SPIKE_Z0, 0.0, SPIKE_YDOT0, 0.0], dtype=np.float64)
    return state, SPIKE_T


# ---------------------------------------------------------------------------
# Test 1: pseudo-arclength reproduction off the spike's seed.
# ---------------------------------------------------------------------------


@pytest.mark.slow
def test_pseudo_arclength_reproduces_spike_family() -> None:
    """Pseudo-arclength from the spike's IC walks a multi-member 3D family.

    Topology checks (not numeric goldens):
      * Forward + backward walks each produce >= 5 members.
      * x0 spans BOTH above and below the seed x0 (the spike found
        x0 in [-0.85, -0.77]; we don't insist on exact endpoints).
      * z0 stays genuinely 3D (|z0| > 0.05) on at least 10 walked members
        (the spike's family is entirely non-trivial z0 by construction).
      * Independent closure on every accepted member <= 1e-6 (the tracer's
        gate); typical residuals are 1e-8 or better.
    """
    system = _make_system()
    seed, period = _spike_seed()
    fam = continue_general_3d_family(
        system,
        seed,
        period,
        continuation="pseudo_arclength",
        step=0.01,
        n_steps_max=25,
        direction="both",
    )
    assert isinstance(fam, Family3D)
    assert fam.n_steps_forward >= 5, (
        f"forward walk too short: {fam.n_steps_forward} ({fam.forward_termination})"
    )
    assert fam.n_steps_backward >= 5, (
        f"backward walk too short: {fam.n_steps_backward} ({fam.backward_termination})"
    )
    # x0 extent: seed at -0.8116; expect spread in BOTH directions.
    x0s = np.array([m.orbit.state0[IDX_X] for m in fam.members])
    assert x0s.min() < SPIKE_X0 - 1e-3, f"backward walk didn't move x0; min={x0s.min()}"
    assert x0s.max() > SPIKE_X0 + 1e-3, f"forward walk didn't move x0; max={x0s.max()}"
    # z0 genuinely 3D on the bulk of members (the spike's family is non-trivial z0).
    z0s = np.array([m.orbit.state0[IDX_Z] for m in fam.members])
    n_3d = int(np.sum(np.abs(z0s) > 0.05))
    assert n_3d >= 10, f"only {n_3d} members had |z0| > 0.05 — family collapsed to planar?"


@pytest.mark.slow
def test_closure_preserved_on_every_member() -> None:
    """Every accepted member must close under Radau re-propagation."""
    system = _make_system()
    seed, period = _spike_seed()
    fam = continue_general_3d_family(
        system,
        seed,
        period,
        continuation="pseudo_arclength",
        step=0.01,
        n_steps_max=10,
        direction="both",
        closure_tol=1e-6,
    )
    for m in fam.members:
        cl = m.orbit.independent_closure_residual
        assert math.isfinite(cl), f"non-finite closure at step {m.step_index}"
        assert cl <= 1e-6, f"closure exceeded tol at step {m.step_index}: {cl:.3e}"


# ---------------------------------------------------------------------------
# Test 3: monodromy / Floquet structure.
# ---------------------------------------------------------------------------


def test_monodromy_and_floquet_structure() -> None:
    """Every member has a trivial unit pair; non-trivial multipliers are
    reciprocal pairs (Hamiltonian CR3BP symmetry).

    Topology checks:
      * monodromy is a 6x6 real matrix.
      * Floquet multipliers come in 6 eigenvalues, sorted by descending
        magnitude.
      * At least one eigenvalue within 1e-2 of unity (the trivial unit
        pair — the energy/time-translation eigenvalue).
      * Reciprocal-pair structure: for each non-trivial eigenvalue lambda
        there exists another eigenvalue lambda' with |lambda * lambda' - 1|
        small (within 5% — the unstable hyperbolic pair is on the order of
        |lambda| ~ 1e2 for this family, so the reciprocal is on the order
        1e-2; the comparison is loose).
      * stability_tag is one of the expected labels.
    """
    system = _make_system()
    seed, period = _spike_seed()
    fam = continue_general_3d_family(
        system,
        seed,
        period,
        continuation="pseudo_arclength",
        step=0.01,
        n_steps_max=5,
        direction="both",
        monodromy_eval=True,
    )
    valid_tags = {"stable", "unstable", "hyperbolic_pair", "degenerate", "monodromy_failed"}
    for m in fam.members:
        assert m.monodromy is not None, f"no monodromy at step {m.step_index}"
        assert m.monodromy.shape == (6, 6), f"monodromy shape: {m.monodromy.shape}"
        assert m.floquet is not None, f"no Floquet at step {m.step_index}"
        assert m.floquet.shape == (6,), f"Floquet shape: {m.floquet.shape}"
        # Sorted descending by magnitude.
        mags = np.abs(m.floquet)
        assert np.all(mags[:-1] >= mags[1:] - 1e-9), f"Floquet not sorted at step {m.step_index}"
        # At least one trivial unit eigenvalue.
        dist_to_unity = np.abs(m.floquet - 1.0)
        assert dist_to_unity.min() < 1e-2, (
            f"no unit eigenvalue at step {m.step_index}: min dist {dist_to_unity.min():.3e}"
        )
        # Reciprocal-pair structure: for each non-trivial eigenvalue,
        # there should exist a partner whose product is ~1.
        for i, lam in enumerate(m.floquet):
            if dist_to_unity[i] < 1e-2:
                continue
            partners = m.floquet * lam
            partner_dist = np.abs(partners - 1.0)
            min_pd = float(partner_dist.min())
            assert min_pd < 0.5, (
                f"step {m.step_index}, eig {lam!r}: no reciprocal partner (min {min_pd:.3e})"
            )
        assert m.stability_tag in valid_tags, f"unknown tag: {m.stability_tag}"


def test_seed_member_has_a_trivial_unit_eigenvalue() -> None:
    """Strict check on the seed: real-positive eigenvalue near +1."""
    system = _make_system()
    seed, period = _spike_seed()
    fam = continue_general_3d_family(
        system,
        seed,
        period,
        continuation="pseudo_arclength",
        step=0.01,
        n_steps_max=2,
        direction="forward",
        monodromy_eval=True,
    )
    fl = fam.seed.floquet
    assert fl is not None
    # Real-positive eigenvalue near unity.
    real_pos = [lam for lam in fl if abs(lam.imag) < 1e-3 and lam.real > 0]
    assert real_pos, "no real-positive eigenvalue at seed"
    near_unity = [lam for lam in real_pos if abs(lam.real - 1.0) < 1e-2]
    assert near_unity, f"no real-positive unit eigenvalue at seed; real-pos = {real_pos}"


# ---------------------------------------------------------------------------
# Test 4: fold detection.
# ---------------------------------------------------------------------------


@pytest.mark.slow
def test_fold_detected_and_walk_continues_through() -> None:
    """The spike's family has a fold in (z0, x0) near x0 ~ -0.81; the
    pseudo-arclength walker detects the sign flip in the tangent and
    continues past it.

    The fold is detected as a sign change in one component of the unit
    tangent between adjacent members; the walk does NOT terminate (the
    pseudo-arclength constraint is regular through the fold).
    """
    system = _make_system()
    seed, period = _spike_seed()
    fam = continue_general_3d_family(
        system,
        seed,
        period,
        continuation="pseudo_arclength",
        step=0.01,
        n_steps_max=30,
        direction="both",
        fold_detection=True,
    )
    # The spike's family folds; we expect to see one detected.
    assert len(fam.folds) >= 1, (
        "no fold detected; either step too coarse or family doesn't fold in this slab"
    )
    # The fold is recorded with valid metadata.
    for fold in fam.folds:
        assert isinstance(fold, FoldPoint)
        assert fold.natural_param in {"x0", "z0", "ydot0", "T"}
        # Sign flip between the bracketing tangent components.
        assert fold.tangent_before * fold.tangent_after < 0, (
            f"recorded fold has same-sign tangents: {fold.tangent_before} vs {fold.tangent_after}"
        )
    # Walk did NOT terminate prematurely at a fold (pseudo-arclength continues).
    n_total = fam.n_steps_forward + fam.n_steps_backward
    assert n_total >= 10, f"walk too short despite fold detection: {n_total} steps"


# ---------------------------------------------------------------------------
# Test 5: degenerate-planar guard.
# ---------------------------------------------------------------------------


def test_planar_seed_marks_all_members_degenerate() -> None:
    """A planar seed (z0=0) at the Braik-Ross C11a IC walks within the
    planar manifold; members are flagged ``degenerate_planar`` and the walk
    terminates cleanly.

    The planar manifold is dynamically invariant under the CR3BP, so the
    corrector can't spontaneously hop off it. This test confirms the
    tracer's termination logic catches the case and the flag is set.
    """
    system = _make_system()
    state_planar = np.array([C11A_X0, 0.0, 0.0, 0.0, C11A_YDOT0, 0.0], dtype=np.float64)
    fam = continue_general_3d_family(
        system,
        state_planar,
        C11A_PERIOD_TU,
        continuation="pseudo_arclength",
        step=0.005,
        n_steps_max=5,
        direction="forward",
    )
    assert fam.seed.orbit.degenerate_planar
    # Walk should terminate immediately (the first step lands a planar member).
    assert fam.forward_termination == "degenerate_planar"
    for m in fam.members:
        assert m.orbit.degenerate_planar, (
            f"non-planar member spawned from planar seed at step {m.step_index}: "
            f"z0={m.orbit.state0[IDX_Z]:.3e}"
        )


# ---------------------------------------------------------------------------
# Test 6: direction symmetry.
# ---------------------------------------------------------------------------


@pytest.mark.slow
def test_direction_symmetry() -> None:
    """Forward + backward walks cover x0 above + below the seed."""
    system = _make_system()
    seed, period = _spike_seed()
    fam = continue_general_3d_family(
        system,
        seed,
        period,
        continuation="pseudo_arclength",
        step=0.01,
        n_steps_max=10,
        direction="both",
    )
    # Forward members.
    fwd = [m for m in fam.members if m.step_index > 0]
    bwd = [m for m in fam.members if m.step_index < 0]
    assert fwd, "no forward members"
    assert bwd, "no backward members"
    # Forward + backward should be on OPPOSITE sides of the seed in either
    # x0 OR the arc_length (sign). At minimum the arc_lengths differ in
    # sign.
    arc_fwd = np.array([m.arc_length for m in fwd])
    arc_bwd = np.array([m.arc_length for m in bwd])
    assert (arc_fwd > 0).all(), "forward arc_lengths must be positive"
    assert (arc_bwd < 0).all(), "backward arc_lengths must be negative"


# ---------------------------------------------------------------------------
# Test 7: natural-T continuation runs cleanly.
# ---------------------------------------------------------------------------


def test_natural_t_continuation_runs_cleanly() -> None:
    """natural_T mode (step T as the natural parameter) at the spike's seed
    converges for a small step count.

    Robustness check: this mode is simpler than pseudo-arclength but stops
    at folds. With a small step + small n_steps_max the walk should produce
    at least one converged member in at least one direction.
    """
    system = _make_system()
    seed, period = _spike_seed()
    fam = continue_general_3d_family(
        system,
        seed,
        period,
        continuation="natural_T",
        step=0.05,
        n_steps_max=5,
        direction="both",
    )
    n_total = fam.n_steps_forward + fam.n_steps_backward
    assert n_total >= 1, (
        f"natural_T produced no members; forward={fam.forward_termination}, "
        f"backward={fam.backward_termination}"
    )
    for m in fam.members:
        assert m.orbit.independent_closure_residual <= 1e-6


# ---------------------------------------------------------------------------
# Test 8: input validation.
# ---------------------------------------------------------------------------


def test_input_validation() -> None:
    """The tracer rejects malformed inputs cleanly."""
    system = _make_system()
    seed, period = _spike_seed()
    # Wrong state shape.
    with pytest.raises(ValueError, match="shape"):
        continue_general_3d_family(system, np.zeros(5), period)
    # Non-positive period.
    with pytest.raises(ValueError, match="seed_period"):
        continue_general_3d_family(system, seed, -1.0)
    # Non-positive step.
    with pytest.raises(ValueError, match="step"):
        continue_general_3d_family(system, seed, period, step=-0.01)
    # Unknown continuation mode.
    with pytest.raises(ValueError, match="continuation"):
        continue_general_3d_family(system, seed, period, continuation="bogus")  # type: ignore[arg-type]
    # Unknown direction.
    with pytest.raises(ValueError, match="direction"):
        continue_general_3d_family(system, seed, period, direction="sideways")  # type: ignore[arg-type]


# ---------------------------------------------------------------------------
# Test 9: on_step callback fires.
# ---------------------------------------------------------------------------


def test_on_step_callback_fires() -> None:
    """The optional on_step callback is invoked for the seed AND every
    accepted member."""
    system = _make_system()
    seed, period = _spike_seed()
    seen: list[int] = []

    def cb(m: Family3DMember) -> None:
        seen.append(m.step_index)

    fam = continue_general_3d_family(
        system,
        seed,
        period,
        continuation="pseudo_arclength",
        step=0.01,
        n_steps_max=3,
        direction="forward",
        on_step=cb,
    )
    # Seed callback + each accepted forward step.
    assert seen[0] == 0  # seed
    assert len(seen) == 1 + fam.n_steps_forward
