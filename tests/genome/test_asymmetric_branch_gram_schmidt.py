"""#347 Phase 2 P2.1 — Gram-Schmidt fix in _select_saddle_center_eigenvector.

The Phase 1 (3, 2) Earth-Moon C32 saddle-center is harmless against the
trivial time-translation eigenvector: the (3, 2) parent's tangent
``ẋ₀ = f(state0, mu)`` is purely (ẋ, ẏ) because the IC has y=z=ẋ=ż=0 by
construction (perpendicular x-axis crossing); the marginal saddle-center
eigenvector at i=124 is purely (z, ż). These two subspaces are orthogonal
by symmetry, so omitting Gram-Schmidt was a no-op for the Phase 1 anchor.

For general Phase 2 sweep targets — different (k1, k2) families, off-axis
ICs, families whose saddle-center direction has in-plane components —
the secondary-pair eigenvector returned by ``np.linalg.eig`` is degenerate
against the trivial pair (clustered eigenvalues) and is *numerically*
rotated to carry a non-zero ẋ₀ component. Perturbing along an
un-orthogonalised eigenvector slides the corrector along the parent
family's own orbit, not off it.

Gates:

  1. Synthetic non-orthogonal eigenvector: a hand-constructed monodromy
     whose returned secondary-pair eigenvector has a deliberate ẋ₀
     component; verify the Gram-Schmidt fix subtracts it cleanly.
  2. The (3, 2) Phase 1 anchor: with the Gram-Schmidt fix in place, the
     i=124 eigenvector is UNCHANGED (in-plane vs out-of-plane subspaces are
     orthogonal). The Phase 1 PASS verdict must survive the fix.
  3. The (3, 2) Phase 1 ``branch_at_saddle_center`` end-to-end: branching
     from i=124 at eps=5e-4 still converges with the Gram-Schmidt fix in
     place (the actual orthogonalisation is a no-op at the (3, 2) anchor,
     so the converged orbit's residual must match Phase 1's < 1e-10 gate).
  4. Degenerate guard: if the candidate eigenvector is PARALLEL to the
     parent tangent, the function returns ``None`` (the saddle-center is
     not cleanly distinguishable from the time-translation symmetry — the
     caller treats as no-pick).
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.genome.asymmetric_branch import (
    _select_saddle_center_eigenvector,
    branch_at_saddle_center,
)
from cyclerfinder.search.bifurcation_detector import monodromy
from cyclerfinder.search.reachable_representatives import braik_ross_system

_REPO_ROOT = Path(__file__).resolve().parent.parent.parent
ARTIFACT_PATH = _REPO_ROOT / "data" / "floquet_phase1_c32_family.jsonl"


def _load_member_at_index(idx: int) -> dict:  # type: ignore[type-arg]
    with ARTIFACT_PATH.open() as fh:
        rows = [json.loads(line) for line in fh]
    members = [r for r in rows if r.get("kind") == "member"]
    result: dict = members[idx]  # type: ignore[type-arg]
    return result


def _state0_for_member(row: dict) -> np.ndarray:  # type: ignore[type-arg]
    return np.array(
        [float(row["x0"]), 0.0, 0.0, 0.0, float(row["ydot0"]), 0.0],
        dtype=np.float64,
    )


def test_gram_schmidt_subtracts_parent_tangent_component() -> None:
    """Synthetic non-orthogonal eigenvector: the fix removes the ẋ₀ component.

    Constructs a synthetic 6x6 monodromy with 4 distinct real eigenvalues + a
    near-+1 eigenvalue whose eigenvector has a deliberately large component
    along a chosen ``parent_tangent``. Without Gram-Schmidt the returned
    eigenvector retains that component; with Gram-Schmidt the returned
    eigenvector is exactly orthogonal to the tangent.
    """
    # Pick eigenvalues that satisfy the secondary-pair classification:
    # - 2 close-to-1 (trivial pair): 1.0 + 1e-6, 1.0 - 1e-6 (the closest)
    # - 2 large |log|lambda|| (primary saddle): 3.0, 1/3.0
    # - 2 secondary near +1 (saddle-center direction): 1.01, 0.985
    # The CLOSER-TO-+1 secondary is 1.01 (|1.01-1|=0.01 < |0.985-1|=0.015), so the
    # detector picks the index-4 eigenvector (the one we attach the parent_tangent
    # contamination to).
    eigvals = np.array([1.0 + 1e-6, 1.0 - 1e-6, 3.0, 1.0 / 3.0, 1.01, 0.985])
    # Construct eigenvectors. The secondary "1.01" eigenvector deliberately has
    # a component along a chosen parent_tangent direction.
    parent_tangent = np.array([0.0, 0.0, 0.0, 0.4, 0.5, 0.0], dtype=np.float64)
    # The dirty secondary eigenvector: (z, ż) subspace + a deliberate
    # parent_tangent component. Choose it well-separated from the other columns
    # to keep the eigenvector matrix well-conditioned.
    raw_secondary = np.array([0.0, 0.0, 0.8, 0.4, 0.5, 0.6], dtype=np.float64)
    raw_secondary /= np.linalg.norm(raw_secondary)
    # Trivial + primary eigenvectors: pick orthogonal canonical directions.
    eigvecs = np.zeros((6, 6), dtype=np.float64)
    eigvecs[:, 0] = np.array([1.0, 0.0, 0.0, 0.0, 0.0, 0.0])  # trivial
    eigvecs[:, 1] = np.array([0.0, 1.0, 0.0, 0.0, 0.0, 0.0])  # trivial
    eigvecs[:, 2] = np.array([0.0, 0.0, 0.0, 1.0, 0.0, 0.0])  # primary
    eigvecs[:, 3] = np.array([0.0, 0.0, 0.0, 0.0, 1.0, 0.0])  # primary
    eigvecs[:, 4] = raw_secondary  # the eigenvector with the dirty parent_tangent component
    eigvecs[:, 5] = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 1.0])  # complementary

    # M = V * diag(eigvals) * V^-1.
    mono_mat = eigvecs @ np.diag(eigvals) @ np.linalg.inv(eigvecs)

    # Without parent_tangent: the returned eigenvector retains the parent_tangent component.
    pick_no_gs = _select_saddle_center_eigenvector(mono_mat)
    assert pick_no_gs is not None
    v_no_gs, lam_no_gs = pick_no_gs
    # The pre-Gram-Schmidt eigenvector should be (proportional to) the raw secondary.
    overlap_no_gs = abs(float(np.dot(v_no_gs, parent_tangent / np.linalg.norm(parent_tangent))))
    assert overlap_no_gs > 0.05, (
        f"sanity: synthetic eigenvector should have a non-trivial parent_tangent component "
        f"(overlap={overlap_no_gs:.3e})"
    )

    # With parent_tangent: the returned eigenvector is orthogonal to parent_tangent.
    pick_gs = _select_saddle_center_eigenvector(mono_mat, parent_tangent=parent_tangent)
    assert pick_gs is not None
    v_gs, lam_gs = pick_gs
    overlap_gs = abs(float(np.dot(v_gs, parent_tangent / np.linalg.norm(parent_tangent))))
    assert overlap_gs < 1e-10, (
        f"post-Gram-Schmidt eigenvector should be orthogonal to parent_tangent "
        f"(overlap={overlap_gs:.3e})"
    )
    # Eigenvalue selection is unchanged by Gram-Schmidt.
    assert abs(lam_gs - lam_no_gs) < 1e-12, "Gram-Schmidt must not change the eigenvalue selection"
    # Eigenvector is unit-normalised.
    assert abs(np.linalg.norm(v_gs) - 1.0) < 1e-12


def test_gram_schmidt_is_noop_at_c32_i124() -> None:
    """The (3, 2) Phase 1 anchor's eigenvector is unchanged by Gram-Schmidt.

    At i=124 the parent tangent ``ẋ₀ = f(state0, mu)`` is purely (ẋ, ẏ): the IC
    has y=z=ẋ=ż=0 by the perpendicular-x-axis-crossing construction. The
    marginal eigenvector is purely (z, ż). These subspaces are orthogonal by
    construction. The Gram-Schmidt step subtracts a quantity proportional to
    the inner product (which is ~0), leaving the eigenvector unchanged.

    Phase 1's verdict (PASS) must survive the introduction of the fix; if the
    no-op were silently corrupting the eigenvector this test would catch it.
    """
    m = _load_member_at_index(124)
    system = braik_ross_system()
    state0 = _state0_for_member(m)
    mono = monodromy(system, state0, float(m["period_TU"]))
    pick_no_gs = _select_saddle_center_eigenvector(mono)
    assert pick_no_gs is not None
    v_no_gs, lam_no_gs = pick_no_gs

    parent_tangent = cr3bp.cr3bp_eom(0.0, state0, system.mu)
    pick_gs = _select_saddle_center_eigenvector(mono, parent_tangent=parent_tangent)
    assert pick_gs is not None
    v_gs, lam_gs = pick_gs

    # Sign-ambiguity-tolerant comparison: |v_gs ± v_no_gs|.
    # The (3, 2) i=124 in-plane components are ~6e-6 at the cluster point
    # (variable-step DOP853 contamination, P372.3 cross-check); Gram-Schmidt
    # against the parent tangent removes that level of in-plane noise floor.
    # That's the WANTED behaviour: we drop in-plane noise from the eigenvector,
    # which makes the perturbation slightly cleaner — but the dominant (z, ż)
    # components are preserved.
    delta_pos = float(np.linalg.norm(v_gs - v_no_gs))
    delta_neg = float(np.linalg.norm(v_gs + v_no_gs))
    delta = min(delta_pos, delta_neg)
    assert delta < 1e-4, (
        f"(3, 2) i=124 eigenvector changed by Gram-Schmidt (delta={delta:.3e}) — "
        f"the in-plane vs out-of-plane subspaces should be near-orthogonal by symmetry; "
        f"~6e-6 in-plane noise floor is expected per P372.3 P-R 2016 finding"
    )
    # The dominant (z, ż) components must be unchanged at the 1e-4 level.
    z_no_gs = float(v_no_gs[2])
    z_gs = float(v_gs[2])
    zdot_no_gs = float(v_no_gs[5])
    zdot_gs = float(v_gs[5])
    # Account for sign ambiguity by taking the sign that minimises the delta.
    sign = 1.0 if (z_no_gs * z_gs) > 0 else -1.0
    assert abs(sign * z_gs - z_no_gs) < 1e-6
    assert abs(sign * zdot_gs - zdot_no_gs) < 1e-6
    assert abs(lam_gs - lam_no_gs) < 1e-12


def test_phase1_branch_at_saddle_center_survives_gram_schmidt() -> None:
    """The Phase 1 P1.4 end-to-end branch (3, 2) i=124 eps=5e-4 still converges.

    With the Gram-Schmidt fix wired into ``branch_at_saddle_center`` (the
    function now always computes parent_tangent and passes it through to
    ``_select_saddle_center_eigenvector``), the Phase 1 PASS gate must
    survive: corrector residual < 1e-10, topology changes from (3, 2) to a
    different (k1', k2').
    """
    m = _load_member_at_index(124)
    system = braik_ross_system()
    state0 = _state0_for_member(m)
    period = float(m["period_TU"])
    result = branch_at_saddle_center(system, state0, period, epsilon=5e-4)
    assert result is not None, (
        "Phase 1 PASS gate regressed: branch_at_saddle_center returned None at i=124 eps=5e-4 "
        "after the Gram-Schmidt fix landed"
    )
    bo = result.branched_orbit
    assert bo.corrector_residual < 1e-10, (
        f"Phase 1 residual gate regressed: corrector residual {bo.corrector_residual:.3e} >= 1e-10"
    )
    assert result.topology_changed, (
        f"Phase 1 topology-change gate regressed: branched topology "
        f"({result.branched_topology.k1}, {result.branched_topology.k2}) == parent (3, 2)"
    )


def test_gram_schmidt_returns_none_for_parallel_eigenvector() -> None:
    """Degenerate guard: when the candidate eigenvector is parallel to ẋ₀, return None.

    If the secondary-pair eigenvector is exactly proportional to the parent
    tangent (a pathological case where the Floquet eigendecomposition has
    aliased the trivial pair into the secondary pair), Gram-Schmidt zeros it
    out and the function returns ``None`` rather than emitting a degenerate
    perturbation direction.

    Construction: build a DIAGONAL monodromy (so eigvec[:, i] is the i-th
    canonical basis vector) and pick the parent_tangent to align EXACTLY with
    one of the canonical basis vectors that maps to the secondary eigenvalue
    we're testing.
    """
    # Diagonal monodromy. Eigenvalue ordering:
    # idx 0: 1.0 + 1e-6 (trivial)
    # idx 1: 1.0 - 1e-6 (trivial)
    # idx 2: 3.0 (primary)
    # idx 3: 1/3.0 (primary)
    # idx 4: 1.01 (secondary, closest to +1)
    # idx 5: 0.985 (secondary, further from +1)
    # With a diagonal matrix np.linalg.eig returns eigvec = identity (each
    # column is a canonical basis vector).
    mono_mat = np.diag([1.0 + 1e-6, 1.0 - 1e-6, 3.0, 1.0 / 3.0, 1.01, 0.985])
    # Parent tangent: aligned with idx 4's canonical eigenvector (e_4 = unit
    # vector in component 4). The Gram-Schmidt will subtract the full
    # eigenvector and leave a zero residual.
    parent_tangent = np.array([0.0, 0.0, 0.0, 0.0, 1.0, 0.0], dtype=np.float64)

    pick = _select_saddle_center_eigenvector(mono_mat, parent_tangent=parent_tangent)
    assert pick is None, (
        "expected None when the secondary eigenvector is parallel to parent tangent"
    )


def test_gram_schmidt_validates_parent_tangent_shape() -> None:
    """Defensive contract: a malformed parent_tangent raises ValueError."""
    import pytest

    mono_mat = np.eye(6, dtype=np.float64)
    with pytest.raises(ValueError, match="parent_tangent"):
        _select_saddle_center_eigenvector(mono_mat, parent_tangent=np.array([1.0, 2.0]))
