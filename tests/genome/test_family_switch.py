"""Family-switching corrector tests (#266 Phase 3).

Two gates:

  1. **Petal-count k on a real switch** -- the family-switch corrector,
     applied to the L2 Southern NRHO at its k=2 bifurcation point, returns a
     2-petal orbit. This is the end-to-end functional gate (the synthetic-
     eigenvector approach is much harder to set up without inventing a
     CR3BP-like dynamics; the real continuation provides a known-truth
     parent + bifurcation pair).

  2. **No-converge returns None** -- a clearly-bogus eigenvector step size
     either re-converges on the parent (rejected by the petal-count gate) or
     fails to converge; in both cases the corrector returns ``None`` rather
     than raising.

The discipline forbids tuning the eigenvector step or the corrector tolerance
to make the petal-count gate pass.
"""

from __future__ import annotations

import numpy as np

from cyclerfinder.genome.family_switch import (
    _select_period_multiplying_eigenvector,
    switch_family,
)
from cyclerfinder.genome.tulip import KOBLICK_2023_TABLE4_PAPER, koblick_system, petal_count
from cyclerfinder.search.bifurcation_detector import BifurcationPoint, floquet_multipliers
from cyclerfinder.search.nrho_continuation import (
    SymmetricNRHO,
    continue_nrho_family,
    correct_symmetric_nrho,
)


def _build_parent_at_bifurcation() -> tuple[SymmetricNRHO, BifurcationPoint]:
    """Helper: continue the Koblick NRHO family until the closest-to-(-1)
    member and return it + the first k=2 BifurcationPoint bracket.

    Seed source: Koblick 2023 AMOSTECH Table 4 paper row 1 (sourced).
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
        bif_tol=1e-1,
        bif_k_max=4,
    )
    # Closest member to lambda = -1 (the period-doubling bifurcation point).
    best_idx = -1
    best_dist = float("inf")
    for i, m in enumerate(branch.members):
        if m.monodromy is None:
            continue
        eigs = floquet_multipliers(m.monodromy)
        d = min(abs(complex(e) + 1.0) for e in eigs)
        if d < best_dist:
            best_dist = d
            best_idx = i
    k2_bifs = [b for b in branch.bifurcations if b.k == 2]
    assert best_idx >= 0 and k2_bifs, (
        f"bifurcation gate setup failed: best_idx={best_idx}, "
        f"k2_bifs={len(k2_bifs)}, branch_members={len(branch.members)}, "
        f"stop={branch.stop_reason.value}"
    )
    return branch.members[best_idx], k2_bifs[0]


# ---------------------------------------------------------------------------
# Gate 1: petal_count on the switched member equals k.
# ---------------------------------------------------------------------------


def test_switch_family_petal_count_is_k_on_real_bifurcation() -> None:
    """family-switch + correct on the real Koblick NRHO k=2 bifurcation lands
    a 2-petal orbit.

    Test seed source: the parent is derived from the Koblick AMOSTECH Table 4
    paper Np=1 row via :func:`continue_nrho_family`. The bifurcation is the
    first k=2 bracket detected on that continuation. The switch should:

      - converge with closure residual below the gate (1e-8),
      - the switched orbit's Jacobi should be approximately the parent's
        (the bifurcation is along a single energy level),
      - and the petal count of the switched orbit should equal k=2.
    """
    sysm = koblick_system()
    parent, bif = _build_parent_at_bifurcation()
    switched = switch_family(
        parent,
        bif,
        sysm,
        k=2,
        eigenvector_step=1e-2,
        tol=1e-9,
    )
    assert switched is not None, (
        f"switch_family returned None on a real k=2 bifurcation; parent "
        f"x0={parent.x0:.6f}, T={parent.T_TU:.4f}"
    )
    assert switched.converged
    # Topology gate: petal count must equal k=2.
    s0 = np.array([switched.x0, 0, switched.z0, 0, switched.ydot0, 0])
    n = petal_count(s0, switched.T_TU, sysm)
    assert n == 2, (
        f"switched orbit petal_count = {n}, expected 2; "
        f"switched T={switched.T_TU:.4f}, parent T={parent.T_TU:.4f}, "
        f"ratio T/T_parent = {switched.T_TU / parent.T_TU:.4f}"
    )
    # The switched period should be ~2x the parent period (period doubling).
    ratio = switched.T_TU / parent.T_TU
    assert 1.95 < ratio < 2.05, (
        f"switched/parent period ratio {ratio:.4f} not near 2 (period doubling)"
    )


# ---------------------------------------------------------------------------
# Gate 2: bogus step returns None (not exception, not lie).
# ---------------------------------------------------------------------------


def test_switch_family_returns_none_on_bogus_step() -> None:
    """A huge eigenvector_step that pushes the IC clearly off-family should
    fail to converge and return ``None``.

    The corrector may still converge on something (the perturbed state can
    land on yet another nearby periodic orbit) -- in that case the petal-
    count gate rejects the result, also returning ``None``. Either path is
    acceptable: the contract is "no exception, no relabeling of a wrong
    orbit as the switched family".
    """
    sysm = koblick_system()
    parent, bif = _build_parent_at_bifurcation()
    # 100x the normal step in state-vector units -- the perturbed IC is far
    # outside the linearisation neighborhood.
    switched = switch_family(
        parent,
        bif,
        sysm,
        k=2,
        eigenvector_step=10.0,
        tol=1e-9,
    )
    # Either no-converge or wrong petal count -- both yield None.
    assert switched is None


# ---------------------------------------------------------------------------
# Defensive: bad k or missing monodromy raises ValueError.
# ---------------------------------------------------------------------------


def test_switch_family_raises_on_bad_k() -> None:
    """k=1 (or k<2) is rejected."""
    import pytest

    sysm = koblick_system()
    parent, bif = _build_parent_at_bifurcation()
    with pytest.raises(ValueError, match="k must be >= 2"):
        switch_family(parent, bif, sysm, k=1)


def test_switch_family_raises_on_missing_monodromy() -> None:
    """parent.monodromy must not be None."""
    import pytest

    sysm = koblick_system()
    parent, bif = _build_parent_at_bifurcation()
    # Strip the monodromy.
    parent_no_mono = SymmetricNRHO(
        x0=parent.x0,
        z0=parent.z0,
        ydot0=parent.ydot0,
        T_TU=parent.T_TU,
        jacobi=parent.jacobi,
        converged=parent.converged,
        closure_residual=parent.closure_residual,
        n_iter=parent.n_iter,
        monodromy=None,
    )
    with pytest.raises(ValueError, match="monodromy is None"):
        switch_family(parent_no_mono, bif, sysm, k=2)


# ---------------------------------------------------------------------------
# Eigenvector selection: real -1 multiplier yields a real eigenvector.
# ---------------------------------------------------------------------------


def test_select_period_multiplying_eigenvector_picks_minus_one() -> None:
    """A synthetic monodromy with eigenvalue exactly -1 yields the -1
    eigenvector (k=2 pick)."""
    # Build a 6x6 with eigenvalues (-1, 0.5, 2, 1, 1, 0.5) -- two reciprocal
    # pairs and a trivial pair. Eigenvector for -1 is e0.
    diag = np.array([-1.0, 0.5, 2.0, 1.0, 1.0, 0.5])
    mat = np.diag(diag)
    pick = _select_period_multiplying_eigenvector(mat, 2)
    assert pick is not None
    v, lam, dist = pick
    assert abs(lam + 1) < 1e-12, f"lam = {lam}, expected -1"
    assert dist < 1e-12
    # The eigenvector should be e0 (the column corresponding to lambda=-1).
    np.testing.assert_allclose(np.abs(v), np.array([1, 0, 0, 0, 0, 0]), atol=1e-12)


def test_select_period_multiplying_eigenvector_returns_none_when_far() -> None:
    """If no eigenvalue is near a primitive k-th root, returns None."""
    mat = np.diag([2.0, 2.0, 2.0, 2.0, 2.0, 2.0])
    pick = _select_period_multiplying_eigenvector(mat, 2)
    assert pick is None
