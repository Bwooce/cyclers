"""Tests for per-point torus stability + manifold generation (#522).

Reuses the sourced #299 Neimark-Sacker bracket fixture already validated by
``tests/genome/test_qp_tori.py`` (not re-deriving a seed here) as the real
CR3BP data source. The structural positive control is the CR3BP monodromy's
RECIPROCAL-PAIR eigenvalue property, a documented mathematical fact of
Hamiltonian/symplectic STMs (see ``search/bifurcation_detector.py``'s own
"reciprocal-pair structure" discussion) -- not a value this module computed.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pytest

import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.genome.qp_tori import QPTorus, correct_qp_torus
from cyclerfinder.genome.qp_torus_manifold import (
    local_stability,
    torus_manifold_grid,
    torus_point_stm,
)

DATA_DIR = Path(__file__).resolve().parents[2] / "data"
SUBFAMILIES_FILE = DATA_DIR / "family_296_3d_subfamilies_299.jsonl"
PARENT_FAMILY_FILE = DATA_DIR / "family_296_3d_em_11.jsonl"

EM_MU = 1.2150584270572e-2
EM_L_KM = 384400.0
EM_T_S = 375699.8


def _em_system() -> cr3bp.CR3BPSystem:
    return cr3bp.CR3BPSystem(mu=EM_MU, primary="earth", secondary="moon", l_km=EM_L_KM, t_s=EM_T_S)


def _load_parent_at_step(step_index: int) -> dict[str, Any]:
    if not PARENT_FAMILY_FILE.exists():
        pytest.skip(f"parent family file not present: {PARENT_FAMILY_FILE}")
    with PARENT_FAMILY_FILE.open() as f:
        for line in f:
            d = json.loads(line)
            if d.get("type") == "header":
                continue
            if d.get("step_index") == step_index:
                return dict(d)
    raise RuntimeError(f"step_index={step_index} not found in {PARENT_FAMILY_FILE}")


def _load_first_neimark_sacker_bracket() -> tuple[dict[str, Any], dict[str, Any]]:
    if not SUBFAMILIES_FILE.exists():
        pytest.skip(f"subfamilies file not present: {SUBFAMILIES_FILE}")
    with SUBFAMILIES_FILE.open() as f:
        for line in f:
            d = json.loads(line)
            if d.get("type") == "header":
                continue
            br = d.get("bracket")
            if br is None:
                continue
            if br.get("classification") == "neimark_sacker":
                parent = _load_parent_at_step(int(br["step_a"]))
                return br, parent
    raise RuntimeError("no Neimark-Sacker bracket found in subfamilies file")


def _sourced_torus() -> QPTorus:
    br, parent = _load_first_neimark_sacker_bracket()
    system = _em_system()
    parent_state = np.asarray(parent["state_nd"], dtype=np.float64)
    parent_period = float(parent["T_TU"])
    k = int(br["k"])
    lam_a = complex(br["eig_a_re"], br["eig_a_im"])
    lam_b = complex(br["eig_b_re"], br["eig_b_im"])
    torus = correct_qp_torus(
        system,
        parent_state,
        parent_period,
        (lam_a, lam_b),
        k=k,
        n_long=16,
        n_trans=2,
        initial_torus_amplitude=5e-4,
        tol=1e-8,
        max_iter=40,
        independent_tol=1e-3,
        notes="qp_torus_manifold_test",
    )
    # Matches tests/genome/test_qp_tori.py::test_sourced_neimark_sacker_smoke's
    # own gates directly -- torus.converged uses the STRICTER tol=1e-8 passed
    # above as the correction target, which this bracket does not always hit
    # exactly, but the actual published-practice gates below are what matter.
    assert torus.invariance_residual < 1e-5
    assert torus.independent_closure_residual < 1e-3
    return torus


def test_torus_point_stm_eigenvalues_are_reciprocal_paired() -> None:
    """The CR3BP STM over any time span is the Jacobian of a Hamiltonian
    flow, hence symplectic -- its eigenvalues occur in reciprocal pairs
    (lambda, 1/lambda). This is a structural property, not something this
    module computed; it holds regardless of which torus point the STM is
    evaluated at.
    """
    torus = _sourced_torus()
    _state, stm = torus_point_stm(torus, 0.7, 1.3)
    eigvals = np.linalg.eigvals(stm)
    # For each eigenvalue, ITS reciprocal (1/lambda) must be present
    # somewhere else in the set -- checked via nearest-match rather than an
    # assumed sort order, since near-unit-circle complex-conjugate pairs
    # (e.g. the trivial energy/phase pair coexisting with a genuine
    # oscillatory stable mode) can tie in |lambda| and interleave under a
    # magnitude sort, which does NOT imply they pair with EACH OTHER.
    n = eigvals.shape[0]
    reciprocals = 1.0 / eigvals
    used = np.zeros(n, dtype=bool)
    for i in range(n):
        dists = np.abs(eigvals - reciprocals[i])
        dists[used] = np.inf
        j = int(np.argmin(dists))
        assert dists[j] < 1e-4, (
            f"eigenvalue {eigvals[i]} has no reciprocal-pair partner in the "
            f"spectrum (closest candidate {eigvals[j]}, distance {dists[j]}) "
            f"-- violates the symplectic STM reciprocal-pair structure"
        )
        used[j] = True


def test_local_stability_eigenvector_continuity_sign_fix() -> None:
    torus = _sourced_torus()
    _state, stm = torus_point_stm(torus, 0.7, 1.3)
    stab = local_stability(_state, stm)
    if stab.vec_u is None:
        pytest.skip(
            "this bracket's torus point has no real unstable direction (locally non-hyperbolic)"
        )
    # Feed the NEGATED vector as "previous" -- continuity must flip it back.
    stab_flipped = local_stability(_state, stm, prev_vec_u=-stab.vec_u)
    assert stab_flipped.vec_u is not None
    assert np.allclose(stab_flipped.vec_u, -stab.vec_u, atol=1e-12)
    # Feeding the SAME-sign vector as "previous" must keep the sign.
    stab_same = local_stability(_state, stm, prev_vec_u=stab.vec_u)
    assert stab_same.vec_u is not None
    assert np.allclose(stab_same.vec_u, stab.vec_u, atol=1e-12)


def test_torus_manifold_grid_runs_and_shapes_are_consistent() -> None:
    """Smoke test: a small grid runs end-to-end without crashing and
    produces consistently-shaped output. Does NOT assert the manifold is
    geometrically meaningful (this sourced bracket's torus is not
    guaranteed a priori to be hyperbolic -- #530 found only ~17% of sampled
    orbits genuinely unstable in an unrelated system) -- that is validated
    separately once a genuinely-unstable case is identified.
    """
    torus = _sourced_torus()
    grid = torus_manifold_grid(
        torus,
        n_long=4,
        n_lat=4,
        branch="unstable",
        eps=1e-6,
        surface_x=1.0 - EM_MU,
        t_max=5.0,
    )
    assert grid.origins.shape == (4, 4, 2)
    assert grid.endpoints.shape == (4, 4, 6)
    assert grid.hyperbolic.shape == (4, 4)
    # Every endpoint is either a valid finite crossing state or all-NaN.
    for i in range(4):
        for j in range(4):
            row = grid.endpoints[i, j, :]
            assert np.all(np.isfinite(row)) or np.all(np.isnan(row))
