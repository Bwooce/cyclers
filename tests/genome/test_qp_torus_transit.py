"""Smoke tests for the empirical transit-branch torus-manifold grid (#548).

Reuses the sourced #299 Neimark-Sacker bracket fixture (via the same loader as
``tests/genome/test_qp_torus_manifold.py``) to build a small, fast QP-torus and
check that :func:`genome.qp_torus_transit.transit_torus_manifold_grid` runs
end-to-end, returns consistently shaped output, and reports a branch-sign field
whose non-transit points are exactly the NaN endpoints. Does NOT assert the
manifold is geometrically meaningful for this particular bracket (its torus is
not guaranteed hyperbolic a priori -- the geometric positive control lives in
the #548 sweep script, not the unit suite).
"""

from __future__ import annotations

import json
import pathlib
from typing import Any

import numpy as np
import pytest

import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.genome.qp_tori import QPTorus, correct_qp_torus
from cyclerfinder.genome.qp_torus_transit import transit_torus_manifold_grid

EM_MU = 0.012153643
_DATA = pathlib.Path(__file__).resolve().parents[2] / "data"
PARENT_FAMILY_FILE = _DATA / "family_296_3d_em_11.jsonl"
SUBFAMILIES_FILE = _DATA / "family_296_3d_subfamilies_299.jsonl"


def _em_system() -> cr3bp.CR3BPSystem:
    return cr3bp.CR3BPSystem(
        mu=EM_MU, primary="Earth", secondary="Moon", l_km=384400.0, t_s=382981.0
    )


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
    torus = correct_qp_torus(
        system,
        parent_state,
        parent_period,
        (complex(br["eig_a_re"], br["eig_a_im"]), complex(br["eig_b_re"], br["eig_b_im"])),
        k=int(br["k"]),
        n_long=16,
        n_trans=2,
        initial_torus_amplitude=5e-4,
        tol=1e-8,
        max_iter=40,
        independent_tol=1e-3,
        notes="qp_torus_transit_test",
    )
    assert torus.invariance_residual < 1e-4
    return torus


def test_transit_grid_shapes_and_sign_consistency() -> None:
    torus = _sourced_torus()
    grid, signs = transit_torus_manifold_grid(
        torus,
        n_long=4,
        n_lat=4,
        branch="unstable",
        surface_x=1.0 - EM_MU,
        t_max=6.0,
        eps=1e-5,
        max_workers=2,
    )
    assert grid.origins.shape == (4, 4, 2)
    assert grid.endpoints.shape == (4, 4, 6)
    assert grid.hyperbolic.shape == (4, 4)
    assert signs.shape == (4, 4)
    # Sign is +/-1 exactly where a finite crossing was recorded, 0 elsewhere.
    for i in range(4):
        for j in range(4):
            has_cross = bool(np.all(np.isfinite(grid.endpoints[i, j, :])))
            if has_cross:
                assert signs[i, j] in (-1, 1)
                # crossing state sits on the section
                assert abs(grid.endpoints[i, j, 0] - (1.0 - EM_MU)) < 1e-6
            else:
                assert signs[i, j] == 0
                assert np.all(np.isnan(grid.endpoints[i, j, :]))


def test_transit_grid_rejects_bad_branch() -> None:
    torus = _sourced_torus()
    with pytest.raises(ValueError, match="branch must be"):
        transit_torus_manifold_grid(
            torus, n_long=4, n_lat=4, branch="bogus", surface_x=0.9, t_max=3.0
        )
