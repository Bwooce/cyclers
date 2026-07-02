"""Wiring smoke test for the #522 linking-number heteroclinic screen.

Uses the same sourced #299 Earth-Moon Neimark-Sacker torus (see
``tests/genome/test_qp_torus_manifold.py``) for both branches -- this is a
SELF-consistency / mechanical-wiring test (stable manifold of the torus
against its own unstable manifold), not a genuine two-torus heteroclinic
validation. The real positive control (two DISTINCT quasi-halo tori
reproducing Owen & Baresi's Earth-Moon Sec 4.1.1 result, 4 connections at
mu=0.012153643, C=3.15) requires sourcing seeds at their specific published
latitudinal frequencies and is tracked separately in
data/OUTSTANDING.md's #522 entry, not asserted here.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pytest

import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.genome.qp_tori import QPTorus, correct_qp_torus
from cyclerfinder.genome.qp_torus_heteroclinic import build_manifold_grids, scan_linking_number

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
        notes="qp_torus_heteroclinic_test",
    )
    assert torus.invariance_residual < 1e-5
    assert torus.independent_closure_residual < 1e-3
    return torus


def test_scan_linking_number_runs_end_to_end() -> None:
    """Mechanical wiring test: builds manifold grids, scans a scanning
    variable, and asserts the pipeline runs without crashing and returns a
    result of the requested shape. Does NOT assert any specific linking
    number sequence (no genuine two-torus connection is set up here).
    """
    torus = _sourced_torus()
    stable_grid, unstable_grid = build_manifold_grids(
        torus,
        torus,
        n_long=4,
        n_lat=4,
        eps=1e-6,
        surface_x=1.0 - EM_MU,
        t_max=5.0,
    )
    d_values = np.linspace(-0.05, 0.05, 5)
    result = scan_linking_number(
        stable_grid,
        unstable_grid,
        scanning_component="z",
        curve_components=("y", "ydot", "zdot"),
        d_values=d_values,
    )
    assert result.d_values.shape == (5,)
    assert result.linking_numbers.shape == (5,)
    # sign_change_locations must not crash and returns floats within range.
    for loc in result.sign_change_locations():
        assert d_values.min() <= loc <= d_values.max()
