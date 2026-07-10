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
import math
from pathlib import Path
from typing import Any

import numpy as np
import pytest

import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.genome.qp_tori import QPTorus, correct_qp_torus
from cyclerfinder.genome.qp_torus_heteroclinic import (
    build_manifold_grids,
    closest_curve_distance,
    scan_linking_number,
)
from cyclerfinder.genome.qp_torus_manifold import ManifoldGrid

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


def test_closest_curve_distance_runs_end_to_end() -> None:
    """Mechanical wiring test (#545) for the deflated-Newton residual: within
    the grids' finite-crossing overlap it returns a finite, non-negative
    distance; outside any valid level curve it returns ``None`` (the same
    "no connection detectable there from this data" convention as
    ``scan_linking_number``). Same self-consistency torus/grid setup as
    ``test_scan_linking_number_runs_end_to_end`` -- a mechanical wiring
    check, not a physical connection claim.
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
    dist = closest_curve_distance(
        stable_grid,
        unstable_grid,
        scanning_component="z",
        curve_components=("y", "ydot", "zdot"),
        d=0.0,
    )
    if dist is not None:
        assert dist >= 0.0
        assert math.isfinite(dist)

    far_outside = closest_curve_distance(
        stable_grid,
        unstable_grid,
        scanning_component="z",
        curve_components=("y", "ydot", "zdot"),
        d=1.0e6,
    )
    assert far_outside is None


# ---------------------------------------------------------------------------
# Synthetic extraction-machinery positive control (#555).
#
# The real-torus scans in #534/#536/#546/#548 produced "linking number
# identically 0", which #553 flagged as ambiguous: `scan_linking_number` also
# emits 0 when `_first_closed_curve` returns None (NaN-heavy grid). These tests
# validate the extraction+linking code path against curves with a KNOWN link,
# and validate the #555 availability instrumentation, WITHOUT depending on any
# expensive CR3BP torus construction -- a genuine 0 must be provably a property
# of curves that were actually extracted.
# ---------------------------------------------------------------------------


def _synthetic_ring_grid(ring_fn: Any, n: int = 40, m: int = 40) -> ManifoldGrid:
    """A synthetic ManifoldGrid whose scanning field (index 3, ``xdot``) is a
    pure function of ``theta_long`` (so a level set is a full ``theta_trans``
    loop) and whose ``(x, y, z)`` trace ``ring_fn(theta_trans)`` along it."""
    tl = np.linspace(0.0, 2.0 * math.pi, n, endpoint=False)
    tt = np.linspace(0.0, 2.0 * math.pi, m, endpoint=False)
    origins = np.zeros((n, m, 2))
    endpoints = np.full((n, m, 6), np.nan)
    for i, a in enumerate(tl):
        for j, b in enumerate(tt):
            x, y, z = ring_fn(b)
            origins[i, j] = (a, b)
            endpoints[i, j] = (x, y, z, math.cos(a), 0.0, 0.0)
    return ManifoldGrid(origins=origins, endpoints=endpoints, hyperbolic=np.ones((n, m), bool))


def test_synthetic_hopf_link_positive_control() -> None:
    """Two interlocking rings (a Hopf link, known linking number +-1) pushed
    through the REAL ``scan_linking_number`` code path must yield a NONZERO
    linking number on a fully-available scan; the same rings pulled apart must
    yield identically 0. This proves the extraction+linking machinery detects a
    genuine link -- so a real-torus 0 is about geometry, not broken code."""
    # Ring A: unit circle in z=0 plane, centroid at the origin.
    grid_a = _synthetic_ring_grid(lambda t: (math.cos(t), math.sin(t), 0.0))
    # Ring B: circle in the xz-plane centered (0.8, 0, 0), r=0.6 -> pierces A's
    # disk at (0.2, 0, 0), which is NOT A's centroid (avoids the documented
    # centroid over-count edge case), Hopf-linking A exactly once.
    grid_b = _synthetic_ring_grid(lambda t: (0.8 + 0.6 * math.cos(t), 0.0, 0.6 * math.sin(t)))
    # Unlinked control: same ring pushed far away.
    grid_far = _synthetic_ring_grid(lambda t: (6.0 + 0.6 * math.cos(t), 0.0, 0.6 * math.sin(t)))
    d_values = np.linspace(-0.85, 0.85, 40)

    linked = scan_linking_number(
        grid_a,
        grid_b,
        scanning_component="xdot",
        curve_components=("x", "y", "z"),
        d_values=d_values,
    )
    unlinked = scan_linking_number(
        grid_a,
        grid_far,
        scanning_component="xdot",
        curve_components=("x", "y", "z"),
        d_values=d_values,
    )

    avail = linked.availability_summary()
    # The whole scan had both curves available (this is what makes the result
    # interpretable -- the exact instrumentation #553 asked for).
    assert avail["both_available"] == avail["n"] == len(d_values)
    # The link is detected: a nonzero linking number somewhere on the scan.
    assert int(np.sum(linked.linking_numbers != 0)) > 0
    assert any(lk != 0 for lk in linked.linking_numbers.tolist())
    # Pulled apart: identically zero, but ALSO fully available -> a genuine
    # "these curves do not link", not a missing-curve artifact.
    assert unlinked.availability_summary()["both_available"] == len(d_values)
    assert int(np.sum(unlinked.linking_numbers != 0)) == 0

    # The metric residual is finite/non-negative and larger when pulled apart.
    dmid = float(d_values[len(d_values) // 2])
    dist_linked = closest_curve_distance(
        grid_a,
        grid_b,
        scanning_component="xdot",
        curve_components=("x", "y", "z"),
        d=dmid,
    )
    dist_far = closest_curve_distance(
        grid_a,
        grid_far,
        scanning_component="xdot",
        curve_components=("x", "y", "z"),
        d=dmid,
    )
    assert dist_linked is not None and dist_far is not None
    assert 0.0 <= dist_linked < dist_far


def test_availability_summary_flags_missing_curves() -> None:
    """``availability_summary`` must distinguish a genuine ``linking == 0``
    (both curves extracted) from a scan where the level lies outside the grid's
    finite range so no curve exists (the ambiguity #553 flagged)."""
    grid_a = _synthetic_ring_grid(lambda t: (math.cos(t), math.sin(t), 0.0))
    grid_b = _synthetic_ring_grid(lambda t: (0.8 + 0.6 * math.cos(t), 0.0, 0.6 * math.sin(t)))
    # The scanning field is cos(theta_long) in [-1, 1]; levels far outside that
    # range extract NO curve on either grid -> both_available must be 0.
    out_of_range = np.linspace(5.0, 9.0, 12)
    res = scan_linking_number(
        grid_a,
        grid_b,
        scanning_component="xdot",
        curve_components=("x", "y", "z"),
        d_values=out_of_range,
    )
    avail = res.availability_summary()
    assert avail["n"] == len(out_of_range)
    assert avail["both_available"] == 0
    assert avail["neither_or_one"] == len(out_of_range)
    # Every linking number is the "no data" 0 -- and now provably so.
    assert int(np.sum(res.linking_numbers != 0)) == 0
