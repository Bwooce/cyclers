"""Task #633: fixed-Saturn-Titan-mu 2D (x0, C, hc) corridor grid search tests.

Covers the genuinely new, reusable pieces of code
``scripts/run_633_titan_corridor_grid.py`` adds:

1. :func:`build_grids`/:func:`flatten_points` -- the grid construction, which
   must match the #629 design read's own spec exactly (x0 count, rho count
   INCLUDING the fine [0.955, 0.995] sub-band, hc set, and the ~16,375 total
   point count), and produce a deterministic, index-addressable flattening
   (required for the chunked --start-idx/--end-idx resumability this task's
   own incremental-progress discipline depends on).
2. :func:`evaluate_point` -- the per-grid-point gate chain. A fast
   regression reproduces one of this task's own timing-pilot points exactly
   (self-consistency, not a sourced golden -- see module docstring
   discipline elsewhere in this repo for that distinction). A second test
   grounds :func:`_direct_gate` against the ALREADY-ADMITTED, catalogue-
   registered Pluto-Charon (3,2) stable cycler (``ross-rt-pc-cycler-32-2026``,
   V2) -- a genuinely sourced correctness check, not circular.
3. :func:`_append_jsonl`/:func:`summarize` -- the incremental-checkpointing
   round trip this task's own process discipline requires
   ([[feedback_incremental_progress_reports]]).
"""

from __future__ import annotations

import json
from typing import Any

import numpy as np
import pytest

import cyclerfinder.core.cr3bp as cr3bp
import cyclerfinder.search.cr3bp_periodic as cp
import scripts.run_633_titan_corridor_grid as run
from cyclerfinder.core.cr3bp import cr3bp_system
from cyclerfinder.search.pluto_charon_kk_sweep import _c_l1

# ---------------------------------------------------------------------------
# 1. Grid construction
# ---------------------------------------------------------------------------


def test_build_grids_matches_629_design_read_spec() -> None:
    x0_grid, rho_grid, hc_list = run.build_grids()

    assert x0_grid[0] == pytest.approx(-0.95)
    assert x0_grid[-1] == pytest.approx(-0.30)
    assert len(x0_grid) == 131  # (0.65 / 0.005) + 1

    assert rho_grid.min() == pytest.approx(0.60)
    assert rho_grid.max() == pytest.approx(0.995)
    # the fine sub-band must be present -- without it, thin-corridor families
    # like (3,3) (anchor rho=0.974) are missed by construction.
    assert np.any(rho_grid >= 0.955)
    assert np.any((rho_grid >= 0.95) & (rho_grid <= 0.995))
    assert len(rho_grid) == 25

    assert hc_list == (1, 3, 5, 7, 9)

    n_total = len(x0_grid) * len(rho_grid) * len(hc_list)
    assert n_total == 16375


def test_flatten_points_is_deterministic_and_index_addressable() -> None:
    x0_grid, rho_grid, hc_list = run.build_grids()
    points = run.flatten_points(x0_grid, rho_grid, hc_list)

    n_total = len(x0_grid) * len(rho_grid) * len(hc_list)
    assert len(points) == n_total

    # idx is sequential and matches the point's own position.
    for i, (idx, _x0, _rho, _hc) in enumerate(points[:50]):
        assert idx == i

    # hc is the outer loop variable (per the module's own documented order):
    # every point in the first len(rho_grid)*len(x0_grid) block shares hc[0].
    block = len(rho_grid) * len(x0_grid)
    assert all(p[3] == hc_list[0] for p in points[:block])
    assert points[block][3] == hc_list[1]

    # re-flattening is byte-for-byte identical (needed for chunk resumability:
    # a --start-idx/--end-idx chunk run days apart must address the same points).
    points2 = run.flatten_points(*run.build_grids())
    assert points == points2


# ---------------------------------------------------------------------------
# 2. evaluate_point -- self-consistency regression + sourced correctness check
# ---------------------------------------------------------------------------


@pytest.mark.slow
def test_evaluate_point_reproduces_timing_pilot_point() -> None:
    """Self-consistency regression: reproduces this task's own pilot output.

    x0=-0.94, rho=0.6, hc=1 at Saturn-Titan mu converged to a genuine (1,0)
    periodic orbit with reaches_secondary=False (observed directly during
    this task's timing pilot, see the #633 OUTSTANDING.md bullet). This is a
    regression on OUR solver's own reproducibility, not a physical claim
    sourced externally -- consistent with this project's existing
    self-consistency test convention (e.g. #627's own short-hop test).
    """
    system = cr3bp_system("Saturn", "Titan")
    c_l1 = _c_l1(system.mu)
    rec = run.evaluate_point(system, 2, -0.94, 0.6, 1, c_l1)

    assert rec["converged"] is True
    assert rec["k1"] == 1
    assert rec["k2"] == 0
    assert rec["reaches_secondary"] is False
    assert rec["target_match"] is None


def test_evaluate_point_negative_radicand_reported_cleanly() -> None:
    """x0=-0.95 (the grid's own edge) at a low rho has no real ydot0 -- must
    be reported as a clean non-convergence, not raise."""
    system = cr3bp_system("Saturn", "Titan")
    c_l1 = _c_l1(system.mu)
    rec = run.evaluate_point(system, 0, -0.95, 0.6, 1, c_l1)
    assert rec["converged"] is False
    assert "error" in rec


def test_direct_gate_recovers_admitted_pluto_charon_32_cycler() -> None:
    """Grounds :func:`_direct_gate` against the catalogue-admitted, V2
    Pluto-Charon (3,2) stable cycler (id ``ross-rt-pc-cycler-32-2026``):
    x0=-0.693198287043369, ydot0=-0.297004785528322,
    C=3.57951501972907, T=11.8334625170346 TU, mu=0.10876473603280369,
    Barden stability_index=3.82e-9 (maximally stable). A genuinely sourced
    (catalogue-registered) correctness check for the exact gate function the
    #633 grid uses to decide "stable + on-topology", not a value this task's
    own code computed circularly.
    """
    pc = cr3bp.CR3BPSystem(
        mu=0.10876473603280369, primary="Pluto", secondary="Charon", l_km=19600.0, t_s=87855.81
    )
    orbit = cp.SymmetricOrbit(
        x0=-0.693198287043369,
        ydot0=-0.297004785528322,
        jacobi=3.57951501972907,
        t_half=11.8334625170346 / 2.0,
        period=11.8334625170346,
        converged=True,
        crossing_residual=0.0,
        n_iter=0,
    )
    gate = run._direct_gate(pc, orbit, 3, 2)
    assert gate["topology"] == (3, 2)
    assert gate["reaches_secondary"] is True
    assert gate["nu_direct"] == pytest.approx(3.82e-9, abs=5e-9)
    assert gate["gate_stability_pass"] is True


# ---------------------------------------------------------------------------
# 3. Incremental checkpointing round trip
# ---------------------------------------------------------------------------


def test_append_jsonl_and_summarize_round_trip(tmp_path) -> None:  # type: ignore[no-untyped-def]
    out_path = tmp_path / "633_test_checkpoint.jsonl"
    recs: list[dict[str, Any]] = [
        {"idx": 0, "converged": True, "k1": 1, "k2": 0, "target_match": None},
        {"idx": 1, "converged": False},
        {
            "idx": 2,
            "converged": True,
            "k1": 1,
            "k2": 1,
            "target_match": "11",
            "gate_stability_pass": True,
            "gate_encounter_pass": False,
        },
    ]
    for r in recs:
        run._append_jsonl(out_path, r)

    lines = out_path.read_text().strip().split("\n")
    assert len(lines) == 3
    for line, expected in zip(lines, recs, strict=True):
        assert json.loads(line) == expected

    # summarize() must not raise on a partial checkpoint and must count the
    # candidate correctly (stdout content isn't asserted -- summarize() is a
    # human-facing report, not an API surface with a return value).
    run.summarize(out_path)


def test_summarize_handles_missing_checkpoint_file(tmp_path) -> None:  # type: ignore[no-untyped-def]
    missing = tmp_path / "does_not_exist.jsonl"
    run.summarize(missing)  # must not raise
