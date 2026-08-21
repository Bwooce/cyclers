"""Tests for the Resonant Atlas pilot Stage A harness (`#859`).

Covers the harness's own logic (`#859`'s own dispatch note: "not slow, not
requiring the full sweep to have run"):

* ``coprime_pairs`` enumeration (correctness, determinism, no duplicates);
* ``is_in_band`` boundary behavior (the |lambda| in [50, 2500] tractable
  band, `#858` Sec. 3.2/7);
* ``build_stage_a_cells`` grid structure (deterministic ordering — the
  ``campaign_runner`` resumability contract depends on it);
* ``stage_a_worker``'s two DETERMINISTIC fast-path outcomes (seed-domain-
  invalid, seed-did-not-converge) -- no real continuation, near-instant;
* ONE real, cheap end-to-end cell (Uranus-Oberon 3:2, ``n_c_steps=1``) to
  confirm the corrector/continuation/classification wiring actually
  produces a sane, non-garbage result (not just "did not crash") -- kept to
  a single continuation step to stay fast (a few seconds, not the smoke
  test's own minutes-long slices);
* checkpoint read/write + resume through the REAL worker over a tiny grid
  (mirrors ``test_campaign_runner.py``'s own resume pattern, but exercises
  this module's own worker rather than a synthetic one).
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import pytest

from cyclerfinder.search.campaign_runner import (
    CampaignRunnerConfig,
    CampaignRunnerRouting,
    run_grid_campaign,
)
from cyclerfinder.search.resonant_atlas_stage_a import (
    STAGE_A_SYSTEMS,
    build_stage_a_cells,
    coprime_pairs,
    is_in_band,
    stage_a_worker,
)

# ---------------------------------------------------------------------------
# coprime_pairs
# ---------------------------------------------------------------------------


def test_coprime_pairs_max_3_matches_hand_count() -> None:
    pairs = coprime_pairs(3)
    assert pairs == [(1, 1), (1, 2), (1, 3), (2, 1), (2, 3), (3, 1), (3, 2)]


def test_coprime_pairs_every_pair_is_actually_coprime() -> None:
    pairs = coprime_pairs(8)
    assert all(math.gcd(p, q) == 1 for p, q in pairs)


def test_coprime_pairs_excludes_non_coprime() -> None:
    pairs = coprime_pairs(8)
    assert (2, 2) not in pairs
    assert (4, 6) not in pairs
    assert (3, 6) not in pairs


def test_coprime_pairs_no_duplicates_and_deterministic() -> None:
    pairs = coprime_pairs(8)
    assert len(pairs) == len(set(pairs))
    assert pairs == coprime_pairs(8)  # same call twice -> identical order


def test_coprime_pairs_count_p_q_le_8() -> None:
    # Hand-cross-checked count (#859's own registration: "all coprime p:q
    # with p,q<=8"): 43 ordered pairs.
    assert len(coprime_pairs(8)) == 43


def test_coprime_pairs_rejects_bad_max_pq() -> None:
    with pytest.raises(ValueError):
        coprime_pairs(0)


# ---------------------------------------------------------------------------
# is_in_band
# ---------------------------------------------------------------------------


def test_is_in_band_boundaries() -> None:
    assert is_in_band(50.0) is True
    assert is_in_band(2500.0) is True
    assert is_in_band(49.999999) is False
    assert is_in_band(2500.000001) is False


def test_is_in_band_interior_and_exterior() -> None:
    assert is_in_band(105.05) is True  # #781's own 4:5-saddle Neptune-Triton eigenvalue
    assert is_in_band(1.0) is False  # near-unit-circle, never tractable
    assert is_in_band(14600.0) is False  # #781's own 4:7-stress FAIL, lambda~=1.46e4


# ---------------------------------------------------------------------------
# build_stage_a_cells
# ---------------------------------------------------------------------------


def test_build_stage_a_cells_grid_size_matches_systems_times_pairs() -> None:
    cells = build_stage_a_cells(STAGE_A_SYSTEMS, max_pq=5)
    assert len(cells) == len(STAGE_A_SYSTEMS) * len(coprime_pairs(5))


def test_build_stage_a_cells_ordering_is_systems_outer_pairs_inner() -> None:
    cells = build_stage_a_cells(STAGE_A_SYSTEMS[:2], max_pq=3)
    n_pairs = len(coprime_pairs(3))
    assert len(cells) == 2 * n_pairs
    assert all(c["system_key"] == STAGE_A_SYSTEMS[0].system_key for c in cells[:n_pairs])
    assert all(c["system_key"] == STAGE_A_SYSTEMS[1].system_key for c in cells[n_pairs:])


def test_build_stage_a_cells_fields_and_defaults() -> None:
    cells = build_stage_a_cells((STAGE_A_SYSTEMS[0],), max_pq=2, n_c_steps=4, d_jacobi=1e-3)
    assert len(cells) == len(coprime_pairs(2))
    for c in cells:
        assert c["primary"] == STAGE_A_SYSTEMS[0].primary
        assert c["secondary"] == STAGE_A_SYSTEMS[0].secondary
        assert c["role"] == STAGE_A_SYSTEMS[0].role
        assert c["x0_sign"] == -1
        assert c["n_c_steps"] == 4
        assert c["d_jacobi"] == 1e-3
        assert math.gcd(c["p"], c["q"]) == 1


def test_build_stage_a_cells_are_json_serializable() -> None:
    """The campaign-runner/parallel_sweep pickle-safety contract, plus this
    project's own JSONL-checkpoint round-trip -- cells must survive both."""
    cells = build_stage_a_cells((STAGE_A_SYSTEMS[0],), max_pq=2)
    for c in cells:
        assert json.loads(json.dumps(c)) == c


# ---------------------------------------------------------------------------
# stage_a_worker: deterministic fast paths (no real continuation)
# ---------------------------------------------------------------------------


def _cell(
    *, system_key: str, primary: str, secondary: str, p: int, q: int, n_c_steps: int = 1
) -> dict[str, object]:
    return {
        "system_key": system_key,
        "primary": primary,
        "secondary": secondary,
        "role": "novel_target",
        "p": p,
        "q": q,
        "x0_sign": -1,
        "n_c_steps": n_c_steps,
        "d_jacobi": 5e-4,
    }


def test_stage_a_worker_seed_domain_invalid_is_a_fast_clean_miss() -> None:
    # (3, 1): a=(1/3)**(2/3)~=0.4807 < 0.5 -> two_body_resonant_seed's own
    # vis-viva sqrt goes negative (module docstring's own documented,
    # deterministic construction limit). Observed this task: instant, no
    # integration attempted.
    cell = _cell(system_key="uranus-oberon", primary="Uranus", secondary="Oberon", p=3, q=1)
    out = stage_a_worker(cell)
    assert out.status == "miss"
    assert out.error == ""
    assert out.payload["reason"] == "seed_domain_invalid"
    assert out.payload["n_members"] == 0
    assert out.payload["n_in_band"] == 0


def test_stage_a_worker_seed_no_converge_is_a_clean_miss_not_an_error() -> None:
    # (1, 1): observed this task to fail the corrector outright (the
    # trivial 1:1 two-body seed sits exactly at the secondary's own
    # synodic period with no room for a perpendicular re-crossing at this
    # module's own seed convention). A real, reportable negative -- not a
    # worker malfunction.
    cell = _cell(system_key="uranus-oberon", primary="Uranus", secondary="Oberon", p=1, q=1)
    out = stage_a_worker(cell)
    assert out.status == "miss"
    assert out.error == ""
    assert out.payload["reason"] == "seed_did_not_converge"
    assert out.payload["n_members"] == 0


def test_stage_a_worker_domain_invalid_covers_every_low_q_over_p_ratio() -> None:
    """Every p:q this task found domain-invalid for p,q<=8 (module
    docstring's own worked cutoff, a=(q/p)**(2/3) < 0.5) is handled as a
    clean miss, not an exception -- across all p, not just one hand-picked
    case."""
    invalid_pairs = [(p, q) for p, q in coprime_pairs(8) if (q / p) ** (2.0 / 3.0) < 0.5]
    assert invalid_pairs  # sanity: the fixture premise is non-vacuous
    for p, q in invalid_pairs:
        cell = _cell(system_key="uranus-oberon", primary="Uranus", secondary="Oberon", p=p, q=q)
        out = stage_a_worker(cell)
        assert out.status == "miss"
        assert out.payload["reason"] == "seed_domain_invalid"


# ---------------------------------------------------------------------------
# stage_a_worker: one real, cheap end-to-end cell
# ---------------------------------------------------------------------------


def test_stage_a_worker_real_cell_recovers_a_sane_family() -> None:
    """Uranus-Oberon 3:2, n_c_steps=1 (seed + one continuation step): the
    positive-control system, the cheapest real-converging ratio observed
    this task (~a few seconds). Confirms the corrector -> half-crossing
    detection -> continue_family -> classification wiring produces a
    genuinely converged, physically sane family member, not just "did not
    raise"."""
    cell = _cell(system_key="uranus-oberon", primary="Uranus", secondary="Oberon", p=3, q=2)
    out = stage_a_worker(cell)
    assert out.error == ""
    assert out.status in ("hit", "miss")  # both are legitimate outcomes
    assert out.payload["n_members"] >= 1  # the seed member itself must have converged
    for m in out.payload["members"]:
        # A converged planar CR3BP orbit's own Barden |lambda| is >= 1 by
        # construction (nu = 1/2(lambda+1/lambda), |nu|>=1 <=> |lambda|>=1
        # for a real eigenvalue pair reciprocal about 1); a value far off
        # that is a sign the classification wiring itself is broken, not a
        # physically real family member.
        assert m["abs_lambda"] >= 1.0 - 1e-6
        assert m["period"] > 0.0


# ---------------------------------------------------------------------------
# Checkpoint read/write + resume through the REAL worker
# ---------------------------------------------------------------------------


def test_stage_a_checkpoint_resume_through_real_worker(tmp_path: Path) -> None:
    cells = build_stage_a_cells((STAGE_A_SYSTEMS[0],), max_pq=2, n_c_steps=1)
    routing = CampaignRunnerRouting(results_path=tmp_path / "results.jsonl")

    # First invocation: cap at 1 batch of 2 cells (partial).
    config_partial = CampaignRunnerConfig(n_workers=1, checkpoint_batch_size=2, max_batches=1)
    stats1 = run_grid_campaign(cells, stage_a_worker, routing=routing, config=config_partial)
    assert stats1.batches_run_this_invocation == 1
    assert stats1.evaluated_total < len(cells)
    lines_after_first = routing.results_path.read_text(encoding="utf-8").splitlines()
    assert len(lines_after_first) == stats1.evaluated_total

    # Second invocation, SAME routing: resumes and finishes the grid.
    config_full = CampaignRunnerConfig(n_workers=1, checkpoint_batch_size=2)
    stats2 = run_grid_campaign(cells, stage_a_worker, routing=routing, config=config_full)
    assert stats2.evaluated_total == len(cells)
    assert stats2.errors == 0

    # Every durable line parses and round-trips this module's own payload schema.
    final_lines = routing.results_path.read_text(encoding="utf-8").splitlines()
    assert len(final_lines) == len(cells)
    seen_indices = set()
    for line in final_lines:
        row = json.loads(line)
        assert row["index"] not in seen_indices  # no duplicate/overwritten cells
        seen_indices.add(row["index"])
        assert row["status"] in ("hit", "miss", "error")
    assert seen_indices == set(range(len(cells)))
