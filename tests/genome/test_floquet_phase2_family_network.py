"""#347 Phase 2 P2.4 — family-network builder tests.

Tests the helper :func:`build_family_network` that converts the Phase 2 sweep
JSONL into a Doedel 2003 Fig. 4-5-style family-network JSONL (nodes + edges).

Gates:

  1. A synthetic sweep input with 1 parent + 1 branch produces:
     1 header + 1 parent node + 1 branch node + 1 edge + 1 footer = 5 rows.
  2. Two branches from the same parent produce: 1 header + 1 parent node +
     2 branch nodes + 2 edges + 1 footer = 7 rows; parent node deduplicated.
  3. branch_no_converge rows are COUNTED but produce no node/edge.
  4. parent_error / branch_error rows are COUNTED but produce no node/edge.
  5. Every node has ``catalogue_status: 'phase2_discovery_candidate'``.
  6. Edge fields trace to the source sweep row (parent_label, eigenvalue,
    eigenvector, residual).
"""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

# Load the builder module from scripts/.
_REPO_ROOT = Path(__file__).resolve().parent.parent.parent
_BUILDER_PATH = _REPO_ROOT / "scripts" / "floquet_phase2_family_network.py"
_spec = importlib.util.spec_from_file_location("phase2_family_network_builder", _BUILDER_PATH)
assert _spec is not None and _spec.loader is not None
_builder = importlib.util.module_from_spec(_spec)
sys.modules["phase2_family_network_builder"] = _builder
_spec.loader.exec_module(_builder)


def _make_branch_record(
    parent_label: str,
    bracket_index: int,
    parent_k: tuple[int, int] = (3, 2),
    branch_k: tuple[int, int] = (3, 3),
) -> dict:  # type: ignore[type-arg]
    """Build a synthetic branch_record matching the sweep driver's schema."""
    return {
        "kind": "branch_record",
        "parent_label": parent_label,
        "parent_k1": parent_k[0],
        "parent_k2": parent_k[1],
        "parent_jacobi_anchor": 3.1294,
        "parent_state0": [-0.284, 0.0, 0.0, 0.0, -2.053, 0.0],
        "parent_period_TU": 17.46,
        "parent_period_days": 75.94,
        "parent_topology": list(parent_k),
        "bracket_index": bracket_index,
        "bracket_param_before": 3.14170,
        "bracket_param_after": 3.14180,
        "bracket_eig_before_real": 0.999,
        "bracket_eig_before_imag": 0.047,
        "bracket_eig_after_real": 1.026,
        "bracket_eig_after_imag": 0.0,
        "branch_state0": [-0.70, -2.91, 1e-22, -2.35, 0.57, 1e-23],
        "branch_period_TU": 23.36,
        "branch_period_days": 101.56,
        "branch_jacobi": 3.797,
        "branch_z0": 1e-22,
        "branch_zdot0": 1e-23,
        "branch_degenerate_planar": True,
        "branch_k1": branch_k[0],
        "branch_k2": branch_k[1],
        "branch_k1_recheck": branch_k[0],
        "branch_k2_recheck": branch_k[1],
        "topology_changed": parent_k != branch_k,
        "corrector_residual": 4.77e-12,
        "independent_closure_residual": 2.59e-11,
        "max_floquet_mag_branched": 1.0,
        "sigma_TU_branched": 2.6e-14,
        "sigma_d_per_day_branched": 6e-15,
        "epsilon_used": 5e-4,
        "eigenvector_sign_used": -1,
        "eigenvalue_used_real": 0.974,
        "eigenvalue_used_imag": 0.0,
        "eigenvector_used": [
            -2.9e-10,
            9.2e-12,
            -0.842,
            2.1e-12,
            1.9e-09,
            0.539,
        ],
        "cycler_candidate_flag": True,
    }


def _write_jsonl(path: Path, rows: list[dict]) -> None:  # type: ignore[type-arg]
    with path.open("w") as fh:
        for row in rows:
            fh.write(json.dumps(row) + "\n")


def test_single_parent_single_branch(tmp_path: Path) -> None:
    """1 parent + 1 branch -> 5 rows (header + parent + branch + edge + footer)."""
    input_path = tmp_path / "sweep.jsonl"
    output_path = tmp_path / "network.jsonl"
    _write_jsonl(input_path, [_make_branch_record("C32", 0)])
    counters = _builder.build_family_network(input_path, output_path)
    assert counters["parents_seen"] == 1
    assert counters["branches_seen"] == 1
    assert counters["nodes_emitted"] == 2
    assert counters["edges_emitted"] == 1
    with output_path.open() as fh:
        rows = [json.loads(line) for line in fh]
    assert len(rows) == 5
    assert rows[0]["kind"] == "header"
    assert rows[1]["kind"] == "node"
    assert rows[1]["node_type"] == "parent"
    assert rows[1]["label"] == "C32"
    assert rows[2]["kind"] == "node"
    assert rows[2]["node_type"] == "branch"
    assert rows[2]["k1"] == 3 and rows[2]["k2"] == 3
    assert rows[3]["kind"] == "edge"
    assert rows[3]["source_node"] == "parent_C32"
    assert rows[3]["target_node"] == "branch_C32_b0"
    assert rows[3]["edge_type"] == "saddle_center_bifurcation"
    assert rows[3]["bifurcation_k"] == 1
    assert rows[4]["kind"] == "footer"


def test_two_branches_same_parent_dedups_parent_node(tmp_path: Path) -> None:
    """A parent emitting two branches has its node emitted ONCE; two edges follow."""
    input_path = tmp_path / "sweep.jsonl"
    output_path = tmp_path / "network.jsonl"
    _write_jsonl(
        input_path,
        [
            _make_branch_record("C32", 0, branch_k=(3, 3)),
            _make_branch_record("C32", 1, branch_k=(2, 0)),
        ],
    )
    counters = _builder.build_family_network(input_path, output_path)
    assert counters["parents_seen"] == 1
    assert counters["branches_seen"] == 2
    assert counters["nodes_emitted"] == 3
    assert counters["edges_emitted"] == 2
    with output_path.open() as fh:
        rows = [json.loads(line) for line in fh]
    # 1 header + 1 parent + 1 branch + 1 edge + 1 branch + 1 edge + 1 footer = 7.
    assert len(rows) == 7


def test_branch_no_converge_does_not_emit_node_or_edge(tmp_path: Path) -> None:
    """A branch_no_converge row is counted but emits no node/edge."""
    input_path = tmp_path / "sweep.jsonl"
    output_path = tmp_path / "network.jsonl"
    _write_jsonl(
        input_path,
        [
            {
                "kind": "branch_no_converge",
                "parent_label": "C21",
                "parent_k1": 2,
                "parent_k2": 1,
                "bracket_index": 0,
                "epsilon": 5e-4,
            },
        ],
    )
    counters = _builder.build_family_network(input_path, output_path)
    assert counters["nodes_emitted"] == 0
    assert counters["edges_emitted"] == 0
    assert counters["branch_no_converge_rows"] == 1


def test_parent_error_does_not_emit_node_or_edge(tmp_path: Path) -> None:
    """A parent_error row is counted but emits no node/edge."""
    input_path = tmp_path / "sweep.jsonl"
    output_path = tmp_path / "network.jsonl"
    _write_jsonl(
        input_path,
        [
            {
                "kind": "parent_error",
                "parent_label": "BAD",
                "stage": "anchor_recovery",
                "error": "synthetic",
            }
        ],
    )
    counters = _builder.build_family_network(input_path, output_path)
    assert counters["nodes_emitted"] == 0
    assert counters["edges_emitted"] == 0
    assert counters["error_rows"] == 1


def test_every_node_has_phase2_discovery_status(tmp_path: Path) -> None:
    """Every emitted node carries catalogue_status='phase2_discovery_candidate'."""
    input_path = tmp_path / "sweep.jsonl"
    output_path = tmp_path / "network.jsonl"
    _write_jsonl(
        input_path,
        [
            _make_branch_record("C32", 0),
            _make_branch_record("C11a", 0),
        ],
    )
    _builder.build_family_network(input_path, output_path)
    with output_path.open() as fh:
        rows = [json.loads(line) for line in fh]
    nodes = [r for r in rows if r["kind"] == "node"]
    assert len(nodes) == 4  # 2 parents + 2 branches
    for n in nodes:
        assert n["catalogue_status"] == "phase2_discovery_candidate"


def test_edge_fields_trace_to_source_row(tmp_path: Path) -> None:
    """Edge eigenvalue/eigenvector/residual come from the source branch_record."""
    input_path = tmp_path / "sweep.jsonl"
    output_path = tmp_path / "network.jsonl"
    src = _make_branch_record("C32", 0)
    _write_jsonl(input_path, [src])
    _builder.build_family_network(input_path, output_path)
    with output_path.open() as fh:
        rows = [json.loads(line) for line in fh]
    edge = next(r for r in rows if r["kind"] == "edge")
    assert edge["eigenvalue_real"] == src["eigenvalue_used_real"]
    assert edge["eigenvalue_imag"] == src["eigenvalue_used_imag"]
    assert edge["eigenvector"] == src["eigenvector_used"]
    assert edge["branch_corrector_residual"] == src["corrector_residual"]
    assert edge["epsilon_used"] == src["epsilon_used"]


def test_empty_sweep_produces_header_footer_only(tmp_path: Path) -> None:
    """A sweep with no branch_records produces just header + footer (no nodes)."""
    input_path = tmp_path / "sweep.jsonl"
    output_path = tmp_path / "network.jsonl"
    _write_jsonl(input_path, [])
    counters = _builder.build_family_network(input_path, output_path)
    assert counters["nodes_emitted"] == 0
    assert counters["edges_emitted"] == 0
    with output_path.open() as fh:
        rows = [json.loads(line) for line in fh]
    assert len(rows) == 2
    assert rows[0]["kind"] == "header"
    assert rows[1]["kind"] == "footer"
