"""#347 Phase 2 P2.4 — family-network artifact builder (Doedel 2003 Fig. 4-5).

Reads the Phase 2 sweep results JSONL ``data/floquet_phase2_sweep_results.jsonl``
and emits a family-network JSONL ``data/floquet_phase2_family_network.jsonl``
with two record types:

  * **node** — one per parent OR branched orbit. Carries the orbit's IC + period +
    topology + Floquet diagnostic. Parents are nodes labelled ``parent_<label>``;
    branched orbits are nodes labelled ``branch_<parent_label>_<bracket_index>``.
  * **edge** — one per (parent, branched_orbit) tuple. Carries the bifurcation
    type (k=1 saddle-center) + the bifurcation parameter (C* bracket) + the
    eigenvalue used.

This is the post-sweep catalogue artifact format recommended in the Doedel 2003
digest (Fig. 4 p.1360 + Fig. 5 p.1361): each family is a node, each bifurcation
is a labelled edge between two nodes. The graph is directed: parent → branch.

Per Phase 0 design doc Section 6 + the user's #347 brief: this is a
DISCOVERY artifact, not a catalogue admission. Each node has a
``catalogue_status: "phase2_discovery_candidate"`` marker; admission to the
project catalogue is Phase 4+ gauntlet work.

Usage:
    uv run python scripts/floquet_phase2_family_network.py \\
        [--input PATH] [--output PATH]
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path
from typing import Any


def _build_node(record: dict[str, Any], kind: str) -> dict[str, Any]:
    """Build a node record. ``kind`` is one of {'parent', 'branch'}."""
    if kind == "parent":
        node_id = f"parent_{record['parent_label']}"
        return {
            "kind": "node",
            "node_id": node_id,
            "node_type": "parent",
            "label": record["parent_label"],
            "k1": int(record["parent_k1"]),
            "k2": int(record["parent_k2"]),
            "jacobi": float(record["parent_jacobi_anchor"]),
            "state0": list(record["parent_state0"]),
            "period_TU": float(record["parent_period_TU"]),
            "period_days": float(record["parent_period_days"]),
            "topology_verified": list(record["parent_topology"]),
            "catalogue_status": "phase2_discovery_candidate",
            "source": "braik_ross_2026_table2",
        }
    elif kind == "branch":
        node_id = f"branch_{record['parent_label']}_b{record['bracket_index']}"
        return {
            "kind": "node",
            "node_id": node_id,
            "node_type": "branch",
            "label": node_id,
            "k1": int(record["branch_k1"]),
            "k2": int(record["branch_k2"]),
            "jacobi": float(record["branch_jacobi"]),
            "state0": list(record["branch_state0"]),
            "period_TU": float(record["branch_period_TU"]),
            "period_days": float(record["branch_period_days"]),
            "z0": float(record["branch_z0"]),
            "zdot0": float(record["branch_zdot0"]),
            "degenerate_planar": bool(record["branch_degenerate_planar"]),
            "max_floquet_mag": float(record["max_floquet_mag_branched"]),
            "sigma_TU": float(record["sigma_TU_branched"]),
            "sigma_d_per_day": float(record["sigma_d_per_day_branched"]),
            "corrector_residual": float(record["corrector_residual"]),
            "independent_closure_residual": float(record["independent_closure_residual"]),
            "topology_changed_from_parent": bool(record["topology_changed"]),
            "topology_recheck": [
                int(record["branch_k1_recheck"]),
                int(record["branch_k2_recheck"]),
            ],
            "cycler_candidate_flag": bool(record["cycler_candidate_flag"]),
            "catalogue_status": "phase2_discovery_candidate",
        }
    else:
        raise ValueError(f"unknown node kind {kind}")


def _build_edge(record: dict[str, Any]) -> dict[str, Any]:
    """Build a directed edge record from parent → branch."""
    parent_node = f"parent_{record['parent_label']}"
    branch_node = f"branch_{record['parent_label']}_b{record['bracket_index']}"
    return {
        "kind": "edge",
        "source_node": parent_node,
        "target_node": branch_node,
        "edge_type": "saddle_center_bifurcation",
        "bifurcation_k": 1,  # k=1: saddle-center / pitchfork at λ=+1
        "bifurcation_param_before": record.get("bracket_param_before"),
        "bifurcation_param_after": record.get("bracket_param_after"),
        "eigenvalue_real": float(record["eigenvalue_used_real"]),
        "eigenvalue_imag": float(record["eigenvalue_used_imag"]),
        "eigenvector": list(record["eigenvector_used"]),
        "epsilon_used": float(record["epsilon_used"]),
        "eigenvector_sign_used": int(record["eigenvector_sign_used"]),
        "branch_corrector_residual": float(record["corrector_residual"]),
    }


def build_family_network(input_path: Path, output_path: Path) -> dict[str, int]:
    """Process the sweep results and emit a family-network JSONL.

    Returns a counters dict: nodes_emitted, edges_emitted, parents_seen,
    branches_seen.
    """
    counters = {
        "nodes_emitted": 0,
        "edges_emitted": 0,
        "parents_seen": 0,
        "branches_seen": 0,
        "branch_no_converge_rows": 0,
        "error_rows": 0,
    }
    parent_nodes_emitted: set[str] = set()
    with input_path.open() as fin:
        rows = [json.loads(line) for line in fin]

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w") as fout:
        # Header
        iso = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
        fout.write(
            json.dumps(
                {
                    "kind": "header",
                    "phase": "347_phase2_p2_4",
                    "iso_timestamp": iso,
                    "input_artifact": str(input_path),
                    "format": "doedel_2003_family_network",
                    "discipline": (
                        "Phase 2 DISCOVERY candidates only; "
                        "catalogue_status=phase2_discovery_candidate on every node. "
                        "Catalogue admission is Phase 4+ gauntlet work."
                    ),
                }
            )
            + "\n"
        )
        # Walk sweep rows; emit nodes + edges.
        for row in rows:
            if row.get("kind") == "branch_record":
                parent_id = f"parent_{row['parent_label']}"
                if parent_id not in parent_nodes_emitted:
                    fout.write(json.dumps(_build_node(row, "parent")) + "\n")
                    parent_nodes_emitted.add(parent_id)
                    counters["parents_seen"] += 1
                    counters["nodes_emitted"] += 1
                fout.write(json.dumps(_build_node(row, "branch")) + "\n")
                counters["branches_seen"] += 1
                counters["nodes_emitted"] += 1
                fout.write(json.dumps(_build_edge(row)) + "\n")
                counters["edges_emitted"] += 1
            elif row.get("kind") == "branch_no_converge":
                counters["branch_no_converge_rows"] += 1
            elif row.get("kind") in ("parent_error", "branch_error"):
                counters["error_rows"] += 1
        # Footer
        fout.write(
            json.dumps(
                {
                    "kind": "footer",
                    "phase": "347_phase2_p2_4",
                    "iso_end": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                    "counters": counters,
                }
            )
            + "\n"
        )
    return counters


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("data/floquet_phase2_sweep_results.jsonl"),
        help="input sweep-results JSONL path",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("data/floquet_phase2_family_network.jsonl"),
        help="output family-network JSONL path",
    )
    args = parser.parse_args()
    if not args.input.exists():
        raise SystemExit(f"input artifact does not exist: {args.input}")
    counters = build_family_network(args.input, args.output)
    print(
        f"family-network built: {counters['nodes_emitted']} nodes "
        f"(parents={counters['parents_seen']}, branches={counters['branches_seen']}), "
        f"{counters['edges_emitted']} edges, "
        f"{counters['branch_no_converge_rows']} no-converge rows, "
        f"{counters['error_rows']} error rows"
    )
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
