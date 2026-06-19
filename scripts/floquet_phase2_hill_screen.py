"""#391 P391.2 — retroactive Hill-fraction screen of the Phase 2 sweep branches.

Reads the #347 Phase 2 family-network nodes
(``data/floquet_phase2_family_network.jsonl``) — every parent + branch — and
runs the :mod:`cyclerfinder.genome.hill_screen` amplitude-vs-Hill-fraction
pre-screen on each, emitting one row per node to
``data/floquet_phase2_hill_screen.jsonl``.

This is a *screening* tool, not a catalogue admission: no catalogue writeback.
It answers the load-bearing #391 question — is ``branch_C11a_b0`` (the
stability-flip (1,1) orbit) PASS (a next #389-style gauntlet candidate) or
V4_DOOMED?

The family-network file is the canonical, cross-checked node list (the sweep
results file's ``branch_record`` rows carry the same states); using the network
file means every node (parents included) gets a Hill classification for context.

Usage:
    uv run python scripts/floquet_phase2_hill_screen.py [--input PATH] [--output PATH]
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path
from typing import Any

import numpy as np

from cyclerfinder.genome.hill_screen import screen_orbit
from cyclerfinder.search.reachable_representatives import TU_DAYS, braik_ross_system

DEFAULT_INPUT = Path("data/floquet_phase2_family_network.jsonl")
DEFAULT_OUTPUT = Path("data/floquet_phase2_hill_screen.jsonl")


def _iter_orbit_nodes(input_path: Path) -> list[dict[str, Any]]:
    """Return the parent + branch node rows (those carrying a state0 + period)."""
    nodes: list[dict[str, Any]] = []
    with input_path.open() as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            row: dict[str, Any] = json.loads(line)
            if row.get("kind") == "node" and "state0" in row and "period_TU" in row:
                nodes.append(row)
    return nodes


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    system = braik_ross_system()
    nodes = _iter_orbit_nodes(args.input)
    iso_start = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w") as fh:
        fh.write(
            json.dumps(
                {
                    "kind": "header",
                    "phase": "391_p391_2",
                    "iso_timestamp": iso_start,
                    "input_artifact": str(args.input),
                    "n_nodes": len(nodes),
                    "l_km": float(system.l_km),
                    "screen_bands": {"PASS": "<0.3", "MARGINAL": "0.3-0.5", "V4_DOOMED": ">0.5"},
                    "discipline": (
                        "#391 retroactive Hill-fraction pre-screen of #347 Phase 2 nodes. "
                        "Threshold rationale: #389 branch_C32_b0 failed V4 at 0.77 R_Hill "
                        "(solar tide ~30% of Earth gravity). SCREENING only — NO catalogue "
                        "writeback. Constants sourced via cyclerfinder.core.constants."
                    ),
                }
            )
            + "\n"
        )
        counts = {"PASS": 0, "MARGINAL": 0, "V4_DOOMED": 0}
        for node in nodes:
            state0 = np.array(node["state0"], dtype=np.float64)
            period = float(node["period_TU"])
            result = screen_orbit(system, state0, period)
            counts[result.classification] += 1
            fh.write(
                json.dumps(
                    {
                        "kind": "hill_screen",
                        "node_id": node.get("node_id"),
                        "label": node.get("label"),
                        "node_type": node.get("node_type"),
                        "topology": [node.get("k1"), node.get("k2")],
                        "jacobi": node.get("jacobi"),
                        "period_TU": period,
                        "period_days": period * TU_DAYS,
                        "max_amplitude_km": result.max_amplitude_km,
                        "earth_sun_hill_radius_km": result.earth_sun_hill_radius_km,
                        "hill_fraction": result.hill_fraction,
                        "solar_tide_to_earth_gravity_ratio": (
                            result.solar_tide_to_earth_gravity_ratio
                        ),
                        "classification": result.classification,
                    }
                )
                + "\n"
            )
            print(
                f"{node.get('label'):16s} type={node.get('node_type'):6s} "
                f"amp={result.max_amplitude_km:11.1f} km  frac={result.hill_fraction:.4f}  "
                f"-> {result.classification}"
            )
        fh.write(
            json.dumps(
                {
                    "kind": "footer",
                    "phase": "391_p391_2",
                    "iso_end": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                    "counts": counts,
                }
            )
            + "\n"
        )
    print(f"\nwrote {args.output} ({counts})")


if __name__ == "__main__":
    main()
