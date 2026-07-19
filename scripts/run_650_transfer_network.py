#!/usr/bin/env python3
"""#650: inter-cycler transfer-compatibility network -- the full pairwise sweep.

Reads ``data/catalogue.yaml``, applies the eligibility rules of
``docs/notes/2026-07-19-650-transfer-network-design.md`` §2 (``orbit_class``
in ``{cycler, quasi_cycler}`` with >=1 usable ``vinf_kms_at_encounters``
entry), enumerates every candidate node pair sharing >=1 usable encounter
body (§2), and for each ``(id_a, id_b, body)`` triple computes the full §3-§5
edge record via ``cyclerfinder.data.transfer_network.compute_edge`` -- the
same-body powered-flyby ``dv_hop_kms`` cost + band, and (gated per §5) the
statistical phase-alignment model over an explicitly-unknown relative phase.

Pure closed-form arithmetic over already-catalogued data -- no integrators,
no network access, no catalogue writeback. Writes the §8 artifact:

* ``data/found/650_transfer_network/edges.jsonl`` -- one record per edge.
* ``data/found/650_transfer_network/summary.json`` -- census: node/pair/edge
  counts, per-band/per-body histograms, cheap-edge count + list,
  ``cheap_dv_phase_indeterminate`` count, connected components + hubs of the
  cheap-edge subgraph, ``phase_status`` counts, data-gap tallies, and the
  parameter conventions echoed verbatim.

Runtime: ~32k candidate pairs, closed-form + phase model on the
dv/moon-gated subset -- minutes (design §9's own estimate; confirmed ~10s
for the edge sweep alone on this developer's machine). No background run,
no checkpointing.

No literature-novelty check, no catalogue writeback -- this is a derived
analysis artifact (the ``#317`` ``data/found/<task>_.../`` precedent), not a
discovery result. Whether a "no cheap edges anywhere" or any other outcome
should be registered in ``data/negative_results.yaml`` is an adjudication
call for the coordinating session, explicitly NOT this script's job (design
§8).
"""

from __future__ import annotations

import json
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

import yaml  # type: ignore[import-untyped]

from cyclerfinder.data.catalog import CATALOGUE_PATH
from cyclerfinder.data.transfer_network import (
    B0_DELTA_VINF_KMS,
    B1_DV_HOP_KMS,
    B2_DV_HOP_KMS,
    CHEAP_EDGE_DUTY_ADJUSTED_P_ALIGN,
    CHEAP_EDGE_WAIT_YEARS,
    DV_PHASE_GATE_KMS,
    N_DELTA0_SAMPLES,
    PHASE_LOCK_P_ALIGN_THRESHOLD,
    SENSITIVITY_WINDOWS_DAYS,
    STATISTICAL_HORIZON_YEARS,
    Edge,
    candidate_pairs,
    compute_edge,
    is_node,
    period_days,
    usable_bodies,
)

_OUT_DIR = Path("data/found/650_transfer_network")


def _load_catalogue_rows() -> list[dict[str, Any]]:
    with open(CATALOGUE_PATH) as fh:
        data = yaml.safe_load(fh)
    assert isinstance(data, list)
    return data


def _connected_components(edges: list[Edge]) -> list[list[str]]:
    """Connected components of the cheap-edge subgraph (union-find)."""
    parent: dict[str, str] = {}

    def find(x: str) -> str:
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a: str, b: str) -> None:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    for e in edges:
        union(e.id_a, e.id_b)

    groups: dict[str, list[str]] = defaultdict(list)
    for node_id in parent:
        groups[find(node_id)].append(node_id)
    return sorted((sorted(members) for members in groups.values()), key=lambda g: (-len(g), g))


def _degree_ranked_hubs(edges: list[Edge], top_n: int = 15) -> list[dict[str, Any]]:
    degree: Counter[str] = Counter()
    for e in edges:
        degree[e.id_a] += 1
        degree[e.id_b] += 1
    ranked = sorted(degree.items(), key=lambda kv: (-kv[1], kv[0]))
    return [{"id": node_id, "degree": deg} for node_id, deg in ranked[:top_n]]


def main() -> None:
    # NOTE (test_scripts_call_preflight.py): this script has no region_id/n_points
    # sweep-region concept to preflight -- it is a fixed-N closed-form pairwise
    # sweep over the already-catalogued node set, the same exemption category as
    # #317/#606/#608/#614/#624/#641/#642/#649 (see that test file's
    # _LEGACY_EXEMPT entry for this script).
    t_start = time.time()
    print(f"[{time.strftime('%Y-%m-%dT%H:%M:%S')}] loading {CATALOGUE_PATH} ...")
    rows = _load_catalogue_rows()
    print(f"  {len(rows)} total catalogue rows")

    nodes = [r for r in rows if is_node(r)]
    print(f"  {len(nodes)} eligible nodes (orbit_class cycler/quasi_cycler, >=1 usable encounter)")

    pairs = candidate_pairs(nodes)
    ts = time.strftime("%Y-%m-%dT%H:%M:%S")
    print(f"[{ts}] {len(pairs)} candidate pairs (sharing >=1 usable body)")

    edges: list[Edge] = []
    band_histogram_overall: Counter[str] = Counter()
    band_histogram_by_body: dict[str, Counter[str]] = defaultdict(Counter)
    phase_status_counts: Counter[str] = Counter()
    cheap_edges: list[Edge] = []
    cheap_dv_phase_indeterminate = 0

    t_sweep = time.time()
    for row_a, row_b in pairs:
        shared = usable_bodies(row_a) & usable_bodies(row_b)
        for body in sorted(shared):
            edge = compute_edge(row_a, row_b, body)
            edges.append(edge)
            band_histogram_overall[edge.band] += 1
            band_histogram_by_body[edge.body][edge.band] += 1
            phase_status_counts[edge.phase.status] += 1
            if edge.cheap_edge:
                cheap_edges.append(edge)
            is_cheap_band = edge.band in ("B0_ballistic_compatible", "B1_cheap")
            if is_cheap_band and edge.phase.status == "phase_locked":
                cheap_dv_phase_indeterminate += 1

    elapsed_sweep = time.time() - t_sweep
    ts = time.strftime("%Y-%m-%dT%H:%M:%S")
    print(f"[{ts}] {len(edges)} edges computed in {elapsed_sweep:.1f}s")

    # Data-gap tallies.
    no_period_data_node_count = sum(1 for r in nodes if period_days(r) is None)
    null_vinf_row_count_eligible = 0
    for r in nodes:
        entries = r.get("vinf_kms_at_encounters") or []
        if any(e.get("vinf_kms") is None for e in entries):
            null_vinf_row_count_eligible += 1
    null_vinf_entries_total_catalogue = sum(
        1
        for r in rows
        for e in (r.get("vinf_kms_at_encounters") or [])
        if e.get("vinf_kms") is None
    )
    null_vinf_row_count_catalogue = sum(
        1
        for r in rows
        if any(e.get("vinf_kms") is None for e in (r.get("vinf_kms_at_encounters") or []))
    )

    components = _connected_components(cheap_edges)
    hubs = _degree_ranked_hubs(cheap_edges)

    _OUT_DIR.mkdir(parents=True, exist_ok=True)
    edges_path = _OUT_DIR / "edges.jsonl"
    with edges_path.open("w") as fh:
        for e in edges:
            fh.write(json.dumps(e.to_json()) + "\n")

    summary = {
        "task": "650",
        "generated_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "elapsed_s": time.time() - t_start,
        "conventions": {
            "b0_delta_vinf_kms": B0_DELTA_VINF_KMS,
            "b1_dv_hop_kms": B1_DV_HOP_KMS,
            "b2_dv_hop_kms": B2_DV_HOP_KMS,
            "dv_phase_gate_kms": DV_PHASE_GATE_KMS,
            "statistical_horizon_years": STATISTICAL_HORIZON_YEARS,
            "n_delta0_samples": N_DELTA0_SAMPLES,
            "sensitivity_windows_days": list(SENSITIVITY_WINDOWS_DAYS),
            "phase_lock_p_align_threshold": PHASE_LOCK_P_ALIGN_THRESHOLD,
            "cheap_edge_wait_years": CHEAP_EDGE_WAIT_YEARS,
            "cheap_edge_duty_adjusted_p_align": CHEAP_EDGE_DUTY_ADJUSTED_P_ALIGN,
            "direction_data": "absent",
            "r_p_floor_convention": "radius_eq_km + safe_alt_km (#426/#427 floor convention)",
            "dv_hop_formula": (
                "|sqrt(v_a^2 + 2*mu/r_p) - sqrt(v_b^2 + 2*mu/r_p)| "
                "(Oberth-optimal same-body powered-flyby handoff, lower bound -- "
                "no V_inf direction data in the catalogue)"
            ),
        },
        "counts": {
            "n_catalogue_rows": len(rows),
            "n_nodes": len(nodes),
            "n_candidate_pairs": len(pairs),
            "n_edges": len(edges),
        },
        "band_histogram_overall": dict(band_histogram_overall),
        "band_histogram_by_body": {b: dict(c) for b, c in band_histogram_by_body.items()},
        "phase_status_counts": dict(phase_status_counts),
        "cheap_edge_count": len(cheap_edges),
        "cheap_edges": [
            {"id_a": e.id_a, "id_b": e.id_b, "body": e.body, "dv_hop_kms": e.dv_hop_kms}
            for e in cheap_edges
        ],
        "cheap_dv_phase_indeterminate_count": cheap_dv_phase_indeterminate,
        "data_gaps": {
            "no_period_data_node_count": no_period_data_node_count,
            "null_vinf_row_count_eligible_nodes": null_vinf_row_count_eligible,
            "null_vinf_row_count_full_catalogue": null_vinf_row_count_catalogue,
            "null_vinf_entries_total_full_catalogue": null_vinf_entries_total_catalogue,
        },
        "cheap_edge_subgraph": {
            "n_components": len(components),
            "component_sizes": [len(c) for c in components],
            "components": components,
            "top_degree_hubs": hubs,
        },
    }
    summary_path = _OUT_DIR / "summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=False) + "\n")

    print(f"[{time.strftime('%Y-%m-%dT%H:%M:%S')}] wrote {edges_path} ({len(edges)} records)")
    print(f"  wrote {summary_path}")
    print()
    print("HEADLINE CENSUS:")
    print(f"  nodes={len(nodes)}  candidate_pairs={len(pairs)}  edges={len(edges)}")
    print(f"  band histogram: {dict(band_histogram_overall)}")
    print(f"  phase_status counts: {dict(phase_status_counts)}")
    print(
        f"  cheap_edge_count={len(cheap_edges)} "
        f"cheap_dv_phase_indeterminate={cheap_dv_phase_indeterminate}"
    )
    print(
        f"  data gaps: no_period_data_node_count={no_period_data_node_count} "
        f"null_vinf_row_count_eligible_nodes={null_vinf_row_count_eligible}"
    )
    sizes = [len(c) for c in components][:10]
    print(f"  cheap-edge subgraph: {len(components)} components, sizes={sizes}...")


if __name__ == "__main__":
    main()
