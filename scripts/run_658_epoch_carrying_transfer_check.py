#!/usr/bin/env python3
"""#658: epoch-locking pilot -- deterministic transfer check on the small
epoch-carrying subset of `#650`'s catalogue.

`#654`'s strategy review (shortlist item 4) hypothesized that 2-3 catalogue
rows -- named as "Jones 2017 VEM-triple" and "Aldrin/Byrnes lineage" --
already carry real, sourced, published epochs, and asked whether pointing
`#650`'s already-built, already-positive-controlled transfer-network
machinery (``cyclerfinder.data.transfer_network``) at JUST those rows could
produce a genuine, non-phase-indeterminate result, since #650's own finding
was that every `orbit_class: cycler` row is `epoch_locked: false` (a
deterministic phase verdict is otherwise impossible).

**Verified against the catalogue directly (not the #654 prose) -- the
premise only partially holds:**

* **Jones 2017 VEM-triple rows: 0 qualify.** Every row whose
  ``first_published.authors`` includes "Jones" (``jones-2017-vem-triple-family``,
  ``jones-2017-vem-emevve-outbound``, ``jones-2017-vem-meevem-inbound``,
  ``vem-emeeve-3syn``, ``hernandez-2017-jovian-ieg-triple-family``) has
  ``epoch_locked: false`` and ``launch_epoch: null`` -- same as every other
  `cycler`-class row. #654's framing of these as epoch-carrying was wrong.
* **"Aldrin/Byrnes lineage" rows: 9 qualify, but they are NOT the steady-state
  `aldrin-classic-em-k1-{outbound,inbound}` rows** (those are also
  `epoch_locked: false`, `cycler`-class, same as Jones). The rows that
  actually carry real ``epoch_locked: true`` + a real ``validity_window``
  are the Rogers/Hughes/Longuski/Aldrin 2012/2015 "establishment" precursor
  trajectories -- ``orbit_class: precursor_mga``, ``n_returns: 1`` one-shot
  V-infinity-leveraging insertion legs, each pointing (``inserts_into``) at a
  different steady-state E-M cycler (Aldrin x2, VISIT-1, VISIT-2, Case-1,
  Case-2, Case-3, S1L1, U0L1). `#650`'s own design doc (§2) deliberately
  EXCLUDES `precursor_mga` from the network's node set, precisely because a
  one-shot insertion leg is not a repeating schedule -- so these 9 rows were
  never in #650's 291-node graph at all. This script does NOT modify that
  exclusion or `transfer_network.is_node`/`ELIGIBLE_ORBIT_CLASSES`; it is an
  additive, narrowly-scoped analysis that calls
  ``transfer_network.compute_edge`` directly on a hand-identified small
  row set, exactly as the design's own epoch-window-intersection branch
  (§5 "Epoch-locked special case") already supports for ANY two rows
  with `epoch_locked=true` + `validity_window`, regardless of `orbit_class`.

Programmatic candidate-selection rule (visible, not hand-picked ids):
``epoch_locked is True AND validity_window.start/end are both non-null AND
"Aldrin" appears in first_published.authors`` -- OR the same rule with
"Jones" (which yields the empty set above). This survives a catalogue edit
better than a hardcoded id list.

For every pair of the 9 Aldrin-lineage rows sharing >=1 usable encounter
body (all 9 share "E"/Earth), calls
``cyclerfinder.data.transfer_network.compute_edge`` UNMODIFIED -- the exact
same dv_hop/band/epoch-window-intersection machinery #650 built and
positive-controlled, no new phase-check math written here.

No catalogue writes. Artifact:
``data/found/658_epoch_carrying_transfer_check/{edges.jsonl,summary.json}``
(the #317/#650 `data/found/<task>_.../` precedent).
"""

from __future__ import annotations

import itertools
import json
import sys
import time
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

import yaml  # type: ignore[import-untyped]

from cyclerfinder.data.catalog import CATALOGUE_PATH
from cyclerfinder.data.transfer_network import (
    Edge,
    compute_edge,
    epoch_window_intersection,
    usable_bodies,
)

_OUT_DIR = Path("data/found/658_epoch_carrying_transfer_check")


def _load_catalogue_rows() -> list[dict[str, Any]]:
    with open(CATALOGUE_PATH) as fh:
        data = yaml.safe_load(fh)
    assert isinstance(data, list)
    return data


def _cites_author(row: dict[str, Any], surname: str) -> bool:
    fp = row.get("first_published") or {}
    authors = fp.get("authors") or []
    return any(surname in a for a in authors)


def _has_real_epoch(row: dict[str, Any]) -> bool:
    if row.get("epoch_locked") is not True:
        return False
    vw = row.get("validity_window") or {}
    return vw.get("start") is not None and vw.get("end") is not None


def epoch_carrying_candidates(rows: list[dict[str, Any]], surname: str) -> list[dict[str, Any]]:
    """Rows citing ``surname`` in ``first_published.authors`` with a real
    ``epoch_locked=true`` + non-null ``validity_window`` (design-doc-§5-style
    "epoch-locked" rows) -- the programmatic form of #654's "already-sourced
    published epochs" framing, checked against the data rather than assumed.
    """
    return [r for r in rows if _cites_author(r, surname) and _has_real_epoch(r)]


def main() -> None:
    # NOTE (test_scripts_call_preflight.py): no region_id/n_points sweep-region
    # concept here -- a fixed, small (<=9-row) hand-identified analysis over
    # already-catalogued data, same exemption category as #650/#317/#606/#608
    # (see that test file's _LEGACY_EXEMPT entry for this script).
    t_start = time.time()
    ts = time.strftime("%Y-%m-%dT%H:%M:%S")
    print(f"[{ts}] loading {CATALOGUE_PATH} ...")
    rows = _load_catalogue_rows()
    print(f"  {len(rows)} total catalogue rows")

    jones_candidates = epoch_carrying_candidates(rows, "Jones")
    aldrin_candidates = epoch_carrying_candidates(rows, "Aldrin")
    print(f"  Jones-lineage epoch-carrying candidates: {len(jones_candidates)}")
    print(f"  Aldrin-lineage epoch-carrying candidates: {len(aldrin_candidates)}")
    for r in aldrin_candidates:
        vw = r["validity_window"]
        inserts = r.get("inserts_into")
        print(f"    {r['id']:40s} inserts_into={inserts!s:35s} {vw['start']}..{vw['end']}")

    candidates = aldrin_candidates  # Jones set is empty; nothing to pair there.
    pairs = list(itertools.combinations(sorted(candidates, key=lambda r: r["id"]), 2))
    ts = time.strftime("%Y-%m-%dT%H:%M:%S")
    print(f"[{ts}] {len(pairs)} candidate pairs among {len(candidates)} epoch-carrying rows")

    edges: list[Edge] = []
    edge_raw_overlap: dict[int, bool] = {}
    overlap_pairs: list[Edge] = []
    raw_overlap_pairs: list[Edge] = []
    for row_a, row_b in pairs:
        shared = usable_bodies(row_a) & usable_bodies(row_b)
        for body in sorted(shared):
            edge = compute_edge(row_a, row_b, body)
            edges.append(edge)
            # Raw calendar-window overlap, computed UNGATED directly via the
            # design's own `epoch_window_intersection` -- #650's own dv_hop
            # compute-gate (design §5: only run the phase model for dv_hop
            # <=1.0 km/s at a heliocentric body) is a compute optimization for
            # the 32k-pair general sweep, not a physical claim; for this
            # 36-pair micro-analysis it costs nothing to check every pair's
            # raw window overlap regardless of dv_hop, so both are reported.
            raw_overlap = epoch_window_intersection(row_a, row_b) is not None
            edge_raw_overlap[id(edge)] = raw_overlap
            if raw_overlap:
                raw_overlap_pairs.append(edge)
            if edge.phase.status == "epoch_window_overlap":
                overlap_pairs.append(edge)

    print(
        f"  {len(edges)} edges computed; {len(raw_overlap_pairs)} pairs have raw calendar-window "
        f"overlap (ungated); {len(overlap_pairs)} go through the design's dv_hop-gated "
        f"'epoch_window_overlap' branch (rest are already dv_hop>1.0, i.e. B2_moderate+ -- "
        f"never cheap_edge regardless of phase)"
    )
    for e in raw_overlap_pairs:
        gated = "not_computed_dv_gated" if e.phase.status != "epoch_window_overlap" else "computed"
        print(
            f"    {e.id_a} <-> {e.id_b} @ {e.body}: dv_hop={e.dv_hop_kms:.3f} km/s "
            f"band={e.band} cheap_edge={e.cheap_edge} phase_compute={gated}"
        )

    b0_edges = [e for e in edges if e.band == "B0_ballistic_compatible"]
    b0_edges_without_overlap = [e for e in b0_edges if e not in raw_overlap_pairs]
    print(f"  {len(b0_edges)} B0_ballistic_compatible edges total (any calendar window)")
    for e in b0_edges_without_overlap:
        print(
            f"    {e.id_a} <-> {e.id_b} @ {e.body}: dv_hop={e.dv_hop_kms:.3f} km/s "
            f"band={e.band} -- ballistically cheap but NO real calendar-window overlap "
            f"(the same 'ΔV-cheap, not real-date-realizable' pattern #650 found generally)"
        )

    _OUT_DIR.mkdir(parents=True, exist_ok=True)
    edges_path = _OUT_DIR / "edges.jsonl"
    with edges_path.open("w") as fh:
        for e in edges:
            fh.write(json.dumps(e.to_json()) + "\n")

    summary = {
        "task": "658",
        "generated_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "elapsed_s": time.time() - t_start,
        "finding": (
            "#654's named candidates (Jones 2017 VEM-triple, Aldrin/Byrnes) are "
            "NOT both epoch-carrying as framed. Jones-lineage rows: 0 qualify "
            "(all epoch_locked=false, launch_epoch=null, same as every other "
            "cycler-class row). Aldrin-lineage rows: 9 qualify, but they are the "
            "Rogers/Hughes/Longuski/Aldrin 2012/2015 one-shot precursor_mga "
            "'establishment' trajectories (n_returns=1), NOT the steady-state "
            "aldrin-classic-em-k1-* cycler rows (also epoch_locked=false). #650's "
            "design deliberately excludes precursor_mga from its node set for "
            "exactly this reason (one-shot, not a repeating schedule)."
        ),
        "candidate_selection_rule": (
            "epoch_locked is True AND validity_window.start/end both non-null AND "
            '"<surname>" in first_published.authors'
        ),
        "counts": {
            "n_catalogue_rows": len(rows),
            "n_jones_lineage_candidates": len(jones_candidates),
            "n_aldrin_lineage_candidates": len(aldrin_candidates),
            "n_candidate_pairs": len(pairs),
            "n_edges": len(edges),
            "n_raw_epoch_window_overlap_pairs": len(raw_overlap_pairs),
            "n_epoch_window_overlap_pairs_dv_gated_computed": len(overlap_pairs),
            "n_b0_ballistic_compatible_edges_any_window": len(b0_edges),
            "n_cheap_edges": sum(1 for e in edges if e.cheap_edge),
        },
        "aldrin_lineage_candidate_ids": [r["id"] for r in aldrin_candidates],
        "raw_epoch_window_overlap_edges": [
            {
                "id_a": e.id_a,
                "id_b": e.id_b,
                "body": e.body,
                "dv_hop_kms": e.dv_hop_kms,
                "delta_vinf_kms": e.delta_vinf_kms,
                "band": e.band,
                "cheap_edge": e.cheap_edge,
                "phase_status": e.phase.status,
            }
            for e in raw_overlap_pairs
        ],
        "b0_edges_without_calendar_overlap": [
            {
                "id_a": e.id_a,
                "id_b": e.id_b,
                "body": e.body,
                "dv_hop_kms": e.dv_hop_kms,
                "delta_vinf_kms": e.delta_vinf_kms,
                "note": (
                    "ballistically cheap (B0) but no real calendar-window overlap -- "
                    "the same 'DeltaV-cheap, not real-date-realizable' pattern #650 "
                    "found for the general heliocentric graph"
                ),
            }
            for e in b0_edges_without_overlap
        ],
        "case3_s1l1_b0_pair_caveat": (
            "case-3-4-3-2-establishment <-> s1l1-4-3-2-establishment is the ONLY "
            "B0_ballistic_compatible pair found (dv_hop=0.010 km/s). This is NOT a "
            "new finding: the catalogue's own notes on s1l1-4-3-2-establishment "
            "already state 'Same LD (12/20/2022) as the Case 3 establishment -- "
            "reflects the similar (a, e, peri, apo) parameter space of these two "
            "cyclers in Rogers Table 1' -- i.e. Rogers 2015's own Table 4 already "
            "documents these as a near-identical near-twin pair (same source launch "
            "date, near-equal Earth v_inf, near-equal orbital elements), not two "
            "independently-scheduled missions. Reported for completeness, not as a "
            "citable transfer-opportunity discovery."
        ),
    }
    summary_path = _OUT_DIR / "summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=False) + "\n")

    ts = time.strftime("%Y-%m-%dT%H:%M:%S")
    print(f"[{ts}] wrote {edges_path} ({len(edges)} records)")
    print(f"  wrote {summary_path}")


if __name__ == "__main__":
    main()
