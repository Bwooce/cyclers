"""#404/#388 Component 4 — DSM multi-arc closure-lane batch driver.

Runs the existing Takao eta-DSM closer (search/dsm_descriptor_seed.close_row_dsm)
over every DESCRIPTOR-BEARING catalogue row (those for which
seed_dsm_chain_from_descriptor returns a seed) on the real DE440 ephemeris, and
records per-row outcome: corrector convergence, the EMERGED per-encounter V_inf vs
the row's SOURCED Russell-table anchor (golden separation — the anchor is never
self-computed), the interior-impulse DSM dV, and the V1 anchor-match verdict.

Honesty contract (orbit-closure-discipline): NO catalogue writeback here. Rows
already at V1+ are re-run as a regression/consistency check; the only V0
descriptor-bearing row (mcconaghy-2006-em-k2, the #404-triaged in-scope #365
negative) is the genuine promotion ATTEMPT. A pass is a PROPOSED V0->V1 held for
session review (the single-arc-degenerate guard is satisfied for it by record: the
#365 single-arc closer already returned NO-CLOSE). n-body V3 confirmation is a
separate downstream step.

Run: uv run python scripts/dsm_closure_batch.py
"""

from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Any

import yaml

from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.search.dsm_descriptor_seed import (
    V1_TOLERANCE_KMS,
    close_row_dsm,
    seed_dsm_chain_from_descriptor,
)

_CATALOGUE = Path("data/catalogue.yaml")


def _n_arc_tofs(row: dict[str, Any]) -> int:
    return sum(1 for a in (row.get("free_return_arcs") or []) if a.get("tof_years") is not None)


def main() -> None:
    rows: list[dict[str, Any]] = yaml.safe_load(_CATALOGUE.read_text())
    seedable: list[dict[str, Any]] = []
    no_close: list[str] = []  # descriptor present (>=2 arcs) but no usable seed (arc doesn't reach)
    for r in rows:
        if seed_dsm_chain_from_descriptor(r) is not None:
            seedable.append(r)
        elif _n_arc_tofs(r) >= 2:
            no_close.append(r["id"])  # the spec's OFF-FAMILY-NO-CLOSE descriptor rows
    print(f"descriptor-bearing (seedable) rows: {len(seedable)} / {len(rows)}", flush=True)
    if no_close:
        print(f"descriptor-present-but-NO-CLOSE (arc does not reach): {no_close}", flush=True)
    print(f"V1 anchor tolerance = {V1_TOLERANCE_KMS} km/s\n", flush=True)

    ephem = Ephemeris("astropy")  # real DE440
    t0 = time.time()
    records: list[dict[str, Any]] = []
    for r in seedable:
        rid = r["id"]
        vlevel = r.get("validation_level", "V0")
        t1 = time.time()
        try:
            res = close_row_dsm(r, ephem)
            err = None
        except Exception as exc:  # report, never abort the batch
            res = None
            err = f"{type(exc).__name__}: {exc}"
        dt = time.time() - t1
        if res is None:
            print(f"[{time.time() - t0:6.0f}s] {rid:24s} [{vlevel}] EXC {err}", flush=True)
            records.append({"id": rid, "validation_level": vlevel, "error": err})
            continue
        # The only V0 seedable row is the genuine promotion attempt.
        promote = (
            vlevel == "V0" and res.converged and res.anchor_match and not res.hyperbolic_impossible
        )
        rec = {
            "id": rid,
            "validation_level": vlevel,
            "converged": res.converged,
            "max_residual_kms": res.max_residual_kms,
            "vinf_anchor_kms": res.vinf_anchor_kms,
            "vinf_per_encounter_kms": list(res.vinf_per_encounter_kms),
            "anchor_match": res.anchor_match,
            "dsm_dv_kms": list(res.dv_dsm_kms),
            "hyperbolic_impossible": res.hyperbolic_impossible,
            "proposed_promotion": "V0->V1" if promote else None,
            "wall_s": round(dt, 1),
        }
        records.append(rec)
        tag = "  *** PROPOSED V0->V1 (held) ***" if promote else ""
        print(
            f"[{time.time() - t0:6.0f}s] {rid:24s} [{vlevel}] "
            f"conv={res.converged!s:5} res={res.max_residual_kms:.3e} "
            f"anchor={res.vinf_anchor_kms:.3f} match={res.anchor_match!s:5} "
            f"dsmdV={sum(res.dv_dsm_kms):.3f}{tag}",
            flush=True,
        )

    # --- runlog (jsonl) ---
    runs = Path("data/runs")
    runs.mkdir(parents=True, exist_ok=True)
    stamp = time.strftime("%Y%m%dT%H%M%S")
    runlog = runs / f"dsm-closure-{stamp}.jsonl"
    with runlog.open("w") as fh:
        for rec in records:
            fh.write(json.dumps(rec) + "\n")

    # --- summary ---
    conv = [r for r in records if r.get("converged")]
    matched = [r for r in conv if r.get("anchor_match")]
    proposed = [r for r in records if r.get("proposed_promotion")]
    print(f"\n=== summary ({time.time() - t0:.0f}s) ===")
    print(f"seedable: {len(seedable)}  converged: {len(conv)}  anchor-matched: {len(matched)}")
    print(f"proposed promotions (HELD for review): {[r['id'] for r in proposed]}")
    print(f"runlog: {runlog}")


if __name__ == "__main__":
    main()
