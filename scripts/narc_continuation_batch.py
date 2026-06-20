"""#388 N-arc real-ephemeris continuation batch driver.

Walks every descriptor-bearing catalogue row through the N-arc real-ephemeris
homotopy-ramp continuation lane
(:mod:`cyclerfinder.search.narc_continuation`) and records, per row, whether a
true-ephemeris (DE440 via ``Ephemeris("astropy")``) closure exists whose emerged
V_inf magnitudes match the row's SOURCED Earth/Mars anchors and whose
intermediate flybys are bend-feasible.

This is the decisive #388 finding: does a genuine in-basin Russell parent reach
DE440 and CLOSE? It NEVER writes back to the catalogue. The only thing it does on
a closing ``mcconaghy-2006-em-k2`` is print a HELD ``PROPOSED V0->V1`` marker
for human review.

Spec: docs (see the #388 N-arc continuation lane). Runlog: data/runs/.
"""

from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Any

from cyclerfinder.core.constants import VINF_CEILING_KMS
from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.data.catalog import load_catalog
from cyclerfinder.search.cycler_assembly import assemble_cycler, descriptor_to_phsi
from cyclerfinder.search.generic_return import RussellModel
from cyclerfinder.search.narc_continuation import (
    candidate_epochs,
    narc_continuation_correct,
    russell_parent_to_ballistic_seed,
)

# The inner corrector takes a single scalar V_inf ceiling. The encounters in the
# Earth-Mars cyclers handled here are E and M, so the binding physical ceiling is
# Mars's (the lowest of the two): use it as the cap. (VINF_CEILING_KMS itself is a
# per-body dict, not a scalar.)
VINF_CAP_KMS = VINF_CEILING_KMS["M"]


def main() -> None:
    m = RussellModel()
    ephem = Ephemeris("astropy")

    t0_wall = time.time()
    rows = load_catalog().entries

    records: list[dict[str, Any]] = []
    n_descriptor = 0
    n_converged = 0
    n_anchor_matched = 0
    n_promotions = 0

    for e in rows:
        phsi = descriptor_to_phsi(e.raw)
        if phsi is None:
            continue
        cyc = assemble_cycler(m, phsi)
        if cyc is None:
            continue
        try:
            seed = russell_parent_to_ballistic_seed(m, cyc, e.raw)
        except ValueError:
            continue

        n_descriptor += 1

        # The parent's exact beginning E-M relative phase is not directly exposed
        # by the idealized Cycler, so v1 uses target_phase=0.0 and relies on the
        # LaunchWindow scan over candidate epochs to find a closing epoch (a
        # documented simplification; the epoch grid still spans Russell's
        # LaunchWindow 1..21).
        target_phase = 0.0
        epochs = candidate_epochs(
            ephem, target_phase, launch_window_synodics=range(1, 22), grid=100
        )

        res = narc_continuation_correct(
            seed, epochs=epochs, final_ephemeris=ephem, vinf_cap=VINF_CAP_KMS
        )

        anchor_e = seed.vinf_anchor_e_kms
        anchor_m = seed.vinf_anchor_m_kms
        best_e = min((abs(v - anchor_e) for v in res.emerged_vinf_kms), default=float("inf"))
        best_m = min((abs(v - anchor_m) for v in res.emerged_vinf_kms), default=float("inf"))
        anchor_match = best_e <= 0.5 and best_m <= 0.5

        vlevel = e.raw.get("validation_level", "V0")

        promote = (
            e.id == "mcconaghy-2006-em-k2" and res.converged and anchor_match and res.bend_feasible
        )

        if res.converged:
            n_converged += 1
        if anchor_match:
            n_anchor_matched += 1
        if promote:
            n_promotions += 1

        wall = time.time() - t0_wall
        vinf_disp = [round(v, 2) for v in res.emerged_vinf_kms]
        print(
            f"[{wall:.0f}s] {e.id:24s} [{vlevel}] conv={res.converged!s:5} "
            f"res={res.max_residual_kms:.3e} vinf={vinf_disp} "
            f"anchorE={anchor_e:.2f} anchorM={anchor_m:.2f} "
            f"match={anchor_match!s:5} bend={res.bend_feasible!s:5}"
            f"{' *** PROPOSED V0->V1 (HELD) ***' if promote else ''}",
            flush=True,
        )

        records.append(
            {
                "id": e.id,
                "validation_level": vlevel,
                "converged": res.converged,
                "max_residual_kms": res.max_residual_kms,
                "emerged_vinf_kms": list(res.emerged_vinf_kms),
                "anchor_e_kms": anchor_e,
                "anchor_m_kms": anchor_m,
                "best_e_err_kms": best_e,
                "best_m_err_kms": best_m,
                "anchor_match": anchor_match,
                "bend_feasible": res.bend_feasible,
                "t0_sec": res.t0_sec,
                "tof_days": list(res.tof_days),
                "winning_epoch_sec": res.winning_epoch_sec,
                "promote_held": promote,
            }
        )

    stamp = time.strftime("%Y%m%dT%H%M%SZ", time.gmtime())
    out_dir = Path("data/runs")
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"narc-continuation-{stamp}.jsonl"
    with out_path.open("w", encoding="utf-8") as fh:
        for rec in records:
            fh.write(json.dumps(rec) + "\n")

    print()
    print(f"runlog: {out_path}")
    print(f"descriptor rows:       {n_descriptor}")
    print(f"converged:             {n_converged}")
    print(f"anchor-matched:        {n_anchor_matched}")
    print(f"proposed-promotions:   {n_promotions} (held)")
    print("HELD — no writeback")


if __name__ == "__main__":
    main()
