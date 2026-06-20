"""Phase-D batch: assemble Russell cyclers from every mappable catalogue row.

Iterates the catalogue, maps each row to a Russell p.h.s.i descriptor
(:func:`descriptor_to_phsi`), and attempts to assemble a concrete generic-return
cycler (:func:`assemble_cycler`). For each assembled candidate it compares the
emerged Earth v_inf against the sourced anchor (campaign gate: within 0.5 km/s)
and reports the AR / TR feasibility gates.

This is a read-only diagnostic: it writes a per-row runlog under ``data/runs/``
and prints a summary. It NEVER writes back to ``data/catalogue.yaml``.
"""

from __future__ import annotations

import json
import time
from pathlib import Path

from cyclerfinder.data.catalog import load_catalog
from cyclerfinder.search.cycler_assembly import assemble_cycler, descriptor_to_phsi
from cyclerfinder.search.generic_return import RussellModel

VINF_MATCH_TOL_KMS = 0.5  # campaign gate on emerged-vs-sourced Earth v_inf


def main() -> None:
    model = RussellModel()
    catalog = load_catalog()

    stamp = time.strftime("%Y%m%dT%H%M%SZ", time.gmtime())
    runlog_path = Path("data/runs") / f"russell-cycler-{stamp}.jsonl"
    runlog_path.parent.mkdir(parents=True, exist_ok=True)

    n_mappable = 0
    n_assembled = 0
    n_vinf_matched = 0
    n_tr_gt1 = 0

    with runlog_path.open("w", encoding="utf-8") as fh:
        for entry in catalog.entries:
            spec = descriptor_to_phsi(entry.raw)
            if spec is None:
                continue
            n_mappable += 1

            vlevel = entry.raw.get("validation_level")
            cyc = assemble_cycler(model, spec)

            if cyc is None:
                print(f"{entry.id} [{vlevel}] NO-ASSEMBLE vinf_e_sourced={spec.vinf_e_kms:.3f}")
                record = {
                    "id": entry.id,
                    "validation_level": vlevel,
                    "assembled": False,
                    "vinf_e_sourced": spec.vinf_e_kms,
                    "vinf_m_sourced": spec.vinf_m_kms,
                    "p": spec.p,
                    "h": spec.h,
                    "s": spec.s,
                    "i": spec.i,
                }
                fh.write(json.dumps(record) + "\n")
                continue

            n_assembled += 1
            vinf_e_match = abs(cyc.vinf_e - spec.vinf_e_kms) <= VINF_MATCH_TOL_KMS
            if vinf_e_match:
                n_vinf_matched += 1
            if cyc.tr > 1.0:
                n_tr_gt1 += 1

            print(
                f"{entry.id} [{vlevel}] assembled={cyc is not None} "
                f"vinf_e_emerged={cyc.vinf_e:.3f} vinf_e_sourced={spec.vinf_e_kms:.3f} "
                f"AR={cyc.ar:.3f} TR={cyc.tr:.3f}"
            )
            record = {
                "id": entry.id,
                "validation_level": vlevel,
                "assembled": True,
                "vinf_e_emerged": cyc.vinf_e,
                "vinf_e_sourced": spec.vinf_e_kms,
                "vinf_e_match": vinf_e_match,
                "vinf_m_sourced": spec.vinf_m_kms,
                "ar": cyc.ar,
                "tr": cyc.tr,
                "n_revs": cyc.generic_return.n_revs,
                "branch": cyc.generic_return.branch,
                "p": cyc.p,
                "h": cyc.h,
                "s": cyc.s,
                "i": cyc.i,
            }
            fh.write(json.dumps(record) + "\n")

    print()
    print(f"runlog: {runlog_path}")
    print(
        f"summary: mappable={n_mappable} assembled={n_assembled} "
        f"vinf_matched(<=0.5km/s)={n_vinf_matched} TR>1={n_tr_gt1}"
    )
    print("HELD - no writeback")


if __name__ == "__main__":
    main()
