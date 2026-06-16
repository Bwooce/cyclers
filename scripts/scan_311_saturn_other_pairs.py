"""#311 Phase 1 Part B -- other Saturnian moon-pair systems.

#285 covered Titan-Rhea-Dione-Tethys (length-3) and Titan-Enceladus-Rhea
(length-5). Part A confirmed Rhea-Dione-Rhea (1,1) at 0.107 km/s is the
genome ceiling for that pair.

Part B asks: do the OTHER regular Saturnian moon-pair combinations
produce nearer-to-gate closures, or a sub-gate (SILVER) closure?

The brief lists five system variants:
  * Titan-Rhea (Titan ~4x Rhea distance from Saturn)
  * Titan-Iapetus (Iapetus is the outermost regular, ~6x Titan distance)
  * Iapetus-Hyperion (the irregular nearest Iapetus)
  * Tethys-Dione-Rhea (3-body, Takubo 2210.14996 tour-design prior)
  * Enceladus-Tethys-Dione (3-body)

Pair-only length-3 cycles enumerate by ``RepeatedMoonTarget`` with
``seq_lengths=(3,)``; 3-body cases use ``seq_lengths=(5,)``. n_rev_grid
covers (k1, k2) in [0, 5]^2 for 2-body and (0..3)^4 for 3-body (matching
#312 Part B's sizing choices to keep length-5 runtime tractable).

NO catalogue writeback. NO novelty claims. Sourced moon mu and SMA from
``src/cyclerfinder/core/satellites.py`` (JPL DE440 + sat441).

Outputs:
  * ``data/scan_311_saturn_titan_rhea.jsonl``
  * ``data/scan_311_saturn_titan_iapetus.jsonl``
  * ``data/scan_311_saturn_iapetus_hyperion.jsonl``
  * ``data/scan_311_saturn_tethys_dione_rhea.jsonl``
  * ``data/scan_311_saturn_enceladus_tethys_dione.jsonl``

Plus a top-level summary index ``data/scan_311_saturn_other_pairs_index.jsonl``.

Run as::

    uv run python scripts/scan_311_saturn_other_pairs.py
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from cyclerfinder.search.saturn_uranus_campaign import run_prioritized_scan  # noqa: E402


def _git_sha() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"], text=True, cwd=ROOT
        ).strip()
    except Exception:
        return "unknown"


# 2-body length-3 pair scans.
# Each scan restricts ``moons`` to a 2-element set so the closed length-3
# enumeration emits exactly the two cycles on that pair (A-B-A and B-A-B).
TITAN_RHEA = dict(
    primary="Saturn",
    moons=("Titan", "Rhea"),
    seq_lengths=(3,),
    n_rev_grid=(0, 1, 2, 3, 4, 5),
)

TITAN_IAPETUS = dict(
    primary="Saturn",
    moons=("Titan", "Iapetus"),
    seq_lengths=(3,),
    n_rev_grid=(0, 1, 2, 3, 4, 5),
)

IAPETUS_HYPERION = dict(
    primary="Saturn",
    moons=("Iapetus", "Hyperion"),
    seq_lengths=(3,),
    n_rev_grid=(0, 1, 2, 3, 4, 5),
)

# 3-body length-5 scans. Takubo 2210.14996 (Saturn tour design) is the
# strongest existence prior for Tethys-Dione-Rhea sub-tours; matching #312
# Part B sizing keeps runtime tractable (4^4 = 256 (k1,k2,k3,k4) cells per
# sequence permutation).
TETHYS_DIONE_RHEA = dict(
    primary="Saturn",
    moons=("Tethys", "Dione", "Rhea"),
    seq_lengths=(5,),
    n_rev_grid=(0, 1, 2, 3),
)

ENCELADUS_TETHYS_DIONE = dict(
    primary="Saturn",
    moons=("Enceladus", "Tethys", "Dione"),
    seq_lengths=(5,),
    n_rev_grid=(0, 1, 2, 3),
)


SYSTEMS = [
    ("titan_rhea", "Titan-Rhea (k=3)", TITAN_RHEA),
    ("titan_iapetus", "Titan-Iapetus (k=3)", TITAN_IAPETUS),
    ("iapetus_hyperion", "Iapetus-Hyperion (k=3)", IAPETUS_HYPERION),
    (
        "tethys_dione_rhea",
        "Tethys-Dione-Rhea (k=5, 3-body; Takubo prior)",
        TETHYS_DIONE_RHEA,
    ),
    (
        "enceladus_tethys_dione",
        "Enceladus-Tethys-Dione (k=5, 3-body)",
        ENCELADUS_TETHYS_DIONE,
    ),
]


def main() -> int:
    sha = _git_sha()
    print(f"[311-B] Other Saturnian moon-pair scans -- sha={sha}", flush=True)
    index_rows: list[dict[str, Any]] = []
    for slug, label, params in SYSTEMS:
        out_path = ROOT / "data" / f"scan_311_saturn_{slug}.jsonl"
        print(f"[311-B] {label} -- out={out_path}", flush=True)
        summary = run_prioritized_scan(
            out_path=out_path,
            gate_residual_kms=0.05,
            phase_samples=12,
            git_sha=sha,
            **params,
        )
        sd = summary.as_dict()
        print(f"[311-B] {label} summary: {sd}", flush=True)

        # Append summary row to the per-system JSONL.
        with out_path.open("a", encoding="utf-8") as fh:
            fh.write(
                json.dumps(
                    {
                        "_meta": True,
                        "kind": "summary",
                        "summary": sd,
                        "git_sha": sha,
                    }
                )
                + "\n"
            )

        index_rows.append(
            {
                "slug": slug,
                "label": label,
                "summary": sd,
                "path": str(out_path.relative_to(ROOT)),
            }
        )

    # Top-level index across all the systems.
    out_index = ROOT / "data" / "scan_311_saturn_other_pairs_index.jsonl"
    with out_index.open("w", encoding="utf-8") as fh:
        fh.write(
            json.dumps(
                {
                    "_meta": True,
                    "kind": "index",
                    "task": "#311 Phase 1 Part B -- Saturn other moon-pair systems",
                    "systems": index_rows,
                    "git_sha": sha,
                }
            )
            + "\n"
        )
    print(f"[311-B] DONE -- index: {out_index}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
