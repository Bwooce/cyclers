"""#312 Phase 1 Part B -- other Uranian moon-pair systems.

#285 covered Titania-Oberon-Umbriel as one band (and Ariel-Umbriel-Titania
as the 3-body length-5 band). Part A re-examined the best near-miss
(Oberon-Titania-Oberon (1,1)) at finer (k1, k2) and phase resolution and
confirmed the residual is the GENOME ceiling, not a phase artifact.

Part B asks: do the OTHER regular Uranian moon-pair combinations produce
nearer-to-gate closures, or are they uniformly worse?

The brief lists five system variants. Three (Titania-Umbriel,
Oberon-Umbriel, Ariel-Titania) are pair-only length-3 cycles --
``RepeatedMoonTarget`` enumerates them by including the relevant moons in
the registry-set. The two 3-body cases (Miranda-Ariel-Umbriel and
Ariel-Umbriel-Titania) need length-5 cycles. #285 Band B already covered
Ariel-Umbriel-Titania length-5; we re-include it here only for direct
comparison at the wider (k1..k4) grid.

NO catalogue writeback. NO novelty claims. Sourced Uranian-moon
mu/SMA values come from the JPL DE440 satellites registry already in
``src/cyclerfinder/core/satellites.py``.

Outputs:
  * ``data/scan_312_uranus_titania_umbriel.jsonl``
  * ``data/scan_312_uranus_oberon_umbriel.jsonl``
  * ``data/scan_312_uranus_ariel_titania.jsonl``
  * ``data/scan_312_uranus_miranda_ariel_umbriel.jsonl``
  * ``data/scan_312_uranus_ariel_umbriel_titania.jsonl``

Plus a top-level summary index ``data/scan_312_uranus_other_pairs_index.jsonl``.

Run as::

    uv run python scripts/scan_312_uranus_other_pairs.py
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


# 2-body length-3 pair scans. Each scan restricts ``moons`` to a 2-element
# set so the closed length-3 enumeration emits exactly the two cycles on
# that pair (A-B-A and B-A-B). The Lambert solver and the phasing logic
# are the same as #285. n_rev_grid covers (k1, k2) in [0, 5]^2 -- 36 cells
# per cycle.
TITANIA_UMBRIEL = dict(
    primary="Uranus",
    moons=("Titania", "Umbriel"),
    seq_lengths=(3,),
    n_rev_grid=(0, 1, 2, 3, 4, 5),
)

OBERON_UMBRIEL = dict(
    primary="Uranus",
    moons=("Oberon", "Umbriel"),
    seq_lengths=(3,),
    n_rev_grid=(0, 1, 2, 3, 4, 5),
)

ARIEL_TITANIA = dict(
    primary="Uranus",
    moons=("Ariel", "Titania"),
    seq_lengths=(3,),
    n_rev_grid=(0, 1, 2, 3, 4, 5),
)

# 3-body length-5 scans. RepeatedMoonTarget enumerates all closed length-5
# cycles on the 3-moon set (the open-vs-closed filter in ``_sequences``
# kills the non-closed permutations). #285 Band B used n_rev_grid (0,1,2);
# we use (0,1,2,3) here for a wider sweep but cap at 3 because length-5
# means 4 legs and (0..3)^4 = 256 (k1,k2,k3,k4) per sequence -- already
# heavier than the 2-body case.
MIRANDA_ARIEL_UMBRIEL = dict(
    primary="Uranus",
    moons=("Miranda", "Ariel", "Umbriel"),
    seq_lengths=(5,),
    n_rev_grid=(0, 1, 2, 3),
)

ARIEL_UMBRIEL_TITANIA = dict(
    primary="Uranus",
    moons=("Ariel", "Umbriel", "Titania"),
    seq_lengths=(5,),
    n_rev_grid=(0, 1, 2, 3),
)


SYSTEMS = [
    ("titania_umbriel", "Titania-Umbriel (k=3)", TITANIA_UMBRIEL),
    ("oberon_umbriel", "Oberon-Umbriel (k=3)", OBERON_UMBRIEL),
    ("ariel_titania", "Ariel-Titania (k=3)", ARIEL_TITANIA),
    (
        "miranda_ariel_umbriel",
        "Miranda-Ariel-Umbriel (k=5, 3-body)",
        MIRANDA_ARIEL_UMBRIEL,
    ),
    (
        "ariel_umbriel_titania",
        "Ariel-Umbriel-Titania (k=5, 3-body)",
        ARIEL_UMBRIEL_TITANIA,
    ),
]


def main() -> int:
    sha = _git_sha()
    print(f"[312-B] Other Uranian moon-pair scans -- sha={sha}", flush=True)
    index_rows: list[dict[str, Any]] = []
    for slug, label, params in SYSTEMS:
        out_path = ROOT / "data" / f"scan_312_uranus_{slug}.jsonl"
        print(f"[312-B] {label} -- out={out_path}", flush=True)
        summary = run_prioritized_scan(
            out_path=out_path,
            gate_residual_kms=0.05,
            phase_samples=12,
            git_sha=sha,
            **params,
        )
        sd = summary.as_dict()
        print(f"[312-B] {label} summary: {sd}", flush=True)

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
    out_index = ROOT / "data" / "scan_312_uranus_other_pairs_index.jsonl"
    with out_index.open("w", encoding="utf-8") as fh:
        fh.write(
            json.dumps(
                {
                    "_meta": True,
                    "kind": "index",
                    "task": "#312 Phase 1 Part B -- Uranus other moon-pair systems",
                    "systems": index_rows,
                    "git_sha": sha,
                }
            )
            + "\n"
        )
    print(f"[312-B] DONE -- index: {out_index}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
