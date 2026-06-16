"""#341 Phase 1 Part B -- other Neptune-system pair / 3-body sweeps.

The brief enumerates four candidate Part B configurations:

  * Triton-Nereid
  * Triton-Proteus (reverse of A)
  * Nereid-Triton
  * Proteus-Triton-Nereid (3-body)

THREE OF FOUR REQUIRE NEREID and Nereid is INTENTIONALLY OMITTED from
``src/cyclerfinder/core/satellites.py`` per the sourced-golden discipline:

    "Nereid OMITTED: JPL SSD lists its GM as 0.0 (mass not determined),
    so per the sourcing discipline we omit it rather than guess from a
    size estimate."
    -- src/cyclerfinder/core/satellites.py lines 175-176

Per ``feedback_golden_tests_sourced_only`` and the task discipline
("Sourced golden for Neptune system mu values: NASA Neptune fact sheet
+ JPL Horizons"), I do NOT fabricate a Nereid GM. The sweeps that
require Nereid cannot be executed within the sourced-golden discipline
and are recorded here as a registered empty region (per
``project_negative_results_registry``).

The fourth configuration (Triton-Proteus) is the REVERSE of Part A's
(Triton, Proteus) moon set. The Part A ``_sweep_one_cycle`` enumerator
already covers both closed length-3 cycles on this pair --
Proteus-Triton-Proteus and Triton-Proteus-Triton -- so it is fully
captured by ``data/scan_341_neptune_proteus_triton_finer.jsonl``. There
is no additional ballistic territory to sweep in Part B as briefed.

This script writes a single ``data/scan_341_neptune_other_pairs.jsonl``
with a meta row + a verdict row documenting the source-discipline
constraint, the registered empty region, and the pointer to Part A's
JSONL for the only executable case. NO data fabrication.

Run as::

    uv run python scripts/scan_341_neptune_other_pairs.py
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from cyclerfinder.core.satellites import PRIMARIES, SATELLITES  # noqa: E402


def _git_sha() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"], text=True, cwd=ROOT
        ).strip()
    except Exception:
        return "unknown"


def _neptune_moons_registered() -> list[str]:
    return sorted(m for m, sat in SATELLITES.items() if sat.primary == "Neptune")


def main() -> int:
    sha = _git_sha()
    out_path = ROOT / "data" / "scan_341_neptune_other_pairs.jsonl"
    out_path.parent.mkdir(parents=True, exist_ok=True)

    registered = _neptune_moons_registered()
    mu_neptune_system = PRIMARIES["Neptune"]

    meta = {
        "_meta": True,
        "task": ("#341 Phase 1 Part B -- other Neptune-system pair / 3-body sweeps"),
        "primary": "Neptune",
        "moons_registered_in_satellites_py": registered,
        "neptune_system_mu_km3_s2": mu_neptune_system,
        "system_mu_source": (
            "JPL DE440 planetary constants (gm_de440), accessed 2026-06-14 "
            "via src/cyclerfinder/core/satellites.py"
        ),
        "git_sha": sha,
    }
    verdict = {
        "_meta": True,
        "kind": "verdict",
        "result": "NON_EXECUTABLE_BY_SOURCED_DISCIPLINE",
        "interpretation": (
            "Three of four Part B configurations (Triton-Nereid, "
            "Nereid-Triton, Proteus-Triton-Nereid) require Nereid, which "
            "is INTENTIONALLY OMITTED from src/cyclerfinder/core/satellites.py "
            "because JPL SSD lists its GM as 0.0 (mass not determined). Per "
            "feedback_golden_tests_sourced_only and the task discipline "
            "(sourced-golden for Neptune system mu values), the Nereid GM "
            "is not fabricated. These sweeps are NOT RUN; they are recorded "
            "here as a registered empty region (per "
            "project_negative_results_registry)."
        ),
        "non_executable_configurations": [
            {
                "label": "Triton-Nereid",
                "missing_body": "Nereid",
                "reason": "JPL SSD GM=0.0 (mass not determined); excluded by satellites.py:175-176",
            },
            {
                "label": "Nereid-Triton",
                "missing_body": "Nereid",
                "reason": "JPL SSD GM=0.0 (mass not determined); excluded by satellites.py:175-176",
            },
            {
                "label": "Proteus-Triton-Nereid",
                "missing_body": "Nereid",
                "reason": "JPL SSD GM=0.0 (mass not determined); excluded by satellites.py:175-176",
            },
        ],
        "covered_by_part_a_configurations": [
            {
                "label": "Triton-Proteus",
                "note": (
                    "Reverse of Part A's (Triton, Proteus) moon set; the "
                    "_sweep_one_cycle enumerator at Part A covers both "
                    "closed length-3 cycles (Proteus-Triton-Proteus AND "
                    "Triton-Proteus-Triton). See "
                    "data/scan_341_neptune_proteus_triton_finer.jsonl, "
                    "rows with sequence=[Triton, Proteus, Triton] (best "
                    "0.2251 km/s at (1,1), phase0=285 deg, rel_off=105 deg, "
                    "tof_scale=2.0)."
                ),
            },
        ],
        "registered_empty_region": {
            "primary": "Neptune",
            "moons_needed_but_unavailable": ["Nereid"],
            "reactivation_condition": (
                "Re-sweep when Nereid GM becomes a sourced quantity "
                "(JPL SSD GM updated above 0.0, or an independent peer-"
                "reviewed determination). Note: Nereid's e=0.75 makes the "
                "circular-coplanar moon-orbit assumption used by the "
                "current genome unsuitable even with a sourced GM; "
                "extending the genome to eccentric moons is its own task."
            ),
            "method_capability_tags": [
                "multi-arc",
                "patched-conic",
                "ballistic",
                "coplanar",
                "circular-moon",
            ],
        },
        "git_sha": sha,
    }

    with out_path.open("w", encoding="utf-8") as fh:
        fh.write(json.dumps(meta) + "\n")
        fh.write(json.dumps(verdict) + "\n")

    print(f"[341-B] Neptune registered moons: {registered}", flush=True)
    print(
        "[341-B] Verdict: NON_EXECUTABLE_BY_SOURCED_DISCIPLINE -- "
        "three of four Part B configurations need Nereid (omitted "
        "from satellites.py due to JPL SSD GM=0.0). The fourth "
        "(Triton-Proteus) is covered by Part A.",
        flush=True,
    )
    print(f"[341-B] wrote {out_path}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
