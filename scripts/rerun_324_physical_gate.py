"""Re-run the physical-sanity flyby gate (#324) over every recent search JSONL.

Per `feedback_bugfix_invalidates_past_searches`: every correctness-affecting
fix mandates a sweep over past results computed under the buggy/missing code.
The #324 physical-sanity gate is bug-fix-equivalent — it would have rejected
candidates that earlier reports admitted (the prompt's worst-case Umbriel
flyby being the type case).

Scope (read-only over the source JSONLs; emits ONE new JSONL):

  * data/scan_285_saturn.jsonl, scan_285_uranus.jsonl
  * data/scan_312_uranus_*.jsonl (all per-pair JSONLs + offset_sweep)
  * data/scan_313_mars_phobos.jsonl, scan_313_mars_deimos.jsonl,
    scan_313_sun_jupiter_europa.jsonl, scan_313_sun_jupiter_io.jsonl
  * data/scan_309_low_thrust_em.jsonl, scan_309_low_thrust_vem.jsonl
  * data/scan_298_galileo_veega.jsonl
  * data/precursor_302_aldrin.jsonl, precursor_302_s1l1.jsonl

For each row that carries (sequence, per-encounter V_inf), the script runs
:func:`cyclerfinder.search.physical_sanity.candidate_passes_physical_gate`
and emits a one-line registry row. Rows without V_inf encounters (e.g.
tulip-orbit IC rows in scan_313 mars-phobos/deimos) are skipped with a
``status: skipped, reason: no_vinf_encounters`` entry so the registry shows
they were considered.

Output: data/rerun_324_physical_gate.jsonl

Discipline: NO catalogue writeback, NO novelty claims, NO modification of
source JSONLs.
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import subprocess
from dataclasses import asdict
from pathlib import Path
from typing import Any

from cyclerfinder.search.physical_sanity import (
    DEFAULT_MIN_USEFUL_BEND_DEG,
    candidate_passes_physical_gate,
)

DATA_DIR = Path(__file__).resolve().parent.parent / "data"

INPUT_FILES = [
    "scan_285_saturn.jsonl",
    "scan_285_uranus.jsonl",
    "scan_312_uranus_3d_probe.jsonl",
    "scan_312_uranus_ariel_titania.jsonl",
    "scan_312_uranus_ariel_umbriel_titania.jsonl",
    "scan_312_uranus_miranda_ariel_umbriel.jsonl",
    "scan_312_uranus_oberon_titania_finer.jsonl",
    "scan_312_uranus_oberon_umbriel.jsonl",
    "scan_312_uranus_other_pairs_index.jsonl",
    "scan_312_uranus_robustness.jsonl",
    "scan_312_uranus_titania_umbriel.jsonl",
    "scan_312_uranus_umbriel_oberon_offset_sweep.jsonl",
    "scan_313_mars_phobos.jsonl",
    "scan_313_mars_deimos.jsonl",
    "scan_313_sun_jupiter_europa.jsonl",
    "scan_313_sun_jupiter_io.jsonl",
    "scan_309_low_thrust_em.jsonl",
    "scan_309_low_thrust_vem.jsonl",
    "scan_298_galileo_veega.jsonl",
    "precursor_302_aldrin.jsonl",
    "precursor_302_s1l1.jsonl",
]

OUTPUT_FILE = DATA_DIR / "rerun_324_physical_gate.jsonl"


# Schema notes (inferred from inspection):
#   scan_285/312 (repeated-moon scans):
#       sequence: list[str], vinf_per_encounter_kms: list[float], verdict: str
#   scan_298 (Galileo VEEGA):
#       sequence: list[str], per_encounter_vinf_kms: list[float]
#       (also vinf_tuple_kms but that's the input shell, not the closed values)
#   scan_309 (low-thrust):
#       sequence: list[str], vinf_per_encounter_kms: list[float]
#   scan_313 mars_phobos/deimos: tulip-orbit ICs; NO vinf encounters → skipped
#   scan_313 sun_jupiter_europa/io: same shape as scan_313 phobos (CR3BP ICs) → skipped
#   precursor_302_aldrin/s1l1:
#       candidate.sequence: list[str], closure.per_encounter_vinf_kms: list[float]


def _extract_seq_vinfs(row: dict[str, Any]) -> tuple[tuple[str, ...], tuple[float, ...]] | None:
    """Return ``(sequence, vinfs)`` if the row carries flyby data; else ``None``.

    Handles the JSONL schema variants in the input files. Skips rows that
    don't carry per-encounter V_inf (e.g. CR3BP tulip ICs).
    """
    # Repeated-moon / low-thrust / 312 style.
    if "sequence" in row and "vinf_per_encounter_kms" in row:
        seq = row["sequence"]
        vinfs = row["vinf_per_encounter_kms"]
        if isinstance(seq, list) and isinstance(vinfs, list) and len(seq) == len(vinfs):
            return tuple(str(b) for b in seq), tuple(float(v) for v in vinfs)
    # 298-style (Galileo VEEGA).
    if "sequence" in row and "per_encounter_vinf_kms" in row:
        seq = row["sequence"]
        vinfs = row["per_encounter_vinf_kms"]
        if isinstance(seq, list) and isinstance(vinfs, list) and len(seq) == len(vinfs):
            return tuple(str(b) for b in seq), tuple(float(v) for v in vinfs)
    # 302 precursor style (candidate.sequence + closure.per_encounter_vinf_kms).
    if "candidate" in row and "closure" in row:
        cand = row["candidate"]
        closure = row["closure"]
        if isinstance(cand, dict) and isinstance(closure, dict):
            seq = cand.get("sequence")
            vinfs = closure.get("per_encounter_vinf_kms")
            if isinstance(seq, list) and isinstance(vinfs, list) and len(seq) == len(vinfs):
                return tuple(str(b) for b in seq), tuple(float(v) for v in vinfs)
    return None


def _row_signature(row: dict[str, Any]) -> str:
    """Best-effort short identifier for a candidate row (for the registry)."""
    for key in ("candidate_id", "cycler_id", "scan_id", "dedup_signature"):
        if row.get(key):
            return str(row[key])
    return ""


def _row_verdict(row: dict[str, Any]) -> str | None:
    """Original verdict if present (SILVER/GOLD/..)."""
    if "verdict" in row:
        return str(row["verdict"])
    return None


def _git_sha() -> str:
    try:
        result = subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"],
            check=True,
            capture_output=True,
            text=True,
            cwd=DATA_DIR.parent,
        )
        return result.stdout.strip()
    except Exception:
        return "unknown"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--min-useful-bend-deg",
        type=float,
        default=DEFAULT_MIN_USEFUL_BEND_DEG,
        help=f"useful-bend floor (deg), default {DEFAULT_MIN_USEFUL_BEND_DEG}",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=OUTPUT_FILE,
        help="output JSONL path (registry of per-row verdicts)",
    )
    args = parser.parse_args(argv)

    sha = _git_sha()
    timestamp = dt.datetime.now(dt.UTC).isoformat()

    meta = {
        "_meta": True,
        "task": "#324 physical-sanity flyby gate re-run",
        "min_useful_bend_deg": args.min_useful_bend_deg,
        "git_sha": sha,
        "timestamp_utc": timestamp,
        "inputs": INPUT_FILES,
        "note": (
            "Re-runs cyclerfinder.search.physical_sanity.candidate_passes_physical_gate "
            "over every prior session JSONL with per-encounter V_inf. Read-only on "
            "the source JSONLs. NO catalogue writeback."
        ),
    }

    totals: dict[str, int] = {
        "rows_scanned": 0,
        "rows_with_vinfs": 0,
        "rows_skipped_no_vinfs": 0,
        "rows_skipped_meta": 0,
        "passed_gate": 0,
        "failed_gate": 0,
        "silver_that_failed_gate": 0,
        "silver_that_passed_gate": 0,
    }
    silver_failures: list[dict[str, Any]] = []

    with args.out.open("w") as out_fh:
        out_fh.write(json.dumps(meta) + "\n")

        for fname in INPUT_FILES:
            path = DATA_DIR / fname
            if not path.exists():
                out_fh.write(
                    json.dumps(
                        {
                            "_file_missing": True,
                            "input_file": fname,
                            "note": "Input JSONL not present; skipped.",
                        }
                    )
                    + "\n"
                )
                continue
            file_rows_scanned = 0
            file_rows_with_vinfs = 0
            file_passed = 0
            file_failed = 0
            file_skipped = 0

            with path.open() as in_fh:
                for line in in_fh:
                    line = line.strip()
                    if not line:
                        continue
                    try:
                        row = json.loads(line)
                    except json.JSONDecodeError:
                        continue
                    totals["rows_scanned"] += 1
                    file_rows_scanned += 1

                    if row.get("_meta"):
                        totals["rows_skipped_meta"] += 1
                        continue

                    extracted = _extract_seq_vinfs(row)
                    if extracted is None:
                        totals["rows_skipped_no_vinfs"] += 1
                        file_skipped += 1
                        continue

                    seq, vinfs = extracted
                    totals["rows_with_vinfs"] += 1
                    file_rows_with_vinfs += 1
                    try:
                        passed, verdicts = candidate_passes_physical_gate(
                            seq, vinfs, min_useful_bend_deg=args.min_useful_bend_deg
                        )
                    except KeyError as exc:
                        # Unknown body — record but do NOT silently treat as pass.
                        err_record: dict[str, Any] = {
                            "input_file": fname,
                            "row_id": _row_signature(row),
                            "sequence": list(seq),
                            "vinfs_kms": list(vinfs),
                            "original_verdict": _row_verdict(row),
                            "physical_gate_passed": None,
                            "error": f"unknown body in sequence: {exc!s}",
                        }
                        out_fh.write(json.dumps(err_record) + "\n")
                        continue

                    original = _row_verdict(row)
                    record: dict[str, Any] = {
                        "input_file": fname,
                        "row_id": _row_signature(row),
                        "sequence": list(seq),
                        "vinfs_kms": list(vinfs),
                        "original_verdict": original,
                        "physical_gate_passed": bool(passed),
                        "max_bend_deg_per_encounter": [v.max_bend_deg for v in verdicts],
                        "is_useful_per_encounter": [v.is_useful for v in verdicts],
                        "verdicts": [asdict(v) for v in verdicts],
                    }
                    out_fh.write(json.dumps(record) + "\n")

                    if passed:
                        totals["passed_gate"] += 1
                        file_passed += 1
                        if original == "SILVER":
                            totals["silver_that_passed_gate"] += 1
                    else:
                        totals["failed_gate"] += 1
                        file_failed += 1
                        if original == "SILVER":
                            totals["silver_that_failed_gate"] += 1
                            silver_failures.append(
                                {
                                    "input_file": fname,
                                    "row_id": _row_signature(row),
                                    "sequence": list(seq),
                                    "vinfs_kms": list(vinfs),
                                    "max_bend_deg_per_encounter": [
                                        v.max_bend_deg for v in verdicts
                                    ],
                                }
                            )

            out_fh.write(
                json.dumps(
                    {
                        "_file_summary": True,
                        "input_file": fname,
                        "rows_scanned": file_rows_scanned,
                        "rows_with_vinfs": file_rows_with_vinfs,
                        "rows_skipped_no_vinfs": file_skipped,
                        "passed_gate": file_passed,
                        "failed_gate": file_failed,
                    }
                )
                + "\n"
            )

        out_fh.write(
            json.dumps(
                {
                    "_summary": True,
                    "totals": totals,
                    "silver_failures": silver_failures,
                }
            )
            + "\n"
        )

    # Final stdout summary.
    print(f"#324 physical-sanity gate re-run complete -> {args.out}")
    print(f"  rows_scanned          = {totals['rows_scanned']}")
    print(f"  rows_with_vinfs       = {totals['rows_with_vinfs']}")
    print(f"  rows_skipped_no_vinfs = {totals['rows_skipped_no_vinfs']}")
    print(f"  passed_gate           = {totals['passed_gate']}")
    print(f"  failed_gate           = {totals['failed_gate']}")
    print(f"  SILVER passed gate    = {totals['silver_that_passed_gate']}")
    print(f"  SILVER failed gate    = {totals['silver_that_failed_gate']}")
    if silver_failures:
        print("  SILVER survivors flagged unphysical:")
        for sf in silver_failures:
            print(
                f"    {sf['input_file']}  {sf['row_id']}  "
                f"seq={sf['sequence']}  vinfs={sf['vinfs_kms']}  "
                f"max_bend_deg={sf['max_bend_deg_per_encounter']}"
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
