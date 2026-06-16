"""Sweep the topology-audit harness across every recent discovery JSONL (#325).

Per the audit-orchestration discipline (#322 lesson):
  * #322 found that the tulip topology gate misclassified PLANAR Np-petal
    orbits as 3D tulips at extreme small-mu regimes.
  * Generalisation: any topology gate that counts a feature without checking
    the orthogonal degree of freedom is a #322-sibling-bug candidate.
  * Defensive sweep: re-verify every past row through an independent checker.
    Discrepancies are flagged, NOT retracted (read-only audit).
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import subprocess
from dataclasses import asdict
from pathlib import Path

from cyclerfinder.search.topology_audit import audit_topology

DEFAULT_SOURCES: tuple[str, ...] = (
    "data/family_296_3d_em_11.jsonl",
    "data/family_296_3d_subfamilies_299.jsonl",
    "data/scan_312_uranus_3d_probe.jsonl",
    "data/scan_312_uranus_ariel_titania.jsonl",
    "data/scan_312_uranus_ariel_umbriel_titania.jsonl",
    "data/scan_312_uranus_miranda_ariel_umbriel.jsonl",
    "data/scan_312_uranus_oberon_titania_finer.jsonl",
    "data/scan_312_uranus_oberon_umbriel.jsonl",
    "data/scan_312_uranus_robustness.jsonl",
    "data/scan_312_uranus_titania_umbriel.jsonl",
    "data/scan_312_uranus_umbriel_oberon_offset_sweep.jsonl",
    "data/scan_313_mars_deimos.jsonl",
    "data/scan_313_mars_phobos.jsonl",
    "data/scan_313_sun_jupiter_europa.jsonl",
    "data/scan_313_sun_jupiter_io.jsonl",
    "data/scan_285_saturn.jsonl",
    "data/scan_285_uranus.jsonl",
    "data/scan_298_galileo_veega.jsonl",
    "data/scan_309_low_thrust_em.jsonl",
    "data/scan_309_low_thrust_vem.jsonl",
    "data/scan_311_saturn_titan_iapetus.jsonl",
    "data/scan_311_saturn_titan_rhea.jsonl",
    "data/precursor_302_aldrin.jsonl",
    "data/precursor_302_s1l1.jsonl",
)


def _git_sha() -> str:
    try:
        return subprocess.check_output(["git", "rev-parse", "--short", "HEAD"], text=True).strip()
    except subprocess.SubprocessError:
        return "unknown"


def main(out_path: str, sources: tuple[str, ...]) -> None:
    repo_root = Path(__file__).resolve().parent.parent
    out_full = repo_root / out_path
    out_full.parent.mkdir(parents=True, exist_ok=True)
    with out_full.open("w") as out_fh:
        header = {
            "_meta": True,
            "task": "#325 Phase 1 Part B/C - topology-gate audit sweep",
            "ts_utc": dt.datetime.now(dt.UTC).isoformat(timespec="seconds"),
            "git_sha": _git_sha(),
            "n_sources": len(sources),
            "checkers": [
                "tulip_or_3d_periodic (#322 is_three_dimensional)",
                "periodic_orbit_closure (Radau cross-check)",
                "neimark_sacker (unit-circle + truly-complex)",
                "floquet_eig (informational only)",
            ],
            "discipline": (
                "Read-only on source JSONLs. NO catalogue writeback. NO novelty "
                "claims. Discrepancies are flagged for genome-side investigation, "
                "not retracted."
            ),
        }
        out_fh.write(json.dumps(header) + "\n")

        n_total = 0
        n_discrep = 0
        n_files_audited = 0
        n_files_missing = 0
        per_file_summary: list[dict] = []

        for src in sources:
            src_path = repo_root / src
            if not src_path.exists():
                n_files_missing += 1
                per_file_summary.append({"source_jsonl": src, "status": "missing"})
                continue
            try:
                findings = audit_topology(str(src_path))
            except (RuntimeError, ValueError) as exc:
                per_file_summary.append(
                    {
                        "source_jsonl": src,
                        "status": "error",
                        "error": str(exc),
                    }
                )
                continue
            n_files_audited += 1
            file_n = len(findings)
            file_disc = sum(1 for f in findings if f.discrepancy)
            n_total += file_n
            n_discrep += file_disc
            per_file_summary.append(
                {
                    "source_jsonl": src,
                    "status": "audited",
                    "n_findings": file_n,
                    "n_discrepancies": file_disc,
                }
            )
            for f in findings:
                payload = asdict(f)
                payload["source_jsonl"] = src
                out_fh.write(json.dumps(payload) + "\n")
            print(f"{src}: {file_n} findings ({file_disc} discrepancies)", flush=True)

        summary = {
            "_meta": True,
            "summary": True,
            "n_files_audited": n_files_audited,
            "n_files_missing": n_files_missing,
            "n_total_findings": n_total,
            "n_discrepancies": n_discrep,
            "per_file": per_file_summary,
        }
        out_fh.write(json.dumps(summary) + "\n")
        print(
            f"\n#325 audit summary: {n_files_audited}/{len(sources)} files audited, "
            f"{n_total} findings, {n_discrep} discrepancies. Output: {out_full}",
            flush=True,
        )


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", default="data/topology_audit_325.jsonl")
    parser.add_argument("--source", action="append", default=None)
    return parser.parse_args()


if __name__ == "__main__":
    args = _parse_args()
    sources = tuple(args.source) if args.source else DEFAULT_SOURCES
    main(args.out, sources)
