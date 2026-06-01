"""Render a cross-reference markdown table of the seed catalogue.

Usage:
    uv run --with pyyaml python scripts/render-catalogue.py        # markdown
    uv run --with pyyaml python scripts/render-catalogue.py --csv  # spreadsheet ingest

Reads ``data/catalogue.yaml`` at the repo root and emits a
cross-reference of all entries to stdout. The committed catalogue is
YAML; this script is the on-demand renderer so no derivative markdown
table needs to be kept in sync.
"""

from __future__ import annotations

import argparse
import csv
import io
import sys
from pathlib import Path
from typing import Any

import yaml


def _legs_of(entry: dict[str, Any]) -> list[dict[str, Any]]:
    """Schema-v3 leg source: prefer trajectory.segments, fall back to legs[].

    Mirrors ``cyclerfinder.data.catalog._segments_as_legs`` — segments reuse
    the ``tof_days`` key, so leg-aware checks work on either form.
    """
    segments = (entry.get("trajectory") or {}).get("segments")
    if segments:
        return list(segments)
    return list(entry.get("legs") or [])


def _fully_defined(entry: dict[str, Any]) -> bool:
    """Raw-dict mirror of ``CatalogueEntry.fully_defined`` (spec §16.6.4).

    The orbit is completely specified: all core fields present and no
    acknowledged known-unknown (``data_gaps[]``).
    """
    if entry.get("data_gaps"):
        return False
    elems = entry.get("orbit_elements") or {}
    if elems.get("a_au") is None or elems.get("e") is None:
        return False
    vinfs = entry.get("vinf_kms_at_encounters") or []
    if not vinfs or any(v.get("vinf_kms") is None for v in vinfs):
        return False
    legs = _legs_of(entry)
    return bool(legs) and all(leg.get("tof_days") is not None for leg in legs)


def _data_completeness(entry: dict[str, Any]) -> str:
    elems = entry.get("orbit_elements") or {}
    vinfs = entry.get("vinf_kms_at_encounters") or []
    legs = _legs_of(entry)
    has_elems = elems.get("a_au") is not None and elems.get("e") is not None
    has_vinfs = bool(vinfs) and any(v.get("vinf_kms") is not None for v in vinfs)
    has_legs = bool(legs) and any(leg.get("tof_days") is not None for leg in legs)
    score = sum([has_elems, has_vinfs, has_legs])
    if score == 3:
        return "full"
    if score >= 1:
        return "partial"
    return "citation-only"


def _max_vinf(entry: dict[str, Any]) -> str:
    vinfs = entry.get("vinf_kms_at_encounters") or []
    vals = [v.get("vinf_kms") for v in vinfs if v.get("vinf_kms") is not None]
    return f"{max(vals):.1f}" if vals else "—"


def _primary_citation(entry: dict[str, Any]) -> str:
    fp = entry.get("first_published") or {}
    authors = fp.get("authors") or []
    year = fp.get("year") or "?"
    first = authors[0].split(",")[0] if authors else "?"
    et_al = " et al." if len(authors) > 1 else ""
    doi = fp.get("doi")
    cite = f"{first}{et_al} {year}"
    if doi:
        cite += f" ({doi})"
    return cite


def _period_str(entry: dict[str, Any]) -> str:
    period = entry.get("period") or {}
    years = period.get("years")
    k = period.get("k")
    if years and k:
        return f"{years} yr (k={k})"
    if years:
        return f"{years} yr"
    return "—"


def render_md(entries: list[dict[str, Any]]) -> str:
    rows = [
        "| # | id | primary | regime | bodies | period "
        "| max V∞ (km/s) | data | defined | citation |",
        "|---|----|---------|--------|--------|--------"
        "|---------------|------|---------|----------|",
    ]
    for i, e in enumerate(entries, 1):
        rows.append(
            f"| {i} | `{e['id']}` | {e.get('primary', 'Sun')} | "
            f"{e.get('trajectory_regime', 'ballistic')} | "
            f"{','.join(e.get('bodies') or [])} | {_period_str(e)} | "
            f"{_max_vinf(e)} | {_data_completeness(e)} | "
            f"{'yes' if _fully_defined(e) else 'no'} | {_primary_citation(e)} |"
        )
    return "\n".join(rows) + "\n"


def render_csv(entries: list[dict[str, Any]]) -> str:
    buf = io.StringIO()
    w = csv.writer(buf)
    w.writerow(
        [
            "id",
            "primary",
            "trajectory_regime",
            "bodies",
            "period_years",
            "period_k",
            "max_vinf_kms",
            "data_completeness",
            "fully_defined",
            "primary_citation",
        ]
    )
    for e in entries:
        period = e.get("period") or {}
        w.writerow(
            [
                e["id"],
                e.get("primary", "Sun"),
                e.get("trajectory_regime", "ballistic"),
                ",".join(e.get("bodies") or []),
                period.get("years"),
                period.get("k"),
                _max_vinf(e),
                _data_completeness(e),
                _fully_defined(e),
                _primary_citation(e),
            ]
        )
    return buf.getvalue()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0] if __doc__ else "")
    parser.add_argument("--csv", action="store_true", help="emit CSV instead of markdown")
    args = parser.parse_args()
    repo_root = Path(__file__).resolve().parent.parent
    catalogue_path = repo_root / "data" / "catalogue.yaml"
    with catalogue_path.open() as f:
        entries = yaml.safe_load(f)
    if not isinstance(entries, list):
        raise SystemExit(f"Expected a YAML list at {catalogue_path}; got {type(entries).__name__}")
    out = render_csv(entries) if args.csv else render_md(entries)
    sys.stdout.write(out)
    n_defined = sum(1 for e in entries if _fully_defined(e))
    sys.stderr.write(f"fully-defined orbits: {n_defined} / {len(entries)}\n")


if __name__ == "__main__":
    main()
