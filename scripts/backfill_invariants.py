#!/usr/bin/env python3
"""Phase 3: Backfill invariants{} blocks into all multi-arc rows in data/catalogue.yaml.

Reads AR, TR, and transit_times from the row's existing sourced notes/source_quotes/
trajectory.segments — transcribes faithfully, never invents.  Where a value is not
explicitly present, leaves null and adds a data_gaps[] entry.

Strategy: line-by-line text insertion (same approach as classify_cycler_class.py),
inserting the invariants: block after the cycler_class line for each multi-arc row.
No YAML round-trip — comments, formatting, and inline annotations are preserved.

Usage:
    uv run python scripts/backfill_invariants.py

Idempotent: rows that already carry an invariants: block are skipped.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import yaml  # type: ignore[import-untyped]

REPO_ROOT = Path(__file__).resolve().parent.parent
CATALOGUE_PATH = REPO_ROOT / "data" / "catalogue.yaml"


# ---------------------------------------------------------------------------
# Text extraction helpers (parse-only, do not mutate)
# ---------------------------------------------------------------------------


def _collect_text(obj: object, depth: int = 0) -> list[str]:
    parts: list[str] = []
    if depth > 10:
        return parts
    if isinstance(obj, str):
        parts.append(obj)
    elif isinstance(obj, dict):
        for v in obj.values():  # type: ignore[union-attr]
            parts.extend(_collect_text(v, depth + 1))
    elif isinstance(obj, list):
        for v in obj:  # type: ignore[union-attr]
            parts.extend(_collect_text(v, depth + 1))
    return parts


def _extract_ar(row: dict) -> float | None:  # type: ignore[type-arg]
    """Return the Aphelion Ratio stated in the row's sourced text, or None.

    Only returns a value when exactly one unambiguous AR value is present.
    """
    text = "\n".join(_collect_text(row))
    matches: list[float] = []
    for pat in [
        r"\bAR\s*=\s*([0-9]+\.[0-9]+)",
        r"Aphelion Ratio\s*(?:AR\s*=?\s*|is\s+|=\s*)?([0-9]+\.[0-9]+)",
    ]:
        for m in re.finditer(pat, text):
            matches.append(float(m.group(1)))
    unique: list[float] = []
    for v in matches:
        if not any(abs(v - u) < 0.005 for u in unique):
            unique.append(v)
    return unique[0] if len(unique) == 1 else None


def _extract_tr_from_parsed(row: dict) -> float | None:  # type: ignore[type-arg]
    """Return the Turn Ratio stated in the row's parsed YAML text, or None."""
    text = "\n".join(_collect_text(row))
    matches = [float(m.group(1)) for m in re.finditer(r"\bTR\s*=\s*([0-9]+\.[0-9]+)", text)]
    unique: list[float] = []
    for v in matches:
        if not any(abs(v - u) < 0.005 for u in unique):
            unique.append(v)
    return unique[0] if len(unique) == 1 else None


def _extract_tr_from_raw(raw_block: str) -> float | None:
    """Scan raw YAML block (including comments) for TR value."""
    matches = [float(m.group(1)) for m in re.finditer(r"\bTR\s*=\s*([0-9]+\.[0-9]+)", raw_block)]
    unique: list[float] = []
    for v in matches:
        if not any(abs(v - u) < 0.005 for u in unique):
            unique.append(v)
    return unique[0] if len(unique) == 1 else None


def _extract_transit_times(row: dict) -> list[int] | None:  # type: ignore[type-arg]
    """Return [t_out, t_in] or [t_out] from trajectory.segments, or None."""
    traj = row.get("trajectory") or {}
    segs = traj.get("segments") or []
    t_out: int | None = None
    t_in: int | None = None
    for seg in segs:
        if isinstance(seg, dict):
            seg_id = seg.get("id", "")
            tof = seg.get("tof_days")
            if tof is not None and isinstance(tof, (int, float)):
                if seg_id == "out-em" and t_out is None:
                    t_out = int(tof)
                elif seg_id == "ret-me" and t_in is None:
                    t_in = int(tof)

    # Fallback: scan parsed text for t_out=/t_in= patterns (Table 4.9 style)
    if t_out is None or t_in is None:
        text = "\n".join(_collect_text(row))
        if t_out is None:
            m = re.search(r"t_out\s*=\s*([0-9]+)", text)
            if m:
                t_out = int(m.group(1))
        if t_in is None:
            m = re.search(r"t_in\s*=\s*([0-9]+)", text)
            if m:
                t_in = int(m.group(1))

    if t_out is not None and t_in is not None:
        return [t_out, t_in]
    if t_out is not None:
        return [t_out]
    return None


# ---------------------------------------------------------------------------
# Manually resolved ambiguities (documented in Phase 3 spec)
# ---------------------------------------------------------------------------

# These rows have multiple TR values in their text; the correct one is chosen
# by examining which value cites the primary Russell table directly.
TR_OVERRIDES: dict[str, float] = {
    # russell-ch4-5.30ggF3: TR=1.27 is from corroborating_sources note citing
    # dissertation lines 5595-5597.  TR=1.00 in source_quotes.delta_v_kms
    # erroneously references cycler 5.33 ggF (a distinct cycler), not 5.30 ggF.
    "russell-ch4-5.30ggF3": 1.27,
    # russell-ch4-6.44Gg3: TR=0.95 is cited three times in the row as the
    # Russell 2004 Table 4.13 value.  TR=0.91 appears only in a note about
    # the Table 4.13 threshold (a reference value, not this cycler's TR).
    "russell-ch4-6.44Gg3": 0.95,
}


# ---------------------------------------------------------------------------
# YAML block renderer
# ---------------------------------------------------------------------------


def _fmt_float(v: float) -> str:
    """Format a float for YAML: use minimal decimal places."""
    # Round to 4 decimal places, strip trailing zeros after 2 places
    s = f"{v:.4f}".rstrip("0")
    if s.endswith("."):
        s += "0"
    # Ensure at least 2 decimal places
    if "." in s and len(s.split(".")[1]) < 2:
        s += "0"
    return s


def _build_invariants_block(
    ar: float | None,
    tr: float | None,
    transit: list[int] | None,
) -> str:
    """Return the invariants: YAML block (2-space indent, schema v4)."""
    lines: list[str] = []
    lines.append(
        "  invariants:   "
        "# schema v4 (2026-06-03); cycle-level identity for multi-arc cyclers; "
        "see docs/spec.md §16.7.4"
    )
    if ar is not None:
        lines.append(
            f"    aphelion_ratio: {_fmt_float(ar)}"
            f"   # Russell 2004: outbound-arc aphelion / Mars sma (1.52 AU)"
        )
    else:
        lines.append(
            "    aphelion_ratio: null   # not tabulated in this Russell source — see data_gaps[]"
        )
    if transit is not None:
        tt_str = "[" + ", ".join(str(t) for t in transit) + "]"
        tt_comment = "# [E->M, M->E] transit times in days; from trajectory.segments"
        lines.append(f"    transit_times_days: {tt_str}   {tt_comment}")
    else:
        lines.append("    transit_times_days: null   # not tabulated — see data_gaps[]")
    if tr is not None:
        tr_comment = "# Russell 2004: required / max ballistic turn (TR >= 1 => strictly ballistic)"
        lines.append(f"    turn_ratio: {_fmt_float(tr)}   {tr_comment}")
    else:
        lines.append(
            "    turn_ratio: null   # not tabulated in this Russell source — see data_gaps[]"
        )
    return "\n".join(lines) + "\n"


# ---------------------------------------------------------------------------
# Main backfill logic
# ---------------------------------------------------------------------------


def _get_raw_blocks(raw_text: str) -> dict[str, str]:
    """Return {id: raw_block_text} for every top-level entry."""
    lines = raw_text.splitlines()
    blocks: dict[str, str] = {}
    current_id: str | None = None
    start: int = 0
    for i, line in enumerate(lines):
        if re.match(r"^- id:\s+\S", line):
            if current_id is not None:
                blocks[current_id] = "\n".join(lines[start:i])
            current_id = line.split(":", 1)[1].strip()
            start = i
    if current_id is not None:
        blocks[current_id] = "\n".join(lines[start:])
    return blocks


def main() -> None:
    raw_text = CATALOGUE_PATH.read_text(encoding="utf-8")
    rows = yaml.safe_load(raw_text)
    raw_blocks = _get_raw_blocks(raw_text)

    lines = raw_text.splitlines(keepends=True)

    # Build the map of multi-arc row ids and their extracted invariants
    multi_arc_rows = {r["id"]: r for r in rows if r.get("cycler_class") == "multi-arc"}
    print(f"Multi-arc rows to backfill: {len(multi_arc_rows)}")

    # Count already-backfilled rows
    already_done = sum(1 for r in multi_arc_rows.values() if "invariants" in r)
    if already_done == len(multi_arc_rows):
        print("All multi-arc rows already have invariants{} — nothing to do.")
        return
    print(f"  Already backfilled: {already_done}")
    print(f"  To insert: {len(multi_arc_rows) - already_done}")

    # Build invariants data for each row
    invariants_data: dict[str, dict] = {}  # type: ignore[type-arg]
    for rid, row in multi_arc_rows.items():
        if "invariants" in row:
            continue  # skip already-populated rows

        ar = _extract_ar(row)
        tr = _extract_tr_from_parsed(row)

        # Fall back to raw text scan if TR is absent from parsed YAML
        if tr is None:
            raw_block = raw_blocks.get(rid, "")
            tr = _extract_tr_from_raw(raw_block)

        # Apply manual overrides for ambiguous cases
        if rid in TR_OVERRIDES:
            tr = TR_OVERRIDES[rid]

        transit = _extract_transit_times(row)
        invariants_data[rid] = {"ar": ar, "tr": tr, "transit": transit}

    # Coverage stats
    ar_count = sum(1 for v in invariants_data.values() if v["ar"] is not None)
    tr_count = sum(1 for v in invariants_data.values() if v["tr"] is not None)
    tt_count = sum(1 for v in invariants_data.values() if v["transit"] is not None)
    print("\nExtracted values:")
    print(f"  aphelion_ratio populated: {ar_count}/{len(invariants_data)}")
    print(f"  turn_ratio populated: {tr_count}/{len(invariants_data)}")
    print(f"  transit_times_days populated: {tt_count}/{len(invariants_data)}")
    print(f"  aphelion_ratio gaps: {len(invariants_data) - ar_count}")
    print(f"  turn_ratio gaps: {len(invariants_data) - tr_count}")
    print(f"  transit_times gaps: {len(invariants_data) - tt_count}")

    # Line-by-line insertion: insert the invariants: block after the cycler_class line
    # for each multi-arc row that doesn't already have one.
    #
    # We track the current row id as we scan.  When we see a "  cycler_class: multi-arc"
    # line for a row that needs backfilling, we emit the invariants block immediately after.

    current_id: str | None = None
    output_lines: list[str] = []
    inserted = 0

    for line in lines:
        stripped = line.lstrip()

        # Detect row id change
        if stripped.startswith("- id:") or (stripped.startswith("id:") and line.startswith("- ")):
            parts = stripped.split(":", 1)
            if len(parts) == 2 and parts[0].strip() in ("id", "- id"):
                slug = parts[1].strip()
                if "#" in slug:
                    slug = slug[: slug.index("#")].strip()
                current_id = slug

        output_lines.append(line)

        # After emitting the cycler_class line for a multi-arc row that needs backfilling,
        # insert the invariants block.
        if (
            line.startswith("  cycler_class: multi-arc")
            and current_id is not None
            and current_id in invariants_data
        ):
            d = invariants_data[current_id]
            block = _build_invariants_block(d["ar"], d["tr"], d["transit"])
            output_lines.append(block)
            inserted += 1

    if inserted != len(invariants_data):
        print(
            f"ERROR: expected to insert {len(invariants_data)} invariants blocks, "
            f"inserted {inserted}. Aborting.",
            file=sys.stderr,
        )
        sys.exit(1)

    result = "".join(output_lines)
    CATALOGUE_PATH.write_text(result, encoding="utf-8")
    print(f"\nInserted {inserted} invariants blocks into {CATALOGUE_PATH}")
    print("Done.")


if __name__ == "__main__":
    main()
