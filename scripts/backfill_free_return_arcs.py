#!/usr/bin/env python3
"""Backfill free_return_arcs into multi-arc rows whose notes carry explicit
Russell leg descriptors (spec §16.7.7, schema v4.1, 2026-06-03).

Source: Russell 2004 dissertation Tables 4.9-4.13 (pp.127-134).
Descriptor notation explained at p.126:
- "the first number in the parenthesis is the time of flight in years"
- "the number following the colon in the full-rev strings represent the
  number of revolutions by the spacecraft" (M:N resonant orbit)
- uppercase letter = designated transit leg
- letters: g/G = generic, h/H = half-rev, f/F = full-rev

Strategy: line-by-line text insertion to preserve all comments, formatting,
and inline annotations. Each targeted entry has a ``  delta_v_kms:`` line;
we insert ``  free_return_arcs:`` block immediately before it.
No YAML round-trip (avoids comment stripping / reflowing).

Usage:
    uv run python scripts/backfill_free_return_arcs.py
"""

from __future__ import annotations

import sys
from pathlib import Path

# ---------------------------------------------------------------------------
# Source data: descriptor tokens per entry ID.
#
# Source for each entry: Russell 2004 dissertation, as noted.
# Entries without an explicit descriptor are intentionally omitted (gap, not error).
# ---------------------------------------------------------------------------

# Each value is a list of descriptor token dicts with:
#   raw_descriptor: verbatim Russell token string
#   arc_type: "generic" | "half-rev" | "full-rev"  (derived from token letter)
#   resonance: M:N string for full-rev, null otherwise
#   tof_years: float (from first param) for g/h; null for f/F (determined by resonance)
#
# Source: Russell 2004 dissertation Tables 4.9-4.12 pp.127-132.

FREE_RETURN_ARCS: dict[str, list[dict]] = {
    # -------------------------------------------------------------------
    # Table 4.9 row 1: 4.99 gG2 (S1L1, v_inf E=4.99, v_inf M=5.10, t=150d)
    # Source: Russell 2004 Table 4.9 p.127, footnote a="documented in Ref.15"
    # -------------------------------------------------------------------
    "russell-ch4-4.991gG2": [
        {
            "arc_type": "generic",
            "resonance": None,
            "tof_years": 1.4612,
            "raw_descriptor": "g(1.4612,526.02,Ll)",
        },
        {
            "arc_type": "generic",
            "resonance": None,
            "tof_years": 2.8096,
            "raw_descriptor": "G(2.8096,651.46,U)",
        },
    ],
    # -------------------------------------------------------------------
    # Table 4.9 row 2: 8.05 gGf2 (new cycler, v_inf E=8.05, t=93d)
    # Source: Russell 2004 Table 4.9 p.127 row 2 (no footnote = undocumented)
    # -------------------------------------------------------------------
    "russell-ch4-8.049gGf2": [
        {
            "arc_type": "generic",
            "resonance": None,
            "tof_years": 1.4951,
            "raw_descriptor": "g(1.4951,538.24,Ll)",
        },
        {
            "arc_type": "generic",
            "resonance": None,
            "tof_years": 1.7757,
            "raw_descriptor": "G(1.7757,279.24,U)",
        },
        {
            "arc_type": "full-rev",
            "resonance": "1:1",
            "tof_years": None,
            "raw_descriptor": "f(1:1,74.468,-180.000)",
        },
    ],
    # -------------------------------------------------------------------
    # Table 4.9 row 4: 9.35 Gg2 (v_inf E=9.35, t=85d)
    # Source: Russell 2004 Table 4.9 p.127 row 4
    # -------------------------------------------------------------------
    "russell-ch4-9.353Gg2": [
        {
            "arc_type": "generic",
            "resonance": None,
            "tof_years": 1.7238,
            "raw_descriptor": "G(1.7238,260.58,U)",
        },
        {
            "arc_type": "generic",
            "resonance": None,
            "tof_years": 2.5469,
            "raw_descriptor": "g(2.5469,916.9,L)",
        },
    ],
    # -------------------------------------------------------------------
    # Table 4.10 low-energy bold cluster: 3.64 gGg3 (v_inf E=3.64, t=175d)
    # Source: Russell 2004 Table 4.10 p.128 notes line 5321.
    # NOTE: The shorthand label says 'gGg' but the explicit descriptor
    # shows f(1:1,...) as the third leg. The descriptor takes precedence.
    # Source: catalogue notes field "Leg descriptor: g(2.4845,894.42,Ll)
    # G(2.9217,691.79,U) f(1:1,82.995,-180.000)"
    # -------------------------------------------------------------------
    "russell-ch4-3.64gGg3": [
        {
            "arc_type": "generic",
            "resonance": None,
            "tof_years": 2.4845,
            "raw_descriptor": "g(2.4845,894.42,Ll)",
        },
        {
            "arc_type": "generic",
            "resonance": None,
            "tof_years": 2.9217,
            "raw_descriptor": "G(2.9217,691.79,U)",
        },
        {
            "arc_type": "full-rev",
            "resonance": "1:1",
            "tof_years": None,
            "raw_descriptor": "f(1:1,82.995,-180.000)",
        },
    ],
    # -------------------------------------------------------------------
    # Table 4.10 bold: 3.78 Gg3 (v_inf E=3.78, favorable turning angles)
    # Source: Russell 2004 Table 4.10 p.128 (bolded)
    # -------------------------------------------------------------------
    "russell-ch4-3.78Gg3": [
        {
            "arc_type": "generic",
            "resonance": None,
            "tof_years": 2.9043,
            "raw_descriptor": "G(2.9043,685.56,U)",
        },
        {
            "arc_type": "generic",
            "resonance": None,
            "tof_years": 3.5018,
            "raw_descriptor": "g(3.5018,1260.65,L)",
        },
    ],
    # -------------------------------------------------------------------
    # Table 4.11 5.30gGf: 5.30 gGf3 (v_inf E=5.30, 3:2 resonant, t=118d)
    # Source: Russell 2004 Table 4.11 p.129. From notes: "Leg descriptor:
    # g(1.4646,527.25,Ll) G(1.9416,338.97,U) f(3:2,82.487,118.851)"
    # -------------------------------------------------------------------
    "russell-ch4-5.30gGf3": [
        {
            "arc_type": "generic",
            "resonance": None,
            "tof_years": 1.4646,
            "raw_descriptor": "g(1.4646,527.25,Ll)",
        },
        {
            "arc_type": "generic",
            "resonance": None,
            "tof_years": 1.9416,
            "raw_descriptor": "G(1.9416,338.97,U)",
        },
        {
            "arc_type": "full-rev",
            "resonance": "3:2",
            "tof_years": None,
            "raw_descriptor": "f(3:2,82.487,118.851)",
        },
    ],
    # -------------------------------------------------------------------
    # Table 4.11 9.94Gg: 9.94 Gg3 (shortest transit 82d, v_inf E=9.94)
    # Source: Russell 2004 Table 4.11 p.129 (bolded, 9.94Gg)
    # -------------------------------------------------------------------
    "russell-ch4-9.94Gg3": [
        {
            "arc_type": "generic",
            "resonance": None,
            "tof_years": 1.7025,
            "raw_descriptor": "G(1.7025,252.9,U)",
        },
        {
            "arc_type": "generic",
            "resonance": None,
            "tof_years": 4.7037,
            "raw_descriptor": "g(4.7037,2053.31,Ls)",
        },
    ],
    # -------------------------------------------------------------------
    # Table 4.12 row 3.66gfF: 3.66 gfF3 (v_inf E=3.66, t=173d)
    # Source: Russell 2004 Table 4.12 p.132 (bolded)
    # -------------------------------------------------------------------
    "russell-ch4-3.66gfF3": [
        {
            "arc_type": "generic",
            "resonance": None,
            "tof_years": 2.4062,
            "raw_descriptor": "g(2.4062,866.21,Ls)",
        },
        {
            "arc_type": "full-rev",
            "resonance": "1:1",
            "tof_years": None,
            "raw_descriptor": "f(1:1,82.955,87.388)",
        },
        {
            "arc_type": "full-rev",
            "resonance": "3:2",
            "tof_years": None,
            "raw_descriptor": "F(3:2,87.27,0.000)",
        },
    ],
    # -------------------------------------------------------------------
    # Table 4.12 row 5.30ggF: 5.30 ggF3 (v_inf E=5.30, 3:2 resonant, t=207d)
    # Source: Russell 2004 Table 4.12 p.132 (bolded, note 64-day difference
    # in inbound/outbound times)
    # -------------------------------------------------------------------
    "russell-ch4-5.30ggF3": [
        {
            "arc_type": "generic",
            "resonance": None,
            "tof_years": 1.4646,
            "raw_descriptor": "g(1.4646,527.25,Ll)",
        },
        {
            "arc_type": "generic",
            "resonance": None,
            "tof_years": 1.9416,
            "raw_descriptor": "g(1.9416,338.97,U)",
        },
        {
            "arc_type": "full-rev",
            "resonance": "3:2",
            "tof_years": None,
            "raw_descriptor": "F(3:2,82.487,180.000)",
        },
    ],
    # -------------------------------------------------------------------
    # Table 4.12 row 5.75ggF: 5.75 ggF3 (v_inf E=5.75, 2:1 resonant, t=111d)
    # Source: Russell 2004 Table 4.12 p.132 (bolded, short transfer time)
    # -------------------------------------------------------------------
    "russell-ch4-5.75ggF3": [
        {
            "arc_type": "generic",
            "resonance": None,
            "tof_years": 1.8987,
            "raw_descriptor": "g(1.8987,323.54,U)",
        },
        {
            "arc_type": "generic",
            "resonance": None,
            "tof_years": 2.5074,
            "raw_descriptor": "g(2.5074,902.67,L)",
        },
        {
            "arc_type": "full-rev",
            "resonance": "2:1",
            "tof_years": None,
            "raw_descriptor": "F(2:1,85.196,0.000)",
        },
    ],
    # -------------------------------------------------------------------
    # Table 4.13 (near-ballistic): 6.44 Gg3 (v_inf E=6.44, very long ToF)
    # Source: Russell 2004 Table 4.13 p.134 (bolded near-ballistic)
    # -------------------------------------------------------------------
    "russell-ch4-6.44Gg3": [
        {
            "arc_type": "generic",
            "resonance": None,
            "tof_years": 2.087,
            "raw_descriptor": "g(2.087,1111.33,L)",
        },
        {
            "arc_type": "generic",
            "resonance": None,
            "tof_years": 4.3191,
            "raw_descriptor": "G(4.3191,1194.88,L)",
        },
    ],
    # -------------------------------------------------------------------
    # mcconaghy-2006-em-k2: same as russell-ch4-4.991gG2 (same underlying cycler)
    # The descriptor is explicitly quoted in the notes field:
    # "leg descriptor 'g(1.4612,526.02,Ll) G(2.8096,651.46,U)'"
    # Source: catalogue notes §, attributed to Russell 2004 Table 4.9 row 1
    # -------------------------------------------------------------------
    "mcconaghy-2006-em-k2": [
        {
            "arc_type": "generic",
            "resonance": None,
            "tof_years": 1.4612,
            "raw_descriptor": "g(1.4612,526.02,Ll)",
        },
        {
            "arc_type": "generic",
            "resonance": None,
            "tof_years": 2.8096,
            "raw_descriptor": "G(2.8096,651.46,U)",
        },
    ],
}

# Entries that are multi-arc but have no extractable descriptor (data gaps, not errors):
# russell-ch4-8.165Gfh-f2  — G(1.7708,277.48,U) partially mentioned, f/h/f not spelled out
# russell-ch4-3.77Gh3      — h(3.5,3,U,+-2.022) is a non-standard 4-param format
# russell-ch4-5.66Gfh3     — no complete descriptor in notes
# All russell-ocampo-* entries — Russell Ch3 tables don't carry leg descriptors
# These are deliberately left null.

ENTRIES_WITH_DESCRIPTORS = set(FREE_RETURN_ARCS.keys())
assert len(ENTRIES_WITH_DESCRIPTORS) == 12


def _format_arc_block(arcs: list[dict], indent: str = "  ") -> str:
    """Format a free_return_arcs YAML block with proper indentation.

    Produces a block like:
        free_return_arcs:   # schema v4.1 (spec §16.7.7); Russell Earth-Earth arc descriptors
          - arc_type: generic
            resonance: null
            tof_years: 1.4612
            raw_descriptor: "g(1.4612,526.02,Ll)"
    """
    lines = [
        f"{indent}free_return_arcs:   "
        "# schema v4.1 (spec §16.7.7); Russell Earth-to-Earth free-return arc descriptors\n"
    ]
    for arc in arcs:
        tof = arc["tof_years"]
        res = arc["resonance"]
        raw = arc["raw_descriptor"]
        arc_type = arc["arc_type"]

        tof_str = str(tof) if tof is not None else "null"
        res_str = f'"{res}"' if res is not None else "null"

        lines.append(f"{indent}  - arc_type: {arc_type}\n")
        lines.append(f"{indent}    resonance: {res_str}\n")
        lines.append(f"{indent}    tof_years: {tof_str}\n")
        lines.append(f'{indent}    raw_descriptor: "{raw}"\n')
    return "".join(lines)


def main() -> None:
    catalogue = Path(__file__).resolve().parent.parent / "data" / "catalogue.yaml"
    text = catalogue.read_text(encoding="utf-8")
    lines = text.splitlines(keepends=True)

    current_id: str | None = None
    output_lines: list[str] = []
    inserted_count = 0

    for line in lines:
        stripped = line.lstrip()

        # Detect entry start: "- id: <slug>"
        if stripped.startswith("- id:") or (stripped.startswith("id:") and line.startswith("- ")):
            parts = stripped.split(":", 1)
            if len(parts) == 2 and parts[0].strip() in ("id", "- id"):
                slug = parts[1].strip()
                if "#" in slug:
                    slug = slug[: slug.index("#")].strip()
                current_id = slug

        # Before emitting delta_v_kms, insert free_return_arcs block if applicable.
        # The delta_v_kms line is always the first top-level field after invariants/cycler_class.
        if (
            line.startswith("  delta_v_kms:")
            and current_id is not None
            and current_id in FREE_RETURN_ARCS
        ):
            arcs = FREE_RETURN_ARCS[current_id]
            block = _format_arc_block(arcs, indent="  ")
            output_lines.append(block)
            inserted_count += 1

        output_lines.append(line)

    expected = len(ENTRIES_WITH_DESCRIPTORS)
    if inserted_count != expected:
        print(
            f"ERROR: expected to insert {expected} blocks, inserted {inserted_count}. "
            f"Check entry IDs and file structure.",
            file=sys.stderr,
        )
        sys.exit(1)

    result = "".join(output_lines)
    catalogue.write_text(result, encoding="utf-8")
    print(f"Backfilled free_return_arcs for {inserted_count} entries in {catalogue}")


if __name__ == "__main__":
    main()
