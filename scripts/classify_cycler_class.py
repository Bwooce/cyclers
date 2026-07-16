#!/usr/bin/env python3
"""Tag every row in data/catalogue.yaml with cycler_class.

Classification rules (source: docs/notes/multi-arc-classification.md §3 and §9):
- non-keplerian: the 6 rows listed in §3 (primary != Sun).
- multi-arc: exactly the 199 ids in MULTI_ARC_ALLOWLIST (§9).
- single-ellipse: all remaining 28 rows.

Strategy: line-by-line text insertion to preserve all comments, formatting,
and inline annotations. Each entry has a ``  model_assumption:`` line; we
insert ``  cycler_class: <value>  # schema v4 (2026-06-03)`` immediately after
it. No YAML round-trip (avoids comment stripping / reflowing).

Usage:
    uv run python scripts/classify_cycler_class.py
"""

from __future__ import annotations

import sys
from pathlib import Path

# ---------------------------------------------------------------------------
# Classification data (from docs/notes/multi-arc-classification.md §3 and §9)
# ---------------------------------------------------------------------------

MULTI_ARC_IDS: frozenset[str] = frozenset(
    [
        # russell-ocampo family (184 rows) - Russell 2004 Ch3, Tables 3.4-3.11
        "russell-ocampo-2.1.1+2-case2",
        "russell-ocampo-2.3.1+1-case3",
        "russell-ocampo-2.5.1+0",
        "russell-ocampo-3.1.1+3",
        "russell-ocampo-3.1.2+1",
        "russell-ocampo-3.3.1+2",
        "russell-ocampo-3.5.1+1",
        "russell-ocampo-3.5.2+0",
        "russell-ocampo-3.7.1+1",
        "russell-ocampo-3.9.1+0",
        "russell-ocampo-4.0.3+1",
        "russell-ocampo-4.1.1-5",
        "russell-ocampo-4.10.1+2",
        "russell-ocampo-4.11.1-2",
        "russell-ocampo-4.12.1+1",
        "russell-ocampo-4.13.1-1",
        "russell-ocampo-4.14.1+0",
        "russell-ocampo-4.14.1-1",
        "russell-ocampo-4.3.1-4",
        "russell-ocampo-4.3.1-5",
        "russell-ocampo-4.5.1-3",
        "russell-ocampo-4.5.1-4",
        "russell-ocampo-4.5.2-2",
        "russell-ocampo-4.5.3-1",
        "russell-ocampo-4.7.1-3",
        "russell-ocampo-4.9.1-2",
        "russell-ocampo-4.9.2-1",
        "russell-ocampo-5.1.1-7",
        "russell-ocampo-5.1.2-3",
        "russell-ocampo-5.1.5-1",
        "russell-ocampo-5.2.1-7",
        "russell-ocampo-5.2.2+2",
        "russell-ocampo-5.2.5+0",
        "russell-ocampo-5.3.1-6",
        "russell-ocampo-5.3.1-7",
        "russell-ocampo-5.3.3-2",
        "russell-ocampo-5.4.1+5",
        "russell-ocampo-5.4.1+6",
        "russell-ocampo-5.4.3+1",
        "russell-ocampo-5.5.1-4",
        "russell-ocampo-5.5.1-5",
        "russell-ocampo-5.5.1-6",
        "russell-ocampo-5.5.2-2",
        "russell-ocampo-5.5.2-3",
        "russell-ocampo-5.5.4-1",
        "russell-ocampo-5.6.1+3",
        "russell-ocampo-5.6.1+4",
        "russell-ocampo-5.6.1+5",
        "russell-ocampo-5.6.2+1",
        "russell-ocampo-5.6.2+2",
        "russell-ocampo-5.6.4+0",
        "russell-ocampo-5.7.1-3",
        "russell-ocampo-5.7.1-4",
        "russell-ocampo-5.7.1-5",
        "russell-ocampo-5.8.1+2",
        "russell-ocampo-5.8.1+3",
        "russell-ocampo-5.8.1+4",
        "russell-ocampo-5.9.1-2",
        "russell-ocampo-5.9.1-3",
        "russell-ocampo-5.9.1-4",
        "russell-ocampo-5.9.2+1",
        "russell-ocampo-5.9.2-2",
        "russell-ocampo-5.9.3-1",
        "russell-ocampo-5.10.1+2",
        "russell-ocampo-5.10.1+3",
        "russell-ocampo-5.10.2+0",
        "russell-ocampo-5.10.2+1",
        "russell-ocampo-5.10.3+0",
        "russell-ocampo-5.11.1-2",
        "russell-ocampo-5.11.1-3",
        "russell-ocampo-5.11.2+1",
        "russell-ocampo-5.12.1+1",
        "russell-ocampo-5.12.1+2",
        "russell-ocampo-5.13.1-2",
        "russell-ocampo-5.13.1-3",
        "russell-ocampo-5.13.2-1",
        "russell-ocampo-5.14.1+1",
        "russell-ocampo-5.14.1+2",
        "russell-ocampo-5.14.2+0",
        "russell-ocampo-5.15.1-1",
        "russell-ocampo-5.15.1-2",
        "russell-ocampo-5.16.1+0",
        "russell-ocampo-5.16.1+1",
        "russell-ocampo-5.17.1-1",
        "russell-ocampo-5.18.1+0",
        "russell-ocampo-6.0.1+6d",
        "russell-ocampo-6.0.1+7c",
        "russell-ocampo-6.0.1+8b",
        "russell-ocampo-6.0.1+9a",
        "russell-ocampo-6.1.2-4",
        "russell-ocampo-6.1.3-3",
        "russell-ocampo-6.1.4-2",
        "russell-ocampo-6.1.6-1",
        "russell-ocampo-6.2.1+6",
        "russell-ocampo-6.2.1+7",
        "russell-ocampo-6.2.1+8",
        "russell-ocampo-6.2.2+2",
        "russell-ocampo-6.2.2+3",
        "russell-ocampo-6.2.3+1",
        "russell-ocampo-6.2.3+2",
        "russell-ocampo-6.2.4+1",
        "russell-ocampo-6.2.6+0",
        "russell-ocampo-6.3.1-9",
        "russell-ocampo-6.3.4+1",
        "russell-ocampo-6.4.1+4",
        "russell-ocampo-6.4.1+5",
        "russell-ocampo-6.4.1+6",
        "russell-ocampo-6.4.1+7",
        "russell-ocampo-6.5.1-6",
        "russell-ocampo-6.5.1-7",
        "russell-ocampo-6.5.1-8",
        "russell-ocampo-6.5.5-1",
        "russell-ocampo-6.6.1+3",
        "russell-ocampo-6.6.1+4",
        "russell-ocampo-6.6.1+5",
        "russell-ocampo-6.6.1-6",
        "russell-ocampo-6.6.2+1",
        "russell-ocampo-6.6.2+2",
        "russell-ocampo-6.6.5+0",
        "russell-ocampo-6.7.1-6",
        "russell-ocampo-6.7.1-7",
        "russell-ocampo-6.7.2+3",
        "russell-ocampo-6.7.3-2",
        "russell-ocampo-6.7.5+0",
        "russell-ocampo-6.8.1+2",
        "russell-ocampo-6.8.1+3",
        "russell-ocampo-6.8.1+4",
        "russell-ocampo-6.8.1+5",
        "russell-ocampo-6.8.1+6",
        "russell-ocampo-6.8.3+1",
        "russell-ocampo-6.9.1-4",
        "russell-ocampo-6.9.1-5",
        "russell-ocampo-6.9.1-6",
        "russell-ocampo-6.9.2-2",
        "russell-ocampo-6.9.2-3",
        "russell-ocampo-6.9.4-1",
        "russell-ocampo-6.10.1+2",
        "russell-ocampo-6.10.1+3",
        "russell-ocampo-6.10.1+4",
        "russell-ocampo-6.10.1+5",
        "russell-ocampo-6.10.2+1",
        "russell-ocampo-6.10.2+2",
        "russell-ocampo-6.10.4+0",
        "russell-ocampo-6.11.1-4",
        "russell-ocampo-6.11.1-5",
        "russell-ocampo-6.11.2+2",
        "russell-ocampo-6.12.1+2",
        "russell-ocampo-6.12.1+3",
        "russell-ocampo-6.12.1+4",
        "russell-ocampo-6.13.1-3",
        "russell-ocampo-6.13.1-4",
        "russell-ocampo-6.13.1-5",
        "russell-ocampo-6.13.1+5",
        "russell-ocampo-6.13.2-2",
        "russell-ocampo-6.13.3-1",
        "russell-ocampo-6.14.1+1",
        "russell-ocampo-6.14.1+2",
        "russell-ocampo-6.14.1+3",
        "russell-ocampo-6.14.2+0",
        "russell-ocampo-6.14.2+1",
        "russell-ocampo-6.14.3+0",
        "russell-ocampo-6.15.1-2",
        "russell-ocampo-6.15.1-3",
        "russell-ocampo-6.15.1-4",
        "russell-ocampo-6.15.1+4",
        "russell-ocampo-6.15.2+1",
        "russell-ocampo-6.16.1+1",
        "russell-ocampo-6.16.1+2",
        "russell-ocampo-6.17.1-2",
        "russell-ocampo-6.17.1-3",
        "russell-ocampo-6.17.1+3",
        "russell-ocampo-6.17.2-1",
        "russell-ocampo-6.18.1+1",
        "russell-ocampo-6.18.1+2",
        "russell-ocampo-6.18.2+0",
        "russell-ocampo-6.19.1-2",
        "russell-ocampo-6.19.1+2",
        "russell-ocampo-6.19.2+0",
        "russell-ocampo-6.20.1-4",
        "russell-ocampo-6.20.1+0",
        "russell-ocampo-6.20.1+1",
        "russell-ocampo-6.21.1-1",
        "russell-ocampo-6.21.1+1",
        "russell-ocampo-6.22.1+0",
        # russell-ch4 family (14 rows) - Russell 2004 Ch4, Tables 4.9-4.13
        "russell-ch4-4.991gG2",
        "russell-ch4-8.049gGf2",
        "russell-ch4-8.165Gfh-f2",
        "russell-ch4-9.353Gg2",
        "russell-ch4-3.64gGg3",
        "russell-ch4-3.77Gh3",
        "russell-ch4-3.78Gg3",
        "russell-ch4-5.30gGf3",
        "russell-ch4-5.66Gfh3",
        "russell-ch4-9.94Gg3",
        "russell-ch4-3.66gfF3",
        "russell-ch4-5.30ggF3",
        "russell-ch4-5.75ggF3",
        "russell-ch4-6.44Gg3",
        # resolved §7 — same physical cycler as russell-ch4-4.991gG2
        "mcconaghy-2006-em-k2",
    ]
)

NON_KEPLERIAN_IDS: frozenset[str] = frozenset(
    [
        "arenstorf-em-figure8-1963",
        "genova-aldrin-2015-em-3petal-cycler",
        "wittal-2022-em-cycler-family",
        "hernandez-2017-jovian-ieg-triple-family",
        "russell-strange-2009-jovian-multimoon-family",
        "russell-strange-2009-saturnian-multimoon-family",
    ]
)

assert len(MULTI_ARC_IDS) == 199
assert len(NON_KEPLERIAN_IDS) == 6


def classify(row_id: str) -> str:
    """Return the cycler_class for a given row id."""
    if row_id in NON_KEPLERIAN_IDS:
        return "non-keplerian"
    if row_id in MULTI_ARC_IDS:
        return "multi-arc"
    return "single-ellipse"


def main() -> None:
    catalogue = Path(__file__).resolve().parent.parent / "data" / "catalogue.yaml"
    text = catalogue.read_text(encoding="utf-8")
    lines = text.splitlines(keepends=True)

    # We scan line-by-line, tracking the current entry id so we know what
    # class to assign when we hit a ``  model_assumption:`` line.
    #
    # The catalogue has one ``- id: <slug>`` line per entry at the top level
    # (the YAML list item marker). We track the most recently seen id.

    current_id: str | None = None
    output_lines: list[str] = []
    tagged_count = 0

    for line in lines:
        # Detect entry start: "- id: <slug>"
        stripped = line.lstrip()
        if stripped.startswith("- id:") or stripped.startswith("id:"):
            # The canonical form is "- id: <slug>" for top-level list items.
            # The form "  id: <slug>" does not occur at entry level (ids are
            # always the list bullet). Grab the slug.
            parts = stripped.split(":", 1)
            if len(parts) == 2 and parts[0].strip() in ("id", "- id"):
                slug = parts[1].strip()
                # Strip inline comment if any
                if "#" in slug:
                    slug = slug[: slug.index("#")].strip()
                current_id = slug

        output_lines.append(line)

        # After emitting the model_assumption line, insert cycler_class.
        # Only do this if the line is the entry-level model_assumption
        # (2-space indent, matching "  model_assumption:").
        if line.startswith("  model_assumption:") and current_id is not None:
            cls = classify(current_id)
            # Mirror the comment style of other schema-v* fields in the catalogue.
            new_line = (
                f"  cycler_class: {cls}"
                f"   # schema v4 (2026-06-03); single-ellipse | multi-arc | non-keplerian;"
                f" see docs/notes/multi-arc-classification.md\n"
            )
            output_lines.append(new_line)
            tagged_count += 1

    if tagged_count != 233:
        print(
            f"ERROR: expected to tag 233 entries, tagged {tagged_count}. Aborting.",
            file=sys.stderr,
        )
        sys.exit(1)

    result = "".join(output_lines)
    catalogue.write_text(result, encoding="utf-8")
    print(f"Tagged {tagged_count} entries in {catalogue}")


if __name__ == "__main__":
    main()
