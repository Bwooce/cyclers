# ruff: noqa: E501  (catalogue YAML and evidence text carry long lines)
# STATUS (2026-10-12): PARKED (owner hold on gc4-19). Dry-run clean on 3ccbf5e7 (with the ev script; ruff, mypy, tests/data, writeback/census/rediscovery exit 0) using window4-opus's crank-sibling wording. NOT yet updated to jovian-1055-opus's exact sec. 14.1 replacement text and orbit_elements.note wording (sent 2026-10-12); swap those in and re-dry-run before use.
"""Catalogue window 4: add the gc4-19 row (#1055; owner decision 2026-10-11, verbatim label
"Candidate-novel V1 row, like the 15"; docs/notes/2026-10-11-1055-degenerate-family-seed-correction.md
sec. 14).

VARIANT (crank-sibling wording, #1055 sec. 14.1, 977bf215; lead order 2026-10-12): the same patch
as apply_gc4_19.py, but the row states that the continuous-gravity limit is a CRANK-SIBLING of the
#973 member (2:1 Ganymede leg V_inf direction 6.04 deg off; epochs, V_inf, turns and periapses within
amendment 10; amendment 14 fails), and orbit_elements.inclination_deg is 7.08 (the chain tool's 2:1
Ganymede resonant leg) instead of 0.0. Apply this OR apply_gc4_19.py, never both (each refuses if
gc4-19 is already present).

Base script: rebuilt 2026-10-12 from the committed record after the session scratchpad was wiped. It reproduces
jovian-nbody-opus's window-4 script (md5 1cf48acb) and its row file (md5 3661644c, now generated
here from the committed draft and checked against that md5) with these changes: the repo root is a
required argument, every check runs before any file is written, a re-run refuses with a message, and
the edited Python files are run through ``ruff format``.

Edits (each anchor must match exactly once; census counts are read and incremented, not typed):
- data/catalogue.yaml: the gc4-19 row, between gc4-10 and gc4-25. The row is
  data/1055_gc4-19/gc4-19_row_draft.yaml with the sec. 14 changes: the owner decision replaces the
  DRAFT comment; the data_gaps "notes" entry names the method controls (A, B, the amendment-8
  negative control) and the amendment-10 member-identity pass.
- src/cyclerfinder/data/validate.py: the (gc4-19, V1) evidence entry (_V1_1025 plus its own line).
- tests/data/test_cycler_class_census.py: the multi-arc allowlist id, its length assert, the
  census multi-arc count.
- tests/data/test_schema_v45_fields.py: the gc4-19 -> V1 map entry (the V2 count is not touched).
- tests/data/test_validation_tier_census.py: "unvalidated" + 1.
- tests/test_catalogue_rediscovery.py: NON_HELIOCENTRIC + 1 (primary Jupiter).

Usage: uv run python scripts/window4/apply_gc4_19_v2.py <repo-root>
"""

from __future__ import annotations

import hashlib
import re
import subprocess
import sys
from pathlib import Path

RID = "ganymede-callisto-two-working-body-cycler-gc4-19-2026"
NEXT = "ganymede-callisto-two-working-body-cycler-gc4-25-2026"
DRAFT = "data/1055_gc4-19/gc4-19_row_draft.yaml"
ROW_MD5 = "7f5abd0ef8487c16db6c124d5e35a5f6"  # the generated crank-sibling row (sec. 14.1 variant)

# The sec. 14 changes to the committed draft (the first three as jovian-nbody-opus applied them,
# the identity clause reworded and the last four added for the sec. 14.1 crank-sibling correction).
ROW_CHANGES: list[tuple[str, str]] = [
    (
        "  our_status: candidate-novel   # DRAFT, OWNER DECISION PENDING (held 2026-10-11 until "
        "GanCal#1 on jup365 reports): proposed as for the 15 #1025 members (a member of gc-1's "
        "class, gc-1's attribution).",
        '  our_status: candidate-novel   # OWNER DECISION 2026-10-11 (#1055): "Candidate-novel V1 '
        "row, like the 15\". A member of gc-1's two-working-body Ganymede-Callisto class, written "
        "candidate-novel with gc-1's attribution; reached in continuous gravity from the #1055 "
        "corrected seed. docs/notes/2026-10-11-1055-degenerate-family-seed-correction.md sec. 14.",
    ),
    (
        "The branch's sigma -> 0 limit is this #973 member (Richardson, within 14 s; #973 gate "
        're-check, fd5ee948)."',
        "The branch's sigma -> 0 limit is a CRANK-SIBLING of this #973 member, in gc4-19's degenerate "
        "patched-conic family (sec. 14.1, correcting sec. 10.1): the same epochs (Richardson, within "
        "14 s; #973 gate re-check, fd5ee948, which corrects dates only and does not carry the crank), "
        "V_inf magnitudes (within 1e-4), turns (within 0.1 deg) and periapses (r_p / sigma within "
        "0.4 %), but on the 2:1 Ganymede resonant leg the V_inf direction differs by 6.04 deg "
        "(crank about 81 deg against the chain tool's 90 deg; pump the same, 41.76 deg). Method "
        "controls: A, B and the amendment-8 negative control satisfied (secs. 8, 12.3). Member "
        "identity: amendment 10 (sec. 10.1) PASS; amendment 14 (V_inf direction < 0.1 deg at every "
        "node, sec. 13.5) FAILS on the 2:1 leg (6.04 deg). No other gc member is a crank-sibling "
        '(sec. 14.1 screen)."',
    ),
    (
        '    source_hint: "docs/notes/2026-10-11-1055-degenerate-family-seed-correction.md secs. 8, '
        '10, 10.1; docs/notes/2026-10-11-1055-gc4-19-973-gate.md"',
        '    source_hint: "docs/notes/2026-10-11-1055-degenerate-family-seed-correction.md secs. 8, '
        "10, 10.1, 12.3, 13.5, 14, 14.1; docs/notes/2026-10-11-1055-gc4-19-973-gate.md; "
        'data/1055_vinf_direction.json"',
    ),
    (
        "    masses (identity by the limit, IAS15 re-fly 1.5e-05 km, no unscheduled minimum",
        "    masses (identity by the limit in epochs, V_inf, turns and periapses; on the 2:1 Ganymede\n"
        "    leg the limit is a crank-sibling, V_inf direction 6.04 deg off, sec. 14.1; IAS15 re-fly\n"
        "    1.5e-05 km, no unscheduled minimum",
    ),
    (
        "    inclination_deg: 0.0\n",
        "    inclination_deg: 7.08   # the chain tool's 2:1 Ganymede resonant leg (V_inf 2.3565 km/s, "
        "pump 41.76 deg, crank 90 deg; #1055 sec. 14.1). The continuous-gravity crank-sibling's leg "
        "is 6.99 deg. The Lambert legs and the Callisto 1:1 leg lie in the moons' plane.\n",
    ),
    (
        'note: "Multi-arc (no single ellipse). Distance from Jupiter',
        'note: "Multi-arc (no single ellipse); not coplanar: the 2:1 Ganymede resonant leg is inclined '
        "7.08 deg to the moons' plane (#1055 sec. 14.1). Distance from Jupiter",
    ),
    (
        '"data/1055_gc4-19/v1_evidence.json", "data/1025_litcheck_973_members.json"]',
        '"data/1055_gc4-19/v1_evidence.json", "data/1055_vinf_direction.json", '
        '"data/1025_litcheck_973_members.json"]',
    ),
]

PY_FILES = [
    "src/cyclerfinder/data/validate.py",
    "tests/data/test_cycler_class_census.py",
    "tests/data/test_schema_v45_fields.py",
    "tests/data/test_validation_tier_census.py",
    "tests/test_catalogue_rediscovery.py",
]


class Patch:
    """Edits held in memory; nothing is written until every edit has been prepared."""

    def __init__(self, root: Path) -> None:
        self.root = root
        self.files: dict[str, str] = {}

    def text(self, rel: str) -> str:
        if rel not in self.files:
            self.files[rel] = (self.root / rel).read_text()
        return self.files[rel]

    def edit(self, rel: str, old: str, new: str) -> None:
        s = self.text(rel)
        n = s.count(old)
        if n != 1:
            sys.exit(
                f"ABORT (nothing written): {rel}: anchor found {n} times (need 1):\n{old[:300]}"
            )
        self.files[rel] = s.replace(old, new)

    def bump(self, rel: str, pattern: str, comment: str) -> tuple[int, int]:
        """Increment the one integer matched by ``pattern`` (one group); put ``comment``
        (formatted with old, new) on the line above, at the same indentation."""
        s = self.text(rel)
        ms = list(re.finditer(pattern, s, flags=re.M))
        if len(ms) != 1:
            sys.exit(f"ABORT (nothing written): {rel}: {pattern!r} matched {len(ms)} times")
        m = ms[0]
        old = int(m.group(1))
        line_start = s.rfind("\n", 0, m.start()) + 1
        indent = re.match(r"[ ]*", s[line_start:]).group(0)  # type: ignore[union-attr]
        line = s[line_start : m.end()]
        new_line = line[: m.start(1) - line_start] + str(old + 1) + line[m.end(1) - line_start :]
        self.files[rel] = (
            s[:line_start]
            + indent
            + comment.format(old=old, new=old + 1)
            + "\n"
            + new_line
            + s[m.end() :]
        )
        return old, old + 1

    def write(self) -> None:
        for rel, s in self.files.items():
            p = self.root / rel
            if p.read_text() != s:
                p.write_text(s)
                print(f"edited {rel}")


def build_row(root: Path) -> str:
    row = (root / DRAFT).read_text()
    for old, new in ROW_CHANGES:
        n = row.count(old)
        if n != 1:
            sys.exit(f"ABORT (nothing written): {DRAFT}: anchor found {n} times:\n{old[:200]}")
        row = row.replace(old, new)
    md5 = hashlib.md5(row.encode()).hexdigest()
    if md5 != ROW_MD5:
        sys.exit(
            f"ABORT (nothing written): the generated row has md5 {md5}, not {ROW_MD5} "
            f"(the committed draft {DRAFT} has changed since window 4 was prepared)"
        )
    if not (row.startswith("\n- id: " + RID + "\n") and row.endswith("\n")):
        sys.exit("ABORT (nothing written): the generated row does not start with its id line")
    return row


def main() -> None:
    if len(sys.argv) != 2:
        sys.exit("usage: uv run python scripts/window4/apply_gc4_19_v2.py <repo-root>")
    root = Path(sys.argv[1]).resolve()
    if not (root / "data" / "catalogue.yaml").is_file():
        sys.exit(f"ABORT: {root} is not a cyclers checkout (no data/catalogue.yaml)")

    # Refuse a second application before preparing anything.
    for rel in ("data/catalogue.yaml", *PY_FILES[:3]):
        if RID in (root / rel).read_text():
            sys.exit(
                f"REFUSED: the gc4-19 patch is already applied ({RID} is in {rel}); nothing written"
            )

    patch = Patch(root)
    row = build_row(root)

    # 1. the catalogue row, between gc4-10 and gc4-25
    patch.edit("data/catalogue.yaml", f"\n\n- id: {NEXT}\n", "\n" + row + f"\n- id: {NEXT}\n")

    # 2. validate.py V1 evidence entry, before gc4-25's
    patch.edit(
        "src/cyclerfinder/data/validate.py",
        f'    ("{NEXT}", "V1"): (\n',
        f'    ("{RID}", "V1"): (\n'
        '        _V1_1025 + "Row: gc4-19 (lamberthub <= 6.3e-12 m/s; DOP853 re-fly 1.5e-04 km). "\n'
        '        "Continuous gravity reached from the #1055 degenerate-family-corrected seed "\n'
        '        "(data/1055_gc4-19/; owner decision 2026-10-11); its limit is a crank-sibling of "\n'
        '        "the #973 member (2:1 Ganymede leg V_inf direction 6.04 deg off; sec. 14.1)."\n'
        "    ),\n"
        f'    ("{NEXT}", "V1"): (\n',
    )

    # 3. multi-arc allowlist: the id, the length assert and its message, the census dict
    census = "tests/data/test_cycler_class_census.py"
    patch.edit(
        census,
        f'        "{NEXT}",\n',
        "        # 2026-10-11 (#1055, owner decision): gc4-19. multi-arc + 1.\n"
        f'        "{RID}",\n'
        f'        "{NEXT}",\n',
    )
    old, new = patch.bump(
        census,
        r"^assert len\(MULTI_ARC_ALLOWLIST\) == (\d+), \($",
        "# 2026-10-12 (window 4; #1055 owner decision 2026-10-11): + gc4-19: {old} -> {new}.",
    )
    patch.edit(census, f"Allowlist must have {old} entries", f"Allowlist must have {new} entries")
    patch.bump(
        census,
        r'^    expected = \{"single-ellipse": \d+, "multi-arc": (\d+), "non-keplerian": \d+\}$',
        "# 2026-10-12 (window 4; #1055 owner decision 2026-10-11): + gc4-19: multi-arc {old} -> "
        "{new}.",
    )

    # 4. validation_level map (V1) and the unvalidated source-pair tier
    patch.edit(
        "tests/data/test_schema_v45_fields.py",
        f'        "{NEXT}": "V1",\n',
        "        # 2026-10-11 (#1055, owner decision): gc4-19, V1 (validate._V1_1025).\n"
        f'        "{RID}": "V1",\n'
        f'        "{NEXT}": "V1",\n',
    )
    patch.bump(
        "tests/data/test_validation_tier_census.py",
        r'^    "unvalidated": (\d+),$',
        "# {old} -> {new} (2026-10-12, window 4; #1055): + gc4-19, same convention as the 15 "
        "#1025 gc rows.",
    )

    # 5. rediscovery exclusion census: Jovian primary
    patch.bump(
        "tests/test_catalogue_rediscovery.py",
        r"^    ExclusionReason\.NON_HELIOCENTRIC: (\d+),$",
        "# {old} -> {new} (2026-10-12, window 4; #1055): + gc4-19 (primary Jupiter).",
    )

    patch.write()
    subprocess.run(["uv", "run", "ruff", "format", *PY_FILES], check=True, cwd=root)
    print(
        "OK: gc4-19 (crank-sibling variant) applied (catalogue rows + 1, candidate-novel + 1, V1 + 1)"
    )


if __name__ == "__main__":
    main()
