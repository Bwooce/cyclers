"""Structural integrity checks for the real ``data/OUTSTANDING.md``.

Guards against a real failure mode found 2026-07-10: a concurrent agent's
string-replace edit clipped the ``- **#549** (...)`` bullet header off its
own paragraph while appending unrelated content nearby, leaving an orphaned
fragment ("genome sweep. The Tier-1 forward-plan line ...") starting
mid-sentence with no header at all. The TASK ALLOCATIONS ledger line still
claimed #549 was allocated, but no bullet in the document body actually
opened with that number any more -- silently desyncing the ledger from the
document body. Caught by hand while adding a follow-on task; this pins it
as a regression test.
"""

from __future__ import annotations

import re
from pathlib import Path

_OUTSTANDING_PATH = Path(__file__).resolve().parents[2] / "data" / "OUTSTANDING.md"
_LEDGER_ENTRY_RE = re.compile(r"#(\d+) for ")

# Lenient bullet-opening-line check: does the bullet's FIRST line mention a
# task number at all, in EITHER the strict "- **#NNN**" preflight-gate style
# or the looser historical "- ✓ Resolved (date) **#NNN** — ..." style. This
# intentionally does not require the strict immediate-bold-close the
# preflight gate's own regex demands (several older, already-resolved
# entries use the looser style) -- it only needs to confirm SOME bullet
# marker + number exists on the line, ruling out a fully orphaned paragraph
# (no leading "- " at all) that lost its header to a bad edit.
_BULLET_OPENING_WITH_NUMBER_RE = re.compile(r"^-\s.*?#(\d+)", re.MULTILINE)


def _ledger_declared_task_numbers() -> set[int]:
    """Task numbers the TASK ALLOCATIONS ledger line(s) claim to have allocated."""
    text = _OUTSTANDING_PATH.read_text(encoding="utf-8")
    ledger_lines = [line for line in text.splitlines() if "TASK ALLOCATIONS" in line]
    assert ledger_lines, "no TASK ALLOCATIONS ledger line found in data/OUTSTANDING.md"
    declared: set[int] = set()
    for line in ledger_lines:
        declared.update(int(m) for m in _LEDGER_ENTRY_RE.findall(line))
    return declared


def _task_numbers_with_intact_bullet_openers() -> set[int]:
    """Task numbers whose FIRST bullet-opening line mentions that number,
    in either the strict or the looser historical style (see module note).
    """
    text = _OUTSTANDING_PATH.read_text(encoding="utf-8")
    return {int(m) for m in _BULLET_OPENING_WITH_NUMBER_RE.findall(text)}


def test_every_ledger_declared_task_has_an_intact_bullet() -> None:
    """Every "#NNN for ..." the ledger declares must resolve to a bullet
    whose OWN opening line mentions that number -- not just be referenced in
    the ledger line itself or somewhere deep in another task's paragraph. A
    mismatch means a bullet header was clipped or never written, leaving an
    orphaned paragraph with no leading "- " marker at all (the exact #549
    failure mode this test pins).
    """
    declared = _ledger_declared_task_numbers()
    actual = _task_numbers_with_intact_bullet_openers()
    missing = sorted(declared - actual)
    assert not missing, (
        f"TASK ALLOCATIONS ledger declares {missing} but no bullet's opening "
        f"line mentions that number -- a bullet header was likely clipped or "
        f"never written. Search for the task number's orphaned paragraph "
        f"text (it will start mid-sentence with no leading '- ' marker) and "
        f"restore its bullet header line."
    )
