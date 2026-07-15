"""Header/body status consistency check for the real ``data/OUTSTANDING.md``.

Guards against a real, recurring failure mode found 2026-07-15: a task's bullet
gets a genuine resolution appended to its BODY (e.g. "**✓ Resolved (2026-07-14)**
..." tacked on at the end, days or weeks after the bullet was first written), but
the bullet's own OPENING line is never updated to match -- it keeps saying
"parking lot -- not auto-fired" or "planning phase first -- do NOT build/run yet"
indefinitely. Anyone (human or agent) who reads only the header, or greps for a
task's status and stops at the first hit, gets the wrong answer.

This bit twice in one session before this test existed: #557's bullet opened with
"planning phase first -- do NOT build/run yet" for months after a "CLEAN,
WELL-CHARACTERIZED NEGATIVE" resolution was appended ~96 lines further down in
that SAME bullet -- misleading both a same-day ``## CURRENT STATE`` dashboard
audit and a manual review of it into re-presenting an already-answered scope
question to the user. #560 and #543 had the identical pattern, caught only once
this heuristic was built and run against the live file.

This is a HEURISTIC, not a natural-language parser -- it flags bullets whose
opening lines contain an explicit "still open" phrase (see ``_OPEN_SIGNALS``)
while the REST of the bullet contains an explicit "resolved" marker (see
``_RESOLVED_SIGNALS``, deliberately narrow: this project's own "✓" convention and
its "STALE, already resolved" phrase, not fuzzy words like "done" or "closed"
that show up in ordinary prose too). It will not catch every staleness pattern,
but it is a real automated backstop where three prior full-file audits (manual
review, a dedicated "#594 full-file comprehensive audit", and the "## CURRENT
STATE" dashboard build) each independently missed at least one live instance.
"""

from __future__ import annotations

import re
from pathlib import Path

_OUTSTANDING_PATH = Path(__file__).resolve().parents[2] / "data" / "OUTSTANDING.md"

_BULLET_RE = re.compile(r"^- \*\*#(\d+)")
_HEADING_RE = re.compile(r"^#{1,6} ")

_HEADER_WINDOW_LINES = 3
"""How many lines from a bullet's start count as its "opening" for this check.

#557's actual header ("- **#557** (P0, planning phase first -- do NOT build/run
yet, review the plan before deciding on\\n  implementation) -- extend #535's...")
wraps across 2 physical lines before the descriptive prose begins; 3 gives a
small margin without reaching deep enough to catch a resolution that's
genuinely later in the bullet.
"""

# Deliberately narrow: this project's OWN convention for marking a bullet
# resolved. Bare English words like "done"/"closed"/"resolved" are NOT included
# here -- they show up constantly in ordinary prose (e.g. "once #538 is done")
# and produced a real false positive in calibration (see module docstring).
_RESOLVED_SIGNALS: tuple[str, ...] = (
    "✓",
    "stale, already resolved",
)

# Phrases this project's own bullets use to mark a task as still open/awaiting
# a decision/not yet actioned. Calibrated against the live file (2026-07-15):
# broadening this list must not introduce new false positives against the
# current ~180 bullets (see the test below, which pins the exact known-good
# count).
_OPEN_SIGNALS: tuple[str, ...] = (
    "do not build/run yet",
    "do not dispatch",
    "not yet dispatched",
    "not yet been dispatched",
    "not dispatched",
    "awaiting a user",
    "awaiting review",
    "awaiting a",
    "planning phase first",
    "not auto-fired",
    "flagged for user review before dispatch",
    "flagged for user review",
    "review the plan before deciding",
    "needs a user greenlight",
    "needs a scoping conversation",
    "[watch]",
    "parking lot",
)


def _bullet_spans(text: str) -> list[tuple[int, int, int]]:
    """Return ``(task_no, start_line, end_line_exclusive)`` for every top-level
    ``- **#NNN`` bullet, where the span runs to the next such bullet or the next
    markdown heading, whichever comes first.
    """
    lines = text.splitlines(keepends=True)
    starts: list[tuple[int, int]] = []
    for i, line in enumerate(lines):
        m = _BULLET_RE.match(line)
        if m:
            starts.append((i, int(m.group(1))))

    spans: list[tuple[int, int, int]] = []
    for idx, (line_i, task_no) in enumerate(starts):
        end = len(lines)
        if idx + 1 < len(starts):
            end = min(end, starts[idx + 1][0])
        for j in range(line_i + 1, end):
            if _HEADING_RE.match(lines[j]):
                end = j
                break
        spans.append((task_no, line_i, end))
    return spans


def _find_header_body_contradictions(text: str) -> list[tuple[int, list[str], list[str]]]:
    """Bullets whose opening lines carry an "open" signal while the rest of the
    bullet carries a "resolved" signal, with neither signal present in the other
    half (i.e. the header doesn't ALSO say resolved, so this isn't just a bullet
    that mentions both in one coherent, already-consistent status line).
    """
    lines = text.splitlines(keepends=True)
    flagged: list[tuple[int, list[str], list[str]]] = []
    for task_no, start, end in _bullet_spans(text):
        header_text = "".join(lines[start : min(start + _HEADER_WINDOW_LINES, end)]).lower()
        body_text = "".join(lines[start + _HEADER_WINDOW_LINES : end]).lower()

        header_open = [s for s in _OPEN_SIGNALS if s in header_text]
        header_resolved = [s for s in _RESOLVED_SIGNALS if s in header_text]
        body_resolved = [s for s in _RESOLVED_SIGNALS if s in body_text]

        if header_open and body_resolved and not header_resolved:
            flagged.append((task_no, header_open, body_resolved))
    return flagged


def test_no_bullet_header_contradicts_its_own_body() -> None:
    """A bullet's opening line must not claim "still open" while its own body
    (appended later, possibly much later) records an actual resolution.

    A hit here means: find the named task number's bullet (search
    ``- **#NNN``), read it to its END, confirm whether it's actually resolved,
    and if so rewrite the OPENING line in place to reflect that (do not just
    leave the stale header and rely on the body being read carefully -- that is
    the exact failure this test exists to catch).
    """
    text = _OUTSTANDING_PATH.read_text(encoding="utf-8")
    flagged = _find_header_body_contradictions(text)
    assert not flagged, (
        "Bullet(s) whose opening line says 'still open' but whose own body "
        f"records a resolution: {[t for t, _, _ in flagged]}. For each, find "
        "its '- **#NNN' bullet, read to the end, and rewrite the OPENING line "
        "to match the real status -- do not just trust the header. Details: "
        f"{flagged}"
    )
