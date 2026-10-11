# A "hold commits, window open" message does not reach a teammate that is mid-turn, so it commits inside the window

- Date: 2026-10-11
- Agent: main (lead)
- Seen before: yes (coordination.md known hazard: teammates read their mailbox only between turns; 2026-10-07 tracked-runlog entries show the same stash race from another direction)

What happened: the catalogue window opened at 13:05 and every agent got "hold all commits until green" at 13:05 and 13:07. jovian-nbody-opus was inside a long turn and committed 8cde2928 at 13:12:16, inside the window. Its pre-commit hook stashed the lead's unstaged patched files while the ratchet suite was reading them; the catalogue chunk had to be re-run to clear the doubt.
Workaround: the lead re-runs the catalogue ratchets after the window's suite and checks the patched files' diff stat and row count after every foreign commit.
Suggested fix: a mechanism, not a message. A lock file (for example `.git/CYCLERS_WINDOW_OPEN` or `data/.window-open`, untracked) that the lead creates when a window opens and removes on green, and a pre-commit hook step that refuses any commit while the file exists unless `CYCLERS_WINDOW_LEAD=1` is set. The refusal text names the lead and the reason. Register as a task.

Disposition (2026-10-11 19:58, lead): Fixed. #1063 lock landed (24e66a48; core.hooksPath = scripts/githooks set 19:46): a commit during a window is refused by the hook. Cause seen four times today before the fix (jovian 13:12; corpus-liang x8, corpus-henon, corpus-gobetz in window 3).
