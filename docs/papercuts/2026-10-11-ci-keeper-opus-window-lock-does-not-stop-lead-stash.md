# The #1063 window lock does not stop the lead's own commit from stashing everyone's edits

- Date: 2026-10-11
- Agent: ci-keeper-opus
- Seen before: yes (2026-10-11-main-long-hook-commit-wipes-concurrent-edits.md, the same event)

What happened: the README, data/README and OUTSTANDING edits I had ready were wiped between 19:08 and 19:23 AEDT while the window 3 commit's hooks ran. The lock in `scripts/githooks/pre-commit` refuses OTHER agents' commits before their stash. But the lead's commit runs with `CYCLERS_WINDOW_LEAD=1`, goes straight into pre-commit, and stashes the whole tree for as long as its hooks run. So with core.hooksPath set, the wipe in the entry above can still happen.
Workaround: rebuilt the edits as a reapply script with exactly-once anchors; reapplied after "window 3 landed" (ff1eb1e7).
Suggested fix: the main entry's fix (2) is the one that closes this: keep the window commit's hooks short (no full ratchet run inside the hook). Until then, a window rule: no edits to tracked files while the lead's window commit is running, announced by the lock reason.
