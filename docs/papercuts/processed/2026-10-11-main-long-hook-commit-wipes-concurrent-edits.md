# A commit whose pre-commit hooks run for 15 minutes wipes every other agent's unstaged edits made meanwhile

- Date: 2026-10-11
- Agent: main (lead)
- Seen before: yes (coordination.md known hazard: pre-commit stashes other agents' unstaged changes; 2026-10-07 tracked-runlog entries; today's in-window commits)

What happened: the window 3 catalogue commit ran the "catalogue ratchets" hook (pytest on tests/data, including the 600-s Kumar manifold test) on a machine shared with the CI runner, so the hook took over 15 minutes and failed once on a load timeout. For that whole time pre-commit held the other agents' unstaged files in a stash; edits they made while it ran were overwritten when the stash was restored. corpus-henon-opus lost four follow-up edits and ci-keeper-opus lost its README, data/README and OUTSTANDING edits; neither was in ~/.cache/pre-commit/patch*. Thirteen agents then had to be told to stop editing tracked files, by message, mid-turn.
Workaround: both agents rebuilt their edits as reapply scripts with exactly-once anchors and wait for "window 3 landed".
Suggested fix: (1) land #1063's core.hooksPath wrapper so a window refuses other commits mechanically; (2) run the catalogue-window commit with the ratchet hook's tests already proven by the window suite instead of re-running them inside the hook (a hook mode that trusts a fresh suite log, or a lighter hook for catalogue.yaml that checks schema and counts only, with the full ratchets in CI); (3) never open a window while a CI run is in progress on the same machine, which the coordination file already says and which was waived today.

Disposition (2026-10-11 19:58, lead): Backlog #1069 (short hooks on the window commit: the catalogue-ratchet tests run in the window suite and CI, not inside a 15-minute hook that holds everyone's edits in a stash) and Promoted (coordination.md: no edits to tracked files while the lead's window commit is running; the lock reason announces it). Owner decision queued on moving the ratchet run out of the hook.
