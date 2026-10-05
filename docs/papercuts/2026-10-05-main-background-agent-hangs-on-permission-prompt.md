# A background agent hangs for hours on a permission prompt while listed as "running"

- Date: 2026-10-05
- Agent: main
- Seen before: no

What happened: auto mode was off for a few minutes; corpus-review-fable's pending Edit/Bash waited for approval it could not get, and it sat "running" from 17:26 to 19:15 with no progress.
Workaround: detected by transcript mtime, stopped and respawned from its on-disk draft.
Suggested fix: harness should fail a background agent's call instead of waiting; leads check transcript mtime (now in the team-lead skill).
