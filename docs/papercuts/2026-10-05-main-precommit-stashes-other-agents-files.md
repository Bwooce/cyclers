# pre-commit stashes and restores other agents' unstaged edits during every commit

- Date: 2026-10-05
- Agent: main
- Seen before: no

What happened: each pathspec commit printed "Stashing unstaged files ... Restored". In a shared checkout that briefly removes other agents' in-progress edits; a concurrent write or a hook crash in that window can lose them.
Workaround: noted in docs/team/coordination.md; agents commit promptly.
Suggested fix: run hooks only on the committed paths without stashing, or use per-agent worktrees.
