# A tracked runlog that a long run appends to makes every other agent's commit fail its hooks

- Date: 2026-10-07
- Agent: earthmoon-opus
- Seen before: related to 2026-10-05-main-precommit-stashes-other-agents-files.md (pre-commit stashes other agents' unstaged edits). This is a new symptom of the same cause.

What happened: my pathspec commit of `data/997_lineage/` and `scripts/run_997_*.py` failed four times in a row. Each time the message was "mypy ... Failed - files were modified by this hook", followed by "Stashed changes conflicted with hook auto-fixes... Rolling back fixes". mypy itself reported "Success: no issues found". At the same time, `data/968_control/runlog.txt` (tracked, modified, not mine) was being appended by another agent's long run. The whole-tree mypy hook takes about a minute. Any write to a tracked file in that window counts as "modified by this hook", so the commit is rejected.
Workaround: wait for the other run to finish, or retry until a quiet window comes.
Suggested fix: keep runlogs that grow during a run untracked (gitignored) and commit a final copy at the end, or have hooks diff only the committed paths.
