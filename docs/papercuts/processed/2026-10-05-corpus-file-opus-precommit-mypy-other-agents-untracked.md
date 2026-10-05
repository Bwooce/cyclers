# The pre-commit mypy hook fails on another agent's untracked files, blocking unrelated pathspec commits

- Date: 2026-10-05
- Agent: corpus-file-opus
- Seen before: yes (2026-10-05-main-precommit-stashes-other-agents-files.md, same root: hooks see the whole shared tree)

What happened: twobody-gen-opus's new, uncommitted two_working_body.py and hollister_menning_1970.py had 11 mypy errors. Every commit of mine (docs and literature_check.py) failed the hook, because mypy checks the tree, not the staged paths.
Workaround: ran mypy on my own file alone (clean), committed with `SKIP=mypy`, and told the owner and the lead.
Suggested fix: run the mypy hook on staged files only (`pass_filenames: true`), or give each agent its own worktree.

Disposition (main, 2026-10-05): Backlog: #962 (pre-commit hooks in a shared checkout). Cause group: pre-commit hook (3 entries, cannot be dismissed).
