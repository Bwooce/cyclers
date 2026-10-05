# Team coordination — cyclers

Current state only; edit in place. Live roster = `ListAgents`; this file is a snapshot.

## Push model
Teammates commit locally with pathspec commits and never push. The lead pushes, after reading each
commit's file list against its author's ownership (`~/.claude/skills/team-lead/safe-push.sh`).

## Git rules (copy into every brief)
1. `git commit -m "msg" -- path1 path2` only. Never `git add -A`, `commit -a`, bare commit after `git add`.
2. A pathspec commit takes the whole file: `git diff <file>` first if others edit it.
3. Never `commit --amend`, `rebase`, `reset`, `stash`, `checkout -- <file>`, `checkout .`. Fix forward.
4. Commit each item as it finishes. No Co-Authored-By / Claude-Session / AI attribution lines.
5. Never create branches. Never delete files you did not create this session.
6. Commit messages start with the task number (`#NNN: ...`).

## Shared files (diff before commit)
- `data/OUTSTANDING.md` (CURRENT STATE block is edited in place)
- `data/catalogue.yaml` (any change: run ALL ratchets, `uv run pytest tests/data tests/search -q`)

## Long-lived teammates
- see `docs/team/lead-log.md` for current dispatches

## Shared resources
- Machine CPU: only one full-suite pytest run at a time (8-way parallel runs collide).

## Known hazards
- The pre-commit hook stashes other agents' unstaged changes during each commit and restores them afterwards. Avoid committing while another agent is in the middle of editing a file that has the same hook scope. If a teammate's edit disappears, look in `~/.cache/pre-commit/patch*`.
- The `git commit` hook needs the venv: `PATH="$PWD/.venv/bin:$PATH" git commit ...`.
