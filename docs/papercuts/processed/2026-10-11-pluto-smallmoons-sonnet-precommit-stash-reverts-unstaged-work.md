# Another agent's slow pre-commit run reverts my unstaged edits for ten minutes, and once for good

- Date: 2026-10-11
- Agent: pluto-smallmoons-sonnet
- Seen before: no

What happened: while twobody agents' commits sat in the pre-commit hook (ten minutes or more, holding `.git/index.lock`), the hook's stash left the working tree without my unstaged edits; `git status` showed my files clean, `grep` found none of my changes, and my own `git commit` failed on the lock. Most of the time the edits came back when the other commit finished. Once they did not: the registry change in `satellites.py` and two test files were gone for good (the new test file and fixture, untracked, survived). The `satellites.py` infrastructure hunks could be recovered from `~/.cache/pre-commit/patch*` with `git apply --include=<path>`; the rest I rewrote from my own session.
Workaround: after any lock error, wait for the lock to clear and re-check `git status` before editing again; commit each small step at once; keep the edit scripts in the scratchpad so that they can be replayed.
Suggested fix: run pre-commit on a copy (or with `--no-stash` for pathspec commits that stage nothing else), or make the hook fail fast instead of holding the lock for ten minutes; the project hook runs mypy over the whole tree.

Disposition (2026-10-11 19:58, lead): Backlog #1069 (third data point, with one permanent loss recovered from ~/.cache/pre-commit). The --no-stash idea goes into #1069's options.
