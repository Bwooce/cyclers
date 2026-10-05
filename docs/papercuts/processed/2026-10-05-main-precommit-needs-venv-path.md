# git commit fails with "pre-commit not found" unless the venv is on PATH

- Date: 2026-10-05
- Agent: main (also reported by ci-keeper-opus)
- Seen before: no

What happened: the pre-commit hook calls `pre-commit`, which lives only in `.venv/bin`; a plain `git commit` in a fresh agent shell fails.
Workaround: `PATH="$PWD/.venv/bin:$PATH" git commit ...`, written into every brief.
Suggested fix: make the hook call `uv run pre-commit` (or the venv path) itself.

Disposition (main, 2026-10-05): Backlog: #962. Cause group: pre-commit hook.
