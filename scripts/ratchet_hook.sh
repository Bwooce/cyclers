#!/usr/bin/env bash
# #1069: pre-commit entry for the catalogue and manifold-connection ratchet hooks.
#   scripts/ratchet_hook.sh <pytest args>
# If the lead's window suite left a record (scripts/suite_record.py) that is ok, under an hour old,
# taken on exactly the tree being committed and covering these targets, skip the pytest re-run and
# say so; otherwise run `uv run pytest <args>` as before. Owner decision 2026-10-12, design in
# docs/notes/2026-10-12-1069-hook-trusts-fresh-suite-log.md.
set -u

repo_root=$(git rev-parse --show-toplevel)
if python3 "$repo_root/scripts/suite_record.py" check "$@"; then
    exit 0
fi
exec uv run pytest "$@"
