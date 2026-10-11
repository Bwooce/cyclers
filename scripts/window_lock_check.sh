#!/usr/bin/env bash
# #1063: refuse a commit while the lead's catalogue window is open.
#
# The lead opens a window with `scripts/window.sh open "<reason>"`, which writes
# the lock file <git common dir>/CYCLERS_WINDOW_OPEN (inside .git, so never
# tracked), and removes it with `scripts/window.sh close` when the ratchet suite
# is green. While the file exists, every commit is refused unless the committer
# sets CYCLERS_WINDOW_LEAD=1 (the lead's own window commit).
#
# Run from scripts/githooks/pre-commit BEFORE pre-commit's stash (so a refused
# commit does not stash and restore other agents' unstaged files while the
# lead's suite reads them), and again as the first hook in
# .pre-commit-config.yaml for clones that have not set core.hooksPath.
set -u

common_dir=$(git rev-parse --git-common-dir 2>/dev/null) || exit 0
lock="$common_dir/CYCLERS_WINDOW_OPEN"

[ -e "$lock" ] || exit 0
[ "${CYCLERS_WINDOW_LEAD:-}" = "1" ] && exit 0

line=$(head -n 1 "$lock" 2>/dev/null)
who=$(printf '%s' "$line" | cut -d' ' -f1)
when=$(printf '%s' "$line" | cut -d' ' -f2)
why=$(printf '%s' "$line" | cut -d' ' -f3-)
{
    echo "COMMIT REFUSED: a catalogue window is open (#1063)."
    echo "  opened by: ${who:-unknown}"
    echo "  opened at: ${when:-unknown}"
    echo "  reason:    ${why:-not given}"
    echo "Keep your edits in the working tree and commit after the lead reports green"
    echo "(the lead removes $lock with 'scripts/window.sh close')."
    echo "Do not delete the lock file yourself; message the lead (\"main\") if you are blocked."
} >&2
exit 1
