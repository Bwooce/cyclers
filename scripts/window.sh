#!/usr/bin/env bash
# #1063: the lead's catalogue-window lock.
#
#   scripts/window.sh open "<reason>"   create the lock (refuses if one is open)
#   scripts/window.sh close             remove the lock (after the suite is green)
#   scripts/window.sh status            show the lock, if any
#
# The lock is <git common dir>/CYCLERS_WINDOW_OPEN, one line:
#   lead <ISO time> <reason>
# While it exists, scripts/window_lock_check.sh refuses every commit unless
# CYCLERS_WINDOW_LEAD=1 is set, so the lead commits the window with
#   CYCLERS_WINDOW_LEAD=1 git commit -m "..." -- <paths>
set -eu

common_dir=$(git rev-parse --git-common-dir)
lock="$common_dir/CYCLERS_WINDOW_OPEN"

case "${1:-}" in
    open)
        reason="${2:-}"
        if [ -z "$reason" ]; then
            echo "usage: scripts/window.sh open \"<reason>\"" >&2
            exit 2
        fi
        if [ -e "$lock" ]; then
            echo "a window is already open: $(head -n 1 "$lock")" >&2
            exit 1
        fi
        printf 'lead %s %s\n' "$(date -Iseconds)" "$reason" > "$lock"
        echo "window OPEN: $(head -n 1 "$lock")"
        ;;
    close)
        if [ ! -e "$lock" ]; then
            echo "no window is open" >&2
            exit 1
        fi
        echo "closing: $(head -n 1 "$lock")"
        rm -f "$lock"
        echo "window CLOSED"
        ;;
    status)
        if [ -e "$lock" ]; then
            echo "OPEN: $(head -n 1 "$lock")"
        else
            echo "no window open"
        fi
        ;;
    *)
        echo "usage: scripts/window.sh open \"<reason>\" | close | status" >&2
        exit 2
        ;;
esac
