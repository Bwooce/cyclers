#!/usr/bin/env bash
# #1070: pathspec commit that reports its own hash.
#
#   scripts/commit.sh "<message>" <path> [<path> ...]
#
# Runs `git commit -m "<message>" -- <paths>`. On success it prints exactly one line,
#   COMMITTED <short hash> <subject>
# read with `git log -1` from the commit this call made, so a reported hash comes from
# tool output and never from memory (seven invented hashes on 2026-10-11). On failure it
# prints `NOT COMMITTED (exit N)` and the last lines of the git and hook output, and exits
# with git's code.
#
# The hash is taken from git's own "[branch hash] subject" line, not from HEAD: another
# agent may commit between our commit and the log read.
#
# Untracked paths are added first (only those paths) so the pathspec commit can take them;
# on a refused commit they are removed from the index again, the files stay on disk.
set -u

if [ "$#" -lt 2 ] || [ -z "$1" ]; then
    echo "usage: scripts/commit.sh \"<message>\" <path> [<path> ...]" >&2
    exit 2
fi
msg=$1
shift
for p in "$@"; do
    case "$p" in
        . | ./ | :/ | '*' | -*)
            echo "NOT COMMITTED: pathspec '$p' is not a file or directory path (commit by path)" >&2
            exit 2
            ;;
    esac
done

added=()
for p in "$@"; do
    if [ -e "$p" ] && ! git ls-files --error-unmatch -- "$p" >/dev/null 2>&1; then
        git add -- "$p" || { echo "NOT COMMITTED: git add failed for $p"; exit 1; }
        added+=("$p")
    fi
done

out=$(git commit -m "$msg" -- "$@" 2>&1)
rc=$?
if [ "$rc" -ne 0 ]; then
    if [ "${#added[@]}" -gt 0 ]; then
        git rm -q --cached -r -- "${added[@]}" >/dev/null 2>&1
    fi
    echo "NOT COMMITTED (exit $rc)"
    printf '%s\n' "$out" | tail -n 15
    exit "$rc"
fi

# "[main 1a2b3c4] subject" or "[main (root-commit) 1a2b3c4] subject"
hash=$(printf '%s\n' "$out" | sed -n 's/^\[[^]]* \([0-9a-f]\{7,\}\)\] .*/\1/p' | head -n 1)
if [ -z "$hash" ]; then
    echo "NOT COMMITTED (git exited 0 but printed no commit line; check git log)"
    printf '%s\n' "$out" | tail -n 15
    exit 1
fi
echo "COMMITTED $(git log -1 --format='%h %s' "$hash")"
