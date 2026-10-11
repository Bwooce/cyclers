# #1069: the ratchet hooks trust a fresh window-suite log (pre-registration, 2026-10-12)

Status: PRE-REGISTRATION, written before any code. Owner decision on `#1069` (2026-10-12, verbatim label):
"Hook trusts a fresh suite log". Implementer: ci-keeper-opus. Specification: the lead's dispatch of
2026-10-12 morning.

## 1. Problem

The "catalogue ratchets" and "manifold connection ratchets" pre-commit steps re-run pytest inside the window
commit's hook. On a loaded machine this took over 15 minutes on 2026-10-11. During that time pre-commit holds
every other agent's unstaged edits in a stash and overwrites edits made meanwhile when it restores. Three agents
lost edits (papercuts `2026-10-11-main-long-hook-commit-wipes-concurrent-edits.md` and
`2026-10-11-ci-keeper-opus-window-lock-does-not-stop-lead-stash.md`). The `#1063` lock does not help, because the
lead's own commit is the long one.

## 2. Design

Two new files and two changed hook entries.

`scripts/suite_record.py` (stdlib only):
- `start`: run before the first suite chunk. It computes the CONTENT HASH (sec. 3) of the working tree and writes
  `<git common dir>/CYCLERS_SUITE_PENDING` (hash, UTC start time).
- `finish --log <file>:<exit> ... --covers <pytest target> ...`: run after the last chunk. It reads every log
  (the tee'd full output of one chunk) and its exit status, collects every line that starts with `FAILED ` or
  `ERROR `, checks that each log reached `[100%]`, and recomputes the content hash. It writes
  `<git common dir>/CYCLERS_SUITE_OK` (JSON: start and finish times, hash at start and at finish, each log's path,
  exit status and completion flag, the failed-test list, the covers list, `ok`, and the reasons when `ok` is
  false). It removes the pending file. `ok` is true only when every exit status is 0, the failed list is empty,
  every log reached `[100%]`, and the hash did not change between start and finish. It prints the verdict.
- `check <pytest target> ...`: run by the hook. Exit 0 (trusted) only when ALL hold:
  1. the record exists, parses, and has `ok: true`;
  2. its start time is less than 3600 s ago (age is measured from the START of the suite);
  3. the content hash of the tree being committed (sec. 3, from `git write-tree` on the index git gives the hook)
     equals the record's hash;
  4. every hook target (its file part, before any `::`) is equal to or under one of the record's covers paths.
  It then prints `ratchets: trusted window-suite log <start time, AEDT>` and exits 0. Otherwise it prints one line
  `ratchets: running pytest (<reason>)` and exits 1.

`scripts/ratchet_hook.sh <pytest args>`: runs `check` with the args that are test targets; on exit 0 it exits 0;
otherwise it runs `uv run pytest <args>` as today.

The two hook entries change their `entry` from `uv run pytest ...` to `scripts/ratchet_hook.sh ...` with the same
arguments. The cycler-network hook is not in the owner decision and is unchanged.

## 3. What the hash covers

A SHA-256 over the sorted `(mode, blob id, path)` lines of every TRACKED file in `src/`, `tests/`, `data/`,
`scripts/`, plus `pyproject.toml`, `uv.lock`, `conftest.py` and `.pre-commit-config.yaml`, EXCEPT files ending in
`.md`, `.txt` or `.log` (notes, OUTSTANDING and run logs; no ratchet test reads them; this exclusion is the one
judgement call in the design).
- At `start` and `finish`: the tree is built in a temporary index copied from the real one, with `git add -u` on
  those paths, so it is the working tree's content of tracked files (what pytest read), and it includes new files
  the lead has staged.
- In `check`: the tree is `git write-tree` on the index the hook sees (for `git commit -- <paths>` this is git's
  temporary index, i.e. exactly the commit being made).
Consequence: if ANY agent has an uncommitted edit to a covered file during the suite, or the commit leaves out a
covered file that was edited, the hashes differ and the hook runs pytest as today. That is the safe direction.

## 4. Tests (`tests/scripts/test_1069_suite_record.py`), expected results stated before the code

In a throwaway git repository with the scripts copied in:
- T1 no record -> `check` exits 1 ("no record").
- T2 good logs, `start`/`finish`, nothing changed -> `ok: true`; `check` on a covered target exits 0 and prints
  "trusted window-suite log".
- T3 a covered file edited after `finish` and staged -> `check` exits 1 (hash mismatch).
- T4 a record older than one hour (start time rewritten) -> `check` exits 1 (stale).
- T5 a log with a `FAILED` line, or a non-zero exit, or no `[100%]` -> `ok: false`; `check` exits 1.
- T6 a covered file edited between `start` and `finish` -> `ok: false` (tree changed during the suite).
- T7 a hook target outside the covers -> `check` exits 1.
- T8 an edit to an excluded `.md` file after `finish` -> `check` exits 0 (excluded by design).
- T9 through a real `git commit -- <path>` with a hook that calls `ratchet_hook.sh` with a fake pytest: a trusted
  record skips the fake pytest (no marker), a mismatched one runs it (marker present).

## 5. What the lead does in a window (goes to coordination.md)

```
scripts/suite_record.py start
<chunk 1 ... | tee c1.log>   (note each exit status)
<chunk 2 ... | tee c2.log>
scripts/suite_record.py finish --log c1.log:<exit1> --log c2.log:<exit2> --covers <every pytest target run>
CYCLERS_WINDOW_LEAD=1 git commit -m "..." -- <paths>     (within the hour, nothing else edited)
```
The covers must list the targets exactly as run, without `-k` or `--deselect`. The window suite must include the
hook's own targets (`tests/verify/test_verdict_tier_census.py` and `tests/test_catalogue_rediscovery.py` are not
in `tests/data tests/search tests/scripts`); otherwise the hook runs those targets itself, as today.
