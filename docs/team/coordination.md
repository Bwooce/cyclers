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
7. Report a commit hash only from `git log -1` AFTER the commit returns; never from memory (four wrong hashes on 2026-10-07).
8. If a commit is refused with "COMMIT REFUSED: a catalogue window is open (#1063)", keep your edits in the working tree, do other work, and commit after the lead reports green. Never delete the lock file or set `CYCLERS_WINDOW_LEAD`.

## Shared files (diff before commit)
- `data/OUTSTANDING.md` (CURRENT STATE block is edited in place)
- `data/catalogue.yaml` (any change: run ALL ratchets, `uv run pytest tests/data tests/search -q`)

## Long-lived teammates
- see `docs/team/lead-log.md` for current dispatches

## Shared resources
- Machine CPU: only ONE full-suite pytest run at a time on Shed-Air, counting the cyclers CI runner (it runs on this Mac). Do a local full run only while no cyclers CI run is in progress, and do not push while a local full run is going. The offgridgate project also has a CI runner here (next-server, playwright), and its runs add load that we cannot schedule. Any pytest-timeout seen under heavy load is suspect until it has been re-run alone.

## Catalogue windows (#1063 lock; replaces the "hold commits" message)
A "hold commits" message does not reach a teammate in the middle of a long turn, so on 2026-10-11 four agents committed inside window 3. The window is now a lock that git enforces:
- The lead opens a window with `scripts/window.sh open "<reason>"`. This writes `.git/CYCLERS_WINDOW_OPEN` (one line, `lead <ISO time> <reason>`; inside `.git`, so never tracked).
- While the lock exists, every commit is refused with a message naming the lead, the open time and the reason. The lead commits the window with `CYCLERS_WINDOW_LEAD=1 git commit -m "..." -- <paths>`.
- After the ratchet suite is green, the lead runs `scripts/window.sh close`, which removes the lock, and reports green. `scripts/window.sh status` shows the state.
- The check runs in `scripts/githooks/pre-commit` BEFORE pre-commit stashes the tree, so a refused commit does not stash and restore files the lead's suite is reading. The `window-lock` entry at the top of `.pre-commit-config.yaml` is a backstop for clones without the setup below; it refuses the commit, but only after pre-commit has stashed.
- Setup for a fresh clone (once): `uv run pre-commit install && git config core.hooksPath scripts/githooks`. Check it with `git config core.hooksPath` (prints `scripts/githooks`).
- Reinstall the pre-commit hook (new venv, new pre-commit version): with `core.hooksPath` set, `pre-commit install` refuses to run, so do `git config --unset core.hooksPath && uv run pre-commit install && git config core.hooksPath scripts/githooks`.
- Test: `tests/scripts/test_1063_window_lock.py`.

## Known hazards
- The pre-commit hook stashes other agents' unstaged changes during each commit and restores them afterwards. Avoid committing while another agent is in the middle of editing a file that has the same hook scope. If a teammate's edit disappears, look in `~/.cache/pre-commit/patch*`.
- The `git commit` hook does not need a PATH prefix. The installed hook runs pre-commit with this clone's own `.venv/bin/python` (its `INSTALL_PYTHON` line), and the local hooks call `uv run` (#962). After moving or re-creating the venv, reinstall the hook (the reinstall steps are under "Catalogue windows").
- Background agents hang forever on a permission prompt (seen 2026-10-05: auto mode was switched off for a few minutes and the Fable reviewer froze for 2 h, still listed as "running"). Liveness check = the agent's transcript mtime under `~/.claude/projects/<proj>/<session>/subagents/`, not the agent list.
- Teammate messages addressed to "team-lead" go to the file mailbox `~/.claude/teams/session-<id>/inboxes/team-lead.json`, which the in-process lead conversation never reads; "main" is queued straight into the lead conversation. Briefs say "message main". Verified 2026-10-05 from delivery receipts. Original note: (2026-10-05: refcheck-sonnet, corpus-review-fable-2, twobody-gen-opus and corpus-file-opus reports all lost). Briefs must tell teammates to message the lead session by its ListAgents name ("main" works; "team-lead" and "cyclers" do not). A lost report can be recovered from the agent's transcript (`subagents/agent-a<name>-*.jsonl`, SendMessage tool_use input).
- Teammates read their own mailbox only between turns: a teammate deep in one long turn does not see messages (the first #938 Fable reviewer has 5 unread lead messages in its inbox file). Put must-know items in the brief; ask for a receipt within the first minute.

## Launch and numerics rules (from the papercut review 2026-10-07)
- Validate a long-compute launch with the EXACT launch command and flags (same script, --n-cycles, epochs, cells), cutting only wall time (fewer restarts, a --sample); a shorter chain or one cell does not cover per-cell branches (two launches lost on 2026-10-06).
- Every new phase of an instrumented script gets its own progress lines before its first pilot.
- Before passing a kwarg the code used to leave at default, check the library's CURRENT default (`inspect.signature`); SciPy 1.16 changed lm's x_scale and a stall was misread as a no-root result.
- An independent checker takes every physical constant (GM, states, periods) from the reference model, never from the object under test, and prints them.
- A running process never appends to a tracked file. Live logs go under `data/*/live/` (gitignored) or `*.live.log` in the run root; a final static copy is committed when the run ends. A growing tracked file makes every other agent's pre-commit hook fail (papercuts 2026-10-05, 2026-10-07 x2).
- Pilot a long run at the load it will see (or scale the pilot rate by the current load average), and size shards for the 10-minute tool limit at the loaded rate. Every enumeration cell gets a coverage validator (every structure once, no duplicate zero lines) like `run_973_enumerate.py check` (papercut 2026-10-07).
- The full ratchet set (`tests/data tests/search tests/scripts`) runs in foreground chunks under 8 minutes each, tee'd to logs, or the launch is handed to the lead. Teammates never background it: the completion notice does not reach them (papercut 2026-10-08, memory 2026-07-21).
- A report that names a commit hash is sent in a LATER tool round than the commit, never in the same batch: the message is composed before the commit output exists, so a batched hash is invented. Seven invented hashes on 2026-10-11 alone (papercut 2026-10-11-main-hash-guessed-in-same-batch-as-commit.md). Paste the hash from the commit's own output or from `git log -1` run after it returned.
- Never kill processes by command-line pattern on this shared machine; kill only PIDs you launched, by PID, or use the harness's task stop (papercut 2026-10-11-main-killed-sibling-pytest-by-pattern.md).

## Papercuts
Friction you work around goes in `docs/papercuts/YYYY-MM-DD-<agent-name>-<slug>.md`, one file per entry (format: the `papercuts` skill; see existing entries). Fill in "Seen before". Commit the one file by path with your work. The lead reviews at 10 unprocessed entries, at handover, or when the user asks; processed entries move to `docs/papercuts/processed/` with a `Disposition:` line.

## Lead address
Teammates message the lead at `main`. Not "team-lead" (unread mailbox) and not the session name (refused).
- The Read tool's `pages` takes one range ("9-12"), not a comma list: "9,12,18" silently returns page 9 only.
