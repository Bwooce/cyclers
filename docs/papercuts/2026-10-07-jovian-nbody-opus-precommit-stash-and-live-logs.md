# Pre-commit stashes other agents' unstaged files; a tracked live log turned that into a blocker

- Date: 2026-10-07
- Agent: jovian-nbody-opus
- Seen before: yes (docs/team/coordination.md, known hazards: pre-commit stash on a shared tree)

What happened: every pathspec commit printed "Stashing unstaged files ... Restored changes": the
pre-commit hook stashes every other agent's uncommitted work for the length of the hook run. My
control script appended its run log every 2 s to a TRACKED file (data/968_control/runlog.txt), so
other agents' hooks saw the tree change under them and mypy failed with "files were modified by
this hook" (earthmoon-opus, four tries). The lead untracked the file (4729072f); that also removed
it from disk, and the static copy was recovered from git history (c3ce5aca).
Workaround: live logs now go to data/968_control/live/ (gitignored by `data/*/live/`); only a final
static copy is committed.
Suggested fix: briefs say "a running process never appends to a tracked file; live logs under
data/*/live/ or the scratch dir". Longer term, a hook mode that does not stash (or runs on the
staged pathspec only) would remove the race for everyone.
