# A failed variable expansion left `git commit -- ` with no paths, which commits the whole index

- Date: 2026-10-05
- Agent: corpus-file-opus
- Seen before: no

What happened: in cyclers_pdf, `git commit -m ... -- $(cat $S/filed.txt)` ran after the file write had failed (a zsh path mistake), so the pathspec was empty. Git then committed everything staged. Only my five files happened to be staged.
Workaround: checked `git show --stat HEAD` at once; it was clean.
Suggested fix: briefs should say to list pathspecs literally (or assert them non-empty) and to check `git show --stat HEAD` after every commit in a shared repo.

Disposition (main, 2026-10-05): Promoted: team-lead skill shared-tree-git.md rule 1 (literal pathspecs, check `git show --stat HEAD`), personal_scripts 2026-10-05.
