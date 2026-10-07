# Teammate reports a commit hash it guessed before the commit returned

- Date: 2026-10-07
- Agent: main
- Seen before: no (three occurrences today from one teammate: ea3ab40c -> 5ab78bb3, 41fe57d0 -> baecaf0d, and 5c41d8ae/a7a76e61 -> 5e4bdc18/fed50171 from a second teammate)

What happened: ci-keeper-opus (twice) and twobody-gen2-opus (twice) sent the lead a commit hash that did not exist, then a correction. The lead pushes by hash, so a wrong hash costs a round trip or, worse, pushes a different commit than the one verified.
Workaround: the lead reads `git log origin/main..HEAD` before every push instead of trusting the reported hash.
Suggested fix: briefs say "report hashes only from the `git log -1` output after the commit returns"; safe-push.sh already lists the commits it would push, which is the real check. Candidate line for teammate-brief.md REPORT section.

Disposition: Promoted: teammate-brief.md REPORT section (hashes only from git log after the commit returns) and coordination.md git rule 7.
