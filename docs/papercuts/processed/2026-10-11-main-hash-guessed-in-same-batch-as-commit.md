# Teammates still report guessed commit hashes when the commit and the report go out in one tool batch

- Date: 2026-10-11
- Agent: main (lead)
- Seen before: yes (five teammate hash corrections on 2026-10-07 and 2026-10-08, processed 2026-10-07 with a brief-template fix; five more today from jovian-nbody-opus x3, ci-keeper-opus x1, twobody-gen2-opus x1 (3ece9c8f -> accc51ce, 37f6fc6d -> 8cde2928, 4fde5df4 -> 20930461, bc0e58fd -> c5dc87f6))

What happened: the brief says "hashes only from `git log -1` after the commit". Agents follow it when the commit and the report are separate turns, and break it whenever they batch `git commit` and `SendMessage` in one tool round: the message is composed before the commit output exists, so the hash is invented. The lead then verifies a hash that does not exist and gets a correction a minute later.
Workaround: the lead verifies every reported hash against `git log` before pushing, and agents self-correct.
Suggested fix: a mechanism, not a rule. (1) The brief and docs/team/coordination.md say "a report that names a hash is sent in a LATER tool round than the commit; never in the same batch". (2) A `scripts/commit.sh` wrapper that commits by pathspec and prints `COMMITTED <hash> <subject>` on its own line, so the hash is in the tool output that the next turn reads. (3) The lead refuses to push a hash it cannot find and asks for the `git log -1` line instead of guessing which commit was meant.

Disposition (2026-10-11 19:58, lead): Promoted (coordination.md rule, 294d6a5e: a report naming a hash is sent in a later tool round than the commit) and Backlog #1070 (scripts/commit.sh wrapper printing COMMITTED <hash> <subject>). Seven invented hashes on 2026-10-11; after the rule, ci-keeper and jovian adopted the separate-round habit.
