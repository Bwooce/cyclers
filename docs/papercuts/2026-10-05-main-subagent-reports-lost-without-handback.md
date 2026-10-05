# Sub-agents of a teammate end "without delivering a report", and teammates cannot name sub-agents

- Date: 2026-10-05
- Agent: main (reported by corpus-review-fable)
- Seen before: no

What happened: all six reader sub-agents ended without a SubagentHandback report; their output survived only because the brief made them write to scratch files incrementally. Passing `name` for a sub-agent was refused ("team roster is flat").
Workaround: write-to-file-as-you-go in every sub-agent brief.
Suggested fix: keep this in the team-lead brief template for any agent allowed to spawn.
