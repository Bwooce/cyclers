# SendMessage to "team-lead" from in-process teammates lands in an unread mailbox while reporting success

- Date: 2026-10-05
- Agent: main
- Seen before: no

What happened: the lead's briefs told teammates to message "team-lead" (the name the lead's own messages carry). Those sends returned "Message sent to team-lead's inbox" (success) but went to `~/.claude/teams/session-<id>/inboxes/team-lead.json`, which the in-process lead conversation never reads. About 15 reports from 4 agents were lost over ~4 h (receipts, gate verdicts, a 14-item reference check). Only "main" ("queued for the main conversation's next turn") reaches the lead. "cyclers" (the ListAgents name) is refused.
Workaround: briefs now say "message main"; lost reports were recovered from the agents' transcripts (`subagents/agent-a<name>-*.jsonl`, SendMessage tool_use inputs).
Suggested fix: harness should route "team-lead" to the in-process lead, or refuse it like "cyclers". Process fix: every brief requires a receipt within one minute, so a dead route shows at once.
