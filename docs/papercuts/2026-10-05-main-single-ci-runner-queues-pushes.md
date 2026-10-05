# One self-hosted CI runner on the busy dev Mac: each push queues a ~50 min run behind the last

- Date: 2026-10-05
- Agent: main (ci-keeper-opus measured it)
- Seen before: no

What happened: a docs-only push waited 52 min in the queue, then ran ~55 min at load average 33-37 on 8 cores, competing with local agents and test runs.
Workaround: batch pushes every 1-2 h or at milestones.
Suggested fix: skip the full suite for docs-only pushes (paths-ignore), or cancel superseded queued runs (concurrency group with cancel-in-progress).
