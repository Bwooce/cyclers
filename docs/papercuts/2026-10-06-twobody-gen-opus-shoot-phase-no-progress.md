# The chain tool's shoot phase ran 9 minutes with no output

- Date: 2026-10-06
- Agent: twobody-gen-opus
- Seen before: yes (same class as the global "instrument long runs" rule)

What happened: the blend phase of scripts/run_942_realeph_chain.py printed one line per lambda step, but the lambda = 1 shoot phase printed only on a successful restart. A 10-cycle GanCal#1 pilot used its whole 9-minute timeout with no output, so I could not tell slow from hung.
Workaround: per-restart result lines (including non-closures with residual, nfev and time), an evaluation heartbeat, and a checkpoint of the shoot start state per epoch (commit 1b8cfbfa).
Suggested fix: when adding a new phase to an instrumented script, give it its own progress lines before the first pilot.
