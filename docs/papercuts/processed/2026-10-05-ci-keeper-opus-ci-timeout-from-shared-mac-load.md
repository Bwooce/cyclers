# CI on the shared dev Mac times out a long test when local runs and other projects load the machine

- Date: 2026-10-05
- Agent: ci-keeper-opus
- Seen before: yes (processed/2026-10-05-main-single-ci-runner-queues-pushes.md, same root: one busy machine runs CI)

What happened: CI run 37291468519 was a docs-only push. It failed with `test_known_close_pair_73c_plateaus_just_outside_guard` at the 600 s pytest-timeout. The test dates from 2026-08-08 and passed in every earlier green run. The run took 63 min against 43-50 min for green runs. At the time the load average was 43-49 on 8 cores: a local full suite (run2), a non-cyclers next-server at 275 % CPU, and agent work. A red CI that a docs-only push caused costs a log read and a re-run every time.
Workaround: read the log, timed the test alone, re-ran CI.
Suggested fix: one full suite at a time on this Mac, CI included (hold pushes during a local full run, or verify through CI instead of locally); `--durations=25` in the CI pytest step so a slowdown is visible before it becomes a timeout.

Disposition: Backlog: #969 (wall-clock budgets) and Promoted: one full suite at a time incl. CI is in docs/team/coordination.md; Fix now: `--durations=25` added to the CI pytest step (lead, this commit).
