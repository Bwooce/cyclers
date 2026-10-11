# Background CI-watcher completion notices arrived about an hour late, so CI results reached the lead late

- Date: 2026-10-11
- Agent: ci-keeper-opus
- Seen before: no (the global rule "completion notices ... can arrive an hour late or never" names the cause; this is its first entry here)

What happened: I watched CI runs with `run_in_background` loops that exit when a run completes, and waited for their completion notices. Run 38107321585 finished at 15:02 AEDT, but its notice reached me only at 16:00, after the lead asked for a status. The lead had seen no report since 14:11. Every later run was cancelled by newer pushes through the concurrency group, so the one green run was the only result, and it was late.
Workaround: switched to the Monitor tool, which streams each status change as an event, and reported the missed results.
Suggested fix: for any result someone is blocked on, use a streaming monitor (one event per status change), not a single completion notice, and check the output file yourself at the expected finish time. Also: hold pushes while waiting for one run, because each push cancels the queued run (the lead now does this).
