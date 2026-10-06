# A validation slice did not exercise the code path that blocked the real launch

- Date: 2026-10-06
- Agent: twobody-gen-opus
- Seen before: no

What happened: run_942_enumerate.py chose preflight_search(task_no=943 if cell is Jovian else 942). The slices validated before launch (vm, vm2, ev and the earlier gc/ge recall runs) either were #942 cells or ran before the preflight call existed. So the #943 branch first ran at the real gc launch, where the task-number guard (filename must match task number) blocked both shards.
Workaround: a separate scripts/run_943_enumerate.py entry point with its own preflight call; each entry point refuses the other task's cells.
Suggested fix: validate a slice of EVERY cell or branch through the real launch command (same script and flags, small --k or --sample) before sending launch commands. A slice of one cell does not cover per-cell branches.
