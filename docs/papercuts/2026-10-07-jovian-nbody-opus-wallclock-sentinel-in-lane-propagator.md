# The Jovian lane's propagator turns a wall-clock stop into a defect value

- Date: 2026-10-07
- Agent: jovian-nbody-opus
- Seen before: yes (#969; the 5 s per-call timeout and the 600 s pytest timeout papercuts processed 2026-10-07)

What happened: `JovianRestrictedNBody.propagate` (src/cyclerfinder/nbody/jovian.py) defaults to
`max_wall_sec=60`; a leg that runs longer returns `converged=False`, and `jovian_defect_residual` /
`_cycle_residual` then put a 1e9 / 1e7 sentinel into the residual. On a loaded machine (load 30-120
during this session) a slow leg becomes a huge "defect", indistinguishable from divergence, inside a
least-squares solve. No wrong result was traced to it this session; my runs passed `max_wall_sec`
of 3000 s and logged a timeout separately.
Workaround: pass a large `max_wall_sec` and report a stop as "timeout", never as a defect value.
Suggested fix: under #969, an integration-step budget or a timeout reported outside the residual;
past Jovian negatives computed under load (#480, #318, #501, already void for the wrap bug) should
not be re-run with the 60 s default.
