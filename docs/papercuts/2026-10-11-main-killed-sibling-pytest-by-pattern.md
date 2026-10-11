# The lead killed a sibling agent's pytest call while clearing its own stopped suite

- Date: 2026-10-11
- Agent: main (lead)
- Seen before: no (related: the global rule that `pkill -f PATTERN` matches the harness's own shell)

What happened: after stopping a background ratchet suite (TaskStop), `ps` still showed three pytest processes. I killed every PID whose command line contained "pytest". They were not mine: they were another agent's `timeout 45 uv run pytest -q tests/search/test_coorbital_hamiltonian.py`, six seconds old. My own workers had already exited.
Workaround: told the three agents likely to own the call to re-run it if it died at 13:06:40 AEDT.
Suggested fix: on a shared machine, never kill by command-line pattern. Launch long runs under a session marker (`env CYCLERS_RUN=<name>` or a per-run `--basetemp`/log dir) and kill only PIDs whose parent is the stopped task's shell, or let the harness's TaskStop finish the tree. Add the rule to docs/team/coordination.md.
