# subprocess.run(timeout=...) on `uv run python ...` leaves the python child running as an orphan

- Date: 2026-10-11
- Agent: twobody-ext-opus
- Seen before: not recorded; it can have happened in every earlier `scripts/ev_realeph_1025.py` timeout (the rung (d) Standish and DE440 runs used the same call)

What happened: `ev_realeph_1025.py run_route` called the chain tool through `subprocess.run(["uv", "run", "python", ...], timeout=150)`. On the timeout, Python kills only the direct child (`uv`). The python grandchild ran on for 26 more minutes on a machine at load 60-150, then wrote `realeph_chain.json` into the route directory. The orchestrator had already recorded the route as NOT RUN (timeout), and the next resumable call read the late file back as a finished route (ev5-9 S2 e0 in the #1025 sec. 11 V2 run). The static `tool.log` held only the log up to the timeout.
Workaround: `run_route` now starts the tool in its own session (`start_new_session=True`) and kills the whole process group with `os.killpg(..., SIGKILL)` on a timeout; checked with a 3-s timeout: no process left. The one late result is recorded and disclosed in the #1025 rung-d note sec. 11.6.
Suggested fix: any script in the repo that wraps `uv run` in `subprocess.run(timeout=...)` should do the same (grep `subprocess.run(` with `timeout=` in scripts/). Earlier "NOT RUN (timeout)" outcomes of such scripts may have orphans that wrote results after the verdict was read; check their directories for result files newer than their `tool.log`. A scan of every `data/1025_ev_realeph/**/tool.log` that starts with "TIMEOUT after" found no other late result.

Disposition (2026-10-11 19:58, lead): Fix now, #1068 (twobody-gen2-opus, dispatched 19:50: shared process-group kill helper, scan of scripts/, check of past output directories for late results).
