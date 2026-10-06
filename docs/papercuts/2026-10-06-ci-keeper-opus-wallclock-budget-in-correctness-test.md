# A wall-clock budget inside a correctness test fails it under load

- Date: 2026-10-06
- Agent: ci-keeper-opus
- Seen before: yes (2026-10-05-ci-keeper-opus-ci-timeout-from-shared-mac-load.md, same root: load on the shared Mac)

What happened: `tests/search/test_656_pc_higher_kk_sweep.py::test_656_grid_seed_search_recovers_admitted_pc_32_seed` calls `_grid_seed_search(..., per_call_timeout=5)`. Each corrector call is cut off by SIGALRM after 5 s of wall-clock time. At load 16-20 a call ran past 5 s, the seed was "not recovered", and the test failed. Alone on a quiet machine it passes in 1.5 s. The same SIGALRM budget is used in the production sweep (`pluto_charon_kk_sweep._run_with_timeout`), so a search run under load can also report a false "not found".
Workaround: re-ran the test alone.
Suggested fix: #969. Bound correctness by an evaluation or step count, not by seconds. Keep a wall-clock limit only as a separate runaway guard that is reported as a timeout, never as "not found". Do not raise the timeout.
