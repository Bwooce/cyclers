# #513 — R52-U recovery: PARTIAL flip (2 of 3 C32-dominance metrics)

**Date:** 2026-07-01. **Status:** R52-U recovered from a SOURCED IC; C32-dominance gate
partially flips. Betweenness residual left open (not root-caused this session).

## What happened
The user spotted that BinBraik's own repo (adopted #495) already contains R52-U's exact
converged IC — sitting unused in `data/golden/braik_ross_2026_em_family_ics.yaml`:
`x0=-0.2719428329684943, ydot0=-2.137466033981196, jacobi=3.1294, period_days=56.436`
(matches Table 2's period exactly). No blind multiple-shooting search was needed — this IC
was added directly to `OFFLINE_SEEDS` in `reachable_representatives.py` and to
`_recover_subset()`'s 13-member list.

## Result
With R52-U present, the 13-node network scorer output:

| metric | winner | value | 2nd place |
|---|---|---|---|
| strength | **C32** | 3.038 | C11a 2.581 |
| harmonic-closeness | **C32** | 3.362 | R52-U 3.074 |
| betweenness | C21 | 0.3788 | C32 0.3485 (close, ~9% gap) |

**C32 now wins 2 of 3 metrics** (strength + harmonic), up from 0 of 3 before R52-U was
recovered (#497's diagnosis). The gap that remains is narrow — C32 is a close second on
betweenness, not far off as it was pre-recovery.

## Test changes
- `test_validation_gate_c32_dominant`: stays `xfail` (all 3 metrics still required), reason
  updated to reflect the partial flip and pin the exact remaining numbers.
- `test_validation_gate_c32_undominant_faithful_negative`: REWRITTEN — its old premise (C32
  wins nothing) is now false. Now asserts the TRUE current state: C32 wins strength+harmonic,
  C21 wins betweenness by a narrow (~9%) margin. Still a faithful-negative-style test (records
  reality, doesn't force a pass) — just recording a different reality now.
- Module docstring updated: 12-node → 13-node (full source-confirmable set, no member excluded).

## Not investigated (honest scope boundary)
The betweenness residual's root cause was NOT determined this session. Candidate causes
(untested): residual proxy-fidelity specific to the R52-U-C21-C32 shortest-path routing, or a
genuine second-order network-topology difference from Braik's exact grid/horizon choices. This
is a smaller, more tractable gap than the "missing node" diagnosis was — worth a future task,
not pursued further here given the session was closing.

## Verification
`uv run ruff check` + `ruff format --check` + `mypy src tests` all clean.
`pytest tests/search/test_reachable_network_gate.py -m ""` (all markers, incl. slow):
**7 passed, 1 xfailed** (the expected, now-narrower gap). No regressions.
