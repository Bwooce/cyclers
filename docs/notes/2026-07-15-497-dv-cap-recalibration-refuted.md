# #497: DV_CAP_MS recalibration does NOT close the C32-betweenness gate

**Date:** 2026-07-15
**Origin:** `data/OUTSTANDING.md`'s #497 entry proposed "set `DV_CAP_MS` ≈ 51 m/s (from #495: C32
wins centrality in 80% of Braik's budget sweep there; our xfail gate uses 409.3 = near-full
connectivity, no betweenness signal) → flip the C32-dominance gate." Investigated directly before
dispatching, since the actual test file (`tests/search/test_reachable_network_gate.py`) showed the
task's premise was already stale: #513's later R52-U node recovery got the gate to 2/3 metrics
(strength, harmonic-closeness) matching Braik's Table 4, leaving ONLY betweenness unresolved (C21
narrowly beats C32, ~9% margin) — and that residual was explicitly still computed at the ORIGINAL
409.3 m/s cap, meaning #497's proposed additional fix (recalibrating the cap) had never actually
been tried in combination with the R52-U-recovered network.

## Method

Reused the xfail gate test's own `_recover_subset()` (13-node source-confirmable representative
set, R52-U included) and reproduced its `_run_network()` construction (`VoxelGrid`,
`build_reachable_set`, `proxy_matrix`), parametrizing only `dv_cap_ms` in `apply_budget_cap()`.
Swept 9 cap values from 51.16 m/s (Braik's own reference cap, per #495/#497b) up to 409.3 m/s
(the current gate's default). Each run takes ~15s — cheap enough to sweep directly rather than
guess.

## Result: clean, decisive negative

```
cap=  51.16  bw_winner=C11a    C32=0.1364  C21=0.3333  C11a=0.3788
cap=  75.00  bw_winner=C21     C32=0.2273  C21=0.3333  C11a=0.2879
cap= 100.00  bw_winner=C21     C32=0.2273  C21=0.3333  C11a=0.2879
cap= 150.00  bw_winner=C21     C32=0.2273  C21=0.3333  C11a=0.2879
cap= 200.00  bw_winner=C21     C32=0.3485  C21=0.3788  C11a=0.3182
cap= 250.00  bw_winner=C21     C32=0.3485  C21=0.3788  C11a=0.3182
cap= 300.00  bw_winner=C21     C32=0.3485  C21=0.3788  C11a=0.3182
cap= 350.00  bw_winner=C21     C32=0.3485  C21=0.3788  C11a=0.3182
cap= 409.30  bw_winner=C21     C32=0.3485  C21=0.3788  C11a=0.3182
```

**C32 never wins betweenness at any swept cap value.** At the tightest, Braik-matching cap
(51.16 m/s) it actually gets WORSE — C32's betweenness drops to 0.1364 (from 0.3485 at loose
caps) and C11a takes over as the winner entirely, rather than closing the gap toward C32. This
refutes #497's proposed fix outright, not just leaves it unconfirmed: recalibrating `DV_CAP_MS`
is not the missing piece.

## Interpretation

This confirms, rather than merely leaves open, the xfail test's own speculation about the
residual cause: "candidate causes: residual proxy-fidelity on that specific path, or a genuine
second-order network-topology difference from Braik's grid/horizon choices" (docstring, prior to
this investigation, framed as untested candidates — now one of them, the cap-value hypothesis
implicit in #497, is explicitly ruled out). The real gap is almost certainly in how this
project's `build_reachable_set`/`proxy_matrix` approximates pairwise transfer cost/connectivity
relative to Braik's own full DVmatrix computation — a genuine methodology difference, not a
tunable parameter.

## Disposition

- No code change (this is a negative-result investigation, not a fix).
- `tests/search/test_reachable_network_gate.py`'s xfail reason updated to record that the
  cap-recalibration hypothesis was tested and refuted, so a future reader doesn't re-propose it.
- `data/OUTSTANDING.md`'s #497 entry corrected from "active follow-on" to closed/refuted.
- The actual remaining root cause (proxy-fidelity vs. Braik's DVmatrix, or grid/horizon
  differences) is NOT resolved here — flagged as the genuine next step if anyone wants to pursue
  closing this specific gate further, but that's a real investigation (comparing the proxy
  reachable-set computation edge-by-edge against Braik's own DVmatrix.csv), not a quick follow-on.
