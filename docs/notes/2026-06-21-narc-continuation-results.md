# 2026-06-21 — #388 N-arc real-ephemeris continuation batch results

**Verdict: terminal negative (held, no writeback).** No descriptor-bearing
catalogue row closes on DE440 within 0.5 km/s of its sourced V∞ anchors with
bend-feasible flybys. Failures occur at an **intermediate ramp rung**, not at the
DE440 endpoint — the genuine in-basin Russell parents never reach lambda=1.

- Script: `scripts/narc_continuation_batch.py`
- Lane: `src/cyclerfinder/search/narc_continuation.py` (#388)
- Runlog: `data/runs/narc-continuation-20260620T190007Z.jsonl`
- Run wall time: ~72 s total (all 4 rows); no timeout reduction was needed.
- **Launch-window setting used: `range(1, 22)`** (Russell §5.3 LaunchWindow
  1..21, grid=100). The fallback `range(1, 6)` reduction was NOT required.
- `target_phase = 0.0` (documented v1 simplification: the parent's exact
  beginning E-M relative phase is not exposed by the idealized Cycler; the
  LaunchWindow scan over 21 candidate epochs is relied on to find a closing
  epoch).
- V∞ cap: `VINF_CEILING_KMS["M"]` ≈ 58.25 km/s (the inner corrector takes a
  single scalar; Mars is the binding ceiling of the E/M encounters here).

## Per-row outcome

| id | V-level | converged | max residual (km/s) | emerged V∞ (km/s) | anchor E / M (km/s) | anchor-match | bend-feasible | winning epoch (s) |
|---|---|---|---|---|---|---|---|---|
| mcconaghy-2006-em-k2  | V0 | False | inf | [] | 4.70 / 5.00  | False | False | 7.237e8 |
| russell-ch4-4.991gG2  | V3 | False | inf | [] | 4.99 / 5.10  | False | False | 7.237e8 |
| russell-ch4-8.049gGf2 | V3 | False | inf | [] | 8.05 / 10.02 | False | False | 7.237e8 |
| russell-ch4-9.353Gg2  | V1 | False | inf | [] | 9.35 / 10.52 | False | False | 7.237e8 |

Every row produced `max_residual_kms = inf` with empty `emerged_vinf_kms`: the
driver's `best` stayed `None`, i.e. **no candidate epoch's homotopy walk produced
any final-rung result**. The reported `t0_sec`/`winning_epoch_sec` is the default
fall-through (`epochs[0]`), not a converged solution.

## Where the ramp breaks (elastic-mode gap)

Direct rung-by-rung probe of `russell-ch4-4.991gG2` (seq `E-E-M-M`, revs
`(1, 0, 1)`, branch `(low, single, low)`, seed ToFs `(533.70, 150.0, 1026.21)` d,
period 1559.86 d), best epoch, 3 ramp rungs:

```
rung 0: ok  res=inf  conv=False
rung 1: ok  res=inf  conv=False
rung 2: RAISED ValueError: Method 'lm' doesn't work when the number of
        residuals is less than the number of variables
```

So the walk **never reaches the DE440 (lambda=1) endpoint**. The two earlier
rungs "do not raise" but already fail to close (`res=inf`); the final intermediate
rung raises a least-squares dimensionality error (residual count < variable
count) which the driver catches and treats as divergence, aborting that epoch.
This is an **intermediate-ramp-rung failure**, characteristic of an elastic-mode
gap in the homotopy ladder rather than a DE440-endpoint reach failure.

## Held verdict

This is the **honest characterized terminal negative for #388**: the N-arc
real-ephemeris continuation lane, driven over every descriptor-bearing catalogue
row across Russell's full LaunchWindow 1..21 epoch scan, finds **no DE440 closure
matching the sourced V∞ anchors with bend-feasible flybys**. Notably the V3
regression rows (`russell-ch4-4.991gG2`, `russell-ch4-8.049gGf2`) — genuine
in-basin Russell parents — do not survive the homotopy ramp either.

No row closes ⇒ no PROPOSED V0→V1. **HELD, no writeback.** `data/catalogue.yaml`
and `validate.py` were not touched.

The failure is structurally located (intermediate ramp rung, least-squares
under-determination at the final intermediate ephemeris), which points the next
iteration at the ramp ladder / corrector free-variable count — not at the DE440
endpoint physics.
