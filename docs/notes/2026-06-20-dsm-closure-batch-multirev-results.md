# #388 — Multi-rev descriptor-seed lane: large improvement, closure still blocked (honest negative)

Date: 2026-06-20. Outcome of wiring multi-revolution Lambert + sourced seed ToFs into
the DSM descriptor-seed closure lane (spec `2026-06-20-multirev-descriptor-seed-design.md`;
commits `8871ba7`, `506a197`, `50b2e10`) and re-running `scripts/dsm_closure_batch.py`
on real DE440. **No catalogue writeback, no promotion.**

## Result — residuals cut 2.5–4×, rev selection correct, but 0/4 converge

| row | level | single-rev res (baseline) | **multi-rev res** | revs | DSM dV (was → now) | match |
|---|---|---|---|---|---|---|
| mcconaghy-2006-em-k2 | V0 | 40.1 | **9.05** | [1,0,1] | 102.6 → 19.6 | True* |
| russell-ch4-4.991gG2 | V3 | 35.9 | **9.18** | [1,0,1] | 95.4 → 20.0 | True* |
| russell-ch4-8.049gGf2 | V3 | 23.0 | **14.13** | [1,0,0] | 51.2 → 30.5 | False |
| russell-ch4-9.353Gg2 | V1 | 76.9 | **29.47** | [1,0,0] | 131.5 → 39.7 | False |

(* `match=True` is coincidental on a non-converged solve — the gate is
`converged AND anchor_match`, which is empty.) Convergence tol = 0.1 km/s.

**The fix works exactly as designed.** The rev selection is correct: the resonant
same-body legs (E→E, M→M) now pick the 1-rev branch (`revs=[1,…,1]`) and the transit
leg picks 0 — the degenerate single-rev collapse is gone. Residuals dropped 2.5–4×
and DSM dV dropped 3–5×. This validates the literature diagnosis (Russell's `i` rev
index; the #153 single-rev-on-a-resonant-leg degeneracy) and the sourced seeding.

**But the rows still do not close.** Residuals floor at 9–29 km/s, far above the
0.1 km/s tolerance. Multi-rev representation was necessary but not sufficient.

## Seedable dropped 9 → 4 (a correctness improvement, not a regression)

The sourced arc→leg mapping refuses to seed a row when it cannot map each resonant
same-body leg to a published arc ToF. The 5 rows that dropped
(`3.64gGg3`, `3.78Gg3`, `5.30gGf3`, `9.94Gg3`, `6.44Gg3`) all have sequence
`E-E-E-M-M` — **three** same-body resonant legs — but the catalogue publishes only
**two** arc ToFs. The old slack heuristic silently fabricated a seed by reusing one
ToF; the new lane correctly returns `None` (no source for the third resonant leg) and
classifies them NO-CLOSE. This is the honest behavior — it does not invent a leg ToF
that the source does not provide.

## Why closure is still blocked (working hypothesis — model jump, not genome)

With the rev count correct and the seed at the published resonant geometry, the
remaining 9–29 km/s residual is most consistent with the **circular-coplanar →
real-ephemeris model jump**, not a genome deficiency:

- Russell 2004 §5.4.1 explicitly warns the circular-coplanar parent **cannot be jumped
  straight to the accurate ephemeris** — "the gap … prohibits an immediate jump." His
  cyclers reach the real model only via **continuation/homotopy** (fidelity raised in
  steps, each solution seeding the next) plus **multiple-shooting with SNOPT elastic
  mode** that keeps minimizing constraint violation even through intermediate
  infeasibility (§5.4.3–5.4.4).
- Our lane jumps directly from the published (circular-coplanar / mean-element J2000)
  descriptor to DE440 and drives it with `scipy.least_squares` (TRF) — **no homotopy,
  no elastic mode**. The residual floor at correct resonant geometry is the signature
  of that missing continuation, the gap Russell names.

This is a HYPOTHESIS consistent with the evidence; it is not proven here. What IS
established: rev count was a real, large obstruction (now removed) and is not the
whole story.

## Status — characterized partial negative; stop, do not loosen

Per the spec's error-handling gate: the validated regression rows still do not
reconverge, so the obstruction is deeper than rev count — **stop; do not loosen the
tolerance to manufacture a closure.** The multi-rev lane is real, tested, committed
infrastructure and a materially better closer (residuals cut 2.5–4×), but it does not
close these rows on DE440. `mcconaghy-2006-em-k2` stays V0. No writeback.

### Next venue (future #388/#307, not done here)
**Continuation homotopy** circular-coplanar → DE440 (Russell §5.4.1/§5.4.4): raise
ephemeris fidelity in steps, re-seeding each from the last, and/or adopt an
elastic-mode-capable solver that survives intermediate infeasibility. That is a larger
build than this wiring change and is the logical next rung if these promotions are
pursued further.

## References
- `docs/superpowers/specs/2026-06-20-multirev-descriptor-seed-design.md`
- `docs/notes/2026-06-20-dsm-closure-batch-results.md` — the single-rev baseline this improves on.
- `docs/notes/2026-06-07-russell-2004-dissertation-method-mining.md` §5.4 — the model-jump / continuation / elastic-mode method.
- runlog: `data/runs/dsm-closure-20260620T221306.jsonl`
