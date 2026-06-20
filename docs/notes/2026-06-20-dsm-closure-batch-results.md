# #404/#388 — DSM multi-arc closure-lane batch: structural genome mismatch (clean negative)

Date: 2026-06-20. Outcome of running the DSM η-leg closure lane (spec
`2026-06-10-dsm-multiarc-closure-lane-design.md`, Components 1–2 + the Component 4
batch driver `scripts/dsm_closure_batch.py`) over every descriptor-bearing
catalogue row on the real DE440 ephemeris. Honest negative; **no catalogue
writeback, no promotion.**

## Result — 0 / 9 converge

| row | level | converged | residual (km/s) | anchor | match | DSM dV (km/s) |
|---|---|---|---|---|---|---|
| mcconaghy-2006-em-k2 | V0 | False | 40.1 | 5.00 | (coinc.) | 102.6 |
| russell-ch4-4.991gG2 | V3 | False | 35.9 | 5.10 | (coinc.) | 95.4 |
| russell-ch4-8.049gGf2 | V3 | False | 23.0 | 10.02 | False | 51.2 |
| russell-ch4-9.353Gg2 | V1 | False | 76.9 | 10.52 | False | 131.5 |
| russell-ch4-3.64gGg3 | V1 | False | 29.5 | 4.59 | False | 92.0 |
| russell-ch4-3.78Gg3 | V1 | False | 20.0 | 4.63 | False | 64.1 |
| russell-ch4-9.94Gg3 | V1 | False | 49.6 | 10.76 | False | 151.2 |
| russell-ch4-5.30ggF3 | V1 | False | 20.0 | 5.44 | (coinc.) | 47.8 |
| russell-ch4-6.44Gg3 | V1 | False | 33.1 | 3.74 | False | 74.8 |

(`russell-ch4-5.30gGf3`, `5.75ggF3`: descriptor present but the coplanar arc does
not reach the body — the spec's OFF-FAMILY-NO-CLOSE case; correctly skipped, no seed.)

**`mcconaghy-2006-em-k2` does NOT promote** (the only V0 in-scope #365 target stays V0).
The two `match=True` cells are coincidental — an emerged V∞ happening to land within
0.5 km/s of the anchor on a *non-converged* solve carries no meaning; the gate
`promote = converged AND anchor_match` is correctly empty.

## Root cause — single-rev Lambert genome cannot express the multi-rev resonant legs

The corrector is sound: `test_dsm_leg.py::...` (line ~174) proves
`dsm_chain_correct` converges to <0.1 km/s on a transfer constructed by design. The
failure is **upstream, in the seed geometry**, and it is structural — not a tuning gap.

Concrete seed dump (`seed_dsm_chain_from_descriptor`):

```
mcconaghy-2006-em-k2  seq ('E','E','M','M')
  tof seeds (days) = [883.0, 143.2, 883.0]   # legs E→E, E→M, M→M
russell-ch4-4.991gG2  seq ('E','E','M','M')
  tof seeds (days) = [857.5, 168.7, 857.5]
```

The sequences are `E-E-M-M`: the **same-body legs** (E→E, M→M) carry ~880-day ToFs,
i.e. **~2.4 heliocentric revolutions**. The DSM-leg genome
(`dsm_leg.dsm_chain_correct`) models each leg as a **single-revolution** Lambert arc
plus an interior impulse. A 2.4-rev same-body return forced through a single-rev
Lambert is a wildly wrong arc → emerged V∞ and DSM impulses of tens-to-hundreds of
km/s → residuals floor at 20–77 km/s. The corrector cannot recover because the
*model itself* cannot represent the leg.

This is confirmed by the regression rows: **validated V3 cyclers do not reconverge
through this lane either.** Their geometry is the same multi-rev resonant free-return
(the `g`/`G` arcs are 1–3-yr same-body resonant returns by construction), so the
single-rev genome mis-expresses them identically. The lane reproduces *nothing*,
known-good or candidate — the tell that the obstruction is the genome, not the seed
tuning or the corrector.

## Conclusion

The descriptor's `free_return_arcs` are **multi-revolution resonant same-body
returns**; the DSM η-leg genome is **single-revolution Lambert legs**. The two do not
match, so the DSM closure lane cannot close these rows on the real ephemeris.
`mcconaghy-2006-em-k2` stays V0; the #365 negative is *not* cleared by this lane.

### What this rules in / out
- **Corrector / DSM machinery**: sound (unit-test-proven). Not the blocker.
- **Descriptor→seed bridge**: produces geometrically correct *resonant* arcs but
  hands them to a genome that can only integrate single-rev Lambert legs.
- **The lane as a #365 unblock**: does NOT deliver. The #404 triage already showed only
  1 of the 4 ocampo negatives is even descriptor-bearing; this run shows that one row
  also fails, for the structural reason above.

### Next venues (future #388, not done here)
1. **Multi-rev Lambert in `dsm_leg`**: allow each leg a revolution count `N_rev` and
   solve the matching multi-rev branch. The descriptor already implies `N_rev` from the
   g-arc ToF / body period. This is the direct fix — the genome gains the DOF the geometry needs.
2. **Close in the resonant-arc genome natively**: `self_seeding.g_arc_branches` already
   *builds* these arcs as resonant returns. Correct/validate there instead of round-tripping
   through the single-rev DSM genome. Avoids the model mismatch entirely.
3. Treat the descriptor-bearing `ch4` rows as already-validated (they are V1/V3 by their
   sourced state) and stop trying to reconverge them through an inadequate lane — the
   batch is then purely a *regression guard* that documents the genome boundary.

No writeback. The single V0 promotion target remains V0 pending a genome that can
represent its legs.
