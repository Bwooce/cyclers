# Russell-12 like-for-like — #137 free-return genome results (2026-06-07)

Recovered run (the #137 agent was killed before reporting; its fix commits
37b9bbd/8fc2ffb survived — this note records the re-run of its committed
campaign). Circular model, like-for-like, `--genome free-return`,
128 epochs × 14 workers, ~5 s wall.

## Acceptance gate
`tests/search/test_russell12_likeforlike_probe.py` — **4 passed** (10.5 s):
the radial-crossing genome makes the sourced geometry representable
(truth-residual ≈ 0 where the old genome measured 3.2–37.5 km/s).

## Per-row outcomes (V∞ DERIVED, never imposed — non-circular evidence)
| row | outcome |
|---|---|
| mcconaghy-2006-em-k2 | CLOSE-AND-MATCH |
| russell-ch4-4.991gG2 | CLOSE-AND-MATCH |
| russell-ch4-8.049gGf2 | CLOSE-AND-MATCH |
| russell-ch4-9.353Gg2 | CLOSE-OFF-ANCHOR |
| russell-ch4-3.64gGg3 | CLOSE-AND-MATCH |
| russell-ch4-3.78Gg3 | CLOSE-AND-MATCH |
| russell-ch4-5.30gGf3 | CLOSE-AND-MATCH |
| russell-ch4-9.94Gg3 | CLOSE-AND-MATCH |
| russell-ch4-3.66gfF3 | CLOSE-MATCH-SYMMETRIC-ONLY |
| russell-ch4-5.30ggF3 | CLOSE-MATCH-SYMMETRIC-ONLY |
| russell-ch4-5.75ggF3 | CLOSE-MATCH-SYMMETRIC-ONLY (vinf Δ 0.04/0.01) |
| russell-ch4-6.44Gg3 | NO-CLOSE (res 0.2666 — near the 0.1 floor) |

**Counts: 7 CLOSE-AND-MATCH, 3 SYMMETRIC-ONLY, 1 OFF-ANCHOR, 1 NO-CLOSE**
(was 0/9-11/1-3 on the lambert genome, #125/#135).

## Caveats
- The runner's `--probe-at-truth` section still prints the OLD lambert-genome
  diagnostics (script wiring: the probe path predates `--genome`); the
  authoritative truth-residual≈0 pin is the committed acceptance test.
- NO validation_level writeback performed: whether like-for-like
  circular-coplanar reproduction of a circular-coplanar source meets §14 V1
  is a flagged USER decision (the §14 V1 text — lamberthub re-solve + Kepler
  re-prop — is fidelity-agnostic and could be applied to these closed
  trajectories as the next step if approved).
- Implication: the Jones elimination chain's negatives were computed on the
  old genome — re-validation on the free-return genome is now warranted
  (#133), as is the near-miss survey.
