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

## Part 1 — §14 V1 writeback (#137 follow-up, 2026-06-07)

USER-approved principle: like-for-like reproduction qualifies for §14 V1 when the
§14 mechanics pass. Applied the LITERAL §14 V1 mechanics to each matched row's
closed free-return geometry, like-for-like on the CIRCULAR ephemeris (a
circular-coplanar reproduction of a circular-coplanar source):

* path (a) — every leg re-solved with `lamberthub` izzo2015 + gooding1990,
  agreement < `V1_TOLERANCE_MPS` (1e-3 m/s);
* path (c) — Kepler forward re-propagation residual < `KEPLER_REPROP_TOL_KM` (1 km);

both reused verbatim from `cyclerfinder.verify.agreement.crosscheck_code_paths`
(§14 V1 is exactly paths a+c; path b is an extra Forge witness, not part of §14).
Module: `src/cyclerfinder/search/free_return_v1.py`; evidence test:
`tests/search/test_free_return_v1_mechanics.py` (10 passed).

**The decisive honesty gate.** A single free-return ellipse forms a CLOSED,
V_inf-continuous E->M->E cycler only when its descending Earth crossing coincides
with where Earth actually is on the circular orbit. The return-Earth phase mismatch
splits the rows cleanly:

| row | outcome | return-Earth Δphase | Mars V_inf continuity | §14 V1 |
|---|---|---|---|---|
| russell-ch4-5.30gGf3 | CLOSE-AND-MATCH | -0.2° | 0.01 km/s | **PASS → V1** |
| russell-ch4-9.94Gg3 | CLOSE-AND-MATCH | -0.5° | 0.04 km/s | **PASS → V1** |
| russell-ch4-5.75ggF3 | SYMMETRIC-ONLY | -6.1° | 0.18 km/s | **PASS → V1** |
| mcconaghy-2006-em-k2 | CLOSE-AND-MATCH | 175.2° | 24.35 km/s | refused |
| russell-ch4-4.991gG2 | CLOSE-AND-MATCH | 176.0° | 24.58 km/s | refused |
| russell-ch4-3.64gGg3 | CLOSE-AND-MATCH | 179.1° | 22.79 km/s | refused |
| russell-ch4-3.78Gg3 | CLOSE-AND-MATCH | 179.2° | 23.52 km/s | refused |
| russell-ch4-3.66gfF3 | SYMMETRIC-ONLY | 175.9° | 21.93 km/s | refused |
| russell-ch4-5.30ggF3 | SYMMETRIC-ONLY | 170.5° | 23.68 km/s | refused |
| russell-ch4-8.049gGf2 | CLOSE-AND-MATCH | n/a | Lambert-singular | refused |

For the three PASS rows the reconstructed E->M->E arc is a genuine closed cycler:
lamberthub agreement ~1e-8 m/s (≪ 1e-3), Kepler reprop ~3-9e-4 km (≪ 1 km), Mars
V_inf continuous to ≤ 0.18 km/s. For the six 175-180° rows the single free-return
ellipse is genuinely MULTI-ARC — its return needs the eliminated phasing loops, so
forcing a Lambert return leg lands a ~24 km/s Mars V_inf discontinuity (a broken,
non-physical trajectory whose lamberthub/Kepler self-consistency would pass
VACUOUSLY). `VINF_CONTINUITY_TOL_KMS = 0.5` (the campaign match tolerance, not
tightened) rejects exactly that. `8.049gGf2` (deep aphelion, 93-d transit) is
single-rev-Lambert-singular on reconstruction. None of the six promote.

The evidence chain per PASS row: SOURCED aphelion + transit seed → free-return
closure → DERIVED V_inf matches sourced (already shown) → §14 V1 mechanics pass on
a closed, V_inf-continuous reconstruction. **Scope (honest):** this is a
circular-coplanar reproduction of a circular-coplanar source — like-for-like; it is
NOT a real-ephemeris (V3) result.

Writeback (same commit): `_LEVEL_EVIDENCE` extended (3 new entries),
`_LEVEL_BY_ID` extended, `scripts/backfill_validation_level.py` applied
(V0→V1 on the 3 rows), the v4.5 census ratchet
(`tests/data/test_schema_v45_fields.py`) and `data/README.md` updated.

**Validation census (the headline): 5 rows are V1** — the Aldrin pair
(outbound + inbound, #125) and the three Russell free-return rows above (#137
Part 1). Everything else is V0 / untagged. No row claims V2+.

## Caveats
- The runner's `--probe-at-truth` section still prints the OLD lambert-genome
  diagnostics (script wiring: the probe path predates `--genome`); the
  authoritative truth-residual≈0 pin is the committed acceptance test.
  (Part 3 fixes this — the probe now honours `--genome`.)
- Implication: the Jones elimination chain's negatives were computed on the
  old genome — re-validation on the free-return genome is now warranted
  (#133), as is the near-miss survey. (Part 2 assesses the extension.)
