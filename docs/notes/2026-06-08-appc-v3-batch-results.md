# App-C V3 batch — scaling the S1L1 closure pipeline to the 9 `russell-ch4-*` V0 parents (#170)

> **CORRECTION (2026-06-23):** the Sun-only continuous-TCM proxy numbers below
> (163.6 m/s #188, 2040.6 m/s #192) were computed under a since-fixed 63 s UTC/TDB
> epoch bug (#198, commit 439d279). Corrected: **114.4 / 2020.7 m/s**. The PARTIAL
> verdict is **unchanged** — it rests on the sourced published App-C total Δv
> (420 / 1678 m/s), not the proxy. See
> `docs/notes/2026-06-23-appc-s1l1-tcm-epoch-rederivation.md`.

**Date:** 2026-06-08
**Task:** #170 — apply the proven S1L1 closure pipeline (tasks #166/#167/#169,
`russell-ch4-4.991gG2` = first V3) to the 9 remaining V0 `russell-ch4-*` parent
cyclers, the first scalable V3 batch.

**Discipline:** `feedback_orbit_closure_discipline.md` — independent n-body
cross-check mandatory; ALL binding constraints in the residual; verify topology vs
source BEFORE building; mine reproduction data first; same-model golden target
(App-C real-eph breathing v∞); per-row hold (no batch-trust); a clean negative is
success, never tolerance-tune.

---

## The 9 target rows (V0/absent, live loader `(e.raw or {}).get('validation_level')`)

```
russell-ch4-3.64gGg3      russell-ch4-3.66gfF3      russell-ch4-3.77Gh3
russell-ch4-3.78Gg3       russell-ch4-5.30ggF3      russell-ch4-5.66Gfh3
russell-ch4-6.44Gg3       russell-ch4-8.049gGf2     russell-ch4-8.165Gfh-f2
```

## STEP 1 — reachability triage (App-C "DATA NECESSARY TO REPRODUCE" block?)

The S1L1 pipeline depends ENTIRELY on a sourced App-C per-leg real-eph state block.
Russell 2004 Appendix C documents the lowest-Δv DE405 solutions for **77 of the 203**
circular-coplanar parent cyclers; the App-C / Table-5.5 shorthand uses the **real-eph
simple-model v∞** as its prefix, NOT the coplanar v∞E the catalogue `russell-ch4-*`
IDs use. The authoritative bridge from a catalogue row to its App-C block is the
**parent cycler number** — and only two of the nine rows carry one (Russell states it
explicitly; the catalogue notes record it):

| catalogue id | source table | App-C parent | App-C block? | App-C total Δv (real-eph) |
|---|---|---|---|---|
| `russell-ch4-8.049gGf2`   | Table 4.9  | **#188** | YES (line 8643) | 0.436091 km/s (420 m/s, Table 5.5) |
| `russell-ch4-8.165Gfh-f2` | Table 4.9  | **#192** | YES (line 8697) | 1.677496 km/s (1678 m/s, Table 5.5) |
| `russell-ch4-3.64gGg3`    | Table 4.10 | — | NO sourced parent # | n/a |
| `russell-ch4-3.66gfF3`    | Table 4.12 | — | NO sourced parent # | n/a |
| `russell-ch4-3.77Gh3`     | Table 4.10 | — | NO sourced parent # | n/a |
| `russell-ch4-3.78Gg3`     | Table 4.10 | — | NO sourced parent # | n/a |
| `russell-ch4-5.30ggF3`    | Table 4.12 | — | NO sourced parent # | n/a |
| `russell-ch4-5.66Gfh3`    | Table 4.11 | — | NO sourced parent # | n/a |
| `russell-ch4-6.44Gg3`     | Table 4.13 | — | NO sourced parent # | n/a |

**Why the 7 are NOT-REACHABLE (no sourced App-C block):** these rows are catalogued
canonical members of Russell's **circular-coplanar** Tables 4.10/4.11/4.12/4.13. Each
coplanar (v∞,ToF,aphel) group has dozens of descriptor variants; the catalogue picked
one canonical member. Russell's App-C reproduces a *parent-number-indexed* subset
under the DIFFERENT real-eph v∞ shorthand. There is NO published 1:1 mapping from a
coplanar Table-4.x descriptor to a specific App-C parent number except where Russell
prints it (8.049gGf2→#188, 8.165Gfh-f2→#192). Some prefixes are coincidentally close
(e.g. catalogue `3.78Gg3` vs App-C `3.784Gg3` #74), but the two numbers are different
model v∞ definitions (3.784Gg3 #74 has real-eph avg v∞E = 4.79, not 3.78); asserting
the correspondence would be an UNSOURCED inference — exactly the "forcing" the
discipline forbids. Per the brief rule 2, recorded NOT-REACHABLE and SKIPPED; not
forced.

This leaves **2 reachable rows** with a sourced App-C block: `8.049gGf2` (#188) and
`8.165Gfh-f2` (#192). Both have App-C total Δv far above the 120 m/s V3 budget on
their OWN published reproduction data (420 / 1678 m/s over 7 cycles), so the *gate
expectation* is PARTIAL (the App-C block already records they are not ballistic in the
real ephemeris) — but the pipeline is run per-row to measure the in-band encounter
behaviour and continuous TCM honestly rather than pre-judging.

---

## STEP 2+ — per-row results

### Reachable rows (2) — INDEPENDENT n-body cross-check + continuous TCM

Both reconstruct cleanly: all 7 Mars encounters per row intercept the TRUE DE440 Mars
at the published per-leg v∞ (true-longitude rendezvous Δlon < 0.001°), measured on an
INDEPENDENT REBOUND/IAS15 Sun-only integrator (Russell's patched-conic cruise model;
the 3-Mars-SOI band ≈0.0116 AU, the SAME band #165/#167 used, NOT loosened). The
mid-leg App-C Δv is applied (these parents are POWERED, not ballistic). The
continuous-from-one-seed horizon TCM (#169 method, Sun-only patched-conic V3 tier)
is the V3 budget figure.

**`russell-ch4-8.049gGf2` (App-C #188) — VERDICT: PARTIAL**
- 7 Mars encounters in-band: miss 3.2e-6 … 1.1e-5 AU (≤ 1,691 km, ≪ 1 SOI); v∞
  matched to ≤ 1e-3 km/s vs published (8.42, 8.52, 10.07, 10.77, 11.24, 11.61,
  11.71 — breathes, real-eph).
- Continuous-from-one-seed TCM = **163.6 m/s** over 7 cycles (23.4 m/s/cycle).
  Published App-C total Δv = 420 m/s (Table 5.5).
- **TCM 163.6 m/s > 120 m/s V3 budget → NOT V3.** Encounter geometry is real; the
  fuel cost is not V3-grade. NO writeback.

**`russell-ch4-8.165Gfh-f2` (App-C #192) — VERDICT: PARTIAL (DRIFT-leaning)**
- 7 Mars encounters in-band: miss 3.2e-6 … 1.9e-5 AU (≤ 2,906 km); v∞ matched to
  ≤ 1e-3 km/s vs published (6.79, 8.36, 9.56, 10.61, 10.87, 11.67, 11.86).
- Continuous-from-one-seed TCM = **2040.6 m/s** over 7 cycles (291.5 m/s/cycle).
  Published App-C total Δv = 1678 m/s (Table 5.5) — one of the worst-performing
  real-eph cyclers in the catalogue (per the catalogue note).
- **TCM 2041 m/s ≫ 120 m/s V3 budget → NOT V3.** Heavily powered. NO writeback.

| id | verdict | measured v∞ vs App-C | continuous TCM | commit |
|---|---|---|---|---|
| `russell-ch4-8.049gGf2`   | PARTIAL | matches (4 dp), breathes 8.4–11.7 | 163.6 m/s (>120 budget) | (this commit) |
| `russell-ch4-8.165Gfh-f2` | PARTIAL | matches (4 dp), breathes 6.8–11.9 | 2040.6 m/s (≫120 budget) | (this commit) |

Gate: `tests/nbody/test_appc_batch_nbody.py` (2 @slow, parametrized — encounters
in-band + TCM-over-budget pinned), `tests/search/test_appc_corrected.py` (7 fast —
construction mechanics). Both green. NO catalogue writeback for either — the
encounters are real, the maintenance is not V3-grade. A clean, honest negative.

### Why PARTIAL is the right call (discipline check)
The danger signal "it closed!" does NOT apply here — the encounter half closed, but
the OMITTED binding constraint is the maintenance budget, which both rows fail on
their OWN published App-C Δv. S1L1 was a 0.000000-km/s ballistic parent (62 m/s
continuous TCM); these two are 420 / 1678 m/s powered parents. The S1L1 pipeline
scales mechanically (same reconstruction recipe, same independent gate, same band),
but the V3 *verdict* is parent-specific: only a near-ballistic App-C parent clears
the 120 m/s budget. No tolerance was tuned; the band and budget are identical to S1L1.

### NOT-REACHABLE (7, no sourced App-C reproduction block) — recorded, skipped
- `russell-ch4-3.64gGg3`  (Table 4.10 coplanar; no sourced App-C parent #)
- `russell-ch4-3.66gfF3`  (Table 4.12 coplanar; no sourced App-C parent #)
- `russell-ch4-3.77Gh3`   (Table 4.10 coplanar; no sourced App-C parent #)
- `russell-ch4-3.78Gg3`   (Table 4.10 coplanar; no sourced App-C parent #)
- `russell-ch4-5.30ggF3`  (Table 4.12 coplanar; no sourced App-C parent #)
- `russell-ch4-5.66Gfh3`  (Table 4.11 coplanar; no sourced App-C parent #)
- `russell-ch4-6.44Gg3`   (Table 4.13 near-ballistic coplanar; no sourced App-C parent #)

---

## Rows done / rows remaining

- Triage complete: 9/9 classified (7 NOT-REACHABLE, 2 reachable).
- Reachable rows processed: `8.049gGf2` (#188) PARTIAL, `8.165Gfh-f2` (#192) PARTIAL.
- **Promoted to V3: 0 of 9.** The 2 reachable rows are powered cyclers whose
  maintenance exceeds the V3 budget; the 7 others have no sourced App-C block.
- V3 census unchanged at **1** (S1L1 `russell-ch4-4.991gG2` remains the only V3).
- Rows remaining: none. Batch complete.

## Summary

| outcome | count | rows |
|---|---|---|
| CONFIRMED → V3 | 0 | — |
| PARTIAL (in-band, TCM over budget) | 2 | `8.049gGf2`, `8.165Gfh-f2` |
| NOT-REACHABLE (no sourced App-C block) | 7 | `3.64gGg3`, `3.66gfF3`, `3.77Gh3`, `3.78Gg3`, `5.30ggF3`, `5.66Gfh3`, `6.44Gg3` |

The S1L1 pipeline scaled mechanically (generic `appc_corrected.py` reconstructs any
App-C block including powered ones), but no V3 promotions resulted: the only two
reachable parents are powered (420 / 1678 m/s), and the rest lack a sourced
reproduction block. Promoting any would have required either loosening the budget or
asserting an unsourced coplanar→App-C parent mapping — both forbidden.
