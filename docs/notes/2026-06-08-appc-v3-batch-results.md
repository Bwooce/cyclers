# App-C V3 batch — scaling the S1L1 closure pipeline to the 9 `russell-ch4-*` V0 parents (#170)

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

(filled in below as each reachable row is processed)

| id | verdict | measured v∞ vs App-C | continuous TCM | commit |
|---|---|---|---|---|

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
- Reachable rows to process: `8.049gGf2` (#188), `8.165Gfh-f2` (#192).
