# The §14 V3 class-split (V3-ballistic / V3-powered) — #175

**Date:** 2026-06-08
**Task:** #175 — add a V3-POWERED tier to the §14 validation ladder (parallel to
the V2 class-split) and apply it honestly to the two powered parents the #170
batch confirmed (`russell-ch4-8.049gGf2` #188, `russell-ch4-8.165Gfh-f2` #192).

**Discipline:** `feedback_orbit_closure_discipline` + `feedback_golden_tests_sourced_only`
— the encounters-confirmed side is sourced-anchor evidence (#170's independent
DE440 reconstruction at the published v∞), the documented ΔV is the sourced
budget, and a clean negative (#192 over budget) is the correct honest outcome.
No tolerance was tuned to pass it.

---

## Why split V3

V3 is split for the same reason V2 was (spec §14, 2026-06-07): **a single ΔV bar
cannot judge a ballistic and a powered cycler honestly.** The original V3 holds a
cycler to the **~120 m/s ballistic-maintenance budget** — the right instrument
for a cycler that is *meant to be near-ballistic*. A cycler that is **powered by
design** (a documented, non-zero operational ΔV in its source) flies its
maneuvers as part of the nominal trajectory; the honest question for it is not
"is it under 120 m/s?" but **"does flying it cost no more than its authors said it
would?"** The two tiers differ **only in the budget reference**, not in the
closure evidence required.

## The two tiers (verbatim from the spec amendment)

- **V3-ballistic** (the original V3): phase-matched to a real launch window;
  ephemeris-mode horizon TCM over 3–5 laps (~20–30 yr) bounded and within the
  **~120 m/s ballistic-maintenance budget**. Type specimen: S1L1
  (`russell-ch4-4.991gG2`), 62 m/s continuous TCM over 7 cycles (52% of budget).

- **V3-powered** (new): the cycler's encounters are **independently confirmed on
  the real ephemeris** (in-band miss + true-longitude rendezvous + per-leg v∞
  match, on an integrator independent of the finding solver) AND the
  continuous-from-one-seed horizon TCM — which includes executing the cycler's
  documented nominal maneuvers — is **≤ the cycler's documented operational ΔV**
  (the published per-cycle/total maintenance budget), within a stated tolerance.

## The budget criterion (stated explicitly + defended)

A row earns V3-powered iff its measured **continuous-from-one-seed horizon TCM**
(the #169 method — a single seed integrated through the whole multi-cycle
horizon, *executing the documented nominal maneuvers as part of the run*, so the
figure is the true cost of flying the published trajectory, not a per-encounter
re-anchor) is

```
continuous_TCM  ≤  documented_ΔV × (1 + τ),   τ = 0.10  (10%)
```

The 10% slack absorbs the patched-conic-vs-source modelling gap and the choice of
horizon length; it is **not** tuned per row. The comparison is deliberately
like-for-like: continuous TCM is the realistic finite-horizon maintenance cost
(spec §12a), and the documented ΔV is the source's own stated budget for the same
cycler. A row whose continuous TCM **exceeds** its own documented budget fails —
it costs more to maintain than the source claims, so it is not a faithful
realisation of the documented cycler, and promoting it would be over-claiming.

Clearing V3-powered is **not** a stronger claim than V3-ballistic — a powered
cycler is simply held to the right bar for its class. The "credible (V3+)" trust
gate in §14 applies to both.

## Application — the #170 App-C batch (honest verdicts)

Both rows had their 7 Mars encounters independently reconstructed in-band on
REBOUND/IAS15 DE440 at the published per-leg v∞ with true-longitude rendezvous
(#170, `docs/notes/2026-06-08-appc-v3-batch-results.md`). The verdict is decided
by the budget comparison alone.

| id | App-C # | continuous TCM | documented App-C ΔV (Table 5.5) | TCM / budget | verdict |
|---|---|---|---|---|---|
| `russell-ch4-8.049gGf2` | #188 | **163.6 m/s** | **420 m/s** | 0.39× | **V3-powered PASS** |
| `russell-ch4-8.165Gfh-f2` | #192 | **2040.6 m/s** | **1678 m/s** | 1.22× | **FAIL (over budget) — NOT promoted** |

- **#188 PASS.** 163.6 m/s is 0.39× the documented 420 m/s — well under the 1.10×
  bar. Note this is *over* the 120 m/s ballistic budget, so #188 is **not**
  V3-ballistic — it is exactly the case the class-split exists for: a powered
  cycler that flies within its documented budget.
- **#192 FAIL.** 2040.6 m/s is 1.22× the documented 1678 m/s — over the 1.10×
  bar. #192 is one of the worst-performing real-eph cyclers in the catalogue;
  tuning τ to admit it would defeat the criterion. This is the correct honest
  negative.

**Promoted: 1 of 2** (`russell-ch4-8.049gGf2` → V3-powered). #192 stays at the V0
floor.

## Encoding (mirrors V2-powered)

The schema enum is `V0`..`V5`; "V3-powered" is **not** a new enum value. As with
V2-powered (which the catalogue carries as `V2`), the powered class lives in the
`_LEVEL_EVIDENCE` evidence text, not the enum:

- `data/catalogue.yaml` — `russell-ch4-8.049gGf2` row: `validation_level: V3`.
- `src/cyclerfinder/data/validate.py` — `_LEVEL_EVIDENCE[("russell-ch4-8.049gGf2","V3")]`
  records the V3-POWERED class + the full evidence chain (#170 independent DE440
  in-band reconstruction + true-longitude rendezvous + TCM-vs-documented-ΔV).
- `scripts/backfill_validation_level.py` — `_LEVEL_BY_ID["russell-ch4-8.049gGf2"] = "V3"`.
- `tests/data/test_schema_v45_fields.py` — census ratchet: V3 count 1 → 2.

## Census after this change

- V1 = 11, V2 = 1, V3 = **2** (S1L1 V3-ballistic + #188 V3-powered), V4 = 0, V5 = 0.
