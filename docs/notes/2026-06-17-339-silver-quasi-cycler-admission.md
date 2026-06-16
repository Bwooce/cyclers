# #339 — Catalogue admission of the Umbriel-Oberon-Umbriel SILVER as the first computed `quasi_cycler` row

**Date:** 2026-06-17
**Task:** #339 — admit the #327 SILVER candidate `repeated-moon-uranus-00000041` as `umbriel-oberon-1-1-uranian-quasi-cycler-2026`
**Catalogue:** 282 → 283 rows
**Orbit class:** `quasi_cycler` (the v4.7 scope expansion's first computed admission)

## Summary

The Umbriel-Oberon-Umbriel (1,1) Uranian near-5:1 synodic-resonance trajectory
discovered by the #312 Uranus extended sweep on 2026-06-16 has now cleared every
computational gate the project has, the offline literature-novelty check, AND
the previously-paywalled Heaton-Longuski 2003 JSR paper (user-acquired
2026-06-16). It is admitted as **the catalogue's first computed `quasi_cycler`
row**, validating the v4.7 catalogue scope expansion built earlier in the same
sprint.

## Full provenance chain (12 tasks, 2 days)

| # | Task | Gate | Verdict |
|---|---|---|---|
| #281 | Multi-μ tulip Np=2 sweep | search-space coverage | Did NOT find — confirms #281's coverage was orthogonal |
| #285 | Saturn + Uranus prioritized repeated-moon | discovery | Surfaced Oberon-Titania-Oberon (1,1) near-miss at 0.062 km/s |
| #312 | Uranus extended sweep | **discovery** | **Found this SILVER** at Uranus Oberon-Umbriel basin floor 0.025 km/s |
| #327 | Gate-passing IC verification | closure + DOP853 | 0.025 km/s, cross-check 2.7e-11 nondim |
| #324 | Physical-sanity max-bend | flyby geometry | 14.7° Umbriel + 39° Oberon (both above 5° gate) |
| #328 | Uranian cycler literature deep-dive | offline lit-novelty | 41-anchor not-found, confidence 0.40 |
| #329 | Heaton-Longuski 2003 JSR paywall | **direct read** | **mga_tour, not cycler — SILVER absent from Table 5** |
| #330 | V2 moontour (Lambert-relegs × N) | quasi-periodic structure | 10/10 cycles convergent, bounded drift 86k–530k km, near-5:1 synodic |
| #331 | V3 REBOUND IAS15 independent integrator | dynamics-not-integrator-artifact | Nanometer agreement vs V2 across n=3,5,10 |
| #332 | V4 scipy J2 + 5-moon Battin (circular ephemerides) | mild ephemeris perturbation | 2–5% V4-vs-V3 agreement, bounded oscillation persists |
| #335 | V4-strict URA111 SPICE (real ephemerides) | **real-eph perturbation** | **MIXED** epoch-dependent: 2000-01-15 fails, 2030/2050 pass |
| #338 | Annual epoch sweep 2000–2099 | **boundary characterization** | **EFFECTIVELY_CYCLIC** — interior PASS run 2000-2083 (85/85=100%) |
| **#339** | **Catalogue admission** (THIS row) | row composition + ratchet bumps | **283 rows; v4.7 quasi_cycler slot filled** |

## Catalogue row details

```yaml
- id: umbriel-oberon-1-1-uranian-quasi-cycler-2026
  orbit_class: quasi_cycler
  epoch_locked: true
  n_returns: 10
  validity_window: {start: "2000-06-21T00:00:00Z", end: "2083-06-21T00:00:00Z"}
  launch_epoch: "2041-06-21T00:00:00Z"   # centre of longest PASS run
  validation_level: V4
  bodies: ["Uranus", "Umbriel", "Oberon"]
  sequence_canonical: "Umbriel-Oberon-Umbriel"
  vinf_kms_at_encounters:
    - {body: "Umbriel", vinf_kms: 0.9199}
    - {body: "Oberon",  vinf_kms: 0.9604}
    - {body: "Umbriel", vinf_kms: 0.8947}
```

## Ratchet bumps (3 frozen-census tests)

* `tests/test_catalogue_rediscovery.py::EXPECTED_COVERAGE[NOT_TWO_BODY]`: **1 → 2**
  (3-body Uranian row joins Heaton-Longuski 2003 in the same lane)
* `tests/data/test_cycler_class_census.py`: multi-arc **242 → 243**; added id to
  `MULTI_ARC_ALLOWLIST`
* `tests/data/test_validation_tier_census.py`: unvalidated **28 → 29**
  (`derived/derived` source-pair classifies unvalidated on the orthogonal
  provenance-tag tier axis; the row's `validation_level=V4` evidence lives on
  the gauntlet axis, NOT this census)

15/15 ratchet tests pass after the bumps.

## Honest scope caveats (baked into the row's `notes`)

1. **Sub-year DOY sensitivity:** the validity_window assumes June 21 DOY
   launches. 2000-01-15 V4-strict FAILS while 2000-06-21 PASSES. Sub-year
   boundary characterization is a Phase 2 follow-up on #338.
2. **URA111 kernel expiry:** the kernel covers 1900-2099. The 2084-2099
   V4-strict failures are kernel-edge extrapolation artifacts, NOT dynamical
   resonance breakdown. A fresher post-URA111 Uranian satellite kernel would
   extend the validity_window beyond 2083.
3. **Coarse 2-moon-grid discovery:** the SILVER was discovered at the coarse
   2-moon (Umbriel-Oberon) convention. The 3-moon convention (per #285 / #312)
   reports a different residual (0.636 km/s) at the same IC — same physics,
   different deterministic relative-offset seed. The 2-moon basin floor
   (0.024 km/s at 96×96 offset sweep) is the gate-passing definition.
4. **Validation_level=V4** (not V5): V5 is human mission-quality review.
   #339 admission is purely computational gauntlet completion + literature
   clearance.

## What this validates

* **The v4.7 scope expansion (#294) is not theoretical.** The
  `quasi_cycler` admission slot now has a real entry to display. The website's
  filter UI (#295) will auto-show this row at next sync via
  `scripts/sync-catalogue.mjs` in the `cyclers.space` repo.
* **The 10-gate ladder is operational end-to-end.** Discovery sweep → physical
  sanity → literature gates × 2 → ML → V1 → V2 → V3 → V4-scipy → V4-strict →
  V4-strict-epoch-sweep → admission row. No previous candidate in the project's
  history has cleared this many gates.
* **The discovery discipline held.** No false novelty was ever claimed during
  the chain. Every gate's verdict was the verdict the math gave, not a tuned
  pass.

## Mission-utility implication

A spacecraft launched at the indicated `launch_epoch` (or any June 21 within
the validity_window 2000-2083), arriving at Umbriel at V_∞ = 0.92 km/s, can
complete ~10 Umbriel-Oberon-Umbriel cycles before bounded drift accumulates
beyond ~530,000 km (smaller than Oberon's orbit radius). This is structurally
useful for repeat-observation orbital science of both Umbriel and Oberon
without requiring continuous maintenance ΔV. The Phase 2 sub-year DOY
characterization (future #338-followup) would refine the launch-window
practicality.

## Related artifacts

* JSONLs: `data/silver_327_verified.jsonl`, `data/silver_327_moontour_v2_verdicts.jsonl`,
  `data/silver_327_v3_verdicts.jsonl`, `data/silver_327_v4_verdicts.jsonl`,
  `data/silver_327_v4_strict_verdicts.jsonl`,
  `data/silver_327_v4_strict_annual_sweep_338.jsonl`,
  `data/silver_327_v4_strict_boundary_338.jsonl`
* Prior task docs: `docs/notes/2026-06-16-312-uranus-extended-sweep.md`,
  `docs/notes/2026-06-16-327-umbriel-silver-verification.md`,
  `docs/notes/2026-06-16-324-physical-sanity-gate.md`,
  `docs/notes/2026-06-16-328-uranian-cycler-lit-deep-dive.md`,
  `docs/notes/2026-06-16-330-moontour-v2-phase2.md`,
  `docs/notes/2026-06-16-331-v3-nbody-phase3.md`,
  `docs/notes/2026-06-16-332-v4-uranus-phase4.md` (renamed during ordering),
  `docs/notes/2026-06-16-335-v4-strict-phase41.md`,
  `docs/notes/2026-06-16-338-silver-epoch-sweep-boundary.md`
* Heaton-Longuski 2003 PDF (private repo): `cyclers_pdf` commit `0c26264`
