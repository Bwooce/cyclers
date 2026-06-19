# #390 — NAIF SPK direct-V∞ extractor (wholesale #345 unblock)

**Date:** 2026-06-19 AET
**Task:** #390 — derive per-encounter V∞ for the #345 classic-mission `mga_tour`
backlog directly from the missions' archived NAIF SPK kernels at the
already-sourced flyby epochs.
**Status:** Extractor BUILT + VALIDATED on Voyager 2; V∞ EXTRACTED for Voyager 1,
Voyager 2, Mariner 10 (the tractable subset). Catalogue stays 302 — **admission
is a parent-reviewed follow-on** (the V∞ is DERIVED by our own code, so it needs
independent review before writeback per `feedback_orbit_closure_discipline`).

## 1. Why this exists

The #345 digest (`docs/notes/2026-06-19-345-voyager-mariner-mission-digests.md`)
established the firm negative pattern: every Voyager / Mariner-10
mission-overview AND dedicated-navigation paper publishes encounter dates +
closest-approach geometry + maneuver ΔV but **NOT** the per-encounter
hyperbolic-excess velocity (V∞) the §14 V0 `mga_tour` standard requires. Even the
best venue (McKinley-Van Allen 1976 Voyager navigation strategy) tabulates TCM
ΔV but no V∞; the true V∞ values live in unobtainable internal JPL documents.

The acquisition route is therefore exhausted. The V∞ can instead be **DERIVED**
from each mission's archived NAIF SPK kernel at the flyby epochs we already have.
The project already validates DE440-derived heliocentric V∞ to <0.2% (Tito 2018,
Aldrin catalogue rows). #390 extends that machinery to **planetocentric** flyby
V∞ using the spacecraft SPK.

## 2. Extractor design

Module: `src/cyclerfinder/verify/mission_spk.py` (validation infrastructure only;
seeds-not-tracks intact). Two pieces:

* **`ensure_mission_spk(filename, base_url)`** — fetches + caches a named mission
  SPK from NAIF into the astropy cache dir, exactly mirroring
  `spice_kernels.ensure_leapseconds_kernel`'s on-demand-fetch pattern (binary
  kernels never committed; network touched only on first call).

* **`vinf_at_flyby(mission_spk, sc_naif_id, flyby_body, epoch_utc)`** —
  1. `furnsh` DE440 (planet barycenter ephemeris) + LSK (UTC→ET) + the mission
     SPK (spacecraft state).
  2. `str2et(epoch_utc)`.
  3. Over a ±3-day window, sample the spacecraft state relative to `flyby_body`
     (`spkezr`, J2000, body-centered) at 25 epochs.
  4. At each epoch compute the vis-viva hyperbolic-excess speed
     `sqrt(max(0, v_rel² − 2·μ/r_rel))`. Near periapsis this is depressed by the
     deep potential; far outside the SOI it converges to the true asymptotic V∞.
  5. Report V∞ = the value at the **outermost** sample (weakest planet pull),
     plus the outer-window mean/std as convergence evidence, plus the
     closest-approach radius (cross-check vs the published flyby geometry).

Body GM + radius come from `cyclerfinder.core.constants.PLANETS` (sourced JPL
DE440 / IAU 2015 values) — **never hard-coded**. Every record carries SPK-kernel
filename + NAIF body ID + epoch provenance.

### SPK sources (NAIF public archive, fetched on demand; not committed)

| Mission | SPK file | NAIF dir | sc body ID | coverage |
|---|---|---|---|---|
| Voyager 1 | `vgr1_jup230.bsp` | `VOYAGER/kernels/spk/` | −31 | 1979-01-31 → 03-18 |
| Voyager 1 | `vgr1_sat337.bsp` | `VOYAGER/kernels/spk/` | −31 | 1980-08-06 → 11-19 |
| Voyager 2 | `vgr2_jup230.bsp` | `VOYAGER/kernels/spk/` | −32 | 1979-06-04 → 07-20 |
| Voyager 2 | `vgr2_sat337.bsp` | `VOYAGER/kernels/spk/` | −32 | 1981-06-07 → 09-21 |
| Mariner 10 | `M10_archive_1.bsp` | `M10/kernels/spk/` | −76 | **1974-03-21 → 04-04 only** |

## 3. Voyager 2 validation result (P390.3)

Frozen-gate test: `tests/verify/test_390_spk_vinf_extractor.py`. Validation
evidence: `data/390_extractor_validation.jsonl`.

The **tight** self-consistency check is the closest-approach radius (which IS
published, so the EXPECTED side is sourced, not computed by us — non-circular per
`feedback_golden_tests_sourced_only`):

| Flyby | Extracted CA | Published CA | Agreement |
|---|---|---|---|
| V2 Jupiter | 10.09 R_J | 10.0 R_J (Kohlhase-Penzo Table IV, JSX) | 0.9% |
| V2 Saturn | 2.68 R_S | ~2.7 R_S (Voyager 2 Saturn flyby) | <1% |

Vis-viva convergence (V∞ stability across the outer window, target <1%):

| Flyby | V∞ | window std | std/mean |
|---|---|---|---|
| V2 Jupiter | 7.639 km/s | 0.031 km/s | 0.41% |
| V2 Saturn | 10.674 km/s | 0.0007 km/s | 0.007% |

Order-of-magnitude sanity: V2 Jupiter approach V∞ ~7.6 km/s and Saturn ~10.7
km/s are both in the widely-cited bands. **Validation PASSES** — the machinery is
proven; the closest-approach geometry the extractor reads matches the published
mission record to <1%, and the V∞ is vis-viva-stable far outside the SOI.

## 4. Extracted V∞ table (P390.4)

`data/390_mission_vinf.jsonl` — one record per (mission, flyby body, epoch):

| Mission | Flyby | Epoch (CA, UTC) | V∞ (km/s) | std (km/s) | CA radius | CA alt (km) |
|---|---|---|---|---|---|---|
| Voyager 1 | Jupiter | 1979-03-05T12:05 | 10.773 | 0.006 | 4.88 R_J | 277,369 |
| Voyager 1 | Saturn  | 1980-11-12T23:46 | 15.167 | 0.030 | 3.06 R_S | 123,875 |
| Voyager 2 | Jupiter | 1979-07-09T22:29 | 7.639  | 0.031 | 10.09 R_J | 650,060 |
| Voyager 2 | Saturn  | 1981-08-26T03:24 | 10.674 | 0.001 | 2.68 R_S | 101,050 |
| Voyager 2 | Uranus  | 1986-01-24T17:58 | 14.732 | 0.00001 | 4.19 R_U | 81,573 |
| Voyager 2 | Neptune | 1989-08-25T03:56 | 16.742 | 0.0004 | 1.18 R_N | 4,507 |
| Mariner 10 | Mercury | 1974-03-29T20:47 | 10.376 | 0.038 | 1.29 R_Me | 713 |

The Voyager 2 Uranus + Neptune rows were added in the #390 admission follow-on
(P-A): same `vinf_at_flyby()` with zero new code, centered on the true CA epochs
found by a fine periapsis scan. Cross-checks (both sourced): Uranus CA 4.192 R_U
/ 81,573 km alt vs the published 81,500 km above cloud tops (Wikipedia "Voyager
2"; NASA "35 Years Ago: Voyager 2 Explores Uranus" 50,700 mi); Neptune CA 1.182
R_N / 29,271 km radius-from-center vs the published ~29,240 km (4,950 km above
the north pole, Wikipedia "Voyager 2") — both agree to <0.2%. Voyager 2 is now
the complete E-J-S-U-N Grand Tour V∞ tuple.

All closest-approach radii match the published flyby geometry (V1 Jupiter 4.88 vs
published 4.89 R_J; V1 Saturn 3.06 vs ~3.09 R_S; Mariner-10 Mercury I 713 km alt
vs the ~703 km achieved value). The Voyager-1 Saturn V∞ of 15.2 km/s is higher
than Voyager 2's because Voyager 1 took the high-energy Titan-targeted approach.

CA altitude for the giant planets is large only because it is `R·(R_b) − R_eq`
with `R_eq` the 1-bar equatorial radius; the radius-in-planet-radii column is the
meaningful published cross-check, and it matches.

## 5. Honest blockers / partial coverage

* **Mariner-10 Venus flyby (1974-02-05) NOT covered** — `M10_archive_1.bsp`
  spans only 1974-03-21 → 04-04. It captures Mercury encounter I (1974-03-29)
  but not the Venus gravity-assist nor the Mercury II/III re-encounters. No other
  Mariner-10 SPK is published in the NAIF `M10/` archive. **The Mercury-I V∞ is
  the only Mariner-10 datum #390 can deliver.** Venus/Mercury-II/III remain
  blocked pending a wider-coverage Mariner-10 SPK (none currently on NAIF).
* All Voyager flybys at giant planets are fully covered.
* Voyager 2 Uranus (1986) / Neptune (1989) SPKs (`vgr2.ura182.bsp`,
  `vgr2_nep097.bsp`) were extracted in the #390 admission follow-on (P-A) with
  the same code — see the §4 table. Coverage `vgr2_nep097.bsp` spans 1988-11-12
  → 1989-10-01 (Neptune CA well inside); `vgr2.ura182.bsp` covers the Uranus
  encounter. Voyager 2 is the full four-giant-planet Grand Tour.

## 6. RECOMMENDED catalogue `mga_tour` rows (NO writeback — review-gated)

These are RECOMMENDATIONS only. The V∞ values are DERIVED by our own code from
the SPK; per `feedback_orbit_closure_discipline` the "it closed!" danger signal
applies to derived numbers feeding admission, so admission is a parent-reviewed
follow-on, NOT performed here. `data/catalogue.yaml` stays at 302.

Common fields for all three: `source: literature`, `trajectory_regime:
ballistic`, `model_assumption: real-ephemeris`, `cycler_class: multi-arc` (two+
distinct heliocentric arcs; NOT a cycler), `orbit_class: mga_tour`,
`epoch_locked: true`, `n_returns: 1`, `orbit_source: derived`, `vinf_source:
derived`, `orbit_fidelity: real-ephemeris`, `vinf_fidelity: real-ephemeris`,
`source_ephemeris: "DE440 + NAIF reconstructed spacecraft SPK"`,
`validation_level: V0` (sourced epochs + SPK-derived V∞; the closest-approach
geometry independently reproduces the published value to <1%, but the V∞ itself
is not in the published literature — V0 is the honest floor).

### 6a. Voyager 2 (J-S, the validated reference row)

```yaml
- id: voyager-2-jupiter-saturn-mga-tour
  name: "Voyager 2 Jupiter-Saturn gravity-assist tour (MJS77, JSX)"
  orbit_class: mga_tour
  epoch_locked: true
  n_returns: 1
  bodies: ["E", "J", "S"]
  sequence_canonical: "E-J-S"
  launch_epoch: "1977-08-20T00:00:00Z"   # Kohlhase-Penzo 1977 Table IV, JSX launch
  validity_window: { start: "1977-08-20T00:00:00Z", end: "1981-08-26T03:24:00Z" }
  vinf_kms_at_encounters:
    - { body: "J", vinf_kms: 7.639, note: "SPK-derived planetocentric V_inf at Jupiter CA 1979-07-09T22:29 UTC (vgr2_jup230.bsp, NAIF -32); vis-viva outer-window std 0.031 km/s; CA 10.09 R_J vs published 10.0 R_J (Kohlhase-Penzo Table IV)." }
    - { body: "S", vinf_kms: 10.674, note: "SPK-derived planetocentric V_inf at Saturn CA 1981-08-26T03:24 UTC (vgr2_sat337.bsp, NAIF -32); vis-viva outer-window std 0.001 km/s; CA 2.68 R_S." }
  source_quotes:
    "vinf_kms_at_encounters.0.vinf_kms": "DERIVED (#390): planetocentric V_inf 7.639 km/s from NAIF SPK vgr2_jup230.bsp (Voyager 2, body -32) at Jupiter closest approach 1979-07-09T22:29 UTC, vis-viva sqrt(v_rel^2 - 2*mu_J/r_rel) converged outside Jupiter SOI (std 0.031 km/s). CA radius 10.09 R_J independently reproduces Kohlhase-Penzo 1977 Table IV (10.0 R_J)."
    "vinf_kms_at_encounters.1.vinf_kms": "DERIVED (#390): planetocentric V_inf 10.674 km/s from NAIF SPK vgr2_sat337.bsp (Voyager 2, body -32) at Saturn closest approach 1981-08-26T03:24 UTC (std 0.001 km/s)."
  first_published:
    authors: ["Kohlhase, C. E.", "Penzo, P. A."]
    year: 1977
    title: "Voyager Mission Description"
    venue: "Space Science Reviews 21(2):77-101"
```

### 6b. Voyager 1 (J-S)

```yaml
- id: voyager-1-jupiter-saturn-mga-tour
  name: "Voyager 1 Jupiter-Saturn gravity-assist tour (MJS77, JST)"
  orbit_class: mga_tour
  epoch_locked: true
  n_returns: 1
  bodies: ["E", "J", "S"]
  sequence_canonical: "E-J-S"
  launch_epoch: "1977-09-05T00:00:00Z"   # Voyager 1 launch (Kohlhase-Penzo 1977)
  validity_window: { start: "1977-09-05T00:00:00Z", end: "1980-11-12T23:46:00Z" }
  vinf_kms_at_encounters:
    - { body: "J", vinf_kms: 10.773, note: "SPK-derived V_inf at Jupiter CA 1979-03-05T12:05 UTC (vgr1_jup230.bsp, NAIF -31); std 0.006 km/s; CA 4.88 R_J vs published 4.89 R_J." }
    - { body: "S", vinf_kms: 15.167, note: "SPK-derived V_inf at Saturn CA 1980-11-12T23:46 UTC (vgr1_sat337.bsp, NAIF -31); std 0.030 km/s; CA 3.06 R_S (high-energy Titan-targeted approach)." }
  source_quotes:
    "vinf_kms_at_encounters.0.vinf_kms": "DERIVED (#390): V_inf 10.773 km/s from NAIF vgr1_jup230.bsp (body -31) at Jupiter CA 1979-03-05T12:05 UTC; CA 4.88 R_J reproduces the published 4.89 R_J."
    "vinf_kms_at_encounters.1.vinf_kms": "DERIVED (#390): V_inf 15.167 km/s from NAIF vgr1_sat337.bsp (body -31) at Saturn CA 1980-11-12T23:46 UTC."
  first_published:
    authors: ["Kohlhase, C. E.", "Penzo, P. A."]
    year: 1977
    title: "Voyager Mission Description"
    venue: "Space Science Reviews 21(2):77-101"
```

### 6c. Mariner 10 (PARTIAL — Mercury I only; Venus + Mercury II/III blocked)

⚠ The natural Mariner-10 row is E-V-Me-Me-Me (Venus gravity assist + three
Mercury encounters). The published SPK covers **only Mercury encounter I**, so a
complete `vinf_kms_at_encounters` tuple cannot be built. Recommendation: **DO NOT
admit a full Mariner-10 mga_tour row yet** — file the single Mercury-I V∞ as a
supporting datum and revisit if a wider Mariner-10 SPK surfaces. If a
Mercury-I-only partial row is wanted:

```yaml
# PARTIAL — Mercury encounter I only (Venus + Mercury II/III not in any NAIF SPK)
- id: mariner-10-mercury-encounter-1   # partial; not a complete tour tuple
  vinf_kms_at_encounters:
    - { body: "Me", vinf_kms: 10.376, note: "SPK-derived V_inf at Mercury CA 1974-03-29T20:47 UTC (M10_archive_1.bsp, NAIF -76); std 0.038 km/s; CA altitude 713 km vs ~703 km achieved (Giberson-Cunningham 1975)." }
  source_quotes:
    "vinf_kms_at_encounters.0.vinf_kms": "DERIVED (#390): V_inf 10.376 km/s from NAIF M10_archive_1.bsp (Mariner 10, body -76) at Mercury encounter-I CA 1974-03-29T20:47 UTC; CA altitude 713 km reproduces the achieved ~703 km."
```

## 7. Net impact

#390 converts the #345 backlog from "blocked on unobtainable V∞" to "mechanical
SPK extraction" for every mission whose reconstructed spacecraft SPK is on NAIF.
The Voyager 1/2 giant-planet flybys are fully delivered and validated; Mariner-10
is partially delivered (Mercury I only — Venus and the Mercury re-encounters are
genuinely absent from the public SPK archive). The extractor is reusable for the
remaining #345 missions (Galileo, Cassini, Pioneer, Juno) and for the Voyager
Uranus/Neptune encounters with no further code.

The #361/#384/#387 acquisition tasks for Voyager/Mariner V∞ are SUPERSEDED by
#390. Catalogue admission of the recommended rows is the reviewed follow-on.
