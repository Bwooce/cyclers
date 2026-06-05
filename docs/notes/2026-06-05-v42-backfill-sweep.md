# Schema v4.2 backfill sweep (Task #98, RESEARCH-ONLY)

Date: 2026-06-05. Hunting sourced backfill for three v4.2 fields:
1. `trajectory.segments[].center` — a segment arc explicitly centred on a body
   other than the trajectory-level center.
2. `trajectory.segments[].tof_days_bounds` — a published TOF *range* (min–max).
3. top-level `source_ephemeris` — paper explicitly names the ephemeris model its
   numbers were computed against.

STRICT PROVENANCE: a finding counts only if explicitly stated in the source AND
attributable to the specific values the catalogue carries from that source.

Known existing backfill (skipped): Aldrin out-em `tof_days_bounds: [161,172]`
(catalogue line 78); `source_ephemeris: "STOUR analytic ephemeris"` on the 2
Aldrin establishment rows (lines 2639, 2758).

**NOTHING was applied. NOTHING committed.** Catalogue / schema / tests untouched.

---

## Counts

| Field | New defensible findings | Notes |
|---|---|---|
| `source_ephemeris` | **4** rows (1 clean + 3 qualified/forward), plus 1 large forward pointer (77 Russell Appendix-C cyclers → DE405, none currently catalogued) | strongest result |
| `segments[].center` | **0** clean. 1 soft candidate (Genova lunar segment) but no explicit per-segment-center value statement | |
| `tof_days_bounds` | **0** new. Only the already-applied Aldrin [161,172] exists in any source | |

Top finding: **Russell 2004 names JPL's DE405** as the ephemeris for its
*accurate-model* (Appendix C) cyclers — but those 77 cyclers are NOT currently
catalogued; the ~201 catalogued Russell rows are all the circular-coplanar
idealized model (Table 3.4 / 4.9 / 5.5) and therefore have NO ephemeris. This is
the cleanest *forward pointer* in the sweep, and an important negative for the
bulk Russell rows.

---

## Stage 1 — existing mining notes

### `2026-06-05-ccsds-odm-502-mining.md`
- Lines 181-194 (CCSDS 502.0-B-3 §6.2.5.15 + §6.2.5.4(d)): CENTER_NAME is a
  per-TRAJ-block keyword; an OCM can carry consecutive blocks centered on SUN
  then EARTH then MARS. **Spec/definitional support for the `segments[].center`
  concept, NOT a cycler value.** No backfill.
- Lines 254-256: CCSDS records "the ephemeris/EOP model the numbers were computed
  against" — external validation of the `source_ephemeris` field intent. No value.

### `multi-arc-classification.md`
- Lines 280-282 (Hollister & Menning 1970): a "direct return orbit" is "a
  sun-centered elliptical orbit" — confirms the 15 EV rows are Sun-centered
  single ellipses. No alternate-center segment; no TOF range. No backfill.

### `2026-06-05-jones-aas17-577-vem-mining.md`
- Lines 192-199: model is "patched-conic on real planetary ephemeris" / "true
  ephemeris" — generic, NOT a named DE/SPICE kernel (see Stage 2 Jones).
- Lines 352, 364: per-row transit-leg ToFs are **point** values (309 d / 259 d;
  268 d / 223 d), not ranges. Architecture-wide "vary from 219 to 309 days"
  (line 185) is a fleet spread across 7 missions, NOT a per-row TOF bound — does
  not qualify as a row `tof_days_bounds`.

### `2026-06-04-agrawal-landau-howe-mining.md`
- All Howe/Agrawal/Landau ToFs are point values or design constraints (e.g. Howe
  up-escalator 151 d; down-escalator 170/191 d internal inconsistency). No TOF
  range, no ephemeris naming, no alt-center. No backfill.

### `s1l1-target-topology-mining.md`, `2026-06-05-endgame-tisserand-mining.md`,
`2026-06-05-vasile-*` — no qualifying center / TOF-range / named-ephemeris
statements tied to catalogue values. (Vasile rows are not in the catalogue.)

---

## Stage 2 — method sections of held PDFs

(Papers cited by author/title/report number only.)

### 1. Sanchez Net et al. 2022, "Cycler Orbits and Solar System Pony Express,"
*J. Spacecraft & Rockets* (AIAA DOI 10.2514/1.A35091). Catalogue rows: 2
(Fig. 2a near-ballistic EEM; Fig. 2b EM/EEM multi-arc).
- §III "Cycler Orbit Selection" / §III.A "Enumeration of Orbits" (p.862–863),
  verbatim: *"A large number of potentially ballistic cycler trajectories can be
  computed in a patched-conic, two-body problem using the **Star** software
  available at Jet Propulsion Laboratory (JPL)."* And p.863: *"Star is a
  broad-search tool that produces a database of trajectories in the
  patched-conics model."*
- The realistic-model optimisation (*"Cosmic optimal control software within the
  Monte suite"*, ref [12]) is described as a **future** step: p.863 / p.1180
  *"the nearly ballistic patched-conic cycler orbits have to be optimized in the
  full model."* It is NOT applied to the Fig. 2 numbers.
- **FINDING (qualified):** the Fig. 2 values our 2 rows carry are
  patched-conic / two-body from the **Star** broad-search tool. The paper names
  **no DE kernel or SPICE** for these numbers. → If the schema accepts a
  tool/model descriptor, `source_ephemeris: "Star (JPL) patched-conic two-body
  broad search"` is defensible. It is NOT a DE-named ephemeris. Recommend leaving
  null unless the schema explicitly wants the model string, because the Fig. 2
  numbers are not computed against a planetary ephemeris in the DE sense.
- No alternate-center segment; no published TOF range. (Fig. 2 panels DO carry
  per-encounter dates = epochs, but that is not one of the three target fields.)

### 2. Jones, Hernandez & Jesick 2017, "Low Excess Speed Triple Cyclers of
Venus, Earth, and Mars," AAS 17-577. Catalogue rows: the VEM family rows
(jones-2017-vem-triple-family etc.).
- Methodology (p.5), verbatim: *"A zero-sphere-of-influence patched conic gravity
  model is used with the **real planetary ephemeris**, and Lambert's problem is
  solved…"* Optimisation (p.7) "in the true ephemeris" adds Sun + all planets +
  Earth's moon.
- **Searched the entire text + references for a named kernel: NO DE405/DE421/
  DE430/DE440, NO SPICE, NO Standish, NO Horizons.** Only the generic phrases
  "real planetary ephemeris" / "true ephemeris" appear.
- **FINDING (qualified):** `source_ephemeris: "real planetary ephemeris
  (unspecified)"` is the most that is sourceable for Jones-derived rows. Explicit
  that numbers are real-ephemeris (not circular-coplanar), but version unnamed.
  Recommend recording only if the schema accepts a non-versioned descriptor.
- No alternate-center segment statement; per-row ToFs are point values (no range).

### 3. Rogers, Hughes, Longuski & Aldrin 2012, "Preliminary Analysis of
Establishing Cycler Trajectories…via V∞ Leveraging," AIAA 2012-4746.
- p.3 ("STOUR Analysis"), verbatim: *"…verified by using the Satellite Tour
  Design Program (STOUR) developed by the Jet Propulsion Laboratory for the
  Galileo mission tour design. **STOUR uses an analytical ephemeris for the
  location of the planets and the patched-conic method**…"*
- Table 4 caption (p.6), verbatim: *"Table 4. STOUR results for the selected
  cyclers **in the analytic ephemeris**."* (Table 3 = genetic-algorithm results
  "in the circular, co-planar model" — NOT ephemeris.)
- Table 4 covers types: Aldrin (4:3(2)−, 3:2(1)−), VISIT-1, VISIT-2, Case 1,
  Case 2, Case 3, S1L1, U0L1. Columns: K:L(M), launch date (epoch), TOF (days),
  periapse/apoapse AU, V∞,launch, ΔV_DSM, V∞,flyby, h_flyby.
- **FINDINGS:**
  - Already applied: 2 Aldrin establishment rows (lines 2639, 2758).
  - **NEW candidate:** row `mcconaghy-2005-em-case1` — its segment `out-em`
    `tof_days: 365` is explicitly "Rogers 2012 Table 4 (STOUR analytic-ephemeris)
    …365 days for the 4:3(2)- Case 1 trajectory" (catalogue line 9644). A
    Table-4 STOUR value lives on this row → `source_ephemeris: "STOUR analytic
    ephemeris"` is defensible. CAVEAT: the row's a/e are Rogers Table 1
    (coplanar); only the segment TOF is Table 4. Mixed-provenance — apply only if
    the field semantics allow "the ephemeris the row's ephemeris-derived numbers
    were computed against."
  - **NEGATIVE:** row `russell-ocampo-2.1.1+2-case2` — text mentions Table 4
    ("171-day…analytic ephemeris", line 804) but the row's catalogued values
    (transit_times [207], AR, TR, turn angles) are all **Russell Table 3.4
    circular-coplanar**, NOT Table 4. Do NOT set source_ephemeris.
  - **NEGATIVE:** row `niehoff-visit1` — Table 4 V∞,flyby=2.834 / launch=2.540
    are quoted in a note (line 1527) but explicitly NOT used as the row's values
    ("those are establishment quantities, not the cycler's own steady V∞").
    Do NOT set source_ephemeris.
- No alternate-center segment. The Aldrin out-em TOF bounds [161,172] (already
  applied) is the only published TOF range in any source this sweep touched
  (161 for 4:3(2)−, 172 for 3:2(1)−; confirmed against Table 4).

### 4. Russell 2004 PhD dissertation (UT Austin). Catalogue rows: ~201
(Table 3.4 / 4.9 / 5.5 idealized members).
- §1.1 (p.2–3), verbatim: planets' positions found by one of two methods; method
  two "**locates the planets using JPL's DE405 Ephemerides. The gravitational
  parameters of the sun and planets are taken from DE405**…" The *idealized*
  catalogue (the "203 noteworthy idealized Earth-Mars free-return cyclers") is the
  **circular-coplanar** simple model (§1.1: "the simplest possible…the
  circular-coplanar solar system").
- Catalogue cross-check: all Russell rows cite Table 3.4 / 4.9 / 5.5 (circular-
  coplanar). "Appendix C" appears only 3× and DE405 appears **0×** in
  catalogue.yaml. Catalogue note (line ~40632) states Appendix C (DE405 ephemeris
  data for the 77 selected cyclers) is **not individually catalogued**.
- **FINDING (forward pointer, NOT a backfill of any current row):** DE405 is the
  named ephemeris **only** for Russell's Appendix-C accurate-model cyclers. The
  ~201 catalogued Russell rows are circular-coplanar → `source_ephemeris` stays
  **null / "circular-coplanar (none)"** for all of them. If/when any of the 77
  Appendix-C ephemeris cyclers are ingested, their `source_ephemeris: "JPL DE405"`
  is firmly sourced (dissertation §1.1).
- No alternate-center segment; no TOF range (idealized model is epoch-free).

### 5. Genova & Aldrin 2015, "Periodic Earth-Moon Orbit … Cycler," AAS-15.
Catalogue row: `genova-aldrin-2015-em-3petal-cycler` (CR3BP, trajectory-level
`center: "Earth"`, line 7514).
- Method (p.2), verbatim: *"The trajectories were designed using AGI's System's
  Tool Kit (STK) Astrogator module, with a high-fidelity force model including
  gravity fields for the Earth (50X50), Moon (50X50), and Sun (4X0)…
  Runge-Kutta 8th/9th order numerical integrator."* **No DE kernel named.**
- **FINDING (qualified):** `source_ephemeris: "STK Astrogator high-fidelity force
  model (Earth 50x50, Moon 50x50, Sun 4x0)"` is a defensible model descriptor for
  this row; it is NOT a named DE ephemeris.
- **segments[].center — SOFT candidate, NOT a clean finding:** this is an
  Earth-centered cycler whose lunar-encounter segment (`petal-emoon-1`, E→Moon,
  line 7532) is physically a Moon-approach arc. The catalogue trajectory `center`
  is "Earth"; a Moon-centered flyby sub-arc would be the textbook
  `segments[].center: "Moon"` case. BUT the paper presents one continuous
  STK/CR3BP trajectory and makes **no explicit per-segment-center statement** with
  attributable values. Per strict-provenance, this does NOT qualify. Flag for
  human judgement only.
- No published TOF range (cadences are point values: 7/10 d Earth, 26/27.5 d Moon).

### 6. Pascarella et al. 2022, "…Pony Express," AAS-22-015.
- Method (p.8) explicitly describes Step 1 STAR patched-conic → Step 2 two-body
  impulsive (planetocentric flyby TPBVP) → Step 3 medium-fidelity heliocentric
  with planetary perturbations → Step 4 high-fidelity low-thrust "ephemeris
  model." Strong per-segment center semantics (planetocentric flyby arcs vs
  heliocentric interplanetary arcs) and an explicit "ephemeris model."
- **NO catalogue rows reference Pascarella / AAS-22-015** (grep: 0 hits). Nothing
  to backfill. Recorded as corroboration of the `segments[].center` concept only.
  (The PDF's later "ephemeris model" likely names a kernel in a part not read,
  but irrelevant absent rows.)

### 7. Patel 2019 (FIT thesis). Catalogue: S1L1 / Aldrin vehicle-design context.
- p.3 only *quotes Rogers' STOUR* numbers (Aldrin 3.449/6.546; S1L1
  2.492/3.657 km/s; direct ΔV 5.011 / 3.796). No independent ephemeris; values
  already attributed to Rogers in the catalogue. No new finding.

### 8. Hollister & Menning 1970 (JSR 7(10)). Catalogue: 15 EV rows.
- Sun-centered single-ellipse "direct return orbits"; Table 3 gives per-encounter
  dates/V_r/turn angles/R_min. **No named ephemeris** (1970 analytic), **no
  alternate-center segment**, **no TOF min-max range**. No backfill.

---

## READY-TO-APPLY backfill list (apply AFTER schema work completes)

Only one row is a genuinely new, defensible application; the rest are qualified
or forward-only. Apply with human review of field semantics.

| Row id | Field | Proposed value | Quote ref | Confidence |
|---|---|---|---|---|
| `mcconaghy-2005-em-case1` | `source_ephemeris` | `"STOUR analytic ephemeris"` | Rogers 2012 Table 4 caption ("STOUR results…in the analytic ephemeris"); the row's `out-em` `tof_days:365` is the Table-4 STOUR value (catalogue line 9644) | Medium — mixed provenance (a/e are Table 1 coplanar; only segment TOF is Table 4). Apply only if field = "ephemeris the row's ephemeris-derived numbers used." |
| `jones-2017-vem-triple-family` (+ any Jones VEM rows) | `source_ephemeris` | `"real planetary ephemeris (unspecified)"` | Jones 2017 p.5 "patched conic…with the real planetary ephemeris" | Low-medium — version unnamed; record only if schema accepts non-versioned descriptor. |
| 2 Sanchez Net Fig. 2 rows | `source_ephemeris` | `"Star (JPL) patched-conic two-body broad search"` | Sanchez Net 2022 §III.A "patched-conic, two-body problem using the Star software" | Low — model descriptor, not an ephemeris; recommend leaving null. |
| `genova-aldrin-2015-em-3petal-cycler` | `source_ephemeris` | `"STK Astrogator high-fidelity force model (Earth/Moon/Sun gravity)"` | Genova & Aldrin 2015 p.2 | Low — force model, not a named DE kernel. |

### FORWARD POINTER (no current row; do not apply now)
- Russell 2004 Appendix-C ephemeris cyclers (77 of 202): when ingested,
  `source_ephemeris: "JPL DE405"` — sourced from dissertation §1.1 ("locates the
  planets using JPL's DE405 Ephemerides…parameters…taken from DE405"). The ~201
  *currently* catalogued Russell rows are circular-coplanar and must remain
  `source_ephemeris: null`.

### Explicit NEGATIVE results
- `segments[].center`: **no clean sourced finding.** Only soft Genova lunar-arc
  candidate (no explicit per-segment-center value statement).
- `tof_days_bounds`: **no new sourced range** beyond the already-applied Aldrin
  [161,172]. Jones "219–309 d" is a fleet spread, not a per-row bound.
- `russell-ocampo-2.1.1+2-case2`, `niehoff-visit1`: Table-4 mentioned in prose
  but row values are NOT Table 4 → do NOT set source_ephemeris.
- All ~201 idealized Russell rows: circular-coplanar → no source_ephemeris.

---

## NOT CHECKED (honest list)
- Vasile/Summerer/De Pascale 2005 & Vasile-Campagnola 2009 method sections —
  font-broken PDFs, and their rows are not in the catalogue (per task guidance).
- Campagnola & Russell 2009 endgame Parts A/B method sections beyond what the
  existing endgame note already mined (no cycler rows trace there).
- Hollister & Rall 1970 NASA-CR, Genova 2016 PADME, Ozimek 2019 LINX,
  Landau-Longuski 2006/2009, Hiraiwa 2026, CCSDS blue book body — not row-bearing
  cycler sources for these three fields (CCSDS already mined as spec support).
- Pascarella AAS-22-015 pp.10+ (its named ephemeris) — skipped: 0 catalogue rows.
- Full bodies of Sanchez Net / Jones / Rogers / Genova beyond the
  method/numerical-setup sections (per task scope: method sections only).

---

## DECISION (2026-06-05, post-sweep review)

The single ready-to-apply candidate (`mcconaghy-2005-em-case1` →
`source_ephemeris: "STOUR analytic ephemeris"`) is **REJECTED** under the
field's spec §16.7.9 semantics ("the ephemeris model the source paper states
its published numbers were computed against"): the row's headline `a/e` are
Rogers Table 1 **coplanar** values; only the `out-em` segment TOF is Table-4
STOUR. A row-level `source_ephemeris` would mislabel the coplanar elements —
the same reasoning that excluded `aldrin-classic-em-k1-outbound` from the
original v4.2 backfill. Net sweep outcome: **zero new applications**; the
negatives (Russell catalogued rows must stay null; Jones/Sanchez name no
kernel) and the DE405 forward pointer (Russell Appendix-C, if ever ingested)
are the durable yield.
