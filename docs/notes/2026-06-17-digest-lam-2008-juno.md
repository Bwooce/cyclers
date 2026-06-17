# Lam-Johannesen-Kowalkowski 2008 Juno deep-read digest (#382 paper 3)

Per #382 task brief 2026-06-17 AET (#345 / #361 admission candidates for `mga_tour`). All 11 pages read; this digest is sourced strictly from the paper itself.

## 1. Header (verbatim)

- **Title:** "Planetary Protection Trajectory Analysis for the Juno Mission"
- **Authors:** Try Lam (Staff Engineer, GN&C, M/S 301-121), Jennie R. Johannesen (Juno Mission Design and Navigation Manager, GN&C, M/S 301-360), Theresa D. Kowalkowski (Senior Engineer, GN&C, M/S 301-150). All AIAA Members; all Jet Propulsion Laboratory, California Institute of Technology, Pasadena, CA 91109.
- **Venue:** AIAA/AAS Astrodynamics Specialist Conference and Exhibit, 18-21 August 2008, Honolulu, Hawaii.
- **Paper ID:** AIAA 2008-7368. DOI form: `10.2514/6.2008-7368`.
- **Page count:** 11 (pp. 1-11).

## 2. What the paper actually is

A **planetary-protection trajectory-analysis paper** focused on quantifying impact probabilities with the Galilean satellites (Io, Europa, Ganymede, Callisto) over a 150-year propagation horizon following the Juno mission's nominal Jupiter orbit. The paper is **NOT primarily a flight-design paper**. The flight-design baseline is referenced as **Kowalkowski-Johannesen-Lam 2008 AIAA-2008-7368** (Ref. 1 p. 11 — same three authors, same conference; the abstract Ref. 1 labels the launch-period-development paper, which is **NOT the same paper as this one despite the same paper number AIAA-2008-7368 in Ref. 1**). 

**Critical typesetting note:** Ref. 1 on p. 11 says "Kowalkowski, T. D., Johannesen, J. R., and Lam, T., 'Launch Period Development of the Juno Mission to Jupiter,' AIAA/AAS Astrodynamics Specialist Conference, AIAA-2008-7368, Aug. 2008." — but the present paper's own ID is also AIAA 2008-7368. This is a **typographical defect in Ref. 1's paper number** (likely should be AIAA-2008-7369 or AIAA-2008-7370). The referenced launch-period paper is a sibling paper, NOT this paper. The flight-design data (per-leg V_∞ table) would be in the launch-period paper, NOT this planetary-protection paper.

The paper has six substantive sections (§I-VI):

1. **§I Introduction** (p. 1) — Juno overview; polar orbit, periapsis altitude ~4500 km, apoapsis 39 R_J = 71,492 km × 39 ≈ 2.79 million km, 11-day period, 0.95°/orbit apsidal rotation.
2. **§II Juno Mission Overview** (pp. 1-3) — Launch August 2011 Atlas V 551 Cape Canaveral; 2 DSMs near aphelion; Earth gravity assist 2013-10-12; Jupiter arrival 2016-08-03 (~5 yr ToF); 78-day capture orbit; PRM 2016-10-19 reduces to 11-day science orbit (actual 10.9725 d). Nominal mission: 33 science orbits + 1 extra. De-orbit on apojove of orbit 33, 2017-10-16. Solar conjunction 2017-10-26.
3. **§III Planetary Protection Requirements** (pp. 3-4) — Europa 1×10⁻⁴ / others 1×10⁻³ impact probability ceiling.
4. **§IV Impact Probabilities with the Galilean Satellites** (pp. 4-8) — Monte-Carlo simulation methodology. Tables 1-6 give Jupiter gravity field, Atreya atmosphere, orbit-lifetime statistics, failure-mode impact probabilities.
5. **§V Method of Reducing Impact Probabilities** (pp. 9-10) — Lower perijove altitude OR change inclination (85° / 95°).
6. **§VI Conclusion** (p. 11) — Reference mission meets requirements.

References 1-9 (p. 11) cover the sibling launch-period paper (Ref. 1), the Planetary Protection memo (Ref. 2), the JPL Juno Planetary Protection Plan (Ref. 3), DE414 / JUP230 / JUP100 ephemerides (Refs. 4-7), and the Juno Navigation Plan (Ref. 9, JPL D-35557, March 2008).

## 3. Per-encounter data extracted

### 3.1 Body sequence (canonical)

E → DSM → DSM → E → J, then 33 Jupiter orbits + de-orbit. The "2+ ΔV-EGA" trajectory archetype (2 DSMs + Earth Gravity Assist) per §II.

### 3.2 Epoch tuple (sourced p. 1-2)

| # | Encounter | Date | Notes | Source |
|---|---|---|---|---|
| 0 | Earth launch | 2011-08 (month only) | Atlas V 551 Cape Canaveral | p. 1, §II |
| 1 | DSM-1 | "near aphelion" (date not given) | First Deep Space Maneuver | p. 1, §II |
| 2 | DSM-2 | "near aphelion" (date not given) | Second Deep Space Maneuver | p. 1, §II |
| 3 | Earth flyby | 2013-10-12 | "Earth on October 12, 2013, a little over 2 years from launch" | p. 1, §II |
| 4 | Jupiter arrival (JOI) | 2016-08-03 | "Juno to arrive at Jupiter on August 3, 2016, roughly 5 years in flight time" | p. 1, §II |
| 5 | PRM (Period Reduction Maneuver) | 2016-10-19 | "Period Reduction Maneuver (PRM) will be performed at the next perijove pass on October 19, 2016" | p. 1, §II |
| 6 | Science PJ-3 (primary science observation start) | 2016-12 (approx, on 2nd 11-day orbit) | "primary science observation begins on the second 11-day orbit (PJ-3)" | p. 2, §II |
| 7 | De-orbit maneuver | 2017-10-16 | "manuever on the last apojove to de-orbit the spacecraft on October 16, 2017" (p. 4) | p. 4, §III |
| 8 | Jupiter atmospheric impact | before 2017-10-26 | "before solar conjunction on October 26, 2017" (p. 4) | p. 4, §III |

### 3.3 V_∞ — NOT PUBLISHED (zero values)

**Critical finding:** the paper does NOT publish V_∞ at any encounter. There is NO Earth-departure V_∞, NO Earth-flyby V_∞, NO Jupiter-arrival V_∞ at JOI. The paper's focus is impact-probability post-JOI, not heliocentric trajectory design.

**Sourced numerical values that ARE published:**

- **Periapsis altitude (Jupiter):** 4500 km above 1-bar pressure level (§I p. 1, §II p. 1).
- **Apoapsis distance:** 39 R_J (§I p. 1). R_J = 71,492 km (equatorial radius taken; explicit). So apoapsis = 2,788,188 km.
- **Capture orbit period:** 78 days (§II p. 1).
- **Science orbit period:** 11 days nominal; actual 10.9725 days (§II p. 1).
- **Number of science orbits:** 33 + 1 extra (§II p. 2; numbered by perijove between apojoves; Fig. 1).
- **Apsidal rotation rate:** 0.95° per orbit (§I p. 1).
- **Latitude of perijove:** 3°N at Jupiter arrival; 34°N at PJ-33 (§II p. 2).
- **Initial ascending-node distance:** ~37 R_J (Fig. 3, after PRM).
- **End-of-mission ascending-node distance:** ~9 R_J (Fig. 3, near Europa).
- **Total mission ΔV (Reference Trajectory):** **1962.6 m/s** (Table 8 p. 10).
- **Atreya atmospheric model:** density at altitudes 0-4000 km (Table 2 p. 6). Surface 0.146 kg/m³ → 4000 km altitude 9.95×10⁻²⁰ kg/m³.
- **Jupiter gravity field (Table 1 p. 5):** J2 = 14696.430 × 10⁻⁶, J3 = -0.640 × 10⁻⁶, J4 = -587.140 × 10⁻⁶, J5 = 0 (not modeled), J6 = 34.250 × 10⁻⁶, C22 = 0.0065 × 10⁻⁶, S22 = -0.0125 × 10⁻⁶. All ± 1-sigma uncertainties given.

### 3.4 Failure-mode trajectory ΔV breakdown — NOT given in this paper

The total mission ΔV of 1962.6 m/s (Table 8) is given for the reference trajectory only. Per-maneuver breakdowns (DSM-1, DSM-2, JOI, PRM, OTMs, de-orbit) are NOT published in this paper. JOI is described qualitatively as the largest single maneuver but no number given.

For the **alternative trajectory cases** (Table 8):
- Low altitude case: 1960.3 m/s total mission ΔV
- 85° initial inclination: 1959.9 m/s total mission ΔV
- 95° initial inclination: 1960.5 m/s total mission ΔV

These are all within ~3 m/s of the reference trajectory — the alternative cases don't meaningfully change the ΔV budget.

## 4. Catalogue admission verdict

### 4.1 V0 admissibility test

Per the Heaton-Longuski 2003 U00-01 admission precedent and the Wolf-Smith 1995 Cassini / Bourke 1971 / Dunne-Burgess 1978 negative-verdict precedents in this digest sweep, the V0 evidence standard for an `mga_tour` row is:

**Per-flyby epoch + V_∞ tuple + body sequence are minimum.**

Lam-Johannesen-Kowalkowski 2008 publishes:

- Per-flyby epoch: PARTIAL. Earth flyby 2013-10-12 (day-precision; no hour); JOI 2016-08-03 (day-precision; no hour); PRM 2016-10-19; de-orbit 2017-10-16; impact before 2017-10-26. Launch is only "August 2011" (month).
- Body sequence: YES (E → E → J, with 2 DSMs and 33 polar science orbits — though the "33 orbits + de-orbit" is intra-Jupiter-orbit, not an mga_tour leg sequence in the catalogue sense).
- V_∞ tuple: **NO**. No V_∞ at Earth flyby, no V_∞ at JOI, no V_∞ at launch. Only the total mission ΔV (1962.6 m/s) is published.

**Verdict: the paper is ONE COLUMN SHORT of the V0 `mga_tour` minimum** (same shortfall as Wolf-Smith 1995, Bourke 1971, Dunne-Burgess 1978 — fourth negative verdict in this digest sweep).

### 4.2 Two layered problems

1. **V_∞ absent.** The paper's scope is planetary-protection (post-JOI), not flight-design (pre-JOI). The flight-design data is in the **sibling paper Kowalkowski-Johannesen-Lam 2008** (Ref. 1 p. 11 — "Launch Period Development of the Juno Mission to Jupiter") which is referenced but not in this PDF and not in cyclers_pdf at commit `fa38aae` for this admission cycle. **The flight-design paper IS the correct source; this PP paper is not.**

2. **Juno orbit_class question.** The catalogue admission asks for `mga_tour`, but Juno is a 1-leg E-E-J trajectory (one Earth gravity assist between two DSMs, then JOI). It is not a multi-body MGA tour in the Cassini / Galileo / Voyager sense. It is more accurately characterized as:
   - **mga_tour with single Earth flyby:** if any pre-JOI gravity assist counts (Juno's single EGA does qualify under a strict reading of "mga_tour" = "any multi-gravity-assist trajectory"). The "2+ ΔV-EGA" trajectory archetype IS a recognized mga subclass (Diehl-Kaplan-Penzo 1983, Sims-Longuski-Staugler 1997 ΔV-EGA papers).
   - **precursor_mga:** if the requirement is multiple flybys per leg.
   - The scope-expanded catalogue (#294 / `project_catalogue_scope_expanded_2026-06-15.md`) admits `mga_tour` rows for single-EGA archetypes like Juno (and Galileo VEEGA, and the Cassini VVEJGA pre-Saturn arc).

### 4.3 Recommendation

**RECOMMEND: do NOT admit Juno as a `mga_tour` row from Lam 2008 alone.** V_∞ data is missing (V0 minimum fails). Instead:

1. **File Lam 2008 as a KNOWN_CORPUS planetary-protection literature anchor** for the Juno mission post-JOI orbit geometry. This is the canonical reference for the 78-day capture / 33-orbit science design and the 0.95°/orbit apsidal rotation.

2. **(Action — acquisition gate)** **HIGH PRIORITY** acquisition: **Kowalkowski-Johannesen-Lam 2008 "Launch Period Development of the Juno Mission to Jupiter"** (Ref. 1 p. 11). This is the sibling paper with the same three authors, presented at the same AIAA/AAS Astrodynamics Specialist Conference, August 2008. **Despite the typo'd paper number** ("AIAA-2008-7368" in Ref. 1 is the same as this paper's number — a clear typesetting defect), this paper exists and should have per-leg V_∞ for the 2+ ΔV-EGA trajectory. The correct paper number is most likely **AIAA-2008-7369** or **AIAA-2008-7367** (adjacent numbers); query the AIAA Aerospace Research Central catalogue. DOI form expected: `10.2514/6.2008-7369` or similar.

3. **If a V_∞ table surfaces** from the launch-period paper, admit a `juno-egadsm-2-plus-2011` row at V0 with the structural template below.

### 4.4 Structural template (for future use, NOT for writeback now)

If/when V_∞ data surfaces from the Kowalkowski 2008 launch-period sibling paper:

```yaml
- id: juno-egadsm-2-plus-2011-jovian-orbiter
  name: "Juno 2+ ΔV-EGA trajectory and Jupiter polar orbiter (Lam-Johannesen-Kowalkowski 2008 reference)"
  source: literature
  trajectory_regime: ballistic-impulsive   # 2 DSMs + EGA + JOI; pre-JOI is impulsive deep-space + ballistic Earth flyby; post-JOI is ballistic with PRM/OTMs
  model_assumption: analytic-ephemeris   # DE414 + JUP230 (Lam 2008 §IV.A)
  cycler_class: multi-arc   # at least: launch-DSM-DSM-EGA-JOI = 1 heliocentric mga arc + 33+ Jupiter-orbit science arcs
  orbit_class: mga_tour   # schema v4.7; single-EGA 2+ΔV-EGA flight archetype (Sims-Longuski-Staugler 1997 family)
  epoch_locked: true   # 2011-08 launch window selected for the 2013-10-12 EGA geometry; not repeatable
  n_returns: 1   # single Jupiter arrival; nominally de-orbited at end of 33 orbits
  validity_window:
    start: "2011-08-05T00:00:00Z"   # Juno actual launch (not given in Lam 2008; standard NASA fact-sheet)
    end:   "2017-10-26T00:00:00Z"   # before solar conjunction per p. 4 §III
  launch_epoch: "2011-08-05T00:00:00Z"
  validation_level: V0
  source_ephemeris: "DE414 + JUP230 (Jacobson 2003) per Lam-Johannesen-Kowalkowski 2008 §IV.A"
  bodies: ["E", "J", "Io", "Europa", "Ganymede", "Callisto"]
  sequence_canonical: "E-DSM-DSM-E-J(+33 orbits)"
  vinf_kms_at_encounters: <NEEDS-MORE-DATA: not published in Lam 2008; expected in Kowalkowski 2008 launch-period paper>
  encounter_epochs:
    - body: "Earth"   # EGA
      epoch: "2013-10-12T00:00:00Z"   # day-precision per Lam 2008 p. 1
      closest_approach_km: <not in Lam 2008>
      source_quote: "gravity assist with Earth on October 12, 2013" (Lam 2008 p. 1 §II)
    - body: "Jupiter"   # JOI
      epoch: "2016-08-03T00:00:00Z"
      closest_approach_km: 4500   # nominal perijove altitude
      source_quote: "Juno to arrive at Jupiter on August 3, 2016" + "closest approach at roughly 4500 km altitude above the Jupiter 1-bar pressure level" (Lam 2008 §I p. 1)
  delta_v_total_kms: 1.9626   # 1962.6 m/s per Table 8 reference trajectory
  jupiter_orbit:
    capture_period_days: 78
    science_period_days_nominal: 11
    science_period_days_actual: 10.9725
    science_orbit_count: 33
    perijove_altitude_km: 4500   # 1-bar altitude
    apojove_distance_RJ: 39
    apojove_distance_km: 2788188   # 39 × 71492
    apsidal_rotation_deg_per_orbit: 0.95
    perijove_latitude_deg_arrival: 3
    perijove_latitude_deg_pj33: 34
    deorbit_epoch: "2017-10-16T00:00:00Z"
  first_published:
    authors: ["Kowalkowski, T. D.", "Johannesen, J. R.", "Lam, T."]
    year: 2008
    title: "Launch Period Development of the Juno Mission to Jupiter"
    venue: "AIAA/AAS Astrodynamics Specialist Conference, Honolulu, HI, August 2008"
    paper_id: "AIAA-2008-7368"   # PER LAM 2008 REF. 1 — but this is the same number as this paper, almost certainly a typo
  corroborating_sources:
    - authors: ["Lam, T.", "Johannesen, J. R.", "Kowalkowski, T. D."]
      year: 2008
      title: "Planetary Protection Trajectory Analysis for the Juno Mission"
      venue: "AIAA/AAS Astrodynamics Specialist Conference, Honolulu, HI, August 2008"
      paper_id: "AIAA-2008-7368"
      doi: "10.2514/6.2008-7368"
```

Lam 2008 would be the `corroborating_source` for Jupiter-orbit geometry; Kowalkowski 2008 would be the `first_published` row for the heliocentric trajectory and V_∞ tuple.

## 5. KNOWN_CORPUS impact

Search the existing literature_check.py for Juno / Jupiter polar orbiter / EGA anchors. If none, adding a Juno anchor is reasonable:

- Name: `"Juno 2+ΔV-EGA trajectory to Jupiter (Lam-Johannesen-Kowalkowski 2008 / Kowalkowski 2008)"`
- body_set: `frozenset({"E", "J"})` — the heliocentric trajectory; the Galilean satellites are AVOIDED targets in this mission, not encountered, so they should NOT be in body_set despite being mentioned in the paper.
- topology_label: `{"mga-tour"}` — single-EGA mga is a valid subclass per the #294 scope expansion.
- Citation: `"Lam, Johannesen & Kowalkowski, 'Planetary Protection Trajectory Analysis for the Juno Mission,' AIAA 2008-7368 (2008) -- Jupiter polar science orbit + planetary-protection impact-probability analysis (33 11-day science orbits, 78-day capture, perijove 4500 km, apojove 39 R_J, apsidal rotation 0.95 deg/orbit, EGA 2013-10-12, JOI 2016-08-03, de-orbit 2017-10-16); Kowalkowski, Johannesen & Lam, 'Launch Period Development of the Juno Mission to Jupiter,' AIAA Aug 2008 (paper number likely AIAA-2008-7369 — Ref. 1 typesetting needs verification)."` DOI: `10.2514/6.2008-7368` for the Lam paper.
- Authors: `("Lam", "Johannesen", "Kowalkowski")`.

This anchor would prevent a future false-novelty claim for "Juno" / "2+ ΔV-EGA" / "Jupiter polar orbiter" searches.

## 6. Errata

Versus the task brief and the paper text:

- **Brief said "extract per-encounter V_∞ + Earth flyby + Jupiter Orbit Insertion (JOI) + DSM data for Juno's '2+ ΔV-EGA' trajectory."** The paper does NOT publish V_∞ at any encounter. DSM dates are not given (only described as "near aphelion"). DSM ΔV magnitudes are not given individually. JOI ΔV is not given individually. **The brief overestimated the paper's content** — this paper is planetary-protection, not flight-design.
- **Brief said "Enable #345 Juno mga_tour catalogue admission."** Admission cannot be enabled from this paper alone. The sibling paper Kowalkowski 2008 is the correct source.
- **Ref. 1 paper-number typo:** Ref. 1 on p. 11 lists "AIAA-2008-7368" for the sibling launch-period paper, which is identical to this paper's own ID. This is a typesetting defect. The correct sibling paper number must be one off (AIAA-2008-7367 / 7369 / 7370). Recommended verification step: search AIAA Aerospace Research Central for "Kowalkowski Johannesen Lam Juno Launch Period" in 2008 conference proceedings to recover the correct paper ID.
- **Atlas V 551 launch vehicle:** The paper says "Atlas V 551" (p. 1 §II); Juno actually launched 2011-08-05 on Atlas V 551 (NASA fact sheet confirms; not in this paper). Date precision in this paper is "August 2011" — month only.
- **Author affiliations are JPL, all three authors at M/S 301-NNN GN&C section** — consistent with JPL Outer Planets Flight Operations / Mission Design.

## 7. Action items

For the parent / #382 / #345:

1. **DO NOT writeback an `mga_tour` row for Juno from Lam 2008 alone.** V_∞ data is missing (V0 minimum fails); DSM and JOI ΔVs not individually published; only the total mission ΔV is given.

2. **(Action — KNOWN_CORPUS amendment)** Open a small follow-up task to add a `Juno 2+ΔV-EGA / Jupiter polar orbiter` literature anchor to `src/cyclerfinder/search/literature_check.py` covering E/J with topology_label `{mga-tour}` and the citation per §5 above. This prevents future false-novelty claims for "Juno" / "Jupiter polar orbiter" / "2+ ΔV-EGA" searches.

3. **(Action — acquisition gate, HIGH PRIORITY)** **Kowalkowski-Johannesen-Lam 2008 "Launch Period Development of the Juno Mission to Jupiter"** — the sibling launch-period paper at the same Aug 2008 AIAA/AAS conference. Resolve the Ref. 1 paper-number typo by querying AIAA ARC; likely AIAA-2008-7367 / 7369 / 7370. This is the V_∞-table-containing flight-design paper.

4. **(Action — typesetting defect documented)** The Ref. 1 typo (same paper number for two different papers) is a clean errata finding for the project's #356 paper-defect-flagger workflow. Catalogue it.

5. **(Action — Juno planetary-protection data)** Lam 2008's Jupiter gravity-field table (Table 1) and Atreya atmospheric model (Table 2) are useful as a sourced Jupiter-system data reference. If the Jupiter-system gravity / atmosphere models are needed for any cyclerfinder N-body validation, this paper is a citable source.

6. **(Action — Sims-Longuski-Staugler 1997 ΔV-EGA family)** Juno's "2+ ΔV-EGA" archetype is documented in Sims-Longuski-Staugler 1997 (the canonical ΔV-EGA paper). If that paper is not yet in cyclers_pdf, it is a natural acquisition to anchor the flight-design archetype.

End of digest paper 3 of 4.
