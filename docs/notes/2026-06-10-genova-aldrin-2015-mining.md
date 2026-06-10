# Mining Note: Genova & Aldrin (2015), AAS 15

**Source:** Anthony L. Genova and Buzz Aldrin, "A Free-Return Earth-Moon Cycler Orbit for an
Interplanetary Cruise Ship," AAS 15-___ (paper number blank in preprint header),
AAS/AIAA Astrodynamics Specialist Conference, Vail CO, August 9–13, 2015.
NTRS accession 20150018049 / NASA report ARC-E-DAA-TN22765.
**No DOI found in paper.** URL: https://ntrs.nasa.gov/citations/20150018049

**Mined:** 2026-06-10 from full PDF (11 pp. + references).

---

## 1. Which orbit the catalogue row corresponds to

The paper presents TWO Earth-Moon cyclers:

| Cycler | Rotating-frame appearance | Lunar resonance | Origin |
|--------|--------------------------|-----------------|--------|
| **3-petal** (primary subject) | 3 lobes in E-M rotating frame | 3:1 | *This paper* — previously undiscovered |
| 4-petal (secondary / transition) | 4 lobes in E-M rotating frame | Arenstorf (1963) | Ref. 14 / Arenstorf 1963 |

The catalogue row `genova-aldrin-2015-em-3petal-cycler` corresponds to the **3-petal cycler**,
the newly discovered orbit that is the central subject of the paper (p.2 heading "3-Petal
Earth-Moon Cycler"; Fig. 3; Fig. 5). The 4-petal Arenstorf cycler appears only as a comparison
and contingency option (p.6 heading "Transition to a 4-Petal Earth-Moon Cycler"). There is no
ambiguity in the row's identity.

---

## 2. Character of the orbit

**Free-return, circumlunar, 3:1 lunar resonant (3 petals in the Earth-Moon rotating frame).**

- The cycler is a periodic circumlunar orbit, **not** a libration-point orbit (p.1 abstract; p.2 body).
- Each cycle consists of: two sub-lunar phasing orbits (holding orbits at perigee) → figure-8
  free-return arc to the lunar farside → return to 3,000 km Earth perigee (Fig. 3, caption p.4).
- Aldrin's 1985 concept "C-2-R" (ref. 22, p.3) theorized this family; the paper shows it exists
  only when **solar gravity is added** to the CR3BP model (p.3):
  > "the addition of solar gravity to the astrodynamics model and a modest ΔV maneuver causes
  > Aldrin's C-2-R theorized cycler's apogee to drop temporarily below lunar distance which
  > enables lunar phasing."
- The cycler **requires station-keeping maneuvers** and is therefore **not a pure CR3BP periodic
  orbit**; it is a solar-perturbed, maneuver-maintained trajectory (p.2, p.3).
- The force model used: AGI STK/Astrogator with Earth 50×50, Moon 50×50, Sun 4×0 gravity
  fields; Runge-Kutta 8th/9th order integrator (p.2).

---

## 3. Encounter cadences and maneuver budget

All values sourced directly from the paper; page/figure citations given.

### Earth encounters
- Every **7 or 10 days** (abstract, p.1; Fig. 4 caption, p.4).
- 7 days between Earth encounters on legs that include the lunar flyby (Fig. 4 caption).
- Just under 10 days between Earth encounters during the holding orbits (Fig. 4 caption).
- Earth perigee altitude: **3,000 km** (p.3; Fig. 3 caption).

### Lunar encounters
- Every **26 days** (abstract, p.1; Fig. 4 caption, p.4).
- Lunar perilune altitude: **3,000 km** on the **lunar farside** (p.3; Fig. 3 caption).
- The cycler flies a free-return, figure-8 "Arenstorf orbit" arc bringing it back to 3,000 km
  perigee at each Earth encounter (p.3).

### Station-keeping ΔV
- Three maneuvers per cycle (one lunar cycle ≈ 26 days):
  - **13 m/s** at perigee (point C in Fig. 3): enables Moon phasing; yields close-Earth encounters
    every 9.5 days in two sub-lunar holding orbits (p.3).
  - **14 m/s** at perigee (point E in Fig. 3): increases apogee back to true lunar distance (p.3).
  - **8 m/s** 0.5 days after E (point ~E+0.5d): corrects out-of-plane drift due to solar
    perturbations (p.3).
  - Total per cycle: ~**35 m/s** (13 + 14 + 8).
- Solar gravity perturbations cause variance: **20 to 62 m/s per cycle** (p.4).
- Abstract states average **39 m/s per month** (p.1 abstract) — this is the average of the
  20–62 m/s range across the full apsidal rotation analysis period of 553 days.

### Full apsidal rotation period
- **553 days** (p.5; Fig. 5 caption): "Trans-lunar injection (TLI) of the cycler assumed on
  July 12, 2019 with the last perigee (following the 22nd lunar flyby) on Jan. 15, 2021, or
  553 days." 21 complete cycles (22 lunar flybys) span one full revolution in the inertial frame.
- This is the natural "super-period" of the 3-petal orbit in the presence of solar perturbations
  (the line of apsides rotates 360° in 553 days).

### Inclination
- Cycler inclination approximately **24 degrees** (p.8: "For launch latitudes much higher than
  the 24 degrees contained in the cycler…").
- In the Moon true-of-date (TOD) frame the inclination is about **174 degrees** (p.8).

---

## 4. Launch / injection epoch used in analysis

- TLI assumed **July 12, 2019** (Fig. 5 caption, p.5).
- First lunar flyby: **July 15, 2019** (Fig. 4 caption, p.4).
- Analysis runs through Jan. 15, 2021 (22nd lunar flyby + final perigee).

These are the *reference trajectory* dates used in the paper's STK/Astrogator run; they are not
defining constraints of the periodic orbit itself. The orbit family repeats with the apsidal
rotation cycle.

---

## 5. Search for CR3BP initial conditions (state_nd, Jacobi constant, period_nd)

**Result: NO usable CR3BP initial state exists in this paper.**

The paper explicitly states that the 3-petal cycler does **not** exist in the pure (unperturbed)
CR3BP (p.3):
> "Although in 1985 Aldrin theorized the existence of a 3:1 resonance free-return Earth-Moon
> cycler designated as C-2-R (Fig. 2, left), such a cycler was not shown to exist in the
> restricted three-body problem: Arenstorf shows a 3:1 resonance orbit but without the required
> close Earth passes (Fig. 2, center)."

The orbit exists **only** when solar gravity (and modest ΔV station-keeping) is added to the
CR3BP model. The paper:
- Never tabulates a Jacobi constant.
- Never tabulates a non-dimensional period.
- Never tabulates a rotating-frame state vector.
- Never tabulates a stability index.
- Uses AGI STK/Astrogator (not a CR3BP corrector) to design the trajectory.

The `orbit_elements.cr3bp.state_nd: null` in the catalogue row is **correct and cannot be
improved from this source.** The `jacobi_constant`, `period_nd`, `stability_index`, `lunit_km`,
and `tunit_s` fields also remain null.

---

## 6. What the paper does provide that is new vs. the previous catalogue entry

The catalogue was previously populated from the NTRS abstract only. The full PDF adds:

| Field | Previously known | Now confirmed / improved |
|-------|-----------------|--------------------------|
| 39 m/s/month average ΔV | From abstract | Confirmed: 39 m/s avg (abstract p.1); detailed 20–62 m/s range per cycle (p.4) |
| 7/10-day Earth cadence | From abstract | Confirmed: 7 d (lunar-flyby legs), <10 d (holding-orbit legs) — Fig. 4 caption |
| 26-day lunar cadence | From abstract | Confirmed: p.1 + Fig. 4 caption |
| 3,000 km perilune / perigee | In catalogue | Confirmed: p.3 + Fig. 3 caption |
| 3:1 resonance, 3-petal shape | In catalogue | Confirmed: p.2–3 |
| mass_ratio μ = 0.01215 | In catalogue | Confirmed: p.1 |
| cycler_class: non-keplerian | In catalogue | Confirmed: solar-perturbed, maneuver-maintained |
| **Apsidal rotation period: 553 days** | NOT in catalogue (null) | NEW from Fig. 5 caption, p.5 |
| **Inclination ~24° (geocentric)** | NOT in catalogue | NEW from p.8 |
| **Inclination ~174° (Moon TOD frame)** | NOT in catalogue | NEW from p.8 |
| **Three maneuver details per cycle (13/14/8 m/s)** | NOT in catalogue | NEW from p.3 |
| **Reference epoch TLI July 12, 2019** | NOT in catalogue | NEW from Fig. 5 caption, p.5 |
| **AAS paper number** | "AAS 15-___" (blank in preprint) | Still blank in preprint — the AAS number was not assigned at preprint submission |
| CR3BP state_nd | null — correctly null | Still null (orbit does not exist in pure CR3BP) |

---

## 7. Proposed update to `genova-aldrin-2015-em-3petal-cycler`

**NO writeback to data/catalogue.yaml per task rules.**

The following is a sourced proposed update for review:

### Fields to add/update (all sourced; page/figure cited):

```yaml
# --- period sub-block: add apsidal_rotation_days ---
period:
  pair: "E-Moon"
  k: null
  years: null
  apsidal_rotation_days: 553   # SOURCED: Fig. 5 caption p.5 — "553 days" for one full apsidal revolution;
                               # 21 complete cycles (22 lunar flybys) starting TLI July 12, 2019.
  note: |
    Apsidal rotation period 553 days (Fig. 5 caption, p.5). Earth encounters every 7 days
    (legs containing lunar flyby) or just under 10 days (holding-orbit legs) (Fig. 4 caption,
    p.4). Lunar encounter every 26 days (abstract, p.1; Fig. 4 caption, p.4).

# --- orbit_elements: add inclination_deg ---
orbit_elements:
  inclination_deg: 24   # SOURCED: p.8 — "24 degrees contained in the cycler" (geocentric frame)
  inclination_moon_tod_deg: 174   # SOURCED: p.8 — "cycler's inclination is about 174 degrees in
                                  # the Moon true-of-date (TOD) frame"
  # all cr3bp sub-fields remain null — not a pure CR3BP orbit (see data_gaps)

# --- trajectory.maneuvers: add three per-cycle maneuver details ---
maneuvers:
  - label: "MAN-C"
    at: "perigee (point C, Fig. 3)"
    dv_ms: 13   # SOURCED: p.3 — "first maneuver (13 m/s of ΔV) is performed at perigee (Fig. 3, C)
                #  which enables phasing with the Moon and yields close-Earth encounters every 9.5 days"
    purpose: "Moon phasing; yields 9.5-day Earth holding orbits"
  - label: "MAN-E"
    at: "perigee (point E, Fig. 3)"
    dv_ms: 14   # SOURCED: p.3 — "second maneuver (14 m/s of ΔV) at perigee (Fig. 3, E) to increase
                #  apogee back to true lunar distance"
    purpose: "Restore apogee to lunar distance"
  - label: "MAN-OOP"
    at: "0.5 days after MAN-E"
    dv_ms: 8    # SOURCED: p.3 — "third maneuver (8 m/s of ΔV) is performed 0.5 days later to correct
                #  for drift out of the lunar plane, likely due to solar perturbations"
    purpose: "Out-of-plane correction (solar perturbation)"
  # Per-cycle total (nominal): 35 m/s (13+14+8); range across apsidal rotation: 20–62 m/s/cycle
  # (SOURCED: p.4 — "Solar gravity perturbations cause variance in the ΔV requirements,
  #  from 20 to 62 m/s per cycle … needed to maintain the presented Earth-Moon cycler").
  # Abstract states average 39 m/s/month (p.1) — consistent with 35 m/s/26-day nominal.

# --- first_published: add AAS paper number note ---
first_published:
  note: |
    The preprint header reads "AAS 15-____" (number left blank at submission). The specific
    AAS-assigned paper number is not recoverable from the preprint PDF. NTRS accession
    20150018049; NASA report ARC-E-DAA-TN22765. No DOI found.
    Conference: AAS/AIAA Astrodynamics Specialist Conference, Vail CO, August 9–13, 2015.
```

---

## 8. CR3BP backfill assessment (#182)

**This source does NOT unblock the #182 genova-aldrin gap.**

The 3-petal orbit is **solar-perturbed and maneuver-maintained by design**. The paper explicitly
states it does not exist in the pure restricted three-body problem (p.3). There is no Jacobi
constant, no rotating-frame state vector, no non-dimensional period, and no stability index
anywhere in the paper. The `state_nd: null` field in the catalogue is correct; the
`jacobi_constant`, `period_nd`, `stability_index`, `lunit_km`, and `tunit_s` fields also remain
null.

The #182 `correct_periodic` corrector cannot be seeded from this paper. The gap would require
either:
1. A separate numerical search in the bicircular problem (Earth-Moon-Sun) for the analogous
   periodic family, or
2. A future paper that publishes rotating-frame initial conditions for this orbit family.

The existing `data_gaps` entries in the catalogue row remain accurate. The `source_hint`
"Genova & Aldrin 2015 (AAS-15) full PDF" can now be closed (the PDF has been read) and
replaced with "orbit does not exist in pure CR3BP; bicircular model required."

---

## 9. Fidelity caveat for the catalogue model_assumption field

The catalogue records `model_assumption: cr3bp` for this row. This is a **known mismatch** that
was flagged in the existing `data_gaps` block. The paper uses a high-fidelity force model
(50×50 Earth, 50×50 Moon, 4×0 Sun gravity; SRP not mentioned for cycler itself). The CR3BP
label was assigned because the orbit family is *related to* the CR3BP resonance structure, not
because the trajectory was computed in a pure CR3BP model. The mining note confirms this
caveat. Reclassification to `model_assumption: analytic-ephemeris` or a new `bicircular` value
would be more accurate but is a schema decision outside this mining task's scope.
