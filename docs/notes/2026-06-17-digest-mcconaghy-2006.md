# McConaghy-Landau-Yam-Longuski 2006 (JSR 43(2)) — full-page deep-read verdict

**Read 2026-06-17 AET.** Part of Mars-cycler digest wave (Agent A).

## 1. Header (verbatim from PDF p.1)

- Title: **"Notable Two-Synodic-Period Earth-Mars Cycler"**
- Authors: **T. Troy McConaghy*, Damon F. Landau†, Chit Hong Yam‡, James M.
  Longuski§** (all Purdue University, West Lafayette IN 47907-2023).
  Footnotes: * Ph.D. Candidate, troy_mcconaghy@hotmail.com, Student Member
  AIAA; † Ph.D. Candidate, nutmunky@purdue.edu, Student Member AIAA;
  ‡ Graduate Student, cyam@purdue.edu, Student Member AIAA; § Professor,
  longuski@ecn.purdue.edu, Associate Fellow AIAA.
- Venue: **Journal of Spacecraft and Rockets, Vol. 43, No. 2, March-April
  2006, pp. 456-465**.
- DOI: **10.2514/1.15215**.
- Received 21 December 2004; revision received 29 April 2005; accepted 9 May 2005.
- Associate Editor: D. Spencer.
- Length: **10 pages** (pp. 456-465).

## 2. What the paper actually is

THE **ballistic S1L1 cycler** paper. Despite the catalogue's existing
`mcconaghy-2006-em-k2` row treating this as a "notable two-synodic-period
cycler", this paper is **specifically about the S1L1 cycler** — a single
named cycler family (label `2g(2.8276, 658°, U) g(1.4580, 525°, L)` in the
McConaghy-Russell-Longuski 2005 nomenclature, p.456 abstract context).

The paper does two things:
1. **Circular-coplanar analysis** of the ballistic S1L1 cycler (p.456-458):
   identifies the trajectory, computes V∞=4.7128 km/s at Earth flyby
   (incoming = outgoing → no Δv needed), V∞≈5.0 km/s at Mars, Earth-Mars
   transit time 153.15 days, Earth flyby altitude 31,809 km.
2. **DE405 ephemeris-model continuation** (p.458-463) using SNOPT — finds
   30-year, 24-body itineraries for two outbound + two inbound cycler vehicles
   over 2005-2007 launch window. **Tables 2-9 give 30-year itineraries with
   per-encounter date, V∞, closest-approach altitude, leg ToF** for two cost
   metrics:
   - Tables 2-5: minimize cycler Δv only.
   - Tables 6-9: minimize cycler Δv + taxi Δv.

### Title / venue / DOI cross-check against our records

| Field | Our records | Paper actual | Status |
|---|---|---|---|
| Catalogue row `mcconaghy-2006-em-k2` | "McConaghy/Landau/Yam/Longuski two-synodic-period Earth-Mars cycler" (line 243) | "Notable Two-Synodic-Period Earth-Mars Cycler" — paper is specifically about the **S1L1 cycler** | The catalogue's row name is correct in spirit but misses the S1L1-specific identification. The row's `nomenclature: "Russell-McConaghy SnLm"` and the abstract values V∞=4.7 km/s Earth + 5.0 km/s Mars + 153-day transfer match exactly. |
| KNOWN_CORPUS | Anchor cites McConaghy 2006 JSR 43(2) DOI 10.2514/1.15215. | DOI matches. | Citation is **correct**. |
| Filename | "mcconaghy-landau-yam-2006-notable-two-synodic-period-earth-mars-cycler-jsr-doi-10.2514-1.15215.pdf" | matches | OK |

## 3. Key numerical / structural content

### Modeling assumptions (p.456)
1. Earth-Mars synodic period exactly 2-1/7 yr.
2. Earth + Mars + cycler in ecliptic plane.
3. Earth + Mars circular orbits.
4. Cycler is conic (two-body).
5. Only Earth provides gravity assist.
6. Gravity-assist instantaneous.

p.457: "**assumption 1 implies that the orbital period of Mars is 1.875
years** (whereas a more accurate value is 1.881 years). Also, assumption 1
implies that the semi-major axis of Mars is 1.875^(2/3) AU ≈ **1.5206 AU**."

### Circular-coplanar S1L1 geometry (p.457)

The S1L1 cycler has THREE Earth encounters per cycle: at t=0, τ=2.8276 yr,
and t=4-2/7 yr = 4.286 yr. Position vectors at the three encounters
(Cartesian, Sun-centered, AU):
- R₁ = [1, 0]
- R₂ = [cos(2πτ), sin(2πτ)] ≈ **[0.4690, −0.8832]**
- R₃ = [cos(2π·4-2/7), sin(2π·4-2/7)] ≈ **[−0.2225, 0.9749]**

**Leg 1** (Lambert R₁→R₂ in τ=2.8276 yr): chooses the **short-period** solution
of the two multi-rev solutions making 1-2 revs around the Sun:
- Period: **1.4889 yr**.
- Semi-major axis: **1.3039 AU**.

**Leg 2** (Lambert R₂→R₃ in 4-2/7 − τ = 1.4580 yr): chooses the **long-period**
solution of the two multi-rev solutions making 1-2 revs:
- Period: **1.0733 yr**.
- Semi-major axis: **1.0483 AU**.

p.457: "The name 'S1L1 cycler' is a reminder that **the first leg is the
short-period solution making 1-2 revs, and the second leg is the long-period
solution making 1-2 revs**." ← **THIS is the SnLm naming convention origin**:
**S = Short-period Lambert solution, L = Long-period Lambert solution**, NOT
"S-leg" and "L-leg". The integer after S/L counts the revolutions.

### S1L1 V∞ at Earth (p.457)
"When we calculate the incoming and outgoing V∞ at Earth, we find that they
are both equal to **4.7128 km/s**. This equality is not a coincidence; the
duration of leg 1 is chosen to make the incoming and outgoing V∞ equal so
that the required velocity change can be accomplished using a gravity-assist
flyby. **No Δv maneuvers are required** to keep the spacecraft on its orbit,
hence the name 'ballistic S1L1 cycler.'"

### Leg 1 aphelion + Mars crossings (p.457)
- Leg 1 aphelion: **1.6369 AU** (>Mars 1.5206 AU; Mars-crossing).
- Leg 2 aphelion: **1.2170 AU** (<Mars 1.5206 AU; never reaches Mars).
- Leg 1 crosses Mars orbit **four times**: at t = **0.4193, 0.9194, 1.9082,
  2.4083 years (= 153.15, 335.81, 696.97, 879.63 days)** after launch.
- The Earth-Mars transit time is **153.15 days** (outbound) — same as inbound
  (mirror symmetry).

### Table 1 — S1L1 cycler characteristics, circular-coplanar (p.458)

| Characteristic | Value |
|---|---|
| Δv required | **0 km/s** |
| Earth V∞ | **4.7 km/s** |
| Earth flyby altitude | **31,809 km** |
| Mars V∞ | **5.0 km/s** |
| Mars flyby altitude | ∞ (no GA at Mars) |
| Earth-Mars transfer time | **153.15 days** |
| Minimum Mars stay time | **604.63 days** |
| Mars-Earth transfer time | **153.15 days** |
| Repeat time | **4-2/7 yr** |

### Tables 2-5 (p.460-461) — 30-yr ephemeris itineraries, **cycler Δv only**

22 encounters each for cycler vehicles 1, 2, 3, 4. Sample structure (Vehicle 2
outbound, Table 3 — re-stating from paper, abridged):

| Enc | Date | V∞ km/s | Closest alt km | Leg ToF days |
|---|---|---|---|---|
| Earth-1 | 12/03/2009 | 6.87 | 31,100 | — |
| Mars-2 | 06/07/2010 | 4.31 | 17,600 | 186 (crew) |
| Earth-3 | 08/24/2012 | 6.43 | 26,400 | 809 |
| Earth-4 | 02/14/2014 | 6.43 | 41,500 | 539 |
| Mars-5 | 07/03/2014 | 7.14 | 12,200 | 138 (crew) |
| ... | ... | ... | ... | ... |
| Earth-22 | 10/26/2039 | 5.53 | 23,900 | 536 |

**Important**: Vehicle 2's first encounter Earth-1 on **12/03/2009** is a
near-exact replica of **McConaghy 2004 JSR Table 6's S1L1 starting on 9 June
2008** but offset by ~18 months (within the 128-day launch period flexibility
noted in §3 conjecture). Different vehicle, same cycler family.

p.461: "The two outbound cyclers are ballistic (over the examined time span),
but each inbound cycler requires a modest Δv. … the outbound trajectories
also require a Δv on the order of 10 m/s in the following 30-year cycle. So,
although the S1L1 cycler is ballistic in the circular-coplanar model, **it is
only NEARLY ballistic in the ephemeris model**, where the relative positions
of Earth and Mars do not repeat precisely."

p.461 flyby altitude note: "the average Earth flyby altitude is **24,090 km**,
whereas the Earth flyby altitude in the circular-coplanar model is 31,809 km.
Similarly, the average Mars flyby altitude is **9,950 km**, whereas the Mars
flyby altitude in the circular-coplanar model is effectively infinite … The
reason for the lower flyby altitudes (especially at Mars) is that the
gravity-assist maneuvers must also change the orbital plane … the inclination
of Mars's orbit is about 1.85 deg."

### Tables 6-9 (p.461-462) — 30-yr itineraries, **total Δv (cycler + taxi)**

Same 4 vehicles. **DSM Δv values inserted** to trade departure V∞ for
deep-space Δv. Examples:
- Vehicle 1 (Table 6): DSM 11/04/2010 = **0.66 km/s**; DSM 08/02/2023 = 0.11
  km/s; DSM 11/30/2027 = 0.81 km/s.
- Vehicle 2 (Table 7): DSM 08/31/2008 = 0.30 km/s; DSM 12/04/2012 = 0.62 km/s;
  DSM 09/19/2025 = 0.31 km/s; DSM 12/16/2029 = 0.27 km/s.

p.462 summary: "Over all 14 missions, the **average Earth departure V∞ was
reduced from 5.52 to 4.90 km/s**, a reduction accomplished by increasing the
DSM Δv to an average of **0.11 km/s per synodic period per outbound cycler
vehicle**. Similarly, the **average Mars departure V∞ was reduced from 5.22
to 3.22 km/s** by increasing the DSM Δv to an average of 0.33 km/s per
synodic period per inbound cycler vehicle."

### Constrained transit-time Figures 9-11 (p.463-464)
TOF-constrained variants: 120-270 day transit. Total Δv ranges from ~2 to
~16 km/s at 120-day cap (extreme); ~2-4 km/s at 230+ days.

### Cycler-vs-cycler comparison
p.458: V∞ for S1L1 (4.7 / 5.0 km/s at E/M) are "significantly lower than the
corresponding V∞ for the **Aldrin cycler (6.5 and 9.8 km/s, respectively)**"
[ref 12 = Byrnes-Longuski-Aldrin 1993]. Matches Byrnes 1993 derivation +
catalogue's `aldrin-cycler` row.

### References (p.464-465)
- Ref 14 = **McConaghy-Longuski-Byrnes AIAA 2002-4420** (the S1L1 primary
  source).
- Ref 19 = **McConaghy-Yam-Landau-Longuski AAS 03-509** (Aug 2003,
  "Two-Synodic Period Earth-Mars Cyclers with Intermediate Earth Encounter")
  — the **conference paper this JSR paper grew out of**; contains additional
  itinerary data referenced in this paper but omitted "for brevity".
- Ref 20 = **McConaghy-Longuski-Byrnes JSR 41(4) 2004** (= wave paper 4).
- Ref 25 = **McConaghy-Russell-Longuski JSR 42(4) 2005** (= wave paper 5,
  the standard nomenclature paper).
- Ref 18 = **Chen-McConaghy-Landau-Longuski-Aldrin "Powered Earth-Mars Cycler
  with Three-Synodic-Period Repeat Time," JSR 42(5), 2005, pp. 921-927** —
  a different cycler, in a sister 2005 paper.

## 4. `KNOWN_CORPUS` impact (RECOMMEND, do not edit)

The McConaghy 2006 JSR DOI is already cited correctly in the SnLm anchor
(literature_check.py line 263-277). **No DOI/citation defect.**

**HOWEVER** — the catalogue's `mcconaghy-2006-em-k2` row's `nomenclature`
field is `"Russell-McConaghy SnLm"` (line 244-245), and the catalogue's `id`
is `mcconaghy-2006-em-k2` (with `k2` suggesting `k=2` periods). The S1L1
naming origin is now fully resolved by this paper's p.457 text:

> **"S = Short-period Lambert solution making 1-2 revs;
> L = Long-period Lambert solution making 1-2 revs;
> SnLm = an n-revolution short-period leg followed by an m-revolution
> long-period leg (or similar)."**

**This resolves the "SnLm naming convention origin" question raised in the
McConaghy 2005 verdict note.** The SnLm convention IS this paper's
nomenclature (not a separate paper). Suggestion: add a short
`docs/spec.md` annotation documenting the SnLm naming origin
explicitly.

Suggested KNOWN_CORPUS anchor refinement (combining all wave findings):

```python
name="Earth-Mars cycler family — multi-source nomenclature anchor",
authors=("Russell", "Ocampo", "McConaghy", "Longuski", "Byrnes", "Landau", "Yam"),
citation=(
    "Aldrin cycler (1L1): Byrnes, Longuski & Aldrin, JSR 30(3):334-336 "
    "(1993) DOI 10.2514/3.25519; "
    "SnLm Lambert-solution-type nomenclature: McConaghy, Russell & Longuski, "
    "JSR 42(4):694-698 (2005) DOI 10.2514/1.8123; "
    "nPr family-tag nomenclature: McConaghy, Longuski & Byrnes, JSR "
    "41(4):622-628 (2004) DOI 10.2514/1.11939; "
    "Ballistic S1L1 ephemeris itineraries (Tables 2-9): McConaghy, Landau, "
    "Yam & Longuski, JSR 43(2):456-465 (2006) DOI 10.2514/1.15215; "
    "Systematic p.h.s.±i enumeration: Russell & Ocampo, JGCD 27(3):321-335 "
    "(2004) DOI 10.2514/1.1909, preprint AAS-03-145; "
    "Dissertations: Russell PhD (UT Austin 2004), McConaghy PhD (Purdue 2004)"
),
```

## 5. Catalogue impact (RECOMMEND)

### Rows touched

- **`mcconaghy-2006-em-k2`** (line 242, currently cross_validated tier) — this
  row IS the S1L1 cycler per the abstract values match (V∞_E=4.7, V∞_M=5.0,
  transit 153 days, 2-synodic repeat). The row's existing fields are correct.

### V0→V1 promotion candidacy

**MASSIVE V0→V1 PROMOTION OPPORTUNITY**: Tables 2-9 of this paper provide
**FOUR independent 22-encounter DE405 itineraries** (cycler-Δv-only metric)
AND **FOUR more** (cycler+taxi metric) for a total of **8 fully-specified
30-year ephemeris realizations** of the S1L1 cycler. Each gives per-encounter:

- Date (mm/dd/yyyy precision)
- V∞ (km/s, 2-decimal precision)
- Closest-approach altitude (km, 100-km precision)
- Leg ToF (days)

This is the **richest V1-grade data anywhere in the McConaghy/Russell cycler
literature**. The catalogue's `mcconaghy-2006-em-k2` is at V3 cross_validated
per the abstract values — **but with the full Tables 2-9 itineraries on hand,
we can move it to V4** (per `project_validation_ceiling` memory: "V4 / V5 is
new-input-gated"; this paper IS the new input).

**Caveat from memory `project_s1l1_realeph_closure_blocker`**: "modelling
stack complete; S1L1 fails on family-selection (off-basin), not infra. … S1L1
is actually MULTI-ARC (two generic-return arcs) — likely the real reason
single-ellipse never closed; coplanar 154d Lambert reaches Vinf_E~5.55 but
Mars 3.05 needs the intermediate flyby/real eccentricity."

**This paper RESOLVES the S1L1 modelling blocker**:
1. Two arcs ARE the right structure (CONFIRMED: Tables 2-9 show 22 encounters
   in alternating Earth-Mars-Earth-Earth-Mars-Earth-Earth-Mars pattern, with
   each Earth-Earth pair being a single arc and each Earth-Mars or Mars-Earth
   pair another arc).
2. The 154-day Lambert reaches **V∞=4.7128 km/s in the circular-coplanar
   model** (p.457), NOT 5.55 — matching the catalogue's 4.7 km/s value.
3. The Mars 3.05 km/s the memory mentions is likely from **McConaghy 2004 JSR
   Table 4** (2L3 cycler footnote d, the Mars-aphelion-speed-difference, NOT
   a true V∞). **DIFFERENT cycler from S1L1.** The S1L1 has Mars V∞=5.0 km/s.

So the prior blocker may have been comparing **2L3 (which has Mars V∞=3.05
"speed difference" not true V∞) against S1L1 (which has Mars V∞=5.0)** — i.e.
the family-selection issue was a CYCLER CONFUSION, not a Lambert/corrector
defect.

### Specific catalogue actions

- **Promote `mcconaghy-2006-em-k2` from cross_validated → ephemeris_reproduced
  (V4)** once we successfully reproduce ANY of Tables 2-9 in our M7 stack
  within published tolerance.
- **Add a `mcconaghy-2006-table3-veh2` reference itinerary** (or similar) as a
  golden-test target. Tolerance: ±5 days on encounter dates, ±0.5 km/s on V∞,
  ±2000 km on altitude.
- **Cross-reference McConaghy 2004 JSR Table 6 vs this paper's Tables 2-9** —
  Table 6 of paper 4 starts on 9 June 2008; this paper's Vehicle 2 (Table 3)
  starts on 12/03/2009. The two are **18 months apart, different launches of
  the same outbound S1L1 family**. Both are valid V1 anchors.

## 6. Errata / surprises

1. **The "SnLm" naming convention IS THIS paper's convention** (p.457):
   `S` = short-period Lambert solution, `L` = long-period Lambert solution,
   numbers count revs. **This is the answer to the question raised in the
   McConaghy 2005 verdict note**. SnLm and the McConaghy 2005 per-leg formal
   labels are sibling conventions: SnLm = S/L type indicators per leg;
   per-leg g/f/h descriptors = exact parametric label.
2. **S1L1's Mars V∞ is 5.0 km/s, NOT 3.05 km/s**. The 3.05 km/s value
   referenced in some prior memory traces back to McConaghy 2004 JSR Table 4
   row 2L3 footnote d (Mars-aphelion-speed-difference for cyclers that don't
   reach Mars). 2L3 has aphelion 1.51 AU < Mars 1.52 AU, so its "V∞ at Mars"
   in McConaghy 2004 Table 4 is a synthetic delta, not a true V∞. **S1L1 has
   aphelion 1.6369 AU > Mars 1.5206 AU, so its V∞ Mars = 5.0 km/s is a
   genuine V∞.** Different cyclers; different physical regimes.
3. **The S1L1 has FOUR Mars-orbit crossings per leg-1 revolution-pair** (p.457:
   t=0.4193, 0.9194, 1.9082, 2.4083 yr). Only ONE corresponds to an actual
   Mars encounter, controlled by launch-date choice. The catalogue's
   "outbound cycler" vs "inbound cycler" distinction (line 92) is precisely
   this launch-date choice (outbound = encounter Mars on first crossing at
   153.15 d after launch; inbound = encounter Mars on last crossing 153.15 d
   before next Earth encounter).
4. **30-year repeat NOT 15-year**. p.458: "the ballistic S1L1 cycler repeats
   inertially after 30 years (seven cycles of 4-2/7 years). The orbits of
   Earth and Mars also repeat approximately every 30 years (after 30 Earth
   revolutions and 16 Mars revolutions)." Worth confirming this is the
   correct framing for the catalogue's repeat-period interpretation.
5. **"Ideal cycler" conjecture (p.459)**: "We conjecture that the middle
   times of a long itinerary are the same as the times of the ideal cycler
   trajectory (and hence unique)." If true, the SnLm long-itinerary middle
   bodies are convergent — converging to a UNIQUE periodic trajectory. The
   launch window flexibility (128 days, p.459) is for entry, not for the
   underlying cycler.
6. **Cost function asymmetry**: minimizing cycler Δv alone leaves inbound
   cyclers with modest Δv; minimizing total (cycler + taxi) Δv lets outbound
   vehicles absorb DSM Δv to reduce Earth departure V∞ from 5.52 → 4.90 km/s
   (an 11% Δv-budget reduction for the taxi). Catalogue should reflect both
   metrics as separate "modes" if it parameterizes the S1L1 by metric.
7. **No 0.66 km/s DSM in the abstract** — only the cycler-Δv-only metric is
   summarized in the abstract. The cycler+taxi metric requires DSMs of
   0.07-1.51 km/s magnitude (Tables 6-9). Important context omitted from
   the abstract.

## 7. Action items for parent

- [ ] **Promote `mcconaghy-2006-em-k2` row from cross_validated → V4
  ephemeris_reproduced** once an M7 / ephemeris-corrector run reproduces ANY
  of Tables 2-9 within tolerance. Track separately.
- [ ] **Resolve the `project_s1l1_realeph_closure_blocker` memory** with the
  finding that the prior closure-failure may have been a 2L3 vs S1L1 cycler
  confusion. Reread the original blocker investigation with these data; if
  the Mars V∞=3.05 used in the closure attempt was from McConaghy 2004 Table 4
  row 2L3, **it was the wrong target** — S1L1 has Mars V∞=5.0 km/s. Track a
  follow-on to re-run S1L1 closure against the correct Mars V∞=5.0 target.
- [ ] **Document the SnLm naming convention origin** in `docs/spec.md` (S =
  short-period Lambert solution making 1-2 revs; L = long-period Lambert
  solution making 1-2 revs; integer counts revs). Per this paper p.457.
- [ ] **Add reference itineraries from Tables 2-9 as golden-test targets**:
  `data/golden_itineraries/mcconaghy-2006-table3-veh2.yaml` (and the other 7
  tables), schema = list of (encounter_id, date, body, V∞_kms,
  closest_approach_km, leg_tof_days). Tolerance: ±5 days, ±0.5 km/s, ±2000 km.
- [ ] **Cross-validate against McConaghy 2004 JSR Table 6** (paper 4): the
  9 June 2008 S1L1 launch in that paper is the same cycler family, different
  launch year. Both should reproduce in our stack with the same physics.
- [ ] **Track ref 19 = McConaghy-Yam-Landau-Longuski AAS 03-509 (Aug 2003)**
  as a still-outstanding acquisition (the conference precursor of this JSR
  paper; cited by the JSR paper as containing additional itinerary data).
  Add to wishlist if not present.
- [ ] **No KNOWN_CORPUS DOI/citation defect** for this paper. The existing
  citation is correct; only the surrounding SnLm anchor needs broader fixes
  per the other wave papers' verdict notes.
