# Byrnes-Longuski-Aldrin 1993 (JSR 30(3)) — full-page deep-read verdict

**Read 2026-06-17 AET.** Part of Mars-cycler digest wave (Agent A).

## 1. Header (verbatim from PDF p.1)

- Title: **"Cycler Orbit Between Earth and Mars"**
- Authors: **Dennis V. Byrnes* (Jet Propulsion Laboratory, Caltech, Pasadena CA
  91109), James M. Longuski† (Purdue University, West Lafayette IN 47907),
  Buzz Aldrin‡ (Starcraft Enterprises, Laguna Beach CA 92651)**.
  Footnotes: * Member Technical Staff, Mission Design Section, Member AIAA;
  † Assistant Professor, School of Aeronautics and Astronautics, Senior Member
  AIAA; ‡ Astronaut, Fellow AIAA.
- Venue: **Journal of Spacecraft and Rockets, Vol. 30, No. 3, May-June 1993,
  pp. 334-336**.
- DOI: **10.2514/3.25519**.
- Received Aug. 17 1991; revision received March 16 1992; accepted March 16 1992.
- Associate Editor: James A. Martin.
- Length: **3 pages** (pp. 334-336).

## 2. What the paper actually is

The **first peer-reviewed publication of the Aldrin cycler**. Aldrin's
original concept appeared as a 1985 SAIC internal report (ref 8: "Cyclic
Trajectory Concepts," Science Applications International Corp., Aerospace
Systems Group, Hermosa Beach CA, Oct. 1985). Friedlander-Niehoff-Byrnes-
Longuski AIAA-86-2009 (ref 6) was the first conference paper analyzing the
concept. This JSR paper is the **archival JSR record**.

The paper does two things:
1. **Circular-coplanar analytic derivation** of the cycler's orbital elements
   (a, e) from the 15/7 synodic-period precession constraint (Δψ = 2π/7 rad =
   51.4°), via a multi-revolution Lambert solve.
2. **Real-ephemeris numerical results** for both outbound and inbound cyclers
   over a complete 15-year cycle 1995-2010, using a "multiconic optimization
   program" (D'Amario refs 10, 11) with a 1000 km minimum Earth-flyby altitude
   constraint.

### Title / venue / DOI cross-check against our records

| Field | Our records | Paper actual | Status |
|---|---|---|---|
| Authors | "Aldrin, Byrnes, McConaghy, Longuski" (KNOWN_CORPUS line 257) | **Byrnes, Longuski, Aldrin** (3-author 1993 paper). McConaghy is NOT on this paper. | The KNOWN_CORPUS Aldrin anchor mixes 1993 + 2002 author lists. Defensible if it's meant to cover both, but worth disambiguating. |
| Title | KNOWN_CORPUS doesn't print a title for the Aldrin anchor | "Cycler Orbit Between Earth and Mars" | — |
| Venue (KNOWN_CORPUS) | "Byrnes, McConaghy & Longuski, AIAA 2002-4420 (Aldrin cycler)" (line 259) | This paper is JSR 30(3) 1993, NOT AIAA 2002-4420. **Different paper**. McConaghy 2002 is the AIAA 2002-4420; THIS is its 1993 progenitor. | Errata-grade gap: the JSR 1993 paper (DOI 10.2514/3.25519) is NOT cited in the KNOWN_CORPUS Aldrin anchor at all. The catalogue's `aldrin-cycler` row DOES cite it explicitly (corroborating_sources, lines 133-147 of catalogue.yaml). |
| Filename | "byrnes-longuski-aldrin-1993-cycler-orbit-earth-mars-jsr-doi-10.2514-3.25519.pdf" | matches | OK |

## 3. Key numerical / structural content

### Circular-coplanar analytic derivation (p.334-335)
- Earth period **1 yr** (e=0.0168 in real world, 0 in this analysis); Mars
  period **1.8808 yr** (e=0.0934, i=1.85° in real world, **approximated as
  1.875 yr** for analytic 15-year resonance) (p.334).
- Synodic period: 15/7 = **2.1429 yr** (p.334).
- Δψ = 2π/7 rad = 51.4° precession per cycler orbit (p.334).
- True-anomaly Earth encounter: **θ = ±25.7°** (= Δψ/2) (p.334-5).
- Approximation: a ≈ (2)^(2/3) = **1.59 AU**, e ≈ **0.387** (p.335, eq above 2).
- **Exact Lambert solution** (Eqs 2-3, p.335): **a = 1.60 AU, e = 0.393, orbit
  period 2.02 yr**. (This is the 1.60/0.393 pair the catalogue's `aldrin-cycler`
  row uses — vs the "spec.md a≈1.659, e≈0.41, peri≈0.98, apo≈2.34" alternative
  from the early SAIC presentation that the catalogue notes as a discrepancy.)
- Eq 4: γ = 7.18° (flight-path angle at Earth).
- V_E = 29.8 km/s (Earth heliocentric), V = 34.9 km/s (s/c heliocentric at
  Earth), **ΔV = 8.73 km/s** (heliocentric).
- Eq 7-8: **V∞ = 6.54 km/s** at Earth flyby; **turn angle 2δ = 83.8°**;
  required perigee r_p = 4640 km (p.335). **Earth radius 6371 km, so the
  required flyby is BELOW THE SURFACE — physically infeasible without a
  maintenance maneuver.**
- Stated maintenance need: **~230 m/s at aphelion to rotate ω enough to make
  the flyby occur at 1000 km altitude above Earth's surface** in the
  circular-coplanar model (p.335).

### Numerical results: outbound cycler 1995-2010 (Table 1, p.336)

**Encounter conditions** for outbound cycler (E1→M2→E3→M4...E15), 19 Nov 1996
launch (Earth-1, V∞=6.19 km/s) through 13 Nov 2011 (Earth-15, V∞=5.81 km/s):

| Encounter | Date | V∞ (km/s) | r_p (planet radii) |
|---|---|---|---|
| Earth-1 | 19 Nov 1996 | 6.19 (launch) | — |
| Mars-2 | 01 May 1997 | 10.69 | 5.8 |
| Earth-3 | 01 Jan 1999 | 5.94 | 1.3 |
| Mars-4 | 28 May 1999 | 11.74 | 29.1 |
| Earth-5 | 08 Feb 2001 | 5.67 | 1.2 |
| Mars-6 | 06 July 2001 | 10.22 | 1.3 |
| **Maneuver** | 13 March 2002 | **0.54 (ΔV₁)** | — |
| Earth-7 | 16 April 2003 | 5.67 | 1.2 |
| Mars-8 | 12 Sept 2003 | 7.28 | 1.3 |
| **Maneuver** | 17 May 2004 | **0.74 (ΔV₂)** | — |
| Earth-9 | 07 July 2005 | 5.87 | 1.2 |
| Mars-10 | 13 Dec 2005 | 6.05 | 3.4 |
| **Maneuver** | 23 July 2006 | **0.45 (ΔV₃)** | — |
| Earth-11 | 06 Sept 2007 | 5.87 | 1.8 |
| Mars-12 | 16 Feb 2008 | 7.43 | 6.4 |
| Earth-13 | 10 Oct 2009 | 5.89 | 1.9 |
| Mars-14 | 28 March 2010 | 8.66 | 5.0 |
| Earth-15 | 13 Nov 2011 | 5.81 | 1.9 |

(Maneuver units are km/s per the table header "Approach V∞ km/s" — but the
text on p.336 says "the optimal trajectory covering this 15-year cycle requires
propulsive maneuvers on only three of the seven orbits. ... the sum of these
three maneuvers is approximately seven times the per-orbit requirement of the
circular coplanar problem" — i.e. 0.54+0.74+0.45 ≈ 1.73 km/s ≈ 7×0.23 km/s ×
1.07. So ΔV₁,₂,₃ ARE km/s, not m/s, in Table 1. Confirms maintenance.)

### Outbound cycler V∞ ranges
- **Earth flybys: 5.67 to 6.19 km/s** (matches the catalogue's Rogers-2012
  Table 4 quote of "V_inf,flyby ranging 5.51-6.55 km/s for the 4:3(2)- and
  3:2(1)- Aldrin trajectories" in `aldrin-cycler` row line 53 note).
- **Mars flybys: 6.05 to 11.74 km/s** (much wider spread). Note: catalogue
  currently quotes Mars V∞=9.7 km/s for Aldrin (line 55), which is the
  middle-of-range Russell-2004 dissertation value, NOT the per-encounter Byrnes
  1993 values.

### Numerical results: inbound cycler 1995-2010 (Table 2, p.336)

| Encounter | Date | V∞ (km/s) | r_p (planet radii) |
|---|---|---|---|
| Earth-1 | 05 June 1995 | 5.88 (launch) | — |
| Mars-2 | 20 Jan 1997 | 8.52 | 5.5 |
| Earth-3 | 09 July 1997 | 5.95 | 1.8 |
| Mars-4 | 07 March 1999 | 7.35 | 9.4 |
| Earth-5 | 17 Aug 1999 | 6.01 | 1.4 |
| **Maneuver** | 28 Sept 2000 | **0.27 (ΔV₁)** | — |
| Mars-6 | 15 May 2001 | 6.60 | 5.2 |
| Earth-7 | 08 Oct 2001 | 5.88 | 1.2 |
| **Maneuver** | 04 Dec 2002 | **1.11 (ΔV₂)** | — |
| Mars-8 | 07 Aug 2003 | 7.30 | 1.3 |
| Earth-9 | 02 Jan 2004 | 5.39 | 1.4 |
| **Maneuver** | 02 Feb 2005 | **0.66 (ΔV₃)** | — |
| Mars-10 | 10 Oct 2005 | 9.96 | 1.3 |
| Earth-11 | 12 March 2006 | 5.48 | 1.5 |
| Mars-12 | 19 Nov 2007 | 11.59 | 8.4 |
| Earth-13 | 16 April 2008 | 5.96 | 1.5 |
| Mars-14 | 13 Dec 2009 | 10.55 | 5.0 |
| Earth-15 | 22 May 2010 | 5.93 | 1.8 |

### Per-cycle total ΔV
- Outbound: 0.54 + 0.74 + 0.45 = **1.73 km/s over 15 years (~3 cycles)**
  ≈ 0.58 km/s per cycle (NOT 230 m/s — text says "approximately seven times
  the per-orbit requirement of the circular coplanar problem" meaning 7×0.23 =
  1.61 km/s, close to the 1.73 sum).
- Inbound: 0.27 + 1.11 + 0.66 = **2.04 km/s over 15 years (~3 cycles)**.

### Earth-Mars transit time
- p.336: **"flight-time variation is not large: 147-170 days for both
  Earth-to-Mars transits on the outbound cycler and for Mars-to-Earth transits
  on the inbound cycler"**.
- This matches the catalogue's `aldrin-cycler` outbound transit_time_days=146-161
  range (line 91 "146-d is circular-coplanar; Rogers et al. 2012 reports 161-172
  d for the K:L(M) variants").

### Cycler orbital elements
- p.335 Eqs (2)-(3): **a = 1.60 AU, e = 0.393, period 2.02 yr**.
  (NOT the 1.659 AU / 0.41 / peri 0.98 / apo 2.34 spec.md alternative — the
  Byrnes 1993 JSR value is unambiguously a=1.60, e=0.393.)

## 4. `KNOWN_CORPUS` impact (RECOMMEND, do not edit)

The current Aldrin anchor (literature_check.py line 256-260):

```
authors=("Aldrin", "Byrnes", "McConaghy", "Longuski"),
citation="Byrnes, McConaghy & Longuski, AIAA 2002-4420 (Aldrin cycler)",
```

**DEFECTS / GAPS**:
1. **The 1993 JSR paper is the PRIMARY published Aldrin cycler paper** but is
   NOT cited in the KNOWN_CORPUS Aldrin anchor. McConaghy is NOT a co-author of
   this paper. AIAA 2002-4420 (McConaghy/Longuski/Byrnes) is a separate, later
   paper. The Aldrin anchor's citation should be **Byrnes, Longuski & Aldrin,
   JSR 30(3):334-336 (1993), DOI 10.2514/3.25519**, with AIAA 2002-4420 as a
   secondary/follow-up citation.
2. The authors tuple should be tightened to include the actual 1993 JSR author
   list **(Byrnes, Longuski, Aldrin)**, distinct from the 2002 AIAA list
   **(McConaghy, Longuski, Byrnes)**. Optionally combine as a (paper-tagged)
   set.

Recommended fix:

```python
authors=("Aldrin", "Byrnes", "Longuski", "McConaghy"),  # combined 1993+2002 author rosters
citation=(
    "Byrnes, D. V., Longuski, J. M. & Aldrin, B., 'Cycler Orbit Between "
    "Earth and Mars,' JSR 30(3):334-336 (1993), DOI 10.2514/3.25519 "
    "(canonical Aldrin cycler paper); McConaghy, T. T., Longuski, J. M. & "
    "Byrnes, D. V., 'Analysis of a Broad Class of Earth-Mars Cycler "
    "Trajectories,' AIAA-2002-4420 (Aldrin-class follow-up)"
),
```

## 5. Catalogue impact (RECOMMEND)

### Rows touched

- **`aldrin-cycler`** (id 1, currently V3 cross_validated) — corroborating_sources
  list at line 135 ALREADY cites Byrnes/Longuski/Aldrin 1993 with DOI
  10.2514/3.25519 (verified by line-grep on `byrnes.*aldrin`). The catalogue
  citation is **correct here**; only the KNOWN_CORPUS anchor needs the fix.

### V0→V1 promotion candidacy

The catalogue's `aldrin-cycler` row is already **V3 cross_validated**. The
question is whether Byrnes 1993 gives us new state to lift specific OTHER
catalogue rows.

- **HONEST NEGATIVE for new V1 promotions.** Byrnes 1993 publishes only:
  - Circular-coplanar (a, e, period) for the heliocentric ellipse — already in
    `aldrin-cycler` (a=1.60, e=0.393).
  - Per-encounter (date, V∞, r_p in planet radii) for the 15-year 1995-2010
    real-ephemeris case — Tables 1 and 2.
  - Three maneuver ΔV values per cycler (km/s, no direction).
- **No per-arc orbital elements** for the individual real-ephemeris arcs.
- **No multi-revolution Lambert state vectors** in machine-precision form.

For our `aldrin-cycler` row to gain a NEW V-level from Byrnes 1993, we'd need
to reproduce the 15-year Table 1/2 V∞ sequence in our ephemeris stack within
some tolerance. That is a **V3→V4 (or higher) promotion** path, not a V0→V1
path — it requires the cycler row to ALREADY be in good shape (which it is at
V3). **Not actionable as a V0→V1 step** in this task.

### Cross-validation opportunity (V3→V4 ahead)

Tables 1 and 2 give a **15-year V∞ sequence** that our M7 / ephemeris-corrector
stack could ATTEMPT to reproduce. The published-vs-our delta would be a
genuine V4-level test:
- Outbound launch V∞ = 6.19 km/s on 19 Nov 1996.
- Mars-2 arrival V∞ = 10.69 km/s on 01 May 1997 (≈163 d transit).
- Net 15-year ΔV outbound = 1.73 km/s, inbound = 2.04 km/s.

This is the "real-ephemeris ground truth" the `aldrin-cycler` row's
`maintenance_dv_kms_per_synodic: 1.52` field's COMPUTED note alludes to.
Byrnes 1993 gives published values to compare against. Worth a tracked task.

## 6. Errata / surprises

1. **AIAA 2002-4420 ≠ Byrnes 1993 JSR.** The KNOWN_CORPUS Aldrin anchor cites
   only the 2002 follow-up, not the 1993 JSR archival paper. This is a citation
   gap, not a numerical defect.
2. **Mars period: 1.8808 yr (real) vs 1.875 yr (analytic 15-year resonance
   approximation)**. The 1.875 value enables exact 15:8 Earth:Mars resonance.
   Russell-Ocampo 2003 uses the same 1.875 yr value for the same reason.
3. **a=1.60, e=0.393, period 2.02 yr** is the Byrnes 1993 value (Eqs 2-3).
   The "spec.md a≈1.659, e≈0.41, peri≈0.98, apo≈2.34" alternative on the
   catalogue's `aldrin-cycler` row note (lines 63) appears to be from a different
   source — likely the original Aldrin 1985 SAIC presentation. **The Byrnes 1993
   value should be the catalogue's preferred value** since it's the
   peer-reviewed publication. The catalogue currently flags this discrepancy
   but doesn't resolve it — Byrnes 1993 is the resolution.
4. **3 maneuvers / 7 orbits, not 7 maneuvers / 7 orbits**. p.336:
   "the optimal trajectory covering this 15-year cycle requires propulsive
   maneuvers on only three of the seven orbits". This is a SURPRISE compared
   to the circular-coplanar 1-maneuver-per-orbit picture. The real-world
   Aldrin cycler is **MUCH cleaner** than the circular-coplanar derivation
   suggests.
5. **Outbound vs inbound asymmetry**: outbound launch V∞ = 6.19 km/s, inbound
   launch V∞ = 5.88 km/s. The two cyclers are NOT exact mirrors in the
   real-ephemeris model (the text on p.334 hedges this: "are essentially mirror
   images of each other"). Catalogue should reflect both.
6. **NO maneuver direction / location data** beyond "near aphelion ~8 months
   post-Mars (outbound) / pre-Mars (inbound)" (p.336). Insufficient for
   independent reproduction of the optimization result.

## 7. Action items for parent

- [ ] **Fix KNOWN_CORPUS Aldrin anchor** — add the 1993 JSR paper with DOI
  10.2514/3.25519 as the **primary** Aldrin cycler citation, keep AIAA 2002-4420
  as secondary follow-up. Already-tracked-elsewhere check: this matches the
  wave digest's "Action: deep-read for KNOWN_CORPUS citation strengthening +
  Aldrin row provenance audit".
- [ ] **Add a `byrnes1993-table1-outbound` and `byrnes1993-table2-inbound`
  source tag** to the provenance vocabulary if not present, and consider
  populating the `aldrin-cycler` row's `vinf_kms_at_encounters` `note` fields
  with the per-encounter Byrnes 1993 ranges (5.67-6.19 km/s at Earth outbound,
  5.39-6.01 km/s at Earth inbound, 6.05-11.74 km/s at Mars).
- [ ] **Resolve the a=1.60 vs a=1.659 discrepancy** flagged in catalogue line 63
  by declaring Byrnes 1993 (a=1.60, e=0.393) as the catalogue-of-record value
  and demoting the a=1.659 alternative to a "see also" footnote. Optionally
  keep the test tolerances ±0.01 AU / ±0.02 e wide enough to absorb both.
- [ ] **V3→V4 ephemeris-reproduction sanity check** — track a follow-on task
  to attempt reproducing Byrnes 1993 Table 1's V∞ sequence with our M7 stack.
  Comparison metric: per-encounter V∞ difference vs Table 1.
- [ ] **Document the "3 maneuvers in 15 years" finding** in the catalogue's
  Aldrin maintenance discussion — it's a structural fact that helps interpret
  the published 1.73 km/s figure as cycle-averaged ~0.58 km/s, NOT
  per-encounter ~0.23 km/s as the circular-coplanar derivation predicts.
- [ ] **NO action on V0→V1 promotions** from this paper alone — Byrnes 1993
  doesn't print per-arc state for any cycler beyond the canonical Aldrin one
  (which is already V3).
