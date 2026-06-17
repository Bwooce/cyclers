# McConaghy-Longuski-Byrnes 2004 (JSR 41(4)) — full-page deep-read verdict

**Read 2026-06-17 AET.** Part of Mars-cycler digest wave (Agent A).

## 1. Header (verbatim from PDF p.1)

- Title: **"Analysis of a Class of Earth-Mars Cycler Trajectories"**
- Authors: **T. Troy McConaghy*, James M. Longuski† (Purdue University, West
  Lafayette IN 47907-2023), and Dennis V. Byrnes‡ (Jet Propulsion Laboratory,
  Caltech, Pasadena CA 91109-8099)**.
  Footnotes: * Ph.D. Candidate, Student Member AIAA; † Professor, Associate
  Fellow AIAA; ‡ Principal Engineer, Deputy Manager, Navigation and Mission
  Design Section, Senior Member AIAA.
- Venue: **Journal of Spacecraft and Rockets, Vol. 41, No. 4, July-August
  2004, pp. 622-628**.
- DOI: **10.2514/1.11939**.
- Received 28 January 2003; revision received 16 June 2003; accepted for
  publication 18 July 2003.
- Associate Editor: C. A. Kluever.
- Length: **7 pages** of content (pp. 622-628), plus a half-page advertisement.

## 2. What the paper actually is

The **nPr nomenclature** paper. This is McConaghy's pre-SnLm naming scheme.
A cycler is denoted `nPr±` where:
- **n** = number of synodic periods before repeat,
- **P** = one of `U` (unique-period, only one cycler exists for that n), `L`
  (long-period), or `S` (short-period),
- **r** = number of complete heliocentric revolutions before repeat,
- **±** indicates the sign of ω, the argument of periapsis (only used for `n=7`
  family).

The method (Section "Methodology", p.622-623):
1. Choose **n** (synodic periods before repeat) → derives the position of Earth
   at start and end (Eqs. 2-3).
2. Solve the multi-revolution Lambert problem from R₁=(1,0) to R₂=
   (cos(2πnS), sin(2πnS)) where S = 2-1/7 years (synodic period).
3. **For each n there are multiple solutions** (Fig. 3 shows 7 solutions for n=1
   in a time-of-flight vs period plot).
4. Identify each by `nPr` where P ∈ {U, L, S} and r is the rev count.

Tables 1, 2, 3 (p.623-624) show GRAPHICAL DRAWINGS of all cyclers for n=1..4,
parameterized by r=1..9. Specifically: Table 1 = unique-period (nU0), Table 2 =
long-period (nLr), Table 3 = short-period (nSr). Many entries are "NS" (No
Solution).

Table 4 (p.625) lists **21 of the most promising cyclers** (1 ≤ n ≤ 6) with
numerical properties (Aphelion radius AU, V∞ at Earth km/s, V∞ at Mars km/s,
shortest transfer time days, required turn angle deg, max possible turn angle
deg). Table 5 (p.626) lists **14 n=7 cyclers** that are ALL ballistic (because
their line of apsides doesn't need to be rotated). Table 6 (p.627) lists a
**22-encounter outbound ballistic S1L1 cycler itinerary using DE405**.

### Title / venue / DOI cross-check against our records

| Field | Our records | Paper actual | Status |
|---|---|---|---|
| Catalogue cross-references | The McConaghy 2004 dissertation note (`docs/notes/2026-06-10-mcconaghy-2004-dissertation-mining.md`) refers to a 2004 PhD dissertation — DIFFERENT from this JSR paper. | This is the JSR paper, the Ph.D. candidate sub-paper. **McConaghy was still a Ph.D. candidate at submission (28 Jan 2003)**, completing the dissertation in 2004. The dissertation expands on this paper. | Two distinct works — both cited in catalogue line 705-718. |
| KNOWN_CORPUS | "Russell-Ocampo / McConaghy Earth-Mars SnLm cyclers" anchor (line 263) cites McConaghy 2006 JSR but NOT this 2004 JSR paper. | This paper uses **nPr** naming, NOT **SnLm** (SnLm was introduced in McConaghy-Russell-Longuski 2005 JSR 42(4)). | Citation gap: McConaghy 2004 JSR is the canonical pre-SnLm cycler analysis paper. |
| Filename | "mcconaghy-longuski-byrnes-2004-analysis-class-earth-mars-cycler-trajectories-jsr-doi-10.2514-1.11939.pdf" | matches title | OK |

## 3. Key numerical / structural content

### Solar-system model
Same as Russell-Ocampo 2003 (circular coplanar, Mars period = 1.875 yr,
Earth = 1.0 yr); see Methodology p.622 ("Assumptions 2 and 3 allow us to set
up a planar coordinate system with the sun at the origin and the Earth on the
positive x axis on the launch date").

### Naming convention (nPr)
- `1U0` = single cycler n=1 unique-period r=0.
- `1L1` (Aldrin cycler) = n=1 long-period r=1, ref 14 footnote (Aldrin 1985).
- `1S1, 1L2, 1S2, 1L3, 1S3` = the other six n=1 cyclers.
- `2L3` = Case 1 cycler analyzed by Byrnes et al. ref 34 (Byrnes 2002
  AIAA-2002-4423).
- n=7 family uses **`7(R_p)r±`** where R_p is perihelion radius in AU and ±
  indicates ω sign.

### Table 4 (p.625) — most-promising cyclers (1 ≤ n ≤ 6)

| Cycler nPr | Aphelion AU | V∞ Earth km/s | V∞ Mars km/s | Shortest transfer days | Req turn deg | Max turn deg | Notes |
|---|---|---|---|---|---|---|---|
| **1L1**ᵃ | 2.23 | 6.54 | 9.75 | 146 | 84 | 72 | Aldrin (ref 14) |
| 2L2 | 2.33 | 10.06 | 11.27 | 158 | 134 | 44 | |
| **2L3**ᵇ | 1.51ᶜ | 5.65 | 3.05ᵈ | 280ᵉ | 135 | 82 | Case 1 cycler (Byrnes et al. ref 34); ᶜ Mars sma = 1.52 AU; ᵈ "Difference between Mars's speed and spacecraft aphelion speed"; ᵉ "Time to transfer from Earth to aphelion" |
| 3L4 | 1.89 | 11.78 | 9.68 | 189 | 167 | 35 | |
| 3L5 | 1.45ᶜ | 7.61 | 2.97ᵈ | 274ᵉ | 167 | 62 | |
| 3S5 | 1.52ᶜ | 12.27 | 5.45ᵈ | 134ᵉ | 167 | 33 | |
| 4S5 | 1.82 | 11.23 | 8.89 | 88 | 167 | 38 | |
| 4S6 | 1.53 | 8.51 | 4.07 | 157 | 89 | 54 | |
| 5S4 | 2.49 | 10.62 | 12.05 | 75 | 134 | 41 | |
| 5S5 | 2.09 | 9.08 | 9.87 | 89 | 134 | 50 | |
| 5S6 | 1.79 | 7.51 | 7.32 | 111 | 135 | 62 | |
| 5S7 | 1.54 | 5.86 | 3.67 | 170 | 135 | 79 | |
| 5S8 | 1.34ᶜ | 4.11 | 0.71ᵈ | 167ᵉ | 136 | 103 | |
| 6S4 | 2.81 | 7.93 | 12.05 | 87 | 83 | 59 | |
| 6S5 | 2.37 | 6.94 | 10.44 | 97 | 84 | 68 | |
| 6S6 | 2.04 | 5.96 | 8.69 | 111 | 84 | 78 | |
| **6S7**ᶠ | 1.78 | 4.99 | 6.66 | 133 | 85ᶠ | 90ᶠ | **BALLISTIC** (req < max) |
| **6S8**ᶠ | 1.57 | 4.02 | 3.90 | 179 | 85ᶠ | 104ᶠ | **BALLISTIC** |
| **6S9**ᶠ | 1.40ᶜ | 3.04 | 1.21ᵈ | 203ᵉ | 86ᶠ | 120ᶠ | **BALLISTIC** but aphelion 1.40 AU < Mars 1.52 AU |

ᶠ "Ballistic cycler: required turn angle is less than maximum possible turn angle."

### Table 5 (p.626) — n=7 cyclers (all ballistic, repeat every 15 years)

These are special because `n=7` means the line of apsides doesn't need rotating
(ΔΨ = (n/7)·360° mod 360° = 0). All are ballistic. Period 15/r years.

| r | Period yr | Aphelion AU | Years between Earth | Years between Mars |
|---|---|---|---|---|
| 1 | 15 | (11.16, 12.16) | 15 | 15 |
| 2 | 7.5 | (6.66, 7.66) | 15 | 7.5 |
| 3 | 5 | (4.85, 5.85) | 5 | 15 |
| 4 | 3.75 | (3.83, 4.83) | 15 | 3.75 |
| 5 | 3 | (3.16, 4.16) | 3 | 15 |
| 6 | 2.5 | (2.68, 3.68) | 5 | 7.5 |
| 7 | 2.143 | (2.32, 3.32) | 15 | 15 |
| 8 | 1.875 | (2.04, 3.04) | 15 | 1.875 |
| 9 | 1.667 | (1.81, 2.81) | 5 | 15 |
| **10**ᵇ | 1.5 | (1.62, 2.62) | 15 | 7.5 | (**VISIT 2 cycler**, Niehoff ref 18) |
| 11 | 1.364 | (1.46, 2.46) | 15 | 15 |
| **12**ᶜ | 1.25 | (1.32, 2.32) | 5 | 3.75 | (**VISIT 1 cycler**, Niehoff ref 18) |
| 13 | 1.154 | (1.20, 2.20) | 15 | 15 |
| 14 | 1.071 | (1.09, 2.09) | 15 | 7.5 |

(Aphelion radius "range corresponds to perihelion range R_p ∈ (0,1) AU.")

### Table 6 (p.627) — Outbound ballistic S1L1 cycler itinerary, DE405

**This is V1-grade reproducible data** — 22 encounter dates, V∞ km/s,
closest-approach altitude km, leg-duration days, in the DE405 ephemeris model.

| Enc | Date | V∞ km/s | Closest approach alt km | Leg duration days |
|---|---|---|---|---|
| Earth 1 | 9 June 2008 | 6.89 (launch) | — | 541 |
| Earth 2 | 3 Dec 2009 | 6.90 | 31,114 | 186 |
| Mars 3 | 8 June 2010 | 4.31 | 17,704 | 809 |
| Earth 4 | 24 Aug 2012 | 6.42 | 26,490 | 540 |
| Earth 5 | 14 Feb 2014 | 6.43 | 41,524 | 139 |
| Mars 6 | 3 July 2014 | 7.14 | 12,179 | 890 |
| Earth 7 | 9 Dec 2016 | 4.01 | 27,726 | 530 |
| Earth 8 | 22 May 2018 | 4.03 | 19,923 | 115 |
| Mars 9 | 15 Sept 2018 | 6.47 | 11,570 | 934 |
| Earth 10 | 6 April 2021 | 4.61 | 22,992 | 532 |
| Earth 11 | 20 Sept 2022 | 4.59 | 14,780 | 223 |
| Mars 12 | 1 May 2023 | 2.77 | 7,593 | 793 |
| Earth 13 | 2 July 2025 | 7.08 | 23,858 | 542 |
| Earth 14 | 26 Dec 2026 | 7.09 | 35,164 | 170 |
| Mars 15 | 14 June 2027 | 5.26 | 13,751 | 830 |
| Earth 16 | 21 Sept 2029 | 5.78 | 26,818 | 537 |
| Earth 17 | 12 March 2031 | 5.78 | 39,044 | 125 |
| Mars 18 | 15 July 2031 | 7.70 | 10,566 | 915 |
| Earth 19 | 15 Jan 2034 | 3.78 | 22,938 | 529 |
| Earth 20 | 20 June 2035 | 3.76 | 9,586 | 138 |
| Mars 21 | 13 Nov 2035 | 4.68 | 15,525 | 907 |
| Earth 22 | 7 May 2038 | 5.55 | — | — |

p.627 text: "We note that the flyby V∞ are all less than 7.7 km/s, the flyby
altitudes are all greater than 7500 km, and the Earth-Mars legs range from
115 to 223 days."

**S1L1 footnote on Table 6**: "Outbound ballistic S1L1 cycler itinerary (using
DE405 ephemerides of Earth and Mars). The positions and velocities of Earth
and Mars were determined using the Jet Propulsion Laboratory's DE405
ephemerides. We note that the flyby V∞ are all less than 7.7 km/s, the flyby
altitudes are all greater than 7500 km, and the Earth-Mars legs range from
115 to 223 days." (p.627 text continuation; the **S1L1 ballistic cycler**
identification cites **ref 36 = McConaghy/Longuski/Byrnes AIAA Paper 2002-4420**).

### Numerical results in body text (p.627)

- Aldrin (1L1) required turn = 84°, available turn at 200 km Earth flyby alt =
  72°. ΔV required. Required perihelion at Earth alt ≈ −1731 km (i.e. below
  surface).
- "Total ΔV required by the outbound cycler is 1.73 km/s (or 247 m/s per
  synodic period, on average), and the total ΔV required by the inbound cycler
  is 2.04 km/s (or 291 m/s per synodic period, on average)." ← **EXACT MATCH
  with Byrnes 1993 JSR Tables 1-2** (1.73 km/s outbound, 2.04 km/s inbound).
- "The flight times between Earth and Mars vary between 147 and 170 days." ←
  Match with Byrnes 1993.
- "The V∞ at Earth and Mars vary from 5.39 to 6.19 km/s, and from 6.05 to
  11.74 km/s, respectively." ← Exact match with Byrnes 1993.

### References

- Ref 13 = Aldrin 1985 SAIC presentation (same source as Byrnes 1993 ref 8
  and Nock-Friedlander 1987 ref 14).
- Ref 14 = Byrnes-Longuski-Aldrin 1993 JSR DOI 10.2514/3.25519 (Byrnes 1993 —
  the paper I just digested).
- Ref 18 = Niehoff (VISIT orbits).
- Ref 34 = Byrnes 2002 AIAA-2002-4423 (Case 1 / 2L3 cycler).
- Ref 35 = Niehoff-Friedlander-McAdams 1991 IAF Paper 91-438.
- Ref 36 = **McConaghy/Longuski/Byrnes "Analysis of a Broad Class of
  Earth-Mars Cycler Trajectories," AIAA Paper 2002-4420 (Aug 2002)** ← S1L1
  primary source.

## 4. `KNOWN_CORPUS` impact (RECOMMEND, do not edit)

The current SnLm anchor (literature_check.py line 263, 269, 276):

```python
authors=("Russell", "Ocampo", "McConaghy", "Landau", "Longuski"),
citation="Russell & Ocampo, J. Spacecraft & Rockets 41(1) 2004; "
"McConaghy et al., J. Spacecraft & Rockets 43(2) 2006",
```

**MISSING CITATION**: McConaghy-Longuski-Byrnes 2004 JSR 41(4) DOI
10.2514/1.11939 — THE paper that introduces the **nPr nomenclature** AND
publishes the S1L1 cycler with DE405 itinerary. This paper PREDATES the SnLm
nomenclature (which appears in McConaghy-Russell-Longuski 2005 JSR 42(4)).
The anchor should cite ALL THREE papers in the McConaghy lineage:

```
McConaghy-Longuski-Byrnes 2004 JSR 41(4) (nPr nomenclature, S1L1 DE405 itinerary)
McConaghy-Russell-Longuski 2005 JSR 42(4) (SnLm standard nomenclature)
McConaghy-Landau-Yam-Longuski 2006 JSR 43(2) (notable 2-synodic S2L1)
```

Plus Russell-Ocampo 2003/2004 (which is JGCD 27(3), NOT JSR 41(1), as
established in the Russell-Ocampo 2003 verdict note).

Plus the McConaghy 2004 PhD dissertation (already cited at line 705).

Recommended consolidated citation text:

```
"Russell & Ocampo, JGCD 27(3):321-335 (2004) DOI 10.2514/1.1909, preprint "
"AAS-03-145; McConaghy, Longuski & Byrnes, JSR 41(4):622-628 (2004) DOI "
"10.2514/1.11939 (nPr family, S1L1 DE405 itinerary Table 6); "
"McConaghy, Russell & Longuski, JSR 42(4) (2005) DOI 10.2514/1.8123 (SnLm "
"standard nomenclature); McConaghy, Landau, Yam & Longuski, JSR 43(2) "
"(2006) DOI 10.2514/1.15215; Russell PhD diss (UT Austin 2004); McConaghy "
"PhD diss (Purdue 2004)"
```

## 5. Catalogue impact (RECOMMEND)

### Rows touched

- **`aldrin-cycler`** — McConaghy 2004 JSR p.627 numerically corroborates
  Byrnes 1993 exactly (1.73 km/s outbound, 2.04 km/s inbound, transit
  147-170 d, V∞ E 5.39-6.19, V∞ M 6.05-11.74). **Add to corroborating_sources.**
- **Any catalogue row with `mcconaghy-2002` nPr-style ID** (Cycler 6S6, 6S7,
  6S8, 6S9 are the McConaghy 2002 cyclers Russell-Ocampo 2003 Table A2 cites).
  Check if these are in catalogue under McConaghy IDs.
- **`mcconaghy-2006-em-k2`** (line 242 in catalogue) — McConaghy 2006 notable
  S2 cycler is the 2P2 cycler in the McConaghy 2004 nPr naming = **2L3** in
  this paper (per Table 4 footnote b: "Case 1 cycler analyzed by Byrnes et
  al."), with V∞ E=5.65, V∞ M=3.05, Aphelion=1.51 AU, transfer 280 d. **NOT
  the same as the catalogue's 4.7 / 5.0 km/s value.** ← **POSSIBLE
  DISCREPANCY** — let me re-check.

Wait — the catalogue's `mcconaghy-2006-em-k2` cites V∞ E = 4.7 km/s, V∞ M =
5.0 km/s (from McConaghy 2006 abstract). McConaghy 2004 Table 4 row 2L3 shows
V∞ E = 5.65, V∞ M = 3.05ᵈ (footnote d says "Difference between Mars's speed
and spacecraft aphelion speed" — so V∞ M is NOT 3.05; it's the speed-difference
at aphelion since the cycler aphelion doesn't reach Mars). The 4.7 / 5.0
McConaghy 2006 values may refer to a **different S2 cycler** — McConaghy 2006
is "Notable two-synodic-period Earth-Mars cycler" — I'll verify when I read
paper 6.

### V0→V1 promotion candidacy

**S1L1 cycler (Table 6, p.627) is a VERY STRONG V1 candidate.** Table 6 gives:
- Reference epoch (Earth 1 launch on 9 June 2008).
- 22-encounter DE405 itinerary with date, V∞, closest-approach altitude, leg
  duration.
- Stated bounds: V∞ < 7.7 km/s for all flybys, all altitudes > 7,500 km,
  Earth-Mars legs 115-223 days.

This is V1-grade reproducible data in a published ephemeris (DE405). Any
catalogue row for **S1L1** could be promoted to V1 by demonstrating that our
M7 stack produces the same 22-encounter sequence within tolerance.

**Catalogue check needed**: is there an `s1l1` row in catalogue.yaml?

(From memory note `project_s1l1_realeph_closure_blocker.md`: S1L1 is the
saga's central study; "modelling stack complete; S1L1 fails on family-selection
(off-basin), not infra." So our M7 has NOT successfully reproduced S1L1 in
real ephemeris. **The McConaghy 2004 Table 6 itinerary is the V1 GROUND TRUTH
we should compare against** if we ever do succeed.)

Other Table 4 rows that are BALLISTIC (req < max turn):
- **1L1** (Aldrin): req 84 > max 72 → NOT ballistic, as expected.
- **6S7**: req 85 < max 90, V∞_E=4.99, V∞_M=6.66, aphelion 1.78 AU. BALLISTIC.
- **6S8**: req 85 < max 104, V∞_E=4.02, V∞_M=3.90, aphelion 1.57 AU. BALLISTIC.
- **6S9**: req 86 < max 120, V∞_E=3.04, V∞_M=1.21ᵈ, aphelion 1.40 AU but
  BELOW Mars orbit (1.52). BALLISTIC but unreachable Mars in this model.

These three (6S7, 6S8, 6S9) match the McConaghy 2002 cyclers Russell-Ocampo
2003 Table A2 footnotes a/b/c/d as Cycler 6S9/6S8/6S7/6S6 — **same family**.

### Other actionable findings

The **n=7 family (Table 5)** is all ballistic and includes the VISIT 1 cycler
(7(R_p)12) and VISIT 2 cycler (7(R_p)10). These appear in the catalogue's
`niehoff-visit-*` rows (if any). Worth a check.

## 6. Errata / surprises

1. **Two distinct McConaghy 2004 works exist**: this JSR paper AND McConaghy's
   PhD dissertation (Purdue 2004). Both date from the same year; the catalogue
   currently distinguishes them (line 705 dissertation, this JSR paper not
   explicitly cited).
2. **nPr ≠ SnLm**. Both nomenclatures cover the same trajectory families but
   were introduced in different papers (McConaghy 2004 JSR for nPr; McConaghy-
   Russell-Longuski 2005 JSR for SnLm). The catalogue uses SnLm-style IDs but
   the **earliest published cycler names (Aldrin = 1L1, Case 1 = 2L3, VISIT 1
   = 7(R_p)12) are from nPr**. McConaghy 2005 introduces the SnLm renaming.
3. **Table 4 footnote d "V∞ at Mars" is special**: for cyclers whose aphelion
   doesn't reach Mars (2L3, 3L5, 3S5, 5S8, 6S9), the value is the DIFFERENCE
   between Mars's speed and the cycler's speed at aphelion — i.e. the
   minimum-energy Δv to reach Mars from aphelion. This is NOT a true V∞.
   Easy to misread.
4. **Table 4 footnote e "Shortest transfer time"** for the same cyclers is the
   transfer time from Earth to aphelion (NOT Earth to Mars), since Mars is
   never reached. Easy to misread.
5. **S1L1 ballistic cycler (Table 6, DE405)** uses the SnLm name `S1L1` —
   confusingly, this paper introduces the nPr names but Table 6 uses the SnLm
   name "S1L1". This is a forward-reference to McConaghy-Russell-Longuski
   2005's renaming convention. The author footnote points to **ref 36 = AIAA
   2002-4420** as the S1L1 source paper.
6. **Table 6 starts launch on 9 June 2008** — fixed epoch. Any V1 reproduction
   MUST use this exact launch date (or document the shift if using a different
   epoch).
7. **Table 6 column "Leg duration"** appears to be the duration of the leg
   ENDING at that encounter — confirmed by Earth 1 → Earth 2 dates (9 June
   2008 → 3 Dec 2009 = ~542 days; Table shows 541 days at Earth 1, so leg
   duration is for the NEXT leg). Actually re-checking: Earth 1 → Earth 2 =
   541 days outbound; Earth 2 → Mars 3 = 186 days; this matches Table 6.
   **Leg duration is for the leg STARTING at that encounter**.

## 7. Action items for parent

- [ ] **Fix `KNOWN_CORPUS` SnLm anchor**: add McConaghy 2004 JSR (DOI
  10.2514/1.11939) as a primary citation alongside the existing 2005/2006
  papers. Also fix the Russell-Ocampo JSR→JGCD error per the Russell-Ocampo
  2003 verdict note (one combined fix).
- [ ] **V1 promotion candidate for S1L1**: McConaghy 2004 Table 6 gives a
  full 22-encounter DE405 itinerary starting 9 June 2008. **This is the V1
  ground truth.** If our M7 / ephemeris stack ever closes S1L1 (currently
  blocked per `project_s1l1_realeph_closure_blocker`), the comparison metric
  is per-encounter (date, V∞, alt) vs Table 6. Track separately.
- [ ] **Cross-check `mcconaghy-2006-em-k2` row** against McConaghy 2004
  Table 4 2L3 entry. The 2L3 entry has Aphelion=1.51 AU (Mars orbit 1.52 — so
  aphelion just inside Mars), V∞ E = 5.65 km/s, V∞ M = 3.05 km/s (NOT a true
  V∞, but the Mars-aphelion speed difference). McConaghy 2006 publishes 4.7
  km/s @ Earth + 5.0 km/s @ Mars for the "notable" S2L1 cycler. **Probably a
  different cycler** — McConaghy 2006 is about S2L1 (= 2P? in nPr), not 2L3.
  Verify in paper 6 digest.
- [ ] **Add a catalogue row note to `aldrin-cycler` corroborating_sources**
  citing McConaghy 2004 JSR p.627 as a numerical re-statement of Byrnes 1993
  Tables 1-2.
- [ ] **n=7 family worth a separate catalogue audit**: Table 5 gives 14
  all-ballistic 15-year cyclers including VISIT 1, VISIT 2. Each has an
  aphelion radius RANGE (corresponds to perihelion range R_p ∈ (0,1) AU).
  Catalogue's `niehoff-visit-*` rows (if any) should reference Table 5.
- [ ] **Track `AIAA 2002-4420 (McConaghy/Longuski/Byrnes)` as a still-
  outstanding acquisition** if not already in our corpus — it's ref 36 here,
  the S1L1 cycler's primary source. The catalogue cites it (line 297) but
  the paper itself may not be in our `papers/` directory.
