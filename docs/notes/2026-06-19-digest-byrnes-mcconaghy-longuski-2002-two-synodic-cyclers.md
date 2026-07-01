# Digest — Byrnes, McConaghy & Longuski 2002 "Analysis of Various Two Synodic Period Earth-Mars Cycler Trajectories"

**Date**: 2026-06-19 AET (Agent dispatched by parent #384; corpus-document-policy digest of 3 newly-filed papers)
**Verdict (TL;DR)**: **CROSS-VALIDATION / PRECEDENCE SOURCE for the S1L1 row's V_inf — but with an important caveat.** This is the **Byrnes-led sibling** to the already-digested McConaghy-Longuski-Byrnes "Analysis of a *Broad Class*" paper (AIAA 2002-4420, `2026-06-17-digest-mcconaghy-2002.md`). Both were presented at the **same conference** (AIAA/AAS Astrodynamics Specialist Conference, Monterey CA, Aug 2002). This paper's **Case 3** ("two synodic period cycler with backflip plus 1-year loop") is explicitly identified by the authors as topologically near-identical to the **S1L1-B cycler** of the broad-class companion paper (p.6: "remarkable similarity to Case 3"). The paper publishes a **full 30-year real-ephemeris V_inf-matched table** for Case 3 (the Table, p.6) and circular-coplanar V_inf for Cases 1/2/3. **CRITICAL FINDING:** this paper does NOT publish the literal S1L1-B V_inf pair; it publishes Case 3 (a 3-Earth-flyby *resonant* variant), whose real-ephemeris V_inf at Earth oscillates ~4–7.5 km/s and at Mars ~3–8 km/s (p.7). So it CORROBORATES the catalogue spec-9 anchor values (V_inf_E 5.65 / V_inf_M 3.05) as lying *inside the observed real-world oscillation envelope*, and supplies an INDEPENDENT (Byrnes-authored) source predating McConaghy 2004 JSR — but it is NOT a clean single-pair literal match to S1L1-B, because Case 3 ≠ S1L1-B (the difference is exactly the one the authors call out: Case 3 has 3 Earth flybys and a *resonant* Earth-Earth transfer, S1L1-B has 2 Earth flybys and a *non-resonant* one, p.6).

---

## 1. Header

- **Title** (verbatim, p.1): "Analysis of Various Two Synodic Period Earth-Mars Cycler Trajectories"
- **Authors** (verbatim, p.1 byline + footnotes):
  - **Dennis V. Byrnes** — Principal Engineer, Deputy Section Manager, JPL (Senior Member AIAA, Member AAS)
  - **T. Troy McConaghy** — Graduate Student, Purdue University (Student Member AIAA, Member AAS)
  - **James M. Longuski** — Professor, Purdue (Associate Fellow AIAA, Member AAS)
- **Affiliations**: Jet Propulsion Laboratory, California Institute of Technology, Pasadena, California 91109-8099; Purdue University, West Lafayette IN 47907-1282
- **Venue**: AIAA/AAS Astrodynamics Specialist Conference, Monterey CA, August 2002 (per title + the paper's own ref 4 citing the companion as "AIAA/AAS Astrodynamics Specialist Conference, Monterey, CA, Aug. 2002")
- **Source**: AIAA/AAS Astrodynamics Specialist Conf., Monterey CA, Aug 2002
- **Page count**: 8 pages, text layer present

### Relationship to the companion paper (NOT the same paper)

This paper and `2026-06-17-digest-mcconaghy-2002.md` are **two distinct papers, same conference, overlapping author set, complementary scope**:

| | This paper (Byrnes-led) | Companion (McConaghy-led) |
|---|---|---|
| Title | "Analysis of **Various Two Synodic Period** Earth-Mars Cycler Trajectories" | "Analysis of a **Broad Class** of Earth-Mars Cycler Trajectories" |
| Lead author | **Byrnes** | **McConaghy** |
| Paper ID | not printed on this scan | AIAA 2002-4420 |
| S1L1 treatment | refers to it via the companion ("S1L1-B cycler", p.6) | **first publishes** S1L1 (Table 6, Fig 8: V_inf_E 4.7, V_inf_M 5.0) |
| Method | 3 specific 2-synodic cases, circular-coplanar + real-ephemeris V_inf-matching | systematic nPr enumeration |

This paper's **ref 4** (p.8) IS the companion: *"McConaghy, T. T., Longuski, J. M., and Byrnes, D. V., 'Analysis of a Broad Class of Earth-Mars Cycler Trajectories,' AIAA/AAS Astrodynamics Specialist Conference, Monterey, CA, Aug. 2002."*

## 2. What the paper actually is

A focused study of **three** specific two-Earth-Mars-synodic-period cycler trajectories (plus brief mention of a fourth, worse, variant), analyzed first in a circular-coplanar model and then (Case 3 only) with real Earth/Mars ephemerides via a Sun-centered point-to-point conic, V_inf-matching ("instantaneous V_inf rotation at flybys") model. Stated goal (Introduction, p.1): improve on the Aldrin Cycler, which makes 3 2/7 revolutions in 4 2/7 years, by finding 2-synodic cyclers with lower V_inf at both Earth and Mars (closer to Hohmann conditions).

Two-synodic cyclers all require **four vehicles** (p.1): two on "Up" trajectories (Earth→Mars Type I, ≤ ~9 months) and two on "Down" trajectories (Mars→Earth Type I, short). Contrast: the Aldrin Cycler needs only two vehicles (p.1).

### Circular-coplanar assumptions (p.2)
1. Earth-Mars synodic period = 2 1/7 years (⇒ assumes Mars period = 1 7/8 yr; true value 1.881 yr noted as discrepancy). Geometry repeats every **15 years**.
2. Earth/Mars/cycler all coplanar (ecliptic).
3. Circular planet orbits.
4. Conic, prograde cycler.
5. Only Earth provides gravity assist.
6. Instantaneous flybys.
(Real Mars used later: a = 1.524 AU, perihelion 1.381 AU, aphelion 1.666 AU — p.2.)

## 3. Extracted cycler data (circular-coplanar model)

### Case 1 — Simple two synodic period cycler (Fig 1, p.2-3)
- Period **P = 1.348 yr**
- Aphelion radius **R_A = 1.51 AU**
- **V_inf at Earth = 5.6 km/s** (text, p.2); turn-angle discussion gives V_inf 5.65 km/s, max Earth rotation only ~82° vs ~135° required (p.2)
- Perihelion ~0.93 AU (p.2)
- Earth→Mars transfer ~225°, "a little over nine months", Type I (Up) / Type VI (Mars-Earth leg) (p.2)
- Makes **3 2/7 revs in 4 2/7 yr** total cycle; encounters Earth 2/7 rev (102.9°) from start (p.2)
- DEFICIENCY in real world: does not reach Mars' mean distance; would need ΔV; reaches Mars at most half the time over 7 cycles (p.2). Stated value: "real basis for variations."

### Case 2 — Two synodic period cycler with "backflip" (Fig 2, p.3-4)
- Adds a 2nd Earth flyby ~6 months / 180° after the first (Uphoff "backflip", ref 3)
- Earth-Mars-Earth leg: **P = 1.325 yr**, aphelion **R_A = 1.45 AU**, **V_inf at Earth = 4.15 km/s** (p.3)
- Earth→Mars Type I or II, Mars→Earth Type V; makes 2 11/14 revs in 3 11/14 yr; backflip year-period orbit re-encounters Earth ~6 mo/180° later completing 3 2/7 revs in 4 2/7 yr (p.3)
- Lower V_inf enables Earth to rotate V_inf up to ~102° (p.4)
- DEFICIENCY: does not reach Mars in circular-coplanar model (only does when Mars near perihelion in real world); reaches Mars at most 2 of 7 cycles without ΔV (p.4)

### Case 3 — Two synodic period cycler with "backflip" PLUS 1-year loop (Fig 3, p.4-5) — **the catalogue-relevant variant**
- Adds a **third** Earth flyby: a 1-year Earth-Earth loop in addition to the backflip
- Earth→Mars Type I, Mars→Earth Type III or IV; first Earth encounter at 1 11/14 revs in 2 11/14 yr (p.4)
- Earth-Mars-Earth leg: **P = 1.484 yr**, aphelion **R_A = 1.65 AU**, **V_inf at Earth = 5.4 km/s** (p.5)
- "aphelion approximately equal to Mars' aphelion ⇒ **always crosses Mars orbit in the real world**" (p.5) — this is why Case 3 is the only fully real-world-viable case
- **Identified by authors as ≈ S1L1-B** (see §4)

### Fourth variant (one or two 1-year loops without backflip) — p.5
- Possible but "much higher V_inf's, less desirable than Cases 1, 2, 3, or the Aldrin Cycler." No tabulated numbers.

## 4. THE S1L1 CROSS-REFERENCE (p.6, "Ongoing and Future Work") — verbatim-anchored

> "The work of McConaghy, et. al. [ref 4], presented at this Conference as well, has identified a cycler denoted as the **S1L1-B cycler**. Studying Figure 8 given in that paper shows remarkable similarity to **Case 3** in Figure 3 above. In fact they share significant similarities and **one important difference. The Earth-Mars-Earth legs are nearly identical, however the S1L1 cycler has only two Earth flybys instead of three and the Earth-Earth transfer is not a resonance.** That is, the transfer time and angle between the two Earth flybys is not an exact multiple or half-multiple of the Earth's period. This gives significantly more flexibility ... A preliminary study underway has already identified a **completely ballistic version of the S1L1-B cycler over a 30-year period.**" (p.6)

**Interpretation for the catalogue:**
- The literal S1L1-B V_inf pair (E 4.7 / M 5.0 km/s) is published in the **companion** paper, not this one.
- THIS paper publishes **Case 3** (3 Earth flybys, resonant E-E transfer), which is a *near-twin* of S1L1-B but NOT identical.
- The paper explicitly states a fully ballistic S1L1-B exists over 30 years (a precedence claim for the catalogue's "ballistic two-synodic-period" designation), but defers the per-member S1L1-B data to "a future paper" (p.6). So **this paper is a precedence/topology corroborator for S1L1-B's existence and ballistic-ness, not an independent numeric V_inf source for the literal S1L1-B pair.**

## 5. Case 3 real-ephemeris V_inf table (the Table, p.6) — INDEPENDENT BYRNES DATA, sourced-only

"Detailed Analysis of Case 3" (p.6) gives a full 7-cycle (30-year) V_inf-matched conic table. Reproduced verbatim (V-inf in km/s; Altitude in km, Mars flybys only):

| Planet | Date | V-inf (km/s) | Altitude (km) |
|---|---|---|---|
| Earth | 8/25/05 22h | 4.362 | |
| Mars  | 2/21/06 3h  | 3.238 | 8000 |
| Earth | 6/5/08 19h  | 7.234 | |
| Earth | 12/5/09 16h | 7.234 | |
| Mars  | 6/2/10 12h  | 4.534 | 19000 |
| Earth | 8/19/12 10h | 6.945 | |
| Earth | 2/18/14 7h  | 6.945 | |
| Mars  | 7/1/14 8h   | 7.361 | 8000 |
| Earth | 11/30/16 15h| 4.622 | |
| Earth | 6/1/18 12h  | 4.622 | |
| Mars  | 9/15/18 10h | 6.675 | 7000 |
| Earth | 3/30/21 23h | 5.121 | |
| Earth | 9/29/22 20h | 5.121 | |
| Mars  | 4/21/23 16h | 2.973 | 10000 |
| Earth | 6/28/25 17h | 7.436 | |
| Earth | 12/28/26 14h| 7.436 | |
| Mars  | 6/10/27 21h | 5.494 | 16000 |
| Earth | 9/15/29 1h  | 6.387 | |
| Earth | 3/16/31 22h | 6.387 | |
| Mars  | 7/15/31 2h  | 7.827 | 7000 |
| Earth | 1/8/34 23h  | 4.151 | |
| Earth | 7/10/35 20h | 4.151 | |
| Mars  | 11/13/35 17h| 4.897 | 9000 |
| Earth | 5/3/38 0h   | 6.12  | |

**Key sourced ranges (p.7):**
- Mars V_inf: "vary between about **3 km/s and 8 km/s**" (table min 2.973, max 7.827); circular-coplanar Case value would be 5.3 km/s.
- Earth V_inf: "vary between about **4 km/s and 7.5 km/s**" (table min 4.151, max 7.436); circular-coplanar 5.4 km/s.
- "Earth-Earth transfers constrained to be exactly **1.5 years** apart" (resonant — the defining Case 3 ≠ S1L1-B distinction).
- Mars flybys all at reasonably high altitude (7000-19000 km); gravity assist controls heliocentric inclination + small energy adjust (p.6-7).

**Note on table fidelity**: the paper calls this a V_inf-matched point-to-point conic with instantaneous flyby rotations, "sufficient accuracy for developing long term trajectory scenarios that can be closely reproduced with fully numerically integrated trajectory models" (p.6). It is NOT a ballistic-closed integrated solution; "Some initial optimal trajectory simulation ... has begun" and "adjustments to the dates ... will be necessary to minimize required deterministic ΔV's. It may be possible, although it has not yet been demonstrated, that a completely ballistic trajectory may be attainable" (p.7). So the Case 3 table is a **V0-grade design table, not a ballistic-verified closure.**

## 6. References in this paper (p.8)
1. Aldrin, B., "Cyclic Trajectory Concepts," SAIC presentation to the Interplanetary Rapid Transit Study Meeting, JPL, Oct. 28, 1985.
2. Byrnes, D.V., Longuski, J.M., and Aldrin, B., "Cycler Orbit Between Earth and Mars," J. Spacecraft & Rockets, Vol. 30, No. 3, May-June 1993, pp. 334-336.  *(already digested: `2026-06-17-digest-byrnes-longuski-aldrin-1993.md`)*
3. Uphoff [incomplete citation in original — "backflip" lunar trajectory concept]
4. McConaghy, T. T., Longuski, J. M., and Byrnes, D. V., "Analysis of a Broad Class of Earth-Mars Cycler Trajectories," AIAA/AAS Astrodynamics Specialist Conference, Monterey, CA, Aug. 2002.  *(already digested: `2026-06-17-digest-mcconaghy-2002.md`)*

## 7. Catalogue / KNOWN_CORPUS relevance

### 7.1 S1L1 row (`s1l1-2syn-em-cpom`, catalogue.yaml L484) — RE-ANCHOR ANALYSIS

The row currently has `vinf_source: mcconaghy-2002` with V_inf_E 5.65 / V_inf_M 3.05 traced ONLY to spec.md §9 (data L511, L513 carry explicit provenance caveats — the values could not be corroborated in any mined source). The companion `mcconaghy-2002` digest established that the **literal S1L1-B pair is E 4.7 / M 5.0 km/s** (broad-class Table 6/Fig 8), which already CONTRADICTS the spec-9 5.65/3.05 anchor.

**What THIS paper adds:**
1. **Precedence corroboration** that S1L1-B is a real, ballistic, 30-year two-synodic cycler (Byrnes' own statement, p.6) — strengthens the existence/ballistic claim independent of McConaghy's lead authorship.
2. **An independent real-ephemeris V_inf envelope** for the near-twin Case 3: Earth ~4–7.5 km/s, Mars ~3–8 km/s (p.7). The spec-9 anchor (E 5.65 / M 3.05) and the companion's S1L1-B pair (E 4.7 / M 5.0) **both fall inside this envelope** — so this paper is best used as an *envelope-consistency* cross-check, not a single-pair literal source.
3. **A defect flag**: it shows the spec-9 3.05 km/s "Mars" value is consistent with the *low end* of the Case 3 Mars oscillation, but Case 3 ≠ S1L1-B, and the companion S1L1-B literal Mars value is 5.0. The 5.65/3.05 spec-9 pair does not match either authored source's literal S1L1-B numbers.

### 7.2 RE-ANCHOR RECOMMENDATION (for the parent — DO NOT WRITE here)

- **DO add this paper to the s1l1 citation chain as a corroborating precedence + V_inf-envelope source** (suggest a new provenance key e.g. `byrnes-mcconaghy-2002` distinct from the existing `mcconaghy-2002`, since they are different papers). It is the Byrnes-authored, conference-precedence companion and the only source publishing a per-encounter real-ephemeris Case 3 table.
- **It is NOT a clean single-pair cross-validation source for the literal S1L1-B V_inf.** The literal S1L1-B pair (E 4.7 / M 5.0) remains the companion broad-class paper's Table 6/Fig 8. This paper publishes Case 3 (the resonant near-twin), not S1L1-B per se.
- **Tier guidance**: the catalogue's spec-9 5.65/3.05 anchor remains *uncorroborated as a literal pair* by either 2002 paper. The honest move is either (a) re-anchor s1l1's literal V_inf to the companion's authored S1L1-B pair (E 4.7 / M 5.0) and cite THIS paper as the envelope/precedence corroborator, demoting the spec-9 values to a documented historical discrepancy; or (b) keep the spec-9 values but add this paper's envelope (E 4–7.5 / M 3–8) as the cross-validation bracket showing the spec-9 pair is *envelope-consistent but not literal-matched*. Recommend (a): it eliminates a long-standing uncorroborated literal anchor and uses two authored 2002 sources. **This is a V0→cross-validated (envelope) move, not a clean literal V0→V1, because no source publishes the exact 5.65/3.05 pair.**
- **Caution**: the Case 3 table is a V_inf-matched design table, NOT a ballistic-verified closure (§5). Do not upgrade s1l1's `trajectory_regime` confidence on the basis of this paper; the "completely ballistic version" is asserted but its data is deferred to an unpublished future paper.

### 7.3 KNOWN_CORPUS
This paper should appear in the literature corpus as a distinct entry (different title/lead author from AIAA 2002-4420). For the literature-novelty gate it confirms the two-synodic-period Earth-Mars cycler family is published prior art (Cases 1/2/3 + S1L1-B), reinforcing that no two-synodic E-M cycler can be claimed novel.

### 7.4 Establishes nothing new for the corrector/optimizer lane
Methodology here is V_inf-matching point-to-point conic + instantaneous-flyby rotation — design-grade, not a BVP corrector. (The corrector/optimizer provenance is the D'Amario-Byrnes-Stanford 1981 paper, digested separately.)
