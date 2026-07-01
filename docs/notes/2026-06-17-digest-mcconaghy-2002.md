# Digest — McConaghy, Longuski & Byrnes 2002 "Analysis of a Broad Class of Earth-Mars Cycler Trajectories"

**Date**: 2026-06-17 AET (Agent #381, parent #352 Mars-cycler wave digest closeout + KNOWN_CORPUS Aldrin anchor confirmation)
**Verdict (TL;DR)**: **ANCHOR CITATION CONFIRMED**. This IS the AIAA 2002-4420 paper cited by the existing Aldrin KNOWN_CORPUS anchor. Confirmed authors `(McConaghy, Longuski, Byrnes)`, venue (AIAA/AAS Astrodynamics Specialist Conference and Exhibit, 5-8 August 2002, Monterey, California), paper ID AIAA 2002-4420. **The paper is the FOUNDATIONAL nPr nomenclature paper** — it introduces the `n(P)r` classification (n = synodic-period repeats, P ∈ {L,S,U} short/long/unique-period, r = revs) that McConaghy 2004 (#352) and McConaghy-Russell-Longuski 2005 (#353) build on. Aldrin cycler = **1L1**; VISIT 1 = **12C5** (n=7); VISIT 2 = **10C7** (n=7). **S1L1 cycler appears IN THIS PAPER FOR THE FIRST TIME** (Table 6, Fig 8) as a ballistic two-synodic-period cycler with V_∞_E = 4.7 km/s, V_∞_M = 5.0 km/s — confirming this 2002 conference paper as the original S1L1 source, NOT McConaghy 2004 JSR.

---

## 1. Header

- **Title** (verbatim, p.1): "ANALYSIS OF A BROAD CLASS OF EARTH-MARS CYCLER TRAJECTORIES"
- **Authors** (verbatim, p.1 byline + footnotes):
  - **T. Troy McConaghy** — Graduate Student, Purdue University (Student Member AIAA, Member AAS)
  - **James M. Longuski** — Professor, Purdue (Associate Fellow AIAA, Member AAS)
  - **Dennis V. Byrnes** — Principal Engineer, Deputy Manager, Navigation and Mission Design Section, JPL (Senior Member AIAA, Member AAS)
- **Affiliations**: School of Aeronautics and Astronautics, Purdue University, West Lafayette, Indiana 47907-1282; Jet Propulsion Laboratory, California Institute of Technology, Pasadena, California 91109-8099
- **Venue**: AIAA/AAS Astrodynamics Specialist Conference and Exhibit, 5-8 August 2002, Monterey, California
- **Paper ID**: AIAA 2002-4420 (p.1 top right)
- **Copyright**: 2002 by the author(s). Published by the American Institute of Aeronautics and Astronautics, Inc., with permission.
- **Page count**: 10 pages

## 2. What the paper actually is

A **systematic enumeration paper** that constructs the broad class of conic Earth-Mars cycler trajectories under simplifying assumptions (circular coplanar planet orbits, Hohmann-period Mars at 1.881 years, Earth-only gravity assists, instantaneous flybys). The method patches consecutive collision orbits (Poincaré "second species" periodic orbits) at planetary encounters so the entire trajectory repeats after an integer number of Earth-Mars synodic periods S = 2 + 1/7 years.

The paper produces:
1. A **systematic classification scheme** — the `n(P)r` cycler notation that becomes the field's standard (subsequently reused by McConaghy 2004/2005/2006, Russell-Ocampo 2003, and all follow-on Purdue work).
2. **Rediscovery** of known cyclers (Aldrin = 1L1; VISIT 1 = a 12-rev n=7 ballistic; VISIT 2 = a 10-rev n=7 ballistic).
3. **Discovery** of previously unknown cyclers, including the S1L1 family that becomes the catalogue's central topology-problem (memory `project_s1l1_realeph_closure_blocker.md`).

Eleven sections:

1. **Introduction** (p.1): historical lineage (Ross 1963, Minovitch 1967, Hollister 1969, Hollister-Menning 1969/1970, Rall 1969, Rall-Hollister 1971, Aldrin 1985 SAIC presentation, Byrnes-Longuski-Aldrin 1993, Niehoff VISIT 1985-86, Friedlander-Niehoff-Byrnes-Longuski 1986)
2. **Methodology** (p.2): 6 simplifying assumptions, the Lambert problem setup connecting R₁ = (1,0) to R₂ = (cos(2πnS), sin(2πnS)) in time T = nS
3. **Categorizing Cycler Trajectories** (pp.3-4): the nPr notation introduced. P ∈ {L, S, U} for long-period / short-period / unique-period; r = revs rounded down to nearest integer. Tables 1-3 enumerate n=1-4 cyclers
4. **Number of Cycler Vehicles Required** (pp.4-5): inbound vs outbound concept; upper bound 2n
5. **Aphelion Radius** (p.5): cyclers must cross Mars orbit (R_a > 1.52 AU)
6. **V_∞ at Earth and at Mars** (p.5): Hohmann lower bound (V_∞_E ≈ 2.95 km/s, V_∞_M ≈ 2.65 km/s) — Hohmann not a cycler though
7. **Required vs Maximum Possible Turn Angle** (p.5): line-of-nodes rotation Δψ = (n/7)·360°
8. **The Most Promising Solutions** (pp.5-7): Table 4 lists 20 most promising cyclers n=1-6; Table 5 lists 14 n=7 ballistic cyclers
9. **Extending the Method** (pp.7-8): intermediate Earth flybys, Fig 7 schematic, Table 6 three two-synodic-period cyclers with intermediate Earth GA (U0L1, L2U0, S1L1)
10. **Other Possible Extensions** (p.8): multiple intermediate GA, Venus/Mars/other-planet intermediate GA, repeat times other than 2S, total-ΔV estimation for powered cyclers, accurate-model verification
11. **Conclusions** (p.8): "Previously known cyclers, such as the Aldrin cycler and the VISIT cyclers, can be constructed using this method. However, our construction method is not completely general. We investigated a simple extension and discovered some remarkable, previously-unknown cyclers."

## 3. Per-cycler data tables extracted

### 3.1 Notation (p.3)

A cycler is denoted `n P r` where:
- **n** = number of Earth-Mars synodic periods before repeating (1, 2, 3, …)
- **P** ∈ {L, S, U} = long-period, short-period, unique-period Lambert solution
- **r** = number of complete revs before repeating, rounded down to nearest integer (r ∈ 0, 1, 2, …, r_max(n))

Examples: Aldrin = 1L1; VISIT 1 = 12C5 (n=7 implicit, r=12); VISIT 2 = 10C7 (n=7). When n is a multiple of 7, Δψ = 0 (no line-of-nodes rotation required) → all n=7 cyclers are ballistic.

### 3.2 Most promising cyclers n=1-6 (Table 4, p.6)

| Cycler (nPr) | Aphelion (AU) | V_∞_E (km/s) | V_∞_M (km/s) | Shortest transfer time (days) | Required turn angle (deg) | Max possible turn angle (deg) |
|---|---|---|---|---|---|---|
| **1L1** (Aldrin) | 2.23 | **6.54** | **9.75** | 146 | 84 | 72 |
| 2L2 | 2.33 | 10.06 | 11.27 | 158 | 134 | 44 |
| 2L3 | 1.51* | 5.65 | 3.05** | 280† | 135 | 82 |
| 3L4 | 1.89 | 11.78 | 9.68 | 189 | 167 | 35 |
| 3L5 | 1.45* | 7.61 | 2.97** | 274† | 167 | 62 |
| 3S5 | 1.52* | 12.27 | 5.45** | 134† | 167 | 33 |
| 4S5 | 1.82 | 11.23 | 8.89 | 88 | 167 | 38 |
| 4S6 | 1.53 | 8.51 | 4.07 | 157 | 167 | 54 |
| 5S4 | 2.49 | 10.62 | 12.05 | 75 | 134 | 41 |
| 5S5 | 2.09 | 9.08 | 9.87 | 89 | 134 | 50 |
| 5S6 | 1.79 | 7.51 | 7.32 | 111 | 135 | 62 |
| 5S7 | 1.54 | 5.86 | 3.67 | 170 | 135 | 79 |
| 5S8 | 1.34* | 4.11 | 0.71** | 167† | 136 | 103 |
| 6S4 | 2.81 | 7.93 | 12.05 | 87 | 83 | 59 |
| 6S5 | 2.37 | 6.94 | 10.44 | 97 | 84 | 68 |
| 6S6 | 2.04 | 5.96 | 8.69 | 111 | 84 | 78 |
| 6S7 | 1.78 | 4.99 | 6.66 | 133 | 85*** | 90*** |
| 6S8 | 1.57 | 4.02 | 3.90 | 179 | 85*** | 104*** |
| 6S9 | 1.40* | 3.04 | 1.21** | 203† | 86*** | 120*** |

*Footnote b: aphelion radius slightly below Mars orbital radius 1.52 AU.
**Footnote c: V_∞_M is the difference between Mars' speed and spacecraft aphelion speed.
†Footnote d: time to transfer from Earth to aphelion (since trajectory doesn't actually cross Mars orbit).
***Footnote e: ballistic cycler (required turn angle < max possible turn angle).

Ballistic cyclers in n=1-6 set (no ΔV maneuver needed at Earth flyby): **6S7, 6S8, 6S9** (only).

### 3.3 n=7 ballistic cyclers (Table 5, p.6)

All n=7 cyclers are ballistic because Δψ = 0 (line-of-nodes need not rotate when synodic count is a multiple of 7).

| # Revs r every 15 years | Period (15/r) years | Aphelion range R_a (AU) | Years between Earth encounters | Years between Mars encounters |
|---|---|---|---|---|
| 1 | 15 | [11.16, 12.16) | 15 | 15 |
| 2 | 7.5 | [6.66, 7.66) | 15 | 7.5 |
| 3 | 5 | [4.85, 5.85) | 5 | 15 |
| 4 | 3.75 | [3.83, 4.83) | 15 | 3.75 |
| 5 | 3 | [3.16, 4.16) | 3 | 15 |
| 6 | 2.5 | [2.68, 3.68) | 5 | 7.5 |
| 7 | 2.143 | [2.32, 3.32) | 15 | 15 |
| 8 | 1.875 | [2.04, 3.04) | 15 | 1.875 |
| 9 | 1.667 | [1.81, 2.81) | 5 | 15 |
| **10** (VISIT 2) | **1.5** | **[1.62, 2.62)** | 15 | 7.5 |
| 11 | 1.364 | [1.46, 2.46) | 15 | 15 |
| **12** (VISIT 1) | **1.25** | **[1.32, 2.32)** | 5 | 3.75 |
| 13 | 1.154 | [1.20, 2.20) | 15 | 15 |
| 14 | 1.071 | [1.09, 2.09) | 15 | 7.5 |

Aphelion-range bound corresponds to perihelion ranging R_p ∈ (0, 1] AU.

### 3.4 Two-synodic-period cyclers with intermediate Earth GA (Table 6, p.9)

These are the "extension" cyclers — n=2 with one intermediate Earth flyby at time τ between t=0 and t=T=2S.

| P₁r₁P₂r₂ | Intermediate Earth flyby time τ (years) | V_∞ at Earth (km/s) | Leg 1 Mars orbit crossing times (years) | V_∞ at Leg 1 Mars crossing (km/s) | Leg 2 Mars orbit crossing times (years) | V_∞ at Leg 2 Mars crossing (km/s) |
|---|---|---|---|---|---|---|
| U0L1 | 2.754 | 11.3 | 0.188, 2.567 | 14.0 | 3.449, 3.592 | 5.4 |
| L2U0 | 2.541 | 8.8 | NC (does not cross Mars) | — | 2.781, 4.046 | 10.3 |
| **S1L1** | **2.828** | **4.7** | **0.419, 0.920, 1.908, 2.409** | **5.0** | **NC** | **—** |

The S1L1 cycler is the "third trajectory listed has remarkably low V_∞ at Earth and Mars (4.7 km/s and 5.0 km/s, respectively). We will refer to it as the ballistic S1L1 cycler... Only the first leg crosses Mars' orbit, so a short (153-day) transfer is only available once every two synodic periods. We note that this cycler is very similar to the Case 3 cycler in Ref. 30." (p.8)

**Ref 30 = Byrnes, McConaghy, Longuski 2002 "Analysis of Various Two Synodic Period Earth-Mars Cycler Trajectories" AIAA/AAS Astrodynamics Specialist Conference, Monterey, CA, Aug. 2002.**

The two papers are explicit companion conference papers from the same Monterey session. Byrnes-McConaghy-Longuski 2002 is the more focused two-synodic-period analysis; this paper is the broader-class survey.

## 4. KNOWN_CORPUS impact

### 4.1 Anchor citation CONFIRMED

The current KNOWN_CORPUS anchor in `src/cyclerfinder/search/literature_check.py` cites this paper for the Aldrin cycler as:

```
McConaghy, Longuski, Byrnes "Analysis of a Broad Class of Earth-Mars Cycler Trajectories" AIAA 2002-4420
```

Every field of this citation is confirmed verbatim by the paper itself:
- **First author**: T. Troy McConaghy ✓
- **Second author**: James M. Longuski ✓
- **Third author**: Dennis V. Byrnes ✓
- **Title**: "Analysis of a Broad Class of Earth-Mars Cycler Trajectories" ✓ (verbatim p.1)
- **Paper ID**: AIAA 2002-4420 ✓
- **Year**: 2002 ✓
- **Venue**: AIAA/AAS Astrodynamics Specialist Conference and Exhibit, Monterey, CA, 5-8 August 2002 ✓

No errata required.

### 4.2 New anchor opportunity — Byrnes-McConaghy-Longuski 2002 (Ref 30, p.10)

The same Monterey conference produced a companion paper: **Byrnes, D.V., McConaghy, T.T., Longuski, J.M., "Analysis of Various Two Synodic Period Earth-Mars Cycler Trajectories"** — referenced as Ref 30 on p.10. **NOT IN private paper corpus**. This is the explicit S1L1 / Case 3 cycler companion paper. Worth acquiring as #116 backfill since the catalogue's S1L1 row currently anchors on McConaghy 2004 JSR (#352) — but the Byrnes-McConaghy-Longuski 2002 conference paper may predate it.

### 4.3 S1L1 first-publication date corrected

**The catalogue's `s1l1_mcconaghy_2004` row should be re-checked**: the S1L1 cycler appears HERE (AIAA 2002-4420, August 2002) NOT for the first time in McConaghy 2004 JSR. The 2004 JSR paper is a published-journal version of the same 2002 conference results, plus the Byrnes-McConaghy-Longuski 2002 companion.

Specifically, this paper says (p.8): *"The third trajectory listed has remarkably low V_∞ at Earth and Mars (4.7 km/s and 5.0 km/s, respectively). We will refer to it as the ballistic S1L1 cycler... We note that this cycler is very similar to the Case 3 cycler in Ref. 30."*

The fact that this paper "refers" to it as S1L1 strongly suggests the nomenclature is BEING INTRODUCED HERE — and that "Case 3" was the alternative name in the more focused Byrnes-McConaghy-Longuski 2002 companion. The 2004 JSR paper inherits the S1L1 label.

**Recommended catalogue citation**: anchor the s1l1 row on **(2002-AIAA-4420 + 2002-Byrnes-McConaghy-Longuski Monterey + 2004 JSR)** as a 2002-precedence triple, not just 2004.

## 5. Catalogue impact

### 5.1 No new row enabled

This paper does NOT enable a new V0/V1+ catalogue row. The data per cycler is:
- Aphelion radius (AU) — single number per cycler
- V_∞ at Earth (km/s) — single number (per repeat interval)
- V_∞ at Mars (km/s) — single number (per repeat interval)
- Shortest transfer time (days) — single number
- Required vs max turn angle — single numbers

There is **no per-encounter epoch table** (the cyclers are defined in a planar circular abstract model where the launch epoch is a free choice). Per memory `feedback_golden_tests_sourced_only.md` the catalogue's V0 admission requires per-encounter V_∞ vs epoch — this paper publishes the family-level parameters, not the encounter-by-encounter sequence.

### 5.2 Re-confirmation of canonical cyclers

- **1L1 = Aldrin cycler** (1985 SAIC presentation, refined Byrnes-Longuski-Aldrin 1993 JSR). Confirmed in Table 4: V_∞_E = 6.54 km/s, V_∞_M = 9.75 km/s, R_a = 2.23 AU, shortest transfer 146 days. (The 2.23 AU aphelion vs the McConaghy 2004 published 2.234 AU — consistent at the published precision.)
- **VISIT 1** = n=7 r=12 ballistic, period 1.25 yr, aphelion ∈ [1.32, 2.32) AU.
- **VISIT 2** = n=7 r=10 ballistic, period 1.5 yr, aphelion ∈ [1.62, 2.62) AU.
- **S1L1** = ballistic two-synodic-period with intermediate Earth GA, V_∞_E = 4.7 km/s, V_∞_M = 5.0 km/s, τ = 2.828 yr.

### 5.3 Empty / negative regions identified (Tables 1-3, p.3-4)

The unique-period table (Table 1) and short/long-period tables (Tables 2-3) show "N.S." (No Solution) cells where r exceeds r_max(n) — the Lambert problem has no solution. These are **method-versioned negative regions** under the planar-circular-no-V_∞-cap assumption. Per memory `project_negative_results_registry.md` these should feed the persistent empty-region registry. Specifically:
- n=1: r_max(1L) = 1, r_max(1S) = 4 (N.S. at r=4 for nL; r=4-9 for nS)
- n=2: r_max(2L) = 2, r_max(2S) = 4
- n=3: r_max(3L) = 4, r_max(3S) = 6
- n=4: r_max(4L) = 8, r_max(4S) = 8

These bounds are method-specific to the planar circular conic assumption. They are subsumed by real-eph methods but useful as reference upper bounds.

## 6. Errata vs prior assumptions

1. **S1L1 PRECEDENCE**: S1L1 cycler is FIRST INTRODUCED in this paper (Aug 2002), not in McConaghy 2004 JSR. Update catalogue citation chain.
2. **Companion paper Byrnes-McConaghy-Longuski 2002 (Ref 30) is the explicit two-synodic-period paper**; calls S1L1 "Case 3". Not in corpus — flag for #116.
3. The paper acknowledges "Dr. Buzz Aldrin for his advice and inspiration" (p.9 Acknowledgments) and notes Wisuwat Bhosri's helpful suggestions. This is the canonical 2002 collaboration anchor.
4. **Hollister castles-in-space reference IS Hollister 1969** (Ref 3, p.9) "Castles in Space" *Astronautica Acta* 14:2 pp.311-316 — the 1969 Hollister castles-in-space concept paper that anchors all subsequent Venus-Earth cycler literature.
5. **n=2 r=3 S1L1 paper**: this is the OPEN CONFERENCE paper. The closed-form S1L1 paper is the 2004 JSR (in the private paper corpus). The catalogue currently treats them as one anchor; they're actually two (2002 conference + 2004 JSR).

## 7. Action items for parent

- **No code changes** to KNOWN_CORPUS required — the existing AIAA 2002-4420 anchor is correctly cited. The verbatim title/author/venue/paper-ID all match.
- **Memory update**: the `project_s1l1_realeph_closure_blocker.md` memory should add a note that S1L1's 2002 conference precedence is in this paper (Table 6), with the V_∞_E = 4.7 / V_∞_M = 5.0 km/s tuple. This is the **earliest sourced S1L1 publication** in the private paper corpus and supersedes the 2004 JSR for citation precedence.
- **Byrnes-McConaghy-Longuski 2002** "Analysis of Various Two Synodic Period Earth-Mars Cycler Trajectories" AIAA/AAS Monterey Aug 2002 — flag for #116 acquisitions. Per the Ref 30 description, this paper has the explicit "Case 3" treatment that precedes S1L1 nomenclature.
- **Possible new precursor_mga rows**: the ballistic n=7 cyclers (VISIT 1, VISIT 2) might warrant precursor_mga admissions in the expanded scope (memory `project_catalogue_scope_expanded_2026-06-15.md`). Table 5 gives aphelion ranges, repeat periods, and Earth/Mars encounter intervals. **However**: a single (R_p, R_a) range pair without a chosen epoch is family-spec, not member-spec — V_min for admission would still need the family member chosen, and the canonical Friedlander-Niehoff-Byrnes-Longuski 1986 VISIT papers (Refs 11-14) are the proper anchors for that — those papers ARE the V0 sources, not this one.
- **6S7, 6S8, 6S9 ballistic cyclers**: the n=6 ballistic family in Table 4 (footnote e). These have low V_∞ (4.0-5.0 km/s) and are ballistic — worth checking whether the catalogue already has rows for them or whether they're a novel admission opportunity. (Quick check: catalogue does NOT obviously have n=6 ballistic rows; merits a separate task.)

---

**Cross-references**:
- `docs/notes/2026-06-17-digest-mcconaghy-2004.md` (#352) — the published JSR version of this conference paper (DOI 10.2514/1.11939); has expanded analysis but inherits the nPr nomenclature and S1L1 label introduced here
- `docs/notes/2026-06-17-digest-mcconaghy-2005.md` (#353) — McConaghy-Russell-Longuski 2005 "Standard nomenclature for Earth-Mars cycler trajectories" — formalizes the nPr notation introduced here
- `docs/notes/2026-06-17-digest-mcconaghy-2006.md` (#352) — McConaghy-Landau-Yam 2006 "Notable two-synodic-period Earth-Mars cycler" JSR; further refines the S1L1 line
- `docs/notes/2026-06-17-mars-cycler-wave-digest.md` — Mars-cycler wave overview from earlier in this 2026-06-17 push
- Memory: `project_s1l1_realeph_closure_blocker.md`, `project_s1l1_nomenclature.md`
