# Digest: Pisarevsky, Kogan & Guelman 2008, "Interplanetary Periodic Trajectories in Two-Planet Systems" (#960)

D. M. Pisarevsky, A. Kogan and M. Guelman, "Interplanetary Periodic Trajectories in Two-Planet Systems",
Journal of Guidance, Control, and Dynamics 31(3):729-739 (May-June 2008), DOI 10.2514/1.30046
(Crossref-confirmed 2026-10-05: title, authors, 31(3), pp. 729-739). Received 26 Jan 2007, accepted 22 Nov
2007 (p.729). Filed as `cyclers_pdf/papers/pisarevsky-kogan-guelman-2008-interplanetary-periodic-trajectories-two-planet-systems-JGCD-31-3-729-doi-10.2514-1.30046.pdf`
(11 pages, text layer, md5 b63447eeb0c5fed73dc85d0d4453a604). Journal page p = PDF page + 728.
All 11 pages were read. The text layer was read in full, and the Table 4 page (p.739) was read as a page image.

Evidence tags: READ = printed in the paper (page cited). COMPUTED = my own check on 2026-10-05 (scratch,
not committed). INFERRED = my reading.

Why it was acquired: `#938` review R1 named it a hard gate for `#942` (R1) and `#943` (X1). It is "the
only published two-planet periodic-trajectory method with flybys at both bodies" that the corpus knew of
by citation. Russell & Strange 2009 cite it as ref. [18].

## 0. Verdict in one paragraph

This paper gives a method and a classification for cyclers in a circular, coplanar, two-planet system with
patched-conic "collisional" arcs. Both planets may provide gravity assists. It extends Hénon's 1968
diagram. All numbers are for Earth-Mars:
- one numeric two-massive-planet cycler. It is spatial, class III, has a period of two synodic periods,
  and is purely ballistic (Table 4).
- graphical candidate sets (Figs. 6, 14).

It does NOT compute Earth-Venus, Venus-Mars or any planet-moon pair. It gives no state vectors, no
ephemeris continuation and no catalogue. So R1 cell (a) is PARTIAL. Cells (b) and (c) and X1 get only a
method antecedent from this paper.

## 1. Model (READ pp.729-731)

- "Restricted four-body problem": a central body, two primaries on circular coplanar orbits and a
  spacecraft. The primaries' masses are small and their spheres of influence are infinitesimal, so the
  model is patched conic (p.730, Sec. II.A).
- Flybys are Poincaré "consecutive collisions": the limit mu -> 0 with the flyby distance going to zero,
  and |V_inf| conserved (p.730, Sec. II.B).
- The turn limit is applied separately, Eq. (15): sin(delta/2) = 1/(1 + (rho + h_p) V_inf^2 / mu_2)
  (p.732). Until it is checked, solutions are only "potential candidates" (p.732, Sec. III.B).
- Arcs are either recursive (both ends at the same planet) or transitional (the ends at different
  planets).
  - Recursive arcs are even (2k pi), odd ((2k-1) pi) or generic (Fig. 1). Even and odd arcs may be out of
    plane. Generic arcs are coplanar (p.730).
  - Transitional arcs are fast, slow, short or long, by perihelion/aphelion passage count (Table 1, p.731).
- Hénon's timing equation is Eq. (4) (p.731). The relative-velocity level lines are Eq. (12) (p.732).
  Their intersections give all collisional arcs at one V_inf (Fig. 5, V_inf = 1 canonical).
- The authors restrict loitering arcs to multiples of pi (p.730).

## 2. One-working-body case: Earth-Mars, Mars massless (READ p.733, Sec. IV.A, Fig. 6)

- Periodicity condition, Eq. (16): 2 tau + T_E sum_j n_E,j / 2 = N T_syn. One generic arc plus any number
  of odd and even arcs. One arc must have aphelion above Mars's orbit.
- "Figure 6 presents all possible solutions for generic arcs for {N, n_E} = {1, 0} and {2, 3}. The figure
  shows that there are 11 and 13 solutions, respectively" (p.733). The diagram bounds are those of Hénon's
  Fig. 4, tau/pi < 5.5 and eta/pi < 6.5 (p.731).
- Known cyclers placed on Fig. 6:
  - near-ballistic: the Aldrin cycler, {2,3}, {3,1} (two full revolutions), {2,1} and {3,5} (two full
    revolutions).
  - purely ballistic: {2,5}, {3,1} and {3,5} (one full revolution), {3,3}, {3,7} and {3,9}.
  - All are "described in more detail in [12]" (Russell & Ocampo 2004) (p.733).
- The Aldrin turn is "about 15% above the maximal [12]" (p.733).
- This section reproduces the Russell-Ocampo set graphically. It is not new.

## 3. Two-working-body case (READ pp.733-739, Sec. IV.B-G)

- Every cycler is one P1-recursive trajectory plus optional P1 arcs (Fig. 7, p.734).
  - Leg 1 is P1 -> P2 (transitional).
  - Leg 2 is optional P2 arcs.
  - Leg 3 is P2 -> P1 (transitional).
  - Leg 4 is optional P1 arcs.
- Equal V_inf at both ends gives a1 = a3 (Eqs. 17-20). If both transitional arcs are coplanar, then also
  e1 = e3 (p.733-734).
- Classes come from pairs of transitional arcs: 21 combinations in five classes (Tables 2-3, p.734).
  - I.1: P1 inner. {fast, long} and {short, slow}, coplanar.
  - I.2: P1 outer. {short, fast} and {slow, long}, coplanar.
  - II: the same family twice, symmetric in both synodic frames.
  - III: one coplanar and one noncoplanar arc, where the noncoplanar transfer angle is (2k-1) pi.
  - IV and V: resonant two-planet systems only.
- Synchronisation:
  - class I: Eqs. (23)-(33). Hénon's diagram is the special case without leg 2.
  - class II: Eqs. (34)-(43), using Lagrange F, G.
  - class III: Eqs. (44)-(52), where cos i3 is given by Eq. (52).
  - The cycler period condition is Eq. (53), t1 + t2 + t3 + t4 = N T_syn. For class I.1 it becomes
    Eq. (55).
- Recursive diagrams, all Earth-Mars:
  - Fig. 10a: class I.1, n2 = m2 = 1.
  - Fig. 10b: class I.2, n2 = m2 = 2.
  - Fig. 12a-d: class II, with n2 = m2 = 0 and 2.
  - Fig. 13: class III, n2 = m2 = 2.
  - The class II diagram has infinitely many curves, "twice as large" as Hénon's (p.737).
- **Fig. 14 (p.738):** "Candidates for periodic trajectories in the Earth-Mars system for class I.1 with
  n2 = m2 = 1". The top labels are {N, n4} = {2,3} and {3,4} (text p.739). These are graphical points
  only. No turn-angle feasibility, no elements and no states are printed.
- **Table 4 and Fig. 15 (p.739):** "an example of a purely ballistic class-III cycler in the system of
  massive Earth and Mars". Its period is two synodic periods, the Earth-to-Mars time is 134 days and leg 2
  is absent. "The gravity of Mars is used to change inclination of the spacecraft's trajectory." The
  printed Mars period is "2.13 years", citing [13] (Russell & Ocampo 2005).

Table 4 as printed (READ p.739, page image):

| Leg | Arc | a, AU | e | i, deg | theta, deg | V_inf, km/s | flyby height, km |
|---|---|---|---|---|---|---|---|
| 1 | 1 | 1.3 | 0.284 | 0 | 93.54 | 5.7 | 7617 |
| 3 | 2 | 1.3 | 0.269 | 5.3 | 540 | 6.2 | 2655 |
| 4 | 3 | 1 | 0 | 12 | 180 | 5.7 | 2834 |
| 4 | 4 | 1 | 0.16 | 7.6 | 360 | 5.7 | 7617 |

"The final flyby hyperbolic elements correspond to the arc beginning" (p.739).

**Consistency check (COMPUTED, Tisserand with v_E = 29.785 km/s and r_M = 1.5237 AU).** Each arc's
(a, e, i) gives the following V_inf at Earth (r = 1) and at Mars:

| Arc | Earth V_inf | Mars V_inf |
|---|---|---|
| 1 | 6.27 | 5.74 |
| 2 | 6.24 | 5.72 |
| 3 | 6.23 | (does not reach Mars) |
| 4 | 6.18 | (does not reach Mars) |

The model is self-consistent: Earth about 6.2 km/s on all four arcs and Mars about 5.7 km/s. But the
printed V_inf column carries 5.7 on arcs 1, 3 and 4, which begin at Earth, and 6.2 on arc 2, which begins
at Mars. The two values look interchanged.

With r_M = 2.13^(2/3) = 1.656 AU, that is, if the printed 2.13 yr were the Mars sidereal period used in
the model, Mars V_inf would be about 3.8 km/s, which matches nothing. 2.13 yr is the Earth-Mars synodic
period (2.135 yr). So the "orbital period of Mars" wording is most likely a slip for the synodic period.

Both points are erratum candidates. They are offered with respect: they are typesetting-class slips that
the table's own elements resolve. The control should use (a, e, i, theta, t1 = 134 d) and not the V_inf
column.

## 4. Gate answers for `#942` (R1) and `#943` (X1)

- **R1 cell (a), Earth-Mars where Mars also gives gravity assists: PARTIAL.**
  - What is published:
    - the class and the method (Secs. IV.B-G).
    - one numeric member, Table 4. It is spatial, class III, has a 2-synodic period and is ballistic, and
      the Mars turn is used for inclination only.
    - the graphical class I.1 candidates of Fig. 14, {N, n4} = {2,3} and {3,4}.
  - What is left:
    - an enumeration with turn-angle gating (the paper says Fig. 14 entries are "candidates").
    - numeric coplanar class I and II members (none printed).
    - multi-generic-return chains, which the paper names as future work (p.739: "including any number of
      generic returns to both planets").
    - ephemeris continuation.
  - The `#938` R1 claim that "every held source states the exclusion" no longer holds once this paper is
    held. Any new Earth-Mars two-working-body member must be collision-checked against Fig. 14, which needs
    digitising and is therefore figure-derived, and against Table 4.
  - Attribution: Pisarevsky, Kogan & Guelman 2008 for the class.
- **R1 cell (b), Earth-Venus with Venus returns: OPEN in this paper.** The paper computes no Earth-Venus
  case. Its method claims "any circular, coplanar, two-planet system" (Table 2 caption context, p.734). So
  the method is published, and cite it as the method antecedent alongside Hollister & Menning 1970.
- **R1 cell (c), Venus-Mars: OPEN in this paper.** It is not computed here. See the AAS 07-118 digest for
  the published one-working-body Venus-Mars member.
- **X1 (Ganymede-Callisto, Ganymede-Europa): OPEN in this paper.** No planet-moon pair is computed. The
  four-body formalism applies to any central body with two small primaries on circular coplanar orbits
  (INFERRED from Sec. II.A), so cite it as the method antecedent.
- **Extended Hénon diagram: what it covers and what it leaves out.**
  - Covered: recursive trajectories in classes I-V for loitering arcs that are multiples of pi, coplanar
    and noncoplanar transitional arcs, and Earth-Mars diagrams only.
  - Left out:
    - generic (non-k pi) loitering arcs at P2 and at P1.
    - more than one generic return per planet (future work, p.739).
    - eccentric or inclined planet orbits.
    - turn-angle screening of the diagram points.
    - any other planet pair.
    - any moon system.

## 5. Positive controls with sourced numbers

1. **Table 4 class III cycler (p.739).** a = 1.3 AU on both transitional arcs, e = 0.284 and 0.269, i = 0
   and 5.3 deg, theta = 93.54 and 540 deg. The Earth arcs are (1, 0, 12 deg, 180 deg) and
   (1, 0.16, 7.6 deg, 360 deg). t1 = 134 d. Period 2 T_syn. Flyby heights 7617, 2655, 2834 and 7617 km.
   - The values are only 2-3 significant figures. Use it as a loose structural control: a = 1.3 AU and
     t1 = 134 d within the rounding.
   - Do not use the V_inf column (see sec. 3).
2. **Solution counts (p.733).** 11 generic-arc solutions at {N, n_E} = {1, 0} and 13 at {2, 3}, inside
   tau/pi < 5.5 and eta/pi < 6.5. This is a discrete control on any Hénon-diagram implementation.
3. **Fig. 6 placements (p.733).** The Russell-Ocampo 2004 cyclers (Aldrin and the {N, n_E} list in sec. 2)
   must fall on the diagram. This duplicates the existing `#913` / Russell-Ocampo controls.

## 6. Citation mining (policy step 4)

The background and related-work references, refs. [1]-[18], pp.739, were checked against
`CORPUS_INDEX.md`, the digest bodies and the filenames in `cyclers_pdf/papers/` on 2026-10-05.

Held:
- [6] Byrnes-Longuski-Aldrin 1993.
- [7] McConaghy-Longuski-Byrnes AIAA 2002-4420.
- [9] Byrnes-McConaghy-Longuski AIAA 2002-4423.
- [11] Hénon 1968, `henon-1968-orbites-interplanetaires-...pdf`.
- [14] Russell & Ocampo 2006.
- [15] Szebehely 1967.
- [12] Russell & Ocampo 2004: the journal version is not held. Its content is covered by
  `russell-ocampo-2003-...-AAS-03-145.pdf` and `russell-2004-dissertation.pdf`.
- [13] Russell & Ocampo 2005 JGCD, Global Search: not held as a journal file. Covered by the Russell 2004
  dissertation (INFERRED, the dissertation chapters).
- [8] McConaghy et al. AAS 03-509: covered by `mcconaghy-2004-...-purdue-phd.pdf` and
  `mcconaghy-landau-yam-2006-...jsr...pdf` (INFERRED from title match).
- [2] and [3], the Menning 1968 MS and Rall 1969 PhD theses: partly covered by
  `hollister-menning-1970-...JSR-7-10.pdf` and `rall-1969-free-fall-periodic-orbits-connecting-earth-mars-mit-scd-thesis-msl-te-34-ntrs-19700017824.pdf` (formerly `hollister-rall-1970-periodic-orbits-NASA-CR.pdf`). CORRECTION 2026-10-05: the second file IS Rall's thesis (MSL report TE-34), so Rall [3] is fully held, not partly covered.
  Menning [2] has also been held since 2026-10-05 (`menning-1968-...-mit-ms-thesis-ocr.pdf`).

Not held, in priority order. The DOI was checked through the Crossref API (title, authors, volume and
pages).

1. Hollister, W. M. (1969), "Periodic Orbits for Interplanetary Flight", J. Spacecraft and Rockets
   6(4):366-369, doi 10.2514/3.29664 (CONFIRMED). Unlocks: the original Earth-Mars and Earth-Venus
   periodic-orbit record (R1 literal collision).
2. Patel, M. R., Longuski, J. M. & Sims, J. A. (1998), "Mars Free Return Trajectories", J. Spacecraft and
   Rockets 35(3):350-354, doi 10.2514/2.3333 (CONFIRMED). Unlocks: Hénon-table free returns relative to
   Earth (the R1 one-body generator).
3. Russell, R. P. & Ocampo, C. A. (2004), "Systematic Method for Constructing Earth-Mars Cyclers Using
   Free-Return Trajectories", JGCD 27(3):321-335, doi 10.2514/1.1011 (CONFIRMED). The journal version of
   held material.
4. Russell, R. P. & Ocampo, C. A. (2005), "Global Search for Idealized Free-Return Earth-Mars Cyclers",
   JGCD 28(2):194-208, doi 10.2514/1.8696 (CONFIRMED). The journal version of held material.
5. Rall, C. S. (1969), "Freefall Periodic Orbits Connecting Earth and Mars", PhD, MIT. No DOI (thesis).
6. Menning, M. D. (1968), "Freefall Periodic Orbits Connecting Earth and Venus", MSc, MIT. No DOI (thesis).
7. Niehoff, J. (1986), "Pathways to Mars: New Trajectory Opportunities", AAS 86-172. No DOI (AAS paper).
8. Niehoff, J., Friedlander, A. & McAdams, J. (1991), "Earth-Mars Transport Cycler Concepts", IAF-91-438.
   No DOI (IAC paper).
9. Bruno, A. D. (1994), "The Restricted 3-Body Problem: Plane Periodic Orbits", de Gruyter, doi
   10.1515/9783110901733 (CONFIRMED). The collision-orbit (second-species) theory behind Eq. (1).
10. Poincaré, H. (1899), "Les Méthodes Nouvelles de la Mécanique Céleste", Tome 3. No DOI (book; public
    domain).
11. Battin, R. H. (1999), "An Introduction to the Mathematics and Methods of Astrodynamics", Revised
    Edition, AIAA, doi 10.2514/4.861543 (CONFIRMED). A textbook, low priority.

## 7. Catalogue and code implications (proposals only; no catalogue writes)

- A candidate catalogue row: the Table 4 class III Earth-Mars cycler. The class is two-working-body and
  spatial, the period is 2 T_syn, and it is ballistic in the ideal model. The fields are low-precision and
  there are no states. V0 at most. The V_inf erratum must be noted.
- A `literature_check.py` KNOWN_CORPUS anchor was added in this task (`pisarevsky-2008-two-planet`):
  body_set {E, M}, topology repeated-moon. Any Earth-Mars candidate in which Mars bends is then flagged
  against this paper.
