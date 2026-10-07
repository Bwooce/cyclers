# Digest: Voyatzis & Kotoulas 2005, "Planar periodic orbits in exterior resonances with Neptune" (Planetary and Space Science 53:1189) (#960 batch 36)

George Voyatzis and Thomas Kotoulas (Department of Physics, University of Thessaloniki), Planetary and Space Science
53(11):1189-1199 (2005), doi:10.1016/j.pss.2005.05.001. Received 10 November 2004, revised 26 April 2005, accepted
15 May 2005. Title, authors, volume, issue and pages agree with Crossref.
- File: `39a3b1ac-voyatzis2005.pdf`, 11 pp., md5 `4bb4af5446c24d3569a9e6221182967a`. Publisher PDF ("ARTICLE IN
  PRESS" headers) with a good text layer.
- **Proposed corpus filename:**
  `cyclers_pdf/papers/voyatzis-kotoulas-2005-planar-periodic-orbits-exterior-resonances-neptune-pss-53-1189-doi-10.1016-j.pss.2005.05.001.pdf`
- **Wanted list:** this paper is row 22 (Tier B, `#971` R14 collision check). R14 in the `#971` review calls it
  "wanted row 52" and the Varin-Bruno 2009 digest calls it "row 56": both are older numberings.
- **How I read it:**
  - Whole text layer. On page images (200 dpi, zoomed crops; Fig. 1 at 400 dpi): eqs. (1)-(3) and the mu value
    (p.1190-1191), the closed-curve passage (p.1192), Fig. 1 and Fig. 2 (p.1192), and every cell of Table 1 (p.1195).
    The text layer agreed with the image in every Table 1 cell.
  - For the R14 question I also read Bruno & Varin JAMM 2007 sec. 4.2 (p.944) and sec. 6.5 (p.953, image).
  - Scripts and outputs in this folder: `check_vk2005_table1.py` -> `check_vk2005_table1.out` (fixed-period
    shooting); `check_vk2005_fixed_h.py` -> `check_vk2005_fixed_h.out` and `check_vk2005_fixed_h_35.py` ->
    `check_vk2005_fixed_h_35.out` (fixed-h shooting for two rows the first method missed).
  - YAML: `voyatzis-kotoulas-2005-tables.yaml` (Table 1, printed values). R14 verdict: `r14-verdict.md`.

## 0. Verdict

**A survey of the symmetric resonant families of the exterior first-, second- and third-order resonances with
Neptune (2/3 to 7/8, 3/5, 5/7, 7/9, 4/7, 5/8, 7/10), in the planar circular problem at mu = 5.178e-5, and of the
families they start in the planar elliptic problem.** It is a Kuiper-belt paper. The orbits are heliocentric
resonant orbits with x0 from about 1 to 3 (units of Neptune's orbital radius), not cyclers.
- **The thing Bruno and Varin cite it for is in Fig. 1:** the families I_6/7, II_7/8 and a new unstable family
  II^u_7/8 "form a closed characteristic curve" (p.1192). In Bruno's notation this is the closed family containing
  the circular piece Id_{-7} (my mapping, sec. 2). The other five curves of Fig. 1 are Id_{-2} ... Id_{-6}, so Fig. 1
  holds all six families Bruno-Varin describe ("p = -2 to -7"), again by my mapping.
- **R14 collision: partial.** The closed p = -7 family is exterior (a > 1) and at the Neptune mu only, with no
  numbers. None of R14's interior targets i_5-i_8 (a < 1, family i) is here. R14 is not pre-empted, but its novelty
  wording must be limited to the interior side. Details in `r14-verdict.md`.
- **Data:** the only table (Table 1) gives (e0, h) and stability codes of 34 bifurcation points (BPs) where
  elliptic-problem families start (6 circular BP0, 28 others). No initial conditions (x0, ydot0) are printed
  anywhere. I checked the six BP0 h values by a Kepler estimate (all agree) and 11 of the 28 others by integration:
  8 agree in h within rounding, 2 miss by one unit in the last printed digit, 1 I could not pin (sec. 3).
- **What it gives the project:**
  - A literature anchor for R14 (prior art for the exterior closed family) and a shape-level positive control:
    a correct closed-family solver at mu = 5.178e-5 must find a closed loop in the 6/7-7/8 band with folds near
    (x0, h) = (1.04, -1.502) and (1.335, -1.474) (read from Fig. 1, about +-0.005).
  - Table 1 rows as cheap controls for a period-fixed symmetric shooter at a small mu (sec. 3).
- **Catalogue implication (PROPOSAL only):** none. No row is a cycler.

## 1. Frame and conventions (stated first; pp.1190-1191, image-checked)

- Planar restricted problem, Sun (mass 1 - mu) and Neptune (mass mu), G = 1, a' = 1, T' = 2 pi. Rotating frame Oxy,
  origin at the barycentre O: r1^2 = (x + mu r)^2 + y^2, r2^2 = (x - 1 + mu r)^2 + y^2, with r = 1 in the circular
  problem. So the Sun is at x = -mu and Neptune at x = 1 - mu (also in the Fig. 5 caption).
- **mu = 5.178 x 10^-5** (p.1191 and p.1198; one value for the whole paper).
- Jacobi constant, eq. (2): h = (xdot^2 + ydot^2 - (x^2 + y^2))/2 - (1 - mu)/r1 - mu/r2. This is -C/2 for the usual
  barycentric C. Near-circular exterior orbits have h about -1.5.
- Symmetric orbits: x(0) = x0, y(0) = 0, xdot(0) = 0, ydot(0) = ydot0; periodicity y(T/2) = xdot(T/2) = 0, eq. (3).
  Families are drawn in the (x0, h) plane. "Multiplicity" = number of same-direction crossings of y = 0 per period.
- Families: C (circular, first kind); for each resonance p/q (q > p, order q - p) two resonant families I (body
  "initially at perihelion") and II ("at aphelion") (p.1191). For second order, family I crosses Ox vertically only
  at x > 0 and family II only at x < 0 (p.1193). Segments between collisions are numbered I^n_{p/q}.
- Elliptic problem: an orbit with period T = 2 k pi in a circular-problem family is a BP; from each BP two families
  start, E^{p/q}_{lp} and E^{p/q}_{la} (Neptune initially at perihelion or aphelion), with e' as the parameter
  (p.1191, p.1195). e0 = e(0) is the eccentricity "that corresponds to the initial conditions" (p.1196); the paper
  does not say which centre or GM. Stability index k = a11 + a22 (circular problem; Hénon 1997) and Broucke's k1, k2
  (elliptic). Newton-Raphson shooting to 1e-13 (circular) or 1e-11 (elliptic), Bulirsch-Stoer integration.

## 2. Content

### 2.1 Circular problem: first-order resonances (sec. 3.1, pp.1191-1193, Figs. 1-3)

- The circular family breaks near each m/(m+1) resonance (the "gap" structure: Guillaume 1974, Hadjidemetriou
  1993). Family I_m/(m+1) joins II_(m+1)/(m+2) smoothly through a short near-circular segment (circles in Fig. 1).
  As x0 decreases, the I branches run into a collision with Neptune at h about -1.5.
- **The exception (p.1192, image):** "A rather different structure is observed for the families I_6/7 and II_7/8
  because the resonant family I_6/7 avoids the collision with Neptune. A new resonant family of unstable periodic
  orbits, denoted by II^u_7/8, is found and the three families (I_6/7, II_7/8 and II^u_7/8) form a closed
  characteristic curve." Fig. 2a: n/n' runs between 6/7 and 7/8 on the loop. Fig. 2b: e from about 0.005 (near the
  x0 = 1.094 circle) to about 0.23 (the fold at x0 about 1.335). The near-circular plateau "decreases rapidly as m
  increases, and disappears for the closed characteristic curve". Both folds sit beside close-encounter marks
  (Fig. 1: the "x" at (1.04, -1.503) and the bold tick at about (1.335, -1.474)).
- **My mapping to Bruno's notation (inference).** Bruno's Id_p is the circular piece with p/(p-1) > N > (p+1)/p
  (JAMM p.944, p.953); for p = -7 that is 7/8 > N > 6/7, bounded by E_{7/8} and E_{6/7}. Fig. 1's right-edge labels
  1-6 sit on six curves. On the 400 dpi render, curve j is I_{j/(j+1)} -> circle -> II_{(j+1)/(j+2)} for j = 2-6
  (circles at x0 = 1.256, 1.183, 1.141, 1.114, 1.094). Curve 1 shows II_2/3 and its circle (1.355); its other
  branch leaves the plot. So labels 1-6 carry Id_{-2} ... Id_{-7}, and label 6 is the closed one. This matches
  Bruno-Varin's sentence "families ... including the pieces Id_p ... for p = -2, -3, ..., -7 were calculated for
  mu = 5.178e-5 ... the last family: it is closed" (JAMM p.953).
- A stray second "I_6/7" label sits at (1.21, -1.517) between the I_3/4 and I_2/3 curves. I cannot place it.
- Full families (Fig. 3, 4/5 and 5/6, p.1193): I^1_m/(m+1) (m <= 6) is unstable, multiplicity one, and ends at a
  Neptune collision; later segments are stable, separated by Neptune collisions (one for 2/3 and 3/4, two for
  m >= 4), and end at a Sun collision. II^1 (m <= 7) ends at a Neptune collision; II^2 ends at a Sun collision
  (m <= 4) or a second Neptune collision (5/6, 6/7). All II families are stable except at close encounters.

### 2.2 Circular problem: second- and third-order resonances (secs. 3.2-3.3, Figs. 4-6)

- Both families bifurcate from C where n/n' = p/q. Family I starts unstable, meets a Neptune close encounter at
  h about -1.5, then continues stable. Second order: I_3/5 reaches a Sun collision; I_5/7 meets a second Neptune
  encounter, then a Sun collision; I_7/9 has three encounters. II_3/5, II_5/7, II_7/9 are stable with one, two and
  three Neptune encounters, and all end at a Sun collision (Fig. 4).
- Third order: I_4/7 ends at its second encounter ("seems to terminate or become strongly chaotic"); I_7/10 the
  same after the third; I_5/8 reaches a Sun collision. II_5/8 and II_7/10 have two encounters, II_4/7 one; II_5/8
  ends at a Sun collision, II_4/7 and II_7/10 at a near "double collision" with Neptune and Sun (Fig. 6b, c).
- No asymmetric families were found in the circular problem (abstract; Conclusions (c)).

### 2.3 Elliptic problem (sec. 4, pp.1195-1198, Table 1, Figs. 7-10)

- First order: BPs only on families II. 4/5, 5/6 have two BPs, 6/7 three. E_la families have e0 rising with e';
  most reach e' about 0.9 before the computation fails. E^{5/6}_1p has an unstable stretch 0.17 < e' < 0.55, whose
  ends are candidate BPs of asymmetric families.
- Second order: BP0 on family C (period q pi, continued with period 2 q pi). E^{3/5}_02 and E^{3/5}_1p join at
  e' about 0.85 (one family from BP0 to BP1); the same holds for 5/7.
- Third order: BP0 on C (period 2 q pi / 3, taken three times). For 4/7, E_1a + E_02 join at e' about 0.13 and
  E_1p + E_2p at e' about 0.8. E^{4/7}_2a stops at e' about 0.296 for no reason the authors can find.
- "We did not find any collision orbits along the families studied in the elliptic RTBP" (p.1198). Unstable orbits
  have stability indices only slightly above 2; Fig. 10b integrates three orbits for 100 My (2 pi = 165 yr).

## 3. Checks (`check_vk2005_table1.out`, `check_vk2005_fixed_h.out`, `check_vk2005_fixed_h_35.out`)

Frame as sec. 1, mu = 5.178e-5, DOP853 rtol = atol = 1e-12. A BP of resonance p/q is taken to have T = 2 q pi
(from the text, p.1191, p.1195-1196). e0 is the Sun-centred osculating eccentricity with GM = 1 - mu (my
assumption; GM = 1 changes it by about 1e-5). With T fixed, the computed h does not depend on the e0 convention.

- **BP0 (6 rows), Kepler estimate, not an integration:** a circular orbit with n = p/q has
  h = -(1 - mu)/(2a) - sqrt((1 - mu) a), O(mu) terms dropped. All six printed h agree within rounding
  (differences -0.0001 to -0.0004).
- **Fixed-period shooting:** seed (x0, ydot0) from the printed (e0, h) at an apsis on the +x or -x axis; Newton on
  (x0, ydot0) for y = xdot = 0 at t = q pi. Where several seeds converge, I list the one closest to the printed
  values. The +x apoapsis and -x periapsis (or -x apoapsis) solutions give the same h: they are the same orbit
  started at its other axis crossing, a useful self-check.

| row | printed e0, h | ours e0, h | difference |
|---|---|---|---|
| 2/3 BP1 | 0.469, -1.393 | 0.4692, -1.39241 | **h +6e-4** (one unit) |
| 3/4 BP1 | 0.329, -1.452 | 0.3291, -1.45203 | agrees |
| 4/5 BP1 | 0.253, -1.473 | 0.2531, -1.47299 | agrees |
| 4/5 BP2 | 0.871, -0.960 | 0.8715, -0.95900 | **h +1.0e-3** (one unit), e0 +5e-4 |
| 5/6 BP1 | 0.205, -1.483 | 0.2052, -1.48278 | agrees |
| 5/6 BP2 | 0.749, -1.146 | (fixed-T shooter missed; fixed-h below) | |
| 6/7 BP1 | 0.172, -1.488 | 0.1720, -1.48816 | agrees |
| 6/7 BP2 | 0.649, -1.253 | 0.6482-0.6484, -1.25267 | h agrees; e0 -6e-4 to -8e-4 |
| 6/7 BP3 | 0.960, -0.743 | 0.9608, -0.74273 | h agrees; e0 +8e-4 |
| 3/5 BP1 | 0.427, -1.428 | 0.4270, -1.42790 | agrees |
| 3/5 BP2 | 0.800, -1.065 | see below | |

- **Fixed-h shooting** (h at the printed value, period free): 5/6 BP2 gives a symmetric orbit with
  T/(2 pi) = 6.00000 at x0 = 1.976456 (aphelion), e0 = 0.7496: the row is confirmed. 3/5 BP2: no apsis state has
  e = 0.800 at h = -1.065 (for e = 0.8 the highest apsis h is -1.0671, at the 3/5 radius; e0 about 0.801 is needed).
  Seeding with e0 = 0.802, the symmetric orbits at h = -1.065 near this point have T/(2 pi) = 4.952 and 5.029
  (+x) and 4.976 and 5.046 (-x), all with e0 = 0.8022. So a T = 10 pi orbit lies close by, but I did not isolate it.
  I record this row as "not reproduced to the printed digits", not as a misprint.
- **Disagreements, not misprints:** 2/3 BP1 (h) and 4/5 BP2 (h, e0) miss by one unit in the third decimal. The
  6/7 BP2, BP3 e0 misses (6e-4 to 8e-4) may come from the e0 convention, which the paper does not define. A likely
  cause of the h misses is how precisely the BPs were located along the family; the paper does not say.
- Rows not run: the fixed-period run was stopped after 3/5 (see the last line of `check_vk2005_table1.out`); 5/7,
  7/9, 4/7, 5/8, 7/10 BP1-BP4 (17 rows) were not checked beyond BP0.
- Nothing in the paper allows a check of the closed family itself: it has no printed numbers.

## 4. Relation to the held corpus

- **Bruno & Varin JAMM 2007 sec. 6.5 and Bruno & Varin SSR 2009 "Other closed families"** (both held) describe this
  paper's Fig. 1. JAMM states that the family i covers the internal annulus, and the external annulus is "a
  denumerable set of families" (p.953). The printed "a > 1" there contradicts "internal annulus" and must mean a < 1
  (a JAMM slip; image-checked).
- **Bruno & Varin 2009b (SSR 43(1):26, families c and i at mu = 5e-5, held):** the interior side at almost the same
  mu. With this paper it gives both sides of P2 at a Neptune-class mu: interior family i (Bruno-Varin) and exterior
  Id_{-2..-7} families (here, figure only).
- **Hadjidemetriou 1993 (CMDA 56:201, held inside the Dvorak-Henrard 1993 volume):** the same circular-to-elliptic
  BP method for interior resonances at Sun-Jupiter; its digest has a tolerance scheme the Table 1 rows could share.
- **Hénon 1997 (held):** the family definitions (first and second kind) and the stability index.

## 5. Text and print slips

- Reference list: Morbidelli, Brown & Levison 2003 is printed as "Earth, Moon, Planets 91, 63-93", the same
  volume and pages as Kotoulas & Hadjidemetriou 2002 in the line above. Malhotra 1996 is printed "111, 540-516"
  (the end page is below the start page). Kotoulas & Voyatzis 2005 (IAU Coll. 197) is printed "Kotoulas, A.".
- Text: "the second order ones 5/7 and 7/9" etc. are rendered as stacked fractions; the text layer garbles 7/10 as
  "10 7". The image is clear in every case.
- Fig. 1: a stray "I_6/7" label (sec. 2.1).
- Sec. 3.1 says first-order resonances m/(m+1) "m = 2, 3, ..., 7" are considered, while the family statements use
  1 <= m <= 6 and 1 <= m <= 7 (they include the 1/2 side of curve 1). Not an error in the results.

## 6. Citation mining (29 references)

Checked with `ls cyclers_pdf/papers | grep -i <author>`, `grep -i <author>` in CORPUS_INDEX.md and in the wanted list
(current numbering). Grep hits were checked one by one: "Duncan" is Duncan, Levison & Lee 1998 (SyMBA), not Duncan,
Levison & Budd 1995; "Poincaré", "Press" and "Roy" hits are other papers.

| reference | status |
|---|---|
| Broucke 1969a, JPL TR 32-1360 | HELD (`broucke-1969-...-jpl-tr-32-1360-...`) |
| Broucke 1969b, AIAA J. 7:1003 | HELD (`broucke-1969b-...`) |
| Hadjidemetriou 1993, CMDA 56:201-219 | HELD as a chapter of `dvorak-henrard-eds-1993-...` |
| Hénon 1997, LNP m52 | HELD (`henon-1997-...`) |
| Kotoulas & Voyatzis 2004, CMDA 88:343-363 | not held; **wanted row 56** |
| Voyatzis, Kotoulas & Hadjidemetriou 2005, CMDA 91:191-202 | not held; **wanted row 56**. The 1/2, 1/3, 1/4 symmetric and asymmetric families at the same mu; the closest companion paper |
| Guillaume 1974, A&A 37:209-218 ("families of symmetric periodic orbits ... when the perturbing mass is small") | not held (Guillaume 1973, 1975, 1975a are held); not wanted. **Candidate** (Tier C): the first-order "gap" theory behind Fig. 1 |
| Kotoulas & Hadjidemetriou 2002, Earth Moon Planets 91:63-93 (2/3 and 3/4 exterior resonances, Neptune) | not held; not wanted. **Candidate** (Tier C): its families are the 2/3 and 3/4 curves here |
| Kotoulas & Voyatzis 2005, IAU Coll. 197, pp.349-354 (3D BPs) | not held; not wanted. Low priority |
| Varadi 1999, AJ 118:2526 (3:2 periodic orbits for various mass ratios) | not held; not wanted. Low priority (one resonance, mass ratio varied: a possible exterior "in mu" source) |
| Hadjidemetriou 1988 CMDA 43:371; 1992 CMDA 53:151 | not held; not wanted. Background |
| Celletti, Chessa, Hadjidemetriou & Valsecchi 2002, CMDA 83:239 | not held; not wanted. Background |
| Haghighipour et al. 2003, ApJ 596:1332; Beaugé et al. 2003, ApJ 593:1124 | not held; not wanted. Extrasolar / asymmetric background |
| Benet et al. 1999; Berry 1978; Poincaré 1892; Press et al. 1992; Roy 1982 | not held; not wanted. Textbooks and general dynamics |
| Malhotra 1996; Morbidelli 1999; Morbidelli et al. 2003; Gallardo & Ferraz-Mello 1998; Duncan et al. 1995; Knežević et al. 1991; Torbett & Smoluchowski 1990; Jewitt 1999; Jewitt & Luu 1993 | not held; not wanted. Kuiper-belt background, not in project scope |

Proposals for the wanted list (PROPOSALS only):
- Close row 22 with this file.
- Update R14's "wanted row 52" and the Varin-Bruno 2009 digest's "row 56" pointers to row 22 (or to "held").
- Add Guillaume 1974 and Kotoulas & Hadjidemetriou 2002 as Tier C. DOIs not checked here.

*Filed as `cyclers_pdf/papers/voyatzis-kotoulas-2005-planar-periodic-orbits-exterior-resonances-neptune-pss-53-1189-doi-10.1016-j.pss.2005.05.001.pdf`. Check scripts, outputs and other files named above are filed beside it as `cyclers_pdf/papers/<pdf stem>-<file name>`. Table transcription: `data/sources/voyatzis-kotoulas-2005-tables.yaml`.*

*Wanted-list row numbers in this digest are the pre-batch-36 numbering; the list was renumbered in batch 36.*
