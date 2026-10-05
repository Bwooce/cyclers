# Digest: Minovitch "Letter/Documents-Gallery" part 1 (JPL TM 312-118 and TM 312-130, 1961) (#966)

Source: an owner-supplied PDF, "Letter/Documents-Gallery", M. A. Minovitch's reproductions of his own
original documents in date order. 72 pages, image-only JPEG scans at 150 ppi, no text layer.
- md5 of the scan as supplied: a19adf60365e39bb83e99a7cfae30f0a. The upload name was `LD-1.pdf`, so this
  is probably part 1 of a multi-part gallery.
- Filed (OCR'd copy) as `cyclers_pdf/papers/minovitch-1961-letters-documents-gallery-part1-jpl-tm-312-118-tm-312-130.pdf`.
- The name is spelled "Minovich" on all three 1961 documents; the later published spelling is "Minovitch".

## 0. Verdicts

- **TR 32-464 is NOT in this file, not even in part.** The file holds a telegram and two 1961 JPL
  technical memos. The wanted-list row for TR 32-464 stays open.
- The file contains **no numerical trajectory tables, no dates of real launch windows, no V-infinity
  values and no worked example.** Both memos are theory and algorithm only. There is nothing here to use as
  a numerical positive control.
- The value of the file is (a) priority and method history: TM 312-130 (23 August 1961) states the
  multi-planet free-fall tour and its patched-conic solution procedure, with the V-infinity-magnitude
  match at each flyby as the closing condition; and (b) two small method points that agree with choices
  the project has already made (sec. 4).

## 1. Method of reading

- Two OCR engines, both run on all 72 pages, rendered at 300 dpi greyscale:
  - Engine 1: tesseract 5 (`--psm 4`).
  - Engine 2: macOS Vision, `VNRecognizeTextRequest`, accurate level, language correction off.
- Agreement: tesseract found 7,004 words of three or more letters and Vision found 6,641. The median
  per-page word-sequence similarity was 0.81, and 36 of 72 pages were at 0.8 or above. The low pages are
  the telegram (p2, tesseract read nothing useful on the orange paper) and the equation-dense pages (p15,
  p33, p37, p40, p41, p53, p66).
- Numeric tokens: tesseract 1,887, Vision 2,706, and 1,364 matched as a multiset. Almost all the numeric
  tokens are equation numbers, subscripts and page numbers, because the memos have no tables. Both engines
  garble the mathematics badly (Vision less so on the typed text). The OCR text is useful only for search.
- **Every number and formula quoted below was read on the page image**, not taken from OCR. The pages
  read on the image are p1-4, p20-22, p25-26, p33, p36, p45-48, p63-64, p66-70 and p72. The remaining
  pages were read from the OCR text of both engines (enough to place them in the structure, not to quote
  formulas).

## 2. Page inventory

| PDF pages | Type | Doc id | Date | Title | Author |
|---|---|---|---|---|---|
| 1 | cover sheet | - | - | "Letter/Documents-Gallery (Most of these reproductions are reproductions of the original documents in chronological order.)" | - |
| 2 | Western Union telegram | LA167 | 14 Apr 1961 | Offer of summer employment at JPL as research engineer at $690.00 a month, made by Dr William G. Melbourne, start date 9 June 1961 | A. E. Locke, JPL employment supervisor, to Michael A. Minovich Jr, Dykstra Hall, Los Angeles |
| 3-23 | JPL technical memo | TM #312-118 | 11 Jul 1961 | "An Alternative Method for the Determination of Elliptic and Hyperbolic Trajectories" | M. Minovich |
| 24-72 | JPL technical memo | TM #312-130 | 23 Aug 1961 | "A Method For Determining Interplanetary Free-Fall Reconnaissance Trajectories" | M. A. Minovich |

- TM 312-118 is memo pages -1- to -21-, complete. Distribution: Section 312 Engineers, J. F. Scott,
  W. Scholey.
- TM 312-130 is memo pages -1- to -47-, plus -12a- (PDF p36) and -27a- (PDF p52). It ends with the
  references on -47- (PDF p72) and the typist mark "MM:ls". It appears complete. Same distribution.
- PDF metadata: created 6 January 2006 by an HP scanner driver.

## 3. Contents

### 3.1 TM 312-118 (11 July 1961): Lambert's problem by the semi-major axis

- Problem: given two position vectors, the focus and the flight time T, find the conic (p3).
- Construction: the second (vacant) focus lies on circles of radius 2a - r1 and 2a - r2 about P and Q.
  There are two candidate vacant foci F* and F~*, so two ellipses. The minimum-energy ellipse has
  2a_m = s = (r1 + r2 + c)/2 (p4).
- Flight time: Lagrange's form, eq. (1)-(4) (p4), with x1 = 1 - s/a and x2 = 1 - (s - c)/a:
  T = f(a) = sqrt(a^3/mu) { sqrt(1 - x2^2) + asin x2 - sqrt(1 - x1^2) - asin x1 } and
  T~ = f~(a) = sqrt(a^3/mu) { pi + sqrt(1 - x2^2) + asin x2 + sqrt(1 - x1^2) + asin x1 }.
  It cites Battin, "The Determination of Round-Trip Planetary Reconnaissance Trajectories", ARS Journal /
  Space Sciences, pp. 550-52 (p4; this is Battin 1959, already on the wanted list, row 52).
- Solution: Newton iteration a_{k+1} = a_k - (F(a_k) - T)/F'(a_k). Convergence is guaranteed by the
  convexity of the lower half of the T(a) curve. The error is E_{k+1} ~ (1/2) E_k^2 |F''/F'| (p22, step v).
- Hyperbolic branches h(a) and h~(a) for T < T~_0, with limits as a goes to 0 and to infinity (p13-19).
  The p20 figure ("Graph of T vs. a") shows the five branches: short-time and long-time elliptical,
  short-time hyperbolic, and long-time hyperbolic for cases 1 and 2.
- Summary procedure, steps (i)-(vii) (p21-23), with eccentricity formulas for each branch.
- The "round trip" on p20 is only the degenerate radial limit c -> 0, T~_0 = (4/3) sqrt(r^3/(2 mu)). It
  is not a planetary round trip.
- No multi-revolution branch. No numbers.

### 3.2 TM 312-130 (23 August 1961): free-fall reconnaissance trajectories to one planet and to N planets

Abstract (p24): "After solving the trajectory problem to one planet and back the more general problem of
determining a free-fall reconnaissance trajectory to N planets before returning to the launch planet will
be solved. No assumptions will be made as to the geometry of the solar system; indeed, it will not matter
how eccentric the planets orbits are or how much their planes of motion differ from each other. ... As far
as the author knows the method and results are new."

Model (p24-25): a three-dimensional patched conic.
- I. Inside the sphere of influence only the planet acts; outside it only the Sun acts.
- II. The planet's velocity is constant during the encounter, equal to its value at closest approach.
- III. The angle of the asymptote equals the angle of the hyperbola at the sphere-of-influence crossing
  (phi_inf = phi_rho*).
- Sphere of influence: rho* = (m/M)^(2/5) c, cited to Tisserand, Traite de Mecanique Celeste, Tome IV,
  p. 198 (1889) (p25, p72).

One-planet round trip (p34-48):
- Only the launch time t_0 and the closest-approach time t_CA are prescribed (p34).
- The outbound leg is a Lambert arc by TM 312-118.
- The return leg is tabulated as a function of the unknown return time t_3: a_3(t_3) from eq. (22)
  (p45). The return may make k complete circuits of the Sun first (eq. 19 and 20, p44). The k = 1 case
  adds 2 pi sqrt(a_3^3/mu_s) to the right side of eq. (22) (p48).
- **Closing condition, eq. (25)-(26) (p46):** because rho_1 = rho_2, the planet-relative speeds are equal
  (v'_1^2 = v'_2^2), so the heliocentric energy change is
  v_2^2 - v_1^2 = 2 V_Q . (V_2 - V_1), that is
  mu_s (1/a_1 - 1/a_3(t_3)) = 2 V_Q . [V_2(t_3) - V_1].
  The unknown t_3 is found by comparing the two tables and picking the entry where they agree.
  This is the V-infinity-magnitude match at a flyby, with the next encounter date as the unknown.
- Flyby geometry (p47-48): v'^2 from (29) as the average of the two sides; the hyperbola semi-major axis
  a_2 = rho* mu_Q / (v'^2 rho* - 2 mu_Q) ~ mu_Q / v'^2 (eq. 30); the eccentricity
  e_2 = sqrt( 2 v'_1 v'_2 / (v'_1 v'_2 - V'_1 . V'_2) ) (eq. 31); and the closest-approach altitude
  d = a_2 (e_2 - 1) - R_Q (p48).
  - I checked these against the standard forms. Eq. (30) is vis-viva at the sphere-of-influence radius.
    Eq. (31) is the same as sin(delta/2) = 1/e, with the turn angle delta = pi - 2 phi.
  - d is the only feasibility gate. There is no explicit turn-angle limit apart from d > 0.
- **Branch logic (p48):** if d < 0, the vehicle "cannot return to the launch planet on a short-time
  elliptical path without first making at least one complete circuit of the sun", so set k = 1. Or try
  the long-time return arc with k = 0.
- **Four branches (p63):** for one launch time and one energy there are "in most cases" four distinct
  round trips with k = 0 (short or long departing arc, times short or long returning arc).
- The memo also computes the position, velocity and time along each arc at N equal steps, for planning
  observations (p51-62).

N-planet tour (p63-68):
- The problem (p63-64): launch at t_02, prescribed closest approach to the first planet at t_1CA, then
  visit N-1 more planets in a prescribed order and return to the launch planet.
- **The example (p64):** "at t_02 the vehicle leaves the 'center' of the earth and makes a closest
  approach to the first planet Venus at time t_1CA. It then proceeds to visit the remaining N-1 planets in
  the following order: Mars, Earth, Saturn, Pluto, Jupiter, Earth." That is E-V-M-E-S-P-J-E. Earth is
  flown by in mid-tour. The sequence is stated, not computed.
- **The algorithm (p66-68), steps (i)-(xxxiii):**
  - Only the first closest-approach date is prescribed. Each later closest-approach date t_jCA is solved
    by eq. (26) at planet j.
  - Branch order per leg: short-time arc with k = 0, then long-time with k = 0, then the first leg is
    switched to long-time, then k = 1, 2, ... (step xxvii).
  - After each flyby, compute d_j. If all d_j > 0, the tour is accepted. At the first d_j < 0, "the next
    best value of t_1CA is calculated and the process is continued" (p68).
  - In effect this is a sequential date-marching shooter over a one-parameter family (t_1CA), for a fixed
    launch date and planet sequence.
- Conclusion (p69), in full: "In conclusion, we notice the remarkable fact that if E is the total
  heliocentric energy of a departing free-fall reconnaissance vehicle to one planet and back, it may be
  possible to send the vehicle on a trajectory which will take it to N-1 more planets before returning to
  its launch planet without any appreciable change in E."
- Appendix (p70-71): successive approximation to refine the patched-conic solution toward the exact one
  (one planet; N > 1 "very similar" and not written out).

Not in either memo:
- No periodic or repeating trajectory, and no cycle that repeats. The tour returns to the launch planet
  once.
- No Earth-Venus or Earth-Mars resonance or period ratios.
- No numbers for any planet.

## 4. Points that bear on the project's methods

1. **Closing condition at a massive flyby.** TM 312-130 eq. (26) closes each flyby on the energy (that is,
   V-infinity magnitude) match, with the next encounter date as the unknown. This is the same unknown and
   residual as the #942 corrector (Hollister & Menning 1970, pp. 1194-1195). Minovitch states it 8 years
   earlier, for non-periodic tours. It does not change the #942 code.
2. **Branch enumeration.** The memo enumerates short-time and long-time arcs and k complete solar
   circuits for every leg, and says four branches exist "in most cases" (p63). This supports the #864
   finding that multi-revolution Lambert legs must be in the search (the June negatives that lacked them
   were void).
3. **Periapsis gate only.** The feasibility test is d > 0 (above the surface). It has no minimum altitude
   margin and no explicit turn limit. Our #888/#937 turn gate is stricter. Do not use this memo as a source
   for a turn limit.
4. **Sphere of influence.** rho* = (m/M)^(2/5) c (Laplace's radius, via Tisserand). Assumption III is
   the usual approximation that the asymptote direction is the direction at the sphere-of-influence
   crossing.

## 5. Ideas stated but not computed (candidate routes)

- The N-planet tour with no appreciable change in heliocentric energy (p69). This is the gravity-assist
  idea, not a cycler. It adds nothing that the #938 routes do not already cover.
- The example E-V-M-E-S-P-J-E (p64). It includes an Earth-Venus-Mars-Earth segment, which is the
  topology of the Earth-Venus-Mars round trips that TR 32-464 later computed (as Menning 1968 cites it).
  The memo has no numbers for it.

## 6. Useful to the project (ranked)

1. **Get the later parts of the gallery (new request, high value, low cost).** This file is named LD-1
   and is in date order, and it stops in August 1961. TR 32-464 (1963) is probably in a later part. Ask
   the owner whether LD-2 and later parts exist. If they do, rerun this pipeline on them. This is the
   shortest path to wanted-list row 1, which gates the pre-Rall Earth-Venus-Mars collision check for
   #938 R1 cells (a) and (c).
2. **#942 / #943 (R1, X1) method provenance.** Cite TM 312-130 eq. (26) (p46) and the date-marching
   algorithm (p66-68) as the earliest statement of the V-infinity-magnitude-match residual with
   encounter dates as unknowns. This is a note-level citation only. No code or gate change.
3. **#864 lesson support.** p48 and p63 (k circuits; four short/long branches) are a primary source for
   "include multi-revolution and long-way legs". Add it as a citation where the #864 lesson is recorded.
   No new task.
4. **Priority date for the novelty record.** 23 August 1961 is the dated JPL memo for the N-planet
   free-fall tour (E-V-M-E-... example). Any project text that dates multi-flyby round trips should use
   this date for the method, and TR 32-464 (1963) for the first numbers.
5. **Not useful:** TM 312-118 is a textbook Lambert solver (single revolution, Newton in a). The project
   already has better solvers. Battin 1959 (cited on p4) is already on the wanted list (row 52).

No positive control, no new cycler topology and no new #938 route comes from this file.

## 7. Files

- OCR text (both engines, per page and concatenated) and the inventory are in the session scratch
  directory `minovitch-ocr-opus/` (`all_tess.txt`, `all_vision.txt`, `inventory.md`, `stats.txt`). They
  are not committed. The OCR'd PDF in the corpus carries a tesseract text layer for search.
