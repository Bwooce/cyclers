# #942 / #943: two-working-body cycler generator (Hollister date corrector), controls and enumeration

Status: BUILD AND CONTROLS DONE; ENUMERATION PENDING the lead's launch (sec. 6). No catalogue writes.
Novelty language: every enumerated cycler is "candidate, pending collision check" until sec. 7 clears it.

## 1. What was built

- `src/cyclerfinder/search/two_working_body.py`: the corrector.
  - Legs:
    - `LambertLeg`: transfers and same-body generic or "symmetric" returns.
    - `ResonantLeg`: an n:m full-revolution return. Its duration is fixed and its V-infinity direction
      is free on a circle.
    - `HalfRevLeg`: an n-pi return. It is analytic and has discrete directions.
  - Unknowns: the start dates of the Lambert legs. Residuals: the V-infinity magnitude difference at
    each massive junction, or the full vector difference at a massless target. This is Hollister &
    Menning 1970, pp.1194-1195. A square system, solved by least squares.
  - Free directions: chosen by minimax of the demanded/available turn ratio over the flybys of a
    block. Ties at inner flybys are broken by spreading the turn evenly. The tie-break cannot change
    the largest ratio, so it cannot change a gate verdict.
  - Every flyby goes through the #888/#937 gate (`verify/turn_gate.py`), on inertial-axis vectors.
    Wrap vectors are rotated in the circular model.
  - Near-180-degree rejection: PROVISIONAL threshold of 175 deg (sec. 8). The largest demanded turn is
    stored, so the number can change without a rerun.
  - Systems:
    - `CircularSystem`: the ideal model.
    - `MeanElementSystem`: Standish-Williams J2000 elements, inclined and elliptic, with optional
      commensurate periods.
  - Independent checks:
    - `encounter_self_consistency`: each leg is re-propagated by its own Kepler step, not the Lambert
      solver, and must land on its body to well inside the sphere of influence (#480 rule).
    - `kepler_step`: an eccentric-anomaly propagator. The core `propagate` fails to converge on some
      ordinary heliocentric arcs; this was reported to the lead with a reproducer.
- `src/cyclerfinder/search/two_working_body_enum.py`:
  - structures (cycle templates).
  - multi-start seeding spread across basins.
  - de-duplication.
  - assessment (gate at the project floor and at H&M's 1.1 radii, extents, sphere of influence).
  - placement of a massless target.
- `src/cyclerfinder/search/hollister_menning_1970.py`: the Table 3 loader, with transcription and print
  fixes (sec. 3.1).
- Scripts:
  - `scripts/run_942_hm_control.py`: the pre-registered H&M control.
  - `scripts/run_942_enumerate.py`: the sharded, resumable driver.
  - `scripts/analyse_942_enumeration.py`: merges physically identical zeros and mirror twins.
- Tests: `tests/search/test_two_working_body.py`.

## 2. Bend capacity ("massive enough to bend", #943 brief)

Available bend (deg) at the project floor. GM and radius are from R-S 2009 Table 2 for the moons and
from the registry for the planets.

| Body | Floor (km) | V_inf 1 | 2 | 3 | 4 | 5 | 6 km/s |
|---|---|---|---|---|---|---|---|
| Europa | 100 | 82.4 | 38.0 | 20.3 | 12.3 | 8.2 | 5.8 |
| Ganymede | 100 | 103.1 | 56.7 | 33.3 | 21.2 | 14.5 | 10.5 |
| Callisto | 200 | 94.4 | 48.1 | 27.1 | 16.9 | 11.4 | 8.1 |

| Body | Floor (km) | V_inf 3 | 4 | 5 | 6 | 7 | 8 km/s |
|---|---|---|---|---|---|---|---|
| Earth | 200 | 121.1 | 104.6 | 90.1 | 77.7 | 67.1 | 58.2 |
| Venus | 300 | 116.5 | 99.2 | 84.4 | 71.9 | 61.4 | 52.7 |
| Mars | 200 | 69.4 | 50.5 | 37.6 | 28.8 | 22.6 | 18.1 |

At the R-S V-infinities (1.7-4.3 km/s), Callisto and Europa bend tens of degrees. Rall's condition holds.

## 3. Positive control 1: Hollister & Menning 1970 Table 3 (pre-registered; FAILED as registered)

### 3.1 Source audit
- Table 3 is the INCLINED-ELLIPTIC 16-yr solution set on real dates (p.1194-1195). It is not
  circular-coplanar. The circular 3.2-yr orbits (Hollister 1969, not held) have no printed numbers.
  The 15 catalogue rows carry `model_assumption: circular-coplanar`; that looks like a mislabel.
  Reported, not edited.
- Transcription errors in `data/sources/hollister-menning-1970-table3.yaml`, checked on 450-dpi page
  images:
  - orbit 1 row 12 is "E 3163" (the YAML has V).
  - orbit 6 row 3 is "995" (the YAML has 993).
- Print errors in the paper (225-d Venus step broken):
  - orbit 2 "5715" (5815).
  - orbit 4 "2573" (2583).
  - orbit 5 "5662" (5652).
  - orbit 6 "2585, 2810" (2360, 2585).
  - orbit 8 "4870" (5870).
  - orbit 5 "977" (inferred 1017: a 265-d first gap gives an 8.9 km/s residual at the printed dates).
- Paper-internal: the text puts the 1.16-radius Earth pass in orbit 5. The table has it in orbit 4.
- 25 printed rows have a (V_r, theta, Rmin) triple that disagrees with r_p = mu/v^2 (1/sin(theta/2) - 1)
  by more than 10 %. In several of them our theta reproduces THEIR Rmin.

### 3.2 Results
Pre-registered criteria are in the header of `scripts/run_942_hm_control.py`:
- dates within 3 d.
- V_r within 0.005 EMOS.
- theta within 3 deg.
- Rmin within 10 %.
- for >= 90 % of encounters.

Stage results:
- Stage F (V_r at the printed dates): failed on all 15 orbits (0.32-0.78 within 0.006 EMOS).
- Stage C, real periods: 0/15. Dates slip about 0.3 d per Venus full revolution (224.70 d against the
  printed 225 d).
- Stage C on H&M's stated "exact periodicity" model, after amendments 1-2 (post hoc, labelled): 4/15
  pass. Earth 5844/16 d and Venus 5844/26 d; anchor JD 2443363.

| Orbit | Earth/Venus returns | Exact zero | Match | Match, source-consistent rows | Note |
|---|---|---|---|---|---|
| 1 | all full-rev | yes | 0.92 | 1.00 | Table 2 1H seed and perturbed seed reach the same zero |
| 2 | Venus FS | yes | 0.92 | 0.92 | |
| 3 | Earth SY x5 | yes | 0.24 | 0.32 | 7 zeros within 35 d; best 0.28 |
| 4 | Earth SY x1 | yes | 0.64 | 0.73 | |
| 5 | Earth SY x2 | yes | 0.44 | 0.48 | |
| 6 | Earth SY x2 | no (fold) | 0.16 | 0.16 | LS minimum 0.0057 EMOS, just above H&M's 0.005 |
| 7 | Earth SY x3 | yes | 0.44 | 0.50 | |
| 8 | Earth SY x4 | yes | 0.76 | 0.76 | |
| 9 | Venus S x1 | no (fold) | 0.92 | 0.96 | LS minimum 0.0015 EMOS, inside H&M's tolerance |
| 10 | Venus S x2 | yes | 0.84 | 0.91 | |
| 11 | Venus S x3 | yes | 0.96 | 0.96 | |
| 12 | Venus S x4 | yes | 0.96 | 0.96 | |
| 13 | Venus S + Earth SY | yes | 0.88 | 0.88 | |
| 14 | | yes | 0.80 | 0.91 | |
| 15 | | yes | 0.88 | 0.96 | |

Reading:
- The family with Earth full-revolution returns (orbits 1, 2, 9-15) reproduces.
- The Earth-symmetric-return family (orbits 3-8, H&M's "sequential modification" set) does not. Our
  model has nearby zeros that differ by about 0.01 EMOS.
- The homotopy from the printed dates folds for orbits 6 and 9.
- The control is FROZEN. A further model change chosen by match would be tuning.
- The open discrepancy goes to the lead and owner.

### 3.3 Amendment 3, 2026-10-05, written and committed BEFORE the re-run

Context: the table was corrected to the print (commit 259efc0d, 30 cells). The re-run with
amendments 1-2 gave 6/15 on the exact-periodicity model and 0/15 on real periods.

Then I read Menning 1968 (MIT S.M. thesis): ch. 1, ch. 3, ch. 4 pp.17-25, ch. 5 pp.29-31, p.33, p.44 and
A-1/A-2. I also read Hollister 1969 (JSR 6(4) pp.366-369) and the H&M 1970 text (pp.1194-1195).

**Which model Table 3 was computed in (sources):**
- H&M 1970 p.1194: "inclined elliptic case"; "the error made by assuming exact periodicity of the solar
  system is of the same order of magnitude as [the patched conic model]".
- Menning p.30: "For the inclined elliptic case, Earth and Venus repeat their absolute orientation to
  within several degrees every 16 years" (his ref. 12, Gillespie & Ross 1966).
- Hollister 1969 p.366: the 32-yr (and 8-yr) resonances hold "to within a few degrees". The general
  case is "eccentric and inclined".
- No source gives the planetary elements or periods used.

Reading: real (non-commensurate) inclined-elliptic planet orbits. The exact periodicity is assumed only
when the 16-yr cycle is closed (last date = first + 5844 d). So our REAL-PERIOD model (stage 1) is the
sourced one, and its result is the honest headline. The exact-periodicity model (stage 2) is our choice.

Unexplained data point, not a source: 129 of the 130 printed Venus full-revolution steps are exactly
225 d (the other is 224 d). A 224.70-d period would print about 30 % 224-d steps, and 224.77 d about
23 %. So their effective Venus full-revolution time was close to 225.0 d. Neither model explains this.
I do NOT adopt it, because it was inferred from the table being matched.

**Every difference between our model and theirs:**

| # | Item | Theirs (cite) | Ours before amendment 3 | Amendment 3 |
|---|---|---|---|---|
| D1 | Planet elements | inclined, elliptic; values not stated (Menning p.30; Hollister 1969 p.366) | Standish & Williams J2000 mean elements, fixed | unchanged (no source) |
| D2 | Periodicity | real orbits, cycle closed at +5844 d (H&M p.1194; Menning p.30) | stage 1 real; stage 2 commensurate 5844/16, 5844/26, anchor JD 2443363 | both still reported; stage 1 is the headline |
| D3 | Full-revolution circle | spacecraft heliocentric speed equals the planet's, \|v_sc\| = \|V_P\| (Menning p.18) | vis-viva with a from the model period (differs only in stage 2, by a few m/s) | \|v_sc\| = \|V_P\| for every 1:1 return |
| D4 | Turn selection | cone-vector rules, secs. 4.21 (one return) and 4.22 (two returns), pp.22-25 | global minimax plus an even-spread tie-break | Menning's rules for blocks of 1-2 full-revolution returns (`menning_block`); same largest turn, different inner turns |
| D5 | Symmetric return | conventional transfer plus one revolution, Eq. (3.6); the root that is not the planet's own orbit (pp.12-15) | 1-rev Lambert, branch chosen by the residual at the printed dates | unchanged; the chosen branch is never the planet-orbit root (that root has V_inf near 0) |
| D6 | Convergence | assumed at summed \|delta V_r\| = 0.005 EMOS (Menning p.33; H&M p.1195) | exact zero (max residual < 1e-9 km/s) | summed \|residual\| at the least-squares minimum <= 0.005 EMOS |
| D7 | Flyby constants (mu, radius, EMOS) | not stated | registry; EMOS 29.785 km/s | unchanged |
| D8 | Patched conic, hyperbola time ignored, planet-centre ellipses | Menning pp.4-5 | same | - |
| D9 | Acceptance floor | 1.1 radii (Menning p.6) | not part of the match test | - |

Pass rule unchanged: >= 0.90 of the 25 encounters within all four pre-registered tolerances (date 3 d,
V_r 0.005 EMOS, theta 3 deg, Rmin 10 %). D6 replaces only the "exact zero" requirement.

Reported side by side:
- raw.
- source-consistent rows only. These are rows whose own printed (V_r, theta, Rmin) agree with the
  flyby formula; orbit 14's three inconsistent rows can never pass.

One re-run per model (real periods = headline; exact periodicity = secondary).

### 3.4 Result of the amendment-3 re-run (2026-10-05, commit a04d4696, corrected table)

Criterion: >= 0.90 of the encounters within the pre-registered tolerances, and summed |dV_r| <= 0.005
EMOS (Menning p.33). Cells give raw / source-consistent-rows match, [summed residual in EMOS].

**HEADLINE, real periods (the sourced model): 0/15 pass raw; 1/15 on source-consistent rows (orbit 14).**

| Orbit | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 | 15 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| raw | 0.76 | 0.64 | 0.52 | 0.00 | 0.36 | 0.00 | 0.48 | 0.48 | 0.84 | 0.84 | 0.72 | 0.60 | 0.84 | 0.80 | 0.68 |
| consistent | 0.76 | 0.64 | 0.52 | 0.00 | 0.39 | 0.00 | 0.48 | 0.48 | 0.84 | 0.84 | 0.72 | 0.60 | 0.84 | 0.91 | 0.71 |

Orbit 6's residual is 0.0056 EMOS and orbit 9's is 0.0037; every other orbit is an exact zero.

**Our variant, exact periodicity: 7/15 pass raw (orbits 1, 2, 9, 11, 12, 13, 15); 8/15 on
source-consistent rows (adds orbit 14).**

| Orbit | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 | 15 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| raw | 1.00 | 0.92 | 0.48 | 0.64 | 0.48 | 0.16 | 0.60 | 0.84 | 0.96 | 0.84 | 0.96 | 0.96 | 1.00 | 0.88 | 0.96 |
| consistent | 1.00 | 0.92 | 0.48 | 0.73 | 0.52 | 0.16 | 0.60 | 0.84 | 0.96 | 0.84 | 0.96 | 0.96 | 1.00 | 1.00 | 1.00 |

Orbit 6's residual is 0.0057 EMOS (fails H&M's tolerance) and orbit 9's is 0.0015 (passes it).

Comparison with earlier runs, exact-periodicity model:
- original pre-registered rule: 4/15.
- after the table correction: 6/15.
- after amendment 3: 7/15.
- On real periods: 0/15 throughout.

The Earth-symmetric-return family (orbits 3-8) fails in both models, with nearby zeros whose V_inf is
about 0.01 EMOS off. The full-revolution family passes only on the exact-periodicity variant.

### 3.5 Reproduction ladder (owner ruling, 2026-10-05): step 1, the planet model (sourced)

Hollister 1969 (JSR 6(4)) states the general-case procedure on p.368 (text layer, read):
- "The periodic orbits of the last section can be used as starting points, but they now take 16 yr
  before repeating exactly."
- Steps: "1) Use the circular, coplanar analysis to establish the 10 t_i for the case of zero
  eccentricity and inclination. 2) Increase the eccentricity and inclination to their actual values."

Its circular model (p.366) is "Earth makes 32 revolutions of the sun in 32 yr ... Venus makes 52
revolutions ... in 32 yr", so the periods are 1 yr and 8/13 yr.

So the inclined-elliptic orbits keep the COMMENSURATE periods and repeat EXACTLY in 16 yr. With 5844 d =
16 yr (the Table 3 footnote), that is Earth 365.25 d and Venus 224.769 d, the same values as our
"exact-periodicity" model.

**Amendment 4:** the exact-periodicity model is the SOURCED model, and its result is the headline:
- 7/15 raw.
- 8/15 on source-consistent rows (sec. 3.4).

This supersedes sec. 3.3's reading. Sec. 3.3 had taken H&M p.1194 "assuming exact periodicity" and
Menning p.30 "within several degrees" to mean real periods with a closed cycle. Hollister's procedure is
the more specific statement. The real-period run (0/15) is kept as the non-sourced variant.

What stays unsourced (D1):
- the planets' phases (mean longitudes) and the orientation of the ellipses. Ours: Standish & Williams
  J2000 elements, with mean longitudes exact at JD 2443363.
- the e and i "actual values". Ours: Standish & Williams.

The 129/130 printed Venus steps of 225 d are still not explained: 224.769 d would print about 23 % 224-d
steps.

### 3.6 Step 2: orbits 3-8 (the Earth-symmetric family) against Menning's method

Checked against Menning chs. 1, 3, 4, 5 and H&M pp.1194-1195:
- Date convention: encounter dates are planet-centre crossings, and hyperbola time is ignored (Menning
  pp.4-5). Ours is the same.
- Symmetric-return root: Eq. (3.6), a conventional transfer plus one revolution. The transfer root is
  the one that is not the planet's own orbit (pp.12-15). Ours is a 1-rev Lambert whose other root has
  V_inf 0.000 (the planet orbit); we use the non-planet root. Same.
- Transfer legs: a conventional (0-rev) Lambert, Eq. (3.1). Ours is the same.
- Turn rule: secs. 4.21-4.22, applied after the dates are solved, so it cannot change the zero.
  Adopted in amendment 3.

No method difference was found. Diagnosis at the printed dates, orbit 3, exact-periodicity model:
- The symmetric legs reproduce the printed V_r to 0.002-0.006 EMOS (they are long, 490 d, and not
  sensitive).
- The short transfers (78-135 d) miss by up to 0.019 EMOS (e.g. V->E 3758->3885: 0.186 against 0.204).
- Orbit 1's transfers (155-223 d) miss by at most 0.012.

INFERRED (not sourced): the Earth-symmetric family has the shortest transfers. Those are the most
sensitive to planet phase, which is exactly the unsourced item D1. A planet-phase difference of a degree
or so between our ephemeris and theirs would act on that family first. Testing this by moving our phase
anchor would be tuning, so it is not done.

### 3.7 Step 3: Hollister 1969 Table 1 V_inf DIRECTIONS (PRE-REGISTERED 2026-10-05, before any computation)

Source: Hollister 1969 p.368, Table 1, "Periodic orbit I". It lists 21 transfer endpoints (LV/AR) with
JD - 2440000, V (EMOS), Ang and Elev (deg). Footnotes, read from the page image:
- "The angle is in the orbital plane clockwise from the circumferential direction."
- "The elevation is positive when above the orbital plane."

The dates are Menning's 1H = Table 3 orbit 1, so the test uses our orbit-1 solution.

Pre-registered conventions:
- "orbital plane" = the plane of the encountered planet's orbit (normal r x v of the planet).
- "circumferential" = the unit vector in that plane perpendicular to the planet's radius vector, in the
  direction of motion.
- "Clockwise" is ambiguous in sign. Convention A measures toward the outward radial; convention B toward
  the inward radial. Both are reported; the control passes if EITHER passes. That is one declared bit
  of freedom.

Pass rule: at >= 90 % of the 21 events,
- |dV| <= 0.005 EMOS,
- |dAng| <= 10 deg (wrapped),
- |dElev| <= 10 deg.

Models:
- headline: the sourced exact-periodicity model, at the amendment-3 converged orbit-1 zero (also
  reported at the printed dates).
- secondary: real periods.

The vectors are the Lambert-leg end V_inf, which carry no free direction. So this tests the transfer
geometry and the ephemeris, independent of the turn rule.

### 3.8 Step 3 result: the direction control PASSES (data/942_hollister1969_table1_directions.json)

| Model / point | Convention A (toward outward radial) | Convention B |
|---|---|---|
| exact periodicity (sourced), converged orbit 1 | **20/21 = 0.95, PASS**; median abs dAng 0.7 deg, dElev 0.5 deg | 0.05 (median dAng 51 deg) |
| real periods, converged orbit 1 | 0.95, PASS; median 1.1 / 0.4 deg | 0.05 |
| exact periodicity, at the printed dates | 0.33 (median 1.7 / 1.7 deg) | 0.00 |

Reading the result:
- The angle convention is not in doubt. A fits to about 1 deg; B misses by about 50 deg.
- The one failing event is AR E 3935 (V 0.1562 against 0.156; angle within 0.2 deg). Its elevation
  is -45 deg against a printed +46: the same magnitude with the opposite sign. It is the only sign
  disagreement in 21 events. This may be a sign slip in the print, offered with respect; it was not
  adjusted.
- Every other converged event agrees to within 4 deg (real periods) or 1.5 deg (exact periodicity).
- The converged zero fits much better than the printed dates. The corrector moves the dates by under
  2.2 d onto Hollister's own solution.

So for orbit I (= Table 3 orbit 1) the transfer GEOMETRY, and not only the magnitudes, reproduces an
independent published solution.

### 3.9 Amendment 5 (INFERRED SOURCE), pre-registered and committed before its runs

Source: Rall, C. S. (1969), MIT Sc.D. thesis TE-34, filed as
`hollister-rall-1970-periodic-orbits-NASA-CR.pdf`.
- p.10: "Hollister and Menning, however, took care of this periodicity problem by modeling the planets'
  orbits as truly periodic." This supports amendment 4.
- p.129: "The ephemerides are based on the mean orbital elements of 1960."
- p.136 (program listing, page image), Earth and Venus:
  - a = 1.0 and 0.723332 AU.
  - e = 0.016726 and 0.006793.
  - PER = 365.25636 and 224.7008 d.
  - GFP (true longitude of perihelion, 1960 equinox) = 102.25253 and 131.00831 deg.
  - 1960 perihelion dates (JD - 2440000) = 2.124962 - 3065 and -27.01776 - 3065.
- INFERRED: that H&M used the same numbers. Rall is Hollister's student in the same group, but no source
  says so.

The model (`hollister_menning_1970.RallElementSystem`):
- the Rall elements.
- GFP precessed +0.5588 deg (50.29 arcsec/yr x 40 yr) to J2000.
- mean longitude exact at the anchor (Rall's elements and real PER), then advanced at the truly
  periodic 365.25 / 224.769 d.
- Venus inclination and node from Standish & Williams J2000 (OUR CHOICE; Rall's setup lists none).
- Earth on the ecliptic.

Anchors, both pre-registered, one run each, both reported, neither picked as better:
- (a) the 1960 element epoch, JD 2436935.
- (b) JD 2443363 (as amendment 4).

At (a), Venus sits about 3.24 deg in longitude from the amendment-4 ephemeris; at (b), about 0.1 deg.

Everything else is unchanged from amendment 3: Menning turn rules, |v_sc| = |V_P|, the 0.005-EMOS
summed-residual convergence, and the four tolerances with the 0.90 pass fraction.

THE HEADLINE STAYS AMENDMENT 4 (sec. 3.5) unless a source ties H&M to Rall's elements. Amendment 5 is
evidence, not a reproduction claim.

Suspected print error (recorded, not adjusted): Hollister 1969 Table 1, AR E 3935, elevation +46. Ours
is -45, with matching speed and angle (sec. 3.8).

### 3.10 Amendment 5 results (commit 240b3b2c) and the end of the ladder

Cells give raw / source-consistent-rows match, [summed residual in EMOS]; P = pass.

| Orbit | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 | 15 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| (a) anchor 1960 | 1.00 P | 0.40 | 0.20 | 0.12 | 0.00 | 0.00 | 0.48 | 0.44 | 0.80 | 0.36 | 0.56 | 0.64 | 0.76 | 0.52/0.59 | 0.76/0.79 |
| (b) anchor JD 2443363 | 1.00 P | 0.96 P | 0.48 | 0.64/0.73 | 0.48/0.52 | 0.24 [0.0057] | 0.68 | 0.84 | 1.00 P [0.0014] | 0.84 | 0.96 P | 0.96 P | 1.00 P | 0.88/1.00 | 0.96/1.00 P |

Totals:
- Anchor (a): 1/15 raw, 1/15 on source-consistent rows.
- Anchor (b): 7/15 raw, 8/15 on source-consistent rows. This is the same set as amendment 4.

Reading:
- Rall's 1960 elements do NOT fix orbits 3-8 at either anchor.
- Anchored at their own 1960 epoch, the truly periodic Venus drifts about 3.2 deg by the 1970s, and
  almost every orbit fails.
- Anchored mid-span they reproduce amendment 4 to within a few hundredths. So Standish and Rall
  elements agree there; the e and perihelion differences do not matter.

The anchor (planet phase) is the one model item no source states, and the result depends strongly on it.
That fits the inference in sec. 3.6, but it does not reproduce orbits 3-8.

LADDER STOPPED here (lead's instruction). Headline (amendment 4, sourced): 7/15 raw and 8/15 on
source-consistent rows; orbits 3-8 not reproduced.

## 4. Positive control 2: recall in the circular model (production enumerator)

E-V, k = 2 (3.197 yr), production settings (n_phase 36, n_split 12, n_refine 40, basin-spread seeds):
- 1H (full-rev Earth, two full-rev Venus): 3 gate-passing zeros, e.g. V_inf E/V 2.99/3.19 and
  6.14/5.28 km/s.
- 2H (full-rev Earth, full-rev plus symmetric at Venus): 3 gate-passing zeros, e.g. 5.60/6.02 km/s.
- 3H (symmetric Earth, two full-rev Venus): 3 gate-passing zeros, e.g. 4.31/5.51 km/s. The symmetric
  return lasts 1.36 yr. H&M's sequential modification inserted a 1.37-yr one (p.1195); this is a
  consistency check only.

## 5. Positive controls 3-6: Russell & Strange (expected values from the papers' tables)

All four were run BLIND through the production enumerator. In each run the published cycler is the only
gate-passing physical cycler.

| Row | Cell | Ours | Published |
|---|---|---|---|
| GanCal#5 | gc1, k=3 | V_inf 3.238 / 3.34; alt 328.0 km; r 821,913-2,390,844 km | 3.24 / 3.34; 328; 821,915-2,390,844 |
| GanCal#1 | gc1, k=3, f(2:1) leg | 3.180 / 3.255; alt 247.3 km | 3.18 / 3.26; 247 |
| GanEur#43 | ge1, k=2 | 1.872 / 3.893; alt 8862 km; r 564,558-1,072,318 | 1.87 / 3.89; 8861; 564,558-1,072,330 |
| GanEur#5 | ge1, k=5 | 1.658 / 2.573; alt 1819.8 km | 1.66 / 2.57; 1819 |
| VenMar#45 | vm, k=2 | 8.220 / 12.965; alt 19,784 km; r_max 341,571,415 | 8.22 / 12.96; 19,784; 341,571,371 |

GanCal#1 contains a resonant f(2:1) leg, so its 247 km altitude also checks the minimax free-direction
convention against a published value. VenMar#45's published minimum distance (108,067,501 km) is the
conic's perihelion, which the leg does not pass; ours is the along-leg minimum.

## 6. Enumeration (PENDING)

Settings and pruning (first pass):
- ev and em: <= 2 returns per block, 0-rev transfers, generic 1-rev.
- vm: <= 3 Venus returns, transfers 0-1 rev, generic 1-2.
- gc and ge: <= 1 return per block, transfers 0-2 rev, generic 1-2.
- k = 1-3 in every cell.
- Catalogue: resonances (1:1, 2:1, 1:2, 3:2, 2:3); half-revs (1, 0, peri/apo), (3, 1, peri/apo);
  generic n-rev, low and high.

## 7. Literal-collision checks (to be completed per candidate)

Against:
- Hollister & Menning 1970 (catalogue rows hollister-menning-1970-ev-orbit-01..15).
- Jones 2017 VEM.
- R-S 2007 and 2009 rows: EurGan, GanEur, GanCal, GanIo, VenMar.
- Pisarevsky 2008: Table 4 class III, and the Fig. 14 class I.1 candidates for Earth-Mars.
- `data/empty_regions.jsonl`.
- `literature_check.py`.

## 8. Open items for the owner

- The H&M control discrepancy (orbits 3-8).
- The near-180 threshold (provisional 175 deg).
- The catalogue `model_assumption` of the 15 H&M rows.
- The two YAML transcription fixes.
