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

Hollister 1969 p.367, circular orbit I (sourced numbers): transfers of 0.485 yr and 0.6 rev, with V_inf
0.107 EMOS at Earth and 0.126 at Venus. Our corrector finds exactly this symmetric zero: both transfers
0.4846 yr, 216 deg, gate-passing. Its V_inf is 0.1008 EMOS at Earth and 0.1075 at Venus. An independent
shooting solve with the Kepler step gives the same numbers.

Readings checked (lead's request, 2026-10-06):
- (a) The inclined-elliptic orbit I does not fit. Hollister's own Table 1 averages 0.166 EMOS at Earth
  and 0.198 at Venus, with ranges 0.154-0.191 and 0.178-0.223.
- (b) No transfer angle at the stated 0.485-yr flight time gives both values. At 0.485 yr: 200 deg
  gives 0.120/0.101, 216 deg gives 0.101/0.108, 230 deg gives 0.099/0.150.
- (c) A scan over angle and flight time fits 0.107/0.126 only at geometries other than the stated one,
  e.g. 189.5 deg / 0.385 yr or 228.5 deg / 0.505 yr.
- (d) Rescaling by another planet's mean orbital speed does not give the pair either. Our Venus value of
  3.20 km/s equals the printed EARTH value of 3.2 km/s, but our Earth value of 3.00 km/s matches neither
  printed number.

Status: NOT REPRODUCED FROM THE STATED GEOMETRY. This is a suspected unit or print issue, offered with
respect. It is recorded only and does not affect the other controls.

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

## 6. Enumeration: PRE-REGISTRATION (owner ruling "search where controls pass", 2026-10-06)

Written and committed before any production run. Code: `scripts/run_942_enumerate.py` (sharded,
resumable, progress and ETA per structure) and `scripts/analyse_942_enumeration.py`.

### 6.1 Common settings (every cell)
- Ideal model: circular coplanar.
  - Heliocentric cells: R-S 2007 Table 2 constants for Venus-Mars; Venus 0.61520 yr for Earth-Venus;
    Mars 1.875 yr for Earth-Mars.
  - Jovian cells: R-S 2009 Table 2 constants.
- Return catalogue per body:
  - resonant n:m in {1:1, 2:1, 1:2, 3:2, 2:3}.
  - half-rev (1, 0, peri/apo) and (3, 1, peri/apo), both mirror signs.
  - generic same-body Lambert legs (revs as stated per cell, both branches).
- Structures: [A-block, A->B, B-block, B->A], one visit per cycle, period k synodic periods.
- Seeds: n_phase 36, n_split 12, n_refine 40, spread across basins (min separation 0.03 T).
  These are the settings at which every recall control in secs. 4-5 passed.
- Zero: max |residual| < 1e-8 km/s (exact zero; candidates are held to exact closure, not to H&M's
  0.005-EMOS tolerance).
- Free directions: minimax of demanded/available turn, with the even-spread tie-break.
- Gate at every massive flyby: #888/#937 demanded-turn gate (`verify/turn_gate.py`, tri-state), at the
  registry floors. Earth 200 km, Venus 300 km, Mars 200 km, Ganymede and Europa 100 km, Callisto 200 km.
- Near-180 rule: a demanded turn >= 175 deg is a rejection (PROVISIONAL; owner ruling #937/#906 gives
  no number). Every zero stores its largest demanded turn.
- PASS (a "gate-passing zero"):
  - exact zero.
  - gate status "pass" at every massive flyby ("indeterminate" is reported separately and is NOT a
    pass).
  - no near-180 demand.
  - independent re-propagation miss < 1 km at every encounter (far inside every SOI).
- Dedupe (analysis script): two zeros are the same physical cycler when the cyclic sequence of massive
  flybys (body, V_inf to 1 m/s, turn to 0.1 deg) agrees up to rotation. Mirror (time-reversed) twins
  are merged and flagged. Split labels of one conic at a massless target collapse automatically.
- A gate-passing physical cycler is only "candidate, pending gauntlet". The gauntlet:
  - collision with Hollister-Menning (15 rows plus 1H-3H), Rall 1969 / Rall-Hollister 1971 (M4-1, M5-1,
    M5-2), R-S 2007/2009 (all rows incl. VenMar#45, EurGan, GanEur, GanCal), Campagnola 2019 GCGC,
    Jones 2017 and the catalogue (`our_status`).
  - `literature_check.py`.
  - an independent cross-check (re-solve by a different code path).
  - SOI self-consistency at every encounter.
  - No catalogue writes.

### 6.2 Cells, in launch order

| # | Cell | Bodies (massive) | k | Returns per block (A,B) | Transfer revs | Generic revs | Structures | In-run recall control |
|---|---|---|---|---|---|---|---|---|
| 1 | vm (R1(c), one working body) | Venus; Mars massless | 1-3 | 3, 0 | 0-1 | 1-2 | 14,044 | VenMar#45 at k = 2 must appear as a gate-passing cycler, else the run is void |
| 2 | vm2 (R1(c), two working bodies) | Venus, Mars | 1-3 | 3, 0 (Mars bends at pass-through only) | 0-1 | 1-2 | to be counted before launch | the vm2 zeros with zero Mars turn must coincide with vm cyclers |
| 3 | gc (X1) | Ganymede, Callisto | 1-3 | 1, 1 | 0-2 | 1-2 | 4,043 | gc1 recall of GanCal#1/#5 (passed, sec. 5) re-run as a slice; GCGC (Campagnola 2019, V_inf 3.5/4.5 km/s) as a geometry check |
| 4 | ge (X1) | Ganymede, Europa | 1-3 | 1, 1 | 0-2 | 1-2 | 5,533 | ge1 recall of GanEur#5/#43 (passed, sec. 5) |
| 5 | ev (R1(b)) | Earth, Venus | 1-3 | 2, 2; EARTH returns full-revolution only (`--resonant-only E`) | 0 | 1 | to be counted | Hollister 1H and 2H topologies must appear |
| 6 | em (R1(a)) | Earth, Mars | 1-3 | 2, 2 | 0 | 1 | to be counted | every candidate compared with Rall's families |

Notes on the cells:
- Earth-symmetric families (orbits 3-8 type) are excluded from ev by `--resonant-only E`, per the owner
  ruling.
- Empty-cell stamps go to `data/empty_regions.jsonl` with these exact settings as the method scope.
- Workers: 2 shards per cell (`--shard 0/2`, `1/2`), launched by the lead.
- Timing pilot for vm: 2.8 s per structure (20-structure k = 3 sample). Serial about 11 h; about 5.5 h
  on 2 workers, longer under load.

### 6.3 Cell 1 result: vm (Venus-Mars, Mars massless), 2026-10-06

Run (launched by the lead, 2 shards): 14,044 structures, 2,046 exact zeros, 344 physical cyclers, 2
gate-passing. Gauntlet: `data/942_cell_vm_gauntlet.json`.

- In-run control: VenMar#45 recovered (LITERAL). The run is valid.
- Bug fix during adjudication: `kepler_step` diverged near aphelion at e about 0.93. All zeros were
  reassessed (`scripts/reassess_942_zeros.py`); no verdict changed.
- Candidate vm-1, a class member, NOT novel. k = 2; one Venus 1-rev generic return plus a 0-rev
  V->M->V conic.
  - V_inf: Venus 41.81 km/s, Mars 20.78 km/s.
  - Venus turn 3.05 deg of 3.26 available (742 km required altitude).
  - Perihelion 0.069 AU (about 15 solar radii), aphelion 1.95 AU.
  - DOP853 re-fly miss 0.009 km.
  - Recorded as a MEMBER OF THE PUBLISHED R-S 2007 CLASS, not a find (lead ruling 2026-10-06). R-S ran
    this exact class in the same ideal model and archived the complete trajectories (AAS 07-118 p.9), so
    it is not an author-excluded class under #875.
  - Not a credible design: it is a sun-grazer.
  - Full state: k = 2, period 667.85 d, key k2|LV>V/1l|LV>M/0s|LM>V/0s. Lambert-leg start dates (d,
    model epoch, all bodies at angle 0 at t = 0): 316.278, 729.731, 915.491.
    - Legs: V-V 1-rev low, 413.5 d, a 1.0074 AU, e 0.9318. V-M-V 0-rev, a 0.8737 AU, e 0.9405.
  - Registry stamp: heliocentric-venus-mars-one-working-body-rs2007-ideal-k1-3-942.

### 6.4 Cell 2 result: vm2 (Venus-Mars, both massive), 2026-10-06

Run: 14,044 structures, 6,744 exact zeros, 2,020 physical cyclers, 3 gate-passing (reassessed with the
fixed `kepler_step`; no verdict changed). Gauntlet: `data/942_cell_vm2_gauntlet.json`.
- VenMar#45 (LITERAL) and vm-1 (sec. 6.3): both have a Mars turn of 0.0, so they are the one-body
  cyclers again.
- **vm2-1: candidate, pending owner adjudication. NOT novel.** Mars bends, so it is two-working-body.

**vm2-1 full state** (R-S 2007 ideal model: circular coplanar, Sun mu 1.3271244e11, Venus period
19,414,153 s, Mars 59,354,429 s, all bodies at angle 0 at t = 0).
- k = 3, period 1001.770 d. Key k3|LV>V/2l|LV>M/0s|LM>V/0s. Lambert-leg starts (d): 52.916705,
  614.929975, 834.808349.

| Leg | Type | Flight (d) | a (AU) | e | Perihelion / aphelion (AU) |
|---|---|---|---|---|---|
| V->V | generic return, 2 rev, low | 562.01 | 0.74537 | 0.17372 | 0.616 / 0.875 |
| V->M | 0-rev transfer | 219.88 | 1.12250 | 0.35896 | 0.720 / 1.525 |
| M->V | 0-rev transfer | 219.88 | 1.12250 | 0.35896 | 0.720 / 1.525 |

On the generic return Venus makes 2.5 revolutions and the spacecraft 2. The two transfers form a
symmetric 440-d round trip.

| Flyby | V_inf (km/s) | Turn (deg) | Available (deg) | Ratio | Required altitude |
|---|---|---|---|---|---|
| Venus (x2) | 6.086 | 69.55 | 70.90 at 300 km | 0.981 | 556 km (1.092 radii) |
| Mars | 4.849 | 16.42 | 39.27 at 200 km | 0.418 | 7,533 km (3.2 radii) |

Checks:
- Exact zero (residual < 1e-8 km/s).
- Independent DOP853 re-fly: miss 0.013 km, V_inf vector error 4e-9 km/s, gate pass on the integrated
  vectors.
- lamberthub (Izzo/Gooding) agreement on all three legs: 1e-7 m/s.
- SOI fraction 2e-8.

**Both acceptance criteria, side by side:**
- PASSES the project's registry floor: Venus 300 km.
- FAILS Rall's and H&M's own acceptance, "flybys ... beyond 1.1 planet radii" (Menning p.6; H&M p.1193),
  i.e. 605 km at Venus. vm2-1 needs 556 km.
- The neighbouring zero of the same structure (V_inf 6.11 / 5.87) fails the gate at ratio 1.07.

**Prior attempt (Rall 1969 Sc.D. thesis, sec. 4.4, pp.84-86, read on the page images):**
- "The approach to obtaining such a periodic orbit was to have the trajectory go from Venus to Mars to
  Venus in a low energy fashion (taking around 400 or 450 days for the round trip and making about one
  revolution of the Sun) while Venus makes about two revolutions of the Sun. Then find appropriate
  direct return trajectories in the vicinity of Venus until the next opportunity for a transfer to Mars
  presents itself."
- "Attempts were made with repeating periods of up to four Mars-Venus synodic periods. In all cases,
  the attempted method either did not converge, or the trajectory intersected the surface of Venus."
- "More promising direct return orbits appeared to be those which travel around the Sun a different
  number of times than does Venus ... In each case, the attempt either failed to converge; or the
  resulting trajectory intersected the surface of Venus."
- "this investigation did not prove that no periodic orbits of the type considered connect Mars and
  Venus; this investigation simply failed to find any."

vm2-1 is in that attempted class: a 440-d V-M-V round trip, plus a Venus return that circles the Sun 2
times while Venus does 2.5, at 3 synodic periods. It closes ballistically only between the two flyby
criteria above.

Other gates:
- R-S 2007: their ideal model has Mars massless (p.2), and the massive target is named as future work
  (p.18). VenMar#45 is the only printed Venus-Mars row. No collision.
- Pisarevsky 2008: the method is general, but its numbers are Earth-Mars only. Its diagrams restrict
  loitering arcs to multiples of pi and name "any number of generic returns" as future work (digest sec.
  4). vm2-1's Venus return is generic (non-k pi), so it is outside their covered diagrams.
- Turner AAS 07-175: unheld.
- Catalogue: only Jones VEM rows on this pair.
- Offline literature_check: "published" via the R-S Venus-Mars body-pair anchor (a flag, not a
  verdict).

Approved next steps: robustness pre-registration (sec. 6.5); registry stamp for the rest of the cell.

### 6.5 vm2-1 robustness: PRE-REGISTRATION (lead approval 2026-10-06), committed before any run

(a) Venus altitude floor, 300 to 700 km in 50-km steps. Re-gate EVERY vm2 zero (all 6,744, exact
    zeros unchanged) at each floor; report the gate-passing physical cyclers per floor, and the floor at
    which vm2-1 stops passing. Its required altitude is 556.3 km, so the prediction is "passes up to
    556 km". The question is whether any other vm2 zero joins or replaces it at other floors.

(b) Mars mass:
- In the patched-conic ideal model the Mars mass does not enter the zero; it enters only the Mars
  available bend. So: (b1) find the smallest Mars GM fraction at which the 16.42-deg Mars turn is
  available at 200 km.
- (b2) The one-body (massless-Mars) member of the SAME structure. From the vm run it is the zero with
  V_inf 6.110 / 5.866 km/s, Mars turn 0, Venus ratio 1.075 (gate FAIL). It is in R-S's archived class.
- (b3) Is vm2-1 continuously connected to it? Use a homotopy that keeps the Mars magnitude match and
  scales the allowed Mars turn, adding the planets' period ratio as the extra free parameter. Report
  whether the path connects; no claim beyond that.

(c) Neighbours:
- A targeted enumeration of the same template, both cells (vm2 and vm):
  - k = 2..5.
  - Venus generic returns of 1-4 revs, both branches.
  - V<->M transfers of 0-1 rev.
  - up to one extra Venus full-rev or half-rev return.
- Same seeds and gate as sec. 6.1. Report every gate-passing physical cycler with a nonzero Mars turn,
  and whether vm2-1 is isolated or part of a k-series.

(d) Real ephemeris: a patched-conic ephemeris continuation in the R-S 2007 style.
- An OPEN chain of n consecutive vm2-1 cycles (n = 7, i.e. about 7 V-M synodic periods = 2337 d, the
  near-commensurability noted by Gillespie & Ross), with every encounter date free and no periodic wrap.
- Start in the R-S circular model; homotopy to fixed mean elements (Standish J2000 Venus/Mars, real
  periods, e and i) in 10 steps.
- Report per step: the convergence residual, the Venus and Mars turn ratios at every encounter of the
  chain, and the epoch dependence (5 start epochs spread over one 32-yr cycle).
- Pass at this rung: every encounter ballistic at the registry floors, at >= 1 epoch.
- Then the #866 V3 lane only if that passes (a separate launch request).

All four are diagnosis. None changes vm2-1's label.

### 6.6 vm2-1 robustness results so far (a, b1, b2)

- (a) Floor sweep, all 6,744 vm2 zeros re-gated (`data/942_vm2_1_robust_a_floor_sweep.json`): vm2-1
  passes at Venus floors 300-550 km and fails from 600 km. No other cycler joins at any floor. VenMar#45
  and vm-1 pass at all floors to 700 km.
- (b1) The 16.42-deg Mars turn needs at least 0.329 of Mars's GM at the 200-km floor.
- (b2) The one-body (massless-Mars) member of the same structure is the vm zero with V_inf
  6.110 / 5.866 km/s and Venus ratio 1.075: it FAILS the gate. So the Mars turn is what makes vm2-1
  feasible. It lowers the Venus demand from 1.075 to 0.981 of capacity.
- (b3), (c) and (d) are pending; (c) and (d) need launches.

### 6.7 Cell 5 result: ev (Earth-Venus, Earth returns full-revolution only), 2026-10-06

Run: 2,770 structures, 13,915 exact zeros, 3,534 physical cyclers, 31 gate-passing (reassessed; no
verdict changed). In 23 of them the transfer skeletons are distinct (same transfer V_inf; they differ
only in the block returns). Gauntlet: `data/942_cell_ev_gauntlet.json`.
- In-run control: Hollister 1H and 2H appear, so the run is valid. 3H is excluded by
  `--resonant-only E`.
- Hollister family (PUBLISHED, Hollister 1969 orbits I/II), 4 skeletons:
  - (k=2; E/V V_inf) 2.99/3.19 (1H; also a variant with a half-rev pair at Venus).
  - 6.14/5.28 (1H).
  - 4.05/7.07 and 5.60/6.02 (2H).
- Earth-massless (Earth turn 0; Venus hosts every return; the R-S one-body architecture at Venus),
  5 skeletons:
  - k=2: 9.07/13.17, 13.85/12.77.
  - k=3: 11.06/5.97, 11.12/6.13, 17.03/13.72.
  - AAS 07-118 ran no Earth-Venus set (#960 gate: R1(b) OPEN).
- Two-working-body, outside Hollister's 3.2-yr itineraries, 14 skeletons:
  - k=2: 3.51/4.41, 4.40/8.82, 4.89/10.36, 6.04/4.16.
  - k=3: 3.34/4.04, 4.30/5.82, 5.06/5.93, 5.58/3.77, 5.92/3.58, 5.99/3.70, 6.24/4.34, 6.26/6.13,
    8.01/10.93, 9.43/5.24.

All are "candidate, pending owner adjudication". NOT novel. Collision context:
- Menning 1968 p.41-42 estimates "a minimum of 1024" Earth-Venus orbits of the FR/SY type, and names
  half-rev and order variations as possible but not computed.
- Pisarevsky 2008: the method is general; its numbers are Earth-Mars only.
- No Earth-Venus catalogue row other than the 15 H&M rows. The DOP853 re-fly miss is < 0.1 km for all 31.

### 6.8 Amendment to (d), before its full run

The first n = 7 attempt failed to converge at lambda = 0.1 at epoch 0 (residual 0.029 km/s), and the
unconverged least-squares calls ran past the 10-minute limit. The single-cycle validation (n = 1) went to
lambda = 1 in 10 steps.

Numerical change only, with the criteria unchanged:
- adaptive lambda step: 0.1, halved on failure, down to 1/640.
- least-squares evaluation cap of 50 x unknowns.

The last converged lambda is reported for each epoch.

### 6.9 vm2-1 robustness (d) result: it does NOT survive toward the real ephemeris

Data: `data/942_vm2_1_robust_d_realeph_chain.json`. A 7-cycle open chain was continued from the R-S circular
model toward Standish J2000 mean elements, at 5 epochs from 2030 to 2052.

| Epoch (JD) | Last converged lambda | Venus ratio there | Last lambda with a gate pass |
|---|---|---|---|
| 2462782.4 | 0.067 (fold) | 1.006 | 0.050 |
| 2465105.9 | 0.061 (fold) | 1.000 (indeterminate) | 0.056 |
| 2467463.2 | 0.067 (fold) | 1.003 | 0.050 |
| 2469784.5 | 0.077 (fold) | 1.013 | 0.0 (fails from the first step) |
| 2472130.0 | 0.125 (run cut at the 10-min limit) | 1.041 | 0.0 |

Findings:
- At every epoch the Venus turn exceeds capacity within the first 5-8 % of the way from circular to the
  real eccentricities and inclinations.
- At four of the five epochs the chain solution itself ends at a fold before lambda = 0.08.
- The single-cycle version (n = 1, epoch 0) continues to lambda = 1 but fails the gate from lambda = 0.2
  (Venus ratio 1.23, Mars 1.47 at lambda = 1).

Verdict at this rung: vm2-1 is an ideal-model-only object. Its 2 % Venus margin does not survive real
orbit eccentricity (Mars e = 0.093) and inclination. The #866 V3 lane is not requested.

### 6.10 vm2-1 final record (lead ruling 2026-10-06)

vm2-1 is an "ideal-model curiosity / negative-adjacent, not a find". It does not persist toward the real
ephemeris: the Venus ratio is over 1 by lambda of about 0.06 (sec. 6.9). Step (b3) was dropped. Step (c),
the neighbours, stays: a neighbour with more margin is the only way this family could survive real
eccentricity. The lead launches it after em.

### 6.11 Cell 3 result: gc (Ganymede-Callisto, both massive), 2026-10-06

Run: 4,043 structures, 2,300 exact zeros, 775 physical cyclers, 3 gate-passing (reassessed; no verdict
changed; 0 errors). Gauntlet: `data/943_cell_gc_gauntlet.json`.
- In-run control: GanCal#5 is LITERAL. GanCal#1 needs two Ganymede returns, so it lies outside this cell's
  `--max-returns 1,1` scope; it was recovered blind by the targeted run of sec. 5.
- GCGC (Campagnola 2019): not present. It has an alternating G1 C2 G3 C4 pattern with V_inf 3.5/4.5; no
  gate-passer matches within 0.5 km/s at both moons. The Lam 2015 13F7 pump-down is non-resonant, not a
  cycler, so there is no collision. Liang 2024 CGCEC is a three-moon tour.
- The #576 symmetric G-C-G closures (V_inf pairs 3.97/2.47, 7.59/3.71, 2.80/3.47, 3.28/4.70, 6.54/5.51,
  1.93/1.49, 7.66/3.67) have no returns and were judged by the capacity gate. None matches a gc
  gate-passer.

| Candidate | gc-1 | gc-2 |
|---|---|---|
| Status | candidate, pending owner adjudication, NOT novel | candidate, pending owner adjudication, NOT novel |
| Key | k3\|LGanymede>Ganymede/1l\|LGanymede>Callisto/0s\|RCallisto/1:1\|LCallisto>Ganymede/0s | k3\|LGanymede>Ganymede/1l\|LGanymede>Callisto/0s\|LCallisto>Callisto/1h\|LCallisto>Ganymede/0s |
| Period | 37.57 d (3 G-C synodic) | 37.57 d |
| V_inf G / C (km/s) | 2.397 / 1.807 | 3.617 / 3.039 |
| Ganymede turn | 29.35 deg of 45.44 (ratio 0.646, 2,437 km) | 19.72 deg of 25.02 (ratio 0.788, 1,024 km) |
| Callisto turns | 2 x 40.15 deg of 54.45 (ratio 0.737, 1,801 km) | 2 x 6.87 deg of 26.56 (ratio 0.259, 9,800 km) |
| Distance from Jupiter (km) | 888,745 to 1,955,850 | 791,455 to 2,337,392 |
| Lambert starts (d) | 0.766855, 11.345133, 31.824260 | 0.838679, 11.684556, 14.289304, 35.803635 |

Leg detail:
- gc-1:
  - G-G 1-rev low, 10.578 d (a 1,142,885 km, e 0.2224).
  - G->C 3.790 d and C->G 6.512 d on one conic (a 1,481,355 km, e 0.3203).
  - Callisto 1:1 full-rev return.
- gc-2:
  - G-G 1-rev low, 10.846 d (e 0.3300).
  - G->C and C->G of 2.605 d each (a 1,572,225 km, e 0.4120).
  - C-C 1-rev high, 21.514 d (e 0.3856).

Checks, both:
- DOP853 re-fly miss < 1e-5 km.
- lamberthub agreement < 3e-9 m/s on every leg.
- Gate pass on the integrated vectors.

Both moons bend in both, so these are two-working-body cyclers. R-S 2007/2009 name the massive target
as future work; GCGC is the only published two-working-body G-C cycler, and it is a different class.

### 6.12 ev shortlist (lead ruling 2026-10-06)

Ranking: worst gate ratio, then the lowest maximum V_inf; anything inside Menning's p.41-42
"variations" is demoted; the best Earth-massless skeleton is included. Full list:
`data/942_cell_ev_shortlist_ranking.json`. One skeleton is a Menning variation (k=2, 3.51/4.41: FR at
Earth, half-rev plus symmetric at Venus) and is demoted to last.

Top 3, for the full ladder:
- ev-A: k=2, E/V 4.89/10.36, worst ratio 0.574, two-working-body. Key k2|LE>V/0s|RV/1:1|LV>V/1h|LV>E/0s.
- ev-B: k=3, 8.01/10.93, worst ratio 0.673, two-working-body. Key k3|RE/1:1|LE>V/0s|LV>V/1l|RV/3:2|LV>E/0s.
- ev-C: k=2, 9.07/13.17, worst ratio 0.746, Earth-massless (R-S architecture at Venus). Key
  k2|LE>V/0s|LV>V/1h|LV>E/0s.

The other 16 non-Hollister skeletons are recorded as "candidate, not laddered".

### 6.13 Ladder for the shortlisted candidates: PRE-REGISTRATION (2026-10-06), committed before running

Candidates: gc-1, gc-2 (sec. 6.11); ev-A, ev-B, ev-C (sec. 6.12).

Rung (a), floor sweep: re-gate all of the cell's zeros at the target body's floor, raised from the
registry value to 3x it in 6 steps; report where each candidate stops passing.

Rung (d), real ephemeris: `scripts/run_942_realeph_chain.py`.
- Generic tool: any cell and key, full-rev and half-rev legs supported.
- The first date is fixed; every other date and the chain period are free; H&M closure B at the end.
- A homotopy from the ideal circular model, rotated into the bodies' mean orbital plane and
  phase-matched, to the real ephemeris:
  - Standish J2000 mean elements for Earth-Venus.
  - NAIF jup365 (SPICE) for the Galilean moons.
- Adaptive lambda step; every interior flyby gated at the registry floors.
- Chains:
  - gc-1, gc-2: n = 10 cycles (376 d).
  - ev-A, ev-C (k = 2): n = 5 cycles (16 yr).
  - ev-B (k = 3): n = 4 cycles (19.2 yr).
- 5 epochs, every 6.4 yr from 2030-01-01.
- RUNG PASS: lambda = 1 reached and every interior flyby passes, at >= 1 epoch. All epochs are reported.
- Validation:
  - the tool reproduces the vm2-1 single-cycle result (worst ratio 1.46 at lambda = 1, against 1.465 from
    the earlier script).
  - gc-1, single cycle, epoch 0: passes at lambda = 1 (worst 0.742).

Rung (c), neighbours: only if (d) passes. A separate launch request.

### 6.14 Ladder rung (d) results (2026-10-06), real-ephemeris chains

Data: `data/942_943_ladder_realeph_chains.json`. Independent re-fly: `scripts/check_942_realeph_chain.py`
(DOP853, rtol 1e-13, from the real-ephemeris states).

| Candidate | Legs | Rung (d) by the tool (5 epochs) | Independent DOP853 re-fly at lambda = 1 | Verdict at this rung |
|---|---|---|---|---|
| gc-2 | Lambert only | PASS at 5/5; worst ratio 0.813-0.823 (Ganymede) | miss <= 2.4e-5 km, V_inf error <= 1.4e-10 km/s | PASS: ballistic over 10 cycles (376 d), patched conic, NAIF jup365, all epochs |
| ev-C | Lambert only | PASS at 5/5; worst 0.869-0.914 (Venus); Earth ratio 0.06-0.11 | miss <= 0.12 km, V_inf error <= 3e-8 km/s | PASS: ballistic over 5 cycles (16 yr), patched conic, Standish mean elements, all epochs |
| gc-1 | includes a Callisto 1:1 full-rev | "PASS" at 5/5 (0.75) | miss 4,300-6,500 km, V_inf error 19-27 m/s | NOT VALID: the full-rev legs do not close |
| ev-A | includes a Venus 1:1 full-rev | "PASS" at 5/5 (0.58) | miss about 52,000 km, V_inf error 17 m/s | NOT VALID: same defect |
| ev-B | includes 1:1 and 3:2 full-revs | 2/5 by the tool | not re-flown | NOT VALID: same defect |

The defect ("it closed!" was the danger signal):
- On the real ephemeris the chain tool timed a full-revolution return at the IDEAL model's body period,
  with |v_sc| = |V_P|. The real body does not return to the same point in that time, so those legs
  do not close.
- The tool's own gate saw only V_inf magnitudes and directions, not the leg closure. The independent
  re-fly caught it.
- Fix needed before gc-1, ev-A or ev-B can be judged: on the real ephemeris a full-rev return must
  become a free-date 1-rev same-body Lambert leg (its circle of directions collapses to discrete
  solutions), or be charged its correction Delta-V. Until then those three are "not judged at rung (d)".

Reproduction checks of the tool:
- the vm2-1 single cycle gives 1.46 (the earlier script gave 1.465).
- the vm2-1 chains (Lambert only) are unaffected by the defect.

Status: gc-2 and ev-C are the strongest results so far. They remain "candidate, pending owner
adjudication", NOT novel. Gauntlet items still open for them:
- a web literature search (only the offline corpus was checked);
- an n-body (V3-lane) check;
- for ev-C, a DE440 run (Standish fixed mean elements are not DE440). DONE 2026-10-06: on DE440
  (astropy backend) ev-C passes at all 5 epochs, worst 0.869-0.914. The dates move by about 0.23 d
  against the Standish run. The DOP853 re-fly miss is <= 0.097 km.

### 6.15 Chain-tool positive controls (PRE-REGISTERED 2026-10-06, before running)

What happened: the chain tool was NOT run on a published member before the sec. 6.14 verdicts (lead's
condition). Only one of our own candidates was used to validate it. That is why the full-rev defect
reached the candidates. Controls, judged by the SAME criterion as the candidates (rung pass by the
tool AND the independent DOP853 re-fly: miss < 1 km, V_inf error < 1e-6 km/s, at every leg):

- C1, VenMar#45 (R-S 2007: "easily converges to ballistic" in a patched-conic ephemeris model, Fig. 9a).
  Run in cell vm2 (Mars massive): R-S restrict masslessness to the ideal model (p.2). In 3-D a massless
  target is over-determined by 2 equations per passage, so the tool now refuses massless cells.
  - 7-cycle chain; Standish elements and DE440; 5 epochs.
  - EXPECTED: pass.
- C2, GanCal#5 (Lambert-only), in cell gc (Callisto massive). 10 cycles, jup365, 5 epochs.
  - Sanity only: R-S printed no ephemeris result for #5.
- C3, Hollister 1H (full-rev returns), cell ev.
  - 5 cycles (16 yr), Standish, 5 epochs.
  - EXPECTED to FAIL the re-fly under the current tool (the defect). This is a NEGATIVE control: it
    shows the re-fly catches the defect.
- After the full-rev fix: C4 = GanCal#1 (R-S 2007: "ballistic over 10 cycles", Fig. 9b; it contains an
  f(2:1) leg) and C3 again. Both must pass by the re-fly before gc-1, ev-A or ev-B are judged.

### 6.16 Chain-tool control results (2026-10-06), data 

| Control | Expected | Result |
|---|---|---|
| C1 VenMar#45, cell vm2, 7 cycles, Standish | pass | PASS at 5/5 epochs; worst ratio 0.31-0.49 (Venus), Mars 0.08-0.25. Re-fly miss <= 0.03 km, V_inf error <= 9e-9 km/s |
| C1 VenMar#45, DE440 | pass | PASS at 5/5; same ratios to 3 decimals. Re-fly miss <= 0.025 km |
| C2 GanCal#5, cell gc, 10 cycles, jup365 | sanity | FAILS at 5/5. The Ganymede ratio rises from 0.94 (ideal) to 1.10-1.23 by lambda = 1 (gate pass up to lambda 0.3-0.4). R-S printed no ephemeris result for #5, so there is no contradiction; it bounds what a 6 % ideal margin survives |
| C3 Hollister 1H, 5 cycles, Standish | negative control | the chain does not converge past lambda = 0.29-0.39 at any epoch |

Reading:
- The Lambert-only path of the tool reproduces a published real-ephemeris result (VenMar#45). That
  supports the gc-2 and ev-C passes (Lambert-only, re-flown).
- The 1H negative control fails earlier than predicted: the run never reached lambda = 1, so the re-fly
  was not the stage that exposed it. The full-rev defect shows up as a fold.

Full-rev fix (next, approved): on the real ephemeris each full-rev leg becomes a shooting leg.
- Unknowns: the departure direction (2 angles) and the flight time.
- Equations: the arrival position equals the body's position (3 equations).
- Its |V_inf| is fixed by the junction.
- Then C4 (GanCal#1) and C3 again must pass by the re-fly before gc-1, ev-A or ev-B are judged.

## 7. Literal-collision checks (to be completed per candidate)

R1(a) gate addition (lead ruling, 2026-10-05): Rall 1969 and Rall & Hollister 1971 (JSR 8(10):1017, doi
10.2514/3.59763) computed Earth-Mars periodic orbits with Mars as a ballistic working body (Mars turns
2.3-13.6 deg, never direct returns at Mars). Any R1(a) candidate is checked against the families M4-1,
M5-1 and M5-2 (circular and eccentric-inclined) before novelty language is used.

Method cross-check: the full-rev circle reproduces Russell & Ocampo 2005 Eqs. (13) and (17) exactly (test,
Fig. 10 setup). Our half-rev radial and transverse speeds are Eqs. (18)-(19) with r1 = r2.

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
