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

### 6.17 Full-rev fix, and the runs that depend on it: PRE-REGISTERED 2026-10-06, before running

Fix (`scripts/run_942_realeph_chain.py`, `chain_eval`): every fixed leg (full-rev, half-rev) is a shooting
leg.
- Unknowns: the departure direction (theta, phi about the body velocity) and the flight time tau.
- |V_inf| is fixed by the junction.
- Equations: the arrival position equals the body's position (3 residuals; 1 km counts as 1e-3 km/s).
- The free-direction minimax is gone on the real ephemeris: the directions are solved for.
- The re-fly checker now re-flies EVERY segment (Lambert and shot) with DOP853.
- Unit check: gc-1, single cycle, epoch 0 converges at lambda = 1 with a re-fly miss of 3e-6 km.

Runs, all with the same settings as sec. 6.13 (5 epochs; rung pass = lambda 1 reached, every interior
flyby passes, and the re-fly miss is < 1 km at every segment):
1. C4 GanCal#1 (published "ballistic over 10 cycles", R-S 2007 Fig. 9b): cell gc, 10 cycles, jup365.
   EXPECTED pass. Dates from the gc zero: 0.991426, 26.049946, 35.954310 d.
2. C3 Hollister 1H (published inclined-elliptic periodic): cell ev, 5 cycles, Standish. EXPECTED pass.
3. Only if C4 and C3 pass: gc-1 (10 cycles), ev-A (5), ev-B (4).
4. Regression: C1 VenMar#45, epoch 0 only (Lambert-only path unchanged).

### 6.18 Amendment (approach change, criteria unchanged) and what the controls now show (2026-10-06)

Approach change (advisor review, recorded before the remaining runs):
- Standish runs use a Keplerian RAMP (`RampedKepler`): e and i ramped from 0 to their actual values
  (Hollister 1969 p.368), with a, mu and the epoch longitude interpolated. A full-rev leg timed at the
  current Keplerian period is then exact at every lambda, and the minimax stays valid.
- Blend runs (jup365, DE440) are time-and-minimax continuation devices only. At lambda = 1 every fixed
  leg is shot, with restarts over its direction, and the closing solution with the smallest worst ratio
  is kept.
- DE440 is not Keplerian: report the closure Delta-V of full-rev candidates there, with no pass/fail.
- Bug fixed: shot fixed-leg parameters were assigned in block order rather than leg order whenever a key
  begins with a fixed leg (1H does). This produced the C3 "residual inf".

Findings that change how rung (d) must be read:
1. ENDPOINT CONTROL PASSES. The published H&M Earth-full-rev orbits 1, 2, 11, 12, 13 and 15, solved on
   the real-period Standish model:
   - all close as 16-yr chains.
   - all pass the gate at the registry floors (worst ratio 0.78-0.89).
   - all re-fly with DOP853 to a miss <= 0.001 km, including the full-rev legs.
   So the gate, the full-rev timing and the re-fly are validated on published members at lambda = 1.
2. THE HOMOTOPY CAN FOLD WHERE A REAL-EPHEMERIS SOLUTION EXISTS. Ramp-mode continuation of our circular
   1H member folds at lambda of about 0.29-0.39, in both the free-period and the fixed-period
   formulation. Yet 1H-family orbits exist at lambda = 1 (finding 1). Our circular 1H member need not be
   the parent of H&M's orbit 1, but either way a fold is NOT a negative. Gate verdicts at intermediate
   lambda describe artificial blended models and carry no physical meaning; only lambda = 1 counts.
3. Consequences:
   - vm2-1: its 7-cycle "negative" (sec. 6.9) is PATH-LIMITED (folds before lambda 0.13). The
     single-cycle lambda = 1 failure stands (Venus 1.23, Mars 1.47). The record must say "single-cycle
     real-ephemeris failure; 7-cycle chain not reached".
   - GanCal#5 (C2) reached lambda = 1 and failed there: a real failure.
   - gc-2 and ev-C reached lambda = 1, passed and re-flew: unaffected.
   - gc-1, ev-A, ev-B: not judged. A lambda = 1 seeding route that does not depend on the homotopy is
     needed (e.g. a direct solve from the ideal dates at several epochs, or pseudo-arclength through
     the fold).

### 6.19 Amendment (solver only) and the C4 GanCal#1 control status (2026-10-06)

Solver changes, recorded before the next recorded run (criteria unchanged: closure max residual
< 1e-6, gate at the registry floors, DOP853 re-fly):
- `x_scale="jac"` restored for the lm shoot. Commit 1b8cfbfa forced 1.0 (SciPy >= 1.16 defaults lm to
  "jac"), and that stalled the solver; the 1-cycle sweep made with it is void (papercut).
- `--shoot-jac sparse`: a column-grouped forward-difference Jacobian. It reproduces the dense result
  exactly on the 1-cycle GanCal#1 sweep (same closures, same ratios) and is about 3x faster at 10 cycles.
- A Gauss-Newton polish of near-closures (max residual between 1e-6 and 1e-2), with steps accepted only
  if they reduce the max residual.
- `--shoot-nfev-per-var 60` for the 10-cycle runs.

C4 GanCal#1 on jup365 (diagnostics, not the pre-registered rung):
- 1 cycle, 5 epochs, 20 restarts: closes at 3 of the 4 epochs that reach lambda = 1. The best closure
  per epoch is gate "indeterminate" at worst ratio 0.968, 0.978 and 4.91 (fail). The first two are inside
  three tidal turn scales of the margin; GanCal#1's ideal-model ratio is already 0.961-0.966, so it is
  marginal in every model. One epoch folds in the blend at lambda 0.45.
- 10 cycles, epoch 1, restart 0: the max residual falls from 1e-3 to 5.5e-6 and then stops (LM ends on
  its own tolerance; the polish cannot improve it, even with central differences). The
  LGanymede>Ganymede/1l legs span 179.74-179.78 deg between their endpoints, but that is NOT the
  cause (correction, same day): the 1-cycle epoch-1 closure has the same leg at 179.75 deg and
  converges to 2.9e-9. The stall grows with chain length: ten weakly determined crank angles, and
  a forward-difference Jacobian whose smallest singular value changes with the step size
  (1e-4 to 6e-3). Cause not yet isolated. So the pre-registered 10-cycle C4 cannot be decided with this
  formulation; it is NOT a failure of the member.

gc-1 (diagnostics, same settings): its Ganymede-Ganymede leg spans about 172 deg.
- 1 cycle, 5 epochs: closes at every epoch, with gate passes at worst ratio 0.726-0.742 (ideal 0.737).
- 10 cycles, epoch 0: restart 0 closes in about 1 s and passes the gate at worst ratio 0.758. The DOP853
  re-fly misses by at most 2.6e-6 km over 40 segments.
No gc-1 verdict until the control question is ruled on (lead's instruction).

### 6.20 gc-2 and the R-S GanCal family: Callisto-mass question, PRE-REGISTERED 2026-10-06 (before running)

Question (lead and corpus-file-opus): is gc-2 a two-working-body relative of GanCal#1/#5 (Callisto
massless), as tested by a Callisto-mass homotopy?

Model fact that settles the patched-conic part without computation: the date residuals (V_inf
magnitude match at each junction) do not contain any body's GM. GM enters only the gate's turn
capacity. So a homotopy in Callisto's GM leaves every solution fixed, including gc-2, GanCal#1 and
GanCal#5; only the gate verdicts change. A massless-Callisto member (GanCal) must also meet the vector
match at Callisto (zero turn). gc-2 turns 6.87 deg at each Callisto encounter, so it is not a GanCal
member at any Callisto mass, and in this model no Callisto-mass path joins them. Checking this
through continuous gravity (CR4BP or n-body) would need a different model: not done here, an owner
choice.

Structure (from `data/943_cell_gc_gauntlet.json`):
- GanCal#5 (LITERAL in-run) = LG>G/1l | LG>C/1h | LC>G/0s, Callisto turn 0.0.
- gc-2 = LG>G/1l | LG>C/0s | LC>C/1h | LC>G/0s.
- gc-2 is GanCal#5's skeleton with GanCal#5's 24.25-d G->C 1-rev leg replaced by a 2.60-d G->C leg plus
  a 21.51-d C->C 1-rev leg: one extra Callisto encounter inserted.

Computations (script `scripts/analyse_943_gc2_gancal_relation.py`, ideal circular gc model):
- H1a: gate verdict of gc-2 and of GanCal#5 against a Callisto GM scale s (log grid 1e-3 to 1; then
  bisection to 1e-4 relative for s*, the smallest s at which gc-2's Callisto encounters are still
  feasible). This is the only effect of the homotopy in this model.
- H1b (descriptive): along GanCal#5's G->C 1h leg (Kepler conic, 0.01-d steps, excluding the last 0.5 d),
  the minimum distance to Callisto and its time after departure. Compare it with Callisto's sphere of
  influence (`sphere_of_influence_km`) and with gc-2's insertion time, 2.60 d after departure. Also
  compare the conic elements (a, e) of the two legs.

Interpretation rules, fixed now:
- If the H1b minimum distance is below Callisto's SOI: "gc-2 has GanCal#5's skeleton with a Callisto
  encounter inserted where GanCal#5's leg passes close to Callisto: structurally related (INFERRED),
  not the same orbit".
- Otherwise: "the inserted encounter has no counterpart on GanCal#5's leg; the relation is the shared
  G-G 1l leg and the period only".
- GanCal#1 (LG>G/1l | RGanymede/2:1 | LG>C/0s | LC>G/0s; three Ganymede encounters and one Callisto
  encounter) has a different encounter count from gc-2 (two and two). No insertion test applies; it is
  compared by V_inf distance only (sec. 6.11 numbers).

### 6.21 Result of 6.20 (2026-10-06), data `data/943_gc2_gancal_relation.json`

Instrument check first: my first run showed no change of the Callisto ratio with s. Cause:
`dataclasses.replace` handed the scaled system the original's flyby-body cache (`_fb`). Fixed in the
script (fresh cache); the ratio now scales as expected (201 at s = 0.001, 0.259 at s = 1).

- H1a: as stated in 6.20, no solution moves. gc-2's Callisto encounters (6.87 deg each) are feasible
  down to s* = 0.2136 of Callisto's GM (ratio 0.695 at s = 0.316, 1.189 at s = 0.178). Below s*, gc-2
  does not exist as a ballistic cycler; it does not become a GanCal member. GanCal#5 (Callisto turn 0)
  passes at every s.
- H1b: GanCal#5's G->C 1-rev leg (24.25 d; a 1,683,097 km, e 0.4205) is nearest Callisto in its interior
  at 170,800 km, 2.90 d after leaving Ganymede (183,305 km at gc-2's insertion time of 2.60 d). That is
  4.5 times Callisto's SOI (37,681 km). Excluding the last 0.5 d, the minimum is at the window edge
  (148,120 km at 23.74 d), on the approach to the final arrival.
- Rule applied (fixed in 6.20): the minimum is above the SOI, so: **the inserted encounter has no
  counterpart on GanCal#5's leg; the relation between gc-2 and GanCal#5 is the shared G-G 1l leg and the
  period only.** Descriptive only: GanCal#5's leg does pass Callisto at about 4.5 SOI near the time gc-2
  meets it. gc-2's own legs are G->C and C->G on one conic (a 1,572,225 km, e 0.4120) and C->C 1-rev
  (a 1,686,915 km, e 0.3856, close in a to GanCal#5's leg but not in e).
- GanCal#1: different encounter count (three Ganymede and one Callisto against two and two); V_inf
  distance only (sec. 6.11: 0.44 at G, 0.22 at C).
- Not tested: a continuation in a continuous-gravity model (CR4BP or n-body), where encounters are not
  fixed in number. That is an owner choice.

### 6.22 Direct lambda = 1 route for the Standish full-rev rows: PRE-REGISTERED 2026-10-06 (before running)

Why: ramp continuation folds for Hollister 1H, yet 1H-family orbits exist at lambda = 1 (sec. 6.18). So
fold-based negatives mean nothing for ev-A or ev-B, and those need a route that does not depend on the
homotopy path. The route itself must first be controlled (advisor point).

Route (`--direct`, `scripts/run_942_realeph_chain.py`): one date solve at lambda = 1 (the real-period
Standish mean-element model, where a full-rev leg timed at the Keplerian period is exact), started from
the ideal-model dates phase-matched at each epoch. No continuation; nothing else changes. Same 5 epochs.

Re-fly for ramp-mode output (`scripts/check_942_realeph_chain.py`): at lambda = 1 the ramped model is
the mean-element system (checked to 1e-5 km over 16 yr), so every leg, full-revs included, is re-flown
with DOP853 against it (`gauntlet_942.cross_check`). The junction |V_inf| mismatch is printed too (it
includes the chain wrap, which closes in magnitude only).

Control D1, Hollister 1H (k2|RE/1:1|LE>V/0s|RV/1:1|RV/1:1|LV>E/0s, x 474.597193, 1100.618380 d), cell ev,
5 cycles:
- PASS if at least one epoch converges (date residual < 1e-6), passes the gate on the interior flybys,
  and re-flies with miss < 1 km and V_inf vector error < 1e-6 km/s. (A published 1H-family orbit exists
  at lambda = 1, sec. 6.18, so a sound route should find one.) The landing solution's V_inf is compared
  with H&M orbit 1 (descriptive).
- If D1 fails at every epoch, the route cannot judge ev-A or ev-B, and they stay "not judged".
Then, only if D1 passes, with the same settings and criteria: ev-A (5 cycles) and ev-B (4 cycles).

### 6.23 GM defect found by D1, and the results of 6.22 (2026-10-06), data `data/942_direct_route_and_gm_fix.json`

Defect (found by the D1 re-fly, fixed in c4ff9a41):
- The heliocentric ideal model's GM is a convention: a 1-AU circle with exactly a 365.25-d period. It is
  3.78e-5 above the real solar GM (132,717,453,060 against 132,712,440,018 km^3/s^2).
- Blend and RampedKepler kept that GM for the spacecraft at every lambda, so at lambda = 1 the bodies
  were real but the Lambert arcs were not. The re-fly checker took the GM from the solver, so it repeated
  the error instead of catching it.
- D1 showed it: the solver closed to 1e-12, but the junctions rebuilt against the mean-element system
  disagreed by 2.7e-4 km/s. Now the GM follows lambda, and the checker takes it from the real system.
- Affected: every ev-cell real-ephemeris result (ev-C, C3, D1). vm2 (1.4e-10) and gc (about 1e-8) are
  negligible.

Results with the fix (criteria as pre-registered):
- **ev-C re-run** (5 cycles):
  - Standish (now ramp mode): PASS at 5/5, worst 0.869-0.914; re-fly miss <= 0.093 km, V_inf error
    <= 2.5e-8 km/s.
  - DE440: PASS at 5/5, the same ratios; re-fly miss <= 0.116 km.
  - Verdict unchanged.
- **D1, Hollister 1H, direct route: PASS** (rule: at least one epoch).
  - Converges at 2 of 5 epochs, gate pass at worst 0.862 and 0.822.
  - Re-fly miss <= 1.4e-4 km, V_inf error <= 3.3e-11 km/s, junction mismatch <= 1.1e-11 km/s.
  - The other 3 epochs do not converge.
  - The landings span V_inf E 3.68-5.58 and V 4.78-6.41 km/s, the range of H&M orbit 1 (E 0.155 EMOS =
    4.62, V 0.179-0.206 EMOS = 5.33-6.13 km/s; descriptive).
- **ev-A** (k2|LE>V/0s|RV/1:1|LV>V/1h|LV>E/0s, 5 cycles, Standish, direct):
  - Converges and passes at 5/5, worst 0.579-0.580 (ideal 0.574).
  - Re-fly miss <= 1.9e-3 km, V_inf error <= 5.4e-10 km/s, including the Venus full-revs.
- **ev-B** (k3|RE/1:1|LE>V/0s|LV>V/1l|RV/3:2|LV>E/0s, 4 cycles, Standish, direct):
  - Passes at 2/5 (worst 0.933 and 0.757); re-fly miss <= 1.6e-2 km.
  - Gate fail at 2 epochs (1.092, 1.034); no convergence at 1.
- Limits:
  - The direct route finds one landing solution per epoch, and other solutions may exist. A
    non-converged epoch is not a negative.
  - Standish fixed mean elements are not DE440. By the 6.18 amendment, ev-A's and ev-B's full-rev legs
    on DE440 get a closure Delta-V report (blend + shoot), not a pass/fail; not yet run.
- Status: ev-A and ev-B are "candidate, pending owner adjudication, NOT novel", like ev-C.

### 6.24 DE440 closure Delta-V for full-rev rows: definition, PRE-REGISTERED 2026-10-06 (before computing)

Amendment 6.18 promised a closure Delta-V for full-rev candidates on DE440 but did not define it. On
DE440 a full-rev leg timed at the Keplerian period, from the minimax direction, returns to its
departure point while the planet has drifted by D (non-Keplerian motion over one period).

Definition (script `scripts/analyse_942_fullrev_closure_dv.py`, input: the blend-mode chain at
lambda = 1, i.e. the date solution with minimax directions):
- For each full-rev leg: D = |planet position at arrival - spacecraft position at arrival| (DOP853
  two-body re-fly), km.
- A mid-course correction at half the leg's flight time retargets the arc onto the planet at the same
  arrival time. It is solved exactly (two-body Newton on the 3 arrival-position equations) and reported
  as dv_mid in m/s.
- The induced junction error: the arrival |V_inf| changes; |change| in m/s is reported.
- Per chain: sum of dv_mid, max dv_mid and max induced |V_inf| change, per epoch.
- Descriptive only, no pass/fail. Reported next to the result of the lambda = 1 shoot (a ballistic
  closure, with its gate verdict, when one is found).
Runs: ev-A (5 cycles) and ev-B (4 cycles), DE440, 5 epochs, `--shoot-jac sparse`, default restarts.

### 6.24a Amendment to 6.24 (before the ev-B numbers were produced)

For ev-B's RV/3:2 leg (3 Venus periods, 2 spacecraft revolutions), a correction at half the leg time
leaves exactly one spacecraft revolution. The position map is then singular, and the Newton solve
failed. The correction point is moved to half a spacecraft revolution before arrival:
remaining time = 0.5 x leg time / sc_revs. For the 1:1 legs (all of ev-A's, and ev-B's RE/1:1) this is
the same point as before, so ev-A's numbers are unchanged.

### 6.25 ev-A and ev-B on DE440 (2026-10-06), data `data/942_evAB_de440_closure_dv.json`

Blend continuation reaches lambda = 1 at 5/5 epochs for ev-A and at 4/5 for ev-B (epoch 1 folds at
lambda 0.70). At lambda = 1 with minimax directions the interior gate is: ev-A 0.580 at 5/5; ev-B 0.936,
1.095, 0.758 and 1.033. The lambda = 1 shoot finds ballistic closures in every case, but all fail the
gate (ev-A 1.12-2.99; ev-B 3.7-6.9): the closing directions are far from the minimax ones. So no
gate-passing ballistic DE440 member was found (WEAK: restarts used the 0.3/0.6-rad perturbations, see 6.34-6.35); the closure Delta-V below is the measure 6.24 defines.

| | ev-A (5 cycles, 5 Venus 1:1 legs) | ev-B (4 cycles, 4 Venus 3:2 + 4 Earth 1:1 legs) |
|---|---|---|
| Full-rev miss | 5,704-8,118 km | Venus 1,310-10,358 km; Earth 28,958-1,336,854 km |
| Mid-course dv, sum per chain | 1.56-2.53 m/s (max 0.91 per leg) | 119-173 m/s (Venus 0.26-0.51 per leg; Earth 1.6-62 per leg) |
| Induced arrival |V_inf| change | <= 1.70 m/s | Venus <= 2.4; Earth up to 200 m/s |

Diagnosis of ev-B's Earth legs (checked, not assumed): on DE440 the Earth (not the Earth-Moon
barycentre) moves with the lunar reflex. Its osculating heliocentric period at the four Earth
departures of epoch 0 is 364.80, 365.77, 364.92 and 365.26 d, against the 365.25-d leg time. The 1:1 rule
|v_sc| = |V_Earth| gives the spacecraft that period. The period errors (-0.45, +0.52, -0.33, +0.01 d)
predict the misses (1.17, 1.34, 0.86 and 0.03 million km), which match the observed misses
(1.18, 1.34, 0.88 and 0.03 million km). So the Earth-leg numbers measure the seed rule on a
non-Keplerian Earth, not the cycler's need. A seed with a period matched to the Earth-Moon barycentre
would be the fair measure; it is not computed here. ev-B's Earth figures are an upper bound for this seed.

Reading:
- ev-A on DE440: a few m/s of mid-course correction over 16 yr from the minimax solution (0.580). Near-
  ballistic at the TCM level, but no ballistic gate-passing member was found.
- ev-B on DE440: Venus legs like ev-A; Earth legs not fairly measured (above).
- Standish (sec. 6.23): ev-A 5/5, ev-B 2/5, both ballistic and re-flown. That verdict stands; DE440
  carries no pass/fail by the 6.18 amendment.

### 6.26 gc-1 rung (d) on jup365, 10 cycles (lead launch 2026-10-06 14:22), data `data/943_gc1_realeph/`

- Run: blend + shoot, `--shoot-nfev-per-var 60 --shoot-jac sparse`, 20 restarts. The code is from
  4735f609, which is before the GM fix c4ff9a41; the ideal Jovian GM is 126,686,535 km^3/s^2.
- Result: the rung passes at 5/5 epochs, worst 0.755-0.761 (ideal 0.737). At every epoch only restart 0
  (the unperturbed minimax seed) closes; the other 19 do not (WEAK: restarts used the 0.3/0.6-rad perturbations, see 6.34-6.35).
- Re-fly (DOP853):
  - With the solver's GM: closure residual 1.1e-9, miss <= 2.6e-6 km (epoch 0, checked).
  - With JUP365's Jupiter-alone GM (126,686,534; 7.9e-9 lower), as the fixed checker now does: miss
    <= 0.19 km at all 5 epochs, inside the 1 km criterion. The 0.19 km is the GM sensitivity of a
    376-d chain, not a closure defect (checked by rebuilding the segments at the run-time GM).
- Status: candidate, pending owner adjudication, NOT novel. The control question for full-rev Jovian
  rows (C4 GanCal#1, sec. 6.19) is still with the lead. This result does not settle it.

### 6.27 Cell 4 (ge) adjudication and its real-ephemeris control: PRE-REGISTERED 2026-10-06, before the control runs

Author: twobody-gen2-opus (takes over from twobody-gen-opus).

Run record, checked from the files:
- 5,533 structures (shards 2,767 + 2,766), 5,437 zeros. The first 568 structures of each shard ran
  before the #963 crash (79193f6e) and the rest after the resume, so the zeros were assessed by two
  code versions. All 5,437 zeros were reassessed at 86a4d5ab (`scripts/reassess_942_zeros.py`):
  0 verdict changes, 0 errors.
- 65 raw gate-pass records (28 of them from the pre-crash part), 1,721 physical cyclers, 22
  gate-passing physical cyclers. Gauntlet (sec. 6.1 criteria): `data/943_cell_ge_gauntlet.json`.

Classification rules for the 22 (the first rule is the vm-1 ruling of sec. 6.3; the second is NEW
and goes to the lead before it is applied):
- A gate-passer with a turn of exactly 0 at one moon is a one-working-body cycler: a member of
  the R-S class (GanEur: Ganymede hosts, Europa massless; EurGan: the reverse). R-S ran both
  Jovian one-body G-E searches in this ideal model (AAS 07-118 Table 1), so these are class
  members, NOT finds.
- PROPOSED: a gate-passer whose turn at one moon is below 1 deg with a ratio below 0.05 is
  "near-one-body" and is demoted in the shortlist (like the Menning variations in sec. 6.12).
- The gauntlet's LITERAL test compares V_inf and period only. In this cell the whole EurGan
  family shares V_inf of about 2.35-2.40 / 4.05-4.10 km/s, so a LITERAL label also needs the same
  encounter count and the same return types as the R-S nomenclature (Table 5).

Real-ephemeris control for the ge cell, chosen from the page images (AAS 07-118 p.15, Fig. 10):
- GanEur#316, "10 cycles in ephemeris model, 40 G. & 10 E. flybys, start=4-24-2019,
  TOF=493.5 days, Delta-v_TOTAL=0 m/s" (Fig. 10(a) title). Nomenclature g h(1.5, 540 deg) g G:
  it has a 3-pi half-rev leg, so it tests the shot fixed-leg path at 10 cycles. It does NOT test a
  full-rev leg; the gc-1 caveat (sec. 6.19) stays for full-rev rows.
- EurGan#131, "10 cycles ..., Delta-v_TOTAL=61 m/s" (Fig. 10(b) title). It is NOT a ballistic
  published result. The R-S 2007 digest line "GanEur#316 and EurGan#131 are ballistic (Fig. 10)" is
  wrong for #131 (reported to the corpus owner). Not a control.
- Also read for context: Fig. 9(b) GanCal#1 is "10 cycles, 30 G. & 10 C. flybys, start=9-27-2013,
  TOF=375.7 days, Delta-v_TOTAL=0 m/s". The C4 runs (sec. 6.19) used the 2030-2056 epochs, not
  R-S's 2013 epoch.

Step 1, recall (`scripts/recall_943_ganeur316.py`, cell ge, k = 7, 200 structures, blind). EXPECTED
(Table 3): V_inf G/E 3.20/3.81 km/s, period 49.4 d, minimum altitude at Ganymede 1,447 km, distance
to Jupiter 592,969-1,496,829 km, transits G->E 7.60 d and E->G 12.23 d, Europa turn 0.0. PASS if a
gate-passing zero matches V_inf to 0.05 km/s, the altitude to 5 % and the distances to 0.1 %.

Step 2, chain control (`scripts/run_942_realeph_chain.py`, cell ge, jup365, blend + shoot,
`--shoot-jac sparse --shoot-nfev-per-var 60`, 20 restarts: the gc-1 settings of sec. 6.26).
- Judged at R-S's own epoch, 2019-04-24 (JD 2458597.5), 10 cycles. PASS = the sec. 6.17 criteria
  (lambda = 1 reached, every interior flyby passes the gate at the registry floors, DOP853 re-fly miss
  < 1 km at every segment). R-S claimed ballistic only at 2019, so the 5 standard epochs (2030-2056)
  are reported but do not decide the control.
- First a 1-cycle slice at the R-S epoch (validation of the settings), then the 10-cycle run (a lead
  launch).
- If it passes: the Jovian fixed-leg path is validated at 10 cycles on a published member, for
  half-rev legs. If it fails or does not decide: the ge candidates with fixed legs are "not judged" at
  rung (d).

### 6.28 Results of 6.27 so far (2026-10-06)

ge classification (22 gate-passing physical cyclers, `data/943_cell_ge_gauntlet.json`; the index is the
row order in that file):
- 16 are one-working-body (turn exactly 0 at one moon): R-S class members, NOT finds.
  - Europa hosts (EurGan type), 9: #0, #4, #5, #7, #13, #14, #17, #19, #20. #17 is EurGan#131
    (LITERAL: V_inf 2.402/4.103, 2 Europa + 1 Ganymede encounters, Europa 2:1 full-rev).
  - Ganymede hosts (GanEur type), 7: #2, #6, #8, #9, #10, #12, #18. #6 is GanEur#43 (LITERAL).
  - So the cell recovers two published rows in-run. GanEur#5 (k = 5) is outside k = 1-3.
- 6 are two-working-body (both moons turn). All six contain a full-rev leg.

| Name | Row | k | Key | V_inf G / E (km/s) | Ganymede turn, ratio | Europa turn, ratio | Worst |
|---|---|---|---|---|---|---|---|
| ge-1 | #1 | 1 | k1\|LG>E/0s\|RE/1:1\|LE>G/0s | 3.734 / 8.184 | 13.08 deg, 0.550 | 2 x 0.98 deg, 0.305 | 0.550 |
| ge-2 | #15 | 3 | k3\|LG>E/1h\|RE/3:2\|LE>G/0s | 1.371 / 1.620 | 23.94 deg, 0.291 | 2 x 30.68 deg, 0.612 | 0.612 |
| ge-3 | #21 | 3 | k3\|LG>G/1h\|LG>E/1h\|RE/1:1\|LE>G/1l | 3.881 / 8.413 | 2 x 19.71 deg, 0.883 | 2 x 0.56 deg, 0.182 | 0.883 |
| ge-4 | #11 | 2 | k2\|LG>G/1h\|LG>E/0s\|RE/1:1\|LE>G/0s | 4.033 / 8.705 | 2 x 12.48 deg, 0.595 | 2 x 0.10 deg, 0.035 | 0.595 |
| ge-5 | #16 | 3 | k3\|RG/1:1\|LG>E/0s\|RE/2:1\|LE>G/0s | 4.052 / 2.379 | 2 x 0.94 deg, 0.045 | 2 x 17.93 deg, 0.609 | 0.609 |
| ge-6 | #3 | 2 | k2\|RG/1:1\|LG>E/0s\|LE>G/0s | 4.049 / 2.347 | 2 x 0.89 deg, 0.042 | 19.22 deg, 0.640 | 0.640 |

(G = Ganymede, E = Europa in the keys.) ge-4, ge-5 and ge-6 are near-one-body by the rule proposed in
6.27 (turn under 1 deg, ratio under 0.05 at one moon). ge-5 is within 0.05 km/s of EurGan#131 in
V_inf, but it has an extra Ganymede full-rev (2 + 2 encounters against R-S's 2 + 1).
- Checks, all 22: DOP853 re-fly miss <= 2.3e-4 km, V_inf vector error <= 4.2e-9 km/s, gate pass on
  the integrated vectors. Offline literature_check: "published" for every row via the Jovian
  body-pair anchor (a flag, not a verdict).
- Collisions, two-working-body rows: no R-S row within 0.3 km/s except ge-5 (above). Clipper 11-F5
  (Buffington et al. 2012) has one switch-flip E-G-G-G-E segment (E 3.89, G 2.75-2.79 km/s), one-shot
  and not repeating: no collision. Lam et al. 2015 and the 21F31 tour (Cangahuala et al. 2025) have no
  G-E chain. Liang 2024 CGCEC (three moons; E 4.5-12.0, G 7.0-10.5 km/s) and the Hernandez 2017 and
  Lynam-Longuski 2011 Io-Europa-Ganymede triples are three-moon structures: no collision. Kumar,
  Anderson & de la Llave 2023 (G-E resonant tori, CR4BP) is a different model: context only.
- Every two-working-body ge row has a full-rev leg. So at rung (d) each carries the gc-1 caveat (no
  Jovian full-rev positive control decided at 10 cycles), whatever GanEur#316 shows.

GanEur#316 recall (`data/943_ganeur316_recall.json`):
- The blind grid at n_split 6 did NOT reach it (seed density). Seeded from R-S's printed leg times
  (`--seeded`), it is found: key k7|LG>G/1h|HG/3,1,p|LG>G/1h|LG>E/1l|LE>G/2l, dates 1.5316071495938,
  21.6589661731763, 31.0545022800922, 38.6557803419799 d. It is the only gate-passer of the 200
  structures.
- Against Table 3: V_inf 3.198/3.813 (3.20/3.81); minimum altitude 1,447.9 km (1,447); r_min 592,973
  km (592,969); Europa turn 0.0. PASS on V_inf, altitude and r_min.
- r_max: ours 1,281,581 km against 1,496,829. Cause, checked: `leg_extent` samples the Lambert legs
  only, and the half-rev's apoapsis (p = r_G, e = 0.28493: 1,070,338 / (1 - e) = 1,496,840 km) is not
  sampled. With it the value matches to 1e-5. So the r_max (and r_min) reported for every zero with a
  full-rev or half-rev leg, in every cell, covers the Lambert legs only. Descriptive values only; no
  gate or zero is affected. Not fixed yet.
- The h leg's tilt: our alpha at V_inf 3.198 is 4.14 deg; R-S print -3.98557. If that number is the
  same tilt angle (not checked), the 4 % difference is unexplained; it does not enter the checks above.

Chain control, 1-cycle slice at R-S's epoch (JD 2458597.5), `data/943_ganeur316_realeph/n1_rs2019/`:
blend reaches lambda = 1 (worst 0.671), the shoot closes (restarts 0 and 2), gate pass, worst 0.751.
DOP853 re-fly: miss 5.7e-6 km, velocity difference 9.1e-11 km/s. The 10-cycle control is a lead launch.

em reassessment check (the em run used code from before 79193f6e): every pass (24) plus 2,000 random
other zeros of the 28,090 were reassessed at HEAD. 3 changed, all "fail" to "no-directions"; no pass
or indeterminate changed. The em gate-pass set stands.

### 6.29 Lead rulings on ge, and C4 at R-S's own epoch: PRE-REGISTERED 2026-10-06, before the 10-cycle run

Lead rulings (2026-10-06):
- The near-one-body rule of 6.27 is APPROVED: ge-4, ge-5 and ge-6 are "near-one-body R-S-class
  relatives", demoted. ge shortlist: ge-1, ge-2, ge-3; ge-2 (V_inf 1.37/1.62) goes first in the
  ladder.
- The R-S 2007 digest correction (EurGan#131 is 61 m/s, not ballistic) goes to corpus-file-opus.
- GanEur#316 10-cycle runs launched by the lead: A at R-S's epoch (n10_rs2019), B at the 5 standard
  epochs (n10_std).
- `leg_extent` is to be fixed in src with a test (full-rev and half-rev legs), and the r_min/r_max of
  the reported candidates re-reported.
- C4 (GanCal#1) at R-S's own epoch: approved.

C4 at R-S's epoch. R-S Fig. 9(b): "10 cycles in ephemeris model, 30 G. & 10 C. flybys, start=9-27-2013,
TOF=375.7 days, Delta-v_TOTAL=0 m/s".
- Key k3|LGanymede>Ganymede/1l|RGanymede/2:1|LGanymede>Callisto/0s|LCallisto>Ganymede/0s, cell gc.
  Dates re-solved to 1.8e-13 km/s: 0.9914256850702445, 26.049946188944165, 35.95431049622312 d
  (V_inf 3.180/3.255, altitude 247.3 km; sec. 5).
- 1-cycle slices at the R-S epoch, done before this pre-registration:
  - Blend continuation folds at lambda 0.58 (JD 2456562.9); in a scan of 8 epochs over 36 d around it,
    it folds at lambda 0.54-0.97 at 7 and reaches lambda 1 at 1 (closure fails the gate, 1.44).
  - `--direct` (one date solve at lambda = 1 from the phase-matched ideal dates, then the shoot; the
    6.22 route): converges at lambda 1 at the 5 epochs from JD 2456550.4 to 2456562.9 and fails at the 3
    later ones. At the R-S epoch (JD 2456562.9) the shoot closes with a gate pass, worst 0.964
    (Ganymede); DOP853 re-fly miss 4.1e-5 km. Two other runs at the same phase-matched JD (other
    restart seeds) closed only with gate fails (4.67) or not at all, so the closing solution depends on
    the restart: one landing solution, not the only one (WEAK: restarts used the 0.3/0.6-rad perturbations, see 6.34-6.35).
- RUN (lead launch): `--direct`, 10 cycles, `--epochs 1 --first-epoch-jd 2456562.5`, `--shoot-jac
  sparse --shoot-nfev-per-var 60 --shoot-restarts 20`.
- PASS = the sec. 6.17 criteria at this epoch: the date solve converges at lambda = 1, the shoot closes
  (max residual < 1e-6) with a gate pass at every interior flyby, and the DOP853 re-fly miss is < 1 km
  at every segment.
- If it passes: the Jovian full-rev path has a decided 10-cycle published control, and the gc-1 caveat
  is withdrawn. If the solve or the shoot does not converge: undecided (as in 6.19), the caveat stays.
  A closure that only fails the gate: the control is NOT passed and the caveat stays; R-S's model is
  not jup365, so this is not a contradiction of R-S.

### 6.30 The GanEur#316 10-cycle launches failed: a half-rev root defect, fixed; `leg_extent` fixed (2026-10-06)

What failed: both lead launches (n10_rs2019, n10_std) stopped at once in `initial_fixed_params`
(`assert res is not None`). My 1-cycle slice (6.28) had passed. The 2-cycle chain already fails, so
the slice did not test the launch (papercut 2026-10-06-twobody-gen2-opus-one-cycle-slice-hid-chain-failure).

Cause (checked block by block):
- `_solve_half_rev_e` returned the FIRST root of the flight-time equation, scanning e up from 0.
- For a (3, 1) half-rev (and any (2k+1, k) half-rev) the body's own circle, e = 0, is also a root.
- In cycle 2 the blended Ganymede radius is 6e-4 above the circular one. That moves the circle's root
  to e = 0.002-0.005, and the scan took it in place of the leg's conic (e = 0.285).
- Result: the minimax chose a tilted near-circle and demanded 87-deg turns (ratio 2.89 from
  lambda = 0.1), then found no directions at lambda >= 0.9.
- In cycle 1 the radius happened to be below the circular one, so the slice passed.

Fixes:
- `two_working_body._solve_half_rev_e`: the circle's root lies before the flight time's minimum over
  e, the leg's after it. When roots lie on both sides, the leg's is kept.
  - Test: GanEur#316's half-rev e = 1 - r_G / r_max(Table 3) = 0.2849 at radius scales 1 -/+ 1e-3.
  - The old code gave 0.0 at a scale of 1 + 1e-6.
- `run_942_realeph_chain.initial_fixed_params`:
  - The bare assert is now an error that names the block, the body, the time, |V_inf| and the leg types.
  - Where the blended model has no minimax solution (a half-rev whose |V_inf| is below its conic's
    minimum there), the shoot is seeded from the untilted conic. This is a seed only; the shoot solves
    the directions.
- `two_working_body_enum.leg_extent` now also samples the fixed legs, from the flyby directions of
  `cycle_flybys`; `assess` passes them. Test against R-S Table 3:
  - GanCal#1: 826,589-2,415,871 km; ours was 2,357,860, now 2,415,872.
  - GanEur#316: 592,969-1,496,829 km; ours was 1,281,581, now 1,496,828.
  - Not in the test: the in-run EurGan#131 now gives 669,299-1,459,266 km, equal to Table 3 to the km
    (it was 1,458,675).

Reassessment after the half-rev fix (it moves only the assessment of zeros with an H leg; the zeros
themselves do not use the root):
- vm: 1,668 H-leg zeros, 0 changes. gc: 1,312, 0. ge: all 5,437 zeros reassessed with both fixes,
  0 verdict changes. `data/943_cell_ge_gauntlet.json` is regenerated; only r_min/r_max moved.
- ev: the first 400 of 8,464 H-leg zeros gave 12 changes, all "fail" to "no-directions" (the old
  code had taken the circle's root at round-off). The run was stopped at the 10-min limit; the rest
  is a lead launch. CORRECTION (6.34): the full ev run under this version lost 5 gate-passing
  cyclers; this version was itself a regression, replaced in 6.34.
- vm2, em, vm2n, vmn: lead launches.

Distances re-reported (ideal model; `data/942_943_extent_rereport.json`). A full-rev leg's extent is
at its minimax direction.

| Candidate | Fixed legs | r_min (old -> new) | r_max (old -> new) |
|---|---|---|---|
| gc-1 | Callisto 1:1 | 888,745 km (same) | 1,955,850 -> 2,294,675 km |
| gc-2 | none | 791,455 km | 2,337,392 km (same) |
| ev-A | Venus 1:1 | 0.5481 -> 0.5116 AU | 1.2429 AU (same) |
| ev-B | Earth 1:1, Venus 3:2 | 0.6555 -> 0.6104 AU | 1.6861 AU (same) |
| ev-C | none | 0.5090 AU | 1.6514 AU (same) |

ge rows with fixed legs moved too (in the regenerated gauntlet):
- ge-1: r_min 296,786 -> 289,471 km.
- ge-2: r_max 1,071,917 -> 1,087,498 km.
- ge-3: r_min 284,027 -> 279,867 km.

Validation with the exact launch-A flags (10 cycles, R-S epoch), after the fixes: the blend reaches
lambda = 1. The seed is used at blocks 28 and 32, and the shoot starts. Restart 0 stalls at a max
residual of 3.0e-5, like C4 in 6.19. The run continues (scratch); its verdict will be reported.

### 6.31 C4 at R-S's 2013 epoch, 10 cycles (lead launch 2026-10-06 ~16:00), data `data/943_c4_rs2013/n10_direct/`

- The direct date solve converges at lambda = 1 (JD 2456562.9); the minimax gate there is pass, worst 0.966.
- Shoot: no closure in 20 restarts. Restart 0 (the unperturbed minimax seed) stalls at a max residual of
  4.0e-6 (criterion 1e-6) after 3,074 evaluations. Four other restarts (WEAK: restarts used the 0.3/0.6-rad perturbations, see 6.34-6.35) stop at 2.8-8.6; the rest fail
  at once.
- Verdict by the 6.29 rule: UNDECIDED. The gc-1 caveat stays. Restart 0's stall has the signature of the
  6.19 long-chain stall (one cycle closes to 1e-9; ten cycles stall between 1e-6 and 1e-5).
- This run used the code from before 1feb8f2b. GanCal#1 has no half-rev leg, so the half-rev fix does
  not touch it.

### 6.32 GanEur#316 launch A (10 cycles, R-S's 2019 epoch; the exact-flag validation run adopted by the lead), data `data/943_ganeur316_realeph/n10_rs2019_v2/`

- Code 1feb8f2b (in the working tree at launch). The blend reaches lambda = 1; the seed fallback fires
  at blocks 28 and 32.
- Shoot: no closure in 20 restarts (WEAK: restarts used the 0.3/0.6-rad perturbations, see 6.34-6.35). Restart 0 stalls at 3.0e-5; restart 5 at 0.57; the others at
  2.8-1.7e3.
- Verdict by the 6.27 rule: the control does NOT decide. The ge candidates with fixed legs (all of
  ge-1, ge-2, ge-3) are "not judged" at rung (d).
- Both Jovian fixed-leg controls (C4@2013 and #316@2019) now end in the same 10-cycle shoot stall.
  They are the reproducers for the stall item (approved by the lead, after em).
- (The lead's failed first launch left `n10_rs2019/` and `n10_rs2019.log`; they are not this run.)

### 6.33 The long-chain shoot stall: diagnosis, method and control, PRE-REGISTERED 2026-10-06 (before the control run)

Lead approval: the stall item, with C4@2013 (GanCal#1, R-S Fig. 9(b): 10 cycles, 0 m/s) as the
positive control. It must not be tuned on the candidates; development uses only C4 and GanEur#316.

Diagnosis (C4@2013 restart-0 stall point, `data/943_c4_rs2013/n10_direct/`):
- The residual is smooth at large steps but has a floor of about 5e-8 (second differences at steps of
  1e-9 to 1e-8 d), from the rounding of absolute times. The shoot's date unknowns are absolute days
  (about 16,600), so the body states carry the rounding of t (about 2.4e-7 s, i.e. 2.3e-6 km of
  Ganymede position). A 14-d full-rev leg amplifies that about 25 times.
- The forward-difference step that LM and the sparse Jacobian use is sqrt(eps) x |y|. For an
  absolute date that is about 2.5e-4 d (21 s), set by the calendar, not by the problem. Over that step
  the curvature term is far larger than the 1e-6 closure level.
- Test (diagnosis only, not a verdict): from the stall point, the same LM, with the dates as offsets
  from the chain's first date, converges from 4.0e-6 to 2.1e-8 in 9 evaluations. That closure FAILS
  the gate (worst 5.11, Ganymede): the stall point had drifted from the minimax directions.

Method (`scripts/run_942_realeph_chain.py --shoot-rel-time`):
- In the shoot (phase 2) the date unknowns are offsets in days from x0. The system is wrapped by
  `RelTime`: every state is taken at t_ref + t, corrected to first order by the rounding error of
  that sum.
- Everything else is unchanged:
  - the dense or sparse forward-difference Jacobian;
  - x_scale "jac", the tolerances, the restarts and the seed;
  - the Gauss-Newton polish;
  - the closure threshold of 1e-6;
  - the gate and the re-fly.
- Outputs are converted back to absolute days.

Control runs (lead launches, with the same flags as 6.29 and 6.32 plus `--shoot-rel-time`):
1. C4@2013: cell gc, 10 cycles, `--direct`, epoch 2456562.5.
2. GanEur#316@2019: cell ge, 10 cycles, blend, epoch 2458597.5.

PASS for each = the 6.17 criteria at that epoch:
- the shoot closes (max residual < 1e-6) at some restart;
- the gate passes at every interior flyby (the best closure by worst ratio);
- the DOP853 re-fly miss is < 1 km at every segment.

Readings, fixed now:
- C4 passes: the stall item is closed. The tool may then judge Jovian fixed-leg rows, each in its own
  pre-registered run.
- Closures, but none passes the gate: the stall is fixed, and the open item becomes "the shoot finds
  closures but not the gate-passing member". The caveat stays.
- No closure: the method fails, and the caveat stays.

### 6.34 Half-rev keys: one key per geometry in every cell (replaces the 1feb8f2b rule), and the C4 result with `--shoot-rel-time` (2026-10-06)

Regression found in 1feb8f2b (by the full ev reassessment):
- 5 of ev's 31 gate-passing physical cyclers were lost, including the Hollister 1H variant
  2.994/3.19 with a half-rev pair at Venus.
- Cause: at e = 0 the flight-time equation is round-off in the ideal model, and its sign differs
  between cells. In the heliocentric cells (f(0) > 0):
  - before 1feb8f2b, a (3,1,peri) key returned the body's tilted circle and the apo key nothing, so
    the leg's own conic (e = 0.285) was NEVER assessed;
  - 1feb8f2b swapped that, so the circle was never assessed.
- In the Jovian cells (f(0) < 0) both were covered (peri = conic, apo = circle). So vm, gc and ge
  were unaffected; ev, em, vm2, vm2n and vmn have the coverage gap.

New rule (`two_working_body.half_rev_conic`):
- |f(0)| <= 1e-12 of the target counts as a root.
- Roots are split into the circle's (e < 0.05, or before the flight time's minimum) and the conic's.
- peri key = the leg's own conic through periapsis.
- apo key = its own apo conic if one exists, else the tilted circle, from whichever equation holds
  the root. The radial sign follows the equation used.
- Tests (identical in ev, em and ge, at the circular radius and 1 ulp either side):
  - (3,1,p) gives 0.28493 (GanEur#316 Table 3);
  - (3,1,a) and (1,0,a) give 0;
  - (1,0,p) gives none.
- The blend radius test (scales 1 -/+ 1e-3) is kept.

Acceptance check, ev H-leg zeros (8,464) reassessed with the new rule, against the set from before
1feb8f2b:
- gate-passing physical cyclers: 31 before, 31 after, lost 0, gained 0. The in-run Hollister 1H and
  2H controls are in the set.
- Physical cyclers in all: 3,534 -> 4,018 (the newly assessed conics). None of them passes the gate.
- Other cells under the new rule: vm (1,668 H zeros), gc (1,312) and ge (2,944): 0 verdict changes.
- em, vm2, vm2n and vmn: H-leg reassessment with this code is a lead launch. The lead's em full
  reassessment (intermediate code) stands for its non-H zeros only.
- Registry consequence: the stamps for vm2 and ev (and later em, vm2n and vmn) were made under the
  old coverage gap. The ev scope is now complete with no change of verdict. vm2 is to be re-stamped
  or annotated after its reassessment.

C4@2013 with `--shoot-rel-time` (the 6.33 control), `data/943_c4_rs2013/n10_direct_rel/`:
- Restart 0 closes in about 2 s (max residual below 1e-6; it stalled at 4.0e-6 before). The DOP853
  re-fly of that closure (checker `--include-failed`) misses by 3.8e-4 km. The closure is real.
- It FAILS the gate: worst 5.106 (Ganymede), min altitude -2,605 km. No other restart closes (WEAK: restarts used the 0.3/0.6-rad perturbations, see 6.34-6.35).
- Reading (fixed in 6.33): the stall is fixed, but the shoot finds a closure that is not the
  gate-passing member. The control is NOT passed; the gc-1 caveat stays.
- Restart defect, found after the run: the restarts perturb every fixed leg's theta by N(0, 0.3) rad
  and phi by N(0, 0.6) rad. That puts the arrivals about 4e6 km off (residual about 4e3). LM's first
  step then leaves the domain (a Lambert or flight-time failure), and the restart ends at nfev 2 with
  an infinite residual. 14 of 20 restarts died that way, so only about 6 were real attempts. Changes
  to the restarts, a chain-length continuation (1 -> 10 cycles from the 1-cycle gate-passing closure)
  or a gate-constrained shoot are method changes. Each needs pre-registration and lead approval.

### 6.35 Chain-length continuation (lead ruling (b) with (a)): PRE-REGISTERED 2026-10-06, before the control runs

Lead ruling: method (b), chain-length continuation, with (a) small-perturbation restarts; (c) only if
(b) fails. The restart defect (6.34) is fixed regardless: `--restart-sigma` defaults to 0.03,0.05 rad
(theta, phi). Statements in this note that rest on the old restarts are marked WEAK.

Method (`scripts/run_942_realeph_chain.py --grow-chain`, with `--shoot-rel-time`):
- Phase 1 (the date solve, `--direct` or blend) is unchanged and gives the seed.
- The shoot runs at 1 cycle first: the seed's first-cycle dates, the period divided by n, and the
  first cycle's minimax fixed-leg parameters.
- The gate-best closure at k cycles seeds k + 1. The last cycle's dates are moved by one cycle
  (period / k) and appended, the period is extended by one cycle, and the last cycle's (theta, phi,
  tau) are copied.
- At every length: restart 0 from the seed, then `--shoot-restarts` - 1 perturbed seeds (sigma
  0.03/0.05 rad). A length without a closure ends the run (no closure).
- Unchanged:
  - the shoot (`shoot_once`: LM, the sparse forward-difference Jacobian, the polish);
  - the 1e-6 closure threshold;
  - the gate at the registry floors; "indeterminate" is not a pass (sec. 6.1);
  - the DOP853 re-fly (< 1 km at every segment).
- Validation slice (3 cycles, 3 restarts, C4@2013, before this pre-registration):
  - k = 1 closes at restart 2 with a gate pass, 0.964;
  - k = 2 and k = 3 close at restart 0, gate "indeterminate", 0.973-0.974.
  - GanCal#1 is marginal in every model (ideal 0.961). An "indeterminate" verdict at 10 cycles is
    therefore a live possibility.

Control runs:
1. C4@2013 (cell gc, `--direct`, 10 cycles, 20 restarts, epoch 2456562.5).
2. Then GanEur#316@2019 (cell ge, blend, 10 cycles, 20 restarts, epoch 2458597.5).

Readings, fixed now:
- C4: PASS = a 10-cycle closure with a gate pass at every interior flyby and a re-fly miss < 1 km.
  Then the full-rev path is validated and gc-1 and the ge shortlist go up the ladder (each
  pre-registered).
- C4 closes at 10 cycles but only "indeterminate" (within the tidal turn scale of capacity): NOT a pass.
  It is reported as "closes, marginal as in the ideal model". The owner decides whether a marginal
  control can validate; the caveat stays until then.
- A gate fail, or no closure at some k: (b) fails, and (c) is next (with a lead ruling).

Also recorded here: GanEur#316 B (5 standard epochs, 10 cycles, absolute-time shoot, the old
restarts; descriptive by 6.27), `data/943_ganeur316_realeph/n10_std/`:
- It closes with a gate pass at 3 of 5 epochs: JD 2462503.8 (0.853), 2464844.7 (0.820) and
  2469519.5 (0.769). Each time restart 0 closes. DOP853 re-fly miss <= 2.4e-5 km.
- At 2467178.6 and 2471853.3 the best restarts stop at 1.0e-3 and 1.3e-3 (WEAK).
- So a published Jovian half-rev member closes ballistically over 10 cycles on jup365 with the
  pre-6.33 shoot at 2030-2056 epochs, but not at R-S's own 2019 epoch (6.32).

### 6.36 Cell 6 result: em (Earth-Mars, both massive), 2026-10-06

Run record:
- 9,901 structures and 28,090 zeros. The run code predates 79193f6e; it has no Traceback.
- Full reassessment at 1feb8f2b (lead launch). Its 25,178 H-leg zeros were then reassessed again under
  the 6.34 half-rev rule (12,860 verdict changes, mostly fail <-> no-directions).
- Merged set: 6,247 physical cyclers, 8 gate-passing (the same 8 as on the original assessment).
- Gauntlet: `data/942_cell_em_gauntlet.json`. Every one re-flies with DOP853 to <= 3.6e-3 km and passes
  the gate on the integrated vectors; SOI fraction <= 6.3e-9.

Positive control: NONE in-run. The pre-registration (6.2) set no in-run recall control for em.
- The nearest Russell-Ocampo catalogue rows (V_inf rounded to 0.1 km/s; no leg structures in the
  catalogue) do not match structurally.
- Our k = 3 one-body 5.151/9.153 is within 0.049 km/s of R-O 3.5.2+0, but R-O list 5 Earth flybys
  (turns 83 x 4 and 24 deg) and ours has 3 (33, 42, 75 deg).
- So the em cell is NOT validated by a published recall. Its results are conditional on that gap.
  A targeted R-O recall (leg structures from Russell 2004) is the missing control.

Classification:
- One-body (Mars turn 0; the R-O Earth-hosted free-return class, Mars massless in their ideal model;
  class members, not finds), 3:
  - k2 7.593/9.865 (NEAR R-O 2.5.1+0 by V_inf);
  - k3 5.151/9.153 (V_inf NEAR R-O 3.5.2+0, structure differs);
  - k3 5.591/9.372.
- Two-working-body (both planets turn), 5. None is near-one-body: Mars turns >= 5.6 deg, ratio >= 0.44.

| Name | k | Key | V_inf E / M (km/s) | Earth worst ratio | Mars turn, ratio | Worst | r (AU) |
|---|---|---|---|---|---|---|---|
| em-1 | 3 | k3\|RE/2:1\|LE>M/0s\|LM>M/1l\|LM>E/0s | 5.333 / 4.713 | 0.708 | 2 x 38.4 deg, 0.939 | 0.939 | 1.000-2.175 |
| em-2 | 3 | k3\|RE/1:1\|RE/1:1\|LE>M/0s\|LM>M/1l\|LM>E/0s | 5.333 / 4.713 | 0.763 | 2 x 38.4 deg, 0.939 | 0.939 | 0.903-1.930 |
| em-3 | 3 | k3\|LE>E/1l\|HE/1,0,a\|LE>M/0s\|LM>M/1l\|LM>E/0s | 4.684 / 4.539 | 0.983 | 2 x 38.4 deg, 0.892 | 0.983 | 0.880-1.914 |
| em-4 | 3 | k3\|LE>E/1l\|RE/1:1\|LE>M/0s\|RM/1:1\|LM>E/0s | 5.490 / 9.799 | 0.988 | 2 x 5.6 deg, 0.442 | 0.988 | 0.816-2.302 |
| em-5 | 3 | k3\|RE/1:1\|LE>E/1h\|LE>M/0s\|RM/1:1\|LM>E/0s | 5.977 / 10.028 | 0.996 | 2 x 5.9 deg, 0.483 | 0.996 | 0.745-2.331 |

- em-1 and em-2 share the transfer and Mars part (a 1-rev generic Mars return, 38.4-deg Mars turns)
  and differ in the Earth block. em-3's Earth half-rev is the tilted-circle geometry (6.34).
- Collision checks:
  - Rall 1969 / Rall & Hollister 1971: their M4-1, M5-1 and M5-2 have k = 4-5, two round trips per
    pattern, Mars swing-bys of 2.3-4.3 deg at about 9.35 km/s, and no direct returns at Mars. No em row
    is a member (k = 3, one round trip, Mars returns).
    - em-4 and em-5 are in Rall's energy regime at Mars (9.8-10.0 km/s, 5.6-5.9-deg Mars turns). They add
      Mars 1:1 returns, so they are Rall-adjacent, not Rall members.
  - Pisarevsky 2008: their diagrams cover loitering arcs that are multiples of pi only. Every em
    two-working-body row has a generic loitering arc (Mars generic for em-1/2/3, Earth generic for
    em-3/4/5), so all lie outside the covered diagrams. Table 4 (E 6.2 / M 5.7, k = 2) is not within
    0.3 km/s of any em row. Fig. 14 (class I.1 graphical points) is not digitised: an open check.
  - Catalogue (273 Earth-Mars rows incl. 222 R-O and 15 Rall): no V_inf match within 0.3 km/s for any
    two-working-body row.
  - Offline literature_check: "published" via the Earth-Mars anchors (a flag).
- Status: em-1..em-5 are "candidate, pending owner adjudication, NOT novel". The margins are thin
  (worst 0.939-0.996).
- Shortlist (6.12 rule; no demotions): em-1, em-2 (0.939), em-3 (0.983).
  - Real-ephemeris viability:
    - em-1 has an Earth 2:1 full-rev, and em-2 two Earth 1:1 full-revs. On DE440, Earth full-revs carry
      the lunar-reflex seed issue of 6.25; on Standish (ramp, exact) they are fine.
    - em-3's Earth leg is a tilted-circle half-rev; that leg type has no real-ephemeris control.
  - With ideal margins of 6 % or less, all three are expected to be fragile (GanCal#5 failed with a
    6 % margin, 6.16).

### 6.37 vm2n and vmn (vm2-1 neighbours, k = 2-5), 2026-10-06

- vmn (Mars massless, cell vm): 8,616 structures, 1,707 physical cyclers, 14 gate-passing, all
  one-body R-S class members. In-run control VenMar#45 is LITERAL. vm-1 is present. No finds.
  `data/942_cell_vmn_gauntlet.json`.
- vm2n (Mars massive, cell vm2): 8,616 structures, 4,476 physical cyclers, 23 gate-passing.
  `data/942_cell_vm2n_gauntlet.json`. H-leg reassessment (6.34 rule): vm2n 0 changes, vmn 0 changes.
- Two-working-body vm2n rows (Mars turns > 0), ranked:

| Name | k | Key | V_inf V / M | Worst (body) | Mars turn | r (AU) |
|---|---|---|---|---|---|---|
| vm2n-1 | 4 | k4\|LV>V/1h\|LV>M/0s\|LM>V/1h | 18.651 / 7.594 | 0.297 (V) | 5.4 deg (0.276) | 0.435-1.524 |
| vm2n-2 | 4 | k4\|RV/2:1\|RV/2:1\|LV>M/0s\|LM>V/0s | 6.221 / 4.864 | 0.421 (M) | 16.5 deg | 0.722-1.574 |
| vm2n-3 | 4 | k4\|LV>V/3h\|LV>M/0s\|LM>V/0s | 26.223 / 10.491 | 0.591 (V) | 5.0 deg | 0.249-1.524 |
| vm2n-4 | 4 | k4\|RV/2:1\|LV>M/1h\|LM>V/0s | 6.516 / 5.039 | 0.880 (M) | 32.7 deg | 0.714-1.573 |
| (k = 5 rows) | 5 | 5 rows | 6.2-6.9 / 4.9-5.1 | 0.908-0.996 | 16-34 deg | |

- vm2-1 (k3, 0.981) reappears.
- Answer to "a vm2-1 neighbour with more margin?": YES, vm2n-2.
  - It has the same V-M-V round trip (V_inf 6.22 / 4.86 against vm2-1's 6.09 / 4.85; Mars turn 16.5
    against 16.4 deg), with the Venus loiter as two 2:1 full-revs (k = 4) in place of the 2-rev generic
    return (k = 3).
  - Venus ratio 0.289 against 0.981; worst 0.421 (Mars).
  - Its one-body twin in vmn (k4, V 6.208 / M 5.871, RV/2:1 x 2) is an R-S class member.
- vm2n-1 and vm2n-3 have Venus V_inf of 18.7 and 26.2 km/s (r_min 0.435 and 0.249 AU): high-energy.
- All are "candidate, pending owner adjudication, NOT novel". Collisions: R-S has no two-body V-M rows;
  the catalogue has only the Jones VEM rows on this pair; Rall's thesis records the failed V-M search
  (6.4).
- Real-ephemeris rung for vm2n-2 (Venus 2:1 full-revs, heliocentric): the Standish direct route (6.22)
  is controlled (D1, the H&M endpoint orbits). That is a launch request with its own pre-registration.

### 6.38 6.33 control 2: GanEur#316@2019 with `--shoot-rel-time` PASSES (lead launch 2026-10-06), data `data/943_ganeur316_realeph/n10_rs2019_rel/`

- 10 cycles, R-S's epoch (JD 2458597.6), blend to lambda = 1, then the shoot (old 0.3/0.6-rad
  restarts; restart 0 is the unperturbed seed).
- Restart 0 closes (max residual < 1e-6) with a gate pass at all 49 interior flybys, worst 0.804.
  DOP853 re-fly: max arrival miss 1.5e-5 km, velocity difference 9.0e-11 km/s. The other 19 restarts do
  not close (WEAK, old restarts).
- Verdict by 6.33: PASS. This is the first decided 10-cycle Jovian control on a published member
  (R-S 2007 Fig. 10(a): ballistic, 10 cycles, start 2019-04-24). It validates the shot fixed-leg path
  for a HALF-REV leg on jup365. The full-rev path (C4) is still not validated.
- With 6.32 (the absolute-time shoot stalled at 3.0e-5 at the same epoch), this confirms the 6.33
  diagnosis on a second published member.

### 6.39 Real-ephemeris rung (d) for vm2n-2 and the em shortlist: PRE-REGISTERED 2026-10-06, before running

Route: the 6.22 direct route on the Standish J2000 mean-element model (`--real mean --direct`, ramp
mode; a full-rev leg timed at the Keplerian period is exact there). Its positive controls passed with
the same criteria:
- D1 (Hollister 1H, direct, 2/5);
- the H&M endpoint orbits 1, 2, 11, 12, 13 and 15 (6.18);
- C1 VenMar#45 (Lambert-only path).
No control covers a tilted-circle half-rev on a non-circular body. em-3's verdict is therefore
"uncontrolled leg type" whatever it shows.

Runs (5 epochs from 2030-01-01, every 6.4 yr; `scripts/run_942_realeph_chain.py`):
- vm2n-2: cell vm2, key k4|RV/2:1|RV/2:1|LV>M/0s|LM>V/0s, dates 282.4400200294492,
  500.8850093560939; 5 cycles (18.3 yr).
- em-1: cell em, k3|RE/2:1|LE>M/0s|LM>M/1l|LM>E/0s, dates 756.5892857142851, 1061.9235532084485,
  2068.7907325058313; 3 cycles (19.2 yr).
- em-2: cell em, k3|RE/1:1|RE/1:1|LE>M/0s|LM>M/1l|LM>E/0s, the same dates; 3 cycles.
- em-3: cell em, k3|LE>E/1l|HE/1,0,a|LE>M/0s|LM>M/1l|LM>E/0s, dates 33.800011351304605,
  748.8785600772824, 1062.738652758267, 2067.9756329569373; 3 cycles.

PASS (the 6.22 criteria):
- at >= 1 epoch the lambda = 1 date solve converges (< 1e-6);
- every interior flyby passes the gate at the registry floors ("indeterminate" is not a pass);
- the DOP853 re-fly against the mean-element system misses by < 1 km with V_inf vector error
  < 1e-6 km/s.

A non-converged epoch is not a negative. All results are on Standish fixed mean elements, not DE440.

### 6.40 Result of 6.39 (2026-10-06), data `data/942_rung_d_vm2n_em/`

| Candidate | Epochs converged (direct) | Gate at the converged landings | Rung (d) |
|---|---|---|---|
| vm2n-2 | 1/5 (JD 2472014.4) | fail: Mars 5.91, Venus 1.47 | NOT PASSED |
| em-1 | 1/5 (JD 2472366.1) | fail: Mars 2.94, Earth 1.07 | NOT PASSED |
| em-2 | 1/5 (same epoch and dates as em-1) | fail: Mars 2.94, Earth 0.77 | NOT PASSED |
| em-3 | 3/5 | fail: worst 3.2-24.3 (Mars) | NOT PASSED (uncontrolled leg type anyway) |

- By 6.22 a non-converged epoch is not a negative, and the direct route finds one landing per epoch.
  So these are "not passed by the controlled route", not proofs of absence.
- At every converged landing the MARS flyby fails. Mars's real eccentricity (0.093) moves the V_inf
  at Mars far from the ideal values, and the thin ideal margins (em: Mars 0.89-0.94) do not survive,
  as with vm2-1 (6.9).
- vm2n-2's ideal Mars margin (0.421) did not help at the one converged landing either.
- A ramp-continuation run (6.18) could land on other solutions; that is not pre-registered and not run.
- Status: vm2n-2 and em-1..3 stay "candidate, pending owner adjudication, NOT novel", recorded as
  ideal-model members with no controlled real-ephemeris pass.

### 6.41 Prior-art results received 2026-10-06 (from the lead and corpus-file-opus)

- ev-C: D. Ross, "A cycler-quartet between Venus and Sol/Terra L1" (unrefereed manuscript; possible
  class collision, see 0ba585b2). Its 2L4 is powered, never reaches Earth (aphelion 0.979 AU; it
  targets Sun-Earth L1), needs a Venus flyby below the surface, and has a Venus V_inf of 3.86 km/s.
  NO collision with ev-C; it is the nearest unrefereed class relative. Search note:
  `docs/notes/2026-10-06-942-evC-prior-art-search.md`.
- gc-1, gc-2, ge-1, ge-2, ge-3: Campagnola et al., ISSFD 2024 paper 19-3 (the conference form of
  JAS 2025, the 21F31 reference trajectory).
  - NO repeating G-C or G-E segment. "Cycler" appears once, as an untried option.
  - The flown tour uses a Callisto petal rotation: Callisto V_inf 4.8 -> 3.6 -> 5.1, Ganymede 5.35 /
    4.30.
  - G-E appears only in the one-shot pump-down (Ganymede 7.65-8.73, Europa 6.30-6.38 km/s, read on
    the page image of Cangahuala 2025 Table 2).
  - Nearest: the petal-rotation Callisto V_inf of 3.6 is 0.56 km/s above gc-2's. No collision.
- gc-1, gc-2: Golubev, Grushevskii, Koryanov & Tuchin 2014 (JCSSI 53(3):445). Its "crossed" G-C-G
  manoeuvres are one-shot (no cycler, no periodic orbit, no G-C V_inf values). No collision.
- Still open (lead log): the web prior-art search for ge-*, em-* and ev-A.

### 6.42 Lead rulings on 6.40, and the em recall control: PRE-REGISTERED 2026-10-06, before running

Ruling on 6.40 (lead): no ramp continuation now. vm2n-2 and em-1..3 are recorded as "ideal-model;
rung (d) not passed at the converged landings, all failing at Mars (e = 0.093); other landings
untested". A ramp continuation stays an open option.

em recall control (approved by the lead). The em cell's scope holds no published Russell-Ocampo row:
- 2.5.1.+0 has three Earth returns;
- Byrnes case 3 (2.3.1.+1) has a 1-rev transfer;
- the rest have more returns or k = 4.
So the rows' own structures are solved with the em cell's model and the production solver
(`scripts/recall_942_ro.py`; every cyclic order of the Earth block, k = 2, n_phase 36, n_split 12,
n_refine 40). Structures are from the McConaghy, Russell & Longuski 2005 Table 2 labels; expected
values are from Russell 2004 Table 3.4. Turn ratio there = max / required, so ours (required / max)
is its inverse.

- 2.5.1.+0 = g(1-11/14 yr, 11/14 rev) f(1:1) h(0.5 yr, 0 rev, tilt 15.081 deg) f(1:1).
  EXPECTED: V_inf E 7.8 / M 9.9; E->M 94 d; three Earth turns of 54 deg; our ratio 1/1.12 = 0.893.
  Gate PASS.
- 2.3.1.+1 (Byrnes case 3) = g(2-11/14 yr, 1-11/14 rev) f(1:1) h(0.5 yr, 0 rev, tilt 10.388 deg).
  EXPECTED: V_inf E 5.4 / M 5.3; E->M 143 d; turns 93, 93 deg; our ratio 1/0.92 = 1.087. Gate FAIL
  (R-O list it as near-ballistic).

CONTROL PASS = a zero of each structure reproduces:
- V_inf within 0.05 km/s (the table's rounding);
- every Earth turn within 1 deg;
- the worst ratio within 0.02 of the expected value;
- 2.5.1.+0 passes the gate.
If it passes, em's results stop being conditional on the control gap. If not, the gap stays and the
cause is reported.

### 6.43 Result of the em recall control (2026-10-06), data `data/942_em_ro_recall.json`

| Row | Our best zero | V_inf E / M | Earth turns (deg) | Worst ratio | Gate | Against the pre-registration |
|---|---|---|---|---|---|---|
| 2.5.1.+0 | k2\|RE/1:1\|HE/1,0,a\|RE/1:1\|LE>M/0s\|LM>E/0s | 7.817 / 9.945 (7.8 / 9.9) | 53.56 x 4 (R-O: 54, 54, 54) | 0.897 (0.893) | PASS | PASS on every criterion |
| 2.3.1.+1 (Byrnes case 3) | k2\|HE/1,0,a\|RE/1:1\|LE>M/0s\|LM>E/1l (and its 0+1 twins) | 5.393 / 5.314 (5.4 / 5.3) | 92.85, 89.53, 38.49 (R-O: 93, 93) | 1.092 (1.087) | fail (as published: near-ballistic) | V_inf and ratio PASS; turns FAIL (89.5 against 93; one extra value) |

Notes:
- 2.5.1.+0 is found as the time-reversed twin. Our Mars-to-Earth leg is 93.8 d, R-O's Earth-to-Mars
  time is 94 d; the analysis merges mirror twins (6.1).
- R-O list one turn fewer than our encounter count in both rows (3 for 4, and 2 for 3). For 2.5.1.+0
  all four of ours equal 53.56 deg, so the mapping is immaterial.
- For case 3 the 1:1 full-rev direction is free (lambda in the McConaghy label). Our minimax picks a
  different split of the turns than R-O's tabulated choice, at the same worst ratio (1.092 against
  1.087).

Reading:
- The published BALLISTIC row (2.5.1.+0) is recovered blind, with V_inf, turns, worst ratio and the
  gate pass all reproduced. This is the em-cell positive control for a gate pass.
- The near-ballistic case 3 reproduces V_inf and the worst ratio, but not R-O's per-flyby turn split.
  The pre-registered all-criteria rule is not met for that row.
- Recommendation to the lead: count the control as passed for the gate-pass question (2.5.1.+0), so
  em-1..5 are no longer conditional on a control gap. Keep the case-3 turn mismatch on record.
- The em cell's own scope (<= 2 returns per block, 0-rev transfers) excludes both rows; they were
  solved by the same code path outside that scope, as GanEur#5 was for ge (sec. 5).

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
