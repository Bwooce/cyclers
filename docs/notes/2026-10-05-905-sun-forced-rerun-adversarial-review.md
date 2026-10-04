# #905 Sun-forced Earth-Moon rerun: adversarial review

Date: 2026-10-05. Subject: `data/found/905_sun_forced_rerun/` (launched at commit `84f09e75`; driver
`src/cyclerfinder/search/sun_forced_905.py`, run script `scripts/run_905_sun_forced_rerun.py`).
Reviewer scripts (scratch, not committed): `verify.py` (independent re-integration), `compare.py` /
`uniq.py` (duplicates, symmetry classes, model pairing), `lb.py` (distance to Leiva & Briozzo 2008),
`litcheck.py` (literature check).

## 1. Verdict

- The 114 eps = 1 branch records stored by the run are 72 distinct orbits: 36 in the bicircular
  model, 36 in the quasi-bicircular (QBCP) model. Two orbits count as the same when they are related
  by a whole-Sun-period shift or by the s1, s2 or s3 reflection, so each of the 72 is one
  symmetry class, not one trajectory. Three pairs of parent rows are duplicates of each other:
  Braik-Ross C11a with Ross-RT 1:1, Braik-Ross C32 with Ross-RT 3:2, and C21 spatial #438 with
  C21 3D corridor #682. So 42 of the 114 records repeat another record.
- All 72 pass the independent checks: closure, Sun phase, surfaces, minimal period, reciprocal
  multipliers and planarity (section 3). None is an artefact. Class (d) = 0.
- Class (a), a reproduction of a published Sun-forced orbit: 0 among the run's eps = 1 orbits. The
  driver's separate `controls/` stage reproduces Oshima 2022 and Leiva & Briozzo 013_t3/t4; that
  was reviewed in the #905 build and is not repeated here.
- Class (b), a completion of a published gap: **2 orbits, both in the QBCP**. These are the
  closures of Leiva & Briozzo's (2008) 5/2 periodic arcs 180A_1_t1/t2 (C32, C = 3.18010) and
  357_t3/t4 (C31, C = 3.05583). In each case both of the paper's arc epochs lie on one of our
  orbits, and the multipliers match the paper (section 4).
- Class (c), a Sun-forced counterpart of a published three-body family not found in the held
  literature: 70 orbits. Here "not found in the held literature" does not mean new to the
  project. The #884 reviewer's corrected-sense table already contains the bicircular C32 5/2
  values to every quoted digit (OUTSTANDING #905), and the 5/2 QBCP orbits were already computed
  as #905 controls. Of the 70, 53 have a periselene inside the lunar Laplace sphere (66,183 km).
  The other 17 are R21-S (a near-circular Earth orbit, not a cycler), R52-S 2/1 and 5/2, and the
  C21 3D 8/3 orbits, with periselenes of 70,000 to 196,000 km.
- **Record to correct:** the 180A_2 claim in OUTSTANDING #905 ("Leiva & Briozzo's 5/2 C32/C31
  members ... COMPLETE", 180A_2 at 17,750 and 22,239 km) is not a completion. Neither 180A_2 orbit
  passes within 1e-2 length units (3,800 km) of the published 180A_2_t1 arc state, at any
  half-Sun-period epoch. Their multipliers (4.6e6 and 1.2e7) do not match the arc either: s1^2 =
  5153^2 = 2.7e7. The two orbits have the same three-body parent as a published arc but are
  different objects, so they are class (c).
- Nothing here is a catalogue row yet (section 8).

## 2. What was checked, and how

- **Independent dynamics.** Each stored eps = 1 orbit was re-integrated with `core/bcr4bp.py`
  (`bcr4bp_eom`, `bcr4bp_stm_eom`) or `core/qbcp.py` (`qbcp_eom`, `qbcp_stm_eom`, PV to PM via
  `state_pv_to_pm`). These are pure-Python fields, not the driver's numba field. The integrator was
  scipy DOP853 with rtol 1e-13 and atol 1e-14. The driver supplied only the node times (from
  `node_offsets` on the stored parent); a wrong time would show up as a large mismatch.
- **Closure.** The test is node-to-node multiple shooting over P = N Tg exactly (Tg = 2 pi /
  omega_S), so the last segment ends with the Sun at its starting phase. The single-shot closure
  from node 0 is also recorded, but at multipliers of 1e3 to 1e15 it measures only how the
  rounding of node 0 is amplified.
- **Passes.** Every local minimum of the distance to the Moon (at 1 - mu) and to the Earth (at
  -mu) was found on the dense output, sampled at 6,000 points per TU and refined by bounded
  minimisation. Lunar surface is 1,737 km; Earth surface is 6,378 km.
- **Minimal period.** The full state was compared with node 0 at tau + j Tg for j = 1 .. N-1. The
  smallest return is 0.28 length units or more for every orbit, so none is a multiple cover.
- **Multipliers.** The monodromy is the product of the core-model segment STMs. A backward
  integration of every segment gives M^-1 independently. The largest |lambda| agrees with the
  largest |eig(M^-1)| to 1e-5 for all 114 records, which confirms the reciprocal pairing for the
  dominant pair. For the subdominant pairs, the forward and backward values disagree by up to a
  factor 3.8 when |lambda_max| > 1e12 (13 records: C32 8/3 and C31 8/3). Those subdominant
  multipliers are numerically meaningless; the driver's note says "unreliable above about 1e6".
- **Planarity.** max |z| is below 1e-19 for every planar orbit. The C21 3D orbits reach |z| of
  0.21 to 0.26.
- **Jacobi-like drift.** The forced models have no integral. The three-body Jacobi constant
  evaluated on the orbit varies over a band of width 0.005 to 0.08 in the bicircular model and
  0.02 to 0.3 in the QBCP. The QBCP band is larger because its coordinates pulsate with the
  Earth-Moon distance (alpha_1 has a 1.4 % harmonic). This is reported, not used as a pass/fail.
- **Duplicates and model pairing.** Positions were sampled on the absolute clock grid (128 per
  Tg) and compared under k Tg shifts and the s1, s2, s3 images. QBCP samples were shifted by
  Tg/2 to align the Sun direction: the bicircular Sun is on +x at t = 0, the QBCP Sun on -x.
- **Published orbits.** Each QBCP orbit with N = 3, 4 or 5 was compared with all 17 usable Leiva &
  Briozzo (2008) Table 2-3 states. Paper frame to project frame is a sign flip of x, y, xdot and
  ydot; the state is taken at QBCP clock t_i (reading B, per the #896 addendum). The comparison
  was made at t_i + k Tg/2, k = 0 .. 2N-1, so an arc's return epoch is included. The other held
  Sun-forced orbits (Oshima 2022, Boudad-Howell-Davis 2020, Leiva & Briozzo 2005, Rosales et al.
  2021/2023, Andreu POL1/POL2, Jorba et al. 2020) were compared by period, geometry and
  multipliers, from the digests.

## 3. Failures and flags found

1. **The QBCP orbits are closed over M T*, not over N Tg (59 of 59 QBCP branches).** The driver's
   own `sun_phase_returns` is False for every QBCP branch and True for every bicircular one, and
   `summary.json` does not show this. For the C31 8/3 QBCP orbit, the period is shorter than
   N Tg by 1.2e-10 TU.
   - Over M T*, the node mismatch is 8e-13 in every segment.
   - Over N Tg, the mismatch is 1.7e-9, and all of it is in the last segment.
   So the driver's numba QBCP field agrees with `core/qbcp.py`, and the defect is only the period
   used. It is negligible physically, but the strict claim "N Sun periods with the Sun back at its
   starting phase" holds to 1e-9 in the QBCP, against 1e-12 in the bicircular model. One Newton
   step at P = N Tg would remove it.
2. **Duplicate parents.** Three pairs of parent rows give the same orbits (section 1), so the
   58 summary rows overstate the census.
   - The duplicates' parent states differ by 1e-5 to 2e-5 at the same period. This is a phase
     shift of the reference point (tau differs by 4e-5).
   - The near-identical C21 3D member-0 parents still give different outcomes: bicircular 2
     versus 1 branches reached; QBCP 0 versus 1. The continuation is fragile at the 1e-5 level of
     its starting point, which is a reproducibility flag.
3. **`after_bifurcation` depends on the walk path, not on the orbit.** It is False for Braik-Ross
   C32 m0 and True for the identical Ross-RT 3:2 m3 orbit (0 against 14 recorded events before
   the member). The flag also fires on in-plane index crossings with m = 1, k = 0, a crossing of
   the index s = 2: the family changes stability and may branch (pitchfork or tangent bifurcation).
   No check in the run shows that the walk stayed on the catalogued family after such a crossing.
   For C32 5/2 and C31 5/2, the parents' Jacobi constants match Leiva & Briozzo Table 1 to 8
   digits, so those parents are right. For the others it is unverified.
4. **Branch points crossed.** 19 of the 72 distinct orbits recorded one or two branch points
   (sign changes of the bordered determinant) between eps = 0 and 1. The eps = 1 orbits are valid,
   but whether they are the continuation of the stated Melnikov zero, or of a branch the path
   switched onto, is not established.
5. **Bicircular/QBCP disagreements.** Most orbits have a counterpart in the other model within
   0.006 to 0.06 length units. That is the expected size, given the QBCP's pulsating coordinates
   and the model difference. The exceptions:
   - C32 5/2 at C = 3.15168: the bicircular orbit with periselene 23,199 km and the QBCP orbit
     with 17,750 km have no counterpart in the other model. The nearest are at 0.12 and 0.56.
   - C21 planar 3/1 at C = 3.12523: the bicircular orbit with periselene 26,066 km (|lambda|
     4.0e3) has no QBCP counterpart. Its QBCP branches folded back below eps = -0.5.
   - C21 3D 5/2 at C = 3.01251: the two QBCP orbits (|lambda| 62 and 76, periselene 58,271 and
     65,938 km) have no bicircular counterpart. Both bicircular branches stopped at eps of about
     1e-4 to 1e-5 (step collapse, or maximum steps with 80 folds). This is unexplained and almost
     certainly numerical: for a simple zero, Rhouma-Chicone guarantees continuation at small eps.
   - C21 3D 8/3: one bicircular class has no QBCP partner. The QBCP walk of the #682 duplicate
     stalled at eps = 1.5e-8, with 150 alternating "folds" of size 1e-8. That is chattering, not
     folds.
   - C32 8/3: the pairs are 0.13 to 0.14 apart. That is plausible at |lambda| about 6e15, but not
     shown to be the same branch.
6. **8/3 members.** Leiva & Briozzo's first-order (quadrupole) condition vanishes identically for
   q = 3 (their Eq. 28). The driver's Melnikov amplitudes for 8/3 are 1e-5 to 4e-4, against 1e-2 to
   1e-1 for 5/2 and 3/2, so the zeros come from the octupole and higher terms. Several 8/3 orbits
   reached eps = 1 anyway. The persistence radius is unknown, and Brown et al.'s higher-order
   functions were not used.
7. **Near-twin orbits.** Many parents give two orbits whose Sun phases differ by pi, with
   periselenes that differ by 2 to 30 km (for example C11 2/1: 24,999 and 25,005 km). They are
   distinct solutions, closed to 1e-13. The likely explanation, my inference: the Sun's
   quadrupole tide is pi-periodic in the Sun angle, and only the octupole (about 1/389 smaller)
   separates the twins. Counting them as two orbits is correct but inflates the impression of
   variety.
8. **Literature check.** `search/literature_check.py` (offline corpus, `check_literature`) returns
   "published" for every signature tried (C32, C31, C11, C21 planar and spatial, R52-S, R21-S),
   each time citing Braik & Ross 2026. The match is on the three-body parent family. The tool has
   no Sun-model dimension and no Leiva & Briozzo, Oshima, Boudad or Brown anchors. Its verdict is
   correct for the parents but cannot tell whether a Sun-forced counterpart is new.
   `is_novelty_claimable` is False for all.

## 4. The two class (b) orbits

| Orbit (QBCP) | Matches arc | Distance at the arc epochs (pos / vel, LU / VU) | Our |lambda_1|, |lambda_2| over 5 Tg | Paper s1^2, s2^2 (s over 5/2 Tg) | Periselene: ours / paper d_M | Perigee: ours / paper d_E |
|---|---|---|---|---|---|---|
| C32, C = 3.18010 (RT-32 m1 b1 = BR-C32 m2 b1) | 180A_1_t1 and _t2 (both arc epochs lie on this one orbit, Tg/2 apart) | 2.2e-4 / 4.0e-4 (85 km, 0.4 m/s) | 2467, 151.7 | 49.9^2 = 2490, 12.4^2 = 154 | 7,348 / 7,371 and 7,365 km | 118,464 / 118,458 km |
| C31, C = 3.05583 (RT-31 m0 b0) | 357_t3 and _t4 | 4.3e-4 / 1.1e-3 (165 km, 1.1 m/s) | 2.89e7, 24.9 | 5366^2 = 2.88e7, 5.2^2 = 27 | 14,715 / 14,713 and 14,716 km | 96,783 / 96,931 km |

The residual distance is what separates an arc (a fixed point of the 5/2 Tg map with the Sun at
phase pi) from the closed 5 Tg orbit next to it. It also includes the paper's mu (0.0121505482)
against ours (0.0121505816).

- The bicircular counterparts are at 0.020 and 0.014 length units: periselene 7,670 km with
  |lambda| 3.2e3, and periselene 14,903 km with |lambda| 2.6e7. They are class (c), since Leiva &
  Briozzo worked only in the QBCP.
- The other Sun phase of each parent (180A_1: periselene 9,933 km, |lambda| 3.9e4; 357:
  13,871 km, 4.6e7) is not in the paper, which reports only one pair of epochs per family. These
  are class (c).
- **The strongest class (b) candidate is the C32 180A_1 closure:** QBCP, 5 Sun periods (147.6 d),
  periselene 7,348 km, perigee 118,464 km, |lambda| 2,467 per period (an e-folding time of
  18.9 d), closure 6e-10 over N Tg.

## 5. Usefulness and distance from a real orbit

- Both models are idealised: a circular or quasi-bicircular Sun, a circular or Fourier-fitted
  Moon, no lunar eccentricity (0.055), no lunar inclination (5.1 deg). Nothing was continued to an
  ephemeris.
- The spread between the bicircular and QBCP versions of the same orbit is 0.006 to 0.14 length
  units (2,000 to 54,000 km). That is a lower bound on how far the real-ephemeris counterpart
  could sit from either one.
- Instability:
  - Doubling times over the cycler-class orbits range from 4.5 d (|lambda| about 6e15 over
    8 Sun periods) to about 13 d (180A_1 closure).
  - The 8/3 C32 and C31 and the 5/2 C11 orbits (|lambda| 1e12 to 6e15) cannot be flown without
    continual control and are of mathematical interest only.
  - The least unstable lunar-pass orbits are the QBCP C21 3D 5/2 pair (|lambda| 62 to 76 per
    147.6 d, doubling about 24 d, periselene 58,000 to 66,000 km). They have no bicircular
    counterpart.
  - Next are R52-S 2/1 (|lambda| 7.6 to 9.6 per 59 d, periselene 70,000 to 82,000 km, outside
    the Laplace sphere) and C21 3D 8/3 (|lambda| 30 to 35, periselene 120,000 km).
- Every perigee is 76,000 km or more, about 12 Earth radii. None of these orbits offers
  low-Earth access. Genova & Aldrin's 3-petal real-ephemeris cycler (3,000 km perigee and
  perilune altitudes, 39 m/s per month of station-keeping) is a different object, and no
  bicircular version of it is held.
- R21-S 1/2 is a near-circular prograde Earth orbit at 187,000 to 193,000 km with period Tg (two
  laps). Its multipliers are on the unit circle or within 1.04 to 1.09. It is not a cycler and has
  no lunar pass inside 194,000 km. It is distinct from Oshima's 1:1 spatial Earth-centred orbits
  (at 1.09 to 1.12 length units), from Leiva & Briozzo 2005 (lunar pass at about 12,000 km) and
  from the Rosales L1/L2 replacement orbits. The held corpus has no source for it.

## 6. Cases that did not reach eps = 1, all accounted for

- **R21-S m0, bicircular:** Melnikov identically zero. The theorem says nothing, and the run did
  not continue it.
- **R21-S m0, QBCP (5/6):** step collapse at eps = 6e-4, and maximum steps at eps = -2e-5.
- **Vaquero 2:1 (5/6, q = 6):**
  - bicircular: step collapse at eps = 3e-4 and 1e-4.
  - QBCP: folds back below eps = -0.5.
  - The first-order function is flat for q > 2.
- **C32 m3 (C = 3.18264, 8/3, parent periselene 4,457 km):** both models stall at eps = 0.004 to
  0.007 (a fold at 0.0086 in one). This is a low-periselene parent with an octupole-order
  Melnikov amplitude.
- **C21 3D 5/2, bicircular:** stalls at eps of about 1e-4 (section 3, item 5). Unexplained.
- **C21 3D 8/3, QBCP (#682 copy):** chatters at eps = 1.5e-8. The #438 copy reached eps = 1 on
  one branch. This is fragility, not non-existence.
- **C21 planar 3/1 m0, QBCP:** two of the four branches fold to eps < -0.5.
- **Casoliva 7:3 (b) and (c):** no catalogue state or perpendicular crossing, so not run. The
  other Casoliva rows have no commensurate member, which is why none appears in `forced/`.

## 7. All eps = 1 orbits

Column guide:
- class: (b) completes a published arc; (c) Sun-forced counterpart not in the held literature;
  "x" means no counterpart in the other model within 0.1 length units.
- cyc: periselene at or below 66,183 km.
- closure: the maximum node mismatch over N Tg, from the core models. The QBCP values of 1e-10 to
  2e-9 are the M T* period effect of section 3, item 1.
- doubling: ln 2 / ln |lambda| times the period, in days.
- min return jTg: smallest return distance of the state at tau + j Tg, j < N.
- bp: branch points recorded during the continuation.
- counterpart: position distance in length units after Sun-direction alignment, with the
  counterpart's periselene.
- nearest published: Leiva & Briozzo 2008 row (Tables 2-3; QBCP only, N = 3 or 5), with the
  position / velocity distance at the best arc epoch.
- "none held at this N": the digests hold no Sun-forced orbit with that period and geometry.
  Oshima's orbits are 1:1 and 12:11 Earth-centred retrograde; Boudad's are NRHOs; Rosales's and
  Andreu's are L1/L2 replacements; Leiva & Briozzo 2005 is a 1-Tsun high-energy orbit.

Run aliases: RT = Ross-RT, BR = Braik-Ross, m = member, b = branch index in the forced file.

| # | model | N/M | parent C | runs (family member branch) | class | cyc | closure | periselene km | perigee km | max abs multiplier | doubling d | min return jTg | max abs z | bp | other-model counterpart, LD | nearest published (L&B 2008), pos/vel |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | bcr4bp | 1/2 | 3.42072 | BR-R21S m1 b1 | (c) | no | 2e-14 | 195395 | 187853 | 1.04e+00 | 555.5 | n/a | 0.00 | 0 | 0.006 (pk 196585) | none held at this N |
| 2 | bcr4bp | 1/2 | 3.42072 | BR-R21S m1 b0 | (c) | no | 2e-14 | 195680 | 188217 | 1.00e+00 | stable | n/a | 0.00 | 0 | 0.006 (pk 194011) | none held at this N |
| 3 | bcr4bp | 2/1 | 3.19716 | BR-R52S m0 b1, BR-R52S m0 b3 | (c) | no | 2e-13 | 70196 | 126833 | 7.60e+00 | 20.2 | 1.76 | 0.00 | 0 | 0.017 (pk 70421) | none held at this N |
| 4 | bcr4bp | 2/1 | 3.19716 | BR-R52S m0 b2 | (c) | no | 2e-13 | 80217 | 129867 | 9.26e+00 | 18.4 | 1.81 | 0.00 | 0 | 0.026 (pk 81578) | none held at this N |
| 5 | bcr4bp | 2/1 | 3.19716 | BR-R52S m0 b0 | (c) | no | 2e-13 | 80243 | 129878 | 9.26e+00 | 18.4 | 1.81 | 0.00 | 0 | 0.024 (pk 81425) | none held at this N |
| 6 | bcr4bp | 5/2 | 3.18010 | BR-C32 m2 b1, RT-32 m1 b1 | (c) | yes | 5e-12 | 7670 | 119571 | 3.15e+03 | 12.7 | 1.74 | 0.00 | 0 | 0.020 (pk 7348) | none held at this N |
| 7 | bcr4bp | 5/2 | 3.18010 | BR-C32 m2 b0, RT-32 m1 b0 | (c) | yes | 7e-13 | 10048 | 117150 | 2.20e+04 | 10.2 | 1.78 | 0.00 | 1 | 0.014 (pk 9933) | none held at this N |
| 8 | bcr4bp | 5/2 | 3.15168 | BR-C32 m1 b0, RT-32 m2 b0 | (c) | yes | 6e-13 | 21806 | 108153 | 9.76e+06 | 6.4 | 1.82 | 0.00 | 1 | 0.015 (pk 22239) | none held at this N |
| 9 | bcr4bp | 5/2 | 3.15168 | BR-C32 m1 b1, RT-32 m2 b1 | (c) x | yes | 6e-13 | 23199 | 109013 | 1.95e+07 | 6.1 | 1.82 | 0.00 | 1 | 0.124 (pk 7348) | none held at this N |
| 10 | bcr4bp | 3/2 | 3.15107 | BR-C11a m1 b1, RT-11 m2 b1 | (c) | yes | 1e-13 | 21080 | 187062 | 6.56e+05 | 4.6 | 1.60 | 0.00 | 0 | 0.020 (pk 22614) | none held at this N |
| 11 | bcr4bp | 3/2 | 3.15107 | BR-C11a m1 b0, RT-11 m2 b0 | (c) | yes | 3e-12 | 25104 | 180265 | 3.86e+07 | 3.5 | 1.62 | 0.00 | 0 | 0.038 (pk 24390) | none held at this N |
| 12 | bcr4bp | 8/3 | 3.12875 | BR-C32 m0 b1, RT-32 m3 b1 | (c) x | yes | 6e-13 | 24496 | 98509 | 6.37e+15 | 4.5 | 0.67 | 0.00 | 1 | 0.134 (pk 23996) | none held at this N |
| 13 | bcr4bp | 8/3 | 3.12875 | BR-C32 m0 b0, RT-32 m3 b0 | (c) x | yes | 1e-12 | 24504 | 98519 | 6.40e+15 | 4.5 | 0.67 | 0.00 | 0 | 0.141 (pk 23996) | none held at this N |
| 14 | bcr4bp | 3/1 | 3.12523 | RT-21 m0 b2 | (c) | yes | 1e-13 | 25103 | 221975 | 5.58e+04 | 5.6 | 0.43 | 0.00 | 0 | 0.034 (pk 25791) | none held at this N |
| 15 | bcr4bp | 3/1 | 3.12523 | RT-21 m0 b0 | (c) | yes | 1e-13 | 25109 | 221974 | 5.58e+04 | 5.6 | 0.43 | 0.00 | 0 | 0.035 (pk 25786) | none held at this N |
| 16 | bcr4bp | 3/1 | 3.12523 | RT-21 m0 b1, RT-21 m0 b3 | (c) x | yes | 1e-13 | 26066 | 215381 | 4.02e+03 | 7.4 | 0.39 | 0.00 | 0 | 0.570 (pk 20927) | none held at this N |
| 17 | bcr4bp | 2/1 | 3.12107 | BR-C11a m2 b0, RT-11 m0 b0 | (c) | yes | 5e-14 | 24999 | 227567 | 3.44e+05 | 3.2 | 1.82 | 0.00 | 0 | 0.020 (pk 25608) | none held at this N |
| 18 | bcr4bp | 2/1 | 3.12107 | BR-C11a m2 b2, RT-11 m0 b2 | (c) | yes | 7e-14 | 25005 | 227581 | 3.44e+05 | 3.2 | 1.82 | 0.00 | 0 | 0.020 (pk 25602) | none held at this N |
| 19 | bcr4bp | 2/1 | 3.12107 | BR-C11a m2 b1, BR-C11a m2 b3, RT-11 m0 b1, RT-11 m0 b3 | (c) | yes | 7e-14 | 25931 | 223037 | 4.24e+05 | 3.2 | 1.82 | 0.00 | 0 | 0.026 (pk 26065) | none held at this N |
| 20 | bcr4bp | 3/2 | 3.09204 | BR-C11a m0 b1, RT-11 m3 b1 | (c) | yes | 4e-13 | 20531 | 138751 | 2.34e+10 | 2.6 | 1.76 | 0.00 | 0 | 0.016 (pk 21171) | none held at this N |
| 21 | bcr4bp | 3/2 | 3.09204 | BR-C11a m0 b0, RT-11 m3 b0 | (c) | yes | 3e-13 | 23396 | 145860 | 1.47e+10 | 2.6 | 1.73 | 0.00 | 0 | 0.008 (pk 23142) | none held at this N |
| 22 | bcr4bp | 5/2 | 3.07692 | RT-11 m1 b1 | (c) | yes | 1e-12 | 18514 | 213882 | 8.89e+13 | 3.2 | 1.54 | 0.00 | 0 | 0.015 (pk 19668) | none held at this N |
| 23 | bcr4bp | 5/2 | 3.07692 | RT-11 m1 b0 | (c) | yes | 6e-14 | 21010 | 236253 | 1.00e+12 | 3.7 | 1.38 | 0.00 | 0 | 0.063 (pk 20323) | none held at this N |
| 24 | bcr4bp | 3/1 | 3.07518 | RT-21 m1 b2 | (c) | yes | 2e-13 | 17699 | 159551 | 4.25e+07 | 3.5 | 0.49 | 0.00 | 0 | 0.013 (pk 18249) | none held at this N |
| 25 | bcr4bp | 3/1 | 3.07518 | RT-21 m1 b0 | (c) | yes | 2e-13 | 17701 | 159548 | 4.25e+07 | 3.5 | 0.49 | 0.00 | 0 | 0.012 (pk 18240) | none held at this N |
| 26 | bcr4bp | 3/1 | 3.07518 | RT-21 m1 b1, RT-21 m1 b3 | (c) | yes | 6e-13 | 21068 | 167437 | 2.44e+07 | 3.6 | 0.55 | 0.00 | 0 | 0.013 (pk 20927) | none held at this N |
| 27 | bcr4bp | 5/2 | 3.05583 | RT-31 m0 b1 | (c) | yes | 5e-12 | 13481 | 96281 | 4.88e+07 | 5.8 | 2.07 | 0.00 | 0 | 0.043 (pk 13871) | none held at this N |
| 28 | bcr4bp | 5/2 | 3.05583 | RT-31 m0 b0 | (c) | yes | 5e-13 | 14903 | 95217 | 2.61e+07 | 6.0 | 2.12 | 0.00 | 0 | 0.014 (pk 14715) | none held at this N |
| 29 | bcr4bp | 8/3 | 3.02516 | C21-3D#682 m0 b0, C21sp#438 m0 b0 | (c) | no | 3e-14 | 120325 | 221146 | 3.00e+01 | 48.2 | 0.42 | 0.21 | 2 | 0.036 (pk 122507) | none held at this N |
| 30 | bcr4bp | 8/3 | 3.02516 | C21-3D#682 m0 b1 | (c) x | no | 3e-14 | 120321 | 221127 | 3.05e+01 | 47.9 | 0.42 | 0.21 | 1 | 0.861 (pk 122507) | none held at this N |
| 31 | bcr4bp | 8/3 | 3.02472 | RT-31 m1 b1 | (c) | yes | 5e-13 | 6664 | 86231 | 1.85e+12 | 5.8 | 0.66 | 0.00 | 2 | 0.051 (pk 6815) | none held at this N |
| 32 | bcr4bp | 8/3 | 3.02472 | RT-31 m1 b0 | (c) | yes | 6e-13 | 6670 | 86238 | 1.85e+12 | 5.8 | 2.16 | 0.00 | 1 | 0.055 (pk 6824) | none held at this N |
| 33 | bcr4bp | 5/2 | 3.02010 | BR-R52S m1 b0 | (c) | no | 4e-13 | 76076 | 90031 | 2.18e+06 | 7.0 | 0.77 | 0.00 | 0 | 0.016 (pk 76582) | none held at this N |
| 34 | bcr4bp | 5/2 | 3.02010 | BR-R52S m1 b1 | (c) | no | 3e-11 | 82484 | 89650 | 1.97e+06 | 7.1 | 0.79 | 0.00 | 0 | 0.014 (pk 82172) | none held at this N |
| 35 | bcr4bp | 8/3 | 2.97116 | BR-R52S m2 b1 | (c) | yes | 6e-13 | 57209 | 76134 | 1.30e+09 | 7.8 | 0.65 | 0.00 | 2 | 0.055 (pk 57378) | none held at this N |
| 36 | bcr4bp | 8/3 | 2.97116 | BR-R52S m2 b0 | (c) | yes | 6e-13 | 57230 | 76148 | 1.30e+09 | 7.8 | 0.64 | 0.00 | 1 | 0.060 (pk 57449) | none held at this N |
| 37 | qbcp | 1/2 | 3.42072 | BR-R21S m1 b0 | (c) | no | 4e-11 | 194011 | 186879 | 1.09e+00 | 238.0 | n/a | 0.00 | 0 | 0.006 (pk 195680) | none held at this N |
| 38 | qbcp | 1/2 | 3.42072 | BR-R21S m1 b1 | (c) | no | 4e-11 | 196585 | 187815 | 1.00e+00 | stable | n/a | 0.00 | 0 | 0.006 (pk 195395) | none held at this N |
| 39 | qbcp | 2/1 | 3.19716 | BR-R52S m0 b1, BR-R52S m0 b3 | (c) | no | 6e-11 | 70421 | 125143 | 7.57e+00 | 20.2 | 2.94 | 0.00 | 0 | 0.017 (pk 70196) | none held at this N |
| 40 | qbcp | 2/1 | 3.19716 | BR-R52S m0 b2 | (c) | no | 1e-11 | 81425 | 133162 | 9.61e+00 | 18.1 | 2.90 | 0.00 | 0 | 0.024 (pk 80243) | none held at this N |
| 41 | qbcp | 2/1 | 3.19716 | BR-R52S m0 b0 | (c) | no | 1e-11 | 81578 | 133376 | 9.62e+00 | 18.1 | 2.90 | 0.00 | 0 | 0.026 (pk 80217) | none held at this N |
| 42 | qbcp | 5/2 | 3.18010 | BR-C32 m2 b1, RT-32 m1 b1 | (b) | yes | 6e-10 | 7348 | 118464 | 2.47e+03 | 13.1 | 1.55 | 0.00 | 0 | 0.020 (pk 7670) | 180A_1_t2: 2.2e-04 / 4.0e-04 |
| 43 | qbcp | 5/2 | 3.18010 | BR-C32 m2 b0, RT-32 m1 b0 | (c) | yes | 6e-10 | 9933 | 119146 | 3.95e+04 | 9.7 | 1.55 | 0.00 | 1 | 0.014 (pk 10048) | 032B_1_t4: 2.5e-02 / 3.6e-02 |
| 44 | qbcp | 5/2 | 3.15168 | BR-C32 m1 b1, RT-32 m2 b1 | (c) x | yes | 3e-10 | 17750 | 107805 | 4.58e+06 | 6.7 | 0.95 | 0.00 | 0 | 0.563 (pk 23199) | 180A_1_t2: 1.3e-02 / 7.9e-02 |
| 45 | qbcp | 5/2 | 3.15168 | BR-C32 m1 b0, RT-32 m2 b0 | (c) | yes | 7e-10 | 22239 | 109914 | 1.24e+07 | 6.3 | 1.66 | 0.00 | 1 | 0.015 (pk 21806) | 032B_1_t4: 1.2e-02 / 1.5e-01 |
| 46 | qbcp | 3/2 | 3.15107 | BR-C11a m1 b1, RT-11 m2 b1 | (c) | yes | 9e-11 | 22614 | 185717 | 3.51e+06 | 4.1 | 1.87 | 0.00 | 0 | 0.020 (pk 21080) | 146A_t3: 2.1e-01 / 4.7e-01 |
| 47 | qbcp | 3/2 | 3.15107 | BR-C11a m1 b0, RT-11 m2 b0 | (c) | yes | 3e-11 | 24390 | 183834 | 2.25e+07 | 3.6 | 1.86 | 0.00 | 0 | 0.038 (pk 25104) | 146A_t4: 8.6e-02 / 1.0e-01 |
| 48 | qbcp | 8/3 | 3.12875 | BR-C32 m0 b1, RT-32 m3 b1 | (c) x | yes | 1e-10 | 23996 | 99635 | 5.70e+15 | 4.5 | 0.73 | 0.00 | 0 | 0.141 (pk 24504) | none held at this N |
| 49 | qbcp | 8/3 | 3.12875 | BR-C32 m0 b0, RT-32 m3 b0 | (c) x | yes | 1e-09 | 23996 | 99540 | 5.68e+15 | 4.5 | 1.54 | 0.00 | 1 | 0.134 (pk 24496) | none held at this N |
| 50 | qbcp | 3/1 | 3.12523 | RT-21 m0 b2 | (c) | yes | 3e-11 | 25786 | 221094 | 7.43e+04 | 5.5 | 0.68 | 0.00 | 0 | 0.035 (pk 25109) | 146A_t3: 1.5e-01 / 5.4e-01 |
| 51 | qbcp | 3/1 | 3.12523 | RT-21 m0 b0 | (c) | yes | 3e-11 | 25791 | 221223 | 7.35e+04 | 5.5 | 0.68 | 0.00 | 0 | 0.034 (pk 25103) | 146A_t3: 1.5e-01 / 5.4e-01 |
| 52 | qbcp | 2/1 | 3.12107 | BR-C11a m2 b3, RT-11 m0 b3 | (c) | yes | 2e-11 | 25602 | 226554 | 4.12e+05 | 3.2 | 2.60 | 0.00 | 0 | 0.020 (pk 25005) | none held at this N |
| 53 | qbcp | 2/1 | 3.12107 | BR-C11a m2 b1, RT-11 m0 b1 | (c) | yes | 2e-11 | 25608 | 226767 | 4.11e+05 | 3.2 | 2.60 | 0.00 | 0 | 0.020 (pk 24999) | none held at this N |
| 54 | qbcp | 2/1 | 3.12107 | BR-C11a m2 b0, BR-C11a m2 b2, RT-11 m0 b0, RT-11 m0 b2 | (c) | yes | 2e-11 | 26065 | 224587 | 3.63e+05 | 3.2 | 2.58 | 0.00 | 0 | 0.026 (pk 25931) | none held at this N |
| 55 | qbcp | 3/2 | 3.09204 | BR-C11a m0 b1, RT-11 m3 b1 | (c) | yes | 3e-11 | 21171 | 137969 | 2.41e+10 | 2.6 | 1.76 | 0.00 | 0 | 0.016 (pk 20531) | 146A_t3: 2.1e-01 / 5.4e-01 |
| 56 | qbcp | 3/2 | 3.09204 | BR-C11a m0 b0, RT-11 m3 b0 | (c) | yes | 2e-11 | 23142 | 145936 | 1.59e+10 | 2.6 | 1.74 | 0.00 | 0 | 0.008 (pk 23396) | 146A_t4: 6.3e-02 / 3.1e-01 |
| 57 | qbcp | 5/2 | 3.07692 | RT-11 m1 b1 | (c) | yes | 9e-11 | 19668 | 215075 | 6.97e+13 | 3.2 | 1.79 | 0.00 | 0 | 0.015 (pk 18514) | 053d_2_t4: 5.5e-02 / 3.7e-01 |
| 58 | qbcp | 5/2 | 3.07692 | RT-11 m1 b0 | (c) | yes | 1e-10 | 20323 | 229232 | 6.45e+12 | 3.5 | 1.77 | 0.00 | 0 | 0.063 (pk 21010) | 180A_1_t2: 8.3e-02 / 4.3e-01 |
| 59 | qbcp | 3/1 | 3.07518 | RT-21 m1 b2 | (c) | yes | 2e-11 | 18240 | 158638 | 4.49e+07 | 3.5 | 0.50 | 0.00 | 0 | 0.012 (pk 17701) | 146A_t3: 7.4e-02 / 4.5e-01 |
| 60 | qbcp | 3/1 | 3.07518 | RT-21 m1 b0 | (c) | yes | 2e-11 | 18249 | 158819 | 4.46e+07 | 3.5 | 0.50 | 0.00 | 0 | 0.013 (pk 17699) | 146A_t3: 7.4e-02 / 4.5e-01 |
| 61 | qbcp | 3/1 | 3.07518 | RT-21 m1 b1, RT-21 m1 b3 | (c) | yes | 2e-11 | 20927 | 168070 | 2.37e+07 | 3.6 | 0.68 | 0.00 | 0 | 0.013 (pk 21068) | 146A_t3: 8.5e-02 / 3.7e-01 |
| 62 | qbcp | 5/2 | 3.05583 | RT-31 m0 b1 | (c) | yes | 1e-09 | 13871 | 95016 | 4.56e+07 | 5.8 | 2.04 | 0.00 | 0 | 0.043 (pk 13481) | 357_t3: 1.2e-01 / 3.6e-01 |
| 63 | qbcp | 5/2 | 3.05583 | RT-31 m0 b0 | (b) | yes | 1e-09 | 14715 | 96783 | 2.89e+07 | 6.0 | 2.03 | 0.00 | 0 | 0.014 (pk 14903) | 357_t3: 4.3e-04 / 1.1e-03 |
| 64 | qbcp | 8/3 | 3.02516 | C21sp#438 m0 b1 | (c) | no | 9e-11 | 122507 | 221046 | 3.51e+01 | 46.0 | 0.42 | 0.21 | 2 | 0.036 (pk 120325) | none held at this N |
| 65 | qbcp | 8/3 | 3.02472 | RT-31 m1 b0 | (c) | yes | 2e-09 | 6815 | 87733 | 1.94e+12 | 5.8 | 1.96 | 0.00 | 1 | 0.051 (pk 6664) | none held at this N |
| 66 | qbcp | 8/3 | 3.02472 | RT-31 m1 b1 | (c) | yes | 1e-10 | 6824 | 87848 | 1.95e+12 | 5.8 | 0.59 | 0.00 | 0 | 0.055 (pk 6670) | none held at this N |
| 67 | qbcp | 5/2 | 3.02010 | BR-R52S m1 b0 | (c) | no | 4e-11 | 76582 | 88915 | 2.18e+06 | 7.0 | 0.78 | 0.00 | 0 | 0.016 (pk 76076) | 053d_2_t4: 7.2e-02 / 4.7e-01 |
| 68 | qbcp | 5/2 | 3.02010 | BR-R52S m1 b1 | (c) | no | 3e-11 | 82172 | 91068 | 1.94e+06 | 7.1 | 0.78 | 0.00 | 0 | 0.014 (pk 82484) | 357_t3: 2.1e-02 / 4.5e-02 |
| 69 | qbcp | 5/2 | 3.01251 | C21-3D#682 m1 b1, C21-3D#682 m1 b3, C21sp#438 m1 b1, C21sp#438 m1 b3 | (c) x | yes | 4e-11 | 58271 | 177332 | 7.58e+01 | 23.6 | 0.71 | 0.26 | 1 | 1.392 (pk 7670) | 180A_1_t1: 1.8e-02 / 1.3e-01 |
| 70 | qbcp | 5/2 | 3.01251 | C21-3D#682 m1 b0, C21-3D#682 m1 b2, C21sp#438 m1 b0, C21sp#438 m1 b2 | (c) x | yes | 4e-11 | 65938 | 188177 | 6.16e+01 | 24.8 | 1.00 | 0.26 | 2 | 1.136 (pk 14903) | 084_2_t3: 3.3e-02 / 8.9e-02 |
| 71 | qbcp | 8/3 | 2.97116 | BR-R52S m2 b0 | (c) | yes | 6e-11 | 57378 | 77566 | 1.39e+09 | 7.8 | 0.52 | 0.00 | 2 | 0.055 (pk 57209) | none held at this N |
| 72 | qbcp | 8/3 | 2.97116 | BR-R52S m2 b1 | (c) | yes | 6e-11 | 57449 | 77673 | 1.40e+09 | 7.8 | 0.52 | 0.00 | 1 | 0.060 (pk 57230) | none held at this N |

## 8. Recommendation

- **No catalogue row from this run yet.**
- The two class (b) orbits (the C32 180A_1 and C31 357 closures in the QBCP) are the only results
  that both survive an independent check and add something to a published record. Each completes
  a numerical failure that Leiva & Briozzo reported (single shooting over 5 Tsun did not converge,
  their p. 239). The completion is modest: the paper already has the arcs, and our orbits sit
  85 to 165 km from them with matching multipliers. If the owner wants rows, these are the
  candidates, labelled as "Sun-forced (QBCP) closure of Leiva & Briozzo 2008 arc 180A_1 / 357",
  not as discoveries.
- Gates still needed for those two, and for any class (c) orbit:
  1. Re-correct the QBCP orbits at P = N Tg exactly (one Newton step) and record the closure.
  2. A literature check that knows the Sun-forced literature. Add Leiva & Briozzo 2005/2008,
     Oshima 2022, Boudad-Howell-Davis 2020 and Brown et al. 2025 as anchors, with a forcing-model
     field, then rerun.
  3. Continue to a real ephemeris (DE440, eccentric and inclined Moon) with the Sun phase set to a
     real epoch. Report the closure or the maintenance cost against Genova & Aldrin's
     39 m/s per month.
  4. For any (c) orbit with an "x" or a branch point: show that the eps path did not switch
     branch, and confirm the parent is on the catalogued family across the after-bifurcation
     crossings.
- Correct OUTSTANDING #905: 180A_2 was not completed (section 1). The C32 5/2 bicircular values
  are a reproduction of the #884 reviewer's table, not a new result. The run summary counts
  duplicated parents three times over and hides the QBCP `sun_phase_returns = False`.
- Investigate the bicircular C21 3D 5/2 stall at eps of about 1e-4, and the QBCP C21 3D 8/3
  chattering. These are numerical failures, not negative results, and must not enter the
  negative-results registry as "does not persist".
