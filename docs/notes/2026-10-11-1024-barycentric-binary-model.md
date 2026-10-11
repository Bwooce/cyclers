# #1024: a barycentric circular-model option for binaries in the two-working-body generator

Status: DONE. Pre-registration secs. 1-7 (with amendments 3a and 7a), committed before the code
and runs; results sec. 8. No catalogue edits. The #998 adjudication (#1036) is a later owner step. #949 (the
binary-star cases) is covered by the #1024 entry but is not run here.

## 1. The defect being fixed

`scripts/run_998_enumerate.py::pluto_system` builds Pluto-Charon as follows (#998 note sec. 1):
- every body is on a circle about the ORIGIN, radius from its period by Kepler III with the system GM
  975.5 km^3/s^2;
- the spacecraft legs are two-body about the origin with that GM.

That puts Charon at 19,596 km, the radius of the RELATIVE (Pluto-Charon) orbit, moving at
sqrt(GM/a) = 0.2231 km/s. In the plu060 kernel Charon is at a mean 17,464 km from the barycentre,
moving at about 0.199 km/s. The mass fraction mu_C = GM_Charon / GM_system = 106.1 / 975.5 = 0.10877,
and 1 - mu_C = 0.8912: the ideal model's Charon is 12 % too far out and 12 % too fast for a body on the
barycentric circle. The small moons are already consistent: they circle the barycentre, and their
radii come from Kepler III with the system GM.

## 2. The model option (exactly)

Opt-in, default OFF. Every existing cell and run is bitwise unchanged with the option off.

With the option ON, for a binary made of a primary (P) and one massive secondary (S, the flyby body),
mass fraction mu_S = GM_S / (GM_P + GM_S):
- **Origin:** the barycentre, inertial.
- **Spacecraft legs:** two-body about the barycentre with the TOTAL GM (the cell's existing `mu`,
  975.5 for Pluto-Charon). Lambert legs, full-rev and half-rev legs, and the re-fly all use it.
- **Secondary (Charon):** a circle about the barycentre of radius (1 - mu_S) a_rel, where a_rel is the
  existing Kepler-III radius (19,596 km), giving 17,464.5 km. Same period (the fitted 6.38722 d),
  same phase convention (angle 0 at t = 0). Its speed is the KINEMATIC 2 pi r / P = 0.1988 km/s, not
  sqrt(GM/r). Charon is not a test particle of the barycentre monopole.
- **Primary (Pluto):** on the opposite circle, radius mu_S a_rel = 2,131.5 km, same period. It is used
  only for the surface screen (sec. 4), as a distance from the MOVING Pluto. It is not a flyby body and
  not a gravitating point in the legs.
- **Unchanged:** the target moons (Styx, Nix, Kerberos, Hydra), with their radii, periods and
  massless status. Also unchanged: Charon's GM, radius, floor and tidal band, the gate and its floors,
  the near-180 rule, the seeds and the pass definitions.

Why the TOTAL GM at the barycentre, and not the primary's GM as the #1024 entry's wording suggests:
- This makes the option exactly the circular-coplanar limit of the model in which #998's
  real-ephemeris check already ran. `scripts/realeph_998.py` uses the system GM 975.5 at the
  barycentre with the kernel's barycentric body states.
- It is also the far-field (monopole) field of the binary, which is what the small moons orbit.
- With the primary's GM alone, the small-moon circles (system-GM Kepler radii) would be inconsistent
  with the legs.

What the option does NOT fix:
- Inside the binary (r comparable to a_rel) the field is not a point mass at the barycentre. Every
  #998 strong cycler has r_min of 1,256-16,067 km, inside Charon's 17,464-km circle, so every one
  spends part of its cycle where the monopole is wrong. This is a limit of any one-centre patched
  conic.
- It is reported as a diagnostic, r_min / a_S and the cycle fraction inside 2 a_S, not a gate.
- The real ephemeris, and eventually an n-body or CR3BP check, stays the arbiter. The option corrects
  Charon's STATE at the flyby, which drives the V_inf, not the field along the legs.

## 3. Code (where and what)

Code audit done before writing this note. Two places in `src/cyclerfinder/search/two_working_body.py`
assume each body's circle is Keplerian about `mu`:
- `resonant_circle`, the 1:1 shortcut `|v_sc| = |V_body|`. For a non-Keplerian body the general
  formula (spacecraft a from the period, `v^2 = mu (2/r - 1/a_sc)`) must be used. With the barycentric
  Charon the shortcut would give 0.199 km/s where the correct return speed is 0.249 km/s.
- `half_rev_vectors`, which assumes the body speed equals sqrt(mu/r) (p = r, transverse speed = W). It
  must REFUSE a non-Keplerian body (a clear error), not return a wrong conic.

Also in that file, `wrap_rotation` (uniform rotation about the barycentre) still holds.
`sphere_of_influence_km`: an uncommitted #1060 edit by another agent reads `a` from `system.bodies`.
If it lands, Charon's SOI in the barycentric cells uses 17,464 km (SOI 7,192 km) instead of 19,596 km
(8,070 km). That is a diagnostic and filter input, reported both ways.

Plan:
- New module `src/cyclerfinder/search/two_working_body_binary.py`:
  - `BarycentricCircularSystem`: a `CircularSystem` with a kinematic-speed secondary and the moving
    primary;
  - `barycentric_binary(system, secondary, mass_fraction)` to build one;
  - `primary_min_distance_km(...)`: the closest approach to the moving primary along every leg,
    full-rev legs included.
- Two small guarded edits in `two_working_body.py`: the `resonant_circle` branch, and the
  `half_rev_vectors` refusal. Both are keyed on an attribute only the new class sets, so the default
  path is bitwise unchanged. They are made only after the #1060 edit has landed or with its owner's
  agreement.
- Tests in `tests/search/` at production scales (km, km^3/s^2):
  - mass fraction 0 reproduces the default states to 1e-12 relative;
  - the Pluto-Charon barycentric radius is within 1 km of the plu060 mean 17,464 km, and the speed
    is 2 pi r / P;
  - the primary is opposite the secondary, at mu_S a_rel, with zero total momentum;
  - a full-rev leg in the barycentric model re-flies (kepler_step, and DOP853) back onto Charon to
    < 1 km;
  - `half_rev_vectors` raises for the barycentric Charon;
  - `primary_min_distance_km` with mass fraction 0 equals `leg_extent`'s r_min.
- Driver `scripts/run_1024_barycentric.py` (declares #1024).

## 3a. Amendment A1 (lead, 2026-10-11), before any edit to `two_working_body.py`

The two edits are restated exactly, now that the lead has confirmed the file is pluto-smallmoons-sonnet's
#1060 edit in progress. They are made only after that commit lands.

1. `resonant_circle`: the 1:1 shortcut `vf2 = |V_body|^2` is used only when the body is Keplerian about
   `mu` (`getattr(system, "is_keplerian", None)` is absent, or returns True). Otherwise the general
   branch runs: `a_sc` from the period, `vf2 = mu (2/r - 1/a_sc)`.
2. `half_rev_vectors`: raises `ValueError` ("half-rev legs need a Keplerian body circle (#1024)") for a
   non-Keplerian body. The continuation EXCLUDES half-rev structures before solving. That is 48 of the
   457 #998 gate-passers; the 17 of the 224 strong are:
  - ps: `k4|HCharon/3,1,p|RCharon/1:2|LCharon>Styx/0s|LStyx>Charon/1l`
  - ps: `k4|HCharon/3,1,p|LCharon>Styx/0s|LStyx>Charon/1l`
  - ps: `k4|HCharon/3,1,p|RCharon/2:3|LCharon>Styx/0s|LStyx>Charon/0s`
  - ps: `k5|HCharon/3,1,p|HCharon/3,1,p|LCharon>Styx/1l|LStyx>Charon/1h`
  - ps: `k5|HCharon/3,1,p|HCharon/3,1,p|LCharon>Styx/0s|LStyx>Charon/1l`
  - ps: `k6|RCharon/3:2|HCharon/3,1,p|LCharon>Styx/1h|LStyx>Charon/1l`
  - ps: `k6|RCharon/2:1|HCharon/3,1,p|LCharon>Styx/1h|LStyx>Charon/1l`
  - ps: `k6|HCharon/3,1,p|LCharon>Styx/1h|LStyx>Charon/1l`
  - ps: `k6|HCharon/3,1,p|LCharon>Styx/1h|LStyx>Charon/0s`
  - pn: `k4|HCharon/3,1,p|HCharon/3,1,p|LCharon>Nix/0s|LNix>Charon/0s`
  - pn: `k5|HCharon/3,1,p|RCharon/2:1|LCharon>Nix/1l|LNix>Charon/0s`
  - pn: `k5|HCharon/3,1,p|LCharon>Nix/1l|LNix>Charon/1h`
  - pn: `k5|HCharon/3,1,p|LCharon>Nix/0s|LNix>Charon/1l`
  - pk: `k3|HCharon/3,1,p|LCharon>Kerberos/0s|LKerberos>Charon/0s`
  - pk: `k5|HCharon/3,1,p|HCharon/3,1,p|LCharon>Kerberos/0s|LKerberos>Charon/0s`
  - ph: `k3|HCharon/3,1,p|LCharon>Hydra/0s|LHydra>Charon/0s`
  - ph: `k4|HCharon/3,1,p|LCharon>Hydra/0s|LHydra>Charon/0s`

3. Charon SOI diagnostic: as in sec. 3. If #1060's `sphere_of_influence_km` reads `a` from
   `system.bodies`, then in the barycentric cells Charon's SOI uses 17,464.5 km (7,192 km) instead of
   19,596 km (8,070 km). Both are reported. The "strong" count is given under each.

REGRESSION TEST (new, pinned):
- For the existing Keplerian systems, `resonant_circle` and `half_rev_vectors` must return
  BIT-IDENTICAL results before and after the edits. The cells are ev, gc, vm and the #998 `pn` cell,
  each over a grid of epochs, V_inf values and the catalogue's resonant and half-rev legs.
- Expected values are recorded from the pre-edit code at its current commit into a test fixture, then
  asserted with `==` (not a tolerance) against the edited code.
- So the #973 and #998 numbers cannot move without the test failing.

## 4. Method: continue each #998 member, do not re-enumerate

For each of the 457 #998 gate-passing cyclers (representative zero from
`data/998_pluto_smallmoons/<cell>_gauntlet.json`):
- **Ramp:** lambda from 0 to 1, with Charon's radius (1 - lambda mu_C) a_rel and the kinematic speed
  at every lambda; the primary circle follows. Steps of 0.1, halved on failure down to 1/640. At each
  step, `correct_dates` (tol 1e-8 km/s) from the previous solution.
- **lambda = 0 control (#998 off path):** the residual at the stored x_days must be < 1e-8 km/s, and
  the assessed V_inf must equal the stored ones to 1e-9 km/s.
- **Branch-jump guard:** a step is rejected (and halved) if any date moves by more than 0.5 d, or any
  node V_inf by more than 0.02 km/s, per 0.1 of lambda (scaled to the step).
- **Outcomes:**
  - CONTINUES: reaches lambda = 1;
  - FOLDS / STOPS: the smallest step fails. A fold is NOT a negative, and is reported with the last
    lambda reached;
  - EXCLUDED: the structure has a half-rev leg (48 of the 457; 17 of the 224 strong). The half-rev
    formula does not apply to a non-Keplerian Charon. Reported, not run.
- **At lambda = 1, re-assessed exactly as #998** (`assess`; `gauntlet_998` logic):
  - exact zero;
  - #888/#937 gate (tri-state with Charon's tidal band);
  - near-180 rule;
  - Kepler-step re-propagation < 1 km;
  - DOP853 re-fly;
  - the #998 strong filters: the Pluto surface as the closest approach to the MOVING Pluto
    > 1,188.3 km, and the Charon periapsis inside its SOI. The SOI is taken as `sphere_of_influence_km`
    gives it at run time, and the 8,070-km value is also reported.
- **Diagnostics:** r_min / a_S, and the cycle fraction inside 2 a_S.

A full barycentric re-enumeration (about 27,000 structures) is offered to the lead as an option, not
run here.

## 5. Controls

1. **Off path, bitwise.** gc-1 and GanCal#5 (the #943 gc gauntlet) and VenMar#45 (the #942 vm
   gauntlet / `control_generator_venmar45.json`) give their stored V_inf to 1e-12 km/s after the code
   change. Plus the lambda = 0 check of sec. 4 on every #998 member.
2. **On path, low mass fraction: the option must reproduce the current results.**
   - VenMar#45 with Venus as the secondary: mu_V = 324,860 / (1.3271244e11 + 324,860) = 2.448e-6;
     v_V = 35.0 km/s; mu v = 8.6e-5 km/s. TOLERANCE, derived: |Delta V_inf| <= 10 mu v_S = 8.6e-4
     km/s against the option-off value, and the published 8.22 / 12.96 within 0.05 (R-S 2007 Table 3
     rounding). Gate status unchanged.
   - gc-1 with Ganymede as the secondary (the gc cell treated as the Jupiter-Ganymede binary; Callisto
     unchanged): mu_G = 9,887.83 / (126,686,535 + 9,887.83) = 7.80e-5; v_G = 10.88 km/s; 10 mu v = 8.5e-3
     km/s. TOLERANCE: |Delta V_inf| <= 8.5e-3 km/s at both moons, gate status unchanged.
   - The Russell-Strange cells use the primary's GM (Jupiter alone, the Sun alone) as `mu`. Reading it
     as the system total is an O(mu) inconsistency, below these tolerances.
   - A failed control stops the run; nothing is loosened.
3. **The real-world test (falsifiable).** Charon's barycentric state matching the kernel proves
   nothing: the option is built to match it. The test is #998 sec. 6 repeated for the 12 best strong
   cyclers. Each member is seeded from its CONTINUED barycentric zero (if it continues), and the defect
   on the ideal dates on plu060 states (before the refit) is compared with #998's off-model defects
   (4.3-32.0 m/s).
   - PREDICTION: the defect falls for at least 9 of the 12, and the median falls by at least 30 %.
   - If not: the 12 % Charon-state error is not what dominates the ideal-to-kernel defect (the small
     moons' eccentricities and inclinations, or the epoch alignment, would be the candidates). The
     option is then not shown to improve the ideal model as a seed generator. It stays opt-in, and
     #998 stays seeded from the default model.
   - Same script logic (`scripts/realeph_998.py`: frame, conjunction epochs, residual), with the
     barycentric x_days. One cycle, five epochs.

## 6. Expected effect on the 224 strong candidates (written before any run)

What moves, by construction:
- Charon's radius falls by 10.9 % (19,596 -> 17,464.5 km), and its speed by 10.9 % (0.2231 -> 0.1988
  km/s, Delta 0.024 km/s).
- Hence the V_inf at Charon, which is 0.16-0.36 km/s today. Expected shifts are of order 0.02-0.06 km/s
  (10-35 %), from the 0.024 km/s speed change plus the geometry change of the Lambert legs. The sign
  varies with the encounter geometry.
- Hence the demanded turns, the available bend (strongly V_inf-dependent at these speeds: about 117 deg
  at 0.16 km/s and 65 deg at 0.36 km/s at the 100-km floor), the ratios (tens of percent either way),
  the required periapsis altitudes, the Lambert-leg dates (hours), r_min and the distance to the
  moving Pluto (up to +-2,132 km against the origin-based r_min).
- The Charon-SOI filter tightens if the #1060 SOI change lands (8,070 -> 7,192 km).

What does not move: the periods, k, the structures, the small-moon circles, the central GM, Charon's GM,
radius and floor, and the gate.

Not predicted: how many of the 224 stay gate-passing or strong. The direction is mixed (V_inf can rise
or fall). The candidates with r_min near Pluto (1,256-4,000 km) are the ones most at risk from the
moving-Pluto screen.

## 7. Runtime and operations

- A 10-member pilot is timed under the current load before any counted run. The projected time per moon
  cell decides the run. If every call fits in under about 8 minutes (shards of members, `--members i/n`)
  I run it; otherwise the lead gets the exact command.
- No backgrounding. Commits by pathspec. Hashes are reported from `git log -1` in a step after the
  commit.

## 7a. Amendment A2 (implementation bugs, found in the first continuation run before any analysis; 2026-10-11 16:00 AEDT)

The first run of `continue` (before commit; its outputs are discarded) showed two implementation bugs.
The criteria are unchanged:
1. **Records keyed by the cycle key.** The #998 gauntlet's representative keys are not unique (two
   merged groups can share a key), so members were lost (ps 183 -> 181, pk 62 -> 56). Records are now
   keyed by the gauntlet index. The `realeph` step matches a #998 member by key AND its ideal worst
   ratio.
2. **Floating-point lambda.** Ten steps of 0.1 sum to 0.9999999999999999. The next step was then
   about 1e-16, and the branch-jump guard, "scaled to the step", allowed only about 1e-16 of motion.
   Any solver noise was rejected, and 81 of 107 "folds" ended at lambda = 1 - 1e-16. The fix:
   - lambda snaps to 1 within 1e-9;
   - the guard's scale is floored at the smallest legal step, 1/640 of lambda. That floor is inside
     the pre-registered "halved down to 1/640", so the guard is unchanged for every legal step.

The continuation is re-run in full after this amendment. The controls are unaffected: they reached
lambda = 1 through the same code, and are re-run anyway.

## 8. Results (2026-10-11, 16:00-16:04 AEDT; code at the amendment-A2 commit)

Data: `data/1024_barycentric/`:
- `controls.json`;
- `<cell>_continue.json`, one record per #998 gate-passer, keyed by gauntlet index;
- `realeph_defect.json`;
- a static log beside each.
Wall time about 4 min in all (1 s per member).

### 8.1 Code and controls

Code:
- `two_working_body.py` edits: commit 80b3420c.
- The bit-identity regression on the Keplerian ev, gc, vm and pn cells (459 rows) passes after the
  edits.
- A full-rev n:m leg from the barycentric Charon re-flies onto Charon to < 1 km (Kepler step and
  DOP853, n:m = 1:1, 2:1, 1:2, 3:2). Positive control of that test: the old 1:1 shortcut, forced on,
  misses by 19,270 km.
- `half_rev_vectors` refuses the barycentric Charon.
- `mypy src tests` clean.

Controls (sec. 5), all OK:
- **Off path, bitwise:** VenMar#45, gc-1 and GanCal#5 give their stored V_inf exactly (|dV_inf| = 0).
  At lambda = 0, every one of the 409 non-excluded #998 members has a residual <= 1.6e-14 km/s and a
  Charon V_inf equal to the stored value to <= 3.1e-16 km/s.
- **On path:**
  - VenMar#45 with Venus as the secondary (mass fraction 2.45e-6): |dV_inf| 1.5e-4 km/s (tolerance
    8.6e-4), published 8.22 / 12.96 within 0.005, gate pass, dates moved by up to 22 s.
  - gc-1 with Ganymede as the secondary (7.80e-5): |dV_inf| 3.6e-4 km/s (tolerance 8.5e-3), gate pass,
    dates moved by up to 37 s.

### 8.2 The 457 #998 gate-passers under the option

| Cell | Gate-passing | Continues to mass fraction 0.1088 | Folds / stops | Excluded (half-rev) | Gate pass at the end | Strong at the end | of which strong in #998 | newly strong |
|---|---|---|---|---|---|---|---|---|
| Styx (ps) | 183 | 140 | 16 | 27 | 139 | 55 | 54 of 108 | 1 |
| Nix (pn) | 149 | 127 | 7 | 15 | 116 | 59 | 41 of 46 | 18 |
| Kerberos (pk) | 62 | 59 | 1 | 2 | 55 | 36 | 36 of 41 | 0 |
| Hydra (ph) | 63 | 57 | 2 | 4 | 52 | 25 | 24 of 29 | 1 |
| all | 457 | 383 | 26 | 48 | 362 | 175 | 155 of 224 | 20 |

- "Strong" is the #998 definition with the Pluto surface measured from the MOVING Pluto. The Charon
  SOI is the run-time value (7,192 km after #1060). With the 8,070-km SOI the strong count is 176:
  one Nix member is the difference.
- **The 224 #998 strong members:**
  - 204 continue, and 203 of them still pass the gate (one is "indeterminate");
  - 3 fold or stop (all Styx; a fold is not a negative);
  - 17 are excluded (half-rev, listed in sec. 3a).
  - Of the 204, 155 stay strong. 48 lose the strong label on the Pluto-surface screen: once Pluto
    moves (up to +-2,131 km against the origin-based r_min), their closest approach to it is below
    1,188.3 km.
- **The 20 newly strong members** were outside #998's strong set on the old Pluto-surface screen (18:
  Nix 17, Hydra 1) or on the SOI screen (2).
- **What moved, for the 204 continued strong members:**
  - Charon V_inf rose by +0.003 to +0.023 km/s (median +0.011).
  - The worst gate ratio changed by -5 % to +18 % (median +11 %).
  - Dates moved by hours.
  - The pre-registered expectation (sec. 6) was V_inf shifts of 0.02-0.06 km/s of mixed sign and ratio
    changes of "tens of percent either way". The observed shifts are smaller and all positive. The
    expectation overstated them.
- **Diagnostic (sec. 2, "does not fix"):** the continued strong members reach r_min / a_C = 0.06-0.83
  (median 0.23), and spend 5-50 % of their Lambert-leg time (median 16 %) inside 2 a_C. Every one of
  them flies deep inside the binary, where the barycentre monopole is not the field. The Pluto-surface
  screen is evaluated in that region, so the 48 "lost" and 20 "gained" labels are no more reliable
  than the monopole there.

### 8.3 The real-world test (sec. 5, control 3): PASSED as predicted

The plu060 defect on the ideal dates, before the refit, for the 12 best strong cyclers (#998 sec. 6),
now seeded from their continued barycentric zeros:
- 11 of the 12 continued. The defect fell for all 11.
- The median of the per-cycler maxima over 5 epochs fell from 7.9 to 2.6 m/s (-67 %).
- Prediction (written before the run): at least 9 of 12 fall and the median falls by at least 30 %.
  MET.
- Per cycler, #998 -> barycentric (max over epochs, m/s):

| Cycler | #998 | Barycentric |
|---|---|---|
| Styx k3 | 7.9 | 3.0 |
| Styx k5 | 4.3 | 1.2 |
| Nix k4 (generic) | 31.1 | 2.4 |
| Nix k4 (1:2 pair) | 12.3 | 3.4 |
| Nix k3 | 10.2 | 2.6 |
| Kerberos k6 (1:2) | 5.4 | 1.6 |
| Kerberos k6 (2:1 pair) | 14.0 | 1.6 |
| Kerberos k6 (2:1, 3:2) | 7.3 | 2.6 |
| Hydra k4 | 6.3 | 3.1 |
| Hydra k6 | 5.2 | 3.1 |
| Hydra k5 | 32.0 | 24.9 |
| Styx k4 | 11.8 | did not continue (stops at lambda 0.978) |

  The Hydra k5 outlier is one epoch (32.0 -> 24.9); its other epochs are 0.5-3.0 m/s.

Reading: the 12 % Charon-state error was the dominant part of #998's ideal-to-kernel defect. With the
option on, the ideal model is a better seed generator for the plu060 step.

### 8.4 A finding about #998 sec. 6.1 (not fixed here)

`scripts/realeph_998.py`'s `EphSystem` uses the kernel Charon states (barycentric, about 0.199 km/s),
but it has no `is_keplerian`. So `resonant_circle` still takes the 1:1 shortcut |v_sc| = |V_Charon|
there. For a spacecraft about the 975.5 monopole that is not the speed of a 1:1 return; the return
needs about 0.249 km/s at Charon's barycentric radius (sec. 3). The test of 8.1 shows the shortcut
misses by about 19,000 km in the barycentric circle.
- 135 of the 457 #998 gate-passers (90 of the 224 strong) have an RCharon/1:1 leg. Their #998 sec. 6.1
  kernel "refit to under 1 m/s" assumed a full-rev return that the conic does not make.
- The 12 best of sec. 6 have no 1:1 leg and are unaffected.
- Fix: `EphSystem.is_keplerian` returns False for Charon. That is a #998 script change, owned by its
  task, then a re-run of sec. 6.1 for those 135. Proposed as a task for the lead to register.

### 8.5 What this means for #998 / #1036

- With the option on, 155 of the 224 strong #998 cyclers stay strong (204 continue). 20 others become
  strong. 48 lose the label to the moving-Pluto screen. 17 half-rev members cannot be run in this
  model. 3 do not continue.
- The option brings the ideal seeds closer to plu060 (sec. 8.3). It does not fix the field inside the
  binary, where every one of them flies.
- #1036's adjudication should use the barycentric continued members (`data/1024_barycentric/`) as the
  ideal-model record, with the 8.4 caveat for the 1:1 members until sec. 6.1 is re-run.
- No catalogue edit; nothing is called novel.
- A full barycentric re-enumeration (about 27,000 structures; it would find members that do not
  continue from the old model) is offered as a lead launch.

## 9. Pre-registration: #1065, full re-enumeration in the barycentric model (2026-10-11 16:20 AEDT, before any counted run)

Sec. 4 continued the #998 members. It cannot find a barycentric member whose old-model parent did not
exist or did not pass. This run enumerates the four #998 cells from the start with the option on.
Driver: `scripts/run_1065_barycentric_enumerate.py` (subcommands enumerate, check, gauntlet, controls).
The lead launches it after window 3's suite.

### 9.1 Cells

| cell | system | bodies | structures |
|---|---|---|---|
| psb | barycentric_binary(pluto_system("Styx"), "Charon", 106.1/975.5) | Charon, Styx | 3,462 |
| pnb | same, Nix | Charon, Nix | 3,346 |
| pkb | same, Kerberos | Charon, Kerberos | 3,146 |
| phb | same, Hydra | Charon, Hydra | 3,025 |
| total | | | 12,979 |

- The systems are those of secs. 2 and 3: legs about the 975.5 monopole at the barycentre; Charon on
  its 17,464.5 km circle at the kinematic speed; Pluto on the opposite 2,131.5 km circle, used only by the
  surface screen. The small moons are on their #998 circles.
- Flags are the #998 ones: k 1-6, max-returns 2,0 (the moon is massless), transfer-revs 0,1,
  generic-revs 1,2, n-phase 24, n-split 10, n-refine 60.
- One change: Charon's catalogue has no half-rev returns (`half_revs=()`). `half_rev_vectors` refuses
  a non-Keplerian body (amendment A1). So the run has 12,979 structures, not the 26,831 of #998. The
  48 half-rev gate-passers of #998 (17 of them strong) are outside this run.
- Counts are from `enumerate --count-only` at the driver commit. All 383 continued keys of sec. 8
  (psb 140, pnb 127, pkb 59, phb 57) are in the structure sets (checked with `structures()` before
  this commit).

### 9.2 Pass definition (fixed before the run)

Per merged physical cycler (`scripts/analyse_942_enumeration.py` grouping), on its best member:
1. Exact zero: date residual below 1e-8 km/s (the enumeration's acceptance).
2. Gate: status "pass" under the #888/#937 demanded-turn gate (near-180 rule included), and the
   independent DOP853 re-fly misses every encounter by less than 1 km (`is_pass`).
3. #973 amendment 2.6/2.7 unscheduled-pass screen (`closest_unscheduled`, with the barycentric
   states): an interior distance minimum inside a body's radius on a fixed-geometry leg is an impact
   (reject); one inside the SOI is model-invalid (reject); on a free-direction full-rev leg it is a
   "direction-dependent" flag and passes.
4. Pluto surface: the minimum distance to the moving Pluto (`primary_min_distance_km`) above
   1,188.3 km.
5. Charon SOI: the largest demanded Charon periapsis inside the SOI of the run (7,192 km, #1060).
   Reported both ways: `strong` with 7,192 km and `strong_8070` with the #998 value of 8,070 km.
- "strong_1065" = 1-5 all true. The direction-dependent flag is carried in the record.
- Deep-dip caveat, a flag that does not gate: `r_min_over_a_charon` and the fraction of Lambert-leg
  time inside 2 a_C. The patched conic about the monopole is weakest inside the binary (sec. 8.5).

### 9.3 Coverage

`check` on each cell's shard directories: every structure done once without error, no duplicated
zero line, settings the same. A cell that fails coverage is re-run before any analysis.

### 9.4 Controls

- Positive: each of the 383 sec. 8 continued members counts as re-found when a zero in the same cell
  has the same key and a Charon V_inf within 1e-4 km/s of the continued value. "Same verdict" = the
  zero's `is_pass` equals the continued member's `pass`. Expected: at least 370 of 383 re-found (the
  phase grid is the #998 grid, and the zeros moved by at most 0.023 km/s in V_inf), and every re-found
  member with the same verdict. FAIL if fewer than 370 are re-found or any re-found member changes
  verdict. A miss is then examined one by one (seed grid against the moved zero) before the run's
  members are used.
- Negative: no zero in the run may sit at a #998 default-model solution: same key and Charon V_inf
  within 1e-4 km/s of the #998 value, over all 457 sec. 8 records. For continued members the
  V_inf moved 0.003-0.023 km/s, 30 times the tolerance. FAIL if any continued member's old zero is
  present (the option was not on). A hit on a non-continued record is reported and examined.
- `controls` writes both lists per cell.

### 9.5 Expected result (written before the run)

- Gate-passing physical cyclers: at least the re-found continued members (370-383), plus members
  that have no old-model parent. No number is predicted for the new members.
- Strong: of order 175 (sec. 8.2, 7,192 km SOI), plus new members, before the unscheduled-pass screen;
  that screen may remove some. No number is predicted for it (it was not run on #998).
- Nothing is called novel. No catalogue edit. #1036's adjudication takes this as the ideal-model
  record for the non-half-rev members.

### 9.6 Runtime and shards

- Pilot (16:10-16:12 AEDT, load average 22-32, 20 random structures, seed 1065): 2.8 s per structure
  wall in pnb, 3.4 s in psb (0.7-1.1 s CPU). #998 at its load: 0.6-1.9 s per structure by cell.
- Estimate: 12,979 structures at about 3.1 s = 11.2 h serial, about 5.6 h on 2 workers at the
  current load; about 2.3 h on 2 workers at the #998 rates.
- Shards: 32 per cell, 128 in all, about 100 structures each, about 5-6 min each at the loaded rate
  (under the 10-min tool limit of the coordination rule). A shard re-run resumes from its
  `structures.jsonl`; `check` catches a duplicated zero line from a kill mid-structure.
- Raw-shard rule (#973 k5): commit `settings.json` and `structures.jsonl` of each shard; archive the
  zeros with md5s.

### 9.7 Exact commands

```
cd /Users/bruce/dev/cyclers
for c in psb pnb pkb phb; do for i in $(seq 0 31); do echo "$c $i"; done; done |
  xargs -n2 -P2 sh -c 'mkdir -p data/1065_barycentric/$0 && uv run python scripts/run_1065_barycentric_enumerate.py enumerate --cell $0 --k 1,2,3,4,5,6 --out data/1065_barycentric/$0/s$1 --max-returns 2,0 --transfer-revs 0,1 --generic-revs 1,2 --n-phase 24 --n-split 10 --n-refine 60 --shard $1/32 --timing-pilot-s 3.4 > data/1065_barycentric/$0/s$1.log 2>&1'
# after all 128 shards, per cell (psb shown). Glob fixed 2026-10-12 (sec. 10.1): s*/, not s*, which also matched the s*.log files
uv run python scripts/run_1065_barycentric_enumerate.py check data/1065_barycentric/psb/s*/
uv run python scripts/run_1065_barycentric_enumerate.py gauntlet --cell psb data/1065_barycentric/psb/s*/ --out data/1065_barycentric/psb_gauntlet.json
uv run python scripts/run_1065_barycentric_enumerate.py controls --cell psb data/1065_barycentric/psb/s*/ --out data/1065_barycentric/psb_controls.json
```

## 10. #1065 results (enumeration: lead launch 2026-10-11 19:49-21:47 AEDT; check, controls, gauntlet: 2026-10-12 09:11-09:25 AEDT)

Data: `data/1065_barycentric/`:
- per shard, `settings.json` and `structures.jsonl` (committed);
- `<cell>_gauntlet.json`, `<cell>_controls.json` and `summary.json`;
- `zeros_md5.txt`: 128 shards, 12,472 zero lines.
The zeros and the shard logs are not committed. The full directory is archived at
`~/dev/references/cyclers-runs/1065/` (24 MB; `diff -rq` clean at archive time).

### 10.1 Logs and coverage

- Every one of the 128 shard logs ends in `DONE ... structure_errors=0 zero_assessment_errors=0`.
  No log contains a Traceback, an exception or a warning (grep over all 128).
- The Tracebacks the lead saw come from the sec. 9.7 command, not from the run. That command's glob
  `data/1065_barycentric/<cell>/s*` also matches the `s*.log` files, and `check` then fails on
  `s0.log/settings.json` (NotADirectoryError). The fix is the glob `s*/`, used for every step below.
  No shard needs a re-run.
- `check` (`s*/`): psb 3,462, pnb 3,346, pkb 3,146 and phb 3,025 structures. Every structure was done
  once, with no error lines, no duplicated zero line and none at the n_refine ceiling. All four cells
  are OK. Total 12,979, as in sec. 9.1.

### 10.2 Controls (sec. 9.4): both FAIL as pre-registered; examined one by one

| Cell | Continued members | Re-found (key + Charon V_inf within 1e-4 km/s) | Same verdict | Old-model V_inf hit |
|---|---|---|---|---|
| psb | 140 | 135 | 135 | 1 (index 103) |
| pnb | 127 | 124 | 124 | 0 |
| pkb | 59 | 58 | 58 | 0 |
| phb | 57 | 50 | 49 | 0 |
| total | 383 | **367** (threshold 370) | **366** (1 change) | **1** |

**Positive control: FAIL** (367 < 370, and one verdict change). **Negative control: FAIL by its
literal criterion** (one continued member's #998 V_inf is present).

Examination. Script: scratchpad `examine2.py`, run after the counted controls; the output is
recorded below. For each of the 16 misses plus phb 26:
1. Is the continued member an exact zero of the run's system?
2. Does `solve_structure` reach it at n_phase 24 (the run), 96 and 384?
3. Basin: from the member's dates perturbed by 0.05, 0.2, 0.5 and 1 d (8 random directions each),
   how often does `correct_dates` return to it?

| Member | 1024 pass | Residual at the member (km/s) | Grid finds it at n_phase 24 / 96 / 384 | Basin returns (of 8) at 0.05 / 0.2 / 0.5 / 1 d |
|---|---|---|---|---|
| psb 51 | yes | 5e-16 | no / no / no | 8 / 7 / 8 / 6 |
| psb 81 | yes | 2e-15 | no / yes / yes | 8 / 8 / 8 / 6 |
| psb 83 | yes | 2e-15 | no / yes / yes | 8 / 8 / 8 / 3 |
| psb 93 | yes | 6e-16 | no / yes / yes | 8 / 8 / 8 / 8 |
| psb 126 | yes | 3e-16 | no / yes / yes | 8 / 8 / 7 / 7 |
| pnb 91 | yes | 4e-16 | no / no / no | 8 / 8 / 7 / 4 |
| pnb 100 | no | 5e-16 | no / no / no | 8 / 8 / 8 / 6 |
| pnb 106 | no | 2e-15 | no / no / no | 8 / 8 / 2 / 0 |
| pkb 7 | no | 6e-15 | no / yes / no | 7 / 5 / 0 / 0 |
| phb 12 | yes | 3e-16 | no / no / no | 8 / 8 / 7 / 5 |
| phb 21 | no | 4e-16 | no / yes / yes | 6 / 1 / 0 / 1 |
| phb 22 | yes | 2e-16 | no / no / no | 4 / 2 / 1 / 0 |
| phb 23 | yes | 6e-17 | no / no / no | 2 / 5 / 1 / 0 |
| phb 24 | yes | 2e-16 | no / no / no | 2 / 4 / 1 / 2 |
| phb 41 | yes | 4e-16 | no / no / no | 8 / 8 / 6 / 2 |
| phb 42 | yes | 3e-16 | no / no / no | 8 / 8 / 3 / 2 |
| phb 26 (the verdict change) | yes | 7e-16 | no / no / no | 8 / 4 / 6 / 4 |

Findings:
- Every missed member is an exact zero of the barycentric system (residual <= 6e-15 km/s) with a real
  basin. The misses are recall misses of the seed grid (24 phases x the simplex split, the 60
  best-residual seeds refined), not model errors. Six are found at a finer phase grid; 11 are not
  found even at 384 phases. The grid's seed choice, not its density, loses them.
- **phb 26, the "verdict change", is a matching artefact.**
  - The zero matched by key + V_inf has the member's Charon and Hydra V_inf to 1e-15, but other
    dates (+2.46 d on the first three) and a failing gate (1.287 against the member's 0.311). It is
    another zero of the same structure.
  - The member itself is not found, so it is a 17th miss.
  - Strict count: 366 of 383 re-found (95.6 %), all with the same verdict.
- **psb 103, the negative-control hit, is a matching artefact too.**
  - The zero within 1e-4 km/s of the #998 V_inf (0.173348 against 0.173356 km/s) has other dates
    ([6.395, 27.725, 30.207] d against #998's [9.113, 23.141, 26.887]) and fails the gate (1.545).
    It is an exact zero of the barycentric system (8e-16), a different solution of the same
    structure.
  - The continued member itself (V_inf 0.19193) is present.
  - No #998 default-model solution is in the run: every zero in the run is a zero of the barycentric
    system by construction.
- Both artefacts come from the key + V_inf match rule: keys are not unique, and one structure carries
  several zeros (papercut 2026-10-11-gauntlet-keys-not-unique). The sec. 9.4 rule should have matched
  on dates (`_same`, modulo the phase period).
- Reading:
  - The run's recall of the known barycentric members is 95.6 %. All 17 missed members exist (the
    sec. 8 continuation holds them).
  - The ideal-model record for #1036 should be the union of the #1065 gate-passers and the 17 missed
    continued members. The #1024 assessment of those 17 stands (`data/1024_barycentric/`).
  - Whether the pre-registered FAIL blocks the use of the run's 97 new members is for the lead.

### 10.3 Gauntlet (sec. 9.2): gate-passers and strong cyclers

| Cell | Zeros | Physical cyclers | Gate-passing (pass def. 1-2) | of these with a #1024 parent / new | Pluto surface fails | Unscheduled screen: pass / direction-dependent / model-invalid / impact | strong (7,192 km SOI) | strong_8070 | **strong_1065** (with the screen) | new and strong_1065 |
|---|---|---|---|---|---|---|---|---|---|---|
| psb | 6,363 | 1,621 | 181 | 128 / 53 | 70 | 157 / 5 / 17 / 2 | 54 | 54 | **50** | 5 |
| pnb | 3,355 | 907 | 128 | 111 / 17 | 32 | 112 / 5 / 10 / 1 | 59 | 60 | **56** | 1 |
| pkb | 1,745 | 525 | 74 | 55 / 19 | 19 | 70 / 4 / 0 / 0 | 43 | 44 | **43** | 7 |
| phb | 1,009 | 336 | 52 | 44 / 8 | 14 | 49 / 2 / 1 / 0 | 22 | 22 | **22** | 1 |
| total | 12,472 | 3,389 | **435** | 338 / **97** | 135 | 388 / 16 / 28 / 3 | 178 | 180 | **171** | **14** |

- "With a #1024 parent" means a group containing a continued member, matched by dates (`_same`). 367
  re-found members fall into 338 groups, because the merge by flyby sequence joins some.
- No strong_1065 cyclers carries the direction-dependent flag: the 16 flagged gate-passers all fail
  another strong filter.
- Deep-dip flag (recorded, not gating) over the 171 strong_1065:
  - r_min / a_C is 0.018-0.83. The lowest, pnb at 306 km from the barycentre, is at least
    1,826 km from Pluto's centre on the opposite circle (638 km above its surface).
  - The Lambert-leg time inside 2 a_C is 5-50 %.
  - The patched conic about the 975.5 monopole is at its weakest there (sec. 8.5).

Against sec. 9.5:
- Gate-passers: expected at least 370-383 re-found plus new members. Got 435, with 97 new.
- Strong: expected of order 175 plus new members before the unscheduled screen. Got 178 before the
  screen (7,192 km), and 171 after.
- The screen removes 7 strong cyclers (4 psb, 3 pnb).

No catalogue edit. Nothing is called novel.

**Lead ruling (2026-10-12), recorded:**
- The proposed reading is accepted. The ideal-model record for #1036 is the union of the #1065
  gate-passers (435; strong_1065 171 after the screen) and the 17 missed continued members, whose
  #1024 assessment stands.
- The pre-registered control FAIL stays on record as written. It is read as a seed-grid recall limit
  (366/383 strict), not a model error.
- The 97 new members may be used. Each carries the recall caveat and the deep-dip flag.
- The recall follow-up (384-phase or basin-aware seeding) is registered as #1072 (backlog).
- The sec. 9.7 glob slip is fixed in place (`s*/`).

### 10.4 For the lead

1. FAIL as pre-registered, examined. Both controls fail. The examination finds recall misses of real
   barycentric zeros (17 of 383) and two key + V_inf matching artefacts; it finds no model error.
   - Proposed reading: use the run with the union rule of 10.2.
   - The 97 new gate-passers (14 strong_1065) are members with no old-model parent.
2. The sec. 9.7 glob is wrong (`s*` matches the logs). The fix is `s*/`. No re-run is needed.
3. Possible follow-up, not proposed for now: seed-grid recall. 11 of the 17 misses are not reached
   even at 384 phases, so the 60-best-seed selection is the limit. A recall test of the generator
   against the #1024 continued set is cheap (scripted above), if anyone wants the generator's recall
   on record.
