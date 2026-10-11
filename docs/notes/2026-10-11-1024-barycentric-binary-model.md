# #1024: a barycentric circular-model option for binaries in the two-working-body generator

Status: PRE-REGISTRATION (secs. 1-7), written and committed before any code or run. Results follow in
later sections. No catalogue edits. The #998 adjudication (#1036) is a later owner step. #949 (the
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
