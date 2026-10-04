# #895: ballistic Titania-Oberon flyby arcs in a URA111 force model, rebuilt in the repository

Task `#895` (registered from the `#890` adversarial review,
`docs/notes/2026-10-04-890-adversarial-review.md` section 6). Question: near the `#890` model orbit
(`docs/notes/2026-10-04-890-titania-oberon-candidate.md`), do manoeuvre-free arcs with alternating
Titania and Oberon flybys exist in a Uranus-centred force model that tracks the URA111 satellite
ephemeris, at several epochs, reproducibly, from code in the repository? Code:
`src/cyclerfinder/search/titania_oberon_realeph_895.py`, driver
`scripts/screen_895_titania_oberon_realeph.py`, tests `tests/search/test_titania_oberon_realeph_895.py`;
outputs `data/found/895_titania_oberon_realeph/`.

Method lineage: the homotopy from mean circular ephemerides to the real ephemeris, with continuity
constraints solved at each step, is the continuation of Bradley & Russell (2014), J. Astronaut. Sci.
61:227 (digest `docs/notes/2026-10-04-digest-bradley-russell-2014-patched-conics-to-full-gravity-continuation.md`),
there applied to open tours from a patched-conic start. Here it starts from the `#890` four-body
orbit, not from a patched conic, so the mass part of their continuation is not used.

The `#890` reviewer's scratch implementation of the same computation was NOT read before this
section was committed and will not be read until this build's first converged real-ephemeris arc
and its numbers are committed (section 3 will record the commit).

## 1. Pre-registration (committed before the force-model positive control and before any arc was computed; not edited afterwards)

Before this commit only the following were done: the kernel's contents and header constants were
listed (body list, coverage, GM table, J2, J4, pole) and the background notes and the `#890` code
were read. No force-model control, no propagation of any moon or spacecraft, no arc.

### 1.1 Force model

Frame and units: centred on Uranus (NAIF 799), non-rotating, J2000 axes; km, km/s; time is TDB
seconds (SPICE ephemeris time). States from the kernels are geometric (no light time).

Kernel: `ura111.bsp` as shipped with the local GMAT R2022a installation
(`cyclerfinder.data.validation.v4_uranus_strict.DEFAULT_URA_PATH`). It holds the satellites 701 to
705 and Uranus 799 relative to the Uranus barycentre (URA111, R. A. Jacobson, released
2014-01-08) and the DE430 segments for the Sun (10), the Uranus barycentre (7) and Earth, so the
Sun relative to Uranus comes from the same file. Leap seconds: the repository's
`src/cyclerfinder/verify/kernels/naif0012.tls` (used only to convert calendar dates).

Constants: ONE set, the one URA111 was fitted with, read from the kernel's comment area ("Bodies
on the File" and "Additional Constants on the File"), hard-coded in the module and checked against
the kernel by a test that skips without the kernel:

| quantity | value | source |
|---|---|---|
| GM Uranus (planet only) | 5,793,951.322279009 km^3/s^2 | URA111 header |
| GM Miranda, Ariel, Umbriel | 4.319516899232100, 83.46344431770477, 85.09338094489388 | URA111 header |
| GM Titania, Oberon | 226.9437003741248, 205.3234302535623 | URA111 header |
| GM Sun | 1.327132332639e11 | URA111 header (GM10) |
| J2, J4 at R = 25,559 km | 3.510685384697763e-3, -3.416639735448987e-5 | URA111 header (J702, J704, RADIUS) |

The registry (`core/satellites.py`) GMs differ in the fourth figure (Titania 226.9) and are not
mixed in. The repository's `v4_uranus.URANUS_J2` (3509.291e-6, French et al. 2024) differs from
the kernel's J2 by 4.0e-4 relative; the kernel's value is used because the positive control is
against the kernel, and the French value is run as a reported variant of the control. The kernel's
J2 is the "3510.68e-6 from memory" of the `#890` review, which this identifies.

Uranus: point mass plus the zonal harmonics J2 and J4 about the pole. Including J4 goes beyond the
brief's "point mass plus J2"; it is declared here as a deviation, and the positive control is
reported with and without it. Pole: the IAU 2009 value in the GMAT PCK
(`SPICEPlanetaryConstantsKernel.tpc`, BODY799_POLE_RA = 257.311 deg, BODY799_POLE_DEC = -15.175 deg,
no rate terms), held constant in J2000 axes. Only the axis enters a zonal field, not its sense. The
URA111 header carries its own pole (ZACPL7 77.30990 deg, ZDEPL7 15.17246 deg, the opposite sense,
with rate terms DACPL7, DDEPL7 whose units the header does not state); its angle from the IAU axis
will be reported, and the control will be repeated with that pole held constant. The rate terms are
not interpreted unless the control discriminates between interpretations.

Moons: Miranda, Ariel, Umbriel, Titania and Oberon as UNSOFTENED point masses at their URA111
positions relative to 799, with the direct term on the spacecraft and the indirect term (the
acceleration of Uranus by that moon, GM_k r_k/|r_k|^3). The Sun likewise, direct and indirect. The
force is never switched off or softened near a body. Positions come from a table sampled from SPICE
every 600 s (positions and velocities) and evaluated by cubic Hermite interpolation; a test bounds
the interpolation error against direct SPICE calls at mid-knots (Miranda, the worst case, below 5 m;
the Sun and the other moons below 1 m).

Left out, with the reason: the reaction of the moons on Uranus's oblate figure (the indirect part of
J2; relative size GM_k/GM_U times the J2 acceleration at the moon, below 1e-9 of the central
attraction); the other planets (Neptune's tide at Uranus is below 1e-3 of the Sun's); the moons'
own figures; Puck and the small moons (GM 0 in the kernel); relativity. The spacecraft is massless.

Integrators: the solver integrates each segment and its 6x6 state transition matrix with an
adaptive DOP853 written for this module (the coefficients of scipy's DOP853, error control on the
six state components, rtol 1e-12, compiled with numba). Verification uses scipy's own DOP853 and
scipy's LSODA (a multistep method), with the same force function but not the module's integrator.

### 1.2 Model-level checks that must pass first (step 1)

* **P1, force-model positive control (pass/fail).** Each moon in turn is the test body: start it at
  its own URA111 state, remove it from the perturber list, add its GM to the central term (the
  relative motion of a massive moon about a massive planet), propagate with the full model (J2, J4,
  the other four moons, the Sun), and compare with URA111. Five start epochs: the five arc epochs of
  1.4. **Pass if every moon at every epoch is within 10 km of URA111 after 30 days and within 40 km
  after 123 days.** (The review's model reached 0.4 to 4 km and 2 to 17 km.) **Discrimination
  (pass/fail):** with J2 = J4 = 0, and separately with the repository's former pair (J2 = 3.34343e-3
  at 25,559 km, J4 = 0), at least Miranda and Ariel must FAIL the 30-day bound at every epoch;
  otherwise the control does not see the field and P1 does not count. Reported, not pass/fail: the
  same control without the Sun, without J4, with the French et al. J2, with the kernel's pole, and
  with the moon's own GM left out of the central term.
* **P2, the `#890` model orbit from independent code (pass/fail; a permanent test that needs no
  kernel).** The same force function and integrator as the arcs, fed by an analytic body provider
  that encodes the `#890` model (planar; Titania on a circle about Uranus at the registry radius;
  Oberon on a circle of its registry radius about the Uranus-Titania barycentre at its `#890` rate;
  central GM = registry system GM minus Titania and Oberon; indirect term = the barycentre's
  acceleration plus Uranus's motion about it; no J2, no Sun, no inner moons), started from the
  refined state in `data/found/890_titania_oberon_candidate/refine.json` converted by this module's
  own frame code. Pass if (i) at the half cycle the `#890` symmetry conditions hold: in the `#890`
  rotating frame |y| <= 10 m and |v_x| <= 1 mm/s; (ii) the periapsis altitudes are 1,976.9 +/- 1 km
  (Titania) and 1,364.2 +/- 1 km (Oberon); (iii) after one cycle the state returns to within 1 km
  and 1 cm/s of the start in the rotating frame (the leading multiplier is 8.4e5, so integration
  noise alone gives tens of metres; 1 km is the bound that is not ill-posed). **Negative control
  (must fail):** with Oberon on a circle about Uranus instead of the barycentre (a 17 km change),
  the one-cycle return error must exceed 100 km, otherwise the test cannot see a model difference
  of that size and P2 does not count.
* **P3, unit tests (pass/fail):** the analytic gravity gradient (point masses, J2, J4) against
  central differences of the acceleration; the zonal accelerations against their closed forms on
  the pole axis and in the equator; the module's integrator against scipy's DOP853 on the same
  right-hand side; its state transition matrix against central differences of propagations; the
  ephemeris table against SPICE (bounds above).

If P1 or P2 fails, no arc result is reported as a result of this model; the failure is reported.

### 1.3 The correction

* Arc layout for N cycles from a start epoch: nodes at fixed epochs t_k = t_c - 2h + k h,
  k = 0 .. 24N + 4, with h = T_cyc/24, T_cyc = 5 synodic periods of the fitted mean motions
  (about 123.19 d, h about 5.13 d), t_c a mean Titania-Oberon conjunction (1.4). The intended
  flybys are Titania near t_c + j T_cyc (j = 0 .. N) and Oberon near t_c + (j + 1/2) T_cyc
  (j = 0 .. N-1): 2N + 1 flybys, two extra nodes (about 10 d) beyond the first and the last.
* Unknowns: all six components of every node state; node epochs are fixed (no time variables).
  Constraints: continuity of position and velocity at each of the 24N + 4 junctions. No manoeuvre
  variable exists anywhere. The system is underdetermined by six; each Newton step is the
  MINIMUM-NORM solution of the linearised continuity equations, in scaled variables (position
  divided by 436,300 km, velocity by 3.64 km/s, roughly Titania's distance and speed).
* Jacobian: the segment state transition matrices from the variational equations (analytic
  gravity gradient), no finite differences.
* Newton stops when every junction is below 1 cm and 1e-8 km/s (0.01 mm/s) as computed by the
  solver's integrator, or after 12 iterations (failure). A step that increases the largest
  junction residual is halved up to four times (failure after that).
* Homotopy, parameter lam from 0 to 1:
  - lam = 0: Titania and Oberon on circular orbits in one plane (Titania's mean orbital plane over
    the fit window, +z along Titania's mean orbital angular momentum, which is checked to be within
    1 deg of the antipode of the IAU north pole), radius = the mean URA111 distance and longitude a
    linear least-squares fit to the URA111 longitude in that plane, both over the fit window
    [D - 30 d, D + 12 T_cyc + 30 d] for target date D; Miranda, Ariel and Umbriel lumped into the
    central GM; J2 = J4 = 0; no Sun.
  - lam: Titania and Oberon at (1 - lam) circle + lam URA111 (Bradley & Russell's linear blend);
    Miranda, Ariel, Umbriel at URA111 positions with GM lam GM_k, central GM
    GM_U + (1 - lam)(sum of the three); J2, J4 and the Sun's GM multiplied by lam. At lam = 1 this is
    exactly the force model of 1.1.
  - Steps: first step 0.1; predictor = secant extrapolation from the last two converged solutions
    (the previous solution for the first step); a failed step is halved, down to a minimum of
    1/512; after a step that converged in 4 iterations or fewer the next step is 1.5 times larger,
    capped at 0.25. At most 60 attempted steps per run. Every step is recorded (lam, step,
    iterations, residual history, flyby summary). The run is checkpointed after every step.
  - **Branch check at every converged step:** the list of approaches within 2 Hill radii of
    Titania or Oberon must be exactly the 2N + 1 intended flybys in the order T, O, T, ..., T; each
    must keep the side it had at lam = 0 (sign of (r - r_moon) . r_moon at periapsis, i.e. outside
    or inside the moon's orbit) and the turn sense (sign of the flyby's angular momentum on the
    reference +z). A step that fails the check is treated as a failed step (halved), and the event
    is recorded. A branch change is never accepted.
* Seed at lam = 0: the `#890` refined orbit, converted from its rotating frame to Uranus-centred
  inertial states with this module's code, lengths scaled by the ratio of the fitted Titania radius
  to the registry radius, time scaled so that its cycle equals T_cyc, placed in the reference plane
  with Titania's direction matching the fitted circle at the node epoch, repeated N times. If
  Newton fails at lam = 0 from this seed, a fallback is allowed and must be reported: a homotopy
  in the circle constants from the `#890` registry circles to the fitted circles.

### 1.4 Start epochs

Five epochs, all decided now: the mean Titania-Oberon conjunction (from the fitted circles of 1.3)
nearest to each target date D = 2030-01-12, 2031-06-13 (the review's two dates), 2035-07-01,
2040-07-01 and 2045-07-01 (00:00 TDB). Spread 15.5 years. Arcs of N = 3 cycles at each.

### 1.5 Pass/fail criteria for a converged arc (all on the lam = 1 arc)

* **(c) No manoeuvre.** Every segment re-integrated serially with scipy DOP853 (rtol 1e-13,
  atol 1e-9 km and 1e-15 km/s) and separately with scipy LSODA (rtol 1e-13, same atol): every
  junction discontinuity below 1 m and 1 mm/s under EACH. LSODA's self-convergence (rtol 1e-12
  against 1e-13) is reported with it; if it exceeds 0.3 m at any junction, LSODA is declared not to
  resolve the problem at that tolerance and that is reported as a failure of the verifier, not a
  pass.
* **(d) Every flyby integrated through.** For flyby j at node index 2 + 12 j, one continuous
  integration from node 2 + 12 j - 6 to node 2 + 12 j + 6 (the middles of the neighbouring legs,
  about 61.6 d; from node 0 for the first flyby and to the last node for the last), with each of
  the two verification integrators. Pass if the end lands within 1 km and 1 cm/s of the stored
  node. Prediction (stated so that the limit is well posed): the landing error expected from the
  junction residuals alone is d_pred = sum over junctions k inside the span of
  Phi(t_end, t_k) delta_k, with delta_k the measured junction jump and Phi the product of the
  solver's segment matrices. The limit is applied only if |d_pred| <= 0.3 km and 0.3 cm/s; if the
  prediction exceeds that, the Newton tolerance is tightened and the arc re-converged before (d) is
  evaluated, and if that is not possible (d) is reported as not evaluable, not as a pass.
  Consistency (pass/fail): the observed landing error must not exceed
  |d_pred| + 3 x (difference between the two integrators' continuous end states) + 10 m; a larger
  error is something the junctions do not explain.
* **(e) Flybys.** Periapsis of each of the 2N + 1 flybys found by root-finding on the range rate of
  the verified DOP853 trajectory relative to the moon's position from SPICE. Pass if the altitude is
  at least 50 km above the surface (Titania 788.9 km, Oberon 761.4 km radius) and the osculating
  moon-centric energy is positive (hyperbolic). Reported for each: date, altitude, speed relative to
  the moon at periapsis, osculating excess speed sqrt(v^2 - 2 GM/r), osculating eccentricity and
  turn, side, and distance out of the reference plane.
* **(f) Other approaches.** Every local minimum of the distance to each of the five moons over the
  arc within 2 Hill radii (r_H = a (GM_k/(3 GM_U))^(1/3), a the mean distance) other than the
  2N + 1 targeted flybys is reported, and the closest approach to each moon and to Uranus. Fail only
  on an approach below 50 km altitude of any moon or below 1,000 km altitude of Uranus.
* **(b) Coverage.** The result "arcs exist at several epochs" holds only if, among the epochs whose
  N = 3 arc passes (c), (d), (e) and (f), at least three exist whose start dates span at least ten
  years. Evaluated on the passing set only.
* **(g) Extension.** From every epoch whose N = 3 arc passes, the whole homotopy is run again with
  N = 6 and N = 12 (same epoch, same fit). If N = 12 fails and N = 6 passes, N = 9 is tried. The
  largest N that converges and passes (c) to (f) is reported, and for each failure where it stopped:
  lam, step size, residual, the branch-check event or the flyby that went below the floor, or the
  verification criterion that failed. Reported whether it succeeds or not; there is no pass/fail.
* **Reported, not pass/fail:**
  - The Sun's effect: from each converged lam = 1 arc, one Newton solve with the Sun's GM set to
    zero, started from that arc; the change in each flyby altitude and the largest node
    displacement. (Separate homotopies are not compared, because with six free parameters two
    converged arcs in the same model can differ.)
  - Comparison with the review's arcs (only after this build's first arc is committed): a variant
    with the review's model (no Sun, no J4, J2 = 3.510685e-3) at the review's two start times and
    node layout, continued from this build's converged arc at that epoch where the layouts agree;
    flyby dates and altitudes compared. Reported only, because the arcs form a six-parameter family
    and the minimum-norm step picks one member depending on the path.

### 1.6 What counts as failure, and that it will be reported

* P1 or P2 failing: the model is not established; any arc is reported only as a computation in an
  unverified model.
* An epoch's homotopy reaching the minimum step without converging, a branch change that cannot be
  avoided, or a converged arc failing (c), (d), (e) or (f): that epoch fails at N = 3, reported with
  the stage and the numbers.
* Fewer than three passing epochs spanning ten years: criterion (b) fails, and the result is
  "ballistic arcs were found at the following epochs only", or "none found", stated as such.
* Failures of the corrector are failures of the method at that epoch, not evidence that no
  ballistic arc exists; they are reported in those words.
* Nothing in this task is a catalogue row or a claim about one.

## 2. Results (run after commit `f3860e91`, which holds section 1)

Sections 2.1 and 2.2 were written and committed before the `#890` reviewer's scratch code was
opened. Everything here was computed by this build unless marked READ or INFERRED.

### 2.1 Step 1: model checks (commit `71630e05`; `control.json`, `p2.json`)

P1 passes at all five epochs: with the full model every moon is within 0.01 to 1.68 km of URA111
after 30 days and within 0.0 to 6.8 km after 123 days (worst: Ariel at E4, 1.68 / 6.8 km). The
control discriminates: without the zonal field Miranda and Ariel are 6,920 to 6,990 km and 2,713 to
2,719 km off at 30 days, and with the repository's former J2 pair 333 to 336 km and 128 to 130 km.
Without the Sun Titania and Oberon are 0.9 to 2.0 km and 0.9 to 3.6 km off at 30 days (5 to 16 km
at 123 days) against 0.05 to 0.17 km and 0.01 to 0.75 km with it: the Sun is visible in the moons'
own motion at this level. Without J4 Miranda degrades from 0.1 to 0.5 km to 2.9 to 3.6 km. The
French et al. J2 is worse than the kernel's for Miranda and Ariel (2.5 to 3.2 km and 0.6 to 1.1 km),
as expected for an ephemeris fitted with the kernel's value. The kernel's own pole changes nothing
material. Leaving the moon's own GM out of the central term costs 25 to 730 km at 30 days (the
control needs it). Interpolation error of the 600 s Hermite table at mid-knots: Miranda 0.31 m,
Ariel 0.05 m, the others 1.4 mm or less. The fitted mean plane of Titania is 0.037 to 0.058 deg from
the antipode of the IAU pole. The mean conjunctions found by this build for the review's two dates
are 2030-01-11 23:57:20 TDB and 2031-06-13 09:12:10 TDB (the review gives 23:57 and 09:13, READ),
an agreement between two independent fits.

P2 passes with the module's integrator at rtol 1e-13: half-cycle |y| = 1.4 cm and |v_x| = 0.0014
mm/s; flyby altitudes 1,976.86 km (Titania) and 1,364.18 km (Oberon); one-cycle return 12.8 m and
0.089 cm/s (scipy DOP853 1e-13: 12.1 m; scipy Radau 1e-12: 0.63 m; these reproduce the `#890`
figures 12 m and 0.68 m). At the solver's rtol 1e-12 the return is 125 m and 0.86 cm/s, which
shows that the 1 cm/s part of P2(iii) was close to ill-posed: the return velocity error tracks the
position error along the unstable direction at about 7e-5 per second, so 1 cm/s corresponds to
about 145 m. scipy LSODA at rtol 1e-13 returns 3.0 km and 21 cm/s after one cycle: over a whole
cycle of this orbit LSODA is far less accurate than DOP853. Negative control: with Oberon circling
Uranus instead of the barycentre the half-cycle |y| is 8.5 km and the trajectory strikes Titania
0.12 km from its centre at 122.62 d; there is no return.

### 2.2 First converged real-ephemeris arc: E1, N = 3 (`arcs/E1_N3*`, `verify_E1_N3.json`)

Route. Newton at lam = 0 from the scaled `#890` seed did NOT converge (step halving exhausted at
iteration 3, residual 2,900 km), so the pre-registered fallback was used: the `#890` orbit placed on
the `#890` circles converged at sigma = 0 in one Newton step (seed junctions 0.56 km, from the
interpolation of the stored orbit), then the constants homotopy sigma = 0 to 1 took 39 attempts
(29 accepted; the first steps had to be cut to 0.003), then lam = 0 to 1 took 7 attempts (6
accepted), every converged step with the same flyby count, order, sides and turn sense.

The arc moved a long way along the family during the constants homotopy, not during the ephemeris
homotopy. Flyby altitudes (T, O, T, O, T, O, T), km: at sigma = 0, 1,977 / 1,364 / ... (the `#890`
orbit); at sigma = 1 (circular model with the URA111 mean motions and radii), 5,454 / 4,073 / 6,519
/ 4,991 / 8,541 / 3,609 / 6,709; at lam = 1, the values below. The ephemeris homotopy moved them by
less than 400 km.

Verification at lam = 1 (criteria of section 1.5):

| criterion | result | numbers |
|---|---|---|
| solver junctions | | 7.0 mm, 1.2e-4 mm/s |
| (c) DOP853 1e-13 / LSODA 1e-13 | PASS | 11.7 mm, 1.1e-4 mm/s / 9.6 mm, 6.8e-5 mm/s; LSODA 1e-12 vs 1e-13: 6.9 cm |
| (d) seven flybys integrated through | PASS | landing within 6 m (DOP853) and 45 m (LSODA), at most 0.027 cm/s; predictions 0 to 5 m |
| (e) altitudes, hyperbolic | PASS | all hyperbolic; altitudes below |
| (f) other approaches | PASS | none within 2 Hill radii; closest Umbriel 180,322 km, Ariel 254,476, Miranda 312,050; Uranus 441,699 to 580,599 km |

| flyby | date (TDB) | altitude km | speed at periapsis km/s | osculating excess km/s | e | side | spacecraft z km |
|---|---|---|---|---|---|---|---|
| Titania | 2030-01-11 20:23 | 5,562 | 0.314 | 0.166 | 1.77 | outside | -850 |
| Oberon | 2030-03-14 04:25 | 4,221 | 0.351 | 0.203 | 2.00 | inside | 796 |
| Titania | 2030-05-14 18:43 | 6,827 | 0.277 | 0.132 | 1.58 | outside | -1,337 |
| Oberon | 2030-07-14 23:55 | 4,946 | 0.328 | 0.188 | 1.99 | inside | 976 |
| Titania | 2030-09-14 18:30 | 8,945 | 0.222 | 0.050 | 1.11 | outside | -948 |
| Oberon | 2030-11-15 01:09 | 3,393 | 0.379 | 0.211 | 1.90 | inside | 634 |
| Titania | 2031-01-15 20:34 | 6,734 | 0.272 | 0.118 | 1.46 | outside | -223 |

So a manoeuvre-free arc with seven alternating flybys exists at E1 in this model, but it is NOT
close to the `#890` orbit's flyby geometry: its flybys are two to five times higher and slower
(osculating excess speed 0.05 to 0.21 km/s against 0.26 to 0.27 km/s), and the fifth flyby, at
8,945 km (about 0.95 Hill radii) with an excess speed of 0.05 km/s, is barely an encounter. The
minimum-norm continuation picked this member of the six-parameter family; another route could
pick a different one. INFERRED: the drift happened while the circles' mean motions changed by
about 2e-4 (the `#890` cycle is 123.162 d, the URA111 one 123.187 d), which the open arc absorbs by
changing its flybys rather than its period.

### 2.3 Epochs E2 to E5, N = 3, pre-registered route

All four converged at lam = 1 and pass (c), (d), (e) and (f) (`verify_E*_N3.json`). E2 converged
at lam = 0 directly from the scaled seed; E3, E4 and E5 needed the constants homotopy, as E1 did.
With E1 the five epochs span 2030-01-11 to 2045-06-23, so criterion (b) passes for this route. The
flyby altitudes differ a great deal between epochs, for the route reason described in 2.2 (numbers
in section 4).

## 3. Deviation D1, declared before it was run: a periodic start at lam = 0

Why. Section 2.2 shows that the pre-registered route reaches lam = 1 but lets the open arc drift
along its six-parameter family while the circles' constants change, so the arcs it returns are not
"near the `#890` orbit" in their flyby geometry, and the drift differs from epoch to epoch. The
cause of the failed direct start is also identified (INFERRED, arithmetic): the seed's velocity
scaling s q = 0.99976 is not a similarity of a Kepler orbit (which needs s^(-1/2) = 1.00002), about
1 m/s of error, which a 5-day node interval on a 9-day e-folding orbit turns into the hundreds of
kilometres seen at the seed's junctions. Neither point changes section 1; both are reported as
faults of the pre-registration.

D1 replaces only the lam = 0 starting arc. Everything after it (lam homotopy, Newton settings,
branch check, criteria (c) to (f), extension (g)) is as pre-registered.

* Periodic shooter at lam = 0: one cycle, 24 nodes at t_c + k h(sigma), h(sigma) = T_cyc(sigma)/24
  with T_cyc(sigma) = 10 pi / (n_T(sigma) - n_O(sigma)) from the blended mean motions, anchored at
  the conjunction t_c; continuity at the 23 interior junctions and the wrap condition
  phi(x_23) = R x_0, R the rotation about the reference +z by n_T(sigma) T_cyc(sigma) (the circular
  model maps to itself under that time shift and rotation). 144 equations, 144 unknowns, the same
  Newton, step-halving and continuation rules as section 1.3, the same branch check on the tiled
  orbit.
* Control at sigma = 0: the unscaled `#890` orbit must satisfy the wrap condition to its
  interpolation noise before any Newton step, and converge in at most two steps.
* Continue sigma from 0 to 1 (the `#890` circles to the URA111-fitted circles). Independent check,
  stated now: at sigma = 1 the periodic orbit's flyby altitudes should be within 50 km of
  1,800 to 1,820 km (Titania) and 1,255 to 1,266 km (Oberon), the range READ from the review
  (section 2 table and arc A's circular end, computed there with URA111 mean motions). A larger
  difference is investigated before anything else is run.
* The periodic orbit tiled N times is the lam = 0 seed of the open arc (same node layout as 1.3);
  Newton at lam = 0, then the lam homotopy, at all five epochs, N = 3, 6 and 12 (9 if 12 fails).
* Reported per arc in addition: the largest absolute difference between a flyby's altitude and the
  periodic orbit's altitude at sigma = 1 ("nearness").
* Every arc is labelled by route ("pre-registered" or "D1") and criterion (b) is evaluated per
  route.
* The Sun-off solve and the review-model comparison (section 1.5) are started from D1 arcs.
* The `#890` reviewer's scratch code is still unread; it will be opened only after the first D1
  arc is committed.

### 3.1 First D1 results (written before the review's scratch code was opened)

Periodic orbit at E1 (`arcs/E1_N1_D1periodic*`). Control at sigma = 0: the unscaled `#890` orbit
satisfied the wrap condition to 12 m (largest junction 18 m, the interpolation of the stored orbit)
and converged in one Newton step, with flybys at 1,976.9 and 1,364.2 km. (A first run failed
because the closed-chain Jacobian subtracted the wrap rotation at every junction instead of only at
the last; fixed before any result was used.) Its multipliers in three dimensions at sigma = 0:
8.38e5 and 1.19e-6, 1.786 and 0.560 (the planar pairs, as in `#890`), and an out-of-plane pair
24.3 and 0.041, so the out-of-plane motion is unstable too. Sigma 0 to 1 took 22 accepted steps
with the flyby sides unchanged; at sigma = 1 the periodic orbit's flybys are at 1,805.9 km (Titania)
and 1,256.7 km (Oberon). The pre-stated check (1,800 to 1,820 and 1,255 to 1,266 km within 50 km)
passes; these agree with the review's figures, computed there by other code.

D1 arc E1, N = 3 (`arcs/E1_N3_D1*`, `verify_E1_N3_D1.json`): Newton at lam = 0 from the tiled
periodic orbit converged, the lam homotopy took 10 accepted steps, every one with the same
signature. At lam = 1 it passes (c) (DOP853 13 mm, LSODA 31 mm, at most 0.0008 mm/s), (d) (landing
within 10 m with DOP853 and 18 m with LSODA, predictions up to 8 m), (e) and (f) (no other approach
within 2 Hill radii; closest Umbriel 170,374 km).

| flyby | date (TDB) | altitude km | speed at periapsis km/s | osculating excess km/s | e |
|---|---|---|---|---|---|
| Titania | 2030-01-11 23:09 | 1,588.6 | 0.515 | 0.273 | 1.78 |
| Oberon | 2030-03-14 14:21 | 1,144.2 | 0.540 | 0.275 | 1.70 |
| Titania | 2030-05-15 04:16 | 1,764.4 | 0.499 | 0.267 | 1.80 |
| Oberon | 2030-07-15 19:06 | 1,082.6 | 0.547 | 0.277 | 1.69 |
| Titania | 2030-09-15 09:32 | 1,752.9 | 0.500 | 0.267 | 1.80 |
| Oberon | 2030-11-16 00:20 | 1,075.8 | 0.548 | 0.278 | 1.69 |
| Titania | 2031-01-16 15:12 | 1,582.7 | 0.516 | 0.274 | 1.78 |

Pre-registered route at N = 6 (all five epochs): all converge and pass (c) to (f)
(`verify_E*_N6.json`); numbers in section 4.

### 3.2 Deviation D2, declared before it was run: Newton tolerance 5 cm

Observed (`arcs/E2_N6_D1.json`, `arcs/E4_N6_D1.json`, `arcs/E1_N12_D1.json`): three D1 runs failed
without any branch event and without any sign of divergence. In each, Newton reached a junction
residual of 1.1 to 2.3 cm and could not go lower (step halving exhausted with the residual flat at
that level), so it never met the pre-registered stop of 1 cm; repeated failures then shrank the
homotopy step below the minimum (E2 at lam = 0.577, E4 at lam = 0.216) or stopped the run at
lam = 0 (E1, N = 12, 2.3 cm). That level is the noise floor of the solver's integrator at
rtol 1e-12 (adaptive step selection makes the segment map non-smooth at the 1e-12 relative level,
and a flyby segment amplifies that by tens). The 1 cm stop was set too close to it; the
pre-registration did not measure the floor first.

D2: Newton stops at 5 cm and 5e-8 km/s instead of 1 cm and 1e-8 km/s. Nothing else changes; every
criterion of section 1.5 is applied unchanged (junctions under 1 m and 1 mm/s with both
verification integrators, flybys integrated through with the prediction check). A 5 cm residual
gives a (d) prediction of order 50 m at the per-leg growth of about 1,000, well inside the 0.3 km
well-posedness bound. D2 is used only for runs that are labelled with it; the failed runs under the
pre-registered tolerance stay recorded as failures of the method at that tolerance.
