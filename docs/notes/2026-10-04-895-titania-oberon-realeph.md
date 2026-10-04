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
