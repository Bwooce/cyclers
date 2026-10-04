# #890: is the Titania-Oberon-Titania closure a real periodic orbit of the four-body model?

Task `#890` (registered from `#888`). Candidate: the one symmetric two-moon closure that passes the
demanded-turn gate (`data/found/888_turn_gate/candidates.jsonl`, first record; note
`docs/notes/2026-10-04-888-demanded-turn-gate.md` section 5). Code: `src/cyclerfinder/search/two_moon_periodic_890.py`,
`scripts/screen_890_titania_oberon_candidate.py`, tests `tests/search/test_two_moon_periodic_890.py`;
outputs `data/found/890_titania_oberon_candidate/`.

## 1. Pre-registration (written and committed before steps 2 to 5 were run)

### 1.1 Model

Planar concentric circular restricted four-body problem, `core/ccr4bp.py` unmodified, through the
vectorised planar right-hand side of `search/ccr4bp_strob_connection.py`. Titania is the base moon
(second primary of the rotating frame), Oberon the periodic perturber. Units: length a_Titania =
436,298 km, time 1/n_Titania. Constants from `core/satellites.py` (GM_Titania 226.9, GM_Oberon 205.3,
radii 788.9 / 761.4 km, a 436,298 / 583,511 km). Mass scale lam multiplies both moon GMs; the planet
GM is `PRIMARIES["Uranus"] - lam (GM_T + GM_O)` so that the total is the registry system GM at every
lam (that GM already contains Miranda, Ariel and Umbriel, which orbit inside the spacecraft's orbit
and act on it approximately as central mass). At lam = 0 the moons move at exactly the patched-conic
mean motions of the `#888` enumeration. Uranus J2, the moons' eccentricities and inclinations, the
other moons' individual gravity and moon-moon forces are NOT in the model (step 5 lists what they
would need).

The model's forcing (Titania-Oberon synodic) period depends on lam through the moon masses in the
mean motions; the cycle is always taken as five model forcing periods, never a hard-coded 123.16 d.

### 1.2 Conventions (reconciliation of the two descriptions)

`#888`'s enumeration puts Titania and Oberon both at longitude 0 at t = 0 ("rel_offset 0"): the
Titania flyby is at t = 0 with both moons on the same side. The leg time is 2.5 synodic periods, so
the relative angle has advanced by 5 pi at the middle (Oberon) flyby: Oberon is then on the opposite
side of Uranus from Titania. The coordinator's "far moon on the opposite side of the symmetry axis,
koff = 180 degrees" is measured at the middle encounter and is the same configuration. In the
rotating frame both moons lie on the x-axis at tau = 0 (Oberon at angle 0) and at tau = T = 2.5 P
(Oberon at angle pi), so the planar model is invariant under (tau, x, y, vx, vy) -> (-tau, x, -y,
-vx, vy) about both instants. A trajectory crossing the x-axis perpendicularly at tau = 0 and at
tau = T is periodic with period 2T = 5P, a fixed point of the five-period stroboscopic map. The
symmetry is checked numerically before it is used (test `test_mirror_time_reversal_symmetry`).

### 1.3 Why the plan's mass continuation from zero is changed

The task plan says to continue "from zero (where the Kepler closure is exact)". It is not exact
there: with massless moons nothing can turn the trajectory, and the patched-conic closure has
velocity jumps of 42 and 58 degrees at the flybys. The patched conic is the SINGULAR limit lam -> 0
in which the flyby periapsis shrinks like lam (at fixed V-infinity and turn, r_p = GM (e - 1)/V^2).
The continuation therefore starts at small positive lam (1e-3), with the flyby hyperbola of the
patched-conic turn scaled to that lam as the guess, where the patched conic is an excellent
approximation (r_p / r_Hill ~ lam^(2/3)). The lam = 0 check is kept, but as a test of the frame
conversion only.

### 1.4 Method

* Guess at each flyby: the moon-centric hyperbola that turns the patched-conic incoming V-infinity
  onto the outgoing one (periapsis along v_in_hat - v_out_hat, periapsis velocity along v_in_hat +
  v_out_hat), placed at the moon at the flyby epoch. No node is ever seeded with a state aimed at a
  moon's centre. Interior nodes: the leg-0 Kepler conic at the same fraction of the leg.
* Symmetric multiple shooting: unknowns (x0, vy0) at the Titania perpendicular crossing tau = 0,
  M = 11 interior 4-state nodes at tau_k = k T/12, and (xT, vyT) at the Oberon perpendicular
  crossing tau = T; segments forward from tau = 0 and between interior nodes, backward from tau = T;
  4M + 4 equations, 4M + 4 unknowns; Newton with the variational equations, step cap and
  backtracking.
* Continuation in lam: 1e-3, 2e-3, 5e-3, 0.01, 0.02, 0.05, 0.1, 0.2, ..., 1.0, halving the step on
  failure down to a minimum step of 1e-4 (relative to lam). The sign of det J is tracked; a sign change
  or a step collapse is reported as a fold at that lam. Also tried: Newton directly at lam = 1 from
  the lam = 1 hyperbola guess.
* Positive control for the corrector: it must converge at lam = 1e-3 from the patched-conic guess,
  and the converged flyby periapses there must agree with the patched-conic hyperbola (r_p and turn)
  to a few percent. If it fails there, a later non-convergence says nothing about the physics.

### 1.5 Step 2 (first decisive test): what is reported

Two propagations of one cycle (5 P) at lam = 1, DOP853 rtol = atol = 1e-13: (a) the literal stored
conic state at the leg-0 apoapsis nearest mid-leg (this conic was aimed at the moons' centres, so
the result partly measures that artefact); (b) the honest patched-conic realisation: the Titania
hyperbola periapsis state at tau = 0. For each: every approach to either moon inside 2 Hill radii
(distance, time, relative speed, moon-centric osculating eccentricity), and the state after one
cycle against the start. Reporting only; no pass/fail. The decision about "in the neighbourhood of
something that closes" is made by step 3.

### 1.6 Pass/fail criteria for an orbit to count as the candidate, a periodic two-moon orbit

All at lam = 1 (physical masses). Lengths: 1 km = 2.29e-6 nondim; 1 cm/s = 2.74e-6 nondim.

* **V1 corrector convergence**: max |residual| of the multiple-shooting system <= 1e-10 (nondim).
* **V2 closure, DOP853**: one continuous DOP853 integration (rtol = atol = 1e-13) from the tau = 0
  node over 5P. Pass if the return differs from the start by <= 1 km in position and <= 1 cm/s in
  velocity. If the monodromy's largest |eigenvalue| Lambda exceeds 1e6, a closure failure may be the
  instability amplifying round-off rather than absence of an orbit; then V2 passes instead if the
  half-cycle symmetric conditions |y(T)| and |vx(T)| meet the same 1 km / 1 cm/s bounds AND the
  full-cycle closure error is <= 10 Lambda (1e-10 + 1e-12) nondim.
* **V3 closure, Radau**: the same with Radau (rtol 1e-12, atol 1e-13, analytic Jacobian); and the
  Radau and DOP853 states at tau = T agree to <= 1 km and <= 1 cm/s.
* **V4 it is a two-moon flyby orbit**: at each of the two targeted encounters per half cycle,
  (a) periapsis altitude >= 50 km (the project floor for Titania and Oberon,
  `data/flyby_altitude_references.yaml`, Heaton-Longuski 2003 design floor); (b) periapsis distance
  <= the Laplace sphere of influence a (m/M)^0.4 (Titania 7,532 km, Oberon 9,678 km); (c) the
  moon-centric osculating eccentricity at periapsis > 1 (a flyby, not a temporary capture).
* **V5 demanded-turn gate on the corrected orbit** (`verify/turn_gate.py`, project floor 50 km), with
  the V-infinity vectors extracted as the moon-relative inertial velocity at the inbound and the
  outbound crossing of the Laplace sphere of influence. Pass if `turn_feasible`. The magnitude
  mismatch is reported (it measures how far from a two-body hyperbola the encounter is). The gate is
  also applied with the osculating asymptotes at periapsis, which only restates V4(a) and is NOT
  counted as independent evidence.
* **V6 no unintended encounter**: no approach to either moon inside 2 Hill radii (Titania 20,545 km,
  Oberon 26,576 km) other than the two targeted flybys per half cycle; no impact anywhere.
* **V7 sensitivity (reported, not pass/fail)**: from the interior node nearest the apoapsis before
  each flyby, perturbations of 1 km (radial, along-track) and 1 m/s (radial, along-track); the
  change of the next flyby's periapsis distance and time, linear (STM) and nonlinear. A maintenance
  estimate per cycle is given only if it can be made honestly.
* **Identity**: the converged orbit is "the candidate" only if it is reached by continuation from the
  lam = 1e-3 patched-conic solution with the same turn sense at both flybys (periapsis on the same
  side of each moon) and the same revolution count per leg (five apoapsis passages between flybys).
  An orbit found otherwise is reported as a different object.

### 1.7 Fail outcomes (any one is a verdict "does not survive", reported with the step where it fails)

* F1: no convergence at lam = 1 and the continuation folds or its step collapses before lam = 1.
* F2: the orbit at lam = 1 has the opposite turn sense at either flyby or a different revolution count.
* F3: a periapsis below the 50 km floor or below the surface.
* F4: a periapsis outside the Laplace sphere of influence, or moon-centric e <= 1 there.
* F5: an unintended approach inside 2 Hill radii.
* F6: V2 or V3 fails.

Step 5 (a first real-ephemeris look) is run only if V1 to V6 pass.

## 2. Results (run after commit 0e5d63ae, which holds section 1)

Outputs in `data/found/890_titania_oberon_candidate/` (`rebuild.json`, `decisive.json`,
`continuation.jsonl`, `verify.json`, `refine.json`, `sensitivity.json`, `variants.json`,
`realeph.json`). Everything below was run unless marked "inference".

### 2.1 Deviations from the pre-registration (stated before the verdict)

1. **Continuation started at lam = 0.01, not 1e-3.** At lam = 1e-3 the flyby periapses are 5 km
   (Titania) and 3 km (Oberon); one multiple-shooting evaluation with the variational equations
   took 399 s (the backward Oberon segment alone exceeded 60 s; the Titania segment needed 25,000
   right-hand-side calls), which does not fit the 8-minute command limit. At lam = 0.01 the patched
   conic is still an excellent approximation (r_p / r_Hill = 0.024 at Titania).
2. **Mass-scaling predictor added.** A first lam = 0.01 -> 0.02 step from the unchanged lam = 0.01
   solution started at residual 0.8 (a flyby at the old periapsis distance turns far more under
   twice the mass). The step was killed and the predictor `rescale_flyby_nodes` (periapsis offset
   proportional to lam, moon-relative velocity unchanged: the patched-conic scaling of section 1.3)
   was added. Every later step then converged in 3 to 5 Newton iterations.
3. **Newton tolerance 1e-10** (the V1 bound) rather than 1e-11: at lam = 0.01 the residual floor
   of the integrations is about 3e-11 to 2e-10 and the 1e-11 tolerance only cycled.
4. **The identity criterion's "five apoapsis passages" was a miscount in section 1.6.** The
   patched-conic leg itself has SIX apoapsis passages (5 full revolutions from just before
   periapsis at Titania to just past apoapsis at Oberon); counted the same way, the lam = 0.01
   and lam = 1 orbits have 6 per leg. The criterion's intent (same revolution count as the conic)
   is met.
5. **Stages `refine` (post hoc) and `realeph`** were added after `verify`; see 2.5 and 2.7.

### 2.2 Rebuild and frame checks (pass)

The stored record is reproduced exactly: V-infinity 0.274732 / 0.268904 / 0.274732 km/s
(difference 0), demanded turns 58.44 deg (Oberon) and 42.16 deg (Titania), required altitudes
2,215 and 4,562 km; leg conic a = 511,105 km, e = 0.14849, r_p = 435,210 km, r_a = 587,001 km,
period 11.039 d. Frame round trip exact to 9e-16 km/s; at lam = 0 the CCR4BP flow reproduces the
Kepler leg to 6e-12 (position, relative) and 3e-12 (velocity) and ends 3.5e-7 km from Oberon; the
model's 2.5 forcing periods equal the stored leg time to 1e-12 at lam = 0; at lam = 1 the forcing
period is 24.6324 d (cycle 123.1622 d) and Oberon is at 180.000 deg at tau = T. The mirror
symmetry holds to 5e-13 with theta0 = 0 and fails (3e-4) with theta0 = 0.3 (control).

### 2.3 First decisive test: the patched-conic state does not stay on the cycle (reported)

(a) The literal stored conic, started at its leg-0 apoapsis nearest mid-leg (t = 27.85 d), passes
Oberon at 4,259 km (0.32 Hill radii, t = 60.83 d), then has no approach inside 2 Hill radii of
either moon for the rest of the cycle; after one cycle it is 1,045,000 km and 3.3 km/s from its
start. (b) The honest realisation (the Titania hyperbola periapsis of the patched-conic turn,
r_p 5,351 km, at tau = 0) never meets Oberon at tau = T: at T it is 967,000 km from the intended
Oberon periapsis; its only later close approach is an unplanned Oberon pass at 155 km from the
centre (inside the moon) at t = 99.4 d. Why (run, `diag`): at lam = 1 the trajectory leaving
the periapsis is already 11,000 km and 98 m/s off the leg conic after one day and 46,000 km and
335 m/s after three, scaling linearly with lam (1.1, 10, 98 m/s at one day for lam 0.01, 0.1,
1). The flyby's 5,000-8,000 km offset from the moon's centre changes the Uranus-centred energy
by a few percent, and over 5.5 revolutions that is a full phase slip. So the patched-conic
state is not on, or near, a closing trajectory in the four-body model; whether a nearby orbit
exists is step 3's question.

### 2.4 Correction (converged)

Symmetric multiple shooting, 11 interior nodes. Positive control at lam = 0.01, from the
patched-conic guess (initial residual 2.2e-3): converged in 14 iterations to 4.7e-12; flyby
periapses 51.9 km at Titania (patched conic 53.5) with osculating turn 42.96 deg (42.16), and
29.6 km at Oberon (29.8) with 58.60 deg (58.44), on the same sides. Continuation 0.01, 0.02,
0.05, 0.1 to 1.0 in steps of 0.1: every step converged (final residuals 1e-11 to 3e-11), the sign
of det J stayed -1 throughout (no fold or branch point detected between ladder points), the
periapsis stayed OUTSIDE Titania (on the far side from Uranus) and INSIDE Oberon (Uranus side),
turn sense unchanged. The periapsis distance grows sub-linearly with lam (Titania 52, 101, 240,
448, 811 ... 2,766 km; Oberon 30, 59, 145, 285, 549 ... 2,126 km) and the turn grows (Titania 43.0
-> 66.2 deg, Oberon 58.6 -> 70.1 deg, osculating at periapsis).

**The lam = 1 orbit** (`verify.json`): Titania flyby periapsis 2,765.8 km from the centre
(altitude 1,977 km, 0.27 Hill radii), moon-relative speed 0.482 km/s, osculating e = 1.83,
osculating V-infinity 0.261 km/s, turn 66.2 deg; Oberon flyby 2,125.6 km (altitude 1,364 km,
0.16 Hill radii), 0.515 km/s, e = 1.74, turn 70.1 deg. V-infinity measured at the Laplace sphere
crossings is larger (0.351 km/s Titania, 0.333 km/s Oberon), and each sphere crossing takes
about a day: the encounters are not two-body hyperbolae, as expected at this V-infinity.
Monodromy eigenvalues 8.39e5, 1.789, 0.559, 1.19e-6 (det 0.9967; the departure from 1 is the
integration error of a matrix with an 8e5 entry): two real reciprocal pairs, i.e. strongly
unstable (factor 8.4e5 per 123-day cycle) with a second, mildly hyperbolic pair.

### 2.5 Verification against section 1.6

| criterion | result | numbers |
|---|---|---|
| V1 residual <= 1e-10 | PASS | 2.6e-11 |
| V2 DOP853 closure <= 1 km and 1 cm/s | **FAIL as pre-registered** | 0.238 km, **1.64 cm/s**; Lambda = 8.39e5 < 1e6, so the fallback clause did not apply |
| V3 Radau closure, and DOP853-Radau agreement at T | **FAIL as pre-registered** (closure); agreement passes | 0.230 km, 1.59 cm/s; agreement at T 9e-6 km, 9e-10 km/s |
| V4 periapsis >= 50 km altitude, inside the Laplace sphere, e > 1 | PASS | Titania 1,977 km alt, 2,766 km < 7,532 km, e 1.83; Oberon 1,364 km alt, 2,126 km < 9,678 km, e 1.74 |
| V5 gate on Laplace-sphere-crossing V-infinity, 50 km floor | PASS | Titania 55.7 / 86.8 deg (ratio 0.64), Oberon 66.0 / 88.1 deg (0.75); in/out magnitudes equal by symmetry |
| V6 no other approach inside 2 Hill radii, no impact | PASS | none; distance to Uranus 435,000 to 588,000 km |
| identity | PASS | continuation from lam = 0.01, same sides and turn sense, 6 apoapses per leg as in the conic |

On V2/V3 (POST HOC, labelled as such; `refine.json`): the closure error is what the instability
does to the corrector's residual (8.4e5 x 2.6e-11 = 2.2e-5 nondimensional is the scale; the
observed error is 4.5e-6), and the half-cycle symmetry conditions that define the orbit are met
to 0.27 m and 2.6e-8 km/s. One Newton step on the full-cycle map with the monodromy moves the
start by less than a millimetre and the one-cycle closure becomes 5.2 m and 0.036 cm/s (DOP853)
and 0.68 m and 0.005 cm/s (Radau), both inside the pre-registered bounds. My pre-registered
fallback clause was conditioned on Lambda > 1e6 where it should have been conditioned on
Lambda times the residual; as written, V2 and V3 fail, and the pass comes only after a post hoc
refinement. I judge the orbit to exist in the model (the post hoc evidence is strong and
independent of the clause), but the reader should know the pre-registered test was not passed as
written.

### 2.6 Sensitivity (V7, reported)

From the apoapsis before each flyby (`sensitivity.json`; nonlinear runs at the full and at a
tenth of the perturbation agree to within 10 percent, so the response is close to linear):

| flyby (apoapsis lead) | 1 km radial | 1 km along-track | 1 m/s radial | 1 m/s along-track |
|---|---|---|---|---|
| Oberon (11.7 d) | -9.1 km, -52 s | +0.3 km, +4 s | -63 km, -7 s | -1,243 km, -2.6 h |
| Titania (5.8 d) | +7.0 km, +33 s | +1.6 km, +11 s | +438 km, +0.7 h | +1,430 km, +1.5 h |

A 1 m/s along-track error moves the next periapsis by more than the Oberon flyby's altitude
margin (1,364 km): the orbit is only flyable with routine targeting, about 1.3 km of periapsis
shift per mm/s of along-track error. A maintenance cost per cycle was NOT estimated (it needs a
navigation-error model and a targeting scheme; inference: with centimetre-per-second
navigation the per-flyby corrections are of order cm/s, as in any multi-flyby tour, but that is
not computed here).

### 2.7 First real-ephemeris look (run, uncorrected; `realeph.json`)

Model: Uranus point mass (system GM minus the five moons), J2 = 3.34343e-3 about Titania's orbit
normal (not the IAU pole: a first-look approximation), and Miranda, Ariel, Umbriel, Titania and
Oberon as point masses from URA111 (with indirect terms). The refined lam = 1 orbit's
Titania-relative periapsis state was placed at Titania at three Titania-Oberon conjunctions
(2030-01-12, 2030-02-05, 2030-03-02), in Titania's instantaneous orbit plane, and propagated
for one cycle without any correction. All three reach Oberon about 0.5 to 0.7 d early, at
1,986 km, 842 km and 382 km from its centre (the last is inside the moon, radius 761 km: an
impact; distances from a 266-s sampling); none returns to Titania (690,000 to 792,000 km away at
the end of the cycle). This is the expected outcome for an orbit with 1.3 km of flyby shift per
mm/s: Titania's real eccentricity alone (0.0011) changes its orbital speed by about 4 m/s.
**What a real-ephemeris correction would need** (not attempted): multiple shooting in the full
force model with nodes at every apoapsis and the flyby periapses free, over the validity window
(several cycles, not one, because the real system is not periodic); the IAU Uranus pole for J2
(and J4); out-of-plane states (inclinations of 0.08 and 0.06 deg give hundreds of km at the
flybys); a deterministic manoeuvre budget per cycle as the measure of whether "ballistic"
survives; and the same encounter criteria (V4-V6) at every flyby.

### 2.8 Variants (step 6, from `data/found/888_turn_gate/extended_uranus.jsonl`)

Not needed for the verdict, recorded for completeness. Among all stored Titania-Oberon symmetric
closures, the next-best worst-case gate ratios are 2.11 (n = 4, i.e. two synodic periods per
leg, 4 revolutions, high branch), 2.53 (n = 3, 3 revolutions, low branch) and 3.46 (n = 2): all
fail. Legs of 3.5 and 4.5 synodic periods (86 and 111 d) are beyond the enumeration's
6 sqrt(P_T P_O) = 65 d limit and were not enumerated.

### 2.9 Verdict

In the planar circular restricted four-body problem with both Titania and Oberon at their
physical masses, there IS a symmetric periodic orbit (period five synodic periods, 123.16 d)
that continues the patched-conic closure without a fold, encounters Titania at 1,977 km altitude
and Oberon at 1,364 km altitude once each per half cycle (both inside the Laplace sphere,
hyperbolic, turn-feasible), and has no other close approach. It is strongly unstable (8.4e5 per
cycle). The pre-registered closure tests V2 and V3 were failed as written (1.6 cm/s against
1 cm/s) and passed only after a post hoc refinement. It is not a trajectory in the real system
without correction: uncorrected real-ephemeris propagations lose the return to Titania in one
cycle (and one hits Oberon). The patched-conic state itself is not near a trajectory (section
2.3); the orbit found differs from it by thousands of km and about 20 degrees of turn at each
flyby. Status: a model-level periodic orbit, not a catalogue row. Still open: the
real-ephemeris correction and its manoeuvre budget, the spec 16.4 literature check (Canales,
Howell and Fantino treat a one-way Titania to Oberon transfer; near-Hohmann Titania-Oberon
cycling is an obvious construction, so novelty is doubtful), and an adversarial review.
