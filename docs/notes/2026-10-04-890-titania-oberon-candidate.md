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
