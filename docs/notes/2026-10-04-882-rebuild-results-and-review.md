# #882 — the rebuilt four-body torus-connection machinery: results, checks, and the review's verdict

**Date:** 2026-10-04. **Module:** `src/cyclerfinder/search/ccr4bp_strob_connection.py` (tests
`tests/search/test_ccr4bp_strob_connection.py`, driver `scripts/screen_882_ccr4bp_strob_connection.py`,
outputs `data/found/882_ccr4bp_strob_connection/`; commits `72c6c92f` to `f3425a87`). **Why it was
built:** `2026-10-03-882-umbriel-torus-row-adversarial-review.md`.

**Bottom line.** The rebuilt machinery produces real model objects: one-trajectory connections
that leave a torus and return to it, phase-consistent by construction. They were checked three
ways (the build's own verification, the coordinator's independent re-integration, an adversarial
review with its own computations). **Nothing is catalogued.** The Uranian object is on HOLD: it is
a genuine dynamical object of the three-moon-free model, but it encounters no moon, sits outside
the catalogue's classes, and the lane still has no positive control against a published
connection.

## 1. Method

Everything is done in the one-forcing-period (stroboscopic) map of the planar concentric circular
restricted four-body problem. A torus is an invariant circle of the map, corrected by a GMOS
scheme and continued in the perturber's mass at fixed rotation number. Stable and unstable
bundles come from the linearised circle map. Manifold points are parameterised by the circle
angle and an offset in a fundamental domain, with a second-order term. A connection is a zero of
a 4-vector in 4 unknowns at one forcing phase. A verification step integrates one trajectory
through all periods.

Design corrections the build agent had to make (the coordinator's design was wrong or incomplete
on each): the fundamental-domain identity had the shift on the wrong side; unperturbed
connections are not isolated (the autonomous problem has a continuous family), so the refinement
is rank-deficient there and needs Levenberg-Marquardt; a trivial zero exists at zero offset;
first-order departures leave a floor of order offset squared, hence the second-order term;
"continuing the unperturbed connection" has to mean continuing the zeros of a Melnikov-type
mismatch function selected from the unperturbed family.

## 2. Results as built

- **R1.** A homoclinic connection of the unstable Jupiter-Europa 3:4 resonant orbit with the
  perturber's mass set to zero (junction counts (18, 10); its time-reversal partner also found).
- **R2.** Connections continued to Ganymede's physical mass: branches 0, 2 and 4. Branch 3 folds
  near 0.17 of physical mass. Branches 1 and 5 never refined.
- **R3a.** The Uranus-Umbriel 1:2 orbit (the parent of the withdrawn row) under Titania: the
  corrector fails at a thousandth of Titania's mass.
- **R3b.** A homoclinic connection of the Umbriel 3:4 torus (seed eccentricity about 0.1,
  multiplier 1.636) at Titania's physical mass, junction counts (20, 19).
- The old `#694` "positive control" seed was an ELLIPTIC orbit (multipliers 0.974 +/- 0.228i)
  that crosses Ganymede's orbit: there was no hyperbolic object for its "manifolds" to belong to.

## 3. Coordinator's independent check

Own integration with the core propagator, own circle interpolation, own distance measure.

- Circles invariant under the core propagator to 5e-13 (Uranus) and 2e-12 (Jupiter) on sample
  nodes.
- Re-integrating each stored junction state backward and forward reproduces the stored history.
  Uranus: starts 6.7e-4 from the torus, reaches 0.452, returns to 3.0e-4 at the junction count and
  2.1e-4 four periods later; never nearer than 0.092 to Umbriel or 0.30 to Titania. Jupiter
  (three branches): start 3e-4 to 8e-4, reach 0.24, return to 2e-5 to 6e-5.
- The 12 gate tests pass, including one that rejects a deliberately phase-inconsistent
  connection.

## 4. Adversarial review (independent agent, own computations)

| Claim | Verdict | Notes |
|---|---|---|
| R1 | ESTABLISHED (numerical) | Re-refining at smaller offsets moves the junction by 5e-8 to 2.4e-7. The Jacobi agreement is by construction and certifies nothing. Not in the test file; never compared with a published connection. |
| R2 | ESTABLISHED, but TWO independent trajectories, not three | Branches 0 and 4 are mirror images of each other (7.5e-8 apart after the symmetry); branch 2 is self-symmetric and distinct. Branch 0 is near a fold: its smallest singular value falls from 4.1e-4 (half mass) to 1.56e-4 (physical) and refinement fails at 1.10 of physical mass. |
| R3a | NOT ESTABLISHED as a non-existence claim | Non-convergence is not diagnostic (the corrector also fails on 5 of 7 hyperbolic orbits that stay inside Titania's orbit). Supporting evidence added: the unperturbed torus passes 0.00285 from Titania's centre, inside its radius of 0.00297, and the residual and Fourier tail scale with mass, the signature of a singular forcing. Wording to use: "no smooth continuation found; the torus contains collision phases". |
| R3b | ESTABLISHED as a model object; stored precision overstated | In bundle coordinates the unstable coordinate shrinks backward at 0.611 to 0.612 per period for 14 periods and the stable one shrinks forward at 0.608 to 0.612 (multiplier 0.6114). Random 1e-6 perturbations at departure miss by 5e-4 to 7e-3: the return is specific to about 1e-8. Six re-refinements at offsets from 1e-3 to 1e-5 agree to about 2e-7. The STORED junction is off by 2.6e-6 (the refinement left the fundamental domain). Continues in mass from 0.9 to 1.5 of Titania's. Not self-symmetric, so its mirror image is a second connection. |

**The most serious weakness.** The quoted junction residual of 1e-10 is the residual of the
offset MODEL, not the accuracy of the trajectory, and the verification's "independent" target is
computed by the same parameterisation. Real accuracy is 1e-7 to 3e-6. And there is still no
positive control against a published connection.

**Other findings.** The unexplained stored diagnostics (`second_order_residual`) are the cubic
coefficients of the manifold expansion, not defects. The eigenvalue selection is right, but the
spectrum contains spurious reciprocal pairs, so the product of the two multipliers being one
discriminates nothing. The changed verification criterion (four extra periods) was a fair
response to an unsatisfiable original, but it is weak: it accepts trajectories perturbed by 1e-8
to 1e-6. The review proposes a pre-registered replacement: (1) the stable coordinate contracting
within 2 percent of the multiplier for at least 8 periods, unstable-to-stable ratio at arrival
below 1e-3, and the mirror test at departure; (2) junction agreement between an offset and a
tenth of it within 1e-5; (3) random 1e-6 perturbations must fail. The stored results pass these;
the stored junctions should be regenerated at an offset of 1e-4 inside the fundamental domain.

## 5. Literature position (review; first four quotes checked against the held papers)

- Kumar, Anderson and de la Llave, Acta Astronautica 2023 (held as the accepted manuscript): the
  Jupiter-Europa-Ganymede four-body problem at physical masses, with Jupiter-Europa 3:4 tori and
  their whiskers. So R2's TORUS and MANIFOLDS are published objects. It reports only mesh
  near-intersections between Jupiter-Ganymede and Jupiter-Europa tori: "it was not possible to
  differentially correct ... to yield an exact, zero delta-v solution"; "The ultimate goal is to
  find true heteroclinic connections in the CCR4BP". No homoclinic.
- The same authors, SIAM J. Applied Dynamical Systems 2025: the planar ELLIPTIC three-body
  problem. It prints a refined 3:4 to 5:6 torus heteroclinic with coordinates
  (x, y, px, py) = (-0.96064, 0.88783, -0.51377, -0.64714), rotation data 1.558039 and 1.030011.
  **This is a published connection with printed numbers: the positive control the lane lacks.**
- Kumar and Anderson, AAS 24-288: "we do not yet look at the heteroclinics ... in the CCR4BP".
- "First refined four-body torus homoclinic" is defensible wording, but the object is the
  expected persistence of a known three-body homoclinic, around a torus that is itself published
  at Jupiter.

## 6. Decisions

- **No catalogue row.** The Uranian object (R3b) is recorded here only. Reasons: it encounters no
  moon (closest to Umbriel 0.092 of Umbriel's orbit radius, about 5.4 Hill radii; 0.30 to
  Titania); its osculating semi-major axis stays in 1.210 to 1.219 and eccentricity in 0.096 to
  0.104, so it is a resonant-phase excursion of one ellipse and belongs to none of the
  catalogue's classes; the model omits Ariel and Oberon, which the review estimates pull about
  25 and 10 percent as hard as Titania, against a transversality that is itself proportional to
  the perturber's mass; and the lane has no published positive control.
- **R3a wording** as in the table above.
- **Next (registered as `#889`):** reproduce the SIAM 2025 printed connection with an elliptic
  right-hand side in this module. Then regenerate the stored junctions at a 1e-4 offset, adopt
  the pre-registered verification criteria, and add R1 and an offset-convergence gate to the
  tests. Only after that does the rest of `#886` (connections at Titan-Rhea and at Uranus) go
  ahead, and only objects that actually encounter a moon are candidates for the catalogue.
