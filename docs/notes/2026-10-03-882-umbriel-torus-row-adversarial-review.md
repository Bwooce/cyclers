# #882 — adversarial review of `umbriel-1-2-torus-homoclinic-uranus-2026`: not labelled; the catalogued connection is not a trajectory of its own model

**Date:** 2026-10-03. **Trigger:** the owner asked whether the row is novel and, if so, to label it
(`#865` residual item 1). An independent adversarial review was run on the most capable model, told
to look for reasons the row is NOT novel. The coordinating session then checked the review's core
finding directly against the code and the stored result. Nothing was computed for this note.

**Outcome: `our_status` stays absent. The coordinator's same-day proposal of `candidate-novel`
(recorded in the `#865` ledger bullet) is withdrawn.** The literature is clear, but the object as
catalogued does not satisfy the equations of the model it is stated in, so there is nothing yet to
attach a novelty label to.

## 1. The defect (verified by the coordinator)

The concentric circular restricted four-body problem is time-dependent: the perturber's angle is
`theta = theta0 + omega * t` (`src/cyclerfinder/core/ccr4bp.py`, the line computing `theta` from
`system.omega_gan * t`). A connection is a single trajectory, so where the outgoing and incoming
arcs are joined the SPACECRAFT STATE and the PERTURBER PHASE must both agree.

- `search/ccr4bp_manifold_globalize.py::manifold_state_at` starts BOTH branches at the same
  absolute time `t0 = theta1_section / omega1` and integrates `+t_flow` for the unstable branch and
  `-t_flow` for the stable one. The junction is therefore reached at time `t0 + t_u` on one arc
  and `t0 - t_s` on the other.
- `search/ccr4bp_heteroclinic_search.py` refines on a residual of four numbers only,
  `(x, y, vx, vy)` of one arc minus the other. The forcing phase is not in the residual and is not
  constrained anywhere else.
- For the catalogued solution (`data/found/701_ccr4bp_umbriel_titania_search/result.json`,
  `best_robust_genuine_connection_corrected`): `t_u = 19.002859`, `t_s = 19.092875`, sum
  `38.095733` TU; forcing period `torus_period_tu = 11.991105` TU; so the two arcs' junction times
  differ by **3.17700 forcing periods**. Titania is **63.7 degrees** apart on the two sides of the
  junction.

The "closure to 5.66e-10 km / 5.1e-14 km/s" is a match of the spacecraft state between two arcs
that live at different perturber phases. Joined, they are not a solution of the CCR4BP. A correct
residual needs the stable arc to arrive on the torus at section phase
`theta1_section + omega1 * (t_u + t_s)` (or, equivalently, a stroboscopic-map formulation in which
both arcs are compared at the same phase, as Kumar, Anderson and de la Llave do).

Supporting evidence in the same result file (verified): the SAME solution (`t_u` about 19.00,
`t_s` about 19.09, `theta2_u` about 3.045, `theta2_s` about 3.237, and its mirror) is returned for
all four combinations of the unstable and stable lobe signs. A real manifold intersection depends
on which side of the torus the arcs leave from. Insensitivity to the lobe sign is what a
self-intersection of the torus's own projection onto `(x, y, vx, vy)` looks like: two torus points
with the same state at different forcing phases. That reading is the reviewer's inference; it was
not tested by computation.

## 2. Further findings by the reviewer

Each is marked with what the coordinator checked.

- **The parent orbit crosses Titania's orbit** (verified: the result file's own `base_orbit_note`
  gives the extent as 1.477-1.698 Umbriel units "bracketing Titania's own a_gan=1.640"). Kumar and
  Anderson (AAS 24-288, p. 17; quote verified in the held paper), on Oberon orbits that intersect
  Titania's orbit: "This immediately precludes the existence of quasiperiodic dynamical
  equivalents of those 4:3 orbits, as they would experience an infinite perturbation from the
  singularity at Titania for any mu3 > 0." The stored torus has `torus_closure_residual` 1.4e-4
  (verified), about 140 times the 1e-6 manifold offset used on top of it. So the torus itself is
  an approximate object, and by the published argument an exact one cannot exist for this parent.
- **The real-ephemeris evidence (`#704`/`#705`) does not test a connection** (checked against the
  module docstring only): it flies the unstable half from the model's departure state and compares
  the result with the model's own target state. The return arc is never flown. "10 of 10 epochs
  recur" shows that the outgoing arc is shadowed at favourable phases, not that a connection
  exists.
- **The Jupiter-Europa-Ganymede "positive control" (`#694`) used the same residual** (not checked
  by the coordinator). If so, it did not validate the method against a published connection in
  the sense intended.
- **V1 does not cover the gap** (agreed): the Radau-against-DOP853 agreement checks each
  integration, not the missing phase constraint.
- **Row wording that overclaims** (agreed): "NEW dynamical species", "closing ... to machine
  precision", "real, recurring, correctable near-connection".

## 3. Literature (no collision found; policy fit is the owner's call)

- Kumar and Anderson, AAS 24-288, is the only located four-body work at Uranus. It treats Oberon
  resonant orbits with Titania as perturber, and says (p. 12, verified) "we do not yet look at the
  heteroclinics between mean motion resonances in the CCR4BP" and (p. 19, verified) "there are
  three other large moons of Uranus - Titania, Umbriel, and Ariel - for which this study should be
  repeated".
- No Umbriel 1:2 orbit, torus or connection was found in any held paper or in eleven live queries
  by the reviewer, or in the coordinator's two web queries and check of the first author's
  publication list earlier the same day.
- Pergola, Geurts, Casaregola and Andrenucci, IEPC-2007-305 (found by the reviewer; now held)
  computes three-body invariant manifolds for each Uranus-moon system including Umbriel, linked by
  electric propulsion in a one-way five-moon tour. It is not a collision (no resonant orbit, no
  torus, no two-moon model by a first text search; digest in progress) but it is prior
  three-body work at Uranus-Umbriel and must be an anchor.
- The 2026 ice-giant review (Liang et al.) does not cite AAS 24-288 at all, so it cannot be used
  as evidence that no dynamical-systems work exists at Uranian moons.
- **Policy fit if the object is ever rebuilt.** Spec 16.4 case (ii) reads "a primary/moon system
  the authors never treated (Uranus, Neptune)". Its examples are primaries, and Kumar and Anderson
  did treat Uranus with this architecture and this perturber. Whether an untreated BASE MOON in a
  treated primary is "a system never treated" or "a member they did not happen to tabulate" is not
  settled by the text. It is the same question as the unenumerated Io-Callisto pair, and the
  coordinator's proposal answered it at pair level without saying so. The owner must rule.

## 4. Errors found in existing project material

- The row's notes and the corpus anchor "Kumar Uranus-Oberon PCRTBP MMR study (2025)" credit
  arXiv:2509.03655 with the Uranus-Titania-Oberon four-body secondary resonances. That paper does
  not mention Titania (verified: zero occurrences in the held file); the content is in AAS 24-288.
- The coordinator's proposal said the connection is "exact only in the idealized model". It is not
  exact there either.
- `#699`'s statement that no Uranian moon-pair paper is held is stale (AAS 24-288 and the Kumar
  paper have been held since `#728`).
- The reviewer also reports that the `#728` digest compares a stroboscopic rotation number of
  5.9956 radians with the integer 6 (it is 0.9542 cycles); not checked by the coordinator.

## 5. What follows

1. The row is not labelled. A notice has been added at the top of its notes; its numerical fields
   are unchanged pending the owner's decision on withdrawal (remove the row until rebuilt) against
   keeping it with the notice and demoting `validation_level`.
2. The website's "Discovered by this project" strip now requires a novelty label, so a discovered
   row with no `our_status` is no longer presented as a discovery.
3. `#882` holds the rebuild: a phase-consistent residual; a positive control that reproduces a
   PUBLISHED four-body connection with the phase constraint in force; a parent orbit that does not
   cross the perturber's orbit; and a real-ephemeris check that flies both arcs.
4. The four torus-connection negatives stamped under `#865` (`#695`, `#696`, `#703`, `#716`) used
   the same search. A true connection needs the four-component state match AND the phase match,
   so the absence of a state match is, on its face, still a valid negative; that argument has not
   been checked against how each search measured its closest approach, and the stamps should say
   so.

## 6. Addendum, same day: the defect is in the whole lane, and it is deeper than the missing phase

Further checks by the coordinator before starting the rebuild.

- **Every stored result has the phase defect.** From the result files: `#694`
  (Jupiter-Europa-Ganymede, the lane's "positive control") joins its arcs 316.8 degrees apart in
  Ganymede's phase (`t_u + t_s` = 1.8799 periods); `#701` (Umbriel) 63.7 degrees; `#703`'s
  closest candidate 286.5 degrees. `#694` was not a comparison with a published connection: it
  was the same search run on the Europa 3:4 torus and judged by the same guards. The lane has
  never had a positive control.
- **The departure offset was below the torus's own error.** The tori were computed with one
  Fourier harmonic in the forcing angle (`n1=1`) and have invariance residuals of 1.1e-4 to
  1.8e-4 for `#694`, `#695` and `#701` (7.9e-7 for `#703`, 2.2e-8 for `#716`), while the manifolds
  were started 1e-6 off the torus. Where the torus error exceeds the offset, the trajectories are
  governed by the torus error and not by the stable or unstable direction. That is why `#701`
  returns the same solution for both signs of the offset. The module docstring argues that an
  offset "two orders of magnitude BELOW that floor" is safe; it is the other way round.
- **The dimension count in `ccr4bp_manifold_globalize.py` is wrong.** In the extended phase space
  `(x, y, vx, vy, forcing phase)` the manifolds of a 2-torus are three-dimensional and meet
  generically along isolated trajectories. Fixing the departure phase and the offset size leaves a
  two-dimensional slice of each, and two such slices generically contain no connection at all.
  The search could not have found one even with the phase constraint added.
- **Consequence for the four negatives.** Section 5 item 4 above said the absence of a state
  match "is, on its face, still a valid negative". That was wrong: the searches looked in slices
  that generically miss every connection. The four stamps (`#695`, `#696`, `#703`, `#716`) are now
  marked METHOD-INVALID in `data/empty_regions.jsonl`, with the original text kept.
- **The rebuild** (`#882`) works in the stroboscopic map: the torus is an invariant circle of the
  one-period map, corrected to about 1e-10; the manifolds are two-dimensional in the
  four-dimensional map space and meet in points; both sides of a junction are at the same forcing
  phase by construction; and a claimed connection must pass a single-trajectory test (one
  integration from the departure point that leaves the torus and returns to it). The positive
  control is the classical homoclinic of the unperturbed resonant orbit, continued in the
  perturber's mass to its physical value.
