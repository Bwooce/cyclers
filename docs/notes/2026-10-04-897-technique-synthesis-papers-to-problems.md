# `#897` Technique synthesis: what the methods in the newly filed papers can do for the project's problems

Date: 2026-10-04 (Sydney). Author: the `#897` synthesis agent. Status: first complete version. Reading and reasoning
only; no code was run beyond catalogue counts and text searches; no existing file was edited.

The owner's criticism this note answers: the thirty-odd papers filed on 2026-10-03 and 2026-10-04 were used to check
constants and models. This note is about their METHODS: what each does, under what stated hypotheses, what the project
already has, and what could be built, in an order that favours results that would survive an adversarial review.

## How to read the evidence tags

Every statement about a paper carries one of these tags.

- `[P p.N]` READ by me in the paper itself (page image or text layer), page N.
- `[R p.N]` READ in the paper by one of six source readers I dispatched for this task (Opus agents with page-cited
  briefs; their reports are in the session scratch area, not in the repository). I did not re-read these pages. The
  readers' own inferences are tagged `[I]` here, not `[R]`.
- `[D name]` READ in a project digest (file under `docs/notes/`), not checked at the source by me.
- `[I]` my inference. `[C]` a small calculation of mine.

Papers are cited as "filed in the private paper corpus as <filename>" where a filename is needed. Nothing in this note
is described as new; where the held papers do not show a result, the wording is "not found in the held papers".

Pages I read at the source myself: Casoliva et al. 2010, pp. 1627-1630 (the section on second-species seeds and
continuation in the mass ratio); Russell & Ocampo 2006, the abstract, the flyby-model paragraph, the results paragraph
and the conclusions (text layer); the project's own notes on `#388`. Everything else is `[R]` or `[D]`.

## 0. The findings that change plans

1. **For one small body, the "second species as a generator" idea is already published at the Earth-Moon mass ratio,
   in a paper the project catalogues from.** Casoliva et al. (JGCD 33(5), 2010) seed second-species p-q resonant cyclers
   at mass ratio 1e-6 from Barrabes & Gomez's matched in/out maps and continue them in the mass ratio to 0.0121529529,
   with a three-step strategy (mass, then Jacobi constant to lift the periselene, then mass again) `[P pp.1627-1630]`.
   Their Table 3 is nine of the project's catalogued rows. The project has never built this strategy; `#780`(d) tried
   only the other one (two-body seed corrected directly at the physical mass), which the authors say fails for
   semi-major axes below about 0.6 lunar distances `[P p.1627]`. So `#899` for one moon is a reproduction with a
   strong positive control at both ends, not an extension. (Sent to the coordinator when found.)
2. **The second-species generator cannot produce the C11, C21, C31, C32 families.** A collision arc at zero mass has
   relative speed squared 3 - C at the secondary, so it needs C < 3 (Font, Nunes & Simo 2002: "3 - C_J > 0 bounded away
   from 0" `[R p.118]`; Bolotin & MacKay 2000: C in (-sqrt 8, 3) `[R p.56]`). The Ross & Roberts-Tsoukkas and Braik &
   Ross planar rows sit at C = 3.13 to 3.18 `[D gomez-olle]`, and their authors build them from L1 Lyapunov tube
   intersections, not from Kepler arcs `[R Ross & Roberts-Tsoukkas AAS 25-621 p.10]`. The check set for a
   second-species generator is the Casoliva rows.
3. **Bradley & Russell's continuation does not address the step where `#388` failed.** Its input is "An initial ZSOI
   trajectory in some ephemeris model" `[R Algorithm 1, PDF p.4]`. `#388` failed while carrying the circular-coplanar
   parent to the ephemeris, which is Russell & Ocampo's (2006) problem. Their published outcome is that most parents
   are NOT ballistic in the real ephemeris: "There are 9, 39, and 74 parent cyclers that have at least one launch date
   that requires a total delta-v of less than 1, 10, and 300 m/s" out of 203 `[P results paragraph]`. The right target
   for the heliocentric rows is the published real-ephemeris cost, not a ballistic closure.
4. **Seventy-seven published real-ephemeris Earth-Mars cycler solutions with full reproduction data exist in a held
   source, and the project has used three.** Russell's 2004 dissertation, Appendix C, prints "the only values required
   to reproduce the solution" for the lowest-cost seven-cycle solution of 77 parents `[D 2026-06-07-russell-2004-member-
   tables-transcription]`; `search/appc_corrected.py` reconstructs any such block; only parents 83, 188 and 192 have
   been run (`#170`).
5. **Leiva & Briozzo's "periodic arcs" are a numerical fallback, by their own words, not evidence that periodic
   orbits do not exist.** They single-shot over five or more Sun periods in double precision and report that the
   half-integer members "proved difficult to continue"; where the iteration failed "we reduced the integration interval
   to [0, tau], obtaining - when continuation now succeeded - periodic arcs in the QBCP instead of periodic orbits"
   `[P p.239, text layer]`. The `#884` reviewer's finding of true periodic orbits for
   the same members is therefore a specific, checkable extension of a published result, with a published selection
   rule to pre-register (section 4, idea iii).
6. **Centre-manifold reduction does not give bounded motion around an unstable planar cycler** (it has no centre
   directions in the plane). What the Jorba school does give is the correct object for a cycler whose period is not
   commensurate with the Sun: an invariant curve of the stroboscopic map `[P Rosales et al. 2021 p.5, text layer; R pp.13-14]`. For
   STABLE cyclers this is a direct test of a conjecture printed by Ross & Roberts-Tsoukkas (section 5.3).
7. **A technique the candidate list did not have:** Russell & Strange's use of repeated free returns to ONE moon to
   rotate the excess velocity in steps. Every symmetric two-moon closure the project enumerated asked one flyby per
   moon to supply the whole turn, and all but one failed the turn gate (`#888`). Chains with resonant returns split
   the turn over several flybys (section 5.1).

## 1. Technique inventory

Each entry: what the method does; its hypotheses and limits as the source states them; what the project has; what is
missing.

### T1. Continuation from a zero-sphere-of-influence tour to full gravity (Bradley & Russell 2014)

- Does: one scalar kappa in [0, 1] scales the flyby bodies' GM, radius and sphere of influence and blends a mean
  Keplerian ephemeris into the real one; each flyby is re-seeded with the hyperbola that keeps the turning angle
  (r_p = kappa r_p,final); multiple shooting with patch points at periapsis plus or minus a fixed time and at leg
  midpoints; continuity and altitude constraints, no objective `[D bradley-russell, R PDF pp.11-16]`.
- Limits as stated: input is a ballistic ZSOI trajectory already on an ephemeris `[R PDF p.4]`. "For a relatively
  small number of encounters (two to four), the method reliably finds a feasible solution ... the maximum number of
  encounters for a Jovian moon tour is typically between five and nine" `[D, PDF pp.23-24]`. "Large values of Delta
  kappa may lead to a solution in a different family of trajectories" `[D, PDF p.13]`. For the heliocentric example,
  "the algorithm may converge directly when kappa_0 = 1" `[D, PDF p.19]`. Examples: a five-encounter interplanetary
  transfer and a seven-encounter Jovian tour of about 86 days, chained to twelve with a patch manoeuvre of at most
  1 m/s. Not stated anywhere: periodic trajectories, multi-year durations, multi-revolution heliocentric legs,
  low excess speeds `[R]`.
- Project has: `search/titania_oberon_realeph_895.py` (multiple shooting with variational Jacobian, ephemeris blend,
  unsoftened moons, J2 and J4, for Uranus); `search/two_moon_periodic_890.py` (mass continuation with a periodicity
  and symmetry condition).
- Missing: a force-model provider for other systems (Sun with DE440 planets; Saturn and Jupiter with satellite
  kernels: only Voyager-era Saturn and Jupiter kernels are cached locally); altitude inequality handling (the `#895`
  route drifted to high flybys without the declared deviation D1).

### T2. Carrying an idealised cycler to the ephemeris as a constrained optimisation (Russell 2004; Russell & Ocampo 2006; Russell & Strange 2009)

- Does: zero-radius flybys on real planet positions; unknowns per leg are the full outgoing excess-velocity vector and
  the leg's start and end times (5n unknowns); position match at arrival; one conditional flyby cost that is zero when
  magnitudes match and the demanded turn fits, and otherwise the minimum impulse; SNOPT elastic mode minimises the
  summed violation; the planet model is ramped circular, then eccentric, then inclined, then one step to the
  ephemeris; 21 launch windows times six step counts, keep the best; seven cycles `[D 2026-06-07-russell-2004-
  continuation-deepdive; R Russell & Ocampo pp.355-362]`.
- Limits as stated: "All flybys are modeled with the assumption that they occur instantaneously and the radius of the
  planet's sphere of influence is zero" `[P]`. Lambert-leg formulations are rejected: fixing the integer structure is
  "a difficult if not insurmountable obstacle when dealing with a gradient-based optimizer" `[P; R p.356]`. Step-count
  results "are inconsistent" `[R p.362]`. Success is governed by the parent's turn ratio: "the general performance of
  an accurate ephemeris cycler is closely related to the turn ratio of the circular-coplanar parent" `[P]`. The step
  to integrated flybys "is expected to be minor compared to the step demonstrated here" `[P conclusions]`.
- For moons (Russell & Strange 2009): the mean-element ramp "struggles in the case of planetary moons" and is replaced
  by "a linear interpolation between the ideal model and ephemeris locations" `[R p.146]`; patched-conic ephemeris
  costs of 0 to 222 m/s over 5 to 10 cycles; high-fidelity results for single cycles only (32, 63, 11, 121 m/s,
  Table 8) `[R pp.152-155]`; "the transition to a high-fidelity force model is substantially more difficult for the
  planetocentric cyclers" `[R p.152]`.
- Project has: the model ramp (`search/continuation.py`, the same schedule); a reduced corrector
  (`search/correct.py::ballistic_correct`: unknowns t0 and leg times, Lambert legs, magnitude residual or a bend
  hinge); the conditional flyby cost as a function (`core/flyby.py`, `search/dsm_leg.py`); the faithful generator of
  the ideal parents (`search/generic_return.py`, `search/cycler_search.py`); the Appendix C reconstruction
  (`search/appc_corrected.py`); a full n-body shooter with nodes at the bodies (`nbody/shooter.py`).
- Missing: the optimiser itself as published (free excess-velocity vectors per leg, non-Lambert legs, conditional
  cost as the objective, elastic handling of infeasibility, seven-cycle seed, window and step-count scan with
  keep-best). The project's own June note says so: "We were not running Russell's actual corrector"
  (`docs/notes/2026-06-21-narc-continuation-results.md`). The next thing tried was an n-body shooter at physical
  planet masses, not this.

### T3. Matched asymptotics for one close passage (Perko; Gomez & Olle 1991 I; Barrabes & Gomez 2002, 2003)

- Does: joins an outer Kepler arc to an inner hyperbola about the small body. Gives the impact parameter mu K0, the
  periapsis r_p = mu (e - 1)/v_inf^2 with sin(turn/2) = 1/e, and the periapsis time shift mu ln mu/|v|^3
  `[D gomez-olle, Part I pp.123, 136, 146]`. Barrabes & Gomez give explicit initial conditions on a sphere of radius
  mu^alpha about the small body for a p-q resonant orbit, periodic "with an error of the order mu^(1 - alpha)"
  `[R 2003 abstract]`; these are the seed formulas Casoliva reproduce as their eqs. 14-18 `[P p.1628]`.
- Limits: existence for mu below an unquantified mu_1; nondegeneracy checked numerically only for k = 1 to 4
  `[D gomez-olle, p.133]`; "the generating orbits (for mu = 0) associated to p-q resonant orbits are bifurcation orbits
  of 1st species-2nd species" `[R Barrabes & Gomez 2003, about p.170]`, where passages scale as mu^nu with nu < 1 and
  one generating orbit splits into two families `[D gomez-olle, Part II pp.151-152]`. First-order validity needs
  mu |ln mu|/v^3 small and r_p well below mu^(1/2): the `#890` chain is outside it `[D gomez-olle, computed there]`.
- Project has: the gate's periapsis relation, which is this theory's first-order periapsis (`verify/turn_gate.py`)
  `[D gomez-olle sec. 8e]`; a two-body state transition matrix; regularised propagation
  (`core/cr3bp_regularized.py`).
- Missing: the Barrabes & Gomez seed (Casoliva eqs. 14-18) as code; an encounter diagnostic that compares a computed
  flyby with the asymptotic prediction; a digest of the two Barrabes & Gomez papers (held, not digested).

### T4. Shadowing theorems for collision chains and the discrete action (Bolotin & MacKay 2000, 2006; Marco & Niederman 1995; MacKay 2005)

- Does: for an autonomous system with one small centre at fixed energy, any chain from a FINITE set of nondegenerate
  collision arcs whose consecutive velocities are neither parallel nor antiparallel is shadowed, for small enough
  mass, by a unique nearby orbit `[R Bolotin & MacKay 2006 Thm 2.1 p.435; 2000 p.50]`. The orbits are found as
  critical points of a discrete action whose variables are points on small spheres around the centre, two scalars per
  junction in the plane `[R 2000 pp.66-67]`.
- Limits as stated: the mass bound is never quantified `[R]`. The direction change must be taken in the rotating
  frame: "we mistakenly studied the direction change in the inertial frame" `[R 2006 footnote 3 p.444]`; MacKay 2005
  footnote a restates it with one exceptional relation `[R p.3]`. A body of finite size is allowed "provided its radius
  is less than c eps" `[R 2000 p.55]`. The orbits are strongly unstable: "Lyapunov exponents of order log eps^-1"
  `[R 2006 Thm 2.2 p.436]`; "the angle change in a close encounter is a smooth function of impact parameter with slope
  of order 1/eps" `[R MacKay 2005 p.2]`. A second small body on a different orbit is not covered by any of these
  papers `[R]`; the slightly elliptic case is in Bolotin 2005 and 2006, not held `[R 2006 p.434]`.
- A remark that bears on the project's leg solvers: "It is fortunate that we want nondegeneracy for fixed Jacobi
  constant rather than fixed time ... because any point P on a Kepler ellipse is equal-time self-conjugate"
  `[R 2000 p.62]`. So whole-revolution returns are degenerate for a fixed-time (Lambert) formulation and nondegenerate
  at fixed energy `[I]`. This is the same fact as the project's finding that full-revolution and half-revolution legs
  are singular Lambert cases (`#793`, `docs/notes/2026-06-23-388-v1-for-fh-legs-wellposedness-verdict.md`).
- Project has: the turn gate, which is the finite-mass, per-junction form of "radius less than c eps" `[I, R]`.
- Missing: a fixed-energy (fixed excess speed) arc solver for same-body returns; the discrete-action Newton step. The
  first is useful; the second is not needed if T6 is built, since T6 does the same job numerically.

### T5. The arc alphabet and the return map (Henon; Font, Nunes & Simo 2002, 2009)

- Does: at one Jacobi constant, enumerate the Kepler arcs that leave and return to the secondary's position (S arcs
  begin and end at different points of the ellipse, T arcs at the same point) `[R 2009 p.144]`: "41 orbits, 18 resonant
  orbits or T-arcs ... and 23 symmetric orbits or S-arcs" at C = 2.8 with p, q up to 4 `[R 2009 p.154]`. Chains are
  realised through a return map on a circle of radius mu^alpha about the secondary `[R 2002 p.121]`.
- Limits: mu(N) "<< N^-10" `[R 2002 eq. 57]`; periodic orbits computed only at mu = 1e-4; strips followed to 1.5e-3,
  where they are lost by an extra encounter or "strongly deformed" so "the analytical theory cannot be applied (mu is
  too large)" `[R 2002 pp.139-141]`; repeating the same T arc is forbidden; one chain "fails to intersect" at 1e-4 but
  "subsists for 1e-6" `[R 2009 p.160]`; extended precision was used `[R 2002 p.137]`; stability parameters of 2.6e4
  to 2.9e13 `[R]`.
- Project has: the eight printed periodic orbits as tests (`#896`, commit `a18cf2bc`).
- Missing: an arc enumerator (pure Kepler; the 41-arc table is the control); symbol labels on the catalogued rows.

### T6. Small-mass seed plus continuation in mass and energy (Casoliva et al. 2008, 2010)

- Does: for each (p, q), 500 trial Jacobi constants in the admissible interval, seeds from T3 at mu = 1e-6 corrected to
  a periodicity error below 1e-10; then continuation in mu to the Earth-Moon value `[P p.1628]`.
- Limits as stated: fixed-period continuation in mu "in most cases ... leads to a periodic orbit that impacts the moon
  (especially for resonant orbits with high p/q), showing that resonant cyclers with the seed (p, q) do not exist
  uniformly in mu in the family"; fixed-C continuation "does not preserve the resonance relation"; "Instead a
  three-step continuation strategy was necessary" `[P p.1629]`; "small values of mu require high integration accuracy
  and a very small continuation step" `[P p.1628]`; a final continuation in C at the physical mass is needed to find
  where the family crosses the wanted p/q, and not every segment does (four segments for 7-3, three crossings)
  `[P p.1629]`. Stated benefit: "more systematic and thorough in finding resonant orbits for mu = mu_M", and the only
  way to obtain tight cyclers `[P p.1630]`.
- Project has: `search/mu_continuation.py` (pseudo-arclength in (x0, C, mu), turns folds), the symmetric
  fixed-Jacobi corrector, the Table 3 rows verbatim (`search/earth_moon_resonant_families.py`), regularised
  propagation.
- Missing: the seed (T3); the three-step driver with a periselene trigger; the final family walk in C with a
  resonance-crossing detector.

### T7. Which periodic orbits of the three-body problem continue under a periodic forcing, and at which phases (Leiva & Briozzo 2008; Brown et al. 2025, two papers; Oshima 2022; Boudad et al. 2020; Rosales et al. 2021)

- Does: a periodic orbit whose period is p/q forcing periods "becomes a periodic orbit of period pT once eps is set
  to be different from zero. This is a consequence of the Implicit Function Theorem" `[P Rosales 2021 p.5, text layer]`; all
  others "gain, generically, the frequency of the perturbation and become ... two dimensional invariant tori"
  `[R pp.13-14]`. The admissible phases are zeros of a first-order work integral: Leiva & Briozzo's
  Im[e^(-2 i phi) sum ...] = 0, tan 2 phi = Im/Re, four phases a quarter Sun period apart, nonzero only for q = 1 or 2
  ("When q [is not] 1, 2 no condition is obtained at first order ... This course will not be pursued") `[R pp.232-233]`;
  Brown et al.'s Melnikov-type function M, whose zeros come in fours (two integrations suffice by their
  Proposition 2), is identically zero for some symmetric orbits and must then be recomputed at higher order
  `[R SIADS p.351, Props. 2, 3]`. In the bicircular problem the lowest nonzero order is j = a - 1 with zero spacing
  (pi/omega)/(j + 1), an empirical rule "we have not derived an analytical proof", checked only for a of 3 or more
  `[R JAS p.32, p.15]`. For symmetric orbits Oshima fixes the phase by symmetry (Sun angle 0 or pi) and obtains two
  families per orbit, four for doubly symmetric ones `[R p.1328]`.
- Limits: Leiva & Briozzo's condition is necessary, not sufficient `[D leiva-briozzo-2008 p.236]`; their homotopy also
  switches on the primaries' non-circular motion, not the Sun's mass alone `[R eq. 5]`; Brown et al. state no theorem
  and no nondegeneracy hypothesis `[R]`; "some of them may not be truly periodic" for high-order resonances
  `[R JAS p.13]`; "the algorithm does not guarantee the convergence into periodic orbits in the BCR4BP" `[R Oshima
  p.1330]`.
- Project has: `search/sun_forced_periodic_884.py` (member selection, symmetric phases, a Melnikov amplitude, multiple
  shooting, continuation in Sun mass), the corrected models with published controls (`#891`, `#892`, `#896`).
- Missing: the Leiva & Briozzo phase formula for orbits without the symmetry (their "Aa" families); the
  orbit-versus-object deduplication of Brown et al. (same states with a shifted phase are one object `[R SIADS
  p.356]`); an explicit arc-versus-orbit criterion; a separate, labelled class for q = 3.

### T8. Following a branch through folds (Oshima 2022; Brown et al.; Rosales et al. 2021; Casoliva)

- Does: pseudo-arclength in the forcing strength, continued past folds and through negative values of the parameter,
  so that the whole solution curve is known: a branch "encounters a fold point near eps = 1 and eventually returns
  back to eps = 0" onto a different three-body orbit; "period-m orbits can be hopeful candidates" `[R Oshima
  pp.1328-1330]`. A family that stops at a turning point at small forcing does not reach the physical system
  `[R Brown JAS p.23]`.
- Project has: pseudo-arclength in several modules, fold detection by sign of a determinant (`#890`), and the
  tangent-cosine guard the `#884` reviewer used.
- Missing: a policy that every continuation in a model parameter is pseudo-arclength with a branch-jump guard. `#884`
  first ran steps of 0.1 and "every branch reached physical mass"; one had jumped (ledger).

### T9. Dynamical substitutes, Floquet reduction and centre manifolds (Jorba, Jorba-Cusco & Rosales 2020; Rosales et al. 2023; Jorba-Cusco, Farres & Jorba 2018)

- Does: about a periodic orbit of a periodically forced Hamiltonian whose period equals the forcing period, a
  symplectic Floquet change and a partial normal form remove the hyperbolic directions and leave the centre manifold,
  where bounded motion is read from Poincare sections `[R Jorba 2020 pp.6-10]`.
- Limits as stated: "Let us suppose that we want to describe the motion around some periodic orbit of period T of a
  T-periodic time-dependent Hamiltonian system" `[R p.4]`: any such orbit, not only a libration-point substitute. But
  the reduction needs elliptic directions to keep, and small divisors limit the radius: "they usually have a very small
  radius of validity" `[R p.10]`; at L2 in the bicircular problem "the radius of convergence of the computed center
  manifold was very small, and it was concluded that this approach was not suitable" `[R Rosales 2021 p.8]`. No paper
  of the four prices station keeping `[R]`.
- Project has: Floquet multipliers and monodromy; nothing of the Lie-series machinery.
- Missing: everything, and section 4 (idea v) argues it should not be built for cyclers.

### T10. Invariant curves of the stroboscopic map (Rosales et al. 2021; Villegas-Pinto et al. 2023; Kumar, Anderson & de la Llave 2021)

- Does: solves W(theta + rho) = F(W(theta)) by Fourier collocation with a phase condition, continuing in rotation
  number or forcing strength; for very unstable cases, multiple shooting over sub-intervals of the period
  `[R Rosales 2021 pp.9-13]`. Kumar et al. add the frame of tangent, stable and unstable bundles and solve by a
  quasi-Newton step of cost N log N `[R CMDA 2021]`. Villegas-Pinto et al. take a resonant orbit to a periodic orbit
  of ONE intermediate model "by design" and then switch on the second perturbation to obtain a two-torus `[D villegas-
  pinto p.8]`.
- Limits: families are Cantorian, "The gaps in this family are due to resonances"; large gaps are crossed by lowering
  the Sun's mass and returning `[R Rosales 2021 p.13]`. Kumar's method starts from "the stable and unstable unit
  eigenvectors of the periodic orbit monodromy matrix" `[R CMDA sec. 4.9]`, so it needs a real saddle pair; Rosales's
  does not. Tori "passing close to the singularity" need many more nodes `[R Kumar CMDA sec. 4.12]`. No ephemeris
  validation in Villegas-Pinto: "analyses ... in a full-ephemeris model would be necessary" `[D p.18]`.
- Project has: `search/ccr4bp_strob_connection.py` and `search/pertbp_strob_889.py` (invariant-circle corrector,
  bundles, manifolds; lessons recorded under `#886`: refinement ladder, Fourier-diagonal step needed for high-energy
  members); `genome/bcr4bp_torus.py` (frame helpers corrected under `#891`).
- Missing: the corrector pointed at the Sun-forced Earth-Moon models for cycler parents; a multiple-shooting
  ("r-invariant curve") variant for unstable parents.

### T11. High-order manifold expansions and the layered search for torus connections (Kumar, Anderson & de la Llave)

- Does: Taylor or Fourier-Taylor parameterisations of stable and unstable manifolds valid 100 to 1000 times farther
  than the linear approximation; connections as a zero of four equations in four unknowns at one forcing phase;
  layers of fundamental domains make the search complete up to a flight time `[R SIAM secs. 5-7]`.
- Limits as stated: "may actually find mostly intersections which correspond to near-misses rather than true
  intersections ... especially when the periodic perturbation is weak" `[R SIAM printed p.26]`; in the four-body
  Jupiter case the correction to a zero-cost connection failed and the result is small-manoeuvre transfers
  `[D, R]`. Connections are transfers between resonances, not repeating orbits `[I]`.
- Project has: the rebuilt lane with one published control reproduced (`#889`).
- Relevance to cyclers: indirect. A connection is an edge of a network, not a catalogue row, unless it is closed into
  a periodic chain (`#790`, blocked on `#872`).

### T12. Newton for long periodic orbits of a map in a Floquet frame (Kumar 2026, arXiv 2601.00149)

- Does: a periodic orbit with q points of the stroboscopic map is corrected with cost proportional to q, "against a
  4q x 4q linear system" `[R Remark 4]`; periodic orbits were continued to the physical perturber mass in cases where
  the tori were not (Ganymede 4:3 with Europa; Oberon 6:5 with Titania) `[D]`.
- Limits: four-dimensional maps; needs the saddle frame; "does not handle fold bifurcations during numerical
  continuation" `[R]`.
- Relevance: the `#884` and `#890` objects are exactly periodic points of a stroboscopic map with q = 2 to 5. For
  those sizes ordinary multiple shooting is enough; the value is as an independent cross-check at larger q.

### T13. Folds in eccentricity as a diagnostic for ephemeris transition (Park & Howell 2024; Singh, Park & Howell 2026)

- Does: continues p:q resonant L2 halo orbits in the lunar eccentricity by pseudo-arclength and finds that "the
  proliferation of ER3BP PO fold bifurcations before reaching e = 0.055 successfully predicts challenges in HFEM
  transition" `[R p.14, verified in the PDF by the reader]`.
- Limits as stated: the metric "is planned for verification across various orbital families" `[R p.18]`: one family
  so far. "the presence of a fold bifurcation alone does not confirm a broken branch" `[D pp.10-11]`. In the
  quasi-bicircular model "no PO demonstrates a limit" for the halo interface region, so the diagnostic depends on the
  intermediate model `[D singh pp.17-18]`. Periodic orbits of the elliptic problem exist only for periods commensurate
  with the lunar period.
- Project has: `core/er3bp.py` with a published control (four of five Neelakantan & Ramanan orbits), `search/
  er3bp_periodic.py`, `genome/er3bp_continuation.py`.
- Missing: pseudo-arclength in e with fold detection for the catalogue's commensurate rows.

### T14. Speed as an independent variable in weak-stability maps (Mako & Salamon 2025)

- Does: scans the initial speed over [0.9 circular, 1.2 escape] at each position in the planar elliptic problem, where
  earlier work fixed the osculating eccentricity and so sampled circular-to-escape speeds only; at fixed position the
  stable speeds are a countable union of open intervals; stability depends on the primaries' true anomaly
  `[D mako-salamon pp.3-12]`.
- Limits: one revolution; Sun-Earth only; no transfer design; no tables `[D]`.
- Project has: `core/wsb.py` (the e_2-parameterised surface), `genome/bct_transfer.py`, the sweeps `#378` and `#681`.
- Missing: the speed axis and the true-anomaly axis (`#908` registers the sampling defect).

### T15. Free-return families on the excess-velocity sphere, with a working body and a passive target (Russell & Strange 2009)

- Does: full-revolution, half-revolution and generic free returns to the flyby body; the target body "is considered
  massless" and is met on a transit leg; the turn at each flyby of the working body is checked against
  sin(delta/2) = mu/(mu + r_p v^2) at a minimum altitude; the total time is an integer number of synodic periods
  `[R pp.144-149]`. Solutions with up to about 50 legs `[R p.145]`.
- Project has: the published rows (30 at V0), the turn gate, which reproduces their minimum altitudes within 1 km for
  9 of 10 cyclers (`#888`), `search/resonant_construct.py` (known defect with repeated bodies, `#480`).
- Missing: a chain builder in which BOTH moons bend the trajectory and same-moon returns are used to spread the turn
  (section 5.1).

### T16. Global search on a graph of (body, epoch) nodes with a flyby-defect cost (Bellome et al. 2023; Oshima 2024)

- Does: Lambert arcs on grids of dates; a node is a pair of consecutive (body, epoch) points; the edge cost is the
  flyby "velocity defect" (zero if the demanded turn fits, otherwise the minimum impulse: the same function as T2's
  conditional cost); dynamic programming keeps the best or the Pareto set per node; same-body resonant legs, "not
  generated by standard Lambert solvers", are added analytically `[D bellome]`. Oshima 2024 builds a graph of
  periapsis-to-periapsis arcs of the three-body flow and runs shortest-path searches, then regularised multiple
  shooting `[D oshima-2024]`.
- Limits: complete on the grid only; one-way transfers; no periodicity condition in either paper `[D]`.
- Project has: `search/campaign_runner.py`, a transfer network between rows (`#650`), Lambert with multiple
  revolutions.
- Missing: the triplet-node graph on the real ephemeris with a closure condition (section 5.4).

### T17. Regularisation and extended precision for close passages and very unstable orbits

- Sources: Levi-Civita in Font, Nunes & Simo and in Marco & Niederman `[R]`; a Levi-Civita form for the elliptic
  problem in Kumar et al. `[R CMDA sec. 6.2]`; "single shooting ... with an extended precision arithmetic of 128 bits"
  for a multiplier of 4.6e8 `[R Jorba-Cusco et al. 2018 p.12]`; Leiva & Briozzo name quadruple precision as the
  untried remedy for their failed continuations `[D p.238]`; four shooting sections for the bicircular L2 substitute
  `[R Rosales 2021 p.5]`.
- Project has: `core/cr3bp_regularized.py`; multiple shooting in several lanes. `#895` extensions "stop where Newton
  stalls at its 1 to 2 cm integration noise floor" (ledger), which is this limit.
- Missing: an extended-precision integrator option for verification of orbits with multipliers above about 1e10
  (the `#884` C11 members reach 1e14).

### T18. Linking numbers for torus connections (Owen & Baresi 2023, 2024)

- Limits as stated: the method "would not be applicable ... in non-autonomous Hamiltonian systems such as the
  bicircular restricted four-body problem" `[R pp.16-18]`; for planar two-degree-of-freedom problems the manifolds on
  a section are curves and direct intersection already works `[I, R]`. It detects libration-orbit transfers, not
  cyclers. `#883` (one control at the published energy) is the only use for it.

### T19. Patched two-body and three-body model for transfers between moons (Canales, Howell & Fantino 2021)

- Does: three-body arcs near each moon out to a sphere defined by an acceleration ratio, osculating conics about the
  planet in between, with analytic conditions for the two conics to intersect (their Theorems 1 and 2); one impulse
  at the intersection; carried to an ephemeris with three impulses. Titania to Oberon is a worked case `[D canales]`.
- Limits: one-way transfers between libration orbits; no return; "It is not a dynamical four-body problem, but rather
  a kinematical blending" `[D]`.
- Relevance: at the `#890` excess speed (0.27 km/s, a twelfth of the moons' orbital speed) the zero-radius patched
  conic is at its weakest and the first-order asymptotics are invalid (T3). This is the published seed model for that
  regime `[I]`.

## 2. Map from techniques to the project's problems

### Blocker 1. Heliocentric rows stuck at the lowest validation level (`#388`)

What `#388` did (project notes read by me): (a) a conic walk with Lambert legs, unknowns t0 and leg times, magnitude
or bend-hinge residuals, a closure tolerance of 0.1 km/s, on four descriptor-bearing rows; it reached DE440 and
settled on a low-energy solution with excess speeds about half the published ones; (b) a full n-body multiple shooter
with nodes at the planets at physical mass, seeded from the constructed ideal parent; three rows stalled at defects of
about 3.5e9 and one converged to a different high-energy family; (c) one low-energy row (4.3.1-5) recovered its
published excess speeds only under a 25-start, epoch-targeted search and failed the single-shot path.

What the sources say this differs in:

| Item | Russell & Ocampo (published) | `#388` (a) | `#388` (b) |
|---|---|---|---|
| Legs | Kepler propagation from a free excess-velocity vector; Lambert rejected `[P]` | Lambert | n-body |
| Unknowns | 5 per leg (vector, two times) `[R p.357]` | t0 and leg times | 6 per node at the planet |
| Flyby | zero radius, conditional impulse allowed `[P; R]` | zero radius, must be ballistic | physical mass from the start, node at the body |
| Objective | minimise summed violation, accept a nonzero cost `[R p.361]` | drive a residual to zero | drive defects to zero |
| Seed | seven cycles of the parent at a phase-matched epoch `[D]` | one cycle chain | constructed parent |
| Scan | 21 windows times 6 step counts, keep best `[R]` | 21 windows, one ladder | 3 epochs |
| Outcome | 9, 39, 74 of 203 under 1, 10, 300 m/s at some window `[P]` | "no row closes" | "no row closes" |

So the wall is, in large part, a difference in the QUESTION. The published answer for most parents is a cost of tens
to thousands of m/s (the Aldrin cycler averages 3,297 m/s over all windows `[R Table 6 p.366]`), and a method that
insists on zero will relax to whatever nearby family has zero. Three techniques apply:

- T2 as published, with Appendix C as the control (proposal P5).
- The 77 Appendix C blocks themselves, which need no optimiser at all (proposal P1).
- T1 after either, to replace zero-radius flybys by integrated ones for the few ballistic solutions (part of P1).

Objects: the 14 `russell-ch4-*` rows, the 200 `russell-ocampo-*` rows (182 with no validation level), `mcconaghy-
2006-em-k2`, `s1l1-2syn-em-cpom`. Progress would be: each Appendix C parent reconstructed, its printed total cost and
transit table reproduced, its encounters passed through the turn gate, and a validation level assigned by the
existing V3-ballistic or V3-powered rule.

T16 is the method for the remaining question `#388` named ("a global / family-targeted method"), and comes last.

### Blocker 2. Moon-system discovery produced closures that matched speeds and not directions

- The turn gate is, by T3 and T4, the theory's own admissibility condition at finite mass; it is necessary, and by
  the same theory not sufficient (the arcs must also be nondegenerate, and the chain's periapsis is fixed by the
  outer arcs, so a closure at any other periapsis is not on the branch `[D gomez-olle sec. 8e]`).
- T15 says how published moon cyclers satisfy it: many flybys of one working body with free returns in between.
  The project's enumerations were two-leg symmetric closures. Objects: Uranian pairs first (Ariel with any partner by
  the community poll's stated priority, as the ledger records under `#885`), then the 30 Russell & Strange rows as the control population.
- T1 is how a gated chain becomes a trajectory (`#890`, `#895` did this once). T19 supplies seeds where the excess
  speed is too low for a zero-radius model.
- Progress: one chain that passes the gate at every encounter, converges in the four-body model with both moons
  massive, and continues to the real ephemeris at three epochs, with the pre-registration discipline of `#895`.
- A separate, safer line: the 30 published Russell & Strange rows carried to integrated flybys (proposal P4).

### Blocker 3. Do Earth-Moon cycler families survive the Sun, eccentricity and the ephemeris

The papers split the question into three cases, and the project has been treating them as one:

| Parent | With one periodic perturbation | Published method | Catalogue rows |
|---|---|---|---|
| Period commensurate with the forcing (p/q, q = 1, 2) | isolated periodic orbits at selected phases | T7, T8 | the `#884` members; for eccentricity, the Casoliva rows (periods 2 pi, 4 pi, 6 pi) |
| Period commensurate at q = 3 or more | weak, higher-order; may behave as a torus | Brown JAS | the 8/3 members |
| Generic period, unstable | a hyperbolic invariant curve | T10 | most C11, C21, C31, C32 members |
| Generic period, stable | an elliptic invariant curve with bounded motion around it | T10 (Rosales form) | the stable Ross & Roberts-Tsoukkas members; the stable Casoliva rows |

With both perturbations (Sun and eccentricity) a commensurate orbit of one becomes a two-torus (Villegas-Pinto), and
only after that is an ephemeris transition attempted. Progress would be, per family, a table of which case applies
and one computed object per case with its lunar distance and stability, in the corrected models.

### Blocker 4. Few results survive adversarial review

Recurring features of the papers' own negative statements, each of which matches a project failure:

- A single non-convergence is not non-existence (Leiva & Briozzo's arcs; `#886` part one; `#882` R3a).
- A large continuation step changes family (Bradley & Russell p.13; `#884` first pass).
- A mesh hit is not a connection (Kumar SIAM p.26; `#882`).
- A control that passes under a wrong model is not a control (`#884`'s "Brown control" passed with the Sun reversed).
- A fixed-time formulation hides degenerate arcs (Bolotin & MacKay p.62; `#793`).

The proposals below are ordered so that the first ones have printed numbers on the expected side.

### Other problems in the ledger

- `#378`, `#681`, `#908` (capture sweeps): T14. `#378` is already stamped method-invalid by `#891`.
- `#790`, `#872` (periodic chains of resonant orbits): T5 gives the symbolic bookkeeping and the rule that a T arc
  may not be repeated; T12 is the cheap corrector for long chains.
- `#886` (four-body torus connections): T11's near-miss warning is the standing risk; no change of plan.
- `#883` (linking numbers): T18; autonomous spatial problem only.
- `#650` (transfer network): Braik & Ross's edges are "energy-preserving instantaneous heading-change" proxies, "a
  strong restriction relative to general impulsive maneuvers" `[R p.16]`, at one Jacobi constant; the project's
  network should say the same of its own edges.
- `#900` (V4 for moon tours): T1 with the `#895` force model is the unsoftened propagator that lane needs.

## 3. Ranked proposals

Probabilities are my judgement. (a) is the chance of a validated upgrade of existing catalogue rows; (b) the chance of
an object worth sending to an adversarial review. Costs assume the existing modules named.

### P1. Appendix C in full: 77 published real-ephemeris Earth-Mars cyclers reconstructed, gated and levelled

- Goal: reproduce every Appendix C solution of Russell (2004) from its printed data and assign each a validation
  level by the existing rule.
- Method applied: Russell's own reconstruction recipe (dissertation p.202, already implemented); for the solutions
  with zero printed cost, Bradley & Russell (2014) to replace zero-radius flybys by integrated ones.
- Inputs in the repository: `search/appc_corrected.py`, `search/s1l1_corrected.py`, the dissertation text layer, the
  parent index in `docs/notes/2026-06-07-russell-2004-member-tables-transcription.md` (35 of 77 listed), the V3 rules
  in spec section 14, `verify/turn_gate.py`.
- To build: a parser for the remaining blocks with a checksum against the printed "Total delta v"; a continuous
  propagation mode (one seed, printed manoeuvres executed, zero-radius turn applied only if the gate passes); a
  heliocentric provider for the `#895` corrector for the ballistic subset (parents 1, 54, 58, 83, 117 and others
  printed as 0.000000).
- Positive control: parents 83, 188 and 192 already reproduced (`#170`); the printed Earth-to-Mars transit table of
  each block is the expected side.
- Pass: per block, every encounter within the documented band and the printed total cost reproduced within 1
  percent; turn gate satisfied at every flyby the block calls ballistic. Level by spec: V3-ballistic or V3-powered
  only when the continuous cost is within the rule. Fail: reported per block, as for parent 192.
- A clean negative would be: blocks that do not reconstruct from their own printed data (a transcription or
  convention finding, worth an erratum note framed with care), or blocks that fail the turn gate (which would be a
  finding about the published solutions under the project's altitude floors).
- Cost: 6 to 10 agent-hours; minutes of compute. Integrated-flyby stage: 10 to 15 more.
- How it could fool us: each leg restarted from its own printed excess velocity matches arrivals without ever
  comparing incoming and outgoing vectors (the `#888` defect). Guard: the turn gate on the printed vector pairs, and
  the continuous mode. Second risk: mapping a block onto the wrong catalogue row (the June note found the shorthand
  prefixes differ between the ideal and real-ephemeris tables). Guard: map by parent number only; where no row
  carries the number, add the parent as its own literature row with `known-reproduction` status.
- (a) 85 percent for at least ten rows at V3 by the existing rule (some as added literature rows). (b) 0: these are
  published solutions.

### P2. The published one-moon generator: Barrabes & Gomez seeds and Casoliva's three-step continuation

- Goal: reproduce Casoliva's Table 3 from seeds at mass ratio 1e-6 without using Table 3's states, then run the same
  generator over resonances they did not tabulate.
- Method applied: Barrabes & Gomez 2002, 2003 (seeds); Casoliva et al. 2010 section IV.C (continuation).
- Inputs: `search/mu_continuation.py`, `search/cr3bp_periodic.py`, `core/cr3bp_regularized.py`,
  `search/earth_moon_resonant_families.py::TABLE3_ROWS`, the 2008 conference paper's table of seeds at 1e-6.
- To build: the seed (their eqs. 14 to 18 `[P p.1628]`); the 500-point scan in C per (p, q); the three-step driver
  with a refined-periselene trigger; the family walk in C at the physical mass with a resonance-crossing detector;
  symbol labels.
- Positive control, two-ended: the 2008 table of seeds at 1e-6 (start) and Table 3 at 0.0121529529 (end), including
  the count "four family segments for the seed 7-3 resonance, but only three intersections" `[P p.1629]`.
- Pre-registered pass: at least seven of the nine clean Table 3 rows reached from seeds, each matching the printed
  C_J, period and stability class; the 7-3 segment count reproduced. Fail: fewer than five.
- Then: resonances with p/q above 2.3 that the paper names but does not tabulate (8-3 is named `[P p.1627]`) and any
  other (p, q) with q up to 4. Label under spec 16.4: `known-class-member` (the authors searched this class with this
  method), not a novelty label.
- A clean negative: the generator reproduces the table and finds no further tight cyclers above the lunar surface;
  that would bound the class at the Earth-Moon mass, which the paper does not state.
- Cost: 15 to 25 agent-hours; a few CPU-hours (the paper reports 3 to 23 hours for 500 values of C in 2010 hardware,
  and 35 to 40 minutes per continuation `[P pp.1628-1630]`).
- How it could fool us: converging onto a Table 3 row from its own printed state (circular); a branch jump in mu that
  lands on a first-species orbit with no lunar passage and the right period; a sampled, not refined, periselene.
  Guards: start files contain only seeds; minimal-period and periselene checks at every step; tangent-cosine guard;
  the seven-of-nine rule set before running.
- (a) 60 percent (the nine rows gain an independent second derivation and symbol labels; no level change).
  (b) 45 percent (further members of the published class as catalogue rows).

### P3. `#884` rerun with the published selection rules, and the arcs question settled

- Goal: decide, with controls, whether the 5/2 members of C32 and C31 have true periodic orbits at full Sun strength
  in the quasi-bicircular model, where Leiva & Briozzo printed arcs.
- Method applied: Leiva & Briozzo 2008 (phases), Brown et al. 2025 (multiplicity, objects versus orbits), Oshima 2022
  (fold following), Rosales et al. 2021 (multiple shooting over sections).
- Inputs: `search/sun_forced_periodic_884.py`, corrected `core/bcr4bp.py` and `core/qbcp.py`, the `#896` tests of
  Leiva & Briozzo's Tables 1 to 5 and Oshima's tables.
- To build: the phase formula for non-symmetric orbits; deduplication by time shift; the explicit criterion below.
- Pre-registered rule (from reader C's proposal, which I adopt): members with period (p/q) Sun periods, q = 1 or 2;
  up to four start phases per member, expecting two distinct objects; a PERIODIC ORBIT closes over p Sun periods with
  the Sun back at its starting phase, by multiple shooting, with a monodromy; closure over the half interval with the
  Sun at pi is recorded as an ARC; pseudo-arclength in the Sun's strength with a tangent guard, each branch recorded
  as reaches-one, folds-back, collides or bifurcates; q = 3 reported apart.
- Positive controls, in order: Leiva & Briozzo Table 3 arc 180A_1_t1 (must be reproduced AS AN ARC); Table 2 orbit
  013_t3 (as an orbit); Oshima Table 4 row 1. A further cross-check `[C, I]`: for an arc whose forcing is nearly pure
  quadrupole, the multiplier over the doubled interval should be about the square of the arc's stability parameter.
  The printed arc values 49.9, 5153 and about 5370 `[R Table 5]` squared are 2.5e3, 2.7e7 and 2.9e7; the `#884`
  reviewer's corrected-sense bicircular multipliers for the same members are 3.1e3 to 2.2e4, 9.8e6 to 1.9e7, and
  2.6e7 to 4.9e7 (ledger). That is the right order for all three, assuming their parameter is close to the multiplier
  when large (their definition was not checked).
- Pass: the three members have periodic orbits in the corrected quasi-bicircular model that reproduce the arcs'
  lunar distances within the papers' printed spread. Fail: any member exists only as an arc under multiple shooting.
- A clean negative: an orbit that closes as an arc and not over the full interval would show that the odd harmonics
  (about 1/389 of the quadrupole) or the model's own terms destroy it, which the published paper could not test.
- Cost: 10 to 15 agent-hours; minutes of compute.
- How it could fool us: calling an arc an orbit; counting time-shifted copies as different; a jump of branch. Guards
  are in the rule. The statement to make if it passes is narrow: "periodic orbits of the quasi-bicircular problem at
  the 5/2 members, where the 2008 paper obtained arcs and attributed the failure to its numerical method".
- (a) 30 percent (no level exists for this evidence; rows gain a sourced note). (b) 65 percent.

### P4. Russell & Strange's moon cyclers carried to integrated flybys with the `#895` corrector

- Goal: for the 30 catalogued Russell & Strange rows, obtain multi-cycle trajectories with every flyby integrated in
  a real-ephemeris force model, with or without small manoeuvres, and compare with the paper's printed costs.
- Method applied: Bradley & Russell 2014 (written by the same group after the 2009 paper asked for "Improving
  strategies for transitioning" `[R p.145]`); their Jovian example has the same bodies as ten of the rows.
- Inputs: `search/titania_oberon_realeph_895.py`, `verify/turn_gate_closures.py` (the rows' ideal chains already
  reproduce the printed minimum altitudes), `core/satellites.py`.
- To build: Saturn and Jupiter force-model providers with a positive control of the `#895` kind (each moon
  propagated as a test body against its kernel; current satellite kernels must be fetched, the cached ones are
  Voyager-era); an ideal-to-ephemeris position interpolation as in the 2009 paper; a per-leg manoeuvre variable so
  that the result can be a cost, not only a yes or no.
- Positive control: Russell & Strange Table 8 single-cycle high-fidelity costs (32, 63, 11, 121 m/s) `[R p.155]`.
  This is a loose control (their tool and force model differ), so pre-register "same order, same ranking".
- Pass: for at least three rows, three or more cycles with all flybys integrated, altitudes inside a band set
  beforehand around the printed ones, two integrators, cost no greater than the printed patched-conic ephemeris cost.
- A clean negative: cost rising with the number of cycles for every row would quantify the paper's own remark that
  the planetocentric transition is "substantially more difficult", per cycler.
- Cost: 25 to 40 agent-hours; hours of compute.
- How it could fool us: the converged trajectory keeps the encounter sequence and loses the cycler (flybys drift up by
  thousands of km, as the registered `#895` route did). Guard: altitude and excess-speed bands fixed before the run;
  a result outside them is reported as "a different trajectory".
- (a) 50 percent (rows move from V0 to a real-ephemeris level under a rule the owner would have to set for moon
  rows). (b) 25 percent (a multi-cycle integrated solution is not printed in the paper).

### P5. Russell & Ocampo's optimiser as published

- Goal: give every constructible heliocentric row its real-ephemeris seven-cycle cost by the published method.
- Inputs: `search/generic_return.py`, `search/cycler_search.py`, `search/continuation.py`, `core/flyby.py`.
- To build: the 5n transcription with the conditional cost, a bounded least-squares or SLSQP stand-in for elastic
  mode, the seven-cycle seed, the scan.
- Positive control: Appendix C (from P1): starting from parent 83's ideal seed the optimiser must return a solution
  of zero cost at the printed window; from parent 188 about 0.436 km/s; and the tier counts 9, 39, 74 as a
  population check if all 203 parents can be constructed.
- Pass: the controls within 10 percent in cost and the right ordering of ten parents. Fail: the optimiser relaxes
  off-family on the controls (then the `#388` diagnosis of family selection stands for this method as well).
- A clean negative is informative: it would say the published tiers depend on details of SNOPT's elastic mode that
  a stand-in does not reproduce.
- Cost: 30 to 50 agent-hours; CPU-hours.
- How it could fool us: a low cost reached by changing the sequence's timing so far that it is no longer the parent
  (the paper notes the optimiser "sometimes departs" from the total duration `[R p.366]`). Guard: report the excess
  speeds and durations against the parent's and reject drift beyond a pre-set band.
- (a) 45 percent. (b) 5 percent.

### P6. A test of the printed persistence conjecture for stable cyclers

- Goal: for the stable Earth-Moon members printed by Ross & Roberts-Tsoukkas and the stable Casoliva rows, compute
  the Sun-forced invariant curve at physical Sun mass and test bounded motion in the bicircular model and then in the
  ephemeris.
- The conjecture, as printed: "near-commensurable stable cyclers should persist under the dominant perturbations
  neglected here, including primary eccentricity and, in the Earth-Moon case, solar gravity. This conjecture can be
  tested in the elliptic and bicircular restricted problems and, ultimately, in full ephemeris models"
  `[P 2026 arXiv version, closing paragraph, text layer]`. The conjecture is about stable members NEAR a
  commensurability with the Sun; the object computed here for each stable member is its invariant curve `[I]`.
- Method applied: Rosales et al. 2021 section 3 (invariant curve with rotation number 2 pi T_Sun/P, continuation in
  Sun strength, gaps crossed by detour); Park & Howell's fold monitor for the eccentricity leg.
- Inputs: the invariant-circle corrector of `search/ccr4bp_strob_connection.py`, `core/bcr4bp.py`, `genome/
  hill_screen.py` (the Hill-radius pre-screen from `#391`), the stable rows' states.
- Positive control: Rosales et al. 2021 Table 3 (a printed torus, rotation number 1.380018549762754 `[R]`) and
  Table 2 multipliers (already a test under `#896`).
- Pass or fail, pre-registered per member: curve converged on a node ladder with invariance error below 1e-9 under an
  independent propagator; then 100 Sun periods in the bicircular model and 10 years in DE440 from three epochs without
  manoeuvre, lunar encounters retained. A member that fails the Hill pre-screen is reported as such, not run.
- A clean negative is a result about a printed conjecture, and the `#389` case (a stable (3,3) branch that escaped
  in the ephemeris) says a negative is likely for the large members.
- Cost: 20 to 30 agent-hours; hours of compute.
- How it could fool us: a curve that "converges" with too few nodes; boundedness over too short a time. Guards: the
  node ladder; a known-unstable member must be seen to leave within the same horizon (negative control).
- (a) 35 percent. (b) 45 percent.

### P7. Casoliva's rows in the elliptic problem, with the fold count

- Goal: continue the Casoliva rows whose periods are whole lunar periods to lunar eccentricity 0.0549 as periodic
  orbits of the elliptic problem, and record folds.
- Method applied: Gomez & Olle 1991 Theorem 11 (symmetric second-species solutions persist for every eccentricity
  only with period 2 k pi and perpendicular crossings with the secondary at pericentre or apocentre `[D gomez-olle
  p.141]`); Park & Howell 2024 (pseudo-arclength in e, fold detection).
- Positive control: Gomez & Olle Part II's printed elliptic families at mass ratio 1e-6; Park & Howell's 53:21 halo
  (fold near e = 0.003 to 0.004) and 13:7 (no fold) `[R p.11, p.14]`.
- Pass: both phases of each row followed to 0.0549 or to a recorded fold. Then one ephemeris multiple-shooting
  attempt per row, to see whether folds and transition difficulty go together for this class (the authors claim this
  for one family only).
- Cost: 15 to 20 agent-hours. (a) 40 percent. (b) 30 percent.
- How it could fool us: the rows' periods are 2 k pi only to about 1e-5; the continuation must start from the family
  member of exact period, and the independent variable is the true anomaly, not time `[D gomez-olle sec. 8f]`.

### P8. Two-moon chains with resonant returns, gated, then continued

- Goal: enumerate Uranian two-moon chains in which each moon may be flown several times in succession on free
  returns, keep those that pass the turn gate at every encounter, and continue survivors as `#890` and `#895` did.
- Method applied: Russell & Strange 2009 (architecture), with both moons bending; Bradley & Russell (continuation).
- To build: a fixed-excess-speed return solver for same-moon legs (T4's remark: not Lambert); the chain enumerator
  with per-encounter self-consistency for repeated bodies (the `#480` defect); the period condition (whole synodic
  periods).
- Positive control: the enumerator restricted to one working moon must regenerate Russell & Strange's Table 2 to 6
  cyclers (the gate already reproduces their altitudes), and restricted to two legs must return exactly `#890`.
- Pass: at least one chain converged in the four-body model at physical masses and in the ephemeris model at three
  epochs. A clean negative, with the control passed, would be a real empty region for the stated chain lengths.
- Cost: 40 to 60 agent-hours. (a) 5 percent. (b) 30 percent.
- How it could fool us: "it closed" in the ideal model. Only integrated trajectories count, and the owner's rule that
  `#563`-class sweeps are not repeated applies: this is a different genome (returns), and the note must say so.

### P9. A global search for zero-radius real-ephemeris cyclers by dynamic programming (ranked last; costs more than the others)

- Method applied: Bellome et al. 2023 with a closure condition added. Control: it must find Appendix C parents 1 and
  83 at their printed windows. Cost: 60 agent-hours or more. (a) 15 percent. (b) 15 percent. Do after P1 and P5.

## 4. The six candidate ideas, tested against the sources

**(i) Bradley & Russell applied to the heliocentric cyclers to break `#388`.** Not sound as stated. The method's
input is a zero-radius trajectory already on an ephemeris `[R]`; its authors' heliocentric example hardly needed the
continuation `[D p.19]`; Russell & Strange say the patched conic "is clearly better in the case of the heliocentric
cyclers" `[R p.152]`; Russell & Ocampo expect the step to integrated flybys to be minor `[P]`. `#388` failed one step
earlier, and partly because it asked for a ballistic closure where the published result is a cost. Untested by the
authors for this use in any case: more than nine encounters in one solve (a seven-cycle cycler has 15 to 49 legs
`[R]`), durations of years, a mean-Keplerian fit over decades, periodicity `[R, I]`. Smallest decisive experiment:
P1's last stage, the continuation applied to Appendix C parent 83 (zero printed cost) one cycle at a time. It is
published for open tours and not, in the held papers, for cyclers; the result would be a validation, not a discovery.

**(ii) Second species as a generator.** Sound for one moon, and already published for the project's objects
(Casoliva 2010), so it is a reproduction with a two-ended control (P2). It generates the Casoliva-type cyclers
(C < 3), not C11 to C32. Which chains survive to finite mass: the theorems give no bound `[R]`; numerically, chains
are lost between 5e-4 and 1.5e-3 `[R]`, fixed-period families "do not exist uniformly in mu" `[P]`, and lunar impact
is the usual end `[P]`; the p-q generating orbits are first-to-second species bifurcation orbits where the passage is
not O(mu) `[R]`, so a turn gate applied at zero mass to those seeds can reject the true generators `[I, reader B2]`.
Gomez & Olle's explorations are at 1e-6 only and describe no corrector `[D]`. For two moons no theorem applies
`[R]`; what carries over is the local flyby analysis and the turn gate; the published architecture is Russell &
Strange's (P8). Smallest decisive experiment: the 7-3 seeds at 1e-6 continued to the three printed 7-3 rows.

**(iii) Phase selection and higher-order Melnikov functions for the `#884` rerun.** Sound, and they fit together:
Leiva & Briozzo's quadrupole condition is Brown et al.'s lowest-order function for q = 1, 2; q = 3 first appears at
the octupole, about 1/389 weaker `[D brown-jas; R]`. The question the papers "leave open" is answered by Leiva &
Briozzo's own text: arcs were adopted because single shooting over five Sun periods failed numerically
`[P p.239, text layer]`.
The project's orbits close over the full interval with multiple shooting (the `#884` table lists these members as
"5/2 (2)"), so there is no contradiction, and the squared stability parameters agree in order `[C]`. Published for
these members: arcs. Not found in the held papers: the periodic orbits. Smallest decisive experiment: P3's first
control (reproduce one arc as an arc, then close it over the doubled interval in the same model).

**(iv) Folds in eccentricity as a predictor.** Sound as a diagnostic; as a predictor it is the authors' claim for
one family, "planned for verification across various orbital families" `[R]`, and Singh et al. show it depends on
the intermediate model `[D]`. It applies only to parents commensurate with the lunar period: the Casoliva rows, not
the Ross & Roberts-Tsoukkas rows, whose periods are 1.64 to 3.09 lunar periods `[D gomez-olle sec. 8f]`. Smallest
decisive experiment: P7 on casoliva-2-1a and casoliva-7-3a.

**(v) Floquet and centre-manifold reduction for bounded motion and a maintenance cost.** Not sound for unstable
planar cyclers: there is no centre part to reduce to, the normal form's radius shrinks with small divisors, and the
authors themselves abandoned it at L2 in the bicircular problem `[R]`. None of the papers gives a maintenance cost;
Leiva & Briozzo promise a control algorithm "in a forthcoming paper" and state only that an error "is multiplied by
a factor ~ s_i at each orbital period" `[R p.244]`. What is sound: the invariant curve as the Sun-forced counterpart
of a generic-period cycler (T10), and, for stable cyclers, bounded motion around it (P6). A maintenance cost for
unstable flyby chains is better taken from the flyby sensitivity itself: the turn changes with impact parameter at a
rate of order 1/mu `[R MacKay p.2]`, which is what `#890` measured (1 m/s moves the periapsis 1,200 to 1,400 km). A
per-flyby targeting budget from the state transition matrix is the honest V2 for such orbits `[I]`; it needs an
owner ruling, as `#890` already notes.

**(vi) Speed as an independent variable in the weak-stability sweeps.** Sound as a sampling correction (`#908`).
Published for the Sun-Earth planar elliptic problem, one revolution, no transfers `[D]`. It widens what the `#378`
and `#681` negatives did not sample; it does not by itself make a repeating capture more plausible. Smallest decisive
experiment: reproduce the true-anomaly thresholds of their Figs. 7 and 8, then rescan a few `core/wsb.py` grid points
over speed. Low expected yield for cyclers.

## 5. Techniques not on the candidate list

1. **Spread the turn over repeated flybys of one moon (Russell & Strange).** The withdrawn Uranian rows demanded 1.8
   to 28 times the available bend in ONE flyby per moon. Published moon cyclers reach large net turns with many
   flybys of the working body on full-revolution and generic returns, up to 54 flybys `[R p.156]`. Chains with returns
   were never enumerated with both moons bending (P8).
2. **Ask for a cost, not a closure (Russell & Ocampo; Bellome).** The conditional flyby cost turns "does it close"
   into "what does it cost", is what the published heliocentric results are stated in, and is what the spec's
   V3-powered rule already accepts (P1, P5).
3. **Test a printed conjecture (Ross & Roberts-Tsoukkas).** The authors state how their claim about solar gravity
   and eccentricity should be tested and have not printed the test. The project holds the models, now with controls
   (P6). Either outcome is a result with a clear published referent.
4. **Global search on the real ephemeris by dynamic programming over (body, epoch) triplets (Bellome).** This is a
   concrete form of the "global / family-targeted method" that the `#388` verdict asked for and never specified; it
   enumerates families on a grid instead of following one seed (P9).
5. **Fixed-energy arcs for same-body returns (Bolotin & MacKay's remark).** Whole and half revolution returns are
   degenerate at fixed time and regular at fixed excess speed. A return solver posed on the excess-velocity sphere
   removes the singular Lambert cases in P8 and in any reconstruction of f and h legs.
6. **Two intermediate models before the ephemeris, one perturbation at a time (Villegas-Pinto).** Choose the member
   commensurate with one forcing so that it stays periodic "by design", then add the other to get a two-torus, then
   go to the ephemeris. For cyclers: Casoliva rows are commensurate with the lunar period (eccentricity first); the
   `#884` members with the Sun (Sun first).
7. **Follow the whole solution curve, including negative forcing (Oshima; Rosales).** It shows which three-body
   orbits a forced orbit connects to and exposes forced orbits with no three-body parent at the same period.
8. **Orbit versus object (Brown et al.).** Count equivalents after removing time shifts. The `#884` counts ("39 of
   56 branches") should be restated this way.
9. **Squared-parameter check between an arc and its doubled orbit** `[C, I]`: a cheap consistency test between a
   published half-interval result and a full-interval computation (used in P3).
10. **Extended precision as a verification instrument** for multipliers above about 1e10 (Jorba-Cusco et al. 2018;
    Font, Nunes & Simo): the `#884` C11 members and any second-species orbit at small mass are in that range.
11. **A patched two-body and three-body seed model for slow encounters (Canales, Howell & Fantino)**: for chains
    like `#890`, where the zero-radius model and the first-order asymptotics are both outside their range.
12. **Symbol labels (Font, Nunes & Simo's (S or R, p, q, s))** as an identity key for flyby-chain rows: cheaper and
    more discriminating for the literature gate than body lists, and it encodes the rule that a T arc is not repeated.

## 6. What the project should stop doing

1. Stop treating "does not close ballistically in the real ephemeris" as a failure of method for heliocentric
   parents. The published result is that 9 of 203 do at some window `[P]`.
2. Stop Lambert-leg, time-only correctors for cycler transition; the published method rejected them for stated
   reasons `[P; R p.356]`. Stop seeding an n-body shooter with nodes at physical-mass planets.
3. Do not dispatch `#898` as registered; replace it by P1 and P5.
4. Stop seeding tight Earth-Moon cyclers from two-body ellipses at the physical mass (`#780`(d) records a fourth
   failure of this); the authors of the rows say it fails and give the remedy `[P p.1627]`.
5. Do not use the second-species generator, or its first-order formulas, on the C11 to C32 families or on `#890`.
6. Stop reading a single failed continuation or a single non-converged torus as non-existence, and stop single
   shooting over five or more forcing periods for multipliers above about 1e3.
7. Do not build the Lie-series centre-manifold machinery for cyclers.
8. Stop symmetric two-leg closure enumeration as a discovery method (already on the do-not list); the reason from the
   sources is in section 5.1.
9. Do not extend the linking-number pipeline to forced models; its authors say it does not apply `[R]`.
10. Stop counting forced branches before deduplication, and stop mixing q = 3 members with q = 1, 2.

## 7. Papers still needed

Held but not digested (digest before P2): Barrabes & Gomez, "Spatial p-q resonant orbits of the RTBP", Celest. Mech.
Dyn. Astron. 84:387-407 (2002); Barrabes & Gomez, "Three-dimensional p-q resonant orbits close to second species
solutions", Celest. Mech. Dyn. Astron. 85:145-174 (2003); Casoliva et al., AIAA 2008-6434 (its table of seeds).

Not held (citations as printed in the held papers unless marked):

- Bruno, A. D., "On periodic flybys of the Moon", Celest. Mech. 24:255-268 (1981). Directly on the one-moon case.
- Henon, M., Generating Families in the Restricted Three-Body Problem, Lecture Notes in Physics Monographs 52,
  Springer (1997); and its second volume (2001; citation from memory, to be verified).
- Henon, M., "Sur les orbites interplanetaires qui rencontrent deux fois la terre", Bull. Astron. (3) 3:377-402
  (1968) (the held papers print the pages three different ways).
- Hitzl, D. L. & Henon, M., "Critical generating orbits for second species periodic solutions of the restricted
  problem", Celest. Mech. 15:421-452 (1977).
- Guillaume, P., "Periodic Symmetric Solutions of the Restricted Problem", Celest. Mech. 8:199-206 (1973) (the only
  cited computation near mass ratio 1e-2); Guillaume, Celest. Mech. 11:449-467 (1975).
- Perko, L. M., Celest. Mech. 14:395-427 (1976); 16:275-290 (1977); 24:155-171 (1981).
- Olle, M., doctoral thesis (1989), cited by Gomez & Olle Part II as the full description of the families.
- Bolotin, S., Celest. Mech. Dyn. Astron. 93:345-373 (2005) and Discrete Contin. Dyn. Syst. 14:235-260 (2006):
  second species in the slightly elliptic problem, the nearest time-dependent theory.
- Lantoine, G. & Russell, R. P., "Near ballistic halo-to-halo transfers between planetary moons", J. Astronaut. Sci.
  58(3):335-363 (2011): continuation from three-body to four-body ephemeris models for moon pairs.
- Rhouma, M. B. H. & Chicone, C., "On the continuation of periodic orbits", Methods Appl. Anal. 7:85-104 (2000)
  (the theorem behind Brown et al.'s resonance condition; citation from memory, to be verified).
- Jorba, A. & Villanueva, J., "On the persistence of lower dimensional invariant tori under quasi-periodic
  perturbations", J. Nonlinear Sci. 7:427-473 (1997) (the theory behind T10; citation from memory, to be verified).
- Already listed under `#909` and `#884`: Peng & Xu 2015; Sanaga & Howell (DOI 10.1007/s42064-024-0250-4); Komachi,
  ASC 2026.

## 8. Methods in papers previously set aside (#919)

Added 2026-10-04 by the `#919` agent. The owner asked whether the papers set aside as out of scope, triaged,
reference-only or background-only might hold techniques for the project's problems. Those verdicts were given on what
the papers are about (a sail, a formation, an asteroid, a neural network). This section judges each on what its method
can do for the problems of section 2. Tags as in the legend at the top: `[P p.N]` READ by me at the source; `[D note]`
READ in a project digest or triage note, not checked at the source by me; `[I]` INFERRED; `[C]` my calculation. No
source readers were used for this section: two were dispatched and returned nothing before this version was written. The
papers were found with the corpus-index search the task named (rows marked out of scope, triaged, reference-only or
background-only), plus the papers the task named (Villegas-Pinto et al. 2023; de la Fuente Marcos 2018) and two
owner-supplied papers filed under `#869` and marked out of scope or low relevance (Zhang et al. 2024 and 2025). Papers
the synthesis already used for their methods (Park & Howell 2024, Mako & Salamon 2025) are not repeated.

### 8.1 Periodic orbits without symmetry and period doubling in eccentricity (Wang, Ye, Xiao & Li 2026)

- Method: resonant halo orbits of the Sun-Mercury elliptic problem found by multiple shooting on the full six-state
  closure `X_{N+1} - X_1 = 0` with midpoint matching, solved with fmincon or fsolve `[D wang-ye-xiao-li sec. 2.3, p5]`.
  Continuation is staged: eccentricity, then sail loading, then cone angle, then the sail azimuth beta `[D sec. 2.4,
  p6]`. Orbits with no symmetry appear only under the azimuth law: the symmetry is broken by the control parameter beta,
  and "all six states are needed" for the initial condition `[D sec. 2.3, p16]`. Along the sail-free 3:1 families,
  monodromy eigenvalues reach -1 at e = -0.0179, 0.0089, 0.04 (L1) and -0.0189, 0.0094, 0.044 (L2), where a negative e
  denotes the aphelion start; the doubled-period 6:2 branches born at the second and third reach e = 0.2056 and the
  first stop converging near 0.037 and 0.038 `[D sec. 2.5, pp7-8]`.
- Bears on: the Casoliva 7:3(b) and 7:3(c) members that `#884` skipped for having no perpendicular crossing
  (`docs/notes/2026-10-04-884-sun-forced-em-cyclers.md`), and P7 (`#917`).
- Hypothesis (1) tested. Wang's orbits lack symmetry because a forcing parameter breaks it in a system that is already
  time-periodic, where periodic orbits are isolated and a six-state closure is a square system `[I]`. The 7:3(b) and
  7:3(c) members are asymmetric orbits of the autonomous three-body problem, where a periodic orbit carries a phase
  freedom and lies on a one-parameter family, so the closure needs a phase condition and a fixed Jacobi constant or period
  `[I]`. The project already has the full-state corrector for that case: `genome/asymmetric_branch.py` routes through
  `search/cr3bp_general_periodic_3d.py::correct_general_periodic_3d` in "full-asymmetric mode" with the period free
  (module docstring, read). Once Sun forcing is on, `#884`'s multiple-shooting problem (`newton_fixed_eps` in
  `search/sun_forced_periodic_884.py`) is already a six-state closure with no symmetry imposed (code, read). The gap that
  kept 7:3(b) and (c) out is in the three-body stage: `walk_family` and `correct_symmetric_fixed_period` are symmetric
  only, and the phase rule for orbits without the symmetry is missing (already listed under T7). Wang's paper adds no
  method for that. Its period-doubling route is already in `genome/family_switch.py` for the autonomous case.
- What it does give: printed eigenvalue crossings in eccentricity, a control for a -1 multiplier monitor along an
  eccentricity continuation. That monitor is an add-on to P7: Casoliva's rows continued in e may period-double before
  they fold, and a doubled branch is a further object of the elliptic problem `[I]`. The control prints three or four
  decimals only and no states `[D sec. 5]`, so it tests the monitor and the continuation, not a closure.
- Verdict: NOT USEFUL as a source of a method the project lacks for asymmetric cyclers (the project has the method; the gap is plumbing
  in `#884`'s three-body walk plus T7's phase rule). USEFUL as a control for a -1 monitor added to P7. No separate proposal.

### 8.2 Bifurcation equations for symmetry breaking (Shu & Lin 2025)

- Method: a Lindstedt-Poincare series about the collinear points of the spatial elliptic problem, to order 7, with
  "bifurcation equations" Delta(eta, e, alpha1, alpha2, alpha3) = 0 whose nonzero roots eta mark pitchfork bifurcations
  `[D shu-lin secs. 0, 6]`. Near the collinear points the eccentricity shifts the circular bifurcation condition at order
  e^2 `[D sec. 7]`. No differential correction, continuation, Floquet analysis or integration is reported, and no table is
  printed `[D sec. 0]`.
- Hypothesis (5) tested. The series is an expansion in the amplitude about an equilibrium. Cyclers are large-amplitude
  orbits with close passages of the small primary, far outside any such expansion `[I]`. A semi-analytical equation for
  where asymmetric families branch from symmetric cycler families would need an expansion about the cycler family, which
  this paper does not attempt. The project already locates such branch points numerically (multiplier +1 crossings,
  `search/bifurcation_detector.py`, then `genome/asymmetric_branch.py`).
- Verdict: NOT USEFUL. The method is local to the libration points by construction.

### 8.3 Two intermediate models and a two-torus; the `#895` arcs (Villegas-Pinto, Baresi, Locoche & Hestroffer 2023)

The method is already in T10 and section 5.6. What `#919` adds is a direct answer to the question `#907` left open:
are the `#895` real-ephemeris Titania-Oberon arcs pieces of an invariant torus that can be computed as one object?

- Method: GMOS invariance condition with a Fourier rotation operator, multiple shooting with `N_0` nodes, phase and
  pseudo-arclength conditions, Newton to 1e-10 `[D villegas-pinto sec. 2.2, pp7-11]`. A resonant orbit is first made
  periodic in one intermediate model "by design" (Sun first for synodic resonances, eccentricity first for sidereal),
  then the second perturbation is switched on to make a two-torus, with `N_1` = 30 to 50 points on the circle (more for
  lower perilunes) `[D pp8-9, p15]`. The orbits treated are near-stable: "remains very close to the near-stable behavior
  presented by their periodic counterparts" `[D sec. 2.3, p24]`. Many other resonant orbits "fail to do so, even after
  being computed in the intermediate models" `[D sec. 2.2, p12]`. No numbers are printed: no states, periods or rotation
  numbers `[D sec. 2.3]`.
- Answer for `#907`, all `[I]`. In the real ephemeris the forcing on the spacecraft has many frequencies: two moon mean
  motions, their eccentricities and apsidal precession under J2 and J4, the other three major moons, the Sun. The arcs lie
  on no finite-dimensional torus of that model that the held papers would let one compute. In an intermediate model with
  ONE added frequency they can. The `#890` orbit is periodic in the planar circular four-body model, with period a whole
  number of Titania-Oberon synodic periods. Adding Oberon's eccentricity, with its apsidal precession under J2, adds one
  incommensurate frequency. By the persistence result Rosales et al. 2021 state for a forced periodic orbit (T7; read
  at the source by the synthesis author, not by me) (the general theorem is Jorba & Villanueva 1997, not held), the periodic orbit should become a two-dimensional
  invariant torus if it is nondegenerate. It is hyperbolic, with multipliers about 1e6 per period (`#890` ledger). Such a
  torus would explain one feature of `#895`: arcs exist at all five epochs tried, 2030 to 2045, and each epoch would be a
  different phase on the torus.
- What it would and would not settle. A torus computed in that intermediate model is one object with a definition, a
  rotation number and a stability. It is still not the ephemeris trajectory. The ephemeris arcs would be compared with
  it, and the residual measures the other perturbations. The torus is strongly hyperbolic, unlike Villegas-Pinto's, so
  single shooting over one torus period cannot work. It needs the multiple-shooting invariant-circle variant that T10
  lists as missing and that P6 (`#916`) also needs `[I]`. Kumar's bundle method (T10) applies here, since there is a real
  saddle pair.
- Verdict: USEFUL, as proposal Q2 below. The positive controls cannot come from this paper, which prints no numbers.

### 8.4 Co-orbital asteroids: clone ensembles and the Arjuna domain (Carruba et al. 2026; de la Fuente Marcos & de la Fuente Marcos 2018)

- Carruba et al. method: each co-orbital asteroid is integrated for 78,000 years with two Bulirsch-Stoer codes (SWIFT-BS
  at tolerance 1e-8 and REBOUND); objects where the two codes disagree are labelled "dubious" (alternating); each object
  is labelled from its argument-of-perihelion time series (circulating, librating, intermittent; Table 1); then 729 clones
  are formed from the nominal elements plus or minus one sigma in each of six elements (3^6) and summarised by the mean
  angle and libration amplitude across clones; Hill-sphere encounters are counted over 10^5 years `[D carruba secs. 1.3,
  1.4, 1.5, 4]`. "There is no new equation, Hamiltonian or semi-analytical model in this paper beyond (1) and the
  Hill-radius expression (2)"; the averaged model of Michel & Thomas (1996) is reproduced as a figure and used
  qualitatively `[D sec. 1.2]`.
- De la Fuente Marcos method: Monte Carlo sampling of orbital elements, uniform within the Arjuna box, to compute the
  minimum distance between orbits, relative speed at perigee and gravitational focusing `[P text layer, sec. 3, "Computing
  the perigee: a Monte Carlo"; D fuente-marcos]`. The synodic periods of the class are of decades or longer `[P sec. 1]`.
- Hypothesis (4) tested. First half: the averaged co-orbital model is not in Carruba et al.; they derive none `[D]`, so
  there is nothing to take. Second half: the clone ensemble as a robustness test. For the project's unstable rows it
  would only restate the multipliers (with growth of 1e3 per flyby pair, a one-sigma clone leaves the chain within a few
  flybys, which the monodromy already says) `[I]`. It adds information in one place: labels of the "bounded versus
  divergent" kind, the criterion by which the quasi-cycler rows `#339` and `#344` pass, which are themselves classifications
  of a long time series as Carruba's are `[I]`. There the robustness of the LABEL across an ensemble spanning the
  printed-digit interval of the source and the integrator tolerance is a cheap, pre-registrable check. The
  two-integrator "dubious" rule is already project practice (`#895` used two integrators; the `#896` Mako & Salamon
  checks ran over the printed rounding interval). Quasi-satellite and horseshoe orbits of one planet repeat encounters
  with that planet only, at synodic periods of decades for the Arjuna class, and quasi-satellites never leave the
  planet's vicinity. They are not cyclers between two bodies, and no method in either paper builds one `[I]`.
- Verdict: Carruba et al., POSSIBLY USEFUL as a robustness check on quasi-cycler labels only. What would decide it: a
  quasi-cycler row whose bounded label is close to its threshold. If none is, the check adds nothing. De la Fuente Marcos,
  NOT USEFUL for the blockers; the method is geometric sampling the project's Tisserand tools already cover. It remains
  target data for the dormant `#308`.

### 8.5 Linking two manifold sets at a common time, phase-free first (Vergaaij 2018)

- Method: unstable and stable manifold sets of the departure and arrival libration orbits, each built in its own
  three-body model, are transformed to one heliocentric inertial frame at a linkage time; the position and velocity
  errors are computed for all n x n pairs, and pairs with errors below 1e-3 and 1e-2 in Sun-Earth units are kept
  `[P p.8, steps 1 to 5]`. The free parameters (departure, link and arrival times, sail angle) are first searched in
  a problem made autonomous by removing the Earth-Moon phase, "because any relative orientation will occur once every
  synodic period" `[P p.9, sec. IV.III]`. Survivors seed a direct pseudospectral solver (PSOPT) `[P p.7]`.
- Hypothesis (7) tested. The step that `#650` lacks is the vector comparison: its edges carry `direction_data: "absent"`,
  and its design already says a CHEAP edge "is a candidate requiring vector-level follow-up before any claim"
  (`docs/notes/2026-07-19-650-transfer-network-design.md`, read). The vector comparison itself is the conditional flyby
  cost that `core/flyby.py` already computes (T2, T16). Vergaaij adds an ORDER of work: compare the two sets with phase
  removed first, then put the phase back only for pairs that pass. For `#650` that means comparing excess-velocity
  DIRECTIONS in the planet's frame for every pair of rows at the shared body, then searching epochs only for the pairs
  whose directions are compatible within the available bend `[I]`. The sail optimiser is not relevant.
- Verdict: POSSIBLY USEFUL, as the shape of `#650`'s registered vector-level follow-up, not as a separate proposal. What
  would decide it: whether the phase-free direction screen removes most CHEAP edges. If it does, the epoch search is
  small.

### 8.6 Statistical trim-manoeuvre cost from a stated navigation model (Rinker, Jacobson & Wood 1976)

- Method `[P pp.510-511]`. Errors at the manoeuvre time are mapped to the impact plane at encounter, z = R[Phi_11 dr +
  Phi_12 dV] (eq. 3). The guidance law y = -A z nulls the crosstrack and out-of-plane errors with the downtrack component
  set to zero (eqs. 6, 7), where y is the velocity change in impact-plane coordinates (eq. 4) and A is built from F = R
  Phi_12 R^T (eqs. 5, 7). The estimated state perturbation is taken as Gaussian with zero mean and covariance X^ (p.510),
  so the correction has covariance Y = A R [Phi_11 Phi_12] X^ [Phi_11 Phi_12]^T R^T A^T (eq. 8). The correction's magnitude has
  closed-form statistics from the two nonzero eigenvalues lambda1 <= lambda2 of Y: mean = (2 lambda2/pi)^(1/2) E[(1 -
  lambda1/lambda2)^(1/2)] / u for the duration at thrust acceleration u, with E the complete elliptic integral of the
  second kind (eq. 15), the variance (eq. 16), and a series for the distribution function G (eq. 17). The covariances are
  propagated through the manoeuvre with the execution-error covariance Q (eqs. 18 to 23). The authors point to Lee &
  Boain (1973) for the corresponding statistics of |Delta V| in a ballistic midcourse correction `[P p.511]` (not held).
- Bears on, `[I]`: the maintenance budget that `#907` and the `#890` and `#895` reviews ask for "from a stated
  navigation model", and the V2 question in section 4(v). The method is per encounter, which is the cadence an orbit with
  growth of 1e3 per flyby pair needs. Applied to a flyby chain, it gives one correction per leg, aimed at the next
  encounter's impact-plane coordinates. With an impulsive correction the same algebra holds: y is a velocity change,
  and u Delta t = |y|. The post-correction covariance (eqs. 21, 23) is then propagated through the flyby to the next
  leg's correction. The project has every STM this needs (`#895`'s variational Jacobian; the two-body STM). The
  navigation model (orbit-determination covariance at each correction time, execution error) is an input to state and
  pre-register, not something the paper or the project can derive.
- Control checked `[C]`. Table 2 prints, for five navigation-coast lengths, the mean and standard deviation of the
  manoeuvre duration and G at 0.4, 0.8 and 1.2 hours `[P p.511]`. Inverting eqs. 15 and 16 for lambda1/u^2 and
  lambda2/u^2 from each row's printed mean and standard deviation, then integrating the density of eq. 13, reproduces all
  15 printed G values within 0.001 (14 exactly to three decimals; the 4-day row gives G(0.4) = 0.615 against 0.614)
  (scratch calculation; scipy `ellipe` takes the parameter m = k^2 = 1 - lambda1/lambda2). This tests a statistics
  implementation against printed numbers. It does not test any navigation model.
- Verdict: USEFUL, as proposal Q1 below.

### 8.7 Recovery margin after a missed manoeuvre (Venigalla, Englander & Scheeres 2020)

- Method: the missed-thrust recovery margin is the longest forced coast from a point on the trajectory after which the
  spacecraft can still reach its terminal manifold; a "virtual swarm" of spacecraft, spawned at points along the nominal
  and each required to satisfy the terminal constraints after its forced coast, is co-optimised with the nominal
  trajectory, adding spacecraft at the worst violations (Algorithm 1) `[P text layer, abstract and "Virtual swarm method"
  section and Algorithm 1, PDF pp.1, 4-5; D marginal-papers-triage #1]`. Transcription is Sims-Flanagan `[D]`.
- Bears on, `[I]`: the operational risk of a flyby chain is a correction that is late or missed. The question "how
  long can the correction before encounter k slip before the next encounter cannot be recovered within a stated delta-v"
  can be answered post hoc from the same STMs as Q1, without co-optimisation. The virtual-swarm co-design belongs to a
  later design stage the catalogue does not reach.
- Verdict: POSSIBLY USEFUL, as a column computed inside Q1 (a slip margin per encounter), not as its own proposal. What
  would decide it: whether Q1's budgets are dominated by the correction's timing. If they are, the margin is worth
  reporting.

### 8.8 Conditional initial guesses for increasing robustness (Sinha & Beeson 2025)

- Method: compares a global search sampling from a fixed distribution with a "conditional global search" that seeds the
  more robust problem from solutions of a less robust one, a sequential ladder `[P abstract; D marginal-papers-triage #2]`.
- Bears on: seeding long chains from shorter ones. The project already does this in the form that matters (`#895`
  extended three cycles to six and twelve from the shorter solutions; P5's seven-cycle seed) `[I]`.
- Verdict: NOT USEFUL beyond current practice; the method is a sequential seeding rule the project already uses.

### 8.9 Teardrop hovering with one impulse per period (Fu, Peng, Gong & Shi 2025)

- Method: a deputy returns to a fixed relative position after one period of the 9:2 near rectilinear halo orbit with one
  impulse per revisit; natural-parameter continuation in the revisit distance with a linear predictor from STM entries
  and a least-squares solve, then fmincon `[D fu-peng-gong-shi secs. 1, 3, pp8-10]`. Impulses of order 1e-4 to 8e-2 m/s
  `[D sec. 1]`.
- Hypothesis (3), this paper: one correction per orbit period works for a near-stable reference. At multipliers of 1e6
  per period, an error that one correction per period must absorb has grown by that factor `[I]`. The cadence is wrong
  for the project's orbits. The continuation is natural-parameter with a fixed step and no fold handling `[D]`, weaker
  than what the project has. For taxi rendezvous with a cycler vehicle: the relative-motion set-up is the right
  kinematics for the last phase of a rendezvous, but the paper solves station keeping near a chief, not an approach from
  a planet. No held set-aside paper gives a taxi rendezvous method.
- Verdict: NOT USEFUL; one impulse per period cannot hold an orbit whose multiplier per period is 1e6, and the
  continuation is weaker than the project's.

### 8.10 Learned guidance for continuous thrust (Blender & Singh 2025; Hu, Yang, Li & Baoyin 2024; Singh & Junkins 2022; Li, Topputo & Baoyin 2019)

- Methods: gradient-boosted trees trained on sigma-point ensembles around an extremal bundle to map the belief state to
  control corrections (Blender & Singh); reinforcement learning per phase with a reachability term in the reward (Hu et
  al.); Gaussian-process prediction of costates from an extremal-field bundle (Singh & Junkins); a supervised network
  flying many-revolution orbit raising (Li et al.) `[D marginal-papers-triage #5, #6; background-papers-read-triage #2,
  #3]`. Not read at the source by me.
- Bears on: closed-loop guidance of a thrusting vehicle. The project's orbits are ballistic between impulsive
  corrections, and the open question is the size of the correction budget, not an on-board policy `[I]`. Blender & Singh's
  sigma-point ensemble is an unscented alternative to Q1's linear covariance. It would serve as Q1's nonlinear
  cross-check, but a direct Monte Carlo of the impulsive chain is simpler and has no training step `[I]`.
- Verdict: NOT USEFUL; the methods learn a control policy for continuous thrust, and the project's question is a budget
  for impulsive corrections.

### 8.11 Orbit determination and filtering (Sarang & Capannolo 2025; Peng & Bai 2021; Caldas & Soares 2024; Zhang, Li, Li, Zhang & Sang 2024)

- Methods: EKF, second-order EKF and UKF compared for angles-only relative navigation of a chaser near a Gateway-like
  target on a quasi-periodic orbit (Sarang & Capannolo); machine learning fused with an EKF to correct orbit predictions
  of tracked objects (Peng & Bai); a survey of machine learning in orbit estimation (Caldas & Soares); angles-only
  initial orbit determination with the Gooding algorithm and kernel-density solution selection (Zhang et al. 2024)
  `[D ml-orbit-estimation-triage; D 2026-10-03-digest-zhang-2024-iod-gooding-algorithm]`. Sarang & Capannolo state their
  filter set-up, including an initial covariance with 50 km per position axis and 1e-3 km/s per velocity axis
  `[P text layer, PDF p.10, "Filter Parameters"]`; I did not read their results.
- Bears on: Q1 needs an orbit-determination covariance at each correction time. None of these papers gives one for
  deep-space radiometric tracking of a vehicle on a flyby chain; Sarang & Capannolo's is for optical relative navigation
  near the Moon `[I]`.
- Verdict: NOT USEFUL for the blockers. Sarang & Capannolo are POSSIBLY USEFUL as a stated, cited covariance for an
  Earth-Moon row's Q1 run, if their results give steady-state errors (not checked).

### 8.12 A learned pre-filter for a cycler search (Ozaki et al. 2022)

- Method: Earth free-return legs (full-revolution, half-revolution, generic; Russell & Ocampo's classification) chained
  by Earth flybys, with an asteroid flyby inserted in each Earth-Earth leg; a network predicts the result of the
  Earth-asteroid-Earth block optimisation as a residual from the free return; beam search over sequences; survivors are
  re-optimised end to end by multiple shooting `[D ml-surrogate-trio-triage sec. 1.1-1.2]`. At least 7e6 samples were
  needed before the validation loss fell two orders of magnitude (Table 5); KKT-preserving "pseudo-asteroids" raise the
  data rate from 8.77 to 47.6 samples per second (Table 4); 10,000 block costs in 10 s against 1,140 s; the whole search
  took about 10 h against about 7 days; single-block error about 0.1 km/s, with errors accumulating along a chain `[D]`.
  The triage gave the prune rate as about 10k of 300k children surviving the screen per parent `[D]`; no false-negative
  rate is given in the triage.
- Hypothesis (6) tested. The economics are those of a combinatorial outer loop (15,340 asteroids times a deep tree) over
  one fixed, expensive block. The project's blockers are formulation, not compute: `#388` asked a ballistic question
  where the published answer is a cost; the moon closures matched speeds and not directions. A faster search over the
  same formulation does not touch either `[I]`. The June investigation found the project's exact cells cheap (31.6 ms to
  0.36 s) and preferred memoisation `[D ml-surrogate-investigation, Verdict A]`. The one place the economics could turn is
  P8: chains with resonant returns are combinatorial in the number of returns per moon and their types, which is
  Ozaki's architecture with a second moon in place of the asteroid `[I]`. A further way it could fool us: a region the
  pre-filter pruned was never solved, so it could not be entered in the negative-results registry as empty. Every "no
  chain found" from a filtered search would be conditional on the network.
- Verdict: POSSIBLY USEFUL, for P8 only. What would decide it: P8's chain count times the exact cost per candidate. If the
  product is under a few CPU-days, the network is not worth training. The Ozaki free-return chain decomposition is the
  Russell & Strange architecture P8 already uses.

### 8.13 Other learned surrogates (Zhang, Acciarini et al. 2026; Zhang, Michelotti et al. 2026; Viavattene & Ceriotti 2021; Leifsson et al. 2022; Wu, Sicard & Gadsden 2024; Silvestrini & Lavagna 2022)

- Methods: networks that approximate low-thrust optimal cost and reachability, with a Lambert solution as the best input
  feature (Zhang, Acciarini et al.); neural porkchop plots for low-thrust rendezvous (Zhang, Michelotti et al.); a ranker
  that prunes multi-target sequences, each survivor re-solved as a full optimal control problem (Viavattene & Ceriotti);
  adaptive sampling driven by network uncertainty (Leifsson et al.); reviews (Wu et al.; Silvestrini & Lavagna)
  `[D ml-surrogate-investigation; ml-surrogate-trio-triage; background-papers-read-triage #1]`.
- Bears on: the same pre-filter question as 8.12. The June note's finding stands. The best surrogate input "is itself a
  Lambert solution", so against a Lambert-bound inner loop the speedup disappears `[D ml-surrogate-investigation]`.
  Its governing rule also stands: a surrogate may prune or propose, never be evidence. The basin classifier that note
  deferred, "Build B", is the only learned tool with a use case. It remains blocked on persisted labels, as stated there.
- Verdict: NOT USEFUL; the methods learn costs of low-thrust optimal control, a baseline much more expensive than the
  project's Lambert and ballistic cells.

### 8.14 Dynamic programming with an error bound (Zhang, Guo, Wu, Baoyin, Li & Topputo 2025)

- Method: a fixed-sequence multi-flyby problem recast as an N-stage decision problem with state (epoch, velocity, mass)
  at each flyby, under the stated Markov property "only the velocity and mass of spacecraft at the flyby epoch influence
  subsequent trajectories"; states discretised; stage costs approximated; Bellman recursion; a bound J[P2] <= J[P0] + N
  eps_max with eps_max the largest single-stage error, estimated by random sampling, not derived `[D
  2026-10-03-digest-zhang-2025-global-optimality-multi-flyby-asteroid, eq. 42]`. A flyby is a position match with no
  velocity change `[D]`.
- Bears on: P9 (`#918`), whose limit in T16 is "complete on the grid only". The bound's algebra is the generic sum of
  stage errors in a Bellman recursion. It does not depend on what a flyby does to the velocity, so it carries over when
  the stage state is (epoch, incoming excess-velocity vector) and the stage cost includes the conditional flyby cost of
  T2 `[I]`. Two things weaken it for P9. First, eps_max is sampled, so the bound is an estimate, not a guarantee `[D]`.
  Second, a gravity assist makes the stage state a three-component excess-velocity vector, not a magnitude, so the grid
  grows with its cube; the paper's own fast case dropped the velocity from the state `[D]`. A periodicity condition links
  the last stage to the first, which DP handles by fixing the start state and repeating over start states `[I]`.
- Verdict: POSSIBLY USEFUL, as a component of P9: it would turn a P9 negative from "not found on this grid" into "no
  solution cheaper than the found one by more than N eps_max, as estimated". What would decide it: whether eps_max for the
  conditional flyby cost on a feasible velocity grid is small compared with the 1 to 10 m/s tiers of Russell & Ocampo.

### 8.15 Shape-based initial trajectories for a sail with gravity assists (Fan et al. 2025)

- Method: each heliocentric cylindrical coordinate fitted by an 8th-order Bezier curve with end coefficients fixed by
  the boundary conditions; free coefficients, flight times, flyby pericentre altitudes and pre-flyby velocities optimised
  by fmincon for minimum time, under thrust and a 70 km minimum flyby altitude `[D background-papers-read-triage #14]`.
- Bears on: nothing ballistic. Shape-based methods replace a thrust history; the project's legs are conics or integrated
  arcs `[I]`.
- Verdict: NOT USEFUL; a shape-based method exists to parametrise a thrust history, which ballistic legs do not have.

### 8.16 Entry corridors for returning crews (Putnam, Braun, Rohrschneider & Dec 2005)

- Method: for a given inertial speed at the atmospheric interface, the aerodynamic corridor is bounded by lift-down and
  lift-up trajectories, and the flyable corridor is narrowed by limits on peak heat rate, peak deceleration and
  integrated heat load `[P text layer, sec. "Entry Corridor Definition", PDF pp.4-5]`.
- Bears on: no blocker. For a cycler architecture it bears on the taxi's return: the entry speed implied by a row's
  Earth excess speed, against the corridor limits, would be a per-row figure of merit `[I]`.
- Verdict: NOT USEFUL for the blockers; POSSIBLY USEFUL as an architecture column if the owner wants one. That is a
  scope decision, not a method gap.

### 8.17 Lecture slides (Rickman, NESC)

- Introductory orbital mechanics for thermal engineers `[D rickman-nesc-slides-triage]`. Verdict: NOT USEFUL; textbook
  content with no method beyond what the project implements.

### 8.18 Proposals from this section

Format and probability definitions as in section 3: (a) a validated upgrade of existing rows; (b) an object worth an
adversarial review.

#### Q1. A per-encounter correction budget from a stated navigation model

- Goal: for any flyby-chain trajectory the project holds (the `#895` arcs, the `#890` orbit, the Russell & Strange
  rows once P4 runs, the Appendix C solutions of P1), compute the expected and 99th-percentile correction delta-v per
  cycle under a stated, pre-registered navigation model, and a slip margin per encounter (8.7).
- Method applied: Rinker, Jacobson & Wood 1976, eqs. 3 to 23, with the impulsive reading (y is the velocity change);
  one correction per leg at a stated lead time before each encounter, aimed at the impact-plane coordinates; covariance
  propagated through each correction (eqs. 21, 23) and through each flyby by the STM.
- Inputs: `search/titania_oberon_realeph_895.py` (variational STMs along the arcs), the two-body STM, `verify/
  turn_gate.py` (the flybys).
- To build: the covariance chain; the statistics of eqs. 15 to 17 (already checked against Table 2 `[C]`); a nonlinear
  Monte Carlo of the same chain; a navigation-model file with sources or labelled as convention (orbit-determination
  covariance at each correction, execution error as in the paper's example: 15 mrad pointing, 10 percent magnitude,
  `[P p.512]`, or a cited modern value).
- Positive controls: (i) Rinker et al. Table 2, all 15 G values from the printed means and standard deviations (passes
  now, scratch calculation); (ii) a negative control: on a stable orbit (a stable Casoliva row) the budget must be small
  and nearly independent of the number of cycles; on the `#895` arcs it must grow with the lead time as the STM predicts.
  No held paper prints a correction budget for a flyby chain, so there is no published end-to-end value to match. Lee &
  Boain (1973) would give the ballistic statistics and is not held.
- Pre-registrable criteria: the navigation model and lead times fixed before the run. Linear and Monte Carlo 99th
  percentiles agree within 20 percent, or the linear result is not reported. The result is reported per cycle with the
  model stated, never as a property of the orbit alone.
- How it could fool us: the linear covariance crossing a flyby, where the turn's slope in impact parameter is of order
  1/mu `[R MacKay p.2, section 4(v)]`. A covariance that is small in kilometres can still be large compared with the
  linear range there. Guard: the nonlinear Monte Carlo, and a check that the 3-sigma impact-plane error stays below a
  stated fraction of the flyby periapsis. Second: a budget computed with an optimistic navigation model. Guard: report
  two models (optimistic and conservative), both stated.
- Cost: 10 to 15 agent-hours; minutes of compute.
- (a) 20 percent: it gives the rows a cost column, but a level change needs the owner's ruling on V2 for unstable orbits
  (section 4(v), `#907`). (b) 35 percent: it is the missing item before `#895` can go to an adversarial review as a row
  candidate.

#### Q2. The `#890` orbit as a two-torus in a model with one added frequency, and the `#895` arcs compared with it

- Goal: compute the invariant two-torus that the `#890` periodic orbit should become when Oberon's eccentricity (with
  its J2 apsidal precession) is added to the planar circular four-body model, and measure how far the `#895` ephemeris
  arcs lie from it.
- Method applied: Villegas-Pinto et al. 2023 (one perturbation at a time); Rosales et al. 2021 (invariant curve of the
  stroboscopic map, multiple shooting over sections for unstable cases); Kumar et al. (bundle frame), as in T10.
- Inputs: `search/two_moon_periodic_890.py`, the `#895` force model (it already has circular-orbit bodies and
  kernel-driven ones), the invariant-circle corrector of `search/pertbp_strob_889.py` (planar elliptic problem only, to
  be generalised), `genome/qp_tori.py`.
- To build: a model with Oberon on a precessing ellipse; the multiple-shooting invariant-circle corrector (shared with
  P6); continuation in Oberon's eccentricity from 0.
- Positive controls: Rosales et al. 2021 Table 3 (printed rotation number, as in P6); the e = 0 limit, where the circle
  must collapse to the `#890` orbit (independently reproduced under `#895`).
- Pre-registrable criteria: invariance error below 1e-9 on a node ladder under an independent propagator; continuation
  reaches Oberon's eccentricity, or the stop is recorded (fold, gap, non-convergence) with the step that failed. Then,
  at each of the five `#895` epochs: the arc's flyby altitudes inside the torus's altitude envelope widened by a band set
  beforehand from the moons' circular-to-kernel position differences.
- How it could fool us: a "torus" converged with too few nodes; a resonance gap read as non-existence (Rosales: gaps are
  Cantorian and crossed by detours); the comparison band chosen after the fact. Guards: node ladder; gap rule from
  Rosales; band in the pre-registration.
- Cost: 20 to 30 agent-hours if P6's corrector exists first, about 10 more otherwise; hours of compute.
- (a) 0 percent (no row exists to upgrade). (b) 30 percent: a computed invariant object in a stated model, with the
  ephemeris arcs as its realisations, is a stronger basis for a `quasi_cycler` row than arcs alone.

#### Add-ons to existing proposals (not separate proposals)

- P7 (`#917`): a -1 multiplier monitor along the eccentricity continuation, controlled by Wang et al.'s printed
  crossings (8.1).
- P9 (`#918`): the stage-error bound of Zhang et al. 2025 (8.14).
- `#650`: the phase-free direction screen before the epoch search (8.5).
- P8: a learned pre-filter only if the chain count makes exact evaluation cost more than a few CPU-days, with the
  registry caveat (8.12).

#### Where they rank among P1 to P9

Q1 sits after P3 and before P4: it is cheap, it has a passed control for its statistics, and it removes a stated blocker
for `#895` and for any moon-tour row that P4 produces. Q2 sits after P6 and before P7: it needs the same corrector as P6
and should be built second, on the corrector P6 has validated.

#### The seven hypotheses, in one line each

1. Wang et al.: rejected as a source of a method the project lacks for asymmetric cyclers (the project has full-state correctors; the gap is
   `#884`'s symmetric three-body walk and T7's phase rule); kept as a control for a period-doubling monitor in P7.
2. Torus for the `#895` arcs: accepted in a qualified form. One torus in the ephemeris model is not the right object.
   A two-torus in a model with one added frequency is computable and testable (Q2).
3. Navigation and guidance papers: accepted for Rinker et al. (Q1, control passed) and, as a column inside Q1,
   Venigalla et al.; rejected for Fu et al. and the learned-guidance papers on method grounds. No held set-aside paper
   gives a taxi rendezvous method.
4. Co-orbital papers: the averaged model is not in Carruba et al.; the clone ensemble is possibly useful for
   quasi-cycler labels only; quasi-satellite orbits are not a two-body repeated-encounter class.
5. Shu & Lin: rejected; the expansion is local to the libration points.
6. Learned pre-filter: possibly useful for P8 only, decided by the chain count; the blockers are formulation, not
   compute.
7. Vergaaij: the vector linkage is already registered for `#650` and `core/flyby.py` has the cost; the phase-free first
   pass is the useful part.

No set-aside paper read for this section contradicts the project's code or a project belief.

## Appendix: what was read, and limits of this note

- Project state: README; spec sections 14, 16.4, 16.5; `data/OUTSTANDING.md` lines 1 to 1540 and the `#388`, `#378`,
  `#135` threads; the `#388` result notes of 2026-06-21 and 2026-06-23; the Russell continuation deep-dive; catalogue
  counts by a short script (392 rows; 240 with no validation level, 99 at V0, 43 at V1, 8 at V2, 2 at V3).
- Digests read in full or in their method and implication sections: Bradley & Russell; Gomez & Olle; Leiva & Briozzo
  2008; Brown et al. JAS; Casoliva (`#725`); Bellome; Li, Qiao & Li; Park & Howell 2026; Oshima 2024; Canales et al.
  The other digests of 2026-10-03 and 2026-10-04 were covered through the six source readers, who were told to read
  the digest and then the paper. I did not read every digest in full myself; where this note rests on a digest
  that I did not open, the tag is `[R]` with the reader's page, or the statement is absent.
- The six reader reports are not in the repository. If any `[R]` statement becomes load-bearing for a build, re-read
  that page first.
- Not done: no computation of any proposal; no check of the Appendix C blocks beyond the June transcription note; no
  check of which satellite kernels would be needed for P4 beyond listing the cached ones.
