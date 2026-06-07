# External novel-algorithms research survey (Task #143)

**Date:** 2026-06-07
**Scope:** docs-only literature survey. No source/test/data changed. Web-research
pass mapping four external algorithm families to live project frontiers. Every
citation below was verified against a search result or primary URL during this
pass; where a detail could not be read from the primary PDF (binary/encoded) it
is flagged as such rather than asserted.

## Frontier context (read first)

- **M-ED** (`data/OUTSTANDING.md` §"In progress"): the real-ephemeris ballistic
  corrector closes a *different, higher-V∞, powered* family than the sourced
  Jones VEM members. #110 (dense scan), #120 (3D inclination), #122
  (vector residual) all REFUTED as the fix — the blocker is the
  single-ellipse-per-leg topology + degenerate-basin selection, not scan density.
- **Free-return results** (`docs/notes/2026-06-07-russell12-freereturn-results.md`):
  the radial-crossing genome makes 8/12 Russell rows CLOSE-AND-MATCH on the
  circular model; 6 rows are now V1. But this is like-for-like circular-coplanar,
  NOT real-ephemeris (V3). ~190 Russell catalogue rows remain gapped.
- **Performance profile** (`docs/notes/2026-06-06-performance-profile.md`): solves
  are ephemeris-bound (state() = 87-89%); a dense scan of 2816 points/row is the
  current broad-search hammer. The maintenance (Aldrin E-M-E) solve is ~16.75 s
  astropy and the per-cycle ΔV is in the tens-of-m/s to km/s regime.

The recurring disease across all of these is the same: **our broad search finds a
degenerate / wrong basin, and dense grids do not cure it.** Three of the four
threads below attack exactly that; the fourth attacks the maintenance-ΔV side.

---

## Thread 1 — Monotonic Basin Hopping (MBH)

### The algorithm (~10 lines)

MBH is a stochastic global-search meta-algorithm for problems with many local
optima clustered into "funnels." Loop:
1. Start from an incumbent decision vector `x*` (best so far).
2. **Perturb**: `x' = x* + r`, where `r` is drawn componentwise from a
   **long-tailed distribution** (Cauchy, or — Englander & Englander 2014 — Pareto;
   the long tail is what lets it hop *between* basins, not just within one).
3. **Local-solve**: run an NLP local optimiser (SNOPT/SQP) from `x'` to the
   nearest feasible local optimum `x_loc`.
4. **Accept/reject**: if `x_loc` is feasible AND better than `x*`, set `x* = x_loc`
   (monotonic — only accept improvements; the "hop" comes from the perturbation, not
   from accepting uphill moves).
5. **Stop**: after a fixed wall-clock / iteration budget, or after N non-improving
   hops. Trajectory transcription is Sims-Flanagan multiple-shooting; flyby sequence
   can be searched by an outer integer GA.

The key insight vs a dense grid: MBH spends almost all effort *refining* with a
gradient solver and uses randomness only to escape the current funnel, so it
samples basins far more efficiently than uniform sampling at equal cost.

### Key citations (acquisition-grade)

- Englander, J. A., & Englander, A. C. (2014). *Tuning Monotonic Basin Hopping:
  Improving the Efficiency of Stochastic Search as Applied to Low-Thrust
  Trajectory Optimization.* 24th International Symposium on Space Flight Dynamics
  (ISSFD), Laurel, MD. PDF: https://www.issfd.org/ISSFD_2014/ISSFD24_Paper_S7-3_Englander.pdf
  — **the canonical perturbation-distribution paper** (Cauchy/Pareto long tails).
- Englander, J. A., Conway, B. A., & Williams, T. (2012). *Automated Mission
  Planning via Evolutionary Algorithms.* / Englander & Conway, *Automated Solution
  of the Low-Thrust Interplanetary Trajectory Problem*, JGCD (open via PMC:
  https://pmc.ncbi.nlm.nih.gov/articles/PMC5837074/ ). The EMTG architecture paper.
- Vavrina, M., Englander, J., et al. (2018). *Parallel Monotonic Basin Hopping
  for Low Thrust Trajectory Optimization.* NTRS 20180002379.
- (Track record on flyby tours: MBH+Sims-Flanagan is the engine behind multiple
  GTOC-winning entries and EMTG production mission studies; the GTOP database
  https://www.esa.int/gsp/ACT/projects/gtop/ is the standard benchmark set.)

### Open-source availability + licence

- **EMTG** (NASA Goddard) — https://github.com/nasa/EMTG , v9.02, C++. Licensed
  under the **NASA Open Source Agreement (NOSA) v1.3** (OSI-approved). Bundles
  Boost/CSPICE/GSL; **the only non-free dependency is SNOPT** (commercial SQP) —
  a real adoption blocker. Determines flyby sequence via integer GA + trajectory
  via NLP solver and MBH. NASA software catalog GSC-18459-1.
- **PyGMO / pykep** (ESA) — https://esa.github.io/pykep/ — ships an MBH
  user-defined-algorithm (`pygmo.mbh`) wrapping an arbitrary local solver, MPL-2.0,
  pure-Python-callable. **Far easier to adopt than EMTG** for our stack (no SNOPT;
  we can wrap our own scipy `least_squares`/SLSQP as the inner solver). pyoptgra
  (GPL-3 / ESCL) provides the OPTGRA SQP if we want a near-linear-constraint solver.

### Fit × effort verdict for OUR codebase

- **What it replaces/augments:** the dense epoch×branch grid in `search/scan.py`
  (`scan_parallel`, currently 2816 pts/row brute force) and the seed ladder in
  `search/seed_ladder.py`. MBH would wrap our existing `ballistic_correct`
  (the local solve) and hop between seeds instead of gridding them.
- **Why it fits the disease:** #110/#120/#122 proved the failure is basin
  *selection*, not density. MBH is the trajectory community's standard answer to
  exactly this. Around the free-return seeds it could plausibly find the
  bend-feasible Jones basin that the grid+corrector never lands in — IF such a
  basin exists for the single-ellipse topology (the #122 evidence suggests it may
  not, in which case MBH would confirm the topology blocker faster/cheaper than a
  denser grid, which is itself valuable negative science).
- **Scope estimate:** MEDIUM. A first MBH loop over `ballistic_correct` with a
  Cauchy perturbation on `(t0, ToFs)` and monotonic accept is ~150-250 LOC +
  tests; no new heavy deps (reuse scipy as inner solver, or optionally pull
  pykep's `pygmo.mbh`). Parallelisable on the existing ProcessPool.
- **Uncertainty:** whether MBH escapes to a *qualitatively different* topology is
  unproven — MBH hops between basins of a fixed transcription, so if the
  single-ellipse transcription has no bend-feasible Jones basin, MBH cannot
  conjure a multi-arc one. It is a basin-selection cure, not a topology cure.

---

## Thread 2 — Pseudo-arclength / family continuation

### The algorithm (~10 lines)

Given ONE converged family member `X0` satisfying `F(X, λ)=0` (F = the
closure/periodicity residual, λ = a family parameter, e.g. Jacobi constant, a
synodic-repeat count, or a geometric descriptor):
1. **Tangent predictor:** compute the null-space tangent `Ẋ` of the augmented
   Jacobian `[∂F/∂X | ∂F/∂λ]` (the direction along the family).
2. **Step:** `X_pred = X0 + Δs · Ẋ` for a pseudo-arclength step `Δs` (parameterise
   by arclength `s`, NOT by λ — this is what lets you turn folds where λ reverses).
3. **Corrector:** Newton-solve the augmented system
   `{F(X,λ)=0 ; (X−X0)·Ẋ + (λ−λ0)·λ̇ − Δs = 0}` (the extra row pins arclength).
4. **Adapt:** grow/shrink `Δs` on corrector convergence; detect **folds**
   (tangent λ̇ → 0) and **bifurcations** (Jacobian rank drop) by monitoring the
   determinant/eigenvalues.
5. Repeat → sweep the whole family from the single seed; no re-search per member.

This is the standard "free-variable / constraint differential-correction +
pseudo-arclength continuation (PAC)" used by the Howell group for periodic-orbit
families, and by Russell & Ocampo specifically to walk Earth-Mars cyclers from
analytic circular-coplanar seeds to full-ephemeris solutions.

### Key citations (acquisition-grade)

- **Russell, R. P., & Ocampo, C. A. (2006). *Optimization of a Broad Class of
  Ephemeris Model Earth-Mars Cyclers.* Journal of Guidance, Control, and Dynamics,
  29(2), 354-367.** — THE directly-on-target paper: a continuation algorithm
  transitions infinitely-continuous circular-coplanar solutions to finite-duration
  accurate-ephemeris solutions; reports nine parent cyclers (incl. Aldrin) with a
  finite-duration member needing **< 1 m/s total maneuver over seven full cycles.**
  This is the exact circular→ephemeris bridge our V1→V3 gap needs.
- Russell, R. P., & Ocampo, C. A. (2004). *Systematic Method for Constructing
  Earth-Mars Cyclers Using Free-Return Trajectories.* JGCD 27(3), 321-335.
  (The free-return-genome lineage already in our catalogue.)
- Howell / McCarthy group on PAC for periodic-orbit families: McCarthy, B. P., &
  Howell, K. C. (2021). *Quasi-Periodic Orbits in the Sun-Earth-Moon BCR4BP*,
  AAS 21-270 (https://engineering.purdue.edu/people/kathleen.howell.1/Publications/Conferences/2021_AAS_McCHow.pdf)
  — canonical free-variable/constraint + PAC formulation, fold/bifurcation handling.
- General method reference: Keller (1977) pseudo-arclength continuation; standard
  predictor-corrector treatment.

### Open-source availability + licence

- No single drop-in for *cycler* continuation. Generic numerical-continuation
  engines exist: **AUTO-07p** (periodic-orbit continuation, GPL), **pyDStool**
  (BSD), **COCO** (MATLAB). The Howell-group astrodynamics tooling is not openly
  released; the method is standard enough to implement directly (the corrector is
  our existing differential corrector + one augmented arclength constraint row).

### Fit × effort verdict for OUR codebase

- **What it replaces/augments:** the per-row re-search. We already have a
  converged member generator (`search/free_return.py` /
  `free_return_v1.py` close individual rows; `ballistic_correct` is the corrector).
  Continuation would let us close ONE descriptor-family member, then **predict +
  correct** along the family (sweeping the ~190 gapped Russell rows that belong to
  the same descriptor family) instead of seeding each from scratch.
- **Why it fits the disease:** continuation stays *attached* to a known-good basin
  by construction (small arclength steps + Newton corrector), so it sidesteps the
  basin-selection failure entirely — it never re-searches, it tracks. As gapped
  rows gain sourced data, family members can be compared against the sourced rows.
- **Scope estimate:** MEDIUM-HIGH. Augmented-Jacobian tangent + arclength
  corrector wrapper around `ballistic_correct` is ~200-300 LOC; fold/bifurcation
  detection and a sensible family-parameter choice (synodic-repeat count vs a
  geometric descriptor) are the design risk. Needs the corrector to expose its
  Jacobian (scipy `least_squares` gives it; check it is accessible).
- **Uncertainty:** which parameter parameterises the Russell descriptor families
  cleanly is non-obvious — the #137 finding that `6.44Gg3`'s aphelion+transit and
  its V∞ describe *different* ellipses warns that the catalogue descriptor axis may
  not be a clean continuation parameter. Continuation along the WRONG parameter
  re-creates the off-anchor problem. Start with the Russell&Ocampo 2006 parameter.

---

## Thread 3 — Primer vector theory

### The algorithm (~10 lines)

Primer vector `p(t)` = the velocity-adjoint (costate) from Lawden's first-order
optimality conditions for impulsive transfers. Necessary conditions for a
ΔV-optimal impulse sequence:
1. `|p(t)| ≤ 1` everywhere, with `|p| = 1` exactly at each impulse time.
2. At each impulse, the thrust direction is `p̂` (impulse parallel to primer).
3. At interior impulses, `dp/dt` is continuous and `d|p|/dt = 0` (the primer
   touches the unit circle tangentially).
4. At the endpoints, `ṗ·p̂` conditions fix whether to add a coast.

**Diagnostic use (Lion & Handelsman 1968):** propagate `p(t)` for a *given*
(non-optimal) trajectory; if `|p(t)| > 1` anywhere, the trajectory is NOT optimal,
and the time/location of the peak tells you **where to add a midcourse impulse**.
A first-order formula predicts the new interior-impulse position and the resulting
ΔV reduction (Jezewski & Rozendaal iterate this to convergence). If `|p| < 1`
everywhere with unit endpoints, the impulses are optimally placed.

### Key citations (acquisition-grade)

- Lawden, D. F. (1963). *Optimal Trajectories for Space Navigation.* Butterworths.
  (Foundational — defines the primer vector.)
- Lion, P. M., & Handelsman, M. (1968). *Primer Vector on Fixed-Time Impulsive
  Trajectories.* AIAA Journal 6(1), 127-132. DOI 10.2514/3.4452. (The
  add-an-impulse / improve-a-non-optimal-trajectory extension — the diagnostic.)
- Jezewski, D. J., & Rozendaal, H. L. (1968). *An efficient method for calculating
  optimal free-space N-impulse trajectories.* AIAA Journal 6(11). (Iterative
  impulse-adding algorithm.) See also Jezewski (1980), *Primer vector theory
  applied to the linear relative-motion equations*, Opt. Control Appl. Methods.
- Guzman, J. J., Mailhe, L. M., Schiff, C., Hughes, S. P., & Folta, D. C. (2002).
  *Primer Vector Optimization: Survey of Theory, New Analysis and Applications.*
  IAC-02-A.6.09, 53rd International Astronautical Congress, Houston. NTRS
  20030032208 — a clean modern survey incl. where the linearised theory *fails*
  (singularities along arcs), directly relevant before we trust it on cyclers.
- Prussing, J. E. (2010). *Primer Vector Theory and Applications*, Ch. 2 in Conway
  (ed.), *Spacecraft Trajectory Optimization*, Cambridge UP (standard textbook).

### Open-source availability + licence

- No widely-used standalone primer-vector library. It is short to implement: it
  needs the state-transition matrix (STM) along each arc (to propagate the
  adjoints) plus the unit-magnitude diagnostic — both within reach of our existing
  Lambert/Kepler propagation. GMAT (Apache-2.0) has related optimal-control
  tooling but not a packaged primer diagnostic.

### Fit × effort verdict for OUR codebase

- **What it augments:** the maintenance / TCM work in `search/maintain.py`
  (`optimise_aldrin_maintenance_dv`, the per-cycle flyby-ΔV sum) and the #134
  horizon-chaining. Primer analysis answers a question we currently do NOT ask:
  **are our per-cycle maneuvers optimally placed/timed?** A `|p|>1` excursion on a
  computed maintenance schedule is direct evidence an extra/relocated TCM would
  cut the total ΔV.
- **Why it fits:** it is a *diagnostic on top of* an existing solution, not a new
  search — low blast radius, and it produces a sourced-optimality witness
  (a published necessary condition), which suits our evidence-chain culture.
- **Scope estimate:** MEDIUM. STM propagation along the maintenance arcs +
  `|p(t)|` evaluation + the add-impulse predictor is ~200-300 LOC; the corrector to
  *re-optimise* with the new impulse reuses `maintain.py`. First deliverable can be
  pure-diagnostic (compute `|p|` profile for the current Aldrin schedule, report
  whether it violates unity) before any re-optimisation.
- **Uncertainty:** the Guzman survey flags that linearised primer theory has
  singularities and can fail on long, multi-rev arcs — cyclers are exactly long
  multi-rev arcs, so validate the STM/adjoint propagation on a short known case
  first. Whether the current ~tens-of-m/s..km/s ΔV is reducible is unknown until
  the `|p|` profile is computed; do that diagnostic before promising a number.

---

## Thread 4 — STOUR / STOUR-LTGA + EMTG broad-search architectures

### The architecture (~10 lines)

STOUR (Satellite Tour Design Program; JPL-origin, extended by Longuski/Purdue) is
a **patched-conic broad-search** pathfinder for multiple-gravity-assist (MGA)
trajectories:
1. Fix a flyby body sequence (e.g. E-V-E-M).
2. **Grid search** over launch date × launch V∞ (and, per leg, transfer revs).
3. For each grid node, solve each leg's Lambert/conic arc and **propagate the
   flyby**: at each intermediate body, enforce V∞-magnitude continuity and compute
   the required bend; **prune** any node whose required bend exceeds the
   gravity-assist limit (r_p < body radius) — this is the feasibility filter.
4. Surviving paths are the candidate trajectories within the search box; tabulate
   their V∞ at each body, dates, and ΔV.
5. **STOUR-LTGA** (Petropoulos & Longuski) swaps the ballistic Lambert legs for a
   2-D exponential-sinusoid *shape-based* low-thrust arc, enabling broad low-thrust
   MGA search at the same grid-search level.

EMTG is the modern descendant: same broad-search spirit but Sims-Flanagan
transcription + MBH (Thread 1) + integer-GA sequence search instead of a fixed
sequence + uniform grid.

### Key citations (acquisition-grade)

- Petropoulos, A. E., & Longuski, J. M. (2004). *Shape-Based Algorithm for the
  Automated Design of Low-Thrust, Gravity Assist Trajectories.* Journal of
  Spacecraft and Rockets, 41(5), 787-796. DOI 10.2514/1.13095. (STOUR-LTGA / the
  exponential-sinusoid shape model.)
- Patel, M. R., Longuski, J. M., et al. — STOUR Earth-Mars cycler / MGA studies
  (the STOUR lineage our sources cite; see the Purdue AAC publications list).
- **Direct cross-check target:** *Fast Mars Free-Returns via Venus Gravity Assist*,
  AIAA 2014-4109, Purdue AAC
  (https://engineering.purdue.edu/AAC/wp-content/uploads/2012/09/FastMarsFreeReturnsviaVenusGravityAssist-AIAA-2014-4109.pdf)
  — a STOUR-style broad search of E-V-M free returns; tabulates V∞ per body and
  launch windows. **This is the same topology family as our Jones VEM frontier
  (#110), so its published V∞ table is an external cross-check for our blocked
  EMEVVE/MEEVEM survey.** (PDF was binary on fetch — acquisition item to pull text.)
- Vasile, M., & De Pascale, P. (2006). *Preliminary Design of Multiple
  Gravity-Assist Trajectories.* J. Spacecraft & Rockets / arXiv 1105.1822 — an
  alternative broad-search pruning architecture worth comparing.

### Open-source availability + licence

- STOUR itself is **not open-source** (JPL/Purdue internal). The open descendants
  are **EMTG** (NOSA-1.3, needs SNOPT) and **pykep/PyGMO** (MPL-2.0 — includes MGA
  and MGA-1DSM problem builders + the GTOP benchmark set), the practical route to a
  STOUR-like broad search in Python.

### Fit × effort verdict for OUR codebase

- **What it augments:** the Forge novelty loop (`scripts/forge_novelty_run.py`,
  `data/OUTSTANDING.md` Forge Phases 4+5) — STOUR's patched-conic-tree + bend-prune
  pruning is a more principled broad-search than our current epoch×branch grid, and
  its published V∞ tables (esp. AIAA 2014-4109) give us a **sourced cross-check**
  for the Jones VEM V∞ values where our corrector finds only the high-V∞ basin.
- **Why it fits:** our `hunt_vem_ballistic.py` already *is* a STOUR-like grid +
  bend-feasibility prune (we report "0 bend-feasible"). Reading the STOUR papers
  tells us (a) whether published E-V-M free returns achieve low V∞ at all (if yes,
  with what topology — likely multi-arc, confirming the M-ED front-runner), and (b)
  exact numbers to validate our pipeline against. The biggest near-term value here
  is **validation/cross-check, not a new engine.**
- **Scope estimate:** LOW for the cross-check (acquire AIAA 2014-4109, compare its
  V∞/date table to our windows — a diligence task, no code). HIGH if we rebuild the
  broad-search engine; not recommended over adopting pykep's MGA builders if we
  want that.
- **Uncertainty:** the published free-return V∞ values may come from a *different
  topology* (multi-arc, real eccentricity) than our single-ellipse corrector — in
  which case the cross-check confirms the topology gap rather than validating our
  numbers. Either outcome is informative.

---

## Ranked recommendation list (what to implement first, and why)

1. **Primer vector DIAGNOSTIC on the existing maintenance schedule (Thread 3).**
   Lowest blast radius (a read-only diagnostic on a solution we already compute),
   directly answers an unasked question (#134: are our TCMs optimally placed?),
   and produces a sourced-optimality witness fitting our evidence culture. If
   `|p|>1`, we have evidence-backed ΔV to recover; if `|p|≤1`, we have a published
   proof our placement is already optimal. Either way it is a clean, citable result.

2. **MBH wrapper around `ballistic_correct` (Thread 1), tried first on the
   free-return seeds.** The trajectory community's standard cure for our exact
   degenerate-basin disease, which #110/#120/#122 proved is a basin-*selection*
   (not density) problem. Medium scope, no heavy deps (reuse scipy inner solver or
   pull `pygmo.mbh`). Honest caveat: it cures basin selection within a fixed
   transcription, so if the Jones bend-feasible basin doesn't exist for the
   single-ellipse topology, MBH will *confirm that faster than a denser grid* —
   still a positive outcome.

3. **STOUR cross-check of the Jones VEM frontier (Thread 4), as diligence.** No
   code: acquire AIAA 2014-4109 (Fast Mars Free-Returns via Venus GA) and compare
   its published per-body V∞ / windows against our #110 EMEVVE/MEEVEM survey.
   Cheap, and it directly tests whether the sourced low-V∞ family is even
   single-ellipse-representable — decisive for whether to invest in multi-arc.

4. **Pseudo-arclength family continuation (Thread 2).** Highest long-term payoff
   (sweep ~190 gapped Russell rows from one closed member; the V1→V3 bridge that
   Russell & Ocampo 2006 literally demonstrate for this exact problem), but highest
   scope and the family-parameter-choice design risk (the #137 `6.44Gg3` warning).
   Do it AFTER #2 establishes a reliable single-member generator and AFTER the
   Russell&Ocampo 2006 parameterisation is in hand.

### The single best near-term win

**The primer-vector diagnostic (item 1).** It is the smallest change, touches only
the maintenance layer (orthogonal to the M-ED/topology fight), needs no new
dependency, and converts the open #134 question ("are our per-cycle maneuvers
optimally placed?") into a citable yes/no with a recoverable-ΔV estimate — exactly
the kind of sourced, low-risk result this codebase rewards.

### Acquisition-list additions (for the #116 channel)

1. **Russell, R. P., & Ocampo, C. A. (2006). Optimization of a Broad Class of
   Ephemeris Model Earth-Mars Cyclers. JGCD 29(2):354-367.** — the
   circular-coplanar→ephemeris continuation method; the V1→V3 bridge for Russell rows.
2. **Englander & Englander (2014), Tuning Monotonic Basin Hopping (ISSFD24
   S7-3).** — open PDF already located; pull for the Cauchy/Pareto perturbation
   spec before implementing MBH.
3. **Fast Mars Free-Returns via Venus Gravity Assist, AIAA 2014-4109 (Purdue
   AAC).** — STOUR-style E-V-M free-return V∞/window tables; the external
   cross-check for the Jones VEM frontier. (PDF located but binary on fetch — needs
   text extraction.)
4. **Guzman, Mailhe, Schiff, Hughes & Folta (2002), Primer Vector Optimization:
   Survey of Theory, New Analysis and Applications, IAC-02-A.6.09 (NTRS
   20030032208).** — modern primer survey incl. failure/singularity cases; read
   before trusting primer theory on long multi-rev cycler arcs. (PDF binary on
   fetch — needs extraction.)
5. **Lion, P. M., & Handelsman, M. (1968), Primer Vector on Fixed-Time Impulsive
   Trajectories, AIAA J. 6(1):127-132, DOI 10.2514/3.4452.** — the add-an-impulse
   diagnostic that the maintenance work would implement.

> Caveat on uncertainty: four primary PDFs (the primer survey, the Mars-free-return
> STOUR paper, and two others) returned binary/encoded content on automated fetch,
> so their *internal numeric tables* were NOT read in this pass — their citations
> are confirmed via search metadata, but the V∞ tables and exact formulae must be
> verified on acquisition before any number from them is quoted as sourced.
