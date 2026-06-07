# Ceriotti 2010 Glasgow PhD — MGA global-optimisation methods (algorithm mine)

Mined 2026-06-07 (Task #143 algorithm pass). Method chapters only (Ch.2 trajectory
models §2.3 position formulation; Ch.3 Incremental Pruning in full; Ch.5 ACO-MGA in
full; Ch.6 conclusions intro). Background/literature (Ch.1, Ch.4 Laplace+BepiColombo
real-case results, appendices A/C/D/E) read selectively or skipped — see §7.

**Source (cite exactly, no file path):**
Ceriotti, M., "Global Optimisation of Multiple Gravity Assist Trajectories," Ph.D.
thesis, Department of Aerospace Engineering, University of Glasgow, 2010. 304 pp.
(Supervisor: M. Vasile.)

> Clean digital typeset; equations, algorithms and tables read unambiguously. Vision
> read of front matter (pp.v, vii-xxii, xxix-xxx), §2.3.4 (pp.58-62), Ch.3 (pp.63-118),
> Ch.5 §5.1.4-5.4 (pp.182-212), Ch.6 intro (p.213).

---

## 1. The two contributions in 3 lines each

**(A) INCREMENTAL PRUNING (Ch.3).** Decompose the MGA decision vector into ordered
levels (one per leg/swing-by); exploit that the total objective is a *sum of
non-negative partial objectives each depending only on levels 1..i* (Eq.3.3-3.4), so
by Bellman's principle of optimality the partial objective is a lower bound for the
whole (p.67). Search level 1 for the *feasible set* (all points with partial objective
below a pruning threshold `f̄_i`, Eq.3.5), bound it by axis-aligned **boxes**, and
search only inside those boxes at level 2 — pruning the rest. Iterate leg by leg with
**back-pruning** (removing earlier-level nodes that become infeasible later). The
feasible-set *finder* at each level is a pluggable optimiser (Multi-Start, modified
MBH Alg.3.1, or modified MACS Alg.3.2); boxing is one of 4 clustering methods.

**(B) ACO-MGA (Ch.5).** Recast the *whole* MGA design (planetary sequence + leg types
*together*) as a planning/scheduling problem on an acyclic tree, and solve with an
Ant-Colony-inspired search (Alg.5.3-5.7). Each leg = 5 discrete vars (1 DSM-magnitude
index, 2 rev counts, 2 binaries) + target planet; ants build a solution leg-by-leg
choosing probabilistically from a per-leg pheromone vector. Two list structures
replace classical ACO's edge-pheromone: a **feasible list** (complete solutions, used
to deposit pheromone weighted by 1/objective) and a **taboo list** (partial solutions
that proved infeasible, used to forbid re-exploration). The phasing (scheduling) of
each leg is solved internally by a 1-D Brent root-find of Δθ(λ)=0.

---

## 2. The algorithms stated precisely

### 2.1 Incremental-pruning formalism (Ch.3.2, pp.65-67)
- Partition `x = [x_{L,1} … x_{L,N_L}]`; cumulative vectors `x_i = [x_{L,1}…x_{L,i}]`,
  domains `D_i = ∏_{k=1}^i D_{L,k}` (so `D_{N_L} = D`).
- **Separable objective** (Eq.3.3): `f(x) = Σ_{i=1}^{N_L} f_i(x_i)`, each `f_i ≥ 0`
  (Eq.3.4) and depending only on levels 1..i. "We will call each function `f_i`
  *partial objective function* for level i. Thus, the total objective is built-up level
  by level, incrementally." (p.66)
- **Pruning threshold + feasible set** (Eq.3.5-3.6):
  `D̄_{L,i} = { x_{L,i} ∈ D_{L,i} | f_i([x_{i-1} x_{L,i}]) ≤ f̄_i }`, then
  `D̄_i = D̄_{i-1} × D̄_{L,i}`. Boolean **pruning criterion** `Φ_i(x_i) = (f_i ≤ f̄_i)`
  (Eq.3.7). Justification verbatim: "According to Bellman's principle of optimality
  [108], if all the partial solutions from 1 to i are optimal, `f_i` is a lower bound
  for `f_j`, when j>i, and for the whole objective function f." (p.67)
- **Pruning criterion need not be the ΔV itself.** For a swing-by level the ΔV is "not
  a good pruning criterion. Instead, by considering the aim of this particular
  swing-by, that is to increase the semimajor axis of the transfer orbit, it is
  possible to prune the second level by requiring that the Venus swing-by increases
  the semimajor axis by at least 27,986,077 km" → `Φ_2 = (Δa > 27986077)` (p.77).
  Third (DSM) level: `Φ_3 = (Δv_2 < 0.5 km/s)`. **Numerical pruning values "set from
  the various data in literature about similar trajectories [41,74]"** (p.77).
- **Complexity** (Ch.3.3.2, p.80): no-pruning grid scan needs `∏_{i=1}^n p_i`
  evaluations; the incremental process needs `Σ` of partial-grid products which can be
  *more* function evals than the full scan in the worst (no-pruning) case, but each
  early-level eval is cheaper, so wall-time can still be lower. (Table 3.22 / Table C.6
  give per-level eval costs: level-1 eval 2.33 ms vs level-2 3.68 ms on a P4 3 GHz.)

### 2.2 The feasible-set finders (Ch.3.4.1, pp.81-85)
- **Multi-Start (MS):** N uniform samples → take best m → run `fmincon` (MATLAB SQP)
  from each, but **stop the local solve as soon as a point falls below the pruning
  threshold** (not at the minimum) — "we are not interested in fully converging" — to
  save evals and keep feasible points *sparse* (spread over the feasible set, not piled
  at the minimum) (p.82). This sparse-not-converged trick is the key to mapping a point
  optimiser onto *set* identification.
- **Modified MBH (Algorithm 3.1, p.83):** standard MBH (perturb in neighbourhood
  `N(x*, ρ_i)` of incumbent, local-solve, keep best) but **stores *all* candidate
  points satisfying `Φ_i`** into `L_feas` (not just the best — basic MBH "would be
  unusable to explore the feasible set once one point is identified", p.84), plus a
  **restart after `n_trials,max` non-improving trials** to avoid covering only one
  disconnected component. `ρ_i = 10%` of the variable range in the tests.
- **Modified MACS (Algorithm 3.2, p.85):** a derivative-free Multi-Agent Collaborative
  Search (population of agents w/ individualistic + communication operators); archives
  all feasible samples; restart partitions `D_i` and picks a sub-domain by a Pareto
  dominance index over `ψ_{D_q} = [−V_{D_q}, n_feas,Dq / n_feas,Di]` (Eq.3.15-3.16).
- **Boxing** (Ch.3.4.2, pp.86-89): 4 methods — (1) envelope each feasible point in a
  fixed box; (2) boxes on a grid; (3) grid + wrap adjacent boxes; (4) **Mean-Shift
  clustering** (single bandwidth parameter, [113]). Method 4 used in the Ch.4 real
  cases; method 3 used in the Ch.3.7 comparison.

### 2.3 A reusable partial-objective for resonant/gravity legs (Eq.3.17, p.107)
For any leg whose downstream swing-by needs an outgoing velocity parallel to the
planet's (to brake or accelerate — "aphelion rising / perihelion lowering gravity
manoeuvres"):
`f_i = β (v_{θ,i}^2 + v_{h,i}^2) / v_{r,i}^2 + Σ_{l=1}^i Δv_l`, with β = 1 km/s.
This maximises the radial component of the incoming relative velocity before the next
swing-by while minimising DSM. Derived empirically from Fig.3.18 (the best EEM
solutions all have near-radial incoming relative velocity). **General to two classes
of MGA transfer**, per the author.

### 2.4 ACO-MGA detail (Ch.5.2, pp.185-194)
- **Solution coding** (§5.2.1, Fig.5.10): vector of `2·n_legs` positive integers;
  odd entries index the target planet (from a per-leg ordered list `q_{P,i}`), even
  entries index a row of the Cartesian-product matrix `G_i` (Eq.5.8) of the 5 discrete
  leg parameters `{m_DSM, n_rev,1, n_rev,2, f_{p/a}, f_{1/2}}` (Eq.5.7). `t_0` and
  launch angle `φ_0` are **pre-fixed and removed from the search** (rationale: scan the
  launch window separately once a per-date solver exists, p.186).
- **Pheromone update** (Eq.5.9): for each candidate planet j at leg i,
  `τ_{P,i,j} ← τ_{P,i,j} + (1/y_l) w_planet` summed over feasible-list solutions whose
  partial sequence matches; selection prob `Pr = τ / Στ`. **Pheromone reset to all-ones
  every time an ant is at leg i** and **is not capped at 1** (differs from standard
  ACO, p.190). `w_planet`/`w_type` are the learning rates; the two-step protocol
  (§5.3, Eq.5.10) runs step 1 with `w=0` (pure random sampling to seed the feasible
  list) then step 2 with `w = ŵ·y_est` (intensify around feasible solutions).
- **Taboo list** (§5.2.2): per-leg matrix of partial solutions found infeasible at leg
  i; "Taboo lists have no equivalent in classical ACO" (p.194). When all transfer types
  at a leg are taboo, the pheromone vector is all-zeros and the ant's branch is dropped
  (Alg.5.6 line 13-14).
- **Phasing solve** (§5.1.4, pp.179-182): each leg's scheduling = find zeros of
  `Δθ(v_0)` and `Δθ(r_ps)` by **Brent's method** [119] from a supplied set of starting
  points; resonant transfers give *multiple* zeros (discontinuous Δθ), so all zeros are
  carried forward as a *tree* of trajectories (Alg.5.2). This is the integrated
  sequence+trajectory step that the two-level incremental approach does separately.

---

## 3. Comparison data (directly informs the #143 roadmap)

### 3.1 Incremental vs all-at-once, and the stochastic-optimiser bake-off

**EEM (Earth-Earth-Mars, single GA, 2 levels; best-known ΔV = 2.908 km/s).** All-at-once
success rate (% of 100 runs reaching ≤ 2.958 km/s), Table 3.15 (p.114):

| Solver | 20k eval | 40k | 80k | 160k |
|---|---|---|---|---|
| DE  (< 2.958) | 0% | 7% | 27% | 27% |
| MBH (< 2.958) | 1% | 5% | 18% | 41% |
| **MS (< 2.958)** | **22%** | **32%** | **52%** | **67%** |
| DIRECT (deterministic, km/s) | 4.317 | 4.317 | 3.822 | 3.809 |
| MCS (deterministic, km/s) | 3.840 | 3.840 | 3.840 | 3.812 |

**EVVMeMe (Earth-Venus-Venus-Mercury-Mercury, 14-D, BepiColombo-class; ESA ref
9.467 km/s).** Success (% of 100 runs beating 9.467 km/s), Table 3.19 (p.117):

| Solver | 200k | 400k | 800k | 1.6M |
|---|---|---|---|---|
| DE  (< ESA) | 16% | 16% | 21% | 16% |
| MBH (< ESA) | 4% | 3% | 7% | 16% |
| MCS (km/s, deterministic) | 14.35 | 13.05 | 13.05 | 12.01 |

**EVM (Earth-Venus-Mars, single GA; best ΔV = 2.9818 km/s, Table 3.7).** % of 200 runs
finding ≤ 3 km/s, Table 3.8 (p.102): DE 6.5/5/7%, MS 2.5/3/3%, PSO 2/2.5/7.5%; DE
beats MCS 100% but only 71-78% (PSO) to 99.5% (DE) beat DIRECT.

**Key verbatim findings for our optimiser choice:**
- "Table 3.15 ... suggests that sophisticated global search methods, such as DE and
  MBH, are not the right choice; in particular, **DE is the worst performing
  algorithm**. The reason is the fast convergence of DE with the selected settings."
  (p.112) "**A simple Multi-Start algorithm, instead, can yield better performance
  provided that the local optimisation algorithm converges fast.**" (p.113)
- "if the objective function is globally non-convex, i.e. presents multiple similar
  funnel structures, **MBH may not be effective** and DE could quickly converge but
  within a single funnel" (p.113). MGA cost surfaces are exactly such multi-funnel
  surfaces (Fig.3.22: many near-optimal solutions with very different vectors).
- DE settings used: F=0.8, Cr=0.75, strategy `best`, pop 90-140. MBH: ρ=10% range.
  **These are the settings under which DE under-performs — a tuning caveat, not a
  blanket result** (Englander & Englander 2014 in our survey tunes MBH's perturbation
  for exactly this reason).

**Incremental pruning payoff (Table 3.16, p.114, EEM, 100 runs):** inclusion of best
solution 100% (both MACS & MBH finders); pruned space 90.43% (MACS) / 93.51% (MBH);
avg boxes 8.64 (MACS) / 15.9 (MBH); avg coverage 88.58% / 72.21%. After pruning, the
*same* DE/MBH/MS re-run on the reduced space have far higher success (Table 3.17):
MS reaches 100% at 40k evals on the box containing the reference. Net wall-time
424,400 s (incremental) vs 478,400 s (all-at-once) at equal evals (p.115).

### 3.2 ACO-MGA vs GA/NSGA-II/PSO (the integrated approach)

**Cassini (EVVEJS, 5 legs; the ESA-ACT "Cassini1" GTOP benchmark).** 100 runs,
Table 5.10 (p.209):

| Evals | Optimiser | Avg best (km/s) | % < 16 km/s | % feasible |
|---|---|---|---|---|
| 4000 | **ACO-MGA** | **16.24** | **44%** | **91%** |
| 4000 | GATBX (GA) | 16.349 | 14% | 25% |
| 4000 | NSGA-II | 20.426 | 5% | 26% |
| 4000 | PSO | 24.93 | 1% | 3% |
| 6000 | **ACO-MGA** | **15.434** | **80%** | **100%** |
| 6000 | GATBX | 16.526 | 17% | 28% |
| 6000 | NSGA-II | 20.122 | 7% | 37% |
| 6000 | PSO | 18.133 | 1% | 14% |

**BepiColombo (3-leg to Mercury, extended r_p bound), Table 5.5 (p.201):** ACO-MGA avg
5.67 km/s, 98% < 6 km/s, 100% feasible; GATBX 8.15 / 26% / 97%; NSGA-II 9.58 / 7% /
100%; PSO 11.32 / 5% / 97%. ACO-MGA finds the EVVMe sequence automatically.

The headline: ACO-MGA's win is overwhelmingly on **feasibility rate** (finding *any*
valid plan) on hard combinatorial sequences — exactly where a GA wastes evals on
infeasible sequences and ACO-MGA's taboo list does not.

---

## 4. Golden-eligible benchmark numbers (best-known objective values, with page refs)

These are *published best-known/reference objective values* for standard MGA
benchmarks. EXPECTED side traces to the thesis (and, where noted, to the ESA-ACT GTOP
database it reproduces). Usable as fixtures **if** we ever wire a comparable MGA
patched-conic searcher — NOT cycler-catalogue rows (see §5).

**Anchor C1 — Cassini1 (EVVEJS) reference solution, Table 5.12 (p.210).** ACO-MGA best
= **6.9686 km/s objective** (`y = v_∞ + β T`, β=1/1000 km/s/d; p.209). Per-component
ACO-MGA vs ESA-ACT reference: v_0 = 3.14 / 3.259 km/s; Δv_1 = 600/480 m/s; Δv_2 =
350/398 m/s; v_∞ = 4.21 / 4.246 km/s; legs T = 168/423/53/596/2290 d (ACO-MGA). The
ESA reference total-ΔV "Cassini1" GTOP value is 4.93 km/s minimum (the well-known GTOP
figure; the thesis quotes the per-component reference, not a single ΔV total here).

**Anchor C2 — ACT "Cassini1" position-formulation optimum, p.61.** "the total Δv,
including launch excess velocity and all the powered swing-bys, resulted to be **4.45
km/s**" (position formulation, no DSM); velocity-formulation analogue 4.58 km/s
(Fig.2.20, p.62). These bracket the known Cassini1 GTOP optimum and are a clean
patched-conic golden if we build the powered/position model.

**Anchor C3 — BepiColombo (EVVMe) reference, Table 5.7 (p.204).** ESOC reference:
v_0 = 3.79 km/s, v_∞ = 5.68 km/s, T = 438/674.1/630.9 d, Δv_1 = 7 m/s, Δv_3 = 11 m/s.
ACO-MGA 2D: v_0 = 3.63, v_∞ = 5.51 (lower because 2D ignores Mercury inclination);
3D-reoptimised matches ESOC within ~1%. (Launch t_0 = 4974.5 d MJD2000 = 15 Aug 2013.)

**Anchor C4 — EVM single-GA best, Table 3.7 (p.101).** Best total ΔV = **2.9818 km/s**
(t_0=4472.013 MJD2000, T_1=172.29 d, γ_1=2.9784, r_p,1=1, α_2=0.5094, T_2=697.61 d).
Bounds in Table 3.6 (p.100): t_0 ∈ [3650, 9128.75], T_1 ∈ [50,400], γ_1 ∈ [−π,π],
r_p,1 ∈ [1,5], α_2 ∈ [0,1], T_2 ∈ [50,700]. Self-contained reproducible fixture.

**Anchor C5 — EEM single-GA best, Table 3.11 (p.104) + Fig.3.21.** Best total ΔV =
**2.908 km/s** (t_0=5271 d, E swing-by at T_1≈495 d, DSM1 = 439 m/s, M arrival T_2=881
d, DSM2 = 506 m/s). Bounds Table 3.10 (p.104).

**Anchor C6 — physical constants used throughout, Table 3.5 (p.99).** μ_Sun =
1.3272e11 km³/s²; per-body μ, mean radius, period, semimajor axis for Me/V/E/M/J/S +
Ganymede/Callisto. AU = 149,597,870.7 km (Nomenclature p.xxix). A constants-table
sanity golden if cross-checking any reproduced number.

> Provenance flag on all of the above: these are **patched-conic, mean-element / SPICE
> de405 analytic ephemeris** results (Fig.2.3 compares the thesis's analytic ephemeris
> to JPL NAIF de405). NOT DE4xx-integrated, NOT cycler geometry. Any fixture must carry
> `source_ephemeris: patched-conic analytic (de405-fit) (Ceriotti 2010)`.

---

## 5. v4.2 backfill check — NO catalogue-eligible rows

- **No Earth-Mars cyclers anywhere in this thesis.** It is method-only MGA *mission*
  design (Cassini→Saturn, BepiColombo→Mercury, Laplace→Ganymede/Callisto). No
  periodic/repeating cycler trajectories, no S/L resonant-interval geometry, no Aldrin
  family. Nothing maps to a catalogue cycler row.
- **center:** all heliocentric (interplanetary) or jovicentric (Laplace moon tour);
  per-leg relative velocity is planet-local. No center ambiguity, but also nothing to
  backfill.
- **tof_days_bounds:** the benchmark tables give per-leg T bounds (Tables 3.6, 3.10,
  3.18, 4.2, 5-series) — relevant only if an MGA fixture is added, not a cycler row.
- **source_ephemeris:** patched-conic analytic ephemeris fitted to JPL de405 (Fig.2.3,
  Fig.2.4 "EphSS"); **mark any reproduced number `de405-fit analytic`, never DE440.**
- **Verdict: zero catalogue rows to add. The value here is algorithmic + a small set of
  MGA-benchmark goldens (§4) gated behind building a comparable searcher.**

---

## 6. Maps-to-our-X verdicts

| Ceriotti construct | Our code / concept | Verdict + scope |
|---|---|---|
| **Incremental pruning** (separable objective Eq.3.3, per-leg pruning criterion Eq.3.5-3.7, feasible-set boxing, back-pruning) | `search/scan.py` dense epoch×branch grid; `search/sequence.py` enumeration; the M4 enumeration + M8 VEM 3-body search | **MAPS — strong fit, this is the headline import.** Our broad search is exactly the "all-at-once dense grid" he beats. His per-leg pruning criterion that is *not the ΔV* but a geometric gate (Δa>threshold, `bend`-feasibility) is precisely our bend-feasibility prune already done post-hoc — promote it to a *level gate* and prune the upstream grid. **MED-HIGH** (needs a separable per-leg objective + box bookkeeping; the corrector already gives per-leg residuals). |
| **Pruning criterion ≠ objective** (Φ_2 = Δa-gain; Φ_3 = Δv<0.5) | `hunt_vem_ballistic.py` "0 bend-feasible" prune (post-hoc) | **MAPS — direct upgrade.** Move the bend-feasibility / Δa-gain test *before* the next leg's grid, not after. Cheapest high-value adoption. **LOW-MED.** |
| **Modified MBH stores the whole feasible *set*, not the best** (Alg.3.1) + restart | #145 MBH wrapper around `ballistic_correct` (just landed) | **MAPS — extends #145.** Our MBH keeps the incumbent; his keeps every `Φ_i`-feasible candidate to *characterise* the basin, enabling family/launch-window output. **LOW** add-on to the existing wrapper (collect-all + restart-after-N). |
| **Stop the local solve at the threshold, not the minimum** (MS, p.82) | inner solver in MBH / scan | **MAPS — efficiency trick.** For *feasible-set* identification (not single optimum) stop `least_squares` early once residual < tol; spreads samples. **LOW.** |
| **MS beats DE & MBH on these MGA problems** (Table 3.15) | #143 roadmap optimiser ranking (MBH ranked #2) | **INFORMS — caveat the roadmap.** A plain Multi-Start with a fast inner solve out-performed DE and MBH here; DE was *worst* (premature convergence at F=0.8/Cr=0.75). Our roadmap should (a) keep a Multi-Start baseline as the control, (b) treat DE as low-priority, (c) tune MBH's ρ rather than assume default — consistent with the Englander tuning paper in our survey. |
| **ACO-MGA integrated sequence+type search** (Alg.5.3-5.7) w/ feasible+taboo lists | `search/sequence.py` enumeration; Forge novelty loop | **PARTIAL — adopt the taboo list, not the ants.** Full ACO-MGA is HIGH scope and aimed at MGA *missions*, not cyclers. But the **taboo list of infeasible partial sequences** (forbid re-expanding a prefix that died at leg i) is a clean win for our sequence enumeration / VEM 3-body search to skip dead prefixes. **MED** for the taboo-list idea alone. |
| **Brent multi-zero phasing → tree of trajectories** (Alg.5.2) | our single-ellipse-per-leg corrector (M-ED blocker) | **INFORMS the multi-arc question.** His leg model carries *all* phasing zeros forward as a tree (resonant legs have multiple). This is structurally how a multi-arc topology gets represented — relevant to the M-ED single-ellipse blocker (memory `project_s1l1_realeph_closure_blocker`), though his is patched-conic not free-fall cycler. |
| Position formulation: polynomial complexity in #DSMs/#swing-bys; per-arc independence (§2.3) | `search/correct.py` multiple-shooting (per Russell mine) | **MAPS conceptually** — same per-arc independence as the Russell Ch.5 corrector; his §2.3 is the patched-conic-with-DSM analogue. Powered-swing-by "super-optimal" caveat (p.61) echoes the Russell Eq.5.5 powered-Δv note. |
| GTOP/Cassini/BepiColombo benchmarks (§4 anchors) | (no MGA test fixtures today) | **DOES NOT MAP YET** — golden-eligible only if a patched-conic MGA searcher is built. Park in acquisition list. |

---

## 7. Single most implementable finding (this thesis)

**Promote our post-hoc per-leg feasibility test into an incremental *pruning gate*
(Ch.3.2, Eq.3.5-3.7), with a non-ΔV geometric criterion.** We already compute exactly
the quantities Ceriotti prunes on — bend-feasibility and Δa-gain per swing-by
(`hunt_vem_ballistic.py` reports "0 bend-feasible"). Today that test runs *after* a leg
is built. Ceriotti's contribution is the *ordering*: search level 1, **box the feasible
survivors, and only grid level 2 inside those boxes** — with back-pruning. Because the
objective is a sum of non-negative per-leg terms (Bellman, p.67), this never discards
the global optimum yet collapses the search volume by ~90% (Table 3.16). For our M8 VEM
3-body search this directly attacks the "dense grids don't cure basin selection"
disease (#110/#120/#122) from the other side: not a better optimiser, but a *smaller,
pre-pruned* space the existing optimiser can actually cover. Scope MED-HIGH; the
LOW-cost first step is just moving the existing bend/Δa gate ahead of the next leg's
grid (no boxing yet).

**Runner-up:** extend the just-landed #145 MBH (Alg.3.1) to *store the whole feasible
set* + restart-after-N, turning it from a single-optimum finder into a launch-window /
family characteriser — and add a Multi-Start control, since MS beat both DE and MBH on
every MGA problem in this thesis (Table 3.15).

---

## 8. Honest "not read / not extractable" list

- **Ch.1 (Introduction, pp.1-27)** — background/literature on MGA missions and classic
  design; skipped per method-focus.
- **Ch.4 (Application to Real Case Studies, pp.123-169)** — full Laplace
  (Ganymede/Callisto resonant tour) and BepiColombo incremental-pruning results, Tables
  4.1-4.31. Read only via the List of Tables/Figures. These are moon-tour & Mercury
  *mission* results (jovicentric / Mercury-resonant), NOT cyclers — no catalogue
  relevance — but if an MGA fixture is ever wanted, the GCGC/GGGG Ganymede-tour ΔV
  numbers live here.
- **Appendix A (impl, MATLAB/C), Appendix C (affine transformation, pp.237-250),
  Appendix D (testing procedure for global-opt algorithms, pp.251-256), Appendix E
  (Tisserand plane, pp.257-258)** — not read. Appendix C's EEM/EEVVMe affine-transform
  result tables (C.2, C.5) and Appendix D's convergence-test definition (Alg.D.1-D.2,
  the formal success-rate methodology behind the 100-run % figures) are the rigour
  behind §3 if those numbers are ever promoted to goldens.
- **§2.1-2.2 (patched-conic arcs, velocity formulation, pp.29-48)** — read only §2.3.4
  discussion; the Lambert/swing-by/launch parameterisation equations (2.x) not
  transcribed (standard patched-conic).
- **Ch.5.1.1-5.1.3 (launch/swing-by/deep-space-leg 2D model eqs, pp.171-178)** — read
  from §5.1.4 on; the upstream 2D-model equations not transcribed.
- **GTOP "Cassini1" canonical minimum-ΔV** (≈4.93 km/s total ΔV, ESA-ACT): the thesis
  quotes per-component reference values (Anchor C1) and a 4.45 km/s position-formulation
  figure (Anchor C2), not a single headline ΔV total — confirm against the live GTOP
  database before quoting 4.93 as sourced-from-Ceriotti.
