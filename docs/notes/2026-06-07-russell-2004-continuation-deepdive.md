# Russell 2004 continuation method — implementable deep-dive (unblocks #143 rank-4)

Deep-dive mine 2026-06-07. Companion to
`2026-06-07-russell-2004-dissertation-method-mining.md` (which covered Ch.5 at
survey level). This note extracts the IMPLEMENTABLE detail that survey skipped:
the continuation parameterisation/schedule, the analytic Jacobian, the SNOPT/
elastic-mode/post-processing settings, the identity of the <1 m/s members, and
the failure modes. Everything below has a page number; all pages are the printed
(document) page numbers of:

**Source (cite exactly, no file path):** Russell, R. P., "Global Search and
Optimization of Free-Fall Cycler Trajectories [Earth-Mars cyclers]," Ph.D.
dissertation, The University of Texas at Austin, 2004 (advisor C. Ocampo). This
is the substantive open equivalent of the paywalled Russell & Ocampo 2006 JGCD
paper.

> Vision read of Ch.5 §5.4 (pp.146-171) and §5.5 Results (pp.172-182): the
> continuation method (§5.4.1, p.146), non-analytic rationale (§5.4.2, p.147),
> gradient/SNOPT (§5.4.3, pp.148-149), multiple-shooting transcription (§5.4.4,
> pp.149-151), powered-SOI Δv (§5.4.5 + Eqs.5.1-5.5, pp.152-155), analytic
> Jacobian (§5.4.6 + Eqs.5.6-5.13, pp.156-164), infeasibility/post-processing
> (§5.4.7 + Eqs.5.14-5.15, pp.165-167), the algorithm + homotopy schedule
> (§5.4.8 + Figs.5.3/5.4 + Table 5.4, pp.168-171), and the results
> (Tables 5.5/5.6, Figs.5.5-5.10, pp.172-182).

---

## 1. The continuation parameterisation — the thing #143 said we need (§5.4.1, §5.4.8)

**What is continued: the SOLAR-SYSTEM MODEL fidelity, NOT the trajectory genome.**
The cycler's per-leg unknown vector (the `5n` of §2) is held as the optimisation
variable; the homotopy parameter is the *planet model* that the position-match
constraints reference. Russell ramps two physical model parameters in sequence —
**eccentricity, then inclination** — and finally takes one step from the resulting
mean-element model to the true ephemeris. He explicitly does **not** continue an
epoch or a blend scalar, and he does **not** ramp `a, Ω, ω, ν` — those are pinned
to their J2000 mean values from the start (see "fix initial values", p.169).

**The homotopy schedule (Figure 5.3, p.169 — verbatim):**
```
0) Circular-coplanar model (mean values for Earth and Mars at J2000 for a, Ω, ω,
   and ν; but e = i = 0)
1) Use continuation to "walk" towards mean eccentricities for Earth and Mars with
   nstep equally spaced "steps"
2) Use continuation to "walk" towards mean inclinations for Earth and Mars with
   nstep equally spaced "steps"
3) Make one final "step" from the now mean-orbital-element model directly to the
   accurate ephemeris
```
So one full transition = `nstep` e-steps + `nstep` i-steps + 1 final ephemeris
step = `2·nstep + 1` SNOPT runs. Each step's converged unknowns seed the next
step (classic imbedding/homotopy, §5.4.1 p.146: "The solution to each sub-problem
becomes the initial guess for each successive sub-problem").

**Step-size control = the `nstep` ladder.** `nstep = 3^(steploop−1)` giving the
ladder `nstep ∈ {1, 3, 9, 27, 81, 243}` (Fig.5.4 step 7.1, p.170; powers of 3
"chosen somewhat arbitrarily... to efficiently sample a full range of values",
p.172). The full algorithm runs the SAME parent/launch-window through ALL SIX
`nstep` values and keeps the lowest-Δv result (Fig.5.4 step 8, p.170), because
"the `nstep` value that leads to the best solution varies depending on the parent
cycler and launch window" (p.171) — i.e. there is **no single good step size**;
robustness comes from trying several and keeping the best. For a single
production run he found **`nstep = 5` is "an excellent compromise between
efficiency and performance"** (p.171), and using the whole ladder buys only ~50%
solution improvement at large cost.

**How the circular-coplanar parent seeds the first ephemeris step (Fig.5.4 steps
1-7, p.170):** set Earth/Mars classic elements to J2000 mean values (Table 5.4),
zero their e and i (step 0 model), propagate the *simple* model forward
`LaunchWindow` synodic periods, INPUT the parent's beginning phase angle and
propagate until that phase angle is achieved, record that as the epoch time, then
INPUT seven cycles of the simple-model parent beginning at the epoch — that
seven-cycle simple-model trajectory is the `5n` initial guess for the homotopy.
LaunchWindow ∈ 1..21 (units of synodic periods after J2000; 21 chosen to observe
three complete cycles of the ~7-synodic-period true repeat — §5.3).

**Why the gap can't be jumped directly (§5.4.1, p.146; §5.4.8 "model changes too
rapidly", p.169):** continuation is mandatory because "the gap between the
circular-coplanar model and the accurate ephemeris model prohibits an immediate
jump." It is "very important to fix the initial values of `a, Ω, ω, ν`... to
their respective mean values at J2000, otherwise the model changes too rapidly"
(p.169) — i.e. only the two small perturbations (e≈0.0167/0.0934, Δi≈1.851°) are
ramped; the large elements are frozen so each homotopy step is a small
perturbation.

**Table 5.4 — Mean Elements at J2000 (p.169), the frozen+target endpoint values:**

| Element | Earth | Mars |
|---|---|---|
| a (AU)   | 1.00000101812E+00 | 1.52367934749E+00 |
| e        | 1.67086171540E-02 | 9.34006199474E-02 |
| i (deg)  | 0.0               | 1.84972647778E+00 |
| Ω (deg)  | 0.0               | 4.95655237028E+01 |
| ω (deg)  | 1.02937348083E+02 | 2.86494710278E+02 |
| ν (deg)  | −2.47089957222E+00| 1.93730406472E+01 |

(These are MEAN elements at J2000, not a DE-integrated ephemeris — see backfill
note §6. The "accurate ephemeris" of step 3 is whatever ephemeris file is plugged
in; the dissertation's force model is patched-conic two-body legs between true
ephemeris planet positions, p.146.)

---

## 2. Multiple-shooting transcription (§5.4.4, §5.4.5; Table 5.3 p.154)

**The `5n` unknown vector (per leg `i`, p.150 + Table 5.3):**
- `v∞+(i)` — 3 components of the *initial planet-referenced* spacecraft velocity (3n)
- `t0(i)` — leg initial time (n)
- `tf(i)` — leg final time (n)
- (the two planet IDs per leg are FIXED, not optimised → not counted)

Each leg is made **completely independent** (no leg depends on the previous leg's
output) to avoid the "highly erratic and non-linear perturbations at the end of
the long trajectory" that single-shooting a multi-decade patched trajectory
produces (p.150). Spacecraft initial position is pinned to the planet's position
at `t0(i)` via the ephemeris.

**Constraint set (Table 5.3, p.154 — verbatim):**

| Parameter | Number | Constraint | Number | Type |
|---|---|---|---|---|
| `v∞+(i)` | 3n | `r_f(i) = r_f,plan(i)` | 3n | nonlinear equality (position match) |
| `t0(i)`  | n  | `t_f(i−1) = t_0(i)`    | n−1 | linear equality (time continuity) |
| `tf(i)`  | n  | `Δv(i) = 0`            | n−1 | nonlinear equality (flyby Δv) |
| **total**| **5n** | | **5n−2** | |

So it is a **non-linear root-find: `5n` unknowns, `3n + 2(n−1)` equality
constraints, no explicit objective** (p.153). The time-continuity constraints are
linear and "for all practical purposes non-existent because SNOPT solves [them]
completely during every [major] iteration" (p.152), which is what keeps the
Jacobian block-diagonal-sparse.

**Powered-SOI flyby → single conditional Δv (§5.4.5, Eqs.5.1-5.5, pp.152-155).**
The ballistic-flyby pair of constraints `|v∞−(i−1)| = |v∞+(i)|` (Eq.5.1) and
`ω_req ≤ ω_avail` (Eq.5.2) are REPLACED by one conditional equality `Δv(i) = 0`
whose value is (Eq.5.5, p.154):
```
Δv0(i) = | v∞+(i) − v∞−(i−1) |                                          if ω_req ≤ ω_avail
       = sqrt( v∞+(i)² + v∞−(i−1)² − 2 v∞+(i) v∞−(i−1) cos(ω_req − ω_avail) )  if ω_req > ω_avail
```
with `ω_req = acos( v∞−(i−1)·v∞+(i) / (v∞−(i−1) v∞+(i)) )` (Eq.5.3) and
`ω_avail = 2 asin( μ_plan / (μ_plan + r_p,min v∞,small²) )` (Eq.5.4). The maneuver
is placed just outside a zero-radius SOI, *before* the flyby if v∞−(i−1) < v∞+(i)
else after; powered burns inside the hyperbola are excluded as too risky (p.155).
**Minimum flyby radii (p.165): `r_p,min,Earth = 6578.0 km`, `r_p,min,Mars =
3598.5 km`** (Case B uses a 200 km periapse altitude, p.156). This is exactly the
continuous-Δv replacement for our binary `bend_feasible` gate flagged last pass.

**Derivatives = ANALYTIC, computed alongside the constraints, no extra calls
(§5.4.6, pp.156-164).** Russell explicitly rejects finite differences (a full
central-difference Jacobian costs `2k` extra constraint calls for `k` unknowns,
p.157) and bypasses SNOPT's automatic-differentiation hook because "the general
structure of the constraints and unknowns is fixed and it is possible to derive
the analytic Jacobian... the run-time benefit of using analytic derivatives is
crucial" over thousands of cases (p.158). The Jacobian is `(5n−2)×(5n)` =
`25n²−10n` entries, almost all zero (block-bidiagonal, structure shown p.160).
Key building blocks:
- The hard partials are the final-state-vs-initial-state and -initial-time
  blocks, obtained from the **two-body state transition matrix Φ(t,t0)** (Eq.5.8,
  p.161) — and for Keplerian motion Φ is **analytic** in classic elements with
  the only non-identity entry `Φ_α6,1 = −3(t−t0)√(μ/a⁵)/2` (p.162). So no numeric
  STM integration is needed.
- To avoid the e→0 / i→0 singularities in the classic-element partials (which
  *every* cycler hits — there is always an ecliptic Earth-Earth generic leg, and
  most half-rev returns have near-zero e, p.162), he switches to a **non-singular
  element set** `β = [a, e·sin(Ω+ω), e·cos(Ω+ω), sin(i/2)·sin(Ω),
  sin(i/2)·cos(Ω), Ω+ω+M]ᵀ` (p.162) and maps via `∂x/∂x0 = (∂x/∂β)(∂β/∂x0)`
  (Eq.5.10, p.163). **This directly informs our family-parameter design: the
  e=i=0 endpoint is a coordinate singularity, so the partials (and any of our own
  STM/derivative work) must use equinoctial-style non-singular elements near the
  circular-coplanar parent** — the #137 6.44Gg3 conditioning warning is the same
  failure surface.
- Partial of Δv constraint vs the v∞ unknowns is given in closed form
  (Eq.5.4.6.3, p.164), with separate branches for ω_req ≤ / > ω_avail matching
  Eq.5.5.

**SNOPT settings (§5.4.3 p.149, §5.4.7 pp.165-167, §5.4.8 p.171):**
- Optimiser: **SNOPT** (Sparse Nonlinear OPTimizer, Stanford Business Software);
  chosen over **VF13 (Harwell)** which "often fails when the constraints become
  infeasible," and over GA/Monte-Carlo (too slow over all cases/dates).
- **Elastic mode**: when constraints can't be met, SNOPT automatically minimises
  "a weighted sum of the absolute values of the constraint violations" (p.149) so
  the continuation can keep "walking" even when favourable solutions vanish on the
  path. Because there is no explicit objective, "a specified weight is irrelevant"
  inside pure elastic mode (p.153) — the composite objective becomes
  `J = Σ ( Δv0(i) + |Δr_x,f(i)| + |Δr_y,f(i)| + |Δr_z,f(i)| )` (Eq.5.14, p.166).
- **Jacobian supplied analytically every major iteration; Hessian intrinsically
  approximated by SNOPT** (p.158).
- **Default SNOPT input tuning is left unchanged** for the production runs (p.171)
  — he tried tuning loops but they "only provide marginal solution improvement."
- **Convergence/tolerances:** the dissertation does not print explicit SNOPT
  feasibility/optimality tolerances. Convergence quality is reported downstream as
  residual Δv (m/s over 7 cycles) — the operative "tolerance" for our purposes is
  the <1 m/s / <10 m/s / <300 m/s tiers in §3. (Honest gap: exact SNOPT
  `Major feasibility tolerance` etc. are not in the text.)

**Post-processing to close residual position (§5.4.7, Eq.5.15, pp.166-167):** the
elastic-mode solution is NOT guaranteed to intercept the planet, so a hybrid is
used: optimise with a weighted objective `J = Σ [ Δv0(i) + W(|Δr_x,f| + |Δr_y,f| +
|Δr_z,f|) ]` (Eq.5.15) with **`W = 1` exactly** ("the serendipitous best value of
1 is indicative of the benefits of using canonical units," p.167 — μ_sun = 1, so
position and velocity magnitudes are comparable). Then **post-process**: for each
leg with leftover position residual, add an intermediate Δv and **globally
optimise its epoch with a 1-D grid search over 100 equally spaced intervals**
(avoiding `nπ` transfer angles), and add those Δv's to the total cycler cost
(Fig.5.4 step 9, p.170; p.167). Final trajectories are scored on **total
accumulated Δv after post-processing**.

---

## 3. The <1 m/s over seven cycles result (§5.5, Tables 5.5/5.6, pp.172-182)

**The headline claim (p.137 abstract; p.176):** of 203 parent cyclers ×
21 launch windows = 4263 cases, **nine parent cyclers have at least one launch
date with total maneuver requirement < 1 m/s over seven full cycles** (also
thirty-nine < 10 m/s, seventy-four < 300 m/s). "Over seven full cycles" = the
seven-synodic-period propagation that is the true-model approximate repeat time
(§5.3); each cycle is the parent's simple-model repeat (1/2/3-synodic).

**What was measured:** *total accumulated Δv* (sum of all powered-flyby Δv from
Eq.5.5 PLUS post-processed intermediate position-continuity Δv) over the full
seven-cycle propagation, for the **minimum-Δv launch window** of that parent.
Reported in **m/s** in Table 5.5 (p.178, "total Δv" column).

**The golden-eligible targets (Table 5.5, p.178 — the members reading 0 m/s at
their best launch window):**
- **6.399G1 (#1) — the Aldrin cycler.** total Δv = **0 m/s**, launch **Aug-2003**
  (specifically Aug 6, 2003, Fig.5.7 p.181), avg E-M transit 143 d, avg v∞
  E−/E+/M−/M+ = 6.02/6.63/9.33/9.33 km/s. One-synodic-period. Text (p.176):
  "a seven cycle propagation with a launch date of August 6, 2003... is found
  that is completely ballistic." Note the previously published impulsive Aldrin
  (Nov 1996 launch) needed 1.73 km/s/7cyc (Byrnes/Longuski, his Ref.13) — so this
  ballistic version is the headline find. **This is the cleanest single golden
  target: a named, dated, 0-m/s, one-synodic-period cycler.**
- **4.991gG2 (#83) — the S1L1 cycler.** total Δv = **0 m/s**, launch Jun-2025,
  avg transit 165 d, v∞ 5.37/5.37/5.48/5.48 km/s. Two-synodic-period;
  "essentially ballistic for all launch dates" (p.176), consistent with his
  Ref.15. **Directly relevant to our S1L1 closure blocker** (see §5).
- Other Table 5.5 rows reading total Δv = 0 m/s at their min-launch window:
  **3.768Gh-3 (#54)** Dec-2024; **3.768Gh+3 (#55)** Jul-2018; **5.658Gfh-3
  (#126)** Feb-2031; **5.658Gfh-3 (#147)** Oct-2037 (reads 6 m/s — see note);
  and several 3-synodic 5.658-prefix members hovering at Δv≈0 (Fig.5.6, p.175).
  The exact nine-member list is not enumerated as a single table; the nine are
  those with a 0/sub-1-m/s entry — #1 and #83 are the unambiguous, named ones.

**Lowest-energy (not lowest-Δv) member, useful as a transit-time anchor:**
**cycler #23 (3.406gGh+rh-3)** has the lowest excess speed of any solution — avg
v∞ 4.12 km/s (Earth) / 4.79 km/s (Mars), avg transit 187 d (p.177).

**Force model / ephemeris for the result:** patched-conic two-body legs between
**true ephemeris planet positions** (p.146); the continuation endpoint is "the
accurate ephemeris" (Fig.5.3 step 3). The dissertation does not name DE405
explicitly — the homotopy *target* is whatever ephemeris file is supplied; the
intermediate models use Table 5.4 J2000 mean elements. **For a V3 rung golden,
target the Aldrin #1 Aug-6-2003 case and reproduce 0 m/s/7cyc against a true
ephemeris (DE40x), accepting that "0" means below the post-processing/SNOPT
residual floor, not literally zero.**

**Avg-Δv-over-all-21-windows ranking (Table 5.6, p.180):** best is **5.658Gfh+f3
(#162) at 24 m/s avg**; twenty parents average < 100 m/s/7cyc. This is the
"launch-flexibility" metric (a parent that is cheap across many windows, not just
one). #1 Aldrin averages 3297 m/s (cheap only at its repeat windows), #83 S1L1
averages 714 m/s — i.e. the Aldrin/S1L1 are *spiky* (great at specific windows,
expensive otherwise), whereas the 5.658-family three-synodic cyclers are *flat*
and cheap everywhere. Performance tracks the circular-coplanar **turn ratio**
(p.180).

---

## 4. Failure modes Russell documents (informs our family-parameter risk)

1. **Continuation can dead-end on intermediate infeasibility (§5.4.1, p.146):**
   "solutions to each sub-problem must exist in order to continue along the
   path... the problem may become infeasible at one or more of the intermediate
   steps." A solution can exist at the true ephemeris yet be unreachable because
   the homotopy path passes through an infeasible region. **Workaround: SNOPT
   elastic mode** — keep minimising violation and "walk" through the infeasible
   patch rather than halt (p.149). (Our `scipy.least_squares` has no elastic mode
   — the robustness gap from last pass; this is the single biggest missing piece.)

2. **Integer-structured (analytic) transfers can't morph under a gradient
   optimiser (§5.4.2, pp.147-148):** full-rev/half-rev/generic analytic transfers
   each need ≥2 integer inputs (revs, fast/slow) that "fix the structure of each
   leg." As dates become free, a short-period 1-rev Mars→Earth transfer "may
   cease to exist," and changing the integer "is a difficult if not insurmountable
   obstacle when dealing with a gradient-based optimizer." **Workaround: the
   non-analytic velocity-vector search** — numerically search for the v∞ vector
   that intercepts the target, "eliminat[ing] the integer programming problem,
   and the solution structure is free to morph." He notes analytic fixed-structure
   "works very well for cyclers with favorable turning angle requirements, such as
   the S1L1 cycler" but is too rigid in general. **This validates our radial-
   crossing (a,e)-genome pivot** and is the canonical written justification to
   cite in `free_return.py`.

3. **VF13 fails on infeasible constraints (§5.4.3, p.149):** the first SQP package
   tried "often fails when the constraints become infeasible." Workaround: switch
   to SNOPT (elastic mode). A pure-gradient method also gives only *local*
   solutions; he accepts this, arguing the ballistic circular-coplanar parent is a
   good basin seed (p.149) — i.e. **family selection must start from a genuinely
   ballistic parent, or the local solve lands off-basin** (our recurring S1L1
   off-basin lesson).

4. **Coordinate singularities at e→0, i→0 in the partials (§5.4.6.1, p.162):**
   classic-element partials blow up at zero e/i, which *every* cycler hits (the
   ecliptic Earth-Earth generic leg; near-zero-e half-rev returns). Workaround:
   the non-singular `β` element set (§2 above). **Direct family-parameter design
   risk: any derivative/STM/continuation we build that touches the circular-
   coplanar endpoint must use equinoctial-style non-singular elements** — this is
   the analytic root of the #137 6.44Gg3 conditioning warning.

5. **Model changes too rapidly if the big elements aren't frozen (§5.4.8,
   p.169):** if `a, Ω, ω, ν` are not pinned to J2000 mean values, "the model
   changes too rapidly" between homotopy steps and the solve fails. Workaround:
   freeze them; ramp only e and i. **Design implication: our continuation driver
   must ramp ONLY the small perturbations (e, i), never the large elements or the
   epoch.**

6. **No universal step size (§5.4.8, p.171):** "the `nstep` value that leads to
   the best solution varies depending on the parent cyler and launch window."
   Workaround: run the `{1,3,9,27,81,243}` ladder and keep the best (or use
   `nstep=5` as a single-shot compromise). **Design implication: the driver needs
   an outer "try several step counts, keep best-Δv" loop, not a fixed schedule.**

---

## 5. Maps-to-our-X + build sketch (unblocks #143 rank-4)

| Russell construct | Our asset | Verdict / contribution |
|---|---|---|
| `5n` multiple-shooting unknowns (v∞ vec + t0 + tf per leg), constraints position-match + time-continuity + Δv=0 | `search/correct.py` free-return corrector | **MAPS** — our corrector is his §5.4.4 method, simplified (we carry magnitudes+ToFs; he carries full 3-vec v∞ + position match). Reuse as the inner solve at each homotopy step. |
| Powered-SOI conditional Δv (Eq.5.5) + r_p,min 6578/3598.5 km | `core/flyby.max_bend` / `bend_feasible` | **PARTIAL → upgrade** — we already compute ω_avail (Eq.5.4); wire Eq.5.5 to emit a continuous Δv instead of a binary gate. ~10 lines. |
| Continuation schedule e-ramp → i-ramp → ephemeris, `nstep=3^(steploop-1)`, freeze a/Ω/ω/ν | `search/seed_ladder.py` (partial) | **THE NEW PIECE** — the *continuation driver* is what #143 says we lack. seed_ladder gets us part way; the e/i homotopy + nstep ladder + keep-best is genuinely new. |
| `optimise_cell_ephemeris` (the cell-level ephemeris solve we want) | — | **= his Fig.5.4 inner block** (steps 7.1-7.4 + post-process step 9). Build `optimise_cell_ephemeris` as: seed from circular-coplanar cell → for nstep in ladder: walk e, walk i, step to ephemeris (each = one corrector run) → keep best → post-process residual Δv. |
| n-body harness (REBOUND/planets-on-rails) | `nbody/shooter.py` | **= his "accurate ephemeris" step-3 force model**, upgraded. Russell's force model is patched-conic two-body legs between true-ephemeris planet positions; our n-body harness IS the higher-fidelity version of his step 3, so it slots in as the homotopy endpoint. |
| Analytic Jacobian via two-body STM + non-singular β elements (Eqs.5.8-5.13) | — | **MED-scope new** — only needed if we want his run-time (thousands of cases). For first cut, finite-diff is acceptable (he says it's "easiest, most common," just slower). Non-singular elements are the must-have, not the analytic partials. |
| SNOPT elastic mode | `scipy.least_squares` | **GAP** — no elastic mode. Mitigations: (a) the homotopy itself reduces infeasibility by keeping steps small; (b) consider `scipy.optimize.minimize(method='SLSQP')` with a violation-penalty objective (Eq.5.14/5.15) as a crude elastic surrogate; (c) the keep-best-over-nstep loop recovers some robustness. |

**Concrete build sketch (the continuation driver, `optimise_cell_ephemeris`):**
```
optimise_cell_ephemeris(parent_cell, launch_window, ephemeris):
  u = seed_5n_from_circular_coplanar(parent_cell, launch_window)   # Fig.5.4 steps 1-7
  best = None
  for nstep in [1,3,9,27,81,243]:                                  # or just [5]
    u_h = u.copy()
    for k in range(1, nstep+1):                                    # step 1: ramp e
        model = blend(e=k/nstep * e_target_J2000, i=0, big=frozen_J2000)
        u_h = correct(u_h, model)        # = correct.py, NON-SINGULAR elements
    for k in range(1, nstep+1):                                    # step 2: ramp i
        model = blend(e=e_target, i=k/nstep * i_target_J2000, big=frozen_J2000)
        u_h = correct(u_h, model)
    u_h = correct(u_h, ephemeris)                                  # step 3: true ephemeris
    dv = total_dv(u_h) + postprocess_position_residual(u_h)        # Eq.5.5 + grid-search Δv
    best = min(best, (dv, u_h))
  return best
```
- `correct(u, model)` = our corrector with the planet model swapped (the only
  change per step is what the position-match constraint references).
- Freeze `a, Ω, ω, ν` at Table 5.4 values throughout; ramp only `e`, `i`.
- Use non-singular `β` elements inside the corrector near the e=i=0 seed
  (failure mode 4) — this is the part most likely to bite if skipped.
- `postprocess_position_residual` = add per-leg intermediate Δv, 1-D grid-search
  its epoch over 100 intervals (Eq.5.15 hybrid, W=1).

**Scope:**
- Wire Eq.5.5 continuous flyby Δv into `correct.py`: **LOW** (~10 lines, geometry
  already present).
- Continuation driver (e/i homotopy + nstep ladder + keep-best + post-process):
  **MED** — it's orchestration over the existing corrector; the new logic is the
  model-blend and the keep-best loop.
- Non-singular element corrector internals: **MED** (needed for the e=i=0 seed to
  not blow up; equinoctial conversion is well-trodden).
- Analytic two-body-STM Jacobian: **HIGH**, and **deferrable** — finite-diff works
  for a first golden, only run-time suffers.

**The single thing that unblocks #143 rank-4:** the continuation driver above. It
is genuinely new (not in our tree), but it is *orchestration* over assets we
already have (`correct.py`, `nbody/shooter.py`, the cell seeds) — MED scope, not
a from-scratch solver. First milestone: reproduce **Aldrin #1, launch Aug 6 2003,
0 m/s / 7 cycles** as a V3 rung golden against a true ephemeris.

---

## 6. v4.2 backfill checks

- **center:** heliocentric (μ_sun = 1 canonical units, p.167); planet-referenced
  v∞ is per-leg local. No center ambiguity.
- **tof_days_bounds:** legs are 15-49 per cycler over 7 cycles; per-leg transit
  times appear in Table 5.5 (avg E-M transit, days) and fully in Appendix C
  (pp.201-245, not transcribed). Aldrin #1 avg transit 143 d; #23 (low-energy)
  187 d; high-energy #188 95 d. Source any added Russell row's tof bounds from
  Appendix C.
- **source_ephemeris:** the continuation *intermediate* models use **mean orbital
  elements at J2000 (Table 5.4)** — NOT a DE ephemeris. The final step targets
  "the accurate ephemeris" (file not named; not stated to be DE405). Force model =
  patched-conic two-body legs between true ephemeris planet positions. Any
  catalogue row tracing to a Russell-2004 *intermediate* solution must carry
  `source_ephemeris: mean-element J2000 (Russell 2004)`; a final-step solution is
  `accurate-ephemeris (patched-conic, Russell 2004)` — do not assert "DE405".

---

## 7. Honest "not extractable / not in text" list

- **Exact SNOPT numeric tolerances** (Major feasibility/optimality, iteration
  limits) are not printed — defaults left unchanged (p.171). Convergence is
  reported only as residual Δv tiers (<1 / <10 / <300 m/s).
- **The full nine-member <1 m/s list** is not given as one table; it is the set of
  Table 5.5 rows with a 0/sub-1-m/s min-launch entry. #1 (Aldrin) and #83 (S1L1)
  are the named, unambiguous, 0-m/s ones; the rest are 3-synodic 5.658/3.768-prefix
  members read off Table 5.5 / Figs.5.5-5.6.
- **The full β-element partial expressions** (∂x/∂α, ∂x/∂β) are not re-written —
  Russell defers to his Refs.62/64 ("require several pages... not re-written
  here," p.163), noting Ref.64 Appendix-2 FORTRAN is correct but Appendix-1
  equations have typos. If we implement the analytic Jacobian we need those refs,
  not this dissertation.
- **Appendix C** (pp.201-245): the 77 full reproducible trajectories (per-leg
  states, the source for catalogue member backfill) — not transcribed this pass.
  The other 4186 solutions are "archived electronically; contact the author."
- The `nstep=5` single-shot runtime ~1 day on a 3.2 GHz CPU; full ladder ~4 days
  on 15×2 GHz machines (p.171) — sizing reference for our own runs.
