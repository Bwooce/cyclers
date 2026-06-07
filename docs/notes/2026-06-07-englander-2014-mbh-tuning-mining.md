# Englander & Englander 2014 — "Tuning Monotonic Basin Hopping" (algorithm mine)

Mined 2026-06-07 (sources the `src/cyclerfinder/search/mbh.py` perturbation spec,
which was an UNSOURCED documented default pending this paper — see
`docs/notes/2026-06-07-mbh-wrapper.md`). Full read of all 33 pages.

**Source (cite exactly, no file path):**
Englander, J. A. and Englander, A. C., "Tuning Monotonic Basin Hopping: Improving
the Efficiency of Stochastic Search as Applied to Low-Thrust Trajectory
Optimization," 24th International Symposium on Space Flight Dynamics (ISSFD24),
2014, paper S7-3. (Jacob A. Englander, NASA GSFC Navigation & Mission Design
Branch; Arnold C. Englander, Deployment Technologies Inc.)

> Clean digital typeset; equations, algorithms, tables (1-4) and figures (1-17)
> read unambiguously.

---

## 1. The one-line headline

For MBH, draw the per-step random-walk perturbation from a **long-tailed
distribution (Cauchy or bi-polar Pareto) centred at zero**, NOT from the classical
uniform or from Gaussian. Long tails make MBH both more **efficient** (better
solution in less time) and more **robust** (efficiency that does not depend on
tuning the step-size or on the problem's boundaries/constraints). Of the two
long-tailed options the authors **recommend bi-polar Pareto** as the single best
default because it is the most robust to the choice of its excursion parameter.

---

## 2. The perturbation distribution spec (the thing mbh.py needed sourced)

### 2.1 What is perturbed and how (the classical baseline, §3.3, pp.6-7)

MBH is a two-step loop (Algorithm 1, p.9): (1) a *global reset* — random point in
the decision space, run the NLP solver, repeat until feasible; (2) a *local hop* —
"a small random perturbation vector is added to `x*`, producing a new `x'`, and
then the NLP solver is run" (p.6). Accept iff the result is feasible AND superior:

> "If the resulting solution is both feasible and superior to `x*`, then it is
> adopted as the new `x*` and the hopping process begins again. Otherwise, MBH
> attempts a new hop from the current `x*` and an 'impatience' counter
> `N_not improve` is incremented." (p.6)

Classical random step (the historical default we were standing in for):

> "In the classical version of MBH, the random step is drawn from a uniform
> probability distribution in `[−σ, σ]`." (p.7)

> "MBH has three parameters - the stopping criterion, the parameter
> `N_not improve`, and the type of random step used to generate the perturbed
> points `x'`." (p.7)

**The contribution (p.7):**

> "The contribution of this paper is the investigation, both experimentally and
> theoretically, of using RVs from distributions other than uniform for MBH,
> specifically Cauchy and Power Law distributions chosen because of their very
> long tails as originally suggested in [3]." (p.7)

### 2.2 The four distributions and their excursion parameters (Table 3, p.14)

The "excursion parameter" is the per-distribution knob (= step-size in classical
MBH). Centre/location is **zero** in all cases — i.e. the perturbation is
**bi-directional / symmetric about zero** ("the RW is equally likely to move in
any direction", p.13). All RVs are **i.i.d.** and applied per-step.

| Distribution | Excursion parameter = | 16-value sweep used (Table 3, p.14) |
|---|---|---|
| Uniform(−stepsize, stepsize) | abs. bound (±stepsize) | 0.01, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14, 0.16, 0.18, 0.20, 0.22, 0.24, 0.26, 0.28, 0.30 |
| Gaussian(mean=0, σ) | std deviation σ | same 0.01 … 0.30 sweep |
| Cauchy(location=0, scale) | scale factor | 0.000125, 0.00025, 0.0005, 0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.008, 0.010, 0.012, 0.014, 0.016, 0.018, 0.020 |
| Bi-polar Pareto(location=0, alpha) | Power-Law exponent α | 1.0025, 1.0050, 1.0075, 1.010, 1.015, 1.020, 1.030, 1.040, 1.050, 1.060, 1.070, 1.080, 1.090, 1.10, 1.11, 1.12 |

Note the **scale ranges differ by distribution by orders of magnitude**: the
useful Cauchy *scale* is ~0.0001-0.02 (≈10-30× SMALLER than the uniform/Gaussian
step), and Pareto is parameterised by an *exponent* α≈1.0-1.12 (NOT a step-size at
all). The long tail does the large excursions; the scale is kept tiny so the
*typical* step stays small. **This is the single most important quantitative point
for us:** a Cauchy `scale` is not interchangeable with a Gaussian `σ`.

> "In the case of a Cauchy distribution, the excursion parameter is what
> probability theory specialists refer to as the scale factor, and in the case of
> a bi-polar Pareto distribution, the excursion parameter is what probability
> theory specialists call the Power Law exponent alpha." (p.13)

> "For each of the four probability distributions, the 'average' value of the RVs
> was set to zero, meaning that range of possible values is symmetric around
> zero." (p.13)

> "the 'average' is undefined because the infinitely log[long] tails of these
> distributions prevents their being integrated … However probability theory
> allows for the range of possible values drawn from Cauchy and bi-polar Pareto
> distributions to be centered on zero by using setting parameter called location
> to zero." (p.13)

### 2.3 The exact RV generators (Table 4, p.14)

`r = uniform(0.0, 1.0)`; `s` = a fair coin flip (equal prob −1.0 / +1.0):

| Distribution | RV generator |
|---|---|
| Uniform | `2(r − 0.5)` |
| Gaussian | `s/(σ√(2π)) · exp(−r²/(2σ²))` |
| Cauchy | `tan(π(r − 0.5))` *(generated from the CDF; scale multiplies this)* |
| Pareto | `(s/ε)·((ε/(ε+r))^(−α))` *("bi-polar" = the `s` coin flip gives both signs)* |

> "The Uniform, Gaussian, and Pareto RVs were generated from their probability
> density functions (PDFs) while the Cauchy RVs were generated from the cumulative
> density function (CDF)." (p.14)

The "bi-polar" in **bi-polar Pareto** is literally the `s` fair-coin sign flip
making the otherwise one-sided Pareto symmetric about zero. (`ε` is a small
positive offset in the generator; not separately tabulated.)

### 2.4 Per-variable vs global scaling

The paper applies **one excursion parameter per run across the whole decision
vector** (it sweeps a single scalar parameter, Table 3). It does NOT prescribe
per-variable scales — but note the EMTG decision variables are pre-normalised by
their own bounds via Table 1 (§3.4, pp.7-8: bounds chosen heuristically so each
variable lives on a comparable normalised range; flight-time bounds are functions
of the bodies' periods `τ_i`). So in EMTG a single global excursion parameter is
effectively per-variable *because the variables are already scaled to their
ranges*. **This is a caveat for us (see §6): our genes are NOT range-normalised, so
a single global scale is wrong for us — confirming why mbh.py already needed the
per-gene relative/absolute split.**

### 2.5 QUANTITATIVE comparison data backing the recommendation

Test problem (§4.1, Table 2, p.10): a single 503-variable / 146-constraint
low-thrust VSI radioisotope-electric mission to Uranus, EESU (Earth-Earth-
Saturn-Uranus) flyby sequence, 40 time-steps/phase, objective = maximize final
mass, best known cost = −0.3007 (3007 kg delivered). 64 cases = 4 distributions ×
16 excursion parameters; each case = 10000 MBH steps, binned into 1000 bins of 10
("a path"), recording best fitness so far per bin (§4.2, pp.13-14).

There are **no success-rate tables** in this paper (unlike Ceriotti's % tables);
the evidence is the best-fitness-vs-step convergence curves, Figures 9-12:

- **Fig. 9 (best paths, p.17):** "The best solution was found by MBH driven by the
  Cauchy distribution, but a very nearly equally as good solution was found in
  much less search time by MBH driven by the bi-polar Pareto distribution. The
  bi-polar Pareto distribution is therefore more efficient in the sense of finding
  a better solution. Both the Cauchy and the bi-polar Pareto MBH found a better
  solution and were more efficient than the Gaussian and Uniform MBH." (p.15)
  — i.e. **best objective: Cauchy ≳ Pareto > Uniform ≈ Gaussian; time-to-good:
  Pareto fastest.**
- **Fig. 10 (worst paths, p.18):** "the bi-polar Pareto MBH finds a better
  solution than the other variants in much less time." (p.17) — i.e. Pareto wins
  even at its *worst* excursion parameter (robustness).
- **Fig. 11 (averaged across 16 excursion params, p.18):** "MBH driven by RVs
  drawn from bi-polar Pareto distribution finds a better solution in less time
  than the other variants." (p.18)
- **Fig. 12 (std-dev across the 16 params, p.19):** "The bi-polar Pareto
  distribution is the most robust in the sense that it's time-series of standard
  deviations across its excursion parameters decays the soonest and steepest."
  (p.19) — i.e. Pareto's performance varies *least* with its knob.

**Recommendation, verbatim (p.31):**

> "We find that of the two long-tailed distributions, bi-polar Pareto (in contrast
> to Cauchy), is the more robust to variations in excursion parameters. Therefore,
> we are inclined to recommend the bi-polar Pareto as the distribution that best
> improves the performance of MBH over the performance achieved by the classical
> use of the uniform distribution." (p.31)

> "The results found in this paper are directly applicable to improving the
> performance of EMTG on low-thrust trajectory optimization problems and have been
> adopted as the default settings for EMTG's MBH+SNOPT optimizer." (p.31)

---

## 3. The acceptance / restart policy (§3.3, Algorithm 1, p.9)

- **Monotonic acceptance:** accept the hop iff feasible AND `f(x*) < f(x_current)`
  — strict improvement only; never an uphill move (Algorithm 1, p.9; p.6). This is
  *exactly* what mbh.py does (`improves = feasible and cand_obj < best_obj`).
- **Two-phase search-then-improve:** (a) global reset loop draws uniform random
  points and runs the NLP until a feasible point is found ("MBH will explore the
  solution space via the global reset operator until a feasible solution is
  found", p.6); (b) then local-hop exploit loop. "The exploitation process
  continues until the algorithm fails to improve the current point
  `Max_not improve` times, and then MBH will switch back to exploring the entire
  solution space." (p.6)
- **Stall / restart criterion:** the impatience counter `N_not improve` increments
  on each non-improving hop; when it reaches `Max_not improve` MBH **resets to a
  fresh global random point** (the inner `while N_not improve < Max_not improve`
  loop in Algorithm 1). mbh.py's `stop_after_stall` *stops* on this condition
  rather than *resetting* — see §6 caveat.
- **Archive:** "Each feasible solution is stored in an archive" (p.6); the run
  returns the best in the archive. mbh.py keeps only the incumbent (Ceriotti's
  modified-MBH "store the whole feasible set" is the richer variant — see
  `2026-06-07-ceriotti-2010-mga-global-opt-mining.md` §6).
- **Stopping:** "MBH is run until either a specified number of iterations (trial
  points attempted) or a maximum CPU time is reached" (p.7). The benchmark used a
  96-hour wall-clock budget (Table 2, p.10).
- **SNOPT freeze guard (impl detail):** "EMTG contains a timer that ends any SNOPT
  run that continues for longer than some threshold time typically set to a few
  minutes." (p.7) — relevant only if an inner solver can hang; ours can't here.

---

## 4. Problem-dependence caveats (honest scope)

**Specific to the experiment / low-thrust Sims-Flanagan NLP:**
- The ONLY benchmark is a single 503-var EESU low-thrust mission (Table 2). The
  authors themselves flag this: "the numerical results presented in this work are
  based on a limited set of test cases … it is possible that the best tunings for
  the various distributions are dependent on the problem chosen." (p.29) → **the
  numeric scale/α VALUES (Table 3 winners) are NOT portable; only the qualitative
  long-tail ordering is.**
- The decision variables are EMTG's: normalised by bounds (Table 1), 500+ dims,
  tightly bounded with plentiful internal constraints (flyby-feasibility Eqs.4-6,
  match-point continuity Eq.3). The robustness argument is *built on* these
  boundary/constraint properties (§5.3).

**Plausibly general (the authors' theory claim):**
- §5 (pp.19-29) is a diffusion-theory argument that the long-tail advantage is
  *not* problem-specific: MBH is a random walk; long-tailed RVs give super-diffusive
  (Lévy-flight) walks (MSD `⟨r²(t)⟩ ∝ t^α, α>1`, p.20) that "jump over" internal
  constraints, whereas uniform/Gaussian walks get anti-correlated ("pinned")
  against boundaries and turn sub-diffusive (§5.3, pp.23-28; Figs.13-16).
- "Overall, this theory section explains our belief that our experimental results
  are not problem specific and should be applicable to other problems of the same
  class (i.e. other trajectory optimization problems, and perhaps also problems in
  other fields)." (p.29) "We expect that this work will have broader applicability
  beyond low-thrust trajectory optimization and even beyond Astrodynamics in
  general." (p.31)
- Prior art the theory rests on: Szu & Hartley 1987 (Cauchy for Fast Simulated
  Annealing), Tsallis & Stariolo 1996 (Generalized SA) — both on *unconstrained*
  multi-modal spaces (p.22). The paper's *new* contribution is extending this to
  the *bounded + internally-constrained* case (p.22), which is the regime our
  cycler correctors also live in. **So the qualitative finding is plausibly
  importable; the magnitudes are not.**

**Verdict for us:** import the *qualitative* prescription (long-tailed,
bi-directional, zero-centred, small scale; Pareto ≳ Cauchy ≫ Gaussian ≳ Uniform).
Do **not** import the numeric Table-3 winners — single low-thrust 500-D benchmark,
no error bars, no second problem.

---

## 5. Cross-check vs Ceriotti 2010 (does tuning reconcile the MBH underperformance?)

Ceriotti 2010 (`2026-06-07-ceriotti-2010-mga-global-opt-mining.md` §3.1) found
plain **Multi-Start beat MBH** on MGA problems at default settings (EEM Table 3.15:
MS 22→67% vs MBH 1→41% over 20k→160k evals), with MBH's `ρ = 10%` of the variable
range and a uniform/neighbourhood step. Ceriotti's own caveat (his p.113):

> "if the objective function is globally non-convex, i.e. presents multiple similar
> funnel structures, MBH may not be effective."

**Does Englander's tuning explain Ceriotti's MBH underperformance? Partial yes,
with one honest tension.**

- **The supporting half (tuning):** Ceriotti's MBH used a **uniform** step at a
  **large** excursion (ρ = 10% of range). Englander's whole result is that a
  uniform step is the *worst* choice and that its efficiency is fragile to the
  step-size (Figs. 9-12; the uniform/Gaussian "no-tail/normal-tail" walks get
  pinned by boundaries and go sub-diffusive, §5.3). A 10%-of-range *uniform* step
  is precisely the configuration Englander shows is dominated by a long-tailed
  step with a *much smaller* scale (Cauchy scale ≈0.0001-0.02). So Ceriotti's MBH
  was running in exactly the un-tuned regime Englander cures → his MBH result is at
  least partly an artifact of an untuned uniform perturbation, and a long-tailed
  re-run would plausibly close some of the gap to Multi-Start. The Ceriotti note
  already anticipated this (its §3.1: "tune MBH's ρ rather than assume default —
  consistent with the Englander tuning paper").

- **The tension half (problem class):** the two papers test *different problem
  shapes*, and Ceriotti's failure mode is the one Englander does NOT claim to fix.
  Ceriotti's MGA cost surfaces are **multi-funnel** (his Fig. 3.22: many
  near-optimal solutions with very different decision vectors), and he argues MBH
  is structurally weak there because hopping from one incumbent only explores *one*
  funnel. Englander's benchmark is a **single** low-thrust mission where MBH is
  already the production method (it is not a multi-funnel-vs-single-funnel
  comparison; there is no Multi-Start baseline in Englander at all). Long tails
  help a walk *cover more ground / jump constraints* within the search — they make
  the global-reset/hop walk more diffusive — but they do **not** add the
  *Multi-Start-style independent restarts* that Ceriotti's multi-funnel surface
  rewards. (Englander's MBH does have global resets after `Max_not improve`
  stalls, which is restart-like; but the paper does not test whether that, tuned,
  beats Multi-Start on a multi-funnel MGA problem.)

**Reconciliation verdict:** *Plausible but unproven.* Englander's tuning very
likely explains a meaningful share of Ceriotti's MBH underperformance — Ceriotti
ran MBH in the exact untuned-uniform-large-ρ regime Englander shows is
sub-optimal, and a long-tailed small-scale step would diffuse better. But it does
**not fully** reconcile the two: Ceriotti's headline reason (MGA surfaces are
multi-funnel, where single-incumbent hopping is structurally inferior to
independent multi-starts) is a *different axis* that Englander neither tests nor
claims to fix. Honest position: tuning narrows the gap; it does not erase
Ceriotti's multi-funnel argument. Keep a Multi-Start control (per the Ceriotti
note) AND switch MBH's perturbation to long-tailed (per this paper) — they are
complementary fixes, not substitutes.

---

## 6. Concrete recommended changes to `src/cyclerfinder/search/mbh.py` (PLAN ONLY)

Current state (read, not edited): mbh.py has `PERTURBATIONS = ("cauchy",
"gaussian", "uniform")`, default `perturbation="cauchy"`, `perturbation_scale=0.05`,
per-gene relative/absolute sizing, monotonic acceptance, `stop_after_stall` that
*stops* (no global reset), incumbent-only (no archive). The module docstring and
`2026-06-07-mbh-wrapper.md` flag the perturbation spec as UNSOURCED pending THIS
paper.

**Diff-plan (each item flagged for evidence strength):**

1. **Add `"pareto"` (bi-polar Pareto) to `PERTURBATIONS` and make it the documented
   recommended default.** [SOURCED — the paper's headline recommendation, p.31.]
   - `PERTURBATIONS = ("pareto", "cauchy", "gaussian", "uniform")`.
   - Generator (Table 4, p.14): `s/ε · (ε/(ε+r))^(−α)` with `s` a fair ±1 coin,
     `r = rng.uniform(0,1)`, small fixed `ε`. Parameterised by exponent `α`, NOT a
     step-size — needs a new param (`perturbation_alpha`, default ≈1.08, the
     mid-range value the paper's MSD demo used, p.21). The existing
     `perturbation_scale`/`absolute_scale` still apply as the per-gene multiplier
     on top of the unit Pareto draw.
   - **Whether to flip the *default* from `"cauchy"` to `"pareto"`:** the paper
     recommends Pareto, but our mbh.py docstring/gates and the existing free-return
     gate were validated on Cauchy/Gaussian. CAVEAT: flipping the default changes
     behaviour for current callers. Recommend: ADD `"pareto"`, document it as the
     paper's recommendation, but keep `"cauchy"` as the code default until a gate
     re-run on our problems confirms Pareto ≥ Cauchy here (their single low-thrust
     benchmark is not our cycler corrector — §4 caveat).

2. **Re-base the default scale to reflect the long-tail/small-scale coupling.**
   [PARTIALLY SOURCED — direction is sourced, the number is not.] The paper's
   Cauchy *scale* sweep is 0.000125-0.02 (Table 3) — 10-30× smaller than our
   current `perturbation_scale=0.05`. Their tiny scale works because the *tail*
   does the big jumps. FLAG: their values are on EMTG's bound-normalised variables,
   not our raw genes, so the number is NOT directly portable — but the qualitative
   lesson (long-tailed ⇒ use a SMALLER typical scale than you would for Gaussian)
   should be recorded in the docstring and the default for `cauchy`/`pareto`
   reconsidered downward, validated by a gate sweep, not adopted blind.

3. **Update the module docstring + `2026-06-07-mbh-wrapper.md` "SPEC CAVEAT"
   block to SOURCED.** [SOURCED.] Replace "NOT yet acquired / documented sensible
   default" with the citation (Englander & Englander 2014, ISSFD24 S7-3) and the
   verbatim recommendation (bi-polar Pareto > Cauchy ≫ Gaussian/Uniform;
   bi-directional zero-centred long-tailed step). This is the audit-trail fix the
   task targets.

4. **Audit trail must record the distribution + its excursion parameter.**
   [SOURCED need — the paper's whole point is the distribution choice matters.]
   `MBHResult` currently echoes only `rng_seed`. ADD the distribution name and its
   excursion parameter (`scale` or `alpha`) to `MBHResult` so a run is
   reconstructable from the result alone (this repo's auditable-provenance
   culture). Minimal: add `perturbation: str` and `perturbation_param: float`
   (or carry the full kwargs) onto `MBHResult`.

5. **Consider restart-on-stall (not just stop-on-stall).** [SOURCED — Algorithm 1,
   p.9.] The paper's `Max_not improve` triggers a **global reset to a fresh random
   point**, then continues; mbh.py's `stop_after_stall` terminates. A faithful MBH
   would re-seed (uniform random in bounds) and keep its best-so-far. This needs
   bounds (mbh.py is currently bounds-free), so it is a larger change — FLAG as a
   follow-up, not a one-liner. The current stop-on-stall is a defensible
   simplification for our bounded-by-construction seeds; document the divergence.

6. **(Optional, cross-paper) keep `gaussian`/`uniform` for ablation only.** Both
   papers agree these are inferior; retain them for the mechanics gate and
   for reproducing the negative baseline, with a docstring note that they are
   NOT recommended for production hops.

**Do NOT import:** the specific Table-3 winning scale/α numbers as production
defaults (single low-thrust 500-D benchmark, no error bars, EMTG-normalised
variables — §4). Import the qualitative spec + the bookkeeping, validate
magnitudes on our own correctors via a gate sweep.

---

## 7. v4.2 backfill check — NO catalogue-eligible rows

Pure algorithm/methods paper. The only trajectory is a hypothetical low-thrust
EESU mission to Uranus (Table 2, Fig. 3) used as an optimiser benchmark — not an
Earth-Mars cycler, no S/L resonant geometry, no Aldrin family, no periodic
geometry. **Zero catalogue rows.** No `center` / `tof_days_bounds` /
`source_ephemeris` backfill applies. (For completeness: EMTG uses SPICE/JPL NAIF
ephemeris, §2.3 / ref [14] — but nothing here is a reproducible cycler number.)
The value of this paper is purely the sourced MBH perturbation spec for mbh.py.
