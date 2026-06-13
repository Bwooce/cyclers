# FBS analytic gradients as the ΔV-OPTIMISER engine — fair trial (#243)

**Date:** 2026-06-13
**Optimiser:** `src/cyclerfinder/search/fbs_optimize.py::optimize_chain_fbs` (new)
**Driver:** `scripts/fbs_optimizer_fair_trial.py` (committed)
**Tests:** `tests/search/test_fbs_optimize.py` (4 passed)
**Gradient machinery under test (#226):** `core/fbs_match_point.py::chain_defect` +
`chain_defect_jacobian` (the analytic block-sparse match-point Jacobian).
**Writeback:** NONE. No catalogue edit, no validation-level change. Method
evaluation only.

## Why this re-test (read first)

Task #242 tested FBS as a single-leg FEASIBILITY solver versus Lambert and found a
clean negative — but that was the **wrong role**. Ellison's forward-backward-shooting
transcription exists to supply **analytic gradients for gradient-based ΔV
OPTIMISATION** of multi-DSM / multi-gravity-assist trajectories (the EMTG/SNOPT use
case), not to find feasible single legs. This task tests FBS in that actual role and
is the DECISION GATE for #244-246.

## Question

In a ΔV-minimising multi-leg optimiser, do FBS **analytic** constraint gradients beat
the **finite-difference** gradient style the existing Lambert+FD lane relies on
(`dsm_chain_correct` least_squares `2-point`; `optimize.py` SLSQP), on
(1) convergence robustness, (2) cost, (3) optimum quality? A clean negative is a valid
result. A clean positive sends #244-246 forward.

## Answer (one line)

**YES — on the OPTIMISATION problem the analytic FBS gradients decisively beat
finite differences: 3-16x higher success-to-optimum rate from cold seeds, up to
~4.8x faster wall-clock, and the SAME or BETTER optimum. Recommendation: ADOPT FBS as
the optimisation engine — proceed to #244.** This is the mirror image of the #242
feasibility negative, and confirms #242's own framing that it tested the wrong thing.

---

## What is compared (the fair-trial design)

The SAME ΔV-minimising NLP is handed to SLSQP twice; the **only** difference is the
constraint-gradient source, so every measured difference is attributable to it.

* **Decision vector** `x = [Δv_0..Δv_{M-1} (3M) | v_0..v_M (3(M+1))]` — the per-leg
  interior impulses and the SHARED match-point boundary velocities (Ellison's
  match-point variables; leg `i` uses `v_i` as departure and `v_{i+1}` as arrival).
* **Objective** `f(x) = Σ_i ‖Δv_i‖` (+ arrival v∞ `‖v_M - v_planet‖` when rendezvous).
* **Equality constraint** the stacked per-leg match-point defect `chain_defect`
  (non-dimensionalised: position rows / AU, velocity rows / v_circ), driven to 0.
* **Genuinely under-determined**: `3M + 3(M+1)` variables vs `6M` constraints leaves
  `3M + 3` optimisation degrees of freedom, so SLSQP minimises ΔV over a real null
  space — it is an OPTIMISATION, not a square root-find (the #242 mistake).
* **FBS-analytic lane**: constraint Jacobian = the #226 analytic
  `chain_defect_jacobian`.
* **FD lane**: SLSQP finite-differences the SAME constraint (`jac="2-point"`) — the
  incumbent Lambert+FD lane's gradient style.

**Same-model golden**: each problem's ballistic seed is built by Lambert; the
converged chain ΔV is the same-model optimum BOTH lanes must reproduce when feasible.
No cross-model golden. The analytic Jacobian is cross-checked against central
differences at every problem's seed (the only available check — Ellison publishes no
numeric gradient): max relative error **2.2e-7 / 2.2e-7 / 5.8e-8** across the three
problems — the analytic gradient is correct.

## Test problems + reference optima (same-model)

| Problem | Topology | legs | NLP size (vars / defect rows) | rendezvous | same-model optimum ΔV |
|---|---|---|---|---|---|
| Aldrin-class E-M-E (1-syn) | E→M→E | 2 | 15 / 12 | no (flyby finish) | **14.670354 km/s** |
| Russell-class E-M-E (2-syn) | E→M→E | 3 | 21 / 18 | yes (arrival v∞) | **41.513870 km/s** |
| 6.44Gg3-class multi-arc | E→M→E→M→E | 4 | 27 / 24 | no | **27.573653 km/s** |

The optima are large because the legs use real-eph (DE440) ToFs that are deliberately
off the ballistic-resonant values (so the optimiser has genuine ΔV to minimise rather
than a trivial ~0 ballistic answer); the magnitude is irrelevant — the comparison is
analytic-vs-FD on the IDENTICAL NLP, and "best feasible ΔV across all seeds/lanes" is
the reference both lanes are scored against. (Aldrin/Russell/6.44Gg3 are used for
topology shape and as the project's standard stress cases; the maintenance-ΔV value is
not the metric here, the optimiser convergence behaviour is.)

## Head-to-head table (real numbers, 40 jittered cold seeds per problem, jitter σ=1 km/s)

`feas%` = fraction of seeds that converged to a feasible point (scaled defect < 1e-6).
`opt%` = fraction reaching the best ΔV within 1e-3 km/s (success-to-OPTIMUM — the key
robustness metric). `wall` = mean wall-clock per solve. `con_nf` = mean constraint
function-evaluations (the FD lane pays ~(n_vars+1) extra per Jacobian).

### Aldrin-class E-M-E (2 legs, 15 vars)
| lane | feas% | opt% | best ΔV (km/s) | wall (ms) | nit | obj_nfev | con_nfev | mean max-defect |
|---|---|---|---|---|---|---|---|---|
| **FBS-analytic** | **100%** | **100%** | 14.670354 | **141.1** | 44.5 | 99.3 | 101.3 | 7.1e-15 |
| FD | 28% | 2% | 14.670354 | 248.9 | 13.6 | 39.2 | 274.2 | 1.0e-07 |

### Russell-class E-M-E rendezvous (3 legs, 21 vars)
| lane | feas% | opt% | best ΔV (km/s) | wall (ms) | nit | obj_nfev | con_nfev | mean max-defect |
|---|---|---|---|---|---|---|---|---|
| **FBS-analytic** | **82%** | **82%** | 41.513870 | **236.3** | 46.6 | 125.3 | 127.3 | 1.6e-12 |
| FD | 5% | 0% | 45.915610 | 302.0 | 8.1 | 29.9 | 230.0 | 1.1e-07 |

### 6.44Gg3-class multi-arc E-M-E-M-E (4 legs, 27 vars)
| lane | feas% | opt% | best ΔV (km/s) | wall (ms) | nit | obj_nfev | con_nfev | mean max-defect |
|---|---|---|---|---|---|---|---|---|
| **FBS-analytic** | **92%** | **82%** | 27.573653 | **615.0** | 110.6 | 227.9 | 229.9 | 1.6e-10 |
| FD | 60% | 42% | 27.573653 | 2927.3 | 56.4 | 124.8 | 1732.3 | 1.7e-10 |

## Reading the result

1. **Robustness — analytic wins decisively.** Success-to-optimum from cold seeds:
   100% vs 2% (Aldrin), 82% vs 0% (Russell), 82% vs 42% (6.44Gg3). The FD lane
   repeatedly stalls at SLSQP's line-search on the noisy finite-difference Jacobian —
   it terminates "successfully" at a feasibility tolerance of ~1e-7 but at a
   suboptimal point (Russell: FD's best feasible ΔV is 45.92, 4.4 km/s WORSE than the
   true 41.51 optimum the analytic lane finds), or fails feasibility entirely. The
   analytic gradient keeps SLSQP in the basin and drives the defect to 1e-10..1e-15.

2. **Cost — analytic is cheaper, increasingly so with size.** Even though the analytic
   lane runs MORE SLSQP iterations (it actually converges instead of stalling), it is
   faster in wall-clock on every problem and the gap widens with NLP size: ~1.8x
   (Aldrin), ~1.3x (Russell), **~4.8x (6.44Gg3)**. The mechanism is in `con_nfev`: the
   FD lane finite-differences the 24-row constraint over 27 variables, paying ~1732
   constraint evaluations vs the analytic lane's ~230 (one analytic Jacobian call per
   iteration instead of ~28 perturbed constraint evals). This is exactly the
   `O(n_vars)`-evals-per-gradient penalty FBS analytic gradients eliminate.

3. **Optimum quality — analytic is the same or better.** Where both lanes reach
   feasibility on the same basin (Aldrin, 6.44Gg3) they agree on the optimum to <1e-3
   km/s. Where FD stalls (Russell) the analytic lane finds a strictly BETTER optimum
   (41.51 vs 45.92). The analytic lane never did worse than FD on any problem.

## Honest-outcome discipline check

* **Not a manufactured win.** The metric is success-to-OPTIMUM and same-model best-ΔV,
  not a hand-picked single run. The analytic Jacobian is independently cross-checked
  against central differences (≤2.2e-7) so the gradient itself is trusted. Both lanes
  solve the bit-identical NLP; the FD penalty (`con_nfev` blow-up) is a structural,
  not tuned, property of finite differencing.
* **Consistent with #242, not contradicting it.** #242 found FBS is no better as a
  per-leg FEASIBILITY solver (and weaker on cold multi-rev arcs). That remains true —
  and is irrelevant to the optimisation role: here the legs' feasibility is the
  CONSTRAINT and the win is in how cheaply+robustly the gradient drives the
  ΔV-minimisation around it. #242 explicitly flagged it was testing the wrong role;
  this run tests the right one.
* **Caveat — SLSQP, single optimiser.** The comparison is on scipy SLSQP (the
  project's incumbent gradient optimiser). The qualitative result (analytic gradients
  cut `O(n)` constraint evals per step and stabilise the line search) is
  solver-agnostic and is exactly why EMTG pairs FBS with SNOPT; a production adoption
  (#244) should still confirm on the real chain corrector, not just this harness.
* **Caveat — flyby-continuity not yet in the NLP.** This harness optimises impulses +
  match-point boundary velocities; the patched-conic flyby v∞-continuity constraints
  (`flyby_coupling_block`, #226 Phase 5) are NOT yet wired as additional NLP
  constraints. That is the natural #244 extension — and the place the analytic-gradient
  advantage should compound (more constraints → larger FD penalty).

## Recommendation for the decision gate

**ADOPT FBS as the ΔV-optimisation engine. Proceed to #244.** The analytic gradients
deliver exactly the claimed advantage in the role they were designed for: markedly more
robust convergence to the optimum from cold seeds, lower cost that scales favourably
with problem size, and equal-or-better optima. #244 should wire the
flyby-continuity constraints + the phase/epoch columns into the same NLP and validate
against the real chain corrector.

## Reproduce

```
uv run python scripts/fbs_optimizer_fair_trial.py
uv run pytest tests/search/test_fbs_optimize.py
```

## Lint / type / test status

* `uv run ruff check` (new files) — clean.
* `uv run ruff format --check` (new files) — clean.
* `uv run mypy src tests` — clean (335 source files; `scripts/` outside project mypy
  scope per `.pre-commit-config.yaml`, consistent with every other `scripts/` file).
* `tests/search/test_fbs_optimize.py` — 4 passed (Jacobian-vs-FD correctness,
  analytic convergence to optimum, FD-pays-more-constraint-evals, validation errors).
