# 2026-06-19 MBH Tuning Plan per Englander 2014

## Honest Priority Framing
MBH is already fully implemented with all three of Englander's recommendations available (Cauchy, Pareto, restart-on-stall). This task is default-distribution tuning (polish), not a capability gap. It won't find new cyclers; at best, it improves yield/robustness on bounded-genome closures. It is worth doing but is **lower-leverage than #399 (SPK admissions, real catalogue growth) or #392 (Floquet re-scope, the discovery frontier)**. This plan is queued and should not preempt those higher-leverage tasks.

## Proposed Analysis 

### 1. The Benchmark Suite
We will test the tuning across a suite spanning easy to hard closures to see where the heavy tail matters. To isolate MBH wrapper mechanics cleanly while honoring the requirement for 3 distinct genomes, we will use three free-return testbeds across 8 unique mis-seeds:
- **Genome 1:** `mcconaghy-2006-em-k2` (S1L1, a=1.30, e=0.257, 4.27 yr synodic cycle)
- **Genome 2:** `niehoff-visit1` (VISIT-1, a=1.17, e=0.193, 14.95 yr synodic cycle)
- **Genome 3:** `niehoff-visit2` (VISIT-2, a=1.52, e=0.34, 14.95 yr synodic cycle)

### 2. Parameter Sweeps
- **Unfreeze variables:** Unfreeze all relevant genes (e.g., `a` and `e` in free-return) and apply appropriate relative/absolute scales.
- **Alpha Sweep:** Sweep the Pareto exponent `alpha` over the paper's recommended range (1.01 to 1.12).
- **Algorithm 1:** Evaluate the default-off restart-on-stall (`restart_bounds`) global reset.

### 3. Primary Metrics & Seed-Sensitivity Guard
Every metric must be reported as `mean ± σ` over multiple RNG seeds, never from a single run.
- **Primary Metric:** Success-rate-across-seeds (robustness). False negatives are the enemy; a config that is faster but misses closures is worse.
- **Secondary Metric:** Time-to-first-recovery.

### 4. Pre-registered Win Criterion
To avoid p-hacking the optimizer config, flipping the default distribution from `cauchy` to `pareto` requires meeting this pre-committed bar:
> **`pareto` becomes default only if it beats `cauchy` on success-rate across ≥8 seeds, on ≥2 of the 3 benchmark genomes, by a margin > 1σ of the seed spread.**

### 5. Deliverables
- `scripts/mbh_englander_tuning.py` to run the above benchmark.
- A comprehensive experimental results note saved to `docs/notes/2026-06-19-mbh-tuning-results.md`.
- If the pre-registered win criterion is met, we will flip the default perturbation in `src/cyclerfinder/search/mbh.py` to `"pareto"`.
