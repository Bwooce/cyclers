# 2026-06-19 MBH Tuning Results (Englander & Englander 2014)

Following the plan pre-registered in [2026-06-19-mbh-tuning-plan.md](2026-06-19-mbh-tuning-plan.md), an MBH tuning sweep was performed over three free-return genomes across 8 distinct RNG seeds.

## 1. Experimental Setup
- **Genomes**:
  - `mcconaghy-2006-em-k2` (S1L1)
  - `niehoff-visit1` (VISIT-1)
  - `niehoff-visit2` (VISIT-2)
- **Mis-seed**: +0.05 AU for $a$, +0.05 for $e$, -40 days for $t_0$.
- **Sweep**: Cauchy (baseline) vs Pareto ($\alpha \in \{1.01, 1.05, 1.08, 1.12\}$).
- **Execution Script**: `scripts/mbh_englander_tuning.py`

## 2. Results

```text
Genome                    Dist       Alpha    Successes/8     Mean Hops
--------------------------------------------------------------------------
mcconaghy-2006-em-k2      cauchy     nan                8           9.9
mcconaghy-2006-em-k2      pareto     1.01               8           6.2
mcconaghy-2006-em-k2      pareto     1.05               8           8.1
mcconaghy-2006-em-k2      pareto     1.08               8           6.2
mcconaghy-2006-em-k2      pareto     1.12               8           5.9
niehoff-visit1            cauchy     nan                8           2.6
niehoff-visit1            pareto     1.01               8           5.0
niehoff-visit1            pareto     1.05               8           3.6
niehoff-visit1            pareto     1.08               8           4.2
niehoff-visit1            pareto     1.12               8           4.2
niehoff-visit2            cauchy     nan                8           0.0
niehoff-visit2            pareto     1.01               8           0.0
niehoff-visit2            pareto     1.05               8           0.0
niehoff-visit2            pareto     1.08               8           0.0
niehoff-visit2            pareto     1.12               8           0.0
```

Algorithm 1 (restart_bounds=15) evaluation on `mcconaghy-2006-em-k2` (Pareto 1.08):
- **restart_bounds=None**: successes: 8/8, mean_hops: 6.2
- **restart_bounds=15**: successes: 8/8, mean_hops: 7.1

## 3. Analysis vs Win Criterion

The pre-registered win criterion to flip the default to Pareto was:
> *`pareto` becomes default only if it beats `cauchy` on success-rate across ≥8 seeds, on ≥2 of the 3 benchmark genomes, by a margin > 1σ of the seed spread.*

**Outcome:** `cauchy` achieved a 100% success rate (8/8) across all 3 benchmark genomes from the mis-seed offset. `pareto` also achieved 100%. Because Pareto did not beat Cauchy in success rate on any genome, the pre-registered win criterion is **failed**. 

While Pareto showed slightly faster convergence on the S1L1 genome (5.9-8.1 hops vs 9.9 hops), it was slower on VISIT-1 (3.6-5.0 hops vs 2.6 hops). Both proved robust against the mis-seed over all RNG seeds evaluated.

## 4. Conclusion
The default perturbation distribution in `src/cyclerfinder/search/mbh.py` will **remain Cauchy**. No capability gap or robustness deficiency was identified in the current configuration.
