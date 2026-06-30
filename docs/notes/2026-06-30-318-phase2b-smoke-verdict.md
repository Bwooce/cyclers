# #318 Phase 2b — Sobol Smoke Verdict (2026-06-30)

**Date:** 2026-06-30. **Status:** Phase-2b smoke run complete. **Verdict:** EMPTY.

## Configuration

- Sequence: `Callisto→Ganymede→Callisto→Europa→Callisto`
- Epoch window: 2033-01-01 to 2035-01-01
- n_revs range: (1, 2)
- ToF seed range: (15.0, 45.0) d/leg
- Sobol samples: 256 (seed=0, scrambled)
- N-body shot top-K: 10
- git: `d1bbf57`

## Results

| Metric | Value |
|---|---|
| Cells submitted | 256 |
| Prefilter succeeded | 256 |
| Prefilter failed | 0 |
| Feasible (surrogate) | 1 |
| N-body shot | 1 |
| Closed (jovian_shoot) | 0 |
| Sweep wall time | 66.5s |
| Shoot wall time | 0.5s |

## Top-10 Feasible Prefilter Survivors (by closure_defect_ms)

| Rank | Epoch | n_revs | branches | defect_ms | cycle_tof_d | min_alt_km |
|---|---|---|---|---|---|---|
| 1 | 2033-10-01 | [1, 1, 2, 1] | ['low', 'low', 'high', 'high'] | 1702.208 | 114.6 | 105.1 |

## N-body Shoot Results (top-10)

| Rank | Epoch | Converged | defect_norm | correction_dv_kms | wall_s |
|---|---|---|---|---|---|
| 1 | 2033-10-01 | False | N/A | N/A | 0.0 |

## Interpretation

0 / 10 shot cells closed under `jovian_shoot(jacobian='stm')`.
The 256-cell Sobol smoke scan in the CGCEC joint manifold (epoch 2033-01-01-2035-01-01, n_revs (1, 2), ToF (15.0, 45.0) d/leg) is **compute-bounded empty**.

This is consistent with the Phase-2 design doc's honest prior:
> *The realistic Phase-2 outcome is a compute-bounded empty-region map, not a discovery*

Registered in `data/empty_regions.jsonl`. Scale-up or surrogate-based sampling
required to improve coverage beyond this smoke result.

## Method Notes

- Positive control (Liang Member D, 2033-09-25): PASSED prefilter.
- Prefilter: `evaluate_joint_cell` (Nelder-Mead patched-conic, JUP365 real-eph).
- N-body shoot: `jovian_shoot(jacobian='stm')` — analytic block-bidiagonal STM Jacobian,
  DOP853+STM co-integration per leg, avoids REBOUND variational-particle gotcha.
- Convergence criterion: Jones AAS 17-577 §2.5 floors (1e-3 km / 1e-6 km/s per component).
- NO catalogue writeback. NO novelty claim.
