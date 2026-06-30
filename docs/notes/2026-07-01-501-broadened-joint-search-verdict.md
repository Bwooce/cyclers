# #501 Broadened Galilean Joint-Search Verdict (2026-07-01)

**Date:** 2026-06-30. **Campaign:** #501 broadened real-eph discovery. **Verdict:** EMPTY (all sequences).

## Configuration

- Epoch window: 2030-01-01 to 2040-01-01
- n_revs range: (1, 3)
- Powered flyby floor: 100.0 km
- Cells per sequence: 512 (Sobol, seed=0, scrambled)
- Top-K shoot per sequence: 5
- Total cells: 3072 (6 sequences x 512)
- git: `d0c8def`

## Positive Control

- Liang Member D (CGCEC, 2033-09-25): **PASSED**

## Per-Sequence Results

| Tag | Sequence | Cells | Prefilter OK | Feasible | Shot | Closed | Sweep s | Shoot s |
|---|---|---|---|---|---|---|---|---|
| EGE | `Europa→Ganymede→Europa` | 512 | 512 | 77 | 5 | 0 | 14.4 | 346.1 |
| GCG | `Ganymede→Callisto→Ganymede` | 512 | 512 | 71 | 5 | 0 | 16.9 | 712.0 |
| EGCE | `Europa→Ganymede→Callisto→Europa` | 512 | 512 | 8 | 5 | 0 | 48.5 | 512.0 |
| IEI | `Io→Europa→Io` | 512 | 512 | 51 | 5 | 0 | 19.2 | 792.3 |
| IEGI | `Io→Europa→Ganymede→Io` | 512 | 512 | 4 | 4 | 0 | 44.7 | 542.7 |
| EGCGE | `Europa→Ganymede→Callisto→Ganymede→Europa` | 512 | 512 | 2 | 2 | 0 | 113.3 | 108.5 |

**Totals:** 3072 cells, 213 feasible, 26 shot, 0 closed.  Total elapsed: 3272s.

## Top Feasible Prefilter Survivors (per sequence, best by defect_ms)

**EGE:**

| Rank | Epoch | n_revs | branches | defect_ms | cycle_tof_d | min_alt_km |
|---|---|---|---|---|---|---|
| 1 | 2034-07-28 | [3, 3] | ['high', 'low'] | 0.000 | 49.1 | 399.6 |
| 2 | 2036-10-26 | [2, 2] | ['high', 'low'] | 0.000 | 24.8 | 1148.9 |
| 3 | 2035-12-29 | [3, 3] | ['low', 'high'] | 0.000 | 67.1 | 611.2 |
| 4 | 2032-07-04 | [3, 3] | ['high', 'low'] | 0.000 | 35.6 | 34935.8 |
| 5 | 2030-06-06 | [2, 1] | ['low', 'high'] | 0.000 | 53.6 | 4864.8 |

**GCG:**

| Rank | Epoch | n_revs | branches | defect_ms | cycle_tof_d | min_alt_km |
|---|---|---|---|---|---|---|
| 1 | 2032-06-19 | [3, 1] | ['high', 'low'] | 0.000 | 62.5 | 3735.6 |
| 2 | 2036-03-03 | [2, 3] | ['high', 'high'] | 0.000 | 65.3 | 1261.8 |
| 3 | 2033-08-17 | [1, 2] | ['high', 'low'] | 0.000 | 45.7 | 537.0 |
| 4 | 2036-06-05 | [3, 1] | ['low', 'high'] | 0.000 | 55.5 | 2517.6 |
| 5 | 2030-12-05 | [1, 3] | ['high', 'low'] | 0.000 | 65.8 | 1915.2 |

**EGCE:**

| Rank | Epoch | n_revs | branches | defect_ms | cycle_tof_d | min_alt_km |
|---|---|---|---|---|---|---|
| 1 | 2031-09-30 | [2, 1, 1] | ['low', 'high', 'low'] | 0.000 | 95.7 | 2000.8 |
| 2 | 2036-02-21 | [3, 2, 3] | ['high', 'high', 'low'] | 0.000 | 93.1 | 137.7 |
| 3 | 2030-06-13 | [3, 1, 1] | ['low', 'low', 'high'] | 0.000 | 53.2 | 1035.1 |
| 4 | 2032-04-20 | [3, 3, 3] | ['low', 'low', 'high'] | 0.000 | 106.1 | 1781.1 |
| 5 | 2036-10-13 | [1, 2, 1] | ['high', 'high', 'low'] | 0.000 | 84.5 | 9720.2 |

**IEI:**

| Rank | Epoch | n_revs | branches | defect_ms | cycle_tof_d | min_alt_km |
|---|---|---|---|---|---|---|
| 1 | 2031-05-28 | [1, 2] | ['high', 'low'] | 0.000 | 31.6 | 737.5 |
| 2 | 2034-07-28 | [3, 3] | ['high', 'low'] | 0.000 | 24.6 | 1183.5 |
| 3 | 2030-10-14 | [1, 3] | ['low', 'high'] | 0.000 | 14.3 | 574.2 |
| 4 | 2032-08-17 | [2, 2] | ['high', 'low'] | 0.000 | 31.8 | 1131.3 |
| 5 | 2034-12-29 | [2, 3] | ['high', 'low'] | 0.000 | 33.5 | 1609.7 |

**IEGI:**

| Rank | Epoch | n_revs | branches | defect_ms | cycle_tof_d | min_alt_km |
|---|---|---|---|---|---|---|
| 1 | 2031-01-15 | [2, 2, 2] | ['low', 'low', 'high'] | 0.000 | 67.2 | 1742.5 |
| 2 | 2030-04-02 | [2, 3, 2] | ['high', 'high', 'low'] | 0.000 | 49.7 | 1183.6 |
| 3 | 2032-09-21 | [2, 3, 2] | ['high', 'low', 'high'] | 0.000 | 77.1 | 146.1 |
| 4 | 2035-04-18 | [2, 1, 1] | ['low', 'high', 'low'] | 0.000 | 58.1 | 100.3 |

**EGCGE:**

| Rank | Epoch | n_revs | branches | defect_ms | cycle_tof_d | min_alt_km |
|---|---|---|---|---|---|---|
| 1 | 2037-09-08 | [1, 1, 1, 1] | ['high', 'low', 'high', 'low'] | 0.000 | 145.1 | 103.4 |
| 2 | 2038-08-21 | [3, 3, 2, 3] | ['low', 'high', 'low', 'low'] | 1099.777 | 145.8 | 502.0 |

## N-body Shoot Results (all sequences)

| Tag | Rank | Epoch | Converged | defect_norm | correction_dv_kms | wall_s |
|---|---|---|---|---|---|---|
| EGE | 1 | 2034-07-28 | False | 455031.904417726 | 13.849456357717811 | 57.7 |
| EGE | 2 | 2036-10-26 | False | 2699.8070905088 | 5.300115102588087 | 45.5 |
| EGE | 3 | 2035-12-29 | False | 95133.26798328604 | 2.3135430442573384 | 68.0 |
| EGE | 4 | 2032-07-04 | False | 750188.3913303658 | 0.5642899548163336 | 125.9 |
| EGE | 5 | 2030-06-06 | False | 299277.87997113453 | 9.909993316847538 | 48.9 |
| GCG | 1 | 2032-06-19 | False | 308691.05178947886 | 3.1338662894272153 | 49.2 |
| GCG | 2 | 2036-03-03 | False | 1738813.5419091538 | 21.980435898534168 | 134.6 |
| GCG | 3 | 2033-08-17 | False | 928909.182516883 | 44.3271662004154 | 93.3 |
| GCG | 4 | 2036-06-05 | False | 163344.74986253088 | 12.809221310853607 | 340.9 |
| GCG | 5 | 2030-12-05 | False | 551458.4143890633 | 6.684260122909293 | 94.0 |
| EGCE | 1 | 2031-09-30 | False | 1976.7868669528668 | 2.054973497501818 | 93.1 |
| EGCE | 2 | 2036-02-21 | False | 1441154.6967421107 | 7.8108692289372925 | 135.8 |
| EGCE | 3 | 2030-06-13 | False | 75616.7429989014 | 1.0916207654375787 | 103.3 |
| EGCE | 4 | 2032-04-20 | False | 1645231.3977945377 | 19.066246639700104 | 86.1 |
| EGCE | 5 | 2036-10-13 | False | 2157823.0121483267 | 1.1195194579168861 | 93.7 |
| IEI | 1 | 2031-05-28 | False | 3172.738578071989 | 1.8512113783553446 | 44.9 |
| IEI | 2 | 2034-07-28 | False | 12144.716349393388 | 3.299166513453076 | 59.3 |
| IEI | 3 | 2030-10-14 | False | 329899.29471982404 | 0.6005095398529243 | 43.0 |
| IEI | 4 | 2032-08-17 | False | 605154.7862116889 | 45.806296524915226 | 593.2 |
| IEI | 5 | 2034-12-29 | False | 488633.5651458122 | 2.693476565658837 | 51.9 |
| IEGI | 1 | 2031-01-15 | False | 130033.16973435217 | 8.92824526717298 | 67.9 |
| IEGI | 2 | 2030-04-02 | False | 1095090.6256884024 | 3.905029817443607 | 243.9 |
| IEGI | 3 | 2032-09-21 | False | 175873.1868805768 | 1.8382459432957852 | 144.1 |
| IEGI | 4 | 2035-04-18 | False | 1507452.6672436711 | 4.502584976274784 | 86.8 |
| EGCGE | 1 | 2037-09-08 | False | 1010027.1225539899 | 26.214180799880175 | 106.7 |
| EGCGE | 2 | 2038-08-21 | False | N/A | N/A | 0.0 |

## Interpretation

0 / 26 shot cells closed across all 6 sequences.

The broadened #501 campaign (EGE, GCG, EGCE, IEI, IEGI, EGCGE) over epoch 2030-01-01-2040-01-01, n_revs (1, 3), is **compute-bounded empty** at the 3072-cell budget.

Each sequence is registered as a compute-bounded empty region in `data/empty_regions.jsonl`.
Future scale-up, surrogate-based importance sampling, or a qualitatively different global method is required to break the empty-region wall.

## Empty-Region Summary

Registered 6 compute-bounded empty regions in `data/empty_regions.jsonl`.

| Tag | Sequence | Epoch Window | n_revs | cells | verdict |
|---|---|---|---|---|---|
| EGE | `Europa→Ganymede→Europa` | 2030-01-01-2040-01-01 | (1, 3) | 512 | EMPTY |
| GCG | `Ganymede→Callisto→Ganymede` | 2030-01-01-2040-01-01 | (1, 3) | 512 | EMPTY |
| EGCE | `Europa→Ganymede→Callisto→Europa` | 2030-01-01-2040-01-01 | (1, 3) | 512 | EMPTY |
| IEI | `Io→Europa→Io` | 2030-01-01-2040-01-01 | (1, 3) | 512 | EMPTY |
| IEGI | `Io→Europa→Ganymede→Io` | 2030-01-01-2040-01-01 | (1, 3) | 512 | EMPTY |
| EGCGE | `Europa→Ganymede→Callisto→Ganymede→Europa` | 2030-01-01-2040-01-01 | (1, 3) | 512 | EMPTY |

## Method Notes

- Positive control (Liang Member D, CGCEC 2033-09-25): PASSED prefilter.
- Prefilter: `evaluate_joint_cell` (Nelder-Mead patched-conic, JUP365 real-eph).
- N-body shoot: `jovian_shoot(jacobian='stm')` -- analytic block-bidiagonal STM Jacobian.
- Convergence: Jones AAS 17-577 s2.5 floors (1e-3 km / 1e-6 km/s per component).
- Sequences swept: EGE, GCG, EGCE, IEI, IEGI, EGCGE.
- NO catalogue writeback. NO novelty claim.
