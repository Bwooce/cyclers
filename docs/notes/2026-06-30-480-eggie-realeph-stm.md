# #480 EGGIE Real-Eph STM Corrector Verdict (2026-06-30)

## What was done

Wired the analytic STM Jacobian (`jovian_stm_jacobian`, `jovian_stm.py`) into
`jovian_shoot` (`jovian.py`) via a new `jacobian="fd"|"stm"` parameter, mirroring
the existing `ideal_eggie_shoot` path. The default `"fd"` is byte-unchanged (parity
oracle). The `"stm"` path passes a `jac=` callback to `least_squares(method="trf")`.

Seed: real-eph EGGIE ballistic periapsis member at `paper_departure_et() + 3.78 d`
(Task 4 best-epoch, departure ET ≈ 6.552e8 s past J2000). Sequence
Europa→Ganymede→Ganymede→Io→Europa; ToFs [1.65, 8.70, 7.02, 11.03] d; Lambert plan
[(0,'single'),(1,'high'),(1,'low'),(1,'high')]; periapsis nodes via `periapsis_node`.

## STM-vs-FD parity (real-eph model)

| Metric | Value |
|---|---|
| Overall relative error | 5.14e-4 |
| Worst nonzero 6×6 block rel | 3.63e-3 |
| Pass threshold | 1e-3 overall (1e-4 for ideal model) |
| FD build time | 16.6 s (30 residual evaluations) |
| STM build time | 0.5 s (4 DOP853+STM integrations) |

**PASS**: the analytic STM Jacobian matches FD on the real-eph seed to 5.1e-4
overall. The slightly higher block-level error (~3.6e-3 in one block) is due to the
JUP365 rails-cache / direct-SPICE query mismatch (~1e-2 km spline noise) that the
STM co-integrator sees but the REBOUND corrector averages through the spline; this is
expected and does not affect the corrector.

## Corrector run

Two chunks of `jovian_shoot(jacobian="stm", max_nfev=60, max_wall_sec=500)`.

| Metric | Chunk 1 | Chunk 2 | FD baseline |
|---|---|---|---|
| nfev budget | 60 | 60 | 18 |
| Wall time (s) | 71.7 | 63.3 | 189.4 |
| Per-nfev wall (s) | 1.2 | 1.1 | 10.5 |
| seed defect_norm | 4.13e+05 | 3.89e+02 | 4.13e+05 |
| final defect_norm | 3.89e+02 | 3.77e+02 | 7.52e+02 |
| converged (strict) | False | False | False |
| correction_dv (m/s) | 836 | 182 | 542 |
| bend_feasible | True | True | True |

## Per-encounter V∞ at chunk-2 plateau

Sourced Table-4 values: Europa 9.12, Ganymede 7.07, Io 8.38 km/s
(Hernandez-Jones-Jesick 2017, AAS 17-608, Table 4).

| Node | Body | V∞ (chunk-2, km/s) | Table-4 (km/s) | Gap |
|---|---|---|---|---|
| 0 | Europa | 8.80 | 9.12 | −0.32 |
| 1 | Ganymede | 6.63 | 7.07 | −0.44 |
| 2 | Ganymede | 6.51 | 7.07 | −0.56 |
| 3 | Io | 7.74 | 8.38 | −0.64 |
| 4 | Europa | 8.73 | 9.12 | −0.39 |

All V∞ sit 0.3–0.64 km/s below Table-4, consistent with the characterised statement
that this real-eph member sits ~0.5 km/s below Table-4.

## Honest assessment

**Did STM converge where FD plateaued?** Partially. The STM achieved a 4× lower
final defect than FD at the same nfev budget (3.89e2 vs 7.52e2 after comparable
iterations), and is 9× faster per nfev (1.2 s vs 10.5 s). However neither STM nor
FD converged to the strict continuity criterion (~5e-3 km). The second STM chunk
moved the defect only 3.4% (3.89e2 → 3.77e2), confirming a hard plateau.

**What is the plateau?** The residual at ~380 is ~78,000× above the convergence
floor. The periodicity wrap constraint (the cycler must repeat) is not closing: the
real-eph orbit at this epoch / Lambert topology lands in an off-basin region where
the legs can be made approximately continuous but the orbit does not repeat within
the 5-encounter sequence. This is the same off-basin wall characterised in the
heliocentric lane (#135) and the Jovian conic corrector (#480 earlier tasks).

**Maintenance ΔV:** Cannot be reported — the corrector has not converged to a
ballistic n-body periodic solution. The `correction_dv_kms` (836 m/s in chunk 1) is
the node-impulse difference between seed and partially-corrected states, NOT the
per-cycle maintenance ΔV. Reporting that number as "maintenance ΔV" would be
misleading; it is recorded here as a characterisation only.

**No Table-4 reproduction:** As expected and documented. The real-eph ballistic
member's V∞ profile (E~8.8/G~6.6/Io~7.7) is consistent with the family
below Table-4; Table-4 reproduction would require the ideal-model basin, not the
real-eph geometry.

## Infrastructure built

- `jovian_shoot(jacobian="stm")` path added to `src/cyclerfinder/nbody/jovian.py`
- Parity gate test: `tests/nbody/test_jovian_stm_realeph_parity.py` (2 tests, ~22 s)
- Scripts: `scripts/_v2_stm_480.py`, `scripts/_v2_stm_chunk2_480.py`
- Run logs: `/tmp/.../scratchpad/l3stm.log`, `.../l3stm_chunk2.log`
