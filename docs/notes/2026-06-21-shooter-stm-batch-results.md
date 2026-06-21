# 2026-06-21 — #388 STM-Jacobian literal-parent shooter batch

**Status: COMPLETE. Verdict = decisive characterized NEGATIVE. HELD — no catalogue
writeback.** All 4 descriptor rows × 3 best-phase epochs were shot to LM
convergence/stall under the now-tractable analytic STM Jacobian. None close to
the sourced family.

## What changed vs the FD-era run

The earlier full-shooter-from-Russell-parent attempt
(`2026-06-21-shooter-russell-results.md`) was **compute-infeasible**: the LM
solver's finite-difference Jacobian is `(6*n_nodes+1)` full-state residuals per
iteration. For the `E-E-M-M` rows (4 nodes → 24 free components) that is 25
residual evaluations × 3 leg propagations each ≈ 75 leg-propagations **per LM
iteration**, ~27 min/iteration. That forced a `K=1` / `max_nfev=2` cap that could
only sample the seed basin, never run to convergence.

This batch used the **analytic block-bidiagonal STM Jacobian**
(`shooter.shoot(jacobian="stm")`, commit `8cfe894`): the multiple-shooting
residual is full-state continuity, so the Jacobian is exactly
`∂c_i/∂node_i = Φ_i` (per-leg STM), `∂c_i/∂node_{i+1} = −I`, assembled from **one
co-integrated variational propagation per leg**. The variational particles carry
the full dynamics incl. the rails perturbers via the gravity-gradient-tensor fix
(commit `8196f55`). Parity vs the FD oracle: rel `5.1e-7`.

The FD-era coverage cap was removed (per `feedback_long_runs_acceptable`): `K=3`
best-phase-error epochs over `LaunchWindow 1..21`, `max_nfev=100`, run to LM
stall, per-(row,epoch) append+fsync checkpoint + `--resume` + a per-evaluation
heartbeat for intermediate visibility (commit `4082f70`).

## STM-vs-FD speedup (measured)

On `mcconaghy-2006-em-k2` (E-E-M-M, ToFs 534/153/1026 d):
- 1 residual = **64 s**; 1 STM Jacobian = **121 s** → **~185 s/LM-iteration**.
- FD equivalent: 25 residuals × 64 s ≈ **1600 s/iteration**.
- **≈ 8.6× faster per iteration** — the lever that made the solve finish. (The
  augmented-leg overhead ~2× a plain leg eats most of the naive `25×` structural
  factor; the net is still an order of magnitude and the difference between
  "infeasible" and "12 shoots in ~2.4 h".)

## Per-row results (best epoch by lowest defect)

Sourced anchors are the per-row E and M V∞ (km/s); a row "matches" only if an
emerged per-encounter V∞ lands within 0.5 km/s of BOTH anchors AND the solve
converges to the SNOPT continuity floor.

| id | V | anchors E/M | conv | best defect | emerged V∞ (km/s) | E/M resid | match | bend |
|----|---|-------------|------|-------------|-------------------|-----------|-------|------|
| mcconaghy-2006-em-k2 | V0 | 4.70 / 5.00 | False | 2.46e9 | 17.6, 19.0, 0.0, 3.0e5 | 4.70 / 5.00 | No | No |
| russell-ch4-4.991gG2 | V3 | 4.99 / 5.10 | False | 2.46e9 | 15.7, 19.5, 0.0, 2.8e5 | 4.99 / 5.10 | No | No |
| russell-ch4-8.049gGf2 | V3 | 8.05 / 10.02 | False | 3.46e9 | 8.4, 33.5, 0.0, 6.3e4 | 0.33 / 1.64 | No | No |
| russell-ch4-9.353Gg2 | V1 | 9.35 / 10.52 | False | **7.4e3** | 30.2, 18.9, 7.6, 16.6 | 1.76 / 2.93 | No | No |

(All 12 (row,epoch) records are in `data/runs/shooter-stm-batch.jsonl`; the table
shows the lowest-defect epoch per row. Each row's three epochs were mutually
consistent.)

## Verdict — decisive characterized NEGATIVE

The full multiple-shooting corrector, seeded from the **literal constructed
Russell parent** and now **run to LM convergence** under the tractable STM
Jacobian, does **not** land the sourced family for any of the 4 SnLm rows. Two
distinct failure modes, both confirming the seed-basin / family-selection wall:

1. **Parent too far from any continuous n-body trajectory (3 of 4 rows:
   mcconaghy-2006-em-k2, 4.991gG2, 8.049gGf2).** Seed defect ~3.5e9; the LM
   reduces it by at most ~30% before stalling (nfev 14–15), and one node runs
   away to 1e5–1e6 km/s. The conic-constructed parent is not a usable seed for
   the full-fidelity corrector on these rows — the basin around it does not
   contain a continuous trajectory the local solver can reach. (The 0.33 E-resid
   on 8.049gGf2 is coincidental on an unconverged stall, not a near-match.)
2. **Geometric closure into the WRONG (high-energy) family (1 of 4 rows:
   9.353Gg2).** Here the corrector DOES drive the defect down 6 orders of
   magnitude (3.5e9 → 7.4e3, all nodes finite, nfev 32–39) — but it relaxes into
   a high-V∞ family member (emerged 7.6–30 km/s vs sourced 9.35/10.52; resid
   1.76/2.93), never the published cycler. This is the textbook #135 signature
   ("closes geometrically but lands OFF-ANCHOR, our V∞ 9–28 vs sourced 3–10"),
   now reproduced under the full STM shoot from the constructed parent.

Both modes are the **seeding / family-selection** obstruction the #135
like-for-like and the S1L1 saga ([[project_s1l1_realeph_closure_blocker]])
already identified — now established with a **converged LM solve** behind it
rather than the FD-era compute-capped probe, and from the **literal Russell
parent** rather than the near-miss survey. The obstruction is NOT solver
performance (the STM made the solve tractable and it still does not close) and
NOT model fidelity (the conic continuation lane reached DE440 and Russell's own
mean-element model and landed the same off-anchor basin,
`2026-06-21-narc-continuation-results.md`). It is the family the local corrector
selects from this seed.

`mcconaghy-2006-em-k2` stays **V0**; the V3/V1 rows are unchanged. No catalogue
writeback (`data/catalogue.yaml` / `validate.py` untouched).

## Where this leaves #388

The real-ephemeris closure thread is now characterized to the bottom on the
direct lane: multi-rev wiring (banked), conic N-arc continuation (reaches DE440,
off-anchor), and the full n-body multiple-shooting corrector from the literal
parent (now tractable via the STM, still off-anchor) all converge on the same
wall. Breaking it needs a **global / family-targeted** method — multi-start MBH
over the seed/direction space, or a homotopy that pins the sourced |V∞| family
through the relaxation — not more local-corrector performance. That is a distinct
research build, not a re-run. The STM Jacobian itself is a durable asset (also
the missing piece #347 flagged for n-body monodromy/eigenvalue work).

## References
- `docs/notes/2026-06-21-shooter-russell-results.md` (the FD-era infeasibility
  this run answered), `docs/notes/2026-06-21-narc-continuation-results.md`
  (conic off-anchor basin).
- `nbody/shooter.py::_stm_jacobian`, `nbody/propagator.py::_install_rails_forces`
  (variational perturber gradient), `scripts/shooter_russell_batch.py`.
- Runlog: `data/runs/shooter-stm-batch.jsonl` (12 records).
- Memory: `feedback_long_runs_acceptable`, `project_dsm_closure_modeljump_blocker`,
  `project_s1l1_realeph_closure_blocker`, `feedback_orbit_closure_discipline`,
  `feedback_never_give_up_reproducing_papers`.
