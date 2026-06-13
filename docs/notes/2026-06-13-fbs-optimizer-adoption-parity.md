# FBS analytic gradients on the REAL corrector — catalogue-wide adoption parity sweep (#244)

**Date:** 2026-06-13
**Driver:** `scripts/fbs_optimizer_adoption_parity.py` (committed)
**Corrector under test:** `src/cyclerfinder/search/dsm_descriptor_seed.py::close_row_dsm`
→ `src/cyclerfinder/search/dsm_leg.py::dsm_chain_correct` (the REAL production chain
corrector, `charge_flyby_continuity=True`), via the `gradient=` opt-in landed in
`4abb1b8` / `c2fb57a`.
**Gradient machinery:** `core/fbs_match_point.py` (`chain_defect` + analytic
`chain_defect_jacobian`, #226) + `search/fbs_optimize_flyby.py` (the flyby-continuity
constraints wired with analytic gradients, `b25f2d1`).
**Writeback:** NONE. The opt-in is off by default; no `data/catalogue.yaml` edit, no
validation-level change. Nothing the catalogue depends on moves. Method evaluation only.

## Why this run (read first)

#243's fair trial proved FBS analytic gradients beat finite differences on an
**ISOLATED** match-point NLP (3–16× success-to-optimum, up to ~4.8× faster). #243 left
two caveats explicit: (1) confirm on the **real** `dsm_chain_correct` lane, not a toy
NLP; (2) wire the patched-conic flyby-continuity constraints (done in `b25f2d1`). This
run closes both. It is the decision gate for **#245 (default flip)** vs **shelve**.

## Question

On the REAL DSM chain corrector, run twice through the bit-identical charged seed and
identical DE440 model (`Ephemeris("astropy")`) — once `gradient="lambert"` (incumbent
Lambert+finite-difference lane) and once `gradient="fbs-analytic"` — does the
FBS-analytic lane reach the Lambert+FD optimum within tolerance (or better) on every
catalogue row that has a multi-leg/DSM formulation? Does the #243 advantage hold here?
A clean negative is a valid result and shelves #245–246.

## Answer (one line)

**The #243 advantage holds DIRECTIONALLY on the real corrector — FBS-analytic reaches a
strictly lower residual on 9/9 rows (median 6.6×, up to 36.6× closer to feasibility) —
but NEITHER lane CONVERGES on ANY of the 9 rows. The binding blocker on the real lane is
the seed/basin (the known multi-arc / off-basin problem), NOT the gradient backbone. So
parity-to-convergence is UNDEFINED, and #245 (the default flip) is NOT yet justified by
a convergence win. Recommendation: do not flip the default on this evidence — HOLD #245
behind a seed-basin fix; keep the FBS lane as the opt-in it already is.**

## Rows swept

All catalogue rows that `seed_dsm_chain_from_descriptor` can seed (a g/G DSM
descriptor) — 9 rows: the McConaghy 2006 Aldrin-family E-E-M-M, four Russell ch.4
2-synodic E-E-M-M rows, and five Russell ch.4 3-synodic E-E-E-M-M rows (incl. 6.44Gg3).
There are no Jones-VEM DSM-seedable rows in the catalogue (none carry a g/G descriptor
the lane can charge), so the Jones set is not in scope of this seedable-rows sweep.

## Jacobian cross-check (caveat 2 — analytic gradient independently trusted)

FD-vs-analytic max relative error of the **flyby-continuity-wired** constraint Jacobian
(`_flyby_constraint_jac_analytic` vs central differences over `_flyby_constraint_vector`)
= **2.71e-11**. The analytic gradient on the real flyby-wired chain is TRUSTED — the
zero-convergence result below is NOT a gradient-correctness artefact.

## Parity table (real numbers, REAL corrector, both lanes identical seed + DE440)

`lcv/fcv` = converged (Lambert / FBS). `*_res` = max match-point residual (km/s; tol =
0.1). `res_x` = lam_res / fbs_res (how many × closer FBS gets). `*_dV` = total DSM ΔV at
the returned (non-converged) point. `*_s` = wall seconds.

| row | seq | lcv | fcv | lam_res | fbs_res | res_x | lam_dV | fbs_dV | lam_s | fbs_s |
|---|---|---|---|---|---|---|---|---|---|---|
| mcconaghy-2006-em-k2 | E-E-M-M | False | False | 40.104 | 7.094 | 5.65 | 102.56 | 14.89 | 4.6 | 21.9 |
| russell-ch4-4.991gG2 | E-E-M-M | False | False | 35.880 | 4.011 | 8.94 | 95.39 | 8.47 | 5.2 | 34.8 |
| russell-ch4-8.049gGf2 | E-E-M-M | False | False | 22.999 | 6.015 | 3.82 | 51.15 | 9.15 | 4.4 | 26.1 |
| russell-ch4-9.353Gg2 | E-E-M-M | False | False | 76.926 | 2.103 | 36.58 | 131.55 | 3.85 | 5.6 | 38.1 |
| russell-ch4-3.64gGg3 | E-E-E-M-M | False | False | 29.531 | 8.796 | 3.36 | 91.99 | 27.20 | 2.5 | 50.8 |
| russell-ch4-3.78Gg3 | E-E-E-M-M | False | False | 19.956 | 4.867 | 4.10 | 64.05 | 16.62 | 0.3 | 30.8 |
| russell-ch4-9.94Gg3 | E-E-E-M-M | False | False | 49.570 | 7.493 | 6.62 | 151.17 | 18.38 | 5.2 | 10.4 |
| russell-ch4-5.30ggF3 | E-E-E-M-M | False | False | 19.955 | 1.998 | 9.99 | 47.82 | 3.29 | 8.5 | 48.5 |
| russell-ch4-6.44Gg3 | E-E-E-M-M | False | False | 33.069 | 3.504 | 9.44 | 74.81 | 11.84 | 5.6 | 24.9 |

Aggregate: 9 rows swept; **lam converged 0, fbs converged 0, both 0, neither 9**;
neither-conv residual signal: **FBS strictly closer on 9/9 rows, median lam/fbs residual
= 6.62×**.

## Disagreement ledger

Because no row converges on either lane, every row is a disagreement under the only
honest definition available here (parity-to-convergence is undefined when neither lane
converges — the script no longer silently scores zero-convergence as "clean parity").
The decision-relevant divergence is the residual gap, tabulated above. All 9 rows
diverge the SAME way: FBS-analytic is markedly closer to feasibility, but both lanes
stop short of the 0.1 km/s convergence tolerance. There is no row where the lanes reach
genuinely different feasible optima, and no row where one converges and the other does
not — because none converge.

## Reading the result (honest)

1. **The #243 advantage is real and reproduces directionally on the real lane.** FBS-
   analytic gets the corrector strictly closer to feasibility on all 9 rows, by a median
   6.6× and up to 36.6× (russell-9.353Gg2: 76.9 → 2.1 km/s residual). The analytic
   gradient does exactly what #243 claimed — it keeps the corrector descending toward the
   match-point manifold far more effectively than the finite-difference Lambert lane.

2. **But it does NOT cross the convergence threshold, because that is not the binding
   constraint.** These multi-leg E-E-M-M / E-E-E-M-M rows are the known off-basin / multi-
   arc cases (cf. the S1L1 saga, MEMORY `project_s1l1_realeph_closure_blocker`): the
   single charged seed lands in the wrong basin, and no gradient backbone fixes a
   wrong-basin seed. A better gradient gets you to the nearest local feasibility
   floor faster — here ~2–9 km/s — but the nearest floor is still not a closed cycler.

3. **Wall-clock: FBS is SLOWER on the real lane, opposite to the isolated #243 trial.**
   #243's NLP harness saw FBS faster (analytic Jacobian replaces O(n) perturbed
   constraint evals). On the real corrector the FBS lane runs a full Ellison match-point
   solve per leg per outer iteration, so it is 4–14× slower in wall-clock (e.g. 0.3 → 31 s
   on russell-3.78Gg3). The #243 cost win does NOT carry over to the production chain
   lane as currently wired — the per-leg FBS solve dominates.

## Honest-outcome discipline check

* **Not a manufactured win, and not a manufactured loss.** The first sweep pass produced
  a misleading "CLEAN PARITY" verdict purely because the disagreement logic gated on
  `lam.converged AND fbs.converged` and silently passed the all-fail case. That logic was
  fixed: zero-convergence now reports "parity UNDEFINED — blocker is seed/basin, not
  gradient backbone", and the residual-improvement signal is surfaced explicitly. The
  numbers above are the corrected, real sweep.
* **Same-model golden discipline held.** Both lanes solve the bit-identical charged seed
  on the identical DE440 model; the only difference is the gradient backbone, so the
  residual gap is attributable to it. No cross-model / published target was used. The
  analytic Jacobian is cross-checked against central differences (2.71e-11).
* **No writeback.** The opt-in stays off by default; the catalogue is untouched.

## Recommendation for the decision gate (#245)

**HOLD / do NOT flip the default (#245) on this evidence.** The justification for #245
was "FBS reaches the Lambert+FD optimum (or better) and is cheaper". On the real
corrector: it is closer-to-feasible on every row (the optimum-quality half holds
directionally) but it is **slower**, and crucially **neither lane actually converges**,
so there is no convergence parity to flip toward — flipping the default would make the
production lane slower with no convergence gain. The FBS lane should remain the **opt-in**
it already is (`gradient="fbs-analytic"`), valuable as a closer-to-feasibility probe and
for any future seed-basin work, but it is not a drop-in default.

The real, decision-changing blocker this run exposes is upstream of the gradient: the
**single-charged-seed basin selection** for multi-arc E-E-M-M / E-E-E-M-M rows. #245–246
should be re-scoped behind a seed/basin fix (multi-start, or the multi-arc seeding the
S1L1 work pointed at), not a default flip. Re-run this exact sweep once a seedable row
actually converges on EITHER lane — that is the first point at which a genuine
convergence-parity comparison (and therefore a default flip) becomes meaningful.

## Reproduce

```
uv run python scripts/fbs_optimizer_adoption_parity.py
```

## Lint / type / test status

* `uv run ruff check scripts/fbs_optimizer_adoption_parity.py` — clean.
* `uv run ruff format --check scripts/fbs_optimizer_adoption_parity.py` — clean.
* `uv run mypy src tests` — clean (338 source files; `scripts/` is outside the project
  mypy scope per `.pre-commit-config.yaml`, consistent with every other `scripts/` file
  including the #243 fair-trial driver).
* No new test module — the script is an evaluation driver; the analytic-gradient
  correctness it relies on is already gated by `tests/search/test_fbs_optimize.py` and
  the flyby-wired Jacobian by its own test module from `b25f2d1`.
