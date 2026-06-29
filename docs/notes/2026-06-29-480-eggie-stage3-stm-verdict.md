# #480 Stage 3 verdict — analytic STM corrector breaks the FD plateau, but EGGIE still does not close

**Date:** 2026-06-29
**Task:** #480 follow-up, Stages 2+3 (analytic-STM Jovian corrector).
**Plan:** `docs/superpowers/plans/2026-06-29-480-eggie-analytic-stm-corrector-plan.md`
**Builds on:** Stage 1 (`d619c44`, analytic state+STM co-integrator), Stage 2
(`9f98bb1`, block-bidiagonal Jacobian), Stage 3 wiring (`193569c`).
**Supersedes the wall diagnosis in:** `docs/notes/2026-06-29-480-eggie-stage2-nbody-verdict.md`.

## One-line verdict — PARTIAL

The analytic block-bidiagonal Jacobian (one STM co-integration per leg, parity-gated
to <=2.4e-6 vs the FD oracle) **breaks the Stage-2 FD-Jacobian plateau** (total |r|
3.4e2 → ~2.0e2) and steers the corrector to a markedly better point: leg *position*
continuity collapses from 2.5–54 km to **0.5–4 km**, the per-encounter V∞ tighten to
within **0.28 km/s of all five Table-4 targets**, and the correction ΔV drops from
4101 m/s to **~1450 m/s**. But the corrector **still does not reach ballistic
closure**: the leg *velocity* continuity plateaus at **~0.06–0.15 km/s** — ~10⁴ above
the 1e-6 km/s floor — and descends only slowly thereafter. So the FD noise floor was
**a** wall (the exact Jacobian clears it), **not the sole wall**. A deeper
seed/discretisation barrier remains.

## Stage 2 gate (PASS) — analytic Jacobian == FD oracle

`jovian_stm_jacobian` vs a finite-difference Jacobian of `jovian_defect_residual`
on the EGGIE periapsis seed (`tests/nbody/test_jovian_stm.py`):

- overall Frobenius rel = **1.4e-6**; worst nonzero 6×6 block = **2.4e-6** (target 1e-4).
- analytic build ~**0.4 s** vs FD oracle ~**15.9 s** (~40×), and one STM/leg vs the
  `6·n_nodes+1` FD re-propagations per LM step.

## Stage 3 numbers (ideal model, `scripts/_stm_close_480.py`, jac="stm")

| quantity | SEED | FD Stage-2 (trf 25) | **STM trf 400** | STM lm 400 |
|---|---|---|---|---|
| total \|r\| | 8.69e5 | 3.39e2 | **1.99e2** | 2.40e2 |
| leg-continuity block | 8.69e5 | 2.82e2 | **1.91e2** | 2.30e2 |
| periodicity-wrap block | 5.99e2 | 1.88e2 | **5.67e1** | 6.93e1 |
| leg position defects | 2.4e4–8.6e5 km | 2.5–54 km | **0.5–4.1 km** | 1.9–9.6 km |
| leg velocity defects | 0.20–10.4 km/s | 0.10–0.19 km/s | **0.057–0.122 km/s** | 0.069–0.149 km/s |
| wrap velocity | 0.47 km/s | 0.19 km/s | **0.057 km/s** | 0.069 km/s |
| correction ΔV | — | 4101 m/s | **1450 m/s** | 1171 m/s |
| converged (leg<floor 4.9e-3) | — | False | **False** | False |
| wall / iters | — | — | 457 s / 400 | 464 s / 400 |

Per-encounter V∞ at the STM trf/400 plateau (km/s, tgt): Europa 9.26 (9.12),
Ganymede 7.29 (7.07), Ganymede 7.00 (7.07), Io 8.66 (8.38), Europa 9.20 (9.12) —
**all within 0.28 km/s of Table 4** (vs Stage-2 FD where Europa drifted to 10.27).
lm and trf agree (both ~2e2, leg vel ~0.06–0.15 km/s).

Provenance of the V∞ targets (sourced, never self-computed): Hernandez-Jones-Jesick
2017 (AAS 17-608) Table 4 — Europa 9.12 / Ganymede 7.07 / Io 8.38 km/s, via
`docs/notes/2026-06-26-digest-hernandez-2017-ieg-triple-cyclers-aas-17-608.md`.

## Descent shape — it is a slow plateau, not pre-quadratic convergence

trf/jac=stm: total drops fast to ~5e2 by nfev≈60 (8.69e5→3.28e3→6.7e2→5.4e2), then
**grinds linearly**: 5.4e2 (nfev 60) → 2.5e2 (nfev ~300) → 2.0e2 (nfev 400) — roughly
halving per ~300 evals, with no late quadratic drop. A true ballistic solution would
show the residual collapsing by squares as the nodes approach a zero-defect arc; this
shallow linear crawl is the signature of converging toward a **non-zero local
minimum** (the ~0.1 km/s velocity-continuity floor), not toward closure. Reaching the
1e-6 km/s floor from here is ~4–5 orders of magnitude away at a linear rate —
infeasible by iteration count alone.

## Interpretation — the residual wall is seed/discretisation, not FD noise

The Stage-2 note hypothesised the wall was the FD-Jacobian noise floor. Stage 3
**falsifies that as the sole cause**: with an exact analytic Jacobian (parity 1e-6)
the corrector still plateaus at ~0.1 km/s velocity continuity. The exact Jacobian
*did* buy a real, large improvement (position 10×, V∞ in-band, ΔV ⅓) — so FD noise
**was** masking part of the basin — but a second, deeper wall sits underneath:

- **Single shooting node per leg over a periapsis-grazing arc (diagnosis pt 1).** Each
  node sits at a flyby periapsis deep in the moon well; one node per leg cannot absorb
  the perijove sensitivity, so the velocity continuity bottoms out at ~0.1 km/s while
  position continuity (less sensitive) closes to ~km. This is a discretisation /
  seed-quality limit, not a solver-gradient limit.

## Status of the #480 reproduction

- **Stage 1 (resonant conic): SUCCESS** (unchanged) — all three V∞ on Table 4.
- **Stage 2 (analytic Jacobian): SUCCESS** — STM Jacobian == FD oracle to 1e-6, ~40×
  cheaper; the documented FD-plateau lever is built and parity-gated.
- **Stage 3 (n-body close): PARTIAL / characterised negative (advanced).** The STM
  corrector breaks the FD plateau and lands a far better, in-band, low-ΔV point, but
  EGGIE is **not** reproduced as a ballistic n-body cycler: a ~0.1 km/s
  velocity-continuity wall remains. No catalogue change; the golden
  (`tests/verify/test_ieg_reproduction_golden.py`) stays skipped.

## Next levers (future work — not this task)

The corrector gradient is no longer the wall; the **discretisation/seed** is. In
priority order:
1. **Sub-arc the legs** (multiple interior shooting nodes per leg) so each defect is
   small and the perijove sensitivity is split — the STM Jacobian already supports
   arbitrary node counts, so this is the highest-leverage next build.
2. **Gravity homotopy** (ramp moon GM 0→full): at GM=0 the Lambert seed closes
   exactly; ramp the flyby gravity on, re-converging with the STM Jacobian each step.
3. Only after an in-band *ballistic ideal* close does the ideal→real homotopy +
   real-eph confirm become worthwhile.

Commits: Stage 2 `9f98bb1`, Stage 3 wiring `193569c`.
Scratch (not committed): `scripts/_stm_parity_480.py`, `scripts/_stm_close_480.py`.
