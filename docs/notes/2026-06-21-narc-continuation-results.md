# 2026-06-21 — #388 N-arc real-ephemeris continuation batch results

**Verdict: genuine DE440 negative (held, no writeback).** The walk now reaches
the DE440 (lambda=1) endpoint for every descriptor-bearing catalogue row, and
**no row closes**: residuals floor at 1.4–2.5 km/s (tol 0.1) and the emerged V∞
magnitudes are far below the sourced anchors. This is a TRUE physics result — the
parent reached DE440 — not the solver-config crash reported in the prior version
of this note.

## Correction to the prior (crashed) run

The earlier batch reported `0 converged` with `max_residual_kms = inf` and empty
`emerged_vinf_kms` for every row. That was a **solver-config / driver bug**, not a
closure result:

- `narc_continuation_correct` re-seeded each homotopy rung with the corrector's
  FULL ToF list, which re-inserts the eliminated slack leg. So the free-var vector
  grew by one leg per rung (4 → 5 → 6 …) until it outnumbered the residuals (5 in
  vector mode for E-E-M-M) and `least_squares(method="lm")` raised
  *"Method 'lm' doesn't work when the number of residuals is less than the number
  of variables"* on rung 2. The driver caught it as divergence and aborted the
  epoch, so **the walk never reached lambda=1**.

Diagnosed counts at the failing rung (`russell-ch4-4.991gG2`, vector mode):
`len(x0)=6`, `n_res=5` → m<n → lm raises. (Rung 0: x0=4, res=5, OK. Rung 1: x0=5,
res=5. Rung 2: x0=6, res=5 → crash.) Root cause = the per-rung ToF-list growth,
surfaced as the lm m<n error.

Fix (commit `cc4a2c4`): pin a fixed slack leg (longest seed leg), feed the
corrector only the free legs, and strip the reconstructed slack leg on every
re-seed so the free-var count is stable; also drive the corrector with
`method="trf"` (handles m≤n) for robustness. The lane now walks all rungs to
DE440 for every row.

## Run configuration

- Script: `scripts/narc_continuation_batch.py`
- Lane: `src/cyclerfinder/search/narc_continuation.py` (#388)
- Runlog: `data/runs/narc-continuation-20260620T191235Z.jsonl`
- Run wall time: ~47 s total (all 4 rows).
- Launch-window setting: `range(1, 22)` (Russell §5.3 LaunchWindow 1..21, grid=100).
- `target_phase = 0.0` (documented v1 simplification: the parent's exact beginning
  E-M relative phase is not exposed by the idealized Cycler; the LaunchWindow scan
  over 21 candidate epochs is relied on to find a closing epoch).
- V∞ cap: `VINF_CEILING_KMS["M"]` (the inner corrector takes a single scalar; Mars
  is the binding ceiling of the E/M encounters here).
- residual_mode = `vector` (bend feasibility steered inside the residual; 5
  residuals for E-E-M-M).

## Per-row outcome (post-fix, reaches DE440)

| id | V-level | converged | max residual (km/s) | emerged V∞ (km/s) | anchor E / M (km/s) | anchor-match | bend-feasible |
|---|---|---|---|---|---|---|---|
| mcconaghy-2006-em-k2  | V0 | False | 2.453 | [1.64, 2.50, 1.25, 0.02] | 4.70 / 5.00  | False | False |
| russell-ch4-4.991gG2  | V3 | False | 2.392 | [2.07, 2.98, 1.20, 0.00] | 4.99 / 5.10  | False | False |
| russell-ch4-8.049gGf2 | V3 | False | 1.369 | [2.46, 3.15, 1.88, 1.29] | 8.05 / 10.02 | False | False |
| russell-ch4-9.353Gg2  | V1 | False | 2.383 | [2.19, 3.14, 1.20, 0.00] | 9.35 / 10.52 | False | False |

Batch summary: descriptor rows 4, converged 0, anchor-matched 0,
proposed-promotions 0 (held).

## Reading the negative

Every row now relaxes onto the true DE440 ephemeris but lands a non-closed
configuration:

- **Residuals floor at 1.4–2.5 km/s** — an order of magnitude above the 0.1 km/s
  closure tolerance. The V∞-continuity / periodicity-closure residual cannot be
  driven to zero from the idealized parent's geometry at these epochs.
- **Emerged V∞ magnitudes (~1.2–3.2 km/s) sit far below the sourced anchors**
  (E 4.7–9.4, M 5.0–10.5). The relaxed chain is a low-energy, non-anchor-matching
  configuration, not the published high-energy cycler.
- The trailing near-zero V∞ on the final encounter (0.0–0.02 for three rows)
  indicates the slack leg often reconstructs to a near-degenerate ToF — the
  homotopy has walked the chain off the parent family rather than onto a DE440
  closure of it.

This matches the standing #388 / S1L1 picture: a single-ellipse-per-leg N-arc
homotopy from an idealized planar parent does not, by itself, relax onto the real
multi-arc cycler at a sourced V∞ anchor. The crash previously masked this; the
result is now genuine.

## Held verdict

**HELD, no writeback.** No row closes (converged False for all 4; anchor-match
False for all 4; bend-feasible False for all 4), so there is no PROPOSED V0→V1 —
including for `mcconaghy-2006-em-k2` and the V3 regression rows
`russell-ch4-4.991gG2` / `russell-ch4-8.049gGf2`. `data/catalogue.yaml` and
`validate.py` were not touched.

This is the decisive #388 result: with the solver crash removed, the N-arc
real-ephemeris continuation lane reaches DE440 for every descriptor-bearing row
and produces a clean, characterized negative — the next iteration is gated on the
seed/homotopy modelling (multi-arc legs, anchor-targeted continuation), not on a
solver bug.

## Post-run analysis — "why can't we replicate Russell?" (2026-06-21)

Three follow-up findings sharpen the diagnosis. The negative is NOT "Russell's
result is unreproducible"; it is a **seeding / family-selection wall**, and the
faithful machinery to test it already exists.

1. **The target model is NOT the cause.** Re-ran the continuation targeting
   Russell's own accurate model — the fully-ramped mean-element-J2000 backend
   `continuation.ramped_ephemeris(1.0, 1.0)` (real e,i; a,Ω,ω,ν at J2000 mean) —
   instead of DE440. Result is the SAME for all three probed rows: res ≈ 2.35–2.46
   km/s, same low-energy off-anchor basin, emerged V∞ ≈ [1.4–1.8, 2.4–2.6, 1.2,
   0.0]. So continuing to DE440 (perturbed) vs Russell's mean-element model makes
   no difference — the divergence from Russell is not the model.

2. **We were not running Russell's actual corrector.** `correct.ballistic_correct`
   is Russell §5.4.4 *simplified*: magnitude-continuity + `[t0, ToFs]` free vars.
   Russell's real solve is **5n multiple-shooting** — per-leg full v∞ *vector* +
   t0 + tf, with position-match + Δv=0 constraints, from a seven-cycle seed at a
   phase-matched epoch, under SNOPT elastic mode. The magnitude-only reduced-DOF
   version admits the degenerate low-energy family member (trailing v∞→0) that
   Russell's full-DOF solve + seed does not.

3. **The full corrector already exists, and this was already diagnosed (#135).**
   `nbody/shooter.py` IS the full-state multiple-shooting corrector (the Jones/SNOPT
   analogue). Its docstring records the binding **#135 verdict**: the full shooter
   "closes geometrically but lands OFF-ANCHOR (V∞ 9-28 vs sourced 3-10 km/s); the
   verdict is **seeding/basin, not solver deficiency**," and it "MUST be seeded from
   the #133 near-miss conic survey, NEVER from a blind scan." The S1L1 saga
   ([[project_s1l1_realeph_closure_blocker]]) confirms even the complete stack fails
   on family selection.

**Conclusion.** The reason we don't land Russell's published cyclers is the **seed
basin / family selection**, not the model and not a missing corrector. The full
shooter exists; the only lever never tried is feeding it the *constructed,
golden-validated Russell ψ generic-return parent* (new this session) instead of the
searched #133 near-miss survey. That is the next build: a seed bridge
(Russell parent v∞-vector nodes → `nbody/shooter.seed_from_conic` → `shoot()` on
the real n-body model). It is the strongest remaining shot at breaking the
seeding/family-selection wall — or the most well-earned confirmation that these
specific cyclers are not landable by our stack even from the literal parent. Spec
to follow: `2026-06-21-shooter-russell-seed-*`.
