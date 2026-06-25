# #333 QP-GMOS 2-Tori Family Continuation — Harvest Note

**Date:** 2026-06-26
**Status:** Phase-2 continuator built + run. Report-only — NO catalogue writeback, NO novelty claim.
**Lineage:** #290 Phase 1 (single GMOS torus genome) → #319 (V0–V2_qp gauntlet) → #320 (first sweep, 2 SILVER Earth-Moon tori) → **#333 Phase 2 (pseudo-arclength FAMILY continuator)**.

## What was built

`src/cyclerfinder/genome/qp_tori_arclength.py` — the GMOS 2-torus pseudo-arclength
**family** continuator, mirroring the proven `er3bp_continuation` arclength walker.
It promotes the Jacobi constant `C_J` into the augmented unknown vector, replaces
the fold-blind amplitude pin with a pseudo-arclength row + energy-tie row (phase
pin kept), computes the SVD null tangent, and corrects onto the augmented system
with a Levenberg-Marquardt damped Gauss-Newton step. Continuation runs in a
**variable-normalized coordinate** `w = z / scale` so a single `ds` is a uniform
fractional step across the badly-scaled mode (~5e-4) vs energy/period (~1–6)
unknowns. Tests: `tests/genome/test_qp_tori_arclength.py` (Tasks 1–8, all green).

Three correctness fixes were required beyond the plan draft, each evidence-driven
(see the commit messages for the diagnostics):

1. **LM-damped corrector.** The augmented GMOS Jacobian is intrinsically
   ill-conditioned (cond ~5e8 on the #290 smoke torus — a soft mode direction).
   A plain Newton solve diverged; LM damping (matching the Phase-1 `least_squares`
   trust-region behaviour) is required.
2. **Variable normalization.** A raw `ds=5e-3` unit-tangent step blew the N=2
   truncation tail to `|c_2|/|c_1| = 0.3` in one step. Normalizing `z` by its
   per-component scale keeps the mode block moving by `ds·|c_1| ~ 2.5e-6`, holding
   the torus shape stable (tail stays ~1%).
3. **`resonance_tol = 1e-4`, `fold_eps = 1e-2`, `mode_truncation_guard = 0.1`.**
   The #290 seed is a thin near-1/4 torus (drift ~3e-4 from 1/4); the plan-draft
   defaults (1e-3 resonance band, raw fold sign-flip, 1e-4 tail ceiling) would
   falsely flag every member as phase-locked / a fold / truncation-invalid.

## Campaign

```
uv run python scripts/run_333_qp_family.py --seed smoke --ds 5e-3 --max-steps 20
```

Seed provenance: the first accepted #299 Neimark-Sacker bracket (k=4, `step_a`)
off the #296 Braik-Ross (1,1) Earth-Moon family — the same proven #290 smoke seed
(`n_trans=2`, `initial_torus_amplitude=5e-4`). Output: `data/family_333_qp_smoke.jsonl`
(one row per member, written incrementally via the `on_step` callback).

### Result (both directions, 20 steps each)

| quantity            | value                                       |
|---------------------|---------------------------------------------|
| members             | 41 (40 stepped + seed)                       |
| terminated_reason   | `max_steps`                                  |
| folds crossed       | 0                                            |
| resonance crossings | 0                                            |
| C_J range           | 3.1278519 … 3.1278526 (span **6.3e-7**)      |
| rho range           | 1.5686995 … 1.5688500 (span **1.5e-4**)      |
| freq_ratio range    | 0.249666 … 0.249690 (all ≥ 3.1e-4 from 1/4) |
| residual_norm (max) | 7.1e-7                                        |
| independent (max)   | 7.6e-5 (< 1e-4 V1_qp floor)                   |
| all irrational      | yes                                          |

Every member also PASSes the read-only V1_qp gauntlet (Fourier-norm ~4e-7,
off-grid invariance ~3.5e-5).

## Finding: the family is near-iso-energetic (a rotation-number family)

The headline structural result: **this Neimark-Sacker torus family near the 1:4
bracket is essentially iso-energetic AND near-iso-amplitude.** `C_J` varies only
~6e-7 and `|c_1|` only ~3e-7 across the whole 41-member walk; the genuinely-moving
slow coordinate is the **rotation number `rho`** (span 1.5e-4). The family is a
thin tube of quasi-periodic tori winding around the parent (1,1) Earth-Moon orbit
at near-constant Jacobi energy — exactly the leading-order Neimark-Sacker picture
(`rho = 1/k + O(amplitude^2)`, `C_J ≈ C_J(parent)` at small amplitude). The
continuation correctly traverses it (41 valid irrational tori, no folds, no
phase-lock), but `C_J` is a near-degenerate family coordinate here; `rho` (or the
pseudo-arclength) is the informative one. This iso-energetic structure is itself
the finding, and it is why the Task-4/8 capability tests assert arclength/`rho`
motion + energy *containment* rather than a wide `C_J`/amplitude spread.

## Was a fold crossed? Did #290 and #320 SILVER Bracket 2 fall on one family?

- **Fold:** No. 0 folds over 41 members both directions. With the analytic
  energy/phase Jacobian rows + LM corrector de-noising fold detection (the Phase-2
  win over the amplitude stub), the absence of a fold here is a real negative for
  this short iso-energetic walk, not FD-noise masking — a longer walk toward
  larger amplitude (where the tube thickens and the truncation guard eventually
  fires) is where a turning point would first appear.
- **#290 ↔ #320 SILVER Bracket 2 connection (design-draft §5 hypothesis):** NOT
  established by this run. The #290 smoke torus sits at `C_J ≈ 3.1279`; SILVER
  Bracket 2 is at `C_J ≈ 3.0320`. The continued family stayed within `C_J ±3e-7`
  of the #290 seed — it did **not** reach the `C_J ≈ 3.03` neighbourhood, because
  the family is iso-energetic in this coordinate. Connecting the two SILVER tori
  to the #290 seed (if they lie on one family at all) needs a continuation that
  actually moves in energy — i.e. seeding off a *different* bracket, or a much
  longer amplitude-driven walk through the family's thickening regime — not this
  near-bifurcation iso-energetic segment. Left as an open question, not claimed.

## No catalogue admission

No member is self-admitted. A QP-torus is never "novel" until it clears
`search/literature_check.py` against the published record, and the QP lit-check
fingerprint (Olikara/Howell anchors) + a V3+_qp tier are separate deferred work.
Any candidate quasi-cycler surfacing from a future energy-moving family walk must
be flagged for human gauntlet review, not written back here.

## Follow-on (corpus acquisition, deferred)

> No fully-tabulated published Earth-Moon QP-torus family is digested. Acquire and
> digest member coordinates from one of Olikara-Scheeres 2010 (AAS 10-179),
> Olikara 2016 (Purdue PhD), Howell-Howell 2014, or Henderson-Howell 2008 to
> upgrade the #333 capability-golden to a sourced reproduction golden. Per
> `feedback_never_give_up_reproducing_papers`, this is a registered acquisition
> task, not an accepted stop.
