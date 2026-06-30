# #319 V3_qp gauntlet built + #320 QP-tori SILVERs PASS V3

**Date:** 2026-06-30. Built the deferred **V3_qp** gauntlet (the #319 Phase-2 / #320
follow-on infra blocker) and ran the two genuinely-fresh #320 QP-tori SILVERs through it.
**Both PASS V3_qp** — the torus invariance signature V1/V2 reported is integrator-
independent (REBOUND IAS15 reproduces the scipy DOP853 result), not a corrector artifact.
No catalogue writeback (V4 + human review still follow).

## Why this (the #320 reframe, 2026-06-30)
The #320 candidate adjudication established that the Saturn/Pluto repeated-moon candidates
are V0-known (Russell-Strange census) and Neptune doesn't close — so the **only genuinely-
fresh AND closing #320 candidates are the Earth-Moon QP-tori SILVERs** (brackets 2 & 10,
k=4), a distinct class (quasi-periodic invariant 2-tori, not the published moon-cycler
family). They passed V1_qp + V2_qp but were blocked on V3_qp, which #319 Phase 1 deferred.
This builds it.

## V3_qp — what it is (`src/cyclerfinder/data/validation/v3_qp.py`)
The QP-torus analogue of the strict-periodic V3's IAS15 cross-check (`v3_3d_periodic`):
re-evaluate the V2_qp invariance equation under an INDEPENDENT integrator.
- V1_qp / V2_qp propagate the invariant-circle samples with `cr3bp.propagate` (scipy
  DOP853). V3_qp propagates with `_ias15_propagate_cr3bp_rotating` (**REBOUND IAS15**,
  rotating-frame callback; LSODA fallback). Same CR3BP model, different integrator family.
- Two gated metrics: (1) the **IAS15 invariance residual** — does `phi^IAS15_{(k+1)t_strob}
  (u(theta))` still lie on `u(theta + (k+1)rho)`, bounded ≤ the V2 drift floor (5e-2); (2)
  the **IAS15-vs-DOP853 terminal disagreement** — bounded ≤ a 1e-3 agreement floor (one
  order under the V2 floor), proving V2's bounded invariance is real, not DOP853 noise.
- Tests `tests/data/test_v3_qp.py` (4 fast: floors-sourced, arg-validation, zero-amplitude
  limit PASS, audit fields; 1 slow: the #299 sourced smoke torus PASS). ruff + mypy clean.

## Result (`scripts/run_320_qp_tori_v3_gauntlet.py` → `data/silver_320_qp_v3_verdicts.jsonl`)
Both SILVER tori regenerated from their #299 bracket parent state, V2 re-confirmed, V3 run:

| bracket | k | V2 max drift (nondim) | V3 IAS15 invariance (nondim) | IAS15-vs-DOP853 disagreement | V3_qp |
|---|---|---|---|---|---|
| 2 | 4 | 3.657e-3 | 3.703e-3 (floor 5e-2) | 9.06e-6 (floor 1e-3) | **PASS** |
| 10 | 4 | 2.795e-2 | 2.573e-2 (floor 5e-2) | 2.53e-7 (floor 1e-3) | **PASS** |

The IAS15 invariance residual tracks V2's DOP853 invariance to ~1% on both, and the two
integrators agree to 9e-6 / 3e-7 nondim (~100-10^5× under the agreement floor). The torus
invariance is a real CR3BP property, integrator-independent.

## Standing / discipline
- **Both #320 QP-tori SILVERs are now V1+V2+V3 confirmed** — the project's only fresh +
  closing #320 candidates have cleared the same-model gauntlet through V3.
- **No catalogue writeback / no admission.** V3 PASS is necessary-not-sufficient. Remaining
  gates before any admission: **V4_qp** (real-physics / higher-fidelity model for a 2-torus —
  the QP analogue of #332's HFEM, still to be defined/built) and the **live lit-novelty +
  human gauntlet**. These are Earth-Moon 3D CR3BP invariant tori seeded from the
  Antoniadou-Voyatzis 2018 family (#299 brackets) — the lit-novelty check (is THIS torus
  family published?) is the next decisive gate, not yet run live.
- The #319 V3_qp infra is now SHIPPED — unblocks V3 for any future QP-torus candidate, not
  just these two.
