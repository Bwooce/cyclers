# #593: does the #592 QBCP alpha_6 fix change past search conclusions?

**Date:** 2026-07-14
**Origin:** `docs/notes/2026-07-14-qbcp-alpha6-scaling-fix.md` (#592) — per
[[feedback_bugfix_invalidates_past_searches]], scope whether any past QBCP-based search conclusion
(the #533/#537/#538/#544 task chain) needs re-checking.

## Step 1: bound the error magnitude

`alpha_6` is a genuine time-varying Fourier coefficient (`evaluate_alphas()`), not fixed at 1 in
general use. Sampled over a full synodic period (2000 points):

```
alpha_6 range: 0.99376 to 1.00843
max |alpha_6 - 1|: 0.00843  (~0.8%)
```

The buggy code scaled the WHOLE Newtonian potential gradient (Earth+Moon+Sun) by `alpha_6` instead
of only the Sun term — so the systematic force-law error in the Earth/Moon terms (the DOMINANT
terms near the primaries, where periodic orbits/tori of interest live) was up to ~0.8%. This is a
non-trivial force-law error for a differential corrector — several orders of magnitude larger than
typical corrector convergence tolerances (1e-6 to 1e-11).

## Step 2: empirical before/after comparison (not argument-from-magnitude alone)

Reconstructed the exact pre-#592 `qbcp_eom`/`qbcp_potential_second_derivatives` inline (byte-diffed
against `git show 9e80969:...` to confirm exactness — not a paraphrase), monkeypatched it into
`cyclerfinder.core.qbcp`, and re-ran the SE-L2 torus correction from `scripts/run_533_qbcp_connection.py`
STEP 1 (the one CLEAN, fast-converging case in the whole task chain — a genuine reproducible positive
control) under both the fixed and reconstructed-buggy code, same seed (the BCR4BP-stage seed doesn't
call `qbcp.py` at all, so it's bit-identical going into both runs).

```
FIXED:  residual=1.538e-05  rho=1.4511582716  mean_state=[-0.19945, 0.06816, ~0, 0.06611, 0.18597, ~0]
BUGGY:  residual=2.537e-05  rho=1.4522961     mean_state=[-0.00566, 0.02338, ~0, 0.02218, 0.00586, ~0]

max |delta| in mean state: 0.194   (huge relative to the ~0.01-0.02 scale of the coordinates themselves)
rho delta:                 -1.14e-3  (small, ~0.08%, consistent with the alpha_6 magnitude)
residual delta:            -1.0e-5   (both converge cleanly by their OWN model's criterion)
```

**Both solutions independently satisfy their own model's invariance-residual criterion (~1e-5,
"clean convergence") — but they are NOT the same solution.** The rotation number `rho` is nearly
identical (consistent with the small alpha_6 magnitude), but the torus's actual mean-state location
differs by an amount comparable to the coordinates' own scale, not a small perturbation of it.

## Interpretation

This is NOT the "small, probably-safe" outcome the O(alpha_6-1)~1e-3 framing in #592's own note
might suggest for a casual read. For the SE-L2 torus case specifically, the corrector's basin of
convergence is sensitive enough to the force-law scaling that the buggy and fixed models land on
visibly different points in phase space from the identical seed. Two honest readings, not
distinguished by this test alone:
- The two models have distinct nearby fixed points at comparable rho, and the corrector's Newton
  step happened to land on different ones depending on the (small) force-law difference; or
- Some residual gauge/representation freedom in the torus parametrization (the same family #544
  diagnosed a phase-pin gauge bug in, separately, already fixed on current main) makes the "mean
  state coefficients" comparison less physically meaningful than it looks, and the underlying
  physical torus is less different than these raw numbers suggest.

Distinguishing these needs more than this scoping pass affords (it would mean characterizing the
local solution landscape near SE-L2, not a re-run). What IS established, cleanly: **the SE-L2 result
quoted throughout #533/#537/#544 as the clean "positive control" baseline was computed under the
buggy code, and does not reproduce the same converged state under the corrected code.**

## Consequences for the #533-#544 task chain

- **#544's SE-L2-vs-EM-L2 contrast is likely still directionally valid but not quantitatively
  trustworthy as stated.** EM-L2's failure is dominated by residuals of 0.34-3.4 (unconverged by 4-5
  orders of magnitude) — utterly swamped by a ~0.2-scale force-law-sensitivity effect, so "EM-L2
  doesn't converge, SE-L2 does" as a qualitative statement almost certainly survives. But the SPECIFIC
  SE-L2 numbers cited (residual 3.093e-05 → 2.537e-05 after the phase-pin fix) were computed under the
  buggy alpha_6 scaling and should not be treated as the correct SE-L2 torus's true state.
- **#544's high-harmonic-content diagnosis for EM-L2** (mode-ratio ≈0.87, needing ~40 Fourier modes)
  is a statement about the orbit's own geometry (how deep it swings toward the Moon), not something
  the alpha_6 scaling bug would materially change — likely robust, but not independently re-verified
  here.
- **#533's own build-time verification ("verified against circular BCR4BP limits and finite
  differences") could not have caught this bug by construction**: a finite-difference check of the
  *same* buggy analytic formula will agree with that formula's own numerical derivative regardless of
  whether the formula encodes correct physics — FD verification catches coding-vs-derivative
  inconsistency, not a wrong physical premise. Worth noting for future QBCP verification passes:
  cross-check analytic formulas against the SOURCE PAPER's equations directly (as #592 eventually did),
  not just against their own finite-difference twin.
- **Any future QBCP work should re-derive its own baselines under the now-fixed code** rather than
  reuse #533/#537/#538/#544's specific numbers as a reference point.

## Recommendation

**Do not blindly re-run the whole #538/#544 multi-hour effort** (that's a much larger undertaking
than this scoping task, and #538/#539/#540 are not currently active work). Instead: the NEXT time
anyone picks up #538/#539/#540 or any QBCP-torus-based search, they should treat #533/#537/#544's
specific quoted numbers as **stale/unverified under the fixed code**, not as an established baseline
— re-derive the SE-L2 (or whatever) reference case fresh under current `main` before building on it.
No code changes made here; this is a scoping conclusion, not a fix.

## Update: reconciling with #544's OWN prior alpha_6 experiment (found after initial closure)

After this note was first written, closer reading of `data/OUTSTANDING.md`'s #544 entry (its own
"FOLLOW-UP-2" section) revealed something this scoping pass had missed: **#544 already tried this
exact alpha_6 fix on 2026-07-10 and explicitly reverted it**, reporting two findings:
1. It moved a QBCP "L1 substitute" periodic orbit ~30% closer to the published POL1 golden
   (2.12e-2 → 1.47e-2, distance metric).
2. It allegedly **regressed** the SE-L2 torus positive control 16× (3e-5 → 4.2e-4) — the opposite
   of what this note's Step 2 found (2.5e-5 → 1.5e-5, a slight *improvement*, not a regression).
3. It flagged an unresolved, still-open concern: whether the in-repo alpha_6 Fourier table (fitted
   by Gimeno-Jorba 2018) assumed the OLD or the Rosales-Jorba Eq. 3 scaling convention — swapping
   conventions without re-fitting could be "two wrongs," a question no numerical test alone settles.

**This is a serious discrepancy the original close-out (below) missed — the fix was applied without
searching for this prior history first, an error in the #592 verification process.** Resolving it
required two more, independent checks:

### Re-check A: is the SE-L2 test's "regression" reproducible with the actual code?

Tested whether an INCONSISTENT partial fix (patch only `qbcp_eom`, leave
`qbcp_potential_second_derivatives` buggy, or vice versa) could explain the discrepancy — e.g. a
mismatched EOM/STM Jacobian degrading Newton convergence even though the underlying physics is more
correct. Result: **`correct_qbcp_torus`/`correct_bcr4bp_torus` never call
`qbcp_potential_second_derivatives` at all** (patching it alone reproduces the fully-fixed result
exactly) — so this hypothesis is ruled out. No combination reproduces anything resembling 4.2e-4.
**The reported SE-L2 regression does not reproduce with the current, faithfully-reconstructed code.**

### Re-check B: independently rebuild the L1-substitute/POL1 test from scratch

Rather than try to replicate the 2026-07-10 session's exact, never-committed, undocumented
construction, built a genuinely independent multi-shooting periodic-orbit corrector (12 segments,
analytic STM via `qbcp_stm_eom`/`bcr4bp_stm_eom`, Newton's method) implementing the same conceptual
pipeline described in their note: (1) the exact CR3BP EM-L1 equilibrium via an independently-derived
quintic solve (cross-checked: force residual 3.3e-16, genuine fixed point; x=0.836915145 vs POL1's
published 0.836914168, delta ~1e-6 — sanity-confirms POL1 really is this L1 point's substitute);
(2) continue mu_sun 0→full within BCR4BP (which reduces exactly to CR3BP at mu_sun=0, per
`bcr4bp_stm_eom`'s own docstring) to a periodic BCR4BP orbit; (3) hand off to a QBCP multi-shooter
with the real time-varying alphas, under both the fixed and reconstructed-buggy `qbcp_eom`.

**Result — CONFIRMS the direction of #544's own finding, independently:**
```
FIXED (#592):  converged=True, resnorm=1.13e-14, distance-to-POL1 = 1.805e-02
BUGGY:         converged=True, resnorm=9.69e-14, distance-to-POL1 = 4.139e-02
#544's own numbers (different method): buggy=2.12e-2, fixed=1.47e-2
```
Both independent constructions (theirs and this one) agree on the DIRECTION and rough MAGNITUDE:
the alpha_6 fix moves the L1 substitute roughly 2-2.3× closer to the published POL1 golden. This is
a real, reproducible, doubly-confirmed improvement — not an artifact of either construction.

### Net resolution

- **POL1 agreement: genuinely improves under the fix, confirmed by two independent methods.** Solid
  evidence the fix is a real physics improvement in this respect.
- **SE-L2 "regression": NOT reproducible** with the current, exactly-verified code — the claim from
  2026-07-10 does not hold up against a careful, byte-checked re-test. Most likely explanation
  (not proven): the 2026-07-10 session's on-the-fly experiment applied the fix in some subtly
  different or inconsistent way than what ended up in the stash/#592 (that session's exact code was
  never committed, only described in prose).
- **The Gimeno-2018-fitting-convention question remains genuinely open** — no numerical test settles
  it, it needs the source paper. But with one of #544's two original objections refuted and the other
  now independently confirmed as a genuine improvement, the balance of evidence now favors #592 as a
  net-positive fix, not a reason for caution.

**Recommendation updated: keep #592 applied.** No revert warranted. The Gimeno-2018 convention
question is worth resolving if the source PDF ever becomes accessible, but is not an actionable
blocker today.

## #593 STATUS: CLOSED (fully reconciled)

Scoping complete, including reconciliation with #544's own prior (reverted) attempt at this exact
fix. Net finding: #592 is a genuine improvement (independently confirmed via the POL1 metric two
ways), the previously-reported SE-L2 regression does not reproduce, and #592 should stay applied.
Full re-verification of the #538/#544 chain's other specific numbers is still deferred to whoever
next picks that work up (per the original recommendation above) — this reconciliation only settles
whether #592 itself is trustworthy, not whether #533/#537/#544's OTHER quoted numbers need
re-deriving (they still do, since they were computed pre-fix).
