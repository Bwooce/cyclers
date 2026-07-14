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

**#593 STATUS: CLOSED** (scoping complete; genuine, non-trivial impact confirmed for at least one
past result; full re-verification of the #538/#544 chain deferred to whoever next picks that work up).
