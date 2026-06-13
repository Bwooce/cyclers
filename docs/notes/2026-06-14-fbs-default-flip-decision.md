# FBS optimizer default-flip decision (#245 / #246)

**Date:** 2026-06-14
**Decision:** **DO NOT flip the default.** Keep the Lambert+FD backbone as the
default; FBS-analytic stays an **opt-in** (`gradient="fbs-analytic"`), with the
new reversible **`gradient="auto"`** (FBS-first, Lambert-fallback) available but
NOT the default. This resolves #245 (rejected on evidence) and #246 (this note is
the documentation).

## Why the flip was on the table
- **#243** (optimizer fair trial): FBS-analytic gradients won decisively *in the
  SLSQP optimizer lane* — ~4.8x faster success-to-optimum from cold seeds, same
  or better optimum, far fewer constraint evals.
- **#244** (catalogue-wide parity): on the 9 hard E-E-M-M / E-E-E-M-M rows FBS got
  strictly closer to feasibility (median 6.6x) but **neither lane converged** —
  parity-to-convergence undefined, so #244 recommended HOLD behind a seed-basin fix.
- **#248** provided that seed-basin fix and reported one row (mcconaghy) closing on
  the FBS lane at 0.0987 km/s where Lambert missed — which *appeared* to lift the
  hold.

## The deciding evidence (post-#248 parity, this task)
Apples-to-apples through the **multi-arc closure corrector** (`multiarc_closure_run`,
identical seeds, 8 multi-starts, both lanes), HEAD 52bdd34:

| row | seq | Lambert res / conv / s | FBS res / conv / s |
|---|---|---|---|
| mcconaghy-2006-em-k2 | E-E-M-M   | 0.3827 / no / 21.7 | 0.3886 / no / 61.7 |
| russell-ch4-4.991gG2 | E-E-M-M   | 0.4042 / no / 21.1 | 0.3963 / no / 57.2 |
| russell-ch4-9.353Gg2 | E-E-M-M   | 0.7135 / no / 11.9 | 0.7186 / no / 54.7 |
| russell-ch4-6.44Gg3  | E-E-E-M-M | 3.2683 / no / 25.8 | 3.2647 / no / 106.8 |

**Reading:** in the closure-corrector lane at equal multi-start budget, FBS reaches
the *same* residual as Lambert (within ~1%), converges no more often (neither
converges), and costs **3-5x more wall-clock**. The #243 advantage does NOT carry
into this lane; it is specific to the SLSQP optimizer (`optimize_chain_fbs`). The
#248 0.0987 single-row close was a deeper search that found a better basin than the
8-start default reaches here; it is not reproduced as an FBS-vs-Lambert *advantage*
at equal budget, so it does not justify a global default flip.

(Caveat, honest: the #248 0.0987 basin is not reached by the 8-start default set on
either lane here — both stall at ~0.38. Whether that deeper basin is FBS-specific
is a #248-domain question, not shown here; it does not change the equal-budget
parity conclusion.)

## What stands
- Default: **Lambert+FD** (unchanged, byte-identical for every existing caller).
- Opt-in: `gradient="fbs-analytic"` (the #226/#244 lane) for the SLSQP optimizer
  where #243's win lives, and `gradient="auto"` (#245, committed 52bdd34) — FBS-first
  with Lambert fallback, never converges worse than Lambert — for callers that want
  the optimizer-lane robustness with a safety net.
- A clean negative is a success: we did not flip a system default on a single-row
  result that the equal-budget parity does not corroborate.
