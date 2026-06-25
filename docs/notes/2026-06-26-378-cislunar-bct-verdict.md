# #378 cislunar BCT substrate — verdict (capability-only; clean negative on the quasi-cycler)

**Date:** 2026-06-26
**Issue:** #378 (cislunar BCT integration substrate)

## What was built (CAPABILITY — delivered)

A complete, sourced, reusable Belbruno weak-stability-boundary (WSB) /
ballistic-capture-transfer (BCT) substrate on the project's incoherent BCR4BP:

| Module | Content | Goldens |
|---|---|---|
| `core/wsb.py` | E_2 (eq 3.6), periapsis sigma (eq 3.9), analytic-W (eq 3.29), C_1 validity (Def 3.22), numerical stability-class (§3.2.1) | parabolic C=±√2 (Lemma 3.34), C_1≈3.184 from L_1 Jacobi, E_2(L)<0 sign (Lemma 3.30) |
| `genome/bct_transfer.py` | BCTTarget/QF-on-W, `construct_bct_backward` (arc II → Sun apoapsis), `correct_bct_forward` (2×2 targeting), `build_hiten_bct`, BCT novelty self-test | Hiten signature band (apoapsis 3.9 LD ±30%, ballistic E_2≤0, ΔV_capture=0) |
| `search/cislunar_bct_search.py` | grid sweep + transfer-vs-quasi-cycler classifier + novelty-emittable gate | — |

All tests in the default (non-slow) suite; ruff + mypy clean.

## Honest assessment: capability vs catalogue object

**This produced a CAPABILITY, not a catalogue-class transport object.** Per the
design's critical-honesty mandate (and Belbruno himself: "WSB orbits are not
cyclers"):

* **Phase 0.2 gate PASS** — the incoherent BCR4BP shapes a Sun-driven Earth
  apoapsis to 4.6–5.6 LD from LEO (Hiten 3.9 LD reachable). R1 did not fire.
* **Hiten apoapsis signature REPRODUCED** — the backward constructor reaches the
  Hiten band; θ₂=1.25 lands at 3.95 LD, a bullseye on Belbruno's published
  3.9 LD, with exact on-W ballistic capture (E_2 = −0.064 < 0, ΔV_capture = 0).
  4 of 6 in-family configs fall in [2.7, 5.1] LD.
* **Clean NEGATIVE on the quasi-cycler** — across the θ₂/e₂/branch sweep, **no
  BCT's return leg re-acquires W**: every constructed transfer classifies as
  `transfer` (precursor_mga capability), **zero `quasi_cycler_candidate`**. This
  is the expected default (design R2/§6) and is consistent with Belbruno
  Theorem 3.58 (capture on W is a *chaotic* process — cuts against clean
  periodicity). No catalogue row. NOTHING self-admitted.

## Two documented model-gap boundaries (the empty-region re-open keys)

1. **Forward-from-LEO on-W capture did not converge** in the incoherent BCR4BP:
   a single forward arc reaches the Moon's vicinity (r_23 → ~2000 km) but
   arrives hyperbolically (E_2 > 0), not ballistically captured. The on-W
   capture is delivered by the *backward* QF (exact, on W by construction); the
   forward stitch is the convergence-engineering piece in the un-acquired
   [39] Belbruno-Miller 1993 and likely needs the coherent model.
2. **Coherent QBCP** (Andreu α-tables, un-digested) is the higher-fidelity model
   a future re-sweep would use.

Both are carried as `method_capability` tags on the empty-region record so a
later, more-capable method knows to re-open this region.

## Disposition

* Capability ships (3 modules + goldens), unconditional.
* Clean negative logged to `data/empty_regions.jsonl`
  (`cislunar-bct-wsb-quasicycler-2026-06-26`) with WSB/BCT re-open keys.
* No catalogue edit; no novelty claim (the Hiten signature is correctly flagged
  non-novel against the Belbruno/Hiten corpus by the BCT self-test).
* No human-gauntlet flag raised: no transport-relevant quasi-cycler candidate
  was found.
