# Ross 5-row V1 -> V2-ballistic assessment (#224)

**Date:** 2026-06-13
**Rows:** `ross-rt-em-cycler-{11,21,31,32,33}-2025` (all V1 since #216)
**Question:** does evidence already produced by the #218/#219 campaign satisfy the
spec §14 **V2-ballistic** gate ("≥3 continuous laps; bounded drift in the dynamic
rotating frame ... evaluated in the row's defining model")?

## What the existing evidence shows

Each of the 5 published representatives was, during the campaign
(`cr3bp-jacobi-continuation-v1`, commits f99ba30/c4d6c57):

1. Re-derived by the rotating-frame DOP853 fixed-Jacobi corrector to the printed
   15-digit C/T precision (the V1-class reproduction, #212).
2. Re-propagated through the **independent inertial REBOUND/IAS15 harness**
   (different code path AND frame; Jacobi-conserving, i.e. the same CR3BP
   defining model in inertial coordinates) over a **5-period span**, passing
   R1 (one-period recurrence ≤ 0.1A), R2 (Jacobi drift ≤ 1e-9 over the span),
   R3 (no departure earlier than the linear prediction).
3. Stability ν reproduced inside Ross's published stable windows (all 5 rows
   |ν| < 1).

Mechanically, item 2 reads as the V2-ballistic gate: ≥3 continuous laps,
bounded drift, defining model, independent integrator.

## Why mechanical satisfaction is NOT enough (the honest-limit problem)

The campaign results note records the gate's honest limit: the corrector lands
members at ~1e-10 crossing residual, so the inertial seed error δ₁ is
~1e-9..1e-8 — and over only 5 periods **even a strongly unstable member
(|λ| ~ 10²..10³) does not amplify that seed past the 3A bound**. The 5-period
instrument is therefore weakly discriminating for residual-zero ICs: it would
grade an unstable orbit "bounded" too. Promoting on that evidence alone would
be exactly the false-confidence pattern the orbit-closure discipline exists to
block ("it stayed bounded!" is only meaningful if it COULD have failed).

Contrast the V2 precedent: Aldrin's V2-powered came from a multi-lap
propagation in which divergence was genuinely possible and maintenance was
genuinely accounted. V2 must mean the same evidentiary strength across rows.

## Recommendation: promote AFTER one cheap discriminating run (per row)

Extend the inertial IAS15 propagation per row to a span where the instrument
has teeth: choose T_span such that a hypothetical unstable multiplier would
amplify the actual seed error past the band, i.e.
`T_span ≥ T · ln(3A/δ₁)/ln(|λ_hypo|)` with a conservative |λ_hypo| (e.g. the
smallest |λ| > 1 the corrector can resolve, or simply 30–100 periods, which for
δ₁ ~ 1e-9 and A ~ 1 nd makes even |λ| ~ 2 fail visibly). For a genuinely
stable member (all 5 rows: |ν| < 1, Ross's named stable windows) drift must
stay in a bounded oscillation band over the full span — that claim is
non-vacuous, matches Ross's own stability-window assertion (a same-model
cross-check of the published claim, not just of our corrector), and runs in
seconds per row.

**Evidence package for the writeback (review-gated, per standing rule):**
- the existing 5-period R1/R2/R3 PASS (campaign artifact),
- the new long-span bounded-drift run (the discriminating instrument),
- the ν-in-published-window reproduction,
- an explicit evidence note that the STABLE verdict derives from Barden ν and
  the long-span run, never from the 5-period gate (the campaign's honest-limit
  distinction, preserved verbatim).

**Do not** promote on the campaign artifacts alone; **do not** claim V3-class
language anywhere (no real-ephemeris content exists for these rows — CR3BP is
the defining model and V2-ballistic is this lane's ceiling until an
ephemeris-model continuation exists).

## Disposition

- Follow-up execution task: run the long-span discriminating propagation for
  the 5 rows + assemble the evidence package + HOLD writeback for user review.
- Consistency check at writeback time: V2 evidence fields must record the
  like-for-like model scope (CR3BP) the same way the V1 scoping convention
  does (spec §14 V2-ballistic clause).
