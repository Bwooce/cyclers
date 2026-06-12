# Ross 5-row long-span discriminating run — V2-ballistic evidence package (#229)

**Date:** 2026-06-13
**Rows:** `ross-rt-em-cycler-{11,21,31,32,33}-2025` (all V1 since #216)
**Spec for this task:** `docs/notes/2026-06-13-ross-v2-assessment.md` (#224)
**Instrument:** `scripts/ross_v2_longspan.py` (method tag `ross-v2-longspan-v1`),
reusing the campaign harness `cr3bp-inertial-rebound-ias15-v1`
(`scripts/cr3bp_family_search.py::inertial_crosscheck`, frame conversion +
Jacobi bookkeeping unchanged) at a 100-period span.
**Pinning test:** `tests/search/test_cr3bp_v2_longspan.py` (slow-marked,
(1,1) row, 40-period still-discriminating span).
**Status:** catalogue writeback **HELD FOR USER REVIEW** — `data/catalogue.yaml`
is untouched by this task. The "Proposed writeback" section below is a
proposal only.

## Model scope (like-for-like — read first)

**CR3BP is the defining model of all five rows**, and every number here is
evaluated in that model: the inertial REBOUND/IAS15 harness integrates the
CR3BP idealisation in inertial coordinates (the two primaries on the exact
circular two-body rail, the spacecraft a massless test particle), then
back-transforms to the rotating frame. This is a *different code path AND
frame* than the rotating-frame DOP853 corrector — the false-consensus
independence the spec §14 V2-ballistic clause wants — but it is **NOT a
real-ephemeris claim**. No bicircular, solar-gravity, lunar-eccentricity or
ephemeris content exists for these rows (the paper's stated future work).
**V2-ballistic is this lane's ceiling** until an ephemeris-model continuation
exists; no V3-class language applies anywhere in this package.

## Why this run, and why the span has teeth (per-row justification)

The campaign's 5-period inertial gate is weakly discriminating for
residual-zero ICs: with seed error δ₁ ~ 1e-9..1e-8 and amplitude A ~ O(1) nd,
even |λ| ~ 10²..10³ cannot amplify δ₁ past the 3A departure band in 5 periods.
The discriminating requirement (assessment note) is

    T_span ≥ T · ln(3A/δ₁) / ln(|λ_hypo|),  with conservative |λ_hypo| = 2,

evaluated with the **actual measured δ₁** per row. The measured δ₁ (table
below) gives N_req = ln(3A/δ₁)/ln 2 ≈ **29.6–32.7 periods** across the five
rows; the run used **100 periods — 3× the requirement for every row**. A
hypothetical |λ| = 2 instability would have amplified the measured seed error
past the 3A band before period ~33 on every row; nothing remotely like that is
observed. In wall-clock terms the spans are 12.3–24.9 years of Earth-Moon
dynamics (100 × T × 4.348377 d/TU).

## Integrator hygiene (recorded settings)

- REBOUND 5.0.0, IAS15, `epsilon = 1e-9` (the campaign's setting, unchanged);
  G = 1 with GM-pair set so the circular rail's mean motion matches 1/t_s
  exactly (self-consistent rail); 128 samples/period (mid-period departures
  cannot hide between recurrence samples).
- Noise floor monitored per row (max of primary-rail deviation and
  centre-of-mass drift): 1.2e-11–3.1e-11 nd over the full 100-period span —
  2+ orders below every measured δ₁, so the recurrence measurements are real
  signal, not integrator noise.
- Member preparation: fixed-Jacobi symmetric corrector
  (`correct_symmetric_fixed_jacobi`, tol 1e-10, DOP853 rtol=atol=1e-12),
  Barden ν from the half-period STM at rtol=atol=1e-12, amplitude A from the
  128-sample one-period sweep (`_gate_metrics`, rtol=atol=1e-12).

## Per-row results (run of 2026-06-13, deterministic; ~8 s total)

| row | family | C (sourced) | T (nd) | ν (Barden) | A (nd) | δ₁ (nd) | N_req | span | max drift 1st half | max drift 2nd half | ratio | Jacobi drift (span) | verdict |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| `ross-rt-em-cycler-11-2025` | (1,1) | 3.151175879508174 | 10.29206923 | −0.00334 | 1.8149 | 4.923e-09 | 30.0 | 100 | 4.687e-07 | 9.244e-07 | 1.97 | 1.040e-10 | **BOUNDED** |
| `ross-rt-em-cycler-21-2025` | (2,1) | 3.129389531088256 | 19.44016043 | +0.05007 | 1.4838 | 5.619e-09 | 29.6 | 100 | 4.698e-07 | 1.103e-06 | 2.35 | 8.560e-11 | **BOUNDED** |
| `ross-rt-em-cycler-31-2025` | (3,1) | 3.161784147013429 | 14.78826794 | +0.01545 | 1.2743 | 1.648e-09 | 31.1 | 100 | 9.695e-07 | 1.985e-06 | 2.05 | 5.593e-11 | **BOUNDED** |
| `ross-rt-em-cycler-32-2025` | (3,2) | 3.182762663084288 | 17.90058012 | −0.01175 | 1.3169 | 5.480e-10 | 32.7 | 100 | 3.994e-06 | 7.744e-06 | 1.94 | 2.541e-10 | **BOUNDED** |
| `ross-rt-em-cycler-33-2025` | (3,3) | 3.177224018696528 | 18.14546057 | +0.06001 | 1.3367 | 3.276e-09 | 30.2 | 100 | 5.347e-07 | 1.066e-06 | 1.99 | 2.704e-10 | **BOUNDED** |

Column notes: δ₁ = measured full-state rotating-frame recurrence after the
first period (the seed error the discrimination formula uses); drift columns =
max per-period recurrence |state(kT) − state(0)| sampled at every integer
period over each half of the span; ratio = 2nd-half max / 1st-half max;
Jacobi drift = max |C(t) − C(0)| over all 12,800 samples. The departure band
3A is 3.82–5.44 nd per row — the observed max drift sits **5–7 orders of
magnitude inside it** on every row.

**Honest reading of the drift shape (bounded band, not perfect closure):** the
per-period recurrence grows ~linearly with k (half-ratio ≈ 2 ≙ a 2× growth
between span halves), the signature of slow along-track *phase* drift from the
finite seed error / T-resolution — exactly what a linearly stable orbit with a
slightly-off IC does. Linear extrapolation of the worst row ((3,2)) reaches
the 3A band only after ~5e5 periods. The exponential alternative is excluded
with enormous margin: |λ| = 2 over the second half alone predicts a half-ratio
of ~2⁵⁰ ≈ 1e15, vs the observed ≤ 2.35. The bounded-band verdict criterion is:
no divergence, every sample inside 3A, AND half-ratio ≤ 10 (separating linear
phase drift from exponential amplification by 14+ orders; band never
loosened).

## Evidence package (the four elements, per the assessment note)

1. **The existing 5-period R1/R2/R3 PASS** (campaign artifact): every row's
   published representative re-propagated in the inertial harness, R1
   one-period recurrence ≤ 0.1A, R2 Jacobi drift ≤ 1e-9, R3 no early
   departure. `docs/notes/2026-06-12-cr3bp-continuation-results.md`
   (REPRODUCTION rows, commits f99ba30/c4d6c57, method
   `cr3bp-jacobi-continuation-v1` + gate `cr3bp-inertial-rebound-ias15-v1`).
2. **The new long-span bounded-drift run (the discriminating instrument)**:
   this note's table — 100 periods (≥3× the per-row N_req from the measured
   δ₁), all five rows BOUNDED, Jacobi drift ≤ 2.7e-10 over the full span,
   noise floor 2+ orders under signal. Re-runnable:
   `uv run python scripts/ross_v2_longspan.py`; pinned by
   `tests/search/test_cr3bp_v2_longspan.py` (slow).
3. **The ν-in-published-window reproduction**: Barden ν per row (table above)
   reproduces inside Ross's published stable windows, |ν| < 1 and ν ~ 0
   midpoint (Table 3 members), as already gated by
   `tests/search/test_cr3bp_ross_families.py::test_ross_family_reproduced`.
4. **The honest-limit distinction, preserved verbatim** from the campaign
   results note (`docs/notes/2026-06-12-cr3bp-continuation-results.md`):

   > **Honest limit of this gate (do not overread a PASS).** The corrector
   > lands every member to a ~1e-10 perpendicular-crossing residual, so the
   > inertial one-period recurrence delta1 is ~1e-9..1e-8 -- and even a
   > strongly unstable member (|nu|~10^2, |lambda|~10^2..10^3) does not
   > amplify that tiny seed past 3A within the 5-period span, so it stays
   > numerically bounded and R3 records departure LATER than predicted (a
   > PASS). The inertial gate therefore confirms each member is a GENUINE
   > periodic orbit through an independent integrator+frame (it re-closes and
   > conserves Jacobi) -- it is the false-consensus *consistency* gate. It
   > does NOT independently certify STABILITY: stability is the Barden nu
   > (reported per member), which IS discriminating here (nu spans -0.0
   > stable to +360 wildly unstable across these branches). Read 'inertial
   > PASS' as 'a real orbit, cross-checked', and the STABLE verdict as 'from
   > nu', never conflated.

   Accordingly: **the STABLE verdict for every row derives from Barden ν plus
   the long-span discriminating run (elements 2+3), never from the 5-period
   gate (element 1)** — element 1 certifies "a real orbit, cross-checked",
   nothing more.

## PROPOSED WRITEBACK (HELD FOR REVIEW)

Proposal only — nothing below has been applied. All five rows passed the
discriminating run, so no row is excluded. Per row
(`ross-rt-em-cycler-{11,21,31,32,33}-2025`), the writeback would change:

1. **`data/catalogue.yaml`** — `validation_level: V1` → `validation_level:
   V2-ballistic`, with the inline comment updated to cite spec §14
   V2-ballistic (≥3 continuous laps, bounded rotating-frame drift, evaluated
   in the row's **defining model = CR3BP** — the like-for-like scope recorded
   the same way the V1 scoping convention does) and to point at this note +
   `tests/search/test_cr3bp_v2_longspan.py`.
2. **`src/cyclerfinder/data/validate.py` `_LEVEL_EVIDENCE`** — re-key each
   row's entry from `(id, "V1")` to `(id, "V2-ballistic")` and extend the
   evidence string with: the 100-period inertial REBOUND/IAS15 bounded-band
   result (max drift ≤ 7.8e-6 nd vs band 3A ≥ 3.8 nd; Jacobi drift ≤ 2.7e-10
   over the span; span ≥ 3× the measured-δ₁ discrimination threshold
   N_req ≈ 30–33), the honest-limit clause (STABLE from Barden ν + long-span
   run, never the 5-period gate), and the CR3BP like-for-like scope sentence
   (no real-ephemeris claim; V2-ballistic is the lane ceiling).
3. **`scripts/backfill_validation_level.py` `_LEVEL_BY_ID`** — the five ids
   move from `"V1"` to `"V2-ballistic"` (kept consistent with 1+2).
4. **No other fields change.** No V3-class language anywhere; `delta_v_kms`
   stays `null` (ballistic); `source_ephemeris` stays `null` (pure CR3BP).

Consistency check at writeback time (assessment note disposition): the V2
evidence fields must record the like-for-like model scope (CR3BP) the same
way the V1 scoping convention does (spec §14 V2-ballistic clause).
