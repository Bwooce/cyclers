# Task #505: Pluto-Charon (3,2) V2-ballistic long-span evidence

**Date:** 2026-07-01
**Row:** `ross-rt-pc-cycler-32-2026` (V1 since #494)
**Instrument:** `scripts/pc_v2_longspan.py` (method tag `pc-v2-longspan-v1`),
reusing the campaign harness `cr3bp-inertial-rebound-ias15-v1`
(`scripts/cr3bp_family_search.py::inertial_crosscheck`) at a 100-period span.
**Evidence test:** `tests/search/test_505_pc_v2_longspan.py` (default suite,
NOT slow-marked; 50-period discriminating span, ~4 s).
**Status:** catalogue writeback **HELD FOR USER ADJUDICATION** — `data/catalogue.yaml`
is untouched by this task.

---

## Model scope (like-for-like — read first)

**CR3BP is the defining model of this row**, and every number here is evaluated
in that model: the inertial REBOUND/IAS15 harness integrates the CR3BP
idealisation in inertial coordinates (the two primaries on the exact circular
two-body rail, the spacecraft a massless test particle), then back-transforms
to the rotating frame.  This is a *different code path AND frame* than the
rotating-frame DOP853 corrector — the false-consensus independence the spec
§14 V2-ballistic clause requires — but it is **NOT a real-ephemeris claim**.
**V2-ballistic is this lane's ceiling** (no bicircular, heliocentric-
perturbation, or ephemeris content applies to this CR3BP row).

---

## Positive-control confirmation

The corrector reproduces the admitted IC:

| Quantity | Catalogue (DERIVED, #494) | This run | Match? |
|---|---|---|---|
| μ | 0.10876473603280369 | same (cr3bp_system) | — |
| C | 3.57951501972907 | 3.57951501972907 | exact |
| x₀ | −0.693198287043369 | −0.693198287043 | <1e-10 |
| ẏ₀ | −0.297004785528322 | −0.297004785528 | <1e-10 |
| T (nd) | 11.8334625170346 | 11.833462517068 | <1e-10 |
| Barden ν | 3.82e-9 (catalogue) | 1.213e-06 | ~0 (maximally stable) |

Closure residual: 9.39e-12 (well within the 1e-6 cross-check tolerance).
Topology confirmed (3,2) prograde by `winding_topology` (repeated in #504).
Independent Radau cross-check: PASS (dJ < 1e-12, per #494 + #504).

---

## Why the span has teeth (discrimination criterion)

The 5-period inertial gate from the campaign is weakly discriminating for
residual-zero ICs: with seed error δ₁ ~ 3.5e-10 nd and amplitude A ~ 1.64 nd,
even |λ| ~ 10⁴ cannot amplify δ₁ past the 3A departure band in 5 periods.
The discriminating requirement (same as #229 EM evidence):

    N_span >= ln(3A / delta_1) / ln(|lambda_hypo|),  |lambda_hypo| = 2

With the measured δ₁ = 3.536e-10 nd and A = 1.639645 nd:

    N_req = ln(4.919 / 3.536e-10) / ln(2) = 33.7 periods

The run used **100 periods — ~3x the per-row requirement**.  A hypothetical
|λ| = 2 instability would have amplified the measured δ₁ past the 3A band
before period ~34; nothing remotely like that is observed.

---

## Long-span result (100-period run of 2026-07-01, ~1.4 s)

| row | C (DERIVED) | T (nd) | ν (Barden) | A (nd) | δ₁ (nd) | N_req | span | max drift 1st half | max drift 2nd half | ratio | Jacobi drift (span) | verdict |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| `ross-rt-pc-cycler-32-2026` | 3.57951501972907 | 11.83346252 | 1.21e-6 | 1.6396 | 3.536e-10 | 33.7 | 100 | 3.475e-09 | 6.619e-09 | 1.90 | 4.012e-10 | **BOUNDED** |

**Harness settings:** REBOUND 5.0.0, IAS15, `epsilon = 1e-9`; G=1 with
GM-pair set so the circular rail's mean motion matches 1/t_s exactly
(self-consistent rail); 128 samples/period.

**Key numbers:**

- δ₁ = 3.536e-10 nd (measured seed error after 1 period)
- Departure band 3A = 4.919 nd
- **max drift over 100 periods: 6.619e-09 nd — 8.9 orders of magnitude inside the band**
- Half-ratio 1.90: consistent with slow linear phase drift (a real |λ|=2 over
  the second half alone predicts a ratio of ~2⁵⁰ ≈ 10¹⁵; observed ≤ 1.90)
- Jacobi drift: 4.012e-10 (< 1e-9 IAS15-hygiene bound)
- Noise floor: 1.068e-11 nd (2 orders below δ₁; measurements are real signal)
- Diverged: False; all 12,800 samples finite

**Reading the drift shape:** the per-period recurrence grows slowly and
linearly (ratio 1.90 ~ 2), the signature of along-track *phase drift* from the
finite seed error — exactly what a linearly stable orbit with a slightly-off
period IC does.  Linear extrapolation of this worst-case rate reaches the 3A
band only after ~7e8 periods (~8 × 10⁸ × 12.03 d ~ 3 × 10¹⁰ years).
Exponential divergence is excluded with 14+ orders of margin.

---

## Evidence package (the four elements, per spec §14 / #229 pattern)

1. **Positive-control confirmation**: corrector closes on the catalogue IC
   (C, x0, T all within 1e-10), Barden ν = 1.21e-6 (~0, maximally stable),
   (3,2) prograde topology confirmed, independent Radau PASS. This task +
   #494 + #504.

2. **Long-span bounded-drift run (the discriminating instrument)**: the 100-period
   table above — span ≥ 3× the per-row N_req from the measured δ₁, BOUNDED,
   Jacobi drift ≤ 4.0e-10 (< 1e-9 bound), noise floor 2 orders under signal,
   drift 8.9 orders inside the band. Re-runnable:
   `uv run python scripts/pc_v2_longspan.py`; pinned by
   `tests/search/test_505_pc_v2_longspan.py` (default suite, NOT slow).

3. **ν-in-published-window**: Barden ν = 1.21e-6, well inside |ν| < 1 (stable);
   the nu=0 midpoint construction (#494) guarantees this is the maximally-stable
   member of the stable window. Confirmed again by `winding_topology` in #504.

4. **Honest-limit distinction (verbatim from #229 pattern)**:

   > The 5-period inertial gate confirms each member is a GENUINE periodic orbit
   > through an independent integrator+frame — it is the false-consensus
   > *consistency* gate. It does NOT independently certify STABILITY: stability
   > derives from Barden ν (reported per member) PLUS the long-span discriminating
   > run.  Read "inertial PASS" as "a real orbit, cross-checked"; read
   > "STABLE" as "from ν + long-span run", never conflated.

---

## V1 → V2-ballistic recommendation

**RECOMMENDATION: PROMOTE `ross-rt-pc-cycler-32-2026` from V1 to V2-ballistic.**

Justification:
- All spec §14 V2-ballistic criteria are satisfied:
  ≥3 continuous laps (100 periods used, no departure), bounded rotating-frame
  drift (max 6.62e-9 nd vs band 4.92 nd), evaluated in the row's defining
  model (CR3BP), independent integrator (REBOUND IAS15, different code path +
  frame from the DOP853 corrector).
- The discriminating criterion N_span ≥ N_req is met with 3× margin.
- The orbit is maximally stable (ν ~ 0), the nu=0 midpoint of the only stable
  (3,2) window found at Pluto-Charon mu — this is the strongest-possible
  stability claim within the CR3BP model.
- The evidence test `test_505_pc_32_longspan_bounded_band` runs in the DEFAULT
  CI suite (no `@pytest.mark.slow`), satisfying the project rule that V-claim
  evidence tests must be CI-gated.

**Scope limitation (unchanged from V1):** V2-ballistic refers to the CR3BP
model only — the same model scope as the V1 claim.  No real-ephemeris
continuation exists for this row (the paper's stated future work).
V2-ballistic is the ceiling of this lane; no V3-class language applies.

**Decision is the user's.** This note + test constitute the evidence package;
`data/catalogue.yaml` is untouched.

---

## Proposed writeback (HELD for adjudication)

If promoted, the following changes would be made:

1. `data/catalogue.yaml` — `validation_level: V1` → `validation_level: V2-ballistic`,
   with inline comment updated to cite spec §14 V2-ballistic (100-period
   REBOUND/IAS15 inertial bounded-band run, max drift 6.62e-9 nd vs band 4.92 nd,
   Jacobi drift 4.0e-10, span 3× N_req=33.7, CR3BP like-for-like scope) and
   pointing at this note + `tests/search/test_505_pc_v2_longspan.py`.

2. No other fields change. `delta_v_kms` stays null (ballistic);
   `source_ephemeris` stays null (pure CR3BP).
