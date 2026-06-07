# Aldrin promotion assessment from the #158 continuation evidence (task #160)

**Type:** assessment / recommendation only. No catalogue writeback, no code change,
no test run beyond cheap inspection. The user decides; this note recommends.

**Evidence read:**
- `docs/notes/2026-06-07-continuation-driver-results.md` (the #158 result)
- `tests/search/test_continuation.py` (the encoded golden gate)
- `docs/notes/2026-06-07-russell-2004-continuation-deepdive.md` (what Russell's
  model actually is)
- `docs/spec.md` §14 (the V-ladder, as amended 2026-06-07) + §16.7.12 (the
  validation_level history)
- `src/cyclerfinder/verify/{gauntlet,propagate,v2_powered}.py` (the gauntlet
  machinery)
- `src/cyclerfinder/search/continuation.py` (what the driver actually emits)
- `src/cyclerfinder/data/validate.py` `_LEVEL_EVIDENCE` (what a writeback must cite)
- `data/catalogue.yaml` rows `aldrin-classic-em-k1-{outbound,inbound}`,
  `russell-ch4-4.991gG2`, and the `russell-ocampo-*` family

---

## 0. The bottom line

**Hold. Do not promote any row off the continuation evidence as it stands.** The
#158 result is a real, sourced, end-to-end validation of the *continuation
driver*, and it reproduces Russell Table 5.5's 0-m/s Aldrin entry in Russell's own
model. But it does **not** satisfy the spec's words for any unearned ladder gate:
it is a corrector closure residual, not the spec's multi-lap bounded
rotating-frame drift (V2-ballistic) and not an n-body / phase-matched horizon-TCM
realisation (V3). The honest counter-arguments (V∞ 1.5 km/s off anchors, t0 26 d
off the window, off-Aldrin elements) further argue against attaching it to the
Aldrin catalogue row without a like-for-like model match. Details below.

---

## 1. Which row's evidence is this? (the row-attribution question)

**The continuation reproduced Russell's Table 5.5 member `6.399G1` (#1), in
Russell's Chapter-5 model and his Aug-6-2003 window.** The deep-dive establishes
that `6.399G1 (#1)` *is* the Aldrin cycler:

> "**6.399G1 (#1) — the Aldrin cycler.** total Δv = **0 m/s**, launch **Aug-2003**
> (specifically Aug 6, 2003, Fig.5.7 p.181)" — `…russell-2004-continuation-deepdive.md` §3

So at the level of *physical cycler identity*, this is the Aldrin family, and the
only Aldrin catalogue rows are `aldrin-classic-em-k1-{outbound,inbound}`. There is
**no `6.399G1` row** in the catalogue (grep: the string appears nowhere in
`data/catalogue.yaml`), and no `russell-ocampo-*` row carries that descriptor —
the `russell-ocampo-*` rows use Russell's Table-3.4/4.x descriptors (`2.1.1+2`,
`3.5.1+1`, …), and the Aldrin row itself is cross-referenced to Table 3.4 cycler
**`1.0.1.-1`**, not the Chapter-5 `6.399G1`. So the candidate target is
`aldrin-classic-em-k1-outbound`.

**But "same cycler" is not "same defining model", and the spec's V2/V3 checks are
model-relative.** The spec is explicit that the like-for-like model scope is what
the gate is evaluated against:

> "evaluated **in the row's defining model** (for a circular-coplanar row that is
> the idealized propagation; the like-for-like model scope is recorded in the
> evidence, same convention as the V1 scoping)" — spec §14, V2-ballistic row

The `aldrin-classic-em-k1-outbound` row's defining model is **circular-coplanar**:

```
model_assumption: circular-coplanar      # catalogue.yaml line 27
orbit_source: rogers-2012-t1             # a=1.60, e=0.393
orbit_fidelity: circular-coplanar
```

The continuation evidence is the *opposite* model: a solution walked all the way
**to the true ephemeris (DE440)**, ending on Russell's Chapter-5 force model
(patched-conic two-body legs between true-ephemeris planet positions). That is a
*higher-fidelity realisation* of the Aldrin family, not a like-for-like
circular-coplanar reproduction of the row's defining model.

**Attribution verdict:**
- The continuation result is **evidence about the Aldrin family** and is most
  naturally attached to `aldrin-classic-em-k1-outbound` (the outbound Aldrin row),
  because that is the only catalogue home for the Aldrin cycler and the result is
  the up-escalator one-synodic free-return Russell labels #1.
- It is **corroboration**, not a like-for-like V-gate, for that row, because its
  model (true ephemeris / Russell Ch.5) is not the row's defining model
  (circular-coplanar). The spec's model-relative phrasing means the same physics in
  a *different* model does not automatically satisfy a gate scoped to the row's
  model.
- It is **not** evidence for any `russell-ocampo-*` row (different descriptors,
  different cyclers) and **not** for `russell-ch4-4.991gG2` (that is the S1L1
  member `4.991gG2 (#83)`, which the same #158 run reports as a *strict-xfail break
  point*, not a closure).

---

## 2. Does the 7-cycle ballistic closure satisfy **V2-ballistic**?

**No — not as the evidence stands, and not without RUNNING the V2-ballistic check.
Residual-bounded corrector closure is not the spec's "bounded drift in the
dynamic rotating frame".**

The spec's V2-ballistic words are:

> "**V2-ballistic** | Multi-lap periodicity (ballistic) | ≥3 continuous laps;
> **bounded** drift in the dynamic rotating frame (tolerant of geometric
> breathing) … | propagate.py" — spec §14

Two distinct things are being conflated and must be separated:

**(a) What the continuation emits.** `continuation.py` returns a
`ContinuationResult.best_final` carrying `max_residual_kms`, `t0_sec`,
`vinf_kms`, `converged`. The "7 cycles" is Russell's *seed construction* — the
deep-dive: "INPUT seven cycles of the simple-model parent beginning at the epoch —
that seven-cycle simple-model trajectory is the `5n` initial guess for the
homotopy" (§1). The reported `0.00158 km/s` is the **final-step corrector
residual** at the true-ephemeris model — a V∞/closure root-find residual, *not* a
propagated lap-to-lap geometric drift. The driver does **not** propagate a
continuous trajectory across laps and does **not** compute a rotating-frame drift.

**(b) What V2-ballistic measures.** `verify/propagate.py` propagates a `Cycler`
"continuously across `n_laps` consecutive laps … and measure[s] the lap-to-lap
drift in the spec §12(c) dynamic rotating frame" with a binding km tolerance
(`DRIFT_TOLERANCE_KM`). It distinguishes *secular accumulation* from *bounded
periodic breathing* over ≥3 laps. This is a different instrument from a corrector
residual.

> "≥3 laps = two intervals — the minimum that distinguishes *secular
> accumulation* (drift that grows lap-over-lap) from *bounded periodic breathing*
> (drift that oscillates within a band)." — spec §14 note

A 0.00158 km/s residual says the corrector *found* a closing free-return at one
true-ephemeris epoch. It says nothing about whether propagating that solution
forward for ≥3 laps stays within a bounded drift band — which is exactly the
V2-ballistic question. They are not the same claim.

**There is a further, structural obstacle specific to Aldrin.** The Aldrin row's
*defining model is circular-coplanar*, and in that model Aldrin is **powered, not
ballistic**:

> `trajectory_regime: powered` … "the 1L1 EARTH (geocentric) flyby cannot supply
> the required ~84 deg turn (max ~72 deg at a 200 km Earth flyby)" — catalogue.yaml
> line 26

So the *like-for-like circular-coplanar* V2 gate for this row is **V2-powered**
(which it already holds), not V2-ballistic. V2-ballistic would have to be evaluated
in a *different* model (the true-ephemeris one where the continuation found a
ballistic closure) — and the spec scopes the gate to the row's defining model. You
cannot earn V2-ballistic for a row whose defining model says the cycler is powered
by pointing at a ballistic closure in a different model; that would be the exact
model-scope violation the §14 amendment was written to prevent.

**What running V2-ballistic would take (if one wanted to evaluate it in the
true-ephemeris model as a *new* scoped claim).** It is **not cheap and I did not
run it**: the continuation emits a corrector solution `(t0, a, e, vinf)`, not a
`Cycler` object, so one would have to (i) build a `Cycler` from the continued
true-ephemeris solution, (ii) propagate it ≥3 laps on the DE440 ephemeris via
`propagate.py`, (iii) measure dynamic-rotating-frame lap-to-lap drift against
`DRIFT_TOLERANCE_KM`, and (iv) decide and record the *model scope* of the new
claim (true-ephemeris, NOT the row's circular-coplanar defining model). Each lap is
a full astropy/DE440 propagation; this is the multi-lap ephemeris machinery, in the
minutes-to-tens-of-minutes class, and — more importantly — it would be evidence
for a *different-model* claim than the row carries. Not run here; it is not a
<10-min slot-in, and the model-scope problem makes it the wrong instrument for a
circular-coplanar row anyway.

**V2-ballistic verdict: NO.** (a) The evidence is a corrector residual, not the
spec's bounded rotating-frame drift; the check has not been run. (b) Even if run,
the Aldrin row's defining model is circular-coplanar/powered, so a ballistic
true-ephemeris closure is the wrong-model instrument for *this* row's V2-ballistic.

---

## 3. Does it bear on **V3**?

**It is genuine, sourced V3-*relevant* progress, but it does not satisfy the
spec's V3 words, and it is not n-body.** V3 reads:

> "**V3** | Ephemeris realisation | phase-matched to a real launch window;
> ephemeris-mode horizon TCM over 3–5 laps (~20–30 yr) bounded and within ΔV
> budget | astropy backend" — spec §14

Point-by-point against the evidence:

- **"phase-matched to a real launch window"** — *partial.* The continuation lands
  at t0 = 2003-07-11, **26 d** from Russell's sourced Aug-6-2003 window (well
  inside the test's 200-day gate). That is genuine phase-matching to a real window
  — the strongest V3-shaped part of the result.
- **"ephemeris-mode horizon TCM over 3–5 laps (~20–30 yr) bounded and within ΔV
  budget"** — **not demonstrated.** The continuation produces a *single closing
  arc* at the true ephemeris with a 0.00158 km/s closure residual. It does **not**
  propagate 3–5 laps over ~20–30 yr and report a *bounded horizon TCM within a ΔV
  budget*. That horizon-TCM accounting is precisely what V3 requires and what V2
  defers to V3:

  > "The full geometric-modulation horizon (~7 laps / ~15 yr …) is **V3's**
  > burden, not V2's. V3's phase-matched horizon-TCM gate (3–5 laps, ~20–30 yr) …"
  > — spec §14 note

  A one-arc closure residual is not a 3–5-lap horizon-TCM budget.

**"ephemeris-positions two-body" vs our n-body harness fidelity rungs.** Russell's
Chapter-5 force model — and therefore the continuation endpoint — is **patched-conic
two-body legs between true-ephemeris planet positions**, NOT n-body integration:

> "the dissertation's force model is patched-conic two-body legs between true
> ephemeris planet positions" — `…continuation-deepdive.md` §1

The spec's **V4** ("High-fidelity external … Tudat/pykep n-body") and our planned
n-body harness (REBOUND, planets-on-rails) are the n-body fidelity layer. The
continuation's "true ephemeris" is *ephemeris planet positions with two-body
spacecraft legs* — the same fidelity class as V3's "astropy backend / ephemeris
realisation", and explicitly the precursor (not the substitute) for the n-body
harness:

> "our n-body harness IS the higher-fidelity version of his step 3" —
> `…continuation-deepdive.md` §5

So: the evidence bears on **V3** (ephemeris realisation, two-body legs), and is
**below** the n-body harness rungs (V4-class). It does not reach into the n-body
fidelity layer at all, and it does not complete V3's horizon-TCM requirement.

**V3 verdict: NO (but it is real V3-relevant progress).** The phase-match to a
real window is V3-shaped and sourced; the bounded 3–5-lap horizon-TCM budget — the
binding V3 deliverable — is absent. The model is ephemeris-positions two-body
(V3-class fidelity), not n-body, so this never reaches V4/n-body.

---

## 4. The honest counter-arguments — do they undermine the claims?

All three are real and were surfaced by the #158 authors themselves; none changes
the V2/V3 verdicts (those already say "no" on stronger structural grounds), but
they reinforce *why a writeback would be premature* and why this is corroboration,
not a like-for-like gate.

1. **Emerged V∞ E=6.08 / M=8.88 vs sourced anchors 6.5 / 9.7 — ~1.5 km/s off.**
   The test's own tolerance is `<= 1.5` and it "acknowledges the genome/model gap
   while still pinning the right family" (test docstring). 1.5 km/s is a *family
   confirmation*, not an identity-grade match. For attaching the result to the
   Aldrin row as *corroboration* this is fine; for any gate that asserts "this IS
   the Aldrin V∞ signature realised on the ephemeris" it is too loose. It argues
   *corroboration*, not *promotion*.

2. **t0 26 d from the sourced Aug-6-2003 window.** Inside the 200-day gate, so the
   phase-match claim survives. But 26 d is a real offset; the continuation found a
   nearby basin, not Russell's exact dated solution. Acceptable as evidence of the
   right launch window; not a basis for asserting an exact reproduction.

3. **e/a differ from the classic Aldrin elements.** The continuation converged
   `(a, e) = 1.5249 AU, 0.3616`; the catalogue Aldrin row is `a=1.60, e=0.393`
   (Rogers 2012). These differ (a by ~0.075 AU, e by ~0.03) — because the
   continued solution is the *true-ephemeris* realisation, not the circular-coplanar
   idealisation the row stores. This is the crux: **it is the same family in a
   different model**, which is exactly why it is corroboration for the
   circular-coplanar row rather than a like-for-like reproduction of it. (Note the
   row already documents an a/e source discrepancy — spec §9's a≈1.659/e≈0.41 vs
   Rogers' 1.60/0.393 — so a third, ephemeris-model number does not "correct" the
   row; it is a different-model value.)

**Counter-argument verdict:** none undermines the (already-negative) V2/V3
verdicts. Collectively they confirm the result is **family-level corroboration of
the Aldrin row in the true-ephemeris model**, not a model-matched gate pass.

---

## 5. Recommendation

**HOLD on any validation-level promotion. Record the continuation result as
corroborating evidence on `aldrin-classic-em-k1-outbound`, not as a level bump.**

Concretely:

- **Do NOT** promote `aldrin-classic-em-k1-outbound` past its current **V2
  (V2-powered)**. The continuation does not produce the V2-ballistic instrument
  (bounded rotating-frame drift over ≥3 laps), and V2-ballistic is the wrong-model
  gate for a circular-coplanar/powered row anyway (§2).
- **Do NOT** stamp **V3**. V3's binding deliverable — a bounded 3–5-lap (~20–30 yr)
  ephemeris-mode horizon TCM within a ΔV budget — is not present; only a single-arc
  closure residual and a phase-match-to-window are (§3).
- **Do NOT** create a `6.399G1` row or attach this to any `russell-ocampo-*` or
  `russell-ch4-4.991gG2` row. There is no such catalogue row and the S1L1 member is
  a documented break point in the same run (§1).

**What it MAY justify (a non-level, additive record), if the user wants it now:**
a corroborating-evidence pointer on `aldrin-classic-em-k1-outbound` noting that the
#158 continuation driver walked the Aldrin family from circular-coplanar to the
true ephemeris (DE440), closing ballistically to 0.00158 km/s at t0 26 d from the
sourced Aug-6-2003 window, reproducing Russell Table 5.5's 0-m/s `6.399G1 (#1)`
entry — *in the true-ephemeris (Russell Ch.5 two-body-legs) model*, explicitly NOT
the row's circular-coplanar defining model, and NOT a §14 level change. This is the
same honesty convention the row already uses for `maintenance_dv_kms_per_synodic`
(computed, not source-attested) and is value-adding without over-claiming.

**The cheap additional verification that WOULD support a future V3 claim** (in
priority order; none is a <10-min slot-in, so not run here):
1. Build a `Cycler` from the continued true-ephemeris solution and run
   `verify/propagate.py` for ≥3 laps on DE440, recording dynamic-rotating-frame
   drift — this is the *true-ephemeris-model* V2-ballistic instrument (a NEW,
   model-scoped claim, not the circular-coplanar row's V2).
2. Extend that to the V3 horizon: 3–5 laps (~20–30 yr) with ephemeris-mode TCM
   accounting and a stated ΔV budget — the actual V3 deliverable.
3. Tighten the V∞ match (currently 1.5 km/s) toward identity grade, or document
   the residual gap as the genome/model limitation it is.

**Exactly what evidence chain a writeback WOULD cite (when earned).** Per
`src/cyclerfinder/data/validate.py`, any level above V0 requires a
`_LEVEL_EVIDENCE[(rid, level)]` entry pointing at recorded in-repo mechanical test
evidence (the over-claim guard). A future V3 entry for
`("aldrin-classic-em-k1-outbound", "V3")` would have to cite:
- the sourced EXPECTED: Russell 2004 dissertation Table 5.5 (p.178), `6.399G1 (#1)`
  Aldrin, total Δv = 0 m/s / 7 cycles, launch Aug 6 2003;
- the in-repo EVIDENCE: a *passing, non-xfail* test that propagates the continued
  solution over 3–5 laps on DE440 and asserts bounded horizon TCM within budget
  (which does **not exist yet** — `test_continuation.py` asserts a single-step
  corrector residual `< TOL_KMS` and window proximity, not a horizon-TCM budget);
- the model scope recorded as **true-ephemeris (patched-conic two-body legs,
  Russell Ch.5)** — explicitly distinct from the row's circular-coplanar defining
  model, and explicitly **not** "DE405" or n-body (per the deep-dive §6 backfill
  rule: a final-step solution is `accurate-ephemeris (patched-conic, Russell
  2004)`).

Until that test exists and passes, the golden-discipline floor (`when in doubt,
V0`; here, hold at the earned V2-powered) applies. The continuation result is an
excellent, sourced, end-to-end validation of the *driver* and strong family-level
corroboration — but the ladder is mechanical-evidence-gated, and the mechanical
evidence for a level bump is not yet on disk.
