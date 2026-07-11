# #568 — Writeback-readiness verdict: #312 + the 5 #566 representatives

**Date:** 2026-07-11
**Author:** Opus (trust-bearing interpretation pass, per #568 recommended model)
**Inputs (no new computation):** `data/scan_567_epoch_robustness.jsonl` (committed `a777c2d`, post-#567(1)+(2)
bugfix code `6c54bba`); the #567 step-(3) singleton-flip diagnosis
(`docs/notes/2026-07-11-567-post-fix-singleton-flip-diagnosis.md`) and its Fable adversarial
corrections (folded into #568's scope in `data/OUTSTANDING.md`).
**Scope:** produce the per-candidate writeback-readiness verdict. Does NOT edit `data/catalogue.yaml`
(separate writeback step, proposed below as #569). Does NOT run new computation. Does NOT re-litigate
the #567 diagnosis (settled + Fable-confirmed: the isolated flips are genuine synodic-phase aliasing of
a real `planet_crossing_infeasible` feasibility boundary, no third solver artifact).

---

## 1. Method — what "duty cycle" means and how it was computed

The #567 diagnosis established that the raw PASS% at any single fixed epoch grid is **not** a wall-clock
`validity_window`: each candidate has a genuine, fast, synodic-phase-dependent feasibility boundary (the
transfer's osculating periapsis crossing Uranus's `r_eq` = 25,559 km), recurring every synodic period.
The meaningful figure is the **feasible synodic DUTY CYCLE** — the fraction of synodic phase for which the
cycle is V4-strict-feasible.

The daily-2000 (366 pt) and daily-2030 (365 pt) sweeps sample launch epoch at 1-day spacing across many
synodic cycles (≈120 cycles for #312's 6.0 d boundary, ≈228 for the ~3.2 d Ariel pairs). Because the
1-day grid is incommensurate with each pair's synodic period, the 731 combined daily samples are
**phase-random** — so the PASS *fraction* converges to the true duty cycle by the law of large numbers
**even though individual cells alias** (the aliasing is exactly what produces the isolated-flip pattern
`detect_isolated_singleton_anomalies` reports; per the #567 diagnosis that flip% is a DIAGNOSTIC of the
aliasing, NOT an instability discount, and is not used against any duty cycle here).

**Fable correction #2 applied — drift-floor rows excluded from the denominator.** The duty cycle is
computed over `planet_crossing_infeasible`-mode FAILs only:

```
duty_cycle = n_PASS / (n_PASS + n_planet_crossing_infeasible)      [daily-2000 + daily-2030 combined]
```

The tiny `converged`-but-drift-floor-exceeded population (a conceptually DIFFERENT failure mechanism, see
§3) is set aside from this fraction and reported separately, not silently folded in.

**Fable correction #1 applied — measured boundary periods, no "zero-offset ⇒ P/2" shortcut.** Each
candidate's synodic boundary period is taken from the MEASURED daily PASS/FAIL series (I re-ran the FFT
myself on the committed data; §4), not from any derived halving rule. Only Ariel-Umbriel actually locks
at synodic/2; the three other `rel_offset=0°` symmetric candidates lock at their FULL synodic period.

---

## 2. Per-candidate feasible synodic duty cycle (daily-2000 + daily-2030, N = 731 each)

| candidate (pair) | rel_offset | PASS | PC-fail | drift-fail | **duty cycle** | raw daily PASS% |
|---|---|---|---|---|---|---|
| **#312** Umbriel-Oberon (catalogued) | 180° | 570 | 151 | **10** | **79.1 %** | 78.0 % |
| Titania-Oberon (row 23, mandatory) | 180° | 544 | 187 | 0 | **74.4 %** | 74.4 % |
| Ariel-Oberon (row 7) | 0° | 519 | 212 | 0 | **71.0 %** | 71.0 % |
| Umbriel-Titania (row 10) | 0° | 501 | 230 | 0 | **68.5 %** | 68.5 % |
| Ariel-Titania (row 5) | 0° | 486 | 245 | 0 | **66.5 %** | 66.5 % |
| Ariel-Umbriel (row 1) | 0° | 451 | 280 | 0 | **61.7 %** | 61.7 % |

- **Every fail across all 5 representatives is `planet_crossing_infeasible`** — zero drift-floor fails,
  zero `lambert_no_solution`, zero `integrator_failure`. Their duty cycle therefore equals their raw
  daily PASS% exactly; the drift-floor correction moves only #312.
- Duty cycles span **61.7 %–79.1 %** — all comfortably in the majority-feasible regime. None is a
  knife-edge (the #567 sub-daily zoom showed the boundary is a smooth, deep periapsis excursion, not a
  jittery near-threshold crossing).

## 3. #312's secondary drift-floor-exceeded population — reported SEPARATELY (Fable correction #2)

Distinct from the synodic planet-crossing boundary, #312 alone shows a small population where **all legs
converge but the V4-vs-V3 drift exceeds the 50,000 km floor** — a threshold on a different continuous
metric (real high-perturbation epochs), confirmed genuine physics in the #567 diagnosis, not a
knife-edge artifact. Enumerated from the scan:

| sweep | epoch | drift vs V3 (km) |
|---|---|---|
| annual | 2029-06-21 | 58,059 |
| daily-2000 | 2000-01-03 | 55,846 |
| daily-2000 | 2000-01-15 | 90,617 |
| daily-2000 | 2000-01-27 | 57,560 |
| daily-2000 | 2000-02-20 | 97,390 |
| daily-2000 | 2000-04-08 | 71,432 |
| daily-2000 | 2000-08-03 | 50,084 |
| daily-2030 | 2030-07-03 | 54,644 |
| daily-2030 | 2030-07-12 | 103,413 |
| daily-2030 | 2030-07-27 | 84,702 |
| daily-2030 | 2030-10-01 | 52,952 |

**11 / 831 rows total (~1.3 %); 10 of those in the 731 daily rows (~1.37 %). Range 50,084–103,413 km.**
This is #312's own characteristic — **none of the 5 representatives exhibit it at all**. It is NOT
folded into #312's 79.1 % duty cycle; it is a separate, minor, sporadic secondary-window caveat.
(The 2000-01-15 row here is exactly the epoch #338 originally flagged as a ~91,000 km FAIL — now
correctly classified as this drift-floor mechanism, not a solver artifact.)

## 4. Synodic boundary period — MEASURED per candidate (Fable correction #1)

FFT dominant period of each candidate's daily `planet_crossing_infeasible` boundary series, computed
directly from the committed scan (2000 / 2030 windows), matching the #567 diagnosis table:

| candidate | measured boundary period | synodic period | relation |
|---|---|---|---|
| Titania-Oberon | **24.4 d** (2000) / 24.3 d (2030) | ~24.6 d | full synodic |
| Umbriel-Titania | **8.0 d** / 7.9 d | ~7.9 d | full synodic |
| #312 Umbriel-Oberon | **6.0 d** / 6.0 d | ~5.99 d | full synodic |
| Ariel-Titania | **3.6 d** / 3.5 d | ~3.55 d | full synodic |
| Ariel-Umbriel | **3.2 d** / 3.2 d | ~6.43 d | **synodic / 2** (only one that halves) |
| Ariel-Oberon | **3.1 d** / 3.1 d | ~3.10 d | full synodic |

Confirms Fable correction #1: the `rel_offset=0°` symmetric structure does **not** imply period-halving
— Ariel-Titania, Ariel-Oberon and Umbriel-Titania share it but lock at their full synodic period;
only Ariel-Umbriel halves. The MEASURED value is used, not a derived rule.

**Sampling note (non-gating):** the three shortest-boundary candidates (Ariel-Oberon 3.1 d, Ariel-Umbriel
3.2 d, Ariel-Titania 3.5 d) sit near the daily grid's Nyquist limit. Individual cells alias (that IS the
flip signal), but the *duty-cycle fraction* is averaged over ~200+ synodic cycles and is well-converged;
the tier verdict below does not turn on the exact percent, so no sub-daily refine is needed to reach it.
If a future mission-window statement needed the exact duty cycle pinned to <1 %, a phase-uniform
sub-synodic grid on those three would tighten it — a refinement, not a blocker.

---

## 5. Per-candidate writeback-readiness verdict

**All six clear at #312's own V4-strict validation level (windowed).** A duty-cycle-limited feasibility
band is **not** a validation defect — it is exactly the epoch-locked, windowed validity the `quasi_cycler`
class is defined for (`epoch_locked=true` + `validity_window`; the taxonomy note
`docs/notes/2026-06-16-catalogue-scope-taxonomy.md` §"two coexisting semantics" uses #312 *itself* as the
worked example, with the moon-pair synodic period as the system-period). A low duty cycle is a
launch-epoch-must-be-timed-to-synodic-phase constraint — the same kind of constraint as an interplanetary
launch window — and belongs in the row's `validity_window`, not as a downgrade.

| candidate | duty cycle | boundary period | verdict |
|---|---|---|---|
| **#312** Umbriel-Oberon | 79.1 % (+1.3 % drift-floor caveat) | 6.0 d | **V4 — already catalogued; UPDATE row** (family reframing + duty cycle + drift-floor caveat; its current `{2000→2083}` window misleadingly implies *continuous* validity) |
| Titania-Oberon | 74.4 % | 24.4 d | **V4-equivalent — READY** (widest boundary → fewest flips; mandatory literature-clearance rep) |
| Ariel-Oberon | 71.0 % | 3.1 d | **V4-equivalent — READY** |
| Umbriel-Titania | 68.5 % | 7.9 d | **V4-equivalent — READY** |
| Ariel-Titania | 66.5 % | 3.5 d | **V4-equivalent — READY** |
| Ariel-Umbriel | 61.7 % | 3.2 d | **V4-equivalent — READY** (lowest duty cycle, still majority-feasible; report the window, do NOT cap) |

Grounds the 5 representatives clear at #312's level:
1. **Same failure mechanism, no artifacts** — every fail is the genuine `planet_crossing_infeasible`
   geometry boundary; the two #559 artifact generators are confirmed gone (`lambert_no_solution=0`,
   `integrator_failure=0` in all 18 groups).
2. **Same gauntlet already passed** — #566 ran all 5 through the full V2→V3→V4-scipy→V4-strict chain at
   the anchor epoch, passing identically to #312.
3. **Wide, smooth window** — 62–79 % duty cycle, none knife-edge (smooth deep periapsis excursion per
   the #567 zoom).
4. **Cleaner than #312 on the drift axis** — zero drift-floor failures vs #312's own ~1.3 %.

**Is the 61.7 % (Ariel-Umbriel) low duty cycle disqualifying?** No. It is a real, reportable
`validity_window` characteristic, not a defect. The candidate is V4-strict-feasible across the majority of
its synodic phase, on a smooth boundary, with the same clean mechanism as #312. Capping it below the
others would misread a legitimate launch-window constraint as an instability — precisely the error the
#567 diagnosis warns against.

**One writeback-MECHANICS qualifier (not a validation gap):** #312's catalogue V4 carries a frozen-gate
provenance convention (per-tier JSONL verdicts + `tests/verify/test_silver_327_v*.py` frozen pytest gates
registered in `src/cyclerfinder/data/validate.py::_LEVEL_EVIDENCE`, spec §16.7.12). The 5 representatives
have the #566 gauntlet JSONL + this #567 duty-cycle evidence, but the frozen per-row gate artifacts must
be **generated** as part of the writeback to match the registration convention. That is mechanical
writeback scope (#569 below), not a shortfall in the validation itself.

---

## 6. Overall writeback recommendation

**Proceed with catalogue writeback (as a separate task, #569 below) of all 5 representatives as validated
`quasi_cycler` family members at #312-equivalent V4 (windowed), each carrying its own measured
`validity_window` (duty cycle + synodic boundary period).**

Each new row should say (tying to #564 §5's reframing, now with REAL per-candidate data replacing the
placeholder):
- `orbit_class: quasi_cycler`, `epoch_locked: true`, `model_assumption: cr3bp`, one row **per moon pair**
  (census-minimum = one representative per non-Miranda pair, per #565 §2/§4).
- The **family relation**: each is a member of the **30-member symmetric-closure family** (six non-Miranda
  Uranian moon-pair directions × exact 0°/180° symmetric closures, #563 enumeration), of which **#312 was
  the first documented member**. #312 is NOT a unique novel point but the first-found representative of a
  structured family.
- `validity_window`: the calendar `{start, end}` span (bounded by URA111 kernel coverage) **plus** the
  measured duty cycle and synodic boundary period from §2/§4 — the system-period-relative semantics the
  taxonomy note requires for a moon-system quasi_cycler.
- Literature clearance: Titania-Oberon is the flagged mandatory literature-clearance obligation (#565);
  re-confirm the `literature_check.py` not-found against Canales/Kumar for it before its row is written
  (structurally adjacent, not a direct topology match, per the #562 finding).

**Update #312's own existing row** — materially out of date on two counts:
1. **Family reframing** — add the "first-documented member of a 30-member symmetric-closure family"
   relation (was written when #312 looked like an isolated novel point).
2. **Duty-cycle / window correction** — its current `validity_window: {2000-06-21 → 2083-06-21}` reads as
   *continuous* 83-year validity. The #567 characterization shows it is a ~79.1 % synodic duty cycle
   (feasible band recurring every 6.0 d) within that calendar span, PLUS a distinct ~1.3 % secondary
   drift-floor-exceeded population (50–103 k km). Both belong in the updated `validity_window` /
   notes block; the existing "sub-year DOY sensitivity flagged" note is now superseded by this concrete
   duty-cycle + drift-floor characterization.

**Do not** change any V4-strict driver code (the #567 diagnosis found no third bug — the aliasing is real
physics), and do not extend to the other 25 same-pair-redundant closures (census-minimum is one rep per
pair).

---

## 7. Headline

| candidate | duty cycle | synodic boundary | drift-floor 2nd pop | verdict |
|---|---|---|---|---|
| #312 Umbriel-Oberon | 79.1 % | 6.0 d | ~1.3 % (50–103 k km) | V4, UPDATE row |
| Titania-Oberon | 74.4 % | 24.4 d | none | V4-equiv, READY |
| Ariel-Oberon | 71.0 % | 3.1 d | none | V4-equiv, READY |
| Umbriel-Titania | 68.5 % | 7.9 d | none | V4-equiv, READY |
| Ariel-Titania | 66.5 % | 3.5 d | none | V4-equiv, READY |
| Ariel-Umbriel | 61.7 % | 3.2 d (synodic/2) | none | V4-equiv, READY |

**Bottom line:** all six are writeback-ready at #312's V4-strict-equivalent level. The duty cycles
(62–79 %) are legitimate windowed-validity characteristics to record, not disqualifiers. Writeback should
proceed as #569, reframing the whole set as a 30-member symmetric-closure family (first-documented member
#312), with per-candidate measured `validity_window` data — after a Fable second-opinion pass, matching
this chain's discipline.
