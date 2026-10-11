# #1007 V2-ballistic campaign for gc-1, gc-2, ev-C, ev-A — PRE-REGISTRATION (2026-10-08), committed before any run

Spec basis: the 2026-10-08 owner amendment to spec §14 V2-ballistic (commit 93170995).

The amendment covers patched-conic multi-arc rows whose period is exactly commensurate in the ideal
model. For these, V2 is judged in the real-ephemeris chain:
- >= 3 continuous ballistic cycles;
- bounded drift;
- the gate passes at every flyby;
- the DOP853 re-fly misses < 1 km;
- at >= 3 of 5 epochs.

A positive and a negative control run first.

## 1. What is judged (the existing rung (d) chains; no new solve for the rows)

| Row | Model | Chain | Epochs |
|---|---|---|---|
| gc-1 | NAIF jup365, blend+shoot, 10 cycles | `data/943_gc1_realeph/e0..e4/realeph_chain.json` | the 5 standard epochs (2030-2056) |
| gc-2 | NAIF jup365, 10 cycles | `data/942_943_ladder_realeph_chains.json` ["gc2"] | 5 |
| ev-C | Standish mean elements, ramp, 5 cycles (16 yr) | `data/942_direct_route_and_gm_fix.json` ["evC_standish_rerun"] | 5 |
| ev-A | Standish mean elements, ramp, 5 cycles | ["evA_direct"] | 5 |

- Heliocentric rows are judged on Standish, which is the rung (d) model for ev.
- DE440 is reported, not judged:
  - ev-C on DE440 is ballistic at 5/5 (6.23);
  - ev-A on DE440 is near-ballistic, 2.18-3.55 m/s per 7 cycles (6.57a), so not ballistic there.
- No extra cycles are needed: every chain already has >= 5 cycles.

## 2. Criteria at one epoch (all must hold)

1. Continuity: the chain's re-evaluated residual at lambda = 1 is < 1e-6, in km/s (|V_inf| magnitude
   match at every junction). Shot fixed legs also carry their arrival miss / 1000 km.
2. Gate: the #888/#937 demanded-turn gate passes at every interior flyby ("indeterminate" is not a
   pass). Uses the chain tool's own gate:
   - ramp: `interior_gate`, minimax full-rev directions;
   - blend+shoot: the shot directions, `shoot_best`.
3. Re-fly: `scripts/check_942_realeph_chain.py` logic.
   - Ramp: `cross_check`, max arrival miss < 1 km and V_inf vector error < 1e-6 km/s.
   - Blend+shoot: segment re-fly with the real GM, max arrival miss < 1 km.
4. Bounded drift (project convention, fixed here before any row's drift is computed):
   - Notation: the chain has n cycles and L Lambert legs per cycle. t_ij is the start date of Lambert
     leg j in cycle i (t_00 = the chain start). P = the chain period / n.
   - Date drift: D_i = max_j |t_ij - t_0j - i P| (the departure from a strictly periodic chain at the
     chain's own mean period).
   - |V_inf| drift: W_i = max_j | |V_inf,dep|_ij - |V_inf,dep|_j(template) |, where the template is the
     row's ideal-model zero evaluated in the ideal model.
   - BAND:
     - max_i D_i <= 0.05 x the template's shortest Lambert-leg flight time;
     - max_i W_i <= 0.10 x the template's smallest |V_inf,dep|.
   - NO GROWTH, for D and W each: the max over the second half of cycles 1..n-1 is <= 2 x the max
     over the first half. Growth below the floors 0.5 d (D) and 0.05 km/s (W) is ignored.

A row passes V2-ballistic if all four criteria hold at >= 3 of its 5 epochs.

## 3. Controls (run first; expected outcomes stated now)

POSITIVE, each must PASS the same criteria at >= 3 of 5 epochs:
- ev: Hollister 1H (Hollister 1969 orbit I; key `k2|RE/1:1|LE>V/0s|RV/1:1|RV/1:1|LV>E/0s`), Standish,
  5 cycles, at the 5 standard epochs.
  - Method: the RAMP CONTINUATION (no `--direct`), the same as the ev-C rows' rung.
  - Expected: unknown. The `--direct` run (6.23 D1) converged at 2/5 only; the continuation has not
    been run.
  - If 1H fails, the lane has no ev positive control. STOP and report (no substitution after seeing the
    result).
- gc: R-S GanEur#316 (R-S 2009 Table 3), NAIF jup365, 10 cycles, at the 5 standard epochs.
  - Chain: the existing `data/943_ganeur316_realeph/n10_std`; plus its R-S epoch (2019) chain
    `n10_rs2019_rel`, reported.
  - Expected: PASS at 3/5 if the drift criterion holds (the chain and gate passed at 3/5).
- R-S GanCal#1 at 2013 (owner-named, R-S Fig. 9(b); `data/943_c4_rs2013/n10_grow`): REPORTED, not a
  pass/fail control.
  - Its gate is "indeterminate" (0.9917) by construction, as in R-S's own ideal model (owner ruling
    6.45), so it cannot pass criterion 2 as written.
  - It is judged on criteria 1, 3 and 4 as the full-revolution-path control.
  - Deviation from the owner's list, stated here; GanEur#316 is the gc control that must pass.

NEGATIVE, each must FAIL. Each takes the same chain with ONE interior Lambert start date shifted, not
re-solved:
- ev: ev-C epoch 0, the 2nd Lambert start + 1.0 d.
- gc: gc-1 epoch 0, the 2nd Lambert start + 0.05 d.

Expected: criterion 1 fails (junction mismatch >> 1e-6). This is the #830 lesson: the instrument must
see a broken chain.

## 4. Process

- `scripts/v2_chain_1007.py`: controls first, then the four rows. Output `data/1007_v2_chain.json`.
- Report per row and epoch: n, the 4 criteria, D and W maxima, growth ratios, verdict.
- If a positive control fails or a negative control passes: STOP and report.
- Row-level V2 promotions go in the next catalogue window, with the ratchets (lead).

## 5. Result of the first control (2026-10-08): the ev POSITIVE CONTROL FAILS -> STOP

Hollister 1H: Standish, 5 cycles, ramp continuation, the 5 standard epochs (`data/1007_v2/hollister1H_ramp/`).
The continuation stalls at every epoch:

| Epoch (JD) | Last converged lambda |
|---|---|
| 2462984.5 | 0.288 |
| 2465314.3 | 0.391 |
| 2467656.4 | 0.375 |
| 2469988.2 | 0.366 |
| 2472323.5 | 0.350 |

- 0/5 epochs give a chain. For comparison, the `--direct` solve (6.23 D1) converged and passed the
  gate at 2/5.
- Per sec. 3: the lane has no ev positive control, so the campaign STOPS and is reported. No
  substitution was made after seeing the result.
- The gc controls, the negative controls and the rows were NOT run.

## 6. gc lane results (2026-10-08; `scripts/v2_chain_1007.py --lane gc`, `data/1007_v2/v2_chain_gc.{json,log}`)

Controls, as pre-registered:
- POSITIVE, R-S GanEur#316, 10 cycles: PASS at 3/5 epochs (2462503.8, 2464844.7, 2469519.5). The other
  two epochs have no stored chain (the 6.4x run did not converge there). The R-S-epoch 2019 chain also
  passes all four criteria.
- NEGATIVE, gc-1 epoch 0 with the 2nd Lambert start + 0.05 d: FAILS as required (junction mismatch
  0.23 km/s, re-fly miss 1.3e4 km).
- GanCal#1 at 2013 (reported): criteria 1, 3 and 4 hold; the gate is "indeterminate" 0.992, as
  expected.

Rows, LITERAL verdict under sec. 2:

| Row | Epochs passing | Failing criteria |
|---|---|---|
| gc-1 | 0/5 | criterion 1 at all 5 (the mixed residual 1.4-1.9e-4); also criterion 4 at 2 epochs (D 0.222 and 0.195 d against the 0.190 band) |
| gc-2 | 2/5 (2464840.6, 2469524.3) | criterion 4 at 3 epochs: the W band (0.340-0.358 km/s against 0.304). D and growth hold everywhere |

gc-2 needed one adapter fix before the run: the ladder file stores 4n dates then the period, not an end
date. The first run read it wrongly ("chain-eval-fail" at 5/5). Corrected before judging; the fix is
format only.

Two findings about the pre-registered instrument, for the lead (NOT applied):
1. Criterion 1 as written ("shot fixed legs also carry their arrival miss / 1000 km") contradicts the
   spec amendment, which reads "|V_inf| continuity at every junction < 1e-6 km/s".
   - gc-1's shot full-rev legs close to 0.16-0.19 km (that is criterion 3's business, < 1 km).
     Scaled by 1/1000 this gives 1.6-1.9e-4 in criterion 1.
   - The junction-only |V_inf| residual is 8.6e-7 at every gc-1 epoch (< 1e-6).
   - Under the spec's wording gc-1 would pass criterion 1 at 5/5, and overall at 3/5 (epochs 0-2;
     epochs 3 and 4 fail the D band).
2. W measures the real-versus-ideal |V_inf| OFFSET as well as drift. gc-2's W is 0.34-0.36 km/s with
   no growth (W growth 1.0-2.0); this is a systematic ephemeris offset, not lap-over-lap drift. The
   growth test already captures drift.

The verdicts above are the pre-registered ones. Any change of criterion is the lead's/owner's.

## 7. Amended rules and re-judgement (2026-10-11), recorded as made AFTER seeing the sec. 6 results

- Criterion 1 (lead ruling 2026-10-08): corrected to the spec amendment's own wording, junction-only
  |V_inf| continuity < 1e-6 km/s. The shot legs' arrival miss belongs to criterion 3. This is a
  post-hoc correction of a pre-registration error against the spec text, which predates the results.
- W (OWNER decision 2026-10-08 evening): "bounded drift" is judged by GROWTH only. The systematic
  real-versus-ideal |V_inf| offset is recorded in the row, not banded.
  - Changed after seeing results; recorded as such.
  - The D band and both growth tests are unchanged.

Run: `scripts/v2_chain_1007.py --lane gc --rules amended` (`data/1007_v2/v2_chain_gc_amended.{json,log}`).
- Controls are unchanged in outcome: GanEur#316 passes 3/5; the broken gc-1 chain still fails
  (junction 0.23 km/s).

| Row | Registered rules | Amended rules | Notes |
|---|---|---|---|
| gc-1 | 0/5 | **3/5 PASS** (epochs 2462511.7, 2464840.9, 2467182.8) | MARGINAL: junction residual 8.6e-7 against 1e-6. Epochs 3-4 fail the D band (0.222 / 0.195 d against 0.190) |
| gc-2 | 2/5 | **5/5 PASS** | Systematic |V_inf| offset to record: W 0.16-0.36 km/s, growth 1.0-2.0 |

V2 promotions for gc-1 and gc-2 go in the next catalogue window (lead).

## 8. ev positive control, take 2: Hollister-Menning 1970 orbit 1 at its own epoch — PRE-REGISTERED (2026-10-11), before running

Lead ruling (a), 2026-10-08. Lesson recorded from sec. 5: a published object with no epoch (Hollister
1969 orbit I) is a poor real-ephemeris control.

- Object: H&M 1970 Table 3 orbit 1.
  - Its five blocks are identical (Earth full-rev, Venus full-rev twice; `block_types`), so it is a
    true repeating k = 2 cycle, as the drift metric needs.
  - Its structure is Hollister 1H: key `k2|RE/1:1|LE>V/0s|RV/1:1|RV/1:1|LV>E/0s`, ideal dates
    474.597193, 1100.618380 d.
  - Chosen over orbit 2 by stage-1 reproduction quality (0.68 of encounters within tolerance against
    0.64). It is not one of the open orbits 3-8.
- Epoch: H&M's own. Table 3 orbit 1's first Earth date is JD 2440441, used as `--first-epoch-jd
  2440441 --epochs 1`. The chain tool phase-matches the ideal cycle near that date.
- Length: 5 cycles = 16 yr = H&M's published period (5 blocks); >= 3 cycles as the amendment requires.
- Model: Standish J2000 mean elements, the ev rows' rung (d) model. H&M used an inclined-elliptic model
  of the same kind.
- Route, in this order: the direct lambda = 1 solve (`--direct`); if it does not converge, the ramp
  continuation (no `--direct`). The ev rows used both routes.
- Criteria: the four of sec. 2 under the sec. 7 amended rules (junction-only continuity; drift by D band
  and growth only).
- PASS rule for this control: all four criteria hold at its single published epoch. H&M published one
  epoch, so the 3-of-5 rule does not apply to it.
- Expected: PASS. H&M computed this orbit as periodic at these dates in an inclined-elliptic model.
- If it FAILS: STOP and report. The ev V2 question then waits (lead ruling).
- The ev NEGATIVE control (sec. 3: ev-C epoch 0 with the 2nd Lambert start + 1 d) runs with it and must
  fail.

### 8.1 Result (2026-10-11 09:21 AEDT): the H&M orbit 1 control FAILS -> STOP

- The chain tool phase-matched the ideal cycle at JD 2440790.5. Its first Lambert start (E->V) is 15 d
  before H&M's printed E->V departure (JD 2440806), so this is the same 1969 window.
- Direct solve (`data/1007_v2/hm_orbit1_direct/`): no convergence (residual 0.42 km/s).
- Ramp continuation (`data/1007_v2/hm_orbit1_ramp/`): stalls at lambda 0.398. At that lambda the chain
  is closed and gate-passing (worst 0.738), but it never reaches lambda = 1.
- This is the same failure mode as Hollister 1H at the 2030-2056 epochs (sec. 5).

Per sec. 8 and the lead's ruling, the ev V2 question WAITS. The negative control was not run.

Diagnosis (not tested): both ev controls fail at the IDEAL -> REAL continuation, not at the criteria.
The full-rev-only 1H/orbit-1 branch of the circular model does not continue to the inclined-elliptic
orbit H&M published. A route that seeds the chain from H&M's PRINTED dates was not pre-registered and
was not run. It would test the published object directly (stage-1 of the H&M reproduction converged
from those dates; `data/942_hm1970_control_summary.json`).

### 8.2 ev positive control, take 3: H&M orbit 1 seeded from its PRINTED dates — PRE-REGISTERED (2026-10-11), before running

Lead ruling 2026-10-11. This tests the published object directly; there is no ideal -> real
continuation.

- Seed: H&M 1970 Table 3 orbit 1 via `hollister_menning_1970.load_table3` (the printed dates with the
  module's documented print-error fixes; orbit 1 has none).
  - Lambert-leg starts per block k = 0..4: E->V at the block's second Earth date (row 5k+1), V->E at
    its third Venus date (row 5k+4).
  - Chain start x0 = the first E->V date (JD 2440806). Chain-period seed 5844 d (H&M's 16-yr repeat).
  - Key `k2|RE/1:1|LE>V/0s|RV/1:1|RV/1:1|LV>E/0s` x 5 cycles (orbit 1's five identical blocks).
- Model: Standish J2000 mean elements as the ev rows' rung (d) uses them: `MeanElementSystem()`, NO
  H&M period overrides. H&M assumed exact 16-yr periodicity; we do not, so the open chain's period is
  free.
- Solve: the chain tool's date residual (`date_chain_residual`, junction |V_inf| magnitudes) at the real
  model directly, Levenberg-Marquardt from the printed seed, no continuation.
- Judged by `scripts/v2_chain_1007.py`'s `judge` (ramp mode) under the sec. 7 amended rules. The template
  is the ideal 1H zero (474.597193, 1100.618380 d).
- PASS: the solve converges (< 1e-6 km/s) and all four criteria hold at this single published epoch.
- Expected: PASS. The solution should lie within a few days of the printed dates, since H&M round to
  one day and assume exact periodicity.
- Meaning:
  - A PASS shows the ev lane and the V2 criteria accept the published object at its epoch. The failures
    in secs. 5 and 8.1 were then the ideal -> real continuation route, not the lane.
  - A FAIL means that, at its own dates, the lane cannot reproduce the published 16-yr orbit as a
    V2-ballistic chain.
  - Either way, no ev-member V2 judgement follows from this alone; the ev V2 question stays with the
    owner.
- The ev NEGATIVE control (ev-C epoch 0, the 2nd Lambert start + 1 d, no re-solve) runs in the same
  script and must FAIL.

### 8.3 Result of take 3 (2026-10-11 09:24 AEDT; `scripts/v2_hm_printed_1007.py`, `data/1007_v2/hm_orbit1_printed.{json,log}`)

POSITIVE control, H&M orbit 1 from its printed dates: **FAIL**, on criterion 4 only.
- The solve converges from the print: residual 1.14 km/s at the printed dates -> 6.0e-13. The largest
  date shift from the print is 2.88 d. The chain period is 5842.37 d (print 5844).
- Criterion 1 (continuity) holds. Criterion 2: the gate passes (worst 0.883). Criterion 3: re-fly
  3.4e-4 km.
- Criterion 4 FAILS:
  - D per cycle 0 / 18.8 / 19.0 / 13.1 / 39.5 d against the 8.8-d band;
  - D growth ratio 2.08 against 2.0.
  - W growth is 1.22, which holds; W's offset 2.2-3.5 km/s is recorded, not banded.

NEGATIVE control (ev-C epoch 0, + 1 d): FAILS as required (junction 0.22 km/s).

Per sec. 8.2, the ev lane cannot reproduce the published 16-yr orbit as a V2-ballistic chain under the
sec. 2/7 drift criterion. STOP; the ev V2 question stays with the owner.

Diagnosis (for the owner; nothing is changed here):
1. D measures GEOMETRIC BREATHING, not drift, for an Earth-Venus chain.
   - D compares each cycle with a strictly periodic chain at the chain's own mean period. The
     Earth-Venus geometry only repeats after 16 yr (5 cycles, H&M's period), so within one 16-yr span
     the cycles must differ.
   - H&M's own published orbit departs up to 39.5 d from 3.2-yr periodicity. That is the planets, not
     instability.
   - The 5-cycle window spans exactly one geometric period, so the growth test (second half against
     first half) also reads breathing as growth.
   - Spec §14 asks for drift "tolerant of geometric breathing". A breathing-tolerant form would compare
     SAME-PHASE cycles (cycle i against cycle i + 5) over >= 2 geometric periods (>= 10 cycles,
     32 yr). The gc lane passed because the Galilean geometry breathes far less.
2. The control's template is the wrong member.
   - H&M orbit 1 has 1H's STRUCTURE but not its V_inf: about 4.6 km/s (H&M V_r 0.155 EMOS), against
     1H's ideal 2.99/3.19.
   - So the ideal 1H zero is not orbit 1's circular counterpart. This also explains why the ideal ->
     real continuation from 1H (sec. 8.1) could not reach it.
   - Under the owner's growth-only W rule this does not change the verdict, but the template should be
     orbit 1's own circular counterpart.

### 8.4 DRAFT — same-phase drift metric. PROPOSED, NOT RUN, AWAITING THE OWNER (written 2026-10-11)

Why: sec. 8.3 showed that criterion 4's D compares each cycle with a strictly periodic chain. For
Earth-Venus it then reads the 16-yr geometric breathing as drift, and the published H&M orbit 1 fails.
Spec §14 asks for drift "tolerant of geometric breathing". This draft makes the criterion compare
cycles at the SAME geometric phase.

**Geometric repeat L (per system and cycle period P).** L is the smallest integer such that, for every
encounter body, L P / T_body lies within 0.02 rev of an integer (T = the body's real mean period).
Values (real periods; computed for this draft):

| Lane | P | L | L P |
|---|---|---|---|
| ev, k = 2 (ev-A, ev-C, H&M orbit 1) | 3.197 yr | 5 | 16.0 yr |
| ev, k = 3 (ev-B) | 4.796 yr | 5 | 24.0 yr |
| gc, k = 3 (gc-1, gc-2) | 37.57 d | 4 | 150.3 d |
| gc, k = 4 | 50.09 d | 1 | 50.1 d |
| gc, k = 5 | 62.62 d | 4 | 250.5 d |
| gc, k = 6 | 75.14 d | 2 | 150.3 d |

**Metric.** The chain must have n >= 2L cycles (Earth-Venus k = 2: >= 10 cycles, 32 yr). For
i = 0 .. n-L-1:
- D'_i = max over Lambert legs j of |t_(i+L)j - t_ij - L P_mean|;
- W'_i = max over j of | |V_inf,dep|_(i+L)j - |V_inf,dep|_ij |.

The chain is compared with ITSELF one geometric period later, so NO ideal-model template is needed.
This also removes sec. 8.3's second problem (the wrong-member template). The systematic real-versus-ideal
offset is still recorded in the row, not banded (the owner's 2026-10-08 decision stands).

Band (project convention, unchanged in spirit):
- max D' <= 0.05 x the chain's shortest Lambert-leg flight time;
- max W' <= 0.10 x the chain's smallest |V_inf,dep|.

No growth: over i, the max of the second half <= 2 x the max of the first half, ignoring growth below
0.5 d (D') and 0.05 km/s (W').

Criteria 1-3 are unchanged: junction-only continuity < 1e-6 km/s; the gate at every flyby; DOP853
re-fly < 1 km. The epoch rule is unchanged: >= 3 of 5 epochs; a published single-epoch control is
judged at its epoch.

**Controls (run first).**
- POSITIVE, ev: H&M orbit 1 from its printed dates (sec. 8.2 route).
  - Two H&M periods: the printed dates plus the same dates + 5844 d, 10 cycles, Standish, a direct LM
    solve.
  - Expected PASS. Only cycle i against i+5 is compared, so the 16-yr breathing cancels. What remains
    is the real non-commensurability (Earth 0.013, Venus 0.012 rev per 16 yr), absorbed by the chain's
    mean period. Its size is about the 2.9-d print-rounding scale against a band of about 8 d.
- POSITIVE, gc: R-S GanEur#316, the existing 10-cycle jup365 chains. GanEur#316's cycle and L must be
  computed first (k may differ from 3; if L > 5, the existing chains are too short and it needs a
  20-cycle run). Expected PASS at >= 3/5, as under the current metric.
- NEGATIVE: the broken-date chain (2nd Lambert start + 1 d), as before; it must fail continuity.
  - LIMIT: no chain with genuine secular drift is in hand, so the drift criterion itself has no
    negative control.
  - A candidate negative: the ev-B epoch-2 or epoch-4 chain. These converge but gate-fail, so they
    test the gate, not drift. I would state this limit in the pre-registration rather than invent a
    synthetic drift.

**Runtime.**
- ev: the H&M control and the ev rows need NEW 10-cycle chains. ev-A and ev-C on Standish are stored at
  5 cycles. That is 1 + 2 x 5 runs of seconds to about a minute each (5-cycle ramp chains took 10-15 s).
- gc: the existing 10-cycle chains suffice for L = 4 (n >= 8).
- Total well under 30 minutes. No long compute.

**The spec §14 wording change it needs** (in the 2026-10-08 amendment's "bounded drift" bullet):
- Current: "bounded drift: the per-cycle encounter-date offsets and |V∞| against the ideal-model
  template oscillate within a band without lap-over-lap growth".
- Proposed: "bounded SAME-PHASE drift: with L the system's geometric repeat (the smallest L for which
  every encounter body's mean-anomaly advance over L cycles is within 0.02 rev of an integer), the chain
  has ≥ 2L cycles, and the encounter dates and |V∞| of cycle i+L against cycle i stay within a band
  without growth. The systematic offset against the ideal model is recorded, not banded."

**Wording options for the owner (added 2026-10-11 after the sec. 8.5 pre-check):**
- (a) As drafted above: the comparisons run over i = 0 .. n-L-1.
- (b) START-CYCLE REFINEMENT: the comparisons run over i = 1 .. n-L-1. Cycle 0 is the open chain's
  start; its first departure leg has no incoming junction constraining |V_inf| and is tied only by the
  magnitude closure at the chain's far end. Under (b) the chain needs >= 2L + 1 cycles; the existing
  10-cycle gc chains (L = 4) still qualify, with 5 comparisons.
  - Spec wording under (b): "... the encounter dates and |V∞| of cycle i+L against cycle i, for every
    i except the open chain's first cycle, stay within a band without growth ...".
  - Effect (sec. 8.5, exploratory): gc-2 goes from 2/5 under (a) to 5/5 under (b). Every gc-2 failure
    under (a) is at i = 0 (W' 0.32-0.34 km/s), and all later comparisons are <= 0.022 km/s. gc-1 is
    5/5 under both. The H&M positive control passes under both. The Earth-Venus rows are in sec. 8.6.

**Effect on the gc rows already promoted to V2 (gc-1, gc-2; window of 2026-10-11).**
- They would be RE-JUDGED under the new wording on their existing 10-cycle chains (L = 4, n = 10,
  6 comparisons each), before the wording takes effect. A row that then fails 3/5 goes back to V1 in the
  next window, with the reason recorded.
- Expected: both HOLD.
  - D' removes the periodic part and keeps any secular part.
  - Their current D (0.03-0.22 d) is already small, and their growth ratios are 1.0-2.0.
  - gc-1's two D-band failures (0.222 / 0.195 d against 0.190) are the kind of accumulated phase that
    D' should not count. They may turn into passes, which would lift gc-1 from 3/5 to 5/5.
  - Not computed: this is an expectation, not a result.

Nothing in this section has been run.

### 8.5 EXPLORATORY pre-check of the sec. 8.4 draft metric (2026-10-11) — exploratory, under the draft metric, awaiting the owner

Not a spec or catalogue change. Run: `scripts/v2_samephase_1007.py` (`data/1007_v2/samephase_exploratory.{json,log}`),
on the EXISTING chains (gc, L = 4, n = 10, 6 comparisons). The H&M chain is new: solved from the
printed dates over two periods (10 cycles, L = 5, 5 comparisons). Continuity is junction-only.

| Chain | Continuity | D' max / band (d) | W' max / band (km/s) | Growth D / W | Draft drift |
|---|---|---|---|---|---|
| gc-1, 5 epochs | 8.6e-7 (pass) | 0.016-0.037 / 0.185 | 0.008-0.023 / 0.175 | 0.72-1.23 / 0.25-0.63 | PASS 5/5 |
| gc-2, JD 2464840.6 and 2469524.3 | 9.3e-8 (pass) | 0.038-0.040 / 0.13 | 0.157-0.160 / 0.30 | <= 1.01 / 0.06-0.07 | PASS |
| gc-2, JD 2462511.3, 2467182.5 and 2471853.7 | 9.2e-8 (pass) | 0.046-0.057 / 0.13 | **0.319-0.342** / 0.30 | <= 1.01 / 0.06 | FAIL (W' band) |
| NEGATIVE: gc-1 e0, + 0.05 d | **0.23 (fail)** | 0.043 / 0.185 | 0.25 / 0.175 | 0.74 / 0.02 | FAIL, as required |
| POSITIVE: H&M orbit 1, printed, 10 cycles | 1.4e-12 (pass) | 2.09 / 7.52 | 0.30 / 0.45 | 0.88 / 1.06 | **PASS** |

- The H&M solve converges within 6.75 d of the printed dates over the 32 yr. The chain period is
  11682.2 d against 2 x 5844.
- The draft metric passes the published orbit that the current metric failed (sec. 8.3).

Readings (exploratory):
- gc-1: 5/5 under the draft metric, against 3/5 (marginal) under the current wording. As sec. 8.4
  expected, the two D-band misses were accumulated phase.
- gc-2: 2/5 under the draft metric as written, against 5/5 under the current wording, so it would
  NOT hold.
  - Every gc-2 W' failure is at i = 0: cycle 4 against cycle 0, with W' 0.16-0.34 km/s. Every later
    comparison is <= 0.022 km/s.
  - The cause is the chain's first departure leg. It is the open chain's START and has no incoming
    junction constraining its |V_inf|. Only its magnitude closure at the far end of the chain ties it,
    so it is an end effect, not drift.
  - gc-1 does not show it, because its chain starts on a shot full-rev block.
  - A refinement for the owner: exclude cycle 0 (the unconstrained start), i.e. compare i = 1 ..
    n-L-1. Under that refinement gc-2 would pass 5/5 (W' <= 0.022). Not applied; for the owner's
    decision with the rest of sec. 8.4.
- The H&M positive control passes with margin. The ev rows are not computed (they need new 10-cycle
  chains). Nothing here changes any row.

### 8.6 EXPLORATORY: the ev rows under the draft metric (2026-10-11) — exploratory, under the draft metric, awaiting the owner

New 10-cycle Standish chains at the 5 standard epochs, same routes as each row's rung (d): ev-A
`--direct`; ev-C and ev-B ramp continuation. Data: `data/1007_v2/ev10_ev{A,B,C}/`.
Judgement: `scripts/v2_samephase_1007.py --ev-rows` (`data/1007_v2/samephase_ev_rows_exploratory.{json,log}`),
L = 5, the gate from the chain tool's minimax, re-fly not computed (exploratory).

| Row | Option (a) | Option (b) | What fails |
|---|---|---|---|
| ev-A | 2/5 | 2/5 | only the GROWTH test, at tiny absolute drift: D' 0.65-1.53 d against an 8.7-d band, W' <= 0.078 against a 0.49-km/s band. Growth ratios 2.03-6.0 are ratios of maxima of 2-3 values each |
| ev-B | 0/5 | 0/5 | at 10 cycles the gate fails or is indeterminate at 3 epochs (worst 1.011, 1.019, 0.9994); no 10-cycle chain at 1 epoch; D' 2.1-14.4 d against a 2.0-d band |
| ev-C | 0/5 | 0/5 | D' 7.6-65.8 d against about 5.8 d; W' 0.57-5.7 km/s. Growth ratios 0.8-1.9, so large but NOT growing |

Readings (exploratory; no row changes):
- ev-A's failures are the growth test's small-sample fragility. With 4-5 comparisons the "halves" hold
  2-3 values, and the drift itself is far inside the band.
- ev-C's 10-cycle chains differ from cycle i to i+5 by up to 66 d, with no growth. Its 5-cycle chains
  were not tested this way. A large constant same-phase offset suggests the 10-cycle continuation
  landed on a chain whose two 16-yr halves differ (another branch?). Not diagnosed.
- ev-B does not hold the gate over 10 cycles (48 yr).
- The H&M positive control passes (sec. 8.5).

### 8.4a OWNER DECISION (2026-10-11), recorded

The lead put the sec. 8.4 question to the owner with options (a) as drafted and (b) the start-cycle
refinement. The owner's answer, verbatim option label: **"Same-phase, skip start cycle"**, i.e.:
- sec. 8.4 with option (b), comparing i = 1 .. n-L-1;
- a pre-registered re-judge of gc-1 and gc-2 on their existing chains before the wording takes
  effect;
- the offset against the ideal model recorded, not banded.

### 8.7 Re-judge of gc-1 and gc-2 under the FINAL wording — PRE-REGISTERED (2026-10-11), before running

The final wording (owner, 8.4a) has four criteria:
1. junction-only |V_inf| continuity < 1e-6 km/s;
2. the demanded-turn gate passes at every flyby;
3. DOP853 re-fly < 1 km;
4. bounded SAME-PHASE drift, option (b):
   - L = 4 for gc k = 3;
   - compare cycle i+L with cycle i for i = 1 .. n-L-1 (5 comparisons on the 10-cycle chains);
   - D' <= 5 % of the chain's shortest Lambert leg and W' <= 10 % of its smallest |V_inf,dep|;
   - no growth (second-half max <= 2 x first-half max; floors 0.5 d, 0.05 km/s);
   - the ideal-model offset recorded, not banded.

A row passes at >= 3 of 5 epochs.

- Chains: the EXISTING ones (gc-1 `data/943_gc1_realeph/e0..e4`; gc-2 the ladder file); no re-solve.
- Tool: `scripts/v2_samephase_1007.py --rejudge-gc`.
  - Criteria 1-3 use `scripts/v2_chain_1007.py`'s `judge` (amended rules), as in sec. 7.
  - Criterion 4 uses `same_phase(start=1)`.
- Controls, first:
  - POSITIVE: H&M orbit 1 from its printed dates, 10 cycles (sec. 8.5 route), judged by the same four
    criteria (L = 5, i = 1..4). It must pass.
  - NEGATIVE: gc-1 epoch 0, the 2nd Lambert start + 0.05 d. It must fail.
- Expected (stated with knowledge of the sec. 8.5 exploratory option-(a) numbers, so not blind):
  - gc-1: 5/5, no longer marginal on drift. Criterion 1 stays at 8.6e-7, the same margin as before.
  - gc-2: 5/5. W' after the start cycle is <= 0.022 km/s.
- If a control misbehaves: STOP and report. If a row fails: it is demoted to V1 in the window, with the
  reason recorded.

### 8.8 Result of the sec. 8.7 re-judge (2026-10-11 11:05 AEDT): the POSITIVE CONTROL FAILS -> STOP

Run: `scripts/v2_samephase_1007.py --rejudge-gc` (`data/1007_v2/rejudge_gc_final.{json,log}`).

- POSITIVE, H&M orbit 1 (10 cycles, L = 5, comparisons i = 1..4): criteria 1-3 pass (continuity, gate,
  re-fly). Criterion 4 FAILS on the |V_inf| growth test alone:
  - W' per comparison 0.092 / 0.070 / 0.301 / 0.270 km/s (band 0.454);
  - second-half max / first-half max = 3.26 > 2, above the 0.05-km/s floor;
  - D' 1.66-2.09 d against a 7.52-d band, growth 0.88 (passes).
- NEGATIVE (gc-1 e0 + 0.05 d): fails continuity and re-fly, as required.
- Per sec. 8.7: STOP. gc-1 and gc-2 were NOT re-judged; their V2 status is unchanged and nothing goes
  to the window.

Diagnosis (for the owner; nothing changed):
- Under option (a) (with i = 0) the same chain gives W' 0.285 / 0.092 / 0.070 / 0.301 / 0.270. The
  first-half max is then 0.285, and the growth ratio is 1.06 (pass).
- The sequence rises, falls and rises again. That looks like oscillation of the same-phase difference
  itself, not monotonic growth, but 4-5 comparisons cannot tell the two apart.
- Skipping the start cycle removes the comparison that happened to carry the first large value, so the
  half-max ratio flips. This is the same small-sample fragility that failed ev-A at 2-3 values per half
  (sec. 8.6).
- So the growth test as worded (the ratio of the maxima of two halves) is not robust at n = 2L to
  2L + 1. Options for the owner:
  - (i) require more comparisons: n >= 3L (ev: 15 cycles, 48 yr), so each half holds >= 4-5 values;
  - (ii) replace the half-max ratio by a fitted linear trend of D'/W' over i, which must be
    insignificant against the band (e.g. projected change over the chain < 0.5 x band);
  - (iii) both.
  - Any choice is a further change after seeing results, so it must be pre-registered and recorded as
    such.
