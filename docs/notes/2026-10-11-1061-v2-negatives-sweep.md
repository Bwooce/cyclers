# #1061 Sweep of past V2 negatives after the V2-ballistic drift-rule change (2026-10-11)

Read-only sweep. No re-judgement was run under this task. No catalogue, `validate.py` or test edit.

## 0. Answer to the owner: "does anything we ruled out in the past pass now we tuned the rule?"

Known now: the rule change has already rescued every case it can rescue as written, and all of
them were re-judged on 2026-10-11 in the `#1007` lane: the Hollister-Menning orbit 1 control, gc-1
(3/5 marginal to 5/5), gc-2 and ev-A (2/5 to 5/5). No other past negative in the repository failed
ONLY on a drift or growth reading inside the lane the new wording covers. Every older V2 negative
outside `#1007` failed on something the new wording does not touch: a V-infinity continuity or closure
residual, the demanded-turn gate (`#888`), Lambert infeasibility, a continuation stall, a genuinely
unstable orbit, a wrong topology, or a powered row judged by a ballistic test. Three items are open and
need a run, not a re-reading: (1) ev-C's 1/5 is recorded as "not diagnosed" in the `#1007` note
(secs. 8.6 and 8.14), so it may be a second branch in the 10-cycle continuation and not real drift;
a cheap comparison of stored chains would settle it. (2) The Russell chapter-4 Earth-Mars free-return
rows were declined V2 in June on a structural ground that `#830` made obsolete, and on a drift number
that `#830` showed was an artifact of the rounded period. They now fall inside the 2026-10-08 lane, but
no real-ephemeris chain exists for them, and `#830`'s blind real-ephemeris check did not converge for
any of them. (3) The Titania-Oberon four-body orbit (`#890`/`#895`) waits for the owner's ruling on
what V2 means for an orbit too unstable to fly three laps. One structural point blocks any Earth-Mars
re-judgement until the owner rules: under the `#1007` repeat definition (0.02 rev) Earth-Mars needs
L = 32-37 cycles, so 2L = 64-74 cycles (158-410 yr). At a looser 0.06 rev, L = 5-11 (the spec's "~7 laps").
The tolerance must be ruled on before an Earth-Mars same-phase test is meaningful.

## 1. Method

Sources swept (grep, then the cited section read in full):
- `docs/notes/*.md`, 118 notes that contain "V2", plus the terms drift, lap-over-lap, rotating-frame,
  breathing, "not promoted", "stays V1", "V2 fail";
- `data/OUTSTANDING.md`;
- `src/cyclerfinder/data/validate.py` (`_LEVEL_EVIDENCE` strings with "NOT V2");
- `data/catalogue.yaml` (read-only: level comments, notes, `data_gaps`);
- tests: `tests/verify`, `tests/search`, `tests/data` (V2 docstrings, xfail markers);
- `docs/spec.md` sec. 14 and its amendments of 2026-06-07 and 2026-10-08;
- the V2 instruments themselves: `verify/propagate.py`, `data/validation/v2_moontour.py`,
  `v2_3d.py`, `v2_bcr4bp.py`, `v2_qp.py`.

Classification rule (applied to the BINDING criterion, so the table can be audited):
- **(A)** The drift or growth reading was the only failing criterion, or the test window was shorter
  than the system's geometric repeat. Everything else (continuity, gate, re-fly, Lambert feasibility,
  continuation) passed or was not at issue.
  - **A-in**: inside the 2026-10-08 lane (multi-arc patched-conic row, exact synodic multiple in the
    ideal model, judged in the real-ephemeris chain). The new wording applies as written.
  - **A-out**: the same misreading is possible, but the case is outside that lane. The amendment says
    "For this class only", so a re-judgement needs an owner extension.
- **(B)** Any other criterion failed: continuity or closure, the demanded-turn gate, re-fly, Lambert
  feasibility, continuation to lambda = 1, model or topology mismatch, genuine dynamical instability,
  or the case never reached V2. A drift failure that comes WITH one of these is (B).
- **(C)** Powered rows (V2-powered, or a powered solution judged by a ballistic test). Out of scope.
- **(D)** Unclear. The table says what would settle it.
- **(N)** Never attempted, or declined before any run. Not a negative; listed in sec. 3 only.

## 2. Case table

| # | Row or object | Date | Source | Model and metric | Recorded reason (quoted) | Class |
|---|---|---|---|---|---|---|
| 1 | `aldrin-classic-em-k1-outbound` (powered solution) | 2026-06-07 | `tests/verify/test_aldrin_v2_v3_campaign.py` docstring (#134); spec sec. 14 amendment 2026-06-07 | DE440; original single V2 gate: rotating-frame repeat over 3 laps, 200,000 km | "~4.14e8 km ... Mars's heliocentric radius breathes because the 2.135 yr cycler period is not commensurate with Mars's 1.881 yr orbit" | **C**. The failure WAS a breathing misread. The June class split resolved it: the row passes V2-powered and is V2. Nothing to re-judge. |
| 2 | same row, V2-ballistic | 2026-06-07 | `docs/notes/2026-06-07-aldrin-promotion-assessment.md` sec. 2 | not run | "V2-ballistic verdict: NO ... the Aldrin row's defining model is circular-coplanar/powered" | **C** |
| 3 | `russell-ch4-5.30gGf3`, `-9.94Gg3`, `-5.75ggF3`, `-9.353Gg2` (the #137 free-return rows) | 2026-06-07 | `tests/search/test_free_return_v2_ballistic.py` docstring; `docs/notes/2026-08-12-829-833-835-830-russell12-bend-gate.md` secs. 5, 7 | circular ephemeris; `propagate.py` lap template on the rounded catalogue period; 50,000 km | "the V2-ballistic gate is structurally inapplicable to a single-arc slice"; measured "~9.4e7-1.2e8 km". `#830` sec. 7: "`test_free_return_v2_ballistic`'s cited drift is a period-rounding artifact" | **D** (top of the list, sec. 4). The recorded ground is obsolete (`#830`), the drift number is an artifact, and the rows now fall inside the 2026-10-08 lane, which has not been run for them. |
| 4 | the 7 admissible russell12 rows (adds `russell-ch4-3.78Gg3`, `-5.30ggF3`, `mcconaghy-2006-em-k2` (V0), `russell-ch4-4.991gG2` (S1L1, V3)) | 2026-08-12 | `#830` note sec. 5; OUTSTANDING `#830`, `#850` | circular, exact synodic period, lap template | "GATE-DEGENERATE ... `stable` here is a statement about PERIOD COMMENSURABILITY"; negative control with a 35 km/s break still `stable=True` | **D**, merged with case 3. Neither a pass nor a negative. |
| 5 | hand-built 2-encounter direct Earth-Mars cycler (test fixture, not a row) | 2026-06-01, re-anchored 2026-06-12 (#214) | `tests/verify/test_propagate.py::test_2syn_em_direct_cycler_does_not_close_over_3_laps_astropy` | astropy, 3 laps, rotating-frame repeat | "drifts 1e6-1e8 km ... The hand-built 2-encounter direct E-M cycler used here is NOT that topology" | **B** (wrong topology) |
| 6 | `russell-ocampo-4.3.1-5` | 2026-06-23 | `docs/notes/2026-06-23-388-russell-ocampo-4.3.1-5-anchor-recovery.md` | real ephemeris, single cycle | "V2 (≥3-cycle ballistic closure): NO — the close is single-cycle and 164 m/s low_maintenance, not ballistic" | **C** (powered, single cycle) |
| 7 | 12 Pluto repeated-moon candidates | 2026-06-15 | `docs/notes/2026-06-15-gauntlet-run-274.md` | circular-coplanar Lambert; "bar is 3.0 x the V1 residual" | table "V2 0 / 12", but every candidate "falls at V1" | **B** (never reached V2) |
| 8 | rows declined on instability, no run: `braik-ross-c11a/-c11b/-c32-cycler-2026`; `casoliva-1-2d`, `-2-1b`, `-7-3b`, `-7-3c`; `vaquero-21-c247`, `-21-c266`, `-31-c254`, `-31-c313`. Run and failed: the `#287` 3D C11a extension and one `#301` doubly-hyperbolic k = 4 member (candidates, not rows) | 2026-06-16 to 2026-08-12 | catalogue comment above `braik-ross-c11a-cycler-2026`; `validate.py` NOT-V2 strings; `docs/notes/2026-08-12-828-vaquero-connection-tier-adjudication.md`; `docs/notes/2026-06-16-306-3d-gauntlet-phase1.md` | CR3BP, `v2_3d.py` | rows: "NOT V2 (an UNSTABLE orbit cannot satisfy V2-ballistic bounded-drift-over->=3-laps)"; the two `#306` candidates: "Both correctly FAIL V2 under naive single-cycle propagation — hyperbolic-instability amplification dominates by cycle 2-3" | **B** (genuine instability in an autonomous model with no breathing). The owner's open question on V2 for unstable orbits (`#890` review item 6, `#907`) is separate. |
| 9 | `#392` low-amplitude (3,2) Earth-Moon branch | 2026-06-19 | `docs/notes/2026-06-19-392-floquet-low-amplitude.md` | CR3BP, uncorrected propagation, 50,000 km | 3 cycles 0.047 km PASS, 5 cycles 71.6 km PASS, "10 cycles ... 369,380 km max drift (FAIL)"; max Floquet 37.23; V4 "STRUCTURAL_FAIL_ALL_EPOCHS" | **B** (instability; V4 kills it independently) |
| 10 | `#327` Umbriel-Oberon-Umbriel SILVER | 2026-06-16; withdrawn by `#888` 2026-10-04 | `docs/notes/2026-06-16-330-moontour-v2-phase2.md`; `docs/notes/2026-10-04-888-turn-gate-wiring.md` sec. 2 | circular-coplanar moons; `v2_moontour` (inertial final-encounter offset against cycle 0) | "The drift OSCILLATES rather than growing monotonically"; "The closure residual grows monotonically (0.025 -> 0.349 km/s)"; `#888`: "all reject, with the turn named" | **B** (closure growth and turn gate). The drift part is the inertial reading (sec. 5.1). |
| 11 | `#344` Titan-Rhea-Titan SILVER | 2026-06-17 | `docs/notes/2026-06-17-344-phase2-stage-b-v2-moontour.md` | `v2_moontour` | "by cycle 1 the closure residual has jumped to 5.44 km/s" | **B** |
| 12 | 10 published moon-tour reproductions (5 Jovian, 5 Saturnian) | 2026-06-26; re-adjudicated 2026-06-26 (#473) | `docs/notes/2026-06-26-470-moontour-admission-wave-v2-verdict.md`; `docs/notes/2026-06-26-473-resonance-lock-moontour-search.md` step 2 | `v2_moontour` with the POWERED `MultiRevLeveragingReleg` backend; inertial drift, 50,000 km | "Every tour blows the 50,000 km inter-cycle rendezvous-drift floor by 3×–72×"; `#473` relative criterion: "2/10 flip" (IEG, EGC) to bounded quasi | **C** (powered). Note: the V2-powered (a)+(b) test was never applied to them, so the 0/10 is not a V2-powered negative either. The turn gate was wired after these runs and was not re-run on them. |
| 13 | `#566` five Uranian two-moon representatives (Titania-Oberon, Ariel-Umbriel, Ariel-Titania, Ariel-Oberon, Umbriel-Titania) | 2026-07-11; withdrawn by `#888` | `docs/notes/2026-07-11-566-opus-adjudication-gauntlet-results.md` secs. 1, 3; `#888` wiring note sec. 2 | `v2_moontour`, inertial drift | drift 2.9e5-8.7e5 km, closure ~1e-14; "V2 drift overflow is the class definition"; `#888`: "V2 (with drift and closure floors opened, so the turn is the only criterion left) ... all reject, with the turn named" | **B** (turn gate). Drift alone would have been A-out; the gate is binding. |
| 14 | `#574` Stage B, 15 Titan-Iapetus branches | 2026-07-12 | `docs/notes/2026-07-12-574-stageB-saturn-gauntlet.md` part 3; `#888` sec. 2 | `v2_moontour` | "per-cycle V_inf-continuity residual that grows to several km/s (0.3-6.4 km/s)"; 3 "Lambert infeasible past cycle 0"; `#888`: "All 17 of the 22 ... are turn-infeasible" | **B** |
| 15 | `#890` / `#895` Titania-Oberon-Titania periodic orbit (planar circular four-body model) | 2026-10-04 | `docs/notes/2026-10-04-890-titania-oberon-candidate.md` (results table); `docs/notes/2026-10-04-890-adversarial-review.md` item 10; `#888` sec. 2; OUTSTANDING `#890`, `#907` | DOP853 closure <= 1 km and 1 cm/s (pre-registered); also `v2_moontour` at the default floor | "V2 DOP853 closure <= 1 km and 1 cm/s: FAIL as pre-registered, 0.238 km, 1.64 cm/s"; review: "They were ill-posed"; `#888`: "At the default 50,000 km floor it fails on drift (the V2 drift is an inertial offset of the closing encounter, about 7e5 km here), not on the turn" | **D**. Settled by the owner's ruling on V2 for an orbit with Floquet multiplier 8.4e5 (`#890` review item 6, `#907`). Outside the 2026-10-08 lane. |
| 16 | Hollister 1969 orbit I (1H), ev positive control | 2026-10-08 | `docs/notes/2026-10-08-1007-v2-ballistic-chain-campaign.md` sec. 5 | Standish, ramp continuation | "The continuation stalls at every epoch" (lambda 0.29-0.39) | **B** (continuation). Control, not a row. |
| 17 | Hollister-Menning 1970 orbit 1, take 2 (own epoch) | 2026-10-11 | `#1007` sec. 8.1 | Standish, direct then ramp | "Direct solve ... no convergence"; "Ramp ... stalls at lambda 0.398" | **B** (continuation). Superseded by the printed-date route. |
| 18 | Hollister-Menning 1970 orbit 1, take 3 (printed dates) | 2026-10-11 | `#1007` sec. 8.3 | Standish, template D band and growth | "Criterion 4 FAILS: D per cycle 0 / 18.8 / 19.0 / 13.1 / 39.5 d against the 8.8-d band" | **A-in, already re-judged**: PASS under the final wording (`#1007` sec. 8.11). |
| 19 | `ganymede-callisto-two-working-body-cycler-gc1-2026` | 2026-10-08 | `#1007` secs. 6, 7 | jup365 chain, template D band | literal 0/5 (criterion 1 pre-registration error), then "3/5 PASS ... MARGINAL ... Epochs 3-4 fail the D band (0.222 / 0.195 d against 0.190)" | **A-in, already re-judged**: 5/5 (sec. 8.12). |
| 20 | `ganymede-callisto-two-working-body-cycler-gc2-2026` | 2026-10-08 | `#1007` secs. 6, 7, 8.5 | jup365 chain, W band; then draft same-phase with i = 0 | "criterion 4 at 3 epochs: the W band"; draft: "Every gc-2 W' failure is at i = 0" | **A-in, already re-judged**: 5/5 (secs. 7 and 8.12). |
| 21 | `earth-venus-two-working-body-cycler-eva-2026` (ev-A) | 2026-10-11 | `#1007` secs. 8.6, 8.14 | Standish, draft same-phase with growth test | "only the GROWTH test, at tiny absolute drift" (2/5) | **A-in, already re-judged**: 5/5 (sec. 8.14). |
| 22 | `earth-venus-two-working-body-cycler-evb-2026` (ev-B) | 2026-10-11 | `#1007` secs. 8.6, 8.14 | Standish 10 cycles, final wording | "gate fails at 3 epochs (1.011, 1.019, 0.999); no 10-cycle chain at 1" | **B** (gate and continuation). Footnote on L in sec. 5.3; it does not change the verdict. |
| 23 | `earth-venus-venus-hosted-cycler-evc-2026` (ev-C) | 2026-10-11 | `#1007` secs. 8.6, 8.14 | Standish 10 cycles (ramp continuation), final wording | 1/5: "D' 15-66 d against about 5.8"; sec. 8.6: "Growth ratios 0.8-1.9, so large but NOT growing ... suggests the 10-cycle continuation landed on a chain whose two 16-yr halves differ (another branch?). Not diagnosed."; sec. 8.14: "Its 5-cycle rung (d) chains passed (6.23) ... Not diagnosed further (another branch at 10 cycles is possible)" | **D**. The source does not support "real drift" as settled. See sec. 4 item 2 for the check. |
| 24 | R-S GanCal#1 at 2013 (reported control) | 2026-10-08 | `#1007` secs. 3, 6 | jup365 chain | gate "indeterminate" 0.992 "by construction, as in R-S's own ideal model" | **B** (gate). Reported, not judged. |

Not cases (recorded so they are not re-counted): the `#330` and `#306` "NOT_RUN_PHASE_1" verdicts; the
`#830` broken-chain negative control; the BCR4BP `non_commensurate_no_strict_period` result in
`tests/data/test_v2_bcr4bp.py`, which is a test fixture (no candidate verdict carries it; checked by
grep of `src scripts tests data docs`).

## 3. Never attempted or declined before any run (N): not negatives

41 rows at V1 have no V2 attempt in the repository. Their `_LEVEL_EVIDENCE` or catalogue comment says so
("NOT V2: no real-ephemeris chain has been run", "no multi-lap bounded-drift run has been performed"):
- two-working-body Ganymede-Callisto members gc4-* (9), gc5-* (3), gc6-* (3): no real-ephemeris chain;
- `russell-ocampo-*` V1 rows (7); `russell-ch4-3.64gGg3`, `russell-ch4-6.44Gg3`;
  `aldrin-classic-em-k1-inbound`; `liang-2024-cgcec-*` (3);
- CR3BP: `casoliva-1-2c`, `-1-2e`, `-2-1a`, `-3-2c`, `-7-3a`; `vaquero-21-c198`, `-21-c246` (`#823`);
  `pc-cycler-51-2026` (`#846`); `ross-rt-mu*` (4); `europa-3-4-crnbp-torus-jupiter-2026`.

Declined by class, not by a run: the Tito-MacCallum 2018 free return ("It cannot earn even V2", a
one-shot mission, no row), the Restrepo-Russell non-cycler families (not ingested), and the `#287` 3D
Aldrin scoping spike (predicted fail, not run).

The 240 rows with no `validation_level` and the V0 rows are outside this sweep: none has a V2 attempt.

## 4. Re-judgement list ((A) and (D) cases), priority order

No (A) case is outstanding: all four A-in cases (18-21) were re-judged on 2026-10-11, and no A-out
case exists (every drift-only reading outside the lane also fails another criterion). The list is therefore the (D) cases.

1. **ev-C (case 23).** Cheapest and closest to a verdict. Check, on data already on disk: compare the
   first five cycles of each stored 10-cycle chain (`data/1007_v2/ev10_evC/`) with the stored 5-cycle
   rung (d) chain at the same epoch (`data/942_direct_route_and_gm_fix.json` key `evC_standish_rerun`).
   If the halves differ (a second branch), re-solve 10 cycles by direct LM seeded from the 5-cycle chain
   repeated one 16-yr period later (the H&M control route of `#1007` sec. 8.5), then re-judge under the
   final wording with the same controls. Pre-register before running. Minutes of compute.
2. **The Russell chapter-4 Earth-Mars free-return rows (cases 3 and 4).** In the 2026-10-08 lane by
   class. Needs: (a) the owner's ruling on the Earth-Mars repeat tolerance (sec. 5.3; at 0.02 rev the
   lane needs 64-74 cycles); (b) a ≥ 3-cycle real-ephemeris patched-conic chain with junction continuity,
   the gate at every flyby including the wrap, re-fly < 1 km, at ≥ 3 of 5 epochs; (c) an Earth-Mars
   positive control (none is named yet). Expected blocker: `#830` sec. 6 ran `#388`'s canonical
   real-ephemeris lane blind and no row converged (the two high-V-infinity rows collapse to about
   19 km/s at Mars). If that holds in the chain tool, the result is (B) on continuation, not a rescue.
   Multi-session work.
3. **Titania-Oberon (case 15).** Waits for the owner's ruling on V2 for an orbit too unstable to fly
   three laps (`#890` review item 6, `#907`). No run until then.

## 5. The original V2-ballistic gate and the geometric repeat (task part 3)

Spec text relied on (sec. 14 table): "≥3 continuous laps; **bounded** drift in the dynamic rotating
frame (tolerant of geometric breathing), evaluated **in the row's defining model** (for a
circular-coplanar row that is the idealized propagation ...)". Sec. 14 note: "≥3 laps = two intervals ...
The full geometric-modulation horizon (~7 laps / ~15 yr for an E–M cycler ...) is **V3's** burden, not
V2's."

### 5.1 Where breathing can occur in the original gate

- **Inside the defining model, there is no breathing for most classes.**
  - CR3BP rows (Earth-Moon, Pluto-Charon, Neptune-Triton, Saturn-moon CR3BP): the model is autonomous,
    so the repeat is one orbit period (L = 1). A 3-lap window cannot misread breathing. Growth there is
    real instability (cases 8, 9).
  - Circular-coplanar rows with a period that is an exact synodic multiple of a body PAIR: the
    configuration repeats exactly every cycle in the rotating frame (L = 1). This is the `#830`
    degeneracy: the gate cannot fail.
  - Circular-coplanar rows with three or more bodies (VEM, Io-Europa-Ganymede, Callisto-Ganymede-Europa):
    the configuration repeats only when all pair periods are commensurate. In-model L can exceed 3.
- **Breathing appears when V2 is run in an eccentric or real-ephemeris model**: the `#134` Aldrin run
  on DE440 (case 1) and the `#1007` lane (cases 18-23).
- **The moontour instrument has its own repeat problem, in the ideal model.** `v2_moontour.py` measures
  `||r_final_k - r_final_0||`, the planet-frame position of the closing encounter against cycle 0
  (docstring step 2c). That is an inertial comparison, not the spec's rotating frame. Unless the cycle
  period is a multiple of the anchor moon's own period, the closing encounter moves around the planet
  each cycle and the "drift" oscillates up to twice the moon's orbit radius. `#888` sec. 2 says so: "the
  V2 drift is an inertial offset of the closing encounter". `#473` met the same effect ("invariant under
  a whole-tour rotation") and moved to a relative criterion. No case in sec. 2 flips because of it: every
  moontour negative also fails closure, the turn gate, or is powered. It is an instrument defect for any
  future moontour V2 claim (outside the 2026-10-08 lane, so an owner matter).

### 5.2 Repeat L per system class (real-ephemeris breathing)

Definition (`#1007` sec. 8.4): L is the smallest integer such that, for every encounter body,
L P / T_body is within TOL rev of an integer (T = the body's real mean sidereal period). Computed for this
note with:
- periods (d): Earth 365.256363, Venus 224.7008, Mars 686.979586, Moon anomalistic 27.55455 and synodic
  29.530589 (the Sun's direction in the Earth-Moon frame), Io 1.769138, Europa 3.551181, Ganymede
  7.154553, Callisto 16.689018, Enceladus 1.370218, Rhea 4.518212, Titan 15.945421, Iapetus 79.3215,
  Umbriel 4.144177, Titania 8.705872, Oberon 13.463234, Triton 5.876854, Charon 6.38723;
- P = k x the real synodic period of the pair, or the row's stated period;
- breathing amplitude about e x a of the body's own orbit.

The "3-lap window" column asks whether laps 0..2 contain a same-phase pair, i.e. L <= 2.

| System class | P | L at 0.02 rev | L at 0.06 rev | L P (0.02) | e a of the bodies (km) | 3-lap window covers the repeat? | Breathing misread possible in a real-ephemeris V2? |
|---|---|---|---|---|---|---|---|
| Earth-Mars k = 1 | 779.9 d | 37 | 7 | 79 yr | Earth 2.5e6, Mars 2.1e7 | no | yes, large (case 1 is the instance) |
| Earth-Mars k = 2 | 1559.9 d | 37 | 11 | 158 yr | as above | no | yes |
| Earth-Mars k = 3 | 2339.8 d | 32 | 5 | 205 yr | as above | no | yes |
| Earth-Venus k = 1 | 583.9 d | 5 | 5 | 8.0 yr | Earth 2.5e6, Venus 7.4e5 | no | yes |
| Earth-Venus k = 2 | 1167.8 d | 5 | 5 | 16.0 yr | as above | no | yes (case 18 is the instance) |
| Earth-Venus k = 3 | 1751.8 d | 49 (5 fails by 1e-4) | 5 | 235 yr | as above | no | yes |
| Venus-Earth-Mars (three bodies) | 779.9 d | none <= 400 | 155 | - | as above | no | yes, even in the circular model |
| Earth-Moon, real ephemeris (lunar anomaly and Sun) | 44.8-84.5 d | 322 or none <= 400 | 29-176 | decades | Moon 2.1e4, plus solar perturbation | no | applies to V3/V4 only: the Earth-Moon rows' V2 is in the CR3BP (L = 1) |
| Ganymede-Callisto k = 3 / 4 / 5 / 6 | 37.6 / 50.1 / 62.6 / 75.1 d | 4 / 1 / 4 / 2 | same | 150 / 50 / 251 / 150 d | Ganymede 1.4e3, Callisto 1.4e4 | k = 4, 6 yes; k = 3, 5 no | yes, small (gc-1's sec. 6 D-band misses) |
| Europa-Ganymede k = 1; Io-Europa-Ganymede | 7.05 d | 1 | 1 | 7 d | Europa 6.3e3 | yes | no |
| Ganymede-Io k = 1 | 2.35 d | 3 | 3 | 7 d | Io 1.7e3 | no | in principle; amplitude small |
| Callisto-Ganymede-Europa (Liang, about 100 d) | 100 d | none <= 400 | none | - | up to 1.4e4 | no | in principle |
| Titan-Enceladus k = 1 / 2 / 3 | 1.5 / 3.0 / 4.5 d | 32 / 16 / 39 | 10 / 5 / 7 | 48-175 d | Titan 3.5e4 | no | yes: Titan's e a is close to the 50,000 km V2 floors |
| Titan-Rhea; Titan-Iapetus | 6.3 d; 20.0 d | 43; 4 | 5; 4 | 271 d; 80 d | Titan 3.5e4, Iapetus 1.0e5 | no | yes |
| Uranian pairs (Umbriel-Oberon k = 1; SILVER 29.88 d; Titania-Oberon 123.2 d) | 6.0 / 29.9 / 123.2 d | 9 / none / none | 9 / 114 / 7 | - | 2e2 to 1e3 | no | small; the moontour inertial offset (sec. 5.1) dominates |
| Neptune-Triton | one Triton period | 1 | 1 | 5.9 d | Triton 6 | yes | no |
| Pluto-Charon (3,2) | 12.03 d | 26 | 9 | 313 d | Charon 1 | no | negligible (1 km) |

Readings:
- For every heliocentric class the 3-lap window is shorter than the repeat. A real-ephemeris V2 with a
  no-growth or template test can misread breathing there; Mars's eccentricity makes Earth-Mars the worst.
- Moon systems breathe little, except Titan (and Iapetus). The 50,000 km moontour floor sits near Titan's
  e a of 35,000 km, so a Saturnian real-ephemeris drift test could misread Titan's breathing.
- The CR3BP rows' original V2 cannot misread breathing (L = 1 in-model).

### 5.3 Two findings about L

1. **L depends strongly on the tolerance for Earth-Mars.** 37 at 0.02 rev; 7 at 0.06 (the spec's "~7 laps
   / ~15 yr" is 7 synodic periods, 0.053 rev off for both planets). Under the `#1007` definition as
   written, a same-phase Earth-Mars V2 needs 2L = 64-74 cycles (k = 3 to k = 1). The owner must rule on the tolerance before
   any Earth-Mars re-judgement (sec. 4 item 2).
2. **Earth-Venus k = 3 (ev-B) sits on the tolerance edge.** At L = 5 the residual is 0.0201 rev for both
   planets, so a strict 0.02 reading gives L = 49, not the L = 5 of the `#1007` sec. 8.4 table. L is
   hard-coded in `scripts/v2_samephase_1007.py`, taken from that table; it is not computed in code. This
   does not change ev-B's verdict, which fails on the gate at 3 epochs and on continuation at 1.

## 6. Scope limits of this sweep

- Grounding is by quoted source text; no number was recomputed except the L table (sec. 5.2), whose
  constants are listed above.
- V0 rows and rows with no `validation_level` were checked only for a V2 attempt (none found).
- The catalogue comments for ev-A, ev-B and ev-C still read "NOT V2" because the `#1007` window has not
  landed yet; this sweep used the `#1007` note as the source of truth for them.
