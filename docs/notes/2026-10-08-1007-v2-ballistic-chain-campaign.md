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
