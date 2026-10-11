# #1025: real-ephemeris rung (d) and #415 dV bands for the 22 clean #973 Earth-Venus members

Status: DONE: Standish, DE440 (amendment A1) and the deferred DE440 shoot (sec. 8). Pre-registration secs. 1-5 (commit 4e18f740);
amendment A1 sec. 5a (64c2736f, before any DE440 result); results sec. 6; the deferred DE440 shoot
launch sec. 7.
Owner decision of 2026-10-08 (relayed by the lead on 2026-10-11): rung (d) before any row decision. No
catalogue rows are written, and nothing is called novel.

Members: ev k = 4: ev4-0, 1, 3, 4, 5, 6, 7, 8, 9, 10, 11; ev k = 5: ev5-0 ... ev5-10. Keys, dates and
values are in `docs/notes/2026-10-07-973-two-working-body-k-extension.md` secs. 13-14 and
`data/973_ev/k4_gauntlet.json`, `k5_gauntlet.json`.

## 1. Tools (shared; not edited)

- Chain tool: `scripts/run_942_realeph_chain.py`, last changed in 6ac27fb5 (the 6.55a `--grow-start`).
  The #1007 work did not change it.
- Re-fly checker: `scripts/check_942_realeph_chain.py`.
- Closure dV: `scripts/analyse_942_fullrev_closure_dv.py` (6.24/6.24a).
- Turn-deficit prices: `scripts/dv_band_942.py::price_chain` and `row`, imported.
- Orchestrator: `scripts/ev_realeph_1025.py` (new). It calls the tools as subprocesses, one epoch per
  call under a `timeout`, and writes `data/1025_ev_realeph/<member>/<route>/e<i>/` plus a static
  runlog.

Models:
- Standish: Standish & Williams J2000 mean elements, ramp mode. At lambda = 1 the full-rev legs are
  exact.
- DE440: astropy DE440, blend + shoot. The DE440 "E" is the Earth, not the Earth-Moon barycentre
  (`core/ephemeris.py`, astropy body "earth"). So the 6.25 lunar-reflex issue applies to every
  member: all 22 have Earth full-rev legs. The Earth legs' closure dV measures the 1:1 seed rule on a
  non-Keplerian Earth and is an upper bound. Earth and Venus legs are reported separately.

Chain length, about 16 yr as for ev-A:
- k = 4: 3 cycles (19.19 yr). 2.5 cycles is not a chain.
- k = 5: 2 cycles (15.99 yr).

Epochs: the tool's 5 standard epochs (first 2030-01-01, span 32 yr, phase-matched within 1.1 synodic
periods of each nominal date), as for ev-A, ev-B, ev-C and 1H.

## 2. Route ladder per epoch (fixed now; no route is added after any result is seen)

Standish:
- S1 `--real mean --direct` (6.22): a date solve at lambda = 1 from the ideal dates.
- S2, only if S1 does not PASS: `--real mean` (6.55a step 1), the ramp-lambda continuation at the full
  length.
- S3, only if S2 stalls before lambda = 1 and the length is >= 3 (k = 4 only):
  `--real mean --grow-chain --grow-start 2` (6.55a step 2).

DE440:
- D1 `--real de440 --shoot-jac sparse` (6.18/6.25): blend continuation, then at lambda = 1 the shoot
  of every fixed leg (20 restarts, perturbations 0.03/0.05 rad).
- D2, only if D1 does not reach lambda = 1: `--real de440 --direct --shoot-jac sparse` (date solve at
  lambda = 1, then the shoot).

PASS at an epoch: a converged chain at lambda = 1 of the full length, with every interior flyby
passing the #888 gate at the registry floors ("indeterminate" is not a pass), and the checker's DOP853
re-fly missing by < 1 km. On Standish the V_inf vector error must also be < 1e-6 km/s (the 6.55
criterion). On DE440 the chain must be a ballistic shoot closure. A pass by any route of the ladder
counts; the route is recorded, and the note states that up to three attempts were made.

Each epoch is reported as one of:
- PASS;
- FAIL: a converged lambda = 1 chain whose best worst ratio is > 1, or "indeterminate";
- NOT REACHED: every route folds, does not converge, or finds no shoot closure (6.18 finding 2: a
  fold is NOT a negative);
- NOT RUN: a timeout (420 s per epoch and route). A timeout is never a fail.

The table gives, for example, "Standish 1 pass / 2 fail / 2 not reached".

## 3. Controls first (stop rule)

Through today's tool, same settings:
- C1, Hollister 1H, route S1, 5 cycles. Recorded (6.23 D1): pass at the 2nd and 3rd epochs
  (JD 2465314.3, 2467656.4).
- C2, ev-A, route S1, 5 cycles. Recorded: pass 5/5.
- C3, ev-B epoch 0, 6.55a step 1, 4 cycles. Recorded: pass, worst 0.933.

If C1 passes at no epoch, the route cannot judge the 22, and the campaign stops and reports. If any
control differs from its record, the campaign also stops and reports.

## 4. dV band (#415, `dv_band_acceptance`)

For each member, model and epoch with a converged lambda = 1 full-length chain (the passing chain if
there is one, otherwise the converged chain with the lowest worst ratio):
- Turn-deficit price at every PRICED flyby (`price_chain`; beyond-bend and Oberth-credited
  periapsis), with the minimax directions. Prices on gate-failing chains are therefore upper-bound
  style (6.57).
- DE440 adds the chain's closure dV (sum of the full-rev mid-course corrections, 6.24a), Earth and
  Venus legs separately.
- Per 7 cycles = per-chain total x 7/N (6.57 basis). Known bias, stated rather than fixed:
  `price_chain` drops the wrap block (it is not gated either). That is 1 of 4 blocks at N = 2 (k = 5)
  and 1 of 6 at N = 3 (k = 4), against 1 of 10 for ev-A. The number of priced flybys per chain is
  reported beside each figure.
- The band uses the beyond-bend price (the periapsis price is also reported).

Reported per member and model:
- every epoch's per-7 value ("n/a" where no lambda = 1 chain exists; never 0);
- the best epoch (the minimum over epochs and routes, which is optimistic);
- the range across epochs, as the window dependence.

## 5. Expected outcomes (written before any member is run)

Anchors: ev-A (ideal worst 0.574) passed Standish at 5/5; ev-B (0.673) at 2/5; 1H at 2/5. ev-A's DE440
closures all failed the gate, and its closure dV was 2-4 m/s per 7 cycles.

| Skeleton (E / V km/s) | Members | Ideal worst ratio | Standish expectation | DE440 expectation |
|---|---|---|---|---|
| 4.462 / 5.197 | ev4-0 | 0.858 | at most 1 pass; mostly fail / not reached | no gate-passing ballistic closure; many not reached |
| 4.595 / 4.529 | ev4-1 | 0.835 | the same | the same |
| 6.806 / 4.017 | ev4-3 | 0.990 | 0 passes likely (no margin) | the same |
| 6.830 / 4.038 | ev4-4, 5, 6, 7 | 0.779-0.840 | at most 1 pass each | the same |
| 6.893 / 3.968 | ev4-8, 9, 10, 11 | 0.939 | 0-1 passes | the same |
| 5.227 / 6.284 | ev5-0, 1 | 0.931 | 0-1 passes | the same |
| 8.627 / 10.788 | ev5-2 | 0.915 | 0-1 passes | the same |
| 10.300 / 5.546 | ev5-3 | 0.622 | 1-3 passes (the ev-B-like margin) | the same |
| 11.699 / 6.339 | ev5-4, 5 | 0.801 | 0-2 passes | the same |
| 13.184 / 7.251 | ev5-6, 7 | 0.622 | 1-3 passes | the same |
| 14.253 / 7.946 | ev5-8 | 0.855 | 0-2 passes | the same |
| 15.068 / 8.497 | ev5-9 | 0.309 | the best chance: 2-5 passes | the same |
| 34.174 / 42.154 | ev5-10 (sun-grazer) | 0.944 | 0-1 passes | the same |

These are expectations, not criteria; members of one skeleton differ in their return blocks and are
judged separately.

Pilots disclosed (scratch, before this note): ev4-0 at epoch 0, 3 cycles.
- S1: converged, gate FAIL, worst 2.727.
- S2: folds near lambda 0.19.
- D1: blend folds near lambda 0.26.
- D2: does not converge.
These are not results. The pre-registered run repeats them under this protocol.

## 5a. AMENDMENT A1 (lead ruling (c), 2026-10-11), committed before any DE440 result is read

Why: the pre-registered D1 shoot costs about 35 s per restart, so about 12 min per epoch and route.
- Measured on ev4-1 epoch 0, the only DE440 call made: it timed out at 420 s after 15 restarts,
  none closing. Its log shows the blend reaching lambda 1 with the minimax gate passing (0.819).
  That is disclosed here. It is kept as "NOT RUN" under D1 and is not used as a result.
- The #942 amendment 6.18 also gives full-rev candidates no DE440 pass/fail, only the closure dV
  (6.24).

DE440 now follows the 6.18/6.24 reading for all 22 members. There is no DE440 pass column.
- Route DB1: `--real de440 --shoot-restarts 0`, the blend continuation to lambda 1 with no shoot.
- Route DB2, only if DB1 does not reach lambda 1: `--real de440 --direct --shoot-restarts 0`.
- Per epoch: whether lambda 1 was reached (y/n); the #888 gate at lambda 1 with the minimax
  directions (status, worst ratio).
- The closure dV of the full-rev legs (`analyse_942_fullrev_closure_dv.py`, 6.24a), Earth legs and
  Venus legs separately. The Earth figures are an UPPER BOUND, because DE440's Earth carries the lunar
  reflex (6.25).
- The turn-deficit price at the minimax directions (`price_chain` on DE440), non-zero only where that
  gate fails.
- Per 7 cycles: (closure + deficit) x 7/N, with the band; plus the band of the Venus legs alone (the
  part the lunar reflex does not inflate).
- Summary per member: "DE440: reached r/5, minimax gate pass g/5", the best-epoch band, and the range.

The D1 shoot (a pass/fail by ballistic closure) is deferred to a lead launch, for the members with
Standish >= 4/5 only: ev4-1, 4, 6, 7; ev5-4, 5, 6, 7, 8, 9 (sec. 7 has the command).

## 6. Results (2026-10-11)

Data: `data/1025_ev_realeph/`:
- `<member>/<route>/e<i>/` (chain, tool log, checker log; DE440 closure_dv.json);
- `<member>/de440_blend.json`;
- `summary.json`;
- the static `runlog.txt`.

Controls (sec. 3), all reproduce their record:
- C1 (Hollister 1H, direct): pass at epochs 1 and 2, worst 0.862 / 0.822, re-fly <= 1.4e-4 km.
- C2 (ev-A): 5/5, worst 0.580, re-fly <= 1.9e-3 km.
- C3 (ev-B epoch 0, 6.55a step 1): pass, worst 0.933, re-fly 0.017 km.

Instrument check (dV): the worst ratio re-derived by `price_chain` equals the chain tool's to 3e-12 at
every priced epoch.

Standish (pass / fail / not reached):
- 5 / 0 / 0: ev4-1, 4, 6, 7; ev5-9.
- 4 / 0 / 1: ev5-4, 5, 6, 7, 8.
- 3 passes: ev5-2 (1 fail, 1 not reached), ev5-3 (2 not reached).
- 2 / 2 / 1: ev4-5.
- 1 / 4 / 0: ev5-10 (the sun-grazer).
- No pass: ev4-0, 8, 9, 10, 11 (0 / 2 / 3); ev5-0, 1 (0 / 1 / 4); ev4-3 (0 / 0 / 5).
- Every pass re-flies (DOP853 against the mean-element system) with a miss below 1e-3 km and a V_inf
  vector error below 1e-6 km/s.
- Against the expectations of sec. 5: the k = 4 skeletons 4.595 / 4.529 (ev4-1) and 6.830 / 4.038
  (ev4-4, 6, 7; ev4-5 2/5) did far better than expected (5/5, against "at most 1 pass").
  - The binding Venus ratio falls on the real eccentricities: ev4-4/6/7 run at 0.774-0.981 against
    the ideal 0.779, and ev4-1 at 0.775-0.826 against 0.835.
  - The 6.893 / 3.968 skeleton (ev4-8 to 11) and 6.806 / 4.017 (ev4-3) got 0, as expected. So did
    ev5-0/1.
  - ev5-9 (ideal 0.309) is 5/5, as expected.

Standish dV band:
- A passing Standish chain closes exactly and has no turn deficit, so its band is "strictly
  ballistic": the figures are about 1e-8 m/s, numerical zero. That is true by construction, the same
  weak control as 6.57.
- The gate-failing chains need 5,800 to 68,000 m/s per 7 cycles of turn-deficit impulse at the
  minimax directions (upper-bound style).

DE440 (amendment A1: blend to lambda 1, no shoot):
- The closure dV is the mid-course correction that puts each full-rev leg back on its planet.
- The Earth legs dominate it. Earth-leg misses are 11,000 to 1.34 million km (Venus legs 1,200-34,000 km): the
  lunar reflex of DE440's Earth against the 1:1 seed rule (6.25). So the "total" column is an UPPER
  BOUND.
- The Venus legs plus the turn deficit are the part the reflex does not inflate.
- At the best epoch, the Venus part is essentially ballistic (< 10 m/s per 7 cycles) for ev4-4, 5, 7;
  ev5-3, 9 and 10. It is low-maintenance (< 300) for the other Standish passers.
- ev5-9 has no Venus full-rev leg; its whole closure dV is Earth's.
- Window dependence: for one member the per-7 total varies by 1.2x to about 10x across the reached epochs
  where the minimax gate passes, and by far more where an epoch's minimax gate fails (column 7).

Table (all 22; best epoch = the minimum over the epochs and routes, which is optimistic):

| Member | Standish pass / fail / not reached (of 5) | Standish passing epochs (worst ratio range) | Standish band, best epoch (per 7 cycles, beyond-bend) | DE440 reached / minimax gate pass (of 5) | DE440 per 7 cycles, best epoch: total (upper bound) / Venus legs + deficit | DE440 range of total per 7 cycles |
|---|---|---|---|---|---|---|
| ev4-0 | 0 / 2 / 3 | - | above-low-maintenance, e0 (4.86e+04 m/s) | 1 / 0 | e4: 67627.4 (above-low-maintenance) / 67346.8 (above-low-maintenance) | 67627.4-67627.4 m/s |
| ev4-1 | 5 / 0 / 0 | e0, e1, e2, e3, e4 (0.775-0.826) | strictly-ballistic, e1 (3.74e-09 m/s) | 5 / 5 | e3: 308.4 (above-low-maintenance) / 23.1 (low-maintenance) | 308.4-453.1 m/s |
| ev4-3 | 0 / 0 / 5 | - | n/a | 0 / 0 | n/a | n/a |
| ev4-4 | 5 / 0 / 0 | e0, e1, e2, e3, e4 (0.774-0.981) | strictly-ballistic, e0 (5.03e-09 m/s) | 5 / 5 | e2: 43.3 (low-maintenance) / 5.3 (essentially-ballistic) | 43.3-51.1 m/s |
| ev4-5 | 2 / 2 / 1 | e1, e4 (0.863-0.921) | strictly-ballistic, e4 (1.39e-08 m/s) | 4 / 2 | e4: 24.3 (low-maintenance) / 9.4 (essentially-ballistic) | 24.3-11486.7 m/s |
| ev4-6 | 5 / 0 / 0 | e0, e1, e2, e3, e4 (0.774-0.981) | strictly-ballistic, e0 (5.02e-09 m/s) | 5 / 5 | e0: 48.0 (low-maintenance) / 10.7 (low-maintenance) | 48.0-473.9 m/s |
| ev4-7 | 5 / 0 / 0 | e0, e1, e2, e3, e4 (0.774-0.981) | strictly-ballistic, e0 (5.02e-09 m/s) | 5 / 5 | e3: 225.4 (low-maintenance) / 8.3 (essentially-ballistic) | 225.4-338.8 m/s |
| ev4-8 | 0 / 2 / 3 | - | above-low-maintenance, e4 (5.78e+03 m/s) | 2 / 0 | e4: 5906.7 (above-low-maintenance) / 5864.4 (above-low-maintenance) | 5906.7-6291.7 m/s |
| ev4-9 | 0 / 2 / 3 | - | above-low-maintenance, e4 (5.78e+03 m/s) | 2 / 0 | e4: 5887.2 (above-low-maintenance) / 5864.4 (above-low-maintenance) | 5887.2-6269.8 m/s |
| ev4-10 | 0 / 2 / 3 | - | above-low-maintenance, e4 (5.78e+03 m/s) | 2 / 0 | e4: 6029.0 (above-low-maintenance) / 5864.4 (above-low-maintenance) | 6029.0-6442.9 m/s |
| ev4-11 | 0 / 2 / 3 | - | above-low-maintenance, e4 (5.78e+03 m/s) | 2 / 0 | e4: 6114.5 (above-low-maintenance) / 5864.4 (above-low-maintenance) | 6114.5-6418.0 m/s |
| ev5-0 | 0 / 1 / 4 | - | above-low-maintenance, e2 (7.14e+03 m/s) | 1 / 0 | e2: 7157.6 (above-low-maintenance) / 7114.2 (above-low-maintenance) | 7157.6-7157.6 m/s |
| ev5-1 | 0 / 1 / 4 | - | above-low-maintenance, e2 (7.14e+03 m/s) | 1 / 0 | e2: 7471.8 (above-low-maintenance) / 7114.2 (above-low-maintenance) | 7471.8-7471.8 m/s |
| ev5-2 | 3 / 1 / 1 | e0, e2, e3 (0.784-0.869) | strictly-ballistic, e0 (1.25e-08 m/s) | 4 / 3 | e3: 83.1 (low-maintenance) / 14.1 (low-maintenance) | 83.1-36320.9 m/s |
| ev5-3 | 3 / 0 / 2 | e0, e1, e4 (0.548-0.567) | strictly-ballistic, e0 (2.5e-09 m/s) | 4 / 4 | e3: 16.6 (low-maintenance) / 7.2 (essentially-ballistic) | 16.6-83.3 m/s |
| ev5-4 | 4 / 0 / 1 | e0, e2, e3, e4 (0.585-0.742) | strictly-ballistic, e3 (3.25e-09 m/s) | 4 / 4 | e0: 42.0 (low-maintenance) / 27.7 (low-maintenance) | 42.0-256.3 m/s |
| ev5-5 | 4 / 0 / 1 | e0, e2, e3, e4 (0.755-0.780) | strictly-ballistic, e3 (3.22e-09 m/s) | 4 / 4 | e2: 160.7 (low-maintenance) / 66.2 (low-maintenance) | 160.7-550.8 m/s |
| ev5-6 | 4 / 0 / 1 | e0, e2, e3, e4 (0.405-0.561) | strictly-ballistic, e0 (7.79e-09 m/s) | 4 / 4 | e1: 37.4 (low-maintenance) / 15.1 (low-maintenance) | 37.4-71.0 m/s |
| ev5-7 | 4 / 0 / 1 | e0, e2, e3, e4 (0.533-0.561) | strictly-ballistic, e0 (6.99e-09 m/s) | 4 / 4 | e4: 299.0 (low-maintenance) / 53.7 (low-maintenance) | 299.0-403.5 m/s |
| ev5-8 | 4 / 0 / 1 | e0, e1, e3, e4 (0.602-0.797) | strictly-ballistic, e0 (2.72e-09 m/s) | 5 / 5 | e4: 57.6 (low-maintenance) / 44.8 (low-maintenance) | 57.6-180.0 m/s |
| ev5-9 | 5 / 0 / 0 | e0, e1, e2, e3, e4 (0.304-0.318) | strictly-ballistic, e0 (5.31e-09 m/s) | 5 / 5 | e1: 32.8 (low-maintenance) / 0.0 (strictly-ballistic) | 32.8-41.3 m/s |
| ev5-10 | 1 / 4 / 0 | e0 (0.940-0.940) | strictly-ballistic, e0 (6.77e-08 m/s) | 5 / 1 | e0: 24.6 (low-maintenance) / 1.0 (essentially-ballistic) | 24.6-24816.9 m/s |

No row is written; nothing is called novel.

## 7. Deferred: the DE440 shoot (lead launch; ruling (c), part (a))

Members (Standish >= 4/5): ev4-1, 4, 6, 7; ev5-4, 5, 6, 7, 8, 9. Route ladder D1 -> D2 of sec. 2, with
the pass criterion of sec. 2 (a ballistic shoot closure at lambda 1 that passes the gate and re-flies
< 1 km).

Command (from the repo root; 2 workers; resumable):

```sh
for m in ev4-1 ev4-4 ev4-6 ev4-7 ev5-4 ev5-5 ev5-6 ev5-7 ev5-8 ev5-9; do echo $m; done | xargs -n1 -P2 sh -c 'uv run python scripts/ev_realeph_1025.py run $0 de --root de440_shoot --timeout 1800'
uv run python scripts/ev_realeph_1025.py shootcheck
```

(Rewritten on 2026-10-11 after the run as an equivalent `for` loop: the original `printf` line had been
written with literal line breaks inside its quotes. The lead's launch used the same member list,
flags and 2 workers.)

- Output: `data/1025_ev_realeph/de440_shoot/<member>/{D1,D2}/e<i>/`.
- Live logs: `tool.live.log`, gitignored by `data/1025_ev_realeph/de440_shoot/.gitignore`. The
  static `tool.log`, `realeph_chain.json` and `check.log` are committed.
- Progress lines go to `data/1025_ev_realeph/runlog.txt`.
- Runtime: about 35 s per shoot restart and 20 restarts, so about 12 min per epoch and route (ev4-1
  e0, measured). 10 members x 5 epochs is about 10 h serial and 5 h on 2 workers.
- `--timeout 1800` is per epoch and route. A timeout writes no result file and is reported "NOT RUN",
  never a fail.
- Resume: re-run the same command. Epochs with a `realeph_chain.json` are skipped; timed-out epochs
  are tried again.
- `shootcheck` prints PASS / FAIL / NOT REACHED / NOT RUN per member and epoch and writes
  `de440_shoot/summary.json`. A fresh agent then reads that file. The sec. 5 expectation stands:
  no gate-passing ballistic closure.

## 8. Results of the DE440 shoot (lead launch 2026-10-11 07:27-09:14 AEDT; commit 2858c349)

Data: `data/1025_ev_realeph/de440_shoot/`. Contents:
- `summary.json` from `shootcheck`;
- per member, route and epoch: the static `tool.log` and `realeph_chain.json`, including every
  restart's residual and the shoot start state.
The figures below were re-derived from those files for this section.

Run:
- The sec. 7 command, 2 workers, `--timeout 1800`. No epoch timed out, so no "NOT RUN".
- 10 members x 5 epochs.
- D1's blend reached lambda 1 in 45 member-epochs. The 5 that stalled went to D2 as the ladder says:
  ev5-4 e1 and ev5-5 e1 (at lambda 0.709), ev5-6 e0 and ev5-7 e0 (0.464), ev5-9 e1 (0.473).
- D2 reached lambda 1 only for ev5-9 e1. The other four direct solves did not converge.
- So 46 full-length shoots of 20 restarts each (920 restarts).
- About half the restarts were not real attempts. 481 of the 920 ended at 3 or fewer function
  evaluations: the 6.34 defect, where the first LM step leaves the domain even at the 0.03/0.05-rad
  perturbations. 439 ran their budget.

**Outcome: 0 PASS, 1 FAIL, 49 NOT REACHED.**
- NOT REACHED here means no shoot closure (max residual >= 1e-6), or no lambda 1 on either route.
- The single FAIL is ev4-4 e0. Restart 0, the unperturbed seed, closed ballistically, but the chain
  passes BELOW the Earth's surface: the Earth turn ratio is 1.763, with a required periapsis altitude
  of -4,953 km. Venus is 0.961.
- The sec. 5 expectation, "no gate-passing ballistic closure", is met.

Best shoot residuals per member: the smallest max-over-junctions residual over all restarts and
epochs; in brackets, the range of the per-epoch best.

| Member | Best residual (epoch, route) | Per-epoch best, range |
|---|---|---|
| ev4-7 | 0.0020 (e2, D1) | 0.0020-0.0065 at all five epochs |
| ev5-9 | 0.0008 (e1, D2); D1 best 0.0032 (e3) | 0.0008-0.030 |
| ev4-6 | 0.0026 (e4) | 0.0026-0.50 (e3 0.50) |
| ev4-4 | 0.0039 (e1) | e1-e4 0.0039-0.012; e0 is the FAIL closure (best non-closing restart 0.031) |
| ev4-1 | 0.073 (e2) | 0.073-0.46 |
| ev5-4 | 0.21 (e2) | 0.21-2.0 |
| ev5-7 | 0.32 (e1) | 0.32-3.7 |
| ev5-8 | 0.41 (e0) | 0.41-15.1 |
| ev5-5 | 0.66 (e2) | 0.66-3.5 |
| ev5-6 | 0.99 (e1) | 0.99-1.2 |

What a near-closure means physically:
- The shoot residual is in km/s. Each junction contributes its V_inf-magnitude mismatch, and each
  fixed-leg arrival contributes its position miss scaled as 1 km = 1e-3 km/s
  (`run_942_realeph_chain.py`, `POS_SCALE_KM`). The reported figure is the WORST single component.
- A residual of 0.002 therefore means: at the worst junction, the V_inf magnitudes differ by about
  2 m/s, or a full-rev leg misses its planet by about 2 km. Every other component is smaller.
- The best members stop with a worst junction defect of 0.8-7 m/s (or a 0.8-7 km miss) at: all five
  epochs of ev4-7, four of ev5-9, two of ev4-6 and one of ev4-4. ev4-4's other three non-FAIL
  epochs stop at 9.5-12 m/s.
- That is the same order as the sec. 6 DE440 closure-dV figures for the Venus legs at their best
  epochs: essentially ballistic or low-maintenance, i.e. a few to tens of m/s per 7 cycles. It does
  not contradict them.
- It is NOT a pass. The pre-registered criterion is a closure (< 1e-6), and none was met. A
  worst-component residual is not a dV, so it cannot be put into a band. And the minimiser stopped
  at a local residual floor, not at a closure.

Siblings share their Venus part:
- ev4-4, ev4-6 and ev4-7 have identical x0_days, identical final_y_dates and identical blend steps at
  every epoch. Their lambda 1 minimax worst ratios are identical: 0.920, 0.946, 0.979, 0.831 and
  0.772 at e0-e4, all Venus-bound.
  - They differ only in the order of their Earth full-rev returns in the shoot's fixed-leg seed:
    2 + 2 yr (ev4-4, RE/2:3 twice), 1 + 3 yr (ev4-6), 3 + 1 yr (ev4-7).
- ev5-4/ev5-5 and ev5-6/ev5-7 also have identical dates at every epoch. Their Earth blocks differ
  (ev5-4 one 3:2 return against ev5-5's 1:1 + 2:1; ev5-6 2:1 + 2:1 against ev5-7's 1:1 + 3:2), so
  their Earth turns and worst ratios differ at some epochs.
- So the 10 shoot members represent **6 distinct real-ephemeris transfer-and-Venus chains**:
  - ev4-1;
  - ev4-4/6/7;
  - ev5-4/5;
  - ev5-6/7;
  - ev5-8;
  - ev5-9.
- Grouped the same way (same transfer V_inf, same Venus block, same Earth-block length), the 22 members
  are **14 distinct chains**:
  - ev4-0;
  - ev4-1;
  - ev4-3;
  - ev4-4/6/7;
  - ev4-5 (its Venus block is in the other order);
  - ev4-8/9/10/11;
  - ev5-0/1;
  - ev5-2;
  - ev5-3;
  - ev5-4/5;
  - ev5-6/7;
  - ev5-8;
  - ev5-9;
  - ev5-10.

Positive control for the DE440 shoot (sec. 3 stop rule): **none exists at these settings.**
- The only earlier DE440 shoots of full-rev cyclers are ev-A and ev-B in #942 6.25
  (`docs/notes/2026-10-05-942-943-two-working-body-generator.md` sec. 6.25):
  - the shoot found ballistic closures at every reached epoch, but all failed the gate (ev-A
    1.12-2.99, ev-B 3.7-6.9);
  - they used the old 0.3/0.6-rad restarts (marked WEAK in 6.34-6.35).
- No DE440 shoot has ever closed a known cycler WITH a gate pass:
  - VenMar#45 (6.16, C1) and ev-C (6.23) pass on DE440, but they have no fixed legs, so no shoot
    runs.
  - The H&M endpoint control of 6.18 is on Standish, not DE440.
  - The Jovian gate-passing shoot closures (GanEur#316, 6.35/6.38) used jup365 and other flags.
- So "no closure" here is CONDITIONAL on the formulation: DE440's Earth with its lunar reflex, a 1:1
  seed timed on that Earth, 20 restarts of which about half die at the first step, and no
  chain-length continuation. It is not evidence that no ballistic DE440 orbit exists near these
  members.
- Registered as `#1057` (sec. 9). Cheapest control first:
  - (1) ev-A on DE440 with today's D1 settings at its 5 epochs: does today's shoot still find the
    closures 6.25 found? About 1 h on 2 workers.
  - (2) Then the same shoot for ev4-7, the best near-closure, with the Earth replaced by the
    Earth-Moon barycentre (a body DE440 carries): does the lunar reflex block the closure? That is a small
    chain-tool option, a code change.
  - Only a gate-passing closure in (2) would give the DE440 shoot a positive control.

No row is written; nothing is called novel.

## 9. PRE-REGISTRATION (2026-10-11, before the run): `#1057` step 1, the DE440 shoot control on ev-A; and a proposed step 0 (restart fix)

### 9.1 Step 1: ev-A on DE440 with today's D1 settings

Object:
- ev-A: k2|LE>V/0s|RV/1:1|LV>V/1h|LV>E/0s, ideal dates 583.944386694386, 984.2215547169349 and
  1576.2577920606109 d (the C2 control's key and dates);
- 5 cycles, 16 yr, as in #942 6.23/6.25;
- the tool's 5 standard epochs.

Settings and ladder: identical to the members' (sec. 2, sec. 7):
- D1 `--real de440 --shoot-jac sparse`, 20 restarts at 0.03/0.05 rad;
- D2 (`--direct`) only where D1's blend does not reach lambda 1;
- the same `classify` and the same checker.

Tool change for this, the only one: `scripts/ev_realeph_1025.py` can now address ev-A, as member
"ev-A" (the C2 key and dates, 5 cycles). Also, `shootcheck` takes `--root` and `--members`.

Command (from the repo root; 2 workers; resumable):

```sh
mkdir -p data/1025_ev_realeph/control_1057
printf '%s\n' 0,1,2 3,4 | xargs -n1 -P2 sh -c 'uv run python scripts/ev_realeph_1025.py run ev-A de --root control_1057 --timeout 1800 --epochs $0'
uv run python scripts/ev_realeph_1025.py shootcheck --root control_1057 --members ev-A
```

- Output: `data/1025_ev_realeph/control_1057/ev-A/{D1,D2}/e<i>/` (static `tool.log`,
  `realeph_chain.json`, `check.log`). Live logs `tool.live.log` are gitignored there.
- `control_1057/summary.json` is written by `shootcheck`. Progress goes to
  `data/1025_ev_realeph/runlog.txt`.
- Expected runtime: the members' shoots averaged 4.7 min each (46 in 107 min on 2 workers). ev-A's
  shoot is about the same size: 5 cycles x (3 Lambert + 1 full-rev leg) gives 30 unknowns, against
  36 for a k = 4 member. So about 15-30 min on 2 workers. The worst case is bounded by the timeout:
  3 epochs x 1800 s = 1.5 h.
- A timeout is "NOT RUN", never a fail; resume by re-running the command.

PASS criterion: identical to the members' (sec. 2):
- a ballistic shoot closure at lambda 1 (max residual < 1e-6);
- the #888 gate passes at the registry floors on every interior flyby ("indeterminate" is not a pass);
- the checker's DOP853 re-fly misses by < 1 km.

What each outcome means for the 22:
- **PASS at >= 1 epoch:** today's DE440 shoot can close a known full-rev Earth-Venus cycler with a gate
  pass. The shoot then becomes a usable judge of the members. Their sec. 8 "no closure" results stand
  as recorded, but are now a formulation-checked negative. They still carry the restart-defect caveat
  of 9.2 until step 0 is run.
- **FAIL (closures that fail the gate, as all the #942 6.25 ev-A closures did) or NOT REACHED (no
  closure) at every epoch:** the route cannot judge a full-rev Earth-Venus chain on DE440.
  - The 22 are then decided on the Standish and DE440-blend evidence alone. The owner summary's
    options (sec. 6 there) stay as written.
  - Step 2 of `#1057` (the Earth-Moon barycentre for the Earth) is the next control.
- **Expectation, before the run:** FAIL or NOT REACHED at all five epochs, by the 6.25 precedent (blend
  to lambda 1 at 5/5; closures found but all gate failures, 1.12-2.99). A PASS would contradict that
  precedent and be checked first, in this order: the checker's re-fly, the gate floors, then the
  epoch.

### 9.2 Proposed step 0 (NOT implemented): the restart defect

The defect:
- 481 of 920 restarts (52 %) ended at 3 or fewer function evaluations with an infinite recorded
  residual.
- Reading of the code (inferred, not instrumented): `_res` returns a constant 1e3 penalty vector where
  `chain_eval` fails (no Lambert solution, or a flight-time failure). The residual recorded at the end
  is from `_safe_eval`, which returns None there, hence "inf".
- A perturbed seed that is itself outside the domain therefore starts LM on a constant function. Its
  Jacobian is zero, so LM stops at once (nfev about 2). The restart was never a real attempt.
- (#942 6.34 described it as LM's first step leaving the domain. The nfev of 2 is more consistent with
  the seed itself being out of the domain. Step 0's own check below tells the two apart.)

Proposed fix (chain tool, `run_restarts`), with the same budget of 20 LM runs per shoot:
- Before restart r > 0 is run, evaluate `chain_eval` at the perturbed seed.
- If it fails, re-draw the perturbation (same distribution, same rng stream), up to 200 draws.
- Only an in-domain seed starts LM. Count the draws and record them.
- If 200 draws all fail, record the restart as "no in-domain seed"; it is not a closure attempt.
- Restart 0, the unperturbed seed, is unchanged.

The fix does NOT change:
- the closure threshold;
- the perturbation scale;
- the number of LM runs;
- the gate or the re-fly.

Pre-registered check of the fix (before any re-run is read):
- On one shoot from sec. 8 (ev4-7 e2), the number of restarts ending at nfev <= 3 must fall to about 0.
- If it does not, the 6.34 reading (the first STEP leaves the domain) is the right one. The cheapest
  fix is then a finite domain penalty that grows with the distance to the domain, so LM backtracks.
  That is a separate amendment.

What it would change:
- About twice as many real closure attempts per shoot. The count of attempts changes for certain.
- Whether the result changes is not known in advance:
  - for ev4-7 and ev5-9 the real attempts that ran already reached residuals of 0.0008-0.002 and
    stopped at a floor;
  - more real attempts raise the chance of a closure but do not guarantee one.
- By the rule that a correctness fix invalidates prior negatives, the 49 NOT REACHED of sec. 8 would
  have to be re-run under the fix before they count as negatives. Cost: about 2x today's shoot time,
  roughly 4-5 h on 2 workers for the ten members.
- Order: step 1 (9.1) first. If step 1 or step 2 shows the shoot can pass a known cycler, step 0 then
  comes before any re-run of the members.

Registered as `#1058` (a chain-tool change in a shared tool, with its own check), from `#1057`.

### 9.3 Result of step 1 (lead launch 2026-10-11 10:22-10:33 AEDT; commit 3ae396f7)

Data: `data/1025_ev_realeph/control_1057/ev-A/D1/e<i>/` (static `tool.log`, `realeph_chain.json`) and
`control_1057/summary.json`. The figures below were re-read from those files.

**ev-A FAILS at all 5 epochs.**
- The blend reached lambda 1 at every epoch, so D2 never ran. The minimax gate there passes at every
  epoch (worst 0.58).
- The shoot then CLOSES ballistically at every epoch: 17, 19, 18, 18 and 20 of 20 restarts converge.
  Most land on one or two distinct closures; the distinct worst ratios per epoch are 4, 3, 2, 2 and 4.
- None of the 20 restarts at any epoch ended at nfev <= 3.
- Every closure fails the #888 gate at VENUS. The best closure per epoch:

| Epoch | JD | Closing restarts | Best closure worst ratio (Venus) | Earth ratio | Required Venus periapsis altitude |
|---|---|---|---|---|---|
| e0 | 2462508.0 | 17 / 20 | 1.148 | 0.355 | -848 km |
| e1 | 2464843.6 | 19 / 20 | 1.571 | 0.487 | -2,961 km |
| e2 | 2467763.3 | 18 / 20 | 1.117 | 0.492 | -627 km |
| e3 | 2469516.4 | 18 / 20 | 2.314 | 0.493 | -4,614 km |
| e4 | 2472436.0 | 20 / 20 | 1.151 | 0.494 | -862 km |

- A negative required altitude is below the surface. At the best epoch (e2) the required periapsis is
  627 km below Venus's surface, i.e. 927 km below the 300-km registry floor.

Against the sec. 9.1 pre-registered meanings:
- The outcome is "FAIL at every epoch": the D1 shoot route CANNOT judge a full-rev Earth-Venus chain on
  DE440.
- The 22 are decided on the Standish and DE440-blend evidence alone. The owner options (a)-(d) stand as
  written.
- The members' sec. 8 NOT REACHED results stay CONDITIONAL. They are not a judged negative.
- This matches the pre-run expectation and the #942 6.25 precedent (closures found, all gate failures,
  1.12-2.99, then with the old restarts). It reproduces that precedent at today's restart settings.

What it says about ev-A itself:
- Its Standish standing is unchanged: 5/5, gate-passing, re-flown (C2, sec. 3).
- Its catalogue row (`earth-venus-two-working-body-cycler-eva-2026`) already says "NEAR-BALLISTIC on
  DE440, no gate-passing ballistic member found". Its `data_gaps` entry still calls the DE440 shoot
  restarts "old 0.3/0.6-rad (WEAK)", and that is now superseded.
- Proposed `data_gaps` text, for a LATER catalogue window, not applied now:
  > "On DE440 (astropy Earth, not the Earth-Moon barycentre) the chain-tool shoot at today's settings
  > (restarts 0.03/0.05 rad; #1057 step 1, 2026-10-11) closes ballistic 5-cycle chains at all 5
  > epochs, but every closure fails the #888 gate at Venus: best epoch JD 2467763.3, worst ratio 1.117,
  > required periapsis 627 km below the surface (927 km below the 300-km floor). The blend's minimax
  > solution at lambda 1 passes the gate (0.58) but does not close; its closure dV is 1.56-2.53 m/s per
  > 16-yr chain (6.25)."
- The near-ballistic statement (a few m/s from the gate-passing minimax solution) still holds. What
  changed is that the restart caveat on the shoot is gone.

Why ev-A had no dead restarts, against 481 of 920 for the members. This is an inference from the
structure; it was not instrumented:
- ev-A has 1 fixed leg per cycle (the Venus 1:1), so 5 fixed legs in its 5-cycle chain.
- The ten members have 2 to 4 fixed legs per cycle: ev5-9 two Earth full-revs; ev4-4 two Earth and
  one Venus; ev4-1 two Earth and two Venus. So up to 12 perturbed fixed-leg directions per chain
  (ev4-1, 3 cycles).
- Several members have two Earth full-revs back to back. Under #1058's reading (a seed outside the
  domain), each extra perturbed leg is another chance for a downstream junction to fail.
- The DE440 Earth (lunar reflex) moves the Earth full-rev arrivals. ev-A has no Earth full-rev at all;
  Earth is only a Lambert junction there.
- #1058's check (it records the failure reason per draw) would test this.

Is step 2 (the Earth-Moon barycentre) still worth running?
- **As a positive control for the shoot: not as specified.**
  - The control fails at VENUS, and ev-A has no Earth full-rev leg, so replacing the Earth by the
    barycentre cannot change ev-A's failure.
  - For the members, the barycentre might let more chains close, but the control shows that a closure
    is not the problem. The shoot finds closures readily, and they fail the gate.
- The failure mode is the one #942 6.48 described on jup365. The blend's minimax solution passes the
  gate, but the free-direction shoot lands on other, gate-failing closures. The fix there was the
  anchored shoot (`--shoot-homotopy 10`, with `--grow-chain --grow-from-one --shoot-rel-time`). It
  continues from the minimax seed to the closure connected to it. It passed both of its jup365 controls
  (6.50: C4@2013 indeterminate 0.9917; GanEur#316 pass 0.741).
- **Proposed instead, as `#1057` step 1b** (pre-registration to be written before any run; lead launch):
  ev-A on DE440 with the 6.48 anchored-shoot flags.
  - Pass criterion identical to the members' (closure, gate pass at the registry floors, re-fly
    < 1 km).
  - A PASS would give the DE440 shoot its first positive control, and the members would then need a
    re-run with the same flags.
  - A FAIL would show that ev-A's gate-passing DE440 neighbourhood has no ballistic closure at this
    chain length.
  - Expected cost: of the order of step 1 (about 10-30 min on 2 workers).
- Step 2 (the barycentre) is kept in `#1057` as a later diagnostic for the members' closure EXISTENCE,
  not as the shoot's control.

### 9.4 PRE-REGISTRATION (2026-10-11, before the run): `#1057` step 1b, ev-A with the anchored shoot

Why: step 1 (9.3) showed the #942 6.48 failure mode on DE440. The blend's minimax solution passes the
gate, but the free-direction shoot lands on other closures, all gate-failing. On jup365 the anchored
shoot fixed this and passed both its controls (6.48 method; 6.50 results: C4@2013 indeterminate 0.9917,
GanEur#316 pass 0.741). This step tests it on DE440 for the first time.

Object: ev-A, as in 9.1 (the C2 key and dates, 5 cycles, the tool's 5 standard epochs).

Route ladder:
- A1 = D1 + the 6.48 flags: `--real de440 --shoot-jac sparse --shoot-homotopy 10 --grow-chain
  --grow-from-one --shoot-rel-time`.
- A2, only where A1's phase-1 blend does not reach lambda 1: A1 + `--direct`.
- Everything else as in the members' runs: 20 restarts at 0.03/0.05 rad, the closure threshold
  1e-6, `classify`, the checker.
- The anchored method, from 6.48:
  - the phase-1 date solve on 1 cycle;
  - at 1 cycle each restart follows a 10-step Newton homotopy from the minimax seed to the closure
    connected to it;
  - the chain then grows 1 -> 5 cycles, each length seeded by the previous gate-best closure;
  - a length without a closure ends the run.

Tool change (the smallest): `scripts/ev_realeph_1025.py` gains the model `dea` (routes A1/A2, i.e.
D1/D2 plus the 6.48 flags), and `shootcheck` also reads A1/A2. No test in `tests/scripts` covers this
script. Regression check: `shootcheck --root control_1057 --members ev-A` still reads FAIL at all 5
epochs, and its `summary.json` is byte-unchanged.

Command (from the repo root; 2 workers; resumable):

```sh
mkdir -p data/1025_ev_realeph/control_1057b
printf '%s\n' 0,1,2 3,4 | xargs -n1 -P2 sh -c 'uv run python scripts/ev_realeph_1025.py run ev-A dea --root control_1057b --timeout 1800 --epochs $0'
uv run python scripts/ev_realeph_1025.py shootcheck --root control_1057b --members ev-A
```

- Output: `data/1025_ev_realeph/control_1057b/ev-A/{A1,A2}/e<i>/` (static `tool.log`,
  `realeph_chain.json`, `check.log`). The live logs are gitignored (the directory's `.gitignore`,
  added with this note).
- Expected runtime:
  - Step 1 took 101-215 s per epoch.
  - The anchored shoot runs its 20 restarts at each length 1, 2, 3, 4, 5, plus the 10-step
    homotopy at length 1, so roughly 3-6x step 1: about 5-20 min per epoch, and 15-60 min for 5
    epochs on 2 workers.
  - Bounded by `--timeout 1800`: at most 3 epochs x 30 min = 1.5 h. A timeout is "NOT RUN", never a
    fail.

PASS criterion: identical to the members' (sec. 2):
- a ballistic shoot closure at lambda 1, of the full length (5 cycles), max residual < 1e-6;
- the #888 gate passes at the registry floors on every interior flyby ("indeterminate" is not a pass);
- the checker's DOP853 re-fly misses by < 1 km.

Meanings, written before the run:
- **PASS at >= 1 epoch:** the anchored shoot is a usable judge of a full-rev Earth-Venus chain on DE440
  (its first DE440 positive control). The members would then need to be run on it too. Costs,
  extrapolated from the members' shoot (4.7 min per shoot average) times the 3-6x of the anchored
  method. Step 1b's measured time replaces these figures:
  - the ten Standish >= 4/5 members: about 4-8 h on 2 workers;
  - one member per distinct chain (the 6 among the ten): about 2.5-5 h;
  - all 22 members: about 9-17 h;
  - the 14 distinct chains of the 22: about 6-11 h.
  - Siblings share the transfer and Venus part but differ in their Earth full-rev legs, so a run per
    chain judges the shared part only. A member-level verdict needs that member run.
- **FAIL or NOT REACHED at every epoch:** the shoot formulation still cannot judge a full-rev
  Earth-Venus chain on DE440. The 22 are decided on the Standish and DE440-blend evidence. The owner
  options (sec. 6 of the owner summary) stand as written.
- **Expectation, before the run:** not stated as a prediction. The jup365 controls passed with this
  method, but the DE440 step 1 closures all failed the gate. Either outcome is plausible.

Ratchet rule:
- This is a NEW route. A pass on it does not change any member's recorded result (sec. 8: 0 PASS / 1
  FAIL / 49 NOT REACHED on route D1).
- Only a re-run of the members on route A1/A2, pre-registered on its own, can.
- Until then the members' DE440 shoot results stay as recorded and conditional.

### 9.5 Result of step 1b (lead launch 2026-10-11 10:38-10:50 AEDT; commit 2ebae89b)

Data: `data/1025_ev_realeph/control_1057b/ev-A/{A1,A2}/e<i>/` and `control_1057b/summary.json`. The
figures below were re-read from the chain files (the `grow` records per length).

**No PASS at any epoch.** e0 and e1 FAIL, e2 NOT RUN, e3 and e4 NOT REACHED.

| Epoch | JD | Phase-1 blend (1 cycle) | Grow-chain closures by length: closing restarts, best worst ratio (required Venus altitude) | Outcome |
|---|---|---|---|---|
| e0 | 2462508.0 | lambda 1, minimax gate pass 0.592 | k1 16/20, PASS 0.708 (3,835 km); k2 18/20, PASS 0.701 (4,074 km); k3 20/20, PASS 0.757 (3,109 km); k4 20/20, fail 1.513 (-2,669 km); k5 19/20, fail 1.148 (-848 km) | FAIL (k = 5 closure fails the gate) |
| e1 | 2464843.6 | lambda 1, pass 0.585 | k1 19/20, PASS 0.763 (2,995 km); k2 20/20, fail 1.303; k3 20/20, fail 1.197; k4 15/20, fail 1.327; k5 17/20, fail 1.571 (-2,961 km) | FAIL |
| e2 | - | - | the tool crashed: a solver step left DE440's date range (jplephem OutOfRangeError, "segment only covers dates 1549-12-31 through 2650-01-25"), exit 1, no result file | NOT RUN |
| e3 | 2469516.4 | stalls at lambda 0.794 (A2 direct: does not converge) | none | NOT REACHED |
| e4 | 2472436.0 | stalls at lambda 0.841 (A2 direct: does not converge) | none | NOT REACHED |

- The k = 5 closures at e0 and e1 are the same as step 1's best closures: worst 1.148 and 1.571, the
  same required altitudes -848 and -2,961 km.
- With `--grow-from-one`, phase 1 runs on 1 cycle. Its blend stalls at e3 and e4, where step 1's
  5-cycle blend reached lambda 1. That is a property of the 1-cycle wrap junction (as in #942 6.55a),
  not of the epoch.

Against the sec. 9.4 meanings:
- No PASS at any epoch, so the anchored shoot cannot judge a full-rev Earth-Venus chain on DE440
  either.
- The 22 are decided on the Standish and DE440-blend evidence. The owner options stand as written.
- The members' D1 results stay as recorded and conditional (the ratchet rule of 9.4).

What the k = 3 pass / k = 4 fail pattern says about ev-A as a DE440 ballistic object:
- At e0 the grow chain finds gate-passing ballistic closures at 1, 2 and 3 cycles: 3 cycles = 6.4 yr,
  worst ratio 0.757. The best closure's lowest required periapsis altitude is 3,109 km at
  k = 3 (3,835 and 4,074 km at k = 1 and 2). The anchoring follows them
  continuously from the gate-passing minimax seed.
- From 4 cycles on, every closure found violates the Venus floor: k = 4 needs a periapsis 2,669 km
  below the surface, and k = 5 848 km below.
- So on DE440 at e0, ev-A is a gate-passing ballistic chain for about 6.4 yr. Over 12.8-16 yr the
  ballistic continuation needs Venus flybys that are physically impossible.
- At e1 the gate-passing ballistic chain lasts only 1 cycle (3.2 yr).
- This is consistent with the "near-ballistic" reading (a few m/s of correction per 16-yr chain from
  the minimax solution, 6.25). It is a reading, not a computed link.
- Limits:
  - the k = 3 closures were not re-flown independently (the checker runs only on the final k = 5
    result, and the per-length states are not stored);
  - one epoch only;
  - 3 cycles is below the 5-cycle (16-yr) length of the pre-registered criterion, so it is not a
    pass.

Does it belong in ev-A's proposed `data_gaps` sentence (9.3)? Yes, as a descriptive clause with its
limits stated. Proposed addition, for a later window, not applied now:
> "With the anchored shoot (#942 6.48 method; #1057 step 1b) the ballistic DE440 chain at JD 2462508.0
> passes the gate for 1-3 cycles (6.4 yr; worst 0.757, lowest required periapsis altitude
> 3,109 km) and violates
> the Venus floor from 4 cycles on (the 3-cycle closures not independently re-flown); at JD
> 2464843.6 it passes for 1 cycle only."

The e2 crash is registered as `#1059`: a date guard in the chain tool, so an out-of-range trial step is
a rejected step, not a crash. A correctness fix: e2 of step 1b is NOT RUN, not a negative, until it is
re-run under the guard.

## 10. Catalogue window patch for the 10 ev rows: PRE-REGISTERED expected ratchet deltas (2026-10-11, before the dry run)

Owner decision 2026-10-11: "(a) 10 rows, Standish >= 4/5" (owner summary sec. 7). The rows are ev4-1,
ev4-4, ev4-6, ev4-7, ev5-4, ev5-5, ev5-6, ev5-7, ev5-8, ev5-9.
- Prerequisite: V1 evidence `54b2401f` (`scripts/writeback_1025_ev.py`, `data/1025_ev_writeback_v1.json`,
  writeback test 31 rows).
- Apply script: in the scratchpad (`window/apply_1025_ev_window.py`). Every edit asserts exactly one
  match.
- The same window applies the ev-A `data_gaps` clause of sec. 9.5. twobody-gen2's gc-1/gc-2 V2 notes
  patch is separate and touches other rows.

Expected deltas:
- catalogue rows: 413 -> 423;
- `our_status: candidate-novel`: 19 -> 29;
- multi-arc (allowlist and census): 311 -> 321; single-ellipse 46 and non-keplerian 56 unchanged;
- V1 map (`test_schema_v45_fields`): +10 entries, all V1. The V2 count is unchanged (10), as are V3;
- validation tier census, "unvalidated": 134 -> 144 (this-project rows with no provenance-registry
  tags, as for ev-B);
- rediscovery, MULTI_ENCOUNTER_SEQUENCE: 225 -> 235 (all 10 are heliocentric with multi-encounter
  sequences, like ev-B); NON_HELIOCENTRIC unchanged (111); other exclusion reasons unchanged.

Any other ratchet that moves in the dry run is a surprise. It is reported, not absorbed.

Addendum (2026-10-11, before the window; lead request): commit 38bee4b5 adds three cross-references to
the row text in `scripts/writeback_1025_ev.py`:
- the #973 note secs. 13-14, as the enumeration, gauntlet and screen source;
- `data/973_ev/hm1970_subperiod_check.json` (#973 secs. 2.2, 6.5), as the structural evidence behind
  "outside those itineraries";
- `docs/notes/2026-06-22-dv-band-definitions.md` and `dv_band_acceptance`, as the threshold source.

What was checked:
- `data/1025_ev_writeback_v1.json` is byte-identical (md5 bce09fd8...).
- The apply script is unchanged (md5 3d89aa71...); it regenerates the row text from the writeback
  script.
- In a scratch worktree, the patched catalogue loads with 423 rows and no duplicate ids.
- ruff is clean on the touched files.
- The writeback test passes: 30 rows, 20 earlier plus the 10. The "31" in 54b2401f's message was a
  miscount.

The pre-registered deltas above are unchanged.
