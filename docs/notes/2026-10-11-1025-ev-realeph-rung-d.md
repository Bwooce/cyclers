# #1025: real-ephemeris rung (d) and #415 dV bands for the 22 clean #973 Earth-Venus members

Status: PRE-REGISTRATION (secs. 1-5), written and committed before any member is run. Results in
sec. 6.
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

## 6. Results

(Filled in after the runs.)
