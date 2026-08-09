# Task #808: real-binary grid-path wrong-topology gate — decision and fix

**Date:** 2026-08-09
**Task:** #808 (registered during #807; see
`docs/notes/2026-08-09-807-pc-33-branch-loss-topology-gate.md` for the
trigger and the pluto_charon fix this extends)
**Decision:** **Outcome 1 — extend #807's clean-negative gate** to
`real_binary_kk_sweep`'s grid paths (`_finalize_grid_candidate` AND the
same-gap SRP grid path `sweep_family_grid_srp`). No principled exemption
exists.

---

## Why extend (the evidence, not the assumption)

The #660 scope decision left the grid path ungated on the theory that
`_grid_seed_search` "already filters the SEED for correct topology; whether
that survives the subsequent C-sweep is an existing, separate question".
Investigation found that separate question was already answered — twice,
empirically — before this task:

1. **#656's (5,1) PC grid sweep is a measured in-module instance of the
   exact failure mode** (data/OUTSTANDING.md, #656 bullet, 2026-07-19):
   `_grid_seed_search` found a genuine prograde (5,1) seed
   (x0=-0.6685146994, C=3.05, T=28.1427 TU), then `c_sweep_find_nu_zero`
   (called with `hc=None`, auto-redetecting the crossing count each step)
   walked off that branch onto the unrelated retrograde (4,0) family and the
   run reported `stable_found=True, topology_ok=False` — hand-diagnosed at
   the time as a known fragility of the shared machinery (#627's own
   docstring: "an auto-redetected crossing index can snap onto a different,
   unrelated branch").
2. **#807 gated pluto_charon's own grid-seeded sweeps.** `sweep_21`/
   `sweep_22` use this very same `_grid_seed_search` +
   `c_sweep_find_nu_zero` machinery, and #807 applied
   `_topology_gated_result` to them. There is no architectural difference
   that would justify the sibling module's grid path behaving differently.
3. **No other safeguard exists.** `c_sweep_find_nu_zero` has zero topology
   awareness: failed/non-converged correction steps just `continue` with the
   stale (x0, T) seed, the C-walk is up to 20 + 60 + brentq re-corrections,
   and nothing bounds it to the seed's branch. The seed-topology filter is a
   filter on the STARTING point only.

Additionally found and closed in the same pass: `sweep_family_grid_srp` had
the identical gap, arguably worse — its gravity-only topology-filtered seed
then undergoes TWO un-re-checked continuations (beta-stepping AND a
C_srp-sweep) — while its anchor-path sibling `sweep_family_srp` has gated
inline since #665. New `_topology_gated_result_srp` mirrors #807's helper
for the SRP data shapes (records recovered topology + beta_nd/phi0).

## What changed (reporting only — per #807's precedent)

`src/cyclerfinder/search/real_binary_kk_sweep.py`:

- `_finalize_grid_candidate` now routes through
  `pluto_charon_kk_sweep._topology_gated_result` (reused, not re-implemented)
  BEFORE the #660 clearance gate: a stable orbit of a different winding
  topology is reported as a clean negative for the TARGET family with the
  recovered topology in `note`, instead of the confusing
  `stable_found=True, topology_ok=False` state. This also makes
  `_gate_clearance`'s own documented contract ("only called on a res with
  stable_found and topology_ok both True") actually hold on this path.
- `sweep_family_grid_srp` gets the equivalent gate via the new
  `_topology_gated_result_srp`.
- The topology CHECK itself (integer winding classification) is untouched —
  no solver, classifier, or tolerance change anywhere. The clean negative
  remains method-conditional as all negatives here are; the `note` preserves
  what WAS recovered, so a #656-(5,1)-style "genuine seed lost by the
  C-sweep" case stays diagnosable from the result itself.

`tests/search/test_808_grid_topology_gate.py` (new, 3 tests, ~7 s): feeds
the finalizers a real reconverged orbit (#659's recorded Antiope (2,2) IC,
test_660's own cheap-reconvergence precedent) under a deliberately
MISLABELED (5,1) target — the gate must produce the clean negative with
"recovered (k1,k2)=(2,2)" in `note` and must NOT have run the clearance
evaluation; plus a correctly-labeled SRP pass-through check at beta_nd=0
(exact gravity-only reduction, proven in test_665). The gravity-only
pass-through case is already pinned by
`test_660_antiope_22_fails_clearance_gate_explicitly` (topology gate passes,
clearance gate still rejects with #659's exact 13.75/3.16 km figures).

## Bug-fix-invalidates-past-searches re-verification (actually run)

The gate changes only how `stable_found=True ∧ topology_ok=False` grid
results are REPORTED. Census of every historical grid-path result:

- **#549** (8 grid probes, (2,1)/(2,2) × 4 systems): all clean negatives
  (seed-not-found or no-stable-window) — no result in the affected class.
- **#657** (round 2, incl. Antiope): the only stable grid hit ever,
  Antiope (2,2), is `topology_ok=True` — gate passes it through to the #660
  clearance gate, which still rejects it with #659's exact figures.
- **#660**: its regressions call `_finalize_grid_candidate` directly —
  re-run, unchanged.
- **#665** (SRP): grid phase was 0/24 at the SEED stage (no gravity-only
  (2,1)/(2,2) seed found in any system) — the gated code is never reached;
  anchor phase used the already-gated `sweep_family_srp`. Verdict unchanged.
- **#656** (not on the required list, but the one real historical instance):
  replayed the recorded (5,1) branch-loss deterministically (reconverge the
  recorded genuine (5,1) seed, run the exact `sweep_family_grid` downstream
  wiring). Reproduced #656's wreckage bit-for-recorded-digit
  (C=3.2243893, x0=-0.587245699, T=18.69835 TU) and the gated finalizer now
  returns `stable_found=False` with note "recovered (k1,k2)=(4,0),
  w1=+4.000 ... clean negative" — exactly matching #656's own hand
  diagnosis. #656's stamped verdict is UNaffected: its (5,1) row was
  adjudicated "UNSETTLED (known search-method gap), not certified-empty",
  and that adjudication lives in the #656 bullet + the
  `pluto-charon-kk-45-cycler-sweep-2026-07-19` empty-regions stamp, which
  the gate does not touch. The full SIGALRM-bounded grid re-run was NOT
  used for this comparison because the self-hosted CI runner was loading
  the machine (load ~28) and per-call timeouts are contention-sensitive;
  the deterministic downstream replay avoids that confound entirely.

**Live runs:** all 8 test files importing `real_binary_kk_sweep`
(`test_549`, `test_657`, `test_660`, `test_665`, `test_656`, `test_627`,
`test_629`, `test_808`): 42/42 pass, exit 0 (~4 min, under the new `-n 6`
setting per #809). `uv run ruff check` + `ruff format --check` clean on the
changed files (repo-wide ruff has 3 pre-existing findings in #799's
concurrent, uncommitted `vaquero_em_cyclers.py`/`screen_799_*` files — not
touched, per concurrent-agent rules). Full `uv run mypy src tests`:
clean, 841 files.

**Same verdicts confirmed: no past negative flips, no past positive
(there were none on these paths) affected.**

## Follow-up registered

- **#810**: #656's own recommended (5,1)-specific re-attempt (hold the
  seed's own half-crossing count fixed through the C-sweep instead of
  `hc=None` auto-redetection, mirroring #504's "sweep upward, hc fixed"
  (3,2) convention) was recommended in #656's bullet but never registered
  as a task; the (5,1) topology remains the one k2<=k1<=5 PC case whose
  negative rests on a known search-method gap. Registered now.
