# `#861`: Oberon positive-control gate for the conjugate-apse + fold-turning seeding fix — GATE FAILS

**Task:** `#861`, registered 2026-08-21 (found during `#860`'s advice, dispatched same day). Build
the two cheap fixes `#860` (Fable) recommended for `#859`'s topology-misidentification failure —
a conjugate-apse seed and pseudo-arclength fold-turning continuation — and gate both against
Uranus-Oberon's own six published unstable resonant families (AAS 24-288) before committing any
compute to a novel system.

**Verdict: GATE FAILS, decisively.** 0 of 6 published families recovered with matching topology,
against a `>=4/6` pass criterion. Per `#860`/`#861`'s own decision rule, the Resonant Atlas
direction (`#789`/`#859`) is **SHELVED**; the tractable path forward runs through `#790`'s
corrector/alphabet build instead (`#858`'s own ranking already put `#790` ahead of `#789` on
expected value).

---

## What was built (both pieces work correctly as engineering, just don't solve the problem)

- **Conjugate-apse seed** (`jovian_resonant_families.py`): the excluded, "numerically hazardous"
  `x0_sign=+1` two-body ellipse's SAFE far-apse point, as an alternative starting guess to
  `two_body_resonant_seed`'s default opposition-phase seed. 10 new tests, all passing.
- **Fold-turning Stage A' worker** (`resonant_atlas_stage_a_prime.py`, `fold_turn_family`): wires
  `cr3bp_jacobi_arclength.py`'s (`#249`) pseudo-arclength continuation in place of
  `continue_family`'s natural-parameter walk, so a family that starts on the stable branch is
  walked through its own fold rather than stalling there. Confirmed working as designed — every
  gate cell terminates honestly at `step_underflow` (an arclength floor, the walk's own natural
  boundary), not a crash or silent truncation.

## The gate result

Driver: `scripts/run_861_oberon_gate.py`. Data: `data/found/861_resonant_seeding_oberon_gate/results.jsonl`
(12 records: 6 published ratios x 2 seed kinds, `--max-steps 100-250`, `--ds-max 0.06`).

| p:q | seed | n_unstable | in published C-range | topology matches | max &#124;λ&#124; |
|---|---|---:|---:|---:|---:|
| 3:4 | opposition | 0 | 0 | 0 | 1.00 |
| 3:4 | conjugate_apse | 0 | 0 | 0 | 1.00 |
| 4:5 | opposition | 26 | 2 | 0 | 1.15 |
| 4:5 | conjugate_apse | 12 | 12 | 0 | 4.56 |
| 5:6 | opposition | 0 | 0 | 0 | 1.00 |
| 5:6 | conjugate_apse | 1 | 1 | 0 | 1.16e4 |
| 4:3 | opposition | 0 | 0 | 0 | 1.00 |
| 4:3 | conjugate_apse | 42 | 42 | 0 | 1.28e5 |
| 5:4 | opposition | 0 | 0 | 0 | 1.00 |
| 5:4 | conjugate_apse | 19 | 19 | 0 | 2.23e4 |
| 6:5 | opposition | 40 | 0 | 0 | 1.11 |
| 6:5 | conjugate_apse | 0 | 0 | 0 | 1.00 |

**Topology matches: 0/12 cells, 0/6 ratios.** The `opposition` seed mostly reproduces `#859`'s
own original finding (flat at the stable branch, `|λ|≈1`, with two mild exceptions — 4:5 and 6:5
opposition both find *some* weakly unstable members, `|λ|` 1.11-1.15, but with zero or negligible
topology match either).

**The `conjugate_apse` seed genuinely fixes the "always lands on the stable branch" problem — and
trades it for a different one.** Four of six ratios (4:5, 5:6, 4:3, 5:4) produce STRONGLY unstable
members, sometimes overwhelmingly so (`|λ|` up to 128,060 at 4:3), and every single one of them
falls inside the target ratio's own published Jacobi-constant range. This is not a near-miss or a
borderline case: it looked, at first pass, like a clean win. It is not.

## Why: a real family-mixing artifact, not a topology-checker bug

Independently spot-checked the strongest case (4:3, conjugate_apse, `|λ|`=128,060) before trusting
the `topology_matches=0` verdict: `period_over_2pi = 32.33` against a target `q=3`, `winding_p_inertial
= 32.33` against a target `p=4` — off by roughly 10x, not an off-by-one or sign-convention slip in
the checker. `closest_secondary_approach_nondim = 0.0353` — far tighter than any of the flat-`|λ|=1`
stable branches (0.16-0.42) and tighter than the weakly-unstable opposition-seed members too.

This is a genuine, single, dominant, extremely unstable high-winding resonance structure that the
conjugate-apse continuation reliably walks into REGARDLESS of which target ratio seeded it — its
own Jacobi-constant range happens to overlap several of the published target bands, which is
exactly the "family mixing" `#776` already documented for Neptune-Triton (two topologically
distinct branches sharing a label at nearby C) — here reproduced as an ARTIFACT OF THE FIX ITSELF,
not of the raw two-body seed. Fixing the "converges to the boring stable branch" failure mode
surfaced a worse one: converges to a real but WRONG unstable branch, with high confidence (large
`|λ|`, tight approach) that could look like a genuine discovery if the topology check weren't run.

## Process note: a checkpoint-file race, caught and fixed

Mid-run, a genuinely orphaned background process (the dispatched agent's own self-backgrounding
attempt died, per this project's now well-documented pattern, but the actual compute survived) and
a second process the coordinating session launched to resume the same driver both started within
about a minute of each other, both read the same "already done" checkpoint state, and both began
appending to `results.jsonl` concurrently. Caught via `ps aux` before any interleaved-write
corruption could occur; the redundant process was killed each time (twice), and the resulting file
had 3 exact-duplicate records (`5:6 conjugate_apse`, `4:3 opposition`, `5:4 opposition`, each
appearing twice with byte-identical content — no corruption, just duplication) which were
deduplicated before this analysis. No cell was ever computed with a different result across the
duplicate pair, so the race did not affect the gate's substantive outcome, only the record count.

## Disposition, per `#860`/`#861`'s own decision rule

**Gate fails on 6/6 ratios (need `<=2` failures to pass at `>=4/6`).** Per the rule stated in both
this task's own registration and `#860`'s Sec. 5: **shelve the Resonant Atlas direction**
(`#789`/`#859`) — do not dispatch the real 3-novel-system Stage A run with this seeding, fixed or
unfixed. `#790` (itinerary enumeration, already ranked ahead of `#789` on expected value by `#858`'s
review, and not blocked by this specific seeding problem since it consumes already-published family
alphabets) is the better use of further effort in this campaign family, per `#858`'s own ranking.

This is not a wasted build: the conjugate-apse seed and the fold-turning Stage A' worker are both
real, tested, working additions to the codebase (10 + existing tests, all green) — they simply
answer "does seeding+continuation reach an unstable branch at all" (yes, now) rather than "does it
reach the RIGHT one" (no). A future attempt at this problem should start from that distinction, not
repeat this one — e.g. a topology-constrained seed search (fix `p:q` explicitly in the seed
construction, not just hope the continuation lands there) would be the natural next idea, but is
explicitly out of scope for this task and not attempted here.

No `data/catalogue.yaml` writeback — this task produces no discoveries, only a methodology
validation (which failed).

## Verification

`tests/search/test_resonant_atlas_stage_a_prime.py` (10 tests, ~11s) + the pre-existing
`tests/search/test_resonant_atlas_stage_a.py` (17 tests) + `tests/search/test_cr3bp_jacobi_arclength.py`
(the reused `#249` fold-turner's own suite) all green. `ruff check .` / `ruff format --check .` /
full `mypy src tests` clean. The gate run itself (12 cells, `data/found/861_resonant_seeding_oberon_gate/`)
is a data artifact, not a pytest-covered claim — its own numbers were independently spot-checked by
the coordinating session directly against the raw per-member records before this note was written.

## Addendum: `mu_continuation` cross-check + a concrete diagnosis of the family-mixing mechanism

A second, independent execution pass of this same task's own required scope (`#861`'s registration
explicitly asks for a `mu_continuation` cross-check against a Neptune-Triton table-verified saddle,
and a `deflated_newton` backstop, before declaring a ratio "genuinely unrecoverable" — neither is
covered by the gate/note sections above). Ran this after independently reproducing the exact same
12-cell gate result byte-for-byte (`scripts/run_861_oberon_gate.py`, identical output), so this
section is additive evidence, not a re-litigation of the FAILS verdict, which stands unchanged.

**mu-continuation cross-check #1 (4:5) — CONFIRMS the true family exists, right at the published
range's own edge.** `neptune_triton_resonant_families.ESM_GATE_ROWS["4:5-saddle"]` (table-verified,
`C=2.987089791658`, `half_crossings=3`, `|lambda|~=105` per `#771`'s own survey) continued in mu
(`mu_continuation.py`, `#249`'s sibling) from Neptune-Triton's `2.089503183689124e-04` down to
Oberon's own `3.54326e-5` (~5.9x descent, `#860` Sec. 4(d)'s own "modest, safe" assessment) —
`TARGET_REACHED` cleanly in 20 steps/8.7s, landing at `C=2.986914`, `|lambda|=22.95`, still
genuinely unstable throughout the descent (script: `scripts/run_861_mu_continuation_crosscheck.py`,
data: `data/found/861_resonant_seeding_oberon_gate/mu_continuation_crosscheck.json`). That landing
point's own period/winding (`period_over_2pi=4.966`, `winding_p_inertial=3.966`) match the 4:5
label to <1% — genuine topology confirmation, independent of both `#861` fixes. Its own C
(2.986914) sits ~0.0045 (0.15%) BELOW the paper's own printed lower bound (2.9914) — a small,
real, boundary-adjacent gap, not inside the strict printed range.

**Targeted follow-up (NOT part of either `#861` fix, diagnostic only):** re-corrected that
mu-continued member at `half_crossings=3` (the SAME crossing index the Neptune-Triton anchor uses)
up to `C=2.9914` (the published range's own lower edge) via a plain natural-parameter step, then
fold-turned it (`cr3bp_jacobi_arclength.continue_in_jacobi`, `half_crossings=3` throughout). Result:
a clean, topology-matching unstable segment, `C in [2.9914, 2.9926]` (9 members), `|lambda|` growing
70.8 -> 82 -> 68 (a genuine saddle magnitude, not a near-unit-circle value), `period_over_2pi` in
[4.894, 4.927] (1.1-2.1% off `q=5`) and `winding_p_inertial` in [3.894, 3.927] (1.4-2.7% off `p=4`)
— BOTH pass this task's own `PERIOD_REVIEWER_REL_TOL=0.03` for every one of these 9 members, with a
real close approach (`closest_secondary_approach_nondim=0.0123`, ~7,177 km). The walk then hits a
genuine topology jump at `C~2.9930` (period discontinuously jumps 4.89 -> 7.20) — the SAME
family-mixing artifact the gate note above describes, caught cleanly by `classify_member`'s own
period/winding check on the far side of the jump. Full data:
`data/found/861_resonant_seeding_oberon_gate/hc3_targeted_recovery_4_5.json`.

**This sharpens the gate note's own "family-mixing artifact" diagnosis into a specific, actionable
mechanism**: the genuine 4:5 saddle is real, sits right at the published range's own lower edge, and
IS reachable by fold-turning — but only from `half_crossings=3`, which neither seed phase's own
auto-detection (`_crossing_index_near_half_period`, run once per seed at its OWN natural C) ever
selects (opposition auto-detects `hc=4`; conjugate_apse auto-detects `hc=1` — see the gate table's
own `half_crossings` field in `results.jsonl`, neither is 3). `#860`'s two fixes (seed phase,
continuation method) are both necessary but not sufficient; the auto-detected CROSSING INDEX is a
third, unaddressed free parameter, and for 4:5 it is simply the wrong one before either fix's own
machinery gets a chance to matter. This does not change the gate's own literal verdict (the gate is
scored on `#861`'s two prescribed fixes exactly as built, and they still fail 0/6 as built) but it
answers the "is this a checker bug or a real dynamics fact" question the gate note itself raises,
and gives a concrete, sourced lead for anyone revisiting this direction later.

**mu-continuation cross-check #2 (4:3) — weaker, ambiguous.** Same procedure from
`ESM_GATE_ROWS["4:3-saddle"]` (`C=3.016635194282`, `half_crossings=2`) lands cleanly
(`TARGET_REACHED`, 30 steps/20.6s) at `C=3.016466` — comfortably INSIDE the published 4:3 range
`[2.9836,3.0279]` — staying unstable throughout (`|lambda|=3.95`, real but modest). Its own topology,
however, does NOT cleanly match either `(p=4,q=3)` or `(p=3,q=4)` under this task's own
period/winding convention, even at the Neptune-Triton SOURCE mu before any continuation
(`period_over_2pi=2.037`, ~32% off `q=3` and ~49% off `q=4`) — plausibly reflecting this exact ESM
row's own already-documented "family-mixing" complexity
(`neptune_triton_resonant_families.py`'s own module docstring: "all four printed `Res43+x+h` rows
are NOT samples of one single-crossing-index continuation branch"). Reported honestly as a weaker,
inconclusive cross-check, not a second clean confirmation. Data:
`data/found/861_resonant_seeding_oberon_gate/mu_continuation_crosscheck_4_3.json`.

**`deflated_newton` backstop — attempted, inconclusive, genuinely slow at these settings.** Ran
scalar deflated Newton (`deflated_newton.enumerate_roots`) on the 4:5 corrector residual at the
published range's own midpoint (`C=3.0035`), sweeping `x0 in [-1.9,1.7]` at `half_crossings in
{1,4}` (the two auto-detected indices). Found several distinct roots, including one genuine,
tightly-converged 4:5-topology match (`period_over_2pi=4.998`, `winding_p_inertial=3.998`,
<0.1% both) — but it is STABLE (`|lambda|=1`), not the paper's unstable saddle; no unstable
topology-matching root turned up in this bounded sweep. A wider hc/seed sweep (5 crossing indices x
20 seeds) hit a >250s wall-clock budget without finishing and was abandoned rather than let run
unbounded — a genuine time-budget limitation of this specific attempt, not a negative result on the
method itself (this is exactly the caveat the module's own docstring gives: narrow saddle basins
"still require the deflated Newton to pass nearby"; `hc=3`, the crossing index the mu-continuation
cross-check shows actually hosts the saddle, was not swept). No script/test artifact kept for this
attempt (exploratory only, not reproducible without re-running); numbers are recorded here for the
record.

**Disposition unchanged**: this addendum's own findings are point-for-point consistent with —
and independently corroborate — the gate note's own SHELVE verdict: the fold-turning/conjugate-apse
combination as built and gated genuinely fails 0/6, and the reason is now understood precisely
(crossing-index mismatch, not "no such family" or "checker bug"). Registered `#862` as a follow-up
lead (crossing-index-aware seeding) for whoever next revisits Resonant Atlas, per this project's own
"register everything" discipline — NOT a recommendation to un-shelve `#789`/`#859` now; that
decision still rests with `#790`'s own higher-ranked expected value per `#858`.
