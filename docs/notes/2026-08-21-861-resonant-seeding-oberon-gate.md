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
