# #826 — validation adjudication of #820's 8 truth-region closures: ZERO promotions, and #820's headline is overstated

**Date:** 2026-08-11
**Task:** `#826` (registered by `#820`, dispatched 2026-08-11) — decide, per row, what
validation-tier gate `probe_at_truth`'s evidence actually satisfies under the `#388`
promotion discipline, and perform any earned catalogue writeback.

## Headline

**No row is promoted. No `catalogue.yaml` edit was made.** Seven of the eight rows already
sit AT OR ABOVE the level `#820`'s registration tentatively proposed for them, so those
"tentatively V1 / V3" labels were never promotion proposals; the eighth
(`mcconaghy-2006-em-k2`) stays V0 on stronger grounds than the anchor caveat `#820`
recorded. Independently re-running `probe_at_truth` reproduced every one of `#820`'s
numbers — **and surfaced two things `#820` did not report**:

1. **Two of the eight "closures" are bend-INFEASIBLE** (`5.30ggF3`, `6.44Gg3`).
   `probe_at_truth` computes its verdict from `converged` alone while the result object
   carries `bend_feasible` / `constraints_satisfied` right beside it. Spec §14 **V0**
   requires `bend <= max` as a hard constraint, so a converged-but-infeasible closure is
   not even V0-admissible. #820's "8 of 12" is honestly **6 admissible + 2 not**.
2. **A new sourced cross-check — Russell's published turn ratio — validates 6 rows and
   falsifies 1.** The measured `min(max_bend / required_bend)` reproduces the row's own
   published `invariants.turn_ratio` to **0.001-0.024** on 7 of 8 rows. The one that fails
   is `5.30ggF3` — *the very row `#820` singled out as its strongest result*.

## 1. The tier discriminator (stated once, applied to all 8)

Read from `docs/spec.md` §14 and the `_LEVEL_EVIDENCE` over-claim registry in
`src/cyclerfinder/data/validate.py` (a row may declare above V0 only with a registry entry —
"a level is only as high as the recorded evidence mechanically supports; when in doubt, V0").
All 8 rows are `orbit_class: cycler` (so the strict-cycler gate applies) and
`cycler_class: multi-arc`.

| Tier | What it requires | Does `probe_at_truth` supply it? |
|---|---|---|
| **V0** | hard constraints (V∞ cap, r_p ≥ r_p_min, **bend ≤ max**), V∞ magnitude preserved per flyby, closure residual ≤ tol (idealized) | Partly — residual + V∞ continuity yes; **bend NOT checked by the probe** (§2) |
| **V1** | every leg re-solved with **lamberthub izzo + gooding**, agreement < 1e-3 m/s, AND full re-propagation with the **Kepler** propagator | **No.** The probe calls `ballistic_correct` only — neither half is run |
| **V2-ballistic** | **≥3 continuous laps**, bounded rotating-frame drift, in the row's defining model | **No.** The probe solves ONE period (§4) |
| **V3** | phase-matched real launch window + ephemeris horizon TCM | **No.** The probe ran `--model circular` |
| **V4 / V5** | independent codebase (GMAT) / human review | Not in play |

So the ceiling of this evidence class is V0, and that is what every row already has or
exceeds. Three further grounds, each independently sufficient to decline promotion, follow
`#388`'s own precedent on `russell-ocampo-4.3.1-5`
(`docs/notes/2026-06-23-388-russell-ocampo-4.3.1-5-anchor-recovery.md`, which declined
promotion for an anchor-recovering converged close):

- **Not the canonical determination path.** `#388` promotes/declines on
  `close_row_dsm(row, Ephemeris("astropy"))` — the single-shot path the gauntlet invokes.
  `probe_at_truth` is a bespoke probe on the **circular** ephemeris.
- **Seeded exactly at truth, not blind.** The corrector is handed the row's own sourced
  geometry. `#820`'s own blind 256-epoch grid found the anchor basin on only 2 of 12 rows.
  `#388`'s language for this is "epoch+seed-fragile".
- **Single-cycle.** Every tier above V1 needs multi-lap evidence.

## 2. The bend finding (independent re-run, `#826`)

Re-ran `probe_at_truth`'s exact steps and kept the `BallisticClosureResult` instead of the
dict the probe returns (which drops `bend_feasible` entirely). All of `#820`'s reported
numbers reproduced exactly — truth residual, solved residual 0.000, ToF drift, emerged V∞.
Additionally:

| row | converged | `bend_feasible` | `vinf_cap_ok` | `constraints_satisfied` |
|---|---|---|---|---|
| mcconaghy-2006-em-k2 | ✓ | ✓ | ✓ | ✓ |
| russell-ch4-4.991gG2 | ✓ | ✓ | ✓ | ✓ |
| russell-ch4-9.353Gg2 | ✓ | ✓ | ✓ | ✓ |
| russell-ch4-3.78Gg3 | ✓ | ✓ | ✓ | ✓ |
| russell-ch4-9.94Gg3 | ✓ | ✓ | ✓ | ✓ |
| **russell-ch4-5.30ggF3** | ✓ | **✗** | ✓ | **✗** |
| russell-ch4-5.75ggF3 | ✓ | ✓ | ✓ | ✓ |
| **russell-ch4-6.44Gg3** | ✓ | **✗** | ✓ | **✗** |

The code defect behind the blind spot: `probe_at_truth`'s `stayed` reads
`solved.converged`, and `run_row`'s grid filter is `[r for r in results if
r.get("converged")]`. Neither consults `constraints_satisfied`. The campaign's entire
verdict vocabulary (CLOSE-AND-MATCH / CLOSE-OFF-ANCHOR) is therefore bend-blind. Registered
as **`#829`**. Per `[[feedback_verify_automated_ghost_guard_booleans]]`, `6.44Gg3`'s
razor-thin **+3.302°** excess was exactly the kind of margin that should trigger inspection
of the gate's own code — which is how this was found.

## 3. Russell's published turn ratio as an independent sourced cross-check

Russell-Ocampo 2003 p.13 (digest `2026-06-17-digest-russell-ocampo-2003.md`):

> "Turn Ratio TR = max physically allowable turn angle / max required turn angle (δ_MAX)";
> TR > 1 ⇒ all flybys physically attainable; max allowable based on **200 km altitude Earth
> flyby**.

Our `_max_bend_deg` uses `PLANETS[body].safe_alt_km = 200.0` km for Earth and Mars — itself
sourced from Russell 2004 p.165 `r_p,min`. So `min` over flybys of
`max_bend / required_bend` is a **like-for-like reproduction of Russell's own TR**, and the
published per-row `invariants.turn_ratio` is a legitimate golden expected side. This
constrains **flyby turn angles**, which the V∞ anchors do not — a genuinely independent
axis of evidence.

| row | published TR | measured TR | Δ | per-node (max/req) |
|---|---|---|---|---|
| mcconaghy-2006-em-k2 | 2.65 | 2.658 | +0.008 | M1:∞ E2:2.658 |
| russell-ch4-4.991gG2 | 2.65 | 2.658 | +0.008 | M1:∞ E2:2.658 |
| russell-ch4-9.353Gg2 | 1.70 | 1.702 | +0.002 | M1:∞ E2:1.702 |
| russell-ch4-3.78Gg3 | 1.81 | 1.786 | −0.024 | M1:∞ E2:1.786 |
| russell-ch4-9.94Gg3 | 1.10 | 1.101 | +0.001 | M1:∞ E2:1.101 |
| **russell-ch4-5.30ggF3** | **1.27** | **0.613** | **−0.657** | M1:∞ **E2:0.613** E3:1.255 |
| russell-ch4-5.75ggF3 | 1.34 | 1.347 | +0.007 | M1:∞ E2:84.369 E3:1.347 |
| russell-ch4-6.44Gg3 | 0.95 | 0.956 | +0.006 | M1:∞ **E2:0.956** |

Identical at 64 and 256 phase epochs.

**This resolves the two bend flags into two completely different verdicts:**

- **`6.44Gg3` is NOT a defect.** Measured TR 0.956 vs published 0.95. This row is from
  Russell **Table 4.13, the NEAR-BALLISTIC table**, tagged "NEAR-BALLISTIC: TR = 0.95 < 1.0"
  in `docs/notes/multi-arc-classification.md`. Reproducing TR < 1 is *evidence the closure
  is on the sourced family* — Russell's own cycler requires a turn marginally beyond the
  200 km-flyby limit and needs a small powered nudge. The closure is faithful; it is simply
  not strictly ballistic, exactly as published.
- **`5.30ggF3` IS off-family.** Its third node reproduces the published 1.27 (measures
  1.255), but its **second Earth node demands 141.244° where only 86.582° is available**.
  Russell's cycler has no such turn. So `#820`'s lowest-truth-residual, sole
  STAYED-AT-TRUTH row is not the sourced geometry — a textbook
  `[[feedback_orbit_closure_discipline]]` "it closed!" catch, invisible to the closure
  residual and to the V∞ anchors, caught only by a third, independent, sourced constraint.
  Registered as **`#835`**.

### Caveat on the Mars node (do not over-read the Mars anchor)

`required_bend` is **0.000°** at the Mars node on all 8 rows. That is a genome artifact: the
designated arc is split at the Mars encounter into legs 0 and 1 of the *same* conic, so the
"flyby" imposes no turn. The Mars node therefore constrains only "the arc passes Mars's
position", and the emerged Mars V∞ is close to being determined by the conic seeded from
Russell's own printed ToF and transfer angle. The **Earth** nodes carry real bend and V∞
continuity constraints. So the Earth-side agreement is the load-bearing evidence; the
Mars-side agreement is nearer an internal consistency check on Russell's printed numbers.
(Consistently, Russell's published TR is set by the Earth node on every row, which our
measurement reproduces.)

## 4. V2-ballistic: REOPENED, not answered

`tests/search/test_free_return_v2_ballistic.py` already ran the V2-ballistic gate on four of
these rows and recorded a no-promotion finding whose stated ground was **structural**:

> "the single ellipse does NOT represent the Earth-to-Earth resonant phasing intervals …
> There is therefore no continuous ≥3-lap trajectory to propagate for these objects — the
> V2-ballistic gate is structurally inapplicable to a single-arc slice of a multi-arc cycler."

**#820 removes exactly that obstacle.** Its re-posed genome tiles the FULL cycler period,
loop legs included — a continuous multi-arc trajectory now exists to propagate. So the
recorded reason for inapplicability is obsolete for the 6 admissible rows, and V2-ballistic
is a live question for the first time.

It is **not answered here**, and the answer is not free: the period constraint uses the
rounded `period.years: 4.27` rather than the exact 2 × 2.1354 = 4.2708 yr, and `#820`'s own
note records that this injects ~1.3 d of slack-leg error. That error is **secular** — on the
order of 3×10⁶ km of rotating-frame drift per lap against a 50,000 km tolerance — so a V2
attempt must first re-pose on the exact synodic period. Registered as **`#830`**.

## 5. `mcconaghy-2006-em-k2` — V0, and the probe is not evidence for it at all

The registration framed this as an anchor-mismatch caveat. The row's own fields make the
verdict sharper. `mcconaghy-2006-em-k2` and `russell-ch4-4.991gG2` share, byte for byte:

- `orbit_source: russell-2004-t49_413`
- `free_return_arcs`: `g(1.4612,526.02,Ll)` + `G(2.8096,651.46,U)`
- `aphelion_au: 1.64`, `inclination_deg: 0.0`, `turn_ratio: 2.65`
- `period` k=2 / 4.27 yr, `sequence_canonical: "E-E-M-M"`

They differ only in `transit_times_days` (153 vs 150) and in the V∞ fields
(`vinf_source: mcconaghy-2006`, 4.7/5.0 vs `vinf_source: russell-2004-t49_413`, 4.99/5.10).
The row's own `orbit_elements.note` already records that Russell tags Table 4.9 row 1 as
parent cycler 4.991gG2 (#83) and states "Also known as the 'S1L1' cycler".

`build_genome` reads `free_return_arcs` — Russell's — and never touches the 4.7/5.0 fields.
So the probe on `mcconaghy-2006-em-k2` **is** the `4.991gG2` probe: both return emerged V∞
**5.008 / 5.107 identical to three decimals**, the same `E-M-E-E` sequence, the same
designated arc, and the same measured TR 2.658. It supplies **zero independent evidence for
this row**.

And the row's own cited anchor is not reproduced: `vinf_source: mcconaghy-2006` gives
4.7 km/s at Earth against an emerged 5.008. The 0.308 km/s gap sits inside the campaign's
`TOL_VINF = 0.5` km/s — **that threshold must not launder this**. It is a campaign screening
tolerance, not a tier gate, and a value published to two significant figures missed by ~6.5%
is not agreement. Validating a row means matching the anchor *that row cites*.

**Verdict: V0 stands**, on the ground that the probe is not independent evidence for this row
and its own anchor is unreproduced — not merely "a caveat noted". The probable-duplicate-rows
question (one physical object carried as two rows at V0 and V3, with census implications) is
registered as **`#831`**; rows were NOT merged here, as that ripples through every frozen
census ratchet.

## 6. Per-row verdicts

| row | current level | `#820` tentative | `#826` verdict | ground |
|---|---|---|---|---|
| `russell-ch4-4.991gG2` | **V3** | V3 | **no change** | already V3 (spec §14's own V3-ballistic type specimen); probe is circular-model single-cycle, strictly weaker |
| `russell-ch4-9.353Gg2` | **V1** | V1 | **no change** | already V1; probe supplies no izzo+gooding / Kepler reprop |
| `russell-ch4-3.78Gg3` | **V1** | V1 | **no change** | already V1; ditto |
| `russell-ch4-9.94Gg3` | **V1** | V1 | **no change** | already V1; ditto (binding TR 1.101, the tightest feasible margin) |
| `russell-ch4-5.75ggF3` | **V1** | V1 | **no change** | already V1; ditto |
| `russell-ch4-5.30ggF3` | **V1** | V1 | **no change; evidence RETRACTED** | closure is bend-infeasible AND off-family at node 2 (TR 0.613 vs published 1.27) — not V0-admissible, cannot support any promotion |
| `russell-ch4-6.44Gg3` | **V1** | V1 | **no change** | already V1; closure is bend-infeasible, but faithfully so (published TR 0.95 < 1, near-ballistic) — so it is a valid *near-ballistic* closure, not a *ballistic* one |
| `mcconaghy-2006-em-k2` | **V0** | V0 | **no change** | probe is the 4.991gG2 closure (identical genome); own 4.7/5.0 anchor unreproduced (§5) |

Existing levels were also checked for **downgrade** risk, since `#820` fixed a real
mis-posing and bug fixes cut both ways (`[[feedback_bugfix_invalidates_past_searches]]`).
The V1 evidence for these rows comes from the `#137` free-return single-ellipse genome and
the `#181` joint-(epoch, ToF) closer, not from `build_genome`.
`tests/search/{test_free_return_v1_mechanics,test_closer_sweep_v1,test_russell12_likeforlike_probe}.py`
were re-run: **46 passed**. No downgrade.

## 7. Note on `#388`'s energy-selectivity pattern

`#388` characterised the wall as energy-selective (low-energy near-Hohmann rows recover
anchor; high-V∞ rows collapse off-anchor). `#820`'s closers include high-V∞ rows
(9.35/10.52, 9.94/10.76) that DO recover anchor, which looks like a counterexample — but it
is not a like-for-like comparison: `#388`'s finding is about the **blind real-ephemeris
`close_row_dsm` lane**, whereas these are **circular-model probes seeded exactly at truth**.
The pattern is untouched. Whether it survives in the real-eph lane under the corrected
posing is a real question, registered as part of **`#830`**.

## 8. Deliverables

- **No `catalogue.yaml` edit.** Zero promotions ⇒ nothing earned a writeback. Deliberately
  no "ADDED EVIDENCE" annotations either: citing `#820`'s closure on rows whose closure is
  inadmissible (`5.30ggF3`) or non-independent (`mcconaghy`) would be misleading.
- **`tests/search/test_russell12_probe_bend_feasibility.py`** (new, 17 tests, not `slow`) —
  records the bend finding and the published-TR cross-check with teeth. Golden discipline:
  the expected side of every TR assertion is the row's own published
  `invariants.turn_ratio`, never a computed value.

## 9. Follow-ups registered

- **`#829`** — `campaign_russell12` is bend-blind: gate/report on `constraints_satisfied` in
  both `probe_at_truth` and `run_row`; re-classify all 12 rows.
- **`#830`** — V2-ballistic re-adjudication on the `#820`-reposed full-period genome
  (re-pose on the exact 2 × synodic period first); includes the `#388` real-eph re-check.
- **`#831`** — probable duplicate rows `mcconaghy-2006-em-k2` / `russell-ch4-4.991gG2`.
- **`#832`** — `turn_ratio`'s catalogue comment is INVERTED vs the sourced definition
  (comment says "required / max ballistic turn", source says "max allowable / max required");
  the measured values confirm the source. Comment-only fix across all rows carrying the field.
- **`#833`** — promote the measured-TR check to a reusable sourced validation instrument
  (it caught what the residual and the V∞ anchors both missed).
- **`#834`** — audit `russell-ch4-8.049gGf2`'s V3 `_LEVEL_EVIDENCE` entry against `#820`'s
  35.22 km/s truth residual under the corrected posing (likely a different real-eph lane, so
  probably not a contradiction — but unverified).
- **`#835`** — diagnose `5.30ggF3`'s spurious 141.244° node-2 turn under the reposed genome.

## 10. Verification

- `uv run ruff check .` / `ruff format --check .` — clean.
- `uv run mypy src tests` — clean (845 source files).
- `uv run pytest tests/data tests/search -q` — **exit 0, zero FAILED/ERROR** (full ratchet,
  not a subset).
- `uv run pytest tests/scripts -q` — exit 0.

A FIRST full-ratchet run, made while a self-hosted CI job was saturating the machine
(6 workers, load average 16-32), exited 1 on two tests —
`test_656_pc_higher_kk_sweep.py::test_656_grid_seed_search_recovers_admitted_pc_32_seed` and
`test_earth_moon_class1_resonant_connections.py::test_find_homoclinic_default_k_range_is_too_narrow_for_this_orbit`.
Both are the documented CPU-contention flake class, not a regression:

- `#809`'s own bullet names the second test explicitly as one of two "deliberately
  near-boundary/margin tests" that flip under contention and "pass cleanly in isolation
  (`-n0`) and on some full runs, fail on others; classic isolated-flip contention signature".
- The first drives `_grid_seed_search(..., per_call_timeout=5)` — a wall-clock budget, so it
  is contention-sensitive by construction.
- Both passed on an isolated `-n0` re-run (exit 0), and the full ratchet re-run above, on the
  now-quiet machine, is clean.
- Neither test touches anything this task changed (a new test file, this note, and an
  `OUTSTANDING.md` bullet).

This is `[[feedback_serialize_verification_runs]]` behaving exactly as recorded; the
re-run — not the isolated re-run alone — is the evidence relied on here.
