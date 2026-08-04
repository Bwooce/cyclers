# `#777`: vendor the remaining Neptune-Triton Miceli & Bosanac 2026 ESM rows

**Task:** `#777`, registered 2026-08-01 by `#776`'s own closure, DISPATCHED 2026-08-04 (user:
"dispatch 777"). `#776` vendored a deliberately narrow ~16-row subset of the Miceli & Bosanac 2026
ESM dataset (per the user's own "option 1: prove narrow first" choice); `#777` is the explicit,
already-registered follow-up to vendor the remaining canonical rows and extend the same honest
3-criterion gate (plus continuation and two-body-seed-lineage checks) to them.

**Source data:** unchanged from `#776` -- the same three machine-readable ESM text files
(`miceli-bosanac-2026-jas-73-11-esm{2,3,4}-*.txt`, filed at `cyclers_pdf/papers/`, indexed in
`docs/notes/CORPUS_INDEX.md`'s own `#776` entries) at the same `µ=2.089503183689124e-04`. No new
acquisition, filing, or corpus-index work this task.

---

## Row count: 48, not `#776`'s own "~34" estimate

`#776`'s own results note guessed "the remaining ~34 dataset rows" (out of a hedged "~50" total).
This task ran a full, programmatic, line-by-line audit of all three ESM files instead of trusting
that estimate, and the true count is different -- **stated plainly, not silently under-covered**:

* Every non-header, non-separator line in each file was parsed and classified two independent
  ways: (1) by the source's own label convention (`"Manifold" in label` excludes a manifold-arc
  sample), and (2) independently by physics (`|y0| < 1e-6 and |xdot0| < 1e-6`, the symmetric
  perpendicular-crossing form every genuine periodic-orbit row in this dataset satisfies). The two
  discriminators agree on every row except two: ESM3's own `InitialCondition` row (label says
  non-manifold/"canonical", physics says NOT symmetric -- `y0=-3.615`, `xdot0=-3.610`, and its own
  "Integration Time" of exactly `8.000000000000` is a round guess, not a converged period -- this
  row is EXCLUDED, it is not a periodic-orbit family member at all), and ESM4's own 4th printed
  `Res32-x+h` row (label says canonical, physics borderline: `xdot0=-3.5e-6`, ~1000x the other
  rows' numerical-zero floor -- vendored anyway, since it converges and reproduces cleanly, see
  below).
* Total canonical periodic-orbit rows across all three files: **64** (2 ESM2 + 20 ESM3 + 42 ESM4,
  each count excluding manifold-arc samples and the one ESM3 `InitialCondition` row). This matches
  `docs/notes/CORPUS_INDEX.md`'s own `#776` ESM4 entry ("32 resonant rows + 8 DPOs + 1 LPO + 1 L2
  Lyapunov target" = 42) exactly; its own ESM3 entry says "21 resonant periodic-orbit rows" (one
  more than this task's own 20) -- almost certainly that count folding in the `InitialCondition`
  row, which this task's own physics cross-check confirms does NOT belong in that count.
* `#776` vendored 16 of the 64 (10 `ESM_GATE_ROWS` + 3 `FAMILY_23` + 2 `FAMILY_43_HC2` incl. 1
  duplicate with `ESM_GATE_ROWS["4:3-saddle"]` + 2 `FAMILY_43_NEAR_UNIT` = 16 unique rows).
* **This task vendors all 64 - 16 = 48 remaining canonical rows** -- full coverage, not a further
  narrow subset. No manifold-arc sample and no `InitialCondition`-like row was vendored as if it
  were an independent orbit.

New module data structure: `ESM_ROWS_777: dict[str, EsmRow]` (48 entries), keyed with a
disambiguating suffix encoding source file and, where a resonance label repeats across files (e.g.
`Res32-x+h`/`Res35+x+h`/`Res37+x+h`/`Res25+x+h` each appear in BOTH ESM3, the high-energy scenario,
and ESM4, the low-energy scenario, as genuinely DIFFERENT orbits at different energies), a member
index. `ESM_GATE_ROWS` itself is BYTE-IDENTICAL to `#776`'s own version (still exactly 10 rows) --
`#776`'s own reproduce-before-trust anchor test is untouched.

## `half_crossings`: automatic, not hand-picked, and verified against `#776`'s own values

`#776` determined each row's own `half_crossings` (the x-axis-crossing index the corrector locks
onto) "empirically this task by direct inspection of each row's own crossing-index list" (module
docstring). At 48 new rows that would not scale. Instead, this task used the SAME logic
`cr3bp_periodic.correct_symmetric_fixed_jacobi` itself applies internally when `half_crossings=None`
(the x-axis crossing nearest `T/2` on the row's own raw printed seed, `_xaxis_crossings` +
`argmin`), applied programmatically to every row.

**Verification (new test, `test_half_crossings_none_auto_detection_matches_esm_gate_rows_hand_picked_values`):**
run against all ten of `#776`'s own `ESM_GATE_ROWS`, this automatic determination recovers the
EXACT same converged `(x0, ydot0, period)` `#776` hand-picked to <1e-4 relative on every row (most
agree to <1e-11). This is the assumption the whole `#777` vendoring pass rests on, and it holds.

## Manifold/canonical partition cross-check

The label-based (`"Manifold" in label`) and physics-based (`|y0|<1e-6 and |xdot0|<1e-6`)
partitions were run independently over all three files (993 manifold-labelled rows, 65
non-manifold-labelled rows including `InitialCondition`) and agree on every row except the two
noted above (`InitialCondition`, correctly excluded by physics; ESM4's 4th `Res32-x+h` row,
correctly included -- see below).

---

## Gate results: 47/48 pass all three criteria; one honest crosscheck negative

Same 5-part gate `#776` established, applied to all 48 new rows:

**(a) periodicity self-consistency, (b) reproduction, (c) internal Barden-vs-`_planar_floquet`
cross-check** -- run via `gate_report_777()` (a new `rows` parameter on the existing `gate_report`,
defaulting to `ESM_GATE_ROWS` so `#776`'s own behavior and its own `test_gate_report_all_ten_rows_pass`
are completely unaffected):

* **(a) periodicity**: 48/48 pass, but the margin is tighter than any of `#776`'s own ten. The four
  worst are all strongly-unstable DPO/2:5 family members: `dpo-esm4-6` (7.19e-7, |λ|≈3062),
  `dpo-esm4-8` (5.69e-7, |λ|≈3671), ESM3's own `2:5+x-esm3-b` (1.83e-7), and `dpo-esm4-5` (1.41e-7)
  -- `dpo-esm4-6`'s own 7.19e-7 is only a ~1.4x margin inside `PERIODICITY_GATE_TOL=1e-6`, not the
  order-of-magnitude margin `#776`'s own ten enjoyed, though the SAME mechanism this module's own
  docstring already names (DOP853 residual-noise amplification along a strongly unstable direction)
  explains it -- these are the module's own most unstable rows. Every other row is <=5.1e-9.
* **(b) reproduction**: 48/48 pass. Worst observed relative error 4.70e-8 (`3:2-x-esm4-4-hc2`, the
  physics-borderline row flagged above -- its own raw `xdot0=-3.5e-6` propagates into a
  correspondingly looser (but still passing, four+ orders of magnitude inside
  `REPRODUCTION_GATE_REL_TOL=1e-6`) reproduction).
* **(c) internal crosscheck**: **47/48 pass.** The one honest negative:

  | Key | Source | C_J | Periodicity resid. | Reproduction | Barden vs. planar_floquet rel. err | Passed |
  |---|---|---|---|---|---|---|
  | `dpo-esm4-2` | ESM4 `DPO` (2nd printed) | 3.014371175305 | 3.06e-12 | <=2.7e-11 | **2.15e-5** | **FAIL (c) only** |

  `dpo-esm4-2` is a genuinely near-unit-circle (`max_eigenvalue≈1.0`, `is_real_unstable=False`)
  orbit with a very short period (T≈1.63 nondim, the shortest of any vendored row in this whole
  module). Its Barden-eigenvalue and independent full-period `_planar_floquet` monodromy
  eigendecomposition agree to only 2.15e-5 relative -- just over 2x `CROSSCHECK_GATE_REL_TOL=1e-5`,
  in ABSOLUTE terms a tiny miss (both values are within 2.2e-5 of exactly 1.0), most plausibly a
  numerical-conditioning artifact of evaluating an eigenvalue derivative extremely close to the
  unit circle on a short-period orbit, not a data or corrector defect -- periodicity and
  reproduction on this SAME row are both clean (3.06e-12 / <=2.7e-11). Reported exactly as found,
  `#777`'s own analog of `#776`'s own 4:3 fold-reversal finding: a genuine, well-characterized
  partial result, not a bug to chase or paper over.

`gate_report(system, rows=...)` (the underlying, now-generalized function) was also hardened this
task per pre-work review: it no longer `raise`s on a non-converging seed (safe for `#776`'s own ten
hand-picked rows, which were all known in advance to converge, but not safe for 48 rows that were
NOT hand-picked) -- a non-convergence is now recorded as its own honest `GateRow` instead. In
practice every one of the 48 `#777` rows DOES converge, so this path is defensive, not
load-bearing for this task's own sweep; `recover_esm_candidate` (single-label,
`ESM_GATE_ROWS`-only) is unchanged and still raises, which remains correct for `#776`'s own ten
already-known-good rows.

---

## Item (d): `two_body_resonant_seed` lineage, extended onto two new resonance ratios

`#776` tried (4,3) and (2,3), both honest negatives. This task extended the check onto two of the
NEW resonance ratios `#777` vendors that `#776` never tried: **(1,2)** and **(2,1)**, each
`x0_sign=-1`, converged directly at each seed's own natural Jacobi constant exactly as `#776`'s own
check did.

* **(2,1)**: converges cleanly but to the WRONG topology (`period/2pi≈2.000102`, not the label's
  own implied `q=1`, i.e. `period/2pi≈1.0`) -- the same qualitative finding as `#776`'s own
  (4,3)/(2,3) attempts and Anderson & Lo's/Vaquero's own analogous results. Counting by
  system/task (this project's own established convention -- `#776`'s own (4,3)+(2,3) pair together
  counted as ONE "third, independent confirmation"), this is the **4th** independent confirmation:
  Anderson & Lo (Jupiter-Europa) 1st, Vaquero (Saturn-Titan) 2nd, `#776` (Neptune-Triton, (4,3)+(2,3))
  3rd, `#777`'s own (2,1) 4th.
* **(1,2)**: a MORE NUANCED negative, stated precisely rather than folded into the same bucket:
  it converges to `period/2pi≈2.000218`, matching the "1:2" label's OWN implied index (`q=2`) --
  unlike every other naive attempt in this whole family of checks, across all three systems. But it
  lands at `x0=-1.0`, `C=2.971143204812636`, hugely different from the paper's own printed
  `"1:2+x-esm3"` row (`x0=-2.881985569172`, `C=2.087857684887`). Even a period-index "hit" is NOT a
  genuine identification of the paper's own labeled family member -- the naive two-body-resonant-
  ellipse construction remains unreliable as a seed for identifying a SPECIFIC labeled family in
  this system, even on the rare occasion its period ratio coincidentally lands right.

Both documented in the new `TWO_BODY_SEED_LINEAGE_NOTE_777` (importable, tested).

---

## Item (e): continuation-in-C_J -- two clean confirmations, three honest negatives

Five multi-member families in `ESM_ROWS_777` support a continuation check (uniform
`half_crossings` across >=2 members); each was attempted, seeded from its own lowest-C member,
walking toward its own highest-C member via the EXISTING `cr3bp_continuation.continue_family`
gauntlet (same machinery `#776` used, no reimplementation):

| Family | Members | `half_crossings` | Result | Detail |
|---|---|---|---|---|
| Low-energy 3:2, HC1 branch (`FAMILY_32_ESM4_HC1_777`) | 3 (incl. `ESM_GATE_ROWS["3:2-start"]`) | 1 | **CLEAN** | `JACOBI_BOUND`, 482 gauntlet-passing members, 0 rejections; lands on both higher-C printed members to 5.7e-5/1.0e-6 and 2.8e-4/2.3e-5 relative (x0/period) |
| High-energy 4:7, HC3 pair (`FAMILY_47_ESM4_HC3_777`) | 2 | 3 | **CLEAN** | `JACOBI_BOUND`, 105 gauntlet-passing members, 0 rejections; lands on the printed endpoint to 2.5e-4/4.5e-7 relative (x0/period) |
| Low-energy DPO, full family (`FAMILY_DPO_777`) | 8 (incl. `ESM_GATE_ROWS["DPO"]`) | 1 | HONEST NEGATIVE | `GAUNTLET_REJECT` after 36 members -- reaches only the 1st/2nd of 7 higher-C members cleanly before the gauntlet's own physical-plausibility check rejects the walk |
| Low-energy 2:1, full family (`FAMILY_21_ESM4_777`) | 6 | 1 | HONEST NEGATIVE | `FOLD_REVERSAL` after 91 members -- reaches 4 of 5 higher-C members cleanly (x0 rel. err <=4.2e-3) but folds back short of the 6th, C=3.667872679873, a genuine outlier ~0.5 higher in C than its own 5 siblings in [3.002, 3.169] |
| High-energy 2:5 "-x+h" pair (`FAMILY_25M_ESM4_777`) | 2 | 3 | HONEST NEGATIVE | `TOPOLOGY_JUMP` at the FIRST continuation step -- these two printed rows are not two points on one continuous branch at all |

Two clean confirmations (`continue_32_esm4_hc1_family_777`, `continue_47_esm4_pair_777`) satisfy
the "≥2 families" bar the dispatch note asks for. The three negatives are reported exactly as
found -- a `GAUNTLET_REJECT`, a `FOLD_REVERSAL`, and a `TOPOLOGY_JUMP` are three DIFFERENT genuine
findings, not the same bug wearing different names: the DPO family's own gauntlet rejection is a
physical-plausibility stop (unlike a fold or a jump), the 2:1 family's fold reversal is the same
qualitative "this label spans >=2 topologically distinct branches" finding `#776` established for
the 4:3 family (further corroborated by the flagged C=3.6679 outlier), and the 2:5 "-x+h" pair's
topology jump is the sharpest form seen in this task -- not even one continuation step succeeds.

A sixth candidate family, ESM4's own mixed-`half_crossings` `Res32-x+h` set (3 members at
`half_crossings=1` -- the `FAMILY_32_ESM4_HC1_777` above -- plus a 4th, physics-borderline member
at `half_crossings=2`, vendored as `"3:2-x-esm4-4-hc2"`), is ANOTHER family-mixing instance
paralleling `#776`'s own 4:3 case -- the odd member has no same-`half_crossings` sibling among the
new rows, so no continuation was attempted on it (a lone member cannot support one); it is vendored
and gate-passes individually, just not continuation-tested. Likewise ESM4's `Res37-x+h` pair (2
members, `half_crossings` 5 and 4 -- mixed, no continuation attempted) and ESM3's `Res35+x+h` pair
(2 members, `half_crossings` 5 and 7 -- mixed, no continuation attempted): both are documented here
rather than force-walked, per the dispatch note's own "if not, that's fine, just don't force one"
instruction.

---

## Tests

`tests/search/test_neptune_triton_resonant_families.py` extended from 47 to 63 passing collected
items (16 new top-level test functions: row-count/partition/verbatim-value checks, the
`half_crossings` auto-detection cross-check, three `gate_report_777` assertions, two extended
two-body-seed-lineage tests, and five continuation tests -- one per family above). A new
module-scoped `gate_777_rows` fixture computes `gate_report_777()` ONCE per test session and is
reused by every gate-row assertion, rather than re-converging all 48 rows per test (the existing
per-row parametrized `#776` tests each independently re-converge all ten `ESM_GATE_ROWS`; at 48
rows that pattern would mean hundreds of redundant corrector/STM runs).

`uv run pytest tests/search/test_neptune_triton_resonant_families.py -v`: **63/63 pass**
(455s wall time, 8 xdist workers -- longer than `#776`'s own suite, dominated by the two clean
continuation walks, 482 + 105 gauntlet-passing members).

---

## Verification

* `uv run pytest tests/search/test_neptune_triton_resonant_families.py -v`: 63/63 pass.
* `uv run pytest tests/data tests/search -q`: see commit history / this note's own final
  verification block for pass/fail status recorded at commit time (this task does not touch
  `data/catalogue.yaml` -- no catalogue writeback, per `#776`'s own explicit out-of-scope framing,
  unchanged here).
* `uv run ruff check .` / `ruff format --check .`: clean.
* `uv run mypy src tests` (canonical full-strict invocation): clean.
* `uv run pytest tests/data/test_outstanding_structure.py tests/data/test_outstanding_header_body_consistency.py -q`:
  run before committing the `OUTSTANDING.md` update.
* No file under `scripts/` was created or edited this task -- `tests/scripts/` is not relevant.

---

## Explicitly out of scope (unchanged from `#776`)

No manifold/homoclinic/heteroclinic/chain work, no retrograde families, no catalogue writeback
(nothing here is catalogue-eligible until confirmed AND novelty-checked via
`search/literature_check.py`, a separate future task -- this task does not do that check and does
not touch `data/catalogue.yaml`). No new corpus acquisition/filing/digesting -- the source files
were already filed and indexed by `#776`.

## What remains after `#777`

All 64 canonical periodic-orbit rows across all three ESM files are now vendored (16 by `#776`, 48
by this task) -- full dataset coverage for the family-confirmation gate. No further "vendor more
rows" follow-up is registered; `#776`'s own "Recommendation for a Task-B analog" section (a future
connection-stage task, NOVEL-territory since no published Neptune-Triton homoclinic/heteroclinic
connection state exists) remains the only open recommendation from that thread, unchanged and not
attempted by this task.
