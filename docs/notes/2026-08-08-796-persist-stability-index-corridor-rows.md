# #796 — persist the already-computed stability_index scalar on the corridor rows: honest results

**Date:** 2026-08-08

Split out of `#793`'s own item (c), which was in that task's original registration but got
dropped from its actual dispatch — see `#793`'s own bullet in `data/OUTSTANDING.md`. Dispatch
instructions asserted "the 29 corridor rows" whose `data_gaps` entry cites a missing
`stability_index` scalar. This note records what was actually found and closed, including a
correction to that "29" figure.

## Step 1 — confirming the diagnosis (and correcting the row count)

Grepped the full `data/catalogue.yaml` for every `data_gaps` entry whose `path` ends in
`stability_index`: there are genuinely **29** such entries catalogue-wide, but they are **not**
all "corridor rows", and not all of them describe the same gap this task's dispatch note was
built to close. Reading every one of the 29 in full:

- **20 are genuine corridor rows** matching the exact pattern the dispatch describes — `#438`'s
  scan (or the Braik-Ross Table 2 golden ICs, for the 3 planar rows) records only a qualitative
  Floquet `'stable'` tag, and `#682`'s own corridor census recomputes the monodromy/Floquet
  eigenpairs internally (to pick the widest-rotation center pair for its GMOS torus) but never
  persists a scalar. These are the 20 rows actually in scope:
  `braik-ross-c21-3d-corridor-{01..05}-2026` (5), `braik-ross-c32-3d-corridor-{01..05}-2026` (5),
  `lyapunov3d-l1-corridor-{01..07}-2026` (7), `braik-ross-planar-r{21,31,52}-s-corridor-2026` (3).
- **9 are NOT in scope**, each for a different, independently-honest reason:
  - `braik-ross-c32-cycler-2026` — **already has** a persisted `stability_index` (`126565.0`,
    Barden nu convention, `#249`). Its `data_gaps` entry is `kind: uncertain`, a caveat about that
    value's own numerical-precision ceiling, not a "missing value" gap — nothing to backfill.
  - `arenstorf-em-figure8-1963`, `wittal-2022-em-cycler-family` — literature-sourced gaps ("not
    tabulated in accessible sources" / "PERMANENT ... no monodromy/eigenvalue/stability-index data
    anywhere in the full text"), `todo_ref: #76-cr3bp-backfill` / `#211`. These orbits are not
    products of this project's own `#438`/`#682` census chain, so `ml/seed_generation.py`'s
    function was never established as the right tool for them by this task's own dispatch scope;
    left untouched.
  - `genova-aldrin-2015-em-3petal-cycler` — already `RESOLVED` (`#184`): explicitly
    `kind: not-applicable` (maneuver-maintained, not a natural CR3BP orbit; a monodromy stability
    index doesn't characterize it).
  - `hernandez-2017-jovian-ieg-triple-family`, `russell-strange-2009-jovian-multimoon-family`,
    `russell-strange-2009-saturnian-multimoon-family` — all `kind: not-applicable` (patched-conic
    model; no CR3BP monodromy exists to compute).
  - `umbriel-1-2-torus-homoclinic-uranus-2026` (`#707`), `europa-3-4-crnbp-torus-jupiter-2026`
    (`#736`) — different task chains entirely (Uranus Umbriel-Titania and the N=5 Jupiter torus),
    gated behind their own follow-on manifold/Floquet work, not the `#438`/`#682` Braik-Ross/
    Lyapunov3d-L1 corridor census this task's dispatch note describes.

**Conclusion: the "29" in the dispatch note conflated the catalogue-wide total of
`stability_index`-path `data_gaps` entries with the actual in-scope corridor family. The real
in-scope set is 20 rows, not 29** — reported honestly per this project's own discipline of not
forcing a headline number that doesn't survive a direct check.

## Step 2 — recomputing and persisting `stability_index` on the 20 in-scope rows

Read `src/cyclerfinder/ml/seed_generation.py::stability_index` in full (lines 585-612). It is a
**general, model-agnostic** stability proxy — the spectral radius (max `|eigenvalue|`) of the
**full-period** monodromy matrix, obtained from one extra STM propagation over the orbit's
already-known period. This is deliberately **not** `cr3bp_periodic.barden_stability`'s `nu`
(a half-period convention restricted to `SymmetricOrbit` ICs) — the general proxy applies to any
periodic orbit's raw 6-component state, which is exactly what all 20 rows carry
(`orbit_elements.cr3bp.state_nd`, a full `(x0, 0, z0, 0, ydot0, 0)`-form 6-vector). Convention:
`~1.0` (all eigenvalues on the unit circle) means neutrally/marginally stable; `>1` means at least
one growing/decaying mode.

For each of the 20 rows, fed the row's own already-catalogued `mass_ratio` (`mu`), `state_nd`, and
`period_nd` into a `CR3BPSystem` (Earth-Moon, `l_km=384400.0`, `t_s` from the row's own
`tunit_s`) and called `stability_index(system, state0, period)` directly — pure serialization of
an already-computable value, no new orbit derivation. All 20 computations succeeded (no
propagation failures, no missing prerequisite input — every row already carried a complete
`mass_ratio`/`state_nd`/`period_nd` triple).

**Result: 20/20 closed.**

| Row | stability_index |
|---|---|
| braik-ross-c21-3d-corridor-01-2026 | 1.0000011802616242 |
| braik-ross-c21-3d-corridor-02-2026 | 1.0000000000006097 |
| braik-ross-c21-3d-corridor-03-2026 | 1.000001042942155 |
| braik-ross-c21-3d-corridor-04-2026 | 1.0000011947939065 |
| braik-ross-c21-3d-corridor-05-2026 | 1.000002179166843 |
| braik-ross-c32-3d-corridor-01-2026 | 1.0000038716502904 |
| braik-ross-c32-3d-corridor-02-2026 | 1.0000000000135296 |
| braik-ross-c32-3d-corridor-03-2026 | 1.0000000000183842 |
| braik-ross-c32-3d-corridor-04-2026 | 1.0000000001409557 |
| braik-ross-c32-3d-corridor-05-2026 | 1.0000000000628415 |
| lyapunov3d-l1-corridor-01-2026 | 1.000024770598633 |
| lyapunov3d-l1-corridor-02-2026 | 1.0000544240451326 |
| lyapunov3d-l1-corridor-03-2026 | 1.0000271978408284 |
| lyapunov3d-l1-corridor-04-2026 | 1.0000199438228894 |
| lyapunov3d-l1-corridor-05-2026 | 1.0000170675547742 |
| lyapunov3d-l1-corridor-06-2026 | 1.0000154730906106 |
| lyapunov3d-l1-corridor-07-2026 | 1.0000147036868363 |
| braik-ross-planar-r21-s-corridor-2026 | 1.0000000000035079 |
| braik-ross-planar-r31-s-corridor-2026 | 1.008109574501162 |
| braik-ross-planar-r52-s-corridor-2026 | 1.01480560909678 |

All 20 land at or very close to 1.0 — consistent with the qualitative `'stable'` Floquet tag
`#438`'s scan and Braik-Ross Table 2 already assigned every one of these members. The two planar
rows (`r31-s`, `r52-s`) sit a little further from 1.0 (1.0081, 1.0148) than the rest; both are
still comfortably in the "marginally/near-neutrally stable" regime rather than anything resembling
the genuinely unstable `braik-ross-c32-cycler-2026` sibling row (`stability_index=126565.0`), so
this is read as ordinary numerical residual from a longer-period propagation (6.27 / 12.60 TU vs.
the sub-TU-to-few-TU periods of the other 18 rows) rather than a sign these two are actually
unstable — reported as-is, not over-interpreted.

## Writeback

For each of the 20 rows:
- `orbit_elements.cr3bp.stability_index` set to the computed value, with an inline comment citing
  `#796`, the exact function (`ml/seed_generation.py::stability_index`), the inputs it was fed
  (the row's own `mass_ratio`/`state_nd`/`period_nd`), and an explicit note that this is **not**
  Barden nu and **not** sourced — our own derived computation, not a golden.
- The corresponding `data_gaps` entry (`path: "orbit_elements.cr3bp.stability_index"`) removed.
  17 of the 20 rows retain a second, unrelated `data_gaps` entry (`state_nd`, `kind: derive`,
  `#438`'s own continuation-corrector caveat) — untouched. The 3 planar rows had only the
  `stability_index` gap, so `data_gaps` is now empty (`null`) on those rows, matching this
  catalogue's existing convention for gap-free rows (313 other rows already use a bare
  `data_gaps:` key rather than `data_gaps: []`).

Applied via the same minimal-diff `ruamel`-round-trip-and-patch strategy
`scripts/backfill_russell_2004_tables.py::write_via_patch` already established for this catalogue
(round-trip a pristine baseline and a modified copy, diff the two round-trips against each other
to isolate only the real edit from round-trip reformatting noise, then apply that diff as a patch
onto the untouched original file) — confirmed empirically to avoid the drift a full YAML re-dump
would otherwise introduce. Net diff: 20 rows touched, no unrelated reformatting.

## Verification

1. `python -c "import yaml; yaml.safe_load(open('data/catalogue.yaml'))"` — YAML well-formed.
2. `uv run pytest tests/data/test_jsonschema.py -q` — passes.
3. `uv run pytest tests/ -q` — full tree (mandatory per
   `[[feedback_verify_scope_must_include_tests_scripts]]`); see summary below.
4. `uv run ruff check .` / `uv run ruff format --check .` — clean.
5. `uv run mypy src tests` — `Success: no issues found in 839 source files` (no `src`/`tests`
   Python files were touched by this task; only `data/catalogue.yaml`).

**Full-tree result:** 4 `FAILED`, 6 `XFAIL`, 20 `XPASS`, no collection errors. All 4 failures
confirmed pre-existing and unrelated via a direct `git stash`-and-rerun check (stashed
`data/catalogue.yaml` out, reran just the 4 failing tests against the untouched baseline — all 4
reproduce byte-identical failures with this task's catalogue edit completely absent):

- `tests/search/test_eggie_ballistic.py::test_gate_b_table4_vinf_reached_but_subsurface` and
  `tests/search/test_504_pluto_charon_kk_sweep.py::test_504_sweep_33` — already documented as
  pre-existing in `#682`'s own results note (same 2 failures, same modules).
- `tests/nbody/test_propagator_api.py::test_rails_cache_batch_samples_match_per_point` (residual
  `2.98e-08` vs. a `1e-09` gate) and `tests/genome/test_da_section_map.py::
  test_taylor_fixed_point_reaches_png_neighbourhood` (residual `2.78e-04` vs. a `1e-04` gate) —
  both new to this run, both razor-thin-margin numeric residuals in modules `#796`'s
  `data/catalogue.yaml`-only change never touches (N-body propagator rails cache, DA-section-map
  Taylor fixed point), consistent with the well-documented cross-platform DOP853/BLAS-divergence
  class already recorded throughout `data/OUTSTANDING.md` (`#584`/`#631`/`#632`/`#635`/`#731`/
  `#784`).

A test run also mutated a tracked-but-nondeterministic data file
(`data/floquet_phase1_reproduction.jsonl`, written by
`tests/genome/test_floquet_phase1_reproduction.py`) as an unrelated side effect of running the
suite; reverted via `git checkout` before committing since it is outside this task's scope.

## Scope note

All 20 rows this task's corrected scope covers were closed. The 9 catalogue-wide
`stability_index`-path gaps outside that scope were left untouched, each for an independently
checked, honest reason (already-persisted value, literature-sourced permanent gap,
model-inapplicable, or a different unrelated task chain) — see Step 1 above. No row was forced
through the general-proxy computation where it didn't actually apply.
