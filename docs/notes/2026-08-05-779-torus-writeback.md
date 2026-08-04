# `#779`: writeback of `#766`'s C=3.0041 homoclinic self-connection as seed-lineage provenance

**Task:** `#779`, implementing `#770`'s user-approved recommendation (option b) without
any deviation: keep `orbit_class: quasi_periodic_torus` and `validation_level: V1`
unchanged on `europa-3-4-crnbp-torus-jupiter-2026`; add ONE new optional, non-gating
sub-block `crnbp_provenance.torus.seed_orbit_homoclinic` recording `#761`/`#766`'s
homoclinic self-connection of the row's cited seed periodic orbit; plus three small
companion edits (`notes`, two `data_gaps` entries). See
`docs/notes/2026-07-31-770-torus-connection-writeback-design.md` for the full argument
this task implements verbatim.

## What was added

`data/catalogue.yaml`, row `europa-3-4-crnbp-torus-jupiter-2026`:

1. `crnbp_provenance.torus.seed_orbit_homoclinic` — a new sub-block, sibling to the
   existing `seed_lineage` and `radau_cross_check` fields, following the `#738`
   `radau_cross_check` precedent (row-side, additive, no schema bump — both
   `crnbp_provenance` and its `torus{}` sub-object carry `additionalProperties: true`,
   confirmed by reading `data/catalogue.schema.json` directly). Carries `object`,
   `model`, `orbit{}`, `connection{}`, `evidence_class`, `caveats`, `method` — the exact
   shape `#770` §4 drafted, prose fields copied essentially verbatim.
2. `data_gaps` `orbit_elements.cr3bp.stability_index`: appended one sentence noting the
   seed orbit's own CR3BP saddle character is now established (`#761`), so the gap is
   specifically the N=5 torus object's own Floquet/manifold work — `kind: unknown` and
   the original `#724` qualifier text are unchanged/verbatim.
3. `data_gaps` `orbit_elements.cr3bp.jacobi_constant`: added a pointer to
   `torus.seed_orbit_homoclinic.orbit.jacobi_constant` for the seed's own exact C — the
   original "no connection exists for this object" sentence is unchanged/verbatim
   (still true, object-scoped).
4. `notes`: appended a "Seed-lineage standing (#761/#766, writeback #770/#779)"
   paragraph summarizing the new evidence and its non-gating status, and extended the
   "Discovery + verification chain" line with
   `-> #750 -> #761 -> #766 -> #770 (design) -> #779 (writeback)`.

**Not touched**: `orbit_class` (`quasi_periodic_torus`), `validation_level` (`V1`),
`epoch_locked`, `n_returns`, `model_assumption`, and every other existing field's
value. No `ccr4bp_provenance.connection{}` block was added (per `#770` §2.2/§4, that
block is scoped to a different pipeline/class and roughly half its fields have no
referent here).

## Full-precision re-derivation (mandatory per the dispatch, not just copying the draft's rounded placeholders)

Ran the producing code fresh this task (`src/cyclerfinder/search/jovian_resonant_families.py:continue_34lo_to_kumar_c`,
`src/cyclerfinder/search/jovian_resonant_connections.py:{build_34lo_kumar_c_node, find_homoclinic, homoclinic_reapproach_check}`,
`cyclerfinder.genome.heteroclinic_cycle.{correct_connection, assemble_cycle, crosscheck_cycle}`),
exactly mirroring `tests/search/test_jovian_resonant_connections.py`'s own known-hit
invocations (`_KNOWN_HIT_766`, `max_time_factor=8.0`, Radau `rtol=atol=1e-11`).

**Orbit block** (`continue_34lo_to_kumar_c` endpoint) — bit-identical to `#766`'s own
recorded values:

```
x0      = -1.385248445624164     (identical to #766's -1.3852484456241640)
ydot0   =  0.5988394002678391    (identical)
period  = 25.312119648766764     (identical)
lambda  = 54.589750588953734     (identical)
```

No correction needed here — the design note's draft was already full precision on
this sub-block (only `connection.t_u_tu` and `connection.radau_crosscheck` were
marked `TBD`).

**Connection block** — re-derived via `correct_connection`, seeded at `#766`'s own
recorded `(tau_u0, tau_s0) = (10.72913431175392, 14.582984371438714)`, `tol=1e-9`
(the test suite's own invocation):

- At `tol=1e-9` the corrector accepts the seed **unmoved** — its residual (`1.97e-10`)
  was already inside `tol`, so zero Newton iterations ran. `tau_u`/`tau_s` in the
  written block are therefore `#766`'s own recorded values carried forward, **not** an
  independent re-solve at this tolerance. Verified directly: re-running at
  `tol=1e-12` DOES move the corrector — `tau_u -> 10.729134321384189`,
  `tau_s -> 14.582984370878199`, `residual -> 2.60e-12` — confirming the seed is a
  genuine near-root (not a no-op artifact), just not independently re-derived at
  `tol=1e-9`.
- `crossing_x`, `crossing_xdot`, `newton_residual`, `ghost_distance` **are** freshly
  integrated from that seed (not echoed):
  `crossing_x=-1.4220714951697728`, `crossing_xdot=-2.016580027963677e-10`,
  `newton_residual=1.9723783076056544e-10`, `ghost_distance=0.03682304954560878` —
  all bit-identical to `#766`'s recorded values (as expected: same seed, same
  deterministic Newton path, same production code).
- **Independent corroboration** (not merely re-affirming the handed-in seed): a
  SEPARATE, un-seeded `find_homoclinic(rank_by_residual=True)` grid-scan (its own
  internal `scan_n=12` seeding, no `tau_u0`/`tau_s0` passed) re-landed on the same
  root from scratch: `tau_u=10.729134301732879`, `tau_s=14.582984363185206`
  (~1e-8 from `#766`'s), `crossing_x=-1.4220714951677023`,
  `residual=1.1606645550029221e-10`, `ghost_distance=0.03682304954353821`,
  Radau crosscheck `2.4064910829886726e-08` — agreeing to 9-10 significant digits on
  `ghost_distance`/`crossing_x` with `#766`'s value, from a genuinely independent
  search path. Both re-derivations are recorded in the catalogue block's own comments
  for full transparency.
- `t_u_tu` (full precision, replacing the design note's `134.269 # TBD` placeholder):
  `134.26929450801416` (`homoclinic_reapproach_check`'s own `t_u`); `|t_s| =
  134.26929354447384`, agreeing to `9.6e-7` absolute (`7e-9` relative) — the expected
  forward/backward integration-direction asymmetry, not a discrepancy.
- `radau_crosscheck` (full precision, replacing the design note's `2.42e-8 # TBD`
  placeholder): `2.4222186781124212e-08` (via `assemble_cycle`/`crosscheck_cycle`,
  `method="Radau"`, `rtol=atol=1e-11`, `max_time_factor=8.0` — the exact test-suite
  invocation). Matches the design note's rounded value.
- `backward_reapproach`/`forward_reapproach` (full precision): `5.638445571711433e-08`
  and `2.2523590122573838e-05` — match the design note's rounded values
  (`5.638e-08`, `2.252e-5`) to displayed precision.
- `mirror_pair_note`: unchanged from the draft (not marked `TBD`, no re-derivation
  needed for this task).

**Verdict: no re-derived value differed from the design note's rounded/placeholder
value by more than expected floating-point noise.** Nothing was silently overwritten;
the one genuinely subtle point (that `tau_u`/`tau_s` in the written block are echoed
rather than independently re-solved at `tol=1e-9`) is documented explicitly in the
catalogue block's own inline comment, plus the independent grid-scan corroboration
that closes that gap.

Also confirmed `tests/search/test_jovian_resonant_connections.py` still collects
27 items (matches the `method:` field's claim) via
`uv run pytest tests/search/test_jovian_resonant_connections.py --collect-only -q`.

## Verification

1. `uv run python -c "import yaml; yaml.safe_load(open('data/catalogue.yaml'))"` —
   clean.
2. `uv run pytest tests/data/test_jsonschema.py -q` — 5/5 pass (schema validation).
3. `uv run pytest tests/ -q` (full tree, mandatory per this project's own
   catalogue-edit rule) — run twice (once before, once after a comment-only fix to
   the new block following an advisor review). Both runs showed the SAME small set
   of failures, and a stash-based baseline check (`git stash push -- data/catalogue.yaml`,
   re-run the failing tests in isolation, `git stash pop`) confirmed **all of them
   pre-exist independent of this change**:
   - 4 tests fail identically with or without the `#779` catalogue edit, in isolation,
     confirming genuine pre-existing failures unrelated to this task:
     `test_da_section_map.py::test_taylor_fixed_point_reaches_png_neighbourhood`,
     `test_propagator_api.py::test_rails_cache_batch_samples_match_per_point`,
     `test_eggie_ballistic.py::test_gate_b_table4_vinf_reached_but_subsurface`,
     `test_504_pluto_charon_kk_sweep.py::test_504_sweep_33`.
   - The remaining failures/errors seen only under full-suite 8-way parallel load
     (varying between the two runs: SPICE/nbody/ephemeris tests, one SPICE
     kernel-pool-count test) all PASSED when re-run in isolation — CPU-contention
     flakes under `pytest-xdist`, consistent with this project's own documented
     class of flake (`feedback_serialize_verification_runs`), sharpened here by a
     concurrent agent actively running Neptune-Triton work in the same tree during
     this session (confirmed via `git status` showing unrelated tracked-file
     modifications appear/disappear around the stash operations — not touched by
     this task).
   - Zero new failures traceable to the `#779` catalogue edit.
4. `uv run ruff check .` — all checks passed. `uv run ruff format --check .` — one
   file flagged, `src/cyclerfinder/core/ccr4bp_titan_hyperion.py` — this is an
   UNTRACKED file from a different, concurrently-running task (present at session
   start, not touched by `#779`); not a `#779` regression.
5. `uv run mypy src tests` (full strict run) — clean, "Success: no issues found in
   829 source files".
6. `git diff data/catalogue.yaml` reviewed directly: 70 insertions / 3 deletions
   (later folded into a slightly larger diff after the comment fixes above), exactly
   the one new sub-block plus the three companion edits — nothing else touched.
   `orbit_class`, `validation_level`, `epoch_locked`, `n_returns`,
   `model_assumption` confirmed unchanged by direct inspection of the diff (no lines
   touching those keys appear).

## Confirmation of scope discipline

- `orbit_class` stays `quasi_periodic_torus`. **Not touched.**
- `validation_level` stays `V1` (the `#738` Radau-cross-check basis). **Not touched.**
- `epoch_locked`, `n_returns`, `model_assumption`: **not touched.**
- No `ccr4bp_provenance.connection{}` block added (explicitly rejected by `#770`
  §2.2/§4).
- No schema change (`data/catalogue.schema.json` untouched — verified
  `additionalProperties: true` on both `crnbp_provenance` and its `torus{}`
  sub-object directly from the schema file before writing the block).
