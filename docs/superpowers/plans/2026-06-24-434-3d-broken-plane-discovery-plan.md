# #434 3D / broken-plane cycler discovery — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: superpowers:subagent-driven-development. Steps use `- [ ]` checkboxes.

**Goal:** Discover novel out-of-plane (z≠0) cycler families by lifting the catalogue's planar Aldrin/Braik-Ross cyclers into z0 and continuing them with the already-built 3D genome, screening survivors against the known halo/NRHO/vertical libration families. The #287 spike already proved a genuine 3D extension of the planar Braik-Ross (1,1) Earth-Moon cycler exists (z0≈−0.241, 92,567 km out of plane, ~80 converged members) — this campaign maps the family systematically.

**Architecture:** The 3D infrastructure is ~90% built and stays UNCHANGED: Phase-1 corrector `correct_general_periodic_3d` (`src/cyclerfinder/search/cr3bp_general_periodic_3d.py:285`, compound corrector+independent-Radau-closure gate) and Phase-2 pseudo-arclength tracer `continue_general_3d_family` (`src/cyclerfinder/search/cr3bp_3d_family_tracer.py:602`). New work: (1) a 3D vertical-Lyapunov/halo seed generator, (2) a z-oscillation topology count to distinguish out-of-plane structure, (3) a discovery sweep driver that lifts planar catalogue cyclers into z0, (4) literature-check + validation widening to adjudicate 3D candidates against the spatial-CR3BP corpus. Report-only first; catalogue admission of any genuine survivor is a separate follow-on (the Umbriel template).

**Honest likelihood:** MEDIUM-HIGH for novel structure, but heavily mixed with rediscovery — halos/NRHOs/vertical-Lyapunov/butterfly are exhaustively published (Howell 1984, Folta-Bosanac-Cox-Howell 2017, JPL 3BP catalog). The genuine frontier is 3D Aldrin-class repeating-encounter cyclers lifted from planar roots, which are NOT systematically tabulated. The literature gate is mandatory before any novelty claim.

**Conventions:** work on main; `uv run` ruff + mypy before commit; no Co-Authored-By; pathspec commits; imports at top; goldens trace to a published source OR a closure self-check (never to a value our own code computed); subagents finish through commit, no self-spawned reviewers.

---

## Task 0 (controller, gate): re-confirmation spike

- [ ] Before any new code: run a 4–6h-equivalent natural-z0 continuation on the Braik-Ross C11a (1,1) Earth-Moon seed using the EXISTING `continue_general_3d_family(..., continuation="natural_x0"|"natural_T")` from a z0-perturbed seed, confirming closure < 1e-9 over a z0 range reproducing the #287 spike (`data/spike_287*.jsonl` if present). This is a GO/NO-GO: if the existing tracer cannot reproduce the spike's ~80 members, STOP and debug the tracer before building Tasks 1–4. Record the confirmation in the eventual verdict.

## Task 1: 3D vertical-Lyapunov / halo seed generator

**Files:** Modify `src/cyclerfinder/search/cr3bp_seed_generator.py`; test `tests/search/test_cr3bp_seed_generator.py`

- [ ] Step 1: Failing test for `lyapunov_seed_3d(system, *, point="L1", amplitude_z=0.02, amplitude_x=1e-3) -> (state0(6,), period)`. Extract the z-direction (vertical) linear mode at the collinear point: the out-of-plane frequency is `omega_z = sqrt(c2)` (Koon-Lo-Marsden-Ross 2011 §2.5; the z-equation decouples linearly as `z'' + c2 z = 0`), giving a small-amplitude vertical IC `[x_L, 0, amplitude_z, 0, vy0, 0]` refined by `correct_general_periodic_3d(..., free_vars=FREE_VARS_SYMMETRIC_TULIP, residual_indices=RESIDUAL_PERPENDICULAR_HALF_PERIOD, is_half_period_residual=True)`. GOLDEN (sourced): the Earth-Moon (Koblick μ) L1 vertical/halo seed converges with non-zero z0 (`abs(state0[2]) > 0.01`, not a planar collapse) and closes (`orbit.converged`), period in the documented vertical-Lyapunov band. Assert closure + non-planarity (NOT an unsourced exact IC).

- [ ] Step 2: Run → FAIL. Step 3: Implement (`omega_z=sqrt(c2)`; reuse `_collinear_c2`/`lagrange_collinear_x`; amplitude ladder if convergence is μ-fragile, mirroring `lyapunov_seed`). Step 4: Run → PASS.
- [ ] Step 5: ruff + mypy; commit `search/#434: 3D vertical-Lyapunov/halo seed generator at arbitrary mu`.

## Task 2: z-oscillation topology classifier

**Files:** Modify `src/cyclerfinder/search/binary_star_search.py` (the home of `winding_topology`); test alongside existing topology tests.

- [ ] Step 1: Failing test for `z_oscillation_count(z_series, *, z_center=0.0) -> int` (counts equatorial-plane crossings / z sign-changes over one period) and a `Topology3D(k1, k2, k_z)` extension (add `k_z` field; keep planar `(k1,k2)` from the existing `winding_topology`). GOLDEN (self-check): a propagated planar orbit (z≡0) returns `k_z == 0`; a synthetic `z_series = sin(2π t)` over one period returns `k_z == 2` (two sign changes) — pins the counting convention without an unsourced value.
- [ ] Step 2: Run → FAIL. Step 3: Implement `z_oscillation_count` (count sign changes of `z_series - z_center`, robust to endpoints) + `topology_3d(mu, state0, period)` that propagates, calls `winding_topology` for `(k1,k2)`, and counts `k_z`. Step 4: Run → PASS.
- [ ] Step 5: ruff + mypy; commit `search/#434: z-oscillation topology count (k_z) for 3D orbits`.

## Task 3: 3D broken-plane discovery sweep driver

**Files:** Create `scripts/scan_434_3d_broken_plane_em.py`

**SPIKE-CONFIRMED CAVEAT (Task 0 finding 2026-06-24):** small-z0 lifts of a planar root COLLAPSE back to the planar manifold (|z0|→<1e-14) — the 3D branch only locks in from a seed already well out-of-plane (the #287 spike's C11a 3D branch locks near z0≈−0.24, reached via a dedicated push). So a naive "perturb planar root by tiny z0 and correct" sweep mostly yields planar re-collapses. Two seed routes that DO lock onto 3D structure:
  (i) **vertical-Lyapunov generator (Task 1)** — `lyapunov_seed_3d` produces genuinely out-of-plane ICs at the collinear point; continue each in z0-amplitude / Jacobi.
  (ii) **z0-amplitude lock per planar root** — for each planar Aldrin/Braik-Ross root, step z0 up from a MODERATE start (e.g. |z0|≥0.05) and let `correct_general_periodic_3d` (full-asymmetric `FREE_VARS_FULL_ASYMMETRIC` if symmetric-tulip collapses) find the 3D branch; record which roots HAVE a 3D extension (lock) vs collapse-only.

- [ ] Build the sweep using BOTH seed routes: (i) generate `lyapunov_seed_3d` families at Earth-Moon L1/L2 and continue them; (ii) for each planar Aldrin/Braik-Ross catalogue root, attempt the z0-amplitude lock from |z0|∈{0.05,0.10,0.15,0.20,0.24} (record lock-vs-collapse). For each converged 3D seed run `continue_general_3d_family(continuation="pseudo_arclength", step=0.01, n_steps_max=100, direction="both")`. Log every member to `data/scan_434_3d_broken_plane_em.jsonl`: `x0,z0,ydot0,T_TU,jacobi,k1,k2,k_z,corrector_residual,independent_closure_residual,floquet_tag,seed_row,seed_route`. Use the #321 joblib parallel substrate. `_print_progress` timestamps + per-seed-family member counts + lock-vs-collapse tally + closure distribution.
- [ ] Smoke: one planar seed, z0∈{−0.05}, n_steps_max=20 — confirm it lifts, converges, continues, writes records.
- [ ] GOLDEN gate (sourced): re-running the C11a (1,1) seed at z0≈−0.05 with pseudo-arclength must yield ≥75 members all with `independent_closure_residual < 1e-9` and consistent planar `(k1,k2)=(1,1)` — reproduces the #287 spike family. Encode as a `@pytest.mark.slow` test in `tests/search/test_scan_434_spike_reproduction.py`.
- [ ] ruff + mypy; commit `search/#434: 3D broken-plane discovery sweep driver + spike-reproduction golden`.

## Task 4: 3D literature-check + validation widening

**Files:** Modify `src/cyclerfinder/search/literature_check.py` (extend `CandidateSignature` with optional `topology_3d: dict | None = None`); create `src/cyclerfinder/genome/known_corpus_3d.py` (sourced spatial-CR3BP anchors: Antoniadou-Voyatzis 2018 spatial resonant families arXiv:1811.09442; Folta-Bosanac-Cox-Howell 2017 L1/L2 halos; Howell 1984 halos — citations + topology + μ + Jacobi-band only, NO fabricated ICs); expand `src/cyclerfinder/data/validation/v1_3d.py` if a same-model V1 closure assertion is missing for 3D.

- [ ] Step 1: Failing test: `CandidateSignature(..., topology_3d={"k1":1,"k2":1,"k_z":0,"max_z_km":92567})` constructs and round-trips; `check_literature` with a 3D signature against the new corpus returns a documented status in {published, likely-rediscovery, not-found}. Verify `known_corpus_3d` anchors all have a citation + DOI/arXiv (no anchor without provenance — golden-sourced discipline). Step 2 FAIL. Step 3 implement (additive `topology_3d` field; spatial matcher = `(k1,k2,k_z)` tuple + Jacobi-band overlap). Step 4 PASS.
- [ ] Step 5: `uv run pytest tests/data tests/search -q` (per catalogue-ratchets discipline — ANY signature/corpus change ripples into frozen-census ratchets). ruff + mypy; commit `search/#434: 3D topology literature-check + spatial-CR3BP known corpus`.

## Task 5: harvest + gauntlet survivors + verdict + registry (controller)

- [ ] Launch `scripts/scan_434_3d_broken_plane_em.py` detached (full sweep, all planar seeds) + harness-tracked waiter. Expect tens of minutes to hours (O(seeds × z0-perturbations × ~100 continuation steps) 3D 6-arc corrections; ~50ms each, parallel).
- [ ] Harvest → rank members by novelty: cluster by `(k1,k2,k_z)`, run `check_literature` (Task 4) + `FalsePosFlagger().score()` on representatives. Separate `not-found` (novel-topology candidates) from `published`/`likely-rediscovery` (halos/NRHO/vertical re-derivations).
- [ ] For each `not-found` candidate family: run the available 3D V0–V5 gauntlet entry (`v1_3d`) + ML flagger; do NOT write catalogue (admission is a separate follow-on — the Umbriel template).
- [ ] Verdict `docs/superpowers/plans/2026-06-24-434-3d-broken-plane-verdict.md`: spike-reproduction confirmation, total converged members, family clusters by `(k1,k2,k_z)`, lit-fresh vs rediscovery counts, any gauntlet-passing novel family. Honest framing: even a "all clusters are rediscoveries of published spatial families" is a registry-grade negative + an infrastructure validation.
- [ ] Method-versioned `data/empty_regions.jsonl` entry (3D broken-plane Earth-Moon Aldrin-lift region; method = z0-lift + pseudo-arclength continuation). `uv run pytest tests/data tests/search -q` before commit.
- [ ] If a genuine catalogue-able novel 3D cycler survives all gates → file a follow-on admission task (do not admit inline).

## Self-review
- Reuses the delivered Phase-1 corrector + Phase-2 tracer unchanged; new code is seed-gen + topology + sweep + lit/validation. ✓
- Task 0 spike gate de-risks before the full build. ✓
- Goldens sourced or closure-self-checks; known_corpus_3d anchors all carry provenance. ✓
- Catalogue-ratchet discipline: Task 4 runs the full tests/data+tests/search suite. ✓
- Report-only; survivor admission is a gated follow-on (closure discipline). ✓
