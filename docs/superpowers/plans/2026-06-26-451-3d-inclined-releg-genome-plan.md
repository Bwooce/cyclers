# 3D-inclined releg genome — IMPLEMENTATION PLAN (#451)

**Date:** 2026-06-26
**Design draft:** `docs/superpowers/specs/2026-06-26-451-3d-inclined-releg-genome-design-draft.md`
**Status:** PLAN — bite-sized TDD tasks. No code written by this doc.

Each task is self-contained: a failing test first, the minimal implementation, a
verify step, and a single pathspec commit. Reuse over rebuild — the **3D Lambert**
(`core/lambert.py`, already vector-based), the **vector flyby turn-cone**
(`core/flyby.py::max_bend`/`bend_angle`/`dv_from_turn_deficit`), the **3D PO
machinery** (#287/#291), and **#449's `Releg` seam** (`search/releg_solver.py`,
`releg_moontour.py`) already exist. The new work is: an **inclined moon
ephemeris**, a **3D V∞-vector turn-cone residual**, and **one composed backend** —
NOT a new optimiser, corrector, or Lambert.

**Dependency:** PREFER sequencing this plan AFTER #449's `Releg` seam merges (the
two share the same two call sites — `discovery_campaign.py:487` and
`v2_moontour.py:267` — and would race otherwise). If #449 is unmerged, Task 4
introduces the minimal seam against the documented contract instead. See design
§5.

Conventions: `uv run pytest <path> -q` to run tests; `uv run ruff check . &&
uv run ruff format --check .` before every commit (`feedback_run_ruff_before_commit`);
**pathspec commits only, NEVER `git add -A`** (sibling agents commit to this same
tree); no `Co-Authored-By`. Prefix multi-step bash with `date -Iseconds`. Any
`catalogue.yaml` touch runs ALL ratchets (`feedback_catalogue_edits_run_all_ratchets`)
— this plan does not touch the catalogue.

---

## Task 0 — (REQUIRED, FIRST) transcribe the inclination-aware golden

**Why first (CRITICAL HONESTY).** Unlike #449 (whose VILM golden already passes
in-repo), #451's golden does NOT yet exist. No inclined-releg code is asserted
before the sourced golden exists (`feedback_golden_tests_sourced_only` — EXPECTED
side traces only to the papers, never to a number our code computes). Both
sources are already in-corpus; this is transcription, not acquisition.

- **Read (no code):**
  - Strange-Russell-Buffington 2007 "Mapping the V∞ Globe" (AAS 07-277), in
    `KNOWN_CORPUS`; pointer `2026-06-17-349-cassini-anchor-topology-label.md`.
    Extract the crank/pump sphere relations (a ballistic flyby preserves `|V∞|`
    and rotates the V∞ vector; the per-flyby rotation is bounded; pumping changes
    `|V∞|`, cranking changes its declination/inclination). **If the cassini note
    is only a pointer (not a full relation transcription), deep-read the AAS
    07-277 PDF body in `cyclers_pdf` for the exact crank-angle relation; flag if
    a deeper acquisition is needed.**
  - Heaton-Longuski 2003 (Uranian, JSR), OCR'd at
    `cyclers_pdf/papers/heaton-longuski-…2.3981.txt`. Extract the bracketed tour
    ΔV / inclination targets (V∞ ≤7.5 km/s arrival, 0.5 km/s Tisserand
    increments, 14 deg inclination reduction via Titania flybys, 0.92 km/s Ariel
    insertion).
- **Write:** `data/golden/strange_russell_vinf_globe.yaml` — invariant relations
  (crank/pump) + Heaton-Longuski bracketed targets, each with a `source:` field
  citing the paper + page/equation. NO numbers computed by our code.
- **Verify:** `uv run python -c "import yaml; yaml.safe_load(open('data/golden/strange_russell_vinf_globe.yaml'))"` parses; every value has a `source`.
- **Commit (pathspec):** `git add data/golden/strange_russell_vinf_globe.yaml && git commit -m "data: #451 inclination-aware golden (Strange-Russell V∞-globe + Heaton-Longuski Uranian)"`

---

## Task 1 — add inclination + node to the satellite registry (regression-preserving)

**Goal:** `SatelliteData` carries `inclination_deg` + `node_deg`, sourced from JPL
SSD satellite orbital tables; default `0.0` keeps every existing row coplanar.

- **Failing test:** `tests/core/test_satellites.py::test_inclination_node_fields_present_and_sourced`
  - Assert `SatelliteData` has `inclination_deg` and `node_deg` attributes.
  - Assert Amalthea's `inclination_deg` is non-zero (the headline inclined moon)
    and matches the JPL SSD value to a stated tolerance.
  - Assert every Galilean (`Io`/`Europa`/`Ganymede`/`Callisto`) and every existing
    moon defaults sanely (near-equatorial moons may be `0.0`).
- **Minimal impl:** add the two fields to the `SatelliteData` dataclass (default
  `0.0`); set Amalthea (and any genuinely-inclined registry moon) from the SSD
  table with a sourcing comment matching the existing `sma_km` sourcing block
  (lines 93-99). All other rows unchanged (default `0.0`).
- **Verify:** `uv run pytest tests/core/test_satellites.py -q`
- **Commit (pathspec):** `git add src/cyclerfinder/core/satellites.py tests/core/test_satellites.py && git commit -m "core: #451 add inclination_deg/node_deg to satellite registry"`

---

## Task 2 — inclined moon ephemeris `_moon_state_3d` (coplanar limit = exact regression)

**Goal:** a small, pure inclined circular ephemeris whose `i=0,Ω=0` slice equals
the existing `discovery_campaign._moon_state` bit-for-bit.

- **Failing test:** `tests/search/test_moon_ephemeris_3d.py`
  - `test_coplanar_limit_matches_moon_state`: for a grid of `(θ0, n, t, a, μ)`,
    assert `_moon_state_3d(θ0, n, t, a, μ, i=0.0, raan=0.0)` equals
    `discovery_campaign._moon_state(θ0, n, t, a, μ)` to machine precision (both
    `pos` and `vel`).
  - `test_inclined_state_out_of_plane`: for `i=30°, Ω=0`, assert the returned
    `pos`/`vel` have a non-zero `z` component of the expected magnitude
    (`|z_pos| = a·sin(i)·sin(θ)` style closed form) — derive EXPECTED from the
    rotation `R_z(Ω)·R_x(i)`, not from the function.
- **Minimal impl:** `src/cyclerfinder/search/moon_ephemeris_3d.py` —
  `_moon_state_3d` builds the coplanar `(x,y,0)` state (copy the `_moon_state`
  math) then applies `R_z(raan) @ R_x(incl)` to both `pos` and `vel`.
- **Verify:** `uv run pytest tests/search/test_moon_ephemeris_3d.py -q`
- **Commit (pathspec):** `git add src/cyclerfinder/search/moon_ephemeris_3d.py tests/search/test_moon_ephemeris_3d.py && git commit -m "search: #451 inclined moon ephemeris (_moon_state_3d, coplanar-limit regression)"`

---

## Task 3 — 3D V∞-vector turn-cone residual (the flyby feasibility generalisation)

**Goal:** a helper that, given the in/out V∞ **vectors** at a flyby, returns the
turn-angle, the achievable `max_bend` cone, and the bend-deficit ΔV — the 3D
generalisation of today's magnitude-only continuity check.

- **Failing test:** `tests/search/test_flyby_turn_cone_3d.py`
  - `test_coplanar_turn_equals_magnitude_check`: for two equal-`|V∞|` vectors in
    `z=0` whose turn is inside the cone, deficit ΔV = 0 (matches today's
    behaviour).
  - `test_out_of_plane_turn_deficit`: for two equal-`|V∞|` vectors with a 40°
    out-of-plane angle exceeding `max_bend`, assert deficit ΔV =
    `dv_from_turn_deficit(|V∞|, δ_req, δ_max)` with `δ_req = angle(in, out)` —
    EXPECTED computed independently from `core.flyby` primitives, not from the
    new helper.
  - `test_unequal_magnitude_flagged`: `||V∞_in|−|V∞_out|| > gate` ⇒ infeasible
    (magnitude continuity still binding).
- **Minimal impl:** a function (in `releg_solver.py` or a small
  `search/flyby_continuity_3d.py`) computing `δ_req = arccos(û_in·û_out)`, the
  cone `max_bend(μ, rp_min, |V∞|)`, and the deficit via
  `dv_from_turn_deficit` (or `dv_powered_flyby_periapsis` for the Oberth-credited
  variant). Reuses `core.flyby` entirely.
- **Verify:** `uv run pytest tests/search/test_flyby_turn_cone_3d.py -q`
- **Commit (pathspec):** `git add src/cyclerfinder/search/flyby_continuity_3d.py tests/search/test_flyby_turn_cone_3d.py && git commit -m "search: #451 3D V∞-vector turn-cone flyby residual"`

---

## Task 4 — `InclinedReleg` backend (composes #449's `Releg` seam)

**Goal:** one new `Releg` backend that runs the (already-3D) Lambert on inclined
moon states and applies the 3D turn-cone residual; `dv_kms = 0` if inside the
cone, else the bend deficit. Composes with #449's `DsmReleg` for the
powered×inclined quadrant.

- **Pre-check:** confirm #449's `releg_solver.py` + `RelegResult` are merged. If
  NOT, define `RelegResult` + a minimal `Releg` protocol here per the design §2
  contract (do not duplicate if #449 is present).
- **Failing test:** `tests/search/test_releg_solver.py::test_inclined_releg_*`
  - `test_inclined_releg_coplanar_limit_matches_ballistic`: with `i=0` moon
    states, `InclinedReleg` reproduces `BallisticReleg`'s V∞ and `dv_kms=0`
    bit-for-bit on a fixed Jovian Io→Europa leg (regression — the `i=0` slice is
    the existing genome).
  - `test_inclined_releg_amalthea_poses_out_of_plane_leg`: a Ganymede→Amalthea
    leg (Amalthea inclined, Task 1) yields V∞ vectors with a non-zero
    out-of-plane component and a turn-cone residual the coplanar path could not
    produce.
- **Minimal impl:** `InclinedReleg` in `releg_solver.py`: place moons via
  `_moon_state_3d`, call `core.lambert.lambert` (3D), compute V∞ vectors, apply
  the Task-3 turn-cone residual, return `RelegResult`.
- **Verify:** `uv run pytest tests/search/test_releg_solver.py -q`
- **Commit (pathspec):** `git add src/cyclerfinder/search/releg_solver.py tests/search/test_releg_solver.py && git commit -m "search: #451 InclinedReleg backend (3D leg, composes #449 seam)"`

---

## Task 5 — strange-russell crank/pump invariant golden test

**Goal:** assert the `InclinedReleg` turn-cone residual satisfies the V∞-globe
crank/pump invariants from the Task-0 golden.

- **Failing test:** `tests/search/test_releg_vinf_globe_golden.py`
  - Load `data/golden/strange_russell_vinf_globe.yaml`.
  - `test_ballistic_flyby_preserves_magnitude_rotates_vector`: a flyby inside the
    cone preserves `|V∞|` and the reported turn ≤ `max_bend` (pump/crank
    invariant).
  - `test_inclination_change_per_flyby_bounded`: the per-flyby inclination change
    the residual permits is ≤ the golden's crank bound.
  - EXPECTED side reads only the golden YAML (sourced to AAS 07-277).
- **Minimal impl:** none beyond Tasks 3-4 (this is a golden assertion over
  existing behaviour); if an invariant fails, fix the residual, not the test.
- **Verify:** `uv run pytest tests/search/test_releg_vinf_globe_golden.py -q`
  (NOT `@pytest.mark.slow` — keep the V-evidence test in the default suite,
  `feedback_delegation_fresh_agent_not_fork`).
- **Commit (pathspec):** `git add tests/search/test_releg_vinf_globe_golden.py && git commit -m "test: #451 V∞-globe crank/pump invariant golden (Strange-Russell)"`

---

## Task 6 — inclined-aware driver wiring + Heaton-Longuski bracket + disjoint-contour negative

**Goal:** route the inclined ephemeris + `InclinedReleg` through #449's
`releg_moontour` driver; prove the Uranian bracket and the honest
disjoint-contour negative.

- **Failing test:** `tests/search/test_releg_moontour_inclined.py`
  - `test_uranian_tour_brackets_heaton_longuski`: a Titania/Oberon/Ariel inclined
    skeleton; per-leg V∞ ≤ 7.5 km/s and the inclination-reduction ΔV brackets the
    Task-0 Heaton-Longuski golden values.
  - `test_uranus_disjoint_contour_unbridgeable_even_inclined`: an Ariel→Umbriel
    leg with disjoint Tisserand contours reports **unbridgeable** even inclined
    (reproduces `uranus-neptune-regular-moon-endgame-vilm-2026-06-23` as an
    inclined re-test) — proves inclination does not fabricate a radial bridge.
  - `test_amalthea_positive_control`: a Jovian skeleton including Amalthea POSES
    the out-of-plane leg and scores it (capability proof — the coplanar genome
    could not).
- **Minimal impl:** extend `releg_moontour.py` to accept the inclined ephemeris
  + `InclinedReleg` and use the 3D turn-cone continuity; VILM-floor prefilter
  unchanged (radial test runs first).
- **Verify:** `uv run pytest tests/search/test_releg_moontour_inclined.py -q`
- **Commit (pathspec):** `git add src/cyclerfinder/search/releg_moontour.py tests/search/test_releg_moontour_inclined.py && git commit -m "search: #451 inclined-aware moontour driver + Uranian bracket + disjoint-contour negative"`

---

## Task 7 — v2_moontour inclined-ephemeris parity (the validation mirror)

**Goal:** the validation-side mirror (`v2_moontour.py`, which imports
`_moon_state`) uses the inclined ephemeris in lockstep with discovery, so the two
seams cannot drift.

- **Failing test:** `tests/data/test_v2_moontour_inclined.py::test_v2_inclined_ephemeris_parity`
  - Assert `v2_moontour`'s per-flyby V∞ on a fixed inclined skeleton equals the
    discovery-side `InclinedReleg` V∞ to machine precision (single source of
    truth: `_moon_state_3d`).
  - Assert the `i=0` slice reproduces today's coplanar `v2_moontour` exactly
    (regression).
- **Minimal impl:** swap `v2_moontour`'s `_moon_state` import/use for
  `_moon_state_3d` (with `i=0,Ω=0` default preserving current behaviour) gated by
  the same inclined-aware flag the driver uses.
- **Verify:** `uv run pytest tests/data tests/search -q` (run BOTH suites — the
  moontour change ripples into census/validation ratchets,
  `feedback_catalogue_edits_run_all_ratchets`).
- **Commit (pathspec):** `git add src/cyclerfinder/data/validation/v2_moontour.py tests/data/test_v2_moontour_inclined.py && git commit -m "data: #451 v2_moontour inclined-ephemeris parity with discovery seam"`

---

## Task 8 — full-suite green + ruff + verdict note

**Goal:** the capability lands green; record the verdict / kill-criterion outcome.

- **Verify:**
  - `uv run ruff check . && uv run ruff format --check .`
  - `uv run pytest tests/core tests/search tests/data -q` (full relevant suite;
    confirm the coplanar regression + inclined goldens + disjoint-contour
    negative all pass).
- **Record (no catalogue writeback in this plan):** a short verdict note
  `docs/superpowers/plans/2026-06-26-451-3d-inclined-releg-verdict.md` stating:
  capability landed (inclined legs now representable + scored); which
  kill-criterion (if any) fired — in particular the §7 "inclination representable
  but immaterial for near-equatorial Galileans" outcome if the Amalthea
  out-of-plane deficit is below the gate; and whether the registry re-stamp /
  follow-on campaign is warranted. If a region is re-tested empty, re-stamp the
  corresponding `empty_regions.jsonl` entry with the subsuming inclined method +
  version (capability-subsumption record), never silent deletion.
- **Commit (pathspec):** `git add docs/superpowers/plans/2026-06-26-451-3d-inclined-releg-verdict.md && git commit -m "docs: #451 3D-inclined releg verdict + kill-criterion outcome"`

---

## Notes on reuse (do NOT rebuild)

- **3D Lambert:** `core/lambert.py` is already vector-based (solves the transfer
  plane from `r1, r2`). Do not write a new Lambert.
- **Flyby turn cone:** `core/flyby.py::max_bend`/`bend_angle`/`dv_from_turn_deficit`/
  `dv_powered_flyby_periapsis` are already V∞-vector-ready. Do not reimplement the
  bend geometry.
- **3D PO machinery / fingerprint:** `cr3bp_general_periodic_3d.py`,
  `cr3bp_3d_family_tracer.py` (#287/#291), `genome/known_corpus_3d.py`,
  `spatial_novelty_prefilter.py` exist for the 3D V0 lit-check. Do not rebuild.
- **Releg seam:** `search/releg_solver.py` + `releg_moontour.py` (#449). Extend,
  do not duplicate.
- **The ONLY genuinely new geometry** is the inclined ephemeris (one rotation)
  and the 3D turn-cone residual (vector angle vs cone) — everything else composes
  existing, tested primitives.
