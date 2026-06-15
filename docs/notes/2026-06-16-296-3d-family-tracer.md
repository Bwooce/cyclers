# #296 Phase 2 — 3D family tracer + continuation driver

Date: 2026-06-16
Issue: #296 (Phase 2 of #291)
Predecessor: Phase 1 corrector at commit `6da8ee1`
              (`src/cyclerfinder/search/cr3bp_general_periodic_3d.py`)

## What Phase 2 delivers

A predictor-corrector continuation driver that consumes the Phase 1
`correct_general_periodic_3d` corrector and walks a 1-parameter family curve
through the full 6D state space:

- **Module:** `src/cyclerfinder/search/cr3bp_3d_family_tracer.py`
  - `continue_general_3d_family(...)` — driver supporting three stepping modes:
    - `pseudo_arclength` (default, robust past folds) — walks the unit null
      vector of the 3x4 Jacobian in the 4-D space `(x0, z0, ydot0, T)`,
      correcting onto `{r(z) = 0, tau . (z - z_pred) = 0}`.
    - `natural_T` — advances T as the natural parameter; corrector solves
      `(x0, z0, ydot0)`.
    - `natural_x0` — advances x0; corrector solves `(z0, ydot0, T)`
      (matches the #287 spike's pattern, known to fold around x0 ~ -0.81).
  - `Family3D` dataclass — list of `Family3DMember`s, fold-turning points,
    termination reasons per direction, walk metadata.
  - `Family3DMember` — wraps Phase 1's `Periodic3DOrbit` with monodromy,
    Floquet multipliers, and a coarse stability tag.
  - `FoldPoint` — recorded sign-flips in tangent components (pseudo-arclength
    walks through; natural-parameter walks stop cleanly at folds).

- **Tests:** `tests/search/test_cr3bp_3d_family_tracer.py`
  - 10 tests, all passing in ~80s. Topology-only goldens (per
    `feedback_golden_tests_sourced_only`): the sourced anchor is the planar
    Braik-Ross C11a member at `data/catalogue.yaml` row
    `braik-ross-c11a-cycler-2026`; the 3D walk's correctness is verified by
    (1) spike-extent reproduction, (2) closure preservation, (3) Floquet
    reciprocal-pair / unit-pair structure (Hamiltonian symmetry), (4) fold
    detection + walk-through, (5) planar-collapse guard, (6) direction
    symmetry, (7) natural-T mode sanity, (8) input validation, (9) callback.

- **Family registry:** `data/family_296_3d_em_11.jsonl` (266 lines = 1
  header + 265 member rows). Each member row carries the full state vector,
  period (TU + days), Jacobi C, the 6x6 monodromy matrix (serialized), the
  6 Floquet multipliers, and the stability tag. NOT a catalogue writeback.

- **Generator script:** `scripts/family_296_3d_em_11_generator.py` —
  reproduces the JSONL from the spike's seed with `step = 0.01`,
  `n_steps_max = 200`, `direction = "both"`.

## Reproduction verdict vs the #287 spike

| Quantity | Spike (#287) | Phase 2 tracer (#296) | Notes |
|---|---|---|---|
| Family members | 164 | **265** | Tracer reaches further; pseudo-arclength walks through the spike's z0-fold |
| `x0` extent | [-0.85, -0.77] | **[-0.96, -0.76]** | Wider on the backward side (toward the L2 region) |
| `z0` extent | [-0.24, 0] | [-0.24, ~0] | Same; the planar manifold bounds z0 from above |
| `T_TU` extent | [9.5, 18.9] | [9.04, 10.29] | Spike walked a different branch (period-multiplying / homoclinic side); Phase 2's pseudo-arclength stays on the local family curve |
| `C` extent | [2.78, 3.10] | [2.92, 3.15] | Phase 2's slab; spike covered a deeper-energy region |
| Closure max (independent Radau) | 1e-8 | **3.7e-9** | Both tight; the tracer's gate is 1e-6 |
| Closure median | — | **1.4e-10** | Round-trip noise floor |
| Folds detected | (z0 fold) | **2** — T-fold (step +11), z0-fold (step -3) | Walk continues past both |
| Stability tag | (not computed) | **all 265 hyperbolic_pair** | Family is structurally unstable; reciprocal real pair off the unit circle |

The Phase 2 tracer's pseudo-arclength walk produces a STRICTLY WIDER family
than the spike's natural-parameter walk in x0 (-0.96 vs -0.85 backward bound)
because it doesn't terminate at the z0-fold near the seed (the spike's
natural-x0 stepping is rank-deficient at the fold; pseudo-arclength's
augmented system stays well-conditioned).

The T-extent difference (Phase 2 covers [9.04, 10.29] TU only; spike covered
[9.5, 18.9] TU) reflects that the spike walked into a **different branch** of
the family — likely a period-doubling or homoclinic structure off the
hyperbolic_pair seed. Phase 3's bifurcation tracking will pick that up.

Forward termination: `degenerate_planar` (the walk lands on the planar
manifold near z0 ~ 0 at step +138; correctly flagged).
Backward termination: `corrector_failed` (the corrector loses convergence
near x0 ~ -0.96; likely an integrator floor at the start of a deeper
hyperbolic regime).

## Discipline notes

- **No catalogue writeback.** The 265-member family is in a JSONL registry,
  not in `data/catalogue.yaml`. Per #291 Phase 1's note: the family is most
  likely a rediscovery of the Antoniadou-Voyatzis 2018 spatial CR3BP
  corpus (KNOWN_CORPUS commit `568d8a4`). Phase 3 will run
  `literature_check.py` to formally classify each member.
- **Independent cross-check on every member** (per
  `feedback_orbit_closure_discipline`): Radau at `rtol = atol = 1e-12`,
  L2 closure asserted ≤ `closure_tol = 1e-6`. Median: 1.4e-10.
- **Sourced-only goldens** in tests (per `feedback_golden_tests_sourced_only`):
  the planar Braik-Ross C11a member is the sourced anchor; all assertions
  on the 3D walk are TOPOLOGICAL (extent, fold presence, closure
  preservation, Floquet structure), NEVER "expected z0 == <number our own
  code produced>".
- **No novelty claims** in commit messages or this doc. Phase 3 owns the
  novelty adjudication.

## Phase 3 next steps

1. **Literature-novelty check.** Run `literature_check.py` (or its 3D-aware
   successor) on each of the 265 family members. Compare against the
   Antoniadou-Voyatzis 2018 spatial CR3BP corpus + any other 3D corpus the
   KNOWN_CORPUS surfaces. EXPECTED: all-rediscovery.
2. **Bifurcation tracking** along the family. The hyperbolic_pair tag
   suggests period-doubling / Neimark-Sacker bifurcations may sit at the
   ends of the family slab (where T sharply changes, around step +11 and
   step -3 in particular). Run the existing
   `bifurcation_detector.detect_period_multiplying` on each member's
   Floquet spectrum, refine brackets via the existing
   `family_switch_corrector` infrastructure. CANDIDATE PAYLOAD: a
   period-doubled sub-family that the spike skipped (its natural-parameter
   walk couldn't track through the fold).
3. **V0-V5 gauntlet adaptation for 3D orbits.** The catalogue's V-pipeline
   currently assumes planar IC + (k1, k2) topology. For 3D orbits the
   V0-V3 fingerprints need extension (3D-aware Tisserand, MOID with
   plane-out-of-ecliptic correction). This is the gate for the family
   becoming catalogue-admissible.
4. **Concrete next-step IC if literature-novelty surfaces a fresh
   sub-family.** Recommend Phase 3 specifically watches for any family
   member whose Floquet spectrum is NOT `hyperbolic_pair` (e.g. a complex
   non-unit pair indicating a Neimark-Sacker quasi-periodic torus). The
   265-row registry can be filtered post-hoc: `stability_tag != "hyperbolic_pair"`
   is the fresh-sub-family fingerprint.

## Files

- `src/cyclerfinder/search/cr3bp_3d_family_tracer.py` — module
- `tests/search/test_cr3bp_3d_family_tracer.py` — tests (10 cases, 80s)
- `data/family_296_3d_em_11.jsonl` — 266-line family registry
- `scripts/family_296_3d_em_11_generator.py` — generator (re-runnable)
- `docs/notes/2026-06-16-296-3d-family-tracer.md` — this doc

## Commit

`genome: 3D family tracer + continuation (#296 Phase 2 of #291)` — pending.
