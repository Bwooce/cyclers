# M-ED · Real-ephemeris multi-arc cycler discovery — DESIGN DRAFT

**Status:** design draft (brainstorming, 2026-06-05). No code, no task-level plan
yet. This is the design pass the roadmap requires before M-ED can be planned
(`docs/superpowers/plans/2026-06-02-multirev-3d-vem-ephemeris-roadmap.md:104-118`,
"Ready for task-level planning: NO — needs design first").

**Goal:** make the catalogue's **203 multi-arc rows**
(`docs/notes/multi-arc-classification.md:371`: 184 ocampo + 14 ch4 + mcconaghy-2006
+ 2 Sanchez Net + 2 Jones VEM) solvable on the real DE440 ephemeris — i.e. find
the closed N-arc chain (resonant E-E intervals patched by gravity-assist flybys)
that each multi-arc descriptor names, anchored to its sourced V∞/ToF multiset.

**One-line thesis:** the missing capability is **not** an optimiser — it is a
**ballistic N-arc differential corrector** with a **descriptor-driven seeding
ladder**. `scripts/correct_s1l1_twoarc.py` already proves the corrector closes
real DE440 E-M-E-E chains; this design generalises it and wires it into the
existing `optimise_cell_ephemeris` / `discover` consumption surface.

---

## 0. What the roadmap got wrong (verified read-only, 2026-06-05)

The roadmap M-ED section was written 2026-06-02 and is **stale** in three places
that materially change the design. Cite these in the plan:

1. **`optimise_cell_ephemeris` is NOT a `NotImplementedError` stub.** It is fully
   implemented (`src/cyclerfinder/search/optimize.py:1309-1516`). It resolves a
   launch epoch by V∞ phase-matching (`_resolve_t0_multi_seed`,
   `optimize.py:1207-1306`), then calls the general
   `optimise_maintenance_dv` (`search/maintain.py:425-442`). The roadmap's
   "implement the stub" framing (`roadmap:109`) is obsolete; the real work is a
   **second solver mode**, not filling a stub.

2. **The shipped objective is maintenance-ΔV (turn-deficit), not ballistic
   closure.** `optimise_maintenance_dv`'s objective is
   `_maintenance_dv_chain` — summed flyby turn-deficit ΔV
   (`maintain.py:353-368`, called from `_objective` at `maintain.py:371-391`). The
   working S1L1 prototype instead drives **V∞-magnitude continuity residuals to
   zero** with `least_squares` (`correct_s1l1_twoarc.py:96-106`,
   `_residuals` → three `|V∞_in| − |V∞_out|` terms). These are **different
   problems**: turn-deficit ΔV is ≈0 wherever a flyby can ballistically bend the
   turn (a wide feasible plateau — `maintain.py:564-566` even documents the "flat
   ΔV plateau"), so minimising it does **not** pin the ballistic family. The
   prototype's residual *does*. This is the central architectural finding.

3. **`data/discover.py:99-106` still believes the stub raises.** The V3 branch of
   `_auto_validate` calls `optimise_cell_ephemeris(...)` expecting
   `NotImplementedError` ("raises until M6b lands", `discover.py:105`). It no
   longer raises — so `enable_v3=True` today runs a real (maintenance-ΔV) solve
   and unconditionally sets `level="V3"`. M-ED must replace this branch with a
   real ballistic-closure gate, not just "flip a flag" (`roadmap:112`).

Net: the infrastructure (epoch resolution, multi-rev Lambert threading,
body-agnostic chain build) is **already built and general**. M-ED's genuinely new
content is (a) a ballistic-closure solver mode, (b) the descriptor→Cell seeding
ladder, and (c) honest validation gates.

---

## 1. The prototype, exactly (what to generalise)

`scripts/correct_s1l1_twoarc.py` is the load-bearing prior art. Anatomy:

- **Free variables** `x = [t0_offset_days, T_EM, T_ME]`
  (`correct_s1l1_twoarc.py:125-133`). The phasing leg `T_EE` is **pinned** by the
  sourced 2-synodic period (`T_EE = PERIOD_DAYS − T_EM − T_ME`,
  `:78`; `PERIOD_DAYS = (1.4612 + 2.8096) * 365.25`, `:40` — the two Russell
  4.991gG arc ToFs). So the period is a *constraint*, not a free variable.
- **Residuals** (square 3×3 system, `:96-106`): three flyby V∞-magnitude
  continuity terms — `|V∞_M_in| − |V∞_M_out|` (Mars conserves energy),
  `|V∞_E1_in| − |V∞_E1_out|` (intermediate Earth flyby),
  `|V∞_E2_in| − |V∞_E0|` (cycle closure at Earth). A flyby only *rotates* V∞, so
  magnitude continuity is the ballistic-closure condition; the required turn angle
  is checked **post-hoc** for feasibility (`_max_bend_deg`/`_bend_deg`,
  `:109-120`, `feas` at `:151`).
- **Branch-pinned Lambert legs** (`:65-67`): E→M `max_revs=0`, M→E `max_revs=1`,
  E→E `max_revs=2`; the specific `(n_revs, branch)` are picked by `_pick`
  (`:48-52`) — i.e. the descriptor topology is hard-wired per leg.
- **Multi-start = epoch scan over branches** (`:162-181`): the OUTER loop is over
  EE branch `(2,"low")`/`(1,"low")` and launch-epoch offsets `range(-180, 1100,
  80)`; the corrector runs `least_squares` (method `"lm"`) from each. Closed iff
  `max_res < 0.1` km/s.
- **Result:** ~20 closed ballistic E-M-E-E cyclers on DE440; sourced 5.65/3.05
  is the **cross-check, never fitted** (`:39`, `:159`). Documented limit: the
  closed family **floors Mars V∞ ≈ 6.4 km/s** (project memory
  `project_s1l1_realeph_closure_blocker.md`; OUTSTANDING `data/OUTSTANDING.md:45`)
  — it closes a *different* family than the (unverified) 5.65/3.05 anchor.

The preceding chain — `construct_s1l1_twoarc_realeph.py` (free epoch scan, no
root-find, `:50-97`), `close_s1l1_realeph.py` (drove the shipped
`optimise_cell_ephemeris` and **collapsed into a degenerate basin**, `:1-15`),
`diagnose_s1l1_twoarc.py` (coplanar θ-scan showing the single-arc model floors
~4.9/4.98, never 5.65/3.05) — is the evidence trail for *why* the corrector beats
the optimiser: the optimiser's flat ΔV plateau lets it drift onto a degenerate
high-V∞ family; the corrector's residual pins the ballistic node conditions.

---

## 2. Architecture

### 2.1 The N-arc model

Generalise the 3-variable prototype to an **arbitrary closed encounter sequence**
`B0-B1-…-Bn` with `B0 == Bn` (e.g. `E-M-E-E` for S1L1, `E-M-E-V-V-E` for the
Jones EMEVVE outbound, `M-E-E-V-E-M` for MEEVEM inbound). The chain is exactly
the `Cell.sequence` already enumerated/loaded (`search/sequence.py:107-111`;
`tests/_catalogue_loader.py` admits these as `CONSTRUCTIBLE_MULTIBODY` /
`MULTI_ENCOUNTER_SEQUENCE`).

**Free variables (recommended parameterisation):** `x = [t0, T_1, …, T_{n-1}]` —
launch epoch (s since J2000) + per-leg ToFs (days), the *same* shape
`optimise_maintenance_dv` / `_build_chain` already consume
(`maintain.py:303-350`). This reuses the entire forward map. The **period
constraint** (sum of leg ToFs == sourced `T = _target_period_sec(cell)`,
`optimize.py:236-275`) is imposed either by:
- (a) **elimination** — pin one leg ToF to `T − Σ(others)`, as the prototype
  pins `T_EE` (`correct_s1l1_twoarc.py:78`); dimension `n−1`. Simplest; recommended.
- (b) **constraint residual** — add `Σ T_j − T` as one more residual row.
  Needed only if no single leg is a clean "slack" leg.

**Residuals (the ballistic-closure system):**
- For every **intermediate** encounter `Bi` (1 ≤ i ≤ n−1): `|V∞_in(Bi)| −
  |V∞_out(Bi)|` (the flyby conserves V∞ magnitude).
- **Closure** at `B0`: `|V∞_in(Bn)| − |V∞_out(B0)|` (periodicity).
- This is the prototype's residual set (`correct_s1l1_twoarc.py:102-106`)
  generalised: `n−1` intermediate + 1 closure = `n` residuals against `n−1` free
  vars (after period elimination) → over-determined by one (the prototype's 3×3 is
  the `n=3` special case where they balance because t0 adds the extra DOF). Use
  `least_squares` (Levenberg–Marquardt / trust-region), not a square solve, so the
  system size is irrelevant.

**Bend feasibility (post-hoc, never in the residual):** at each flyby check
`required_turn ≤ max_turn(V∞, r_p_safe)` (`correct_s1l1_twoarc.py:109-120`;
the same `max_bend` / turn-deficit machinery lives in `maintain.py:280-295`). A
ballistically-closed chain whose flyby exceeds `max_turn` is a *powered* cycler,
not ballistic — surface it, do not fit toward it.

### 2.2 Why a corrector, not the existing optimiser

The shipped `optimise_cell_ephemeris` is the **right consumption surface** but the
**wrong objective** for ballistic discovery (§0 finding 2). The design keeps both:

- **maintenance-ΔV mode** (existing) — answers "what is the per-cycle station-
  keeping ΔV of this geometry?" (the M7 TCM-budget question). Correct as-is.
- **ballistic-closure mode** (new) — answers "is there a *ballistic* member of
  this family on the real ephemeris, and where?" This is the M-ED discovery
  question for the 203 multi-arc rows, ~all of which are ballistic/near-ballistic
  (Russell free-return, Sanchez Net ≤10 m/s, Jones <200 m/s).

The corrector becomes a **new mode of `optimise_cell_ephemeris`** (a `mode=`
keyword, `"maintenance"` default vs `"ballistic"`), reusing `_target_period_sec`,
`_ephemeris_tof_seed_and_bounds`, `_resolve_t0_multi_seed`, and the
`OptimisationResult` shape. The ballistic mode runs the §2.1 corrector instead of
`optimise_maintenance_dv`, then reports the same result fields (`converged` =
`max_res < tol`; `constraints_satisfied` = converged ∧ V∞-cap ∧ bend-feasible;
`closure_residual_kms` = `max_res`, which here is a *real* residual unlike the
maintenance-ΔV proxy at `optimize.py:1504-1506`).

### 2.3 Multi-start strategy — the degenerate-basin lesson

The single most important lesson from the prior art (project memory
`project_s1l1_realeph_closure_blocker.md`; `close_s1l1_realeph.py:1-15`):
**`optimise_cell_ephemeris` collapsed into a degenerate basin and never reached
the family the direct construction found.** Seeding matters more than the
optimiser. Therefore:

- Multi-start is an **epoch × branch scan** (prototype `:162-181`), not random
  restarts on a flat objective.
- The **ToF seed is the descriptor's sourced arc ToFs**, asymmetric, never
  equispaced — this is exactly why `_resolve_t0_multi_seed` exists and takes a
  `seed_days` override (`optimize.py:1414-1425` widens bounds to `[0.05,0.95]·T`
  for asymmetric seeds; the docstring at `:1408-1413` names the S1L1 "154 d
  outbound + long return" degenerate-basin failure as the motivation).
- Per-leg `(n_revs, branch)` come from the **descriptor**, not enumeration (§3).

---

## 3. Seeding ladder

The corrector is only as good as its seed. The ladder, cheapest/most-sourced first:

1. **Descriptor genome → topology + TOF seed (primary).** Spec §16.7.7
   `free_return_arcs[]` (`docs/spec.md:975-1016`): each arc carries `arc_type`
   (generic/half-rev/full-rev → maps to leg `n_revs`/`branch`), `tof_years` (the
   E-E leg ToF for g/h arcs), and `resonance` (M:N for full-rev arcs, which *sets*
   the ToF). 12 catalogue rows carry explicit descriptors (`docs/spec.md:1018`).
   The parser maps a descriptor string (e.g. `g(1.4612,…) G(2.8096,…)` →
   two generic arcs, ToFs 1.4612 yr + 2.8096 yr) onto:
   - per-leg `n_revs`/`branch` for the `Cell` (`Cell.per_leg_revs/per_leg_branch`,
     `sequence.py:109-110`), and
   - the asymmetric `tof_seed_days` for `_resolve_t0_multi_seed`.
   This is the "Descriptor as genome" design principle already recorded in spec
   (`docs/spec.md:1074-1079`).

2. **Catalogue anchors (sourced V∞/ToF multisets).** The epoch resolver needs
   per-body V∞ targets (`vinf_targets_kms`, `optimize.py:1355-1358`). Sources:
   - Jones VEM multisets — full per-encounter V∞ for EMEVVE/MEEVEM
     (`docs/notes/2026-06-05-jones-aas17-577-vem-mining.md:345-367`, Tables 2/3;
     transit ToFs 309/259 d and 268/223 d) — these are **real M-ED targets**
     (`tests/_catalogue_loader.py:69-71`).
   - Sanchez Net Cycler 1/2 event V∞ + ToFs (`s1l1-target-topology-mining.md:209-219`;
     ME transit ~136 d / ~1026 d) — already cross-validated by the prototype's
     regime (DE440, 2030–2034 window, near-ballistic).
   - Russell Table 4.9–4.13 arc ToFs for the 14 ch4 rows.
   - **NOT** the S1L1 5.65/3.05 pair — it is **unverified-provenance** (untraceable
     in Patel/McConaghy, `s1l1-target-topology-mining.md:88,240`; flagged
     `unverified` in commit `b9dd64b`). Use it only as a non-gating diagnostic.

3. **Coplanar / idealized warm starts.** `optimise_cell_idealized` (circular
   model, `optimize.py:955`) + the catalogue warm-start mechanism
   (`docs/superpowers/specs/2026-06-01-catalogue-seeded-warm-start-design.md`)
   give a cheap 2D guess for the interior epochs; `diagnose_s1l1_twoarc.py` shows
   the coplanar θ-scan brackets the transfer-angle family even when it cannot
   reach the real V∞. Use the coplanar solution's leg ToFs as the asymmetric seed
   when no descriptor exists.

4. **Epoch scans.** Last resort / completeness: scan `t0` across ±N synodic
   windows around the row's `priority_date` (prototype `:164`;
   `_resolve_t0_multi_seed` already centres a ±10 yr window, `optimize.py:1293`).

The ladder degrades gracefully: descriptor → anchor → coplanar → blind scan, each
rung used only when the rung above is absent for that row.

---

## 4. Where it lives

- **`search/correct.py` (new module): the N-arc ballistic corrector.** Pure,
  body/length-agnostic, depends only on `core/lambert`, `core/ephemeris`,
  `core/constants` (the prototype's imports). Public:
  `ballistic_correct(sequence, per_leg_revs, per_leg_branch, t0_seed, tof_seed,
  period_sec, ephem, *, vinf_cap, rp_factors) -> BallisticClosureResult`. This is
  the prototype's `_solve` (`correct_s1l1_twoarc.py:123-152`) lifted to N arcs.
- **`search/descriptor.py` (new module): the genome parser.** Russell descriptor
  string → `(per_leg_revs, per_leg_branch, tof_seed_days)`. Reads
  `free_return_arcs[]` (spec §16.7.7). Isolated so the corrector stays catalogue-
  agnostic.
- **`optimise_cell_ephemeris` (extend, do not fork):** add `mode="maintenance" |
  "ballistic"` (`optimize.py:1309`). `"ballistic"` reuses epoch resolution +
  seeds, calls `search/correct.py`, maps to `OptimisationResult`. Default unchanged
  → every existing caller/test byte-identical.
- **`generalise `solve_powered_periodic_cycler` off the Aldrin lock`:** `bvp.py:92`
  pins E-M-E via `optimise_aldrin_maintenance_dv` (`bvp.py:61,149`). M-ED does
  **not** need this for ballistic discovery (it goes through `optimise_cell_ephemeris`),
  but the roadmap lists it (`roadmap:110`); recommend deferring it to M7 (powered
  cyclers) as an explicit non-goal here.
- **`data/discover.py`:** (a) fix the stale V3 branch (`discover.py:99-106`) to
  call the ballistic mode and gate on real closure, not assume an exception; (b)
  `discover(..., optimiser="ephemeris")` already routes to
  `optimise_cell_ephemeris` (`discover.py:179-193`) — add the `mode` passthrough.
  **What `discover()` gains:** for the first time it can record a *ballistic
  multi-arc closure* on real DE440 against the `analytic-ephemeris` signature pool
  (`discover.py:218`) and surface a sourced-anchor match for the 203 rows.
- **`find_cyclers`:** unchanged in M-ED. It is circular-only by contract
  (`optimize.py:1592-1598`); ephemeris discovery is the `discover()` path. Note
  but do not touch.

---

## 5. Validation discipline (golden rules)

Golden discipline (project memory `feedback_golden_tests_sourced_only.md`): the
EXPECTED side of every gate traces to a published source; never a value our own
code computed. Applied to M-ED:

- **Gating anchors (sourced, flip xfails green):**
  - **Jones VEM multisets** — `tests/test_vem_rediscovery.py:220-258`
    (`test_emeeve_idealized_optimiser_converges_feasible`, currently `xfail`,
    "Flipped by M-ED"). Ballistic-mode convergence of the EMEVVE/MEEVEM cells to
    the sourced per-encounter V∞ (within a stated tol) is the cleanest M-ED gate:
    the multisets are fully source-traced (`jones-...-mining.md:393-416`).
  - **Sanchez Net regime cross-validation** — already done by the prototype
    (DE440, 2030–2034, near-ballistic E-M-E-E); promote to a test asserting the
    corrector closes ≥1 chain in that regime with V∞ ≤ cap.
  - **Russell ch4 arc ToFs** — period round-trip already gated
    (`test_vem_rediscovery.py:132-168`); M-ED adds geometry closure for the rows
    whose descriptors are complete.
- **Stays xfail / NOT a gate:**
  - **S1L1 5.65/3.05** — unverified provenance (§3.2); the *solver* need not wait
    (it closes the family at Mars-V∞≈6.4), but no green test may assert 5.65/3.05.
    Keep `EXPECTED_SKIPS["s1l1-2syn-em-cpom"]` (`real_closure.py:222-236`) and the
    M5 `test_2syn_em_rediscovers_5_65_kms_earth` xfail until/unless McConaghy 2006
    (4.7/5.0, `data/OUTSTANDING.md:294-295`) is ingested as a traceable anchor.
    The roadmap's "S1L1 5.65/3.05 rediscovered" success criterion (`roadmap:116`)
    is **superseded** — re-target it to the Jones VEM gate.
  - **`mcconaghy-2006-em-k2`** — incomplete leg data (`real_closure.py:237-243`);
    same physical cycler as S1L1 (`multi-arc-classification.md:621-634`).
- **Eventual McConaghy 4.7/5.0** — when ingested with full provenance, becomes the
  S1L1 gate (the corrector must then find a *4.7/5.0* family member, a different
  basin from the 6.4 it floors at today — open research, §7).
- **Census ratchet:** flipping `MULTI_ENCOUNTER_SEQUENCE` /
  `CONSTRUCTIBLE_MULTIBODY` rows to constructible must keep the coverage-audit
  census invariant intact (`roadmap:93`, `tests/_catalogue_loader.py`).

---

## 6. Dependencies + scope cuts

### 6.1 What 2D-on-real-ephemeris already gives (the key clarification)

**The corrector is coplanar-*position*-based on real ephemerides — and that is
already 3D-sufficient for ballistic V∞ closure.** Precisely: the prototype reads
full 3D DE440 states (`ephem.state` returns 3D `r,v`,
`core/ephemeris.py:129`, `:144-151`) and Lambert is **3D-native** (takes 3D
`r1,r2`, `roadmap:76`). The V∞ vectors `v_sc − v_planet`
(`correct_s1l1_twoarc.py:83-84`) are genuine 3D vectors; their *magnitudes*
(the residuals) and the bend angles (`_bend_deg`, full 3D dot product,
`:117-120`) are correct in 3D **without any frame change**. The "coplanar"
limitation the roadmap worries about (`roadmap:69-84`) lives only in the
**rotating-frame periodicity closure** (`frames.py` `synodic_omega` assumes the
ecliptic is the orbit plane) and in the **idealized circular model**. The M-ED
ballistic corrector closes on **V∞-magnitude continuity at nodes**, which needs no
rotating frame at all. So:

- **M-ED does NOT strictly need M-3D.** The ballistic node-closure residual is
  frame-free and already 3D-correct on DE440 states. This is a scope cut the
  roadmap's critical path (`M-L → (M-3D ∥ M-N) → M-ED`, `roadmap:7`) overstates.
- **M-ED DOES need M-L (multi-rev Lambert).** The prototype calls
  `lambert(..., max_revs=1/2)` (`correct_s1l1_twoarc.py:66-67`); the resonant E-E
  intervals are multi-rev by construction. **Verify M-L's status before
  planning** — the prototype runs today, so either M-L landed or the single-rev
  path suffices for the cases tried; the plan must confirm
  (`core/lambert.py` max_revs handling).
- **M-N gives the loader admission** (`CONSTRUCTIBLE_MULTIBODY` already ships,
  §0; `MULTI_ENCOUNTER_SEQUENCE` E-E-M-M rows still excluded,
  `tests/_catalogue_loader.py:240-241`). M-ED needs the multi-encounter rows
  surfaced as `Cell`s with correct per-leg revs/branch — that is the
  descriptor parser (§4) plus the loader generalisation (M-N). Treat M-N's
  structural-inference piece as a co-dependency, buildable inside M-ED via the
  descriptor parser.

### 6.2 Explicit non-goals

- **Rotating-frame multi-lap closure / drift** (M-3D's `REAL_DRIFT_TOLERANCE_KM`,
  `roadmap:77`). M-ED closes node-to-node ballistically; multi-lap drift
  verification is a *separate* validation axis (the existing V2 gate,
  `discover.py:91-97`), left as-is.
- **Powered cyclers** (the Aldrin BVP, `bvp.py`; turn-deficit > max). M-ED reports
  bend-infeasible chains honestly but does not solve the powered TCM problem.
- **M7 per-family TCM-budget machinery** (`roadmap:111`). The maintenance-ΔV mode
  already computes per-cycle ΔV (`maintain.py:_maintenance_dv_chain`); a *budget*
  (per-family ceiling) is a reporting/threshold layer, deferrable. M-ED's gate is
  ballistic closure, not a ΔV budget.
- **Inclination-aware idealized pre-filter** — keep the circular model 2D
  (`roadmap:78`); it is only a warm-start source (§3.3).
- **Blind (anchorless) discovery.** Epoch resolution requires V∞ targets
  (`optimize.py:1356-1358`); M-ED rediscovers *sourced* families. True blind
  search is the Forge Phase 4 (`the-forge-pipeline.md:29-31`), which can consume
  the M-ED corrector + the VILM ΔV-min pruning bound
  (`the-forge-pipeline.md:84-88`) as its feasibility screen.

---

## 7. Architecture options + recommendation

### Option A — Corrector-as-polisher behind the existing optimiser
Run `optimise_cell_ephemeris` (maintenance-ΔV) to get a geometry, then polish with
the ballistic corrector.
- **Pro:** minimal new surface; reuses everything.
- **Con:** the optimiser lands the *degenerate basin* first
  (`close_s1l1_realeph.py:1-15`); polishing from there inherits the wrong family.
  This is the documented failure mode. **Rejected.**

### Option B — Corrector-first construction (recommended)
Seed from the descriptor genome (§3.1), build the chain, run the ballistic
corrector directly (the prototype's structure), report through
`OptimisationResult`. The maintenance-ΔV optimiser becomes a *post-hoc* cost
annotator on the closed geometry, not the family selector.
- **Pro:** matches the only thing proven to work; seeding (not optimisation) drives
  family selection, honouring the degenerate-basin lesson; frame-free 3D-correct;
  smallest physics risk.
- **Con:** needs the descriptor parser + a new corrector module; for rows with no
  descriptor it falls back to coplanar/epoch-scan seeds (graceful, §3).
- **Recommendation: adopt B.**

### Option C — Hybrid: corrector for selection, optimiser for refinement
Corrector closes the ballistic node conditions (family selection); then a bounded
maintenance-ΔV / multi-lap pass refines within the closed basin and produces the
TCM cost.
- **Pro:** B's robustness + a single artefact carrying both closure residual and
  ΔV cost; natural fit for the Forge gauntlet's multi-axis verdict
  (`the-forge-pipeline.md:21-27`).
- **Con:** more moving parts; the refinement must be basin-locked (bounded around
  the corrector solution) or it re-introduces the drift.
- **Recommendation: build B first; layer C's refinement only if the TCM cost is
  needed by a consumer (Forge / M7).** B is the MVP; C is the polish.

### Phased build sketch (not a task plan)
1. **Corrector core** (`search/correct.py`) — lift `_solve` to N arcs; unit-test
   against the prototype's S1L1 closure (the ~20 closed chains are a regression
   fixture; their V∞ are OUR computation → not golden, used only for
   non-regression of the *solver*, not as published anchors).
2. **Descriptor parser** (`search/descriptor.py`) — `free_return_arcs[]` →
   revs/branch/tof_seed; gate on the 12 sourced-descriptor rows.
3. **Mode wiring** — `optimise_cell_ephemeris(mode="ballistic")`; default-inert.
4. **Seeding ladder** — wire descriptor → anchor → coplanar → scan.
5. **Discover + gates** — fix `discover.py` V3 branch; flip the Jones VEM xfail;
   add the Sanchez-regime gate; ratchet the census.
6. **(Optional, Option C)** basin-locked maintenance-ΔV refinement + TCM cost.

---

## 8. Open questions for the user

1. **Objective confirmation.** Do you agree the M-ED objective is **ballistic V∞
   continuity closure** (the prototype), not the shipped maintenance-ΔV
   minimisation — i.e. the corrector is a genuinely new solver mode, not a re-tune
   of `optimise_cell_ephemeris`?
2. **S1L1 success criterion.** The roadmap's M-ED gate is "S1L1 5.65/3.05
   rediscovered" (`roadmap:116`), but 5.65/3.05 is unverified-provenance and the
   corrector floors Mars-V∞≈6.4. Do you accept **re-targeting** the headline gate
   to the **Jones VEM multisets** (sourced, EMEVVE/MEEVEM), leaving S1L1 as a
   non-gating diagnostic until McConaghy 4.7/5.0 is ingested?
3. **M-3D dependency.** Do you accept the §6.1 claim that M-ED's ballistic
   node-closure is **frame-free and already 3D-correct on DE440**, so M-ED does
   **not** block on M-3D (only the rotating-frame *multi-lap drift* validation
   does)? This removes M-3D from M-ED's critical path.
4. **M-L verification.** The prototype already calls `lambert(max_revs=1/2)` and
   runs — has M-L landed, or should the plan's first task be to confirm the
   multi-rev path is real before relying on it?
5. **Scope of the first cut.** Ship **Option B** (corrector-first, ballistic-only,
   no TCM budget) as the M-ED MVP and defer Option C's ΔV refinement + M7 TCM
   budget to a follow-up — agreed?
6. **Descriptor coverage.** Only 12 rows carry complete `free_return_arcs[]`
   descriptors and ~all `russell-ocampo-*` (184 rows) are gapped
   (`docs/spec.md:1018-1022`). For the gapped rows the seed falls back to
   coplanar/AR-TR. Is **descriptor-bearing rows first, gapped rows via fallback
   later** an acceptable staging, or must all 203 close in one pass?

---

## Approval (2026-06-05)

User-approved with all recommendations accepted: (Q1) the M-ED objective is
ballistic V∞-continuity closure — a new solver mode, not a maintenance-ΔV
re-tune; (Q2) headline gate re-targeted to the sourced Jones AAS 17-577 VEM
multisets, S1L1 demoted to non-gating diagnostic pending McConaghy 2006;
(Q3) M-ED does NOT block on M-3D (critical path updated); (Q4) M-L confirmed
landed; (Q5) Option B corrector-first ballistic-only MVP, TCM budget deferred;
(Q6) descriptor-bearing rows first, gapped rows via coplanar/AR-TR fallback
later. Next step: implementation plan.
