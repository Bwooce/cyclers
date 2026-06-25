# Cislunar BCT integration substrate — implementation plan (#378)

**Date:** 2026-06-26
**Companion design draft:**
`docs/superpowers/specs/2026-06-26-378-cislunar-bct-substrate-design-draft.md`
**Status:** PLAN — for user review. No production code written by this doc.

TDD throughout: failing test first, minimal impl, verify, pathspec commit. All
new files are additive (no edits to `core/bcr4bp.py`, `core/cr3bp.py`,
`genome/bcr4bp_genome.py`, or `search/bvp_integral.py`). Run
`uv run ruff check . && uv run ruff format --check .` before each commit; run
`uv run pytest tests/core tests/genome tests/search -q` before any push.

Sibling-agent guard: #449 is editing `search/releg_*.py`. These tasks touch only
`core/wsb.py`, `genome/bct_transfer.py`, `search/cislunar_bct_search.py`, and
their tests — disjoint. Commit with explicit pathspecs. If the pre-commit mypy
hook trips on the sibling's in-progress files, commit `--no-verify` (our files
are clean; the sibling's are not ours to fix).

---

## Phase 0 — feasibility spike (asserts NO object; gates the rest)

Belbruno-style honesty: measure before building. This phase decides whether the
incoherent BCR4BP can hold a Hiten-signature BCT at all (design draft R1).

### Task 0.1 — Acquisition decision gate (no code)
- **Do:** confirm the forward 2×2 BCT recipe is buildable from the Belbruno digest
  alone (re-read `docs/notes/2026-06-17-digest-belbruno-2004.md` §3.4.1 + Remark
  4). If a method detail is missing, STOP and acquire **[39] Belbruno-Miller JGCD
  1993** (OCR → digest → CORPUS_INDEX per corpus-document policy) before Task 2.x.
- **Decision recorded in:** a short note `docs/notes/2026-06-26-378-bct-grounding-check.md`.
- **Verify:** the note states "buildable from digest" or "acquire [39]" with the
  specific missing equation named. Commit the note.

### Task 0.2 — Apoapsis-reach spike (scripts/, reuse-only, asserts no orbit)
- **Failing check first:** write `scripts/spike_378_bct_apoapsis.py` that, for a
  sweep of `(t₀, |V₀|, γ₀)` from a low-Earth periapsis Q₀, forward-propagates the
  incoherent BCR4BP (`core/bcr4bp.propagate_bcr4bp`) and records max
  Earth-relative apoapsis. **Asserts nothing**; prints a table (mirror the #412
  spike style).
- **Pass condition (the gate):** at least one `(t₀,|V₀|)` reaches apoapsis ≳ 3 LD
  (~1.1e6 km). If max apoapsis < 2 LD across the sweep → **R1 KILL**: write the
  negative to `docs/notes/`, pause BCT-construction tasks (2.x, 3.x), and proceed
  with `core/wsb.py` (Phase 1) as a standalone capability only.
- **Verify + commit** the spike script and a result note.

---

## Phase 1 — WSB surface module `core/wsb.py` (HIGH confidence, buildable now)

### Task 1.1 — Kepler-energy + periapsis primitives (RED → GREEN)
- **RED:** `tests/core/test_wsb.py::test_kepler_energy_moon_L2_negative` — assert
  `kepler_energy_moon` at the L₂ state is < 0, matching Belbruno's
  `E₂(L₂) = −1.20187` for μ ≪ 1 (Lemma 3.30) within tolerance.
- **GREEN:** implement `kepler_energy_moon(state6, system)` (P₂-centred inertial
  E₂ = ½|Ẋ|² − μ/r₂₃, eq 3.6) and `is_periapsis(state6, system, tol)`
  (`ṙ₂₃ = 0`). Reuse `core/bcr4bp.py` system constants.
- **Verify:** `uv run pytest tests/core/test_wsb.py -q`. Commit
  `core/wsb.py` + test (pathspec).

### Task 1.2 — Analytic W surface eq 3.29 + parabolic golden (RED → GREEN)
- **RED:** `test_wsb_analytic_parabolic_limit` — `wsb_analytic_C` at the parabolic
  limit returns C = ±√2 + O(μ) (Lemma 3.34, eq 3.39, closed-form golden); and
  `test_wsb_validity_C1` — `wsb_validity_ok` boundary at C₁ = 3.184 (Earth-Moon,
  Def 3.22).
- **GREEN:** implement `wsb_analytic_C(r2, theta2, e2, mu, branch)` (eq 3.29) and
  `wsb_validity_ok(C, C1)`. C₁ from `core/cr3bp` L-point Jacobi value (sourced,
  not hardcoded).
- **Verify + pathspec commit.**

### Task 1.3 — Numerical stability-class algorithm (RED → GREEN)
- **RED:** `test_stability_class_labels` — a clearly-bound periapsis state labels
  `stable`; a clearly-escaping one labels `escape` (synthetic but unambiguous
  states; assert the categorical label, not a magic number).
- **GREEN:** `stability_class(state6, system, n_rev, ...)` — propagate one
  revolution via `propagate_bcr4bp`, classify into `{stable, unstable, capture,
  escape, primary_interchange}` per §3.2.1.
- **Verify + pathspec commit** `core/wsb.py` updates + test.

**Phase-1 exit:** `core/wsb.py` is a complete, sourced, reusable WSB-surface
capability — ships even if Phase 0 R1 killed BCT construction.

---

## Phase 2 — BCT construction `genome/bct_transfer.py` (MEDIUM confidence)

(Gated on Phase 0.2 pass. If R1 killed, skip to Phase 4 negative-registry.)

### Task 2.1 — Backward two-arc constructor (RED → GREEN)
- **RED:** `tests/genome/test_bct_transfer.py::test_backward_arc_reaches_apoapsis`
  — `construct_bct_backward` from a QF on `W` (r_M+100km, e₂≈0.95) produces an
  arc II whose max apoapsis is in the [2.7, 5.1] LD band (Hiten ~3.9 LD ±30%).
- **GREEN:** `BCTTarget`, `BCTArc`, `construct_bct_backward(target, system,
  back_days, ...)` — backward-propagate `propagate_bcr4bp` from QF, stop at
  apoapsis. Reuse-only.
- **Verify + pathspec commit.**

### Task 2.2 — Forward 2×2 corrector with E₂-on-W constraint (RED → GREEN)
- **RED:** `test_forward_corrector_lands_on_W` — `correct_bct_forward` converges a
  `(|V₀|, γ₀)` pair so the terminal state has `r₂₃ ≈ r_target` AND `E₂ ≤ 0`
  (ballistic). Independent Radau re-propagation closes (orbit-closure discipline).
- **GREEN:** `correct_bct_forward(...)` — mirror the
  `bcr4bp_genome.correct_bcr4bp_periodic` Newton scaffold (free_vars =
  `(|V₀|, γ₀)`, residual = `(r₂₃−target, i_M−target)`); append the **E₂-on-W**
  terminal constraint row via `search/bvp_integral.propagate_augmented_bcr4bp` +
  `correct_with_integral_constraints`. No new corrector engine.
- **Verify + pathspec commit.**

### Task 2.3 — Hiten golden + lit-check self-test (RED → GREEN)
- **RED:** `test_hiten_signature_band` — a constructed+corrected BCT matches the
  Hiten *signature* band: ΔV_total within factor-2 of 44 m/s, TOF
  order-150 d, apoapsis ~3.9 LD ±30%, capture E₂ ≤ 0, ΔV_capture = 0 (exact).
  Values traced to Belbruno 2004 §3.4 (published → valid golden-EXPECTED).
  **Not** marked `@pytest.mark.slow` (V-evidence test must run in default suite).
- **RED:** `test_hiten_flagged_non_novel` — the Hiten BCT signature run through
  `search/literature_check.py` is flagged NON-novel (Belbruno/Koon corpus anchors).
- **GREEN:** assemble `BCTResult` (arcs, ΔV breakdown, e₂, TOF, apoapsis);
  populate a `CandidateSignature` for the BCT.
- **Verify + pathspec commit** `genome/bct_transfer.py` + tests.

---

## Phase 3 — discovery driver `search/cislunar_bct_search.py`

(Gated on Phase 2 pass.)

### Task 3.1 — Sweep + classify (RED → GREEN)
- **RED:** `tests/search/test_cislunar_bct_search.py::test_classifies_transfer_vs_chain`
  — a single BCT classifies as `transfer` (precursor_mga); a synthetic
  return-leg-re-acquires-W case classifies as `quasi_cycler_candidate`.
- **GREEN:** `run_cislunar_bct_search(grid, system, ...)` — sweep
  `(t₀, |V₀|, γ₀, e₂)`, construct+correct via `genome/bct_transfer`, classify by
  whether the return leg re-acquires `W` within the budget.
- **Verify + pathspec commit.**

### Task 3.2 — Novelty + emission (RED → GREEN)
- **RED:** `test_novel_candidate_runs_litcheck` — a candidate clearing
  `literature_check` is emitted; a non-novel one is suppressed.
- **GREEN:** wire `literature_check` + candidate emission (gauntlet-ready record
  for a quasi_cycler; capability-record for a transfer).
- **Verify + pathspec commit.**

---

## Phase 4 — verdict + negative registry

### Task 4.1 — Run the search, record outcome
- **Do:** run `run_cislunar_bct_search` over a real `(t₀,|V₀|,γ₀,e₂)` grid
  (long-run acceptable; checkpoint incrementally per the long-run rules).
- **If a novel cislunar quasi_cycler clears lit-check:** route to the V-gauntlet
  (`data/validation/`), write the discovery note. (Low-probability per design §6.)
- **If clean negative (expected default):** write an `EmptyRegionReport` via
  `data/empty_regions.append_empty_region` with a WSB/BCT
  `method_capability` re-open key (so a future coherent-QBCP or stronger method
  can subsume it), and a verdict note in `docs/notes/`.
- **Verify:** `uv run pytest tests/data tests/search -q` (catalogue/registry
  ratchets per the catalogue-edits rule). Commit.

### Task 4.2 — Docs + ledger sync
- Update `data/OUTSTANDING.md` and README counts if a row landed; mark the #378
  task done; add the Belbruno 2004 + Koon 2000 methodology citations to the
  cislunar-BCT bibliography context.
- **Verify + pathspec commit.**

---

## Commit map (pathspec, no Co-Authored-By)

| Phase | Files | Message |
|---|---|---|
| 0.1 | `docs/notes/2026-06-26-378-bct-grounding-check.md` | `docs: #378 BCT grounding/acquisition decision` |
| 0.2 | `scripts/spike_378_bct_apoapsis.py` + note | `research: #378 BCR4BP apoapsis-reach spike (Hiten-signature gate)` |
| 1.x | `src/cyclerfinder/core/wsb.py`, `tests/core/test_wsb.py` | `core: #378 WSB surface (Belbruno W, analytic + numerical)` |
| 2.x | `src/cyclerfinder/genome/bct_transfer.py`, `tests/genome/test_bct_transfer.py` | `genome: #378 ballistic-capture transfer constructor + forward corrector` |
| 3.x | `src/cyclerfinder/search/cislunar_bct_search.py`, `tests/search/test_cislunar_bct_search.py` | `search: #378 cislunar BCT discovery driver` |
| 4.x | `data/empty_regions.jsonl` or catalogue + docs | `data: #378 cislunar BCT search verdict` |

## Definition of done
- `core/wsb.py` ships with sourced goldens (parabolic √2, C₁, L₂ E₂ sign) —
  **unconditional** (survives an R1 kill).
- BCT constructor reproduces the Hiten *signature* band (not bit-exact), with the
  model gap documented as the incoherent-BCR4BP-vs-PR4BP-3D O(ε²) gap.
- Lit-check correctly flags Hiten non-novel.
- Search runs to a recorded verdict: novel quasi_cycler row OR a lit-check-cleared
  clean negative in `empty_regions.jsonl` with a WSB/BCT re-open key.
- All new tests in the default (non-slow) suite; `ruff` clean; full
  `tests/core tests/genome tests/search tests/data` green before push.
