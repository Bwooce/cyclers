# Implementation plan — DA/HOTM global multi-rev enumeration lane (#450)

**Design draft:** `docs/superpowers/specs/2026-06-25-450-da-hotm-enumeration-lane-design-draft.md`
**Status:** PLAN — for review. No code written by this doc.
**Philosophy:** strict TDD (failing test → minimal impl → verify → commit), each
task self-contained, exact file paths, follows existing patterns
(`cr3bp_general_periodic.py` corrector seam; `reachable_impulsive.py` sampling
precedent). Run `uv run ruff check . && uv run ruff format --check .` before every
commit. Catalogue/registry touches run the full ratchet suite
(`uv run pytest tests/data tests/search -q`). Work on `main`, pathspec commits, no
`Co-Authored-By`.

The plan ships the **sampling backend + full Png' validation first** (proves the
lane, re-opens the dead region) and defers the DA-acceleration backend to a final,
gated task. This front-loads the capability win and the published validation
target, per the design draft §0 / §10 recommendation.

> **USER DECISION 2026-06-25 (#450 backend): pure-Python truncated Taylor-map
> fallback (design-draft §8.1 option b — NO MOSEK, no DACEyPy).** The chosen
> *deliverable* backend is `DASectionMap` implemented as a pure-Python truncated
> Taylor-map (single-rev section map expanded to a chosen order about a section
> reference, composed to `Pⁿ` by truncated polynomial composition) with a
> non-commercial fixed-point root-finder (e.g. `scipy.optimize` / polynomial
> homotopy — never MOSEK/cvxpy-commercial). Build order is UNCHANGED and
> de-risks this: build the cheap `SamplingSectionMap` FIRST (Tasks 3–6) as the
> *validation oracle* + Png' recovery proof, THEN build the pure-Python
> Taylor-map `DASectionMap` to the SAME `SectionMap` interface (the former Task 7
> "DA acceleration, deferred" is now PROMOTED to an in-scope deliverable) and
> validate it against BOTH the sourced Png' golden AND the sampling oracle
> (assert the two backends agree on the recovered fixed points to a tolerance).
> Task 0's decision record should reflect this choice rather than "ship sampling,
> defer DA". Do NOT vendor DACEyPy or require MOSEK at any point.

---

## Task 0 — DA backend decision record (acquisition gate, no code)

**Goal:** resolve design-draft §8.1 / §10.1 before any DA work; the validation
path does not depend on this.

1. Probe the env once (record results): `daceypy`, `dace`, `pyaudi`, MOSEK,
   `cvxpy` availability (all absent as of 2026-06-25).
2. Write `docs/notes/2026-06-25-450-da-backend-decision.md` recording the three
   options (vendor DACEyPy+MOSEK / pure-Python Taylor fallback / ship-on-sampling)
   with the recommendation: **sampling backend first**, DA acceleration deferred to
   Task 7 behind a measured cost.
3. **Commit:** `docs: #450 DA backend decision record (ship sampling first)`.

*No production code. This is the honest acquisition flag made concrete.*

---

## Task 1 — Golden file: sourced Png' family ICs

**Test first:** `tests/data/test_png_golden.py` — load `data/golden/png_hybrid_family.yaml`,
assert it contains the 6 members (P5g', P7g'-I, P7g'-II, P7g', P9g', P3g') each
with `x0, xdot0, ydot0, period, jacobi, n` and a `source: "arXiv:2509.12671"`
provenance field; assert P5g' values match the mining-note transcription. **RED**
(file absent).

**Minimal impl:** create `data/golden/png_hybrid_family.yaml` transcribing
Tables 3 + 5 of arXiv:2509.12671 (values already in
`2026-06-13-high-order-transfer-map-2509.12671-mining.md`). EXPECTED side traces
to the paper only (`feedback_golden_tests_sourced_only`).

**Verify:** `uv run pytest tests/data/test_png_golden.py -q` GREEN.
**Commit:** `data: #450 sourced Png' family golden (arXiv:2509.12671)`.

---

## Task 2 — Closure regression lock on the existing corrector

**Test first:** `tests/search/test_png_closure_golden.py` — for each Png' member
in the golden, feed `(x0, xdot0, jacobi, n)` to
`cr3bp_general_periodic.correct_general_periodic` at `cr3bp_system("Earth","Moon")`;
assert `converged`, residual ≤1e-11, recovered period matches golden to ≥8 sig
figs, recovered `x0/xdot0` drift <1e-7. **RED** if any member fails to close (also
flushes out per-member `half_crossings`/`n` mapping issues early). Keep in the
default suite (not `@slow`) — P5g' closes fast.

**Minimal impl:** thin test helper only (the corrector exists). Resolve the
crossing-index (`half_crossings`) → revolution-count `n` mapping per member; record
it in the golden.

**Verify:** GREEN for all 6 members (P5g' already proven to 3.45e-12).
**Commit:** `test: #450 Png' closure regression lock via correct_general_periodic`.

---

## Task 3 — `SectionMap` backend interface + sampling backend

**Test first:** `tests/genome/test_da_hotm_backend.py` —
(a) `SamplingSectionMap(system, c_target).single_rev(s)` returns the first y=0
return state for a section point `s=(x, xdot)` (ẏ from `ydot0_from_jacobi`),
cross-checked against a direct `core.cr3bp.propagate` to ~1e-10;
(b) `compose(s, n)` equals n sequential single-rev returns;
(c) a known DRO section point at C≈3.0002 returns near itself (residual small).
**RED**.

**Minimal impl:** `src/cyclerfinder/genome/da_hotm_backend.py` —
abstract `SectionMap` (methods `single_rev`, `compose(s, n)`, `residual(s, n)`),
and `SamplingSectionMap` built on `core.cr3bp.propagate` +
`cr3bp_periodic.ydot0_from_jacobi` + `jacobi_constant`. Mirror the structure of
`reachable_impulsive.py` (float propagator, no DA).

**Verify:** GREEN. **Commit:** `genome: #450 SectionMap interface + sampling backend`.

---

## Task 4 — Fixed-point enumerator over the (x, ẋ) domain

**Test first:** `tests/genome/test_da_hotm_enumerator.py` —
(a) on a coarse grid over a box containing the n=1 EM DRO at C=3.0002, the
enumerator emits ≥1 candidate within tol of the published DRO section point
`(x0≈0.885, …)`;
(b) candidates are de-duplicated (one representative per residual basin);
(c) returns coarse ICs as `(x0, xdot0, c_target, n)`. **RED**.

**Minimal impl:** `src/cyclerfinder/genome/da_hotm_enumerator.py` —
`enumerate_fixed_points(backend, domain_box, n, *, residual_tol, grid)`:
sample `residual(s, n)` over the box, isolate sub-tol cells, cluster (dedup), emit
coarse ICs. Backend-agnostic (consumes the `SectionMap` interface).

**Verify:** GREEN. **Commit:** `genome: #450 multi-rev fixed-point enumerator`.

---

## Task 5 — The decisive lane-recovery test: enumerate Png' P5g'

**Test first:** `tests/search/test_png_lane_recovery.py` (the design-draft §6
primary proof) — run `enumerate_fixed_points` on EM, C=3.00022, n=5, over the
published domain box; assert it emits a coarse candidate within tol of P5g'
`(x0≈0.8074, xdot0≈−0.0956)`; then feed that *enumerator-emitted* candidate (NOT
the golden IC) to `correct_general_periodic` and assert closure to the golden
P5g' period/residual. Mark `@slow` only if runtime forces it — prefer keeping the
single-member P5g' recovery in the default suite (evidence test). **RED**.

**Minimal impl:** none beyond Tasks 3–4 if the enumerator is correct; this is the
integration assertion that the *global sweep* surfaces the family. Add the domain
box + tolerance constants sourced from the paper's first planar case.

**Verify:** GREEN — the lane recovers a family seed-local continuation cannot.
**Commit:** `test: #450 lane-recovery proof — enumerator surfaces Png' P5g'`.

---

## Task 6 — Driver: filter cascade + novelty routing

**Test first:** `tests/search/test_da_hotm_enumeration.py` —
(a) base-family triangulation: on the untuned EM C≈3.0002 box, the driver recovers
the DRO + L1/L2 Lyapunov members and routes them to the **reproduction** bucket
(via `spatial_novelty_prefilter` / JPL-oracle / catalogue signature dedup);
(b) P5g' is routed to the **novel-PO** bucket (not reproduction), with
`check_literature` invoked (offline, necessary-not-sufficient);
(c) the result ledger carries provenance, residual, period, stability. **RED**.

**Minimal impl:** `src/cyclerfinder/search/da_hotm_enumeration.py` —
`run_enumeration(system, c_band, n_range, domain_box)`: loop `(C,n)` → enumerate →
section-residual gate → dedup → prefilter/JPL/catalogue screen →
`correct_general_periodic` → classify `{reproduction|known-family|novel-PO|
novel-cycler-candidate}`. Reuse existing modules; do NOT modify the corrector,
prefilter, or lit-check.

**Verify:** GREEN, including the base-family triangulation (the anti-circularity
guard). Run `uv run pytest tests/search tests/genome -q`.
**Commit:** `search: #450 DA/HOTM enumeration driver + novelty routing`.

---

## Task 7 — Negative-registry capability-subsumption stamp

**Test first:** extend the registry-consistency test (or add
`tests/data/test_empty_regions_da_hotm.py`) asserting that once the lane ships,
`cr3bp-em-cj3.00-dro-lyapunov-band-newfamily-2026-06-13` carries a `reverification`
/ method-version entry recording the subsuming method `da-hotm-enumeration-v1`
(does NOT delete the entry; capability-subsumption rule). **RED**.

**Minimal impl:** add the `reverification` stamp to that entry (and the three
Saturn-moon Lyapunov entries) in `data/empty_regions.jsonl` recording the new
method+git_sha and the Png'-recovery outcome (re-opened → Png' recovered as a PO,
no cycler). **Run the FULL ratchet suite** (`uv run pytest tests/data tests/search -q`)
per `feedback_catalogue_edits_run_all_ratchets`.

**Verify:** GREEN, full suite. **Commit:**
`data: #450 stamp C≈3.0 + Saturn Lyapunov negatives with DA/HOTM subsumption`.

---

## Task 8 (DEFERRED, gated) — DA-accelerated backend

**Gate:** only after Task 0's decision selects a DA path AND a measured sampling
sweep is too slow for the target production bands.

**Test first:** parity test — `DASectionMap.single_rev` / `compose` must agree with
`SamplingSectionMap` to ~1e-9 on the EM DRO and P5g' section points (the backend is
swappable iff it gives the same geometry). **RED** (backend absent).

**Minimal impl:** `DASectionMap` in `da_hotm_backend.py` using the chosen DA library
(DACEyPy Taylor map + ADS, or pure-Python truncated Taylor) behind the SAME
`SectionMap` interface. The enumerator/driver are unchanged (interface seam).

**Verify:** parity GREEN; re-run Tasks 4–6 tests with the DA backend selected;
measure speedup vs sampling. **Commit:** `genome: #450 DA Taylor-map SectionMap backend`.

---

## Sequencing & checkpoints

- **Tasks 0–2** (decision + golden + closure lock) are pure-validation scaffolding;
  ship first, low risk.
- **Tasks 3–5** build and *prove* the sampling lane (Task 5 is the decisive
  capability proof). After Task 5 the lane is demonstrably working.
- **Task 6** makes it a production driver with novelty routing.
- **Task 7** records the subsumption in the negative registry (the strategic
  payoff: dead region re-opened).
- **Task 8** is the optional DA-acceleration upgrade, fully gated and parity-tested.

Each task commits independently; the lane is usable and validated after Task 6
even if Task 8 is never funded.
