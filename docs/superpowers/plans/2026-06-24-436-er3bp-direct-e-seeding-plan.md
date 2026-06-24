# #436 ER3BP direct-e>0 seeding + branch-switching — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: superpowers:subagent-driven-development. Steps use `- [ ]` checkboxes.

**Goal:** Hunt a novel ER3BP *e>0-only* periodic family (one with NO CR3BP limit) via two mechanisms the #432/#435 continuation campaigns could not exercise: (A) direct e>0 seeding from CR3BP-independent ICs + reverse-continuation to e=0 to prove no-CR3BP-ancestor, and (B) an ER3BP branch-switcher (built but contingent — #432/#435 flagged zero bifurcations, so it has no current target).

**Architecture:** Reuse the #293 ER3BP corrector (`correct_er3bp_periodic`) and the #432 continuation (`continue_er3bp_family_in_e`) unchanged. New: a CR3BP-INDEPENDENT seed source (a coarse symmetric `(x0, ydot0)` grid at fixed period-multiple, converged directly at target e — NOT extrapolated from a CR3BP orbit, which would by construction reverse-continue back to its CR3BP ancestor), a reverse-to-e=0 classifier, a campaign runner, and an ER3BP analogue of `branch_at_saddle_center` (fixed-period: the ER3BP period in true anomaly is locked to a multiple of 2π, so the CR3BP free-period 7-var corrector cannot be reused). Report-only; no catalogue writeback.

**Critical design correction (vs scout proposal):** seeds MUST be independent of any CR3BP family. A seed extrapolated from a CR3BP orbit always has a CR3BP limit, so it can only ever classify as "cr3bp_continuous" — it can never be an e>0-only discovery. The genuine mechanism is a blind symmetric-IC grid converged at target e.

**Honest likelihood:** LOW. The published ER3BP record documents no rich class of e>0-only no-CR3BP-limit families; a blind grid against a 2-var symmetric corrector is low-yield. Expect a clean, registry-grade negative. The value is closing the last untested ER3BP discovery mechanism.

**Conventions:** work on main; `uv run` ruff + mypy before commit; no Co-Authored-By; pathspec commits; imports at top; subagents finish through commit and do not self-spawn reviewers.

---

## Task 1: CR3BP-independent direct-e>0 seed grid

**Files:** Create `src/cyclerfinder/search/er3bp_direct_seeding.py`, `tests/search/test_er3bp_direct_seeding.py`

- [ ] Step 1: Write failing test. `direct_e_seed_grid(system, x0_range, ydot0_range, n_x, n_ydot, period_f, is_half_period_residual=True) -> list[DirectEr3bpSeed]` returns `n_x * n_ydot` seeds, each a `DirectEr3bpSeed(label, system, state0(6,), period_f, is_half_period_residual, target_e, source)` with `state0 = [x0, 0, 0, 0, ydot0, 0]` (symmetric x-axis-crossing IC), `target_e == system.e`, and `source` a non-empty provenance string naming the grid (NOT a CR3BP family). Assert grid covers the corners.

```python
def test_direct_e_seed_grid_corners():
    sys = ER3BPSystem(mu=0.012155, e=0.0549, primary_name="E", secondary_name="M")
    seeds = direct_e_seed_grid(sys, (0.1, 0.9), (-3.5, 3.5), n_x=3, n_ydot=3, period_f=2*np.pi)
    assert len(seeds) == 9
    assert all(s.state0.shape == (6,) for s in seeds)
    assert all(s.state0[1] == 0 and s.state0[3] == 0 and s.state0[5] == 0 for s in seeds)
    assert all(s.target_e == 0.0549 for s in seeds)
    assert seeds[0].state0[0] == 0.1 and seeds[-1].state0[0] == 0.9
    assert "grid" in seeds[0].source.lower() and "broucke" not in seeds[0].source.lower()
```

- [ ] Step 2: Run → FAIL (import).
- [ ] Step 3: Implement `DirectEr3bpSeed` dataclass + `direct_e_seed_grid` (np.linspace over both ranges, build symmetric ICs). `source = f"#436 CR3BP-independent symmetric-IC grid x0∈{x0_range} ydot0∈{ydot0_range}"`.
- [ ] Step 4: Run → PASS.
- [ ] Step 5: ruff + mypy; commit `search/#436: CR3BP-independent direct-e>0 seed grid`.

## Task 2: forward-converge + reverse-to-zero classifier

**Files:** Modify `src/cyclerfinder/search/er3bp_direct_seeding.py`; test in same test file.

The corrector API (verbatim): `correct_er3bp_periodic(system, state_guess, period_f, *, free_vars=(IDX_X, IDX_YDOT), residual_indices=(IDX_Y, IDX_XDOT), is_half_period_residual=True, tol=1e-10, max_iter=60, ...) -> ER3BPPeriodicOrbit` (fields incl. `.state0`, `.corrector_residual`, `.converged`). Continuation: `continue_er3bp_family_in_e(sys_base, seed_state, period_f, e_target, n_steps, *, is_half_period_residual=True, tol=1e-10) -> list[ER3BPPeriodicOrbit]` (raises `ContinuationError` on step failure). Import `IDX_*` from the corrector module.

- [ ] Step 1: Write failing test for `converge_direct_seed(seed, *, tol=1e-10) -> ER3BPPeriodicOrbit | None` (returns converged orbit at `seed.system.e`, else None) and `classify_no_cr3bp_limit(converged_orbit, system, *, n_steps=30, death_floor=1e-3) -> dict`. The classifier reverse-continues from `system.e` toward `e_target=0.0`; returns `{"status": "e_only_candidate", "death_e": <e>}` if continuation dies at `death_e > death_floor`, `{"status": "cr3bp_continuous"}` if it reaches ~0, `{"status": "inconclusive", ...}` if it dies below `death_floor` (numerical noise). GOLDEN (sourced, self-consistent): seed the Broucke-1969 Earth-Moon 7P IC converged at e=0.0549 → it MUST classify `cr3bp_continuous` (it has a CR3BP limit by construction — this validates the classifier's negative branch against a known CR3BP-continuous family).

```python
def test_broucke_classifies_cr3bp_continuous():
    sys = ER3BPSystem(mu=0.012155, e=0.0549, primary_name="E", secondary_name="M")
    seed = DirectEr3bpSeed("broucke-check", sys,
        np.array([0.1520965,0,0,0,3.1608994,0]), 2*np.pi, True, 0.0549, "test")
    orb = converge_direct_seed(seed)
    assert orb is not None
    res = classify_no_cr3bp_limit(orb, sys, n_steps=20)
    assert res["status"] == "cr3bp_continuous"
```

- [ ] Step 2: Run → FAIL.
- [ ] Step 3: Implement both. `converge_direct_seed` wraps the corrector in try/except, returns None on non-convergence. `classify_no_cr3bp_limit` calls `continue_er3bp_family_in_e(..., e_target=0.0, ...)` in try/except: on success read final `.e` (≈0 → cr3bp_continuous); on `ContinuationError` parse/estimate the death eccentricity by re-running with a small `n_steps` ladder or catching the last-converged e (simplest: re-run continuation and track the last successful step's e via a helper that returns partial history — add `continue_er3bp_family_in_e_partial` returning `(orbits, died_at_e | None)` to er3bp_continuation.py if needed, OR wrap by catching ContinuationError and reading the e from its message is fragile — prefer adding a partial-history variant). Implementer: choose the cleanest robust route and note it.
- [ ] Step 4: Run → PASS.
- [ ] Step 5: ruff + mypy; commit `search/#436: forward-converge + reverse-to-e0 no-CR3BP-limit classifier`.

## Task 3: ER3BP branch-switcher (fixed-period analogue)

**Files:** Create `src/cyclerfinder/genome/er3bp_branching.py`, `tests/genome/test_er3bp_branching.py`

Reuse `_select_saddle_center_eigenvector(mono, *, parent_tangent)` from `asymmetric_branch.py` and `er3bp_monodromy(state0, period_f, system)` from `er3bp_floquet.py`. The ER3BP period is FIXED (true-anomaly multiple of 2π), so re-convergence uses `correct_er3bp_periodic` with the SAME `period_f` (NOT the CR3BP free-period corrector). Parent tangent = `er3bp_eom(0.0, state0, mu, e)` (confirm signature in core/er3bp.py).

- [ ] Step 1: Write failing test. `branch_at_saddle_center_er3bp(system, parent_state0, parent_period_f, *, epsilon=1e-3, tol=1e-10, ...) -> tuple[ER3BPPeriodicOrbit | None, dict]`. Perturbs ±epsilon along the marginal eigenvector (with the eps/sign ladder pattern from `branch_at_saddle_center`), re-converges at fixed `period_f`. Test: on a NON-bifurcating parent (Broucke e=0.0549, which #432 showed stays hyperbolic) the function returns either `(None, {...})` or a converged orbit that is the SAME family (cross-check state0 ≈ parent within tol) — i.e. it must not fabricate a spurious branch. Assert it does not raise and returns the documented tuple type.
- [ ] Step 2: Run → FAIL.
- [ ] Step 3: Implement per the scout's Section 4.2 design (monodromy → eigenvector pick → eps/sign ladder → `correct_er3bp_periodic` at fixed period). Compound gate: `corrector_residual < tol` AND independent closure.
- [ ] Step 4: Run → PASS.
- [ ] Step 5: ruff + mypy; commit `search/#436: ER3BP fixed-period saddle-center branch-switcher (contingent infra)`.

## Task 4: campaign runner + verdict + registry (controller-run)

**Files:** Create `scripts/run_436_direct_er3bp.py`

- [ ] Build the runner: for systems [Earth-Moon e=0.0549, Sun-Mercury 0.206, Sun-Mars 0.093, Sun-Pluto 0.249], build a `direct_e_seed_grid` (start ~12×12 over a physically sensible `x0∈(0.1,0.95)`, `ydot0∈(-4,4)` band at full period 2π and half-period π), `converge_direct_seed` each, `classify_no_cr3bp_limit` survivors, adjudicate any `e_only_candidate` via the offline KNOWN_CORPUS `check_literature` (reuse run_435's `offline_search`). Write `data/er3bp_direct_436.jsonl`. Print per-system: grid size, n converged, n cr3bp_continuous, n e_only_candidate, n inconclusive. For any bifurcation/e_only_candidate, also try `branch_at_saddle_center_er3bp`. `_print_progress` with timestamps.
- [ ] Smoke: one system, 3×3 grid — confirm it converges some, classifies, writes records without crash.
- [ ] ruff + mypy; commit `search/#436: direct-e>0 ER3BP discovery runner`.
- [ ] Launch detached + harness-tracked waiter (controller).
- [ ] Harvest → verdict `docs/superpowers/plans/2026-06-24-436-direct-er3bp-verdict.md` (per-system counts; any e_only_candidate + lit status; honest framing — a clean "no e>0-only no-CR3BP-limit family found at this grid" is a registry-grade negative that, with #432/#435, characterises the ER3BP discovery frontier across all three mechanisms: continuation-from-CR3BP, direct-e>0, and branch-switching).
- [ ] Method-versioned `data/empty_regions.jsonl` entry (capability-extends the #432/#435 ER3BP negatives with the direct-seeding mechanism). `uv run pytest tests/data -q -k "empty_region or registry"` before commit.

## Self-review
- Seeds are CR3BP-INDEPENDENT (the scout's extrapolation flaw corrected). ✓
- Golden is sourced/self-consistent (Broucke → cr3bp_continuous validates the classifier; no unsourced exact IC asserted). ✓
- Branch-switcher built but honestly flagged contingent (no current bifurcation target). ✓
- Report-only; no catalogue writeback. ✓
