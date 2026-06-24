# #435 High-e Sun-planet ER3BP seeds — Plan

> REQUIRED SUB-SKILL on execution: superpowers:subagent-driven-development. Design approved 2026-06-24 (parallel with #433).

**Goal:** Generate CR3BP seed orbits at high-e Sun-planet μ (Mercury/Mars/Pluto) from scratch and feed them to the #432 ER3BP discovery pipeline, to probe the high-departure ER3BP regime that #432's Earth-Moon-only run could not reach.

**Architecture:** New `src/cyclerfinder/search/cr3bp_seed_generator.py` builds a converged CR3BP periodic IC at arbitrary μ from an analytic linear seed (collinear-point in-plane mode for a Lyapunov; retrograde co-orbital seed for a DRO) refined by the existing `correct_symmetric_fixed_jacobi`. A thin extension to `scripts/run_432_er3bp_discovery.py` Phase B consumes these as `Er3bpSeed`s and runs the already-built #432 pipeline (continuation + Floquet classification + adjudication) at the real Sun-planet eccentricities. Report-only; no catalogue writeback.

**Honest scope caveat (state in the verdict):** Lyapunov = libration family; DRO = co-orbital. These are *adjacent* to, not strictly, cycler-class — #435 probes the high-e ER3BP-periodic frontier (where the published ER3BP literature actually lives), answering "do high-e Sun-planet families bifurcate into novel ER3BP structure?".

**Confirmed API (use verbatim):**
- `cr3bp_system("Sun","Mercury"|"Mars"|"Pluto") -> CR3BPSystem` (.mu). Sun-Mercury μ=1.660e-7, Sun-Mars μ=3.227e-7, Sun-Pluto μ=7.350e-9.
- `lagrange_collinear_x(mu: float, point: str) -> float` in `src/cyclerfinder/search/reachable_representatives.py:395` (collinear L-point x; point in {"L1","L2","L3"}).
- `correct_symmetric_fixed_jacobi(system, x0_guess, jacobi, period_guess, *, ydot0_sign=-1.0, half_crossings=None, ...) -> SymmetricOrbit` in `src/cyclerfinder/search/cr3bp_periodic.py:204` (needs a good seed; refines a perpendicular-crossing symmetric orbit at fixed Jacobi).
- `jacobi_constant(state6, mu)` in `src/cyclerfinder/core/cr3bp.py:35`.
- `Er3bpSeed(label, system: ER3BPSystem, state0, period_f, is_half_period_residual, target_e, source)` + the #432 pipeline (`continue_and_monitor`, `adjudicate_trace`) in `src/cyclerfinder/search/er3bp_discovery.py`.

**Conventions:** work on main; `uv run` ruff+mypy before commit; no Co-Authored-By; pathspec commits; imports at top; subagents must finish through commit and not spawn their own reviewers.

---

## Task 1: Linear-seed CR3BP Lyapunov generator + golden

**Files:** Create `src/cyclerfinder/search/cr3bp_seed_generator.py`, `tests/search/test_cr3bp_seed_generator.py`

The planar Lyapunov linear seed at a collinear point `x_L`: the linearized in-plane dynamics give `Uxx = 1 + 2*c2`, `Uyy = 1 - c2` where `c2 = mu/|x_L-(1-mu)|^3 + (1-mu)/|x_L+mu|^3` (the standard collinear second-derivative coefficient). The in-plane oscillation frequency `omega = sqrt((2 - c2 + sqrt(9*c2^2 - 8*c2))/2)`; the small-amplitude planar Lyapunov IC is `state0 = [x_L - Ax, 0, 0, 0, -Ax*omega*tau, 0]` with `tau = -(omega^2 + Uxx)/(2*omega)` (the linear y-velocity/x-amplitude ratio), period guess `T = 2*pi/omega`. Refine with `correct_symmetric_fixed_jacobi(system, x0_guess=x_L-Ax, jacobi=jacobi_constant(state0,mu), period_guess=T/2, half_crossings=1)`.

- [ ] Step 1: failing test — `lyapunov_seed(system, point="L1", amplitude=1e-3) -> (state0: (6,), period: float)` returns a state that, propagated for `period`, returns to itself (closure residual < 1e-8) on the **Earth-Moon** system (the golden), and at the published EM L1 Lyapunov small-amplitude period (~2.69 TU for tiny amplitude — assert period in (2.0, 3.5) TU as the linear-regime band, and assert the corrected orbit's perpendicular re-crossing closes < 1e-9). (Sourced anchor: the collinear linear-Lyapunov frequency is textbook — Szebehely 1967 §; the test pins closure + linear-period band, not an unsourced exact IC.)
- [ ] Step 2: run → FAIL (import).
- [ ] Step 3: implement `lyapunov_seed` (linear seed + `correct_symmetric_fixed_jacobi` refine; return the corrected `(state0, full_period)`). Adapt to the real `SymmetricOrbit` return fields (read cr3bp_periodic.py:204).
- [ ] Step 4: run → PASS.
- [ ] Step 5: add `dro_seed(system, amplitude)` — retrograde near-secondary seed `state0=[(1-mu)+r, 0, 0, 0, -v_circ_retro, 0]` with `r=amplitude`, `v_circ_retro = sqrt(mu/r)+...` linear guess, refined by the same corrector; test it closes at Earth-Moon. (If DRO refinement proves finicky at this μ, mark `dro_seed` xfail with a sourced reason and proceed Lyapunov-only — note it; Lyapunov is the validated path.)
- [ ] Step 6: run → PASS (or xfail-noted for DRO).
- [ ] Step 7: ruff+mypy; commit `search/#435: CR3BP Lyapunov/DRO linear-seed generator at arbitrary mu`.

## Task 2: wire high-e Sun-planet seeds into the #432 runner

**Files:** Modify `scripts/run_432_er3bp_discovery.py` (or a new `scripts/run_435_high_e_er3bp.py` to avoid disturbing #432's committed runner — PREFER a new script).

- [ ] Create `scripts/run_435_high_e_er3bp.py` mirroring run_432's structure (offline_search + _print_progress). For each system in [("Sun","Mercury",0.206),("Sun","Mars",0.093),("Sun","Pluto",0.249)]: build the CR3BP system, generate a Lyapunov seed (and DRO if available) via Task 1, wrap as `Er3bpSeed(system=ER3BPSystem.from_cr3bp(cr3bp, e=target_e)... state0, period_f=full_period, is_half_period_residual=True, target_e=target_e)`, run `continue_and_monitor(seed, n_steps=60)` + `adjudicate_trace(...offline_search)`. Write `data/er3bp_discovery_435_highE.jsonl`. Print outcome + dv-free breakdown.
- [ ] Smoke (tiny: one system, n_steps=20) confirming a seed generates + continues without crash.
- [ ] ruff; commit `search/#435: high-e Sun-planet ER3BP discovery runner`.

## Task 3: deliverable run + verdict + registry (controller)

- [ ] Launch `scripts/run_435_high_e_er3bp.py` detached + harness-tracked waiter.
- [ ] Harvest → verdict note `docs/superpowers/plans/2026-06-24-435-high-e-er3bp-verdict.md` (per-system survives/dies/bifurcates; whether the high-departure regime yields a bifurcation = candidate novel family; honest caveat re Lyapunov/DRO not strictly cycler-class).
- [ ] Method-versioned `data/empty_regions.jsonl` entry (capability-extends the #432 EM-only negative to high-e Sun-planet). Run `uv run pytest tests/data -q -k "empty_region or registry"` before commit.

## Self-review
- Golden sourced (collinear linear-Lyapunov closure; no unsourced exact IC asserted). ✓
- Reuses #432 pipeline unchanged; new code is the seed generator + a runner. ✓
- Report-only, no catalogue writeback. ✓
- Caveat (Lyapunov/DRO adjacency to cycler-class) stated in verdict. ✓
