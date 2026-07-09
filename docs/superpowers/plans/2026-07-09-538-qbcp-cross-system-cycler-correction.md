# #538 QBCP Cross-System Periodic Orbit (Cycler) Correction — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Chain the forward (SE-L2→EM-L2) and reverse (EM-L2→SE-L2) QBCP torus connections into a 4-segment boundary-value problem and converge it onto a mathematically exact, synodic-periodic closed orbit (a candidate cislunar cycler) under the genuine time-periodic QBCP model — or, failing that, a rigorously documented clean negative.

**Architecture:** New script `scripts/run_538_qbcp_cycler.py` reusing `core/qbcp.py` (EOM+STM), `genome/qbcp_torus.py` (`QBCPTorus`, `correct_qbcp_torus`, `evaluate_qbcp_torus`) and the #537 candidate (`scripts/run_533_qbcp_connection.py`) as a **starting seed only, not a converged boundary condition** (see Context below — this matters). The corrector is a single **well-posed least-squares Newton solve** over all four torus-phase pairs `θ_i=(θ_long,i, θ_trans,i)`, the four epochs `t_i`, and the two open segment durations `(τ_f, τ_r)`, minimizing a residual that includes **position AND velocity** at both crossing sections plus the synodic-periodicity condition — unlike #537's diagnostic solve, which left velocity and one degree of freedom uncontrolled.

**Tech Stack:** Python 3, numpy, scipy.integrate.solve_ivp (DOP853 primary, Radau independent cross-check), scipy.optimize.least_squares, pytest, uv, ruff. Prior art: `genome/heteroclinic_cycle.py` (#314) and `genome/cross_system_cycle.py` (#405)'s Newton+coarse-scan pattern; #531's `crosscheck_cycle` discipline (independent-integrator re-verification before trusting a "closed" result).

---

## Context — what #537 actually handed us (read before Task 1)

`run_533_qbcp_connection.py` (the #537 script) sweeps SE-L2-torus unstable-manifold ×
EM-L2-torus stable-manifold crossings at section `x=2.0`, coarse-filters candidate pairs,
then refines the top 5 via `least_squares` on a **3-equation residual**
`(Δy, Δz, Δt mod T_s)` against **4 unknowns**
`(θ_long_u, θ_trans_u, θ_long_s, θ_trans_s)`. Two consequences:

1. **The solve is rank-deficient (3 eqs / 4 unknowns).** It converges onto *some* point on
   a 1-parameter family of near-solutions, not a unique connection. The reported
   12,034 km position gap and 911 m/s velocity gap are **printed post-hoc diagnostics from a
   separate propagation**, not quantities the optimizer ever drove toward zero.
2. **Velocity was never in the residual.** A 911 m/s leftover gap is exactly what you'd
   expect from a corrector that only ever asked for position+time agreement.

**Do not treat #537's four `θ_i` values as a tight seed that merely needs polishing.**
Use them as a **basin seed** for a properly-posed solve (this plan's Task 2), and be
prepared for the Newton iteration to walk noticeably away from the #537 numbers — that is
expected and correct, not a sign of a broken seed. This directly motivates why the plan
below (per this project's S1L1-saga lesson that "it closed!" on an under-constrained
residual is exactly the failure mode that produces false positives) tries
**ballistic-first** and only falls back to **powered** if ballistic genuinely fails to
converge in-basin.

`scripts/run_522_coherent_connection.py` and `scripts/search_coherent_connections.py` are
earlier BCR4BP-only precursors to `run_533`'s genuine-QBCP approach — **superseded, not
alternatives**. Leave them in place (git history documents the progression) but do not
build on them; #538 builds only on `run_533` + `core/qbcp.py` + `genome/qbcp_torus.py`.

---

## Conventions every task must follow
- Run from repo root `/home/bruce/dev/cyclers` with `uv run` (uv venv; never pip).
- **Pre-commit mandatory; never `git commit --no-verify`.** Commit with explicit pathspecs.
- **Sourced-golden discipline:** any EXPECTED numeric value in a test traces to a
  publication or an internal-consistency identity, never to our own computed output.
- Type-annotate everything (mypy strict, matching `core/qbcp.py`/`genome/qbcp_torus.py`).
- Expensive convergence/search steps get `@pytest.mark.slow`.
- **Every `scripts/run_*.py` MUST call `preflight_search()` near the top of `main()`** —
  see Task 0. This is not optional flavor; `tests/scripts/test_scripts_call_preflight.py`
  is an AST ratchet that fails the build otherwise.

---

## Task 0 (blocking, do first): fix the currently-broken preflight ratchet

**Finding (2026-07-09):** `tests/scripts/test_scripts_call_preflight.py::
test_non_exempt_scripts_call_preflight_search` **currently fails on `main`** —
`run_534_torus_connection.py` and `run_536_jupiter_europa_connection.py` (already
committed) call no `preflight_search()` and are not in `_LEGACY_EXEMPT`. The two
uncommitted scripts `run_522_coherent_connection.py` and `run_533_qbcp_connection.py`
have the same gap. `scripts/run_538_qbcp_cycler.py` must not add a fifth offender.

**Files:**
- Modify: `scripts/run_534_torus_connection.py`, `scripts/run_536_jupiter_europa_connection.py`,
  `scripts/run_522_coherent_connection.py`, `scripts/run_533_qbcp_connection.py`
- Reference pattern: `scripts/run_515_cross_system_3d_search.py:21,44-55,141-147`

- [ ] **Step 1:** Reproduce the failure: `uv run pytest tests/scripts/test_scripts_call_preflight.py -v` — expect the assertion listing all four filenames.
- [ ] **Step 2:** For each of the four scripts, add near the top of `main()`:
  ```python
  from cyclerfinder.data.preflight import PreflightBlockedError, preflight_search

  preflight_search(
      task_no=<536|534|522|533>,
      region_id="<slug, e.g. qbcp-se-l2-em-l2-connection-2026-07-08>",
      method=MethodCapability(
          genome="<one-line description of the connection genome>",
          corrector="<one-line corrector description>",
          capability_tags=frozenset({"qbcp", "torus", "heteroclinic", ...}),
          git_sha="working-tree",
      ),
      script_path=pathlib.Path(__file__),
      n_points=<actual sweep size>,
  )
  ```
  Wrap the sweep body in `try/except PreflightBlockedError` only if the script is meant to
  be re-runnable safely; otherwise let it raise (that IS the gate working).
- [ ] **Step 3:** Re-run the ratchet: `uv run pytest tests/scripts/test_scripts_call_preflight.py -v` — expect PASS, 2/2.
- [ ] **Step 4:** Run the full data/search ratchet to confirm no regression: `uv run pytest tests/data tests/search tests/scripts -q`.
- [ ] **Step 5:** Commit with pathspecs (the four modified scripts only).

**Recommended model: Haiku (or Sonnet at low effort).** Purely mechanical — a fixed
template applied to four files, with the AST ratchet itself as the deterministic gate that
catches any mistake. No judgment calls; nothing here should burn Opus budget.

---

## Task 1: Well-posed multi-segment residual formulation

**The crux of this task.** Design (on paper / in a design comment, before writing the
solver) the free-parameter vector and residual vector so the system is square or
over-determined, never under-determined like #537's:

- **Unknowns (12):** `θ_0, θ_1, θ_2, θ_3` (4 torus phase-pairs, 8 scalars), `t_0` (departure
  epoch), `τ_f, τ_r` (segment 1 and 3 durations), plus one extra scalar released only if
  needed to keep the system square (e.g. an EM-torus stable-manifold epsilon or the exact
  `x`-crossing value, documented explicitly — do not silently drop the velocity residual to
  make the count work).
- **Residual (≥12):** position (3) + velocity (3) match at the segment-1→2 crossing;
  position (3) + velocity (3) match at the segment-3→4 crossing; the synodic-periodicity
  scalar `(t_3+τ_{return}-t_0) mod T_S = 0`; and the torus-phase return condition
  `θ_3 → θ_0` (as many components as needed to close the loop — do not under-count this
  the way #537 did).

**Files:**
- Create: `scripts/run_538_qbcp_cycler.py` (design/residual-shape comment block only in
  this task; implementation in Task 2)
- Test: `tests/scripts/test_run_538_residual_shape.py` — a cheap unit test asserting
  `len(residual) >= len(unknowns)` for the chosen parameterization (prevents silently
  regressing to an underdetermined system in a future edit).

- [ ] **Step 1:** Write `_residual_shape_ok()` returning `(n_unknowns, n_residuals)` from the chosen parameterization; assert `n_residuals >= n_unknowns` in a test.
- [ ] **Step 2:** Run the test, confirm it fails (function doesn't exist yet).
- [ ] **Step 3:** Implement the function + constants.
- [ ] **Step 4:** Run the test, confirm PASS.
- [ ] **Step 5:** Commit.

**Recommended model: Opus (high effort).** This is the judgment call that #537 got wrong
(a rank-deficient residual that "converged" without actually closing anything). Getting
the parameterization right is trust-bearing numerical-methods judgment, not mechanical
coding — exactly the case the project's model-tiering policy reserves for Opus.

---

## Task 2: Ballistic corrector (attempt first)

**Files:**
- Modify: `scripts/run_538_qbcp_cycler.py`

- [ ] **Step 1:** Seed the 12 unknowns from #537's four `θ_i`, its `t_0`/crossing time (3.66 TU), and equal-split guesses for `τ_f`/`τ_r`.
- [ ] **Step 2:** Implement `_full_residual(params) -> np.ndarray` propagating all 4 segments via `qbcp.propagate_qbcp_pv` (or `qbcp_torus.evaluate_qbcp_torus`) and returning the Task-1 residual vector (position+velocity at both crossings + periodicity + phase-return).
- [ ] **Step 3:** Run `scipy.optimize.least_squares(_full_residual, x0, ...)` (ballistic: no ΔV terms).
- [ ] **Step 4:** Log per-iteration residual norm to a runlog file (`data/runlogs/run_538_qbcp_cycler.jsonl`, append+flush) — per the project's "instrument long runs" convention; this solve may take many iterations from a walked-away seed.
- [ ] **Step 5:** If converged residual norm < 1e-8: record success, proceed to Task 4 (cross-check). If not converged after a reasonable iteration budget (document the budget and why): proceed to Task 3 (powered fallback) rather than loosening the tolerance.
- [ ] **Step 6:** Commit (script + runlog schema, not the runlog data itself if gitignored — check `.gitignore`).

**Recommended model: Sonnet (medium-high effort).** Task 1 has already made the hard
judgment call; this is implementation against a spec-complete residual design, backed by
a deterministic numeric-tolerance gate (`< 1e-8`). Exactly the "spec-complete TDD behind a
strong gate" case for Sonnet.

---

## Task 3: Powered fallback (only if Task 2 does not converge ballistically)

**Files:**
- Modify: `scripts/run_538_qbcp_cycler.py`

- [ ] **Step 1:** Add a ΔV unknown at the `x=2.0` crossing (magnitude + direction, or a
  full 3-vector) and relax the velocity-match residual there to `v_after - v_before - ΔV = 0`.
- [ ] **Step 2:** Re-run `least_squares`, now minimizing `|ΔV|` as an explicit objective
  term (or a two-stage solve: converge position/velocity-with-ΔV-free first, then minimize
  `|ΔV|` holding the closure constraints active).
- [ ] **Step 3:** Log and compare `|ΔV|` against existing catalogue `dv_band` tiers (see
  `verify/dv_band_acceptance.py`) to sanity-check whether a powered result would even be
  catalogue-eligible before investing further.

**Recommended model: Opus (high effort).** Powered fallback changes the objective
structure (an optimization, not just a root-find) and requires judging whether a
converged-but-expensive ΔV is a meaningful result or a sign the connection doesn't
actually exist ballistically — a verdict call, not mechanical implementation.

---

## Task 4: Independent cross-check (mandatory before calling anything "closed")

Per the project's orbit-closure discipline (the S1L1 saga's reusable lesson: "it closed!"
is the danger signal, not the finish line): a converged residual alone does not certify a
real periodic orbit. Mirror #531's `crosscheck_cycle` pattern.

**Files:**
- Create: `tests/scripts/test_run_538_crosscheck.py`

- [ ] **Step 1:** Re-propagate the converged solution with an independent integrator
  (`method="Radau"` instead of `DOP853`) end-to-end through all 4 segments; assert the
  final-state drift vs the DOP853 solution is below a stated tolerance (document the
  tolerance and why, not a round number pulled from nowhere).
- [ ] **Step 2:** Re-derive the periodicity check from the raw propagated states (not from
  the optimizer's own residual output) — an independent recomputation, not a re-read of the
  same number.
- [ ] **Step 3:** If cross-check fails: the result is NOT a confirmed closed orbit — record
  as an honest negative (Task 5's negative-registry branch), do not soften the tolerance to
  make it pass.

**Recommended model: Sonnet (medium effort) for the harness; Opus (high effort) for the
final verdict** — i.e., have Sonnet build the cross-check machinery against a clear spec,
but the human/Opus call on "does this cross-check result mean we have a real orbit" should
not be delegated to a cheaper tier, since a wrong verdict here becomes wrong science that
nothing downstream mechanically catches.

---

## Task 5: Utility analysis + writeback

**Files:**
- Modify: `scripts/run_538_qbcp_cycler.py`
- Modify (only if confirmed): `data/catalogue.yaml`, `data/OUTSTANDING.md`
- Modify (only if honest negative): `data/negative_results.yaml` /
  `data/empty_regions.jsonl` (per `should_sweep()` / #521's registry)

- [ ] **Step 1:** Compute close-approach distance to Earth and Moon (Hill-radius units)
  along the converged (and cross-checked) trajectory.
- [ ] **Step 2:** Compute the monodromy matrix (4-segment composite STM) and its
  eigenvalues; classify stability.
- [ ] **Step 3a (if Task 4 confirms closure):** Draft a catalogue row (`source: discovered`,
  `first_published: cyclerfinder 2026`) but **hold it at the appropriate low validation
  tier (V0/V1) pending an independent Fable-model adversarial review** of the closure
  claim before any writeback — this project's own established pattern (the 2026-07-02
  independent Fable audits) for keeping a fresh, differently-biased model from rubber-
  stamping a result the primary implementer is invested in.
- [ ] **Step 3b (if Task 4 falsifies closure, or Tasks 2+3 both fail to converge):** Record
  a clean, method-versioned negative in `data/negative_results.yaml` and
  `data/empty_regions.jsonl` with the specific `region_id` used in Task 0, per the "clean
  negative is success" convention. This is a fully acceptable, valuable outcome.
- [ ] **Step 4:** Update the `#538` entry in `data/OUTSTANDING.md` with whichever outcome
  landed, and mark it `✓ Resolved (date)`.

**Recommended model: Sonnet** for the mechanical utility calculations and writeback
formatting; **Opus** for drafting the novelty/catalogue-eligibility framing (mirrors the
literature-novelty-check gate: nothing is called "novel" until it clears that check, and
the framing of what got found is a judgment call); **Fable** for the mandatory independent
adversarial review pass in Step 3a before any catalogue writeback.

---

## Verification Plan

### Automated Tests
- `uv run pytest tests/scripts/test_scripts_call_preflight.py -v` — must pass before any
  other task lands (Task 0).
- `uv run pytest tests/scripts/test_run_538_residual_shape.py -v` — residual well-posedness
  regression guard (Task 1).
- Run the corrector: `uv run python scripts/run_538_qbcp_cycler.py`
- Verify the **cross-checked** (not just DOP853-only) matching residual is driven below
  `1e-8` (Task 4), not the single-integrator residual alone.
- `uv run pytest tests/scripts/test_run_538_crosscheck.py -v`
- Full ratchet: `uv run pytest tests/data tests/search tests/scripts -q` (must stay green;
  this is a catalogue.yaml-touching change if Task 5 writes back a row — run the full
  ratchet, not a hand-picked subset).
- Style: `uv run ruff check scripts/run_538_qbcp_cycler.py` /
  `uv run ruff format --check scripts/run_538_qbcp_cycler.py`

### Honest-negative acceptance criterion
A clean negative (Task 5, Step 3b) is a fully valid, complete outcome for this plan. Do
**not** loosen tolerances, relax the residual definition from Task 1, or drop the Task 4
cross-check in order to force a "closed" result — that is the exact failure pattern this
plan's Context section documents #537 as having narrowly avoided being mistaken for.
