# #314 Heteroclinic-Cycle Framework — Design

**Date:** 2026-06-19. **Status:** design approved (brainstorm complete) → ready for writing-plans.

**Goal:** Discover and certify *closed heteroclinic cycles* — chains O₁→O₂→…→O₁ of
transversal invariant-manifold connections among equal-energy unstable orbits — as
**periodic-up-to-rotation** mass-transport loops. This is a new closure definition
distinct from the strict state(T)=state(0) periodicity every existing genome assumes.

**Regime (this phase):** planar CR3BP only. Sun-Jupiter-Oterma (μ=0.0009537, C=3.03)
is the validation system. The novel cross-system search is a separate, blocked task
(#405); the BCR4BP / SE↔EM extension is explicitly **out of scope here**.

**Sourced golden:** Wilczak & Zgliczyński, *Comm. Math. Phys.* (arXiv:math/0201278 Part I,
math/0401146 Part II) — the computer-assisted proof of the **closed** L1↔L2 Lyapunov
homo+heteroclinic cycle in the Sun-Jupiter-Oterma PCR3BP, with ~19-digit interval-enclosed
initial conditions on the section {y=0}. This is the only published candidate that is a
*true closed cycle* (both directions + homoclinic returns) with reproducible ICs. KLMR
Chaos 2000 (DOI 10.1063/1.166509) is the connection-only conceptual fallback (same μ, C;
no IC table). The golden dataset is transcribed by task #403 into
`data/golden/wz_oterma_heteroclinic.yaml`.

---

## Architecture

One new module, building on existing primitives. No edits to reused code.

### File structure
- **CREATE `src/cyclerfinder/genome/heteroclinic_cycle.py`** — the whole framework.
  Placement rationale: `genome/` is the home of modules that *define a new closure and
  corrector for a discoverable orbit class* (`qp_tori.py`, `tulip.py`, `bcr4bp_genome.py`).
  `search/` holds *scorers/sweeps over known orbits* (`resonance_network.py`,
  `lobe_overlap_scorer.py`). #314 certifies a new closure object → `genome/`. Its nearest
  sibling, the QP-tori genome (#290), lives in `genome/`. It *imports from* `search/` the
  same way `qp_tori` imports from `core/` — reuse-of ≠ belongs-with.
- **CREATE `tests/genome/test_heteroclinic_cycle.py`** — sourced-golden tests importing
  EXPECTED values from `data/golden/wz_oterma_heteroclinic.yaml` (#403). Location follows the
  established convention: genome-module tests live in `tests/genome/` (e.g.
  `tests/genome/test_qp_tori.py`).

### Reused primitives (no edits)
- `core/cr3bp.propagate(system, state6, t, *, with_stm, rtol, atol, stm_mode)` → `CR3BPArc`
  (STM/monodromy). `core/cr3bp.cr3bp_eom(t, state6, mu)` for independent re-propagation.
- `search/cr3bp_periodic.correct_symmetric_fixed_jacobi(...)` — generate a Lyapunov orbit
  at a target Jacobi (a Lyapunov orbit is a symmetric planar periodic orbit about a
  libration point). `ydot0_from_jacobi(x0, jacobi, mu, sign=...)` — IC seeding.
- `search/cr3bp_periodic.barden_stability(...)` — monodromy eigenstructure (Floquet
  multipliers + stable/unstable eigenvectors) for manifold seeding.
- `search/cr3bp_periodic.crosscheck_periodic(...)` — independent cross-check discipline.
- The manifold-seeding **pattern** from `search/resonance_network._planar_floquet` and
  `compute_floquet_manifold` (eigenvector perturbation → forward/backward integration →
  Poincaré-event detection). We adapt the section to {y=0} rather than the perilune section.

The pattern mirrors the QP-tori genome: a non-strict-closure object with a Newton corrector,
an independent cross-check, and explicit `converged` flags.

---

## Components

### 1. The connection primitive (core corrector)
A connection Wu(A) → Ws(B) is the classical KLMR section-intersection.

- **Equal-energy guard:** both orbits must share the same Jacobi C (necessary for a
  heteroclinic connection to exist). Mismatch → raise `ValueError`.
- Parameterize Wu(A) by phase τ_u around orbit A; perturb the IC by ε·v_unstable (from
  `barden_stability`); integrate **forward** to the k_u-th crossing of section Σ={y=0};
  record the section point (x, ẋ).
- Parameterize Ws(B) by τ_s; perturb by ε·v_stable; integrate **backward** to the k_s-th
  crossing of Σ; record (x, ẋ).
- **Residual** = 2-D section gap (Δx, Δẋ) between the two crossings. **Free variables** =
  (τ_u, τ_s). Newton with backtracking line-search on ‖residual‖ → transversal intersection.
- The Newton Jacobian is obtained by finite-differencing the section-crossing map over
  (τ_u, τ_s) (2×2) — cheap, robust; analytic STM-based sensitivity is a later optimisation,
  not needed for correctness.

`HeteroclinicConnection` dataclass: `(orbit_from, orbit_to, jacobi, tau_u, tau_s, k_u, k_s,
crossing_xv, residual, converged, n_iter, notes)`.

### 2. Manifold → section machinery
`manifold_section_crossings(system, orbit, *, direction, branch, section_y=0.0, n_samples,
eps, max_time, rtol, atol)` → an array of section crossings (x, ẋ) parameterized by phase.
`direction` ∈ {unstable (forward), stable (backward)}; `branch` ∈ {+1, −1} (which side of
the orbit). Event detection on y=0 with a sign filter on ẏ (Σ⁺/Σ⁻, matching W-Z's split).
Bounded by `max_time`; a sample that never reaches Σ is dropped with a count logged (no
silent truncation — per the negative-results / no-silent-caps discipline).

### 3. Cycle assembly & closure
`assemble_cycle(system, orbit_chain, jacobi, *, tol, …)` → `HeteroclinicCycle`. Certifies a
connection for each consecutive pair in `orbit_chain` (with the chain wrapping back to O₁)
and verifies the chain **returns to O₁**.

- **Closure = periodic-up-to-rotation:** the itinerary returns to the *same orbit* O₁
  (recurrence), with phase along O₁ free — explicitly NOT strict state(T)=state(0). This is
  the new closure definition. For W-Z Oterma: O₁=L1 Lyapunov, O₂=L2 Lyapunov; the
  connections L1→L2 and L2→L1 form the 2-cycle (homoclinic O→O legs are degenerate cycles
  of length 1 and also supported).
- `closed` is True iff every leg converged AND every orbit is at the shared `jacobi`.

`HeteroclinicCycle` dataclass: `(orbits, connections, jacobi, closed, max_leg_residual,
independent_residual, symbol_sequence, notes)`.

### 4. Independent cross-check
Re-propagate every converged leg with the Radau integrator (vs the DOP853 default used in
the corrector) and confirm the section crossing agrees to `independent_tol`. Stored in
`HeteroclinicCycle.independent_residual`. Mandatory before any "closed" claim
(orbit-closure-discipline: independent cross-check is non-negotiable).

---

## Validation ladder (sourced golden, #403)

All EXPECTED values trace to the W-Z publication, never to our own computation.

1. **Conventions self-check.** Reproduce the two Lyapunov fixed points at W-Z's printed
   (x, 0) on {y=0} via `correct_symmetric_fixed_jacobi` at C=3.03. EXPECTED:
   x*₁=0.9208034913207400196, x*₂=1.081929486841799903 (W-Z Part I). Confirms our μ/C/section
   conventions match the paper before any connection is attempted.
2. **Per-leg connection.** Certify L1→L2 and L2→L1 at the tabulated crossing coordinates.
   EXPECTED = W-Z crossing points (from `data/golden/wz_oterma_heteroclinic.yaml`).
3. **Chain closure.** `assemble_cycle([L1, L2], jacobi=C_oterma)` returns `closed=True`,
   `max_leg_residual < tol`.
4. **Independent cross-check.** Radau-vs-DOP853 endpoint agreement `< independent_tol`.

Test file: `tests/search/test_heteroclinic_cycle.py`. The golden YAML is a #403 deliverable;
the #314 corrector build and the #403 transcription run in parallel — the validation tests
(items 2–4) are written against the golden but only pass once #403 lands. Item 1 needs only
the two fixed points (already captured by the scouting pass) and can be written immediately.

---

## Error handling & edge cases

- **Energy mismatch** between chain orbits → `ValueError` (connections cannot exist off
  equal-C).
- **Manifold misses Σ** within `max_time` → that sample dropped, count logged; if NO sample
  reaches Σ, the connection returns `converged=False` with a diagnostic (no hang — the
  lesson from #380's pathological Jacobi-drift quadrature).
- **Non-transversal / no intersection** → `converged=False`, never a fabricated closure.
  A clean negative is an acceptable outcome (a chain that does not close is reported as such).

---

## Scope guard (YAGNI)

This phase delivers the **planar CR3BP corrector + cycle assembler + W-Z validation only**.

- The novel cross-system closed-cycle search is **#405** (blocked on this), and the
  literature confirms no published cross-system closed cycle exists → it is a
  clean-negative-allowed novel search, not a reproduction.
- Lobe/flux *rate* scoring stays in the existing #267/#278 scorers — not rebuilt here.
- BCR4BP / ER3BP / 3D extensions are future phases, not in this design.

## Build sequence (for writing-plans)

1. `manifold_section_crossings` + the {y=0} event machinery + a unit test on a single
   Lyapunov orbit's manifold reaching the section.
2. Conventions self-check test (validation-ladder item 1) — Lyapunov fixed points vs W-Z.
3. `correct_connection` (the 2-D Newton primitive) + per-leg golden test (item 2).
4. `assemble_cycle` + closure definition + chain-closure test (item 3).
5. Independent Radau cross-check (item 4).
6. Public API surface review + docstrings citing W-Z.
