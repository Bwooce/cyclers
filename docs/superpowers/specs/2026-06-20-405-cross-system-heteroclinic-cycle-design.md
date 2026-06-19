# #405 Cross-System Heteroclinic-Cycle Search (Phase A) — Design

**Date:** 2026-06-20. **Status:** design approved (brainstorm complete) → ready for writing-plans.

**Goal:** Search for a CLOSED cross-system heteroclinic cycle — a repeating orbit that
transits between the Sun-Earth (SE) and Earth-Moon (EM) CR3BP regimes via libration-point
manifold tubes and returns to its start (periodic-up-to-rotation) — using a patched-CR3BP
("Shoot the Moon") model. Build the inter-system connection matcher, validate it, then run
a bounded closure search. A closed cycle OR a documented clean negative is an acceptable
outcome (the object is genuinely novel — no published closed cross-system cycler exists,
per the #316 survey).

**Model (decided):** STAGED. Phase A (this spec) = patched-CR3BP: two coupled `core/cr3bp`
instances (SE, EM) matched at a Poincaré section in a common frame. Phase B (deferred,
separate task) = refine any found cycle in the BCR4BP (#292). BCR4BP / broad multi-libration
campaigns are explicitly OUT of scope here.

**Prior art / why novel:** `docs/notes/2026-06-17-316-cross-system-cycler-framework.md`
(18-paper survey) — published work covers two adjacent axes: (1) one-shot SE↔EM heteroclinic
TRANSFERS (KLMR 2001 "Shoot the Moon"; Canalias-Gómez-Marcote-Masdemont 2006; van der
Weg-Vasile 2015), never iterated into a cycler; (2) BCR4BP synodic-resonant periodic orbits
that close in both frames but stay EM-Hill-bounded (Howell school). The union — a repeating
orbit whose support spans the EM Hill sphere AND the SE-L1/L2 region within one period — is
unmined. So the *connection primitive* has a published anchor; *closure* is the novel search.

**Builds on:** #314 `genome/heteroclinic_cycle.py` (validated against Wilczak-Zgliczyński):
reuses `LyapunovNode`, `_seed_on_manifold`, `_section_crossing`, and the 2-D Newton +
coarse-scan connection-corrector pattern. Reuses `core/cr3bp` (propagate + STM + Jacobi).

---

## Validation strategy (two-tier — per the Phase-0 scout)

The genuinely-new inter-system code has no tabulated published golden (no source prints a
machine-precision cross-system patch state vector). So:

1. **Sourced numeric tier** (EXPECTED traces to publication):
   - SE-L1/L2 Lyapunov orbits + single-system L1↔L2 connection reproduce the
     Canalias-Gómez-Marcote-Masdemont per-family Jacobi values (confirmed **C = 3.000863625**
     for the SE L1→L2 heteroclinic family bifurcation). Pins SE μ / Jacobi convention / section.
   - EM-side per-system pieces reuse the #314 W-Z-validated path.
   - Reproduce a KLMR-class one-way SE-L2→EM-L2 ballistic connection: the patch ΔV must be
     low-energy (published figures are ~tens of m/s), confirming the matcher finds a real
     near-ballistic connection.
2. **Internal-consistency tier** (falsifiable identities for the new bridge, where no golden
   exists):
   - **Frame round-trip identity:** SE-rotating → inertial → EM-rotating → inertial →
     SE-rotating returns the state to < 1e-10 (pure transform test, no dynamics).
   - **Ballistic-ΔV continuity:** at a converged patch, the SE-unstable and EM-stable states
     agree in INERTIAL position AND velocity to the manifold ε-tolerance.
   - **Energy bookkeeping:** each leg conserves its own CR3BP Jacobi to integrator tolerance;
     the patch ΔV (inertial velocity gap) is reported and bounded.

Golden source (Canalias C-values) is acquired by a separate #403-style task (Canalias 2007
PhD thesis, UPC, tdx.cat/handle/10803/5927; free PDF). The frame-transform + KLMR-convention
work can proceed in parallel; only the sourced-numeric validation tests gate on it.

Jacobi-convention caution (the #314 lesson): when transcribing Canalias C-values, verify our
`jacobi_constant` convention against theirs (the W-Z value carried a μ(1−μ)/2 term ours omits;
record any offset as `C_ours = C_pub − offset`).

---

## Architecture

One new module. No edits to reused code (`core/cr3bp`, `genome/heteroclinic_cycle`).

### File structure
- **CREATE `src/cyclerfinder/genome/cross_system_cycle.py`** — the whole Phase-A framework.
  Placement: `genome/` (consistent with #314 — it defines a new closure object + corrector).
- **CREATE `tests/genome/test_cross_system_cycle.py`** — two-tier validation tests.
- **APPEND** to `data/negative_results.yaml` if the closure search returns a clean negative
  (method-versioned, with re-sweep condition = "BCR4BP Phase B / wider grid").
- Golden (from the acquisition task): `data/golden/canalias_se_em_connection.yaml` (SE
  Lyapunov ICs + family C-values + any tabulated section data). Tests read EXPECTED from it.

### Reused signatures (no re-implementation)
- `core/cr3bp.CR3BPSystem(mu, primary, secondary, l_km, t_s)`, `propagate(...)`,
  `cr3bp_eom(t, s, mu)`, `jacobi_constant(s, mu)`.
- `genome/heteroclinic_cycle.LyapunovNode.from_libration(system, *, x0_guess, jacobi,
  period_guess, label, ydot0_sign=...)`, `_seed_on_manifold(...)`, `_section_crossing(...)`.

---

## Components

### 1. Inter-system frame transform (the crux new code)
- `FrameTransform` (or module functions `se_to_inertial`, `inertial_to_em`, etc.) mapping a
  6-state between the SE-rotating frame, a common Earth-centered inertial frame, and the
  EM-rotating frame, parameterized by the relative frame phase `theta` (the angle between the
  Sun-Earth line and the Earth-Moon line at the patch epoch) and the systems' length/time
  scalings. Earth is the SE secondary (position 1−μ_SE) and the EM primary (position −μ_EM);
  the transform handles this re-centering + the rotation-rate ratio.
- The relative angular rate is the EM synodic rate in the SE-rotating frame; `theta` is the
  free epoch parameter the connection corrector solves over.
- TEST: round-trip identity to < 1e-10 over random states and `theta` values.

### 2. Cross-system connection corrector
- `CrossConnection` dataclass: `(orbit_from_system, orbit_to_system, orbit_from, orbit_to,
  c_se, c_em, theta, tau_u, tau_s, k_u, k_s, branch_u, branch_s, patch_state_inertial,
  patch_dv_kms, residual, converged, n_iter, notes)`.
- `correct_cross_connection(se_system, em_system, orbit_from, orbit_to, *, section, ...)`:
  propagate the `orbit_from` unstable manifold (in its own system) to the patch Poincaré
  section expressed in the common frame; transform; propagate the `orbit_to` stable manifold
  (backward) to the same section. **Residual = inertial POSITION gap at the patch** (the
  corrector drives the two manifolds to the same point in the common frame); the inertial
  VELOCITY gap there is reported as the patch ΔV. A true (ballistic) heteroclinic connection
  has ΔV→0; the corrector's free variables `(tau_u, tau_s, theta)` close the position match,
  and ΔV→0 is sought by the closure search's energy freedom (a zero-ΔV connection exists only
  at compatible SE/EM energies — that compatibility is what the energy grid in §3 searches
  over). For the KLMR one-way reproduction, ΔV is small-but-nonzero (~tens of m/s) and is the
  validated quantity. Newton with FD-Jacobian + coarse-scan seeding (the #314 pattern). The
  exact free-variable/residual balance (position-on-section is 2 components vs 3 unknowns →
  the extra freedom drives ΔV down) is worked out in the plan.
- Never raises for "no connection" → `converged=False` + diagnostic note. `ValueError` only
  for malformed inputs (e.g. incompatible systems).

### 3. Closure search
- `search_cross_cycle(se_system, em_system, *, energy_grid, resonance_ratios, ...)`:
  for each (C_SE, C_EM, p:q) in a BOUNDED grid, attempt to chain an EM→SE connection and a
  SE→EM return connection into a loop that returns to the start orbit (periodic-up-to-rotation)
  with `theta` advancing commensurately (the ~19-yr Metonic / 235:19 closure from #316 is the
  natural first commensurability to test). Returns a list of `CrossCycle` results (closed or not).
- `CrossCycle` dataclass: `(connections, c_se, c_em, resonance, theta_closure_residual, closed,
  max_leg_residual, independent_residual, notes)`.
- A search that closes nothing over the grid → a registered clean negative (method-versioned).

### 4. Independent cross-check
- Reuse the #314 discipline: re-propagate any found cycle's legs with Radau (vs DOP853) and
  the frame transform, store max disagreement in `independent_residual`. Mandatory before any
  "closed" claim.

---

## Error handling & edge cases
- **Manifold misses the patch section** within the bounded horizon → connection `converged=False`
  with a diagnostic; never hangs (bounded integration, the #380/#314 discipline).
- **Non-commensurate theta closure** → cycle `closed=False` (the spatial connections can each
  close while the epoch fails to return — report it explicitly, don't fake closure).
- **No connection / no closed cycle** over the grid → clean negative registered in
  `negative_results.yaml` with method + re-sweep condition. This is an ACCEPTABLE outcome
  (novel search; per feedback_never_give_up_reproducing_papers, a clean negative on a NOVEL
  search is fine — this is novel, not a reproduction).
- **Convention mismatch** on Canalias C-values → surface + record the offset, don't silently
  fudge tolerances.

## Scope guard (YAGNI)
- Patched-CR3BP only. BCR4BP refinement = deferred Phase B (separate task).
- Bounded energy×resonance grid, SE/EM L1+L2 nodes only as needed — NOT a wide
  multi-libration campaign.
- No catalogue row admission in this phase: a found cycle would be a novel discovery requiring
  the full V0–V5 gauntlet + literature-novelty gate (`search/literature_check.py`) before any
  catalogue claim — that is downstream, not Phase A.

## Build sequence (for writing-plans)
1. Frame transform + round-trip identity test.
2. Per-system Lyapunov node reuse + Canalias SE C-value check (gated on the golden; write the
   test against the golden YAML, skip if not yet acquired — the #314 Task-8 pattern).
3. Cross-system connection corrector + KLMR one-way SE→EM reproduction (ΔV ~ tens of m/s) +
   the ballistic-ΔV-continuity and energy-bookkeeping identities.
4. Closure search over the bounded grid + commensurability handling.
5. Independent Radau cross-check; clean-negative registry path.
