# The Forge — Discovery + Cross‑Validation Pipeline (Phased Plan)

> **For agentic workers:** This is a milestone roadmap (like the `docs/phases/mN` plans). Each phase is executed as its own superpowers:subagent‑driven‑development run with a detailed task breakdown written at execution time. Work on `main` — do NOT branch. Design spec: `docs/superpowers/specs/2026-06-03-the-forge-discovery-pipeline-design.md`.

**Goal:** Implement the Forge — enumerate candidate cyclers, solve construction‑first, run each through a 4‑axis cross‑validation gauntlet, and emit a tiered‑confidence ledger with novel discoveries held at "machine‑confirmed, unsourced" until a human gate.

**Dependency:** Phase 0 = the **data‑validation‑hardening plan** (`2026-06-03-data-validation-hardening.md`) — it builds Axis C (provenance/fidelity/corroboration) and Axis D (falsification/adversarial), which the gauntlet consumes. Do that first.

**Sequencing rule:** each phase ends green (full suite + ruff + mypy), commits locally per task, golden‑discipline throughout. Phases are independently valuable — ship and stop at any boundary.

---

## Phase 0 — Foundation (the data‑validation‑hardening plan)
Provenance + fidelity tags, tiered validation gate, physical invariants, corroboration scoring, falsification + lamberthub guards, loader schema validation. **Prereq for Axes C/D.**

## Phase 1 — Axis B: fidelity‑ladder gate
**What:** a `fidelity_persistence(cell, anchors)` check that solves a candidate at coplanar → analytic‑ephemeris → real‑DE440 and asserts V∞ either persists within tolerance or shifts in the *documented* direction (e.g. Mars V∞ drops with eccentric Mars).
**First tasks:** (1) `solve_at_fidelity(cell, fidelity)` thin wrapper over the existing optimisers/ephemerides; (2) test that S1L1 coplanar 4.9/5.0 → real‑eph shifts toward 5.65/3.05 (the documented, *predictable* fidelity shift — turning the old "bug" into an asserted gate); (3) the persistence classifier.
**Unlocks:** systematic prevention of the cross‑fidelity confusion class.

## Phase 2 — Axis A: code‑path agreement
**What:** for a solved cycler, require ≥2 independent code paths to agree on V∞/geometry: (a) in‑house Lambert vs lamberthub izzo+gooding (have it), (b) resonance‑construction vs free‑optimiser, (c) forward Kepler re‑propagation of the built cycler (not the Lambert that made it).
**First tasks:** (1) `crosscheck_code_paths(cycler) -> AgreementReport`; (2) test that construction and optimiser agree on Aldrin/S1L1 within tol; (3) wire the Kepler re‑propagation residual as a gate.

## Phase 3 — Gauntlet orchestrator + ledger tier
**What:** `run_gauntlet(candidate) -> ValidationVerdict{tier, axis_results, confidence, provenance}` combining Axes A–D into the GOLD/SILVER/BRONZE/REJECTED tiers; write verdict + audit trail to the ledger; add the **validation‑tier census ratchet**.
**First tasks:** (1) the verdict dataclass + combiner; (2) ledger schema extension for verdicts; (3) ratchet test frozen at current counts.

## Phase 4 — Novelty discovery loop
**What:** enumerate cells NOT matching any catalogue signature, solve construction‑first, run the gauntlet; route GOLD/SILVER to the human gate. Loop‑until‑dry on the deepening frontier.
**First tasks:** (1) `discover_novel(bodies, k_max, …)` over `feasible_cells` filtered to non‑matching signatures; (2) novelty classifier vs `canonical_signature`/`match`; (3) a bounded end‑to‑end test that finds the *known* 2‑syn E‑M as a "rediscovery" (sanity) and reports any novel candidate as SILVER‑pending‑human.

## Phase 5 — Adversarial orchestration (workflow) + human gate
**What:** a Workflow that fans candidates out, runs the gauntlet in parallel, applies **per‑finding adversarial panels** (N independent verifiers; majority‑refute kills), loops‑until‑dry, and emits a human‑review queue for novel‑confirmed candidates. Optional stretch: a GMAT high‑fidelity sign‑off hook (V4).
**First tasks:** (1) the workflow script (fan‑out + adversarial‑panel pattern); (2) the human‑review queue artifact; (3) a falsification end‑to‑end (inject a bogus candidate, assert it's REJECTED).

---

## Self‑Review
- **Independence enforced at every tier** (the spec's governing principle): GOLD needs an independent *source*; SILVER explicitly lacks one and is capped pending human review.
- **Builds on what exists** (~70%): enumeration, construction, both optimisers, 3D, lamberthub crosscheck, ledger/match. New work is mostly *gates and orchestration*, not new physics.
- **Golden‑clean:** physical invariants need no source; fidelity gate asserts *documented* shifts; falsification proves teeth; no tolerance loosening; novel discoveries never auto‑promoted.
- **Risk:** Phase 4 novelty yield is unknown — the loop may surface few/no novel ballistic E‑M cyclers (the classic space is well‑mined); the value then shifts to VEM/3D/planet‑centric scopes where the Forge's machinery transfers directly.
