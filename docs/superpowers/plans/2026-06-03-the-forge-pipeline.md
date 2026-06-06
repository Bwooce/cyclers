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
**First tasks:** (1) `solve_at_fidelity(cell, fidelity)` thin wrapper over the existing optimisers/ephemerides; (2) ~~test that S1L1 coplanar 4.9/5.0 → real‑eph shifts toward 5.65/3.05~~ **R1 (2026-06-05): SUPERSEDED — 5.65/3.05 is flagged `unverified-provenance` (untraceable, 2026-06-04) and S1L1 is multi-arc; that gate would violate golden discipline.** Anchor the ladder instead on (a) the **Aldrin outbound TOF shift**: coplanar 146 d vs Rogers 2012 STOUR analytic-ephemeris range **[161, 172] d** — both sides sourced and machine-readable since v4.2 (`segments[].tof_days_bounds` on `aldrin-classic-em-k1-outbound/out-em`), the cleanest documented fidelity shift we hold; and (b) the **Jones AAS 17-577 member rows** (sourced analytic-ephemeris V∞ multisets) as analytic-tier anchors once M-ED can solve VEM cells; (3) the persistence classifier.
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

## Revision R1 (2026-06-05 evening) — post-M8-Core/v4.3 refresh

State changes since this plan was written (binding where they contradict the
phase text above):

1. **Phase 1 anchor replaced** (edited inline above): the S1L1 5.65/3.05
   fidelity target is unverified-provenance; the ladder now anchors on the
   Aldrin TOF shift (146 d coplanar vs sourced STOUR bounds [161,172]) and
   the Jones analytic-ephemeris V∞ multisets.
2. **Phase 4's VEM frontier is now reachable.** The original Self-Review risk
   ("the value then shifts to VEM/3D/planet-centric scopes") has partly
   landed: M8-Core (2026-06-05, `933e75b`..`eb851a2`) gives the optimiser
   `Cell.period_basis` beat dispatch, the same-body Tisserand bypass, and
   `CONSTRUCTIBLE_MULTIBODY` admission — so `discover_novel` can enumerate
   VEM cells from day one. Sanity targets: the two Jones member rows
   (rediscovery), and AAS 17-577 **Table 1's 10-itinerary enumeration** as
   the VEM search-space spec (mining note 2026-06-05). The ballistic
   convergence xfail in `tests/test_vem_rediscovery.py` is the hand-off
   the Phase 4 loop should aim to flip.
3. **Novelty matching must respect v4.3 supersession links** (schema 4.3,
   2026-06-05): a candidate matching a row that carries `superseded_by`
   (e.g. the UNREALIZED `vem-emeeve-3syn`) must surface the supersession
   chain in the verdict, not report a clean "known" match against an
   invalidated premise. Add this to the Phase 3/4 match semantics.
4. **Phase 0 partial overlap — reconcile at execution.** Since the
   data-validation-hardening plan was written, the v4–v4.3 schema work
   shipped two-layer loader validation (its "loader schema validation" item)
   and some of its Task-2 invariants live in `validate_schema_invariants`.
   Re-scope Phase 0 against the live code before executing; in particular
   its `period.years ≈ k × synodic(pair)` invariant must handle **beat
   tokens** (`period.pair: "VEM-syn"` is not a body pair — M8 plan R1
   delta 3) or it will false-fail all four VEM rows.
5. **Powered screen (note only):** Sims-Flanagan machinery (#37, phases 2–5
   in flight) will provide a thrust-bounded maintenance evaluator — a
   candidate future gauntlet input (powered-feasibility screen), not a
   current dependency.

## Sourced algorithm candidates (added 2026-06-05, from the design-reference mining)

Enhancements with published provenance, slotted by phase. None are
prerequisites; each is a drop-in upgrade where named. Sources:
`docs/notes/2026-06-05-endgame-tisserand-mining.md` (Campagnola & Russell,
JGCD 33(2) 2010, DOIs `10.2514/1.44258`/`10.2514/1.44290`, read via AAS
09-224/09-227 preprints) and the JPL STAR tool survey (spec/prior-art notes).

1. **VILM ΔV-minimum quadrature** (Endgame Pt 1, Eq. 13: `ΔV_min = ∫ V∞/Γ dV∞`)
   — closed-form *lower bound* on any V∞-leveraging tour leg. Use in Phase 4
   as an admissible feasibility screen: prune candidate cells/sequences whose
   bound already exceeds the ΔV budget, before any optimiser run. Cheap,
   sourced, and conservative (never prunes a feasible candidate).
2. **Branch & bound / DP over the leveraging graph** (Endgame Pt 1) — their
   explicit search procedure; the DP variant enumerates the full ToF-vs-ΔV
   Pareto front. Candidate replacement for naive enumeration in Phase 4's
   sequence proposal once VILM nodes exist.
3. **Tisserand–Poincaré graph** (Endgame Pt 2) — CR3BP generalisation of
   `tisserand.py` (r_a–r_p axes; T = 3 − v∞² identity our
   `vinf_to_tisserand` already implements; patch-point equation = multibody
   twin of `linkable`; T > 3 unlocks ballistic intermoon transfers).
   The on-ramp for extending Phase 4 to planet-centric scope (#76).
4. **Dijkstra/A\* over the (body, V∞) Tisserand graph** — sequence proposal
   (recorded in spec §16.7.8 item 4; pairs naturally with candidate 1 as the
   heuristic).
5. **STAR-style triplet-ID search** (JPL STAR survey) — polynomial-complexity
   broad search via leg-triplet matching at shared bodies; candidate for
   Phase 4 if enumeration cost becomes the bottleneck.

**Candidate gate anchors for future modules** (golden-discipline compliant —
all literature-sourced, transcribed verbatim in the mining notes; none can
gate *current* modules):
- VILM module: Endgame Pt 1 Tables 1–3 (~25 Jovian/Saturnian intermoon
  transfers, ΔV_min/max + V̄∞ thresholds).
- T-P module: Endgame Pt 2 Table 1 (Europa/Titan insertion Δv max/min;
  two suspected source typos flagged in the note — do not anchor those cells).
- Earth-Moon compute (#76-adjacent): Hiraiwa et al. 2026 (arXiv:2602.17444)
  Table 3 literature values (Hohmann 3954 m/s; Sweetser min 3726 m/s,
  traceable to Topputo 2013); L1 Lyapunov C_J = 3.16.
- Low-thrust MGA (#37): Vasile & Campagnola JBIS Tables 3–4 — BLOCKED on
  re-transcription from a clean copy (digits read off rasters; see
  `docs/notes/2026-06-05-vasile-hiraiwa-scan.md`).

## Completion notes — Phases 4 + 5 (2026-06-06)

**Phases 4 and 5 are SHIPPED.** Modules:
`src/cyclerfinder/data/discover_novel.py` (the construction-first novelty loop +
the `BallisticClosureResult -> Cycler` bridge + supersession-aware verdict
routing), `src/cyclerfinder/verify/adversarial.py` (the N-verifier adversarial
panel: falsification probe, independent re-closure, perturbed-seed robustness;
majority-refute kills), `src/cyclerfinder/data/review_queue.py` (the NON-catalogue
human-review queue, `data/review_queue.jsonl`), and the orchestrator
`scripts/forge_novelty_run.py` (fan-out + per-finding panel + human gate).

**Workflow-tool choice:** no `.claude/workflows/` convention exists in this repo,
so the Python orchestrator alternative the plan sanctions was written
(`scripts/forge_novelty_run.py`), documented in its module docstring.

**Sanity-e2e (Phase 4 task 2):** the demonstrated DE440 E-M-E-E Sanchez-regime
closure flows end-to-end (bridge -> signature -> Axis-A agreement -> gauntlet)
and routes to SILVER (machine-confirmed, unsourced) — green
(`tests/data/test_discover_novel.py::test_sanchez_regime_closure_routes_to_silver_e2e`).

**First real run (Phase 4 task 3, 16-core, tight 2030 Sanchez window, 16 epochs x
4 E-M-E-E topologies):** 12 distinct closed families. **10 REJECTED**
(bend-infeasible — physically inadmissible, auto-falsified at the verdict).
**2 SILVER** (bend-feasible, panel-survived, queued): peak V_inf ~13.0 and ~12.1
km/s, both with `match=novel`. **The prototype's closed family does NOT match any
catalogue row** — its V_inf regime (~9.75/13.0 km/s) sits far above the sourced
Sanchez/S1L1 anchors (~3-6 km/s), the documented "S1L1 Mars-V_inf ~6.4 floor
generalised". These are the project's first machine-confirmed novelty
candidates, held at SILVER pending human review (never auto-promoted).

**VEM frontier (recorded reason, not scanned):** the #110 dense VEM scan produced
ZERO bend-feasible closures (floors 17.9/18.5 km/s); VEM single-ellipse-per-leg
novelty is empirically nil. `discover_novel` defaults to the E-M multi-arc space
where the corrector demonstrably closes bend-feasible chains; it does NOT burn
budget on VEM.

**Falsification e2e (Phase 5 task 3):** a fabricated candidate (perturbed V_inf,
impossible bend, non-converged) is majority-refuted by the panel AND
auto-REJECTED at the verdict — proven in
`tests/verify/test_adversarial.py::test_panel_refutes_fabricated_impossible_bend`
and `tests/data/test_discover_novel.py::test_bend_infeasible_closure_is_rejected_not_silver`.

**Censuses (Phase 5 task 4):** UNCHANGED — by design. Novelty findings never
create catalogue rows (golden discipline), so the catalogue
`EXPECTED_COVERAGE` / validation-tier / cycler_class ratchets are unaffected.
SILVER candidates live only in the non-catalogue review queue.

## Self‑Review
- **Independence enforced at every tier** (the spec's governing principle): GOLD needs an independent *source*; SILVER explicitly lacks one and is capped pending human review.
- **Builds on what exists** (~70%): enumeration, construction, both optimisers, 3D, lamberthub crosscheck, ledger/match. New work is mostly *gates and orchestration*, not new physics.
- **Golden‑clean:** physical invariants need no source; fidelity gate asserts *documented* shifts; falsification proves teeth; no tolerance loosening; novel discoveries never auto‑promoted.
- **Risk:** Phase 4 novelty yield is unknown — the loop may surface few/no novel ballistic E‑M cyclers (the classic space is well‑mined); the value then shifts to VEM/3D/planet‑centric scopes where the Forge's machinery transfers directly.
