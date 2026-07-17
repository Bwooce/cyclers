# #623 — Strategic review: open-item priorities + novel capability combinations (2026-07-17)

Analysis-only pass (no code, no catalogue writes, no dispatches), mirroring `#605`'s format.
Inputs: the reconciled `## CURRENT STATE` dashboard, the full `#605`-`#622` ledger range and each
task's dedicated bullet, the eight new capability modules built across that arc, the
`empty_regions.jsonl` registry (82 entries), and `docs/notes/CORPUS_INDEX.md`.

Context that frames everything below: `#605` diagnosed the program's recurring failure mode as the
**family-selection/basin wall**, and the `#606`-`#620` arc answered it with a genuinely new
capability class — seedless, integration-free spectral correctors (CR3BP periodic `#606`, QBCP
periodic `#611`, CR3BP torus `#612`, QBCP torus `#617`/`#618`, open QBCP arc `#620`) — plus a
generative seed model (`#608`), an interval-certification POC (`#610`), and family-mapping tools
(`#613`/`#615`). That arc's discovery-side outcomes were almost all **well-characterized
negatives** (SE↔EM twice-independently blocked; Phobos-Deimos, small-body multi-moon, Uranian
3-moon, EM-mu niching all clean negatives). The capability inventory is now much richer than the
set of targets it has been pointed at; this review looks for the highest-value untried pairings.

---

## Part A — Priority read on the CURRENT STATE open items

### #542 (learned-seed warm-start; "unfamiliar basin" claim untested) — worth ONE more bounded pilot, then close either way

`#608` answered the feasibility question (~12.25x physically-sane convergence lift over blind
seeding, reproduced bit-for-bit) and `#614` then killed both of `#608`'s own flagged improvement
paths (family-tagging, nonlinear encoder) with a clean statistical negative at this corpus size.
What remains open is exactly one claim: can a model trained on this corpus propose useful seeds in
a *genuinely unfamiliar* basin? Everything tested so far is in-distribution Earth-Moon CR3BP. The
honest read: this is no longer a capability-building question — it is a single falsifiable
transfer-generalization question, and it has a cheap, decisive test (see B1 below, the top
shortlist item). If that pilot shows the lift collapses off-distribution, `#542` should be closed
as answered-negative rather than kept half-open indefinitely; if the lift survives, `#542` becomes
a real discovery lever worth productionizing. Either outcome ends the limbo. Do not invest in
bigger/deeper models first — `#614` already showed model sophistication is not the binding
constraint at 54k examples.

### #520 (aborted 8,640-point 3D sweep) — do NOT revive; supersede it once #622 lands

Disagreeing with the standing "would need a timed pilot + incremental checkpointing first" framing:
that redesign advice is process hygiene, not a reason the sweep should exist. `#520`'s unique
content was coverage of the 4-libration-pair cross-system matrix; since then `#516`/`#517` swept
2 of the 3 untested pairs (both empty, registered), `#411`'s pair is a fully characterized 1-DOF
phase-closure negative, and `#622` (running now) covers the last cell (EM-L1↔SE-L1) with exactly
the narrow, checkpointed format `#520` lacked. What would remain of `#520` after `#622` is only
"the same 4 pairs at a denser grid" — precisely the "NOT more revolutions / finer grids" anti-
pattern `#411`'s own entry warns against — and the higher-fidelity QBCP version of the same
physical question has meanwhile ended in a two-method-independent negative (`#619` manifold
conditioning, `#620` ghost minima). Recommendation: when `#622` closes, mark `#520`
CLOSED/SUPERSEDED (by `#411`+`#515`-`#517`+`#622`) in the dashboard and ledger, and never rebuild
it. The checkpointing lesson is already memorialized (DELTA F2, incremental-progress-reports
memory) and doesn't need a vehicle.

### Ledger-staleness housekeeping (#582/#583/#585/#589 one-liners) — done in this pass

This was cheap (~4 in-place edits) and the failure mode is documented as real, not theoretical: a
stale `#516`/`#517` one-liner propagated into a same-day planning recommendation before being
caught, and the dashboard's own `#557` entry once reproduced a stale header instead of catching it.
All four named one-liners are now corrected in place in the `TASK ALLOCATIONS` paragraph (marked
"one-liner corrected 2026-07-17"), and the dashboard's "Ambiguous" section updated to record that.
Residual honest caveat: no claim these four were the only stale lines — the "rough index only"
caution on the rest of the ledger paragraph stands. A full hand-audit of every one-liner is NOT
recommended as its own task (the per-bullet entries + dashboard are the status sources; the
marginal value of polishing ~180 one-liners is low).

---

## Part B — Ranked shortlist: new combinations and research

### B1. Cross-μ transfer pilot for #608's generative model — the direct test of #542's open claim
**What:** Take `#608`'s already-built model (trained on the Earth-Moon μ=0.01215 corpus,
`src/cyclerfinder/ml/orbit_generative.py`) and evaluate its generate-then-refine lift at a μ it has
never seen — e.g. μ=0.001 (the Ross-RT golden anchor, so positive controls exist) and Sun-Earth
μ≈3.0e-6 — refining with the same `cr3bp_periodic.correct_periodic` at the new μ and comparing
converged-AND-physically-sane rate vs. the uniform-bounding-box baseline at that μ, exactly
`#608`'s own evaluation protocol.
**Why plausible:** the genome (state0, period, Jacobi) is nondimensional and family geometry varies
smoothly in μ; a density model that captures family-shaped structure rather than memorized points
could retain part of its lift. Why it might fail: `#614` showed the model lands in "gaps between
families" — off-distribution μ may amplify exactly that.
**Cheapest pilot:** a few hours; ~100 samples per μ, reuse `run_608`'s harness with a μ parameter.
No new model training.
**Confidence:** MEDIUM that the pilot is decisive (either direction), LOW-MEDIUM that the lift
transfers. High value-per-cost because it converts `#542` from "perpetually not-closed" to a
closed question either way.

### B2. Generalize #610's interval bend-gate certificate to every bend-gate-limited registry entry
**What:** `#610` certified `sup(bend) < 5°` for the Proteus sub-gate via rigorous interval
arithmetic (`scripts/certify_610_proteus_bend_interval.py`, mpmath.iv). At least six other
`empty_regions.jsonl` negatives share the *identical* closed-form failure mode — the `#324`
physical max-bend gate with an undersized moon GM: `#607` (Sylvia/Eugenia/Kleopatra/Elektra, all
four "mass-limited, same failure mode"), `#609` (Phobos+Deimos, bend tops out 0.0159° vs 5°),
`#571`'s two "analytically-empty" Saturn entries, and the Amalthea/Triton precedents. Each entry
already records the survivor V∞ ranges needed to define the certification boxes.
**Why plausible:** pure parameter substitution into an existing, tested, reviewed script — the
Bate-Mueller-White formula is the same; only GM, r_p floor, and the V∞ box change. This is
mechanical (Sonnet/Haiku-tier), not research.
**Cheapest pilot:** generalize the script into a small parameterized helper + one certification run
per entry + a registry field recording "bend sub-gate certified over box [..]". A day of work total.
**Confidence:** HIGH on feasibility. Value is real but bounded: it upgrades the program's main
product (clean negatives) from grid-conditional to continuum-strength on the bend axis, exactly the
`[[project_negative_results_registry]]` ambition — but it discovers nothing new. Honest scope
limit, inherited from `#610`: the Lambert/residual half of those searches stays uncertified
(bounding a multi-rev branch-selecting Lambert solve rigorously is still unsolved).

### B3. Seed #620's arc collocation from #619's best ballistic near-miss — the arc's own last untried refinement
**What:** `#620`'s verdict was precise: collocation can VERIFY a near-real seed but cannot SEARCH
without one, and its bullet explicitly names the one untried refinement — reconstruct `#619`'s
best ~166,000 km / ~400 m/s ballistic near-miss (norm-0.855 local minimum, Radau-confirmed real)
via the manifold machinery, and feed it to `correct_qbcp_arc_connection` as the seed, testing
whether global collocation can refine a shooting near-miss into a real connection. Optionally also
try `#620`'s other flagged variant (over-collocation, m >> order, to suppress ghosts).
**Why plausible (and why probably not):** it is the only seed in existence that is anywhere near
the real basin, and refinement is exactly the regime `#620` proved the machinery handles (both
positive controls). Against it: two independent methods floored at O(1e5 km), and 166,000 km is
~40% of an Earth-Moon distance — "near-miss" is generous; the ghost-minima structure may simply
recapture the solve.
**Cheapest pilot:** one focused session — rebuild `#619`'s best solution (its inputs are fully
recorded), convert to collocation node values, run at orders 40/60 with the existing
`independent_closure_check` as sole arbiter.
**Confidence:** LOW that it closes (~10-15%), but the payoff asymmetry justifies it: a closure
would be the first genuine cross-system SE↔EM cislunar cycler (catalogue-grade, after the mandated
Fable review); a failure completes the `#538`-`#620` arc's negative with its last named alternative
exhausted, making the negative maximally citable. Either way the arc stops carrying an "untried
refinement" asterisk. NOTE: this deliberately reuses the manifold machinery `#620` avoided — that
is fine; the point of `#620`'s avoidance (don't depend on the direction during SEARCH) doesn't
apply to seeding a refinement.

### B4. μ-continue the Ross-RT (k1,k2) ballistic-cycler families DOWN to planet-moon μ values
**What:** `#494` extended the Ross-Roberts-Tsoukkas 2026 stable ballistic prograde (k1,k2)-cycler
families UP in μ (binaries: 0.1085 Pluto-Charon through 0.5), closing `#315`. Nobody has continued
them DOWN below the paper's own μ=0.001 floor to real planet-moon systems: Saturn-Titan (2.4e-4),
Neptune-Triton (2.1e-4), Jupiter-Ganymede (7.8e-5), Jupiter-Europa (2.5e-5), Uranus-Titania
(~4e-5). The solver is μ-agnostic (`#494`'s own reframing) and the Table-I golden + Phase-0
positive control already exist.
**Why plausible:** cheap continuation with proven machinery; stable ballistic single-moon cyclers
at named outer-planet moons would be catalogue-class quasi_cycler/cycler rows, and the published
record (Russell-Strange 2009, Liang 2024, RRT's own EM-focused papers) covers flyby-driven
multi-moon cyclers, not this stable-ballistic class at these systems.
**Honest caveat (the reason this isn't ranked higher):** as μ→0, symmetric stable prograde
periodic orbits shade into classical Poincaré first-kind/resonant families — existence at small μ
is close to theorem-guaranteed, so the *finding* risks being "a known classical object, relabeled."
The novelty gate (`search/literature_check.py`, mandatory per
`[[feedback_literature_novelty_check_baseline]]`) and the stability-index/geometry specifics
(does the perimoon pass give useful encounter distances?) are where any real contribution would
live. Expect census, not a novel species — which is still the program's stated post-pivot job.
**Cheapest pilot:** continue the two EM representatives (k=(1,1), (3,3)) down to Titan μ and check
(a) family survives with |s|<1 stability, (b) perimoon geometry is encounter-relevant, before any
full sweep.
**Confidence:** HIGH on feasibility, LOW-MEDIUM on catalogue-worthy novelty.

### B5. Literature/methods: nothing genuinely new found — one standing reserve, honestly parked
A real pass over `CORPUS_INDEX.md` (354 lines, all sections) looking for un-adopted METHODS (not
IC tables — `#621` closed that angle) found no new actionable method that the `#605`-`#620` arc
hasn't already absorbed or that earlier triages haven't already, correctly, rejected
(Ozaki-2022 surrogate: its screening role is doubly foreclosed by `#317`'s measured
gate-efficiency + weak-signal results; the lobe-dynamics and Keplerian-map items were mined and
their applicable pieces extracted long ago; `#604` killed the map-chaining idea on V∞
incompatibility). The one genuine reserve: **Wilczak-Zgliczyński's computer-assisted-proof
machinery** (both Parts in corpus, digested, and this codebase already reproduces their Oterma
heteroclinic golden to ~1e-10) is the only credible route to *theorem-grade dynamical* claims —
the natural growth path if `#610`-style certification is ever to cover propagation-dependent
negatives (e.g., the SE↔EM non-connection) rather than closed-form sub-gates. That is a
multi-week, research-grade build (rigorous integration/covering relations; note non-existence
over a region is a harder shape than W-Z's existence proofs) with no current forcing function —
park it, but record that the preconditions here (golden reproduction + digests in hand) are
unusually good if a high-stakes claim ever needs it. Minor: the Litteri et al. CMDA 138:25 journal
version (`#608`'s external anchor) is still paywalled/unacquired — a small acquisition item, not a
method gap.

### Explicitly rejected combinations (checked, don't re-propose)
- **#608 seeds → #620 arc search (the obvious pairing):** does not hold up. `#608`'s model is
  trained on CR3BP periodic-orbit *genomes*; `#620` needs a near-real *QBCP transfer-arc
  trajectory* seed — an object class with zero training examples anywhere (none exist; that
  absence is the finding). The ghost-minima wall is basin-structure, not seed-polish; only a
  near-real arc helps, and B3 uses the only one there is.
- **QBCP-fidelity re-pass of the Galilean/Uranian/Saturnian sweeps:** does not hold up, twice
  over. (1) `core/qbcp.py` is hardwired Sun-Earth-Moon (Gimeno-Jorba Fourier-fit alphas); a
  per-system quasi-bicircular refit (Andreu's method) is a real multi-week build per system.
  (2) More fundamentally, those systems' negatives/closures are not fidelity-limited: `#607`/
  `#609`/`#599` fail on moon GM (physics a better model cannot change), `#575` on inclination-
  extension geometry, and the Galilean positives are novelty-blocked (published R-S/Liang
  territory), while `#619`'s conditioning finding predicts the strongly-unstable libration-point
  route hits the same wall there. No target in those systems is waiting on model fidelity.
- **Bigger/fancier generative models at current corpus size:** `#614` already ran this exact
  comparison (family-tags, nonlinear autoencoder) — clean negative at n=54k; don't re-litigate.

---

## Recommended dispatch order (for the coordinating session)
1. **B1** (cross-μ `#608` pilot) — hours, decisive for `#542` either way. Sonnet.
2. **B2** (interval bend-gate generalization) — ~a day, mechanical, registry-strengthening. Sonnet.
3. **B3** (`#619`-near-miss-seeded `#620` refinement) — one session, low odds/high asymmetry. Opus
   (closure adjudication is trust-bearing).
4. **B4** (RRT small-μ continuation pilot at Titan) — cheap pilot first; full sweep only if the
   pilot passes both its gates. Sonnet for the continuation, Opus/Fable for any novelty verdict.
5. Housekeeping already done in this pass (ledger one-liners); `#520` supersession waits on `#622`.
