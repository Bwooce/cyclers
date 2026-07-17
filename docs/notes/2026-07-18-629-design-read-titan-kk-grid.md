# #629 design read — real-moon-μ Ross-RT (k1,k2) cyclers: 2D grid vs homotopy (Fable, 2026-07-18)

Mandatory pre-dispatch design gate for task #629 (per the #629 bullet in `data/OUTSTANDING.md`:
"get a Fable/Opus design read on whether a 2D grid or a homotopy/continuation-in-a-different-
parameter approach is the better fit, and a rough cost estimate, before committing"). Analysis
only — no production code, no catalogue edits. All numbers below were computed with the repo's own
`pluto_charon_kk_sweep._c_l1` / `cr3bp_system("Saturn","Titan")` (run 2026-07-18, session log),
not copied from prose.

## Headline finding: #627's negative is under specific, quantified suspicion of being a
## parameter-scaling artifact — the design question changes shape

Before choosing between "2D grid" and "homotopy," the design read found something that reframes
both: **#627's two failure modes are each quantitatively consistent with a single cause — the
C-tracking walk's absolute `c_margin` does not scale with the corridor it is trying to track,
which shrinks as μ^(2/3).** This does not prove the family survives to Titan μ (a real bifurcation
could still kill it), but it does mean #627's pilot never actually searched the corridor where the
family should live, so "the admissible region is vanishingly small" is NOT established.

### The scaling evidence

The natural energy scale near the L1 neck is `C_L1(μ) − 3 ≈ 3^(4/3) μ^(2/3)` (Hill scaling).
Computed with `_c_l1`:

| μ | C_L1 | C_L1 − 3 |
|---|------|----------|
| 0.012151 (EM anchor) | 3.188341 | 0.188341 |
| 0.001 (RRT floor)    | 3.039949 | 0.039949 |
| 2.36695e-4 (Titan)   | 3.015770 | 0.015770 |

The two sourced (1,1) Table-I anchors sit at a nearly **invariant scaled energy**
`ρ ≡ (C − 3)/(C_L1(μ) − 3)`:

| anchor | C | ρ = (C−3)/(C_L1−3) | C_L1 − C |
|--------|---|--------------------|----------|
| (1,1) @ μ=0.012151 | 3.151176 | **0.8027** | 0.037165 |
| (1,1) @ μ=0.001    | 3.031606 | **0.7912** | 0.008343 |
| (3,3) @ μ=0.012151 | 3.183379 | **0.9737** | 0.004962 |
| (3,2) @ μ=0.1      | 3.573368 | 0.9605     | 0.023586 |

ρ is nearly constant along the (1,1) family across a 12x μ range (0.803 → 0.791). Extrapolating,
the (1,1) family at Titan μ should live near **C ≈ 3 + 0.80×0.015770 = 3.0126** (corridor roughly
C ∈ [3.009, 3.015]); the (3,3) family near **C ≈ 3.0154** (an extremely thin corridor,
C_L1 − C ≈ 4e-4, since its ρ ≈ 0.974).

### What #627's walk actually did, in these coordinates

- **(1,1)**: `c_margin = 0.02` clamps C ≤ C_L1(μ) − 0.02 at every step. Already at the μ=0.001
  anchor itself, C_L1 − C_anchor = 0.00834 < 0.02 — so the walk forced C *below the anchor's own
  natural energy from the very first step*, and by Titan μ had forced C = 2.995770, i.e.
  **ρ = −0.27** — below C = 3 entirely (all zero-velocity necks open, outside the co-orbital/
  ballistic-capture band ρ ∈ (0,1) where every known member of this family lives). Landing on a
  Saturn-only (1,0) branch there is the *expected* outcome, not evidence about the family.
- **(3,3)**: its anchor sits only 0.00496 below C_L1(0.012151); the tuned `c_margin = 0.005`
  equals that entire gap, and as μ falls the clamp drags ρ from 0.974 toward 0.68 — ripping the
  walk off a corridor whose half-width is ~1e-3 in ρ-units. Outright convergence failure partway
  down is again the expected outcome.
- **The 81-point follow-up C-sweep did span the corridor** (C ∈ [2.9558, 3.0358] ⊃ [3.009,
  3.0158]) **but on the wrong branch**: `scan_c_family_at_mu` continues the corrector from the
  landed (1,0) member's x0 = −0.3397 at hc=1 with a period-jump filter — it traced the Saturn-only
  branch through the corridor's C-range, it did not search fresh x0 there. So "0/81 on-target" is
  conditional on the (1,0) seed branch, exactly per
  [[feedback_verify_gauntlet_with_positive_control]]'s "no X found is conditional on search
  formulation."

Also relevant: the fold at the μ=0.012151 (1,1) anchor (#627's dμ=1e-7 branch jump) is a
*separate, real* finding and is already worked around (start from the μ=0.001 anchor); nothing
below re-litigates it.

## Q1 — Is a 2D grid in (x0, C) at Titan μ well-posed? YES, but only in scaled coordinates

Concretely: the residual `r(x0; C, μ) = xdot(t_half) = 0` at fixed (C, μ, hc) has isolated roots
in x0, so the family is a curve in the (x0, C) plane — a grid + corrector-basin finds curves fine
(this is exactly what `pluto_charon_kk_sweep._grid_seed_search` / `sweep_21`/`sweep_22` already do
at Pluto-Charon; the machinery EXISTS, it has just never been pointed at a small-μ target with
correctly scaled bounds).

- **C bounds**: NOT a naive [3.0, C_L1] span at EM-style dc — the target lives at
  ρ ∈ ~[0.6, 0.99], i.e. C ∈ [3.0095, 3.0156] at Titan, a band of width ~6e-3. A naive
  EM-carried-over grid (dc = 0.01) would put 0–1 points in it. Grid in ρ (or equivalently
  C = 3 + ρ(C_L1−3)) with dρ ≈ 0.02 (~20–25 C values), plus a finer dρ ≈ 0.005 sub-band near
  ρ ∈ [0.95, 0.995] for (3,3)-like high-ρ families.
- **x0 bounds**: the weakest-constrained axis. The (1,1) anchors drift −0.768 → −0.647 as μ falls
  0.0121 → 0.001 with no clean two-point scaling exponent (Hill-unit distance from the moon's
  orbit radius *grows*: 1.0 → 3.5 μ^(1/3) units), so use a generous span x0 ∈ [−0.95, −0.30] at
  Δx0 ≈ 0.005 (~130 points) rather than trusting an extrapolation. hc is a third, enumerated axis
  (hc ∈ {1, 3, 5, 7, 9}, per the anchors' hc=1/7 and #504's hc=6 precedent — include even values
  too if budget allows).
- **Is the admissible region vanishingly small?** In absolute C-units it shrinks as μ^(2/3) — but
  so does the natural coordinate; in ρ-units the (1,1) corridor is a fat target (Δρ ≈ 0.2–0.4 wide
  band around 0.8, if the family survives). The genuine thin case is (3,3) (Δρ ~ 0.01–0.03 near
  0.974); flag that its grid needs the fine high-ρ sub-band or it will be missed. #627's
  "shrinks below C_L1" finding is evidence about the *walk's clamping*, not about corridor
  measure.
- **What the grid cannot settle by construction**: whether a *stable* (|ν|<1) window exists on any
  found branch — that's the existing `c_sweep_find_nu_zero` refinement per found branch, exactly
  as `sweep_family_grid` already chains it.

## Q2 — Homotopy in a different parameter?

- **Continuation in C at fixed Titan μ from a known near-C_L1(Titan) member**: no such on-topology
  member exists (that absence is the problem being solved), so this is not available as a
  standalone approach — it is the *refinement stage* after a grid seed (already how
  `sweep_family_grid` works).
- **Hybrid (μ-walk to an intermediate μ, then fixed-μ C-homotopy)**: partially relocates the
  problem — any μ-walk must thread the pinching corridor for the *entire* final stretch, whereas a
  fixed-target-μ search only needs the corridor nonempty *at* Titan μ. But the walk done in the
  right coordinates is cheap and worth one honest retry first:
- **The scaled-margin rerun (recommended Phase A)**: re-run #627's own
  `mu_step_to_system_tracking_c_l1` walks with `c_margin` proportional to the corridor —
  `c_margin = α (C_L1(μ_next) − 3)`, α ≈ 0.1–0.2 (equivalently: track constant ρ) — instead of an
  absolute constant. This is a ~one-parameter change to an existing function plus a rerun
  (~10–15 min compute total, per #627's own 449 s + 166 s timings). It directly addresses BOTH of
  #627's specific failure modes as diagnosed above, rather than generic continuation advice: the
  (1,1) topology loss (walk forced to ρ<0) and the (3,3) convergence loss (margin ≥ its whole
  corridor). The fold at the Rep-2 anchor stays sidestepped by the μ=0.001 waypoint. If a real
  bifurcation kills the family anyway, this rerun fails *in the corridor*, which is then a
  well-characterized negative instead of a suspect one. `continue_in_mu`'s pseudo-arclength (which
  can turn folds, unused by #627's pilot) is a fallback if the natural-parameter walk hits a
  genuine fold mid-range — not needed up front.

## Q3 — Theory-first instead of numerics? Partially — and the cheap part is done above

- The corpus has **no** small-μ asymptotics asset for this family: no Hénon *Generating Families*,
  no Poincaré first/second-species treatment, no co-orbital/quasi-satellite RTBP paper
  (CORPUS_INDEX grep; nearest neighbors are Belbruno 2004 (ballistic capture, adjacent) and the
  Antoniadou-Libert planetary-MMR continuation papers (different model regime, μ=0.001 planetary)).
- The classically available cheap prediction is the **Hill-limit scaling** used above, and it has
  now been extracted from the project's own two sourced anchors: ρ ≈ const along the family. It
  predicts (a) *where* to look (C ≈ 3.0126 for (1,1); C ≈ 3.0154 for (3,3)); (b) *existence
  plausibility*: μ = 2.4e-4 is not qualitatively different from the paper's own 1e-3 floor in
  Hill-scaled units (the limit is μ-independent), so persistence is expected rather than
  surprising. What theory does NOT give: the *stability* window (|ν|<1) at Titan μ — Poincaré
  first-kind theory covers secondary-avoiding orbits; k2≥1 moon-loop members at C>3 are
  co-orbital/ballistic-capture-class objects whose stable windows are known only numerically even
  classically. So the search still does the deciding work — but targeted, not blind. Acquiring the
  classical corpus (Hénon vol. 1/2, Benest, quasi-satellite literature) is optional context, NOT a
  gate; the mandatory `literature_check.py` pass on any hit MUST however cover the co-orbital/
  quasi-satellite/ballistic-capture families explicitly (Hénon family f, Benest, Sidorenko/
  Pousse-Robutel-Vienne QS theory, Titan ballistic-capture and Titan QS mission literature,
  Russell-Strange 2009) — that is where a "known classical object, relabeled" verdict would come
  from.

## Q4 — Cost estimate (honest, order-of-magnitude)

- **Phase A — scaled-margin rerun of #627's own walks**: ~1 h total (tiny param change to
  `mu_step_to_system_tracking_c_l1` call/signature + rerun both representatives, ~15 min compute).
  Decisive about the artifact question either way. Sonnet.
- **Phase B — fixed-μ corridor grid at Titan** (if A doesn't already land on-topology members, or
  to map the corridor even if it does): script reusing `_grid_seed_search`-pattern + existing
  gates (`winding_topology`, `barden_stability`, `c_sweep_find_nu_zero`, `find_perimoon_passage`,
  `flyby_is_useful`, Radau crosscheck, `preflight_search`). ~130 x0 × ~25 C × ~5 hc ≈ 16k
  corrector calls, SIGALRM-bounded 4 s worst-case, ~0.2–0.5 s typical → **~2–8 h single-core;
  overnight with margin**. Build ~half a day (Sonnet) since every component exists. Remember
  `tests/scripts` in the ratchet run ([[feedback_verify_scope_must_include_tests_scripts]]).
- **Phase C — adjudication + mandatory literature check on any hit**: Opus/Fable, ~1 session.
  Writeback held until confirmed per [[feedback_orbit_closure_discipline]].
- **NOT recommended now**: multi-system sweep (Triton/Ganymede/Europa/Titania — defer until the
  Titan corridor answer exists); building analytic-seed machinery; a Hill-problem solver; corpus
  acquisition as a precondition.

## Q5 — Worth building at all right now? GO-with-caveats, for two separable reasons

1. **Registry integrity (the stronger reason, novelty-independent)**: #627's clean negative — a
   registered method-versioned result — is now under specific quantified suspicion of being a
   false negative caused by an un-scaled parameter. The project's own discipline
   ([[feedback_bugfix_invalidates_past_searches]], in spirit: mis-parameterized solvers are
   false-negative generators) makes the ~1-hour Phase A rerun near-obligatory regardless of
   discovery value. Either outcome improves the record: rescue → the negative was an artifact;
   corridor-failure → the negative upgrades from "suspect" to "well-characterized."
2. **Discovery value (honest read)**: even a full success is most likely a "known classical
   object, relabeled" (co-orbital/QS/ballistic-capture territory at Titan is well-plowed) — the
   #629 bullet's own framing stands, do not re-inflate. But #627's live literature check found the
   RRT paper's own claimed range stops at Sun-Jupiter μ (~9.5e-4) and no prior work at real-moon
   μ, so a confirmed stable, encounter-relevant (k1,k2) ballistic Titan member would be a
   legitimate catalogue-grade census row (probable quasi_cycler class) even if the *species* is
   classical. At ~1 Sonnet session + overnight compute, expected value clears the bar. Anything
   bigger does not, yet.

**Recommendation: GO, phased A→B→C as above; NO-GO on the multi-system sweep and on any
new-machinery build until the Titan corridor answer is in.** Recommended models: Sonnet for
Phases A/B (deterministic gates catch its mistakes), Opus/Fable for Phase C verdicts.
