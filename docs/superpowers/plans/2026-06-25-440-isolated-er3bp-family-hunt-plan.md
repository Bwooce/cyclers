# #440 — Direct hunt for ISOLATED ER3BP families (Antoniadou & Libert 2018)

> Status 2026-06-25: SCOPED + GATED **GO**. Two spikes passed (model-match
> scoping + CRTBP 3/2 golden gate). This plan records the validated path. Build
> is multi-day; gate each phase.

**Goal:** find a genuine *e>0-only* periodic family — one with NO circular limit —
which Antoniadou & Libert 2018 (CMDA, DOI 10.1007/s10569-018-9834-8) prove exist
("isolated, not continued from bifurcation points," high-e, stable). Our
#432/#435/#436 continuation-from-circular methods are structurally blind to these
(see `2026-06-25-er3bp-negatives-erratum-antoniadou-libert.md`); this is the one
untested mechanism that could actually yield a novel family.

## De-risking already done (do NOT redo)

1. **Model-match (scoping spike):** Antoniadou's planetary-MMR ER3BP is the SAME
   dynamical system as our `er3bp.py`, different coordinates — their giant
   eccentricity e₂ = our pulsating-frame `e`. **No new EOM needed.** Their
   isolated families have **NO tabulated state-vector ICs** (purely graphical:
   `(x,e₁)`, `(e₁,e₂)`, DS-maps), so seeds must be digitized from figures or
   re-derived via their Poincaré + continuation.
2. **CRTBP 3/2 golden gate: GO.** `correct_symmetric_fixed_jacobi` converges the
   interior 3/2 resonant PO at µ=0.001 to machine precision. **Validated recipe:**
   - seed `x0 = -mu + a1`, `ydot0 = sqrt((1-mu)/a1) - a1` (inertial circular −
     synodic rotation), `jacobi = jacobi_constant(seed, mu)`.
   - call `correct_symmetric_fixed_jacobi(sys, x0, jacobi, period_guess=4*pi,
     ydot0_sign=+1, half_crossings=1, tol=1e-11)`.
   - converged member: a≈0.7596 (on 3:2, n=1.51), T≈11.93, closure 1.6e-15.
   - **GOTCHA:** free-period `correct_periodic` is a FAMILY-SELECTION TRAP here
     (drifts to 1:2 exterior / L4-L5). Use the symmetric perpendicular-crossing
     corrector. T₀=4π is the *linearized bifurcation* period, not the PO period —
     target resonance/crossing structure, let T float.

## Dependency: #437 is the linchpin

The isolated families are reached by continuation that crosses FOLDS in e (both
Peng 2017 and Martínez-Cacho 2025 confirm naive continuation fails at folds).
**#437 (fold-aware pseudo-arclength continuation) is required for #440 Phase 2**
and is now golden-specified: Peng Tables 2/3 (Sun-Mercury ME-Halo ICs +
eigenvalues) and Martínez-Cacho App. B + the e=0.0324 fold (Sun-Mars 2:3) are
same-model regression goldens. Build #437 first; #440 Phase 2 consumes it.

## Phased build

### Phase 1 — e=0 resonant-PO converger, 5 MMRs (independent of #437; ~1 day)
`src/cyclerfinder/search/er3bp_isolated_seeds.py` (new): generalise the validated
3/2 recipe to all five MMRs at µ=0.001 — a1 = 0.763143 (3/2), 0.5428 (5/2),
0.4807 (3/1), 0.3968 (4/1), 0.3419 (5/1); period_guess = kT₀ per the paper.
**Golden:** each converges a symmetric resonant PO landing at its a1 within the
finite-amplitude tolerance (~0.5%) on the correct n=p/q ratio. (Sourced: the a1
gap locations are published; assert the resonance ratio + closure, not an
unsourced IC.)

### Phase 2 — continue to e>0 + isolated-family detection (needs #437; ~3–5 days/MMR)
Continue each e=0 member into our pulsating `e` (= their e₂) with the #437
fold-aware continuator, in the symmetric configs ((0,π) for 3/2 & 5/2; (π,0) for
3/1, 4/1, 5/1). For an ISOLATED family: directly seed at high e (digitised from
the paper's `(e₁,e₂)` figures), converge, then attempt continuation toward e=0
and CONFIRM it dies before e=0 (no circular limit). Distinguish from
connected-family members (which reach e=0). Floquet-classify; literature-check.

### Phase 3 — verdict + (if a genuine isolated family) gauntlet + registry
Report-only first. A reproduced isolated family in the µ=0.001 planetary model is
a CAPABILITY win (we can find the class). The novel-cycler prize is whether the
**Earth-Moon cycler** class has isolated e>0 members — a separate seeding pass at
Earth-Moon µ once the planetary-model method is validated. Hard-scrutinise any
"isolated" claim (the #436 lesson: re-verify the death is a real no-circular-limit
boundary, not a continuation artifact).

## Honest likelihood
MEDIUM. The phenomenon provably exists (Antoniadou & Libert), so this is not a
blind hunt — but (a) it's in the planetary MMR model first, not Earth-Moon
cyclers; (b) seeds must be digitised from figures (precision-limited); (c)
whether the Earth-Moon CYCLER class has isolated members is genuinely open. The
capability (find isolated families at all) is high-confidence given the GO gate;
the novel-Earth-Moon-cycler prize is the speculative part.

## Conventions
main; `uv run` ruff+mypy; no Co-Authored-By; pathspec; report-only until gauntlet;
golden EXPECTED sides sourced (a1 gaps, resonance ratios) never code-computed.
