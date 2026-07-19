# #653 — W-Z proof-machinery scoping for the #646 SE↔EM non-connection (2026-07-19)

Decision-support scoping only (no code, no proof build, no catalogue changes), per the `#653`
registration. Question: does `#646`'s strengthened SE↔EM non-connection negative warrant
dispatching `#636`'s multi-week Wilczak-Zgliczyński (W-Z) computer-assisted-proof build?

Inputs actually read (not summarized from memory): the full `#636`/`#646`/`#619`/`#610`/`#625`
bullets and the `## CURRENT STATE` `#538`-arc entry in `data/OUTSTANDING.md`;
`docs/notes/2026-06-19-digest-wilczak-zgliczynski-oterma-heteroclinic.md`;
Wilczak & Zgliczyński Part I (arXiv:math/0201278) directly — Sections 1-3 (setting, Poincaré
sections, h-sets/covering relations), Section 8 (numerical proof parameters, Tables 4-9),
Section 9 (the `pcr3bp` program / CAPD), and Section 10 (concluding remarks — load-bearing, see
below); `data/golden/wz_oterma_heteroclinic.yaml` + `tests/genome/test_heteroclinic_cycle.py`
(to establish what `#403` actually built).

**Recommendation up front: NO-GO** on pointing W-Z machinery at `#646`'s negative. Keep `#636`
parked, but with its trigger re-scoped: this machinery's natural use in this project is
certifying a future **positive** (a found connection), which is exactly the shape W-Z solved —
not certifying a void, which is not.

---

## 1. The honest formal proposition

### What #646 actually established (numerically)

In the QBCP (time-periodic Sun-Earth-Moon model, post-`#592` α₆ fix), with the GMOS SE-L2 torus
(invariance residual 1.54e-5) and the `#618` pseudospectral EM-L2 torus (C≈3.13, n1=28, n2=9,
invariance rms 9.474e-4), using a CLV-verified, perturbation-stable EM-L2 unstable-manifold
direction: **multistart local optimization** of the 12-unknown/18-residual two-leg ballistic
connection functional floors at a ~166,016 km leg-2 gap (390.9 m/s) in every sampled basin,
integrator-independent (DOP853 vs Radau < 0.001 km). Independently corroborated by `#620`/`#626`
collocation (best Radau-confirmed gap 351,467 km over 24 seeds × 3 orders).

Note what this is **not**, even numerically: it is not a global statement. Multistart sampled
finitely many basins; the crossing scan covered specific section windows and TOF ranges; the
tori sit at one energy/rotation-number pair each; the manifold displacement ε was fixed small.
Three independent methods agreeing on a large floor is strong *evidence* of non-existence in
the searched regime — it is not, and does not claim to be, coverage of a domain.

### The translation to a provable proposition

A theorem cannot be about "no SE↔EM cislunar cycler exists." The strongest honest formal
proposition a rigorous computation could target is a **bounded reachability-avoidance claim**:

> **P.** Fix the QBCP vector field exactly as implemented (μ, μ_S, α-coefficients per
> Rosales-Jorba 2023). Fix explicit compact sets:
> U = a slab around the SE-L2 torus's unstable-manifold departure section, parameterized by
> torus phases θ ∈ T² and offset ε ∈ [ε₁, ε₂];
> S = a closed neighborhood of the EM-L2 torus, with radius large enough to *provably* contain
> the relevant piece of its local stable manifold;
> T_max = an explicit time-of-flight horizon.
> **Claim:** no QBCP trajectory starting in U enters S within time ≤ T_max (and symmetrically
> for the return leg). **Corollary:** no ballistic two-leg cycler with legs in these windows.

Everything outside that box is *not* covered by P: other torus energies/rotation numbers along
both families, longer TOFs, powered connections, multi-leg/lunar-flyby topologies, and — most
importantly — the real ephemeris (P is a theorem about the QBCP, an approximate model; `#592`
found a real EOM bug in this very model in July 2026, so a pre-`#592` "theorem" would have been
rigorous about the wrong equations). The domain quantifiers are not decoration: they are most
of the claim, and they make even a successful theorem dramatically narrower than the informal
"no connection exists" the phrase "theorem-grade negative" suggests.

There is also a hidden prerequisite inside P's *statement*: "the EM-L2 torus" and "its stable
manifold" must be well-defined mathematical objects, i.e., the invariant tori themselves must
be **rigorously validated** (a-posteriori KAM-style existence + enclosure), not just numerically
corrected. See feasibility.

## 2. Feasibility — the shape mismatch is fundamental, not incidental

### What the existing W-Z machinery actually is

Part I proves *existence*: two hyperbolic Lyapunov orbits and homo/heteroclinic connections
between them, in the **planar CR3BP** — autonomous, 2 DOF, with the Jacobi first integral, so
the whole proof lives on a 2D Poincaré section (x, ẋ) of a 3D energy manifold (Part I, Sections
1-2). The proof is a **finite chain of covering relations between 2D h-sets placed along
already-numerically-known orbits** (Sections 3-6), verified with rigorous interval enclosures of
section maps via the CAPD C++ package's C⁰/C¹-Lohner algorithms (Sections 8-9; whole proof
< 32-40 min on a 1.1 GHz Celeron). The topology does the heavy lifting precisely so that the
rigorous numerics only ever has to propagate a few dozen small boxes along one known trajectory
(Tables 4-9: grids of 1-330 subdivisions per edge, in 2D).

Two passages in Part I are load-bearing for this scoping:

1. **Section 1 (p.4) and Section 10 (p.45):** the authors explicitly did *not* rigorously
   compute stable/unstable manifolds even for their own existence proof — "a rigorous
   computation of stable and unstable manifolds for our problem appears to be very difficult
   (requires very extensive C¹-computations). Hence developing tools which avoid a direct
   computation of invariant manifolds is of interest." That is in the *easiest imaginable
   setting* (planar, autonomous, one integral, hyperbolic periodic orbits, 2D sections).
2. **Section 10's stated extension scope:** "this method can be applied to prove a symbolic
   dynamics in any system for which numerical simulations indicate an existence of some kind of
   hyperbolic behavior" — i.e., the method's growth path is *more existence proofs along known
   orbits*. Nothing in Part I or Part II offers a non-existence mode, and the covering-relation
   formalism has none: a failed covering check is silence, not a theorem.

### Why non-existence is a different (and much harder) technique

- **Existence:** exhibit *one* finite chain of boxes along *one* known orbit. Cost is local and
  linear in chain length. Hyperbolicity *helps* (expansion realigns h-sets).
- **Non-existence over a region:** show the rigorous outer enclosure of the reachable set of the
  *entire* departure domain U avoids S for *all* t ≤ T_max. That is validated global
  reachability (set-oriented/GAIO-style with interval enclosures) — a known technique class, but
  with exponential cost in domain dimension and in hyperbolicity, because here expansion
  *hurts*: every box's image inflates by ~|λ_u| per period unless subdivided, and this system's
  one-period amplification at EM-L2 is the measured **~2-3e4** (`#619`), with strong
  non-normality on top (σ_max/|λ_u| up to ~21, `#646`) — interval box growth follows the
  singular values, i.e., is *worse* than the Lyapunov numbers suggest.

Concretely, P's domain is (θ₁, θ₂, ε, TOF) per leg — 4D per leg, two coupled legs (the actual
corrector has 12 unknowns) — in a **6D, time-periodic, no-first-integral** phase space (the QBCP
has no conserved Jacobi constant, so the one classically *rigorous-for-free* non-existence
argument in this domain — closed Hill regions / zero-velocity barriers, as in the PCR3BP of
W-Z's Figure 1 — is unavailable). Each branch-and-bound leaf requires a validated enclosure of a
weeks-to-months model-time propagation with near-Moon passages (`#619` found these arcs stiff
enough to need a collision-guarded propagator even in floating point). W-Z's proof propagated
tens of small 2D boxes for short section-to-section times; P needs on the order of 10⁶-10⁹
long-time enclosures over a 4-8D domain through a 2e4×/period saddle. This is not "the Oterma
machinery, extended" — it is a different, frontier-grade computation that happens to share an
integrator.

### The prerequisite that likely kills it regardless

To state P about "the stable manifold of the EM-L2 torus," the torus must first be rigorously
validated (Figueras-Haro-Luque-style a-posteriori KAM). Such validations need defect residuals
near machine precision; `#618`'s torus sits at rms **9.5e-4** — five-plus orders of magnitude
away, and that residual was itself the hard-won result of crossing the `#544` wall. Validating a
partially-hyperbolic 2D torus in a time-periodic 6D model at this conditioning is an open
research question, not an engineering task. The fallback — inflate S until soft arguments cover
the stable manifold without a validated torus — weakens P further and still leaves the
reachability computation intact.

### Infrastructure reality check

The `#636` bullet's "unusually good preconditions" claim is about *model familiarity*, not proof
infrastructure. `#403`'s golden (`data/golden/wz_oterma_heteroclinic.yaml`,
`tests/genome/test_heteroclinic_cycle.py`) is a **floating-point** reproduction of W-Z's
published crossing coordinates with this project's ordinary corrector/propagator — no interval
arithmetic, no Lohner algorithm, no h-set verification. The codebase's entire rigorous stack is
`mpmath.iv` applied to **closed-form** expressions (`#610`/`#625` bend-gate certificates), and
`#610` itself reported that rigorously bounding even a multi-rev Lambert *solve* is unsolved.
There is no rigorous ODE integrator in the codebase; W-Z's proof rests on CAPD (C++). A build
would start with either CAPD bindings or a from-scratch validated integrator — weeks of work
before the first covering relation, let alone the reachability computation.

### Feasibility verdict

**A rigorous non-existence result of shape P is not a tractable extension of the machinery this
project has or of the machinery W-Z published.** It requires (a) validated quasi-periodic tori
at residuals ~5 orders beyond what is achievable, (b) a validated global reachability
computation whose cost is driven by the very 2e4×/period hyperbolicity that made this system
numerically hard, in a model with no first integral to shrink the domain, and (c) a rigorous
integration stack the codebase does not have. Each is research-frontier alone; P needs all
three.

## 3. Effort/risk estimate (conditional — recorded for completeness, not endorsement)

- **Order of magnitude: months (3-6+), research-grade**, with — my honest estimate — **>50%
  risk of no theorem at the end** (enclosure blow-up through the saddle/near-Moon passages, or
  the torus-validation prerequisite failing outright at rms 9.5e-4).
- **New infrastructure beyond the Oterma reproduction** (which contributes nothing rigorous):
  CAPD bindings or a validated Taylor/Lohner integrator for the QBCP (weeks); rigorous
  interval-QBCP coefficients (the α-series); a torus-validation layer (open research); a
  branch-and-bound reachability driver with checkpointing (per
  [[feedback_long_runs_acceptable]], the compute itself could be days-to-weeks — that is the
  cheap part).
- **Biggest risk to the estimate:** the torus-validation prerequisite is not schedulable — it is
  an open question, and if it fails, the fallback proposition (inflated-S reachability) may not
  be worth proving.
- **The one cheap rigorous stepping stone** that exists: certify a *local* lower bound ("over
  this small box around `#619`'s converged basin, the gap exceeds 10⁵ km") — one box, two
  validated legs. But it adds almost nothing over the existing DOP853-vs-Radau cross-check and
  still requires the full validated-integrator build first. Not worth it on its own.

## 4. GO / NO-GO

**NO-GO.** Reasoning, made explicit:

1. **Shape mismatch (decisive).** `#636`'s own caveat, confirmed here against the primary
   source: W-Z machinery proves existence; non-existence over a propagation-dependent region is
   a different, far harder technique (validated global reachability), not an extension. W-Z
   themselves avoided rigorous manifold computation as too hard *for existence, in a 2D-section
   autonomous problem* (Part I, Section 10).
2. **The prerequisite is likely infeasible.** No validated tori at rms 9.5e-4; without them the
   theorem's subject doesn't formally exist.
3. **Low marginal epistemic value.** The negative already has a 166,000 km floor — 8+ orders
   above closure tolerance — confirmed by three methodologically independent attacks and an
   independent integrator. Rigor upgrades marginal calls, not 166,000 km ones. And even a
   successful P would remain conditional on the QBCP model (per `#592`, the model itself is the
   weakest link — a theorem cannot fix model error), and on domain bounds that exclude most of
   what "no connection" colloquially means. The negative-results registry is *method-versioned
   by design* ([[project_negative_results_registry]]); "empty is conditional on method" is the
   registry's stated epistemic status, and P would still be conditional — just expensively so.
4. **Program context.** Per [[project_capability_frontier_complete]] /
   [[project_validation_ceiling]], the program is in discovery/census mode; nothing downstream
   is gated on upgrading this negative. [[feedback_speculative_high_effort_required]] rightly
   forbids rejecting multi-week builds on cost alone — this NO-GO rests on the technical shape
   mismatch and the low marginal value, not on cost.

**What would flip this to GO** (recommend recording on `#636`'s bullet as its re-scoped
trigger):

- **A found connection needing certification.** If any future search *closes* a high-stakes
  connection (here or elsewhere), a W-Z covering-relations *existence* certificate along the
  found orbit is exactly the published shape — finite chain, known trajectory, hyperbolicity
  helping — and the realistic build is then weeks (CAPD-based, autonomous-model cases first),
  not months. That is the machinery's natural role in this project.
- **A publication goal in its own right.** If a computer-assisted-proof paper becomes a project
  deliverable, the calculus changes — but that is a program decision, not a `#646` trigger.
- **A marginal negative.** If some future negative floors at kilometers, not 166,000 km, where
  rigor could genuinely flip the verdict, revisit — though the shape-mismatch and
  torus-validation obstacles in Section 2 would still apply and must be re-costed, not assumed
  away.
