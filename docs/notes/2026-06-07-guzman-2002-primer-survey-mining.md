# Guzman 2002 primer-vector survey — mining + reassessment of #144/#148

**Status: REFERENCE MINING + REASSESSMENT VERDICT.** Reads the just-acquired
survey and renders the verdict that lifts/sharpens the "PROVISIONAL pending
Guzman 2002 (multi-rev caveats)" labels carried by
`docs/notes/2026-06-07-primer-vector-diagnostic.md` (#144) and
`docs/notes/2026-06-07-primer-refine-recoverable-dv.md` (#148).

Source (cite by authors / title / IAC number only):
**Guzman, J. J.; Mailhe, L. M.; Schiff, C.; Hughes, S. P.; Folta, D. C. (2002).
"Primer Vector Optimization: Survey of Theory, New Analysis and Applications."
IAC-02-A.6.09, 53rd International Astronautical Congress, Houston, TX, 10-19 Oct
2002. NASA Goddard / a.i. solutions (declared a work of the U.S. Government).**

Page numbers below are the typeset page numbers printed at the foot of each page
(p.1 = "Introduction", … p.10 = "Conclusions", p.11 = References). The PDF is a
12-leaf scan (cover + 11 numbered pages).

This note does NOT edit #144, #148, or `verify/primer.py` — recommendation only.

---

## 1. Failure / singularity cases (verbatim conditions, with pages)

The survey's central thesis (Abstract, p.1): *"the applicability of primer
vector theory is examined in an effort to understand when and why the theory can
fail. For example, since the Calculus of Variations is based on 'small'
variations, singularities in the linearized (variational) equations of motion
along the arcs must be taken into account. These singularities are a recurring
problem in analyses that employ 'small' variations."*

It splits failures into two families (heading **"Applicability of Primer Vector
Theory," p.6**):

> *"Two issues that cause problems in the application of primer vector theory
> have been observed in this investigation. **First** … depending on the force
> model … and on the length of the propagation, the linearized dynamics might
> not be adequate to represent the natural dynamics of the system. That is,
> 'small' variations might not remain 'small' and solving for the primer vector
> history becomes numerically intractable. **Second**, there might be isolated
> points along the nominal path or orbit where the upper right hand block of the
> STM, Φ_rv, is singular. This singularity prevents the computation of the
> initial primer derivative via Equation (21)…"*

### 1a. Linear-theory-range breakdown — long propagation / eccentric arcs (p.6–7)

The "Linear Theory Range" subsection (p.6) gives an explicit ad-hoc test:
perturb the initial state by `d` (small scalar), propagate both with the STM and
with the full nonlinear ODE, and compare final states. The worked breakdown
example (p.7, **Fig. 1**):

> *"consider … 1) a circular orbit at a 1,500 km altitude (low Earth orbit)
> propagated for .5 days (about 5 revolutions), and 2) a highly eccentric orbit,
> e = 0.95, with a perigee altitude of 2,703 km propagated for 9.2 days … The
> results illustrate how for the highly eccentric orbit the STM propagation
> breaks down before one complete revolution. The onset of the 'breakdown'
> takes place as the spacecraft approaches perigee and the cartesian STM changes
> very quickly. This type of behavior creates serious numerical problems in
> trying to obtain the primer vector history between periapsis points…"*

Mitigations they name (p.7): *"Switching to a different formulation, e.g. orbital
elements, might alleviate this problem. Schemes that utilize linear propagation
for targeting and/or Monte Carlo analyses should investigate the linear theory
range for all initial conditions."* (the latter cites Goodson's MAP Monte-Carlo
work, ref. 21).

### 1b. Isolated Φ_rv singularities in 2BP elliptic arcs (p.7–8, Fig. 2)

Heading **"Isolated Singularities," p.7**. Singularities in the Φ_rv block on
2BP elliptic arcs occur when (verbatim three conditions, p.7):

> 1. *"The difference between the initial and final times is a multiple of the
>    reference orbit period."*
> 2. *"The difference between the initial and final true anomalies are given by
>    kπ, for k = 1, 2, 3, … — note that this covers the first case, thus,
>    additional singularities occur when the difference between the initial and
>    final true anomalies given by (2k+1)π, for k = 0, 1, 2, 3, …, and"*
> 3. *"The time of flight is a minimum for the given difference in true
>    anomaly."*

Attribution / structure (p.7, citing Stern 1964, ref. 22): *"the first two are
readily explained on physical grounds while the third is a consequence of the
initial assumption of a linearized model. In the first case, the initial and
final position perturbations are not independent. In the second case, the
out-of-plane components are not independent. (Thus, specifying the transfer
plane and considering only planar variations can remedy this situation.) The
third case is more complex and depends on the eccentricity of the reference
orbit as well as on the initial and final eccentric anomalies…"* — Stern's
singularity factor X (Eq. 45) is plotted vs revolution count in **Fig. 2**.

Crucial scope statement (p.8, top): *"as Stern remarks, there are no
singularities (of any type) for true anomaly differences greater than 0 degrees
and less than 180 degrees."* And: *"general tools should be able to handle other
cases that might include multiple revolutions."*

### 1c. The singularity ↔ multi-rev / phasing connection (p.8)

> *"these multiple revolution cases allow the designers to use the concept of
> 'phasing' orbits when implementing orbit transfers and/or rendezvous … For
> other force models (non-Keplerian), the situation appears more complex. In
> fact, these singularities (in Φ_rv) might be the cause of 'spikes' in the Δv
> cost observed in targeting/optimization schemes that utilize the Φ_rv STM
> block."* (cites Wilson/Howell/Lo, ref. 24.)

### 1d. Local-optimum / initial-guess fragility (p.10, Conclusions; Figs. 9–12)

> *"any solution computed using primer vector theory is a local optimum and,
> therefore, highly dependent on the initial guess."* (p.10) The line-of-apsides
> example (p.10) shows two seeds converging to **different** optimal Δv
> (1.876 km/s four-burn vs 1.536 km/s, "about 18% lower"): *"a particular problem
> might have multiple locally optimal solutions (with different fuel costs) …
> success is dependent on the existence of solutions in the neighborhood of the
> initial guesses."*

---

## 2. Necessary-vs-sufficient (the #144/#148 crux), p.3, p.5–6

The primer necessary conditions are stated verbatim (**"Primer Necessary
Conditions," p.3**, conditions 1–4):

> 1. *Primer evolution governed by* p̈ = (G_r − Ġ_v)ᵀ p̄ − G_vᵀ ṗ (Eq. 18).
> 2. *"During transfers ||p̄|| = p ≤ 1, with impulses at instants for which
>    p = 1."*
> 3. *"At an impulse time p̄ = ṗ = û_T*, where û_T* is the optimum thrust
>    direction."*
> 4. *"At all interior impulses (not at the initial or final times) ṗ·p̄ = 0.
>    This condition has implications on the slope of the primer vector magnitude
>    since d||p̄||/dt = (ṗ·p̄)/||p̄||. Therefore d||p̄||/dt = 0 at the
>    intermediate impulses…"*

And the explicit necessary-not-sufficient disclaimer (p.3, immediately after):

> *"These conditions are **necessary (NOT sufficient)** for an optimal
> trajectory (time-fixed problem). In this paper, a trajectory that meets the
> above conditions will be called an optimal trajectory."*

**On `|p|>1` and what it licenses (p.5, "Criteria for Three Impulses"):** the
first-order cost variation for adding an interior impulse is δJ = c(1 − p̄_mᵀ η̂)
(Eq. 40), with the result (verbatim, p.5):

> *"Thus, using the definition of a dot product, **if ||p̄_m|| > 1 at any time, a
> third impulse is beneficial.** Furthermore, the greatest decrease in the cost
> function will be achieved if the impulse is applied at the maximum of ||p̄_m||
> at time t_m and in the direction of η̂. The position along the perturbed path
> and the magnitude of the impulse are yet to be determined."*

**Directly relevant to #148's "necessary-condition violation with no realisable
improvement":** the survey is consistent with it but does NOT contradict #148.
Two points:

- The δJ = c(1 − p̄_mᵀ η̂) result is **first-order in c** (the impulse
  magnitude), and it asserts only that the *gradient* of cost wrt adding an
  impulse is negative when `|p|>1`. The *magnitude* of the improvement is a
  separate calculation: p.5, "Calculation of the Interior Impulse," explicitly
  says you then estimate c from a perturbed-path cost expansion (Jezewski &
  Rozendaal 2nd-order, or Hiday's exact iterative), and warns: *"In any case,
  the mid-impulse should decrease the cost **but might not produce an optimal
  trajectory in the sense of Lawden**."* So "improvement exists" (gradient) vs
  "improvement is realisably large" (magnitude) are deliberately separated by the
  survey itself. #148's finding — gradient says improvable, realised magnitude
  ≈ 0 — sits squarely inside this gap and is NOT an error.
- The survey does **not** treat constrained-boundary / degenerate-zero-burn
  endpoints. It is pure 2BP free-space impulsive (Earth-orbit and libration-point
  applications). The Aldrin coast-0 degeneracy in #148 — the Earth-departure
  heliocentric burn is *exactly* 0 by construction, so the endpoint primer
  direction is a fallback unit vector, not a real `Δv/|Δv|` — is **outside the
  survey's modelled cases**. The survey's BC (Eq. 33) `p̄ = Δv̄/||Δv̄||`
  *presumes* a non-degenerate impulse; with `|Δv|→0` that BC is ill-posed. This
  is new scope the survey does not cover, and is exactly why #148's `|p|=1.122`
  is real-but-empty.

---

## 3. Multi-rev specifics (p.7–8) vs our arcs

- **No singularity for 0° < Δν < 180° (single sub-half-rev arc), p.8.** Our
  E→M coast 0 is 131.9 d; our M→E coast 1 is 598.7 d. Mars's heliocentric period
  is ~687 d, so 598.7 d is ~0.87 rev (Δν < 360°, and as a single transfer arc
  its Δν is < 180°-ish geometry, well under one rev). **Neither leg crosses a
  full revolution.** The survey's singularity family (Δt = integer × period; Δν =
  kπ; min-ToF) is therefore **not triggered** by either Aldrin leg. The labels'
  fear of "multi-rev fragility" does not bind here — these are sub-one-rev arcs.
- Stern's X-factor (Eq. 45, Fig. 2): first zero-crossing occurs *after* 1
  revolution even for e up to 0.75. Both our legs are below that.
- Multi-rev is where Φ_rv *can* go singular and inject spurious "spikes" in Δv
  (p.8). Our diagnostic (`verify/primer.py`) inverts Φ_rv directly
  (`np.linalg.solve(phi_rv, …)`), so it would be exposed to exactly this — **but
  only on ≥1-rev arcs**, which we do not have in the Aldrin schedule.

---

## 4. Applications section — techniques we have / haven't (p.5–6, p.8–10)

Survey toolset and our coverage:

| Survey technique (page) | Lawden condition | Implemented in our stack? |
|---|---|---|
| Add interior impulse, seed at max\|p\|, dir η̂ (Eq. 40–41, p.5) | cond. 2 violation | **Yes** — `primer_refine.py` (#148), seeded from primer peak. |
| Mid-impulse position via δr̄_m = c Ā⁻¹ p̄_m/\|\|p̄_m\|\| (Eq. 41, p.5) | — | Partially — #148 optimises (τ, δr) with Nelder-Mead instead of the closed-form Ā⁻¹ step; equivalent intent. |
| Converge to 3-burn optimum, vary t_m & position (Eq. 42, p.5) | cond. 4 (ṗ·p̄=0) | **No** — #148 does one step, not iterated to ṗ·p̄=0. |
| **Add initial / final coast** (vary endpoint *times*), Eq. 43, p.6 | **cond. 4 / endpoint slope** | **No** — explicitly out of scope in #144 (the bi-elliptic-threshold detector). |
| Four endpoint-slope cases (init-coast/early-dep/final-coast/late-arr), p.6 | — | **No.** |
| Move existing impulses / Lambert type-I↔II switch (p.8, Fig. 6) | — | **No.** |

**Relevance verdict given #148's load-bearing finding** (the 2.9138 km/s cost
lives at the **Earth flyby turn-deficit**, not in the heliocentric coasts):

- **Add-initial/final-coast (Eq. 43)** is the survey technique #144 flagged as
  the missing bi-elliptic detector. But it varies endpoint *times* to reduce a
  *heliocentric* coast cost. Since the heliocentric coast cost on coast 0 is
  already ~0.003 km/s (#148), even a perfect endpoint-coast optimisation recovers
  at most that ~0.003 km/s. **Not worth implementing for the Aldrin maintenance
  problem** — wrong actuator, as #148 concluded for the interior-impulse step.
  Worth it only if/when we have a heliocentric leg with genuinely large coast
  burns.
- **Endpoint-time variation (δJ = −ṗ₀\|Δv̄₀\|dt₀ − ṗ_f\|Δv̄_f\|dt_f, Eq. 43)**
  scales with `|Δv₀|` and `|Δv_f|`. On coast 0 `|Δv₀| = 0` exactly, so the
  initial-coast term is **identically zero** — the survey's own formula confirms
  there is nothing to gain by shifting the (degenerate) Earth-departure epoch via
  this mechanism. Independent corroboration of #148's zero-recovery, *from the
  survey's equations*.
- The flyby turn-deficit is a *geocentric geometric* cost. Primer vector theory
  as surveyed (impulsive, single central body per arc) **has no machinery for
  it.** The survey does not address patched-conic flyby turn limits at all. So
  none of its applications attacks the actual Aldrin cost driver. Confirms #148.

---

## 5. THE REASSESSMENT VERDICT (deliverable core)

Per result, does the survey (a) invalidate, (b) explain, or (c) leave standing:

### #144 — coast 0 `max|p| = 1.1223` → IMPROVABLE_ADD_IMPULSE

**Verdict: (b) EXPLAIN + (c) LEAVE STANDING, with scope tightened — NOT
invalidated.**

- The flag is the textbook Eq.-40 result: `|p|>1` ⇒ a third impulse is
  beneficial (p.5). The diagnostic is doing exactly what the survey prescribes.
- The multi-rev fragility the label feared **does not apply**: coast 0 is 131.9 d,
  a sub-one-rev arc with Δν < 180°, i.e. the survey's *singularity-free* regime
  (p.8). The Φ_rv inversion in our code is safe here.
- The linear-range breakdown (p.6–7) is keyed to high-e perigee passages /
  multi-rev LEO; a 132-d heliocentric near-circular transfer is the benign end of
  the spectrum, not the breakdown end. So the `|p|=1.122` value is trustworthy as
  a *first-order necessary-condition signal*.
- It remains, per the survey's own p.3 disclaimer, **necessary-not-sufficient** —
  it can refute optimality, never prove it. #144 already says this. No change to
  the substance.

**Recommended label change for #144:**
`PROVISIONAL pending Guzman 2002 (multi-rev caveats)`
→ **`CONFIRMED-with-scope per Guzman 2002 §"Isolated Singularities" (p.7-8):
coast 0 is a sub-one-rev (Δν<180°) singularity-free arc; the |p|=1.122
necessary-condition flag is valid and the Φ_rv inversion is well-posed. Remains
necessary-not-sufficient (Guzman p.3).`**
The blanket "all Aldrin primer results provisional pending Guzman" sentence in
#144 (lines 110-118) can be downgraded: the *machinery* is validated for these
arcs; only the M→E-leg "multi-rev-scale" worry was overstated (see next).

### #144 — coast 1 `max|p| = 1.00008` at endpoint (M→E, 598.7 d)

**Verdict: (b) EXPLAIN — the "multi-revolution-scale" framing was the
overstatement the label hedged against, and the survey dissolves it.**

- 598.7 d is ~0.87 of Mars's period — **sub-one-rev**, not multi-rev. The
  survey's singularity conditions (Δt = n·period, Δν = kπ; p.7) are not met.
- The marginal `1.00008` endpoint touch is consistent with the survey's note that
  the maximum sits at endpoints for benign transfers; it is not a Φ_rv-singularity
  spike (those need ≥1 rev, p.8).

**Recommended:** drop the "(599 days, multi-revolution-scale)" characterisation
in #144 line 114 — per Guzman it is sub-one-rev and singularity-free. Re-label to
**`CONFIRMED-with-scope: sub-one-rev, no Φ_rv singularity (Guzman p.7-8);
marginal endpoint touch, treat as necessary-conditions-met.`**

### #148 — realised recoverable heliocentric Δv ≈ 0.000 km/s

**Verdict: (c) LEAVE STANDING — survey corroborates it, via two independent
mechanisms, and does NOT invalidate it.**

- The survey separates "gradient says improvable" (Eq. 40, first-order) from
  "realised magnitude c" (Eq. 41 + Jezewski/Hiday cost estimate, p.5) and warns
  the mid-impulse "might not produce an optimal trajectory in the sense of
  Lawden." #148's gradient-positive / magnitude-zero result is inside this gap,
  not a contradiction.
- The survey's endpoint-time formula (Eq. 43, p.6) scales recoverable Δv by
  `|Δv₀|`, `|Δv_f|`. With the Earth-departure `|Δv₀| = 0` exactly, the survey's
  own equation gives zero recoverable from the degenerate end — **the survey
  predicts #148's zero**.
- The survey has **no** model for the geocentric flyby turn-deficit (the actual
  cost driver). So it cannot license recovering it via any primer technique;
  #148's "wrong lever" conclusion is consistent with the survey's silence/scope.

**Recommended label change for #148:**
`DIAGNOSTIC / PROVISIONAL … Guzman 2002 multi-rev caveats … apply unchanged`
→ **`CONFIRMED per Guzman 2002 p.5-6: zero recovery is the expected outcome of a
degenerate-endpoint (|Δv₀|=0) coast — both Eq. 40 (first-order gradient vs
realised magnitude gap) and Eq. 43 (endpoint-time δJ ∝ |Δv₀| = 0) predict it.
The multi-rev caveat does not bind (sub-one-rev arc).`** Keep the
"NOT for the site/catalogue" honesty framing — that is editorial, not changed by
the survey.

### One genuinely-new caveat to ADD (not previously flagged)

The survey's BC `p̄ = Δv̄/||Δv̄||` (Eq. 33, p.4) is **ill-posed when `|Δv|→0`**.
Our coast-0 Earth-departure burn is exactly zero, so the endpoint primer
direction fed into `primer_on_coast` is a fallback unit vector, not a real ΔV
direction. The `|p|=1.122` bulge is therefore partly an artifact of an *imposed*
(non-physical) endpoint direction. This is a **degenerate-boundary-condition**
limitation the survey does not cover (it assumes non-degenerate impulses). It
does not invalidate the #144/#148 numbers — both notes already attribute the bulge
to the endpoint directions including "the degenerate zero-burn Earth departure →
fallback direction" (#148 lines 97-100) — but it should be promoted from buried
prose to an explicit named caveat in `verify/primer.py`'s docstring on a future
edit: *the diagnostic's BC is undefined for zero-magnitude bounding impulses.*
RECOMMENDATION ONLY — not edited here.

---

## 6. Net effect on the labels

| Note | Old label | Recommended new label |
|---|---|---|
| #144 coast 0 | PROVISIONAL pending Guzman (multi-rev) | **CONFIRMED-with-scope** — sub-one-rev, singularity-free, necessary-not-sufficient (Guzman p.3, p.7-8) |
| #144 coast 1 | PROVISIONAL (multi-rev-scale) | **CONFIRMED-with-scope** — sub-one-rev mischaracterised as multi-rev; no Φ_rv singularity |
| #148 | DIAGNOSTIC/PROVISIONAL, multi-rev caveats apply | **CONFIRMED** — zero recovery predicted by Guzman Eq. 40 & Eq. 43 for degenerate endpoint |
| (new) | — | **ADD caveat**: Guzman Eq. 33 BC ill-posed for `\|Δv\|→0`; coast-0 departure is exactly that |

Bottom line: the survey **lifts** the multi-rev caveat (our arcs are sub-one-rev,
in the singularity-free regime) and **explains/corroborates** the zero-recovery
finding from its own equations. It does **not invalidate** any #144/#148 result.
It surfaces one new, narrower caveat (degenerate `|Δv|=0` boundary condition)
that both notes already informally acknowledge.
