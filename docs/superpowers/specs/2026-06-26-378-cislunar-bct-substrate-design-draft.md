# Cislunar BCT integration substrate — Belbruno WSB / ballistic-capture into the BCR4BP discovery stack (DESIGN DRAFT, #378)

**Date:** 2026-06-26
**Status:** DESIGN-DRAFT — for user review. No production code written by this doc.
**Issue:** #378. Frontier: **cislunar BCR4BP** (Sun-Earth-Moon coupled) — one of
the genuinely-fertile, under-mined venues now that the strict CR3BP/ER3BP
heliocentric and EM-libration envelopes are largely mined out.

This draft answers ONE question: **what is the smallest, well-bounded capability
that adds Belbruno's weak-stability-boundary (WSB) / ballistic-capture-transfer
(BCT) construction to the existing BCR4BP machinery (#292/#380/#412) so the
discovery stack can reach low-energy cislunar transport objects the strict
CR3BP/ER3BP lanes structurally cannot — and what, honestly, does it return:
a catalogue-class transport object (cycler / quasi-cycler) or a capability
(transfer-construction tool) only?**

---

## 0. TL;DR

- **Goal:** make the discovery stack able to **construct, correct, and classify a
  ballistic-capture transfer (BCT)** in the Sun-Earth-Moon BCR4BP, anchored on
  Belbruno's weak-stability-boundary surface `W` (Belbruno 2004, Def 3.12 /
  eq 3.9 / Lemma 3.21 eq 3.29). The substrate is two new pure modules
  (`core/wsb.py`, `genome/bct_transfer.py`) plus a discovery driver; everything
  downstream reuses the existing BCR4BP propagator, the #380 augmented-quadrature
  corrector seam, the gauntlet, and the lit-check.
- **What it targets:** **low-energy cislunar transfer chains** terminating in
  ballistic capture at the Moon (zero capture ΔV by definition of `W`). The
  in-scope object is the **Hiten-class exterior WSB transfer** (Earth phasing arc
  → ~1.5e6 km Sun-shaped apoapsis → fall-back to ballistic capture at the Moon)
  and the **1987 interior WSB transfer** (L₂-neck transit). The dynamical
  mechanism is the Sun perturbation shaping a high apoapsis from which the
  Moon's L₁/L₂ Lyapunov-manifold network funnels the spacecraft to `W` — exactly
  the geometric vehicle the #412 reach spike found the EM-libration family
  **lacks**.
- **Why strict CR3BP/ER3BP can't reach it:** a BCT is intrinsically a
  **four-body** object — the Sun's gravity is the energy-removal agent during
  fall-back (Belbruno 2004 §3.4.2; Koon 2000 decouples it into Sun-Earth +
  Earth-Moon CR3BPs *precisely because one CR3BP cannot hold the whole arc*).
  Strict CR3BP has no Sun term; ER3BP adds primary eccentricity, not a third
  perturber. The #412 spike confirmed empirically that Sun-perturbing an
  EM-libration CR3BP orbit does **not** stretch its reach to the SE-L region
  (best 1.16 LD vs 3.9 LD target) — the reach must come from a Sun-scale apoapsis,
  which only the BCR4BP-with-WSB construction supplies.
- **Honest catalogue verdict (see §8):** **capability-first, with a conditional,
  non-obvious path to a quasi-cycler row.** A single BCT is a one-shot transfer,
  NOT a cycler — Belbruno 2004 §"KNOWN_CORPUS impact" states this outright ("WSB
  orbits are not cyclers"). The catalogue-class object, if any exists, is a
  **periodically-repeating cislunar capture/escape chain** (a WSB-anchored
  quasi-cycler), which this substrate makes *searchable* but does not guarantee
  exists. Plan Phase 0 is a feasibility spike that asserts no object; the
  catalogue-row claim is gated behind it.
- **Golden target:** the **Hiten / MUSES-A exterior WSB transfer** —
  ΔV_total ≈ 44 m/s (14 + 30), TOF ≈ 150 d, Q_a apoapsis ≈ 1.5e6 km (~3.9 LD),
  capture e₂ ≈ 0.95 at r_M + 100 km, ΔV_capture = 0 — fully sourced in
  Belbruno 2004 §3.4 (digested `docs/notes/2026-06-17-digest-belbruno-2004.md`)
  and cross-checked against the Koon-Lo-Marsden-Ross "Shoot the Moon" manifold
  mechanism (digested `docs/notes/2026-06-21-digest-koon-2000-shoot-the-moon.md`).
- **Grounding honesty:** the **WSB surface** (`W = J⁻¹(C) ∩ Σ ∩ σ`, the analytic
  approximation eq 3.29, the validity domain C < C₁, the numerical stability-class
  algorithm) is **fully digested and buildable now**. The **BCT *targeting*
  recipe** (backward-integration two-arc match; forward 2×2 differential
  correction `(|V₀|, γ₀) → (r₂₃, i_M)`) is digested at *recipe* granularity
  (Belbruno 2004 §3.4.1 + Remark 4) but the convergence-detail papers it cites
  ([33] Belbruno-Carrico 2000, [34] Belbruno-Humble-Coil 1997, [39] JGCD 1993)
  are **NOT in-repo**. Plan Phase 0 task 1 is an **acquisition/digest decision**:
  the forward 2×2 method maps onto the existing #380/BCR4BP Newton corrector and
  is buildable from the digest, but if it stalls, [39] JGCD 1993 is the named
  source to acquire. See §9.

---

## 1. What cislunar transport objects this targets, and why the strict lanes can't reach them

### 1.1 The objects

| Object | Definition | Source anchor | Catalogue class (if any) |
|---|---|---|---|
| **Exterior WSB transfer (Hiten-class)** | Earth-phasing arc I + Sun-shaped ~1.5e6 km apoapsis arc II → ballistic capture at Moon | Belbruno 2004 §3.4, Fig 3.14 | precursor_mga (single-leg transfer; NOT a cycler) |
| **Interior WSB transfer (1987 / SMART-1 template)** | L₂-neck transit orbit, 16 d, ~75,000 km out of EM plane | Belbruno 2004 §3.3.2, ref [25] | precursor_mga |
| **WSB-anchored cislunar quasi-cycler** | a *periodically-repeating* capture↔escape chain whose return leg re-acquires `W` (the only catalogue-class candidate) | NOT directly in Belbruno; the open search question | quasi_cycler — **conditional, unproven** |

The first two are **transfers** (capability outputs). The third is the only
genuine catalogue-row candidate, and it is **not asserted to exist** — it is the
hypothesis the substrate makes testable.

### 1.2 Why strict CR3BP/ER3BP structurally cannot reach them

1. **A BCT is a four-body object.** The Sun's gravity removes orbital energy
   during the fall-back from the 1.5e6 km apoapsis — "prior to capture P₃ moves
   in approximate parallel formation with M for about a week" (Belbruno 2004
   §3.4, p 149). Strict Earth-Moon CR3BP has **no Sun term**, so the
   energy-removal mechanism is absent: an EM-CR3BP orbit conserves Jacobi C and
   cannot ballistically capture from a Sun-scale apoapsis.
2. **ER3BP adds the wrong degree of freedom.** ER3BP perturbs the primaries'
   *mutual eccentricity*; it does not add a third gravitating body. The WSB
   surface `W` is defined by primary-*interchange* (P₃ passes between the Sun-Earth
   and Earth-Moon regimes), which needs the third perturber.
3. **The #412 reach spike is the empirical proof.** It Sun-perturbed four
   EM-L₁-Lyapunov amplitudes via `continue_bcr4bp_family_in_musun` and found the
   reach **flat under μ_sun** (best 1.16 LD) and the family **breaking before full
   Sun strength** (9–12 of 40 steps). Its explicit conclusion: *"a Sun-Earth-scale
   seed … is the geometrically correct vehicle, but that is a from-scratch build,
   not a continuation of the EM family."* **This substrate is that from-scratch
   build** — the WSB surface is the Sun-Earth-scale anchor the spike asked for.
4. **No existing negative re-opens here.** `data/empty_regions.jsonl` contains
   only heliocentric/Jovian patched-conic ballistic negatives (28 `ballistic`
   tags, all `ephem_model: circular`, centres Jupiter/Sun). There is **no
   cislunar / WSB / BCR4BP negative** whose re-open key is a WSB/BCT capability —
   so this is **net-new region**, not a subsumption re-open. (Confirmed by scan,
   §3 data flow.)

---

## 2. The WSB / BCT method (what is being built)

### 2.1 Weak stability boundary `W` (Belbruno 2004 Def 3.12, eq 3.9)

`W = J⁻¹(C) ∩ Σ ∩ σ`, the subset of the Jacobi surface at constant C that is
**simultaneously** (Σ) two-body-bound to the Moon (Kepler energy
`E₂ = ½|Ẋ|² − μ/r₂₃ ≤ 0`, P₂-centred inertial coords) and (σ) at periapsis with
respect to the Moon (`ṙ₂₃ = 0`). A point on `W` is a state from which the
spacecraft is *just barely* gravitationally captured — the boundary between
"falls into a Moon orbit" and "escapes to the Earth region."

Two computable representations, both digested:
- **Analytic approximation** (Lemma 3.21, eq 3.29):
  `C = −r₂(±2√(μ(1+e₂)/r₂) + r₂) + μ(1−e₂)/r₂ + A(r₂, θ₂)` on the
  `(r₂, θ₂, e₂)` periapsis surface, valid for **C < C₁** (Def 3.22; for
  Earth-Moon C₁ ≈ 3.184, `W` exists across C ∈ [2.22, C₁]).
- **Numerical stability-class algorithm** (§3.2.1): from a periapsis state,
  propagate one revolution and label it `{stable, unstable, capture, escape,
  primary-interchange}`; `W` is where the label flips. This is the *ground-truth*
  definition the analytic one approximates.

A "ballistic capture" (Def 3.11) is simply `E₂ ≤ 0` at the target — instantaneous
two-body-bound energy, **not** permanent capture. The capture ΔV is therefore
**zero by construction** when you target a point on `W`.

### 2.2 BCT construction + correction (Belbruno 2004 §3.4.1)

- **Backward-integration recipe** (Fig 3.14): pick QF on `W` at r_M + 100 km,
  e₂ ≈ 0.95; integrate the four-body EOM **backward ~45 d** to a Sun-shaped
  apoapsis Q_a (~1.5e6 km) = **arc II**; integrate **forward ~100 d** from the
  Earth-orbit periapsis Q₀ with a small burn ΔV₀ = **arc I**; **match at Q_a**
  (residual ΔV_a). Total ΔV = ΔV₀ + ΔV_a + 0. "Unwieldy" per Belbruno.
- **Forward 2×2 differential correction** (Remark 4, the modern method): control
  `(|V₀|, γ₀)` at Q₀, target `(r₂₃, i_M)` at QF, prescribe `(t₀, r₁₃, i_E)`. A
  Newton iteration with a 2×2 finite-difference (or STM-derived) Jacobian. **This
  maps directly onto the existing BCR4BP single-shooting corrector** (§4).
- **Four-body CR3BP-pair decomposition** (§3.4.2): when full ephemeris is too
  expensive, `CR4BP ≈ CR3B-EM if r₁₃ ≤ ρ else CR3B-ES`, with
  ρ = d(E,M) + 0.368·d(E,M). A coarse model-reduction; the project's standard
  *incoherent BCR4BP* (`core/bcr4bp.py`) is a strictly higher-fidelity stand-in
  and is the substrate's default model.

### 2.3 How it plugs into the #380 corrector

The forward 2×2 method is a single-shooting Newton solve on a residual
`r(u) = [r₂₃(T) − r_target, i_M(T) − i_target]` with `u = (|V₀|, γ₀)`. The #380
BVP integral-constraint corrector (`search/bvp_integral.py`,
`propagate_augmented_bcr4bp`) already provides **augmented-STM quadrature in the
BCR4BP** — so the **Kepler-energy `E₂(T)` and its sensitivity `∂E₂/∂u`** become an
*integral/terminal* constraint row appended to the Newton system, letting the
corrector target "land exactly on `W` (E₂ = 0)" as a closure constraint rather
than a post-hoc check. The point-constraint rows `(r₂₃, i_M)` reuse the existing
BCR4BP STM (`bcr4bp_stm_eom`). **No new corrector is built** — the BCT corrector
is a thin assembler of constraint rows over the #380 + `bcr4bp_genome` seam.

---

## 3. Architecture, new files, seams, data flow

### 3.1 New files (all additive, pure modules)

| File | Responsibility | Reuses |
|---|---|---|
| `src/cyclerfinder/core/wsb.py` | WSB surface: `kepler_energy_moon(state, system)`, `is_periapsis(state)`, `wsb_analytic_C(r2, theta2, e2, mu, branch)` (eq 3.29), `wsb_validity_ok(C, C1)` (Def 3.22), `stability_class(state, system, ...)` (numerical §3.2.1 one-rev label) | `core/bcr4bp.py` (propagate, system constants, `sun_commensurate_period`); `core/cr3bp.py` (Jacobi, L-point C values) |
| `src/cyclerfinder/genome/bct_transfer.py` | `BCTTarget` (QF on W), `BCTArc`, `construct_bct_backward(...)`, `correct_bct_forward(...)` (the 2×2 `(|V₀|,γ₀)→(r₂₃,i_M)` Newton + the E₂-on-W constraint row), `BCTResult` (arcs, ΔV breakdown, capture e₂, TOF, apoapsis) | `genome/bcr4bp_genome.correct_bcr4bp_periodic` Newton scaffold; `search/bvp_integral.propagate_augmented_bcr4bp`; `core/wsb.py` |
| `src/cyclerfinder/search/cislunar_bct_search.py` | discovery driver: sweep `(t₀, |V₀|, γ₀, e₂)`, construct+correct BCTs, classify (transfer vs return-leg-re-acquires-W quasi-cycler), emit candidates | `genome/bct_transfer.py`; `search/literature_check.py`; `data/empty_regions.py` |

### 3.2 Seam into the existing BCR4BP genome

The substrate is a **strict superset consumer** of `core/bcr4bp.py` and
`genome/bcr4bp_genome.py`; it does not modify them. Exact attach points
(verified):
- propagation: `core/bcr4bp.py:330 propagate_bcr4bp`, `:286 bcr4bp_stm_eom`,
  `:147 andreu_default`, `:385 sun_commensurate_period`.
- corrector scaffold to mirror: `genome/bcr4bp_genome.py:200
  correct_bcr4bp_periodic` (free_vars / residual_indices Newton pattern,
  Radau independent re-propagation closure check).
- augmented-quadrature for the E₂-on-W constraint row:
  `search/bvp_integral.py` `propagate_augmented_bcr4bp` +
  `correct_with_integral_constraints` (#380 Steps 1-3, committed).
- system constants: `genome/bcr4bp_systems.py:234 build_bcr4bp_system`.

### 3.3 Seam into gauntlet + lit-check

- **Lit-check:** any BCT/quasi-cycler candidate runs through
  `search/literature_check.py` (`CandidateSignature` → `build_queries` →
  corpus + web). Belbruno 2004 and Koon 2000 are the corpus anchors; a candidate
  matching the Hiten signature is **correctly flagged non-novel** (this is the
  golden self-test). A net-new repeating cislunar chain that clears lit-check is
  the only novelty claim the substrate could make.
- **Gauntlet:** a *transfer* is precursor_mga-class; a *quasi-cycler* enters the
  V-gauntlet (`data/validation/`) only if Phase 0 finds a repeating object.
  Per the expanded scope (four-class), a precursor must `inserts_into` a cycler —
  so a standalone BCT is a **capability artifact, not a catalogue row** unless it
  feeds a cislunar cycler. The substrate writes to `data/empty_regions.py` on a
  clean negative (the expected default outcome of Phase 0).

### 3.4 Data flow

```
Belbruno W surface (core/wsb.py)
        │  QF on W (r_M+100km, e2≈0.95)
        ▼
construct_bct_backward ──arc II──►  Q_a (~1.5e6 km, BCR4BP propagate)
        │
correct_bct_forward (2×2 (|V0|,γ0)→(r23,i_M); E2-on-W row via #380)
        │  BCTResult {arcs, ΔV breakdown, e2, TOF, apoapsis}
        ▼
cislunar_bct_search ── classify ──► transfer  ──► precursor_mga (capability)
        │                          └ return-leg re-acquires W? ──► quasi_cycler candidate
        ▼
literature_check ──► novel? ──► gauntlet (V0..) | empty_regions.jsonl (clean negative)
```

---

## 4. Validation / golden target

**Primary golden — Hiten / MUSES-A exterior WSB transfer** (Belbruno 2004 §3.4,
fully digested; values are *published*, hence valid golden-EXPECTED sides per the
sourced-golden rule):

| Quantity | Sourced value | Tolerance (design) |
|---|---|---|
| ΔV_total | 44 m/s (14 + 30) | within ~factor-2 band (model is incoherent BCR4BP vs Belbruno's PR4BP-3D ephemeris) |
| TOF | ~150 d | order-of-magnitude (≫ 5 d Hohmann) |
| Q_a apoapsis | ~1.5e6 km (~3.9 LD) | ±30% |
| capture e₂ | ≈ 0.95 at r_M + 100 km | E₂ ≤ 0 (ballistic) is the binding check |
| ΔV_capture | 0 (ballistic by construction) | exact (definitional) |

**Secondary cross-check — Koon 2000 "Shoot the Moon"** manifold mechanism:
the BCT's arc II should thread the Moon-system L₁/L₂ Lyapunov stable-manifold
tube (Belbruno §3.4.2: at 7 d pre-capture C(E,M) = 3.17466, just below
C₁ = 3.184, the L₁ neck "just opens"). This is a *mechanism* cross-check, not a
numeric golden.

**Unit goldens (no reproduction needed, analytic):**
1. `kepler_energy_moon` sign agrees with the L₂ value Belbruno computes
   (`E₂(L₂) = −1.20187 < 0` for μ ≪ 1, §3.3 Lemma 3.30).
2. `wsb_analytic_C` at the parabolic limit returns C = ±√2 + O(μ) (Lemma 3.34,
   eq 3.39) — a closed-form golden.
3. `wsb_validity_ok` boundary at C₁ = 3.184 (Earth-Moon) matches Def 3.22.
4. The incoherent BCR4BP reduces to CR3BP at μ_sun → 0 (already tested in
   `core/bcr4bp.py`); WSB built on it inherits the reduction.

The **model gap is explicit and honest**: the project's *incoherent* BCR4BP is
not Belbruno's PR4BP-3D-with-DE403. The golden is therefore a **band/ballpark**
reproduction (the ΔV/TOF/apoapsis *signature*), not a bit-exact match — exactly
as `bcr4bp_genome` already documents the O(ε²) POL1 gap. Asserting bit-exact 44
m/s would be fabrication; asserting the *signature* is honest and sourced.

---

## 5. Risks + kill-criteria

| # | Risk | Mitigation / kill |
|---|---|---|
| R1 | **Incoherent BCR4BP too coarse to hold a BCT** (the Sun shaping is O(ε²); the apoapsis geometry may not reproduce) | Phase 0 spike measures whether arc II reaches ~1.5e6 km apoapsis at all under the incoherent model. **KILL** if apoapsis stays < 2 LD across the (t₀, |V₀|) sweep — then the coherent QBCP (un-digested α-tables) is a hard prerequisite and #378 blocks on that acquisition. |
| R2 | **No catalogue-class object** — every BCT is a single transfer, return leg never re-acquires `W` periodically | This is the *expected* default. The substrate still delivers the WSB/BCT **capability** + a lit-check-cleared **clean negative** in `empty_regions.jsonl`. Honest verdict, not a failure. |
| R3 | **Forward 2×2 corrector won't converge from the digest** ([33]/[34]/[39] not in-repo) | Plan Phase 0 task 1 acquisition gate; the backward two-arc method (§3.4.1, fully digested) is the fallback constructor. If both stall → acquire [39] JGCD 1993. |
| R4 | **Lit-check false-novelty** on a re-derived Hiten transfer | The Hiten golden self-test (must flag non-novel) is the guard; ships in Phase 1. |
| R5 | **Concurrent edits** — sibling #449 touches `search/releg_*.py` | Disjoint files (`core/wsb.py`, `genome/bct_transfer.py`, `search/cislunar_bct_search.py`); pathspec commits only. |

**Overall kill-criterion:** if Phase 0 (R1) shows the incoherent BCR4BP cannot
even produce a Hiten-signature apoapsis, the WSB *surface* module
(`core/wsb.py`) still ships as a reusable capability, but the BCT *construction*
half blocks on the coherent-QBCP α-table acquisition and #378 is paused with a
documented prerequisite — **not** force-fit on a model that can't hold the object.

---

## 6. Honest assessment — catalogue-class object vs capability-only

**This is capability-first.** Belbruno 2004 itself is unambiguous (digest
§"KNOWN_CORPUS impact"): *"WSB orbits are not cyclers — they are single-leg
low-energy lunar-transfer trajectories."* A BCT delivers:

- **Definitely:** a reusable **WSB-surface + BCT-construction capability** (the
  `core/wsb.py` + `genome/bct_transfer.py` modules), a Hiten-signature
  reproduction (validation evidence), and the Sun-Earth-scale vehicle the #412
  spike identified as missing.
- **Conditionally (unproven, the real search question):** a **WSB-anchored
  cislunar quasi-cycler** — a repeating capture↔escape chain. This is *not in
  Belbruno* and is *not asserted to exist*; the substrate makes it **searchable**.
  Its catalogue-row potential is **genuinely uncertain and probably low**: the
  chaotic hyperbolic-network result (Theorem 3.58) says capture on `W` is a
  *chaotic* process, which cuts against clean periodicity. The expanded scope's
  `precursor_mga` class is the realistic home for a BCT (a precursor that
  `inserts_into` a cislunar cycler) — but only if such a cislunar cycler is found
  to insert into.

**Bottom line for the user:** fund this as a **frontier capability + targeted
quasi-cycler spike**, not as a guaranteed catalogue-row generator. The honest
expected yield is: 1 new capability, 1 Hiten reproduction, and (most likely) 1
lit-check-cleared clean negative on the repeating-chain hypothesis — with a small
chance of a genuinely novel cislunar quasi-cycler row.

---

## 7. Confidence: buildable-from-digested-material vs needs-new-corpus

- **WSB surface (`core/wsb.py`): HIGH confidence, buildable now.** Def 3.11 /
  3.12, eq 3.9, eq 3.29, Def 3.22, the parabolic C = ±√2 golden, and the
  numerical stability-class algorithm are all in the Belbruno digest at
  equation/algorithm granularity.
- **BCT construction (`genome/bct_transfer.py`): MEDIUM confidence.** The
  backward two-arc recipe and the forward 2×2 method are digested at *recipe*
  granularity (the control/target variables, the heuristics, the ΔV bookkeeping)
  but the convergence-engineering papers ([33] Belbruno-Carrico 2000, [34]
  Belbruno-Humble-Coil 1997, [39] JGCD 1993) are **not in-repo**. The 2×2 Newton
  maps cleanly onto the existing BCR4BP corrector, so a competent build is very
  likely — but a named acquisition fallback ([39]) exists if it stalls.
- **Repeating-chain quasi-cycler search: LOW confidence of a positive result,
  HIGH confidence the search is well-posed.** The substrate can ask the question
  rigorously; whether nature answers "yes" is unknown and, per Theorem 3.58,
  probably "no clean periodic object" — hence the clean-negative default.

---

## 8. Prerequisite corpus acquisitions

**None hard-required to start** (WSB surface + backward BCT are buildable from the
digest). **Conditional acquisitions, named so the plan's Phase 0 can pull them on
a stall:**

1. **[39] Belbruno & Miller, JGCD 16(4):770-775, 1993** — the published Hiten BCT
   detail (operational targeting). Acquire if the forward 2×2 corrector won't
   converge from the digest (R3).
2. **Coherent QBCP α-tables (Gimeno-Jorba 2018 Table 4)** — already flagged as a
   documented future step in `core/bcr4bp.py`'s header. Acquire **only if** Phase 0
   R1 kill fires (incoherent model can't hold the Hiten apoapsis).
3. **Parker 2007 PhD thesis (low-energy transfers)** — operational depth on
   specific lunar missions; the digest already names it as "the next layer if
   #347 Phase 3 fires." Nice-to-have, not blocking.

Per the corpus-document policy, any acquired PDF gets OCR → digest → CORPUS_INDEX
before it's "done"; the plan's acquisition task carries that checklist.
