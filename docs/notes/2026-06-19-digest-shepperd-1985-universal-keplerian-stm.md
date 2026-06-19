# Digest — Shepperd 1985, Universal Keplerian State Transition Matrix

**Citation:** Shepperd, S.W. (1985), "Universal Keplerian State Transition Matrix," *Celestial Mechanics* **35**, 129–144. DOI/Bibcode: 1985CeMec..35..129S. (Charles Stark Draper Laboratory; supported by NASA under contract NAS9-16023.)

**Corpus file:** `/home/bruce/dev/cyclers_pdf/papers/shepperd-1985-universal-keplerian-state-transition-matrix-celest-mech-35.pdf` (16 pp; born-digital text layer, no OCR.)

**Status:** digested · text-layer · methodology paper (read fully).

---

## 1. What the paper is / contribution

A completely general, closed-form method to compute the **Keplerian (two-body) state transition matrix (STM)** in Goodyear's universal variables, valid for all conic types — elliptic, parabolic, hyperbolic — through a single unified algorithm (Abstract; §1 Introduction, p.129). It is a refinement/repackaging of Goodyear's 1965 *Astron. J.* 70 formulation [ref 5] and the 1966 NASA CR-522 report [ref 6], which Shepperd calls "the definitive works on the subject" (§1, p.129).

Shepperd's three specific enhancements over Goodyear (§1, p.129):
1. A new scheme for **Kepler's problem** (the necessary first step) using a new independent variable that needs only **one transcendental-function evaluation** per iteration.
2. That same transcendental function is chosen so **no additional complex function evaluations are needed** to build the STM itself.
3. The transcendental function is evaluated as a **Gaussian continued fraction** (top-down / Gautschi recurrence) for robust convergence.

The whole thing is unified so the Kepler solve and the STM assembly share intermediate quantities; the paper handles **multi-revolution** transfer-time intervals correctly (§5, p.134–135).

---

## 2. Core method

### 2.1 Universal Kepler propagation (Goodyear/Sundman base) — §2, p.130

Two-body ODE `d²r/dt² + (μ/r³) r = 0` (eq. 1). Sundman transformation `dt/dw = r` (eq. 2, ref [12]) gives a closed-form solution as a function of parameter `w`, universally applicable to all conic types.

Setup scalars (eqs. 3–5, p.130):
- `r₀ = |r₀|` (eq. 3)
- `ν₀ = r₀·v₀` (eq. 4)
- `β = 2μ/r₀ − (v₀·v₀)` (eq. 5) — **constant of motion, equals twice the negative orbital energy**; sign of β selects ellipse (β>0) / parabola (β=0) / hyperbola (β<0).

Range/time relations (eqs. 6–7):
- `r = r₀U₀ + ν₀U₁ + μU₂` (eq. 6)
- `t = r₀U₁ + ν₀U₂ + μU₃` (eq. 7) — this is **Kepler's equation** in universal form; for the STM problem `t` is given and we solve eq. 7 for `w`.

**f and g functions** (eqs. 8–9, p.130) — the Lagrange coefficients giving the propagated state:
- `r = f·r₀ + g·v₀`, with `f = 1 − μU₂/r₀`, `g = r₀U₁ + ν₀U₂` (eq. 8)
- `v = F·r₀ + G·v₀`, with `F = −μU₁/(r r₀)`, `G = 1 − μU₂/r` (eq. 9)

`w` relates to the eccentric-anomaly difference for ellipses by `w = (E − E₀)/√β` (eq. 10, p.131).

### 2.2 Universal functions Uₙ(w, β) — §4, p.132–133

The `Uₙ(w, β)` (Stumpff/universal functions [refs 10,11], popularized by Herrick [7], Pitkin [9], Battin [1]) defined by the series `Uₙ = Σ_{k≥0} (−β)^k w^{n+2k}/(n+2k)!` (eq. 20), with nested/Horner forms (eqs. 21–22). They reduce to trig/hyperbolic functions: `U₀ = cos(√β w)` / `cosh(√−β w)`, `U₁ = sin(√β w)/√β` / `sinh(√−β w)/√−β` (eqs. 23–26) — but the **trig forms are NOT used computationally** (round-off + singularities); the universal functions are used directly (text, p.133). Recurrence `Uₙ + β U_{n+2} = wⁿ/n!` (eq. 27). Derivative identities: `∂U₀/∂w = −βU₁` (eq. 28), `∂Uₙ/∂w = U_{n−1}` for n>0 (eq. 29), `∂Uₙ/∂β = ½[nU_{n+2} − wU_{n+1}]` (eq. 30). Note: do **not** climb to higher Uₙ via eq. 27 because it divides by β (zero/near-zero for parabolic) (p.133).

### 2.3 Recommended independent variable — §5, p.134

Iterating directly on `w` needs ≥2 transcendental evals and risks series-summation trouble for extreme hyperbolas. Shepperd instead iterates on
- `u = U₁(w/4, β) / U₀(w/4, β)` (eq. 31), bounded `−1/√|β| < u ≤ 1/√|β|`,

and an intermediate parameter
- `q = βu²/(1 + βu²)` (eq. 32), bounded `q ≤ ½`.

Low-order U's become simple in u, q (half-angle relations, eqs. 33–37): `U₀(w/2) = 1 − 2q`, `U₁(w/2) = 2(1−q)u`, then double-angle `U₀(w) = 2U₀²(w/2) − 1`, `U₁(w) = 2U₀(w/2)U₁(w/2)`, `U₂(w) = 2U₁²(w/2)`. **Newton–Raphson slope is simply `dt/du = 4(1−q)r` (eq. 43)** — a well-behaved monotonic function; the zero-slope trouble region `q=1` is excluded by construction (the u-range forces `q ≤ ½`).

**Multi-revolution handling (§5, p.134–135):** the raw computed t for ellipses lies in `2ν₀/β − P/2 < t ≤ 2ν₀/β + P/2` (eq. 38) with period `P = 2πμβ^{−3/2}` (eq. 39). To target a transfer time T outside that window, add n periods to the secular terms only: `U₃ ← U₃ + 2nπβ^{−3/2}` (eq. 40), `U ← U + 2nπβ^{−5/2}` (eq. 41), with `n = largest integer ≤ (1/P)[T + P/2 − 2ν₀/β]` (eq. 42).

### 2.4 The STM assembly — §3, p.131–132

The STM linearly propagates state deviations holding **transfer time t constant** (eq. 11, p.131):
```
[Δr; Δv] = Φ [Δr₀; Δv₀],   Φ = [[Φ11 Φ12],[Φ21 Φ22]] = [[∂r/∂r₀, ∂r/∂v₀],[∂v/∂r₀, ∂v/∂v₀]]   (eq. 12)
```
The four 3×3 sub-blocks (Goodyear form, eqs. 13–16, p.131):
- `Φ11 = f·I + [r v] · [[M21 M22],[M31 M32]] · [r₀ v₀]ᵀ` (eq. 13)
- `Φ12 = g·I + [r v] · [[M22 M23],[M32 M33]] · [r₀ v₀]ᵀ` (eq. 14)
- `Φ21 = F·I − [r v] · [[M11 M12],[M21 M22]] · [r₀ v₀]ᵀ` (eq. 15)
- `Φ22 = G·I − [r v] · [[M12 M13],[M22 M23]] · [r₀ v₀]ᵀ` (eq. 16)

Here `[r v]` is the 3×2 matrix of the propagated/initial position–velocity vectors, `I` the 3×3 identity, and the 2×2 coefficient blocks are **overlapping sub-matrices of a single 3×3 "M" matrix** (eq. 17, p.132; reprinted as eq. A.41, p.143). M's entries are built entirely from `f, F, g, G, r, r₀`, the universal functions `U₀, U₁, U₂`, and the scalar `W`:
- `W = gU₂ + 3μU` (eq. 18)
- `U = (U₂U₃ + wU₄ − 3U₅)/3` (eq. 19)

**Key structural point (p.132):** beyond the Kepler solution, only the **nine M coefficients** are additional, and `W`/`U` are the only new quantities — all other M entries come from the Kepler solve (eqs. 3–9). So the full 6×6 STM costs essentially one extra transcendental function `U` on top of Kepler. M has no physical interpretation; it is purely a compact bookkeeping device.

### 2.5 The two transcendentals as continued fractions — §6–7, p.135–140

`U₃` and `U` are the two transcendental functions; both expressed via the **Gauss continued fraction** `G(a,b,c,x) = F(a,b+1,c+1,x)/F(a,b,c,x)` of ratios of hypergeometric series (eqs. 44–45, 68–69):
- `U₃(w,β) = (4/3) U₁³(w/2,β) · G₃(q)`, with `G₃(q) = G(3,0,3/2,q)` (eqs. 46–47, p.136; Gauss [ref 3], universal-function generalization due to Battin [ref 1]).
- `U(w,β) = (16/15) U₁⁵(w/2,β) · G₅(q)`, with `G₅(q) = F(5,1,7/2,q) = G(5,0,5/2,q)` (eqs. 48–59, p.136–137).

Going from hypergeometric series to continued fraction widens convergence from `|q|<1` to `q<1` (Wall [13]); and the deliberate `q ≤ ½` range avoids the `0·∞` indeterminacy near q→1 and the slow-convergence region q>½ (p.136). **Important asymmetry:** compute `U₃` from `U` and not the reverse — `U₃ = βU + (1/3)U₁U₂` (eq. 60, p.137); the reverse divides by β (zero for parabolas, p.137). Since `U` is not needed inside the time loop, it can be evaluated once on the converged q (p.137).

The continued fraction is evaluated by **Gautschi's top-down recurrence** (§7, eqs. 61–67, p.138), reformatted for digital computers as the integer-coefficient recursion eqs. 72–85 (p.139). For `G₃`/`G₅` the params `k, ℓ, d, n` are all integers (p.139). Caution: isolated zero-denominator convergents occur only for q>½ (e.g. G₃ at q=5/6, 7/8; G₅ at q=7/10, 3/4) and are never reached given the `q ≤ ½` range (p.140).

### 2.6 Appendix A — the runnable algorithm (p.140–143)

Inputs: transfer time T, initial state `(r₀, v₀)`. Order: constants `r₀, ν₀, β` (A.1–3) → init `u=0`, `ΔU=0`, period/multi-rev terms for β>0 (A.4–8) → **Kepler iteration loop** computing `q, U₀(w/2), U₁(w/2), U, U₀, U₁, U₂, U₃, r, t` and Newton step `u ← u − (t−T)/(4(1−q)r)` (A.9–19; G via the continued fraction A.20–33 with seed `k=−9, ℓ=3, d=15, n=0`) → **Kepler solution** `f, g, F, G, r, v` (A.34–39) → **transition matrix** `W = gU₂ + 3μU` (A.40), M (A.41), and Φ11..Φ22 / Φ (A.42–46). Outputs: extrapolated `(r, v)` and the STM Φ. (Note: the Appendix omits convergence-criteria/boundary-keeping detail by the author's own statement, p.140.)

---

## 3. Why it matters to THIS project

This is the **canonical closed-form two-body STM** used inside differential-correction / shooting inner loops. Wherever a corrector needs `∂(final state)/∂(initial state)` for a Keplerian arc, Shepperd's Φ gives it analytically — no numerical variational-equation integration, no finite-difference STM. It is the universal-variables alternative to integrating the variational equations alongside the trajectory.

Relevant cyclerfinder touch-points:
- **Lambert / multi-arc Lambert legs** (`src/cyclerfinder/core/lambert.py`, `src/cyclerfinder/search/dsm_leg.py`): a Keplerian-arc corrector that targets arrival position/velocity can use Shepperd Φ for the sensitivity blocks instead of finite differences — faster and exact, and conic-agnostic (the same code path serves hyperbolic departure/arrival v∞ arcs and the elliptic transfer ellipse). The S1L1 multi-arc closure work (two generic-return arcs, per MEMORY) is exactly the kind of multi-segment shooting problem where an analytic per-arc STM chains cleanly.
- **Multi-revolution targeting:** eqs. 38–42 give the correct period bookkeeping for elliptic legs spanning >1 rev — directly applicable to resonant S/L Earth-to-Earth intervals.
- **Relation to Pellegrini–Russell 2016 (task #372):** Pellegrini–Russell give a *fixed-path* STM (variational along an already-computed trajectory). Shepperd is the *closed-form* two-body STM — the natural cheap baseline for the Keplerian (patched-conic) layer of a corrector, before any CR3BP/ER3BP refinement where the analytic two-body STM no longer holds and a Pellegrini–Russell-style numerical STM is required. The two are complementary: Shepperd for the conic legs, numerical/variational STM for the multi-body arcs.

---

## 4. Limits / caveats

- **Two-body only.** This is a Keplerian (inverse-square, single attractor) STM. It is NOT a CR3BP/ER3BP/BCR4BP STM — no third-body, no rotating-frame, no perturbation terms. For cyclerfinder's multi-body corrector arcs it serves only the patched-conic / Lambert layer; multi-body arcs still need numerical variational STMs.
- **Singularity discipline is the whole point but must be respected:** trig forms of U₀/U₁ (eqs. 23–26) are visualization-only — using them reintroduces round-off and singularities (p.133). Never compute higher Uₙ via eq. 27 (÷β) or U from U₃ (÷β); both blow up at the parabolic limit β=0 (p.133, p.137).
- The `q ≤ ½` independent-variable range is load-bearing: it avoids the Newton zero-slope point (q=1, eq. 32/43) and the continued-fraction zero-denominator points (q>½, p.140). An implementation that strays outside this range loses those guarantees.
- **Appendix A is a skeleton.** The author explicitly omits convergence-criteria and boundary-keeping ("keeping the Newton–Raphson adjustment from stepping outside the boundaries of the problem"), which are "crucial to turning a set of equations into a reliable working computer program" (p.140). The note's eqs. are a first step toward code, not production-ready.
- Does NOT cover: non-time-constrained STMs (constraints other than time are only *conjectured* to also need just 9 extra coefficients, p.132), perturbed motion, or partials w.r.t. μ.
