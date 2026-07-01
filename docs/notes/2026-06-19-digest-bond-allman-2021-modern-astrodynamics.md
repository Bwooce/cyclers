# Digest — Bond & Allman, *Modern Astrodynamics: Fundamentals and Perturbation Methods*

**Full citation:** Bond, V.R. & Allman, M.C., *Modern Astrodynamics: Fundamentals and Perturbation Methods*, Princeton University Press, 1996 (2021 reissue), DOI 10.1515/9780691223902. ISBN 0-691-04459-7. 264 pp. (LoC TL1050.B66 1996). First author Victor R. Bond (NASA-JSC Mission Planning & Analysis Division); developed from graduate courses at the University of Houston–Clear Lake.

**Source:** doi:10.1515/9780691223902 (born-digital, has a text layer — no OCR needed).

**Scope:** chapter-summary (TOC + 5 deep chapters/appendices). Deep-read: **Ch 5** (f/g universal-variable, pp. 58–80), **Ch 6** (Two-point BVP / Lambert, pp. 81–90), **Ch 8** (Perturbation theory / VOP, pp. 117–146), **Ch 9** (Special perturbation methods: regularization / Sundman / Sperling-Burdet, pp. 147–183), and **Appendix E** (Stumpff functions, pp. 236–238). Sampled: Ch 1–4, 7, 10, 11, App A–D, F. The book's intellectual center of gravity is Bond's own JSC lineage (Sperling 1961 → Burdet 1968 → Bond-Fraietta, Bond-Gottlieb), so the regularization/Sperling-Burdet material (Ch 9) is the unique contribution and got the deepest read.

---

## Full chapter list (from TOC, pp. v–viii)

**Part I — Fundamentals**
1. **Background** (p. 3) — notation, units, time systems (UT/ephemeris time). Boilerplate.
2. **The Two-Body Problem** (p. 12) — Newton's laws; the six classical integrals (area **c**, energy *h*, Laplacian **P**); the abort problem; Kepler's equation derived. Standard but careful integral-by-integral treatment.
3. **Kepler's Laws** (p. 31) — three laws + planetary-mass determination. Standard.
4. **Methods of Computation** (p. 40) — position/velocity from integrals, from true & eccentric anomaly; Kepler-equation solution (Algorithm 1, 2); orbital elements incl. Delaunay & Poincaré sets. Standard.
5. **The f and g Functions** (p. 58) — *[DEEP]* f/g in ΔE and in the **universal variable** x (Battin-style S(z), C(z)); Kepler's equation solution (Alg. 3, 4, 5); f/g as Taylor series in time with Moulton convergence radius.
6. **Two-Point Boundary Value Problems** (p. 81) — *[DEEP]* Lambert's problem (eccentric-anomaly form and universal-variable form, Alg. 6); the Linear Terminal Velocity Constraint (LTVCON) targeting problem.
7. **Applications** (p. 91) — interplanetary trajectories (heliocentric + planetocentric phases, V-infinity Alg. 7, gravity-turn flyby); shuttle ascent; anytime-deorbit; relative motion (Hill/CW-type). Sampled — practical, JSC-flavored.
8. **Perturbation Theory** (p. 117) — *[DEEP]* general VOP framework; Poisson's method; Lagrange VOP; perturbed two-body element rates (da, dc, di, dΩ, Laplace vector) via Poisson's method.
9. **Special Perturbation Methods** (p. 147) — *[DEEP]* propagation-error growth of Cowell/Encke; **regularization** (Sundman transform, embedding); the **Sperling-Burdet** linearized/regularized element formulation; numerical examples (oblate Earth+Moon, **stable libration point**, continuous radial thrust).
10. **Runge-Kutta Methods** (p. 184) — error classification; fixed-step RK1/RK2; variable-step RKF, RKF4(5) (Alg. 8) with Fehlberg coefficients. Sampled — used as the integrator for the Ch 9 examples.
11. **Types of Perturbations** (p. 198) — *third-body* perturbations (n-body, Battin's f(q) reformulation for the interior problem, Legendre expansion for the exterior problem); planetary geopotential (Legendre/zonal-sectorial-tesseral, recursive nonsingular Gottlieb/Mueller algorithm). Sampled.

**Appendices:** A Coordinate transforms; B Hyperbolic motion; C Conic sections; D Transfer-angle resolution; **E Stumpff functions** *[DEEP]*; F Orbit geometry (ĉ, P̂ as functions of Ω,i,ω). References p. 243; Index p. 247.

---

## Deep-read chapter summaries

### Chapter 5 — The f and g Functions (universal variable) — pp. 58–80

**Content.** §5.3 develops f/g in change-of-eccentric-anomaly ΔE (Alg. 3, pp. 70–71), needing only the semi-major axis *a* (no initial *E*). §5.4 (pp. 71–77) is the **universal-variable** formulation following Battin: defines the truncated-series functions S(z) (eq. 5.49) and C(z) (eq. 5.50), with closed forms `S(z)=(√z−sin√z)/(√z)³`, `C(z)=(1−cos√z)/z` for z>0 and hyperbolic analogues for z<0 (p. 71–72). Crucially it states **S(z)=c₃(z), C(z)=c₂(z)** — the universal functions are Stumpff functions (p. 72). Auxiliary variable `ΔE=√α₀·x`, with `α₀=1/a` (>0 elliptic, <0 hyperbolic, =0 parabolic) (p. 72). Final universal f/g set (p. 74, eq. 5.53–5.56):
- `f = 1 − (x²/r₀)·C(α₀x²)`
- `g = t − t₀ − (1/√μ)·x³·S(α₀x²)`
- `ḟ = (√μ/(r·r₀))·[α₀x³S(α₀x²) − x]`
- `ġ = 1 − (x²/r)·C(α₀x²)`

These are shown to satisfy `r̈ + (μ/r³)r = 0` directly with **no assumption on the sign of α₀** (p. 76) — one solver for ellipse/parabola/hyperbola. Universal Kepler equation (eq. 5.57) and the iteration are Alg. 4 (pp. 76–77) and Alg. 5 (p. 77). §5.5 gives f/g as Taylor series in time (Bond [13]) with recurrence relations (eq. 5.69–5.72, p. 79) and Moulton's convergence-radius ρ formulas for ellipse/parabola/hyperbola (p. 80).

**Relevance to cyclerfinder.** This is the canonical reference for our universal-variable Kepler propagator and corrector front-end. The explicit α₀-sign-agnostic f/g set (eq. 5.53–5.56) is exactly what a single-code-path two-body propagator should implement; useful as a same-model golden cross-check for `core/` two-body routines. The S(z)=c₃, C(z)=c₂ identification (p. 72) ties our Stumpff/universal code to the Sundman material in Ch 9 and the **tulip-orbit Sundman work** — same function family, one is the BVP face and the other the regularized-propagation face.

### Chapter 6 — Two-Point Boundary Value Problems / Lambert — pp. 81–90

**Content.** §6.2 (pp. 81–82) sets up Lambert: given r₀, r, and time-of-flight Δt, find ṙ₀ (then ṙ at any t via Ch 4). Elliptic eccentric-anomaly form (eq. 6.1–6.4) gives three equations in three unknowns (a, p, ΔE) iterated on ΔE (Alg. on p. 83). §6.3 (pp. 83–88) is the **universal-variable Lambert** following Battin (who credits Deyst/MIT IL); BMW call the same thing the "Gauss problem." Defines `y = (r·r₀/p)(1−cosΔφ)` (eq. 6.12), constant `A = sinΔφ·√(r·r₀/(1−cosΔφ))` (eq. 6.13), `z=α₀x²`; iteration eq. 6.14–6.19 with the clean final f/g/ġ forms `f=1−y/r₀`, `g=A√(y/μ)`, `ġ=1−y/r` (eq. 6.20–6.22) and `ṙ₀=(r−r₀f)/g` (eq. 6.23). Full procedure = **Algorithm 6** (pp. 86–88); transfer-angle Δφ branch selection deferred to Appendix D. §6.4 (pp. 88–90) is **LTVCON** (Linear Terminal Velocity CONstraint): an unusual, JSC-internal targeting algorithm (Dave Long / Lineberry / Shepperd, CSDL/JSC memos 1974) that replaces the time-of-flight constraint with a *linear* constraint between terminal radial and horizontal velocity components `ṙ = c₁ + c₂·v_H`, reducing to a quadratic `A·v_H² − 2B·v_H − C = 0` (eq. 6.32).

**Relevance to cyclerfinder.** Direct reference for `core/lambert.py` (currently modified in working tree). The universal-variable single-revolution Lambert (Alg. 6) is the standard our Lambert legs should reproduce; the eccentric-anomaly cross-form (§6.2) gives an independent same-model check. The Δφ-branch caution (Appendix D) is exactly the kind of multi-revolution / transfer-direction subtlety relevant to the **DSM-leg search** (`search/dsm_leg.py`). LTVCON is a niche targeting variant — not currently in our pipeline, but flagged as a possible alternative terminal-constraint formulation if a future leg needs a velocity-component constraint rather than a fixed TOF (loose relevance to the **#380 BVP integral-constraint corrector**, which similarly swaps the natural Lambert constraint for a different closure condition).

### Chapter 8 — Perturbation Theory (VOP framework) — pp. 117–146

**Content.** §8.1 (pp. 117–118) frames the three numerical options for `r̈+(μ/r³)r=F`: **Cowell** (integrate Cartesian directly — simple, diverges fast), **Encke** (integrate the deviation η=r−ρ from a reference two-body solution, needs periodic *rectification*), and the smarter route — pick new dependent variables (the integrals) and derive their ODEs = **variation of parameters**. §8.2 derives **Poisson's method**: for integrals σ_k of the unperturbed problem, along a perturbed solution `σ̇_k = (∂σ_k/∂ẋ)·g` (eq. 8.21). §8.3 derives the **Lagrange VOP**: `ċ = [∂x/∂c]⁻¹ g` (eq. 8.33), requiring the Jacobian `∂x/∂c` to be invertible. Both verified on the harmonic oscillator. §8.4 (p. 126) summarizes the two-body integrals (area, energy, Laplacian, time of pericenter) with the two dependency relations `c·ε=0` and `p=c²/μ`. §8.6 (pp. 127–134) applies Poisson's method to derive element rates under a perturbation F decomposed in the (r̂, φ̂, ĉ) frame:
- `ȧ = (2/μ)a²ṙ·F` (eq. 8.39) — energy/SMA affected only by along-velocity F.
- `ċ = r×F` (eq. 8.41); `ċ = r·φ̂·F` (magnitude, eq. 8.44) — angular momentum affected only by the *horizontal* component.
- `di/dt = (r/c)cos(ω+φ)·ĉ·F` (eq. 8.49) — inclination affected only by the **normal** component.
- `dΩ/dt = (r/(c·sin i))·sin(ω+φ)·ĉ·F` (eq. 8.51) — node also only by the normal component.
- Laplace-vector rate `μ ε̇ = 2(ṙ·F)r − (r·ṙ)F − (r·F)ṙ` (eq. 8.54).

§8.7 works partially-solved cases (conservative potentials, oblate planet, **time-dependent potential giving the Jacobi-integral constant** `h + V − ω·(r×ṙ) = const` when `∂V/∂t = −ω·r×∂V/∂r`, §8.7.3 — this seeds the Ch 9 Jacobi-element work; tethered satellite; drag).

**Relevance to cyclerfinder.** The Poisson/Lagrange VOP machinery is the textbook backing for any element-rate or averaged-dynamics reasoning we do on perturbed legs. The §8.7.3 result that a *rotating-frame* time-dependent third-body potential admits the **Jacobi integral** is the bridge to CR3BP work (`core/cr3bp.py`, modified in tree) — Bond literally builds the Jacobi constant as a VOP element. The clean "which F-component moves which element" rules (eq. 8.39/8.44/8.49) are a sanity-check lens for **#314 heteroclinic transport** intuition (normal vs in-plane forcing). The variational-equation framing also underpins **#372 STM**: Ch 9 (below) gives the explicit STM/propagation-error matrix.

### Chapter 9 — Special Perturbation Methods: regularization & Sperling-Burdet — pp. 147–183

This is the book's signature chapter and the deepest dive.

**§9.1 Propagation error (pp. 147–148).** The error between two neighboring Cowell solutions obeys the **variational equation** `d(PE)/dt = A·PE` (the STM equation), with `PE = (δr, δṙ)` and, for F=0,
```
A = [ 0   I ]        C = −(μ/r³)[ I − r̂ r̂ᵀ ]
    [ C   0 ]
```
(p. 147–148). The solution `PE(Δt)=e^{AΔt}PE(0)`; because **one eigenvalue of A is real and positive, √(2μ/r³)**, the Cowell/Encke propagation error *always grows even with no perturbation* (p. 148). This is the explicit motivation for regularized elements that have *zero* propagation-error growth when F=0.

**§9.2–9.3 Regularization (pp. 148–157).** Following Sperling (1961) and Burdet (1968). Regularization = removing singularities from the *differential equations* (not their solution). Two ingredients:
1. **Sundman transformation** `dt/ds = r` (eq. 9.12) — change independent variable from time *t* to fictitious time *s*, taming the r→0 singularity.
2. **Embedding** the Laplace vector and Keplerian energy as constants into the equation.

Result (eq. 9.21): `r″ + α_t·r = −με + r²F`, where `α_t = −2h` and `(·)′ = d/ds`. This is a **forced harmonic oscillator** in fictitious time — linearized and regularized. Companion scalar distance equation `r″ + α_t·r = μ + r·r·F` (eq. 9.24). §9.3.1–9.3.2 swaps the Keplerian energy for the **Jacobi integral** `α_J = α_t + 2V(r,t) − 2ω·c` (eq. on p. 155) following the Bond-Gottlieb 1989 modification — directly relevant since the Jacobi integral is the CR3BP conserved quantity. The two-body solution is written in **Stumpff functions** `r = r₀c₀ + r₀′ s c₁ − με s²c₂` (eq. 9.33), with `t = t₀ + r₀sc₁ + r₀′s²c₂ + μs³c₃` being **Kepler's equation in Stumpff form** (eq. 9.41).

**§9.3.4–9.3.5 The δ and γ elements (pp. 158–160).** Because the Laplace vector is undefined at e=0, replace it with `δ = −(α_t r₀ + με)` (eq. 9.42) — δ is the initial value of the regularized state `x=−α_t r − με`. Likewise `γ = μ − α_t r₀` (eq. 9.45) is the initial value of the scalar `y=μ−α_t r`. These are non-singular at zero eccentricity.

**§9.4–9.6 Sperling-Burdet element system (pp. 160–177).** The full VOP treatment yields a **system of 15 (effectively 14, via the integral γ+aα_J=μ) first-order ODEs** in fictitious time for the spatial elements (α=r₀, β=r₀′, δ), temporal elements (r₀, r₀′, t₀, γ), and the Jacobian/axial elements (α_J, σ=ω·c). Full summary on p. 177. **Key property (p. 178): when perturbations are zero (F=0, Q=0, P=0) the system reduces to z′=0 — the elements are exactly constant, with zero propagation-error growth** — the opposite of Cowell/Encke.

**§9.7 Numerical examples (pp. 178–183) — sourced validation data:**
- **Oblate Earth + Moon** (e=0.95, J₂ + idealized lunar third-body), Table 9.2: at ~50 revolutions (288.13 days) Sperling-Burdet RSS error **0.318 km at 62 steps/rev** vs Cowell 42.5 km at 240 steps/rev — two orders of magnitude better per-step. Constants: μ=398601 km³/s², a_e=6371.22 km, J₂=1.08265e-3, GM_moon=4902.66 km³/s², ρ_EM=384400 km, Ω_moon=2.665315780887e-6 rad/s.
- **Stable libration point in an idealized Earth-Moon system** (Table 9.3): satellite placed at the **L4/L5 triangular Lagrangian point** (60° ahead, same distance and circular rate as the Moon) treated as a *perturbed two-body* problem in an inertial frame; ω²=(GE+GM)/ρ³, T=27.28459145 days; after 1000 lunar periods position error stays at the **metre level (Δx₁ ≤ 5.6 m, Δx₂ ≤ 2.9 m)**. r₀ = ½ρ(î+√3 ĵ), ṙ₀ = ½ωρ(ĵ−√3 î). Because the libration point is stable there is no propagation-error growth, so residual error is *pure truncation*.
- **Continuous radial thrust** (Tsien's analytic solution, Table 9.4): circular orbit r₀=6800 km thrust to escape; at t=12000 s Sperling-Burdet matches the analytic e=1 escape radius 30682.724 km to **0.13 m** at ~105 s/step.

**Relevance to cyclerfinder.** This chapter is the strongest single justification for our **tulip-orbit Sundman work**: Bond is the authoritative source for `dt/ds=r` regularization, the Stumpff-function solution form, and *why* it beats Cowell (zero propagation-error growth, eq. on p. 148 / p. 178). The **L4/L5 stable-libration-point example (Table 9.3)** is a near-perfect sourced cross-check target: it is exactly a CR3BP triangular-point trajectory recast as perturbed two-body, with published 15-digit-style constants and metre-level published residuals — a candidate same-model golden target for `core/cr3bp.py` and a sanity check against the Earth-Moon mu-double-count bug noted in recent commits. The **Jacobi-integral-as-element** development (§9.3.1, Bond-Gottlieb 1989) is the cleanest textbook link between perturbed-two-body machinery and the CR3BP integral we use. The propagation-error variational matrix (§9.1) is a literal statement of the **#372 STM**: `A=[[0,I],[C,0]]`, `C=−(μ/r³)(I−r̂r̂ᵀ)`. The regularized BVP form `r″+α_t r = −με + r²F` (eq. 9.21) is conceptually adjacent to the **#380 BVP integral-constraint corrector** — it shows how embedding integral constants converts a singular BVP into a well-conditioned oscillator BVP, the same "fold the constraint into the equation" move #380 makes.

### Appendix E — Stumpff Functions — pp. 236–238

**Content.** Defines the Stumpff functions by the single series `c_n(z) = Σ_{k≥0} (−1)^k z^k/(2k+n)!` (p. 236), with z=α_J s². One series represents trig (z>0), hyperbolic (z<0), and parabolic (z=0, c_n(0)=1/n!) cases. Closed forms for c₀, c₁, c₂ (p. 236). Core recurrence `c_n(z) + z·c_{n+2}(z) = 1/n!` (eq. E.1). Derivative formulas in fictitious time `s·c_n′ = c_{n−1} − n·c_n` (eq. E.4) and `c_n′ = α_J·s(n c_{n+2} − c_{n+1})` (eq. E.5). Integration identity `∫ s^k c_k(ρs²) ds = s^{k+1} c_{k+1}(ρs²)` (eq. E.6–E.9). Identities E.10–E.14 (e.g. `c₀²+zc₁²=1`).

**Relevance to cyclerfinder.** These are the exact recurrence/derivative/integral identities the Sperling-Burdet system (Ch 9) and the universal-variable propagator (Ch 5) are built on. Any Stumpff implementation in our `core/` (universal Kepler, Sundman/tulip propagation) should match E.1/E.4/E.5/E.6 verbatim — these are the canonical golden identities for unit tests of the Stumpff layer.

---

## What this book uniquely adds vs our existing digests

Against the digests we already hold — **Szebehely 1967** (CR3BP theory/equilibria/zero-velocity), **Gurfil 2007** (modern astrodynamics methods survey), **Parker 2007** (low-energy/manifold transfers), **Belbruno 2004** (WSB/ballistic capture) — Bond-Allman adds:

1. **The definitive regularization + Sperling-Burdet treatment (Ch 9).** None of our existing digests give the Sundman-transform/embedding derivation, the explicit δ/γ non-singular elements, the 15-element Sperling-Burdet ODE system, or the *proof that regularized elements have zero propagation-error growth when F=0* (§9.1). This is Bond's home turf (Sperling 1961 → Burdet 1968 → Bond-Fraietta → Bond-Gottlieb 1989) and is unmatched elsewhere in the corpus. Directly feeds the **tulip Sundman work**.
2. **Stumpff functions as a first-class, self-contained reference (App E + §5.4).** Szebehely/Parker assume universal variables; Bond derives the series, recurrences, derivatives, and the S=c₃/C=c₂ identity from scratch — a citable golden-identity source for our Stumpff/universal code.
3. **A CR3BP triangular-libration-point trajectory validated as perturbed-two-body with published metre-level residuals (Table 9.3).** Szebehely gives the *theory* of L4/L5; Bond gives a *numerically validated propagation* with sourced constants — a far better same-model golden target than abstract theory.
4. **The Poisson/Lagrange VOP element-rate "which-component-moves-which-element" rules (Ch 8)** and the **propagation-error variational matrix (§9.1)** as explicit STM source — more pedagogically explicit than Gurfil's survey.
5. **Battin's f(q) third-body reformulation for the interior problem and the recursive nonsingular geopotential algorithm (Ch 11)** — practical computational nonsingular forms not in the CR3BP-focused digests.

What it does **not** add: no dynamical-systems-theory manifold/invariant-set material (that is Szebehely/Parker/Belbruno territory), no low-energy/WSB transfer design, no Floquet/periodic-orbit families. It is a *two-body + perturbation-of-two-body* book; CR3BP appears only as a perturbed-two-body example (the L4/L5 case), never as the primary model. So it complements rather than overlaps our CR3BP/manifold digests.

---

## Sourced-discipline note

Every numeric/equation claim above carries a page or equation citation. The Table 9.2/9.3/9.4 values and the §9.7 model constants are published (Bond-Allman pp. 179–183, themselves citing Stiefel-Scheifele [64] and Tsien); they are admissible as EXPECTED-side golden targets only insofar as they are the *authors'* published outputs, not values our code produced. The L4/L5 case (Table 9.3) is the most promising candidate for a same-model golden cross-check of `core/cr3bp.py`, subject to the standard orbit-closure discipline (independent cross-check, verify topology vs source, hold writeback till confirmed).
