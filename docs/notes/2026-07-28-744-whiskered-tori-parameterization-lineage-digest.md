# #744 — Fontich/de la Llave/Sire, Jorba/González/de la Llave/Villanueva,
# Huguet/de la Llave/Sire, Capiński/Gidea/de la Llave, Mireles James/Murray,
# Pérez-Palau/Masdemont/Gómez, Fernández-Mora/Haro/Mondelo: `#730` §5 items
# 32-35, 37-39 (whiskered-tori/parameterization-method lineage, batch 2)

**Task:** `#744`. Seven PDFs filed in the private corpus this session by the
coordinating session, all confirmed text-layer (no OCR needed — page-by-page
`pdftotext` char-count verified far above the 10-char/page floor by the
coordinating session before dispatch; re-confirmed here only by inspecting
`pdftotext -layout` output quality, not re-run from scratch). This is the
second `#730` §5 batch in the Haro/de la Llave/Fontich/Sire/Huguet/Capiński/
Gidea whiskered-tori/parameterization-method lineage; the first batch
(`#741`, items 6-10 of §2's separate top-10 ranking) covered Iuliano,
Calleja/del-Castillo-Negrete/Martínez-del-Río/Olvera, Calleja & de la Llave,
Cabré/Fontich/de la Llave (Part III), and Gonzalez & Mireles James — read in
full for house style before writing this digest (structure only, not
content, per the coordinating session's own framing). `#730` §5 item 41
(Frauenfelder & Moreno 2023, GIT-quotients symplectic-bifurcation paper) is
explicitly **not** in this task's scope — a separate agent is handling it in
parallel; not touched here.

## 0. Acquisition, filing, text-layer check

All seven filed in the private `cyclers_pdf` corpus
(`/Users/bruce/dev/cyclers_pdf/papers/`, separate repo, never committed to
the public `cyclers` repo):

| # | Paper | Filename | Pages | Text layer |
|---|---|---|---|---|
| 1 | Fontich, de la Llave & Sire 2009, JDE 246(8) | `fontich-delallave-sire-2009-...-arxiv-0903.0311.pdf` | 77 | Native, clean (`pdftotext -layout` extracted cleanly through the full document; one recoverable `xref` warning from `pdftotext`, auto-reconstructed, no content loss) |
| 2 | Jorba, González, de la Llave & Villanueva 2005 (short comm.) | `jorba-gonzalez-delallave-villanueva-2005-...-upcommons.pdf` | 8 | Native, clean |
| 3 | Huguet, de la Llave & Sire 2012, DCDS-A 32(4) | `huguet-delallave-sire-2012-...-arxiv-1004.5231.pdf` | 58 | Native, clean |
| 4 | Capiński, Gidea & de la Llave 2016, Nonlinearity 30(1) | `capinski-gidea-delallave-2016-...-arxiv-1510.00591.pdf` | 34 | Native, clean |
| 5 | Mireles James & Murray 2017, IJBC 27(14) | `mirelesjames-murray-2017-...-arxiv-1706.03345.pdf` | 36 | Native, clean |
| 6 | Pérez-Palau, Masdemont & Gómez 2015, CMDA 123(3) | `perezpalau-masdemont-gomez-2015-...-preprint.pdf` | 25 | Native, clean |
| 7 | Fernández-Mora, Haro & Mondelo 2024, SIADS 23(1) | `fernandezmora-haro-mondelo-2024-...-preprint.pdf` | 45 | Native, clean |

All seven were read in substantial detail (introduction, definitions,
theorem statements, worked examples, and reference lists); the three
longest and most proof-heavy (papers 1 and 3, both Fontich/de la Llave/Sire
& Huguet/de la Llave/Sire theory-and-algorithm texts) had their dense
mid-document estimate/convergence-proof sections sampled rather than read
line-by-line, per the same "chapter-summary scope" convention `#741` used.

---

## 1. Fontich, de la Llave & Sire 2009 — "Construction of Invariant Whiskered
Tori by a Parameterization Method. Part I: Maps and Flows in Finite
Dimensions," *J. Differential Equations* 246(8):3136–3213

DOI `10.1016/j.jde.2009.01.037` (CrossRef-confirmed per the task brief).
arXiv:0903.0311.

### 1.1 The a posteriori theorem

**Theorem 3.11** (the paper's main result, maps case): given a Diophantine
frequency ω ∈ D(κ,ν), an exact symplectic analytic map F, and an
**approximate** embedding K₀ satisfying two non-degeneracy conditions
(spectral condition of Def. 3.4, non-singular averages of Q₀/A₀), if the
initial invariance error E₀ = F∘K₀ − K₀∘T_ω is small enough in two explicit
norm bounds (Eq. 18), then there exists a **true** invariant torus K∞ nearby,
with an explicit error bound ‖K∞−K₀‖ ≤ Cκ²δ⁻²ᵛ‖E₀‖. This is the *a
posteriori* format: no closeness to an integrable system, no action-angle
coordinates, and no sequence-of-canonical-transformations proof — instead a
Newton-type iteration adding successive corrections. Theorem 3.14 gives
local uniqueness (up to translation on the torus); Theorem 3.15 deduces the
flow case from the map case via a time-one-map argument.

**The vanishing lemma** (§8.8, Lemma 8.26 for flows; the analogous Lemma 4.9
is used for maps) is the technical device that removes an auxiliary
parameter λ introduced into the invariance equation to make the KAM
iteration solvable on the center subspace. It bounds |λ−λ*| ≤ C‖E‖ᵨ by a
Cartan-formula argument on the exact-symplectic vector field, so that once
the iteration converges (E→0) one recovers λ∞=0 and a genuine (not merely
"translated") invariant torus of the *original*, unmodified system. This is
the mechanism the task brief flagged as load-bearing for the already-corpus
Kumar/Anderson/de la Llave 2022 CMDA whiskered-tori paper's own convergence
argument.

**Primary vs. secondary tori** (§7.4, directly answering the task's
extraction request): primary KAM tori are homotopic to Tˡ×{0}, i.e.
non-contractible, and are the tori persisting from the integrable system.
Secondary tori are homotopic to a *lower*-dimensional torus Tˡ⁻ᵏ×{0} — they
are contractible, are **not present in the integrable system at all**
(generated by the perturbation itself, near resonances — "islands" in 2D
maps), and their existence is not covered by standard perturbative KAM
theory. The paper's method handles both uniformly because it never uses
action-angle coordinates, where primary/secondary tori would require
genuinely different treatments. Proposition 7.2 (bootstrap of regularity)
separately shows that a Cᵣ invariant torus (r > 4ν) of an analytic map is
automatically analytic — licensing the "these tori aren't just numerically
smooth, they're actually analytic" reading used downstream.

### 1.2 Cross-check: consistent with `#741`'s lineage tracing, no new gap

No direct citation of this paper in `src/cyclerfinder/`; `genome/qp_tori.py`
remains on the separate Olikara & Scheeres GMOS lineage (per `#741`'s own
already-established finding). This paper's own reference list is almost
entirely classical KAM/small-divisor background (Zehnder, Moser, Rüssmann,
Kolmogorov-era) already covered by the project's own low-priority-background
convention (§8.3 below); no numeric worked examples of its own (purely
theoretical, consistent with it being the rigorous-existence half of a pair
completed algorithmically by paper 3 below).

---

## 2. Jorba, González, de la Llave & Villanueva 2005 (short communication)
— "KAM Theory Without Action-Angle Coordinates"

**Version caveat (per the task brief):** this 8-page document is a short
communication/preprint, **not** the full 41-page published paper (de la
Llave, González, Jorba & Villanueva, "KAM theory without action-angle
**variables**" [note: published title says "variables," this document's
title says "coordinates" — same work], *Nonlinearity* 18(2):855–895 (2005),
DOI `10.1088/0951-7715/18/2/020`). The document itself explicitly defers to
"full results will be available soon in [6]" (ref [6] is the in-preparation
full paper). Theorem numbering here (just "Theorem 1," "Lemma 1") is
this short document's own, self-contained numbering — **not** a subset of
the full paper's numbering, and the depth here (one core theorem, one
supporting lemma, a list of seven unproven "generalizations... work in
progress") is materially abbreviated relative to what a 41-page journal
paper would contain. Treat this digest's coverage as bounded to what this
short version actually contains.

### 2.1 The core idea (present in this short version)

Classical KAM requires four technicalities the authors want to avoid: (i)
action-angle coordinates, (ii) a perturbative (near-integrable) setting,
(iii) Lagrangian tori, (iv) transformation theory (sequences of canonical
transformations). The paper's replacement: represent an approximately
invariant torus directly as a **parameterization** T(θ) in the system's own
(e.g. Cartesian) coordinates — not necessarily Lagrangian — and define a
"quasi-torus" via the invariance-defect equation L_ω T(θ) = J∇H(T(θ)) +
R(θ). **Lemma 1** is the key structural fact: the antisymmetric matrix
S = Tᵀ J T (measuring failure of T to be Lagrangian) is bounded by
‖S‖ ≤ c‖R‖/(δ^{ν+1}) under a Diophantine condition on ω — i.e. a quasi-torus
with small residual R is automatically *close to Lagrangian*, which lets the
authors decompose the Newton correction ΔT = T·a + J·T⁻¹·b into two scalar
unknowns (a,b) solved by ordinary (non-degenerate) cohomological equations
— never invoking a change of coordinates. **Theorem 1** is the resulting a
posteriori KAM statement: given a non-degenerate approximate torus T with
small enough residual R, there is a true ω-torus T* nearby, with an explicit
error bound. This is the same Newton-without-transformation-theory idea
Fontich/de la Llave/Sire (paper 1 above) and Huguet/de la Llave/Sire (paper
3 below) build on and cite directly ([LGJV05], their own abbreviation for
this exact paper).

### 2.2 What's present vs. deferred

Present: the adapted-Newton-method derivation (§"The adapted Newton
method"), Lemma 1, Theorem 1, and a bulleted list of seven claimed
generalizations (Lindstedt series around a torus, KAM "stickiness"/
Nekhoroshev-type results, exponentially small measure estimates, finite
differentiability, weaker non-degeneracy conditions, lower-dimensional tori,
time-quasi-periodic/general-symplectic-manifold extensions) — stated as
"part of these results have been achieved, and other ones are work in
progress," with **no proofs, no explicit estimates, and no worked numerical
examples** for any of the seven. No celestial-mechanics example is computed
in this short version.

### 2.3 Cross-check

No direct citation in this project's code; the theorem/lemma structure here
is confirmed (by direct comparison) to be the ancestor both `genome/qp_tori.py`'s
different (GMOS) lineage and the Kumar-lineage whiskered-tori papers'
"KAM-without-action-angle-variables" citation both ultimately point back to
— consistent with `#741`'s and `#733`'s prior tracing, no new discrepancy.

---

## 3. Huguet, de la Llave & Sire 2012 — "Computation of Whiskered Invariant
Tori and Their Associated Manifolds: New Fast Algorithms," *DCDS-A*
32(4):1309–1353

DOI `10.3934/dcds.2012.32.1309` (CrossRef-confirmed). arXiv:1004.5231.

### 3.1 The algorithms (the paper is purely algorithmic — no rigorous proofs
of its own, no numeric worked examples of its own)

Explicitly the *algorithmic* companion to the *rigorous* results of paper 1
above ([FdlLS09a/b], i.e. Fontich/de la Llave/Sire) and paper 2's full
version ([LGJV05]) — "we will only discuss the algorithmic issues... [the
rigorous papers] contain estimates." Four algorithms are presented:

1. **Lagrangian KAM tori** and **whiskered KAM tori**: Newton method on the
   invariance equation F(K(θ)) − K(θ+ω) = E(θ), decomposed via the
   invariant-subspace splitting (stable/unstable/center) into three coupled
   sub-equations, each solved with **O(N) storage and O(N log N) operations**
   per Newton step (N = number of Fourier coefficients discretizing the
   torus) — an explicit complexity claim, using the FFT. Truncation error is
   O(exp(−CN^{1/d})), d = torus dimension.
2. **Invariant splittings** (stable/unstable/center bundles): computed via an
   equation for the invariant *projections* directly (not the standard graph
   transform), with an accelerated variant achieving superexponential
   convergence via fast cohomology-equation solvers (Appendix A).
3. **Stable/unstable manifolds**: based on the parameterization method
   ([CFL03a/CFL05], i.e. the Cabré/Fontich/de la Llave series already
   acquired under `#741`), but replacing their original *contractive*
   iteration with a **Newton** iteration, again O(N log N) per step —
   explicitly "the algorithms to compute the stable and unstable manifolds
   had not been previously discussed" (i.e. this paper's own novel
   contribution relative to [CFL03a/CFL05]).
4. Applies uniformly to both **primary and secondary** tori (per the same
   distinction extracted from paper 1, §1.1 above), and to whiskered tori.

No numeric worked examples/timing/accuracy tables are given anywhere in this
paper — the text explicitly defers "implementation details and the results
of several runs" to a companion manuscript ([HdlLS09], never published as
far as this pass could confirm), consistent across a search of the whole
document (confirmed by grep: zero occurrences of "Table" with numeric
content, zero timing figures).

### 3.2 Cross-check

No direct citation in this project's code. Confirms (does not contradict)
`#741`'s and `#733`'s established finding that this whole parameterization-
method/whiskered-tori algorithmic lineage remains correctly un-adopted —
`genome/qp_tori.py` implements the separate Olikara & Scheeres GMOS
collocation/shooting lineage, not this Newton-on-Fourier-coefficients
approach.

---

## 4. Capiński, Gidea & de la Llave 2016 — "Arnold Diffusion in the Planar
Elliptic Restricted Three-Body Problem: Mechanism and Numerical
Verification," *Nonlinearity* 30(1):329

DOI `10.1088/1361-6544/30/1/329` (CrossRef-confirmed). arXiv:1510.00591.

### 4.1 Precisely what this paper proves (and what it does NOT)

**This is not a whiskered-tori/KAM-persistence paper — the paper's own text
says so explicitly: "KAM plays a very minor role."** The mechanism (from a
companion paper [Gidea, de la Llave & Seara 2014, arXiv:1405.0866, itself a
genuine new citation-mining gap — see §8]) is built on **normally hyperbolic
invariant manifolds (NHIMs) and the scattering map**, not on persistence of
KAM tori under perturbation:

- Start from a NHIM Λ (in the PER3BP application: a normally hyperbolic
  manifold foliated by the Lyapunov-orbit family near L2) whose stable and
  unstable manifolds intersect transversally.
- A transverse homoclinic intersection defines a **scattering map**: given a
  point on the intersection, map the asymptotic-past orbit to the
  asymptotic-future orbit. This is computed via Melnikov-type integrals in
  the perturbative regime.
- **Theorem 17** (the paper's main applied result): if a finite,
  explicit, numerically-checkable set of conditions hold (transversality,
  existence of "homoclinic channels," and a sign condition on the derivative
  of the scattering-map-composed-with-itself, Eqs. 42-43), then for
  sufficiently small perturbation ε there exist heteroclinic/diffusing
  orbits connecting arbitrarily different energy levels — with the *size*
  of the energy change independent of ε (unlike prior "micro-diffusion"
  results the paper explicitly contrasts itself against, e.g. its own
  ref [11]).
- Explicitly does **not** require verifying twist conditions, Aubry-Mather
  theory, or KAM non-degeneracy conditions — the whole point of the method
  is to sidestep those.
- Interesting nuance the paper itself states (Remark 19, §1.2): for the
  actual Jupiter-Sun PER3BP application, the authors note they *can* also
  verify the KAM twist condition for a range of energies, and where KAM
  applies it *rules out* diffusion along the inner (torus) dynamics — so
  diffusion in that regime is forced to occur via the alternative
  (jumping/scattering-map) route the paper's mechanism targets. This is the
  one place KAM tori enter at all: as evidence that the *other* diffusion
  channel is closed off, not as the mechanism being used.

Numerically: applied to the Jupiter-Sun PER3BP (planar elliptic, not
circular), around Lyapunov orbits near L2, with error estimates the authors
describe as "not completely rigorous" but accurate to 10⁻⁶–10⁻⁸ relative
error (standard numerical-analysis checks, not a computer-assisted proof —
the authors explicitly flag CAPD-based rigorous validation as future work).
No single clean sourced (μ, energy-jump) numeric table is given as a
reusable positive control; the paper's contribution here is the mechanism
and its verification methodology, not a tabulated result.

### 4.2 Precision on the task's own framing

The task brief asked whether this is "KAM-persistence... underlying
torus-existence claims" — **on direct read, this characterization needs
correcting**: this paper's own text is explicit that KAM/torus-persistence
is *not* the diffusion mechanism (NHIM+scattering map is), and plays only
the minor, negative role described above (ruling out the competing inner
channel). If the already-corpus Kumar/Anderson/de la Llave 2022 CMDA
paper's own §3 cites this Capiński/Gidea/de la Llave paper for a genuine
KAM-persistence step, that citation is drawing on this paper's incidental
Remark 19 observation, not its main theorem — worth flagging for anyone
auditing that CMDA paper's own citation precision, though this digest does
not re-verify the CMDA paper's text (out of this task's scope).

### 4.3 Cross-check

No direct citation in this project's code. The scattering-map/NHIM
machinery here is a genuinely different tool from this project's own
manifold/torus code (which does not implement scattering maps at all) —
flagged as a possible future capability (a scattering-map-based
heteroclinic/diffusion search would be structurally different from the
project's current mesh-intersection-based CCR4BP manifold search in
`kumar-anderson-delallave-2023-...-acta-astro-211...`), not a bug or gap.

---

## 5. Mireles James & Murray 2017 — "Chebyshev-Taylor Parameterization of
Stable/Unstable Manifolds for Periodic Orbits: Implementation and
Applications," *IJBC* 27(14):1730050

DOI `10.1142/S0218127417300506` (CrossRef-confirmed). arXiv:1706.03345.

### 5.1 What "Chebyshev-Taylor" means, precisely

The standard prior approach (the paper's own baseline, e.g. Cabré/Fontich/de
la Llave and the already-acquired Gonzalez & Mireles James 2017 SIADS paper)
is **Fourier-Taylor**: represent the periodic-orbit direction with a Fourier
series (natural since it's periodic) and the manifold's transverse/fold
direction with a Taylor series, solving the resulting homological equations
recursively. Fourier coefficients decay slowly for long or "complicated"
(high-harmonic-content) periodic orbits, forcing more and more coefficients
— for long orbits this becomes impractical. This paper's replacement:
**Chebyshev spectral methods for the periodic-orbit direction too**,
treating the periodic orbit as a sequence of coupled **boundary value
problems on smaller (non-periodic) sub-domains** rather than one global
Fourier expansion. Chebyshev series share Fourier's efficiency advantages
(differentiation is a — tridiagonal, not diagonal — operation in the
transform domain; a fast cosine transform substitutes for the FFT for
nonlinearities) while natively handling non-periodic BVPs, which lets the
decay rate of the coefficients on *each* sub-domain be controlled
independently of the orbit's total harmonic complexity. The method is thus
"Chebyshev-Taylor" in both directions used together — Chebyshev for
period/BVP-subdomain structure, Taylor for the transverse manifold-fold
direction — not a mechanical swap of one basis for another.

### 5.2 Worked examples (sourced)

- **Lorenz system** (not celestial mechanics): periodic orbits labeled by
  symbolic itineraries (e.g. "ABB," "AABBB"); Table 1 gives a genuine
  numeric accuracy/order table — conjugacy-error vs. Taylor order N (20 to
  100) and Chebyshev order K (50 to 250), errors ranging 9.19×10⁻⁴ (N=20)
  down to 7.2×10⁻¹⁰ (N=100) at short test time, with the norm of the last
  Taylor coefficient tracked as an independent convergence diagnostic
  (down to 6.6×10⁻²⁰ at N=100/K=250).
- **Planar CRTBP** (Earth-Moon, μ=0.0123, Jacobi constant 3.17): L1/L2
  Lyapunov orbits, their stable/unstable manifolds, and a computed
  heteroclinic connecting orbit between them (Figs. 10-13) — qualitative
  (figure-based), no digitizable numeric table.
- **Planar equilateral CRFBP** ("Circular Restricted Four Body Problem" in
  the paper's own terminology — **three equal-triangle co-orbiting
  primaries**, structurally a **different** four-body model from this
  project's own `core/ccr4bp.py`, which is the Blazevski-Ocampo-lineage
  *hierarchical/nested* CCR4BP, not the equilateral-triangle CRFBP — see
  §5.3 below): a Lyapunov orbit about the L8 libration point (near the
  heaviest of the three primaries), computed with D=4, N=60, and a reported
  Floquet exponent λ ≈ ±0.0538, explicitly noted by the authors as "much
  closer to zero than in any previous example" in the paper — the weakest
  hyperbolicity case tested, a useful edge-case reference if this project
  ever needs a near-degenerate-hyperbolicity positive control.
- A general complexity note: a degree-N=10 Taylor polynomial in D=6
  variables has 286 coefficients — illustrative of why full polynomial
  (non-recursive) approaches become impractical, motivating the recursive
  homological-equation approach used throughout.

### 5.3 Cross-check: model-mismatch flag for this project's own CCR4BP

**Important, precise finding**: this paper's own "CR4BP"/"CRFBP" is the
classical **equilateral-triangle** restricted four-body problem (three
primaries at the vertices of a Lagrange central configuration, all
co-orbiting at the same period) — a genuinely different four-body model
from this project's own `src/cyclerfinder/core/ccr4bp.py` (which follows
the Blazevski & Ocampo 2012 **concentric/hierarchical** CCR4BP: Jupiter +
Europa as the base CR3BP, Ganymede as a *periodic time-dependent
perturbation*, not a third co-equal primary). Any future adoption of this
paper's Chebyshev-Taylor manifold method for this project's own CCR4BP
periodic orbits would need to re-derive the method against the
hierarchical model's equations of motion — the paper's own CRFBP worked
example (§5.2 above) is not a literal positive control for
`core/ccr4bp.py`, only a same-family (four-body, planar) sanity check of the
*numerical method's* accuracy/order behavior. No direct citation of this
paper exists in the project's code; the method remains an un-adopted future
capability, consistent with the project's "digest ≠ adoption" convention —
not a bug.

---

## 6. Pérez-Palau, Masdemont & Gómez 2015 — "Tools to Detect Structures in
Dynamical Systems Using Jet Transport," *CMDA* 123(3):239–262

DOI `10.1007/s10569-015-9634-3` (CrossRef-confirmed per the task brief).

### 6.1 What "jet transport" is, precisely, and how it differs from the
Cauchy-Green/LCS approach

**Jet transport** (aka flow expansion / box propagation; originated by Berz
& Makino for particle-accelerator physics) replaces a single numerically
integrated initial condition x₀ with a **truncated polynomial**
P(ξ) = x₀ + ξ parameterizing a whole neighborhood, and propagates that
polynomial through the numerical integrator using **polynomial arithmetic**
at every step (the integrator's usual floating-point operations become
polynomial-algebra operations — addition, product, composition — on
truncated Taylor series in ξ). The independent (ξ⁰) term recovers the
ordinary trajectory of x₀; each successive-order term in ξ is exactly the
corresponding order of the variational equations, obtained **without ever
writing those variational ODEs down explicitly** (avoiding the standard
approach's need to symbolically derive and separately integrate 1st/2nd/...
order variational equations).

The classical **Lagrangian Coherent Structures (LCS)** approach (Haller)
instead computes the **Cauchy-Green strain tensor** (transpose of the state
transition matrix times itself) at each grid point and extracts its
eigenvalues (giving the Finite-Time Lyapunov Exponent field) — this only
uses the *first-order* variational information (the STM) and requires an
explicit eigen-decomposition at every sampled point. This paper's
alternative: read the maximum-contraction/expansion behavior **directly
off the jet's own high-order Taylor coefficients**, without ever forming or
diagonalizing a Cauchy-Green tensor. The flagship indicator (§3.1,
"maximal initial boxes") is explicit and simple:
ξ_max = min_{|k|=n} (ε_jet / |a_{m,k}|)^{1/k}, where a_{m,k} are the Taylor
coefficients of the propagated jet — a point near a hyperbolic structure
requires a smaller initial box to hold a given target accuracy ε_jet than a
point far from one, giving a scalar field directly analogous to the FTLE
field but computed from a single jet-propagation pass rather than a
tensor eigenvalue problem. Two further indicators (based on directly
extracting maximum-contraction/expansion *directions* and their rates from
consecutive jet orders) are also given but require order ≥2 jets (vs. order
1 sufficing for the box-size indicator alone).

### 6.2 Worked examples (sourced, though schematic)

Simple and periodically-perturbed pendulum (Taylor-integration order 25,
jet order 5); planar CR3BP and ER3BP with a **toy** mass ratio μ=0.1 (not
Earth-Moon or Sun-Earth), Jacobi-constant-style energy levels
E₀=Ẽ₀=−1.8, jet order 2, final time T=2 (adimensional), ε_jet=10⁻⁶ target
precision — all figure-based (spatial maps of ξ_max over a grid), no
digitizable numeric accuracy/timing table comparable to paper 5's Table 1.
The paper reports (qualitatively, Fig. 6) that going from jet order 1 to
order 5 improves box-size resolution near hyperbolic structures by ~10⁻⁵,
with larger gains farther from them.

### 6.3 Cross-check

No direct citation of this paper, jet transport, differential algebra, or
COSY-INFINITY-style automatic-domain-splitting techniques anywhere in
`src/cyclerfinder/`. This project's own manifold/torus search code uses
explicit variational-equation propagation (STM integration), not jet
transport — a genuinely different, currently-unused technique for the same
underlying goal (locating hyperbolic structures / invariant-manifold
neighborhoods). Flagged as a future-capability candidate (particularly for
LCS-style non-autonomous-system structure detection, which this project's
current STM-based tools do not directly address), not a bug.

---

## 7. Fernández-Mora, Haro & Mondelo 2024 — "Flow Map Parameterization
Methods for Invariant Tori in Quasi-Periodic Hamiltonian Systems," *SIAM J.
Appl. Dyn. Syst.* 23(1):127–166

DOI `10.1137/23M1561257`. **Corrections to the master list's old §5 item
39 entry, both confirmed directly from this PDF's own title page/preprint
header (per the task brief):** (a) correct surname is
**"Fernández-Mora"**, not "Fernández" — used throughout this digest; (b) not
an unpublished/undated preprint — the PDF is headed "This is a preprint of:
... SIAM J. Appl. Dyn. Syst., vol. 23(1), 127–166, 2024. DOI:
[10.1137/23M1561257]," an author-confirmed preprint of the exact published
paper, content should track the published version closely (distinct from
item 2/paper 2 above, which genuinely is an abbreviated precursor
communication, not a preprint of the full paper).

### 7.1 Method

Generalizes the **flow map parameterization method** — previously built for
autonomous/periodic Hamiltonians in the direct predecessor paper Haro &
Mondelo 2021 ([HM21], *CNSNS* 101:105859, a genuinely new citation-mining
gap, see §8) — to **non-autonomous, quasi-periodically-forced** Hamiltonian
systems. Key technical machinery, introduced specifically to make this
generalization work: **fiberwise isotropic tori**, **fiberwise symplectic
deformations**, and their **moment maps**, which together produce "magic
cancellations" licensing solvability of the small-divisor (cohomological)
equations in this genuinely time-dependent setting — without resorting to
the standard trick of artificially extending phase space with fictitious
angle variables to make the system autonomous (which the authors explicitly
avoid because it inflates the dimension and degrades efficiency). The flow
map itself (rather than a fixed-time Poincaré/stroboscopic map) is used to
reduce the torus dimension by one relative to the ambient phase space,
following [GM01, HM21]. Complexity: the same O(N) storage / O(N log N)
per-Newton-step scaling as the rest of this citation lineage, for N Fourier
coefficients/grid points.

### 7.2 Worked example and sourced positive-control data

Sun-Earth **Elliptic** Restricted Three-Body Problem (ERTBP): μ =
3.040357143×10⁻⁶, e = 0.01671123. The authors compute **3-dimensional**
(spatial!) partially-hyperbolic invariant tori around the L1 point:
starting from 2D tori in the L1 center manifold of the (autonomous) CRTBP —
77 tori around vertical Lyapunov orbits, rotation number ρ ∈
[0.03565, 0.0961], extended by continuation in flying time to **8971**
CRTBP tori total — then lifting/continuing each to the ERTBP's added
(true-anomaly) external frequency, of which **4457** successfully reached
the full Sun-Earth eccentricity e=0.01671123. Numerical settings: initial
external-phase grid N₂=16, up to 1024 Fourier coefficients per phase,
multiple-shooting order m=4, Newton tolerance K=10⁻⁹; invariance-equation
residuals as small as **10⁻¹⁵** were achieved for some tori. This is a
concrete, sourced, independently-checkable dataset for anyone wanting to
validate an implementation of this method or a close relative.

### 7.3 The task's core cross-check question: does this generalize to a
spatial CCR4BP?

**Yes, plausibly cleanly, and more directly than a naive reading suggests.**
This project's own `src/cyclerfinder/core/ccr4bp.py` (module docstring,
lines 1-40) is explicit that its CCR4BP is **TIME-PERIODIC in the
Jupiter-Europa synodic rotating frame** — Ganymede's position, viewed in
that frame, evolves at the single beat frequency (n₃−n₂) between the two
moons' mean motions, structurally the same **single-external-frequency
non-autonomous** class as the Sun-Earth **ERTBP** this paper's own worked
example targets (the ERTBP's non-autonomy is likewise a single external
frequency — the true-anomaly rate of the primaries' elliptic orbit). That
is: this paper's demonstrated case (spatial, single-frequency,
non-autonomous, quasi-periodic Hamiltonian torus computation around a
libration point) is a **closer structural match** to a hypothetical spatial
CCR4BP than the "quasi-periodic" framing in the task brief and master list
might suggest — no need to first generalize to genuinely incommensurate
multi-frequency forcing, since a single-perturber CCR4BP is already exactly
the single-frequency case this paper solves.

What would actually be needed to adopt it for a **spatial** CCR4BP:
1. `core/ccr4bp.py`'s own equations of motion are currently **planar only**
   (confirmed by direct grep: Ganymede is forced onto z=0, "Ganymede is
   planar (z=0); indirect z-term vanishes" per the module's own code
   comments) — a genuinely 3D CCR4BP right-hand side would need deriving
   first, independent of which torus-computation method is used downstream.
2. This paper's own worked example is a **3D torus around a spatial
   libration point** (L1, center×center×saddle spectrum) — the CCR4BP's own
   Jupiter-Europa collinear points have the analogous 4D-center/2D-hyperbolic
   spectral structure in the planar case; a spatial extension would need the
   full 6D-center-manifold spectral decomposition this paper's §4.2 works
   through for the ERTBP's L1 (2 elliptic pairs + 1 hyperbolic pair, vs. the
   planar case's 1 elliptic pair + 1 hyperbolic pair).
3. The "magic cancellations"/fiberwise-symplectic-deformation machinery
   (Appendices A-C) is proved for a **general** quasi-periodic non-autonomous
   Hamiltonian and does not appear to depend on any ERTBP-specific structure
   — i.e. this part of the method should transfer directly to a spatial
   CCR4BP's equations of motion once those are written down, without
   requiring new theory.

In short: the theoretical machinery generalizes without modification; the
concrete blocker is that this project's own CCR4BP model needs a spatial
(non-planar) equations-of-motion derivation first, which is a pre-existing,
independently-known gap (not something this paper's citation-mining
surfaces as new) rather than a limitation of this paper's method itself.

### 7.4 Cross-check: code

No direct citation of this paper or of flow-map parameterization methods
anywhere in `src/cyclerfinder/`. Consistent with the "digest ≠ adoption"
convention — flagged as a genuine, well-specified future-capability route
for a spatial CCR4BP extension (not a bug, and not yet started).

---

## 8. Mandatory citation-mining pass — consolidated across all seven papers

Each paper's own introduction/background and reference list was read;
candidates were cross-checked directly against `docs/notes/CORPUS_INDEX.md`
and `docs/notes/2026-07-27-730-acquisition-backlog-master-list.md` (grepped,
not inferred from memory). Given how densely cross-cited this whole
Haro/de la Llave/Fontich/Sire/Huguet/Capiński/Gidea lineage is, the large
majority of citations across all seven papers duplicate gaps `#730`, `#733`,
`#728`, and `#741` already flagged — those are consolidated in §8.2, not
re-flagged individually. **No acquisition was performed for any candidate
below** — flagged only, per policy.

### 8.1 Genuinely new gaps (not found in `CORPUS_INDEX.md` or the master list)

| Candidate | Priority | Why relevant | Source paper(s) |
|---|---|---|---|
| Haro, À. & Mondelo, J.M., "Flow map parameterization methods for invariant tori in Hamiltonian systems," *Commun. Nonlinear Sci. Numer. Simulat.* 101:105859 (2021) | **HIGH** | The direct, explicitly-cited method-predecessor ([HM21]) that paper 7 (Fernández-Mora/Haro/Mondelo 2024) generalizes from autonomous/periodic to quasi-periodic Hamiltonians — cited >10 times in paper 7's own text; the natural "read this first" companion for anyone implementing paper 7's method, and possibly directly sufficient on its own for the CCR4BP's *planar* case (already periodic, no need for the quasi-periodic generalization). | Paper 7 |
| Gidea, M., de la Llave, R. & Seara, T., "A general mechanism of diffusion in Hamiltonian systems: Qualitative results," arXiv:1405.0866 (2014) | **HIGH** | The actual foundational paper ([28]) whose diffusion mechanism (NHIM + scattering map + shadowing, no KAM/twist conditions required) paper 4 (Capiński/Gidea/de la Llave 2016) applies to the PER3BP — paper 4 is explicitly an *application* of this paper's general qualitative theorem, not a self-contained derivation. | Paper 4 |
| Xue, J., "Arnold diffusion in a restricted planar four-body problem," *Nonlinearity* 27(12):2887 (2014) | **HIGH** | Directly a restricted **four-body** problem Arnold-diffusion result — the only four-body-native item surfaced by this whole citation-mining pass, distinct from (and possibly a useful cross-check against) this project's own CCR4BP work, though via a different (equilateral-type) 4BP model, not confirmed which. | Paper 4 |
| Delshams, A., Gidea, M. & Roldán, P., "Arnold's mechanism of diffusion in the spatial circular restricted three-body problem: A semi-analytical argument," *Physica D* 334:29–48 (2016) | Medium-high | A **spatial (3D)** CR3BP Arnold-diffusion mechanism (scattering map + weak-twist argument for out-of-plane amplitude growth) — directly relevant given this project's own interest in spatial extensions (see §7.3 above); explicitly contrasted against paper 4's own (planar, elliptic) mechanism in paper 4's own introduction. | Paper 4 |
| Cabré, X., Fontich, E. & de la Llave, R., "The parameterization method for invariant manifolds. I. Manifolds associated to non-resonant subspaces," *Indiana Univ. Math. J.* 52(2):283–328 (2003) | Medium | Part I of the Cabré/Fontich/de la Llave series whose Part III is already acquired (`#741`); repeatedly cited as [CFL03a]/[CFdlL03a] across nearly every paper in this whole citation cluster (papers 1, 3, 5, 7 here) but never itself flagged as an individual backlog item before this pass — only mentioned in prose. | Papers 1, 3, 5, 7 |
| Cabré, X., Fontich, E. & de la Llave, R., "The parameterization method for invariant manifolds. II. Regularity with respect to parameters," *Indiana Univ. Math. J.* 52(2):329–360 (2003) | Medium | Companion to the above, same status (cited repeatedly, never individually flagged before). | Papers 1, 3, 5, 7 |
| Delshams, A., Gidea, M. & Roldán, P., "Transition map and shadowing lemma for normally hyperbolic invariant manifolds," *Discrete Contin. Dyn. Syst.* 33(3):1089–1112 (2013) | Medium | The specific shadowing-lemma tool underlying paper 4's diffusion-orbit construction methodology. | Paper 4 |
| Wilczak, D. & Zgliczyński, P., "Heteroclinic connections between periodic orbits in planar restricted circular three-body problem — a computer assisted proof," *Comm. Math. Phys.* 234(1):37–75 (2003) | Medium | Rigorous (computer-assisted) heteroclinic-connection existence proofs in the PCR3BP — a different validation standard (rigorous vs. numerical) from most of this project's own corpus, directly on the project's own search domain (heteroclinic connections between periodic orbits). | Paper 4 |
| Galante, J. & Kaloshin, V., "Destruction of invariant curves in the restricted circular planar three-body problem by using comparison of action," *Duke Math. J.* 159(2):275–327 (2011) | Medium | PCR3BP KAM-curve-destruction result; directly relevant to a "how far does the KAM torus region actually extend" question in the same energy regime this project's own CCR4BP/CRTBP torus searches operate in. | Paper 4 |
| Canadell, M. & Haro, À., "Computation of Quasi-Periodic Normally Hyperbolic Invariant Tori," Parts I and II, *SIAM J. Appl. Dyn. Syst.* (2017) [exact volume/pages not independently resolved this pass] | Medium | Directly NHIT (normally hyperbolic invariant tori) computation, the exact object class both this project's own `qp_tori.py`/`variational_*_torus.py` and this whole citation lineage's algorithms target — repeatedly cited across papers 5 and 7 as [CH17a]/[CH17b], never yet individually flagged. | Papers 5, 7 |
| Gawlik, E.S., Marsden, J.E., Du Toit, P.C. & Campagnola, S., "Lagrangian coherent structures in the planar elliptic restricted three-body problem," *Celest. Mech. Dyn. Astron.* 103(3):227–249 (2009) | Medium | Directly PER3BP LCS — a genuine cross-check candidate for paper 6's own jet-transport-vs.-CG-tensor LCS comparison, in exactly the elliptic-RTBP setting paper 4 also uses. | Paper 6 |
| Short, C. & Howell, K., "Lagrangian coherent structures in various map representations for application to multi-body gravitational regimes," *Acta Astronautica* 94(2):592–607 (2014) | Medium | LCS applied specifically to multi-body (not just 3-body) gravitational regimes — directly relevant to this project's own multi-body (CCR4BP/CRNBP) search domain. | Paper 6 |
| Luque, A. & Villanueva, J., "A KAM theorem without action-angle variables for elliptic lower dimensional tori," *Nonlinearity* 24(4):1033–1080 (2011) | Low-medium | Direct descendant of paper 2's own de la Llave/González/Jorba/Villanueva lineage, extended to elliptic (not just hyperbolic/whiskered) lower-dimensional tori — a plausible complementary tool if the project ever needs elliptic-normal-direction tori rather than purely hyperbolic ones. | Paper 7 |

### 8.2 Already-tracked, not double-counted

Confirmed by direct grep of `CORPUS_INDEX.md`/master-list, the following
recurring citations across all seven papers are **already** flagged
(`#730` §5 items 31/33/36, `#733`'s citation-mining, `#741`'s
citation-mining) or already **acquired** — listed only to confirm no new
gap: Haro & de la Llave 2006a/b/2007 (§5 item 31 + already-traced ancestry);
de la Llave/González/Jorba/Villanueva 2005 full paper (§5 item 33, this
task's own paper 2 IS the short precursor of exactly this item); Zhang & de
la Llave 2018 (§5 item 36); Cabré/Fontich/de la Llave 2005 Part III
(**ACQUIRED**, `#741`); Haro/Canadell/Figueras/Luque/Mondelo 2016 book
(**ACQUIRED**, `#733`); Delshams/Masdemont/Roldán 2008 scattering-map paper
(flagged by `#733`); Canalias/Delshams/Masdemont/Roldán 2006 scattering-map
paper (flagged by `#733`); Jorba & Olmedo 2009 (flagged by `#733`); Calleja
& de la Llave 2009 "Fast numerical computation of quasi-periodic
equilibrium states..." (flagged by `#741`); Broucke 1968 (now independently
re-surfaced a 5th time across the whole backlog — by paper 4's own
reference list here — the single most-corroborated unacquired gap in the
entire project backlog); Llibre/Martínez/Simó 1985 (already in master list
§4 item 28); Masdemont 2005 "High-order expansions..." (flagged by `#741`);
Burgos-García/Lessard/Mireles-James 2019 CR4BP rigorous-numerics paper
(flagged HIGH by `#741` — distinct from, but same-lineage as, several
equilateral-CRFBP citations papers 5/6 also surface, see §8.3).

### 8.3 Low-priority background / large clusters not itemized individually

Per the project's own established convention (background/textbook items,
or large same-topic clusters, flagged once as a bucket rather than itemized
DOI-by-DOI): (a) the extensive classical-KAM/small-divisor background common
to papers 1, 2, and 3 (Zehnder 1975/76, Moser 1962/66/67, Rüssmann,
Kolmogorov 1954, Arnol'd 1963, Poincaré, and similar) — none
celestial-mechanics-specific, none independently DOI-resolved this pass; (b)
a substantial **equilateral-triangle CRFBP periodic-orbit/relative-equilibria
census cluster** surfaced by paper 5's own reference list (Baltagiannis &
Papadakis 2011a/b, Barros & Leandro 2014, Leandro 2006, Papadakis 2016a/b,
Burgos-García & Delgado 2013, Burgos-García & Gidea 2015, Cheng & She 2017,
and Burgos-García 2017 "private communication") — genuinely on-topic
(restricted four-body problem) but for the **equilateral**, not this
project's hierarchical/concentric, four-body model (per §5.3's model-mismatch
finding above); not itemized individually given the size of the cluster and
the model mismatch, but flagged as a cluster worth a dedicated pass if this
project ever pursues the equilateral-CRFBP variant specifically; (c) Delshams/
Kaloshin/de la Rosa/Seara 2015 "Global instability in the elliptic restricted
three body problem" (unpublished/2015, cited by paper 4) and Fejoz/Guardia/
Kaloshin/Roldán 2011 "Kirkwood gaps and diffusion along mean motion
resonances in the restricted planar three-body problem" (arXiv preprint,
cited by paper 4) — both genuinely relevant Arnold-diffusion/PCR3BP papers,
Low-medium priority, not independently pursued this pass given the density
of higher-priority items already in §8.1.

---

## 9. Summary: no code fixes required; one significant precision correction

Across all seven papers: (1) paper 1's a posteriori theorem, vanishing
lemma, and primary/secondary-tori distinction are extracted and confirmed
consistent with the already-established Kumar-lineage ancestry — no code
implication (theory remains correctly un-adopted, `genome/qp_tori.py` on
the separate GMOS lineage). (2) Paper 2 is confirmed to be the abbreviated
short-communication precursor, not the full paper — flagged explicitly so
future readers don't assume full coverage; no code implication. (3) Paper
3's algorithms (O(N log N) Newton-on-Fourier-coefficients for tori,
splittings, and manifolds) are purely algorithmic with no numeric worked
examples of its own; consistent with the un-adopted-parameterization-method
finding, no new implication. (4) **Paper 4 required a real correction to
the task brief's own framing**: this paper's diffusion mechanism is NHIM +
scattering map, explicitly **not** KAM-torus persistence ("KAM plays a very
minor role" per the paper's own text) — recorded here for the record,
relevant to anyone citing this paper as a torus-persistence source
downstream. (5) Paper 5's Chebyshev-Taylor method is confirmed as a genuine
second, independent high-order manifold-parameterization technique, but its
own worked CR4BP example is the **equilateral**, not this project's
hierarchical, four-body model — a model-mismatch worth remembering before
treating that example as a literal `core/ccr4bp.py` positive control. (6)
Paper 6's jet-transport indicator is confirmed as a genuinely different
(non-Cauchy-Green-tensor) route to the same LCS-style hyperbolic-structure
detection this project's manifold code does not currently implement — future
capability, not adopted. (7) Paper 7's flow-map parameterization method is
confirmed to plausibly generalize to a spatial CCR4BP with no new theory
required beyond what's already proved for the general quasi-periodic
Hamiltonian case — the actual blocker is this project's own CCR4BP model
being planar-only (`core/ccr4bp.py`, independently confirmed by grep), a
pre-existing, already-known gap, not something this paper's own content
newly reveals. No `CorpusAnchor`, `literature_check.py`, or search-code
changes are warranted by any of the seven papers.
