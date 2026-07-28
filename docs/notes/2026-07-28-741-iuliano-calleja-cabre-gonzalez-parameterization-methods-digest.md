# #741 — Iuliano 2016 thesis, Calleja et al. 2021, Calleja & de la Llave 2010,
# Cabré/Fontich/de la Llave 2005, Gonzalez & Mireles James 2017: the last five
# `#730` top-10 backlog items (§2 items 6-10)

**Task:** `#741`, items 6-10 of `#730`'s own §2 ranking, all five PDFs
user-supplied directly (page-1-verified by the coordinating session before
dispatch). Two of the five (items 8 and 9) required a re-supply round: the
master list's original "standard pattern, not independently re-verified" DOI
guesses for both turned out to be **wrong** — the first uploads at those
guessed DOIs resolved to unrelated papers (a Matthes & Toscani kinetic-theory
paper at the guessed Calleja & de la Llave DOI; an Azevêdo & Ontaneda paper at
the guessed Cabré/Fontich/de la Llave DOI). The coordinating session then
found the correct DOIs via fresh WebSearch and the user re-supplied the
correct PDFs, both confirmed correct in this pass. Sibling tasks `#732`/`#733`
processed the other five of `#730`'s original top-10 items (see
`2026-07-27-732-blazevski-negri-baresi-foundational-papers-digest.md` and
`2026-07-27-733-olikara-thesis-haro-parameterization-book-digest.md`); this
digest is the final tranche.

## 0. Acquisition, filing, text-layer check

All five filed in the private `cyclers_pdf` corpus
(`/Users/bruce/dev/cyclers_pdf/papers/`, separate repo, never committed to the
public `cyclers` repo):

| Paper | Filename | Pages | Text layer |
|---|---|---|---|
| Iuliano 2016 MS thesis | `iuliano-2016-solution-crnbp-planetary-systems-ms-thesis-calpoly.pdf` | 73 | Native, clean (18,030 words via `pdftotext`, no OCR needed) |
| Calleja/del-Castillo-Negrete/Martínez-del-Río/Olvera 2021 | `calleja-delcastillonegrete-martinezdelrio-olvera-2021-new-method-periodic-orbits-symplectic-maps-cnsns-99-105838-doi-10.1016-j.cnsns.2021.105838.pdf` | 17 | Native, clean (11,783 words) |
| Calleja & de la Llave 2010 | `calleja-delallave-2010-numerically-accessible-criterion-breakdown-quasiperiodic-nonlinearity-23-2029-doi-10.1088-0951-7715-23-9-001.pdf` | 31 | Native, clean (16,522 words) |
| Cabré, Fontich & de la Llave 2005 | `cabre-fontich-delallave-2005-parameterization-method-invariant-manifolds-III-overview-applications-jde-218-444-doi-10.1016-j.jde.2004.12.003.pdf` | 72 | Native, clean (30,042 words) |
| Gonzalez & Mireles James 2017 | `gonzalez-mirelesjames-2017-highorder-parameterization-stable-unstable-manifolds-long-periodic-orbits-maps-siads-16-1748-doi-10.1137-16M1090041.pdf` | 48 | Native, clean (21,595 words) |

No OCR needed for any of the five — all are digitally-typeset originals with
full extractable text layers throughout (word counts checked across every
page, not just page 1). Read applying the "chapter-summary scope" convention
for the two longest documents: for the Iuliano thesis (73pp), the EOM/model
chapter (Ch. 2) was deep-read in full, the results chapter (Ch. 4) sampled at
start and end plus its one numeric table, Ch. 1/3/5/6 skimmed via their own
section headers; for the Cabré/Fontich/de la Llave paper (72pp, a dense
theory/proof text), the abstract/introduction and the theorem statements for
Parts I-III were deep-read, the applications sections sampled, and the
proof-heavy middle sections skimmed rather than read line-by-line. The other
three (17, 31, and 48 pages) were read in substantially full detail.

---

## 1. Iuliano 2016 — "A Solution to the Circular Restricted N Body Problem in
Planetary Systems," MS thesis, Cal Poly San Luis Obispo (June 2016)

**Full citation:** Iuliano, Jay R. MS Thesis (Aerospace Engineering),
California Polytechnic State University, San Luis Obispo, June 2016. 73pp.
Committee: Kira Abercromby (chair), Eric Mehiel, Ian Johnson, Paul Choboter.
No DOI; freely available at `digitalcommons.calpoly.edu/theses/1612`.

### 1.1 What the thesis does

Six chapters: Ch. 1 introduces the CR3BP background; Ch. 2 derives a
generalized "circular restricted N-body problem" (CRNBP) — one dominant
"primary" (planet) plus N-1 "secondary" perturbers (its moons) — and its
Jacobi-constant analog; Ch. 3 methodology; Ch. 4 results in the Jovian
system (stable regions, halo orbits/manifolds, low-energy transfers); Ch. 5
future work; Ch. 6 conclusion. The thesis validates its CRNBP model's core
simplifying assumption (treating the system barycenter as coincident with the
primary's own center) against high-fidelity N-body propagation, finding it
holds well below ~1% total-secondary/primary mass ratio (Jupiter+Galileans is
~0.02%) but breaks down near 1-2%.

### 1.2 The equation-defect cross-check (the load-bearing finding)

This project's own `src/cyclerfinder/core/crnbp.py` (lines ~90-110) carries a
code comment crediting the already-acquired Negri & Prado 2022 CRNBP paper
with correcting an "important...to correctly approximate the indirect
effect" gap relative to "Iuliano's now-superseded equations." This digest
independently locates and verifies that gap directly from Iuliano's own
numbered equations rather than taking Negri & Prado's characterization on
faith.

Iuliano's Ch. 2.1 places the primary m₁ at the coordinate origin (not the
standard CR3BP barycentric convention — he explicitly flags this as "the
important differences to note" from his own Ch. 1 background) and derives
(Eqs. 2.3-2.5) the rotating-frame spacecraft acceleration as a **direct**
term (gravity from every body) plus a **flat, single-level indirect** term:

```
ẍₛ − 2ẏₛ − xₛ = −Σᵢ μᵢ(xₛ−xᵢ)/r³ᵢₛ  −  Σⱼ μⱼ xⱼ/r³₁ⱼ
```

(and analogously for y, z), where the indirect sum Σⱼ (j=2..N) is one term
per secondary using only that secondary's own primary-relative position
r₁ⱼ — **no nested/cross term** coupling one secondary's position to
another's. Negri & Prado's corrected per-perturber term is structurally

```
a_j = −μⱼ·[ direct + single-body-indirect + Σ_{k≠j} μₖ·(r_j−r_k)/|r_j−r_k|³ ]
```

i.e. Iuliano's own formula minus exactly the `Σ_{k≠j}` inner coupling sum.
This confirms, from Iuliano's own text, the precise structural defect
`crnbp.py`'s comment attributes to him.

**Important existing-code context, not new to this task:** `crnbp.py`'s own
comment block already proves (and this project's own tests already assert)
that this missing coupling term is dynamically **inert** — it cancels to
exactly zero for the net spacecraft acceleration regardless of N, so
Iuliano's omission, while a real formal/structural gap, would not have
changed any of his own numerical results. This digest's contribution is
confirming the textual grounding of that already-resolved finding, not
reopening it.

A second, weaker, self-admitted caveat in the thesis: origin-at-primary
(rather than true barycenter) is explicitly tested by Iuliano himself and
shown to degrade above ~1-2% mass ratio — a distinct, milder approximation
from the missing coupling term.

### 1.3 Positive-control data

None usable as a rigorous positive control — all "results" in Ch. 4 are
purely graphical (Figs. 2.2-2.11, 4.1-4.28), no digitizable numeric
periodic-orbit/Jacobi-constant tables. The one numeric table (4.1, "Canonical
Conversions") gives Jupiter-Galilean and Sun-Jupiter canonical mass/distance/
time units — useful only for unit-normalization sanity checks.

### 1.4 Note on a separate, unrelated "Iuliano" citation

An already-acquired digest
(`2026-07-27-722-baresi-owen-scheeres-tri-circular-problem-digest.md`, lines
~320-346) separately discusses an unverifiable "Iuliano & Gomes 2019,
*Astrophys. Space Sci.*" journal citation — per `#730`'s own §2 item 6
resolution note, that citation could not be independently verified and this
2016 Cal Poly thesis is the only confirmed-findable paper matching the
"erroneous N+1-body formulation Negri & Prado corrected" description. Do not
conflate the two; this digest concerns only the confirmed, acquired 2016
thesis.

---

## 2. Calleja, del-Castillo-Negrete, Martínez-del-Río & Olvera 2021 — "A new
method to compute periodic orbits in general symplectic maps," *Commun.
Nonlinear Sci. Numer. Simulat.* 99:105838

DOI `10.1016/j.cnsns.2021.105838`, received 5 Mar 2020, accepted 24 Mar 2021.

### 2.1 Method

Targets periodic orbits (rational rotation number p/q) of **general**
area-preserving twist maps that are **not reversible** — for reversible maps
the standard trick restricts the search to 1-D symmetry lines, unavailable
here, forcing a genuinely higher-dimensional root search. The proposed
"compound method" is two-stage: (1) a modified parameterization method
(built on de la Llave et al.'s KAM-without-action-angle-variables machinery),
handling the periodic-orbit case's zero-denominator resonant Fourier modes,
giving an empirically-good (not proven-convergent) seed; (2) a Newton-Gauss
method — explicitly **credited to Haro and collaborators, not original to
this paper** — a block-bidiagonal-with-corner sparse linear solve exploiting
the periodicity-closure structure, contrasted with the naive single-Newton
iteration on the q-fold map DTᵍ, which the paper itself states "is known to
be highly unstable in two or higher dimensions" above period ~10³.

**Stated complexity/scaling:** parameterization step is O(N log N) via FFT
(N = harmonics, empirically N*~2q-4q); Newton-Gauss step storage is O(q)
(12q+4 + 4q memory locations). Combined method demonstrated robust to periods
up to **~10⁷** (quadruple precision), versus ≲10³ for naive direct Newton.

**Benchmark tables (Tables 1-2, own positive control):** Chirikov-Taylor
standard map (κ=0.9600) and a rational harmonic map (κ,α,β)=(1.715,3.0,0.4),
Fibonacci-approximant periodic orbits from 5/8 to 6765/10946: e.g. for
987/1597, error E=6.5275×10⁻⁴⁶, residue R=−3.7914×10⁻¹⁰. Figs. 6-7: residues
converge to the universal critical value **|R|=0.2554** (Greene/MacKay/
Shenker-Kadanoff renormalization prediction) at criticality for both test
maps, Fibonacci orbits N=11 through N=30.

### 2.2 Cross-check: is this "the only directly-competing prior-art method"?

An already-acquired paper, Kumar 2026 "fast-multishooting-periodic-orbits-
symplectic-maps-floquet" (arXiv:2601.00149), frames this Calleja 2021 paper
as its own sole directly-competing prior art. **Partially fair, with a
nuance worth recording**: this paper's genuine, explicitly-stated novelty is
the parameterization-method *seeding* stage that removes the need for a
symmetry-line-derived analytic initial guess; the high-period robustness that
makes it a strong baseline comes from the Newton-Gauss refinement stage,
which the paper itself attributes to prior (Haro et al.) work, not to itself.
A reader wanting to independently audit Kumar 2026's own novelty/speedup
claims should check Kumar's baseline against both halves separately. Also:
this paper's scope is discrete symplectic **maps** (Poincaré-map level);
non-autonomous/time-periodic map extensions are explicitly left to future
work (§6) — worth checking if Kumar 2026 targets continuous-flow periodic
orbits directly, which would be a scope difference beyond a pure speed
comparison.

---

## 3. Calleja & de la Llave 2010 — "A numerically accessible criterion for
the breakdown of quasi-periodic solutions and its rigorous justification,"
*Nonlinearity* 23(9):2029-2058

**Corrected DOI `10.1088/0951-7715/23/9/001`** (the master list's original
guess, `...003`, was wrong — see §0 above and the master-list correction in
this task). Received 31 Aug 2009, published 29 Jul 2010.

### 3.1 The Sobolev-seminorm breakdown criterion

An already-acquired paper (Kumar/Anderson/de la Llave 2023, Acta Astronautica
211, digested at `2026-07-27-727-...-digest.md`) uses "Calleja & de la
Llave['s] Sobolev-seminorm divergence criterion (H²/H³ norms → ∞)" to
distinguish genuine torus breakdown from numerical failure. This paper is the
source, and is entirely theoretical (**no numeric worked examples of its
own** — confirmed by reading all 31 pages). The core pieces:

- **Def. 4.1** (Sobolev norm): Hᵐ on Tⁿ→ℝᵈ, ‖u‖²ₘ = Σ_{k∈ℤⁿ}(1+|k|²)ᵐ|û_k|² —
  H² and H³ (as cited by the Kumar 2023 digest) are just specific integer
  choices of m in this one scale.
- **Meta-theorem 2.2**: the domain of persistence of a Sobolev-regular
  quasi-periodic solution under parameter continuation shrinks only as the
  norm grows; if the norm blows up as λ→λ₀, no persistence past λ₀ — the
  "boundary of KAM tori" detection principle.
- **Theorem 5.3** (a posteriori, symplectic twist-map case): given a small
  enough invariance-equation error, a true nearby solution exists,
  quantitatively bounded — the practical basis for the whole diagnostic.
- **Theorem 5.8** (bootstrap of regularity): a Sobolev-class solution with
  m large enough is automatically analytic — this is what licenses reading
  "Sobolev-norm blow-up" as true analyticity breakdown rather than a mere
  numerical artifact of the corrector.
- Analogous pair (Thm 6.1 + 6.8) for 1-D statistical-mechanics models.
- The paper explicitly does **not** assert a universal blow-up exponent
  itself — §4.1 only notes descriptively that norms blow up "according to a
  power law" near breakdown, citing renormalization-group folklore from
  companion papers, not a result proved here.

**Reuse potential for this project:** directly the kind of "genuine
breakdown vs. numerical-failure" diagnostic that would generalize the
project's own `project_388_wall_energy_selective` pattern (an existing
energy-selective near-anchor-collapse wall) if ever needed for CCR4BP/CRNBP
continuation — flagged as a live, reusable tool, not merely background
theory.

### 3.2 Positive-control caveat

The actual numeric breakdown values (critical ε, rotation numbers, norm-growth
curves) live in two **companion numerical papers**, cited repeatedly here and
not yet in this project's corpus (flagged in §6 below): Calleja & de la Llave
2009a ("Computation of the breakdown of analyticity in statistical mechanics
models," mp_arc preprint 09-56) and 2009b ("Fast numerical computation of
quasi-periodic equilibrium states in 1d statistical mechanics, including
twist maps," *Nonlinearity* 22:1311-1336). If a literal numeric positive
control is ever needed for a reused Sobolev-divergence diagnostic, those two
— not this 2010 paper — are the acquisition targets.

---

## 4. Cabré, Fontich & de la Llave 2005 — "The parameterization method for
invariant manifolds III: overview and applications," *J. Differential
Equations* 218(2):444-515

**Corrected DOI `10.1016/j.jde.2004.12.003`** (the master list's original
guess, `...10.029`, was wrong — see §0 above and the master-list correction).
Received 26 Jul 2004, online 10 Feb 2005.

### 4.1 What it is, and the lineage cross-check

This is explicitly framed (p. 445) as a tutorial overview of the
parameterization method the same three authors introduced across a
three-part series: Parts I ("...manifolds associated to non-resonant
subspaces," *Indiana Univ. Math. J.* 52:283-328, 2003) and II ("...regularity
with respect to parameters," ibid. 52:329-360, 2003) are companion papers by
the same authors; this 2005 JDE paper **is** "Part III" of that series,
covering §4-6 (analytic 1-D stable manifolds, C⁰ stable manifold theorem) and
§7-10 (Cʳ manifolds, non-resonant invariant manifolds for maps and ODEs). It
treats **hyperbolic (non-resonant spectral-subspace) invariant manifolds**,
not quasi-periodic whiskered tori directly — no dedicated whiskered-tori
section exists in this paper itself.

Its own bibliography (p. 514) cites [HdlL03a]/[HdlL03b], the 2003 preprint
precursors of the already-traced Haro & de la Llave 2006ab/2007 papers that
`#733`'s digest (`2026-07-27-733-olikara-thesis-haro-parameterization-book-
digest.md`) already established as the core ancestry for the Kumar-lineage
whiskered-tori/torus machinery. This paper is therefore a direct,
contemporaneous sibling in the same citation graph (shared author de la
Llave), not merely topically adjacent — confirming (with the specific
citation link pinned down) the "recurring x3" flag the master list's own §2
item 9 and multiple Kumar-lineage digests already carried.

### 4.2 Code cross-check verdict: CONSISTENT, no discrepancy

Grepped `src/cyclerfinder/` for "parameteriz|cabre|fontich|haro|whisker": no
direct citation of Cabré/Fontich/de la Llave anywhere in code.
`genome/qp_tori.py` (the project's actual quasi-periodic-torus
implementation) is built on the Olikara & Scheeres 2010/2012 GMOS shooting
scheme (per its own module docstring) — a different, adjacent lineage, not
the Haro/de la Llave/Cabré parameterization-method series. `core/ccr4bp.py`
cites only the Kumar-Anderson-de la Llave-Gunter consumer papers. This is
consistent with the already-flagged `digest ≠ adoption` project pattern: the
Cabré/Fontich/de la Llave-based whiskered-tori parameterization method has
not yet been adopted into any code path — tracked as a future-capability item
(unchanged from the master list's own framing), not a bug or gap requiring a
fix.

---

## 5. Gonzalez & Mireles James 2017 — "High-Order Parameterization of
Stable/Unstable Manifolds for Long Periodic Orbits of Maps," *SIAM J. Applied
Dynamical Systems* 16(3):1748-1795

DOI `10.1137/16M1090041`, received 19 Aug 2016, published 28 Sep 2017.
Companion MATLAB code archived at the authors' webpage (ref [34]).

### 5.1 Method

Computes high-order polynomial approximations of stable/unstable manifolds
attached to **long periodic orbits of maps** by extending multiple-shooting
(the classical periodic-orbit-correction scheme) to the manifold invariance
equation itself, attaching one local chart to *each* point along the orbit
and solving a coupled system of N invariance equations using only Df at each
point — **never** the N-fold composed map fᴺ. The paper's own explicit,
repeated self-description is **"composition-free"** — literally the opposite
of "Taylor-composition with the map."

### 5.2 Cross-check: correcting, not confirming, Kumar 2025's framing

An already-acquired paper, Kumar 2025 "multishooting-parameterization-
invariant-manifolds-heteroclinics-poincare-maps" (arXiv:2509.03655, digested
at `2026-07-27-728-kumar-2025-multishooting-poincare-manifolds-digest.md`),
frames this Gonzalez & Mireles James 2017 paper's own method as "direct
Taylor-composition with the Poincaré map" that its own multi-shooting
approach avoids. **This is not an accurate description of what Gonzalez &
Mireles James actually propose** — it describes the *naive baseline* their
own paper explicitly names and measures itself against, and beats
(§4.2/Fig. 7, Hénon-map period-2 test): the composition-free multiple-shooting
method computes both manifolds to order 50 in 0.004s, vs. 0.035s for one
manifold via naive fᴺ-composition (~5x slower per manifold, ~10x for both);
floating-point cost for the naive approach is "more than 300 times as
expensive" at N=2 already; naive accuracy is "an order of magnitude worse."
The paper's own conclusion: "regardless of the period of the orbit, we deal
only with the nonlinearity of the original map... [naive composition] face[s]
exponential growth of the complexity."

**The real, narrower gap that legitimately motivates Kumar 2025**, stated
precisely so this isn't read as "Kumar's own novelty claim is baseless":
every one of Gonzalez & Mireles James's worked numerical examples (Hénon map,
Lomelí 3D volume-preserving map, standard/Chirikov-Taylor map) is a
closed-form algebraic map — none is a genuine Poincaré return map defined
implicitly by continuous-time ODE flow to a section (e.g. CR3BP). Computing a
flow-defined Poincaré map's own high-order Taylor jet at each orbit point
requires handling an implicit, state-dependent return time — algebraically
awkward in a way it isn't for a closed-form map. This variable-return-time
complication, not "map-composition cost," is the actual gap Kumar 2025's
fixed-time flow-map/adapted-frame construction is built to sidestep. Also:
Kumar's own paper title ("multi-shooting parameterization method") is
literally Gonzalez & Mireles James's own coined framing (multiple shooting,
composition-free), generalized from single-orbit manifolds to Poincaré-map
periodic-orbit manifolds with m>1 crossings — an extension of, not a
departure from, this paper's own contribution.

### 5.3 Positive-control data (sourced, page-cited)

Hénon map (a=1.4, b=0.3), period-2 orbit p0=(−0.475800051175056,
0.292740015352517): unstable eigenvalue of Df² λ̄=−3.010100667740269,
eigenvectors given to 15 digits; K=50 unit-length-eigenvector computation:
0.004s, a posteriori error ε=1.33×10⁻¹⁵ (500 sample points). Rescaled
(×10, K=50/60): same 0.004s, ε≈10⁻⁷/10⁻¹². Rescaled (×22, K=110): ε=6.9×10⁻⁴,
described as "near the limit of what can be done...in double precision."
Period-95 orbit (K=50, tol=10⁻¹⁴). Lomelí 3D volume-preserving map
(a=0.5,b=−0.5,c=1,τ=1.333333333,α=0.3444444444): 1D/2D manifolds to order 25,
error <10⁻¹⁴. Standard map (a=2.1): period-4/25-orbit manifolds to order 200.
These give a reader a concrete, independently-checkable basis to reproduce
this paper's own runtime/accuracy claims (and thereby audit Kumar 2025's
comparative claims against them, per §5.2 above).

---

## 6. Mandatory citation-mining pass — consolidated across all five papers

Each paper's own background/related-work/reference-list section was read;
candidates were cross-checked against `docs/notes/CORPUS_INDEX.md` (grepped
directly, not inferred). **No acquisition was performed for any candidate
below** — flagged only, per policy.

### 6.1 Genuinely new gaps (not found in `CORPUS_INDEX.md`)

| Candidate | Priority | Why relevant |
|---|---|---|
| Burgos-García, J., Lessard, J.-P. & Mireles-James, J.D., "Spatial periodic orbits in the equilateral circular restricted four-body problem: computer-assisted proofs of existence," *Celestial Mech. Dyn. Astron.* 131(1):2 (2019) | **HIGH** | Directly CR4BP; rigorous (computer-assisted) existence proofs for spatial periodic orbits — a genuinely different validation standard (rigorous vs. numerical) from anything else in this project's own CCR4BP/CRNBP corpus. Cited by both Calleja 2021 [46/47] and independently surfaces from the Gonzalez & Mireles James 2017 lineage. |
| Castelli, R., Lessard, J.-P. & Mireles-James, J.D., "Rigorous numerics for the isolation and continuation of periodic orbits in the equilateral circular restricted four-body problem (II): a posteriori analysis and error bounds," *J. Dyn. Diff. Eq.* 30(4):1525-51 (2017) | **HIGH** | Companion/predecessor to the above — same CR4BP rigorous-numerics program, "Part II" implying a Part I also worth locating. |
| Hungria, A., Lessard, J.-P. & Mireles-James, J.D., "Rigorous numerics for analytic solutions of differential equations: the radii polynomial approach," *Math. Comput.* 85(299):1427-59 (2016) | Medium | The underlying rigorous-numerics (radii-polynomial) method the two CR4BP papers above build on. |
| Gonzalez, A., "Accurate High Order Computation of Invariant Manifolds for Long Periodic Orbits of Maps and Equilibrium States of PDE," PhD thesis, Florida Atlantic University (2020) | Medium | Same first-author/advisor lineage as the just-acquired Gonzalez & Mireles James 2017 SIADS paper (item 10 of this task) — plausibly the full thesis this paper's own results are drawn from, likely with additional worked examples/PDE extension. |
| Calleja, R. & de la Llave, R., "Computation of the breakdown of analyticity in statistical mechanics models," mp_arc preprint 09-56 (2009) | Medium | The actual numeric positive-control companion to the just-acquired Calleja & de la Llave 2010 Sobolev-seminorm breakdown-criterion paper (item 8) — see §3.2 above. |
| Calleja, R. & de la Llave, R., "Fast numerical computation of quasi-periodic equilibrium states in 1d statistical mechanics, including twist maps," *Nonlinearity* 22:1311-1336 (2009) | Medium | Companion to the above, same purpose. |
| Capella, A., "A parameterization of the invariant manifolds associated to periodic orbits in the RTBP" (cited by Cabré/Fontich/de la Llave 2005 as "in preparation," 2004) | Low-medium | Directly CR3BP-relevant if it was ever published — publication status not checked this pass; verify before acquiring. |
| Tantardini, M., Fantino, E., Ren, Y., Pergola, P., Gómez, G. & Masdemont, J.J., "Spacecraft trajectories to the L3 point of the Sun-Earth three-body problem," *Celest. Mech. Dyn. Astron.* 108:215-232 (2010) | Low-medium | CR3BP L3-point trajectory design; not core to current search axes but genuinely uncorpused. |
| Masdemont, J.J., "High-order expansions of invariant manifolds of libration point orbits with applications to mission design," *Dyn. Syst.* 20:59-113 (2005) | Low-medium | Classic libration-point manifold-expansion reference; possibly superseded in practice by already-acquired parameterization-method papers, but not independently confirmed as redundant. |
| Koon, W.S., Lo, M.W., Marsden, J.E. & Ross, S.D., "Low energy transfer to the moon," *Celest. Mech. Dyn. Astron.* 81:63-73 (2001) | Low | **Probable duplicate/precursor** of the already-acquired `koon-lo-marsden-ross-2000-shoot-the-moon.pdf`/KLMR 2000 heteroclinic-connections papers (same author group, closely adjacent topic/year) — not independently re-verified as genuinely distinct this pass; check before acquiring. |

### 6.2 Already-tracked, not double-counted

The following recurring citations across all five papers are **already**
flagged in `#730`'s own master list §5 (items 31-34) or already **acquired**
under `#733` — listed here only to confirm no new gap, not as new flags:
Fontich, de la Llave & Sire 2009 (§5 item 32); de la Llave, González, Jorba &
Villanueva 2005 (§5 item 33); Huguet, de la Llave & Sire 2012 (§5 item 34);
Haro & de la Llave 2003a/2003b preprints (superseded by the already-tracked
2006/2007 published versions); Haro, Canadell, Figueras, Luque & Mondelo
2016 book (**already ACQUIRED under `#733`**). The classic Koon/Lo/Marsden/
Ross/Gómez CR3BP-manifold papers cited by Gonzalez & Mireles James 2017
(Gómez et al. 2004 spatial RTBP, the 1999/2000/2001/2002 KLMR papers) are
**already in corpus** (`CORPUS_INDEX.md` lines 82-93) — confirmed by direct
grep, not re-flagged.

### 6.3 Low-priority background (not itemized individually)

Per the project's own §9 convention (background/textbook items flagged once,
low priority, not worth a full DOI-resolution pass): the extensive
renormalization-group/KAM-small-divisor/twist-map theory cited by Calleja &
de la Llave 2010 and Calleja et al. 2021 (Zehnder 1975/76 generalized
implicit-function theorems, Moser 1966/67/86, Kolmogorov 1954, Nash 1956,
Greene 1979's residue method, MacKay 1982/83/89/92, Shenker-Kadanoff 1982,
Chirikov 1979, Berretti-Chierchia/Berretti-Marmi complex-rotation-number
papers, Falcolini-de la Llave 1992, Mather 1982-87 nonexistence criteria,
del Castillo-Negrete's own nontwist-map series, Lichtenberg & Lieberman's
*Regular and Chaotic Dynamics* textbook, and similar) — none of these are
celestial-mechanics-specific and none were independently DOI-resolved this
pass.

---

## 7. Summary: no code fixes required

Across all five cross-checks: (1) Iuliano's equation defect is confirmed
textually but was already established as dynamically inert by `crnbp.py`'s
own existing proof — no action needed. (2) Calleja 2021's characterization as
Kumar 2026's "only directly-competing prior art" is fair with a documented
nuance (the Newton-Gauss refinement stage is itself not original to Calleja
2021) — informational, no code implication. (3) Calleja & de la Llave 2010's
Sobolev-breakdown criterion is confirmed as a genuine, reusable future
diagnostic, not yet needed. (4) Cabré/Fontich/de la Llave 2005's lineage is
confirmed consistent with the already-traced Haro/de la Llave ancestry, and
its parameterization method remains correctly un-adopted (tracked as future
capability, not a bug). (5) Gonzalez & Mireles James 2017's method is
**mischaracterized** by an already-acquired digest's summary of Kumar 2025
(the "Taylor-composition" framing describes Kumar's paper's own claimed
baseline-beaten-by-2017, not the 2017 paper's actual composition-free
contribution) — this is corrected here for the record; it does not affect
any code path (Kumar 2025's own manifold code, if/when built, is unaffected
by which prior paper is more precisely the baseline). No `CorpusAnchor`,
`literature_check.py`, or search-code changes are warranted by any of the
five.
