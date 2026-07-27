# Digest: Kumar 2026, "A new fast multiple-shooting method for computing periodic orbits in symplectic
maps leveraging simultaneous Floquet vector computation to avoid large linear systems"

**Paper:** arXiv:2601.00149 (preprint v1). Full title as above (confirmed against the arXiv abstract
page). **Note on arXiv ID:** the `2601.xxxxx` prefix does denote a January-2026 submission and the page
DOES resolve (confirmed live: `[v1] Wed, 1 Jan 2026`, `arXiv:2601.00149v1 [math.DS] 1 Jan 2026` printed
on the PDF itself) — this is a genuine, currently-indexed preprint, not a placeholder or 404.
**Author:** Bhanu Kumar (Dept. of Mathematics, University of Michigan). **Solo-authored** — like the
companion paper 1 of this task, NOT co-authored with Anderson/de la Llave as the dispatch's title guess
assumed; confirmed directly against the arXiv author field.
**Subjects:** math.DS (primary), astro-ph.EP, nlin.CD. 2020 MSC: 37C25, 37C27, 37J12, 70H12. 33 pages
(PDF page count; downloaded file reports 10 physical pages at the `file` tool level due to dense
multi-column-style LaTeX layout — text extraction confirms the full 33-page content is present).
**Filed:** `kumar-2026-fast-multishooting-periodic-orbits-symplectic-maps-floquet-arxiv-2601.00149.pdf`
(private `cyclers_pdf` repo).
**Acquired/digested:** 2026-07-27 (`#728`).
**Text layer:** confirmed via `pdffonts` (embedded/subsetted Type-1 CM + Concrete/custom fonts) and
`pdftotext -layout` (clean extraction). No OCR needed.

## 1. What this paper actually presents

**Setting:** a 4D symplectic map family F_ε (ε ≥ 0 a perturbation parameter) whose unperturbed map F₀
has a 2D cylindrical normally hyperbolic invariant manifold (NHIM) foliated by 1D invariant tori — the
canonical example being the **stroboscopic map of a 2.5-DOF periodically-perturbed 2-DOF Hamiltonian
flow** (e.g. the CCR4BP as a periodic perturbation of the PCR3BP; §3.2.1 derives this explicitly, citing
Blazevski & Ocampo 2012 for the CCR4BP equations of motion — **not in this project's corpus**, already
flagged in the `#727` digest, see §4). Tori with rational rotation number p/q are foliated by
**subharmonic periodic orbits (SPOs)** — q-iteration periodic orbits of F₀ — of which generically only a
discrete subset persist as ε > 0 (per subharmonic Melnikov theory or symmetry arguments).

**The problem it solves:** given a persisting F₀-SPO (q points), compute the corresponding Fε-SPO
**and simultaneously its Floquet directions/multipliers** (the linearized stable/unstable/center
structure at every orbit point), efficiently enough to matter when q is large (order 100+, as needed for
studying torus breakdown via secondary-resonance overlap — the paper's own motivating application, ref
[3], Kumar/Anderson/de la Llave 2023 AAS-23-397, already in this project's corpus, digested `#688`).

**Why this is hard the traditional way (§1, explicit complexity argument):** classical multiple shooting
for a q-point SPO in 4D solves a dense 4q×4q Newton linear system — O(q³) per iteration — and a
SEPARATE post-hoc eigendecomposition of the resulting 4q×4q monodromy matrix for stability — also O(q³).
For q~100 this becomes expensive and numerically fragile (the paper notes eigenvector extraction from a
large, highly non-normal monodromy matrix is often inaccurate at large q — directly resonant with this
project's own `#619`/`#646` finding that a one-shot eigendecomposition is unreliable for strongly
hyperbolic long-period objects, see §3 below).

**This paper's actual novel contribution (§4, confirmed by direct comparison — the ONLY genuinely new
idea here, everything else is infrastructure to support it):** adapt Kumar/Anderson/de la Llave's own
invariant-TORUS parameterization method (ref [1], the 2022 CMDA "whiskered tori" paper — **also not in
this project's corpus**, see §4) to solve for the SPO points **and** an adapted Floquet frame P_ε(k) with
a near-diagonal Floquet matrix Λ_ε(k) (Eq. 14-15) **simultaneously**, in one coupled but *decoupled-by-
construction* quasi-Newton system. Because Λ_ε(k) is near-diagonal (block form: a 2x2 tangent-to-NHIM
block with one off-diagonal term T, plus separate scalar stable/unstable Floquet multipliers λ_s(k),
λ_u(k)), each quasi-Newton correction step reduces to solving a sequence of **scalar or small linear
"cohomological equations"** of the recurring form `λ_a(k)u(k) - λ_b(k)u(k+1 mod q) = b(k)` (§4.3.1-4.3.3),
solvable either by an explicit closed-form summation (Eq. 19) or by contraction-mapping fixed-point
iteration (guaranteed convergent whenever `|λ_a/λ_b| ≠ 1` in the appropriate direction) — **O(q) work per
correction step, not O(q³)**, and Floquet directions/multipliers fall out of the same solve rather than
requiring a separate eigendecomposition. The paper states this explicitly (Abstract, §7): "the first-ever
[application] of the [torus] parameterization method framework to directly compute periodic orbit points
themselves," improving efficiency over "prior multi-shooting methods for SPOs" — the only cited prior
work in this exact niche is Calleja, del-Castillo-Negrete, del Río & Olvera 2021 (ref [10], **not in
corpus**, see §4), which the paper says only found a 1D curve containing the SPO via a parameterization-
style step, not the actual orbit points, still requiring a traditional large-matrix multi-shooting step
afterward.

A second, dependent contribution (§5) is a **recursive Taylor-series parameterization** of the SPO's weak
stable/unstable **separatrices**, initialized directly from the simultaneously-computed Floquet
directions — this reuses the SAME jet-transport/automatic-differentiation machinery (TaylorSeries.jl /
TaylorIntegration.jl / OrdinaryDiffEq.jl, DP8 integrator) as paper 1 of this task and the same
"cohomological equation" solvers.

## 2. Key numerical results (sourced, with section/page citations; no invented numbers)

Applied in **§6** to two already-published Kumar/Anderson-lineage CCR4BP studies (both already in this
project's corpus, digested):

- **Jupiter-Ganymede-Europa CCR4BP (§6.3, ref [3] = AAS-23-397, `#688`):** Ganymede 4:3 MMR unstable
  orbit family, Europa as perturber (µ₃ = ε, physical value 2.5265×10⁻⁵). SPOs continued for ω/2π ratios
  11/34, 34/105, 23/71, 35/108, 12/37, 25/77 (from [3]) plus, newly reported here, **37/114 and 45/139**
  (not in the original AAS-23-397 paper — a genuinely new result of this manuscript, though presented
  only as additional plotted points, no separate numeric table). Tolerance 10⁻⁷ on the invariance
  equations, continuation step Δε = 5×10⁻⁷ to 10⁻⁶. **All of these SPOs were successfully continued to
  the full physical ε, in contrast with the corresponding tori (none of which survive to ε = 2.5265×10⁻⁵
  for ω < 2.04047)** — this is the physical point of the whole exercise: SPOs persist where tori don't,
  and separatrices computed from their Floquet directions (Fig. 5) confirm consecutive SPOs' separatrices
  intersect, explaining the torus-destruction mechanism. One specific SPO example: the 23/71 SPO's λ₁-λ₂
  Floquet pair flips hyperbolic↔elliptic between ε=3×10⁻⁶ and ε=4×10⁻⁶ (p.24) — handled without accuracy
  loss by the simultaneous method.
- **Uranus-Titania-Oberon CCR4BP (§6.3, ref [27] = Kumar & Anderson 2024 AAS-24-288, "survey of Oberon
  MMR unstable orbit properties" — **not in corpus**, see §4):** Oberon 6:5 MMR family, Titania as
  perturber (µ₃ = ε, physical value 3.9168×10⁻⁵, µ = 3.5433×10⁻⁵ for Uranus-Oberon). SPOs computed for
  ω/2π = 25/69, 21/58, 17/47, 30/83, 13/36, 22/61, 9/25 — same qualitative torus-destruction-via-
  separatrix-intersection mechanism as the Jovian case (Fig. 7).
- **Performance (§6.2, 2021 Apple M1 Pro laptop, Julia):** SPO continuation from the PCR3BP through to
  the full CCR4BP generally took **< 2 seconds per SPO**. Separatrix computation used Taylor truncation
  order d = 20, same jet-transport toolchain as paper 1.

No table analogous to paper 1's Table 1 (fundamental-domain ratio) appears in this paper — its headline
efficiency claim (O(q) vs O(q³)) is an **asymptotic complexity argument** (§1, explicit), not an
empirically tabulated speedup number; no q~100 wall-clock comparison against a literal dense-matrix
baseline is reported. This is a real gap in the paper's own evidence for the claimed practical speedup
(worth noting honestly rather than treating the O(q) vs O(q³) claim as independently benchmarked here).

## 3. Relevance assessment to this project's own codebase

**Grep evidence:** this project has an active CCR4BP model+search stack (`core/ccr4bp.py` and five
per-moon-pair variants, `search/ccr4bp_whisker.py`, `search/ccr4bp_manifold_globalize.py`,
`search/variational_ccr4bp_torus.py`, `search/ccr4bp_heteroclinic_search.py`,
`search/ccr4bp_chained_transfer.py`) grounded directly in this same Kumar/Anderson/de la Llave CCR4BP
paper lineage (the Umbriel-Titania torus-homoclinic connection, per
`[[project_novel_findings_status]]`, is this project's own confirmed novel CCR4BP finding). All of the
existing CCR4BP manifold/whisker machinery I inspected (`ccr4bp_whisker.py`'s
`manifold_direction_segmented_clv`, `ccr4bp_manifold_globalize.py`) computes manifolds of **quasi-
periodic 2-TORI** (via `variational_ccr4bp_torus.py`'s seedless pseudospectral corrector, the CCR4BP
generalization of the `#606`-`#620` seedless-corrector arc), extracting a linear stable/unstable
direction via a segment-anchored discrete-QR/covariant-Lyapunov-vector (CLV) eigendecomposition and then
globalizing by nonlinear-flow propagation from a small offset. **I found no existing module in this
project that computes SUBHARMONIC PERIODIC ORBITS of a CCR4BP stroboscopic map at all** — this
project's periodic-orbit correctors (`cr3bp_periodic.py`, `cr3bp_multiple_shooting.py`) target CR3BP
periodic orbits, not CCR4BP stroboscopic-map SPOs, and its CCR4BP work has gone straight to tori, never
building the intermediate SPO object this paper's method computes.

**Is this a genuine capability gap worth building?** Conditionally yes, but scoped narrowly: this
paper's method matters specifically when (a) the CCR4BP torus of interest has a **resonant/rational**
rotation number under the stroboscopic map (i.e. it does NOT persist as a torus at all, only as isolated
SPOs — exactly the secondary-resonance-overlap regime the already-corpused AAS-23-397/AAS-24-288 papers
document) and (b) q is large enough (order tens-hundreds) that a naive dense multiple-shooting Newton
system would be a real O(q³) cost. This project's existing CCR4BP torus tooling (built for the
Umbriel-Titania discovery and the broader `#686`-`#727` arc) implicitly assumes the torus PERSISTS — it
has no machinery for the case where it doesn't, which this paper's own §6.3 example shows is common
(none of the ω < 2.04047 tori in the Jupiter-Ganymede-Europa case survive to physical ε). **If this
project's discovery program ever wants to systematically search the secondary-resonance/SPO regime of a
CCR4BP system it already models** (e.g. extending the existing Umbriel-Titania or Jupiter-moon CCR4BP
work into the torus-breakdown zone rather than only the persisting-torus zone), this paper's decoupled
O(q) simultaneous SPO+Floquet method — not a dense-matrix multiple shooting reimplementation — is the
right tool, and is a well-scoped, medium-effort future build (the cohomological-equation solvers §4.3.1-
4.3.4 are the same contraction-mapping machinery paper 1 of this task also uses, so building both papers'
methods together would share most of the implementation).

**Cross-check against the documented `#606`-`#646` seedless/variational-corrector arc:** as with paper 1,
**no direct fix for the EM-L2/SE-L2 basin-wall or manifold-conditioning negatives** — those concern
quasi-periodic TORI in the CR3BP/QBCP (not CCR4BP SPOs), and this paper's object class (subharmonic
periodic orbits of a *perturbed* map, specifically requiring a genuine perturbation parameter ε and an
already-known unperturbed-map persisting orbit to continue from) does not match the EM-L2/SE-L2 problem
setup (a single-system torus continuation/connection search, not a perturbative-family SPO continuation).
One methodologically interesting resonance worth naming honestly: this paper's §1 complaint about O(q³)
eigendecomposition being "inaccurate when q is large due to the orbit's instability" is the SAME failure
mode `#619` diagnosed for its one-shot STM eigendecomposition (fixed by `#646`'s segment-anchored CLV
approach) — but this paper solves it with a completely different mechanism (decoupled algebraic
cohomological equations from a near-diagonal Floquet ansatz, not segment-wise QR re-orthonormalization).
This is a genuinely alternative numerical strategy for the same class of problem (robust Floquet
extraction on long/unstable orbits), worth keeping in mind if `#646`'s CLV approach is ever found to have
its own limits on some future long-period target — but it is NOT a demonstrated fix for anything
currently open, since it was never tried against that specific problem.

**Verdict: grounding + a genuine, narrowly-scoped future capability gap (CCR4BP SPO/secondary-resonance
computation, currently entirely absent from this project), not a fix for any currently-open wall.**

## 4. Citation-mining pass (mandatory; cross-checked against `docs/notes/CORPUS_INDEX.md`, grep run
2026-07-27, and the papers/ directory listing directly)

Full reference list ([1]-[31]) read. Flagged gaps, in priority order (several overlap with paper 1 of
this task's own gap list, since both papers share the same core theoretical lineage — noted where
duplicate):

1. **Kumar, Anderson & de la Llave 2022**, "Rapid and accurate methods for computing whiskered tori and
   their manifolds in periodically perturbed planar circular restricted 3-body problems," *CMDA*
   134(1):3, DOI `10.1007/s10569-021-10057-1` [ref 1]. **HIGH priority — the single most load-bearing
   citation in THIS paper** (explicitly: "heavily inspired by the parameterization method of [1]...much
   of the discussion and proofs in this section follow a very similar structure to those of that paper").
   Same gap independently flagged by paper 1 of this task (its ref [15]) and by the `#727` digest
   (possible-distinct-candidate, ref [21] there) — **third independent hit**, now the highest-confidence
   acquisition-gap finding across this whole Kumar-lineage citation-mining effort.
2. **Calleja, del-Castillo-Negrete, del Río & Olvera 2021**, "A new method to compute periodic orbits in
   general symplectic maps," *Commun. Nonlinear Sci. Numer. Simul.* 99:105838 [ref 10]. **HIGH
   priority** — this is the ONLY directly-competing prior-art method this paper positions itself against
   ("the only prior work applying related tools to periodic orbits... only used a parameterization
   method-style algorithm to find a 1D curve containing the long periodic orbit... followed by a
   traditional large matrix-based multiple shooting scheme"); essential for verifying this paper's own
   novelty claim independently rather than taking its word for it.
3. **Haro, Canadell, Figueras, Luque & Mondelo 2016**, *The Parameterization Method for Invariant
   Manifolds*, Applied Mathematical Sciences vol. 195, Springer [ref 5]. Medium priority — the
   foundational textbook, same gap flagged in paper 1's digest and the `#727` digest (now three
   independent hits across this Kumar-lineage citation-mining pass).
4. **Blazevski & Ocampo 2012**, "Periodic orbits in the concentric circular restricted four-body problem
   and their invariant manifolds," *Physica D* 241(13):1158-1167 [ref 17]. Already flagged HIGH priority
   in the `#727` digest (the CCR4BP model-definition paper this project's own `core/ccr4bp.py` traces
   back through the Kumar papers) — this paper cites it for the identical purpose (Eq. 9 derivation);
   fourth-ish independent confirmation this is worth acquiring.
5. **Cabré, Fontich & de la Llave 2005**, "The parameterization method for invariant manifolds III,"
   *J. Differential Equations* 218(2):444-515 [ref 11]. Low-medium priority — same gap as paper 1's [2],
   core theory reference.
6. **Kumar & Anderson 2024**, "A survey of Oberon mean motion resonant unstable orbit properties and
   connections for Uranian tours," AAS/AIAA Astrodynamics Specialist Conference, AAS 24-288 [ref 27].
   Medium priority — same gap as paper 1's [13] (matching title, now with the specific AAS number);
   source of this paper's own Uranus-Titania-Oberon §6.3 example. Conference-only, no arXiv ID found.
7. **De la Llave, González, Jorba & Villanueva 2005**, "KAM theory without action-angle variables,"
   *Nonlinearity* 18(2):855-895 [ref 6]. Low-medium priority — KAM-persistence theory underlying the
   torus-persistence argument in §2.2/§3.3.
8. **Fontich, de la Llave & Sire 2009**, "Construction of invariant whiskered tori by a parameterization
   method. Part I," *J. Differential Equations* 246(8):3136-3213 [ref 21]. Low-medium priority — another
   core parameterization-method theory paper in the same uncorpused cluster as items 1, 3, 5, 7.
9. **Haro & de la Llave 2006**, "A parameterization method for the computation of invariant tori and
   their whiskers in quasi-periodic maps: Numerical algorithms," *Discrete & Continuous Dynamical
   Systems - B* 6(6):1261-1300 [ref 7]. Low-medium priority — same cluster.

**Already in corpus, no action:** Kumar/Anderson/de la Llave 2023 AAS-23-397 (ref [3], digested `#688`),
Chirikov 1960 (ref [4], general theory, previously assessed as not needing separate flagging).
Not domain-specific / low relevance, not flagged: Fenichel 1971 and Hirsch-Pugh-Shub 1977 (NHIM
persistence theory, standard references), Guckenheimer-Holmes 1983 and Guckenheimer-LaMar 2007 (general
dynamical-systems texts), Horn-Johnson 1985 (matrix analysis textbook), Treschev-Zubelevich 1998, Celletti
2010 (already effectively covered — cited identically in paper 1, assessed there as standard background),
Chicone 2006 (contraction mapping theorem reference, standard), the TaylorSeries.jl/TaylorIntegration.jl/
DifferentialEquations.jl software citations (tool references, not domain results), Pérez-Palau/Masdemont/
Gómez 2015 jet-transport paper [ref 24] (same as paper 1's [22], flagged there — not re-flagged here to
avoid duplication), Rasotto et al. 2016 and Berz-Makino 1998 (differential-algebra/verified-integration
tool references, not domain-specific), Spreen/Howell/Davis 2017 NRHO conference paper [ref 31] (mentioned
only in passing as a future-extension motivator, low priority, not flagged).

## Summary for the dispatching session

Method: a quasi-Newton multi-shooting scheme that computes subharmonic periodic orbits of a perturbed 4D
symplectic (e.g. CCR4BP stroboscopic) map SIMULTANEOUSLY with their Floquet directions/multipliers, by
adapting the torus parameterization method's near-diagonal Floquet-matrix ansatz — O(q) per correction
step vs. the traditional O(q³) dense-matrix multiple shooting + separate eigendecomposition. Demonstrated
on already-corpused Jovian (Ganymede 4:3 MMR / Europa perturber) and Uranian (Oberon 6:5 MMR / Titania
perturber) CCR4BP secondary-resonance studies, with SPO continuation < 2 s/orbit and separatrix
computation reusing the same jet-transport toolchain as paper 1. Relevance: this project has an active
CCR4BP model+search stack but NO subharmonic-periodic-orbit machinery at all — it currently only builds
CCR4BP tori (which fail to exist exactly in the secondary-resonance regime this paper targets). A genuine,
narrowly-scoped future capability gap, not a fix for the documented `#606`-`#646` torus-manifold wall
(different object class, different problem setup). Citation-mining reconfirms (now via a THIRD
independent hit across this session's two digests plus the prior `#727` digest) that the foundational
Kumar/Anderson/de la Llave 2022 CMDA whiskered-tori paper and the Haro et al. 2016 parameterization-method
textbook are the two highest-priority acquisition gaps underlying this entire corpus thread.
