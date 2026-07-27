# Digest: Kumar 2025, "Multi-shooting parameterization methods for invariant manifolds and
heteroclinics of 2 DOF Hamiltonian Poincaré maps, with applications to celestial resonant dynamics"

**Paper:** arXiv:2509.03655 (submitted manuscript, single version v1). Full title as above (confirmed
against the arXiv abstract page, not just the task's paraphrase).
**Author:** Bhanu Kumar (Dept. of Mathematics, University of Michigan). **Solo-authored** — NOT
co-authored with Anderson/de la Llave as this task's dispatch assumed; the dispatch's author guess was
wrong, confirmed directly against the arXiv listing (`https://arxiv.org/abs/2509.03655`, author field
`Bhanu Kumar` only). Anderson and de la Llave are cited extensively (this is a direct sequel to
Kumar/Anderson/de la Llave's joint work) but are not authors of this specific paper.
**Submitted:** 3 Sep 2025. MSC 37M21, 37C27, 37C29, 70M20. 37 pages (35 in the downloaded PDF's own
page count; front matter/references pagination differs slightly from arXiv's reported count — not
significant).
**Filed:** `kumar-2025-multishooting-parameterization-invariant-manifolds-heteroclinics-poincare-maps-arxiv-2509.03655.pdf`
(private `cyclers_pdf` repo).
**Acquired/digested:** 2026-07-27 (`#728`).
**Text layer:** confirmed via `pdffonts` (embedded/subsetted Type-1 CM fonts, standard LaTeX/arXiv
output) and `pdftotext -layout` (clean extraction, verified against the visible page 1 abstract text).
No OCR needed.

## 1. What this paper actually presents

This is the **direct sequel** to Kumar, Anderson & de la Llave 2021 (*Commun. Nonlinear Sci. Numer.
Simul.* 97:105691, "High-order resonant orbit manifold expansions..." — cited as ref. [7], **NOT
currently in this project's corpus**, see §4 below) generalized from single-Poincaré-section-crossing
periodic orbits to **multiple-crossing** ones.

**The problem it solves:** given a 2-DOF Hamiltonian system (concretely: the planar CR3BP) and an
unstable periodic orbit γ that crosses a chosen Poincaré section Σ at *m > 1* points (not just once),
compute accurate, high-order (Taylor-series) parameterizations of γ's stable/unstable manifold curves
on Σ — i.e., treat γ as a period-*m* periodic orbit of the Poincaré return map P and parameterize its
manifolds — then use those parameterizations to find heteroclinic connections between two such orbits.

**Why m > 1 sections matter (the paper's actual motivating gap, §1-§2):** the classical choice of
Poincaré section for PCRTBP resonant-orbit work (fixed-x or fixed-y, e.g. Anderson & Lo 2010/2011) lets
unstable resonant orbits cross only once, but has poor transversality — it produces discontinuous
manifold-curve jumps where the flow is tangent to the section. A section transverse everywhere in the
region of interest (osculating true-anomaly ν = 0 "periapse" or ν = π "apoapse") fixes the
transversality problem but forces m:n resonant orbits to cross m times. Prior parameterization-method
machinery for maps (Gonzalez & Mireles James 2017, ref [9], **not in corpus**) required composing the
Poincaré map itself with Taylor series — expensive/awkward since P is defined by ODE propagation, not
closed-form algebra. **This paper's actual novel contribution** is avoiding that composition: it solves
an *intermediate* invariance equation using **fixed-time flow maps** Φ_τ(k) (the per-point first-return
times) instead of the Poincaré map P directly, via an "adapted frame" construction (§4.1, generalizing
Kumar/Anderson/de la Llave's own 2022 torus-adapted-frame machinery, ref [15], **also not in corpus**,
see §4) that near-decouples the order-by-order Taylor recursion into scalar/vector cohomological
equations solvable by a simple contraction-mapping fixed-point iteration (§4.1.1-4.1.2). Manifold curves
are computed to degree ~20-25 via automatic differentiation + jet transport (Julia,
`TaylorSeries.jl`/`TaylorIntegration.jl`/`OrdinaryDiffEq.jl`, DP8 integrator). Heteroclinic connections
between two orbits' manifolds are then found via a "layers" restriction (bounding the connection search
to a finite number of Poincaré-map applications, a rough time-of-flight proxy) + line-segment
intersection (GPU-parallelizable) + bisection refinement (§5).

## 2. Key numerical result (sourced, with page/section citation)

**Table 1 (§6.2.1, p.28), reproduced verbatim (fundamental-domain-size ratio, degree-20 series vs.
linear/eigenvector truncation, computed over Jacobi constants C=3.00 to 3.01 in steps of 0.0001, for
Uranus-Oberon exterior/interior MMR manifolds from the companion paper [13]):**

| MMR | Minimum | Maximum | Mean | Median |
|---|---|---|---|---|
| 3:4 (exterior) | 179.03 | 2562.78 | 808.37 | 525.38 |
| 4:5 (exterior) | 158.49 | 3841.14 | 1320.07 | 784.00 |
| 5:6 (exterior) | 160.97 | 4550.00 | 1533.86 | 496.83 |
| 4:3 (interior) | 295.18 | 1283.00 | 566.65 | 492.42 |
| 5:4 (interior) | 237.96 | 4268.00 | 1465.14 | 814.88 |
| 6:5 (interior) | 215.96 | 4248.82 | 1448.11 | 908.30 |

Text states the lowest single ratio across all cases is 158.49 and the **mean improvement ratio across
all computed manifolds is 1190.37** (p.28-29) — i.e. the degree-20 nonlinear parameterization's
fundamental domain of validity is on average ~1190x wider than the linear (single-order/eigenvector)
truncation's domain. Domains of the degree-20 series were "mostly on the order of 0.1, with a few on the
order of 0.01"; linear-truncation domains were "generally on the order of 10⁻⁴ to 10⁻⁵" (p.28). This
directly means: less numerical integration is required to globalize the manifold, and — implicitly, not
separately quantified in the paper — a wider, better-conditioned neighborhood around the base orbit
before nonlinear effects invalidate the local approximation.

**Performance (p.26, 2019 Mac laptop, Intel i9 8-core + AMD Radeon GPU):** adapted-frame + degree-20
manifold + fundamental-domain computation < 5 s/manifold; globalization/visualization step < 30 s
(often < 15 s); GPU-parallelized (OpenCL.jl) line-segment intersection search < 0.1 s; per-intersection
bisection refinement 0.1-2 s.

**Demonstrated applications (§6, reviewing the author's own prior conference/journal papers [12]-[14],
all EXCEPT [13] already in this project's corpus):** Uranus-Oberon 3:4/4:5/5:6 exterior and 4:3/5:4/6:5
interior MMR manifolds (Figs. 2-3, from ref [13], the AAS/AIAA 2024 Oberon survey — **not in corpus**,
see §4); Earth-Moon 4:1/3:1/2:1 interior MMR manifolds and a computed 3:1→2:1 heteroclinic connection
trajectory at C=3.05 (Figs. 4-6, from refs [12] Rawat/Kumar/Rosengren/Ross 2025 arXiv:2505.10138 — the
same paper already digested in this corpus as `2026-06-30-digest-rawat2026-cislunar-mmr.md`, `#503` — and
[14] Kumar/Rawat/Rosengren/Ross 2024 IAC, already digested `#597`).

## 3. Relevance assessment to this project's own codebase

**Grep evidence gathered before writing this verdict** (per this task's mandate): `cr3bp_periodic.py`
(single-shooting STM Newton corrector for periodic orbits), `cr3bp_multiple_shooting.py` (segmented STM
Newton, #687, built for numerically-ill-conditioned long arcs — a periodicity corrector, not a manifold
tool), `qp_torus_manifold.py` (torus manifold direction via a single-point STM/Floquet eigendecomposition
— linear only), `ccr4bp_whisker.py`/`ccr4bp_manifold_globalize.py` (CCR4BP torus manifold direction via
segment-anchored discrete-QR/CLV eigendecomposition — still a LINEAR direction, just more robustly
extracted — then nonlinear-flow globalization), and **`resonance_network.py` (#267)**.

**The single most concrete, directly-grounded finding of this digest:** `resonance_network.py` — this
project's own "Tier 3" MMR heteroclinic-network prioritizer, explicitly built on Kumar, Rawat, Rosengren
& Ross 2025 (arXiv:2509.12675, already in corpus) — computes Earth-Moon **4:1/3:1/2:1 unstable resonant
periodic orbit manifolds** (`_stable_unstable_manifold_arc`, `unstable_eigenvector` at
`src/cyclerfinder/search/resonance_network.py:300-320,440-469`) via the **exact linear/eigenvector
method** this paper's Table 1 quantifies as the inferior baseline: a monodromy-matrix eigenvector
(`np.linalg.eig` on the one-period 4x4 planar STM) plus a fixed `epsilon` perturbation, forward/backward
propagated by the true nonlinear flow to a perigee/perilune Poincaré section
(`_perilune_event`/`_perigee_event`). This is precisely the object class (Earth-Moon interior MMR
periodic orbits, periapse/apoapse Poincaré section, manifold-tube overlap for heteroclinic detection)
this very paper's §6.3 demonstrates its own method on — and the same object class its OWN cited sibling
papers (already in this project's corpus) use to seed `resonance_network.py`'s `_RESONANT_SEEDS`. If this
project ever wants to compute its own 4:1/3:1/2:1 (or other MMR) manifolds directly, rather than scoring
proximity to literature-sourced Jacobi constants and coarse eigenvector perturbations,
**this paper's method is a concrete, well-specified upgrade path**: replacing `resonance_network.py`'s
degree-1 (linear) manifold representation with this paper's degree-~20 adapted-frame parameterization
would (per Table 1, same MMR-orbit class) plausibly widen the valid perturbation-epsilon domain by
~2-3 orders of magnitude and reduce dependence on `_perilune_event`'s numerical-integration-heavy
globalization for accuracy near the base orbit. This is a genuine, scoped future build candidate, not
speculative hand-waving — it targets a module that already exists and already implements the weaker
version of the same idea.

**Cross-check against the documented `#606`-`#620` seedless/variational-corrector arc (per this task's
mandate — read via `data/OUTSTANDING.md` #606-#646):** **no direct bearing.** That arc's wall
(`#619`/`#646`: EM-L2 unstable-manifold-direction extraction failing/fragile, later fixed by a
segment-anchored CLV extractor, but the resulting connection search still hit a ~166,000 km closure
floor even with a verified-correct direction) concerns the manifold of a **quasi-periodic 2-TORUS**
(a GMOS object from `genome/qp_tori.py`/`search/variational_qbcp_torus.py`), not a periodic orbit. This
paper's method is formally scoped to **periodic orbits of 2-DOF Poincaré maps** — it does not address
torus manifolds at all (the torus-manifold analog is Kumar/Anderson/de la Llave's separate 2022 CMDA
paper, ref [15], cited repeatedly here but itself **not in this project's corpus** — see §4). It would be
speculative and unverified to claim this paper "would have helped" #619/#646; the honest assessment is
that it is the wrong object class for that specific wall. What IS true, and worth flagging: #646's own
diagnosis was that the closure failure persisted even after the manifold DIRECTION was fixed — i.e. the
residual problem may be more than direction accuracy (a genuine non-existence, or higher-order effects
beyond the linear direction the CLV extractor still only reports). This paper's nonlinear domain-of-
validity argument is at least suggestive that a fully nonlinear (not just direction-corrected-linear)
torus/orbit manifold representation could still matter for such walls in general — but confirming that
would require Kumar/Anderson/de la Llave's actual TORUS parameterization paper (uncorpused, see below),
not this periodic-orbit-specific one, and is speculative, not demonstrated.

**Verdict: grounding + one concrete, scoped upgrade opportunity (`resonance_network.py`'s manifold
fidelity), not a fix for any currently-open wall.** Honest, not inflated: this paper does not solve a
documented blocker; it identifies a real quality gap in an already-built, already-corpus-grounded module.

## 4. Citation-mining pass (mandatory; cross-checked against `docs/notes/CORPUS_INDEX.md`, grep
run 2026-07-27, and the papers/ directory listing directly)

Full reference list ([1]-[29]) read. Flagged gaps, in priority order:

1. **Kumar, Anderson & de la Llave 2022**, "Rapid and accurate methods for computing whiskered tori and
   their manifolds in periodically perturbed planar circular restricted 3-body problems," *Celestial
   Mechanics and Dynamical Astronomy* 134(1):3, DOI `10.1007/s10569-021-10057-1` [ref 15]. **HIGH
   priority.** This is the actual TORUS-parameterization-method paper this entire manuscript's
   adapted-frame machinery is generalized FROM (explicitly: "adapted from a similar procedure developed
   in [15] for invariant tori"), and the paper `qp_torus_manifold.py` cites Kumar/Anderson/de la Llave
   2025 for torus stability TYPE but not this 2022 paper's actual nonlinear manifold-parameterization
   method. **Already independently flagged as a "possible-distinct-candidate" (not yet acquired) in the
   `#727` digest** (`2026-07-27-727-...-digest.md` §5, item [21]) — this is a second, independent hit on
   the same gap, raising its priority.
2. **Kumar, Anderson & de la Llave 2021**, "High-order resonant orbit manifold expansions for mission
   design in the planar circular restricted 3-body problem," *Commun. Nonlinear Sci. Numer. Simul.*
   97:105691, DOI `10.1016/j.cnsns.2021.105691` [ref 7]. **HIGH priority** — this is the direct
   single-intersection PREDECESSOR paper this manuscript generalizes; foundational to understanding what
   is actually new here (the m>1 generalization) vs. inherited.
3. **Gonzalez & Mireles James 2017**, "High-order parameterization of stable/unstable manifolds for long
   periodic orbits of maps," *SIAM J. Applied Dynamical Systems* 16(3):1748-1795 [ref 9]. Medium-high
   priority — the prior-art method (Taylor-composition with the map directly) this paper explicitly
   improves on/avoids; useful for understanding the actual novelty claim's baseline.
4. **Kumar & Anderson 2024**, "A Survey of Oberon Mean Motion Resonant Unstable Orbit Properties and
   Connections for Uranian Tours," AAS/AIAA Astrodynamics Specialist Conference [ref 13]. Medium
   priority — the source of this paper's own headline Uranus-Oberon figures (Figs. 2-3, Table 1's data).
   Conference-only, no arXiv ID found — no e-print source path per this project's acquisition policy, so
   flagging only, not a quick-win acquisition.
5. **Haro, Canadell, Figueras, Luque & Mondelo 2016**, *The Parameterization Method for Invariant
   Manifolds: From Rigorous Results to Effective Computations*, Applied Mathematical Sciences vol. 195,
   Springer [ref 3]. Medium priority — the foundational parameterization-method textbook cited by nearly
   every Kumar/Anderson/de la Llave paper in this corpus; **already flagged once** in the `#727` digest
   (medium priority, ref [18] there) — third independent citation of this exact gap across this project's
   Kumar-lineage digests (also implicitly load-bearing for paper 2 of this task, see its own digest).
6. **Cabré, Fontich & de la Llave 2005**, "The parameterization method for invariant manifolds III:
   overview and applications," *Journal of Differential Equations* 218(2):444-515 [ref 2]. Low-medium
   priority — core parameterization-method theory reference, not yet in corpus.

**Already in corpus, no action:** [12] Rawat/Kumar/Rosengren/Ross 2025 (arXiv:2505.10138, digested
`#503`), [14] Kumar/Rawat/Rosengren/Ross 2024 IAC (digested `#597`), [28] Kumar/Anderson/de la Llave 2025
SIAM ADS 24:219 GPU connections (digested, whiskered-tori-connections.md). General-theory/textbook
references not domain-specific to this project (Morbidelli, Bate/Mueller/White, Celletti, Chicone,
Thirring, Chirikov x2) not flagged as acquisition gaps — standard background texts, not novel results.

## Summary for the dispatching session

Method: multi-shooting-style adapted-frame parameterization avoiding Poincaré-map/Taylor-series
composition, generalizing single-crossing manifold parameterization to multi-crossing periodic orbits;
demonstrated ~1190x mean fundamental-domain improvement over linear eigenvector approximations on
Uranus-Oberon and Earth-Moon MMR manifolds (Table 1, p.28). Relevance: no fix for the documented
`#606`-`#646` torus-manifold wall (wrong object class — periodic orbits, not tori), but a concrete,
scoped upgrade candidate for `resonance_network.py`'s own linear-eigenvector manifold construction on
the SAME Earth-Moon MMR orbit class the module already targets. Citation-mining surfaced a coherent,
repeatedly-cited cluster of foundational parameterization-method papers (Haro et al. 2016 book, the 2022
CMDA whiskered-tori paper, the 2021 CNSNS single-intersection predecessor, Cabré/Fontich/de la Llave
2005) that underlie multiple already-corpused Kumar-lineage application papers but are themselves absent
from the corpus — worth a dedicated follow-on acquisition pass given how load-bearing they are across
this whole research thread.
