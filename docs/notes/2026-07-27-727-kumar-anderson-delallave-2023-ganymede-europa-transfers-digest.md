# Digest: Kumar, Anderson & de la Llave 2023 (Acta Astronautica 211:76-87)

**Paper:** "Transfers between Jupiter-Ganymede and Jupiter-Europa Resonant Tori In a Concentric
Circular Restricted 4-Body Model"
**Venue:** Acta Astronautica 211 (2023), pp. 76-87, DOI `10.1016/j.actaastro.2023.05.040`. Open-access
author manuscript, Elsevier user license (ScienceDirect PII S0094576523002813). A previous version was
presented at IAC 2022 (IAC-22-C1,8,4,x73382).
**Authors:** Bhanu Kumar (Georgia Tech), Rodney L. Anderson (JPL/Caltech), Rafael de la Llave (Georgia
Tech).
**Filed:** `kumar-anderson-delallave-2023-transfers-ganymede-europa-resonant-tori-ccr4bp-gpu-manifold-intersections-acta-astro-211-doi-10.1016-j.actaastro.2023.05.040.pdf`
(private `cyclers_pdf` repo, commit `231a9d1`).
**Acquired/digested:** 2026-07-27 (`#727`), user-supplied PDF. Native text-layer PDF (LaTeX-produced,
embedded/subsetted Type-1 + TrueType fonts confirmed via `pdffonts`); no OCR needed.

**Why in corpus:** this is Kumar et al. 2021's (`kumar-anderson-delallave-gunter-2021-...-AAS-21-651`,
digested `#688`) own executed follow-up — the paper that actually searches for transfers between the
Jupiter-Ganymede and Jupiter-Europa resonant tori that 2021 paper computed. It was already flagged as
an acquisition candidate in `2026-07-27-724-final-confirmation-n5-torus-novelty.md` §3.3 ("This paper
is NOT in corpus... flagged as an acquisition candidate"). Companion/sibling Kumar papers already in
corpus: AAS 21-651 (2021, model + tori, `#688`), SIADS 24:219 / arXiv:2109.14814 (2025, whiskered-torus
GPU connection machinery in the *planar CR3BP* single-perturber case, digested
`2026-07-03-...whiskered-tori-connections.md`), AAS 23-397 / arXiv:2309.06073 (2023, secondary-resonance
overlap inside the JG 4:3 family, digested `#688`). This new paper is the CCR4BP-specific application of
the GPU connection-search machinery to the actual Ganymede-Europa transfer problem.

## 1. The GPU-assisted near-intersection method — how it works, what it's good for

Builds directly on the authors' own prior planar-CR3BP GPU method (arXiv:2109.14814, already digested)
but adapted/optimized specifically for the CCR4BP's higher-dimensional (4D stroboscopic-map, vs 2D for
plain CR3BP) manifold representation:

- Each unstable resonant torus's stable/unstable manifold is represented as a parameterized function
  `W(θ,s): T×R → R^4` (Fourier-Taylor series in θ, computed via the quasi-Newton method of their 2022
  CMDA paper [21] extended in this paper to **period maps** as well as stroboscopic maps — §3.4's new
  contribution, needed because period-map tori are cheaper to discretize in some regimes and give a
  cross-check against stroboscopic-map continuation stalls).
- The manifold is discretized on a grid of (θ,s) into a mesh of **quadrilaterals**, each split into 2
  triangles — a mesh in the full 4D `(x,y,px,py)` stroboscopic-map phase space.
- Two manifolds' meshes (an unstable W1^u and a stable W2^s) are tested for intersection by literally
  checking triangle-pairs for a 4D intersection (a 4x4 linear solve + inequality check) — conceptually
  simple, but with tens of thousands of triangles per mesh this is billions of triangle pairs
  (9.7e9 quad pairs in their own benchmark), intractable brute-force.
- **The actual novel contribution of this paper (§4.2.1, "the major improvement... during this study"):**
  a **uniform-grid spatial partition** over the 4D phase space (inspired by k-d trees, borrowed from
  computer-graphics collision detection) that only tests triangle pairs sharing a common grid cell —
  O(A+B) construction instead of O(A·B) all-pairs, cutting a benchmark run from requiring a JPL DGX
  V100 GPU cluster down to 10-26 seconds on a **consumer laptop GPU**. Followed by a cheap axis-aligned
  bounding-box test, then a Möller triangle-triangle test, then only the (now small, "tens of thousands
  at most") survivors get the exact 4x4 intersection solve.
- Restricted further to physically-relevant "layer" pairs (Un vs Sn, Sn-1) per their prior paper's proof
  that all heteroclinic connections must appear in adjacent layers — no other layer pairs need checking.
- **What it is good for:** finding *approximate* (near-) intersections of two given, already-computed
  torus manifolds' meshes at scale — i.e., a systematic, exhaustive geometric intersection search over
  the FULL discretized manifold surfaces, not a sampled grid search over a few free parameters (phase,
  time). §4.4 explicitly discusses (but does not execute) generalizing the method to a **spatial**
  (out-of-plane, 3D) CCR4BP — cuboid/tetrahedra meshes in a 6D stroboscopic phase space — this is about
  adding INCLINATION, not additional perturbing bodies (see §2 below for why this matters to `#724`).

## 2. Key result: Jupiter-Ganymede 4:3 <-> Jupiter-Europa 3:4

- Motivation for the resonance choice: mission design intent is Europa approach via the resonance
  sequence 3:4 -> 5:6 -> 1:1 (L2 Lyapunov) established by prior Anderson-group work [3,6]; this paper
  searches for what feeds INTO the 3:4 Jupiter-Europa resonance FROM Jupiter-Ganymede resonant orbits
  (a genuinely 4-body, different-moon connection — contrasted explicitly in their intro with all prior
  work, which only connected orbits resonant with the *same* moon).
- 3:2 Jupiter-Ganymede (naively expected useful, since Europa/Ganymede are near 2:1 Laplace-resonant)
  was RULED OUT: per their own 2021 paper, most of the unstable/useful part of that family collides with
  Europa's orbit and cannot be continued into CCR4BP tori.
- Tried **7:5 Jupiter-Ganymede** vs **3:4 Jupiter-Europa** first: only λ_u up to 1.4296 computable (slow
  manifold expansion); mesh search (globalized to layer 15) found 91 mesh intersections, but only 2 had
  exact-manifold-point velocity differences < 0.02 (nondimensional; 1 unit = 10,880 m/s, i.e. < 218 m/s),
  only 1 below 0.015 (< 163 m/s); minimum time-of-flight 24 Ganymede periods (~171 days).
- **Then tried 4:3 Jupiter-Ganymede (λ_u = 2.424) vs the same 3:4 Jupiter-Europa torus (λ_u = 8.56)**:
  227 mesh intersections (globalized only to layer 12), 45 with velocity diff < 0.02 units (< 218 m/s),
  **12 with velocity diff < 0.015 units (< 163 m/s)**. First mesh intersections at layer 9, TOF **18
  Ganymede periods ≈ 128 days**. This pair is the paper's headline "most promising candidate," explicitly
  flagged for further investigation.
- **Important negative caveat, directly usable for future positive-control framing:** even at the
  closest near-intersections found, **differential correction to an exact zero-ΔV manifold intersection
  FAILED** — the two manifolds are close but the semi-major-axis/eccentricity ("slow variable") behavior
  differs enough between them (JG 4:3 manifold varies little in a/e; JE 3:4 manifold varies a lot) that
  the authors conclude the manifolds are locally "nearly parallel," not truly intersecting nearby. So
  this is a genuine ΔV transfer opportunity (order 150-220 m/s scale, ~128-171 day TOF), NOT a ballistic
  heteroclinic connection — despite an exhaustive mesh search of the actual computed manifolds.
- §5.3 benchmarks against two prior, non-CCR4BP studies of the same Ganymede->Europa endgame problem:
  Anderson et al.'s Europa Lander optimization study [6] (ΔV "a little less than 150 m/s," TOF ~40 days)
  and Anderson 2021's patched-PCRTBP search [4] (ΔV down to ~55 m/s, TOF up to >200 days). This paper's
  own minimum ΔV is HIGHER than both and its minimum TOF is higher than [6]'s — the authors are explicit
  that no optimization was attempted here; the point was demonstrating the search methodology and
  identifying which resonance pair merits further work, not delivering an optimized transfer.
- §6 (a genuinely interesting side result, not central to our two cross-check questions): more-unstable
  4:3 Jupiter-Ganymede tori (λ_u > 2.458) could NOT be continued to the physical Europa mass in EITHER
  the stroboscopic-map or period-map formulation, both stopping at the same µ3 — ruled out as a numerical
  artifact via Calleja & de la Llave [11]'s Sobolev-seminorm divergence criterion (H^2/H^3 norms -> ∞ at
  the stopping µ3). Diagnosed as internal secondary-resonance gaps (Chirikov-standard-map-style, at
  rational rotation numbers 26π/40, 28π/43, cf. the sibling AAS-23-397 paper's related 11/34, 12/37
  findings, already in corpus) overlapping and destroying the more-unstable tori — a genuinely 4-body
  torus-breakdown mechanism, future work proposed (not yet executed) is computing the resonance-gap
  replacement orbits via Hamiltonian perturbation theory + multiple shooting.

## 3. Cross-check #1 — does this paper's method help `#715`'s stalled search?

**`#715`'s setup** (`src/cyclerfinder/search/ccr4bp_chained_transfer.py`): following Aryan & Fitzgerald
2024 (AAS 24-103, digested `#710`), a Callisto-L1 UNSTABLE manifold and a (separate CCR4BP system)
Europa-L2 STABLE manifold, each carrying Ganymede as perturber, chained via a Ganymede RENDEZVOUS —
i.e. both legs are checked for closest approach to Ganymede's own physical position, not for a mutual
manifold intersection in a shared phase space. Search was a coarse `(theta2, t)` phase/time GRID scan of
an already-globalized manifold tube; it converged the torus/manifold machinery cleanly at the paper's
own quoted Jacobi constants but found no close Ganymede encounter within its (modest) grid and
propagation window. `#715`'s own report explicitly flagged that "a finer/adaptive phase search or the
paper's own Poincaré-map methodology would be needed to chase it further."

**Answer: the method's core PRINCIPLE is exactly the class of upgrade `#715` called for, but it is not
a drop-in fix for `#715`'s specific problem geometry.**

- What transfers directly: replacing point/grid SAMPLING of a manifold with an exhaustive, spatially-
  partitioned MESH intersection test over the full discretized manifold surface is a strictly more
  systematic and complete search than a coarse `(theta2, t)` grid — it is precisely the "finer/adaptive"
  search #715 flagged as missing, and the spatial-partitioning trick (§4.2.1) is what makes it
  computationally tractable on ordinary hardware rather than requiring the finer grid #715 didn't have
  budget for. If `#715`'s two manifold tubes were meshed and mesh-intersected this way (rather than
  scanned for Ganymede proximity on a coarse grid), a genuine near-encounter — if one exists within the
  meshed range — would very likely be found where the grid scan missed it.
- What does NOT transfer directly: this paper's method finds intersections of two manifolds **in the
  SAME 4D CCR4BP phase space** (both W1^u, W2^s live in the same Jupiter-Ganymede-Europa system's
  coordinates — indeed §2.1's frame-transformation machinery exists precisely so that tori computed in
  one moon's synodic frame can be re-expressed in the other's for this comparison). `#715`'s problem is
  structurally different: it is a CROSS-SYSTEM rendezvous — an unstable manifold in the
  Callisto/Ganymede system reaching Ganymede's PHYSICAL position, separately from a stable manifold in
  the Europa/Ganymede system also reaching Ganymede's physical position — closer in spirit to the
  "patching two CRTBPs together" approach this paper's own intro (§1) explicitly contrasts itself
  against ("Past attempts... have involved various approximations, such as patched-conic models or
  patching two CRTBP models together. However... it would be more appropriate to use a restricted
  4-body model instead"). The paper never demonstrates or claims its mesh-intersection method applies to
  a physical-moon-rendezvous chain across two different CCR4BP instances; its own machinery is scoped to
  same-system torus-manifold pairs only.
- Practical implication: applying this paper's method to `#715` would require re-casting the problem —
  e.g., treating "closest approach to Ganymede" as a manifold-intersection-with-a-target-region test
  (mesh the manifold, mesh/represent Ganymede's own trajectory or a small tube around it, then apply the
  same spatial-partition + bounding-box + exact-test pipeline) rather than literal reuse of "intersect
  two resonant-orbit manifolds." This is a plausible, well-motivated NEXT step for `#715` (worth a future
  task), not a paper that already solves it.
- Secondary caveat worth carrying forward: even in this paper's own home turf (same-system torus-torus
  intersection), the mesh search found only NEAR-intersections requiring nonzero ΔV, not exact
  heteroclinic connections — for tori that even the authors judge among the most promising pairs
  available. This tempers expectations for what a mesh-based upgrade to `#715` would deliver: likely a
  ΔV transfer opportunity bound, not necessarily an exact ballistic connection, even if the search were
  made exhaustive.

## 4. Cross-check #2 — does this paper touch N=5 / Io in any way relevant to `#724`'s novelty claim?

**`#724`'s exact claim** (re-read in full from `2026-07-27-724-final-confirmation-n5-torus-novelty.md`):
first computed quasi-periodic invariant-torus substitute of the Kumar et al. 2021 Jupiter-Europa 3:4
resonant orbit (EXTERIOR to Europa), continued to the physical Io mass in a Laplace-locked N=5
(Jupiter-Io-Europa-Ganymede) restricted five-body model, with Ganymede at its PHYSICAL (non-idealized)
synodic rate. `#724`'s own §3.3 already anticipated this exact paper as a "flagged... acquisition
candidate" and pre-judged from its abstract that it is "still N=4 by its own title/abstract"; that prior
judgment is now checked against the FULL text, not just the abstract.

**Answer, from the complete full-text read: no. This paper's scope is entirely N=4 (Jupiter + Ganymede +
Europa + spacecraft). It never mentions Io, the Laplace resonance chain beyond Europa-Ganymede's own 2:1
(cited once, §1, purely to explain why 3:2 Jupiter-Ganymede might naively seem useful for Europa
transfers — Ref. [9], Barnes' "Laplace Resonance" encyclopedia entry, not a computation), or any
extension of the model to a third perturbing moon. `#724`'s claim is unaffected.**

Specific textual grounding (this determines the answer, not inference from silence):

- **§2, Eq. (2)-(3):** the model equations of motion are written for exactly three masses `m1, m2, m3`
  (Jupiter, Europa/Ganymede in either synodic frame) plus the massless spacecraft — a literal N=4
  system, matching Kumar et al. 2021's own model, not extended.
- **§2.2 ("A note on planar model validity"):** discusses only Europa's and Ganymede's (0.462°, 0.207°)
  inclinations relative to Jupiter's equator, justifying the PLANAR assumption — no mention of adding
  Io or any other body.
- **§4.4 ("Some Notes on Potential Generalizations to Higher Dimensional Models"):** the ONE place in
  the paper that discusses extending the model at all. It is explicit that the generalization direction
  contemplated is **spatial** (out-of-plane, 3D tori in a 7D extended phase space `R^6 x T`, i.e. adding
  z, pz, and out-of-plane inclination) — NOT additional perturbing bodies. Direct quote: "there may be
  certain situations where the use of a spatial model is desirable, such as for highly inclined,
  out-of-plane spacecraft orbits." Cuboid/tetrahedra mesh generalizations are proposed for THIS spatial
  extension only. No sentence anywhere proposes or anticipates adding a fourth perturbing moon (Io or
  otherwise) to reach an N=5 model.
- **§7 (Conclusions):** future work listed is (a) computing the resonance-gap-replacement periodic
  orbits/librational tori for the unstable 4:3 Jupiter-Ganymede family (§6's breakdown finding), and
  (b) finding TRUE (zero-ΔV) heteroclinic connections in the (still N=4) CCR4BP. Neither mentions Io,
  N=5, or any additional-moon extension.
- **References [1]-[32]:** none is an N=5 / multi-perturber paper (see §5 citation-mining pass below);
  the closest thing, ref. [9] (Barnes, "Laplace Resonance"), is cited only to note that Europa/Ganymede
  are near a 2:1 commensurability, which is what makes the (N=4) CCR4BP itself periodically forced in
  the first place — not a step toward N=5.

**Conclusion for `#724`:** this paper is confirmed, on a full independent read (not just the abstract),
to be entirely within the N=4 CCR4BP scope. It does not compute, sweep, or even speculatively discuss
extending the Jupiter-Europa 3:4 resonant torus (or any Kumar-class resonant torus) to additional
forcing from Io or any other third moon. `#724`'s "no N=5/multi-perturber extension attempted or
announced" finding for the Kumar lineage is independently reconfirmed by this pass; no correction to
`#724`'s novelty claim is warranted by this paper.

## 5. Mandatory citation-mining pass

Full introduction/background (§1) and reference list ([1]-[32]) read and cross-checked against
`docs/notes/CORPUS_INDEX.md` (2026-07-27 pass) and the `cyclers_pdf/papers/` filename listing directly
(broader net than the index, in case a file predates full indexing).

**Already in corpus** (no action needed):
- [22] Kumar, Anderson, de la Llave, Gunter 2021, AAS 21-651 — `kumar-anderson-delallave-gunter-2021-...`
  (digested `#688`).
- [20] Kumar, Anderson, de la Llave 2021 (AAS 21-349) — this is the companion GPU-connections paper;
  its JOURNAL version (SIADS 24:219, arXiv:2109.14814) is already in corpus and digested
  (`2026-07-03-...whiskered-tori-connections.md`); the AAS-21-349 conference version itself is not
  separately filed but the substantively-identical journal version is — no acquisition gap.
- [21] Kumar, Anderson, de la Llave 2022, CMDA 134(1):3 — the "rapid and accurate methods for computing
  whiskered tori" journal paper; this appears to be the SAME lineage as arXiv:2109.14814/SIADS 24:219
  above (need to double check if it's a distinct CMDA paper vs. the SIADS one — flagging as
  **possible-distinct-candidate**, medium priority: title differs ("Rapid and Accurate Methods for
  Computing Whiskered Tori..." CMDA 134(1):3, 2022) from the SIADS one ("Using GPUs and the
  Parameterization Method..." SIADS 24:219) enough that these may be two different papers by the same
  trio, only one of which (SIADS) is currently in corpus. Not acquired this pass — flagged only.
- [28] Russell & Strange 2009, JGCD 32(1):143-157 — `russell-strange-2009-cycler-trajectories-planetary-moon-systems-...` (digested, THE canonical moon-cycler census, already load-bearing in this
  project).
- [24] Möller 1997 "A fast triangle-triangle intersection test" and [14] Ericson 2005 "Real-time
  collision detection" and [17] Figueiredo et al. 2002 collision-detection survey — pure computer-
  graphics/CS technique references, NOT astrodynamics-domain-overlapping; correctly out of scope for
  this project's corpus, not flagged.
- [32] Weisstein, "Farey sequence" (MathWorld) — general-math reference tool, not flagged.
- [12],[13] Chirikov 1960/1971 — foundational resonance-overlap theory, already well-represented
  conceptually across this project's existing Chirikov-criterion usage (e.g. the sibling AAS-23-397
  digest); not flagged as a new acquisition need (general theory, not new domain content).
- [9] Barnes 2011, "Laplace Resonance" (Springer encyclopedia entry) — a definitional reference, not a
  research result; not flagged.
- [15] Evans et al. 2018, "Monte: The next generation of mission design and navigation software" — JPL
  software-tool description, not flagged (tool reference, not a transferable result/method for our own
  pipeline).
- [19] Hatfield & Rinderle 2001, CATO user's guide — JPL internal tool documentation, not flagged.

**Genuinely new candidates, topically overlapping this project's search domain — flagged, NOT
acquired** (per task instructions):

1. **Anderson, Campagnola, Koh, McElrath & Woollands 2021**, "Endgame design for Europa lander:
   Ganymede to Europa approach," *Journal of the Astronautical Sciences* 68(1):96-119 [ref 6]. **HIGH
   priority** — this is the Europa Lander study whose ΔV (~150 m/s) and TOF (~40 days) this very paper
   benchmarks itself against in §5.3; it is the closest thing to a sourced, published positive-control
   number for the Ganymede->Europa endgame transfer problem this project's CCR4BP work could aim to
   reproduce or compare against.
2. **Anderson, R.L. 2021**, "Tour design using resonant-orbit invariant manifolds in patched circular
   restricted three-body problems," JGCD 44(1):106-119 [ref 4]. Medium-high priority — the "patched
   PCRTBP" comparison point (ΔV ~55 m/s, TOF >200 days) also benchmarked in §5.3; directly relevant
   methodological alternative to CCR4BP tour design.
3. **Anderson, R.L. 2015**, "Approaching moons from resonance via invariant manifolds," JGCD
   38(6):1097-1109 [ref 3]. Medium priority — establishes the 3:4->5:6->1:1 Jupiter-Europa resonance
   sequence this paper's whole search is motivated by.
4. **Anderson, Campagnola & Lantoine 2016**, "Broad search for unstable resonant orbits in the planar
   circular restricted three-body problem," CMDA 124(2):177-199 [ref 7]. Medium priority — a systematic
   PCRTBP resonant-orbit census, directly relevant to this project's own resonant-orbit search
   machinery.
5. **Anderson & Lo 2010/2011** (two items) — "Dynamical systems analysis of planetary flybys and
   approach: Planar Europa orbiter," JGCD 33(6):1899-1912 [ref 8]; "A dynamical systems analysis of
   resonant flybys: Ballistic case," JAS 58 (2011) [ref 1]; "Flyby design using heteroclinic and
   homoclinic connections of unstable resonant orbits," Adv. Astronautical Sci. 140 (2011) [ref 2].
   Medium priority — foundational Anderson/Lo resonant-flyby lineage repeatedly cited across this whole
   Kumar/Jovian thread but not yet directly in corpus (only downstream Kumar/Baresi papers that build on
   it are).
6. **Blazevski & Ocampo 2012**, "Periodic orbits in the concentric circular restricted four-body problem
   and their invariant manifolds," Physica D: Nonlinear Phenomena 241(13):1158-1167 [ref 10]. **HIGH
   priority** — this is THE original CCR4BP model-definition paper; Kumar et al.'s own Eq. (2)-(3) cite
   it directly for the equations-of-motion derivation. Our own `src/cyclerfinder/core/ccr4bp.py` traces
   its model definition through the Kumar papers back to this one, which is not itself in corpus — a
   genuine grounding-chain gap.
7. **Calleja & de la Llave 2010**, "A numerically accessible criterion for the breakdown of
   quasi-periodic solutions and its rigorous justification," Nonlinearity 23(9):2029-2058 [ref 11].
   Medium priority — the Sobolev-seminorm torus-breakdown diagnostic actually used in this paper's own
   §6; methodologically reusable if this project ever needs to distinguish "numerical failure" from
   "genuine torus breakdown" in its own continuation work (a distinction `project_388_wall_energy_selective`-type situations have run into before).
8. **Haro, Canadell, Figueras, Luque & Mondelo 2016**, *The Parameterization Method for Invariant
   Manifolds: From Rigorous Results to Effective Computations*, Applied Mathematical Sciences vol. 195,
   Springer [ref 18]. Medium priority — the foundational parameterization-method reference the entire
   Kumar torus-computation lineage extends; a textbook, so may be harder to acquire than a paper.
9. **Fernández, Haro & Mondelo 2022**, "Flow map parameterization methods for invariant tori in
   quasi-periodic Hamiltonian systems" (arXiv preprint, cited as in-prep/2022) [ref 16]. Low-medium
   priority — an alternative torus-computation technique mentioned as a possible generalization route
   for spatial models (§4.4).
10. **Marchand, Howell & Wilson 2007**, "An improved corrections process for constrained trajectory
    design in the n-body problem," JSR 44(4):884-897 [ref 23]. Low-medium priority — one of three
    differential-correction/patching techniques (with COSMIC/Monte and CATO) mentioned as the route from
    a CCR4BP-model trajectory to a full-ephemeris one (§5.2 remark); relevant to any future full-eph
    validation of a CCR4BP-derived transfer.
11. **Olikara, Z.P. 2016**, "Computation of quasi-periodic tori and heteroclinic connections in
    astrodynamics using collocation techniques," PhD thesis, University of Colorado Boulder [ref 26].
    Medium priority — an alternative (collocation-based, vs. Kumar's parameterization-method) approach
    to the same quasi-periodic-torus computation problem; a useful independent cross-check method if
    this project ever wants a second torus-computation implementation to validate against.
12. **Parker & Anderson 2014**, *Low-Energy Lunar Trajectory Design*, JPL Deep Space Communications and
    Navigation Series vol. 12, Wiley [ref 27]. Low priority (textbook, Earth-Moon-focused, tangential to
    the Jovian thread) — mentioned only for its multiple-shooting periodic-orbit computation scheme.
13. **Sweetser, Maddock, Johannesen, Bell, Penzo, Wolf, Williams, Matousek & Weinstein 1997**,
    "Trajectory design for a Europa orbiter mission: A plethora of astrodynamic challenges" [ref 29].
    Low-medium priority — an early (pre-resonance-hopping-era) Europa mission design study, historical
    context for the patched-conic approaches this whole CCR4BP lineage explicitly improves on.
14. **Vaquero & Senent 2018**, "Poincaré: A multi-body, multi-system trajectory design tool," 7th
    ICATT, Oberpfaffenhofen [ref 30]. Low priority — a JPL tool-description paper; cited only as the
    inspiration (k-d trees) for this paper's own spatial-partitioning idea, not a result to reproduce.
15. **Vaquero Escribano, M. 2013**, "Spacecraft Transfer Trajectory Design Exploiting Resonant Orbits In
    Multi-Body Environments," PhD thesis, Purdue University [ref 31]. Medium priority — this is the
    Vaquero & Howell Saturn-Titan-Hyperion resonant-manifold work cited in this paper's own §1 intro as
    a parallel-system precedent; directly relevant to any future Saturnian-system extension of this
    project's CCR4BP work.

None of the above were acquired this pass (per task scope: flag only). The two HIGH-priority items
(Anderson/Campagnola/Koh/McElrath/Woollands 2021 Europa Lander endgame study, and Blazevski & Ocampo
2012 CCR4BP model-definition paper) are the strongest candidates for a follow-on acquisition task.

## Summary answers (for the dispatching session)

- **Does this paper's GPU method help `#715`'s stalled search?** Its PRINCIPLE — exhaustive, spatially-
  partitioned mesh-intersection search over full manifold surfaces instead of coarse phase-grid sampling
  — is exactly the class of upgrade `#715`'s own report called for, and would very likely outperform a
  coarse grid if applied. But it is not literally reusable: this paper solves same-CCR4BP-system,
  resonant-torus-vs-resonant-torus manifold intersections, not `#715`'s cross-system,
  manifold-vs-physical-moon-position rendezvous chain. Reformulating `#715`'s problem to fit this
  paper's mesh-intersection machinery is a plausible, well-motivated future task, not something this
  paper already delivers. Also worth carrying forward: even on its own best-case pair (JG 4:3 / JE 3:4),
  the paper's exhaustive search still found only near- (nonzero-ΔV), not exact, connections.
- **Does this paper touch N=5/Io in a way bearing on `#724`'s novelty claim?** No. Confirmed by full-text
  read (not just the abstract): the paper is entirely N=4 (Jupiter+Ganymede+Europa+spacecraft); its only
  discussed generalization direction is to a SPATIAL (3D, out-of-plane) N=4 model, not to additional
  perturbing bodies. No mention of Io, N=5, or any third-moon extension anywhere in the text, figures, or
  references. `#724`'s novelty claim is unaffected and requires no correction from this paper.
